#' @title MAP-EM Quantification
#' @description
#' Runs MAP-EM for all spots in parallel using BiocParallel.
#' Optimizations:
#'   - Only processes active ECs (count > 0) per spot
#'   - Vectorized E-step using rowsum (no R for-loops)
#'   - Large ECs (overlimit) use LR prior directly, skip EM iteration
#'   - Sparse result collection: only stores non-zero counts
#' @name em_quantification
NULL


#' Run MAP-EM for all spots
#'
#' @param ec_counts SpotECCounts.
#' @param prior_result List from .build_prior().
#' @param eff_len Named numeric vector.
#' @param max_iter Integer.
#' @param tol Numeric.
#' @param n_cores Integer.
#' @param verbose Logical.
#' @return list:
#'   count_mat: dgCMatrix (transcripts × spots) — EM expected counts
#'   unique_mat: dgCMatrix (transcripts × spots) — unique counts
#'   multi_mat: dgCMatrix (transcripts × spots) — multi-mapping counts
#'   iter_info: data.table (spot_id, n_iter, converged)
#' @keywords internal
.run_em_all_spots <- function(ec_counts, prior_result, eff_len,
                               max_iter = 200L, tol = 1e-6,
                               n_cores = 1L, verbose = TRUE) {

  tx_names   <- ec_counts$transcript_names
  spot_ids   <- ec_counts$spot_ids
  n_tx       <- length(tx_names)
  n_spots    <- length(spot_ids)
  ec_idx     <- ec_counts$ec_index
  spot_ec    <- ec_counts$spot_ec_counts
  alpha_mat  <- prior_result$alpha_mat
  pi_lr      <- prior_result$pi_lr
  overlimit  <- ec_counts$overlimit_ecs

  # Pre-parse EC → transcript indices (integer vectors, fast lookup)
  ec_iso_idx <- lapply(ec_idx$ec_key, function(k) {
    as.integer(strsplit(k, "|", fixed=TRUE)[[1]])
  })
  names(ec_iso_idx) <- ec_idx$ec_id

  # eff_len aligned to transcript_names order
  eff_len_vec <- as.numeric(eff_len[tx_names])
  eff_len_vec[is.na(eff_len_vec)] <- 50.0
  eff_len_vec <- pmax(eff_len_vec, 50.0)

  # Overlimit EC set (use pi_lr directly)
  overlimit_set <- as.character(overlimit)

  if (verbose) message(sprintf(
    "  Running EM: %d spots × %d transcripts, %d cores",
    n_spots, n_tx, n_cores))

  # --- Choose parallel backend ---
  bp <- .make_bp_param(n_cores)

  # --- EM per spot ---
  spot_results <- BiocParallel::bplapply(
    seq_along(spot_ids),
    function(si) {
      sid     <- spot_ids[si]
      ec_row  <- spot_ec[si, ]
      active  <- which(ec_row > 0)
      if (length(active) == 0) {
        alpha_s <- as.numeric(alpha_mat[, si])
        rho     <- alpha_s / pmax(sum(alpha_s), 1e-300)
        return(list(sid=sid, t_idx=integer(0),
                    counts=numeric(0), n_iter=0L, converged=TRUE))
      }

      counts_active <- as.numeric(ec_row[active])
      ec_ids_active <- colnames(spot_ec)[active]

      # Alpha for this spot
      alpha_s <- as.numeric(alpha_mat[, si])
      alpha_s[is.na(alpha_s)] <- 0.01

      .em_one_spot(
        ec_ids_active   = ec_ids_active,
        counts_active   = counts_active,
        ec_iso_idx      = ec_iso_idx,
        n_tx            = n_tx,
        eff_len_vec     = eff_len_vec,
        alpha_s         = alpha_s,
        pi_lr           = pi_lr,
        overlimit_set   = overlimit_set,
        max_iter        = max_iter,
        tol             = tol,
        sid             = sid
      )
    },
    BPPARAM = bp
  )

  if (verbose) message("  EM complete. Assembling result matrices...")

  # --- Collect results into sparse matrices ---
  # Prepare index vectors for sparseMatrix construction
  all_t_i   <- integer(0)
  all_s_j   <- integer(0)
  all_count <- numeric(0)
  iter_rows <- vector("list", n_spots)

  for (si in seq_along(spot_results)) {
    res <- spot_results[[si]]
    if (length(res$t_idx) > 0) {
      all_t_i   <- c(all_t_i,   res$t_idx)
      all_s_j   <- c(all_s_j,   rep(si, length(res$t_idx)))
      all_count <- c(all_count, res$counts)
    }
    iter_rows[[si]] <- data.table::data.table(
      spot_id   = res$sid,
      n_iter    = res$n_iter,
      converged = res$converged
    )
  }

  count_mat <- Matrix::sparseMatrix(
    i    = all_t_i,
    j    = all_s_j,
    x    = all_count,
    dims = c(n_tx, n_spots),
    dimnames = list(tx_names, spot_ids)
  )

  iter_info <- data.table::rbindlist(iter_rows)

  # Check convergence
  n_not_conv <- sum(!iter_info$converged)
  if (n_not_conv > 0) {
    pct <- 100 * n_not_conv / n_spots
    if (pct > 10) warning(sprintf(
      "%.1f%% of spots (%d) did not converge. Consider increasing max_iter.",
      pct, n_not_conv))
    else message(sprintf("  %d spots (%.1f%%) did not converge.", n_not_conv, pct))
  }

  # Unique and multi count matrices (from ec_counts directly)
  unique_mat <- t(ec_counts$unique_counts)  # transcripts × spots

  # multi_count = total - unique
  # Total per transcript per spot from EC distribution (one E-step pass)
  total_mat <- .compute_total_counts(ec_counts)
  multi_mat <- total_mat - unique_mat
  multi_mat@x <- pmax(multi_mat@x, 0)

  list(
    count_mat  = count_mat,
    unique_mat = unique_mat,
    multi_mat  = multi_mat,
    total_mat  = total_mat,
    iter_info  = iter_info
  )
}


#' EM for a single spot — vectorized, no for-loops
#' @keywords internal
.em_one_spot <- function(ec_ids_active, counts_active, ec_iso_idx,
                           n_tx, eff_len_vec, alpha_s, pi_lr,
                           overlimit_set, max_iter, tol, sid) {

  n_active <- length(ec_ids_active)

  # Split into regular and overlimit ECs
  is_overlimit   <- ec_ids_active %in% overlimit_set
  reg_ids        <- ec_ids_active[!is_overlimit]
  reg_counts     <- counts_active[!is_overlimit]
  over_ids       <- ec_ids_active[is_overlimit]
  over_counts    <- counts_active[is_overlimit]

  # Pre-fetch isoform indices for regular ECs
  reg_iso_idx    <- lapply(reg_ids, function(eid) ec_iso_idx[[eid]])
  n_reg          <- length(reg_ids)

  # Build EC expansion vectors for vectorized E-step
  # ec_t_vec: transcript index for each (EC, transcript) pair
  # ec_e_vec: EC index for each (EC, transcript) pair
  # ec_cnt_vec: count for each EC
  if (n_reg > 0) {
    ec_sizes   <- sapply(reg_iso_idx, length)
    ec_t_vec   <- unlist(reg_iso_idx)
    # Guard: only keep indices within valid range
    valid_pos  <- ec_t_vec >= 1L & ec_t_vec <= n_tx
    ec_e_vec   <- rep(seq_len(n_reg), times = ec_sizes)[valid_pos]
    ec_cnt_vec <- rep(reg_counts,    times = ec_sizes)[valid_pos]
    ec_t_vec   <- ec_t_vec[valid_pos]
    t_idx_chr  <- as.character(seq_len(n_tx))
    n_reg      <- length(unique(ec_e_vec))
  }

  # Initialise rho from alpha
  alpha_safe <- pmax(alpha_s, 1e-9)
  rho <- alpha_safe / sum(alpha_safe)

  # Handle overlimit ECs: distribute using pi_lr directly
  # These counts are added once to the M-step delta (fixed, not iterated)
  delta_overlimit <- numeric(n_tx)
  if (length(over_ids) > 0) {
    for (oi in seq_along(over_ids)) {
      k_idx  <- ec_iso_idx[[over_ids[oi]]]
      valid  <- k_idx >= 1L & k_idx <= n_tx
      if (!any(valid)) next
      k_val  <- k_idx[valid]
      pi_sub <- pi_lr[k_val]
      pi_tot <- sum(pi_sub)
      if (pi_tot < 1e-300) next
      delta_overlimit[k_val] <-
        delta_overlimit[k_val] + over_counts[oi] * pi_sub / pi_tot
    }
  }

  # EM iterations (only for regular ECs)
  converged <- TRUE
  n_iter    <- 0L

  if (n_reg > 0) {
    converged <- FALSE
    for (iter in seq_len(max_iter)) {
      n_iter  <- iter
      rho_safe <- pmax(rho, 1e-300)

      # E-step (vectorized)
      unnorm   <- rho_safe[ec_t_vec] / eff_len_vec[ec_t_vec]
      ec_sums  <- rowsum(unnorm, ec_e_vec, reorder=FALSE)[, 1L] + 1e-300
      weights  <- unnorm / ec_sums[ec_e_vec]
      fracs    <- weights * ec_cnt_vec

      # M-step
      rs_tmp   <- rowsum(fracs, ec_t_vec, reorder=FALSE)
      delta_reg <- numeric(n_tx)
      hit       <- intersect(t_idx_chr, rownames(rs_tmp))
      if (length(hit)>0) delta_reg[as.integer(hit)] <- rs_tmp[hit, 1L]
      delta_reg[is.na(delta_reg)] <- 0

      delta     <- delta_reg + delta_overlimit
      rho_new   <- pmax(delta + alpha_s - 1.0, 0.0)
      rho_sum   <- sum(rho_new)
      if (rho_sum < 1e-300) break
      rho_new   <- rho_new / rho_sum

      # Convergence check
      diff <- sum(abs(rho_new - rho)) / sum(rho)
      rho  <- rho_new
      if (diff < tol) { converged <- TRUE; break }
    }
  } else if (length(over_ids) > 0) {
    # Only overlimit ECs: one M-step with overlimit delta
    rho_new <- pmax(delta_overlimit + alpha_s - 1.0, 0.0)
    rho_sum <- sum(rho_new)
    if (rho_sum > 1e-300) rho <- rho_new / rho_sum
    n_iter  <- 1L
  }

  # Final counts = last E-step
  if (n_reg > 0) {
    rho_safe <- pmax(rho, 1e-300)
    unnorm   <- rho_safe[ec_t_vec] / eff_len_vec[ec_t_vec]
    ec_sums  <- rowsum(unnorm, ec_e_vec, reorder=FALSE)[, 1L] + 1e-300
    weights  <- unnorm / ec_sums[ec_e_vec]
    fracs    <- weights * ec_cnt_vec
    rs_tmp2      <- rowsum(fracs, ec_t_vec, reorder=FALSE)
    final_counts <- numeric(n_tx)
    hit2         <- intersect(t_idx_chr, rownames(rs_tmp2))
    if (length(hit2)>0) final_counts[as.integer(hit2)] <- rs_tmp2[hit2, 1L]
    final_counts[is.na(final_counts)] <- 0
    final_counts <- final_counts + delta_overlimit
  } else {
    final_counts <- delta_overlimit
  }

  # Only return non-zero entries (sparse)
  nz <- which(final_counts > 1e-6)
  list(
    sid       = sid,
    t_idx     = nz,
    counts    = final_counts[nz],
    n_iter    = n_iter,
    converged = converged
  )
}


#' Compute total read count per transcript per spot from EC
#' (one E-step pass using current uniform weights)
#' @keywords internal
.compute_total_counts <- function(ec_counts) {
  ec_idx  <- ec_counts$ec_index
  spot_ec <- ec_counts$spot_ec_counts
  tx_names <- ec_counts$transcript_names
  spot_ids <- ec_counts$spot_ids
  n_tx    <- length(tx_names)
  n_spots <- length(spot_ids)

  # For each EC, distribute count equally among its transcripts
  i_idx <- integer(0); j_idx <- integer(0); vals <- numeric(0)
  ec_tot_per_spot <- spot_ec  # spots × ECs

  for (ei in seq_len(nrow(ec_idx))) {
    ec_id  <- ec_idx$ec_id[ei]
    k_idx  <- as.integer(strsplit(ec_idx$ec_key[ei], "|", fixed=TRUE)[[1]])
    valid  <- k_idx >= 1L & k_idx <= n_tx
    if (!any(valid)) next
    k_val  <- k_idx[valid]
    frac   <- 1.0 / length(k_val)
    # Per-spot counts for this EC
    ec_col <- as.numeric(spot_ec[, ei])
    nz_s   <- which(ec_col > 0)
    if (length(nz_s) == 0) next
    for (ki in k_val) {
      i_idx <- c(i_idx, rep(ki, length(nz_s)))
      j_idx <- c(j_idx, nz_s)
      vals  <- c(vals, ec_col[nz_s] * frac)
    }
  }

  if (length(i_idx) == 0) {
    return(Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0),
                                 dims=c(n_tx, n_spots),
                                 dimnames=list(tx_names, spot_ids)))
  }
  Matrix::sparseMatrix(
    i=i_idx, j=j_idx, x=vals,
    dims=c(n_tx, n_spots),
    dimnames=list(tx_names, spot_ids)
  )
}


#' Create BiocParallel backend
#' @keywords internal
.make_bp_param <- function(n_cores) {
  if (n_cores <= 1L) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = n_cores, type = "SOCK")
  } else {
    BiocParallel::MulticoreParam(workers = n_cores)
  }
}
