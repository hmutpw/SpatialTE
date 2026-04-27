#' @title Dirichlet Prior Construction
#' @name prior
NULL


#' Build spatial neighbor graph from spot coordinates
#' @keywords internal
.build_neighbor_graph <- function(spot_coords, neighbor_radius = NULL) {
  spots <- spot_coords$spot_id
  x     <- spot_coords$x
  y     <- spot_coords$y
  n     <- length(spots)

  # Estimate radius from median nearest-neighbor distance
  dists <- as.matrix(dist(cbind(x, y)))
  diag(dists) <- Inf
  nn_dists <- apply(dists, 1, min)
  med_dist <- stats::median(nn_dists)
  radius   <- if (!is.null(neighbor_radius)) neighbor_radius else
              1.5 * med_dist

  message(sprintf("  Spatial: median NN dist=%.1f, neighbor radius=%.1f",
                  med_dist, radius))

  neighbors <- vector("list", n)
  dists2    <- as.matrix(dist(cbind(x, y)))
  for (i in seq_len(n)) {
    nbr_idx <- which(dists2[i,] > 0 & dists2[i,] <= radius)
    neighbors[[i]] <- list(
      idx   = nbr_idx,
      dist  = dists2[i, nbr_idx],
      ids   = spots[nbr_idx]
    )
  }
  names(neighbors) <- spots
  list(neighbors = neighbors, radius = radius, median_dist = med_dist)
}


#' Calibrate gamma_base using reliable transcripts
#'
#' Selects a calibration set of transcripts that are reliably quantified
#' by both LR and SR, then estimates gamma_base as the ratio of
#' LR TPM to SR unique-read-based TPM.
#'
#' @param lr_prior data.table from .load_lr_prior().
#' @param ec_counts SpotECCounts object.
#' @param eff_len Named numeric vector.
#' @return Numeric scalar: calibrated gamma_base.
#' @keywords internal
.calibrate_gamma <- function(lr_prior, ec_counts, eff_len) {
  message("  Calibrating gamma_base...")

  # SR unique counts per transcript (summed across all spots)
  uniq_mat    <- ec_counts$unique_counts
  sr_uniq_tot <- Matrix::colSums(uniq_mat)

  # SR total counts per transcript
  ec_idx <- ec_counts$ec_index
  tx_names <- ec_counts$transcript_names

  # Total SR reads per transcript across all spots
  # = sum over all ECs containing this transcript × ec_count
  # We compute this from the spot_ec_counts matrix
  # For each EC, distribute count to its transcripts
  spot_ec <- ec_counts$spot_ec_counts
  ec_tot  <- Matrix::colSums(spot_ec)  # total count per EC across all spots

  sr_total_tot <- numeric(length(tx_names))
  names(sr_total_tot) <- tx_names

  for (i in seq_len(nrow(ec_idx))) {
    k_idx <- as.integer(strsplit(ec_idx$ec_key[i], "|", fixed=TRUE)[[1]])
    valid <- k_idx >= 1L & k_idx <= length(tx_names)
    if (!any(valid)) next
    sr_total_tot[k_idx[valid]] <-
      sr_total_tot[k_idx[valid]] + ec_tot[i]
  }

  sr_uniq_frac <- ifelse(sr_total_tot > 0,
                          sr_uniq_tot / sr_total_tot, 0)

  # LR combined depth (LR + SR)
  lr_dt  <- lr_prior[match(tx_names, lr_prior$transcript_id)]
  lr_cnt <- lr_dt$em_count
  lr_cer <- lr_dt$certainty
  combined_depth <- lr_cnt + sr_total_tot

  # Calibration set criteria
  el_vals <- eff_len[tx_names]
  calib_mask <- !is.na(lr_cer)       &
                lr_cer > 0.8         &
                sr_uniq_frac > 0.5   &
                combined_depth > 50  &
                !is.na(el_vals)      &
                el_vals >= 500       &
                el_vals <= 2000

  n_calib <- sum(calib_mask, na.rm=TRUE)
  if (n_calib < 20L) {
    warning(sprintf(
      "Calibration set too small (%d transcripts, need ≥20). Using gamma_base=10.",
      n_calib))
    return(10.0)
  }

  # gamma_base = median(LR_count / SR_unique_count × correction)
  # Both normalised by eff_len for comparability
  lr_tpm_calib <- (lr_cnt[calib_mask] / pmax(el_vals[calib_mask], 50)) /
                  sum(lr_cnt / pmax(el_vals, 50), na.rm=TRUE) * 1e6
  sr_tpm_calib <- (sr_uniq_tot[calib_mask] / pmax(el_vals[calib_mask], 50)) /
                  sum(sr_uniq_tot / pmax(el_vals, 50), na.rm=TRUE) * 1e6

  ratio <- lr_tpm_calib / pmax(sr_tpm_calib, 1e-9)
  gamma_base <- stats::median(ratio, na.rm=TRUE)
  gamma_base <- pmax(pmin(gamma_base, 100.0), 1.0)  # clamp to [1, 100]

  message(sprintf("  gamma_base = %.2f (calibration set: %d transcripts)",
                  gamma_base, n_calib))
  gamma_base
}


#' Build Dirichlet prior alpha matrix
#'
#' Computes per-spot, per-transcript alpha values.
#' alpha_k_s = gamma(k,s) × pi_k_LR
#'
#' @param lr_prior data.table from .load_lr_prior().
#' @param ec_counts SpotECCounts.
#' @param eff_len Named numeric vector.
#' @param gamma_base Numeric.
#' @param neighbor_graph List or NULL (spatial mode).
#' @param alpha_local Numeric (spatial weight, 0-1).
#' @param certainty_threshold Numeric.
#' @param gamma_low_certainty_factor Numeric.
#' @return list:
#'   alpha_mat: dgCMatrix (transcripts × spots)
#'   gamma_mat: dgCMatrix (transcripts × spots)
#'   sr_reliability: named numeric vector (per transcript, global)
#' @keywords internal
.build_prior <- function(lr_prior, ec_counts, eff_len,
                          gamma_base,
                          neighbor_graph            = NULL,
                          alpha_local               = 0.7,
                          certainty_threshold       = 0.1,
                          gamma_low_certainty_factor = 1.0) {

  message("  Building Dirichlet prior...")
  tx_names  <- ec_counts$transcript_names
  spot_ids  <- ec_counts$spot_ids
  n_tx      <- length(tx_names)
  n_spots   <- length(spot_ids)

  # --- LR quantities per transcript ---
  lr_dt     <- lr_prior[match(tx_names, lr_prior$transcript_id)]
  pi_lr     <- lr_dt$pi_lr
  lr_rel    <- lr_dt$lr_reliability
  lr_cer    <- lr_dt$certainty
  lr_cnt    <- lr_dt$em_count
  names(pi_lr) <- names(lr_rel) <- names(lr_cer) <- names(lr_cnt) <- tx_names

  # --- SR reliability (global, all spots summed) ---
  uniq_mat   <- ec_counts$unique_counts    # spots × transcripts
  ec_idx     <- ec_counts$ec_index
  spot_ec    <- ec_counts$spot_ec_counts

  sr_uniq_tot <- as.numeric(Matrix::colSums(uniq_mat))
  names(sr_uniq_tot) <- tx_names

  # SR total per transcript
  ec_tot <- as.numeric(Matrix::colSums(spot_ec))
  sr_total_tot <- numeric(n_tx)
  for (i in seq_len(nrow(ec_idx))) {
    k_idx <- as.integer(strsplit(ec_idx$ec_key[i], "|", fixed=TRUE)[[1]])
    valid <- k_idx >= 1L & k_idx <= n_tx
    if (!any(valid)) next
    sr_total_tot[k_idx[valid]] <- sr_total_tot[k_idx[valid]] + ec_tot[i]
  }
  sr_rel_global <- ifelse(sr_total_tot > 0,
                           sr_uniq_tot / sr_total_tot, 0)
  names(sr_rel_global) <- tx_names

  # --- Depth factor ---
  combined_depth <- lr_cnt + sr_total_tot
  max_depth      <- max(combined_depth, na.rm=TRUE)
  depth_factor   <- log(combined_depth + 1) / log(max_depth + 1)
  depth_factor[is.na(depth_factor)] <- 0

  # --- Spot total UMI counts ---
  spot_n <- as.numeric(Matrix::rowSums(uniq_mat))
  names(spot_n) <- spot_ids
  n_median <- max(stats::quantile(spot_n[spot_n > 0], 0.5), 1)

  # --- Per-spot SR reliability (spatial mode: incorporate neighbors) ---
  # In SC mode: sr_rel_s = sr_rel_global (same for all spots)
  # In spatial mode: adjust using neighbor unique counts
  if (!is.null(neighbor_graph)) {
    sr_rel_mat <- .compute_spatial_sr_reliability(
      ec_counts, neighbor_graph, alpha_local, tx_names, spot_ids)
  } else {
    # All spots share global SR reliability
    sr_rel_mat <- NULL  # handled below
  }

  # --- Compute gamma and alpha ---
  # We build sparse matrices: only store non-trivial entries
  # gamma(k,s) = gamma_base × lr_rel(k) / sr_rel(k,s) × depth(k) × spot_factor(s)
  # For efficiency: compute gamma per transcript (global part) first
  # then scale by spot factor

  # LR/SR ratio (with low-certainty guard)
  lr_sr_ratio <- ifelse(
    lr_cer < certainty_threshold,
    1.0,  # guard: don't trust the ratio
    lr_rel / pmax(sr_rel_global, 1e-9)
  )
  lr_sr_ratio <- pmax(pmin(lr_sr_ratio, 100.0), 0.01)
  if (!is.null(gamma_low_certainty_factor) && gamma_low_certainty_factor != 1.0) {
    low_cert <- lr_cer < certainty_threshold
    lr_sr_ratio[low_cert] <- lr_sr_ratio[low_cert] * gamma_low_certainty_factor
  }

  # Global gamma per transcript (before spot factor)
  gamma_tx <- gamma_base * lr_sr_ratio * depth_factor
  gamma_tx[is.na(gamma_tx)] <- gamma_base * 0.1

  # Per-spot spot factor: 1 / (1 + N_s / N_median)
  spot_factor <- 1.0 / (1.0 + spot_n / n_median)

  # Build alpha matrix (transcripts × spots) as sparse
  # alpha_k_s = gamma(k,s) × pi_lr(k)
  # For most spots, pi_lr is non-zero for all transcripts
  # but alpha will be very small for absent transcripts
  # → Use a floor: only store alpha where pi_lr > 1e-9

  # Build as dense then convert (manageable for n_tx × n_spots if n_tx~100k)
  # For 100k × 10k = 1e9 elements: too large for dense
  # Strategy: compute alpha only for transcripts with pi_lr > threshold
  # and treat the rest as the floor value (0.01)

  pi_threshold <- 1e-7  # below this, use alpha_floor
  alpha_floor  <- 0.01

  active_tx <- which(pi_lr > pi_threshold)
  message(sprintf("  Active transcripts for prior: %d / %d",
                  length(active_tx), n_tx))

  # Build sparse alpha matrix
  # For each active transcript k, alpha_k_s = gamma_tx[k] × spot_factor[s] × pi_lr[k]
  # = outer product of gamma_tx[active] × pi_lr[active] with spot_factor
  i_idx <- rep(active_tx, each = n_spots)
  j_idx <- rep(seq_len(n_spots), times = length(active_tx))
  gamma_k <- gamma_tx[active_tx]
  pi_k    <- pi_lr[active_tx]

  alpha_vals <- rep(gamma_k * pi_k, each=n_spots) * rep(spot_factor, times=length(active_tx))
  alpha_vals <- pmax(alpha_vals, alpha_floor)

  alpha_mat <- Matrix::sparseMatrix(
    i    = i_idx,
    j    = j_idx,
    x    = alpha_vals,
    dims = c(n_tx, n_spots),
    dimnames = list(tx_names, spot_ids)
  )

  message(sprintf("  Prior built: gamma_base=%.2f, median alpha=%.4f",
                  gamma_base, stats::median(alpha_vals)))

  list(
    alpha_mat        = alpha_mat,
    gamma_tx         = gamma_tx,
    pi_lr            = pi_lr,
    lr_reliability   = setNames(lr_rel, tx_names),
    sr_rel_global    = sr_rel_global,
    depth_factor     = depth_factor,
    spot_factor      = spot_factor
  )
}


#' Compute per-spot SR reliability incorporating spatial neighbors
#' @keywords internal
.compute_spatial_sr_reliability <- function(ec_counts, neighbor_graph,
                                              alpha_local, tx_names, spot_ids) {
  uniq_mat   <- ec_counts$unique_counts   # spots × transcripts
  n_tx       <- length(tx_names)
  n_spots    <- length(spot_ids)
  spot_pos   <- setNames(seq_along(spot_ids), spot_ids)

  # Total counts per transcript per spot (for denominator)
  # We use the same approach as global but per-spot
  # For efficiency: only compute for spots that have neighbors
  result <- Matrix::Matrix(0, nrow=n_spots, ncol=n_tx,
                            dimnames=list(spot_ids, tx_names), sparse=TRUE)

  neighbors <- neighbor_graph$neighbors
  for (si in seq_along(spot_ids)) {
    sid <- spot_ids[si]
    nbr <- neighbors[[sid]]
    if (is.null(nbr) || length(nbr$idx) == 0) next

    own_uniq  <- as.numeric(uniq_mat[si, ])
    nbr_idx   <- nbr$idx
    nbr_dists <- nbr$dist
    wts       <- 1.0 / pmax(nbr_dists, 1e-9)
    wts       <- wts / sum(wts)

    nbr_uniq  <- Matrix::colSums(
      sweep(uniq_mat[nbr_idx, , drop=FALSE], 1, wts, "*"))

    combined_uniq <- alpha_local * own_uniq + (1 - alpha_local) * as.numeric(nbr_uniq)
    # Store as row in result
    nz <- which(combined_uniq > 0)
    if (length(nz) > 0)
      result[si, nz] <- combined_uniq[nz]
  }
  result
}
