#' @title SpatialTE Result Objects
#' @name objects
NULL


#' Build SpatialTEResult from EM outputs
#'
#' Attempts to create a SpatialExperiment object (preferred),
#' falls back to SingleCellExperiment, then to a simple S3 list.
#'
#' @keywords internal
.build_result <- function(em_result, ec_counts, prior_result,
                           gtf_meta, eff_len, spot_coords,
                           lr_prior, gamma_base, metadata,
                           min_eff_len, n_threshold_pct = 10.0) {

  tx_names  <- ec_counts$transcript_names
  spot_ids  <- ec_counts$spot_ids
  n_tx      <- length(tx_names)
  n_spots   <- length(spot_ids)

  # --- counts and tpm matrices (transcripts × spots) ---
  count_mat  <- em_result$count_mat
  unique_mat <- em_result$unique_mat
  multi_mat  <- em_result$multi_mat
  total_mat  <- em_result$total_mat

  # TPM
  el_vec <- pmax(eff_len[tx_names], min_eff_len)
  rpk    <- count_mat / el_vec   # divide each row by its eff_len
  col_sums <- Matrix::colSums(rpk) + 1e-300
  tpm_mat  <- t(t(rpk) / col_sums) * 1e6

  # --- Confidence score ---
  # unique_frac = unique / total per (transcript, spot)
  total_safe <- total_mat
  total_safe@x <- pmax(total_safe@x, 1e-9)
  uniq_frac_mat <- unique_mat / total_safe

  # Spot N for confidence
  spot_n  <- Matrix::colSums(total_mat)
  n_thr   <- max(stats::quantile(spot_n[spot_n > 0], 0.10), 1)
  depth_factor <- pmin(spot_n / n_thr, 1.0)

  # lr_certainty per transcript
  lr_cer_vec <- lr_prior$certainty[match(tx_names, lr_prior$transcript_id)]
  lr_cer_vec[is.na(lr_cer_vec)] <- 0.0

  # confidence = uniq_frac × depth_factor × lr_certainty
  # depth_factor varies per spot, lr_certainty per transcript
  # Use sparse multiplication
  conf_mat <- uniq_frac_mat
  # Scale columns by depth_factor
  conf_mat <- t(t(conf_mat) * depth_factor)
  # Scale rows by lr_certainty
  conf_mat <- conf_mat * lr_cer_vec
  # Low eff_len penalty
  low_eff_idx <- which(el_vec <= min_eff_len * 1.1)
  if (length(low_eff_idx) > 0) conf_mat[low_eff_idx, ] <- 0

  # --- rowData (transcript metadata) ---
  row_meta <- data.frame(
    transcript_id        = tx_names,
    gene_id              = gtf_meta$gene_id[match(tx_names, gtf_meta$transcript_id)],
    is_novel             = gtf_meta$is_novel[match(tx_names, gtf_meta$transcript_id)],
    tx_length            = gtf_meta$tx_length[match(tx_names, gtf_meta$transcript_id)],
    eff_len              = round(el_vec, 2),
    low_eff_len          = el_vec <= min_eff_len * 1.1,
    lr_certainty         = round(lr_cer_vec, 4),
    lr_reliability       = round(prior_result$lr_reliability[
                                  match(tx_names, names(prior_result$lr_reliability))], 4),
    pi_lr                = round(prior_result$pi_lr[
                                  match(tx_names, names(prior_result$pi_lr))], 6),
    stringsAsFactors = FALSE,
    row.names = tx_names
  )

  # --- colData (spot metadata) ---
  total_reads   <- as.numeric(Matrix::colSums(total_mat))
  uniq_reads    <- as.numeric(Matrix::colSums(unique_mat))
  uniq_frac_sp  <- ifelse(total_reads > 0, uniq_reads / total_reads, 0)
  n_active_tx   <- as.numeric(Matrix::colSums(count_mat > 0.01))
  med_conf      <- apply(conf_mat, 2, function(x) {
    vals <- x[x > 0]; if (length(vals)==0) 0 else stats::median(vals)
  })

  col_meta <- data.frame(
    spot_id               = spot_ids,
    total_reads           = total_reads,
    unique_reads          = uniq_reads,
    unique_frac           = round(uniq_frac_sp, 4),
    n_active_transcripts  = n_active_tx,
    median_confidence     = round(med_conf, 4),
    n_iter                = em_result$iter_info$n_iter,
    converged             = em_result$iter_info$converged,
    stringsAsFactors = FALSE,
    row.names = spot_ids
  )

  # Add spatial coordinates if available
  if (!is.null(spot_coords)) {
    coord_dt <- spot_coords[match(spot_ids, spot_coords$spot_id)]
    col_meta$x <- coord_dt$x
    col_meta$y <- coord_dt$y
  }

  # --- Bulk consistency check ---
  bulk_check <- .compute_bulk_consistency(
    count_mat = count_mat,
    total_reads_per_spot = total_reads,
    pi_lr = prior_result$pi_lr,
    tx_names = tx_names
  )
  row_meta$bulk_consistency_ratio <- round(bulk_check[tx_names], 4)

  # --- Assemble result object ---
  mode <- if (!is.null(spot_coords)) "spatial" else "single_cell"

  assay_list <- list(
    counts     = count_mat,
    tpm        = tpm_mat,
    unique     = unique_mat,
    multi      = multi_mat,
    total      = total_mat,
    confidence = conf_mat
  )

  result <- .try_build_spe(
    assay_list  = assay_list,
    row_meta    = row_meta,
    col_meta    = col_meta,
    spot_coords = spot_coords,
    mode        = mode,
    metadata    = c(metadata, list(
      gamma_base = gamma_base,
      n_tx       = n_tx,
      n_spots    = n_spots,
      mode       = mode
    ))
  )
  result
}


#' Try to build SpatialExperiment → SCE → S3 fallback
#' @keywords internal
.try_build_spe <- function(assay_list, row_meta, col_meta,
                             spot_coords, mode, metadata) {

  # Option 1: SpatialExperiment
  if (requireNamespace("SpatialExperiment", quietly=TRUE) &&
      requireNamespace("SummarizedExperiment", quietly=TRUE) &&
      mode == "spatial" && !is.null(spot_coords)) {
    tryCatch({
      coord_mat <- as.matrix(col_meta[, c("x","y"), drop=FALSE])
      rownames(coord_mat) <- col_meta$spot_id
      spe <- SpatialExperiment::SpatialExperiment(
        assays           = assay_list,
        rowData          = S4Vectors::DataFrame(row_meta),
        colData          = S4Vectors::DataFrame(col_meta),
        spatialCoords    = coord_mat
      )
      metadata[["spatialte_result"]] <- TRUE
      S4Vectors::metadata(spe) <- metadata
      message("  Result: SpatialExperiment object created.")
      # Wrap in S3 envelope (same as SCE) to keep S3 dispatch working
      result <- structure(list(sce=spe, metadata=metadata,
                                n_tx=nrow(row_meta), n_spots=nrow(col_meta),
                                mode=mode),
                           class="SpatialTEResult")
      return(result)
    }, error = function(e) {
      message("  SpatialExperiment creation failed: ", e$message)
    })
  }

  # Option 2: SingleCellExperiment
  if (requireNamespace("SingleCellExperiment", quietly=TRUE) &&
      requireNamespace("SummarizedExperiment", quietly=TRUE)) {
    tryCatch({
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays  = assay_list,
        rowData = S4Vectors::DataFrame(row_meta),
        colData = S4Vectors::DataFrame(col_meta)
      )
      metadata[["spatialte_result"]] <- TRUE
      S4Vectors::metadata(sce) <- metadata
      message("  Result: SingleCellExperiment object created.")
      # Wrap in S3 envelope so SpatialTEResult S3 methods dispatch correctly
      result <- structure(list(sce=sce, metadata=metadata,
                                n_tx=nrow(row_meta), n_spots=nrow(col_meta),
                                mode=mode),
                           class="SpatialTEResult")
      return(result)
    }, error = function(e) {
      message("  SingleCellExperiment creation failed: ", e$message)
    })
  }

  # Option 3: Simple S3 list
  message("  Result: S3 list object (install SpatialExperiment or SingleCellExperiment for richer object).")
  structure(
    list(
      assays     = assay_list,
      rowData    = row_meta,
      colData    = col_meta,
      metadata   = metadata,
      n_tx       = nrow(row_meta),
      n_spots    = nrow(col_meta),
      mode       = mode
    ),
    class = "SpatialTEResult"
  )
}


#' @export
print.SpatialTEResult <- function(x, ...) {
  cat("SpatialTEResult\n")
  cat("  Mode       :", x[["mode"]] %||% "unknown", "\n")
  cat("  Transcripts:", x[["n_tx"]], "\n")
  cat("  Spots/cells:", x[["n_spots"]], "\n")
  cd <- SpatialTE:::.get_coldata(x)
  if (!is.null(cd) && "converged" %in% names(cd))
    cat("  Converged  :", sum(cd$converged), "/", nrow(cd), "spots\n")
  if (!is.null(cd) && "unique_frac" %in% names(cd))
    cat("  Median unique frac:", round(stats::median(cd$unique_frac,na.rm=TRUE), 3), "\n")
  invisible(x)
}


#' Compute bulk consistency ratio per transcript
#' @keywords internal
.compute_bulk_consistency <- function(count_mat, total_reads_per_spot,
                                       pi_lr, tx_names) {
  total_sum <- sum(total_reads_per_spot)
  if (total_sum < 1e-9) return(setNames(rep(NA_real_, length(tx_names)), tx_names))
  wt <- total_reads_per_spot / total_sum
  # rho_k_s ≈ count_k_s / total_s (approximate rho)
  total_per_spot <- pmax(total_reads_per_spot, 1e-9)
  rho_mat <- t(t(count_mat) / total_per_spot)
  bulk_check <- as.numeric(rho_mat %*% wt)
  names(bulk_check) <- tx_names
  pi_lr_aligned <- pi_lr[match(tx_names, names(pi_lr))]
  pi_lr_aligned[is.na(pi_lr_aligned)] <- 1e-9
  ratio <- bulk_check / pmax(pi_lr_aligned, 1e-9)
  setNames(ratio, tx_names)
}

#' @export
summary.SpatialTEResult <- function(object, ...) {
  print(object)
  rd <- SpatialTE:::.get_rowdata(object)
  cd <- SpatialTE:::.get_coldata(object)
  if (is.null(rd) || is.null(cd)) return(invisible(object))
  cat("\nTranscript summary:\n")
  cat("  Novel  :", sum(rd$is_novel, na.rm=TRUE), "\n")
  cat("  Low eff_len:", sum(rd$low_eff_len, na.rm=TRUE), "\n")
  cat("\nSpot summary:\n")
  cat("  Median total reads    :", round(median(cd$total_reads), 0), "\n")
  cat("  Median unique frac    :", round(median(cd$unique_frac), 3), "\n")
  cat("  Median active tx/spot :", round(median(cd$n_active_transcripts), 0), "\n")
  invisible(object)
}
