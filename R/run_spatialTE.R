#' @title Main Entry Points: run_spatialTE and run_scTE
#' @name run_spatialTE
NULL


#' Isoform quantification for spatial transcriptomics data
#'
#' Runs the full SpatialTE pipeline:
#' BAM streaming → effective length estimation → EC construction →
#' Dirichlet prior → MAP-EM → SpatialTEResult.
#'
#' @param pe_bam Character or NULL. PE toTranscriptome BAM path.
#' @param se_bam Character or NULL. SE toTranscriptome BAM path.
#'   At least one of pe_bam / se_bam must be provided.
#' @param gtf_file Character. IsoQuant transcript_models.gtf[.gz].
#' @param lr_input IsoEMResult, data.frame, or file path.
#'   Must contain columns: transcript_id, em_count.
#'   Optional column: certainty.
#' @param spot_coords Character path or data.frame.
#'   Spot coordinates for spatial prior smoothing.
#'   Formats: Visium CSV, Visium HD parquet/CSV, Patho-DBiT,
#'   or generic (spot_id | x | y).
#' @param coords_format Character. "auto", "visium", "visium_hd",
#'   "patho_dbit", or "generic".
#' @param sample_id Character or NULL. Sample name (single-sample mode).
#' @param sample_id_file Character or NULL. Two-column file (read_id | sample_id).
#'   For multi-sample joint STAR runs. Mutually exclusive with sample_id.
#' @param bc_umi_table Character path or data.frame or NULL.
#'   Three columns: read_id | barcode | umi.
#' @param bc_umi_pattern Character regex or NULL.
#'   Pattern to extract barcode and UMI from read ID.
#' @param cb_tag Character. BAM cell barcode tag (auto-detected if NULL).
#' @param ub_tag Character. BAM UMI tag (auto-detected if NULL).
#' @param cb_whitelist Character vector or file path or NULL.
#'   Valid barcodes filter.
#' @param ec_unit Character. "barcode_umi" (default) or "read".
#' @param max_ec_size Integer. ECs larger than this use LR prior directly (default 200).
#' @param frag_dist_source Character. "pe" (default), "se", or "both".
#'   Source for fragment length distribution estimation.
#' @param min_eff_len Numeric. Floor for effective length (default 50).
#' @param gamma_base Numeric or NULL. Prior strength. NULL = auto-calibrate.
#' @param lambda Numeric. Spot-specific prior smoothing (default 10).
#' @param certainty_threshold Numeric. LR/SR ratio guard threshold (default 0.1).
#' @param gamma_low_certainty_factor Numeric. Gamma scaling for low-certainty
#'   transcripts (default 1.0 = no scaling).
#' @param alpha_local Numeric. Weight for own vs neighbour unique reads (default 0.7).
#' @param neighbor_radius Numeric or NULL. Spatial neighbour radius in coord units.
#'   NULL = auto (1.5 × median nearest-neighbour distance).
#' @param max_iter Integer. EM maximum iterations (default 200).
#' @param tol Numeric. EM convergence tolerance (default 1e-6).
#' @param min_mapq Integer. Minimum MAPQ filter (default 0).
#' @param enforce_bulk_consistency Logical. Apply post-EM bulk correction (default FALSE).
#' @param store_ec Logical. Store EC data for write_sharing() (default TRUE).
#' @param n_cores Integer. Parallel workers (default 1).
#' @param verbose Logical (default TRUE).
#' @return SpatialTEResult (SpatialExperiment, SingleCellExperiment, or S3 list).
#' @export
#' @examples
#' \dontrun{
#' result <- run_spatialTE(
#'   pe_bam      = "spatial_pe.toTranscriptome.bam",
#'   gtf_file    = "transcript_models.gtf",
#'   lr_input    = "isoquant_counts.tsv.gz",
#'   spot_coords = "tissue_positions.csv",
#'   sample_id   = "sample_A"
#' )
#' print(result)
#' write_spatialTE(result, outdir = "results/")
#' }
run_spatialTE <- function(
    pe_bam                    = NULL,
    se_bam                    = NULL,
    gtf_file,
    lr_input,
    spot_coords               = NULL,
    coords_format             = "auto",
    sample_id                 = NULL,
    sample_id_file            = NULL,
    bc_umi_table              = NULL,
    bc_umi_pattern            = NULL,
    cb_tag                    = "CB",
    ub_tag                    = "UB",
    cb_whitelist              = NULL,
    ec_unit                   = "barcode_umi",
    max_ec_size               = 200L,
    frag_dist_source          = "pe",
    min_eff_len               = 50,
    gamma_base                = NULL,
    lambda                    = 10,
    certainty_threshold       = 0.1,
    gamma_low_certainty_factor = 1.0,
    alpha_local               = 0.7,
    neighbor_radius           = NULL,
    max_iter                  = 200L,
    tol                       = 1e-6,
    min_mapq                  = 0L,
    enforce_bulk_consistency  = FALSE,
    store_ec                  = TRUE,
    n_cores                   = 1L,
    verbose                   = TRUE
) {
  .run_core(
    pe_bam=pe_bam, se_bam=se_bam, gtf_file=gtf_file,
    lr_input=lr_input, spot_coords=spot_coords,
    coords_format=coords_format,
    sample_id=sample_id, sample_id_file=sample_id_file,
    bc_umi_table=bc_umi_table, bc_umi_pattern=bc_umi_pattern,
    cb_tag=cb_tag, ub_tag=ub_tag, cb_whitelist=cb_whitelist,
    ec_unit=ec_unit, max_ec_size=max_ec_size,
    frag_dist_source=frag_dist_source, min_eff_len=min_eff_len,
    gamma_base=gamma_base, lambda=lambda,
    certainty_threshold=certainty_threshold,
    gamma_low_certainty_factor=gamma_low_certainty_factor,
    alpha_local=alpha_local, neighbor_radius=neighbor_radius,
    max_iter=max_iter, tol=tol, min_mapq=min_mapq,
    enforce_bulk_consistency=enforce_bulk_consistency,
    store_ec=store_ec, n_cores=n_cores, verbose=verbose,
    use_spatial=TRUE
  )
}


#' Isoform quantification for single-cell RNA-seq data
#'
#' Same as \code{run_spatialTE()} but without spatial coordinate support.
#' No spatial prior smoothing is applied.
#'
#' @inheritParams run_spatialTE
#' @return SpatialTEResult (SingleCellExperiment or S3 list).
#' @export
#' @examples
#' \dontrun{
#' result <- run_scTE(
#'   pe_bam   = "sc_pe.toTranscriptome.bam",
#'   gtf_file = "transcript_models.gtf",
#'   lr_input = "isoquant_counts.tsv.gz",
#'   sample_id = "sample_A"
#' )
#' }
run_scTE <- function(
    pe_bam                    = NULL,
    se_bam                    = NULL,
    gtf_file,
    lr_input,
    sample_id                 = NULL,
    sample_id_file            = NULL,
    bc_umi_table              = NULL,
    bc_umi_pattern            = NULL,
    cb_tag                    = "CB",
    ub_tag                    = "UB",
    cb_whitelist              = NULL,
    ec_unit                   = "barcode_umi",
    max_ec_size               = 200L,
    frag_dist_source          = "pe",
    min_eff_len               = 50,
    gamma_base                = NULL,
    lambda                    = 10,
    certainty_threshold       = 0.1,
    gamma_low_certainty_factor = 1.0,
    max_iter                  = 200L,
    tol                       = 1e-6,
    min_mapq                  = 0L,
    enforce_bulk_consistency  = FALSE,
    store_ec                  = TRUE,
    n_cores                   = 1L,
    verbose                   = TRUE
) {
  .run_core(
    pe_bam=pe_bam, se_bam=se_bam, gtf_file=gtf_file,
    lr_input=lr_input, spot_coords=NULL, coords_format="auto",
    sample_id=sample_id, sample_id_file=sample_id_file,
    bc_umi_table=bc_umi_table, bc_umi_pattern=bc_umi_pattern,
    cb_tag=cb_tag, ub_tag=ub_tag, cb_whitelist=cb_whitelist,
    ec_unit=ec_unit, max_ec_size=max_ec_size,
    frag_dist_source=frag_dist_source, min_eff_len=min_eff_len,
    gamma_base=gamma_base, lambda=lambda,
    certainty_threshold=certainty_threshold,
    gamma_low_certainty_factor=gamma_low_certainty_factor,
    alpha_local=0.7, neighbor_radius=NULL,
    max_iter=max_iter, tol=tol, min_mapq=min_mapq,
    enforce_bulk_consistency=enforce_bulk_consistency,
    store_ec=store_ec, n_cores=n_cores, verbose=verbose,
    use_spatial=FALSE
  )
}


# ============================================================
# Internal core pipeline
# ============================================================

#' @keywords internal
.run_core <- function(pe_bam, se_bam, gtf_file, lr_input,
                       spot_coords, coords_format,
                       sample_id, sample_id_file,
                       bc_umi_table, bc_umi_pattern, cb_tag, ub_tag,
                       cb_whitelist, ec_unit, max_ec_size,
                       frag_dist_source, min_eff_len,
                       gamma_base, lambda,
                       certainty_threshold, gamma_low_certainty_factor,
                       alpha_local, neighbor_radius,
                       max_iter, tol, min_mapq,
                       enforce_bulk_consistency, store_ec,
                       n_cores, verbose, use_spatial) {

  if (verbose) message("=== SpatialTE pipeline ===")
  t0 <- proc.time()

  # --- Validate ---
  .validate_inputs(pe_bam, se_bam, gtf_file, lr_input,
                    sample_id, sample_id_file)

  # --- [1] Parse GTF (once) ---
  if (verbose) message("[1/6] Parsing GTF...")
  gtf_meta <- .parse_gtf(gtf_file)
  tx_names_init <- gtf_meta$transcript_id

  # --- [2] Load LR prior ---
  if (verbose) message("[2/6] Loading long-read prior...")
  lr_prior <- .load_lr_prior(lr_input, gtf_meta)

  # --- [3] Load coordinates (spatial mode) ---
  neighbor_graph <- NULL
  spot_coords_dt <- NULL
  if (use_spatial && !is.null(spot_coords)) {
    if (verbose) message("[3/6] Loading spatial coordinates...")
    spot_coords_dt <- .load_coords(spot_coords, coords_format)
    neighbor_graph <- .build_neighbor_graph(spot_coords_dt, neighbor_radius)
  } else {
    if (verbose) message("[3/6] No coordinates — single-cell mode.")
  }

  # --- [4] Stream BAMs (once): g(x), f(d), EC raw ---
  if (verbose) message("[4/6] Streaming BAM(s)...")
  stream_res <- .stream_bams(
    pe_bam           = pe_bam,
    se_bam           = se_bam,
    transcript_names = tx_names_init,
    gtf_meta         = gtf_meta,
    bc_umi_table     = bc_umi_table,
    bc_umi_pattern   = bc_umi_pattern,
    cb_tag           = cb_tag,
    ub_tag           = ub_tag,
    cb_whitelist     = cb_whitelist,
    frag_dist_source = frag_dist_source,
    min_mapq         = min_mapq
  )

  # Compute eff_len
  eff_len <- .compute_eff_len(
    gtf_meta     = gtf_meta,
    g_cov        = stream_res$g_cov,
    frag_len_tab = stream_res$frag_len_tab,
    min_eff_len  = min_eff_len
  )

  # Build ECs
  ec_list <- .build_ec(
    ec_raw           = stream_res$ec_raw,
    transcript_names = stream_res$transcript_names,
    sample_id_file   = sample_id_file,
    sample_id        = sample_id,
    ec_unit          = ec_unit,
    max_ec_size      = max_ec_size
  )

  # --- [5] Prior + EM per sample ---
  if (verbose) message("[5/6] Building prior and running EM...")
  sample_results <- lapply(names(ec_list), function(sid) {
    if (verbose) message(sprintf("  Sample: %s", sid))
    ec_counts <- ec_list[[sid]]

    # Calibrate gamma
    gbase <- if (is.null(gamma_base)) {
      .calibrate_gamma(lr_prior, ec_counts, eff_len)
    } else gamma_base

    # Build prior
    prior_res <- .build_prior(
      lr_prior                  = lr_prior,
      ec_counts                 = ec_counts,
      eff_len                   = eff_len,
      gamma_base                = gbase,
      neighbor_graph            = if (use_spatial) neighbor_graph else NULL,
      alpha_local               = alpha_local,
      certainty_threshold       = certainty_threshold,
      gamma_low_certainty_factor = gamma_low_certainty_factor
    )

    # EM
    em_res <- .run_em_all_spots(
      ec_counts   = ec_counts,
      prior_result = prior_res,
      eff_len     = eff_len,
      max_iter    = max_iter,
      tol         = tol,
      n_cores     = n_cores,
      verbose     = verbose
    )

    # Optional: enforce bulk consistency
    if (enforce_bulk_consistency) {
      em_res <- .apply_bulk_correction(em_res, prior_res$pi_lr,
                                        ec_counts$transcript_names)
    }

    # Build result object
    .build_result(
      em_result    = em_res,
      ec_counts    = ec_counts,
      prior_result = prior_res,
      gtf_meta     = gtf_meta,
      eff_len      = eff_len,
      spot_coords  = spot_coords_dt,
      lr_prior     = lr_prior,
      gamma_base   = gbase,
      metadata     = list(
        sample_id        = sid,
        gtf_file         = gtf_file,
        ec_unit          = ec_unit,
        max_ec_size      = max_ec_size,
        frag_dist_source = frag_dist_source,
        max_iter         = max_iter,
        tol              = tol,
        enforce_bulk_consistency = enforce_bulk_consistency,
        ec_data          = if (store_ec) ec_counts else NULL
      ),
      min_eff_len = min_eff_len
    )
  })
  names(sample_results) <- names(ec_list)

  elapsed <- (proc.time() - t0)[["elapsed"]]
  if (verbose) message(sprintf("[6/6] Done. Total time: %.1f seconds", elapsed))

  if (length(sample_results) == 1L) {
    return(sample_results[[1]])
  }

  # Multi-sample: return list
  structure(
    list(
      samples     = sample_results,
      sample_meta = data.frame(
        sample_id = names(sample_results),
        stringsAsFactors = FALSE
      )
    ),
    class = "SpatialTEDataset"
  )
}


#' Apply post-EM bulk consistency correction
#' @keywords internal
.apply_bulk_correction <- function(em_res, pi_lr, tx_names) {
  count_mat    <- em_res$count_mat
  total_per_sp <- as.numeric(Matrix::colSums(em_res$total_mat))
  total_sum    <- sum(total_per_sp)
  if (total_sum < 1e-9) return(em_res)
  wt <- total_per_sp / total_sum
  bulk_check <- as.numeric(count_mat %*% wt) /
                pmax(as.numeric(Matrix::rowSums(count_mat)) / total_sum, 1e-9)
  pi_aligned <- pi_lr[match(tx_names, names(pi_lr))]
  pi_aligned[is.na(pi_aligned)] <- 1e-9
  correction <- pi_aligned / pmax(bulk_check, 1e-9)
  correction <- pmax(pmin(correction, 10), 0.1)
  em_res$count_mat <- count_mat * correction
  em_res
}


#' @export
print.SpatialTEDataset <- function(x, ...) {
  cat("SpatialTEDataset\n")
  cat("  Samples:", length(x$samples), "\n")
  cat("  IDs    :", paste(names(x$samples), collapse=", "), "\n")
  invisible(x)
}
