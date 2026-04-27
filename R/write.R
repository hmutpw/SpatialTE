#' @title Output Functions and Downstream Compatibility
#' @name write_methods
NULL

# ============================================================
# Helper: extract matrices from SPE/SCE or S3
# ============================================================

# SpatialTEResult is always S3.
# It may contain a 'sce' slot (SingleCellExperiment) or raw assay lists.

.has_sce <- function(x) !is.null(x[["sce"]])

.get_assay <- function(x, name) {
  if (.has_sce(x)) SummarizedExperiment::assay(x$sce, name)
  else x[["assays"]][[name]]
}
.get_rowdata <- function(x) {
  if (.has_sce(x)) as.data.frame(SummarizedExperiment::rowData(x$sce))
  else x[["rowData"]]
}
.get_coldata <- function(x) {
  if (.has_sce(x)) as.data.frame(SummarizedExperiment::colData(x$sce))
  else x[["colData"]]
}
.get_metadata <- function(x) {
  if (.has_sce(x)) S4Vectors::metadata(x$sce)
  else x[["metadata"]]
}


# ============================================================
# write_spatialTE — master output function
# ============================================================

#' Write all SpatialTE outputs to a directory
#' @param result SpatialTEResult or SpatialTEDataset.
#' @param outdir Character. Output directory.
#' @param compress Logical. Use .tsv.gz (default FALSE).
#' @param write_sharing Logical. Write sharing table (default TRUE, needs store_ec=TRUE).
#' @param min_sharing Numeric. Min sharing fraction (default 0.5).
#' @param max_ec_size_sharing Integer. Max EC size for pairwise expansion (default 50).
#' @export
write_spatialTE <- function(result, outdir, compress=FALSE,
                             write_sharing=TRUE,
                             min_sharing=0.5, max_ec_size_sharing=50L) {
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  ext <- if (compress) ".tsv.gz" else ".tsv"
  write_counts(result,       file.path(outdir, paste0("counts",        ext)))
  write_count_matrix(result, file.path(outdir, paste0("count_matrix",  ext)))
  write_tpm_matrix(result,   file.path(outdir, paste0("tpm_matrix",    ext)))
  write_eff_len(result,      file.path(outdir, paste0("efflen_table",  ext)))
  write_qc(result,           file.path(outdir, paste0("qc_summary",    ext)))
  if (write_sharing) {
    ec_data <- .get_metadata(result)$ec_data
    if (!is.null(ec_data)) {
      write_sharing(result, file.path(outdir, paste0("sharing_table", ext)),
                    min_sharing=min_sharing,
                    max_ec_size_sharing=max_ec_size_sharing)
    } else {
      message("Sharing table skipped (store_ec=FALSE or EC data unavailable).")
    }
  }
  message("Outputs written to: ", normalizePath(outdir))
  invisible(result)
}


# ============================================================
# write_counts
# ============================================================

#' Write long-format count table
#' @param result SpatialTEResult.
#' @param file Character output path.
#' @export
write_counts <- function(result, file) UseMethod("write_counts")

#' @export
write_counts.SpatialTEResult <- function(result, file) {
  count_mat  <- .get_assay(result, "counts")
  unique_mat <- .get_assay(result, "unique")
  multi_mat  <- .get_assay(result, "multi")
  total_mat  <- .get_assay(result, "total")
  tpm_mat    <- .get_assay(result, "tpm")
  conf_mat   <- .get_assay(result, "confidence")
  rd         <- .get_rowdata(result)
  cd         <- .get_coldata(result)

  spot_ids <- rownames(cd)
  tx_ids   <- rownames(rd)

  # Convert to long format efficiently
  # Only output non-zero entries
  count_dt <- .sparse_to_long(count_mat,  "em_count")
  uniq_dt  <- .sparse_to_long(unique_mat, "unique_count")
  multi_dt <- .sparse_to_long(multi_mat,  "multi_count")
  total_dt <- .sparse_to_long(total_mat,  "total_count")
  tpm_dt   <- .sparse_to_long(tpm_mat,    "tpm")
  conf_dt  <- .sparse_to_long(conf_mat,   "confidence")

  # Join all on transcript_id × spot_id
  base <- count_dt
  for (sub in list(uniq_dt, multi_dt, total_dt, tpm_dt, conf_dt)) {
    base <- merge(base, sub, by=c("transcript_id","spot_id"), all.x=TRUE)
  }
  base[is.na(base)] <- 0

  # Add transcript metadata
  rd_sel <- rd[, intersect(c("gene_id","is_novel","tx_length","eff_len",
                               "low_eff_len","lr_certainty",
                               "bulk_consistency_ratio"), names(rd)), drop=FALSE]
  rd_sel$transcript_id <- rownames(rd_sel)
  base <- merge(base, rd_sel, by="transcript_id", all.x=TRUE)

  # Add sample_id
  sid <- .get_metadata(result)$sample_id %||% "sample"
  base$sample_id <- sid

  # Reorder columns
  col_order <- c("spot_id","transcript_id","sample_id","gene_id",
                  "is_novel","tx_length","em_count","unique_count",
                  "multi_count","total_count","eff_len","tpm",
                  "certainty","multimapping_rate","confidence",
                  "lr_certainty","low_eff_len","bulk_consistency_ratio")
  # certainty = unique / total
  base$certainty <- ifelse(base$total_count > 0,
                            base$unique_count / base$total_count, 0)
  base$multimapping_rate <- 1 - base$certainty
  present_cols <- col_order[col_order %in% names(base)]
  base <- base[, ..present_cols]

  .write_file(data.table::as.data.table(base), file)
  message("Counts written to: ", file, " (", nrow(base), " rows)")
  invisible(result)
}

#' Convert sparse matrix to long data.table
#' @keywords internal
.sparse_to_long <- function(mat, value_name) {
  idx <- Matrix::which(mat != 0, arr.ind=TRUE)
  if (nrow(idx) == 0) return(data.table::data.table(
    transcript_id=character(0), spot_id=character(0)))
  data.table::data.table(
    transcript_id = rownames(mat)[idx[,1]],
    spot_id       = colnames(mat)[idx[,2]],
    value         = mat[idx]
  ) |> data.table::setnames("value", value_name)
}


# ============================================================
# write_count_matrix / write_tpm_matrix
# ============================================================

#' Write wide-format count matrix (transcripts × spots)
#' @param result SpatialTEResult.
#' @param file Character output path.
#' @export
write_count_matrix <- function(result, file) UseMethod("write_count_matrix")

#' @export
write_count_matrix.SpatialTEResult <- function(result, file) {
  mat <- .get_assay(result, "counts")
  dt  <- data.table::as.data.table(as.matrix(mat), keep.rownames="transcript_id")
  .write_file(dt, file)
  message("Count matrix: ", nrow(dt), " × ", ncol(dt)-1, " → ", file)
  invisible(result)
}

#' Write wide-format TPM matrix
#' @param result SpatialTEResult.
#' @param file Character output path.
#' @export
write_tpm_matrix <- function(result, file) UseMethod("write_tpm_matrix")

#' @export
write_tpm_matrix.SpatialTEResult <- function(result, file) {
  mat <- .get_assay(result, "tpm")
  dt  <- data.table::as.data.table(as.matrix(mat), keep.rownames="transcript_id")
  .write_file(dt, file)
  message("TPM matrix: ", nrow(dt), " × ", ncol(dt)-1, " → ", file)
  invisible(result)
}


# ============================================================
# write_eff_len
# ============================================================

#' Write effective length table
#' @export
write_eff_len <- function(result, file) UseMethod("write_eff_len")

#' @export
write_eff_len.SpatialTEResult <- function(result, file) {
  rd <- .get_rowdata(result)
  dt <- data.table::data.table(
    transcript_id = rownames(rd),
    tx_length     = rd$tx_length,
    eff_len       = rd$eff_len,
    low_eff_len   = rd$low_eff_len
  )
  .write_file(dt, file)
  message("Eff len table: ", nrow(dt), " transcripts → ", file)
  invisible(result)
}


# ============================================================
# write_sharing (lazy, needs store_ec)
# ============================================================

#' Write transcript sharing table
#' @param result SpatialTEResult.
#' @param file Character output path.
#' @param min_sharing Numeric (default 0.5).
#' @param max_ec_size_sharing Integer (default 50).
#' @export
write_sharing <- function(result, file,
                           min_sharing=0.5,
                           max_ec_size_sharing=50L) {
  UseMethod("write_sharing")
}

#' @export
write_sharing.SpatialTEResult <- function(result, file,
                                           min_sharing=0.5,
                                           max_ec_size_sharing=50L) {
  ec_data <- .get_metadata(result)$ec_data
  if (is.null(ec_data)) {
    message("Sharing table unavailable (store_ec=FALSE).")
    return(invisible(result))
  }
  sh <- .compute_sharing(ec_data, min_sharing, max_ec_size_sharing)
  .write_file(sh, file)
  message("Sharing table: ", nrow(sh), " pairs → ", file)
  invisible(result)
}

#' Compute transcript sharing from EC data
#' @keywords internal
.compute_sharing <- function(ec_counts, min_sharing, max_ec_size) {
  ec_idx  <- ec_counts$ec_index
  spot_ec <- ec_counts$spot_ec_counts
  tx_names <- ec_counts$transcript_names
  ec_tot  <- as.numeric(Matrix::colSums(spot_ec))
  total_per_tx <- numeric(length(tx_names))
  names(total_per_tx) <- tx_names

  shared_map <- list()

  for (i in seq_len(nrow(ec_idx))) {
    if (ec_idx$ec_size[i] <= 1L) next
    k_idx <- as.integer(strsplit(ec_idx$ec_key[i],"|",fixed=TRUE)[[1]])
    valid <- k_idx >= 1L & k_idx <= length(tx_names)
    k_idx <- k_idx[valid]
    if (length(k_idx) < 2L) next
    cnt <- ec_tot[i]
    for (ki in k_idx) total_per_tx[ki] <- total_per_tx[ki] + cnt

    if (ec_idx$ec_size[i] > max_ec_size) next  # skip pairwise for large ECs

    pairs <- utils::combn(sort(k_idx), 2)
    for (ci in seq_len(ncol(pairs))) {
      pk <- paste0(pairs[1,ci], "_", pairs[2,ci])
      shared_map[[pk]] <- (shared_map[[pk]] %||% 0) + cnt
    }
  }

  if (length(shared_map) == 0)
    return(data.table::data.table(
      transcript_1=character(0), transcript_2=character(0),
      shared_reads=numeric(0), sharing_fraction=numeric(0),
      recommendation=character(0)))

  rows <- lapply(names(shared_map), function(pk) {
    parts <- strsplit(pk, "_")[[1]]
    k1 <- as.integer(parts[1]); k2 <- as.integer(parts[2])
    sr  <- shared_map[[pk]]
    mn  <- min(total_per_tx[k1], total_per_tx[k2])
    sf  <- if (mn > 0) sr/mn else 0
    if (sf < min_sharing) return(NULL)
    rec <- if (sf >= 0.8) "high_merge_recommended" else
           if (sf >= 0.5) "moderate_caution" else "low_independent_ok"
    list(transcript_1=tx_names[k1], transcript_2=tx_names[k2],
         shared_reads=sr, sharing_fraction=round(sf,4), recommendation=rec)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.table::data.table(
    transcript_1=character(0), transcript_2=character(0),
    shared_reads=numeric(0), sharing_fraction=numeric(0),
    recommendation=character(0)))
  dt <- data.table::rbindlist(rows)
  data.table::setorder(dt, -sharing_fraction)
  dt
}


# ============================================================
# write_qc
# ============================================================

#' Write QC summary
#' @export
write_qc <- function(result, file) UseMethod("write_qc")

#' @export
write_qc.SpatialTEResult <- function(result, file) {
  cd  <- .get_rowdata(result) |> invisible()
  cd2 <- .get_coldata(result)
  dt  <- data.table::as.data.table(cd2, keep.rownames="spot_id")
  # Add bulk consistency summary
  rd  <- .get_rowdata(result)
  if ("bulk_consistency_ratio" %in% names(rd)) {
    n_high <- sum(rd$bulk_consistency_ratio > 2, na.rm=TRUE)
    n_low  <- sum(rd$bulk_consistency_ratio < 0.5, na.rm=TRUE)
    message(sprintf("  Bulk consistency: %d high (>2×), %d low (<0.5×) transcripts",
                    n_high, n_low))
  }
  .write_file(dt, file)
  message("QC summary: ", nrow(dt), " spots → ", file)
  invisible(result)
}


# ============================================================
# Accessor methods
# ============================================================

#' @export
get_counts <- function(x, ...) UseMethod("get_counts")
#' @export
get_counts.SpatialTEResult <- function(x, ...) .get_assay(x, "counts")

#' @export
get_qc <- function(x, ...) UseMethod("get_qc")
#' @export
get_qc.SpatialTEResult <- function(x, ...) .get_coldata(x)

#' @export
get_ec <- function(x, ...) UseMethod("get_ec")
#' @export
get_ec.SpatialTEResult <- function(x, ...) .get_metadata(x)$ec_data

#' @export
get_sharing <- function(x, min_sharing=0.5, max_ec_size_sharing=50L, ...)
  UseMethod("get_sharing")
#' @export
get_sharing.SpatialTEResult <- function(x, min_sharing=0.5,
                                         max_ec_size_sharing=50L, ...) {
  ec_data <- get_ec(x)
  if (is.null(ec_data)) { message("No EC data."); return(NULL) }
  .compute_sharing(ec_data, min_sharing, max_ec_size_sharing)
}


# ============================================================
# Downstream compatibility
# ============================================================

#' @export
as_count_matrix <- function(x, ...) UseMethod("as_count_matrix")
#' @export
as_count_matrix.SpatialTEResult <- function(x, round_counts=TRUE, ...) {
  mat <- .get_assay(x, "counts")
  if (round_counts) round(mat) else mat
}

#' @export
as_tpm_matrix <- function(x, ...) UseMethod("as_tpm_matrix")
#' @export
as_tpm_matrix.SpatialTEResult <- function(x, ...) .get_assay(x, "tpm")

#' @export
as_sce <- function(x, ...) UseMethod("as_sce")
#' @export
as_sce.SpatialTEResult <- function(x, ...) {
  # If already wrapped around an SCE, return it
  if (.has_sce(x)) return(x$sce)
  if (!requireNamespace("SingleCellExperiment", quietly=TRUE))
    stop("Package 'SingleCellExperiment' required.")
  SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = as_count_matrix(x), tpm = as_tpm_matrix(x)),
    rowData = S4Vectors::DataFrame(.get_rowdata(x)),
    colData = S4Vectors::DataFrame(.get_coldata(x))
  )
}

#' @export
as_seurat <- function(x, ...) UseMethod("as_seurat")
#' @export
as_seurat.SpatialTEResult <- function(x, ...) {
  if (!requireNamespace("Seurat", quietly=TRUE))
    stop("Package 'Seurat' required.")
  mat <- as_count_matrix(x, round_counts=TRUE)
  Seurat::CreateSeuratObject(counts=mat, ...)
}


# ============================================================
# QC Plots
# ============================================================

#' Plot spatial QC metrics
#' @param result SpatialTEResult with spatial coordinates.
#' @param metric Column from colData to plot (default "unique_frac").
#' @export
plot_qc_spatial <- function(result, metric="unique_frac") {
  if (!requireNamespace("ggplot2", quietly=TRUE))
    stop("ggplot2 required for plotting.")
  cd <- .get_coldata(result)
  if (!all(c("x","y") %in% names(cd)))
    stop("No spatial coordinates found. Use plot_certainty() for non-spatial data.")
  if (!metric %in% names(cd))
    stop("Metric '", metric, "' not found in colData.")
  ggplot2::ggplot(cd, ggplot2::aes(x=x, y=y, colour=.data[[metric]])) +
    ggplot2::geom_point(size=0.8) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title=paste("Spatial QC:", metric),
                  colour=metric)
}

#' @export
plot_certainty <- function(x, ...) UseMethod("plot_certainty")
#' @export
plot_certainty.SpatialTEResult <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly=TRUE))
    stop("ggplot2 required.")
  rd <- .get_rowdata(x)
  ggplot2::ggplot(rd, ggplot2::aes(x=lr_certainty)) +
    ggplot2::geom_histogram(bins=50, fill="#1D9E75", colour="white") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="LR Certainty Distribution",
                  x="Certainty", y="Count")
}

#' @export
plot_ec_size <- function(x, ...) UseMethod("plot_ec_size")
#' @export
plot_ec_size.SpatialTEResult <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly=TRUE))
    stop("ggplot2 required.")
  ec_data <- get_ec(x)
  if (is.null(ec_data)) { message("No EC data."); return(invisible(NULL)) }
  ec_idx <- ec_data$ec_index
  ggplot2::ggplot(ec_idx, ggplot2::aes(x=pmin(ec_size,20L))) +
    ggplot2::geom_bar(fill="#7F77DD") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="EC Size Distribution",
                  x="EC size (capped at 20)", y="Number of ECs")
}
