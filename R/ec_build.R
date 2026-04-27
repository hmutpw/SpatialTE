#' @title Equivalence Class Construction
#' @description
#' Builds SpotECCounts from raw read-to-transcript assignments.
#' Uses integer encoding and data.table .GRP for fast EC key computation.
#' @name ec_build
NULL


#' Build equivalence classes from raw read assignments
#'
#' @param ec_raw data.table(spot_key, t_idx, nh) from .stream_bams().
#' @param transcript_names Character vector (integer-indexed).
#' @param sample_id_file Character or NULL. Maps spot barcode to sample_id.
#' @param sample_id Character. Used when sample_id_file is NULL.
#' @param ec_unit "barcode_umi" or "read".
#' @param max_ec_size Integer. ECs larger than this use LR prior directly.
#' @return Named list of SpotECCounts objects (one per sample).
#' @keywords internal
.build_ec <- function(ec_raw, transcript_names,
                       sample_id_file = NULL,
                       sample_id      = NULL,
                       ec_unit        = "barcode_umi",
                       max_ec_size    = 200L) {

  if (nrow(ec_raw) == 0) stop("No reads found. Check BAM and BC/UMI settings.")

  # --- Sample assignment ---
  if (!is.null(sample_id_file)) {
    sid_dt <- .read_file(sample_id_file, header=FALSE,
                          col.names=c("qname","sample_id"))
    # spot_key = barcode_umi; qname field maps back
    # For sample_id_file: match on barcode part of spot_key
    # spot_key = paste0(barcode, "_", umi), qname = original read_id
    # We stored spot_key directly; need barcode for join
    # Approach: store barcode separately in ec_raw
    # At this stage ec_raw only has spot_key — we need sample_id
    # The sample_id_file maps read_id → sample_id
    # But we discarded read_id earlier for memory efficiency
    # Solution: add barcode to ec_raw in stream_bams
    # For now: if sample_id_file provided, split spot_key to get barcode
    ec_raw[, barcode := sub("_[^_]+$", "", spot_key)]
    sid_by_bc <- unique(sid_dt[, .(barcode=qname, sample_id)])
    ec_raw <- merge(ec_raw, sid_by_bc, by="barcode", all.x=FALSE)
    samples <- unique(ec_raw$sample_id)
  } else {
    ec_raw[, sample_id := sample_id]
    samples <- sample_id
  }

  # --- Build EC per sample ---
  result <- lapply(samples, function(sid) {
    message(sprintf("  Building ECs for sample: %s", sid))
    sub <- ec_raw[sample_id == sid]
    .build_ec_one_sample(sub, transcript_names, sid, ec_unit, max_ec_size)
  })
  names(result) <- samples
  result
}


#' Build EC for a single sample
#' @keywords internal
.build_ec_one_sample <- function(dt, transcript_names, sample_id,
                                  ec_unit, max_ec_size) {
  if (ec_unit == "barcode_umi") {
    # Group by spot_key (= barcode_umi): collect all t_idx, form EC
    data.table::setorder(dt, spot_key, t_idx)
    read_ec <- dt[, .(
      ec_key   = paste(unique(sort(t_idx)), collapse="|"),
      spot_id  = sub("_[^_]+$", "", spot_key[1])  # barcode = spot_id
    ), by = spot_key]
  } else {
    # ec_unit == "read": each alignment record is a separate unit
    data.table::setorder(dt, spot_key, t_idx)
    read_ec <- dt[, .(
      ec_key  = paste(unique(sort(t_idx)), collapse="|"),
      spot_id = sub("_[^_]+$", "", spot_key[1])
    ), by = .(spot_key, t_idx)][
      , .(ec_key  = paste(unique(sort(t_idx)), collapse="|"),
          spot_id = sub("_[^_]+$", "", spot_key[1])),
      by = spot_key]
  }

  # --- EC index ---
  read_ec[, ec_id := .GRP, by = ec_key]
  ec_counts_dt <- read_ec[, .(ec_count = .N), by = .(spot_id, ec_id)]

  # EC metadata
  ec_meta <- unique(read_ec[, .(ec_id, ec_key)])
  ec_meta[, ec_size := sapply(strsplit(ec_key, "|", fixed=TRUE), length)]

  # EC type
  ec_meta[, ec_type := data.table::fcase(
    ec_size == 1L,                       "unique",
    ec_size > 1L & ec_size <= 5L,        "small_multi",
    ec_size > 5L & ec_size <= max_ec_size, "large_multi",
    ec_size > max_ec_size,               "overlimit"
  )]

  # Flag overlimit ECs (use prior directly in EM)
  overlimit_ecs <- ec_meta$ec_id[ec_meta$ec_type == "overlimit"]
  n_overlimit   <- length(overlimit_ecs)
  if (n_overlimit > 0) message(sprintf(
    "    %d ECs with size > %d will use LR prior directly.",
    n_overlimit, max_ec_size))

  # --- Unique counts (EC size == 1) per spot per transcript ---
  unique_ec_ids <- ec_meta$ec_id[ec_meta$ec_size == 1L]
  uniq_counts_dt <- ec_counts_dt[ec_id %in% unique_ec_ids]

  # Map ec_id → t_idx for unique ECs
  uniq_t <- ec_meta[ec_size == 1L, .(ec_id,
    t_idx = as.integer(ec_key))]
  uniq_counts_dt <- merge(uniq_counts_dt, uniq_t, by="ec_id")

  # --- Spot index ---
  all_spots <- sort(unique(ec_counts_dt$spot_id))
  n_spots   <- length(all_spots)
  n_tx      <- length(transcript_names)
  n_ecs     <- nrow(ec_meta)
  spot_pos  <- setNames(seq_along(all_spots), all_spots)
  ec_pos    <- setNames(seq_len(n_ecs), as.character(ec_meta$ec_id))

  # --- spot × EC sparse count matrix ---
  ec_counts_dt[, spot_i := spot_pos[spot_id]]
  ec_counts_dt[, ec_j   := ec_pos[as.character(ec_id)]]
  spot_ec_mat <- Matrix::sparseMatrix(
    i    = ec_counts_dt$spot_i,
    j    = ec_counts_dt$ec_j,
    x    = as.numeric(ec_counts_dt$ec_count),
    dims = c(n_spots, n_ecs),
    dimnames = list(all_spots, as.character(ec_meta$ec_id))
  )

  # --- spot × transcript unique count matrix ---
  if (nrow(uniq_counts_dt) > 0) {
    uniq_counts_dt[, spot_i := spot_pos[spot_id]]
    uniq_mat <- Matrix::sparseMatrix(
      i    = uniq_counts_dt$spot_i,
      j    = uniq_counts_dt$t_idx,
      x    = as.numeric(uniq_counts_dt$ec_count),
      dims = c(n_spots, n_tx),
      dimnames = list(all_spots, transcript_names)
    )
  } else {
    uniq_mat <- Matrix::sparseMatrix(
      i=integer(0), j=integer(0), x=numeric(0),
      dims=c(n_spots, n_tx),
      dimnames=list(all_spots, transcript_names))
  }

  message(sprintf(
    "    ECs: %d total | %d unique | %d small_multi | %d large_multi | %d overlimit",
    n_ecs,
    sum(ec_meta$ec_type=="unique"),
    sum(ec_meta$ec_type=="small_multi"),
    sum(ec_meta$ec_type=="large_multi"),
    n_overlimit
  ))
  message(sprintf("    Spots: %d | Unique assign rate: %.1f%%",
    n_spots,
    100 * sum(spot_ec_mat[, ec_meta$ec_id[ec_meta$ec_size==1L], drop=FALSE]) /
      max(sum(spot_ec_mat), 1)))

  structure(
    list(
      ec_index       = ec_meta,
      spot_ec_counts = spot_ec_mat,
      unique_counts  = uniq_mat,
      spot_ids       = all_spots,
      ec_ids         = as.character(ec_meta$ec_id),
      transcript_names = transcript_names,
      sample_id      = sample_id,
      overlimit_ecs  = overlimit_ecs
    ),
    class = "SpotECCounts"
  )
}


#' @export
print.SpotECCounts <- function(x, ...) {
  cat("SpotECCounts\n")
  cat("  Sample     :", x$sample_id, "\n")
  cat("  Spots      :", length(x$spot_ids), "\n")
  cat("  Transcripts:", length(x$transcript_names), "\n")
  cat("  ECs        :", nrow(x$ec_index), "\n")
  cat("    unique   :", sum(x$ec_index$ec_type=="unique"), "\n")
  cat("    overlimit:", length(x$overlimit_ecs), "\n")
  invisible(x)
}
