#' @title Internal Utility Functions
#' @name utils
NULL

# ============================================================
# File I/O
# ============================================================

#' Read delimited file with .gz support
#' @keywords internal
.read_file <- function(path, ...) {
  if (!file.exists(path)) stop("File not found: ", path)
  data.table::fread(path, ...)
}

#' Write data.table with optional .gz compression
#' @keywords internal
.write_file <- function(dt, path, ...) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(
    dt, path, sep = "\t", quote = FALSE,
    compress = if (grepl("\\.gz$", path)) "gzip" else "none",
    ...
  )
}

#' Check if file is gzipped
#' @keywords internal
.is_gz <- function(path) grepl("\\.gz$", path, ignore.case = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a


# ============================================================
# GTF parsing
# ============================================================

#' Parse IsoQuant transcript_models.gtf
#'
#' Extracts per-transcript metadata: gene_id, is_novel, tx_length
#' (sum of exon widths, NOT genomic span).
#' Uses rtracklayer if available, otherwise falls back to base R.
#'
#' Novel detection rule: transcript IDs not starting with a known
#' reference prefix (ENS, NM_, NR_, XM_, XR_) are novel.
#'
#' @param gtf_path Character. Path to GTF file (supports .gz).
#' @return data.table: transcript_id, gene_id, is_novel, tx_length.
#' @keywords internal
.parse_gtf <- function(gtf_path) {
  if (!file.exists(gtf_path)) stop("GTF file not found: ", gtf_path)
  message("  Parsing GTF: ", basename(gtf_path))
  if (requireNamespace("rtracklayer", quietly = TRUE))
    return(.parse_gtf_rtracklayer(gtf_path))
  .parse_gtf_base(gtf_path)
}

#' @keywords internal
.parse_gtf_rtracklayer <- function(gtf_path) {
  gr <- tryCatch(
    rtracklayer::import(gtf_path),
    error = function(e) {
      message("  rtracklayer failed, falling back to base parser.")
      NULL
    }
  )
  if (is.null(gr)) return(.parse_gtf_base(gtf_path))
  df     <- as.data.frame(gr)
  tx_col <- intersect(c("transcript_id","transcriptId"), names(df))[1]
  gn_col <- intersect(c("gene_id","geneId"), names(df))[1]
  if (is.na(tx_col)) stop("transcript_id not found in GTF.")

  # Exon-based tx_length
  exon_df <- df[df$type == "exon" & !is.na(df[[tx_col]]), ]
  if (nrow(exon_df) > 0) {
    exon_dt  <- data.table::data.table(
      transcript_id = as.character(exon_df[[tx_col]]),
      width         = as.integer(exon_df$width))
    len_dt <- exon_dt[, .(tx_length = sum(width, na.rm=TRUE)), by=transcript_id]
  } else {
    message("  Warning: no exon features; using transcript span.")
    tx_df  <- df[df$type == "transcript" & !is.na(df[[tx_col]]), ]
    len_dt <- data.table::data.table(
      transcript_id = as.character(tx_df[[tx_col]]),
      tx_length     = as.integer(tx_df$width))
  }
  tx_df  <- df[df$type == "transcript" & !is.na(df[[tx_col]]), ]
  meta   <- data.table::data.table(
    transcript_id = as.character(tx_df[[tx_col]]),
    gene_id = if (!is.na(gn_col)) as.character(tx_df[[gn_col]]) else NA_character_)
  meta   <- meta[!duplicated(transcript_id)]
  dt     <- merge(meta, len_dt, by="transcript_id", all.x=TRUE)
  dt[, is_novel := .is_novel_tx(transcript_id)]
  message(sprintf("  GTF: %d transcripts (%d novel)", nrow(dt), sum(dt$is_novel)))
  dt
}

#' @keywords internal
.parse_gtf_base <- function(gtf_path) {
  if (.is_gz(gtf_path)) {
    con <- gzfile(gtf_path, "rt"); lines <- readLines(con); close(con)
  } else {
    lines <- readLines(gtf_path)
  }
  lines      <- lines[!startsWith(lines, "#")]
  tx_lines   <- lines[grepl("\ttranscript\t", lines, fixed=TRUE)]
  exon_lines <- lines[grepl("\texon\t",       lines, fixed=TRUE)]
  if (length(tx_lines) == 0) {
    warning("No transcript entries found.")
    return(.empty_gtf_dt())
  }

  .attr <- function(s, key) {
    m <- regmatches(s, regexpr(paste0(key,'\\s+"([^"]+)"'), s, perl=TRUE))
    if (!length(m) || m=="") return(NA_character_)
    sub(paste0(key,'\\s+"([^"]+)"'), "\\1", m, perl=TRUE)
  }

  # Transcript metadata
  tx_rows <- lapply(tx_lines, function(l) {
    f <- strsplit(l,"\t",fixed=TRUE)[[1]]
    if (length(f)<9) return(NULL)
    list(transcript_id=.attr(f[9],"transcript_id"),
         gene_id=.attr(f[9],"gene_id"))
  })
  tx_rows <- Filter(Negate(is.null), tx_rows)
  tx_dt   <- data.table::rbindlist(tx_rows)[!duplicated(transcript_id)]

  # Exon-based tx_length
  if (length(exon_lines) > 0) {
    ex_rows <- lapply(exon_lines, function(l) {
      f <- strsplit(l,"\t",fixed=TRUE)[[1]]
      if (length(f)<9) return(NULL)
      s <- suppressWarnings(as.integer(f[4]))
      e <- suppressWarnings(as.integer(f[5]))
      tx_id <- .attr(f[9],"transcript_id")
      if (is.na(tx_id)||is.na(s)||is.na(e)) return(NULL)
      list(transcript_id=tx_id, exon_len=abs(e-s)+1L)
    })
    ex_rows <- Filter(Negate(is.null), ex_rows)
    if (length(ex_rows)>0) {
      ex_dt  <- data.table::rbindlist(ex_rows)
      len_dt <- ex_dt[, .(tx_length=sum(exon_len,na.rm=TRUE)), by=transcript_id]
    } else len_dt <- data.table::data.table(transcript_id=character(0), tx_length=integer(0))
  } else {
    message("  Warning: no exon entries; using transcript span.")
    sp_rows <- lapply(tx_lines, function(l) {
      f <- strsplit(l,"\t",fixed=TRUE)[[1]]
      if (length(f)<9) return(NULL)
      s <- suppressWarnings(as.integer(f[4]))
      e <- suppressWarnings(as.integer(f[5]))
      list(transcript_id=.attr(f[9],"transcript_id"),
           tx_length=if(!is.na(s)&&!is.na(e)) abs(e-s)+1L else NA_integer_)
    })
    len_dt <- data.table::rbindlist(Filter(Negate(is.null), sp_rows))
  }

  dt <- merge(tx_dt, len_dt, by="transcript_id", all.x=TRUE)
  dt[, is_novel := .is_novel_tx(transcript_id)]
  message(sprintf("  GTF: %d transcripts (%d novel)", nrow(dt), sum(dt$is_novel)))
  dt
}

#' Identify novel transcripts by ID prefix
#' @keywords internal
.is_novel_tx <- function(ids) {
  !grepl("^(ENS|NM_|NR_|XM_|XR_)", ids, perl=TRUE)
}

#' Empty GTF data.table
#' @keywords internal
.empty_gtf_dt <- function() {
  data.table::data.table(
    transcript_id=character(0), gene_id=character(0),
    is_novel=logical(0), tx_length=integer(0))
}


# ============================================================
# Long-read prior loading
# ============================================================

#' Load and validate long-read prior
#'
#' Accepts an IsoEMResult object, data.frame, or file path.
#' Aligns transcript IDs with the GTF, fills missing transcripts
#' with uniform (very weak) prior.
#'
#' @param lr_input IsoEMResult, data.frame, or character file path.
#' @param gtf_meta data.table from .parse_gtf().
#' @return data.table: transcript_id, pi_lr, em_count, certainty,
#'   lr_reliability.
#' @keywords internal
.load_lr_prior <- function(lr_input, gtf_meta) {
  message("  Loading long-read prior...")

  # --- Detect input type ---
  if (methods::is(lr_input, "IsoEMResult")) {
    if (!requireNamespace("IsoEM", quietly=TRUE))
      stop("IsoEM package required to use IsoEMResult input.")
    prior_raw <- IsoEM::as_spatial_prior(lr_input)
    dt <- data.table::data.table(
      transcript_id = names(prior_raw$pi_lr),
      em_count      = as.numeric(prior_raw$em_count),
      certainty     = as.numeric(prior_raw$certainty)
    )
  } else if (is.data.frame(lr_input) || data.table::is.data.table(lr_input)) {
    dt <- data.table::as.data.table(lr_input)
  } else if (is.character(lr_input) && file.exists(lr_input)) {
    dt <- .read_file(lr_input)
  } else {
    stop("lr_input must be an IsoEMResult, data.frame, or file path.")
  }

  # --- Column checks ---
  if (!"transcript_id" %in% names(dt))
    stop("lr_input must contain a 'transcript_id' column.")
  if (!"em_count" %in% names(dt))
    stop("lr_input must contain an 'em_count' column.")
  if (!"certainty" %in% names(dt)) {
    warning("'certainty' column missing; defaulting to 1.0.")
    dt[, certainty := 1.0]
  }

  dt[, transcript_id := as.character(transcript_id)]
  dt[, em_count      := as.numeric(em_count)]
  dt[, certainty     := as.numeric(certainty)]

  # --- Align with GTF ---
  gtf_ids <- gtf_meta$transcript_id
  n_before <- nrow(dt)
  # Transcripts in lr_input but not in GTF → warn and drop
  not_in_gtf <- dt$transcript_id[!(dt$transcript_id %in% gtf_ids)]
  if (length(not_in_gtf) > 0) {
    warning(sprintf(
      "%d transcripts in lr_input not found in GTF — dropped.",
      length(not_in_gtf)))
    dt <- dt[transcript_id %in% gtf_ids]
  }
  # Transcripts in GTF but not in lr_input → add with uniform prior
  missing_ids <- gtf_ids[!(gtf_ids %in% dt$transcript_id)]
  if (length(missing_ids) > 0) {
    message(sprintf(
      "  %d GTF transcripts absent from lr_input — assigned uniform prior.",
      length(missing_ids)))
    fill_dt <- data.table::data.table(
      transcript_id = missing_ids,
      em_count      = 0.01,
      certainty     = 0.0
    )
    dt <- data.table::rbindlist(list(dt, fill_dt), fill=TRUE)
  }

  # --- Compute pi_lr ---
  total <- sum(dt$em_count, na.rm=TRUE)
  dt[, pi_lr := em_count / max(total, 1e-9)]

  # --- LR reliability ---
  max_count <- max(dt$em_count, na.rm=TRUE)
  dt[, lr_reliability := (log(em_count + 1) / log(max_count + 1)) * certainty]

  message(sprintf(
    "  LR prior: %d transcripts | pi_lr range [%.2e, %.4f]",
    nrow(dt), min(dt$pi_lr), max(dt$pi_lr)))
  dt
}


# ============================================================
# Coordinate loading
# ============================================================

#' Load spot coordinates
#'
#' @param coords_path Character or data.frame.
#' @param format Character: "auto", "visium", "visium_hd", "patho_dbit", "generic".
#' @return data.table: spot_id, x, y.
#' @keywords internal
.load_coords <- function(coords_path, format = "auto") {
  if (is.data.frame(coords_path) || data.table::is.data.table(coords_path)) {
    dt <- data.table::as.data.table(coords_path)
  } else {
    if (!file.exists(coords_path))
      stop("spot_coords file not found: ", coords_path)
    # .parquet support
    if (grepl("\\.parquet$", coords_path, ignore.case=TRUE)) {
      if (!requireNamespace("arrow", quietly=TRUE))
        stop("Package 'arrow' required for .parquet coords file.")
      dt <- data.table::as.data.table(arrow::read_parquet(coords_path))
    } else {
      dt <- .read_file(coords_path)
    }
  }

  # --- Auto-detect format ---
  if (format == "auto") {
    if (all(c("barcode","pxl_col_in_fullres","pxl_row_in_fullres") %in% names(dt)))
      format <- "visium"
    else if (all(c("barcode","array_row","array_col") %in% names(dt)))
      format <- "visium"
    else
      format <- "generic"
  }

  if (format == "visium" || format == "visium_hd") {
    # Visium: barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
    id_col <- intersect(c("barcode","spot_id"), names(dt))[1]
    x_col  <- intersect(c("pxl_col_in_fullres","array_col","x"), names(dt))[1]
    y_col  <- intersect(c("pxl_row_in_fullres","array_row","y"), names(dt))[1]
    dt <- dt[, .(spot_id = as.character(get(id_col)),
                 x       = as.numeric(get(x_col)),
                 y       = as.numeric(get(y_col)))]
    if ("in_tissue" %in% names(dt)) dt <- dt[in_tissue == 1]
  } else {
    # Generic: requires spot_id, x, y
    need <- c("spot_id","x","y")
    miss <- need[!need %in% names(dt)]
    if (length(miss) > 0)
      stop("Generic coord file missing columns: ", paste(miss, collapse=", "))
    dt <- dt[, .(spot_id=as.character(spot_id), x=as.numeric(x), y=as.numeric(y))]
  }

  message(sprintf("  Coordinates: %d spots", nrow(dt)))
  dt
}


# ============================================================
# Validation helpers
# ============================================================

#' Validate main inputs
#' @keywords internal
.validate_inputs <- function(pe_bam, se_bam, gtf_file, lr_input,
                               sample_id, sample_id_file) {
  # Check sample_id first (fast, no file I/O)
  if (is.null(sample_id) && is.null(sample_id_file))
    stop("Provide either sample_id or sample_id_file.")
  if (!is.null(sample_id) && !is.null(sample_id_file))
    stop("Provide only one of sample_id or sample_id_file, not both.")
  # Check BAM presence
  if (is.null(pe_bam) && is.null(se_bam))
    stop("Provide at least one of pe_bam or se_bam.")
  if (!is.null(pe_bam) && !file.exists(pe_bam))
    stop("pe_bam not found: ", pe_bam)
  if (!is.null(se_bam) && !file.exists(se_bam))
    stop("se_bam not found: ", se_bam)
  if (!file.exists(gtf_file))
    stop("gtf_file not found: ", gtf_file)
  if (!is.null(sample_id_file) && !file.exists(sample_id_file))
    stop("sample_id_file not found: ", sample_id_file)
  invisible(TRUE)
}
