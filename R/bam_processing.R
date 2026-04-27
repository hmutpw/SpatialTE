#' @title BAM Processing: Stream Once for g(x), f(d), and EC Data
#' @description
#' Reads PE and/or SE toTranscriptome BAMs a single time using Rsamtools,
#' simultaneously collecting:
#'   - Coverage profiles for g(x) estimation (NH=1 reads only)
#'   - Read length distribution for f(d) estimation (NH=1 reads only)
#'   - Raw read-to-transcript assignments for EC construction (all reads)
#' @name bam_processing
NULL

# ============================================================
# BC/UMI extraction helpers
# ============================================================

#' Detect CB and UB tags from BAM header and first reads
#' @keywords internal
.detect_bc_umi_tags <- function(bam_path) {
  bf    <- Rsamtools::BamFile(bam_path)
  hdr   <- Rsamtools::scanBamHeader(bf)[[1]]
  # Check @CO lines for tag hints
  co    <- unlist(hdr$text[names(hdr$text) == "@CO"])
  cb_hint <- if (any(grepl("CB:", co))) "CB" else NULL
  ub_hint <- if (any(grepl("UB:", co))) "UB" else NULL

  # Scan first 1000 reads
  param <- Rsamtools::ScanBamParam(
    what = c("qname","flag"),
    tag  = c("CB","UB","CR","UR","XC","XM")
  )
  chunk <- suppressWarnings(
    Rsamtools::scanBam(bf, param = param)[[1]]
  )
  tags_present <- names(chunk$tag)[sapply(chunk$tag, function(x) any(!is.na(x)))]

  cb_tag <- if ("CB" %in% tags_present) "CB" else
            if ("CR" %in% tags_present) "CR" else
            if ("XC" %in% tags_present) "XC" else NULL
  ub_tag <- if ("UB" %in% tags_present) "UB" else
            if ("UR" %in% tags_present) "UR" else
            if ("XM" %in% tags_present) "XM" else NULL

  if (is.null(cb_tag)) stop(
    "No barcode tag (CB/CR/XC) found in BAM. ",
    "Use cb_tag parameter or provide bc_umi_table.")
  if (is.null(ub_tag)) stop(
    "No UMI tag (UB/UR/XM) found in BAM. ",
    "Use ub_tag parameter or provide bc_umi_table.")

  message(sprintf("  Auto-detected BAM tags: CB=%s, UB=%s", cb_tag, ub_tag))
  list(cb_tag = cb_tag, ub_tag = ub_tag)
}


#' Parse BC/UMI from read ID using a regex pattern
#' @keywords internal
.parse_bc_umi_from_readid <- function(qnames, pattern) {
  m <- regmatches(qnames, regexec(pattern, qnames, perl=TRUE))
  bc  <- sapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_)
  umi <- sapply(m, function(x) if (length(x) >= 3) x[3] else NA_character_)
  list(barcode = bc, umi = umi)
}


# ============================================================
# Core: stream BAM(s) once
# ============================================================

#' Stream one or two toTranscriptome BAMs collecting all needed data
#'
#' One pass over each BAM:
#'   NH=1 reads → update coverage bins (g_cov) and read-length table (frag_len)
#'   All reads  → update EC raw data table
#'
#' @param pe_bam Character or NULL.
#' @param se_bam Character or NULL.
#' @param transcript_names Character vector indexed by integer (t_idx).
#' @param gtf_meta data.table from .parse_gtf() (for length bin assignment).
#' @param bc_umi_table data.table (read_id|barcode|umi) or NULL.
#' @param bc_umi_pattern Character regex or NULL.
#' @param cb_tag Character.
#' @param ub_tag Character.
#' @param cb_whitelist Character vector or NULL.
#' @param frag_dist_source "pe", "se", or "both".
#' @param min_mapq Integer.
#' @param chunk_size Integer. Reads per Rsamtools chunk (default 500000).
#' @return list:
#'   g_cov: list(short, medium, long, vlong) — cumulative coverage arrays
#'   frag_len_tab: named integer vector (read_length → count)
#'   ec_raw: data.table(spot_key, t_idx, nh)
#'   transcript_names: character vector (may be extended)
#' @keywords internal
.stream_bams <- function(pe_bam, se_bam,
                          transcript_names,
                          gtf_meta,
                          bc_umi_table    = NULL,
                          bc_umi_pattern  = NULL,
                          cb_tag          = "CB",
                          ub_tag          = "UB",
                          cb_whitelist    = NULL,
                          frag_dist_source = "pe",
                          min_mapq        = 0L,
                          chunk_size      = 500000L) {

  n_bins <- 100L
  g_cov  <- list(
    short  = numeric(n_bins),
    medium = numeric(n_bins),
    long   = numeric(n_bins),
    vlong  = numeric(n_bins)
  )
  frag_len_tab <- integer(0)
  ec_raw_list  <- list()

  # Build transcript → length-bin lookup
  tx_bin_map <- .build_tx_bin_map(transcript_names, gtf_meta)

  # Build bc_umi lookup table if provided
  if (!is.null(bc_umi_table)) {
    bc_umi_dt <- data.table::as.data.table(bc_umi_table)
    data.table::setnames(bc_umi_dt,
      old = names(bc_umi_dt)[1:3],
      new = c("qname","barcode","umi"))
    data.table::setkeyv(bc_umi_dt, "qname")
  } else {
    bc_umi_dt <- NULL
  }

  # CB whitelist set for fast lookup
  wl_set <- if (!is.null(cb_whitelist)) {
    if (is.character(cb_whitelist) && length(cb_whitelist)==1 &&
        file.exists(cb_whitelist)) {
      as.character(data.table::fread(cb_whitelist, header=FALSE)[[1]])
    } else {
      as.character(cb_whitelist)
    }
  } else NULL

  bams_to_process <- list()
  if (!is.null(pe_bam)) bams_to_process[["pe"]] <- pe_bam
  if (!is.null(se_bam)) bams_to_process[["se"]] <- se_bam

  for (bam_type in names(bams_to_process)) {
    bam_path <- bams_to_process[[bam_type]]
    message(sprintf("  Streaming %s BAM: %s", toupper(bam_type),
                    basename(bam_path)))

    # Auto-detect tags if needed
    if (is.null(bc_umi_table) && is.null(bc_umi_pattern)) {
      detected <- tryCatch(
        .detect_bc_umi_tags(bam_path),
        error = function(e) list(cb_tag=cb_tag, ub_tag=ub_tag)
      )
      cb_tag_use <- detected$cb_tag
      ub_tag_use <- detected$ub_tag
    } else {
      cb_tag_use <- cb_tag
      ub_tag_use <- ub_tag
    }

    bf    <- Rsamtools::BamFile(bam_path, yieldSize = chunk_size)
    param <- Rsamtools::ScanBamParam(
      what      = c("qname","rname","pos","qwidth","mapq","flag"),
      tag       = unique(c(cb_tag_use, ub_tag_use, "NH")),
      mapqFilter = min_mapq
    )

    open(bf)
    chunk_idx <- 0L
    repeat {
      chunk <- Rsamtools::scanBam(bf, param = param)[[1]]
      if (length(chunk$qname) == 0) break
      chunk_idx <- chunk_idx + 1L

      # Convert to data.table for fast operations
      dt <- data.table::data.table(
        qname  = chunk$qname,
        rname  = as.character(chunk$rname),
        pos    = chunk$pos,
        qwidth = chunk$qwidth,
        mapq   = chunk$mapq,
        flag   = chunk$flag
      )

      # NH tag
      nh_vec <- chunk$tag[["NH"]]
      if (is.null(nh_vec)) nh_vec <- rep(1L, nrow(dt))
      dt[, nh := as.integer(nh_vec)]

      # Filter unmapped
      dt <- dt[!is.na(rname) & !is.na(pos)]
      if (nrow(dt) == 0) next

      # --- BC/UMI extraction ---
      if (!is.null(bc_umi_dt)) {
        # Method 1: lookup table
        dt_bc <- bc_umi_dt[.(dt$qname), nomatch=NULL]
        dt[, barcode := dt_bc$barcode[match(qname, dt_bc$qname)]]
        dt[, umi     := dt_bc$umi[match(qname, dt_bc$qname)]]
      } else if (!is.null(bc_umi_pattern)) {
        # Method 2: regex from read ID
        parsed <- .parse_bc_umi_from_readid(dt$qname, bc_umi_pattern)
        dt[, barcode := parsed$barcode]
        dt[, umi     := parsed$umi]
      } else {
        # Method 3: BAM tags
        cb_vec <- chunk$tag[[cb_tag_use]]
        ub_vec <- chunk$tag[[ub_tag_use]]
        dt[, barcode := if (!is.null(cb_vec)) as.character(cb_vec) else NA_character_]
        dt[, umi     := if (!is.null(ub_vec)) as.character(ub_vec) else NA_character_]
      }

      # Filter missing BC/UMI
      dt <- dt[!is.na(barcode) & barcode != "" & !is.na(umi) & umi != ""]
      if (nrow(dt) == 0) next

      # Apply CB whitelist
      if (!is.null(wl_set)) dt <- dt[barcode %in% wl_set]
      if (nrow(dt) == 0) next

      # Unique spot key = barcode_umi
      dt[, spot_key := paste0(barcode, "_", umi)]

      # Transcript integer index
      rname_levels <- levels(factor(dt$rname))
      new_tx <- rname_levels[!(rname_levels %in% transcript_names)]
      if (length(new_tx) > 0) {
        transcript_names <- c(transcript_names, new_tx)
        tx_bin_map <- .build_tx_bin_map(transcript_names, gtf_meta)
      }
      dt[, t_idx := match(rname, transcript_names)]
      dt <- dt[!is.na(t_idx)]
      if (nrow(dt) == 0) next

      # -------------------------------------------------------
      # g(x) and f(d) from NH=1 reads only
      # -------------------------------------------------------
      use_for_bias <- (bam_type == "pe" && frag_dist_source %in% c("pe","both")) ||
                      (bam_type == "se" && frag_dist_source %in% c("se","both"))

      if (use_for_bias) {
        nh1 <- dt[nh == 1L]
        if (nrow(nh1) > 0) {
          # f(d): query_length distribution
          ql_tab <- table(nh1$qwidth)
          for (nm in names(ql_tab)) {
            idx <- as.integer(nm)
            if (is.na(frag_len_tab[as.character(idx)])) {
              frag_len_tab[as.character(idx)] <- 0L
            }
            frag_len_tab[as.character(idx)] <-
              frag_len_tab[as.character(idx)] + as.integer(ql_tab[[nm]])
          }

          # g(x): coverage per length bin
          for (bin_name in c("short","medium","long","vlong")) {
            bin_t_idx <- which(tx_bin_map == bin_name)
            sub <- nh1[t_idx %in% bin_t_idx]
            if (nrow(sub) == 0) next
            # Normalise position to [0, 1) and bin into n_bins
            # tx_length from gtf_meta
            sub_lengths <- gtf_meta$tx_length[match(
              transcript_names[sub$t_idx], gtf_meta$transcript_id)]
            valid <- !is.na(sub_lengths) & sub_lengths > 0 & !is.na(sub$pos)
            if (!any(valid)) next
            sub2   <- sub[valid]
            tx_len <- sub_lengths[valid]
            norm_pos <- pmax(pmin((sub2$pos - 1L) / tx_len, 0.9999), 0)
            bin_pos  <- pmin(floor(norm_pos * n_bins) + 1L, n_bins)
            cov_add  <- tabulate(bin_pos, nbins = n_bins)
            g_cov[[bin_name]] <- g_cov[[bin_name]] + as.numeric(cov_add)
          }
        }
      }

      # -------------------------------------------------------
      # EC raw data: all reads
      # -------------------------------------------------------
      ec_raw_list[[length(ec_raw_list)+1]] <- dt[, .(spot_key, t_idx, nh)]
    }
    close(bf)
    message(sprintf("    %s BAM done (%d chunks)", toupper(bam_type), chunk_idx))
  }

  # Combine EC raw
  ec_raw <- if (length(ec_raw_list) > 0)
    data.table::rbindlist(ec_raw_list)
  else
    data.table::data.table(spot_key=character(0),
                            t_idx=integer(0), nh=integer(0))

  list(
    g_cov            = g_cov,
    frag_len_tab     = frag_len_tab,
    ec_raw           = ec_raw,
    transcript_names = transcript_names
  )
}


# ============================================================
# Effective length computation
# ============================================================

#' Compute effective lengths from g(x) and f(d)
#'
#' eff_len(k) = Σ_{x=0}^{l_k} g(x; bin_k) × F(l_k - x)
#' where F is the CDF of f(d).
#'
#' @param gtf_meta data.table.
#' @param g_cov list of 4 numeric(100) coverage arrays.
#' @param frag_len_tab named integer vector (read_length → count).
#' @param min_eff_len Numeric. Floor for effective length.
#' @return named numeric vector: transcript_id → eff_len.
#' @keywords internal
.compute_eff_len <- function(gtf_meta, g_cov, frag_len_tab, min_eff_len = 50) {
  message("  Computing effective lengths...")
  n_bins <- 100L

  # Normalise g(x) per bin (mean = 1)
  g_norm <- lapply(g_cov, function(cv) {
    if (sum(cv) < 1e-9) return(rep(1.0, n_bins))
    cv / mean(cv)
  })

  # Normalise f(d) to a probability distribution
  if (length(frag_len_tab) == 0) {
    message("  Warning: no reads for f(d) estimation. Using uniform f(d).")
    f_d <- setNames(rep(1.0/100, 100), as.character(1:100))
  } else {
    total_reads <- sum(frag_len_tab)
    f_d <- frag_len_tab / total_reads
  }
  # CDF of f(d)
  d_vals  <- as.integer(names(f_d))
  max_d   <- max(d_vals)
  f_vec   <- numeric(max_d)
  f_vec[d_vals] <- f_d
  F_cdf   <- cumsum(f_vec)   # F_cdf[d] = P(fragment_length <= d)

  # Compute eff_len per transcript
  eff_lens <- numeric(nrow(gtf_meta))
  names(eff_lens) <- gtf_meta$transcript_id

  for (i in seq_len(nrow(gtf_meta))) {
    l_k     <- gtf_meta$tx_length[i]
    if (is.na(l_k) || l_k <= 0) { eff_lens[i] <- min_eff_len; next }
    bin_k   <- .tx_length_bin(l_k)
    g_curve <- g_norm[[bin_k]]

    # Interpolate g(x) from 100 bins to l_k positions
    if (l_k <= n_bins) {
      # Upsample: repeat each bin proportionally
      x_norm <- seq(0, 1, length.out = l_k + 1L)[-(l_k+1L)]
      g_interp <- g_curve[pmax(1L, floor(x_norm * n_bins) + 1L)]
    } else {
      # Downsample via linear interpolation
      x_old   <- seq(0, 1, length.out = n_bins)
      x_new   <- seq(0, 1, length.out = l_k)
      g_interp <- approx(x_old, g_curve, xout = x_new)$y
    }

    # eff_len = Σ_x g(x) × F(l_k - x)
    remaining <- l_k - seq_len(l_k) + 1L  # l_k - x for x = 1..l_k
    f_vals    <- ifelse(remaining <= 0, 0,
                  ifelse(remaining > max_d, 1.0, F_cdf[remaining]))
    eff_lens[i] <- sum(g_interp * f_vals, na.rm = TRUE)
  }

  eff_lens <- pmax(eff_lens, min_eff_len)

  pct_low <- mean(eff_lens <= min_eff_len * 1.1) * 100
  if (pct_low > 5) message(sprintf(
    "  Warning: %.1f%% of transcripts have eff_len near floor (%d bp).",
    pct_low, min_eff_len))

  eff_lens
}

#' Assign transcript to length bin
#' @keywords internal
.tx_length_bin <- function(len) {
  if      (is.na(len) || len <  500L) "short"
  else if (len < 1500L) "medium"
  else if (len < 3000L) "long"
  else                  "vlong"
}

#' Build transcript → length bin lookup
#' @keywords internal
.build_tx_bin_map <- function(tx_names, gtf_meta) {
  lens <- gtf_meta$tx_length[match(tx_names, gtf_meta$transcript_id)]
  sapply(lens, .tx_length_bin)
}
