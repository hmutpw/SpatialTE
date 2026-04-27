library(testthat)
library(SpatialTE)
library(data.table)
library(Matrix)

gtf_f      <- system.file("extdata", "toy_gtf.gtf",          package="SpatialTE")
gtf_gz_f   <- system.file("extdata", "toy_gtf.gtf.gz",       package="SpatialTE")
lr_f       <- system.file("extdata", "toy_lr_counts.tsv",    package="SpatialTE")
lr_gz_f    <- system.file("extdata", "toy_lr_counts.tsv.gz", package="SpatialTE")
coords_f   <- system.file("extdata", "toy_coords.csv",       package="SpatialTE")

# ============================================================
# GTF parsing
# ============================================================

test_that(".parse_gtf returns correct columns", {
  dt <- SpatialTE:::.parse_gtf(gtf_f)
  expect_true(all(c("transcript_id","gene_id","is_novel","tx_length") %in% names(dt)))
})

test_that(".parse_gtf detects novel transcripts correctly", {
  dt <- SpatialTE:::.parse_gtf(gtf_f)
  enst_rows  <- dt[grepl("^ENST", transcript_id)]
  novel_rows <- dt[grepl("^novel", transcript_id)]
  expect_true(all(!enst_rows$is_novel))
  expect_true(all(novel_rows$is_novel))
})

test_that(".parse_gtf tx_length is exon sum not genomic span", {
  dt <- SpatialTE:::.parse_gtf(gtf_f)
  # Each toy transcript has 2 exons with a 100 bp intron gap
  # So tx_length < genomic span for all multi-exon transcripts
  expect_true(all(dt$tx_length > 0, na.rm=TRUE))
  # Spot check: ENST00000000001 should have tx_length = 1200 (= tx_lens[1])
  r1 <- dt[transcript_id == "ENST00000000001"]
  expect_equal(r1$tx_length, 1200L)
})

test_that(".parse_gtf works with .gz input", {
  dt1 <- SpatialTE:::.parse_gtf(gtf_f)
  dt2 <- SpatialTE:::.parse_gtf(gtf_gz_f)
  expect_equal(nrow(dt1), nrow(dt2))
  expect_equal(dt1$transcript_id, dt2$transcript_id)
})

test_that(".parse_gtf returns 10 transcripts for toy data", {
  dt <- SpatialTE:::.parse_gtf(gtf_f)
  expect_equal(nrow(dt), 10L)
})

# ============================================================
# LR prior loading
# ============================================================

test_that(".load_lr_prior works with file path", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  pr  <- SpatialTE:::.load_lr_prior(lr_f, gtf)
  expect_true(all(c("transcript_id","pi_lr","em_count","certainty","lr_reliability") %in% names(pr)))
})

test_that(".load_lr_prior pi_lr sums to 1", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  pr  <- SpatialTE:::.load_lr_prior(lr_f, gtf)
  expect_equal(sum(pr$pi_lr), 1, tolerance=1e-6)
})

test_that(".load_lr_prior works with .gz input", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  pr1 <- SpatialTE:::.load_lr_prior(lr_f,    gtf)
  pr2 <- SpatialTE:::.load_lr_prior(lr_gz_f, gtf)
  expect_equal(pr1$pi_lr, pr2$pi_lr)
})

test_that(".load_lr_prior works with data.frame input", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  df  <- as.data.frame(fread(lr_f))
  pr  <- SpatialTE:::.load_lr_prior(df, gtf)
  expect_true(nrow(pr) > 0)
  expect_equal(sum(pr$pi_lr), 1, tolerance=1e-6)
})

test_that(".load_lr_prior fills missing certainty with 1.0", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  df  <- fread(lr_f)[, .(transcript_id, em_count)]
  expect_warning(
    pr <- SpatialTE:::.load_lr_prior(df, gtf),
    "certainty"
  )
  expect_true(all(pr$certainty == 1.0))
})

test_that(".load_lr_prior fills transcripts missing from lr_input with uniform prior", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  df  <- fread(lr_f)[1:5]   # only first 5 transcripts
  pr  <- SpatialTE:::.load_lr_prior(df, gtf)
  expect_equal(nrow(pr), 10L)   # all 10 GTF transcripts present
  missing_ids <- setdiff(gtf$transcript_id, df$transcript_id)
  miss_pr <- pr[transcript_id %in% missing_ids]
  expect_true(all(miss_pr$certainty == 0.0))
})

test_that(".load_lr_prior errors on missing transcript_id column", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  df  <- data.frame(em_count=1:5)
  expect_error(SpatialTE:::.load_lr_prior(df, gtf), "transcript_id")
})

test_that(".load_lr_prior errors on missing em_count column", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  df  <- data.frame(transcript_id=letters[1:5])
  expect_error(SpatialTE:::.load_lr_prior(df, gtf), "em_count")
})

# ============================================================
# Coordinate loading
# ============================================================

test_that(".load_coords loads generic coords", {
  dt <- SpatialTE:::.load_coords(coords_f, format="generic")
  expect_true(all(c("spot_id","x","y") %in% names(dt)))
  expect_equal(nrow(dt), 9L)
})

test_that(".load_coords auto-detects generic format", {
  dt <- SpatialTE:::.load_coords(coords_f, format="auto")
  expect_equal(nrow(dt), 9L)
})

test_that(".load_coords works with data.frame input", {
  df <- data.frame(spot_id=paste0("S",1:5), x=1:5, y=1:5)
  dt <- SpatialTE:::.load_coords(df)
  expect_equal(nrow(dt), 5L)
})

# ============================================================
# Neighbor graph
# ============================================================

test_that(".build_neighbor_graph returns correct structure", {
  coords <- SpatialTE:::.load_coords(coords_f)
  ng     <- SpatialTE:::.build_neighbor_graph(coords)
  expect_true(all(c("neighbors","radius","median_dist") %in% names(ng)))
  expect_equal(length(ng$neighbors), 9L)
})

test_that(".build_neighbor_graph custom radius includes expected neighbors", {
  coords <- SpatialTE:::.load_coords(coords_f)
  # Spots are on a 100-unit grid; radius=150 should find diagonal neighbors
  ng <- SpatialTE:::.build_neighbor_graph(coords, neighbor_radius=150)
  # Corner spots should have 3 neighbors (adjacent only)
  n1 <- length(ng$neighbors[["SPOT_001"]]$idx)
  expect_true(n1 >= 2L)
})

# ============================================================
# Effective length computation
# ============================================================

test_that(".compute_eff_len returns positive values", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  # Minimal g_cov and f_d
  g_cov <- list(
    short  = rep(1.0, 100),
    medium = rep(1.0, 100),
    long   = rep(1.0, 100),
    vlong  = rep(1.0, 100)
  )
  frag_len_tab <- c("90"=500L, "100"=800L, "120"=600L, "150"=300L)
  el <- SpatialTE:::.compute_eff_len(gtf, g_cov, frag_len_tab, min_eff_len=50)
  expect_true(all(el >= 50))
  expect_equal(length(el), nrow(gtf))
  expect_equal(names(el), gtf$transcript_id)
})

test_that(".compute_eff_len respects min_eff_len floor", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  # Very short read length → most eff_len at floor
  g_cov <- list(short=rep(1,100), medium=rep(1,100),
                long=rep(1,100), vlong=rep(1,100))
  frag_len_tab <- c("500"=100L)  # all reads length 500
  el <- SpatialTE:::.compute_eff_len(gtf, g_cov, frag_len_tab, min_eff_len=50)
  expect_true(all(el >= 50))
})

test_that(".compute_eff_len handles empty frag_len_tab with warning", {
  gtf <- SpatialTE:::.parse_gtf(gtf_f)
  g_cov <- list(short=rep(1,100), medium=rep(1,100),
                long=rep(1,100), vlong=rep(1,100))
  expect_message(
    el <- SpatialTE:::.compute_eff_len(gtf, g_cov, integer(0)),
    "Warning"
  )
  expect_true(all(el >= 50))
})

# ============================================================
# EC construction
# ============================================================

.make_toy_ec_raw <- function(n_spots=5L, n_reads_per_spot=20L) {
  set.seed(123)
  tx_names <- c(paste0("ENST", sprintf("%011d", 1:8)), "novel.1", "novel.2")
  spot_ids <- paste0("SPOT_", sprintf("%03d", seq_len(n_spots)))
  rows <- list()
  for (si in seq_along(spot_ids)) {
    for (r in seq_len(n_reads_per_spot)) {
      bc  <- spot_ids[si]
      umi <- paste0("UMI", sprintf("%04d", r))
      key <- paste0(bc, "_", umi)
      # 70% unique, 30% multi-mapping
      if (runif(1) < 0.7) {
        t_idx <- sample(seq_along(tx_names), 1L)
      } else {
        t_idx <- sample(seq_along(tx_names), sample(2:3, 1L))
      }
      for (ti in t_idx) {
        rows[[length(rows)+1]] <- list(spot_key=key, t_idx=ti, nh=length(t_idx))
      }
    }
  }
  list(ec_raw=data.table::rbindlist(rows), tx_names=tx_names)
}

test_that(".build_ec_one_sample returns SpotECCounts", {
  toy  <- .make_toy_ec_raw()
  ec_r <- toy$ec_raw
  ec_r[, sample_id := "sA"]
  sub  <- ec_r
  res  <- SpatialTE:::.build_ec_one_sample(sub, toy$tx_names, "sA",
                                             "barcode_umi", 200L)
  expect_s3_class(res, "SpotECCounts")
  expect_equal(res$sample_id, "sA")
})

test_that("SpotECCounts has correct dimensions", {
  toy <- .make_toy_ec_raw(n_spots=5L)
  sub <- toy$ec_raw; sub[, sample_id := "sA"]
  res <- SpatialTE:::.build_ec_one_sample(sub, toy$tx_names, "sA",
                                           "barcode_umi", 200L)
  expect_equal(ncol(res$unique_counts), length(toy$tx_names))
  expect_true(nrow(res$spot_ec_counts) > 0)
})

test_that("unique_counts only contains EC-size-1 entries", {
  toy <- .make_toy_ec_raw(n_spots=3L, n_reads_per_spot=50L)
  sub <- toy$ec_raw; sub[, sample_id := "sA"]
  res <- SpatialTE:::.build_ec_one_sample(sub, toy$tx_names, "sA",
                                           "barcode_umi", 200L)
  uniq_ecs <- res$ec_index$ec_id[res$ec_index$ec_size == 1L]
  # unique_counts total should equal sum of unique EC counts
  uniq_total <- sum(res$spot_ec_counts[, as.character(uniq_ecs), drop=FALSE])
  expect_equal(sum(res$unique_counts), uniq_total, tolerance=0.01)
})

test_that("overlimit ECs are flagged but not discarded", {
  # Create an EC with > max_ec_size transcripts
  ec_raw <- data.table(
    spot_key = rep("SPOT_001_UMI0001", 10L),
    t_idx    = 1:10,
    nh       = 10L
  )
  tx_names <- paste0("TX", 1:15)
  res <- SpatialTE:::.build_ec_one_sample(ec_raw, tx_names, "sA",
                                           "barcode_umi", max_ec_size=5L)
  expect_true(length(res$overlimit_ecs) > 0)
  # The overlimit EC should still appear in spot_ec_counts
  expect_true(sum(res$spot_ec_counts) > 0)
})

# ============================================================
# Prior construction
# ============================================================

test_that(".build_prior returns alpha matrix", {
  gtf    <- SpatialTE:::.parse_gtf(gtf_f)
  lr_pr  <- SpatialTE:::.load_lr_prior(lr_f, gtf)
  toy    <- .make_toy_ec_raw()
  sub    <- toy$ec_raw; sub[, sample_id := "sA"]
  ec_c   <- SpatialTE:::.build_ec_one_sample(sub, toy$tx_names, "sA",
                                              "barcode_umi", 200L)
  el     <- SpatialTE:::.compute_eff_len(
    gtf,
    list(short=rep(1,100),medium=rep(1,100),long=rep(1,100),vlong=rep(1,100)),
    c("90"=500L,"120"=800L)
  )
  pr     <- SpatialTE:::.build_prior(lr_pr, ec_c, el, gamma_base=10.0)
  expect_true(!is.null(pr$alpha_mat))
  expect_true(methods::is(pr$alpha_mat, "sparseMatrix"))
  expect_true(all(pr$alpha_mat@x > 0))
})

test_that("gamma calibration returns numeric scalar", {
  gtf   <- SpatialTE:::.parse_gtf(gtf_f)
  lr_pr <- SpatialTE:::.load_lr_prior(lr_f, gtf)
  toy   <- .make_toy_ec_raw()
  sub   <- toy$ec_raw; sub[, sample_id := "sA"]
  ec_c  <- SpatialTE:::.build_ec_one_sample(sub, toy$tx_names, "sA",
                                             "barcode_umi", 200L)
  el    <- SpatialTE:::.compute_eff_len(
    gtf,
    list(short=rep(1,100),medium=rep(1,100),long=rep(1,100),vlong=rep(1,100)),
    c("100"=1000L)
  )
  gbase <- SpatialTE:::.calibrate_gamma(lr_pr, ec_c, el)
  expect_true(is.numeric(gbase))
  expect_length(gbase, 1L)
  expect_true(gbase >= 1.0 && gbase <= 100.0)
})

# ============================================================
# MAP-EM
# ============================================================

test_that(".em_one_spot returns named list with expected elements", {
  toy    <- .make_toy_ec_raw(n_spots=1L, n_reads_per_spot=30L)
  sub    <- toy$ec_raw; sub[, sample_id := "sA"]
  ec_c   <- SpatialTE:::.build_ec_one_sample(sub, toy$tx_names, "sA",
                                              "barcode_umi", 200L)
  n_tx   <- length(toy$tx_names)
  ec_idx <- ec_c$ec_index
  ec_iso_idx <- lapply(ec_idx$ec_key, function(k)
    as.integer(strsplit(k,"|",fixed=TRUE)[[1]]))
  names(ec_iso_idx) <- ec_idx$ec_id

  si      <- 1L
  ec_row  <- ec_c$spot_ec_counts[si,]
  active  <- which(ec_row > 0)
  alpha_s <- rep(0.1, n_tx)
  el      <- rep(500.0, n_tx)
  pi_lr   <- setNames(rep(1/n_tx, n_tx), seq_len(n_tx))

  res <- SpatialTE:::.em_one_spot(
    ec_ids_active = colnames(ec_c$spot_ec_counts)[active],
    counts_active = as.numeric(ec_row[active]),
    ec_iso_idx    = ec_iso_idx,
    n_tx          = n_tx,
    eff_len_vec   = el,
    alpha_s       = alpha_s,
    pi_lr         = pi_lr,
    overlimit_set = character(0),
    max_iter      = 100L,
    tol           = 1e-4,
    sid           = "SPOT_001"
  )
  expect_true(all(c("sid","t_idx","counts","n_iter","converged") %in% names(res)))
  expect_true(all(res$counts >= 0))
})

test_that("EM counts are non-negative", {
  gtf   <- SpatialTE:::.parse_gtf(gtf_f)
  lr_pr <- SpatialTE:::.load_lr_prior(lr_f, gtf)
  toy   <- .make_toy_ec_raw(n_spots=5L, n_reads_per_spot=40L)
  sub   <- toy$ec_raw; sub[, sample_id := "sA"]
  ec_c  <- SpatialTE:::.build_ec_one_sample(sub, toy$tx_names, "sA",
                                             "barcode_umi", 200L)
  el    <- SpatialTE:::.compute_eff_len(
    gtf,
    list(short=rep(1,100),medium=rep(1,100),long=rep(1,100),vlong=rep(1,100)),
    c("100"=1000L)
  )
  pr   <- SpatialTE:::.build_prior(lr_pr, ec_c, el, gamma_base=5.0)
  em_r <- SpatialTE:::.run_em_all_spots(ec_c, pr, el,
                                         max_iter=50L, tol=1e-4, n_cores=1L)
  expect_true(all(em_r$count_mat@x >= 0))
})

# ============================================================
# Result object construction
# ============================================================

.make_toy_result <- function() {
  gtf    <- SpatialTE:::.parse_gtf(gtf_f)
  lr_pr  <- SpatialTE:::.load_lr_prior(lr_f, gtf)
  toy    <- .make_toy_ec_raw(n_spots=9L, n_reads_per_spot=50L)
  sub    <- toy$ec_raw; sub[, sample_id := "sA"]
  ec_c   <- SpatialTE:::.build_ec_one_sample(sub, toy$tx_names, "sA",
                                              "barcode_umi", 200L)
  el     <- SpatialTE:::.compute_eff_len(
    gtf,
    list(short=rep(1,100),medium=rep(1,100),long=rep(1,100),vlong=rep(1,100)),
    c("100"=1000L)
  )
  pr     <- SpatialTE:::.build_prior(lr_pr, ec_c, el, gamma_base=5.0)
  em_r   <- SpatialTE:::.run_em_all_spots(ec_c, pr, el,
                                           max_iter=50L, tol=1e-4, n_cores=1L)
  result <- SpatialTE:::.build_result(
    em_result    = em_r,
    ec_counts    = ec_c,
    prior_result = pr,
    gtf_meta     = gtf,
    eff_len      = el,
    spot_coords  = NULL,
    lr_prior     = lr_pr,
    gamma_base   = 5.0,
    metadata     = list(sample_id="sA", ec_data=ec_c),
    min_eff_len  = 50
  )
  result
}

test_that("SpatialTEResult is created without error", {
  expect_no_error(r <- .make_toy_result())
  expect_s3_class(r, "SpatialTEResult")
})

test_that("print.SpatialTEResult works", {
  r <- .make_toy_result()
  expect_output(print(r), "SpatialTEResult")
})

test_that("summary.SpatialTEResult works", {
  r <- .make_toy_result()
  expect_output(summary(r), "SpatialTEResult")
})

# ============================================================
# Output functions
# ============================================================

test_that("write_counts creates file with expected columns", {
  r    <- .make_toy_result()
  tmpf <- tempfile(fileext=".tsv")
  write_counts(r, tmpf)
  expect_true(file.exists(tmpf))
  dt <- fread(tmpf)
  expect_true(all(c("spot_id","transcript_id","em_count","tpm") %in% names(dt)))
  unlink(tmpf)
})

test_that("write_counts supports .gz output", {
  r    <- .make_toy_result()
  tmpf <- tempfile(fileext=".tsv.gz")
  write_counts(r, tmpf)
  expect_true(file.exists(tmpf))
  dt <- fread(tmpf)
  expect_true(nrow(dt) > 0)
  unlink(tmpf)
})

test_that("write_count_matrix creates wide matrix", {
  r    <- .make_toy_result()
  tmpf <- tempfile(fileext=".tsv")
  write_count_matrix(r, tmpf)
  expect_true(file.exists(tmpf))
  dt <- fread(tmpf)
  expect_true("transcript_id" %in% names(dt))
  unlink(tmpf)
})

test_that("write_eff_len creates file with eff_len column", {
  r    <- .make_toy_result()
  tmpf <- tempfile(fileext=".tsv")
  write_eff_len(r, tmpf)
  dt <- fread(tmpf)
  expect_true(all(c("transcript_id","tx_length","eff_len","low_eff_len") %in% names(dt)))
  expect_true(all(dt$eff_len >= 50))
  unlink(tmpf)
})

test_that("write_qc creates file with spot QC", {
  r    <- .make_toy_result()
  tmpf <- tempfile(fileext=".tsv")
  write_qc(r, tmpf)
  dt <- fread(tmpf)
  expect_true(all(c("total_reads","unique_frac","converged") %in% names(dt)))
  unlink(tmpf)
})

test_that("write_spatialTE creates all output files", {
  r      <- .make_toy_result()
  tmpdir <- tempfile()
  write_spatialTE(r, tmpdir, compress=FALSE, write_sharing=FALSE)
  expect_true(file.exists(file.path(tmpdir, "counts.tsv")))
  expect_true(file.exists(file.path(tmpdir, "count_matrix.tsv")))
  expect_true(file.exists(file.path(tmpdir, "efflen_table.tsv")))
  expect_true(file.exists(file.path(tmpdir, "qc_summary.tsv")))
  unlink(tmpdir, recursive=TRUE)
})

test_that("write_sharing works when store_ec=TRUE", {
  r    <- .make_toy_result()
  tmpf <- tempfile(fileext=".tsv")
  write_sharing(r, tmpf, min_sharing=0.0)
  expect_true(file.exists(tmpf))
  unlink(tmpf)
})

# ============================================================
# Downstream compatibility
# ============================================================

test_that("as_count_matrix returns matrix", {
  r   <- .make_toy_result()
  mat <- as_count_matrix(r)
  expect_true(is.matrix(mat) || methods::is(mat, "Matrix"))
})

test_that("as_tpm_matrix returns matrix with positive values", {
  r   <- .make_toy_result()
  mat <- as_tpm_matrix(r)
  expect_true(all(mat@x >= 0))
})

test_that("as_sce returns SCE when available", {
  if (!requireNamespace("SingleCellExperiment", quietly=TRUE)) skip("SCE not available")
  r   <- .make_toy_result()
  sce <- as_sce(r)
  expect_true(methods::is(sce, "SingleCellExperiment"))
})

test_that("get_counts returns count matrix", {
  r   <- .make_toy_result()
  mat <- get_counts(r)
  expect_true(methods::is(mat, "Matrix") || is.matrix(mat))
})

test_that("get_ec returns SpotECCounts when store_ec=TRUE", {
  r  <- .make_toy_result()
  ec <- get_ec(r)
  expect_false(is.null(ec))
  expect_s3_class(ec, "SpotECCounts")
})

test_that("get_sharing returns data.table of pairs", {
  r  <- .make_toy_result()
  sh <- get_sharing(r, min_sharing=0.0)
  expect_true(is.data.frame(sh) || data.table::is.data.table(sh))
})

# ============================================================
# Input validation
# ============================================================

test_that("run_spatialTE errors when no BAM provided", {
  expect_error(
    run_spatialTE(gtf_file=gtf_f, lr_input=lr_f, sample_id="s1"),
    "pe_bam.*se_bam"
  )
})

test_that("run_spatialTE errors when no sample_id provided", {
  expect_error(
    run_spatialTE(pe_bam="x.bam", gtf_file=gtf_f, lr_input=lr_f),
    "sample_id"
  )
})

test_that("run_spatialTE errors when both sample_id and sample_id_file given", {
  expect_error(
    run_spatialTE(pe_bam="x.bam", gtf_file=gtf_f, lr_input=lr_f,
                  sample_id="s1", sample_id_file="f.txt"),
    "only one"
  )
})

test_that("run_scTE errors when no BAM provided", {
  expect_error(
    run_scTE(gtf_file=gtf_f, lr_input=lr_f, sample_id="s1"),
    "pe_bam.*se_bam"
  )
})

# ============================================================
# Utility functions
# ============================================================

test_that(".is_novel_tx correctly classifies IDs", {
  ids    <- c("ENST001","NM_001","novel.1","transcript_001","NR_001","XM_001")
  result <- SpatialTE:::.is_novel_tx(ids)
  expect_equal(result, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that(".tx_length_bin returns correct bins", {
  expect_equal(SpatialTE:::.tx_length_bin(300L),  "short")
  expect_equal(SpatialTE:::.tx_length_bin(1000L), "medium")
  expect_equal(SpatialTE:::.tx_length_bin(2000L), "long")
  expect_equal(SpatialTE:::.tx_length_bin(5000L), "vlong")
})

test_that(".read_file handles .gz", {
  dt <- SpatialTE:::.read_file(lr_gz_f)
  expect_true(nrow(dt) > 0)
})

test_that(".write_file writes and can be read back", {
  tmpf <- tempfile(fileext=".tsv")
  dt   <- data.table(a=1:3, b=letters[1:3])
  SpatialTE:::.write_file(dt, tmpf)
  dt2  <- fread(tmpf)
  expect_equal(nrow(dt2), 3L)
  unlink(tmpf)
})

test_that(".write_file writes .gz", {
  tmpf <- tempfile(fileext=".tsv.gz")
  dt   <- data.table(x=1:5)
  SpatialTE:::.write_file(dt, tmpf)
  dt2  <- fread(tmpf)
  expect_equal(nrow(dt2), 5L)
  unlink(tmpf)
})

# ============================================================
# Toy data files existence
# ============================================================

test_that("all expected toy data files exist", {
  expected <- c(
    "toy_gtf.gtf", "toy_gtf.gtf.gz",
    "toy_lr_counts.tsv", "toy_lr_counts.tsv.gz",
    "toy_lr_counts_multi.tsv",
    "toy_coords.csv", "toy_coords_visium.csv",
    "toy_ec_raw.tsv", "toy_ec_raw.tsv.gz",
    "toy_coverage_bias.tsv", "toy_frag_length.tsv",
    "toy_sample_ids.txt", "toy_sample_ids.txt.gz"
  )
  for (f in expected) {
    path <- system.file("extdata", f, package="SpatialTE")
    expect_true(nchar(path) > 0 && file.exists(path),
                label=paste("missing:", f))
  }
})

# ============================================================
# Visium coordinate loading
# ============================================================

test_that(".load_coords auto-detects Visium format", {
  vis_f <- system.file("extdata","toy_coords_visium.csv", package="SpatialTE")
  dt    <- SpatialTE:::.load_coords(vis_f, format="auto")
  expect_true(all(c("spot_id","x","y") %in% names(dt)))
  expect_equal(nrow(dt), 16L)
})

test_that(".load_coords explicit visium format", {
  vis_f <- system.file("extdata","toy_coords_visium.csv", package="SpatialTE")
  dt    <- SpatialTE:::.load_coords(vis_f, format="visium")
  expect_equal(nrow(dt), 16L)
  expect_true(all(dt$x > 0))
})

# ============================================================
# EC raw data loading and building
# ============================================================

test_that("toy_ec_raw can be loaded and built into SpotECCounts", {
  ec_raw_f <- system.file("extdata","toy_ec_raw.tsv", package="SpatialTE")
  ec_raw   <- fread(ec_raw_f)
  expect_true(all(c("spot_key","t_idx","nh") %in% names(ec_raw)))
  expect_true(nrow(ec_raw) > 0)
  ec_raw[, sample_id := "sA"]
  tx_names <- c(paste0("ENST",sprintf("%011d",1:8)),"novel.1","novel.2")
  res <- SpatialTE:::.build_ec_one_sample(ec_raw, tx_names, "sA",
                                           "barcode_umi", 200L)
  expect_s3_class(res, "SpotECCounts")
  expect_true(length(res$spot_ids) == 9L)
})

test_that("toy_ec_raw.gz loads correctly", {
  ec_raw_gz_f <- system.file("extdata","toy_ec_raw.tsv.gz", package="SpatialTE")
  ec_raw <- fread(ec_raw_gz_f)
  expect_true(nrow(ec_raw) > 0)
})

# ============================================================
# Coverage bias and fragment length profiles
# ============================================================

test_that("toy_coverage_bias has 4 bin columns", {
  cov_f <- system.file("extdata","toy_coverage_bias.tsv", package="SpatialTE")
  dt    <- fread(cov_f)
  expect_true(all(c("bin","short","medium","long","vlong") %in% names(dt)))
  expect_equal(nrow(dt), 100L)
  # g(x) should have mean ~1 per bin
  for (col in c("short","medium","long","vlong")) {
    expect_equal(mean(dt[[col]]), 1.0, tolerance=0.01,
                 label=paste("mean of", col))
  }
})

test_that("toy_frag_length has expected columns", {
  frag_f <- system.file("extdata","toy_frag_length.tsv", package="SpatialTE")
  dt     <- fread(frag_f)
  expect_true(all(c("read_length","count") %in% names(dt)))
  expect_true(nrow(dt) > 0)
  expect_true(all(dt$count >= 0))
  # Peak should be around 100-110 bp
  peak_len <- dt$read_length[which.max(dt$count)]
  expect_true(peak_len >= 95 && peak_len <= 115)
})

# ============================================================
# Multi-sample LR prior
# ============================================================

test_that("multi-sample LR prior loads correctly", {
  lr_multi_f <- system.file("extdata","toy_lr_counts_multi.tsv", package="SpatialTE")
  gtf        <- SpatialTE:::.parse_gtf(gtf_f)
  # Use first sample only
  lr_multi   <- fread(lr_multi_f)[sample_id == "LR_sampleA"]
  pr         <- SpatialTE:::.load_lr_prior(lr_multi, gtf)
  expect_equal(sum(pr$pi_lr), 1, tolerance=1e-6)
})

# ============================================================
# Sample ID file
# ============================================================

test_that("toy_sample_ids.txt has 2 columns", {
  sid_f <- system.file("extdata","toy_sample_ids.txt", package="SpatialTE")
  dt    <- fread(sid_f, header=FALSE, col.names=c("read_id","sample_id"))
  expect_equal(ncol(dt), 2L)
  expect_true(all(c("sample_A","sample_B","sample_C") %in% unique(dt$sample_id)))
})

test_that("toy_sample_ids.txt.gz loads correctly", {
  sid_gz_f <- system.file("extdata","toy_sample_ids.txt.gz", package="SpatialTE")
  dt       <- fread(sid_gz_f, header=FALSE)
  expect_equal(ncol(dt), 2L)
  expect_true(nrow(dt) > 0)
})

# ============================================================
# End-to-end pipeline using toy EC data (no BAM needed)
# ============================================================

test_that("full pipeline runs with toy EC data", {
  gtf      <- SpatialTE:::.parse_gtf(gtf_f)
  lr_pr    <- SpatialTE:::.load_lr_prior(lr_f, gtf)
  ec_raw_f <- system.file("extdata","toy_ec_raw.tsv", package="SpatialTE")
  ec_raw   <- fread(ec_raw_f)
  ec_raw[, sample_id := "sA"]

  # Build EC
  ec_c <- SpatialTE:::.build_ec_one_sample(
    ec_raw, c(paste0("ENST",sprintf("%011d",1:8)),"novel.1","novel.2"),
    "sA", "barcode_umi", 200L)
  expect_s3_class(ec_c, "SpotECCounts")

  # Compute eff_len using toy coverage profile
  cov_f  <- system.file("extdata","toy_coverage_bias.tsv", package="SpatialTE")
  frag_f <- system.file("extdata","toy_frag_length.tsv", package="SpatialTE")
  cov_dt  <- fread(cov_f)
  frag_dt <- fread(frag_f)
  g_cov   <- list(short=cov_dt$short, medium=cov_dt$medium,
                   long=cov_dt$long,   vlong=cov_dt$vlong)
  frag_tab <- setNames(frag_dt$count, as.character(frag_dt$read_length))
  el <- SpatialTE:::.compute_eff_len(gtf, g_cov, frag_tab)
  expect_true(all(el >= 50))

  # Build prior + run EM
  pr   <- SpatialTE:::.build_prior(lr_pr, ec_c, el, gamma_base=10.0)
  em_r <- SpatialTE:::.run_em_all_spots(ec_c, pr, el,
                                         max_iter=100L, tol=1e-4, n_cores=1L)
  expect_true(all(em_r$count_mat@x >= 0))

  # Build result
  r <- SpatialTE:::.build_result(
    em_r, ec_c, pr, gtf, el, NULL, lr_pr, 10.0,
    list(sample_id="sA", ec_data=ec_c), 50)
  expect_s3_class(r, "SpatialTEResult")
  expect_output(print(r), "SpatialTEResult")

  # Write outputs
  tmpdir <- tempfile()
  write_spatialTE(r, tmpdir, compress=FALSE, write_sharing=TRUE)
  expect_true(file.exists(file.path(tmpdir,"counts.tsv")))
  expect_true(file.exists(file.path(tmpdir,"count_matrix.tsv")))
  expect_true(file.exists(file.path(tmpdir,"tpm_matrix.tsv")))
  expect_true(file.exists(file.path(tmpdir,"efflen_table.tsv")))
  expect_true(file.exists(file.path(tmpdir,"qc_summary.tsv")))

  # Verify counts file content
  counts_dt <- fread(file.path(tmpdir,"counts.tsv"))
  expect_true(all(c("spot_id","transcript_id","em_count","tpm") %in% names(counts_dt)))
  expect_true(all(counts_dt$em_count >= 0))
  expect_true(all(counts_dt$tpm >= 0))

  unlink(tmpdir, recursive=TRUE)
})

test_that("spatial pipeline with toy EC data and coordinates", {
  gtf      <- SpatialTE:::.parse_gtf(gtf_f)
  lr_pr    <- SpatialTE:::.load_lr_prior(lr_f, gtf)
  ec_raw_f <- system.file("extdata","toy_ec_raw.tsv", package="SpatialTE")
  ec_raw   <- fread(ec_raw_f); ec_raw[, sample_id := "sA"]
  tx_names <- c(paste0("ENST",sprintf("%011d",1:8)),"novel.1","novel.2")
  ec_c     <- SpatialTE:::.build_ec_one_sample(ec_raw, tx_names,
                                                "sA","barcode_umi",200L)
  cov_f    <- system.file("extdata","toy_coverage_bias.tsv",package="SpatialTE")
  frag_f   <- system.file("extdata","toy_frag_length.tsv",package="SpatialTE")
  cov_dt   <- fread(cov_f); frag_dt <- fread(frag_f)
  g_cov    <- list(short=cov_dt$short,medium=cov_dt$medium,
                    long=cov_dt$long,vlong=cov_dt$vlong)
  frag_tab <- setNames(frag_dt$count,as.character(frag_dt$read_length))
  el       <- SpatialTE:::.compute_eff_len(gtf,g_cov,frag_tab)

  # Load coordinates and build neighbor graph
  coords_f2 <- system.file("extdata","toy_coords.csv",package="SpatialTE")
  coords    <- SpatialTE:::.load_coords(coords_f2)
  ng        <- SpatialTE:::.build_neighbor_graph(coords)
  expect_true(!is.null(ng$neighbors))

  # Build spatial prior
  pr <- SpatialTE:::.build_prior(lr_pr,ec_c,el,gamma_base=8.0,
                                  neighbor_graph=ng,alpha_local=0.7)
  expect_true(!is.null(pr$alpha_mat))

  # Run EM and build result with coords
  em_r <- SpatialTE:::.run_em_all_spots(ec_c,pr,el,max_iter=50L,tol=1e-3)
  r    <- SpatialTE:::.build_result(em_r,ec_c,pr,gtf,el,coords,lr_pr,8.0,
                                     list(sample_id="sA",ec_data=ec_c),50)
  cd   <- SpatialTE:::.get_coldata(r)
  expect_true(all(c("x","y") %in% names(cd)))
})
