suppressPackageStartupMessages({
  library(data.table)
})
set.seed(42)

outdir <- "/home/claude/SpatialTE/inst/extdata"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# ── Transcripts ───────────────────────────────────────────────
tx_ids   <- c(paste0("ENST", sprintf("%011d", 1:8)), "novel.1", "novel.2")
gene_ids <- c(rep("ENSG001",3), rep("ENSG002",3), rep("ENSG003",2), "ENSG004","ENSG004")
tx_lens  <- c(1200L,800L,600L,2000L,1500L,900L,3500L,400L,1100L,700L)
true_props <- c(0.25,0.12,0.06,0.18,0.09,0.05,0.10,0.07,0.05,0.03)
true_props <- true_props / sum(true_props)

# ── GTF ──────────────────────────────────────────────────────
gtf_lines <- c("##format: GTF")
for (i in seq_along(tx_ids)) {
  is_nov   <- grepl("^novel", tx_ids[i])
  tx_type  <- if (is_nov) "novel_not_in_catalog" else "protein_coding"
  start    <- 1000L * i
  e1_len   <- as.integer(tx_lens[i] * 0.6)
  e2_len   <- tx_lens[i] - e1_len
  e1_end   <- start + e1_len - 1L
  e2_start <- e1_end + 101L
  e2_end   <- e2_start + e2_len - 1L
  tx_end   <- e2_end
  attr_tx  <- sprintf(
    'transcript_id "%s"; gene_id "%s"; transcript_type "%s";',
    tx_ids[i], gene_ids[i], tx_type)
  attr_ex  <- sprintf('transcript_id "%s"; gene_id "%s";', tx_ids[i], gene_ids[i])
  gtf_lines <- c(gtf_lines,
    sprintf("chr1\tSpatialTE_toy\ttranscript\t%d\t%d\t.\t+\t.\t%s",
            start, tx_end, attr_tx),
    sprintf("chr1\tSpatialTE_toy\texon\t%d\t%d\t.\t+\t.\t%s",
            start, e1_end, attr_ex),
    sprintf("chr1\tSpatialTE_toy\texon\t%d\t%d\t.\t+\t.\t%s",
            e2_start, e2_end, attr_ex))
}
writeLines(gtf_lines, file.path(outdir, "toy_gtf.gtf"))
con <- gzfile(file.path(outdir, "toy_gtf.gtf.gz"), "wt")
writeLines(gtf_lines, con); close(con)

# ── LR prior ─────────────────────────────────────────────────
cert <- c(0.90,0.85,0.80,0.88,0.82,0.78,0.72,0.91,0.55,0.45)
lr_dt <- data.table(
  transcript_id     = tx_ids,
  sample_id         = "LR_sample",
  em_count          = round(true_props * 50000, 1),
  unique_count      = round(true_props * 35000, 0),
  total_count       = round(true_props * 50000, 0),
  certainty         = cert,
  multimapping_rate = 1 - cert,
  tpm               = round(true_props * 1e6, 1),
  gene_id           = gene_ids,
  is_novel          = grepl("^novel", tx_ids),
  tx_length         = tx_lens
)
fwrite(lr_dt, file.path(outdir, "toy_lr_counts.tsv"), sep="\t")
fwrite(lr_dt, file.path(outdir, "toy_lr_counts.tsv.gz"), sep="\t", compress="gzip")

# ── Spot coordinates (Patho-DBiT generic format) ─────────────
n_spots <- 9L  # 3x3 grid
spots   <- paste0("SPOT_", sprintf("%03d", seq_len(n_spots)))
coords  <- data.table(
  spot_id = spots,
  x       = rep(c(100,200,300), 3),
  y       = rep(c(100,200,300), each=3)
)
fwrite(coords, file.path(outdir, "toy_coords.csv"), sep=",")

cat("Toy data generated successfully.\n")
cat("  toy_gtf.gtf[.gz]\n")
cat("  toy_lr_counts.tsv[.gz]\n")
cat("  toy_coords.csv\n")
