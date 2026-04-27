#!/usr/bin/env Rscript
# Generates all toy/example datasets for SpatialTE
suppressPackageStartupMessages({ library(data.table) })
set.seed(42)
outdir <- "inst/extdata"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
cat("Generating toy data...\n")

tx_ids   <- c(paste0("ENST",sprintf("%011d",1:8)),"novel.1","novel.2")
gene_ids <- c(rep("ENSG001",3),rep("ENSG002",3),rep("ENSG003",2),"ENSG004","ENSG004")
tx_lens  <- c(1200L,800L,600L,2000L,1500L,900L,3500L,400L,1100L,700L)
true_props <- c(0.25,0.12,0.06,0.18,0.09,0.05,0.10,0.07,0.05,0.03)
true_props <- true_props/sum(true_props)
cert_vals  <- c(0.90,0.85,0.80,0.88,0.82,0.78,0.72,0.91,0.55,0.45)
