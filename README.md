# SpatialTE

**Isoform-level quantification for FFPE spatial transcriptomics and single-cell RNA-seq**

[![R-CMD-check](https://github.com/hmutpw/SpatialTE/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/hmutpw/SpatialTE/actions/workflows/R-CMD-check.yml)

SpatialTE performs transcript-level quantification in individual spots or cells of spatial transcriptomics and scRNA-seq data. It integrates paired bulk long-read RNA-seq as a Dirichlet prior to guide isoform assignment in data-sparse conditions (FFPE degradation, shallow coverage, multi-mapping reads). SpatialTE is the companion tool to [IsoEM](https://github.com/hmutpw/IsoEM), which handles the long-read quantification step.

---

## Why SpatialTE?

Standard spatial/scRNA-seq workflows quantify at the **gene level**. SpatialTE goes further:

- **Isoform resolution**: distinguishes transcripts from the same gene, including novel IsoQuant-assembled isoforms
- **FFPE correction**: estimates FFPE-specific effective lengths from 5'→3' coverage bias, correcting for RNA degradation
- **Long-read prior**: uses bulk long-read quantification (from IsoEM) as a Bayesian prior to stabilise estimates in data-sparse spots
- **Spatial smoothing**: incorporates neighbouring spot information to improve reliability in low-UMI regions
- **TE-aware**: correctly handles highly repetitive transposable element loci through equivalence-class compression

---

## Workflow overview

```
Bulk long-read RNA-seq          Spatial / scRNA-seq
    (IsoQuant + IsoEM)              (STAR → toTranscriptome BAM)
           │                                │
           ▼                                ▼
    LR prior (pi_k, certainty)    BAM streaming (one pass):
                                    • f(d): fragment length dist.
                                    • g(x): coverage bias curves
                                    • EC: equivalence classes
           │                                │
           └──────────────┬─────────────────┘
                          ▼
              Dirichlet prior construction
              (gamma auto-calibration, spatial smoothing)
                          │
                          ▼
                  MAP-EM per spot/cell
              (vectorised, BiocParallel)
                          │
                          ▼
              SpatialTEResult
          (SCE/SPE compatible object)
```

---

## Installation

```r
# Install from GitHub
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("hmutpw/SpatialTE")

# Install a specific version
remotes::install_github("hmutpw/SpatialTE@v0.1.0")
```

**Bioconductor dependencies** are installed automatically:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Rsamtools", "BiocParallel",
                        "SingleCellExperiment", "SummarizedExperiment"))
```

**Optional but recommended:**

```r
BiocManager::install("SpatialExperiment")   # for spatial objects
BiocManager::install("rtracklayer")          # faster GTF parsing
```

---

## Quick start

### Spatial transcriptomics (Patho-DBiT / Visium)

```r
library(SpatialTE)

result <- run_spatialTE(
  pe_bam      = "spatial_pe.Aligned.toTranscriptome.out.bam",
  gtf_file    = "transcript_models.gtf",
  lr_input    = "isoquant_em_counts.tsv.gz",   # from IsoEM
  spot_coords = "tissue_positions.csv",          # spot coordinates
  sample_id   = "sample_A"
)

print(result)
# SpatialTEResult
#   Mode       : spatial
#   Transcripts: 98421
#   Spots/cells: 4823
#   Converged  : 4823 / 4823 spots
#   Median unique frac: 0.612

write_spatialTE(result, outdir = "results/sample_A/")
```

### Single-cell RNA-seq

```r
result <- run_scTE(
  pe_bam    = "sc_pe.Aligned.toTranscriptome.out.bam",
  gtf_file  = "transcript_models.gtf",
  lr_input  = "isoquant_em_counts.tsv.gz",
  sample_id = "sample_A"
)
```

### With IsoEM result object

```r
library(IsoEM)
library(SpatialTE)

# Run IsoEM on long-read data
lr_result <- run_isoem(
  counts_file = "transcript_model_counts.tsv.gz",
  gtf_file    = "transcript_models.gtf",
  sample_id   = "bulk_LR"
)

# Pass directly to SpatialTE
result <- run_spatialTE(
  pe_bam    = "spatial.bam",
  gtf_file  = "transcript_models.gtf",
  lr_input  = lr_result,          # IsoEMResult object
  spot_coords = "coords.csv",
  sample_id   = "sample_A"
)
```

---

## Input files

### BAM files

SpatialTE requires **STAR toTranscriptome BAM** files:

```bash
STAR \
  --runMode alignReads \
  --genomeDir /path/to/star_index \
  --readFilesIn R2.fastq.gz \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM \
  --outSAMmultNmax 200 \          # keep multi-mapping reads
  --soloType CB_UMI_Simple \
  --soloCBwhitelist barcodes.txt \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12
```

- Provide `pe_bam`, `se_bam`, or both (at least one required)
- `--outSAMmultNmax 200` is important to retain multi-mapping reads for TE loci

### GTF file

The IsoQuant `transcript_models.gtf` from the same long-read run used to generate the LR prior. This ensures transcript IDs match between inputs.

### Long-read prior

Accepts three formats:

| Format | Description |
|--------|-------------|
| `IsoEMResult` | Output of `IsoEM::run_isoem()` |
| `data.frame` | Must contain `transcript_id` and `em_count` columns; `certainty` optional |
| File path | `.tsv` or `.tsv.gz` with same columns as above |

### Spot coordinates

| Platform | File | Format |
|----------|------|--------|
| Visium | `tissue_positions_list.csv` | `barcode, in_tissue, array_row, array_col, pxl_row, pxl_col` |
| Visium HD | `tissue_positions.parquet` or `.csv` | Same as Visium |
| Patho-DBiT | User-provided | `spot_id, x, y` (generic) |
| Other | User-provided | `spot_id, x, y` (generic) |

### Barcode/UMI extraction

Three methods, tried in priority order:

```r
# Method 1: Provide a mapping table (most flexible)
run_spatialTE(..., bc_umi_table = "read_id_bc_umi.txt")
# File format: 3 columns, no header: read_id | barcode | umi

# Method 2: Regex pattern from read ID
run_spatialTE(..., bc_umi_pattern = "^([A-Z]{16})([A-Z]{12})")

# Method 3: BAM tags (default, auto-detected)
run_spatialTE(..., cb_tag = "CB", ub_tag = "UB")
```

---

## Key parameters

### `run_spatialTE()` / `run_scTE()`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pe_bam` | NULL | PE toTranscriptome BAM |
| `se_bam` | NULL | SE toTranscriptome BAM (at least one required) |
| `gtf_file` | — | IsoQuant GTF[.gz] |
| `lr_input` | — | Long-read prior (IsoEMResult / data.frame / file) |
| `spot_coords` | NULL | Spatial coordinates (`run_spatialTE` only) |
| `coords_format` | "auto" | "visium" / "visium_hd" / "patho_dbit" / "generic" |
| `sample_id` | NULL | Sample name — single-sample mode |
| `sample_id_file` | NULL | Two-column file: read_id \| sample_id — multi-sample mode |
| `ec_unit` | "barcode_umi" | EC grouping unit: "barcode_umi" or "read" |
| `max_ec_size` | 200 | ECs > this use LR prior directly (no EM iteration) |
| `frag_dist_source` | "pe" | Source for f(d): "pe" / "se" / "both" |
| `gamma_base` | NULL | Prior strength (NULL = auto-calibrate) |
| `lambda` | 10 | Spot prior smoothing parameter |
| `certainty_threshold` | 0.1 | LR/SR ratio guard for low-certainty transcripts |
| `alpha_local` | 0.7 | Own vs neighbour weight in spatial prior |
| `neighbor_radius` | NULL | NULL = 1.5 × median NN distance |
| `max_iter` | 200 | EM max iterations |
| `tol` | 1e-6 | EM convergence tolerance |
| `min_mapq` | 0 | Minimum BAM MAPQ filter |
| `cb_whitelist` | NULL | Valid barcode whitelist file or vector |
| `enforce_bulk_consistency` | FALSE | Post-EM bulk correction |
| `store_ec` | TRUE | Retain EC data for `write_sharing()` |
| `n_cores` | 1 | Parallel workers (Linux: fork; Windows: SOCK) |

---

## Output files

```
write_spatialTE(result, outdir = "results/")
```

| File | Description |
|------|-------------|
| `counts.tsv[.gz]` | Long-format quantification (one row per spot × transcript) |
| `count_matrix.tsv[.gz]` | Wide-format: transcripts × spots |
| `tpm_matrix.tsv[.gz]` | Wide-format TPM matrix |
| `efflen_table.tsv[.gz]` | Per-transcript effective lengths |
| `sharing_table.tsv[.gz]` | Transcript pairs sharing reads |
| `qc_summary.tsv[.gz]` | Per-spot QC statistics |

### counts.tsv columns

| Column | Description |
|--------|-------------|
| `spot_id` | Spot/cell barcode |
| `transcript_id` | Transcript identifier |
| `sample_id` | Sample name |
| `gene_id` | Gene identifier (from GTF) |
| `is_novel` | IsoQuant-assembled novel transcript |
| `tx_length` | Transcript length (exon sum, bp) |
| `em_count` | EM-estimated count |
| `unique_count` | Reads from unique-mapping ECs |
| `multi_count` | Reads from multi-mapping ECs |
| `total_count` | Total reads (unique + multi) |
| `eff_len` | FFPE-corrected effective length (bp) |
| `tpm` | Transcripts Per Million |
| `certainty` | `unique_count / total_count` |
| `multimapping_rate` | `1 - certainty` |
| `confidence` | Composite reliability score (0–1) |
| `lr_certainty` | LR unique assignment rate |
| `low_eff_len` | TRUE if eff_len near floor |
| `bulk_consistency_ratio` | Spot-aggregate vs bulk LR ratio |

---

## Multi-sample mode

Multi-sample mode is for data where **multiple samples were jointly processed through STAR and IsoQuant** (sharing the same transcript models).

```r
# Step 1: generate sample ID file from BAM or FASTQ
# (use the helpers from the IsoEM package)
IsoEM::extract_readids_from_bam(
  sample_ids  = c("sample_A", "sample_B", "sample_C"),
  bam_files   = c("sA.bam", "sB.bam", "sC.bam"),
  output_file = "sample_readids.txt"
)

# Step 2: run SpatialTE with merged BAM
result <- run_spatialTE(
  pe_bam         = "merged.Aligned.toTranscriptome.out.bam",
  gtf_file       = "transcript_models.gtf",
  lr_input       = "isoquant_em_counts.tsv.gz",
  spot_coords    = "coords.csv",
  sample_id_file = "sample_readids.txt",    # maps reads to samples
  n_cores        = 8
)
# Returns a SpatialTEDataset (list of SpatialTEResult)
```

---

## Downstream analysis

### Seurat

```r
library(Seurat)
sobj <- as_seurat(result)

# SpatialTE counts are isoform-level — aggregate to gene level if needed
gene_counts <- rowsum(as.matrix(as_count_matrix(result)),
                       rowData(as_sce(result))$gene_id)
```

### SingleCellExperiment

```r
sce <- as_sce(result)
# Access SpatialTE assays
counts(sce)      # EM counts
assay(sce,"tpm") # TPM

# Filter by confidence
high_conf <- assay(sce, "confidence") > 0.5
```

### DESeq2 / edgeR

```r
library(DESeq2)
mat <- as_count_matrix(result, round_counts = TRUE)
cd  <- data.frame(condition = c("tumor","normal"),
                   row.names = colnames(mat))
dds <- DESeqDataSetFromMatrix(mat, cd, ~condition)
```

### QC visualisation

```r
# Spatial map of unique assignment rate
plot_qc_spatial(result, metric = "unique_frac")

# Certainty distribution across transcripts
plot_certainty(result)

# EC size distribution
plot_ec_size(result)
```

---

## Understanding confidence scores

```
confidence_k_s = unique_frac_k_s × depth_factor_s × lr_certainty_k

unique_frac_k_s : fraction of reads that uniquely map to transcript k in spot s
depth_factor_s  : min(N_s / N_10pct, 1), where N_10pct is 10th percentile of spot UMIs
lr_certainty_k  : fraction of long-reads uniquely assigned to transcript k
```

| Score | Interpretation |
|-------|----------------|
| > 0.7 | High confidence |
| 0.4–0.7 | Moderate — some ambiguity |
| < 0.4 | Low — use with caution |
| 0 (low_eff_len) | Effective length near floor — unreliable |

---

## Transcript sharing table

The sharing table identifies isoform pairs that cannot be reliably distinguished in short-read data:

```r
sh <- get_sharing(result, min_sharing = 0.5)
head(sh)
#    transcript_1      transcript_2  shared_reads  sharing_fraction  recommendation
# 1  ENST00000001     ENST00000002         1842             0.91    high_merge_recommended
# 2  ENST00000003     novel.1               541             0.63    moderate_caution
```

Transcripts with `recommendation == "high_merge_recommended"` (sharing > 80%) should be summed for differential expression analysis.

---

## Bulk consistency QC

After EM, SpatialTE computes a per-transcript bulk consistency ratio:

```
bulk_consistency_ratio_k = (Σ_s rho_k_s × N_s / Σ_s N_s) / pi_k_LR
```

- **~1.0**: spatial estimate agrees with bulk long-read
- **> 2**: transcript is enriched in specific spots (expected for cell-type-specific isoforms)
- **< 0.5**: systematic under-estimation (possible degradation or cell-type mismatch)

```r
rd <- rowData(as_sce(result))
# Flag transcripts with suspicious bulk consistency
flagged <- rd[abs(log2(rd$bulk_consistency_ratio)) > 2, ]
```

---

## Frequently asked questions

**Q: Do I need IsoEM to use SpatialTE?**

No. SpatialTE accepts any data frame with `transcript_id` and `em_count` columns. IsoEM is the recommended upstream tool but is not required.

**Q: My spatial data is not FFPE — does FFPE correction matter?**

The FFPE correction (`g(x)` and `f(d)`) is data-driven. If your sample has uniform coverage, `g(x)` will be approximately flat and the correction will have no effect. FFPE correction is always safe to apply.

**Q: What is `--outSAMmultNmax 200` and why do I need it?**

By default STAR reports only 10 multi-mapping alignments per read. Many TE loci have hundreds of copies in the genome. `--outSAMmultNmax 200` increases this limit so that reads from repetitive elements are properly represented in equivalence classes.

**Q: What's the difference between `run_spatialTE` and `run_scTE`?**

`run_spatialTE` uses spot coordinates to build a spatial neighbour graph, which improves prior estimation by incorporating information from neighbouring spots. `run_scTE` works the same way without coordinates — suitable for dissociated single-cell data.

**Q: Can I use STARsolo output directly?**

Yes. STARsolo with `--quantMode TranscriptomeSAM` produces a `Aligned.toTranscriptome.out.bam` file that SpatialTE accepts directly. The CB and UB tags are automatically detected.

**Q: How do I handle very low UMI spots?**

Low-UMI spots are handled automatically. `gamma(k,s)` is inversely proportional to spot UMI depth: sparse spots get stronger LR prior influence, which stabilises estimates. Spots with zero UMIs return zero counts.

---

## Citation

If you use SpatialTE, please cite:

> Tan P. (2025). *SpatialTE: Isoform-level quantification for FFPE spatial transcriptomics*. https://github.com/hmutpw/SpatialTE

---

## License

MIT © Puwen Tan
