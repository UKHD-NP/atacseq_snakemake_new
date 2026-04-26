# ATAC-seq Snakemake Pipeline

A modular Snakemake pipeline for paired-end ATAC-seq from raw FASTQ to peaks, signal tracks, and QC reports.

This repository currently uses:
- `workflow/Snakefile`
- `config/config.yml`
- per-rule Conda environments in `workflow/envs/`

> **WARNING:** This pipeline covers **upstream data analysis only** (QC → alignment → peak calling → QC per-sample counts).
> The per-sample count output (`featurecounts/{sample}.readCountInPeaks.txt`) is **not** ready for differential accessibility (DA) analysis with DESeq2 or edgeR.
> For downstream differential accessibility analysis you must:
> 1. Generate a **consensus peak set** across all samples.
> 2. Re-quantify reads against the consensus peaks to create a **unified count matrix**.
> 3. Run differential accessibility analysis on that count matrix.

## Table of Contents

- [Pipeline Summary](#pipeline-summary)
- [Workflow DAG](#workflow-dag)
- [Rule Quick Reference](#rule-quick-reference)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Files](#input-files)
  - [1. Samplesheet](#1-samplesheet)
  - [2. Reference](#2-reference)
- [Local Run](#local-run)
- [Running on HPC with LSF](#running-on-hpc-with-lsf)
  - [Node roles](#node-roles)
  - [Step 1 — Set up Snakemake environment](#step-1---set-up-snakemake-environment)
  - [Step 2 — Clone the pipeline](#step-2---clone-the-pipeline)
  - [Step 3 — Edit configuration](#step-3---edit-configuration)
  - [Step 4 — Update conda-prefix](#step-4---update-conda-prefix-in-the-lsf-profile)
  - [Step 5 — Validate with a dry-run](#step-5---validate-with-a-dry-run)
  - [Step 6 — Submit to HPC](#step-6---submit-to-hpc)
  - [Monitoring jobs](#monitoring-jobs)
- [Key Configuration](#key-configuration)
  - [Why these default params were chosen](#why-these-default-params-were-chosen)
  - [How to choose macs3_gsize](#how-to-choose-macs3_gsize)
- [BAM Filtering Criteria](#re-mark-duplicates-and-bam-filtering-criteria)
  - [1) SAMtools core filter](#1-samtools-core-filter-from-configyml)
  - [2) BAMTools extra filter](#2-bamtools-extra-filter-optional-but-enabled-when-bamtools-exists)
  - [3) Pysam pair/orphan cleanup](#3-pysam-pairorphan-cleanup-optional-but-enabled-when-python3pysam-exists)
- [Main Outputs](#main-outputs)
- [QC Visual Guide](#qc-visual-guide-from-resources)
- [ENCODE QC Benchmarks](#encode-qc-benchmarks)
- [Space-saving behavior](#space-saving-behavior)
- [Troubleshooting](#troubleshooting)
- [Acknowledgments](#acknowledgments)
- [References](#references)
- [License](#license)

---

## Pipeline Summary

Main workflow (per sample):
1. Merge lanes by `sample_id` (from samplesheet)
2. Raw FastQC
3. Trimming (`trim_galore` default, or `fastp`)
4. Alignment (`bowtie2` default, or `bwa` = BWA-MEM2)
5. Re-mark duplicates (Picard MarkDuplicates; optional by config)
6. BAM filtering to produce `*.filtered.bam`
   - always applies SAM flag/MAPQ filtering
   - optionally applies blacklist exclusion
   - optionally excludes mitochondrial reads via `ref.keep_mito`
   - can preserve the source BAM before filtering
7. BAM stats on filtered BAM (`samtools stats/flagstat/idxstats`)
8. Signal tracks
   - Scaled bedGraph + bigWig from filtered BAM
   - ATAC-shifted BAM + shifted RPGC bigWig
9. Peak calling (MACS3 with Tn5-shifted BED)
10. Peak QC summary plots (`plot_macs_qc.r`)
11. FRiP (two methods: bedtools intersect + featureCounts log)
12. Peak annotation (HOMER + summary)
13. deepTools matrix/profile/heatmap/fingerprint/bamPEFragmentSize
14. NFR analysis — nucleosome-free vs mononucleosomal bigWigs + TSS profile/heatmap
15. ataqv JSON + mkarv HTML report + TSS enrichment / NFR metrics table for MultiQC
16. ATACseqQC — PT score, NFR score, TSSE score + QC plots (narrow peaks only)
17. MultiQC
18. Cleanup temporary/intermediate FASTQ files

## Workflow DAG

Pipeline flow (per sample):

```text
samplesheet + reference prep
        |
        v
merge_raw_fastqs
        |
        +--> fastqc_raw
        |
        v
trimming (fastp | trim_galore)
        |
        +--> fastqc_trimmed (trim_galore mode)
        |
        v
alignment (bwa-mem2 | bowtie2)
        |
        v
    sort_bam
        |
        v
mark_duplicates (optional)
        |
        v
bam_filter  --->  filtered.bam / filtered.bam.bai
        |                    |
        |                    +--> align_stats (samtools stats/flagstat/idxstats)
        |                    +--> bedtools_genomecov -> bedGraphToBigWig -> bigWig
        |                    +--> shift_bam (alignmentSieve --ATACshift) -> shifted.bam + shifted.bigWig
        |                    +--> macs3_callpeak_tn5 (narrow: filtered.bam -> bamtobed -> awk Tn5 shift -> MACS3 BED mode)
        |                               |            (broad:  filtered.bam -> MACS3 BAMPE mode)
        |                               |    
        |                               +--> frip_score (filtered.bam + peaks -> bedtools intersect + featureCounts log -> MultiQC TSVs)
        |                               +--> annotate_peaks (optional)
        |                               +--> deeptools (optional)
        |                               |    +--> computeMatrix / plotProfile / plotHeatmap (shifted.bigWig for narrow peaks; bigWig for broad peaks)
        |                               |    +--> plotFingerprint / bamPEFragmentSize (filtered.bam)
        |                               +--> nfr (optional; when deeptools + narrow peaks)
        |                               |    +--> alignmentSieve → nfr.bigWig + mono.bigWig
        |                               |    +--> computeMatrix / plotProfile / plotHeatmap (nfr vs mono)
        |                               +--> ataqv (optional; filtered.bam)
        |                               |    +--> ataqv_score → ataqv_score.tsv (TSS enrichment, NFR ratio)
        |                               +--> atacseqqc (optional; narrow peaks only)
        |                               |    +--> fragSizeDist → fragsize_dist.png
        |                               |    +--> PTscore → pt_score.png + atacseqqc_score.tsv
        |                               |    +--> NFRscore → nfr_score.png + atacseqqc_score.tsv
        |                               |    +--> TSSEscore → tsse.png + atacseqqc_score.tsv
        |
        v
     multiqc
        |
        v
    delete_tmp
```

## Rule Quick Reference

All per-sample rules in execution order. Toggle columns: `t` = threads, `MB` = `mem_mb`.

| # | Step | Rule(s) | Module | Tool | t | MB | Input → Key Output |
|---|------|---------|--------|------|:-:|:--:|-------------------|
| 1 | Merge lanes | `merge_raw_fastqs` | common.smk | shell (cat / symlink) | 1 | — | raw FASTQ lanes → `raw_merged/*_merged_{1,2}.fastq.gz` |
| 2 | FastQC raw | `fastqc_raw` | qc.smk | FastQC | 2 | 4 096 | merged FASTQ → `fastqc_raw/*_raw_{1,2}_fastqc.{html,zip}` |
| 3 | Trim | `trim_galore` **or** `fastp` | trim.smk | Trim Galore / fastp | 12 | 36 864 / 16 384 | merged FASTQ → `trim/*_trimmed_{1,2}.fastq.gz` |
| 4 | FastQC trimmed | `fastqc_trimmed` | qc.smk | FastQC | 2 | 4 096 | trimmed FASTQ → `trim/*_trimmed_{1,2}_fastqc.{html,zip}` |
| 5 | BWA-MEM2 index | `bwa_mem2_index` | align.smk | bwa-mem2 index | 12 | 65 536 | FASTA → `.amb / .ann / .bwt.2bit.64 / .pac / .0123` |
| 6 | Bowtie2 index | `bowtie2_index` | align.smk | bowtie2-build | 12 | 8 192 | FASTA → `*.{1,2,3,4,rev.1,rev.2}.bt2` |
| 7 | Align | `bwa_mem2_align` **or** `bowtie2_align` | align.smk | BWA-MEM2 / Bowtie2 + samtools view | 12 / 26 | 49 152 / 16 384 | trimmed FASTQ + index → `bam/*.unsorted.bam` |
| 8 | Sort BAM | `sort_bam` | align.smk | samtools sort + index | 6 | 36 864 | `.unsorted.bam` → `bam/*.bam + *.bam.bai` |
| 9 | Mark duplicates | `mark_duplicates` | mark_duplicates.smk | Picard MarkDuplicates | 2 | 49 152 | `.bam` → `bam/*.markdup.sorted.bam + *.MarkDuplicates.metrics.txt` |
| 10 | Samtools stats pre-filter | `samtools_stats_pre_filter` | align_stats.smk | samtools stats / flagstat / idxstats | 1 | 1 024 | pre-filter BAM → `bam/*.pre_filter.bam.{stats,flagstat,idxstats}` |
| 11 | BAM filter | `bam_filter` | bam_filter.smk | samtools view + bamtools filter + pysam bampe_rm_orphan | 6 | 36 864 | `.markdup.sorted.bam` + include_regions → `bam/*.filtered.bam + *.bai` |
| 12 | Samtools stats | `samtools_stats` | align_stats.smk | samtools stats / flagstat / idxstats | 1 | 1 024 | `.filtered.bam` → `bam/*.filtered.bam.{stats,flagstat,idxstats}` |
| 13 | Picard metrics | `picard_collect_multiple_metrics` | align_stats.smk | Picard CollectMultipleMetrics | 1 | 16 384 | `.filtered.bam` → alignment_summary / insert_size / base_dist / quality_* metrics |
| 14 | BigWig (unshifted) | `bedtools_genomecov` → `ucsc_bedgraphtobigwig` | bam_to_bigwig.smk | bedtools genomecov + bedGraphToBigWig | 2 / 1 | 40 960 / 1 024 | `.filtered.bam` → `bigwig/*.bedGraph` → `bigwig/*.bigWig` |
| 15 | Shift BAM | `shift_bam` | shift_bam.smk | deepTools alignmentSieve + bamCoverage (RPGC) | 8 | 40 960 | `.filtered.bam` → `bam/*.shifted.bam` + `bigwig/*.shifted.bigWig` |
| 16 | MACS3 peak calling | `macs3_callpeak_tn5` | call_peaks.smk | bedtools bamtobed + awk Tn5-shift + MACS3 (narrow: BED mode; broad: BAMPE mode) | 2 | 8 192 | `.filtered.bam` → `peaks/*.tn5_shifted.bed + *_peaks.peak + *_peaks.xls` |
| 17 | MACS3 peak QC | `macs3_peak_qc_plot` | call_peaks.smk | R (plot_macs_qc.r) | 2 | — | `*_peaks.peak` → `peaks/*.macs_peakqc.summary.txt + *.plots.pdf` |
| 18 | featureCounts | `featurecounts_in_peaks` | quant.smk | featureCounts / Subread (SAF, paired, unstranded) | 1 | 2 048 | `.filtered.bam` + peaks SAF → `featurecounts/*.readCountInPeaks.txt` |
| 19 | FRiP score | `frip_score` | frip_score.smk | bedtools intersect (+ featureCounts log parse) | 1 | 2 048 | `.filtered.bam` + peaks + flagstat → `peaks/*.FRiP.txt + *_peaks.{FRiP,count}_mqc.tsv` |
| 20 | HOMER annotation | `homer_annotate_peaks` | annotate_peaks.smk | HOMER annotatePeaks.pl + R (plot_homer_annotatepeaks.r) | 1 | 10 240 | peaks + FASTA + GTF → `annotation/*_peaks.annotatePeaks.txt + *.summary.txt` |
| 21 | deepTools | `deeptools_compute_matrix_scale_regions` / `deeptools_compute_matrix_reference_point` / `deeptools_plot_profile*` / `deeptools_plot_heatmap` / `deeptools_plot_fingerprint` / `deeptools_fragment_size_distribution` | deeptools.smk | deepTools computeMatrix / plotProfile / plotHeatmap / plotFingerprint / bamPEFragmentSize | 2–12 | 6 144–20 480 | `.shifted.bigWig` (narrow) or `.bigWig` (broad) + `.filtered.bam` → `deeptools/` matrices, profiles, heatmaps, fingerprint, fragment-size plots |
| 22 | NFR analysis | `nfr_bigwig_nfr` / `nfr_bigwig_mono` / `nfr_compute_matrix` / `nfr_plot_profile` / `nfr_plot_heatmap` | nfr.smk | deepTools alignmentSieve + bamCoverage + computeMatrix / plotProfile / plotHeatmap | 2–12 | 4 096–20 480 | `.shifted.bam` → `nfr/*.nfr.bigWig + *.mono.bigWig` + NFR-vs-mono TSS profile/heatmap |
| 23 | ataqv | `ataqv` / `ataqv_mkarv` / `ataqv_score` | ataqv.smk | ataqv + mkarv + Python (extract_ataqv_score.py) | 1 | 6 144 / 1 024 / 256 | `.filtered.bam` + peaks + TSS + autosomes → `ataqv/*.ataqv.json + *.mkarv_html/index.html + *.ataqv_score.tsv` |
| 24 | ATACseqQC | `atacseqqc_score` | atacseqqc.smk | ATACseqQC R pkg (calc_pt_score.R) | 1 | 16 384 | `.shifted.bam` + BED → `atacseqqc/*.atacseqqc_score.tsv + fragsize/pt_score/nfr_score/tsse .png` |

**Module toggle summary** (config.yml):

| Module | Config key | Default | Requires |
|--------|-----------|:-------:|---------|
| Trimming | `trimming.enabled` | true | — |
| Mark duplicates | `markduplicates.enabled` | true | — |
| BAM filter | *(always on)* | always | — |
| Peak calling | `call_peaks.enabled` | true | bam_filter |
| Peak QC plot | `call_peaks.macs3_peak_qc_plot` | true | call_peaks |
| Annotation | `annotate_peaks.enabled` | true | call_peaks |
| featureCounts + FRiP | *(auto with call_peaks)* | — | call_peaks |
| Shift BAM + shifted bigWig | *(auto for narrow peaks)* | — | narrow peaks |
| deepTools | `deeptools.enabled` | true | call_peaks |
| NFR analysis | *(auto)* | — | deeptools + narrow peaks |
| ataqv | `ataqv.enabled` | true | call_peaks |
| ATACseqQC | `atacseqqc.enabled` | true | call_peaks + narrow peaks |

---

## Requirements

- Linux
- Snakemake ≥ 8 (in a dedicated controller environment)
- Conda/Mamba

> **HPC users:** skip this section and follow [Running on HPC with LSF](#running-on-hpc-with-lsf) instead, which covers environment setup outside your home directory.

For local use, create a minimal controller environment:

```bash
mamba create -n atacseq_snakemake -c conda-forge -c bioconda snakemake
mamba activate atacseq_snakemake
```

Each rule uses its own isolated Conda environment defined in `workflow/envs/*.yml`.
Pass `--use-conda` on every Snakemake invocation so these per-rule envs are built and activated automatically.

## Installation

```bash
git clone https://github.com/UKHD-NP/atacseq_snakemake_new.git
cd atacseq_snakemake_new
```

## Input Files

### 1. Samplesheet

`config/config.yml` key: `samples_csv`

Required columns:
- `sample_id`
- `fq1`
- `fq2`
- `outdir`

Example:

```csv
sample_id,fq1,fq2,outdir
S1,/data/S1_L001_R1.fastq.gz,/data/S1_L001_R2.fastq.gz,results/S1
S1,/data/S1_L002_R1.fastq.gz,/data/S1_L002_R2.fastq.gz,results/S1
S2,/data/S2_R1.fastq.gz,/data/S2_R2.fastq.gz,results/S2
```

Notes:
- Repeated `sample_id` rows are treated as lanes and merged.
- Even with a single lane/sample, the pipeline still runs `merge_raw_fastqs` (symlink mode), so merged FASTQ filenames still use the `*_merged_*` suffix.
- All rows for the same `sample_id` must have the same `outdir`.

### 2. Reference

Set in `config/config.yml -> ref`.

- `assembly`: `hg19`, `hg38`, or `custom`
- For `custom`, provide:
  - `ref.fasta`
  - `ref.gtf`
  - `ref.blacklist`

Pipeline stages references into `references/{assembly}/` and derives:
- `.fai`, `.sizes` (`chromsizes`), `.autosomes.txt`, `.tss.bed`, `.include_regions.bed`
- `ref.bed` is generated from GTF by `prepare_genome` using the dedicated env `workflow/envs/gtf2bed.yml` (Perl + gzip/unzip) for portable HPC runs.

## Local Run

> For cluster execution on HPC, see [Running on HPC with LSF](#running-on-hpc-with-lsf) below.
> The commands here are for single-machine (local) execution only.

**Step 1 — Dry-run first (always).**
Resolves the full DAG and prints every rule that would run — without executing anything:

```bash
snakemake -s workflow/Snakefile --use-conda -n
```

**Step 2 — Optionally verify with the bundled test dataset.**
Runs the full pipeline end-to-end on small test data:

```bash
snakemake -s workflow/Snakefile \
    --configfile config/config_test.yml \
    --use-conda --conda-frontend mamba \
    --cores all
```

**Step 3 — Run with your real config.**

```bash
# Normal run
snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --cores 24

# Rerun only failed/incomplete jobs after fixing an error
snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --cores 24 --rerun-incomplete
```

> `config/config.yml` is loaded automatically by the Snakefile as the default configfile.
> Pass `--configfile path/to/other.yml` only when you want to override it (e.g. for a test config).

## Running on HPC with LSF

This setup uses **IBM Spectrum LSF**.
A ready-made LSF profile is provided at `workflow/profiles/lsf/config.yaml`.

### Node roles

| Node | Purpose | Allowed |
|------|---------|---------|
| `<worker-node>` | Dev, install, testing | ✅ Software install, small runs |
| `<submit-node>` | Job submission only | ✅ Run Snakemake (lightweight), ❌ Processing |
| Cluster nodes | Computation | Jobs submitted automatically via `bsub` |

### Step 1 - Set up Snakemake environment

> **Do this on `<worker-node>`, not on `<submit-node>`.**
> Worker nodes allow software installation. Submission hosts do not.

```bash
ssh YOUR_USERNAME@<worker-node>
```

**Configure conda channels.**
Some HPC clusters ban the `defaults` (Anaconda) channel due to licensing restrictions.
You may need to explicitly restrict to `conda-forge` and `bioconda`:

```bash
cat > ~/.condarc << 'EOF'
channels:
  - conda-forge
  - bioconda
EOF
```

**Load Mamba and initialise your shell.**
This adds `mamba`/`conda` to your `PATH` permanently via `~/.bashrc`:

```bash
module load Mamba   # adjust module name to your site
mamba init bash
source ~/.bashrc   # apply changes to the current shell without re-logging in
```

**Create the Snakemake controller environment outside your home directory.**
Home quota on HPC systems is often limited. Conda environments can easily exceed this — install them on group storage:

```bash
# Set your working directory on group storage (adjust path as needed)
YOUR_WORKDIR="/path/to/group/storage/YOUR_USERNAME"
mkdir -p ${YOUR_WORKDIR}/conda_envs

# Create the controller environment with Snakemake + the LSF executor plugin
mamba create -p ${YOUR_WORKDIR}/conda_envs/snakemake \
    -c conda-forge -c bioconda \
    python=3.11 \
    "numpy<1.25" \
    snakemake=8.30.0 \
    snakemake-executor-plugin-lsf \
    -y

# Activate the new environment
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake
```

> `snakemake-executor-plugin-lsf` translates Snakemake rule resources (`mem_mb`, `runtime`, `threads`) into `bsub` submission flags automatically — no manual `bsub` scripting needed.

### Step 2 - Clone the pipeline

```bash
cd ${YOUR_WORKDIR}
git clone https://github.com/UKHD-NP/atacseq_snakemake_new.git
cd atacseq_snakemake_new
```

### Step 3 - Edit configuration

Open `config/config.yml` and set at minimum:
- `samples_csv`: path to your samplesheet CSV
- `ref.assembly`: `hg19`, `hg38`, or `custom`
- Output directories (via the `outdir` column in the samplesheet)
- Enable/disable optional modules (`deeptools`, `ataqv`, `annotate_peaks`, etc.)

See the [Key Configuration](#key-configuration) section below for all options and defaults.

### Step 4 - Update `conda-prefix` in the LSF profile

`conda-prefix` tells Snakemake where to build and cache the per-rule conda environments (from `workflow/envs/*.yml`).
All rule environments combined take roughly **5–15 GB** and must live outside your home directory.

Update the placeholder path to your actual working directory:

```bash
sed -i "s|/path/to/group/storage/conda_envs|${YOUR_WORKDIR}/conda_envs|g" \
    workflow/profiles/lsf/config.yaml

# Confirm the replacement was applied correctly
grep "conda-prefix" workflow/profiles/lsf/config.yaml
```

> **Note:** Add the following line to your `~/.bashrc` (once, then `source ~/.bashrc`).
> LSF enforces memory limits per-job, so this variable tells the LSF plugin to
> submit the full `mem_mb` value as a per-job request instead of dividing it per slot:
>
> ```bash
> export SNAKEMAKE_LSF_MEMFMT=perjob
> ```

### Step 5 - Validate with a dry-run

Resolves the full DAG and prints every rule that would run — **without executing or submitting any jobs**.
Always do this before submitting to the cluster to catch config errors, missing inputs, or unexpected rule counts.

```bash
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake
cd ${YOUR_WORKDIR}/atacseq_snakemake_new

# Dry-run: prints all rules, checks all inputs, submits nothing
snakemake -s workflow/Snakefile --use-conda -n
```

Confirm that the printed rule count and sample names match expectations before proceeding to Step 6.

> For local testing with the bundled test dataset, see the [Local Run](#local-run) section.

### Step 6 - Submit to HPC

> **Do this on `<submit-node>`**, not on `<worker-node>`.
> Snakemake must run on a submission host to dispatch jobs via `bsub`.

Use `screen` so the Snakemake controller process survives SSH disconnects:

```bash
ssh YOUR_USERNAME@<submit-node>

# Create a named screen session — it keeps running after SSH disconnect
screen -S <session_name>

# Set your working directory (same value as used in Step 1)
YOUR_WORKDIR="/path/to/group/storage/YOUR_USERNAME"

# Activate the Snakemake controller environment
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

# Move into the pipeline directory
cd ${YOUR_WORKDIR}/atacseq_snakemake_new

# Launch the pipeline — Snakemake submits each rule as a separate bsub job automatically.
# The config/config.yml is loaded automatically from the Snakefile; no --configfile needed.
# Concurrency is controlled by `jobs:` in workflow/profiles/lsf/config.yaml.
snakemake --profile workflow/profiles/lsf
```

To rerun only failed/incomplete jobs after fixing an error:

```bash
snakemake --profile workflow/profiles/lsf --rerun-incomplete
```

To rerun with the test dataset config:

```bash
snakemake --profile workflow/profiles/lsf --rerun-incomplete --configfile config/config_test.yml
```

Force rerun examples:

```bash
# Force one rule for all matching jobs (e.g. rerun all trim_galore jobs)
snakemake --profile workflow/profiles/lsf --forcerun trim_galore

# Force specific output files (target-level force)
snakemake --profile workflow/profiles/lsf --force \
  test_data/results/SAMPLE_ID/trim/SAMPLE_ID_trimmed_1.fastq.gz \
  test_data/results/SAMPLE_ID/trim/SAMPLE_ID_trimmed_2.fastq.gz

# Force all jobs in the DAG to rerun from scratch
snakemake --profile workflow/profiles/lsf --forceall
```

| `screen` command | Action |
|-----------------|--------|
| `screen -S <session_name>` | Start new named session |
| `Ctrl+A`, then `D` | Detach - session keeps running after SSH disconnect |
| `screen -ls` | List all active sessions |
| `screen -r <session_name>` | Re-attach to session |
| `screen -S <session_name> -X quit` | Kill the named session |

### Monitoring jobs

| `bjobs` command | Action |
|-----------------|--------|
| `bjobs -w` | List all running/pending jobs |
| `bjobs -w -r` | Running only |
| `bjobs -w -p` | Pending only |
| `bjobs -l JOB_ID` | Detailed info for one job |

## Key Configuration

Below is the default-style config block with practical explanation:

```yaml
# Path to samplesheet CSV (columns: sample_id, fq1, fq2, outdir)
samples_csv: "samplesheet.csv"

ref:
  assembly: "custom"            # hg19 / hg38 / custom
  fasta: "genome.fa"
  gtf: "genes.gtf"
  blacklist: "blacklist.v3.bed"
  bwa_index: ""                 # optional prebuilt BWA prefix; auto-generated if empty
  bowtie2_index: ""             # optional prebuilt Bowtie2 prefix; auto-generated if empty
  mito_name: "chrM"               # MUST match FASTA mito contig exactly (e.g. MT/chrM/M)
  keep_mito: false              # false = exclude mitochondrial reads from include_regions
  
trimming:
  enabled: true
  delete_trimming: true         # delete trimmed FASTQs after pipeline completes
  tool: "trim_galore"           # fastp / trim_galore
  trim_galore_params: "--nextseq 25 --length 36"
  fastp_params: "--cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --trim_poly_g --length_required 36"

align:
  tool: "bowtie2"               # bwa-mem2 / bowtie2
  bowtie2_params: "--very-sensitive --no-discordant -X 2000"
  bwa_params: "-I 0,2000"

bam_filter:
  params: "-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400 -q 30"
  apply_canonical_chromosomes: false   # hg19/hg38: set true. Non-human or custom genomes: set false (see below)
  apply_blacklist: true        # true = exclude blacklist regions; false = keep them
  keep_input_bam: false         # true = preserve BAM before bam_filter; false = delete to save space

markduplicates:
  enabled: true

call_peaks:
  enabled: true
  peak_type: "narrow"           # narrow / broad
  macs3_gsize: "2701495711"     # effective genome size (preferred); if empty, auto-sum from chromsizes
  macs3_narrow_params: "--trackline --shift -75 --extsize 150 --keep-dup all --nomodel --call-summits -q 0.01"
  macs3_broad_params: "--trackline --keep-dup all --nomodel --broad --broad-cutoff 0.1"
  macs3_peak_qc_plot: true      # run plot_macs_qc.r to generate peak QC summary/plots
  frip_overlap_fraction: 0.2
  frip_threshold: 20

annotate_peaks:
  enabled: true

deeptools:
  enabled: true

# Optional: tune NFR/mononucleosomal fragment size boundaries.
# NFR analysis runs automatically when deeptools + narrow peaks are enabled.
nfr:
  nfr_max_fragment: 150      # fragments ≤ this bp → NFR bigWig
  mono_min_fragment: 151     # fragments ≥ this bp → mono bigWig
  mono_max_fragment: 300     # fragments ≤ this bp → mono bigWig

ataqv:
  enabled: true

atacseqqc:
  enabled: true

multiqc:
  config: "workflow/scripts/multiqc_config.yml"

latency-wait: 60
```

Explanation by block:
- `samples_csv`: input table for sample discovery and lane merging.
- `ref`: reference genome/annotation source. For `custom`, `fasta` and `gtf` are required. `blacklist` is only required when `bam_filter.apply_blacklist: true`. `bwa_index` / `bowtie2_index` can be left empty to auto-build.
- `ref.mito_name`: **critical** — must exactly match the mitochondrial contig name in your FASTA (often `MT`, `chrM`, or `M`). Used by `bam_filter` to build `include_regions` and by `ataqv`.
- `ref.keep_mito`: set `true` to retain mitochondrial reads in `include_regions`; `false` (default) excludes them.
- `ref.autosome_pattern` *(optional)*: awk regex for the autosome list passed to `ataqv --autosomal-reference-file`. Default covers **hg19/hg38**: `^chr([1-9]|1[0-9]|2[0-2])$`. Override for other genomes, e.g. `^([0-9]+)$` for ENSEMBL naming or `^(I|II|III|IV|V)$` for `ce11`.
- `trimming`: choose one trimming engine and pass tool-specific options.
- `align`: choose aligner and set aligner-specific CLI parameters.
- `bam_filter.params`: SAMtools core filter flags; see [BAM Filtering](#re-mark-duplicates-and-bam-filtering-criteria) for full breakdown.
- `bam_filter.apply_canonical_chromosomes`: controls whether reads are restricted to standard chromosomes before blacklist and mito filtering.
  - `true` — filter chromosomes using `canonical_chroms_pattern` (see below). Removes noise from unplaced contigs (`chrUn_*`, `*_random`, EBV, decoy sequences) that inflate peak-calling background.
  - `false` — keep all contigs present in the FASTA. Safe for any genome without configuration.
  - `ref.keep_mito` still controls whether the mitochondrial contig is included in the final region set regardless of this flag.
- `bam_filter.canonical_chroms_pattern` *(optional)*: awk regex applied when `apply_canonical_chromosomes=true`. Default covers **hg19/hg38**: `^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$`. Override for other genomes:

  | Genome | Pattern |
  |--------|---------|
  | hg19 / hg38 (default) | `^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$` |
  | Mouse mm10/mm39, rat rn7, zebrafish danRer11 | `^chr([0-9]+|X|Y|M)$` |
  | ENSEMBL naming (no `chr` prefix, e.g. `1`, `X`, `MT`) | `^([0-9]+|X|Y|MT)$` |
  | C. elegans ce11 (I–V, X, MtDNA) | `^(I{1,3}|IV|V{1,3}|X|MtDNA)$` |
- `bam_filter.apply_blacklist`: set `true` (default) to exclude blacklist intervals from `include_regions`; set `false` to keep blacklist regions while still respecting `ref.keep_mito`.
- `bam_filter.keep_input_bam`: set `true` to preserve the BAM entering `bam_filter` (`*.markdup.sorted.bam` when duplicate marking is enabled, otherwise `*.bam`).
- `markduplicates.enabled`: run Picard MarkDuplicates before filtering. When disabled, duplicates are not flagged and `-F 0x0400` in `bam_filter.params` has no effect.
- `trimming.delete_trimming`: when `true`, trimmed FASTQ files are deleted after the pipeline completes.
- `call_peaks.peak_type`: `narrow` uses filtered BAM → `bamtobed` → awk Tn5 shift (+4 forward / -5 reverse) → MACS3 BED mode; `broad` uses filtered BAM directly in MACS3 BAMPE mode.
- `call_peaks.macs3_peak_qc_plot`: when `true`, runs `plot_macs_qc.r` to produce `*.macs_peakqc.summary.txt` and `*.macs_peakqc.plots.pdf`.
- `call_peaks.frip_overlap_fraction`: minimum read-peak overlap fraction for FRiP (passed to both `bedtools intersect -f` and featureCounts `--fracOverlap`).
- `call_peaks.frip_threshold`: FRiP percentage threshold for quality label in `*.FRiP.txt`; samples at or above this value are labelled `good`, below is `bad` (default: 20%).
- `annotate_peaks.enabled`: run HOMER `annotatePeaks` and summary plotting.
- `deeptools.enabled`: run computeMatrix/plotProfile/plotHeatmap/plotFingerprint modules. Requires `call_peaks.enabled=true`. For narrow peaks, computeMatrix uses the Tn5-shifted bigWig (`shifted.bigWig`); for broad peaks, it uses the unshifted bigWig (`bigWig`).
- `nfr` (optional block): tune fragment size cutoffs for nucleosome-free region score (NFR) vs mononucleosomal analysis. Omit to use defaults (NFR ≤150 bp, mono 151–300 bp).
- `ataqv.enabled`: run ATAC-specific QC (`ataqv`) and render interactive HTML (`mkarv`), plus extract short mononucleosomal ratio and TSS enrichment score (TSSE)  score to `ataqv_score.tsv` for MultiQC. Requires `call_peaks.enabled=true`.
- `atacseqqc.enabled`: run ATACseqQC R package to compute Promoter/transcript body score (PT), per-TSS nucleosome-free region score (NFR), and TSS enrichment score (TSSE) and produce QC plots. Outputs `atacseqqc_score.tsv` for MultiQC. Requires `call_peaks.enabled=true` and `call_peaks.peak_type=narrow`.
- `multiqc.config`: path to MultiQC config used by this pipeline.
- `latency-wait`: useful on slow/network filesystems to avoid false missing-output errors.

### Why these default params were chosen

- `trimming.tool: trim_galore` with `--nextseq 25 --length 36`:
  chosen for two-color Illumina runs (poly-G prone) and to remove very short reads that are usually uninformative for peak calling.
- `align.tool: bowtie2` with `--very-sensitive --no-discordant -X 2000`:
  chosen to maximize paired-end sensitivity while constraining improbable pair structure for ATAC fragment lengths.
- `bam_filter.params: ... -q 30`:
  chosen as a relatively strict default to retain high-confidence alignments for downstream peak calling and signal tracks.
- `call_peaks` defaults (`--shift -75 --extsize 150 --nomodel`):
  for narrow peaks, reads are Tn5-shifted inline by awk (+4 on forward strand, -5 on reverse strand) before being passed to MACS3 BED mode. `--shift -75` then pulls each cut-site tag a further 75 bp upstream so that the 150 bp extension lands symmetrically around the insertion site. Positive shift would offset windows away from the cut site.

### How to choose `macs3_gsize`

- `macs3_gsize` is passed to MACS3 as `--gsize`.
- Prefer effective genome size values from deepTools documentation:  
  `https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html`
- If `macs3_gsize` is empty, this pipeline falls back to summing chromosome sizes from `ref.chromsizes`.
- Common values (from deepTools table):
  - `hg19`: `2864785220`
  - `hg38`: `2913022398`

## Re-mark Duplicates and BAM Filtering Criteria

After alignment and coordinate sorting, duplicates are marked with Picard (`mark_duplicates`, when enabled), then `bam_filter` creates `*.filtered.bam`.

Important: `bam_filter.params` is the **SAMtools core filter only**.  
Additional filters are applied by BAMTools and Pysam (if available in the runtime env).

Why this default filtering strategy is used:
- It prioritizes specificity for ATAC peak detection.
- It reduces noisy alignments before MACS3, FRiP, and bigWig generation.
- It is a practical strict default (`-q 30`); for low-depth data you can relax to `-q 20`.

Practical examples:
- Keep mitochondrial filtering, but disable blacklist filtering:
  `ref.keep_mito: false` and `bam_filter.apply_blacklist: false`
- Preserve the BAM before filtering for debugging or alternate peak-calling runs:
  `bam_filter.keep_input_bam: true`

### 1) SAMtools core filter (from `config.yml`)

Current default:

`-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400 -q 30`

Meaning:
- `-F 0x004`: remove unmapped reads
- `-F 0x0008`: remove reads whose mate is unmapped
- `-f 0x001`: keep only reads flagged as paired
- `-F 0x0100`: remove secondary alignments
- `-F 0x0400`: remove reads marked as duplicates
- `-q 30`: keep reads with MAPQ >= 30 (remove lower-confidence multimappers/ambiguous mappings)

`bam_filter` also uses `-L ref.include_regions` with `samtools view`:
- this always constrains BAMs to the generated include regions
- canonical chromosomes are enforced there when `bam_filter.apply_canonical_chromosomes: true`
- blacklist intervals are excluded there when `bam_filter.apply_blacklist: true`
- mitochondrial contig is excluded there when `ref.keep_mito: false`
- if `bam_filter.apply_blacklist: false` and `ref.keep_mito: false`, `chrM`/`MT` is still filtered out as long as `ref.mito_name` matches the FASTA

### 2) BAMTools extra filter (optional but enabled when `bamtools` exists)

From `workflow/scripts/bamtools_filter_pe.json`, extra constraints are:
- mismatches `NM <= 4`
- remove soft-clipped reads (`CIGAR` containing `S`) — **note:** this rejects any read with even 1 bp soft-clip, which can be overly aggressive; see `resources/pipeline_comparison.md` Section 7 for a recommended fix
- keep insert size in `[-2000, 2000]`

### 3) Pysam pair/orphan cleanup (optional but enabled when `python3+pysam` exists)

From `workflow/scripts/bampe_rm_orphan.py` (`--only_fr_pairs` mode):
- remove singleton/orphan reads
- keep only read pairs on the same chromosome
- keep only FR-oriented proper pairs
- remove pairs where one mate fails the pair criteria

## Main Outputs

Per sample under `<outdir>`:

- BAM
  - `bam/{sample}.filtered.bam`
  - `bam/{sample}.filtered.bam.bai`
  - `bam/{sample}.shifted.bam` (narrow peaks only)
  - `bam/{sample}.shifted.bam.bai` (narrow peaks only)
  - `bam/{sample}.markdup.sorted.bam` / `bam/{sample}.markdup.sorted.bam.bai` or `bam/{sample}.bam` / `bam/{sample}.bam.bai` may also be retained when `bam_filter.keep_input_bam: true`
- BAM stats
  - `bam/{sample}.pre_filter.bam.stats`
  - `bam/{sample}.pre_filter.bam.flagstat`
  - `bam/{sample}.pre_filter.bam.idxstats`
  - `bam/{sample}.filtered.bam.stats`
  - `bam/{sample}.filtered.bam.flagstat`
  - `bam/{sample}.filtered.bam.idxstats`
  - `bam/{sample}.markdup.sorted.MarkDuplicates.metrics.txt`
  - `bam/{sample}.CollectMultipleMetrics.alignment_summary_metrics`
  - `bam/{sample}.CollectMultipleMetrics.base_distribution_by_cycle.pdf`
  - `bam/{sample}.CollectMultipleMetrics.base_distribution_by_cycle_metrics`
  - `bam/{sample}.CollectMultipleMetrics.insert_size_histogram.pdf`
  - `bam/{sample}.CollectMultipleMetrics.insert_size_metrics`
  - `bam/{sample}.CollectMultipleMetrics.quality_by_cycle.pdf`
  - `bam/{sample}.CollectMultipleMetrics.quality_by_cycle_metrics`
  - `bam/{sample}.CollectMultipleMetrics.quality_distribution.pdf`
  - `bam/{sample}.CollectMultipleMetrics.quality_distribution_metrics`
- Signal tracks
  - `bigwig/{sample}.bedGraph`
  - `bigwig/{sample}.scale_factor.txt`
  - `bigwig/{sample}.bigWig`
  - `bigwig/{sample}.shifted.bigWig` (narrow peaks only)
- Peaks / FRiP
  - `peaks/{sample}.tn5_shifted.bed` (narrow peaks only)
  - `peaks/{sample}_peaks.peak`
  - `peaks/{sample}_peaks.xls`
  - `peaks/{sample}.FRiP.txt`
  - `peaks/{sample}_peaks.FRiP_mqc.tsv`
  - `peaks/{sample}_peaks.count_mqc.tsv`
  - `peaks/{sample}.macs_peakqc.summary.txt`
  - `peaks/{sample}.macs_peakqc.plots.pdf`
  - `featurecounts/{sample}.readCountInPeaks.txt`
  - `featurecounts/{sample}.readCountInPeaks.txt.summary`
- Annotation
  - `annotation/{sample}_peaks.annotatePeaks.txt`
  - `annotation/{sample}.macs_annotatePeaks.summary.txt`
- Read QC / trimming
  - `fastqc_raw/{sample}_raw_1_fastqc.html` + `.zip`
  - `fastqc_raw/{sample}_raw_2_fastqc.html` + `.zip`
  - `trim/{sample}_1.fastq.gz_trimming_report.txt` and `trim/{sample}_2.fastq.gz_trimming_report.txt` (Trim Galore mode)
  - `trim/{sample}_trimmed_1.fastq.gz` and `trim/{sample}_trimmed_2.fastq.gz` (unless deleted during cleanup)
  - `trim/{sample}_trimmed_1_fastqc.html` + `.zip` and `trim/{sample}_trimmed_2_fastqc.html` + `.zip` (Trim Galore mode)
  - `trim/{sample}.fastp.html` and `trim/{sample}.fastp.json` (fastp mode)
- deepTools
  - `deeptools/{sample}.scale_regions.computeMatrix.gz`
  - `deeptools/{sample}.scale_regions.computeMatrix.tab`
  - `deeptools/{sample}.scale_regions.plotProfile.pdf` + `.tab`
  - `deeptools/{sample}.reference_point.computeMatrix.gz`
  - `deeptools/{sample}.reference_point.computeMatrix.tab`
  - `deeptools/{sample}.reference_point.plotProfile.pdf` + `.tab`
  - `deeptools/{sample}.reference_point.plotHeatmap.pdf` + `.tab`
  - `deeptools/{sample}.plotFingerprint.pdf`
  - `deeptools/{sample}.plotFingerprint.raw_counts.txt`
  - `deeptools/{sample}.plotFingerprint.qcmetrics.txt`
  - `deeptools/{sample}.fragment_size_distribution.pdf`
  - `deeptools/{sample}.fragment_size.raw_lengths.txt`
  - `deeptools/{sample}.fragment_size.qcmetrics.txt`
- NFR analysis
  - `nfr/{sample}.nfr.bigWig` (fragments ≤150 bp)
  - `nfr/{sample}.mono.bigWig` (fragments 151–300 bp)
  - `nfr/{sample}.nfr_vs_mono.computeMatrix.gz`
  - `nfr/{sample}.nfr_vs_mono.computeMatrix.tab`
  - `nfr/{sample}.nfr_vs_mono.plotProfile.pdf` + `.tab`
  - `nfr/{sample}.nfr_vs_mono.plotHeatmap.pdf` + `.tab`
- ataqv
  - `ataqv/{sample}.ataqv.json`
  - `ataqv/{sample}.mkarv_html/index.html`
  - `ataqv/{sample}.ataqv_score.tsv` (NFR ratio + TSSE score for MultiQC)
- ATACseqQC
  - `atacseqqc/{sample}.fragsize_dist.png`
  - `atacseqqc/{sample}.pt_score.png`
  - `atacseqqc/{sample}.nfr_score.png`
  - `atacseqqc/{sample}.tsse.png`
  - `atacseqqc/{sample}.atacseqqc_score.tsv` (PT/NFR/TSSE scores for MultiQC)
- MultiQC
  - `multiqc/{sample}.multiqc.html`
- Cleanup
  - `logs/{sample}.deletion.log`

## QC Visual Guide (from `resources/`)

### 1. FastQC base composition

![FastQC base content](resources/fastqc_base_content.png)

How to read:
- Raw ATAC-seq data can show end-bias and poly-G artifacts (especially on two-color chemistry platforms).
- After trimming, base-content curves should be more stable and less biased at read ends.

### 2. Adapter-content profile

![FastQC adapter content](resources/fastqc_adapter_content.png)

How to read:
- Before trimming, adapter contamination can be high at read tails.
- After trimming, adapter signal should drop strongly (ideally near zero across most cycles).

### 3. Fragment-size distribution

Produced by: `bamPEFragmentSize` → `deeptools/{sample}.fragment_size_distribution.pdf`

![ATAC fragment size examples](resources/fragmentSize_distribution_examples.svg)

How to read:
- Good ATAC libraries usually show a nucleosome-free peak (short fragments) plus mono-/di-nucleosome periodic peaks.
- A flat/noisy pattern without clear peaks often indicates lower signal quality.

### 4. Fingerprint / library complexity trend

![Lorenz / fingerprint examples](resources/lorenz_curve_examples.svg)

How to read:
- Curves help assess enrichment and library complexity.
- Better enrichment typically separates signal from background more clearly.

## ENCODE QC Benchmarks

Sources: [ENCODE ATAC-seq Standards](https://www.encodeproject.org/atac-seq/) · [ataqv](https://github.com/ParkerLab/ataqv) · [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/)

### Quick reference

| Metric | Source | Tool / Output | Target | Acceptable |
|--------|--------|--------------|--------|------------|
| Alignment rate | ENCODE | Bowtie2 / BWA-MEM2 log | >95% | ≥80% |
| Duplication rate | general practice | `*.MarkDuplicates.metrics.txt` | <20% | <30% |
| FRiP score | ENCODE | `*_peaks.FRiP_mqc.tsv` | ≥0.3 | ≥0.2 |
| NFR ratio (short:mono) | ataqv | `ataqv/*.ataqv_score.tsv` | >2 | — |
| TSSE score | ataqv / ENCODE (hg38) | `ataqv/*.ataqv_score.tsv` | ≥7 | ≥5 |
| PT score (2^mean) | ATACseqQC / pipeline | `atacseqqc/*.atacseqqc_score.tsv` | ≥10 | ≥5 |
| NFR score (mean) | ATACseqQC | `atacseqqc/*.atacseqqc_score.tsv` | >0 | — |
| TSSE score | ATACseqQC / ENCODE | `atacseqqc/*.atacseqqc_score.tsv` | ≥7 | ≥5 |

---

### FRiP Score — Fraction of Reads in Peaks

**Source:** ENCODE ATAC-seq Standards

**What it measures:** The proportion of all mapped reads that fall within called peak regions.

```
FRiP = reads_overlapping_peaks / total_mapped_reads
```

**Why it matters:** A high FRiP means most of your sequencing reads captured genuine open chromatin sites, rather than noisy background. Low FRiP suggests poor enrichment, over-amplification, or degraded nuclei.

**Thresholds (from ENCODE):**
- `≥0.3` — strong enrichment; ENCODE recommended target
- `≥0.2` — ENCODE minimum passing standard
- `<0.2` — poor enrichment; inspect fragment size distribution and TSS enrichment

**In this pipeline:** two FRiP estimates are reported per sample — one from `bedtools intersect` and one from `featureCounts`. Both appear in MultiQC. The quality label (`good` / `bad`) in `*.FRiP.txt` uses `call_peaks.frip_threshold` (default 20%).

---

### NFR Ratio — Short-to-Mononucleosomal Fragment Ratio (ataqv)

**Source:** ataqv (`short_mononucleosomal_ratio` field); no official ENCODE threshold — use comparatively across samples

**What it measures:** The ratio of sub-nucleosomal (TF-bound) fragments to mononucleosomal fragments, as computed by ataqv.

```
NFR ratio = count(fragments ≤ 100 bp)   [hqaa_tf_count]
            ────────────────────────────────────────────
            count(fragments 180–300 bp)  [hqaa_mononucleosomal_count]
```

This is a **ratio**, not a percentage of all reads. A ratio of 5 means there are 5 short TF-bound fragments for every 1 mononucleosomal fragment.

**Why it matters:** A healthy ATAC-seq library should be dominated by sub-nucleosomal insertions (nucleosome-free regions). A low ratio means the library is enriched for mononucleosomal fragments, indicating the nuclei had poor chromatin accessibility or Tn5 over-digestion.

**Thresholds:** No official ENCODE cutoff exists for this metric. As a rough community benchmark:
- `>2` — generally considered adequate (more short NFR fragments than mononucleosomal)
- `<1` — likely poor library; mononucleosomal fragments dominate

> Compare values across your own samples. A consistent drop within a batch is more informative than any single absolute threshold.

**In this pipeline:** extracted from ataqv JSON as `short_mononucleosomal_ratio`, reported in `ataqv/{sample}.ataqv_score.tsv`, shown in MultiQC.

---

### TSSE Score — TSS Enrichment Score (ataqv)

**Source:** ataqv (Parker Lab); thresholds from ENCODE ATAC-seq Standards for GRCh38 + RefSeq TSS annotation

**What it measures:** Signal enrichment at Transcription Start Sites (TSS) relative to the flanking background. Calculated by ataqv.

```
TSS enrichment = mean signal in TSS ±150 bp window
                 ─────────────────────────────────
                 mean signal in flanking regions (1400–2000 bp from TSS)
```

**Why it matters:** Tn5 transposase preferentially inserts at open chromatin concentrated at active promoters. A high TSS enrichment confirms the experiment captured nucleosome-free promoter regions. A flat score indicates high background, poor Tn5 enrichment, or degraded chromatin.

**Thresholds** (GRCh38 RefSeq TSS — scale down for non-human genomes):
- `≥7` — ENCODE target; excellent signal-to-noise
- `5–7` — acceptable; usable but noisier peak calls
- `<5` — poor; re-check nuclei isolation and Tn5 titration

> These thresholds are annotation-dependent. For mouse (mm10) or custom genomes, expect lower absolute values. Compare across your own samples rather than using human thresholds as hard cutoffs.

**In this pipeline:** reported in `ataqv/{sample}.ataqv_score.tsv`, shown in MultiQC.

---

### PT Score — Promoter/Transcript Body Ratio

**Source:** ATACseqQC R package (Ou et al., 2018, *Genome Biology*); thresholds are pipeline-defined heuristics (ATACseqQC does not publish fixed cutoffs)

**What it measures:** Whether Tn5 insertions (5′ read ends from shifted BAM) are enriched at promoters relative to gene bodies. Calculated by the ATACseqQC R package.

```
For each transcript:
  promoter_window = [TSS−2000, TSS+500]   (strand-aware)
  body_window     = next 2500 bp downstream of promoter

  PT score (log2) = log2(mean_5prime_density_in_promoter + ε)
                  − log2(mean_5prime_density_in_body + ε)

Final: mean and median PT score across all transcripts
```

The mean PT score is in log2 scale; the equivalent linear ratio is **2^PT_score_mean**.

**Why it matters:** ATAC-seq signal is concentrated at promoters. A high PT score confirms most signal comes from promoter-proximal nucleosome-free regions. Low PT scores indicate high background, poor Tn5 enrichment, or a ChIP-like signal profile.

**Thresholds (linear scale, 2^mean_PT) — pipeline-defined heuristic:**
- `≥10` — strong enrichment; typical for high-quality libraries
- `≥5` — PASS (pipeline threshold in `calc_pt_score.R`)
- `<5` — FAIL flag written to log; ATACseqQC itself does not publish a fixed cutoff

> The ≥5 cutoff is set in this pipeline's `calc_pt_score.R`, not by the ATACseqQC package or its vignette. Compare across your own samples.

**In this pipeline:** calculated on shifted BAM (`shifted.bam`), reported alongside NFR score in `atacseqqc/{sample}.atacseqqc_score.tsv`. A `[INFO] QC: PASS` or `[WARNING] QC: FAIL` line is written to the rule log. Scatter plot saved to `atacseqqc/{sample}.pt_score.png`. Only available for narrow peak mode.

---

### NFR Score (ATACseqQC) — Per-TSS Nucleosome-Free Region Score

**Source:** ATACseqQC R package (`NFRscore()`, Ou et al., 2018, *Genome Biology*)

**What it measures:** Whether the 100 bp window centred on each TSS is more accessible than its flanking nucleosome positions. Computed per TSS and summarised as mean/median across all TSS.

```
For each TSS (400 bp window, strand-aware):
  n1  = upstream 150 bp    (nucleosome flank)
  nf  = middle  100 bp     (nucleosome-free region)
  n2  = downstream 150 bp  (nucleosome flank)

NFR score = log2(nf) − log2((n1 + n2) / 2)
```

A positive score means the NFR window has more Tn5 insertion signal than the average nucleosome flank. Higher values indicate stronger nucleosome depletion at TSS.

**Why it matters:** Complements the PT score. PT score uses a broad 5 kb promoter vs gene-body window; NFR score zooms into the 400 bp TSS window and directly quantifies nucleosome eviction at the TSS itself.

**Thresholds:** ATACseqQC does not publish fixed cutoffs. A positive mean NFR score (> 0) indicates the expected TSS accessibility pattern; compare across samples within the same experiment.

**In this pipeline:** computed alongside PT score on shifted BAM; mean and median NFR score are reported as additional columns in `atacseqqc/{sample}.atacseqqc_score.tsv`. Scatter plot saved to `atacseqqc/{sample}.nfr_score.png`. Only available for narrow peak mode.

---

### TSSE Score — TSS Enrichment Score (ATACseqQC)

**Source:** ATACseqQC R package (`TSSEscore()`, Ou et al., 2018); definition from [ENCODE data standards](https://www.encodeproject.org/data-standards/terms/#enrichment)

**What it measures:** Aggregate read enrichment at TSS relative to flanking background — the same concept as ataqv's TSS enrichment score but computed independently by ATACseqQC.

```
For each TSS (±1000 bp window, 100 bp steps):
  per-step score = depth at step / mean depth at 100 bp end flanks

TSSE = max(mean(per-step score across all TSS))
```

**Why it matters:** An independent TSS enrichment estimate that can be compared to ataqv's `tss_enrichment_score`. Discordance between the two may indicate annotation or BAM handling differences.

**Thresholds** (GRCh38 RefSeq, from ENCODE):
- `≥7` — ENCODE target
- `5–7` — acceptable
- `<5` — poor; same interpretation as ataqv TSS enrichment

**In this pipeline:** computed alongside PT score and NFR score on shifted BAM; reported as `TSSE_score` column in `atacseqqc/{sample}.atacseqqc_score.tsv`. Plot with ENCODE threshold lines saved to `atacseqqc/{sample}.tsse.png`. Only available for narrow peak mode.

---

### Duplication Rate — Library Complexity

**Source:** General bioinformatics practice (not a direct ENCODE threshold — ENCODE uses NRF/PBC1/PBC2 which require a separate counting step not implemented in this pipeline)

**What it measures:** The fraction of reads flagged as PCR/optical duplicates by Picard MarkDuplicates.

```
Duplication rate = duplicate reads / total mapped reads
```

**Why it matters:** High duplication indicates over-amplification or low-complexity library — most reads are copies of the same fragment rather than independent Tn5 insertions. This artificially inflates peak signal.

**Thresholds (community practice):**
- `<20%` — good library complexity
- `20–30%` — acceptable; consider using more input material next time
- `>30%` — poor complexity; reduce PCR cycles or increase input

**In this pipeline:** reported in `bam/{sample}.markdup.sorted.MarkDuplicates.metrics.txt` (`PERCENT_DUPLICATION` column), parsed automatically by MultiQC.

## Space-saving behavior

Pipeline removes some intermediates to reduce storage, for example:
- unsorted BAM after sort
- pre-filter BAM after filtering
- merged/trimmed FASTQ files in cleanup step after MultiQC

## Troubleshooting

- `MissingInputException` in MultiQC:
  - check module toggles and corresponding outputs
  - run dry-run first
- Rule env issues:
  - ensure `--use-conda`
  - delete broken env under `.snakemake/conda/` and rerun
- Large runs:
  - increase `--cores`
  - tune per-rule params in config (aligner/deeptools/featureCounts)
- Job killed / out of memory on HPC:
  - check the LSF job log (`bpeek JOB_ID` or `bhist -l JOB_ID`) to confirm out-of-memory (OOM) as the cause
  - **quick fix:** add or increase `mem_mb` for the failing rule in `workflow/profiles/lsf/config.yaml` under `set-resources` — this overrides the rule default without touching the code
  - **permanent fix:** if the rule's default in `workflow/modules/<rule>.smk` under `resources:` is too low, increase `mem_mb` there so the default itself is correct for all runs

## Acknowledgments

A huge thank you to Dr. Isabell Bludau, Dr.med.Abigail Suwala, Dr. Paul Kerbs, Quynh Nhu Nguyen and Temesvari-Nagy Levente from Heidelberg University Hospital and the German Cancer Research Center (DKFZ) for their support, feedback, and contributions to this pipeline.

## References

Key resources and prior work this pipeline draws from:

1. Niu Y. ATAC-seq data analysis: from FASTQ to peaks. Published March 20, 2019. https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/

2. Patel H, Espinosa-Carrasco J, Langer B, Ewels P, et al. nf-core/atacseq [v2.1.2]. Zenodo; 2022. https://nf-co.re/atacseq/2.1.2/

3. Yuan B. ATAC-seq Data Analysis. Presented at: BaRC Hot Topics; April 4, 2024; Whitehead Institute for Biomedical Research. http://barc.wi.mit.edu/education/hot_topics/ATACseq_2024/ATACseq2024_4slidesPerPage.pdf

4. Ou J, Liu H, Yu J, et al. ATACseqQC: a Bioconductor package for post-alignment quality assessment of ATAC-seq data. BMC Genomics. 2018;19(1):169. Published 2018 Mar 1. doi:10.1186/s12864-018-4559-3

## License

Follow the repository [MIT License](MIT_License.md) and tool licenses used in `workflow/envs/`.
