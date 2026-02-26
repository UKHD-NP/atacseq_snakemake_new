# ATAC-seq Snakemake Pipeline

A modular Snakemake pipeline for paired-end ATAC-seq from raw FASTQ to peaks, signal tracks, and QC reports.

This repository currently uses:
- `workflow/Snakefile`
- `config/config.yml`
- per-rule Conda environments in `workflow/envs/`

## Table of Contents

- [Pipeline Summary](#pipeline-summary)
- [Workflow DAG](#workflow-dag)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Files](#input-files)
  - [1. Samplesheet](#1-samplesheet)
  - [2. Reference](#2-reference)
- [Local Run](#local-run)
- [Running on DKFZ HPC (LSF)](#running-on-dkfz-hpc-lsf)
  - [Node roles at DKFZ](#node-roles-at-dkfz)
  - [Step 1 — Set up Snakemake environment](#step-1---set-up-snakemake-environment-on-odcf-worker01)
  - [Step 2 — Clone the pipeline](#step-2---clone-the-pipeline)
  - [Step 3 — Edit configuration](#step-3---edit-configuration)
  - [Step 4 — Update conda-prefix](#step-4---update-conda-prefix-in-the-lsf-profile)
  - [Step 5 — Validate with a dry-run](#step-5---validate-with-a-dry-run)
  - [Step 6 — Submit to HPC](#step-6---submit-to-hpc)
  - [Monitoring jobs](#monitoring-jobs)
- [Key Configuration](#key-configuration)
  - [Why these default params were chosen](#why-these-default-params-were-chosen)
  - [How to choose macs3_gsize](#how-to-choose-macs3_gsize)
- [Module Dependencies](#module-dependencies-important)
- [BAM Filtering Criteria](#re-mark-duplicates-and-bam-filtering-criteria)
  - [1) SAMtools core filter](#1-samtools-core-filter-from-configyml)
  - [2) BAMTools extra filter](#2-bamtools-extra-filter-optional-but-enabled-when-bamtools-exists)
  - [3) Pysam pair/orphan cleanup](#3-pysam-pairorphan-cleanup-optional-but-enabled-when-python3pysam-exists)
- [Main Outputs](#main-outputs)
- [MultiQC Content](#multiqc-content)
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
7. BAM stats on filtered BAM (`samtools stats/flagstat/idxstats`)
8. Signal tracks
   - Scaled bedGraph + bigWig from filtered BAM
   - ATAC-shifted BAM + shifted RPGC bigWig
9. Peak calling (MACS3 with Tn5-shifted BED)
10. Peak QC summary plots (`plot_macs2_qc.r`)
11. FRiP + MultiQC-ready FRiP / peak-count TSV
12. Peak annotation (HOMER + summary)
13. featureCounts in peaks (SAF)
14. deepTools matrix/profile/heatmap/fingerprint/bamPEFragmentSize
15. ataqv JSON + mkarv HTML report
16. MultiQC
17. Cleanup temporary/intermediate FASTQ files

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
        |                               |
        |                               +--> macs3_callpeak_tn5 (narrow: shifted.bam -> bamtobed -> MACS3 BED mode)
        |                               |                       (broad:  filtered.bam -> MACS3 BAMPE mode)
        |                               |    |
        |                               |    +--> frip_score (filtered.bam + peaks → MultiQC TSVs)
        |                               |    +--> annotate_peaks (optional)
        |                               |    +--> featurecounts_in_peaks (optional; narrow=shifted.bam, broad=filtered.bam)
        |                               |    +--> ataqv (optional; filtered.bam)
        |                               |
        |                               +--> deeptools (optional)
        |                                        +--> computeMatrix / plotProfile / plotHeatmap (shifted.bigWig)
        |                                        +--> plotFingerprint / bamPEFragmentSize (filtered.bam)
        |
        v
multiqc
        |
        v
delete_tmp
```

## Requirements

- Linux
- Snakemake ≥ 8 (in a dedicated controller environment)
- Conda/Mamba

> **DKFZ HPC users:** skip this section and follow [Running on DKFZ HPC (LSF)](#running-on-dkfz-hpc-lsf) instead, which covers environment setup outside your home directory.

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

> For cluster execution on DKFZ HPC, see [Running on DKFZ HPC (LSF)](#running-on-dkfz-hpc-lsf) below.
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
> Pass `--configfile path/to/other.yml` only when you want to override it.

## Running on DKFZ HPC (LSF)

The DKFZ cluster uses **IBM Spectrum LSF**.
A ready-made LSF profile is provided at `workflow/profiles/lsf/config.yaml`.

### Node roles at DKFZ

| Node | Purpose | Allowed |
|------|---------|---------|
| `odcf-worker01/02` | Dev, install, testing | ✅ Software install, small runs |
| `bsub01/02` | Job submission only | ✅ Run Snakemake (lightweight), ❌ Processing |
| Cluster nodes | Computation | Jobs submitted automatically via `bsub` |

### Step 1 - Set up Snakemake environment (on odcf-worker01)

> **Do this on `odcf-worker01`, not on `bsub01`.**
> Worker nodes (`odcf-worker01/02`) allow software installation. Submission hosts (`bsub01/02`) do not.

```bash
ssh YOUR_USERNAME@odcf-worker01.dkfz.de
```

**Configure conda channels.**
The DKFZ cluster bans the `defaults` (Anaconda) channel due to licensing restrictions.
You must explicitly restrict to `conda-forge` and `bioconda`:

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
module load Mamba/24.11.2-1
mamba init bash
source ~/.bashrc   # apply changes to the current shell without re-logging in
```

**Create the Snakemake controller environment outside your home directory.**
Home quota at DKFZ is only 20 GB. Conda environments can easily exceed this - install them on group storage:

```bash
# Set your working directory on group storage (adjust group/username as needed)
YOUR_WORKDIR="/omics/groups/OE0146/internal/YOUR_USERNAME"
mkdir -p ${YOUR_WORKDIR}/conda_envs

# Create the controller environment with Snakemake + the LSF executor plugin
mamba create -p ${YOUR_WORKDIR}/conda_envs/snakemake \
    -c conda-forge -c bioconda \
    snakemake \
    snakemake-executor-plugin-lsf \
    -y

# Activate the new environment
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

# Pin numpy/pandas to versions tested with this pipeline's helper scripts
python -m pip install "snakemake==8.*" "snakemake-executor-plugin-lsf" "numpy==1.26.4" "pandas==2.2.3"

# Verify that all three packages are importable and print their versions
python -c "import snakemake, numpy, pandas; print(snakemake.__version__, numpy.__version__, pandas.__version__)"
```

> `snakemake-executor-plugin-lsf` translates Snakemake rule resources (`mem_mb`, `runtime`, `threads`) into `bsub` submission flags automatically - no manual `bsub` scripting needed.

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
sed -i "s|/omics/odcf/analysis/YOUR_GROUP/conda_envs|${YOUR_WORKDIR}/conda_envs|g" \
    workflow/profiles/lsf/config.yaml

# Confirm the replacement was applied correctly
grep "conda-prefix" workflow/profiles/lsf/config.yaml
```

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

> **Do this on `bsub01` or `bsub02`**, not on `odcf-worker01`.
> Snakemake must run on a submission host to dispatch jobs via `bsub`.

Use `screen` so the Snakemake controller process survives SSH disconnects:

```bash
ssh YOUR_USERNAME@bsub01.lsf.dkfz.de

# Create a named screen session - it keeps running after SSH disconnect
screen -S atacseq

# Set your working directory (same value as used in Step 1)
YOUR_WORKDIR="/omics/groups/OE0146/internal/YOUR_USERNAME"

# Activate the Snakemake controller environment
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

# Move into the pipeline directory
cd ${YOUR_WORKDIR}/atacseq_snakemake_new

# Launch the pipeline - Snakemake submits each rule as a separate bsub job automatically.
# The config/config.yml is loaded automatically from the Snakefile; no --configfile needed.
# -j 100 allows up to 100 concurrent cluster jobs.
snakemake --profile workflow/profiles/lsf -j 100
```

To rerun only failed/incomplete jobs after fixing an error:

```bash
snakemake --profile workflow/profiles/lsf -j 100 --rerun-incomplete
```

To rerun with the test dataset config:

```bash
snakemake --profile workflow/profiles/lsf -j 100 --rerun-incomplete --configfile config/config_test.yml
```

| `screen` command | Action |
|-----------------|--------|
| `screen -S atacseq` | Start new named session |
| `Ctrl+A`, then `D` | Detach - session keeps running after SSH disconnect |
| `screen -ls` | List all active sessions |
| `screen -r atacseq` | Re-attach to session |

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
samples_csv: "test_data/samplesheet_test.csv"

ref:
  assembly: "custom"            # hg19 / hg38 / custom
  fasta: "test_data/references/genome.fa"
  gtf: "test_data/references/genes.gtf"
  blacklist: "test_data/references/ce11-blacklist.v2.bed"
  bwa_index: ""                 # optional prebuilt BWA prefix; auto-generated if empty
  bowtie2_index: ""             # optional prebuilt Bowtie2 prefix; auto-generated if empty
  mito_name: "MT"               # MUST match FASTA mito contig exactly (e.g. MT/chrM/M)
  keep_mito: false              # keep mitochondrial reads in include_regions or not

trimming:
  enabled: true
  tool: "trim_galore"           # fastp / trim_galore
  trim_galore_params: "--nextseq 20 --length 36"
  fastp_params: "--cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --trim_poly_g --length_required 36"

align:
  tool: "bowtie2"               # bwa / bowtie2
  bowtie2_params: "--very-sensitive --no-discordant -p 2 -X 2000"
  bwa_params: "-I 0,2000"

bam_filter:
  enabled: true
  params: "-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400 -q 30"

markduplicates:
  enabled: true

deeptools:
  enabled: true

call_peaks:
  enabled: true
  peak_type: "narrow"           # narrow / broad
  macs3_gsize: "2701262066"     # effective genome size (preferred); if empty, auto-sum from chromsizes
  macs3_narrow_params: "--trackline --shift -75 --extsize 150 --keep-dup all --nomodel --call-summits -q 0.01"
  macs3_broad_params: "--trackline --keep-dup all --nomodel --broad --broad-cutoff 0.1"
  macs3_peak_qc_plot: true      # run plot_macs_qc.r to generate peak QC summary/plots
  frip_overlap_fraction: 0.2

annotate_peaks:
  enabled: true

feature_counts:
  enabled: true
  use_shifted_bam: true         # true = count from shifted.bam, false = filtered.bam

ataqv:
  enabled: true

multiqc:
  config: "workflow/scripts/multiqc_config.yml"

latency-wait: 60
```

Explanation by block:
- `samples_csv`: input table for sample discovery and lane merging.
- `ref`: reference genome/annotation source. For `custom`, paths to `fasta`, `gtf`, and `blacklist` are required. `bwa_index` / `bowtie2_index` can be left empty to auto-build.
- `ref.mito_name`: **critical** — must exactly match the mitochondrial contig name in your FASTA (often `MT`, `chrM`, or `M`). Used by `bam_filter` to build `include_regions` and by `ataqv`.
- `ref.keep_mito`: set `true` to retain mitochondrial reads in `include_regions`; `false` (default) excludes them.
- `trimming`: choose one trimming engine and pass tool-specific options.
- `align`: choose aligner and set aligner-specific CLI parameters.
- `bam_filter.enabled`: disabling this skips all downstream modules that depend on `filtered.bam` (peaks, deeptools, align_stats). See [Module Dependencies](#module-dependencies-important).
- `bam_filter.params`: SAMtools core filter flags; see [BAM Filtering](#re-mark-duplicates-and-bam-filtering-criteria) for full breakdown.
- `markduplicates.enabled`: run Picard MarkDuplicates before filtering. When disabled, duplicates are not flagged and `-F 0x0400` in `bam_filter.params` has no effect.
- `deeptools.enabled`: run computeMatrix/plotProfile/plotHeatmap/plotFingerprint modules.
- `call_peaks.peak_type`: `narrow` uses Tn5-shifted BAM → BED → MACS3 BED mode; `broad` uses filtered BAM in MACS3 BAMPE mode.
- `call_peaks.macs3_peak_qc_plot`: when `true`, runs `plot_macs_qc.r` to produce `*.macs_peakqc.summary.txt` and `*.macs_peakqc.plots.pdf`.
- `call_peaks.frip_overlap_fraction`: minimum overlap fraction for FRiP calculation (passed to bedtools intersect).
- `annotate_peaks.enabled`: run HOMER `annotatePeaks` and summary plotting.
- `feature_counts.use_shifted_bam`: `true` = count reads from `shifted.bam` (recommended for narrow peaks); `false` = use `filtered.bam`.
- `ataqv.enabled`: run ATAC-specific QC (`ataqv`) and render interactive HTML (`mkarv`). Requires `call_peaks.enabled=true`.
- `multiqc.config`: path to MultiQC config used by this pipeline.
- `latency-wait`: useful on slow/network filesystems to avoid false missing-output errors.

### Why these default params were chosen

- `trimming.tool: trim_galore` with `--nextseq 20 --length 36`:
  chosen for two-color Illumina runs (poly-G prone) and to remove very short reads that are usually uninformative for peak calling.
- `align.tool: bowtie2` with `--very-sensitive --no-discordant -X 2000`:
  chosen to maximize paired-end sensitivity while constraining improbable pair structure for ATAC fragment lengths.
- `bam_filter.params: ... -q 30`:
  chosen as a relatively strict default to retain high-confidence alignments for downstream peak calling and signal tracks.
- `call_peaks` defaults (`--shift -75 --extsize 150 --nomodel`):
  negative shift centers 150 bp windows on Tn5 cut sites (reads are already Tn5-shifted by `alignmentSieve --ATACshift`; `--shift -75` pulls each cut-site tag upstream by 75 bp so the 150 bp extension lands symmetrically around the insertion site). Positive shift would offset windows 75 bp away from the cut site.
- `feature_counts.use_shifted_bam: true`:
  chosen to count reads against peaks using the same ATAC-shifted representation used by deepTools signal modules.

### How to choose `macs3_gsize`

- `macs3_gsize` is passed to MACS3 as `--gsize`.
- Prefer effective genome size values from deepTools documentation:  
  `https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html`
- If `macs3_gsize` is empty, this pipeline falls back to summing chromosome sizes from `ref.chromsizes`.
- Common values (from deepTools table):
  - `hg19`: `2864785220`
  - `hg38`: `2913022398`

## Module Dependencies (important)

- `call_peaks` requires `bam_filter.enabled=true`.
- `deeptools` requires `bam_filter.enabled=true`.
- `align_stats` requires `bam_filter.enabled=true`.
- `annotate_peaks` requires `call_peaks.enabled=true`.
- `feature_counts` requires `call_peaks.enabled=true`.
- `ataqv` requires `call_peaks.enabled=true`.

Disabling `bam_filter` effectively disables all downstream analysis — only raw QC (FastQC, trimming) will run.

## Re-mark Duplicates and BAM Filtering Criteria

After alignment and coordinate sorting, duplicates are marked with Picard (`mark_duplicates`, when enabled), then `bam_filter` creates `*.filtered.bam`.

Important: `bam_filter.params` is the **SAMtools core filter only**.  
Additional filters are applied by BAMTools and Pysam (if available in the runtime env).

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
- this keeps only mappable regions outside blacklist intervals
- mitochondrial contig can be excluded there when `ref.keep_mito: false`

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

Why this default filtering strategy is used:
- It prioritizes specificity for ATAC peak detection.
- It reduces noisy alignments before MACS3, FRiP, and bigWig generation.
- It is a practical strict default (`-q 30`); for low-depth data you can relax to `-q 20`.

## Main Outputs

Per sample under `<outdir>`:

- BAM / stats
  - `bam/{sample}.filtered.bam`
  - `bam/{sample}.filtered.bam.bai`
  - `bam/{sample}.filtered.bam.stats`
  - `bam/{sample}.filtered.bam.flagstat`
  - `bam/{sample}.filtered.bam.idxstats`
- Peaks / FRiP
  - `peaks/{sample}_peaks.peak`
  - `peaks/{sample}_peaks.xls`
  - `peaks/{sample}.FRiP.txt`
  - `peaks/{sample}_peaks.FRiP_mqc.tsv`
  - `peaks/{sample}_peaks.count_mqc.tsv`
  - `peaks/{sample}.macs_peakqc.summary.txt`
  - `peaks/{sample}.macs_peakqc.plots.pdf`
- Annotation
  - `annotation/{sample}_peaks.annotatePeaks.txt`
  - `annotation/{sample}.macs_annotatePeaks.summary.txt`
- Peak counts
  - `featurecounts/{sample}_peaks.saf`
  - `featurecounts/{sample}.readCountInPeaks.txt`
  - `featurecounts/{sample}.readCountInPeaks.txt.summary`
- Signal tracks
  - `bigwig/{sample}.bigWig`
  - `bam/{sample}.shifted.bam`
  - `bam/{sample}.shifted.bam.bai`
  - `bigwig/{sample}.shifted.bigWig`
- deepTools
  - `deeptools/*.computeMatrix.*`
  - `deeptools/*.plotProfile.*`
  - `deeptools/*.plotHeatmap.*`
  - `deeptools/*.plotFingerprint.*`
  - `deeptools/{sample}.fragment_size_distribution.pdf`
  - `deeptools/{sample}.fragment_size.raw_lengths.txt`
  - `deeptools/{sample}.fragment_size.qcmetrics.txt`
- ataqv
  - `ataqv/{sample}.ataqv.json`
  - `ataqv/{sample}.mkarv_html/index.html`
- MultiQC
  - `multiqc/{sample}.multiqc.html`

## MultiQC Content

`workflow/scripts/multiqc_config.yml` is configured for:
- Raw and trimmed read QC
- Samtools filtered-BAM stats
- Picard MarkDuplicates
- deepTools QC metrics
- ataqv JSON
- FRiP table (`*_peaks.FRiP_mqc.tsv`)
- Peak count table (`*_peaks.count_mqc.tsv`)
- HOMER annotation summary
- MACS peak QC summary

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

Key metrics from ENCODE (practical targets and acceptable ranges):
- Alignment rate: target `>95%` (acceptable `>80%`)
- FRiP score: target `>0.3` (acceptable `>0.2`)
- TSS enrichment: target `>10` (acceptable `>5`)
- Library complexity (NRF): target `>0.9`
- PBC1: target `>0.9`
- PBC2: target `>3`
- Nucleosome-free regions: should be detectable in called peaks

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

## Acknowledgments

A huge thank you to Dr. Isabell Bludau, Dr.med.Abigail Suwala, Dr. Paul Kerbs, Quynh Nhu Nguyen and Temesvari-Nagy Levente from Heidelberg University Hospital and the German Cancer Research Center (DKFZ) for their support, feedback, and contributions to this pipeline.

## References

Key resources and prior work this pipeline draws from:

1. Niu Y. ATAC-seq data analysis: from FASTQ to peaks. Published March 20, 2019. https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/

2. Patel H, Espinosa-Carrasco J, Langer B, Ewels P, et al. nf-core/atacseq [v2.1.2]. Zenodo; 2022. https://nf-co.re/atacseq/2.1.2/

3. Yuan B. ATAC-seq Data Analysis. Presented at: BaRC Hot Topics; April 4, 2024; Whitehead Institute for Biomedical Research. http://barc.wi.mit.edu/education/hot_topics/ATACseq_2024/ATACseq2024_4slidesPerPage.pdf

## License

Follow the repository [MIT License](MIT_License.md) and tool licenses used in `workflow/envs/`.
