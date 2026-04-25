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
13. featureCounts in peaks (SAF)
14. deepTools matrix/profile/heatmap/fingerprint/bamPEFragmentSize
15. NFR analysis — nucleosome-free vs mononucleosomal bigWigs + TSS profile/heatmap
16. ataqv JSON + mkarv HTML report + TSS enrichment / NFR metrics table for MultiQC
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
        |                               +--> frip_score (filtered.bam + peaks → MultiQC TSVs)
        |                               +--> annotate_peaks (optional)
        |                               +--> featurecounts_in_peaks (optional; always filtered.bam)
        |                               +--> ataqv (optional; filtered.bam)
        |                               |    +--> atac_qc_metrics → atac_qc_mqc.tsv (TSS enrichment, NFR ratio)
        |                               +--> deeptools (optional)
        |                               |    +--> computeMatrix / plotProfile / plotHeatmap (shifted.bigWig for narrow peaks; bigWig for broad peaks)
        |                               |    +--> plotFingerprint / bamPEFragmentSize (filtered.bam)
        |                               |    +--> pt_score (promoter/transcript body ratio → pt_score_mqc.tsv)
        |                               +--> nfr (optional; when deeptools + narrow peaks)
        |                                    +--> alignmentSieve → nfr.bigWig + mono.bigWig
        |                                    +--> computeMatrix / plotProfile / plotHeatmap (nfr vs mono)
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

deeptools:
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

ataqv:
  enabled: true

# Optional: tune NFR/mononucleosomal fragment size boundaries.
# NFR analysis runs automatically when deeptools + narrow peaks are enabled.
# nfr:
#   nfr_max_fragment: 150      # fragments ≤ this bp → NFR bigWig
#   mono_min_fragment: 151     # fragments ≥ this bp → mono bigWig
#   mono_max_fragment: 300     # fragments ≤ this bp → mono bigWig

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
- `deeptools.enabled`: run computeMatrix/plotProfile/plotHeatmap/plotFingerprint/PT score modules. Requires `call_peaks.enabled=true`. For narrow peaks, computeMatrix uses the Tn5-shifted bigWig (`shifted.bigWig`); for broad peaks, it uses the unshifted bigWig (`bigWig`).
- `call_peaks.peak_type`: `narrow` uses filtered BAM → `bamtobed` → awk Tn5 shift (+4 forward / -5 reverse) → MACS3 BED mode; `broad` uses filtered BAM directly in MACS3 BAMPE mode.
- `call_peaks.macs3_peak_qc_plot`: when `true`, runs `plot_macs_qc.r` to produce `*.macs_peakqc.summary.txt` and `*.macs_peakqc.plots.pdf`.
- `call_peaks.frip_overlap_fraction`: minimum read-peak overlap fraction for FRiP (passed to both `bedtools intersect -f` and featureCounts `--fracOverlap`).
- `call_peaks.frip_threshold`: FRiP percentage threshold for quality label in `*.FRiP.txt`; samples at or above this value are labelled `good`, below is `bad` (default: 20%).
- `annotate_peaks.enabled`: run HOMER `annotatePeaks` and summary plotting.
- `ataqv.enabled`: run ATAC-specific QC (`ataqv`) and render interactive HTML (`mkarv`), plus extract TSS enrichment score and NFR ratio to `*.atac_qc_mqc.tsv` for MultiQC. Requires `call_peaks.enabled=true`.
- `nfr` (optional block): tune fragment size cutoffs for NFR vs mononucleosomal analysis. Omit to use defaults (NFR ≤150 bp, mono 151–300 bp).
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

## Module Dependencies (IMPORTANT)

Some modules only run when their upstream module is also enabled:

- `align_stats` requires `bam_filter` (always runs — `bam_filter` has no toggle).
- `call_peaks` requires `bam_filter`.
- `annotate_peaks` requires `call_peaks.enabled=true`.
- `featureCounts` + `frip_score` run automatically when `call_peaks.enabled=true` (no separate toggle).
- `ataqv` requires `call_peaks.enabled=true`.
- `deeptools` requires `call_peaks.enabled=true`. Works with both `narrow` and `broad` peak types — computeMatrix automatically uses the Tn5-shifted bigWig for narrow peaks and the unshifted bigWig for broad peaks.
- `nfr` analysis runs automatically when `deeptools.enabled=true` AND `call_peaks.enabled=true` AND `call_peaks.peak_type=narrow`. Fragment size cutoffs can be tuned via the optional `nfr:` config block (defaults: NFR ≤150 bp, mono 151–300 bp).
- `atac_qc_metrics` (TSS enrichment + NFR ratio table for MultiQC) requires `ataqv.enabled=true`.

> **Note:** `bam_filter` itself has no `enabled` toggle — it always runs. Disabling `markduplicates` is safe but leaves duplicates unflagged (the `-F 0x0400` flag in `bam_filter.params` would then have no effect).

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

Why this default filtering strategy is used:
- It prioritizes specificity for ATAC peak detection.
- It reduces noisy alignments before MACS3, FRiP, and bigWig generation.
- It is a practical strict default (`-q 30`); for low-depth data you can relax to `-q 20`.

Practical examples:
- Keep mitochondrial filtering, but disable blacklist filtering:
  `ref.keep_mito: false` and `bam_filter.apply_blacklist: false`
- Preserve the BAM before filtering for debugging or alternate peak-calling runs:
  `bam_filter.keep_input_bam: true`

## Main Outputs

Per sample under `<outdir>`:

- BAM / stats
  - `bam/{sample}.filtered.bam`
  - `bam/{sample}.filtered.bam.bai`
  - `bam/{sample}.filtered.bam.stats`
  - `bam/{sample}.filtered.bam.flagstat`
  - `bam/{sample}.filtered.bam.idxstats`
  - `bam/{sample}.markdup.sorted.bam` or `bam/{sample}.bam` may also be retained when `bam_filter.keep_input_bam: true`
  - `bam/{sample}.CollectMultipleMetrics.alignment_summary_metrics`
  - `bam/{sample}.CollectMultipleMetrics.base_distribution_by_cycle.pdf`
  - `bam/{sample}.CollectMultipleMetrics.base_distribution_by_cycle_metrics`
  - `bam/{sample}.CollectMultipleMetrics.insert_size_histogram.pdf`
  - `bam/{sample}.CollectMultipleMetrics.insert_size_metrics`
  - `bam/{sample}.CollectMultipleMetrics.quality_by_cycle.pdf`
  - `bam/{sample}.CollectMultipleMetrics.quality_by_cycle_metrics`
  - `bam/{sample}.CollectMultipleMetrics.quality_distribution.pdf`
  - `bam/{sample}.CollectMultipleMetrics.quality_distribution_metrics`
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
  - `deeptools/{sample}.pt_score_mqc.tsv` (PT score: promoter/transcript body signal ratio)
- ataqv
  - `ataqv/{sample}.ataqv.json`
  - `ataqv/{sample}.mkarv_html/index.html`
  - `ataqv/{sample}.atac_qc_mqc.tsv` (TSS enrichment score + NFR ratio for MultiQC)
- NFR analysis
  - `nfr/{sample}.nfr.bigWig` (fragments ≤150 bp)
  - `nfr/{sample}.mono.bigWig` (fragments 151–300 bp)
  - `nfr/{sample}.nfr_vs_mono.computeMatrix.mat.gz`
  - `nfr/{sample}.nfr_vs_mono.plotProfile.pdf` + `.tab`
  - `nfr/{sample}.nfr_vs_mono.plotHeatmap.pdf`
- MultiQC
  - `multiqc/{sample}.multiqc.html`

## MultiQC Content

`workflow/scripts/multiqc_config.yml` is configured for:
- Raw and trimmed read QC
- Samtools filtered-BAM stats
- Picard MarkDuplicates
- Picard CollectMultipleMetrics
  - alignment summary
  - insert size
  - base distribution by cycle
  - quality by cycle
  - quality distribution
- deepTools QC metrics (fingerprint, profiles, fragment size)
- deepTools PT score (promoter/transcript body signal ratio, `*.pt_score_mqc.tsv`)
- deepTools NFR vs mononucleosomal TSS profile (from `nfr/` analysis)
- ataqv JSON
- TSS Enrichment Score + NFR ratio per sample (`*.atac_qc_mqc.tsv`)
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

### Quick reference

| Metric | Tool / Output | Target | Acceptable |
|--------|--------------|--------|------------|
| Alignment rate | Bowtie2 / BWA-MEM2 log | >95% | >80% |
| Duplication rate | `*.MarkDuplicates.metrics.txt` | <20% | <30% |
| FRiP score | `*_peaks.FRiP_mqc.tsv` | >0.3 | >0.2 |
| TSS enrichment | `*.atac_qc_mqc.tsv` (ataqv) | >10 | >5 |
| PT score (2^mean) | `*.pt_score_mqc.tsv` | ≥10 | ≥5 |
| NFR ratio | `*.atac_qc_mqc.tsv` (ataqv) | >48% | >40% |

---

### FRiP Score — Fraction of Reads in Peaks

**What it measures:** The proportion of all uniquely mapped reads that fall within called peak regions.

```
FRiP = reads_overlapping_peaks / total_mapped_reads
```

**Why it matters:** A high FRiP means most of your sequencing reads captured genuine open chromatin sites, rather than noisy background. Low FRiP suggests poor enrichment, over-amplification, or degraded nuclei.

**Thresholds:**
- `>0.3` — ENCODE target; well-enriched library
- `0.2–0.3` — acceptable; borderline, check other metrics
- `<0.2` — poor enrichment; inspect fragment size distribution and TSS enrichment

**In this pipeline:** two FRiP estimates are reported per sample — one from `bedtools intersect` and one from `featureCounts`. Both appear in MultiQC.

---

### TSS Enrichment Score

**What it measures:** Signal enrichment at Transcription Start Sites (TSS) relative to the flanking background. Calculated by ataqv.

```
TSS enrichment = mean signal in TSS ±150 bp window
                 ─────────────────────────────────
                 mean signal in flanking regions (±1900 bp from TSS)
```

**Why it matters:** Tn5 transposase preferentially inserts at open chromatin, which in ATAC-seq is concentrated at active promoters. A high TSS enrichment score confirms the experiment successfully captured nucleosome-free promoter regions. A flat or low score means the library has high background, poor Tn5 enrichment, or degraded/fixed chromatin.

**Thresholds:**
- `>10` — ENCODE target; excellent signal-to-noise
- `5–10` — acceptable; usable but noisier peak calls
- `<5` — poor; re-check nuclei isolation and Tn5 titration

**In this pipeline:** reported in `ataqv/{sample}.atac_qc_mqc.tsv`, shown in MultiQC.

---

### PT Score — Promoter/Transcript Body ratio

**What it measures:** Whether Tn5 insertions (5′ read ends) are enriched at promoters relative to gene bodies. Calculated by the ATACseqQC R package.

```
For each gene:
  promoter_window  = [TSS−2000, TSS+500]   (strand-aware)
  body_window      = next 2500 bp downstream of promoter

  PT score (log2) = log2(mean_5prime_density_in_promoter + ε)
                  − log2(mean_5prime_density_in_body + ε)

Final: mean and median PT score across all genes
```

The mean PT score is reported in log2 scale; the equivalent linear ratio is **2^PT_score_mean**.

**Why it matters:** In ATAC-seq, open chromatin is concentrated at promoters. A high PT score confirms that most signal comes from promoter-proximal nucleosome-free regions rather than evenly distributed gene body background. Low PT scores can indicate high background, poor Tn5 enrichment, or a ChIP-like signal profile.

**Thresholds (linear scale, i.e. 2^mean_PT):**
- `≥10` — excellent; strong promoter-over-body enrichment
- `5–10` — acceptable; ATACseqQC default PASS threshold is `≥5`
- `<5` — borderline/FAIL; interpret with other metrics

**In this pipeline:** reported in `deeptools/{sample}.pt_score_mqc.tsv`, shown in MultiQC. A `[INFO] QC: PASS` or `[WARNING] QC: FAIL` line is written to the rule log.

---

### Duplication Rate — Library Complexity

**What it measures:** The fraction of reads that are PCR/optical duplicates, identified by Picard MarkDuplicates.

```
Duplication rate = duplicate reads / total mapped reads
```

**Why it matters:** High duplication indicates the library was over-amplified or sequenced to saturation — most reads are copies of the same original fragment rather than independent Tn5 insertions. This reduces the effective complexity of the library and inflates peak signal artificially.

**Thresholds:**
- `<20%` — good; library has sufficient complexity
- `20–30%` — acceptable; borderline, check FRiP and fragment size distribution
- `>30%` — high; PCR cycles may need to be reduced in the next experiment

**In this pipeline:** reported in `bam/{sample}.markdup.sorted.MarkDuplicates.metrics.txt` (`PERCENT_DUPLICATION` column), parsed automatically by MultiQC.

> **Note:** NRF, PBC1, and PBC2 (ENCODE complexity metrics) are not calculated by this pipeline. They require a separate read-position counting step not included here. Duplication rate from Picard provides an equivalent practical assessment of library complexity.

---

### NFR Ratio — Nucleosome-Free Region Fraction

**What it measures:** The fraction of fragments that are short (≤150 bp), corresponding to sub-nucleosomal / nucleosome-free insertions.

**Why it matters:** A healthy ATAC-seq library should show a prominent short-fragment peak (<200 bp) in the fragment size distribution, corresponding to Tn5 insertions in nucleosome-free open chromatin. The NFR ratio quantifies how much of the library is in this fraction. Libraries with high mononucleosomal or dinucleosomal contamination will show lower NFR ratios and broader, noisier peaks.

**Thresholds:**
- `>48%` — ENCODE target
- `40–48%` — acceptable
- `<40%` — poor NFR enrichment; check Tn5 titration and nuclei quality

**In this pipeline:** calculated by ataqv, reported in `ataqv/{sample}.atac_qc_mqc.tsv`, shown in MultiQC.

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

## License

Follow the repository [MIT License](MIT_License.md) and tool licenses used in `workflow/envs/`.
