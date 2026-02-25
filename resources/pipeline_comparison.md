# Comparison: OLD vs NEW ATAC-seq Pipeline

## 1. Architecture

| | OLD | NEW |
|---|---|---|
| Structure | Single monolithic Snakefile | 20 modular `.smk` files |
| Configuration | Hardcoded in rules | YAML-driven, each module toggleable on/off |
| References | Manual setup | Auto-staged for hg19/hg38; fully `custom` mode supported |

---

## 2. Step-by-step comparison

| Step | OLD | NEW | Key difference |
|------|-----|-----|----------------|
| **Trimming** | cutadapt | **trim_galore (default)** or fastp | trim_galore is now default (NextSeq poly-G aware); fastp available as alternative |
| **FastQC trimmed** | Not present | `fastqc_trimmed` run only in `trim_galore` mode | fastp mode produces its own HTML/JSON report instead |
| **Alignment** | bwa-mem2 only | **bowtie2 (default)** or bwa-mem2 | bowtie2 is now default (`--very-sensitive --no-discordant -X 2000`); both aligners add RG tags |
| **Filtering** | samtools `-e '([NM]<=4) && sclen<15'` + `-F 0x100 -f 3` | bamtools JSON filter + `-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400 -q 30` | NEW adds orphan removal, uses include-regions BED instead of `grep chrM` |
| **Mark dup** | GATK `MarkDuplicatesWithMateCigar` (removes dups) | Picard `MarkDuplicates` (marks only; removed at filter step via `-F 0x0400`) | Different order: OLD filters then marks; NEW marks then filters (more correct) |
| **Tn5 shift** | 1 massive rule: BAM→BED→shift→BEDPE→shift again→blacklist→RPM→bigwig | 2 separate mechanisms: (1) awk shift BED for peak calling, (2) `alignmentSieve --ATACshift` → `shifted.bam` for BAM-based signal | NEW cleanly separates peak calling from signal track generation |
| **BigWig** | 1 RPM track from shifted data | 2 tracks: RPGC (`shifted.bigWig`) + unscaled bedGraph→bigWig (`bigWig`) | RPGC is more appropriate for cross-sample comparison |
| **Peak calling** | MACS3 on shifted BED (narrow only) | MACS3, **narrow or broad** mode configurable | NEW adds broad mode (`filtered.bam` BAMPE), peak QC plots, FRiP TSVs |
| **Peak annotation** | Not present | HOMER `annotatePeaks` + summary plot | NEW |
| **featureCounts** | Not present | `featurecounts_in_peaks` (SAF-based); configurable shifted vs filtered BAM | NEW; `use_shifted_bam: true` by default |
| **ataqv** | Not present | ataqv JSON + mkarv interactive HTML | NEW |
| **MultiQC** | Not present | Full MultiQC report aggregating all QC modules | NEW |
| **Cleanup** | Not present | `delete_tmp` removes merged FASTQ, pre-filter BAM, unsorted BAM | NEW |

---

## 3. Features in OLD but MISSING in NEW

| Feature | Notes |
|---------|-------|
| **Fragment size distribution** (`bamPEFragmentSize`) | No dedicated rule in NEW; partially covered by ataqv and deepTools |
| **PBC calculation** (PCR Bottleneck Coefficient) | Completely absent — recommended insertion point: after `mark_duplicates`, before `bam_filter` (see Section 8) |
| **Merged Lorenz curve plots** | NEW only has per-sample fingerprint; no cross-sample merge |
| **FastQC on BAMs** | Not present in NEW |
| **Soft-clip filter** (`sclen < 15`) | Not in NEW — bamtools JSON uses `!cigar *S*` (any soft-clip rejected, too aggressive; see Section 7) |
| **BEDPE-based shifting** | NEW only shifts single-end BED for MACS3; BAM shifting via `--ATACshift` |

---

## 4. Features in NEW but MISSING in OLD

| Feature | Notes |
|---------|-------|
| **Multi-lane FASTQ merging** | Multiple sequencing lanes per `sample_id` supported via samplesheet |
| **FastQC on raw reads** | Pre-trimming QC (`fastqc_raw`) |
| **FastQC on trimmed reads** | Post-trimming QC (`fastqc_trimmed`; trim_galore mode only) |
| **Bowtie2 support** | Alternative aligner, now the default |
| **HOMER peak annotation** | Annotates peaks to genomic features with summary plot |
| **featureCounts quantification** | Read counting in peaks (SAF format); configurable BAM source |
| **ataqv** | ATAC-seq-specific QC + interactive mkarv HTML report |
| **deepTools heatmap / plotProfile** | TSS-centered heatmap, gene-body profiles (shifted.bigWig) |
| **deepTools plotFingerprint** | Lorenz curve from filtered.bam |
| **Auto include-regions** | Genome minus blacklist, optionally minus mito (`ref.keep_mito`) |
| **Mito name configurable** | `ref.mito_name` — must match FASTA contig exactly (MT / chrM / M) |
| **Broad peak mode** | `call_peaks.peak_type: broad` uses filtered.bam in BAMPE mode |
| **Disk cleanup** | Auto-deletes unsorted BAM, pre-filter BAM, merged/trimmed FASTQ after MultiQC |
| **FRiP TSVs for MultiQC** | `*_peaks.FRiP_mqc.tsv` + `*_peaks.count_mqc.tsv` picked up by MultiQC |

---

## 5. Pipeline order

**OLD:** Trim → Align → **Filter → Mark dup** → Shift+BigWig → Peak call → QC

**NEW:** Trim → Align → Sort → **Mark dup → Filter** → Shift BAM → BigWig → Peak call → Annotation → QC → MultiQC → Cleanup

> **Important:** OLD filters before marking duplicates. NEW marks duplicates first, then filters.
> The NEW order is more correct because duplicate detection needs to see all reads (including potential duplicates) to accurately identify them.

---

## 6. Detailed parameter comparison

| Parameter | OLD | NEW |
|-----------|-----|-----|
| Default trimmer | cutadapt | **trim_galore** (`--nextseq 20 --length 36`) |
| Default aligner | bwa-mem2 | **bowtie2** (`--very-sensitive --no-discordant -p 2 -X 2000`) |
| MAPQ threshold | Configurable `-q MAPQ` | `-q 30` (configurable via `bam_filter.params`) |
| SAM flag filters | `-F 0x100 -f 3 -F 0x0008` | `-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400 -q 30` |
| NM mismatch filter | `[NM] <= 4` (samtools expression) | `NM:<=4` via bamtools JSON |
| Soft-clip filter | `sclen < 15` (allows small soft-clips) | `!cigar *S*` (rejects ANY soft-clip — too aggressive; see Section 7) |
| Insert size filter | Not filtered | `insertSize >= -2000 && <= 2000` via bamtools JSON |
| Tn5 shift values | +4 (forward) / -5 (reverse) | +4/-5 for BED (peak calling); `--ATACshift` for BAM (signal tracks) |
| MACS3 narrow params | Not specified | `--shift 75 --extsize 150 --keep-dup all --nomodel --call-summits -q 0.01` |
| MACS3 broad params | Not present | `--keep-dup all --nomodel --broad --broad-cutoff 0.1` |
| Effective genome size | Manually configured | Configurable via `macs3_gsize`; auto-sum from chromsizes if empty |
| FRiP overlap threshold | Not present | `frip_overlap_fraction: 0.2` |
| featureCounts BAM source | Not present | `feature_counts.use_shifted_bam: true` (shifted.bam by default) |
| Normalization | RPM | RPGC (shifted.bigWig) + unscaled bedGraph→bigWig |
| Duplicate handling | GATK MarkDuplicatesWithMateCigar (removes) | Picard MarkDuplicates (marks only; `-F 0x0400` removes downstream) |

---

## 7. BAM filtering comparison

### Filter-by-filter breakdown

| Filter | OLD (`samtools -e`) | NEW (`samtools flags` + `bamtools JSON`) | Winner |
|--------|---------------------|------------------------------------------|--------|
| Unmapped reads | Not explicitly filtered | `-F 0x004` (remove unmapped) | NEW |
| Mate unmapped | `-F 0x0008` | `-F 0x0008` | Same |
| Proper pair | `-f 3` (mapped + proper pair) | `-f 0x001` (paired only) + `bampe_rm_orphan.py --only_fr_pairs` (same-chr, FR orientation) | NEW (stricter in practice) |
| Secondary alignments | `-F 0x100` | `-F 0x0100` | Same |
| Duplicates | Not filtered (handled separately) | `-F 0x0400` (remove dup-flagged reads) | NEW |
| MAPQ | `-q MAPQ` (configurable) | `-q 30` (configurable via `bam_filter.params`) | NEW (stricter default) |
| NM mismatch | `[NM] <= 4` (samtools expression) | `NM:<=4` via bamtools JSON | Same (threshold identical) |
| Soft-clip | `sclen < 15` (allows small soft-clips) | `!cigar *S*` (rejects ANY soft-clip) | OLD (NEW is too aggressive) |
| Insert size | Not filtered | `insertSize >= -2000 && <= 2000` via bamtools | NEW |
| MT / blacklist | `grep -v chrM` | `-L include_regions` BED (controlled by `ref.keep_mito` + `ref.mito_name`) | NEW (more flexible) |
| Orphan reads | Not handled | `bampe_rm_orphan.py --only_fr_pairs` | NEW |

### How NEW pipeline handles proper-pair filtering

The OLD pipeline uses `-f 3` (samtools proper-pair flag) which trusts the aligner's definition.

The NEW pipeline takes a more thorough approach via `bampe_rm_orphan.py --only_fr_pairs`:
1. Removes singleton/orphan reads (no mate in BAM)
2. Keeps only pairs on the **same chromosome** (`pair1.tid == pair2.tid`)
3. Keeps only **FR orientation** (read1 forward + read2 reverse, or vice versa, with correct relative positioning)
4. Removes all improper pairs (wrong orientation, interchromosomal, etc.)

Combined with the bamtools insert-size filter (`-2000 to 2000`), this is actually **stricter** than the aligner's `-f 0x002` flag because it verifies orientation and positioning at the read level.

### Remaining issue in NEW pipeline

**Soft-clip filter too aggressive:** `!cigar *S*` in the bamtools JSON rejects any read with even 1 bp of soft-clipping. Most aligners produce some soft-clipping at read ends. This can discard 10–30% of good reads. The OLD approach (`sclen < 15`) is more reasonable — it allows small soft-clips while removing reads with excessive soft-clipping.

### Recommended fix for soft-clip filter

**bamtools JSON** (remove overly aggressive soft-clip filter):
```json
{
    "filters": [
        { "id": "insert_min", "insertSize": ">=-2000" },
        { "id": "insert_max", "insertSize": "<=2000" },
        { "id": "mismatch", "tag": "NM:<=4" }
    ],
    "rule": "insert_min & insert_max & mismatch"
}
```

Two valid approaches for soft-clip filtering:
- **Option A**: Drop it entirely — NM + MAPQ filters already catch most problematic alignments.
- **Option B**: Add `samtools view -e 'sclen < 15'` as an additional pipe step to match OLD behaviour.

---

## 8. PBC Calculation — where to add in NEW pipeline

### Why PBC is not replaced by ataqv

| Metric | What it measures | Tool |
|--------|-----------------|------|
| PBC1, PBC2 | Library complexity / PCR over-amplification (coordinate-level read stacking) | Must be computed separately |
| Duplicate rate | Flag-based duplicate fraction | Picard MarkDuplicates (already in NEW) |
| TSS enrichment | Chromatin accessibility quality | ataqv |
| Fragment size distribution | Nucleosome occupancy | ataqv / deepTools |

**ataqv does NOT compute PBC.** They answer different questions:
- PBC = how much of the library is PCR artifact?
- ataqv = how good is the ATAC signal?

### Recommended insertion point

```
Trim → Align → Sort → Mark dup → [PBC HERE] → Filter → ...
```

**Input:** marked BAM (output of Picard MarkDuplicates — duplicates flagged but NOT yet removed)
**Output:** `{sample}.pbc.txt`

Rationale:
- Computing BEFORE duplicate removal: all reads including duplicates are present for position-stacking counts.
- Computing AFTER duplicate removal (`-F 0x0400`): M1/M2 counts will be wrong — information is lost.
- Computing BEFORE MarkDup: duplicates not yet identified — cannot separate PCR from natural overlaps.

### PBC computation (shell)

From the **marked BAM**, apply basic read filters (unmapped/secondary/MAPQ) but NOT `-F 0x0400`:

```bash
samtools view -F 0x004 -F 0x0008 -F 0x0100 -f 0x001 -q 30 {sample}.marked.bam \
  | awk '{print $3"\t"$4}' \
  | sort \
  | uniq -c \
  | awk '
    BEGIN { Mt=0; M1=0; M2=0 }
    {
      Mt += $1
      if ($1 == 1) M1++
      if ($1 == 2) M2++
    }
    END {
      printf "PBC1\t%.4f\n", M1/Mt
      printf "PBC2\t%.4f\n", (M2 > 0 ? M1/M2 : "NA")
      printf "NRF\t%.4f\n",  (NR > 0 ? NR/Mt : "NA")
    }' > {sample}.pbc.txt
```

### ENCODE thresholds

| Metric | Severe bottleneck | Moderate | Ideal |
|--------|------------------|----------|-------|
| PBC1 | < 0.5 | 0.5–0.8 | > 0.9 |
| PBC2 | < 1 | 1–3 | > 3 |
| NRF | < 0.5 | 0.5–0.8 | > 0.9 |

---

## 9. Running on HPC (DKFZ LSF cluster)

### OLD vs NEW — HPC support comparison

| | OLD | NEW |
|---|---|---|
| Cluster scheduler | Manual `bsub` scripts or none | `snakemake-executor-plugin-lsf` (automatic) |
| LSF profile | Not provided | `workflow/profiles/lsf/config.yaml` |
| Resource declaration | Not per-rule | Per-rule: `mem_mb`, `runtime`, `threads` → auto-translated to `bsub` flags |
| Conda environments | Single shared env (or none) | Per-rule isolated envs under `workflow/envs/*.yml` |
| Conda prefix | In home directory | Must be set outside home (home quota = 20 GB; rule envs take 5–15 GB total) |
| Session persistence | Not documented | `screen` required on `bsub01` to survive SSH disconnects |

### Node roles at DKFZ

| Node | Purpose | Allowed |
|------|---------|---------|
| `odcf-worker01/02` | Dev, install, testing | ✅ Software install, small test runs |
| `bsub01/02` | Job submission only | ✅ Run Snakemake controller (lightweight); ❌ No data processing |
| Cluster nodes | Computation | Jobs submitted automatically via `bsub` |

### Quick setup (NEW pipeline)

**Step 1 — Configure conda and set up env on `odcf-worker01`:**

> Do this on `odcf-worker01`, not on `bsub01`. Worker nodes allow software installation.

```bash
ssh YOUR_USERNAME@odcf-worker01.dkfz.de

# Required: DKFZ bans the defaults (Anaconda) channel due to licensing
cat > ~/.condarc << 'EOF'
channels:
  - conda-forge
  - bioconda
EOF

# Load Mamba and initialise shell permanently
module load Mamba/24.11.2-1
mamba init bash
source ~/.bashrc

# Create controller env outside home (home quota = 20 GB)
YOUR_WORKDIR="/omics/groups/OE0146/internal/YOUR_USERNAME"
mkdir -p ${YOUR_WORKDIR}/conda_envs

mamba create -p ${YOUR_WORKDIR}/conda_envs/snakemake \
    -c conda-forge -c bioconda \
    snakemake snakemake-executor-plugin-lsf -y

mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

# Pin numpy/pandas to versions tested with the pipeline's helper scripts
python -m pip install "snakemake==8.*" "snakemake-executor-plugin-lsf" "numpy==1.26.4" "pandas==2.2.3"

# Verify
python -c "import snakemake, numpy, pandas; print(snakemake.__version__, numpy.__version__, pandas.__version__)"
```

**Step 2 — Clone and configure:**

```bash
cd ${YOUR_WORKDIR}
git clone https://github.com/UKHD-NP/atacseq_snakemake_new.git
cd atacseq_snakemake_new

# Edit config: set samples_csv, ref.assembly, enable/disable modules
# See README Key Configuration section for all options
nano config/config.yml

# Update conda-prefix in LSF profile to point outside home
sed -i "s|/omics/odcf/analysis/YOUR_GROUP/conda_envs|${YOUR_WORKDIR}/conda_envs|g" \
    workflow/profiles/lsf/config.yaml

# Confirm the replacement was applied
grep "conda-prefix" workflow/profiles/lsf/config.yaml
```

**Step 3 — Dry-run from `bsub01` (validate before submitting):**

> Switch to `bsub01` for all Snakemake operations that interact with LSF.

```bash
ssh YOUR_USERNAME@bsub01.lsf.dkfz.de

YOUR_WORKDIR="/omics/groups/OE0146/internal/YOUR_USERNAME"
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake
cd ${YOUR_WORKDIR}/atacseq_snakemake_new

# Dry-run: resolves full DAG, prints all rules, submits nothing
snakemake -s workflow/Snakefile --configfile config/config.yml --use-conda -n
```

**Step 4 — Run in a persistent screen session:**

```bash
# Start named screen session — survives SSH disconnect
screen -S atacseq

# Launch pipeline — Snakemake submits each rule as a separate bsub job
# -j 100: allow up to 100 concurrent cluster jobs
snakemake --profile workflow/profiles/lsf -j 100
```

| screen command | Action |
|---------------|--------|
| `screen -S atacseq` | Start new named session |
| `Ctrl+A`, then `D` | Detach (session keeps running) |
| `screen -ls` | List all active sessions |
| `screen -r atacseq` | Re-attach to session |

**Monitor LSF jobs:**

```bash
bjobs -w           # all running/pending jobs
bjobs -w -r        # running only
bjobs -w -p        # pending only
bjobs -l JOB_ID    # detailed info for one job
```

### Important notes for DKFZ

- `~/.condarc` must list only `conda-forge` and `bioconda` — `defaults` (Anaconda) is banned at DKFZ.
- `conda-prefix` in `workflow/profiles/lsf/config.yaml` must point **outside** home (`/omics/...`). All rule envs together take 5–15 GB.
- Always launch Snakemake from `bsub01/bsub02`, never from `odcf-worker01` — worker nodes cannot submit LSF jobs.
- `ref.mito_name` in `config.yml` must exactly match the mitochondrial contig name in your FASTA (e.g. `MT`, `chrM`, or `M`).
