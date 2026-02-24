# Comparison: OLD vs NEW ATAC-seq Pipeline

## 1. Architecture

| | OLD | NEW |
|---|---|---|
| Structure | Single monolithic Snakefile | 20 modular `.smk` files |
| Configuration | Hardcoded in rules | YAML-driven, each module toggleable on/off |
| References | Manual setup | Auto-download hg19/hg38 |

---

## 2. Step-by-step comparison

| Step | OLD | NEW | Key difference |
|------|-----|-----|----------------|
| **Trimming** | cutadapt | fastp (default) or trim_galore | fastp auto-detects adapters, produces HTML/JSON reports |
| **Alignment** | bwa-mem2 only | bwa-mem2 (default) or bowtie2 | NEW adds RG tags, supports 2 aligners |
| **Filtering** | samtools `-e '([NM]<=4) && sclen<15'` + `-F 0x100 -f 3` | bamtools JSON filter + `-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400` | NEW adds orphan removal, uses include-regions BED instead of grep chrM |
| **Mark dup** | GATK `MarkDuplicatesWithMateCigar` (removes dups) | Picard `MarkDuplicates` (marks only, removed at filter step) | Different order: OLD filters then marks; NEW marks then filters (more correct) |
| **Tn5 shift** | 1 massive rule: BAM->BED->shift->BEDPE->shift again->blacklist->RPM->bigwig | 2 separate mechanisms: (1) awk shift BED for peak calling, (2) deepTools `alignmentSieve --ATACshift` for BAM | NEW cleanly separates concerns |
| **BigWig** | 1 RPM track from shifted data | 2 tracks: RPGC (shifted) + RPM (unshifted) | RPGC is more appropriate for cross-sample comparison |
| **Peak calling** | MACS3 on shifted BED | MACS3, supports narrow + broad | NEW adds peak QC plots |
| **Normalization** | RPM only | RPGC (shifted) + RPM (unshifted) | RPGC is the better method |

---

## 3. Features in OLD but MISSING in NEW

| Feature | Notes |
|---------|-------|
| **Fragment size distribution** (`bamPEFragmentSize`) | No equivalent rule in NEW |
| **PBC calculation** (PCR Bottleneck Coefficient) | Completely absent — recommended insertion point: after `mark_duplicates`, before duplicate removal (see Section 8) |
| **Merged Lorenz curve plots** | NEW only has per-sample, no merge |
| **FastQC on BAMs** | Not present in NEW |
| **Soft-clip filter** (`sclen < 15`) | Not in NEW |
| **BEDPE-based shifting** | NEW only shifts single-end BED |

---

## 4. Features in NEW but MISSING in OLD

| Feature | Notes |
|---------|-------|
| **Multi-lane FASTQ merging** | Supports multiple sequencing lanes |
| **FastQC on raw reads** | Pre-trimming QC |
| **Bowtie2 support** | Alternative aligner |
| **HOMER peak annotation** | Annotates peaks to genomic features |
| **featureCounts quantification** | Read counting in peaks |
| **ataqv** | ATAC-seq specific QC + interactive HTML report |
| **deepTools heatmap/plotProfile** | TSS-centered heatmap, gene-body profiles |
| **Auto reference download** | hg19/hg38 |
| **Auto include-regions** | Genome minus blacklist minus mito |
| **Disk cleanup** | Auto-deletes intermediate files |

---

## 5. Pipeline order

**OLD:** Trim -> Align -> **Filter -> Mark dup** -> Shift+BigWig -> Peak call -> QC

**NEW:** Trim -> Align -> Sort -> **Mark dup -> Filter** -> Shift BAM -> BigWig -> Peak call -> QC -> MultiQC -> Cleanup

> **Important:** OLD filters before marking duplicates. NEW marks duplicates first, then filters.
> The NEW order is more correct because duplicate detection needs to see all reads to accurately identify duplicates.

---

## 6. Detailed parameter comparison

| Parameter | OLD | NEW |
|-----------|-----|-----|
| MAPQ threshold | Configurable `-q MAPQ` | `-q 30` in config (code fallback `-q 20`), configurable via `bam_filter.params` |
| SAM flag filters | `-F 0x100 -f 3 -F 0x0008` | `-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400 -q 30` |
| NM mismatch filter | `[NM] <= 4` (samtools expression) | Via bamtools JSON script |
| Soft-clip filter | `sclen < 15` (samtools expression) | Not present |
| Tn5 shift values | +4 (forward) / -5 (reverse) | +4/-5 for BED; `--ATACshift` for BAM (deepTools) |
| MACS3 params | Not specified in code | `--shift 75 --extsize 150 --keep-dup all --nomodel --call-summits -q 0.01` (narrow) |
| Normalization | RPM | RPGC (shifted) + RPM (unshifted) |
| Duplicate handling | GATK MarkDuplicatesWithMateCigar (removes) | Picard MarkDuplicates (marks only, `-F 0x0400` removes downstream) |
| Effective genome size | Manually configured | Auto-computed from chromsizes or user-configured |

---

## 7. BAM filtering comparison

### Filter-by-filter breakdown

| Filter | OLD (`samtools -e`) | NEW (`samtools flags` + `bamtools JSON`) | Winner |
|--------|---------------------|------------------------------------------|--------|
| Unmapped reads | Not explicitly filtered | `-F 0x004` (remove unmapped) | NEW |
| Mate unmapped | `-F 0x0008` | `-F 0x0008` | Same |
| Proper pair | `-f 3` (mapped + proper pair) | `-f 0x001` (paired only) + `bampe_rm_orphan.py --only_fr_pairs` (same-chr, FR orientation, correct positioning) | NEW (stricter in practice) |
| Secondary alignments | `-F 0x100` | `-F 0x0100` | Same |
| Duplicates | Not filtered (handled separately) | `-F 0x0400` (remove dup-flagged reads) | NEW |
| MAPQ | `-q MAPQ` (configurable) | `-q 30` (config default, code fallback `-q 20`) | NEW (stricter) |
| NM mismatch | `[NM] <= 4` (samtools expression) | `NM:<=4` via bamtools JSON | Same (threshold identical) |
| Soft-clip | `sclen < 15` (allows small soft-clips) | `!cigar *S*` (rejects ANY soft-clip) | OLD (NEW is too aggressive) |
| Insert size | Not filtered | `insertSize >= -2000 && <= 2000` via bamtools | NEW |
| MT / blacklist | `grep -v chrM` | `-L include_regions` BED | NEW (more flexible) |
| Orphan reads | Not handled | `bampe_rm_orphan.py --only_fr_pairs` | NEW |

### How NEW pipeline handles proper-pair filtering

The OLD pipeline uses `-f 3` (samtools proper-pair flag) which trusts the aligner's definition.

The NEW pipeline takes a different, more thorough approach via `bampe_rm_orphan.py --only_fr_pairs`:
1. Removes singleton/orphan reads (no mate in BAM)
2. Keeps only pairs on the **same chromosome** (`pair1.tid == pair2.tid`)
3. Keeps only **FR orientation** (read1 forward + read2 reverse, or vice versa, with correct positioning)
4. Removes all improper pairs (wrong orientation, interchromosomal, etc.)

Combined with bamtools insert-size filter (`-2000 to 2000`), this is actually **stricter** than
the aligner's `-f 0x002` flag because it verifies orientation and positioning at the read level.

### Remaining issue in NEW pipeline

1. **Soft-clip filter too aggressive**: `!cigar *S*` in bamtools JSON rejects any read with even
   1bp soft-clip. Most aligners produce some soft-clipping at read ends. This can discard 10-30%
   of good reads. The OLD approach (`sclen < 15`) is more reasonable -- allows small soft-clips
   while removing reads with excessive soft-clipping.

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

Two valid approaches for soft-clip:
- **Option A**: Drop it entirely -- NM + MAPQ filters already catch most problematic alignments.
- **Option B**: Add `samtools view -e 'sclen < 15'` as an additional pipe step to match OLD behavior.

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
Trim -> Align -> Sort -> Mark dup -> [PBC HERE] -> Filter -> ...
```

**Input:** marked BAM (output of Picard MarkDuplicates, duplicates flagged but NOT yet removed)
**Output:** `{sample}.pbc.txt`

Rationale:
- Computing BEFORE duplicate removal: still have all reads including duplicates to count stacking at each position.
- Computing AFTER duplicate removal (`-F 0x0400`): information is lost — M1/M2 counts will be wrong.
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
| Conda prefix | In home directory | Must be set outside home (home quota = 20 GB) |
| Session persistence | Not documented | `screen` recommended on `bsub01` |

### Node roles at DKFZ

| Node | Purpose | Allowed |
|------|---------|---------|
| `odcf-worker01/02` | Dev, install, testing | Software install, small runs |
| `bsub01/02` | Job submission only | Run Snakemake (lightweight); **no** processing |
| Cluster nodes | Computation | Jobs submitted automatically via `bsub` |

### Quick setup (NEW pipeline)

**Step 1 — Set up env on `odcf-worker01`:**

```bash
module load Mamba/24.11.2-1
mamba init bash && source ~/.bashrc

YOUR_WORKDIR="/omics/groups/OE0146/internal/YOUR_USERNAME"
mkdir -p ${YOUR_WORKDIR}/conda_envs

mamba create -p ${YOUR_WORKDIR}/conda_envs/snakemake \
    -c conda-forge -c bioconda \
    snakemake snakemake-executor-plugin-lsf -y
```

**Step 2 — Clone and configure:**

```bash
cd ${YOUR_WORKDIR}
git clone https://github.com/UKHD-NPS/atacseq_snakemake.git
cd atacseq_snakemake

# Update conda-prefix in LSF profile to your workdir
sed -i "s|/omics/odcf/analysis/YOUR_GROUP/conda_envs|${YOUR_WORKDIR}/conda_envs|g" \
    workflow/profiles/lsf/config.yaml
```

**Step 3 — Dry-run (validate from bsub01):**

```bash
ssh YOUR_USERNAME@bsub01.lsf.dkfz.de
module load Mamba/24.11.2-1
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

snakemake -s workflow/Snakefile --configfile config/config.yml --use-conda -n
```

**Step 4 — Run in a persistent screen session:**

```bash
screen -S atacseq
snakemake --profile workflow/profiles/lsf -j 100
```

| screen command | Action |
|---------------|--------|
| `screen -S atacseq` | Start new named session |
| `Ctrl+A`, then `D` | Detach (keeps running) |
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

- `~/.condarc` must list only `conda-forge` and `bioconda` — `defaults` (Anaconda) channel is banned at DKFZ due to licensing.
- `conda-prefix` in `workflow/profiles/lsf/config.yaml` must point **outside** home (`/omics/...`). All rule envs together take 5–15 GB.
- Always launch Snakemake from `bsub01/bsub02`, never from `odcf-worker01` (worker nodes cannot submit LSF jobs).
