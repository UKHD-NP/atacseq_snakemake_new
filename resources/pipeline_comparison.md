# Comparison: OLD vs NEW ATAC-seq Pipeline

## 1. Architecture

| | OLD | NEW |
|---|---|---|
| Structure | Single monolithic Snakefile | 20 modular `.smk` files |
| Configuration | Hardcoded in rules | YAML-driven, each module toggleable on/off |
| References | Manual setup | Auto-staged for hg19/hg38; fully `custom` mode supported |
| Samplesheet columns | `sample`, `replicate`, `fq1`, `fq2`; no lane merging | `sample_id`, `fq1`, `fq2`, `outdir` → supports multi-lane FASTQ merging |

---

## 2. Step-by-step comparison

| Step | OLD | NEW | Key difference |
|------|-----|-----|----------------|
| **Trimming** | cutadapt with FIXED Nextera adapter `CTGTCTCTTATACACATCT` | fastp or trim_galore | fastp uses automatic adapter detection (`--detect_adapter_for_pe`); OLD used fixed hardcoded Nextera adapter |
| **FastQC raw** | **Not present** | `fastqc_raw` module | Completely absent in OLD |
| **FastQC trimmed** | Present | Present (trim_galore only; fastp has built-in HTML/JSON report) | trim_galore produces FastQC ZIPs; fastp produces JSON/HTML |
| **Alignment** | bwa-mem2 only | **bowtie2 (default)** or bwa-mem2 | bowtie2 is now default (`--very-sensitive --no-discordant -X 2000`); both aligners add RG tags |
| **Pre-filter stats** | **Not present** | `samtools_stats_pre_filter`: stats/flagstat/idxstats on markdup BAM | NEW — goes into MultiQC |
| **Filtering** | samtools `-q 20 -F 0x100 -e '([NM] <= 4) && sclen < 15' -f 3 -F 0x0008` then `grep -v chrM` via idxstats | samtools `-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400 -q 30` + bamtools JSON filter + pysam | NEW adds orphan removal, uses include-regions BED instead of `grep chrM` |
| **Post-filter stats** | **Not present** | `samtools_stats`: stats/flagstat/idxstats on filtered BAM | NEW — goes into MultiQC |
| **Picard CollectMultipleMetrics** | **Not present** | Alignment summary, base distribution, insert size, quality by cycle, quality distribution on filtered BAM | NEW — insert size and alignment metrics go into MultiQC |
| **GATK FixMateInformation** | After initial filter, before MarkDup: `gatk FixMateInformation --ADD_MATE_CIGAR true` | **Not present** | OLD required this to ensure mate CIGAR tags are present for MarkDuplicatesWithMateCigar |
| **Mark dup** | `gatk MarkDuplicatesWithMateCigar` (removes dups via `--REMOVE_DUPLICATES true`) | `picard MarkDuplicates` (marks only; removed at filter step via `samtools view -F 0x0400`) | Different tools; OLD used mate-CIGAR-aware marking |
| **Tn5 shift** | filtered BAM → `bedtools bamtobed` BED → awk shift (+4/-5) | **Peaks:** filtered BAM → `bedtools bamtobed` → awk shift (+4/-5) → BED (same approach as OLD); **Signal tracks:** `alignmentSieve --ATACshift` → `shifted.bam` → RPGC bigWig | Peak calling uses inline awk shift; deepTools signal tracks use alignmentSieve shifted BAM |
| **BigWig (unshifted)** | 1 RPM track from shifted data | RPM-scaled bedGraph→bigWig (`{sample_id}.bigWig`) from filtered BAM | Scale factor = 1,000,000 / mapped reads; saved to `{sample_id}.scale_factor.txt` |
| **BigWig (shifted)** | Not present separately | RPGC `{sample_id}.shifted.bigWig` from shifted BAM via bamCoverage | RPGC is more appropriate for cross-sample comparison |
| **Fragment size classes** | Fragment size distribution plot only | deepTools `bamPEFragmentSize` distribution + NFR/mono/di/tri **read counts** (`fragment_counts_mqc.tsv`) | NEW adds per-class counts for MultiQC |
| **Peak calling** | MACS3 on shifted BED (narrow only) | MACS3 **narrow** on shifted BED (same approach as OLD) or **broad** on unshifted BAM with BAMPE mode | NEW adds broad mode, peak QC plots |
| **TSS enrichment** | Custom `tss3.py` (Greenleaf lab script, `--greenleaf_norm --bins 400 --bp_edge 2000`) | deepTools `computeMatrix` / `plotProfile` / `plotHeatmap` (TSS-centered) + ATACseqQC `TSSEscore` + ataqv TSSE | Different approaches; OLD used a custom Python script |
| **Lorenz curve** | Present | Present (`plotFingerprint`) | Both use the same tool |
| **PBC** | `calculate_pbc` → `merge_pbc_results` (bedtools bedpe method, M0/M1/M2 awk) | **Not present** | Completely absent in NEW — see Section 8 for recommended insertion |
| **Peak annotation** | **Not present** | HOMER `annotatePeaks` + summary plot | NEW |
| **featureCounts** | count reads in peaks from unshifted BAM | count reads in peaks from unshifted BAM (`--fracOverlap 0.2`) | NEW with extra parameters |
| **ataqv** | **Not present** | ataqv JSON + mkarv interactive HTML | NEW |
| **MultiQC** | Separate per-module: `multiQC_trimmed_fastq` + `multiQC_BAMs` (not unified) | Single unified MultiQC report aggregating all QC modules | NEW provides a single aggregated report |
| **FastQC on BAMs** | `fastQC_BAMs` rule (FastQC on final filtered BAM) | **Not present** | OLD ran FastQC on the final BAM; NEW uses samtools stats + Picard instead |
| **Cleanup** | **Not present** | `delete_tmp` removes merged raw FASTQs, trimmed FASTQs (configurable via `trimming.delete_trimming`) | NEW |

---

## 3. Features in OLD but MISSING in NEW

| Feature | Notes |
|---------|-------|
| **PBC calculation** (PCR Bottleneck Coefficient) | Completely absent — recommended insertion point: after `mark_duplicates`, before `bam_filter` (see Section 8). OLD used bedtools bedpe method (`M0/M1/M2` awk pipeline). |
| **FastQC on BAMs** | OLD ran FastQC on the final filtered BAM (`fastQC_BAMs` rule); NOT present in NEW |
| **Soft-clip filter** (`sclen < 15`) | Not in NEW — bamtools JSON uses `!cigar *S*` (any soft-clip rejected, too aggressive; see Section 7) |
| **GATK FixMateInformation** | OLD ran `gatk FixMateInformation --ADD_MATE_CIGAR true` between initial filter and MarkDuplicates to ensure mate CIGAR tags are present; entirely absent in NEW |
| **TSS enrichment** (custom) | OLD used custom Greenleaf-lab `tss3.py` script with `--greenleaf_norm`; NEW uses deepTools `computeMatrix`/`plotProfile` instead |
| **Per-module MultiQC** (split reports) | OLD ran separate MultiQC per data type (`multiQC_trimmed_fastq`, `multiQC_BAMs`); NEW merges everything into one unified report |

---

## 4. Features in NEW but MISSING in OLD

| Feature | Notes |
|---------|-------|
| **Multi-lane FASTQ merging** | Multiple sequencing lanes per `sample_id` supported via samplesheet (`lane` column); OLD used `sample` + `replicate` with no lane merging |
| **FastQC on raw reads** | Pre-trimming QC (`fastqc_raw`) |
| **Bowtie2 support** | Alternative aligner, now the default |
| **Pre-filter BAM stats** | `samtools_stats_pre_filter`: stats/flagstat/idxstats on post-markdup BAM before filtering; goes to MultiQC |
| **Post-filter BAM stats** | `samtools_stats`: stats/flagstat/idxstats on filtered BAM; goes to MultiQC |
| **Picard CollectMultipleMetrics** | Alignment summary, insert size, base distribution by cycle, quality by cycle, quality distribution on filtered BAM |
| **NFR/mono bigWigs** | Fragment-size-filtered RPGC bigWig tracks: NFR (≤150 bp) and mononucleosomal (151–300 bp) from shifted BAM |
| **NFR vs mono TSS profile** | deepTools computeMatrix + plotProfile + plotHeatmap comparing NFR vs mono signal around TSS |
| **Fragment size class counts** | NFR/mono/di/tri read counts per sample → `fragment_counts_mqc.tsv` for MultiQC |
| **HOMER peak annotation** | Annotates peaks to genomic features with summary plot |
| **ataqv** | ATAC-seq-specific QC + interactive mkarv HTML report |
| **deepTools heatmap / plotProfile** | TSS-centered heatmap, gene-body profiles (shifted.bigWig) |
| **Auto include-regions** | Genome minus blacklist, optionally minus mito (`ref.keep_mito`) |
| **Mito name configurable** | `ref.mito_name` — must match FASTA contig exactly (MT / chrM / M) |
| **Broad peak mode** | `call_peaks.peak_type: broad` uses filtered.bam in BAMPE mode |
| **Disk cleanup** | Auto-deletes raw merged FASTQs; optionally trimmed FASTQs (`trimming.delete_trimming: true/false`) |
| **Dual FRiP reporting** | FRiP reported two ways: (1) `bedtools intersect` (reads-in-peaks / total mapped) and (2) featureCounts log ("Successfully assigned %"); both in `*.FRiP.txt` and `*_peaks.FRiP_mqc.tsv` for MultiQC |

---

## 5. Pipeline order

**OLD:** Trim → Align (`bwa-mem2` + `fixmate`) → Filter (samtools + grep chrM) → **FixMateInformation** → **MarkDuplicatesWithMateCigar** → Shift BED (`bamtobed` → awk Tn5 shift) + BigWig (RPM) → Peak call → FRiP (featureCounts) → deepTools → MultiQC (per-module)

**NEW (current):** Trim → Align (`bowtie2`/`bwa-mem2`) → Sort → **MarkDuplicates** → **Pre-filter stats** → Filter (samtools + bamtools + pysam) → **Post-filter stats** → **Picard CollectMultipleMetrics** → Signal tracks: `alignmentSieve` shifted.bam + RPGC shifted.bigWig; Shift BED (`bamtobed` → awk Tn5 shift); RPM unshifted.bigWig → **NFR/mono bigWigs** → **NFR vs mono computeMatrix/plotProfile/plotHeatmap** → Peak call → Annotation → FRiP (bedtools + featureCounts) → **Fragment size class counts** → deepTools (plotFingerprint + bamPEFragmentSize + gene body + TSS) → ataqv → ATACseqQC → MultiQC (unified) → Cleanup

---

## 6. Detailed parameter comparison

| Parameter | OLD | NEW |
|-----------|-----|-----|
| Default trimmer | cutadapt with FIXED adapter `CTGTCTCTTATACACATCT` (hardcoded Nextera); `--nextseq-trim=20 --minimum-length 20` | **trim_galore** with automatic adapter detection (`--nextseq 25 --length 36`) |
| Default aligner | bwa-mem2 (piped with `samtools fixmate -m` during alignment) | **bowtie2** (`--very-sensitive --no-discordant -p 2 -X 2000`) |
| MAPQ threshold | `-q 20` (configurable via `MAPQ_threshold`) | `-q 30` (configurable via `bam_filter.params`) |
| SAM flag filters | `-q 20 -F 0x100 -e '([NM] <= 4) && sclen < 15' -f 3 -F 0x0008` | `-q 30 -F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400` |
| NM mismatch filter | `[NM] <= 4` (samtools `-e` expression, inline with flag filters) | `NM:<=4` via bamtools JSON |
| Soft-clip filter | `sclen < 15` (inline samtools expression; allows small soft-clips) | `!cigar *S*` (rejects ANY soft-clip — too aggressive; see Section 7) |
| Insert size filter | Not filtered at samtools step | `insertSize >= -2000 && <= 2000` via bamtools JSON |
| Fragment length filter | `minFragmentLength 0, maxFragmentLength 2000` applied at Tn5 shift step | Not explicitly filtered |
| MT / chrM removal | `grep -v chrM` | `-L include_regions` BED (controlled by `ref.keep_mito` + `ref.mito_name`) |
| Mate CIGAR repair | `gatk FixMateInformation --ADD_MATE_CIGAR true` (between filter and MarkDup) | Not present |
| Tn5 shift values | +4 (forward) / -5 (reverse) via bedtools + awk | **Peaks:** +4 (forward) / -5 (reverse) via bedtools + awk (same as OLD); **Signal tracks:** `alignmentSieve --ATACshift` |
| MACS3 narrow params | `--shift -75 --extsize 150 --keep-dup all --nomodel --call-summits --nolambda -q 0.05` | `--shift -75 --extsize 150 --keep-dup all --nomodel --call-summits -q 0.01` |
| MACS3 broad params | Not present | `--keep-dup all --nomodel --broad --broad-cutoff 0.1` |
| Effective genome size | Configurable via `effective_genomeSize` | Configurable via `macs3_gsize`; auto-sum from chromsizes if empty |
| FRiP overlap threshold | Not present | `frip_overlap_fraction: 0.2` (bedtools + featureCounts) |
| featureCounts BAM source | Unshifted filtered BAM | Unshifted filtered BAM; `--fracOverlap 0.2` |
| BigWig (unshifted) | RPM (1 `shifted.bigWig` per sample) | RPM-scaled from filtered BAM → `{sample_id}.bigWig`; scale factor saved to `{sample_id}.scale_factor.txt` |
| BigWig (shifted) | Not present separately | RPGC from shifted BAM → `{sample_id}.shifted.bigWig` (bin size 10 bp) |
| NFR/mono bigWigs | Not present | RPGC from alignmentSieve-filtered shifted BAM → `{sample_id}.nfr.bigWig`, `{sample_id}.mono.bigWig` (bin size 10 bp) |
| Fragment size classes | Not present | NFR ≤150 bp, mono 151–300 bp, di 301–500 bp, tri 501–700 bp; configurable via `nfr.*` config keys |
| Duplicate handling | `gatk MarkDuplicatesWithMateCigar` (removes if `--REMOVE_DUPLICATES true`) | `picard MarkDuplicates` (marks only; `samtools view -F 0x0400` removes downstream at filter step) |

---

## 7. BAM filtering comparison

### Filter-by-filter breakdown

| Filter | OLD (`samtools -e` + `grep -v chrM`) | NEW (`samtools flags` + `bamtools JSON` + `pysam`) | Winner |
|--------|---------------------|------------------------------------------|--------|
| Unmapped reads | Not explicitly filtered (relies on proper-pair flag) | `-F 0x004` (remove unmapped) | NEW |
| Mate unmapped | `-F 0x0008` | `-F 0x0008` | Same |
| Proper pair | `-f 3` (mapped + proper pair; trusts aligner) | `-f 0x001` (paired only) + `bampe_rm_orphan.py --only_fr_pairs` (same-chr, FR orientation) | NEW (stricter in practice) |
| Secondary alignments | `-F 0x100` | `-F 0x0100` | Same |
| Duplicates | `gatk MarkDuplicatesWithMateCigar --REMOVE_DUPLICATES true` (removes immediately) | `picard MarkDuplicates` marks; `samtools view -F 0x0400` removes at filter step | Same outcome |
| MAPQ | `-q 20` (configurable via `MAPQ_threshold`) | `-q 30` (configurable via `bam_filter.params`) | NEW (stricter default) |
| NM mismatch | `[NM] <= 4` (inline samtools `-e` expression) | `NM:<=4` via bamtools JSON | Same (threshold identical) |
| Soft-clip | `sclen < 15` (inline; allows small soft-clips) | `!cigar *S*` (rejects ANY soft-clip) | OLD (NEW is too aggressive) |
| Insert size | Not filtered at samtools step | `insertSize >= -2000 && <= 2000` via bamtools | NEW |
| Fragment length | `minFragmentLength 0, maxFragmentLength 2000` applied at Tn5 shift step (post-filter) | Not explicitly filtered | OLD |
| MT / chrM removal | `grep -v chrM` (hardcoded contig name) | `-L include_regions` BED (controlled by `ref.keep_mito` + `ref.mito_name`) | NEW (more flexible, contig-name agnostic) |
| Mate CIGAR repair | `gatk FixMateInformation --ADD_MATE_CIGAR true` (after filter, before MarkDup) | Not present | OLD (required for MarkDuplicatesWithMateCigar) |
| Orphan reads | Not explicitly handled (`-f 3` relies on aligner flags) | `bampe_rm_orphan.py --only_fr_pairs` (removes singletons + wrong-orientation pairs) | NEW |

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

## 9. TSSE: ATACseqQC vs ataqv — why they differ

| | ATACseqQC | ataqv |
|---|---|---|
| **Input BAM** | `shifted.bam` (+4/−5 bp Tn5 shift) | `filtered.bam` (non-shifted) |
| **Read counting** | Full read coverage | 5' end of reads only (Tn5 cut site proxy) |
| **Window** | ±1000 bp around TSS, 20 bins × 100 bp | ±1500 bp around TSS |
| **Smoothing** | LOESS, then take max | None |
| **Normalization** | Max enrichment / mean flanking signal | Flanking region average |

The Tn5 shift (+4/−5 bp) accounts for <0.1% of the window size and is **not** the cause of the difference. The difference is algorithmic. ataqv TSSE is closer to the original ENCODE definition (5' ends only); ATACseqQC TSSE uses full read coverage with LOESS smoothing and will systematically read 1–2 points higher.

**Both inputs are correct:** ataqv is designed for non-shifted filtered BAM; ATACseqQC requires shifted BAM. Do not pass shifted BAM to ataqv.

---

## 10. Running on HPC (DKFZ LSF cluster)

### OLD vs NEW — HPC support comparison

| | OLD | NEW |
|---|---|---|
| Cluster scheduler | Manual `bsub` scripts or none | `snakemake-executor-plugin-lsf` (automatic) |
| LSF profile | Not provided | `workflow/profiles/lsf/config.yaml` |
| Resource declaration | Not per-rule | Per-rule: `mem_mb`, `runtime`, `threads` → auto-translated to `bsub` flags |
| Conda environments | Single shared env (or none) | Per-rule isolated envs under `workflow/envs/*.yml` |
| Conda prefix | In home directory | Must be set outside home (home quota = 20 GB; rule envs take 5–15 GB total) |
| Session persistence | Not documented | `screen` required on `bsub01` to survive SSH disconnects |
