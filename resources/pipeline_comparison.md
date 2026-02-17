# Comparison: OLD vs NEW ATAC-seq Pipeline

## 1. Architecture

| | OLD | NEW |
|---|---|---|
| Structure | Single monolithic Snakefile | 20 modular `.smk` files |
| Configuration | Hardcoded in rules | YAML-driven, each module toggleable on/off |
| References | Manual setup | Auto-download hg19/hg38/m39 |

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
| **PBC calculation** (PCR Bottleneck Coefficient) | Completely absent |
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
| **Auto reference download** | hg19/hg38/m39 |
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
