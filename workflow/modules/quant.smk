FEATURE_COUNTS_CFG = config.get("feature_counts", {}) if isinstance(config.get("feature_counts", {}), dict) else {}
FEATURE_COUNTS_USE_SHIFTED_BAM = as_bool(FEATURE_COUNTS_CFG.get("use_shifted_bam", True), default=True)


def get_bam_for_peak_counting(wildcards):
    """Use filtered BAM for broad peaks (BAMPE); shifted or filtered BAM for narrow peaks (configurable)."""
    if CALL_PEAKS_PEAK_TYPE == "narrow" and not FEATURE_COUNTS_USE_SHIFTED_BAM:
        return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.filtered.bam")
    return get_bam_for_callpeaks(wildcards)


rule featurecounts_in_peaks:
    # Count paired-end reads overlapping called peaks (SAF) with featureCounts.
    input:
        peaks = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        bam = get_bam_for_peak_counting,
    output:
        saf = os.path.join("{outdir}", "featurecounts", "{sample_id}_peaks.saf"),
        counts = os.path.join("{outdir}", "featurecounts", "{sample_id}.readCountInPeaks.txt"),
        summary = os.path.join("{outdir}", "featurecounts", "{sample_id}.readCountInPeaks.txt.summary"),
    params:
        frac_overlap = float(CALL_PEAKS_CFG.get("frip_overlap_fraction", 0.20)),
        strand = 0, # unstranded
    conda:
        os.path.join(workflow.basedir, "envs", "subread.yml")
    message:
        "{wildcards.sample_id}: Running featureCounts on peaks"
    threads: 8
    log:
        os.path.join("{outdir}", "logs", "featurecounts", "{sample_id}.featureCounts.log"),
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.featurecounts_in_peaks.benchmark.txt")
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.counts}")"
        mkdir -p "$(dirname "{log}")"

        awk 'BEGIN{{FS=OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' \
            "{input.peaks}" > "{output.saf}" 2> "{log}" || {{
            echo "[ERROR] Peak-to-SAF conversion failed." >> "{log}"
            exit 1
        }}

        featureCounts \
            -F SAF \
            -O \
            --fracOverlap {params.frac_overlap} \
            -p \
            -T {threads} \
            -a "{output.saf}" \
            -s {params.strand} \
            -o "{output.counts}" \
            "{input.bam}" >> "{log}" 2>&1 || {{
            echo "[ERROR] featureCounts failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.counts}" ] || [ ! -s "{output.summary}" ]; then
            echo "[ERROR] featureCounts outputs missing or empty." >> "{log}"
            exit 1
        fi
        """
