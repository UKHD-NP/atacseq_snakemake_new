FEATURE_COUNTS_CFG = config.get("feature_counts", {}) if isinstance(config.get("feature_counts", {}), dict) else {}
FEATURE_COUNTS_USE_SHIFTED_BAM = as_bool(FEATURE_COUNTS_CFG.get("use_shifted_bam", True), default=True)


def get_bam_for_peak_counting(wildcards):
    """Use shifted BAM for featureCounts if enabled, otherwise filtered BAM."""
    if FEATURE_COUNTS_USE_SHIFTED_BAM:
        return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.shifted.bam")
    return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.filtered.bam")


rule featurecounts_in_peaks:
    # Count paired-end reads overlapping called peaks (SAF) with featureCounts.
    params:
        frac_overlap = float(CALL_PEAKS_CFG.get("frip_overlap_fraction", 0.20)),
        strand = 0,
    input:
        peaks = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        bam = get_bam_for_peak_counting,
    output:
        saf = os.path.join("{outdir}", "featurecounts", "{sample_id}_peaks.saf"),
        counts = os.path.join("{outdir}", "featurecounts", "{sample_id}.readCountInPeaks.txt"),
        summary = os.path.join("{outdir}", "featurecounts", "{sample_id}.readCountInPeaks.txt.summary"),
    log:
        os.path.join("{outdir}", "logs", "featurecounts", "{sample_id}.featureCounts.log"),
    conda:
        os.path.join(workflow.basedir, "envs", "subread.yml")
    threads: 8
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.featurecounts_in_peaks.benchmark.txt")
    message:
        "{wildcards.sample_id}: Running featureCounts on peaks"
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
