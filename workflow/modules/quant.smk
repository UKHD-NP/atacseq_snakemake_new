rule featurecounts_in_peaks:
    # Count paired-end reads overlapping called peaks (SAF) with featureCounts.
    input:
        peaks = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
    output:
        saf = temp(os.path.join("{outdir}", "featurecounts", "{sample_id}_peaks.saf")),
        counts = os.path.join("{outdir}", "featurecounts", "{sample_id}.readCountInPeaks.txt"),
        summary = os.path.join("{outdir}", "featurecounts", "{sample_id}.readCountInPeaks.txt.summary"),
    params:
        strand = 0,  # unstranded
        frac_overlap = float(CALL_PEAKS_CFG.get("frip_overlap_fraction", 0.20)),
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
            -p \
            --fracOverlap {params.frac_overlap} \
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
