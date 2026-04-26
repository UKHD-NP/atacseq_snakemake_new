rule frip_score:
    # FRiP score reported two ways:
    #   1. Bedtools: reads-in-peaks / total mapped reads (from flagstat).
    #   2. featureCounts: "Successfully assigned alignments" % from featureCounts log.
    # fc_summary is input only to ensure featurecounts_in_peaks runs first;
    # the actual log is referenced via params.fc_log.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        peaks = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        flagstat = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.flagstat"),
        fc_summary = os.path.join("{outdir}", "featurecounts", "{sample_id}.readCountInPeaks.txt.summary")
    output:
        frip = os.path.join("{outdir}", "peaks", "{sample_id}.FRiP.txt"),
        peak_count_mqc = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.count_mqc.tsv"),
        frip_mqc = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.FRiP_mqc.tsv")
    params:
        overlap = float(CALL_PEAKS_CFG.get("frip_overlap_fraction", 0.20)),
        frip_threshold = int(CALL_PEAKS_CFG.get("frip_threshold", 20)),
        fc_log = os.path.join("{outdir}", "logs", "featurecounts", "{sample_id}.featureCounts.log"),
    conda:
        os.path.join(workflow.basedir, "envs", "bedtools.yml")
    message:
        "{wildcards.sample_id}: Calculating FRiP score"
    threads: 1
    resources:
        mem_mb = 2048
    log:
        os.path.join("{outdir}", "logs", "frip", "{sample_id}.frip.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.frip_score.benchmark.txt")
    shell:
        """
        set -euo pipefail
        export LC_NUMERIC=C
        mkdir -p "$(dirname "{output.frip}")"
        mkdir -p "$(dirname "{log}")"

        # 1. Bedtools FRiP
        READS_IN_PEAKS=$(bedtools intersect \
            -a "{input.bam}" \
            -b "{input.peaks}" \
            -bed \
            -c \
            -f {params.overlap} \
            2>> "{log}" | awk 'BEGIN{{sum=0}} $NF > 0 {{sum++}} END{{print sum}}') || {{
            echo "[ERROR] bedtools intersect failed." >> "{log}"
            exit 1
        }}

        TOTAL_READS=$(grep 'mapped (' "{input.flagstat}" \
            | grep -v "primary" \
            | awk '{{print $1}}')
        TOTAL_READS=${{TOTAL_READS:-0}}

        if [ "$TOTAL_READS" -gt 0 ]; then
            FRIP_BED=$(awk -v r="$READS_IN_PEAKS" -v t="$TOTAL_READS" \
                'BEGIN{{printf "%.2f", r/t*100}}')
        else
            FRIP_BED="0.00"
        fi

        # 2. featureCounts FRiP (parse log)
        FRIP_FC=$(grep 'Successfully assigned alignments' "{params.fc_log}" \
            | sed -e 's/.*(//' | sed 's/%.*$//' || echo "NA")

        # 3. Quality label
        FRIP_INT=$(printf "%.0f" "$FRIP_BED")
        LABEL=$(if [ "$FRIP_INT" -ge {params.frip_threshold} ]; then echo 'good'; else echo 'bad'; fi)

        # 4. Combined FRiP report
        printf 'Sample\\tFRiP_bedtools\\tFRiP_featureCounts\\tquality\\n' > "{output.frip}"
        printf '%s\\t%s%%\\t%s%%\\t%s\\n' \
            "{wildcards.sample_id}" "$FRIP_BED" "$FRIP_FC" "$LABEL" >> "{output.frip}"

        echo "[INFO] FRiP bedtools: ${{FRIP_BED}}% | featureCounts: ${{FRIP_FC}}% | label: $LABEL (threshold={params.frip_threshold}%)" >> "{log}"

        # 5. Peak count MultiQC TSV
        PEAK_COUNT=$(wc -l < "{input.peaks}" | awk '{{print $1}}')
        printf 'Sample\\tpeak_count\\n' > "{output.peak_count_mqc}"
        printf '%s_peaks\\t%s\\n' "{wildcards.sample_id}" "$PEAK_COUNT" >> "{output.peak_count_mqc}"

        # 6. FRiP MultiQC TSV (both methods)
        printf 'Sample\\tFRiP_bedtools\\tFRiP_featureCounts\\n' > "{output.frip_mqc}"
        printf '%s\\t%s\\t%s\\n' "{wildcards.sample_id}" "$FRIP_BED" "$FRIP_FC" >> "{output.frip_mqc}"

        if [ ! -s "{output.frip}" ] || [ ! -s "{output.peak_count_mqc}" ] || [ ! -s "{output.frip_mqc}" ]; then
            echo "[ERROR] FRiP outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """
