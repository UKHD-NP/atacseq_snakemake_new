rule frip_score:
    # FRiP = reads in peaks / mapped reads.
    # Reuses the flagstat produced by samtools_stats to avoid redundant computation.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        peaks = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        flagstat = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.flagstat")
    output:
        frip = os.path.join("{outdir}", "peaks", "{sample_id}.FRiP.txt"),
        peak_count_mqc = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.count_mqc.tsv"),
        frip_mqc = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.FRiP_mqc.tsv")
    params:
        overlap = float(CALL_PEAKS_CFG.get("frip_overlap_fraction", 0.20))
    conda:
        os.path.join(workflow.basedir, "envs", "bedtools.yml")
    message:
        "{wildcards.sample_id}: Calculating FRiP score"
    threads: 8
    log:
        os.path.join("{outdir}", "logs", "frip", "{sample_id}.frip.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.frip_score.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.frip}")"
        mkdir -p "$(dirname "{log}")"

        READS_IN_PEAKS=$(bedtools intersect \
            -a "{input.bam}" \
            -b "{input.peaks}" \
            -bed \
            -c \
            -f {params.overlap} \
            | awk 'BEGIN{{sum=0}} $NF > 0 {{sum++}} END{{print sum}}') || {{
            echo "[ERROR] Failed to compute reads in peaks." >> "{log}"
            exit 1
        }}

        grep 'mapped (' "{input.flagstat}" \
            | grep -v "primary" \
            | awk -v a="$READS_IN_PEAKS" -v sample="{wildcards.sample_id}" 'BEGIN{{OFS="\\t"}}{{if ($1 > 0) {{print sample, a/$1}} else {{print sample, 0}}}}' \
            > "{output.frip}" 2>> "{log}" || {{
            echo "[ERROR] Failed to write FRiP output." >> "{log}"
            exit 1
        }}

        PEAK_COUNT=$(wc -l < "{input.peaks}" | awk '{{print $1}}')
        {{
            printf "sample\\tpeak_count\\n"
            printf "%s_peaks\\t%s\\n" "{wildcards.sample_id}" "$PEAK_COUNT"
        }} > "{output.peak_count_mqc}" 2>> "{log}" || {{
            echo "[ERROR] Failed to write peak count MultiQC TSV." >> "{log}"
            exit 1
        }}

        {{
            printf "sample\\tfrip\\n"
            cat "{output.frip}"
        }} > "{output.frip_mqc}" 2>> "{log}" || {{
            echo "[ERROR] Failed to write FRiP MultiQC TSV." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.frip}" ] || [ ! -s "{output.peak_count_mqc}" ] || [ ! -s "{output.frip_mqc}" ]; then
            echo "[ERROR] FRiP outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """
