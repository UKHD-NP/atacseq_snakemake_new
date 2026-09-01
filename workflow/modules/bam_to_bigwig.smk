rule bamcoverage_bigwig:
    # Build CPM-normalized bigWig from filtered BAM using deepTools bamCoverage.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.bai")
    output:
        bigwig = os.path.join("{outdir}", "bigwig", "{sample_id}.bigWig")
    params:
        bin_size = 10
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Building CPM-normalized bigWig with bamCoverage - [Source: Filtered BAM, NOT SHIFTED]"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 4096 + (attempt - 1) * 2048,
        runtime = lambda wildcards, attempt: attempt * 480
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.bamCoverage.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bamCoverage.benchmark.txt")
    shell:
        """
        ulimit -Sn $(ulimit -Hn) 2>/dev/null || ulimit -n 65536 2>/dev/null || true
        echo "[INFO] File descriptor limit: $(ulimit -n)" >> "{log}"

        mkdir -p "$(dirname "{output.bigwig}")"
        mkdir -p "$(dirname "{log}")"

        echo "[INFO] bamCoverage start: $(date)" >> "{log}"

        bamCoverage \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing CPM \
            --bam "{input.bam}" \
            -o "{output.bigwig}" \
            2>> "{log}" || {{
            echo "[ERROR] bamCoverage failed: $(date)" >> "{log}"
            exit 1
        }}

        echo "[INFO] bamCoverage done: $(date)" >> "{log}"

        if [ ! -s "{output.bigwig}" ]; then
            echo "[ERROR] bigWig output is missing or empty: $(date)" >> "{log}"
            exit 1
        fi
        """

