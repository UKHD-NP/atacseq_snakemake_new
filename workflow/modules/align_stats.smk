rule samtools_stats:
    # Generate comprehensive statistics for BAM files
    input:
        bam  =  os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        fasta  =  config['ref']['fasta']
    output:
        bam_stats  =  os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.stats"),
        bam_flagstat  =  os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.flagstat"),
        bam_idxstats  =  os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.idxstats")
    message:
        "{wildcards.sample_id}: Running Samtools statistics"
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.samtools_stats.log")
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    threads: 4
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.samtools_stats.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.bam_stats})
        mkdir -p $(dirname {log})

        # Generate comprehensive BAM statistics
        samtools stats --threads {threads} -r {input.fasta} {input.bam} > {output.bam_stats} 2>> {log} || {{ echo "[ERROR] samtools stats failed." >> {log}; exit 1; }}

        # Create flagstat summary
        samtools flagstat --threads 2 {input.bam} > {output.bam_flagstat} 2>> {log} || {{ echo "[ERROR] samtools flagstat failed." >> {log}; exit 1; }}

        # Generate chromosome-level read mapping statistics
        samtools idxstats {input.bam} > {output.bam_idxstats} 2>> {log} || {{ echo "[ERROR] samtools idxstats failed." >> {log}; exit 1; }}
        """


rule picard_collect_multiple_metrics:
    # Collect Picard metrics on filtered BAM.
    params:
        output_prefix = lambda wildcards: os.path.join(
            wildcards.outdir,
            "bam",
            f"{wildcards.sample_id}.CollectMultipleMetrics"
        ),
        tmp_dir = lambda wildcards: os.path.join(wildcards.outdir, "bam", "tmp"),
        xmx = "4915M"
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        fasta = config["ref"]["fasta"]
    output:
        alignment_summary = os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.alignment_summary_metrics"
        )
    message:
        "{wildcards.sample_id}: Running Picard CollectMultipleMetrics"
    log:
        os.path.join("{outdir}", "logs", "picard", "{sample_id}.collect_multiple_metrics.log")
    conda:
        os.path.join(workflow.basedir, "envs", "picard_markduplicates.yml")
    threads: 1
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.collect_multiple_metrics.benchmark.txt")
    shell:
        """
        mkdir -p "{params.tmp_dir}"
        mkdir -p "$(dirname "{log}")"

        picard \
            -Xmx{params.xmx} \
            CollectMultipleMetrics \
            --VALIDATION_STRINGENCY LENIENT \
            --TMP_DIR "{params.tmp_dir}" \
            --INPUT "{input.bam}" \
            --OUTPUT "{params.output_prefix}" \
            --REFERENCE_SEQUENCE "{input.fasta}" \
            > "{log}" 2>&1 || {{ echo "[ERROR] Picard CollectMultipleMetrics failed." >> "{log}"; exit 1; }}

        if [ ! -s "{output.alignment_summary}" ]; then
            echo "[ERROR] Missing Picard output: {output.alignment_summary}" >> "{log}"
            exit 1
        fi

        rm -rf "{params.tmp_dir}"
        """
