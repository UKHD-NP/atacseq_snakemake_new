rule samtools_stats_pre_filter:
    # Generate comprehensive statistics for BAM files before filtering
    input:
        bam   = get_bam_for_filter,
        fasta = config['ref']['fasta']
    output:
        bam_stats    = os.path.join("{outdir}", "bam", "{sample_id}.pre_filter.bam.stats"),
        bam_flagstat = os.path.join("{outdir}", "bam", "{sample_id}.pre_filter.bam.flagstat"),
        bam_idxstats = os.path.join("{outdir}", "bam", "{sample_id}.pre_filter.bam.idxstats")
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "{wildcards.sample_id}: Running Samtools statistics (pre-filter)"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.samtools_stats_pre_filter.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.samtools_stats_pre_filter.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.bam_stats})
        mkdir -p $(dirname {log})

        samtools stats --threads {threads} --reference {input.fasta} {input.bam} > {output.bam_stats} 2>> {log} || {{ echo "[ERROR] samtools stats failed." >> {log}; exit 1; }}

        samtools flagstat --threads {threads} {input.bam} > {output.bam_flagstat} 2>> {log} || {{ echo "[ERROR] samtools flagstat failed." >> {log}; exit 1; }}

        samtools idxstats --threads {threads} {input.bam} > {output.bam_idxstats} 2>> {log} || {{ echo "[ERROR] samtools idxstats failed." >> {log}; exit 1; }}
        """

rule samtools_stats:
    # Generate comprehensive statistics for BAM files
    input:
        bam  =  os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        fasta  =  config['ref']['fasta']
    output:
        bam_stats  =  os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.stats"),
        bam_flagstat  =  os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.flagstat"),
        bam_idxstats  =  os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.idxstats")
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "{wildcards.sample_id}: Running Samtools statistics"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.samtools_stats.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.samtools_stats.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.bam_stats})
        mkdir -p $(dirname {log})

        # Generate comprehensive BAM statistics
        samtools stats --threads {threads} --reference {input.fasta} {input.bam} > {output.bam_stats} 2>> {log} || {{ echo "[ERROR] samtools stats failed." >> {log}; exit 1; }}

        # Create flagstat summary
        samtools flagstat --threads {threads} {input.bam} > {output.bam_flagstat} 2>> {log} || {{ echo "[ERROR] samtools flagstat failed." >> {log}; exit 1; }}

        # Generate chromosome-level read mapping statistics
        samtools idxstats --threads {threads} {input.bam} > {output.bam_idxstats} 2>> {log} || {{ echo "[ERROR] samtools idxstats failed." >> {log}; exit 1; }}
        """

rule picard_collect_multiple_metrics:
    # Collect Picard metrics on filtered BAM.
    input:
        bam=os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        fasta=config["ref"]["fasta"]
    output:
        alignment_summary=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.alignment_summary_metrics"
        ),
        base_distribution_by_cycle_pdf=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.base_distribution_by_cycle.pdf"
        ),
        base_distribution_by_cycle_metrics=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.base_distribution_by_cycle_metrics"
        ),
        insert_size_histogram_pdf=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.insert_size_histogram.pdf"
        ),
        insert_size=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.insert_size_metrics"
        ),
        quality_by_cycle_pdf=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.quality_by_cycle.pdf"
        ),
        quality_by_cycle_metrics=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.quality_by_cycle_metrics"
        ),
        quality_distribution_pdf=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.quality_distribution.pdf"
        ),
        quality_distribution_metrics=os.path.join(
            "{outdir}",
            "bam",
            "{sample_id}.CollectMultipleMetrics.quality_distribution_metrics"
        )
    params:
        prefix=lambda wildcards: os.path.join(
            wildcards.outdir,
            "bam",
            f"{wildcards.sample_id}.CollectMultipleMetrics"
        ),
        tmp_dir = os.path.join("{outdir}", "bam", "tmp"),
    conda:
        os.path.join(workflow.basedir, "envs", "picard_markduplicates.yml")
    message:
        "{wildcards.sample_id}: Running Picard CollectMultipleMetrics"
    threads: 1
    resources:
        mem_mb = 16384,
        jvm_mem_mb = 14336  # mem_mb minus ~2 GB JVM overhead (metaspace, GC, native)
    log:
        os.path.join("{outdir}", "logs", "picard", "{sample_id}.collect_multiple_metrics.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.collect_multiple_metrics.benchmark.txt")
    shell:
        """
        mkdir -p "{params.tmp_dir}"
        mkdir -p "$(dirname "{log}")"

        picard -Xmx{resources.jvm_mem_mb}m CollectMultipleMetrics \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR="{params.tmp_dir}" \
            INPUT="{input.bam}" \
            OUTPUT="{params.prefix}" \
            REFERENCE_SEQUENCE="{input.fasta}" \
            > "{log}" 2>&1 || {{ echo "[ERROR] Picard CollectMultipleMetrics failed." >> "{log}"; exit 1; }}

        if [ -s "{output.alignment_summary}" ]; then
            rm -rf "{params.tmp_dir}"
        else
            echo "[ERROR] Missing Picard output: {output.alignment_summary}" >> "{log}"
            exit 1
        fi
        """
