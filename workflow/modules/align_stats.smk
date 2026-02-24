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
    threads: 2
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
