# Rule for FastQC on raw reads (before trimming)
rule fastqc_raw:
    # Run FastQC on raw FASTQ files before trimming
    input:
        get_paired_fq
    output:
        html1 = os.path.join("{outdir}", "fastqc_raw", "{sample_id}_raw_1_fastqc.html"),
        html2 = os.path.join("{outdir}", "fastqc_raw", "{sample_id}_raw_2_fastqc.html"),
        zip1 = os.path.join("{outdir}", "fastqc_raw", "{sample_id}_raw_1_fastqc.zip"),
        zip2 = os.path.join("{outdir}", "fastqc_raw", "{sample_id}_raw_2_fastqc.zip")
    params:
        outdir = lambda wildcards: os.path.join(wildcards.outdir, "fastqc_raw"),
        temp_dir = lambda wildcards: os.path.join(wildcards.outdir, "fastqc_raw", "temp")
    conda:
        os.path.join(workflow.basedir, "envs", "trim_galore.yml")  # FastQC is included in trim_galore env
    message:
        "{wildcards.sample_id}: Running FastQC on raw reads"
    threads: 2
    resources:
        mem_mb = 4096
    log:
        os.path.join("{outdir}", "logs", "fastqc_raw", "{sample_id}.fastqc_raw.log")
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.temp_dir}

        # Create symlinks with 'raw' in the name so FastQC includes it in sample name
        ln -sf $(readlink -f {input[0]}) {params.temp_dir}/{wildcards.sample_id}_raw_1.fastq.gz
        ln -sf $(readlink -f {input[1]}) {params.temp_dir}/{wildcards.sample_id}_raw_2.fastq.gz

        # Run FastQC on symlinked files
        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            {params.temp_dir}/{wildcards.sample_id}_raw_1.fastq.gz \
            {params.temp_dir}/{wildcards.sample_id}_raw_2.fastq.gz \
            &> {log} || {{ echo "[ERROR] FastQC (raw) failed." >> {log}; exit 1; }}

        # Clean up temp files and directory
        rm -f {params.temp_dir}/{wildcards.sample_id}_raw_1.fastq.gz
        rm -f {params.temp_dir}/{wildcards.sample_id}_raw_2.fastq.gz
        rmdir {params.temp_dir} || true
        """


# Rule for FastQC on trimmed reads (after trimming) - only for trim_galore
rule fastqc_trimmed:
    # Run FastQC on trimmed FASTQ files after trim_galore
    input:
        fq1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2.fastq.gz")
    output:
        html1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1_fastqc.html"),
        html2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2_fastqc.html"),
        zip1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1_fastqc.zip"),
        zip2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2_fastqc.zip")
    params:
        outdir = lambda wildcards: os.path.join(wildcards.outdir, "trim")
    conda:
        os.path.join(workflow.basedir, "envs", "trim_galore.yml")
    message:
        "{wildcards.sample_id}: Running FastQC on trimmed reads"
    threads: 2
    resources:
        mem_mb = 4096
    log:
        os.path.join("{outdir}", "logs", "fastqc_trimmed", "{sample_id}.fastqc_trimmed.log")
    shell:
        """
        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            {input.fq1} {input.fq2} \
            &> {log} || {{ echo "[ERROR] FastQC (trimmed) failed." >> {log}; exit 1; }}
        """


