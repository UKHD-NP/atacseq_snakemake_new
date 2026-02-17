# Merge raw FASTQ lanes for duplicated sample IDs.
# For samples with one lane, this is a simple copy-through via cat.
rule merge_raw_fastqs:
    input:
        fq1 = get_raw_lane_fq1,
        fq2 = get_raw_lane_fq2
    output:
        fq1 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_2.fastq.gz")
    log:
        os.path.join("{outdir}", "logs", "raw_merge", "{sample_id}.merge_raw_fastqs.log")
    message:
        "{wildcards.sample_id}: Merging raw FASTQ lanes"
    shell:
        """
        mkdir -p $(dirname {output.fq1})
        mkdir -p $(dirname {log})
        rm -f {output.fq1} {output.fq2}

        N_R1=$(echo {input.fq1} | wc -w)
        N_R2=$(echo {input.fq2} | wc -w)

        if [ "$N_R1" -eq 1 ] && [ "$N_R2" -eq 1 ]; then
            # Single-lane sample: use symlink to avoid duplicate storage
            ln -sf "$(readlink -f {input.fq1})" {output.fq1}
            ln -sf "$(readlink -f {input.fq2})" {output.fq2}
            echo "[INFO] Mode: symlink (single lane)" > {log}
        else
            # Multi-lane sample: concatenate lanes in listed order
            cat {input.fq1} > {output.fq1}
            cat {input.fq2} > {output.fq2}
            echo "[INFO] Mode: merge (multi-lane)" > {log}
        fi

        if [ ! -s "{output.fq1}" ] || [ ! -s "{output.fq2}" ]; then
            echo "[ERROR] Merged FASTQ output is empty." >> {log}
            exit 1
        fi

        echo "[INFO] Merged R1 inputs: {input.fq1}" >> {log}
        echo "[INFO] Merged R2 inputs: {input.fq2}" >> {log}
        """

# Rule for FastQC on raw reads (before trimming)
rule fastqc_raw:
    # Run FastQC on raw FASTQ files before trimming
    params:
        outdir = lambda wildcards: os.path.join(wildcards.outdir, "fastqc_raw"),
        temp_dir = lambda wildcards: os.path.join(wildcards.outdir, "fastqc_raw", "temp")
    input:
        get_paired_fq
    output:
        html1 = os.path.join("{outdir}", "fastqc_raw", "{sample_id}_raw_1_fastqc.html"),
        html2 = os.path.join("{outdir}", "fastqc_raw", "{sample_id}_raw_2_fastqc.html"),
        zip1 = os.path.join("{outdir}", "fastqc_raw", "{sample_id}_raw_1_fastqc.zip"),
        zip2 = os.path.join("{outdir}", "fastqc_raw", "{sample_id}_raw_2_fastqc.zip")
    log:
        os.path.join("{outdir}", "logs", "fastqc_raw", "{sample_id}.fastqc_raw.log")
    conda:
        os.path.join(workflow.basedir, "envs", "trim_galore.yml")  # FastQC is included in trim_galore env
    threads: 2
    resources:
        mem_mb = 1024
    message:
        "{wildcards.sample_id}: Running FastQC on raw reads"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.temp_dir}

        # Create symlinks with 'raw' in the name so FastQC includes it in sample name
        ln -sf $(readlink -f {input[0]}) {params.temp_dir}/{wildcards.sample_id}_raw_1.fastq.gz
        ln -sf $(readlink -f {input[1]}) {params.temp_dir}/{wildcards.sample_id}_raw_2.fastq.gz

        # Run FastQC on symlinked files
        fastqc \
            --quiet \
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
    params:
        outdir = lambda wildcards: os.path.join(wildcards.outdir, "trim")
    input:
        fq1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2.fastq.gz")
    output:
        html1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1_fastqc.html"),
        html2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2_fastqc.html"),
        zip1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1_fastqc.zip"),
        zip2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2_fastqc.zip")
    log:
        os.path.join("{outdir}", "logs", "fastqc_trimmed", "{sample_id}.fastqc_trimmed.log")
    conda:
        os.path.join(workflow.basedir, "envs", "trim_galore.yml")
    threads: 2
    resources:
        mem_mb = 1024
    message:
        "{wildcards.sample_id}: Running FastQC on trimmed reads"
    shell:
        """
        fastqc \
            --quiet \
            --threads {threads} \
            --outdir {params.outdir} \
            {input.fq1} {input.fq2} \
            &> {log} || {{ echo "[ERROR] FastQC (trimmed) failed." >> {log}; exit 1; }}
        """

