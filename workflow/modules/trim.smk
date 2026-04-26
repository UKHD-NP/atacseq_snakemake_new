TRIMMING_CFG = config.get("trimming", {}) if isinstance(config.get("trimming", {}), dict) else {}


TRIM_TOOL = str(TRIMMING_CFG.get("tool", "fastp")).strip().lower()


if TRIM_TOOL not in ["fastp", "trim_galore"]:
    raise ValueError(f"Invalid trimming tool: {TRIM_TOOL}. Use 'fastp' or 'trim_galore'.")


if TRIM_TOOL == "trim_galore":
    ruleorder: trim_galore > fastp
else:
    ruleorder: fastp > trim_galore


# Rule for Adapter and quality trimming
rule fastp:
    # Perform adapter trimming and quality control on FASTQ files
    input:
        get_paired_fq
    output:
        out1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1.fastq.gz"),
        out2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2.fastq.gz"),
        json = os.path.join("{outdir}", "trim", "{sample_id}.fastp.json"),
        html = os.path.join("{outdir}", "trim", "{sample_id}.fastp.html"),
        unpaired1 = os.path.join("{outdir}", "trim", "{sample_id}_1.fail.fastq.gz"),
        unpaired2 = os.path.join("{outdir}", "trim", "{sample_id}_2.fail.fastq.gz")
    params:
        fastp_params = TRIMMING_CFG.get("fastp_params", ""),
        # Compression level for output files (1-9, higher is more compressed but slower)
        compression_level = 4
    conda:
        os.path.join(workflow.basedir, "envs", "fastp.yml")
    message:
        "{wildcards.sample_id}: Trimming and performing quality control on paired-end FASTQ files"
    threads: 12
    resources:
        mem_mb = 16384
    log:
        os.path.join("{outdir}", "logs", "fastp", "{sample_id}.fastp.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "fastp.{sample_id}.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.out1})
        mkdir -p $(dirname {log})

        fastp \
            --in1 {input[0]} \
            --in2 {input[1]} \
            --out1 {output.out1} \
            --out2 {output.out2} \
            --json {output.json} \
            --html {output.html} \
            --unpaired1 {output.unpaired1} \
            --unpaired2 {output.unpaired2} \
            --thread {threads} \
            --compression {params.compression_level} \
            --detect_adapter_for_pe \
            --report_title "{wildcards.sample_id} fastp report" \
            {params.fastp_params} \
            > {log} 2>&1 || {{ echo "[ERROR] fastp failed." >> {log}; exit 1; }}
        """


rule trim_galore:
    # Trim reads using trim_galore (alternative to fastp)
    input:
        get_paired_fq
    output:
        out1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1.fastq.gz"),
        out2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2.fastq.gz"),
        report1 = os.path.join("{outdir}", "trim", "{sample_id}_1.fastq.gz_trimming_report.txt"),
        report2 = os.path.join("{outdir}", "trim", "{sample_id}_2.fastq.gz_trimming_report.txt")
    params:
        trim_galore_params = TRIMMING_CFG.get("trim_galore_params", ""),
        outdir = lambda wildcards: os.path.join(wildcards.outdir, "trim")
    conda:
        os.path.join(workflow.basedir, "envs", "trim_galore.yml")
    message:
        "{wildcards.sample_id}: Trimming reads with trim_galore"
    threads: 12
    resources:
        mem_mb = 36864
    log:
        os.path.join("{outdir}", "logs", "trim_galore", "{sample_id}.trim.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.trim_galore.benchmark.txt")
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        trim_galore \
            --paired \
            --gzip \
            --cores {threads} \
            --output_dir {params.outdir} \
            --basename {wildcards.sample_id} \
            {params.trim_galore_params} \
            "{input[0]}" "{input[1]}" \
            &> {log} || {{ echo "[ERROR] trim_galore failed." >> {log}; exit 1; }}

        # Rename outputs to match declared output names
        mv "{params.outdir}/{wildcards.sample_id}_val_1.fq.gz" "{output.out1}" || {{ echo "[ERROR] Missing trim_galore R1 output." >> {log}; exit 1; }}
        mv "{params.outdir}/{wildcards.sample_id}_val_2.fq.gz" "{output.out2}" || {{ echo "[ERROR] Missing trim_galore R2 output." >> {log}; exit 1; }}

        # Rename trimming reports (generated from input filename which contains '_merged_')
        R1_BASENAME=$(basename "{input[0]}")
        R2_BASENAME=$(basename "{input[1]}")
        mv "{params.outdir}/${{R1_BASENAME}}_trimming_report.txt" "{output.report1}" || {{ echo "[ERROR] Missing trim_galore R1 trimming report." >> {log}; exit 1; }}
        mv "{params.outdir}/${{R2_BASENAME}}_trimming_report.txt" "{output.report2}" || {{ echo "[ERROR] Missing trim_galore R2 trimming report." >> {log}; exit 1; }}
        """
