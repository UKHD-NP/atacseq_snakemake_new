# =============================================================================
# ALIGN CONFIG
# =============================================================================

config.setdefault("ref", {})

ALIGN_CFG = config.get("align", {}) if isinstance(config.get("align", {}), dict) else {}
ALIGN_TOOL = str(ALIGN_CFG.get("tool", "bwa")).strip().lower()

if ALIGN_TOOL not in ["bwa", "bowtie2"]:
    raise ValueError(
        f"Invalid align.tool: '{ALIGN_TOOL}'. Allowed values: 'bwa' or 'bowtie2'."
    )

# Resolve ambiguity: both aligners produce the same BAM path
# (Snakemake needs this if both rules can create the same output)
if ALIGN_TOOL == "bowtie2":
    ruleorder: bowtie2_align > bwa_mem2_align
else:
    ruleorder: bwa_mem2_align > bowtie2_align


# =============================================================================
# ALIGN PARAMS
# =============================================================================

BWA_PARAMS = str(ALIGN_CFG.get("bwa_params", "")).strip()
BOWTIE2_PARAMS = str(ALIGN_CFG.get("bowtie2_params", "")).strip()


# =============================================================================
# INDEX PREFIXES
# =============================================================================

# ---- BWA index prefix ----
BWA_INDEX_PREFIX = str(config["ref"].get("bwa_index", "")).strip()
if not BWA_INDEX_PREFIX:
    BWA_INDEX_PREFIX = os.path.join(ref_dir, str(config["ref"]["assembly"]).strip())
config["ref"]["bwa_index"] = BWA_INDEX_PREFIX


# ---- Bowtie2 index prefix ----
BOWTIE2_INDEX_PREFIX = str(config["ref"].get("bowtie2_index", "")).strip()
if not BOWTIE2_INDEX_PREFIX:
    BOWTIE2_INDEX_PREFIX = os.path.join(
        ref_dir, str(config["ref"]["assembly"]).strip() + ".bowtie2"
    )
config["ref"]["bowtie2_index"] = BOWTIE2_INDEX_PREFIX


# =============================================================================
# BWA
# =============================================================================

rule bwa_mem2_index:
    # Build BWA-MEM2 index for the configured reference FASTA.
    input:
        fasta = config["ref"]["fasta"]
    output:
        idx_amb  = BWA_INDEX_PREFIX + ".amb",
        idx_ann  = BWA_INDEX_PREFIX + ".ann",
        idx_bwt  = BWA_INDEX_PREFIX + ".bwt.2bit.64",
        idx_pac  = BWA_INDEX_PREFIX + ".pac",
        idx_0123 = BWA_INDEX_PREFIX + ".0123"
    params:
        prefix = BWA_INDEX_PREFIX
    conda:
        os.path.join(workflow.basedir, "envs", "bwa-mem2.yml")
    message:
        "Building BWA-MEM2 index"
    threads: 8
    log:
        os.path.join(ref_dir, "bwa_mem2_index.log")
    benchmark:
        os.path.join(ref_dir, "bwa_mem2_index.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.idx_bwt})
        mkdir -p $(dirname {log})

        bwa-mem2 index -p {params.prefix} {input.fasta} > {log} 2>&1 || {{
            echo "[ERROR] bwa-mem2 index failed." >> {log}
            exit 1
        }}
        """


rule bwa_mem2_align:
    # Align paired-end reads with BWA-MEM2 and emit unsorted BAM.
    input:
        fq       = get_paired_trimmed_fq,
        idx_amb  = BWA_INDEX_PREFIX + ".amb",
        idx_ann  = BWA_INDEX_PREFIX + ".ann",
        idx_bwt  = BWA_INDEX_PREFIX + ".bwt.2bit.64",
        idx_pac  = BWA_INDEX_PREFIX + ".pac",
        idx_0123 = BWA_INDEX_PREFIX + ".0123"
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.unsorted.bam")
    params:
        bwa_params   = BWA_PARAMS,
        index_prefix = BWA_INDEX_PREFIX,
        rg = lambda wc: (
            f"@RG\\tID:{wc.sample_id}\\tSM:{wc.sample_id}"
            f"\\tPL:ILLUMINA\\tLB:{wc.sample_id}\\tPU:1"
        )
    conda:
        os.path.join(workflow.basedir, "envs", "bwa-mem2.yml")
    message:
        "{wildcards.sample_id}: Aligning with bwa-mem2"
    threads: 12
    resources:
        mem_mb = 49152  # ~4 GB per thread
    log:
        os.path.join("{outdir}", "logs", "bwa_mem2", "{sample_id}.align.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bwa_mem2_align.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})

        set -euo pipefail
        bwa-mem2 mem \
            {params.bwa_params} \
            -M \
            -R '{params.rg}' \
            -t {threads} \
            {params.index_prefix} \
            {input.fq[0]} {input.fq[1]} \
            | samtools view -bhS -O BAM --threads {threads} -o {output.bam} - \
            2>> {log} || {{
                echo "[ERROR] bwa-mem2 alignment failed." >> {log}
                exit 1
            }}

        if [ ! -s {output.bam} ]; then
            echo "[ERROR] Missing BWA output BAM: {output.bam}" >> {log}
            exit 1
        fi
        """


# =============================================================================
# BOWTIE2
# =============================================================================

rule bowtie2_index:
    # Build Bowtie2 index for the configured reference FASTA.
    input:
        fasta = config["ref"]["fasta"]
    output:
        idx_1    = BOWTIE2_INDEX_PREFIX + ".1.bt2",
        idx_2    = BOWTIE2_INDEX_PREFIX + ".2.bt2",
        idx_3    = BOWTIE2_INDEX_PREFIX + ".3.bt2",
        idx_4    = BOWTIE2_INDEX_PREFIX + ".4.bt2",
        idx_rev1 = BOWTIE2_INDEX_PREFIX + ".rev.1.bt2",
        idx_rev2 = BOWTIE2_INDEX_PREFIX + ".rev.2.bt2"
    params:
        prefix = BOWTIE2_INDEX_PREFIX
    conda:
        os.path.join(workflow.basedir, "envs", "bowtie2.yml")
    message:
        "Building Bowtie2 index"
    threads: 8
    log:
        os.path.join(ref_dir, "bowtie2_index.log")
    benchmark:
        os.path.join(ref_dir, "bowtie2_index.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.idx_1})
        mkdir -p $(dirname {log})

        bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1 || {{
            echo "[ERROR] bowtie2-build failed." >> {log}
            exit 1
        }}
        """


rule bowtie2_align:
    # Align paired-end reads with Bowtie2 and emit unsorted BAM.
    input:
        fq       = get_paired_trimmed_fq,
        idx_1    = BOWTIE2_INDEX_PREFIX + ".1.bt2",
        idx_2    = BOWTIE2_INDEX_PREFIX + ".2.bt2",
        idx_3    = BOWTIE2_INDEX_PREFIX + ".3.bt2",
        idx_4    = BOWTIE2_INDEX_PREFIX + ".4.bt2",
        idx_rev1 = BOWTIE2_INDEX_PREFIX + ".rev.1.bt2",
        idx_rev2 = BOWTIE2_INDEX_PREFIX + ".rev.2.bt2"
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.unsorted.bam")
    params:
        bowtie2_params = BOWTIE2_PARAMS,
        index_prefix   = BOWTIE2_INDEX_PREFIX,
        rg = lambda wc: (
            f"--rg-id {wc.sample_id} "
            f"--rg SM:{wc.sample_id} "
            f"--rg LB:{wc.sample_id} "
            f"--rg PU:1 "
            f"--rg PL:ILLUMINA"
        )
    conda:
        os.path.join(workflow.basedir, "envs", "bowtie2.yml")
    message:
        "{wildcards.sample_id}: Aligning with bowtie2"
    threads: 12
    resources:
        mem_mb = 24576  # ~2 GB per thread
    log:
        os.path.join("{outdir}", "logs", "bowtie2", "{sample_id}.align.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bowtie2_align.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})

        set -euo pipefail
        bowtie2 \
            {params.bowtie2_params} \
            --threads {threads} \
            -x {params.index_prefix} \
            -1 {input.fq[0]} \
            -2 {input.fq[1]} \
            {params.rg} \
            | samtools view -bhS -O BAM --threads {threads} -o {output.bam} - \
            2>> {log} || {{
                echo "[ERROR] bowtie2 alignment failed." >> {log}
                exit 1
            }}

        if [ ! -s {output.bam} ]; then
            echo "[ERROR] Missing Bowtie2 output BAM: {output.bam}" >> {log}
            exit 1
        fi
        """


rule sort_bam:
    # Sort BAM file by coordinates and index.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.unsorted.bam")
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.bam.bai")
    params:
        tempdir           = os.path.join("{outdir}", "bam"),
        memory_per_thread = "4G"
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "{wildcards.sample_id}: Sorting BAM by coordinates"
    threads: 10
    resources:
        mem_mb = 40960  # 4 GB per thread (matches memory_per_thread param)
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.sort.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bam_sort.benchmark.txt")
    shell:
        """
        mkdir -p {params.tempdir}
        mkdir -p $(dirname {log})

        samtools sort \
            --write-index \
            -m {params.memory_per_thread} \
            -T {params.tempdir}/{wildcards.sample_id}.sorted \
            -@ {threads} \
            -o {output.bam}##idx##{output.bai} \
            {input.bam} > {log} 2>&1 || {{
                echo "[ERROR] BAM sorting failed." >> {log}
                exit 1
            }}

        if [ -s {output.bam} ] && [ -s {output.bai} ]; then
            echo "[INFO] Removing unsorted BAM to save space." >> {log}
            rm -f {input.bam}
        else
            echo "[ERROR] Sorted BAM or BAI missing." >> {log}
            exit 1
        fi
        """
