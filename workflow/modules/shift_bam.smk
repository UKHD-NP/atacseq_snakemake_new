rule shift_bam:
    # Shift ATAC cut sites using deepTools alignmentSieve.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.bai"),
        chromsizes = config["ref"]["chromsizes"]
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam.bai")
    params:
        tmp_dir = os.path.join("{outdir}", "bam", "tmp"),
        memory_per_thread = "3G",
        chunk_length = 100000000
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: ATAC-shifting BAM with alignmentSieve"
    threads: 26
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 73728,
        runtime = lambda wildcards, attempt: attempt * 2880
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.alignmentSieve.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.alignmentSieve.benchmark.txt")
    shell:
        """
        ulimit -Sn $(ulimit -Hn) 2>/dev/null || ulimit -n 65536 2>/dev/null || true
        echo "[INFO] File descriptor limit: $(ulimit -n)" >> "{log}"

        mkdir -p "{params.tmp_dir}"
        mkdir -p "$(dirname "{output.bam}")"
        mkdir -p "$(dirname "{log}")"
        export TMPDIR="{params.tmp_dir}"

        # Step 1: alignmentSieve on full BAM
        echo "[INFO] alignmentSieve start: $(date)" >> "{log}"

        alignmentSieve \
            --numberOfProcessors {threads} \
            --genomeChunkLength {params.chunk_length} \
            --ATACshift \
            --bam "{input.bam}" \
            -o "{params.tmp_dir}/shifted.unsorted.bam" \
            2>> "{log}" || {{
            echo "[ERROR] alignmentSieve failed: $(date)" >> "{log}"
            exit 1
        }}

        if [ ! -s "{params.tmp_dir}/shifted.unsorted.bam" ]; then
            echo "[ERROR] alignmentSieve output is missing or empty: $(date)" >> "{log}"
            exit 1
        fi

        echo "[INFO] alignmentSieve done: $(date)" >> "{log}"

        # Step 2: sort + index
        echo "[INFO] samtools sort start: $(date)" >> "{log}"

        samtools sort \
            --write-index \
            -m "{params.memory_per_thread}" \
            -T "{params.tmp_dir}/sort" \
            -@ {threads} \
            -o "{output.bam}##idx##{output.bai}" \
            "{params.tmp_dir}/shifted.unsorted.bam" \
            2>> "{log}" || {{
            echo "[ERROR] samtools sort failed: $(date)" >> "{log}"
            exit 1
        }}

        echo "[INFO] samtools sort done: $(date)" >> "{log}"

        rm -rf "{params.tmp_dir}"

        if [ ! -s "{output.bam}" ] || [ ! -s "{output.bai}" ]; then
            echo "[ERROR] Sorted BAM or BAI missing: $(date)" >> "{log}"
            exit 1
        fi
        """


rule shifted_bam_to_bigwig:
    # Build shifted bigWig from the shifted BAM.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam.bai"),
        chromsizes = config["ref"]["chromsizes"]
    output:
        bigwig = os.path.join("{outdir}", "bigwig", "{sample_id}.shifted.bigWig")
    params:
        gsize = get_gsize,
        bin_size = 10
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Building shifted bigWig with bamCoverage"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 36864,
        runtime = lambda wildcards, attempt: attempt * 240
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.shifted_bamCoverage.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.shifted_bamCoverage.benchmark.txt")
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
            --normalizeUsing RPGC \
            --effectiveGenomeSize "{params.gsize}" \
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
