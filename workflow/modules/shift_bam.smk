rule shift_bam:
    # Shift ATAC cut sites using deepTools alignmentSieve.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.bai"),
        chromsizes = config["ref"]["chromsizes"]
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam.bai"),
        bigwig = os.path.join("{outdir}", "bigwig", "{sample_id}.shifted.bigWig")
    params:
        gsize = get_gsize,
        bin_size = 10,
        tempdir = os.path.join("{outdir}", "bam", "tmp"),
        memory_per_thread = "2G",
        unsorted_bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.unsorted.bam")
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: ATAC-shifting BAM with alignmentSieve"
    threads: 16
    resources:
        mem_mb = 40960
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.alignmentSieve.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.alignmentSieve.benchmark.txt")
    shell:
        """
        ulimit -Sn $(ulimit -Hn) 2>/dev/null || ulimit -n 65536 2>/dev/null || true
        echo "[INFO] File descriptor limit: $(ulimit -n)" >> "{log}"

        mkdir -p "$(dirname "{output.bam}")"
        mkdir -p "$(dirname "{output.bigwig}")"
        mkdir -p "{params.tempdir}"
        mkdir -p "$(dirname "{log}")"

        # Step 1: alignmentSieve → unsorted BAM
        echo "[INFO] alignmentSieve start: $(date)" >> "{log}"
        alignmentSieve \
            --numberOfProcessors {threads} \
            --ATACshift \
            --verbose \
            --bam "{input.bam}" \
            -o "{params.unsorted_bam}" \
            2>> "{log}" || {{
            echo "[ERROR] alignmentSieve failed." >> "{log}"
            exit 1
        }}
        echo "[INFO] alignmentSieve done: $(date)" >> "{log}"

        if [ ! -s "{params.unsorted_bam}" ]; then
            echo "[ERROR] Unsorted BAM missing or empty." >> "{log}"
            exit 1
        fi

        # Step 2: samtools sort
        echo "[INFO] samtools sort start: $(date)" >> "{log}"
        samtools sort \
            --write-index \
            -m "{params.memory_per_thread}" \
            -T "{params.tempdir}/{wildcards.sample_id}.shifted" \
            -@ {threads} \
            -o "{output.bam}##idx##{output.bai}" \
            "{params.unsorted_bam}" \
            2>> "{log}" || {{
            echo "[ERROR] samtools sort failed." >> "{log}"
            exit 1
        }}
        echo "[INFO] samtools sort done: $(date)" >> "{log}"

        rm -f "{params.unsorted_bam}"

        if [ ! -s "{output.bam}" ] || [ ! -s "{output.bai}" ]; then
            echo "[ERROR] Sorted BAM or BAI missing." >> "{log}"
            exit 1
        fi

        # Step 3: bamCoverage
        echo "[INFO] bamCoverage start: $(date)" >> "{log}"
        bamCoverage \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing RPGC \
            --effectiveGenomeSize "{params.gsize}" \
            --bam "{output.bam}" \
            -o "{output.bigwig}" \
            2>> "{log}" || {{
            echo "[ERROR] bamCoverage failed." >> "{log}"
            exit 1
        }}
        echo "[INFO] Step 3: bamCoverage done: $(date)" >> "{log}"

        if [ ! -s "{output.bigwig}" ]; then
            echo "[ERROR] BigWig missing or empty." >> "{log}"
            exit 1
        fi
        """