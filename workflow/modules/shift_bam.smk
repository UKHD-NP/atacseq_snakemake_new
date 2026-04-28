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
        tempdir = os.path.join("{outdir}", "bam"),
        tmp_bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.tmp.bam"),
        memory_per_thread = "1G"
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
        mkdir -p "$(dirname "{log}")"

        alignmentSieve \
            --numberOfProcessors {threads} \
            --ATACshift \
            --bam "{input.bam}" \
            -o "{params.tmp_bam}" >> "{log}" 2>&1 || {{
            echo "[ERROR] alignmentSieve failed." >> "{log}"
            exit 1
        }}

        samtools sort \
            --write-index \
            -m "{params.memory_per_thread}" \
            -T "{params.tempdir}/{wildcards.sample_id}.shifted" \
            -@ 6 \
            -o "{output.bam}##idx##{output.bai}" \
            "{params.tmp_bam}" >> "{log}" 2>&1 || {{
                echo "[ERROR] BAM sorting failed." >> "{log}"
                exit 1
            }}

        if [ -s "{output.bam}" ] && [ -s "{output.bai}" ]; then
            echo "[INFO] Removing unsorted BAM to save space." >> "{log}"
            rm -f "{params.tmp_bam}"
        else
            echo "[ERROR] Sorted BAM or BAI missing." >> "{log}"
            exit 1
        fi

        bamCoverage \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing RPGC \
            --effectiveGenomeSize "{params.gsize}" \
            --bam "{output.bam}" \
            -o "{output.bigwig}" >> "{log}" 2>&1 || {{
            echo "[ERROR] bamCoverage failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.bigwig}" ]; then
            echo "[ERROR] Shifted bigWig output is missing or empty: {output.bigwig}" >> "{log}"
            exit 1
        fi
        """