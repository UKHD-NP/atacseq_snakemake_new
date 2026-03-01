rule mark_duplicates:
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.bam"),
        fasta = config['ref']['fasta']
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.markdup.sorted.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.markdup.sorted.bam.bai"),
        metrics = os.path.join("{outdir}", "bam", "{sample_id}.markdup.sorted.MarkDuplicates.metrics.txt")
    params:
        tmp_dir = os.path.join("{outdir}", "bam", "tmp"),
        max_records_in_ram = 500000,
        duplicate_scoring_strategy = "SUM_OF_BASE_QUALITIES",
        sorting_collection_size_ratio = 0.25
    conda:
        os.path.join(workflow.basedir, "envs", "picard_markduplicates.yml")
    message:
        "{wildcards.sample_id}: Marking duplicate reads"
    threads: 2
    resources:
        mem_mb = 48000
    log:
        os.path.join("{outdir}", "logs", "markdup", "{sample_id}.markduplicates.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.markdup.benchmark.txt")
    shell:
        """
        mkdir -p "{params.tmp_dir}"
        mkdir -p "$(dirname "{log}")"

        picard -Xmx{resources.mem_mb}m MarkDuplicates \
            SORTING_COLLECTION_SIZE_RATIO={params.sorting_collection_size_ratio} \
            ASSUME_SORTED=true \
            REMOVE_DUPLICATES=false \
            MAX_RECORDS_IN_RAM={params.max_records_in_ram} \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR="{params.tmp_dir}" \
            INPUT="{input.bam}" \
            OUTPUT="{output.bam}" \
            REFERENCE_SEQUENCE="{input.fasta}" \
            METRICS_FILE="{output.metrics}" \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
            DUPLICATE_SCORING_STRATEGY={params.duplicate_scoring_strategy} \
            PROGRAM_RECORD_ID=MarkDuplicates-ATACseq \
            > "{log}" 2>&1 || {{ echo "[ERROR] Picard MarkDuplicates failed." >> "{log}"; exit 1; }}

        if [ -s "{output.bam}" ]; then
            echo "[INFO] Indexing marked-duplicate BAM with samtools..." >> "{log}"
            samtools index -@ {threads} "{output.bam}" "{output.bai}" >> "{log}" 2>&1 || {{
                echo "[ERROR] samtools index failed." >> "{log}"
                exit 1
            }}
        fi

        if [ -s "{output.bam}" ] && [ -s "{output.bai}" ]; then
            echo "[INFO] Removing input BAM to save space." >> "{log}"
            rm -f "{input.bam}" "{input.bam}.bai"
            rm -rf "{params.tmp_dir}"
        else
            echo "[ERROR] MarkDuplicates or BAM indexing failed." >> "{log}"
            exit 1;
        fi
        """
