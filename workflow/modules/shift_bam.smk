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
        tempdir = lambda wildcards: os.path.join(wildcards.outdir, "bam", f"tmp_{wildcards.sample_id}"),
        memory_per_thread = "4G"
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

        mkdir -p "{params.tempdir}"
        mkdir -p "$(dirname "{output.bigwig}")"
        mkdir -p "$(dirname "{log}")"
        export TMPDIR="{params.tempdir}"
 
        # Step 1: split by chromosome, shift in parallel
        echo "[INFO] alignmentSieve start: $(date)" >> "{log}"

        CHROMS=$(samtools idxstats "{input.bam}" | awk '$3 > 0 {{print $1}}')

        if [ -z "$CHROMS" ]; then
            echo "[ERROR] No chromosomes with mapped reads found in {input.bam}" >> "{log}"
            exit 1
        fi

        CHROM_COUNT=$(printf "%s\n" "$CHROMS" | grep -c .)
        echo "[INFO] Chromosomes with mapped reads ($CHROM_COUNT): $(printf "%s\n" "$CHROMS" | tr '\n' ' ')" >> "{log}"

        shift_chr() {{
            CHR=$1
            BAM=$2
            OUTDIR=$3
            LOG=$4

            echo "[INFO] START shift chromosome $CHR: $(date)" >> "$LOG"

            samtools view -b "$BAM" "$CHR" > "$OUTDIR/$CHR.bam" 2>> "$LOG" || {{
                echo "[ERROR] samtools view failed for $CHR: $(date)" >> "$LOG"
                exit 1
            }}

            samtools index "$OUTDIR/$CHR.bam" 2>> "$LOG" || {{
                echo "[ERROR] samtools index failed for $CHR: $(date)" >> "$LOG"
                exit 1
            }}

            alignmentSieve \
                --numberOfProcessors 1 \
                --ATACshift \
                --bam "$OUTDIR/$CHR.bam" \
                -o "$OUTDIR/$CHR.shifted.bam" \
                2>> "$LOG" || {{
                echo "[ERROR] alignmentSieve failed for $CHR: $(date)" >> "$LOG"
                exit 1
            }}

            if [ ! -s "$OUTDIR/$CHR.shifted.bam" ]; then
                echo "[ERROR] Missing shifted BAM for $CHR after alignmentSieve: $(date)" >> "$LOG"
                exit 1
            fi

            rm -f "$OUTDIR/$CHR.bam" "$OUTDIR/$CHR.bam.bai"
            echo "[INFO] DONE shift chromosome $CHR: $(date)" >> "$LOG"
        }}
        export -f shift_chr

        printf "%s\n" "$CHROMS" | xargs -r -P {threads} -I{{}} \
            bash -c 'shift_chr "$@"' _ {{}} "{input.bam}" "{params.tempdir}" "{log}" || {{
            echo "[ERROR] One or more chromosome shifts failed: $(date)" >> "{log}"
            exit 1
        }}

        # Check missing bam
        MISSING=0
        for CHR in $CHROMS; do
            if [ ! -s "{params.tempdir}/$CHR.shifted.bam" ]; then
                echo "[ERROR] Missing shifted BAM for $CHR" >> "{log}"
                MISSING=1
            fi
        done
        if [ "$MISSING" -ne 0 ]; then
            echo "[ERROR] Missing one or more chromosome shifted BAM files." >> "{log}"
            exit 1
        fi

        echo "[INFO] alignmentSieve done: $(date)" >> "{log}"

        # Step 2: merge + sort
        echo "[INFO] samtools merge start: $(date)" >> "{log}"
        samtools merge -f -@ 6 \
            "{params.tempdir}"/merged.bam \
            "{params.tempdir}"/*.shifted.bam \
            2>> "{log}" || {{
            echo "[ERROR] samtools merge failed." >> "{log}"
            exit 1
        }}

        echo "[INFO] samtools sort start: $(date)" >> "{log}"
        samtools sort \
            --write-index \
            -m "{params.memory_per_thread}" \
            -T "{params.tempdir}"/sort \
            -@ 6 \
            -o "{output.bam}##idx##{output.bai}" \
            "{params.tempdir}"/merged.bam \
            2>> "{log}" || {{
            echo "[ERROR] samtools sort failed." >> "{log}"
            exit 1
        }}
        echo "[INFO] samtools sort done: $(date)" >> "{log}"

        rm -rf "{params.tempdir}"

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
        echo "[INFO] bamCoverage done: $(date)" >> "{log}"

        if [ ! -s "{output.bigwig}" ]; then
            echo "[ERROR] Shifted bigWig output is missing or empty: {output.bigwig}" >> "{log}"
            exit 1
        fi
        """
