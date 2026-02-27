def get_bam_for_filter(wildcards):
    """BAM source for filtering depends on markduplicates setting."""
    if is_enabled("markduplicates"):
        return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.markdup.sorted.bam")
    return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.bam")


bam_filter_cfg = config.get("bam_filter", {})
default_bam_filter_params = "-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -F 0x0400 -q 30"


rule bam_filter:
    input:
        bam = get_bam_for_filter,
        include_regions = config["ref"]["include_regions"]
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.bai")
    params:
        bam_filter_params = bam_filter_cfg.get("params", default_bam_filter_params),
        bamtools_script = os.path.join(workflow.basedir, "scripts", "bamtools_filter_pe.json"),
        bampe_rm_orphan_script = os.path.join(workflow.basedir, "scripts", "bampe_rm_orphan.py"),
        tempdir = os.path.join("{outdir}", "bam", "tmp"),
        memory_per_thread = "4G"
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "{wildcards.sample_id}: Filtering BAM"
    threads: 12
    resources:
        mem_mb = 49152  # 4 GB per thread (matches memory_per_thread param)
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.bam_filter.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bam_filter.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.bam}")"
        mkdir -p "$(dirname "{log}")"
        mkdir -p "{params.tempdir}"

        TMP_FILTERED="{params.tempdir}/{wildcards.sample_id}.filtered.tmp.bam"
        TMP_NAME_SORTED="{params.tempdir}/{wildcards.sample_id}.filtered.name_sorted.tmp.bam"
        TMP_CLEANED="{params.tempdir}/{wildcards.sample_id}.filtered.cleaned.name_sorted.tmp.bam"

        set -euo pipefail

        if command -v bamtools >/dev/null 2>&1 && [ -f "{params.bamtools_script}" ]; then
            {{
                samtools view \
                    {params.bam_filter_params} \
                    -L "{input.include_regions}" \
                    -b "{input.bam}" \
                    | bamtools filter \
                        -out "$TMP_FILTERED" \
                        -script "{params.bamtools_script}"
            }} > "{log}" 2>&1 || {{
                echo "[ERROR] BAM filtering (samtools+bamtools) failed." >> "{log}"
                exit 1
            }}
        else
            echo "[WARNING] bamtools not available or script not found, running samtools-only filter." > "{log}"
            samtools view \
                {params.bam_filter_params} \
                -L "{input.include_regions}" \
                -b "{input.bam}" \
                > "$TMP_FILTERED" 2>> "{log}" || {{
                echo "[ERROR] BAM filtering (samtools-only) failed." >> "{log}"
                exit 1
            }}
        fi

        if [ ! -s "$TMP_FILTERED" ]; then
            echo "[ERROR] Intermediate BAM is missing or empty: $TMP_FILTERED" >> "{log}"
            exit 1
        fi

        if [ -f "{params.bampe_rm_orphan_script}" ] && command -v python3 >/dev/null 2>&1 && python3 -c "import pysam" >/dev/null 2>&1; then
            samtools sort \
                -n \
                -m {params.memory_per_thread} \
                -T "{params.tempdir}/{wildcards.sample_id}.name_sort_tmp" \
                -@ {threads} \
                -o "$TMP_NAME_SORTED" \
                "$TMP_FILTERED" >> "{log}" 2>&1 || {{
                echo "[ERROR] Name sort before orphan removal failed." >> "{log}"
                exit 1
            }}

            python3 "{params.bampe_rm_orphan_script}" \
                "$TMP_NAME_SORTED" \
                "$TMP_CLEANED" \
                --only_fr_pairs >> "{log}" 2>&1 || {{
                echo "[ERROR] bampe_rm_orphan.py failed." >> "{log}"
                exit 1
            }}

            FINAL_INPUT="$TMP_CLEANED"
        else
            echo "[WARNING] Skip orphan cleanup (missing script/python3/pysam)." >> "{log}"
            FINAL_INPUT="$TMP_FILTERED"
        fi

        samtools sort \
            --write-index \
            -m {params.memory_per_thread} \
            -T "{params.tempdir}/{wildcards.sample_id}.filtered" \
            -@ {threads} \
            -o "{output.bam}##idx##{output.bai}" \
            "$FINAL_INPUT" >> "{log}" 2>&1 || {{
            echo "[ERROR] Final coordinate sort/index failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.bam}" ] || [ ! -s "{output.bai}" ]; then
            echo "[ERROR] Filtered BAM or BAI is missing." >> "{log}"
            exit 1
        fi

        echo "[INFO] Removing intermediate and input BAM files to save space." >> "{log}"
        rm -f "$TMP_FILTERED" "$TMP_NAME_SORTED" "$TMP_CLEANED" "{input.bam}" "{input.bam}.bai"
        rm -rf "{params.tempdir}"
        """
