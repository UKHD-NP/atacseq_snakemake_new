ATAQV_DESCRIPTION = "NA"
ATAQV_IGNORE_READ_GROUPS = True
ATAQV_MITO_NAME = str(config.get("ref", {}).get("mito_name", "chrM")).strip() or "chrM"


if is_enabled("ataqv") and not is_enabled("call_peaks"):
    fatal("ataqv requires `call_peaks.enabled: true` (peak file input).")


rule ataqv:
    # Build ataqv JSON metrics for ATAC-seq QC.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.bai"),
        peaks = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        tss = config["ref"]["tss"],
        autosomes = config["ref"]["autosomes"]
    output:
        json = os.path.join("{outdir}", "ataqv", "{sample_id}.ataqv.json")
    params:
        ignore_read_groups = "--ignore-read-groups" if ATAQV_IGNORE_READ_GROUPS else "",
        mito_name = ATAQV_MITO_NAME,
        description = ATAQV_DESCRIPTION,
        sample_name = lambda wildcards: f"{wildcards.sample_id}.filtered"
    conda:
        os.path.join(workflow.basedir, "envs", "ataqv.yml")
    message:
        "{wildcards.sample_id}: Running ataqv"
    threads: 1
    resources:
        mem_mb = 6144
    log:
        os.path.join("{outdir}", "logs", "ataqv", "{sample_id}.ataqv.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.ataqv.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.json}")"
        mkdir -p "$(dirname "{log}")"

        ataqv \
            {params.ignore_read_groups} \
            --mitochondrial-reference-name "{params.mito_name}" \
            --peak-file "{input.peaks}" \
            --tss-file "{input.tss}" \
            --autosomal-reference-file "{input.autosomes}" \
            --metrics-file "{output.json}" \
            --threads {threads} \
            --name "{params.sample_name}" \
            "{params.description}" \
            "{input.bam}" > "{log}" 2>&1 || {{
            echo "[ERROR] ataqv failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.json}" ]; then
            echo "[ERROR] ataqv output is missing or empty: {output.json}" >> "{log}"
            exit 1
        fi
        """


rule ataqv_mkarv:
    # Render ataqv interactive HTML report with mkarv.
    input:
        json = os.path.join("{outdir}", "ataqv", "{sample_id}.ataqv.json")
    output:
        index = os.path.join("{outdir}", "ataqv", "{sample_id}.mkarv_html", "index.html")
    params:
        html_dir = os.path.join("{outdir}", "ataqv", "{sample_id}.mkarv_html"),
        jsons_dir = os.path.join("{outdir}", "ataqv", "{sample_id}.mkarv_jsons"),
        json_link = os.path.join("{outdir}", "ataqv", "{sample_id}.mkarv_jsons", "{sample_id}.ataqv.json")
    conda:
        os.path.join(workflow.basedir, "envs", "ataqv.yml")
    message:
        "{wildcards.sample_id}: Building ataqv HTML report with mkarv"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join("{outdir}", "logs", "ataqv", "{sample_id}.mkarv.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.mkarv.benchmark.txt")
    shell:
        """
        mkdir -p "{params.html_dir}"
        mkdir -p "{params.jsons_dir}"
        mkdir -p "$(dirname "{log}")"

        ln -sf "$(readlink -f "{input.json}")" "{params.json_link}"

        mkarv \
            --concurrency {threads} \
            --force \
            "{params.html_dir}" \
            "{params.jsons_dir}"/* > "{log}" 2>&1 || {{
            echo "[ERROR] mkarv failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.index}" ]; then
            echo "[ERROR] mkarv output is missing or empty: {output.index}" >> "{log}"
            exit 1
        fi
        """


rule ataqv_score:
    # Extract TSS enrichment score and NFR ratio from ataqv JSON for MultiQC.
    input:
        ataqv_json = os.path.join("{outdir}", "ataqv", "{sample_id}.ataqv.json")
    output:
        mqc_tsv = os.path.join("{outdir}", "ataqv", "{sample_id}.ataqv_score.tsv")
    conda:
        os.path.join(workflow.basedir, "envs", "ataqv.yml")
    message:
        "{wildcards.sample_id}: Extracting NFR score and TSSE score via ataqv JSON for MultiQC"
    threads: 1
    resources:
        mem_mb = 256
    log:
        os.path.join("{outdir}", "logs", "ataqv", "{sample_id}.ataqv_score.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.ataqv_score.benchmark.txt")
    script:
        os.path.join(workflow.basedir, "scripts", "extract_ataqv_score.py")
