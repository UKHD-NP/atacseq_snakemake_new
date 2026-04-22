NFR_CFG = config.get("nfr", {})
NFR_MAX_FRAGMENT = int(NFR_CFG.get("nfr_max_fragment", 150))
MONO_MIN_FRAGMENT = int(NFR_CFG.get("mono_min_fragment", 151))
MONO_MAX_FRAGMENT = int(NFR_CFG.get("mono_max_fragment", 300))


rule nfr_bigwig_nfr:
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam.bai"),
        chromsizes = config["ref"]["chromsizes"]
    output:
        bigwig = os.path.join("{outdir}", "nfr", "{sample_id}.nfr.bigWig")
    params:
        gsize = get_gsize,
        max_fragment = NFR_MAX_FRAGMENT,
        bin_size = 10,
        tmp_bam = os.path.join("{outdir}", "nfr", "{sample_id}.nfr.tmp.bam"),
        tmp_sorted = os.path.join("{outdir}", "nfr", "{sample_id}.nfr.sorted.bam"),
        memory_per_thread = "1G"
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Creating NFR bigWig (fragments ≤ {params.max_fragment} bp)"
    threads: 8
    resources:
        mem_mb = 16384
    log:
        os.path.join("{outdir}", "logs", "nfr", "{sample_id}.nfr.bigwig.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.nfr.bigwig.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.bigwig}")"
        mkdir -p "$(dirname "{log}")"

        alignmentSieve \
            --numberOfProcessors {threads} \
            --maxFragmentLength {params.max_fragment} \
            --bam "{input.bam}" \
            -o "{params.tmp_bam}" > "{log}" 2>&1 || {{
            echo "[ERROR] alignmentSieve (NFR) failed." >> "{log}"
            exit 1
        }}

        samtools sort \
            --write-index \
            -m "{params.memory_per_thread}" \
            -@ {threads} \
            -o "{params.tmp_sorted}##idx##{params.tmp_sorted}.bai" \
            "{params.tmp_bam}" >> "{log}" 2>&1 || {{
            echo "[ERROR] samtools sort (NFR) failed." >> "{log}"
            exit 1
        }}
        rm -f "{params.tmp_bam}"

        bamCoverage \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing RPGC \
            --effectiveGenomeSize "{params.gsize}" \
            --bam "{params.tmp_sorted}" \
            -o "{output.bigwig}" >> "{log}" 2>&1 || {{
            echo "[ERROR] bamCoverage (NFR) failed." >> "{log}"
            exit 1
        }}
        rm -f "{params.tmp_sorted}" "{params.tmp_sorted}.bai"

        if [ ! -s "{output.bigwig}" ]; then
            echo "[ERROR] NFR bigWig missing or empty." >> "{log}"
            exit 1
        fi
        """


rule nfr_bigwig_mono:
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam.bai"),
        chromsizes = config["ref"]["chromsizes"]
    output:
        bigwig = os.path.join("{outdir}", "nfr", "{sample_id}.mono.bigWig")
    params:
        gsize = get_gsize,
        min_fragment = MONO_MIN_FRAGMENT,
        max_fragment = MONO_MAX_FRAGMENT,
        bin_size = 10,
        tmp_bam = os.path.join("{outdir}", "nfr", "{sample_id}.mono.tmp.bam"),
        tmp_sorted = os.path.join("{outdir}", "nfr", "{sample_id}.mono.sorted.bam"),
        memory_per_thread = "1G"
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Creating mononucleosomal bigWig ({params.min_fragment}–{params.max_fragment} bp)"
    threads: 8
    resources:
        mem_mb = 16384
    log:
        os.path.join("{outdir}", "logs", "nfr", "{sample_id}.mono.bigwig.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.mono.bigwig.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.bigwig}")"
        mkdir -p "$(dirname "{log}")"

        alignmentSieve \
            --numberOfProcessors {threads} \
            --minFragmentLength {params.min_fragment} \
            --maxFragmentLength {params.max_fragment} \
            --bam "{input.bam}" \
            -o "{params.tmp_bam}" > "{log}" 2>&1 || {{
            echo "[ERROR] alignmentSieve (mono) failed." >> "{log}"
            exit 1
        }}

        samtools sort \
            --write-index \
            -m "{params.memory_per_thread}" \
            -@ {threads} \
            -o "{params.tmp_sorted}##idx##{params.tmp_sorted}.bai" \
            "{params.tmp_bam}" >> "{log}" 2>&1 || {{
            echo "[ERROR] samtools sort (mono) failed." >> "{log}"
            exit 1
        }}
        rm -f "{params.tmp_bam}"

        bamCoverage \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing RPGC \
            --effectiveGenomeSize "{params.gsize}" \
            --bam "{params.tmp_sorted}" \
            -o "{output.bigwig}" >> "{log}" 2>&1 || {{
            echo "[ERROR] bamCoverage (mono) failed." >> "{log}"
            exit 1
        }}
        rm -f "{params.tmp_sorted}" "{params.tmp_sorted}.bai"

        if [ ! -s "{output.bigwig}" ]; then
            echo "[ERROR] Mono bigWig missing or empty." >> "{log}"
            exit 1
        fi
        """


rule nfr_compute_matrix:
    input:
        nfr_bigwig = os.path.join("{outdir}", "nfr", "{sample_id}.nfr.bigWig"),
        mono_bigwig = os.path.join("{outdir}", "nfr", "{sample_id}.mono.bigWig"),
        tss = config["ref"]["tss"]
    output:
        matrix = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.mat.gz"),
        values = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.vals.mat.tab")
    params:
        upstream = 2000,
        downstream = 2000,
        nfr_label = lambda wildcards: f"NFR (<{NFR_MAX_FRAGMENT}bp)",
        mono_label = lambda wildcards: f"Mono ({MONO_MIN_FRAGMENT}-{MONO_MAX_FRAGMENT}bp)"
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Computing NFR vs mono matrix around TSS"
    threads: 12
    resources:
        mem_mb = 6144
    log:
        os.path.join("{outdir}", "logs", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.nfr_vs_mono.computeMatrix.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.matrix}")"
        mkdir -p "$(dirname "{log}")"

        computeMatrix \
            reference-point \
            --missingDataAsZero \
            --skipZeros \
            --upstream {params.upstream} \
            --downstream {params.downstream} \
            --samplesLabel "{params.nfr_label}" "{params.mono_label}" \
            --regionsFileName "{input.tss}" \
            --scoreFileName "{input.nfr_bigwig}" "{input.mono_bigwig}" \
            --outFileName "{output.matrix}" \
            --outFileNameMatrix "{output.values}" \
            --numberOfProcessors {threads} > "{log}" 2>&1 || {{
            echo "[ERROR] computeMatrix (NFR vs mono) failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.matrix}" ] || [ ! -s "{output.values}" ]; then
            echo "[ERROR] NFR computeMatrix outputs missing or empty." >> "{log}"
            exit 1
        fi
        """


rule nfr_plot_profile:
    input:
        matrix = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.mat.gz")
    output:
        plot  = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.plotProfile.pdf"),
        table = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.plotProfile.tab")
    params:
        colors = "black red"
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Plotting NFR vs mono profile around TSS"
    threads: 2
    resources:
        mem_mb = 4096
    log:
        os.path.join("{outdir}", "logs", "nfr", "{sample_id}.nfr_vs_mono.plotProfile.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.nfr_vs_mono.plotProfile.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.plot}")"
        mkdir -p "$(dirname "{log}")"

        plotProfile \
            --matrixFile "{input.matrix}" \
            --outFileName "{output.plot}" \
            --outFileNameData "{output.table}" \
            --colors {params.colors} \
            --plotType lines > "{log}" 2>&1 || {{
            echo "[ERROR] plotProfile (NFR vs mono) failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.plot}" ] || [ ! -s "{output.table}" ]; then
            echo "[ERROR] plotProfile (NFR vs mono) outputs missing or empty." >> "{log}"
            exit 1
        fi
        """


rule nfr_plot_heatmap:
    input:
        matrix = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.mat.gz")
    output:
        plot  = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.plotHeatmap.pdf"),
        table = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.plotHeatmap.mat.tab")
    params:
        colors = "autumn_r autumn_r"
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Plotting NFR vs mono heatmap around TSS"
    threads: 2
    resources:
        mem_mb = 20480
    log:
        os.path.join("{outdir}", "logs", "nfr", "{sample_id}.nfr_vs_mono.plotHeatmap.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.nfr_vs_mono.plotHeatmap.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.plot}")"
        mkdir -p "$(dirname "{log}")"

        plotHeatmap \
            --matrixFile "{input.matrix}" \
            --outFileName "{output.plot}" \
            --outFileNameMatrix "{output.table}" \
            --colorMap {params.colors} > "{log}" 2>&1 || {{
            echo "[ERROR] plotHeatmap (NFR vs mono) failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.plot}" ] || [ ! -s "{output.table}" ]; then
            echo "[ERROR] plotHeatmap (NFR vs mono) outputs missing or empty." >> "{log}"
            exit 1
        fi
        """
