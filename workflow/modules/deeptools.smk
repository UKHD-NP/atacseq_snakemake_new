rule deeptools_compute_matrix_scale_regions:
    input:
        regions = config["ref"]["bed"],
        bigwig = os.path.join("{outdir}", "bigwig", "{sample_id}.shifted.bigWig")
    output:
        matrix = os.path.join("{outdir}", "deeptools", "{sample_id}.scale_regions.computeMatrix.mat.gz"),
        values = os.path.join("{outdir}", "deeptools", "{sample_id}.scale_regions.computeMatrix.vals.mat.tab")
    params:
        region_body_length = 1000,
        flank_length = 3000
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Running deepTools computeMatrix (scale-regions)"
    threads: 12
    resources:
        mem_mb = 24576
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.computeMatrix.scale_regions.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.computeMatrix.scale_regions.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.matrix}")"
        mkdir -p "$(dirname "{log}")"

        computeMatrix \
            scale-regions \
            --regionBodyLength {params.region_body_length} \
            --beforeRegionStartLength {params.flank_length} \
            --afterRegionStartLength {params.flank_length} \
            --missingDataAsZero \
            --skipZeros \
            --smartLabels \
            --regionsFileName "{input.regions}" \
            --scoreFileName "{input.bigwig}" \
            --outFileName "{output.matrix}" \
            --outFileNameMatrix "{output.values}" \
            --numberOfProcessors {threads} > "{log}" 2>&1 || {{
            echo "[ERROR] computeMatrix (scale-regions) failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.matrix}" ] || [ ! -s "{output.values}" ]; then
            echo "[ERROR] computeMatrix outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """


rule deeptools_compute_matrix_reference_point:
    input:
        regions = config["ref"]["tss"],
        bigwig = os.path.join("{outdir}", "bigwig", "{sample_id}.shifted.bigWig")
    output:
        matrix = os.path.join("{outdir}", "deeptools", "{sample_id}.reference_point.computeMatrix.mat.gz"),
        values = os.path.join("{outdir}", "deeptools", "{sample_id}.reference_point.computeMatrix.vals.mat.tab")
    params:
        upstream = 3000,
        downstream = 3000
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Running deepTools computeMatrix (reference-point)"
    threads: 12
    resources:
        mem_mb = 24576
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.computeMatrix.reference_point.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.computeMatrix.reference_point.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.matrix}")"
        mkdir -p "$(dirname "{log}")"

        computeMatrix \
            reference-point \
            --missingDataAsZero \
            --skipZeros \
            --smartLabels \
            --upstream {params.upstream} \
            --downstream {params.downstream} \
            --regionsFileName "{input.regions}" \
            --scoreFileName "{input.bigwig}" \
            --outFileName "{output.matrix}" \
            --outFileNameMatrix "{output.values}" \
            --numberOfProcessors {threads} > "{log}" 2>&1 || {{
            echo "[ERROR] computeMatrix (reference-point) failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.matrix}" ] || [ ! -s "{output.values}" ]; then
            echo "[ERROR] computeMatrix outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """


rule deeptools_plot_profile:
    input:
        matrix = os.path.join("{outdir}", "deeptools", "{sample_id}.scale_regions.computeMatrix.mat.gz")
    output:
        plot = os.path.join("{outdir}", "deeptools", "{sample_id}.scale_regions.plotProfile.pdf"),
        table = os.path.join("{outdir}", "deeptools", "{sample_id}.scale_regions.plotProfile.tab")
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Running deepTools plotProfile"
    threads: 12
    resources:
        mem_mb = 8192
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.plotProfile.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.plotProfile.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.plot}")"
        mkdir -p "$(dirname "{log}")"

        plotProfile \
            --matrixFile "{input.matrix}" \
            --outFileName "{output.plot}" \
            --outFileNameData "{output.table}" > "{log}" 2>&1 || {{
            echo "[ERROR] plotProfile failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.plot}" ] || [ ! -s "{output.table}" ]; then
            echo "[ERROR] plotProfile outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """


rule deeptools_plot_heatmap:
    input:
        matrix = os.path.join("{outdir}", "deeptools", "{sample_id}.reference_point.computeMatrix.mat.gz")
    output:
        plot = os.path.join("{outdir}", "deeptools", "{sample_id}.reference_point.plotHeatmap.pdf"),
        table = os.path.join("{outdir}", "deeptools", "{sample_id}.reference_point.plotHeatmap.mat.tab")
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Running deepTools plotHeatmap"
    threads: 12
    resources:
        mem_mb = 8192
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.plotHeatmap.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.plotHeatmap.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.plot}")"
        mkdir -p "$(dirname "{log}")"

        plotHeatmap \
            --matrixFile "{input.matrix}" \
            --outFileName "{output.plot}" \
            --outFileNameMatrix "{output.table}" > "{log}" 2>&1 || {{
            echo "[ERROR] plotHeatmap failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.plot}" ] || [ ! -s "{output.table}" ]; then
            echo "[ERROR] plotHeatmap outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """


rule deeptools_plot_fingerprint:
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam")
    output:
        plot = os.path.join("{outdir}", "deeptools", "{sample_id}.plotFingerprint.pdf"),
        raw_counts = os.path.join("{outdir}", "deeptools", "{sample_id}.plotFingerprint.raw.txt"),
        qc_metrics = os.path.join("{outdir}", "deeptools", "{sample_id}.plotFingerprint.qcmetrics.txt")
    params:
        label = lambda wildcards: wildcards.sample_id,
        num_samples = 500000
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Running deepTools plotFingerprint"
    threads: 12
    resources:
        mem_mb = 16384
    log:
        os.path.join("{outdir}", "logs", "deeptools", "{sample_id}.plotFingerprint.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.plotFingerprint.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.plot}")"
        mkdir -p "$(dirname "{log}")"

        plotFingerprint \
            --skipZeros \
            --numberOfSamples {params.num_samples} \
            --labels "{params.label}" \
            --bamfiles "{input.bam}" \
            --plotFile "{output.plot}" \
            --outRawCounts "{output.raw_counts}" \
            --outQualityMetrics "{output.qc_metrics}" \
            --numberOfProcessors {threads} > "{log}" 2>&1 || {{
            echo "[ERROR] plotFingerprint failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.plot}" ] || [ ! -s "{output.raw_counts}" ] || [ ! -s "{output.qc_metrics}" ]; then
            echo "[ERROR] plotFingerprint outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """
