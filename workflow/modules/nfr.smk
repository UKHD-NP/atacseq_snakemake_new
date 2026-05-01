NFR_CFG = config.get("nfr", {})
NFR_MAX_FRAGMENT  = int(NFR_CFG.get("nfr_max_fragment",  150))
MONO_MIN_FRAGMENT = int(NFR_CFG.get("mono_min_fragment", 151))
MONO_MAX_FRAGMENT = int(NFR_CFG.get("mono_max_fragment", 300))
DI_MIN_FRAGMENT   = int(NFR_CFG.get("di_min_fragment",   301))
DI_MAX_FRAGMENT   = int(NFR_CFG.get("di_max_fragment",   500))
TRI_MIN_FRAGMENT  = int(NFR_CFG.get("tri_min_fragment",  501))
TRI_MAX_FRAGMENT  = int(NFR_CFG.get("tri_max_fragment",  700))


rule nfr_fragment_counts:
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam.bai")
    output:
        mqc = os.path.join("{outdir}", "nfr", "{sample_id}.fragment_counts_mqc.tsv")
    params:
        nfr_max  = NFR_MAX_FRAGMENT,
        mono_min = MONO_MIN_FRAGMENT,
        mono_max = MONO_MAX_FRAGMENT,
        di_min   = DI_MIN_FRAGMENT,
        di_max   = DI_MAX_FRAGMENT,
        tri_min  = TRI_MIN_FRAGMENT,
        tri_max  = TRI_MAX_FRAGMENT
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "{wildcards.sample_id}: Counting NFR, mono, di, and trinucleosomal reads"
    threads: 2
    resources:
        mem_mb = 2048,
        runtime = lambda wildcards, attempt: attempt * 240
    log:
        os.path.join("{outdir}", "logs", "nfr", "{sample_id}.fragment_counts.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.fragment_counts.benchmark.txt")
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.mqc}")"
        mkdir -p "$(dirname "{log}")"

        samtools view \
            --threads {threads} \
            -f 66 -F 4 \
            "{input.bam}" \
        | awk \
            -v nfr_max={params.nfr_max} \
            -v mono_min={params.mono_min} \
            -v mono_max={params.mono_max} \
            -v di_min={params.di_min} \
            -v di_max={params.di_max} \
            -v tri_min={params.tri_min} \
            -v tri_max={params.tri_max} \
            -v sample="{wildcards.sample_id}" \
            -v outfile="{output.mqc}" \
            -v log="{log}" \
        '
        $9 > 0 {{
            total++
            if      ($9 <= nfr_max)                      nfr++
            else if ($9 >= mono_min && $9 <= mono_max)   mono++
            else if ($9 >= di_min   && $9 <= di_max)     di++
            else if ($9 >= tri_min  && $9 <= tri_max)    tri++
        }}
        END {{
            nfr_pct  = (total > 0) ? nfr  / total * 100 : 0
            mono_pct = (total > 0) ? mono / total * 100 : 0
            di_pct   = (total > 0) ? di   / total * 100 : 0
            tri_pct  = (total > 0) ? tri  / total * 100 : 0
            print "Sample\tnfr_reads\tmono_reads\tdi_reads\ttri_reads\tnfr_pct\tmono_pct\tdi_pct\ttri_pct" > outfile
            printf "%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n",
                sample, nfr, mono, di, tri, nfr_pct, mono_pct, di_pct, tri_pct >> outfile
            printf "[INFO] NFR: %d (%.2f%%) | Mono: %d (%.2f%%) | Di: %d (%.2f%%) | Tri: %d (%.2f%%) | Total: %d\n",
                nfr, nfr_pct, mono, mono_pct, di, di_pct, tri, tri_pct, total >> log
        }}
        ' 2>> "{log}"

        if [ ! -s "{output.mqc}" ]; then
            echo "[ERROR] fragment_counts output missing or empty." >> "{log}"
            exit 1
        fi
        """


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
        mem_mb = 16384,
        runtime = lambda wildcards, attempt: attempt * 480
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
        mem_mb = 16384,
        runtime = lambda wildcards, attempt: attempt * 480
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
        matrix = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.gz"),
        values = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.tab")
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
        mem_mb = 6144,
        runtime = lambda wildcards, attempt: attempt * 240
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
        matrix = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.gz")
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
        mem_mb = 4096,
        runtime = lambda wildcards, attempt: attempt * 120
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
        matrix = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.computeMatrix.gz")
    output:
        plot  = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.plotHeatmap.pdf"),
        table = os.path.join("{outdir}", "nfr", "{sample_id}.nfr_vs_mono.plotHeatmap.tab")
    params:
        colors = "autumn_r autumn_r"
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Plotting NFR vs mono heatmap around TSS"
    threads: 2
    resources:
        mem_mb = 20480,
        runtime = lambda wildcards, attempt: attempt * 240
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
