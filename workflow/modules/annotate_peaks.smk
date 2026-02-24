rule homer_annotate_peaks:
    # Annotate called peaks with HOMER annotatePeaks.pl.
    input:
        peaks = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        fasta = config["ref"]["fasta"],
        gtf   = config["ref"]["gtf"]
    output:
        annotation = os.path.join("{outdir}", "annotation", "{sample_id}_peaks.annotatePeaks.txt"),
        summary = os.path.join("{outdir}", "annotation", "{sample_id}.macs_annotatePeaks.summary.txt")
    params:
        plot_homer = os.path.join(workflow.basedir, "scripts", "plot_homer_annotatepeaks.r"),
    conda:
        os.path.join(workflow.basedir, "envs", "homer.yml")
    message:
        "{wildcards.sample_id}: Annotating peaks with HOMER"
    threads: 8
    log:
        os.path.join("{outdir}", "logs", "homer", "{sample_id}.annotatePeaks.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.homer_annotatePeaks.benchmark.txt")
    shell:
        """
        set -euo pipefail

        mkdir -p "$(dirname "{output.annotation}")"
        mkdir -p "$(dirname "{log}")"

        annotatePeaks.pl \
            "{input.peaks}" \
            "{input.fasta}" \
            -gid \
            -gtf "{input.gtf}" \
            -cpu {threads} \
            > "{output.annotation}" 2> "{log}" || {{
            echo "[ERROR] annotatePeaks.pl failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.annotation}" ]; then
            echo "[ERROR] HOMER annotation output missing or empty: {output.annotation}" >> "{log}"
            exit 1
        fi

        Rscript "{params.plot_homer}" \
            -i "{output.annotation}" \
            -s "{wildcards.sample_id}" \
            -p "{wildcards.sample_id}.macs_annotatePeaks" \
            -o "{wildcards.outdir}/annotation" \
            >> "{log}" 2>&1 || {{
            echo "[ERROR] plot_homer_annotatepeaks.r failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.summary}" ]; then
            echo "[ERROR] Missing HOMER summary output: {output.summary}" >> "{log}"
            exit 1
        fi
        """
