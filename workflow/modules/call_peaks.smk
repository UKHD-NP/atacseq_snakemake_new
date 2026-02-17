CALL_PEAKS_CFG = config.get("call_peaks", {}) if isinstance(config.get("call_peaks", {}), dict) else {}
CALL_PEAKS_PEAK_TYPE = str(CALL_PEAKS_CFG.get("peak_type", "narrow")).strip().lower()
if CALL_PEAKS_PEAK_TYPE not in {"narrow", "broad"}:
    CALL_PEAKS_PEAK_TYPE = "narrow"

CALL_PEAKS_MACS3_NARROW_PARAMS = cfg_str(
    CALL_PEAKS_CFG,
    "macs3_narrow_params",
    "--shift 75 --extsize 150 --keep-dup all --nomodel --call-summits -q 0.01",
)

CALL_PEAKS_MACS3_BROAD_PARAMS = cfg_str(
    CALL_PEAKS_CFG,
    "macs3_broad_params",
    "--shift 75 --extsize 150 --keep-dup all --nomodel --broad --broad-cutoff 0.1",
)
CALL_PEAKS_PEAKQC_INPUT_EXT = "broadPeak" if CALL_PEAKS_PEAK_TYPE == "broad" else "narrowPeak"


def get_macs3_params():
    """Resolve MACS3 parameter string for selected peak type."""
    return CALL_PEAKS_MACS3_BROAD_PARAMS if CALL_PEAKS_PEAK_TYPE == "broad" else CALL_PEAKS_MACS3_NARROW_PARAMS


def get_gsize(wildcards, input):
    """Resolve effective genome size from call_peaks.macs3_gsize or chromsizes."""
    call_peaks_gsize = cfg_str(CALL_PEAKS_CFG, "macs3_gsize", "")
    return get_effective_genome_size(input.chromsizes, call_peaks_gsize)


rule macs3_callpeak_tn5:
    # Convert BAM to BED, shift reads to Tn5 cut sites, and call peaks with MACS3.
    params:
        bed = os.path.join("{outdir}", "peaks", "{sample_id}.bed"),
        prefix = lambda wildcards: os.path.join(wildcards.outdir, "peaks", wildcards.sample_id),
        gsize = get_gsize,
        macs3_params = get_macs3_params(),
        raw_peak = lambda wildcards: os.path.join(
            wildcards.outdir,
            "peaks",
            f"{wildcards.sample_id}_peaks.{'broadPeak' if CALL_PEAKS_PEAK_TYPE == 'broad' else 'narrowPeak'}",
        ),
        peak_sanitized_tmp = lambda wildcards: os.path.join(wildcards.outdir, "peaks", f"{wildcards.sample_id}_peaks.tmp.peak")
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        chromsizes = config["ref"]["chromsizes"]
    output:
        tn5_bed = os.path.join("{outdir}", "peaks", "{sample_id}.tn5_shifted.bed"),
        peak = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        peaks_xls = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.xls")
    log:
        os.path.join("{outdir}", "logs", "macs3", "{sample_id}.callpeak.log")
    conda:
        os.path.join(workflow.basedir, "envs", "peakcalling.yml")
    threads: 12
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.macs3_callpeak_tn5.benchmark.txt")
    message:
        "{wildcards.sample_id}: Calling peaks with Tn5-shifted BED and MACS3"
    shell:
        """
        mkdir -p "$(dirname "{output.peak}")"
        mkdir -p "$(dirname "{log}")"

        bedtools bamtobed -i "{input.bam}" > "{params.bed}" 2> "{log}" || {{
            echo "[ERROR] bedtools bamtobed failed." >> "{log}"
            exit 1
        }}

        awk -F $'\\t' 'BEGIN{{OFS=FS}}{{if ($6 == "+"){{$2 = $2 + 4}} else if ($6 == "-"){{$3 = $3 - 5}} print $0}}' "{params.bed}" > "{output.tn5_bed}" 2>> "{log}" || {{
            echo "[ERROR] Tn5 shifting with awk failed." >> "{log}"
            exit 1
        }}

        macs3 callpeak \
            -t "{output.tn5_bed}" \
            -n "{params.prefix}" \
            -f BED \
            --gsize "{params.gsize}" \
            {params.macs3_params} >> "{log}" 2>&1 || {{
            echo "[ERROR] macs3 callpeak failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{params.raw_peak}" ]; then
            echo "[ERROR] MACS3 raw peak output missing: {params.raw_peak}" >> "{log}"
            exit 1
        fi

        awk 'BEGIN{{FS=OFS="\\t"}} !/^track/ && !/^browser/ && !/^#/' "{params.raw_peak}" > "{params.peak_sanitized_tmp}"
        mv "{params.peak_sanitized_tmp}" "{output.peak}"

        if [ ! -s "{output.peak}" ] || [ ! -s "{output.peaks_xls}" ]; then
            echo "[ERROR] MACS3 outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """


rule macs3_peak_qc_plot:
    # Build MACS peak QC summary table and distribution plots.
    params:
        script = os.path.join(workflow.basedir, "scripts", "plot_macs2_qc.r"),
        outdir = lambda wildcards: os.path.join(wildcards.outdir, "peaks"),
        outprefix = lambda wildcards: f"{wildcards.sample_id}.macs_peakqc",
        tmp_peak = lambda wildcards: os.path.join(
            wildcards.outdir,
            "peaks",
            f"{wildcards.sample_id}.macs_peakqc.input.{CALL_PEAKS_PEAKQC_INPUT_EXT}",
        ),
    input:
        peak = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
    output:
        summary = os.path.join("{outdir}", "peaks", "{sample_id}.macs_peakqc.summary.txt"),
        plots = os.path.join("{outdir}", "peaks", "{sample_id}.macs_peakqc.plots.pdf"),
    log:
        os.path.join("{outdir}", "logs", "macs3", "{sample_id}.peak_qc.log"),
    conda:
        os.path.join(workflow.basedir, "envs", "peakcalling.yml")
    threads: 4
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.macs3_peak_qc.benchmark.txt")
    message:
        "{wildcards.sample_id}: Plotting MACS peak QC"
    shell:
        """
        set -euo pipefail

        mkdir -p "{params.outdir}"
        mkdir -p "$(dirname "{log}")"

        cp "{input.peak}" "{params.tmp_peak}"

        Rscript "{params.script}" \
            -i "{params.tmp_peak}" \
            -s "{wildcards.sample_id}" \
            -o "{params.outdir}" \
            -p "{params.outprefix}" > "{log}" 2>&1 || {{
            echo "[ERROR] plot_macs2_qc.r failed." >> "{log}"
            rm -f "{params.tmp_peak}"
            exit 1
        }}

        rm -f "{params.tmp_peak}"

        if [ ! -s "{output.summary}" ] || [ ! -s "{output.plots}" ]; then
            echo "[ERROR] MACS peak QC outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """
