CALL_PEAKS_CFG = config.get("call_peaks", {}) if isinstance(config.get("call_peaks", {}), dict) else {}
CALL_PEAKS_PEAK_TYPE = str(CALL_PEAKS_CFG.get("peak_type", "narrow")).strip().lower()
if CALL_PEAKS_PEAK_TYPE not in {"narrow", "broad"}:
    CALL_PEAKS_PEAK_TYPE = "narrow"

CALL_PEAKS_MACS3_NARROW_PARAMS = cfg_str(
    CALL_PEAKS_CFG,
    "macs3_narrow_params",
    "--shift 75 --extsize 150 --keep-dup all --nomodel --call-summits -q 0.01",
)

# Broad peaks always use BAMPE (-f BAMPE, no --shift/--extsize needed).
CALL_PEAKS_MACS3_BROAD_PARAMS = cfg_str(
    CALL_PEAKS_CFG,
    "macs3_broad_params",
    "--keep-dup all --nomodel --broad --broad-cutoff 0.1",
)


def get_macs3_params():
    """Resolve MACS3 parameter string for selected peak type."""
    return CALL_PEAKS_MACS3_BROAD_PARAMS if CALL_PEAKS_PEAK_TYPE == "broad" else CALL_PEAKS_MACS3_NARROW_PARAMS


def get_gsize(wildcards, input):
    """Resolve effective genome size from call_peaks.macs3_gsize or chromsizes."""
    call_peaks_gsize = cfg_str(CALL_PEAKS_CFG, "macs3_gsize", "")
    return get_effective_genome_size(input.chromsizes, call_peaks_gsize)


def get_bam_for_callpeaks(wildcards):
    """
    Broad peaks:  filtered BAM with -f BAMPE (fragment-level; shifting not applicable).
    Narrow peaks: Tn5-shifted BAM (alignmentSieve).
    """
    if CALL_PEAKS_PEAK_TYPE == "broad":
        return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.filtered.bam")
    return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.shifted.bam")


def get_bai_for_callpeaks(wildcards):
    if CALL_PEAKS_PEAK_TYPE == "broad":
        return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.filtered.bam.bai")
    return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.shifted.bam.bai")


rule macs3_callpeak_tn5:
    # Call peaks with MACS3.
    # Narrow peaks: Tn5-shifted BAM (alignmentSieve) → bamtobed → macs3 BED mode.
    # Broad peaks:  filtered BAM → macs3 BAMPE fragment mode.
    input:
        bam = get_bam_for_callpeaks,
        bai = get_bai_for_callpeaks,
        chromsizes = config["ref"]["chromsizes"]
    output:
        tn5_bed = os.path.join("{outdir}", "peaks", "{sample_id}.tn5_shifted.bed"),
        peak = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
        peaks_xls = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.xls")
    params:
        prefix = lambda wildcards: os.path.join(wildcards.outdir, "peaks", wildcards.sample_id),
        gsize = get_gsize,
        macs3_params = get_macs3_params(),
        use_bampe = CALL_PEAKS_PEAK_TYPE == "broad",
        raw_peak = lambda wildcards: os.path.join(
            wildcards.outdir,
            "peaks",
            f"{wildcards.sample_id}_peaks.{'broadPeak' if CALL_PEAKS_PEAK_TYPE == 'broad' else 'narrowPeak'}",
        ),
        peak_sanitized_tmp = lambda wildcards: os.path.join(wildcards.outdir, "peaks", f"{wildcards.sample_id}_peaks.tmp.peak")
    conda:
        os.path.join(workflow.basedir, "envs", "peakcalling.yml")
    message:
        "{wildcards.sample_id}: Calling peaks with MACS3"
    threads: 12
    resources:
        mem_mb = 16384  # MACS3 peak calling typically needs 8-16 GB
    log:
        os.path.join("{outdir}", "logs", "macs3", "{sample_id}.callpeak.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.macs3_callpeak_tn5.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.peak}")"
        mkdir -p "$(dirname "{log}")"

        if [ "{params.use_bampe}" = "True" ]; then
            # BAMPE mode: call peaks directly from paired-end filtered BAM (broad peaks, fragment-level).
            touch "{output.tn5_bed}"

            macs3 callpeak \
                -t "{input.bam}" \
                -n "{params.prefix}" \
                -f BAMPE \
                --gsize "{params.gsize}" \
                {params.macs3_params} >> "{log}" 2>&1 || {{
                echo "[ERROR] macs3 callpeak (BAMPE) failed." >> "{log}"
                exit 1
            }}
        else
            # BED mode: BAM is already Tn5-shifted by alignmentSieve; convert to BED directly.
            bedtools bamtobed -i "{input.bam}" > "{output.tn5_bed}" 2> "{log}" || {{
                echo "[ERROR] bedtools bamtobed failed." >> "{log}"
                exit 1
            }}

            macs3 callpeak \
                -t "{output.tn5_bed}" \
                -n "{params.prefix}" \
                -f BED \
                --gsize "{params.gsize}" \
                {params.macs3_params} >> "{log}" 2>&1 || {{
                echo "[ERROR] macs3 callpeak (BED) failed." >> "{log}"
                exit 1
            }}
        fi

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
    input:
        peak = os.path.join("{outdir}", "peaks", "{sample_id}_peaks.peak"),
    output:
        summary = os.path.join("{outdir}", "peaks", "{sample_id}.macs_peakqc.summary.txt"),
        plots = os.path.join("{outdir}", "peaks", "{sample_id}.macs_peakqc.plots.pdf"),
    params:
        script = os.path.join(workflow.basedir, "scripts", "plot_macs_qc.r"),
        outprefix = "{sample_id}.macs_peakqc",
    conda:
        os.path.join(workflow.basedir, "envs", "peakcalling.yml")
    message:
        "{wildcards.sample_id}: Plotting MACS peak QC"
    threads: 4
    log:
        os.path.join("{outdir}", "logs", "macs3", "{sample_id}.peak_qc.log"),
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.macs3_peak_qc.benchmark.txt")
    shell:
        """
        set -euo pipefail

        mkdir -p "{wildcards.outdir}/peaks"
        mkdir -p "$(dirname "{log}")"

        Rscript "{params.script}" \
            -i "{input.peak}" \
            -s "{wildcards.sample_id}" \
            -o "{wildcards.outdir}/peaks" \
            -p "{params.outprefix}" > "{log}" 2>&1 || {{
            echo "[ERROR] plot_macs2_qc.r failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.summary}" ] || [ ! -s "{output.plots}" ]; then
            echo "[ERROR] MACS peak QC outputs are missing or empty." >> "{log}"
            exit 1
        fi
        """
