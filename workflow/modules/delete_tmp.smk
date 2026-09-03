# shifted.bam only exists for narrow peaks with shift_bam on (see common.smk).
# Default: delete it once every consumer is done (regenerable from filtered.bam).
DELETE_SHIFTED_BAM_ON = (
    CALL_PEAKS_PEAK_TYPE == "narrow"
    and is_enabled("shift_bam", default=True)
    and as_bool(config.get("shift_bam", {}).get("delete_shifted_bam", True), default=True)
)
DELETE_TRIMMING_ON = is_enabled("trimming") and as_bool(
    config.get("trimming", {}).get("delete_trimming", True)
)


def _shifted_bam_consumer_outputs(wildcards):
    """Outputs that must exist before shifted.bam can be safely deleted - every rule
    that reads shifted.bam directly (shifted bigWig, NFR, ATACseqQC). Empty unless
    delete_shifted_bam is on, so no extra DAG edges are added in the default case."""
    if not DELETE_SHIFTED_BAM_ON:
        return []
    outdir = wildcards.outdir
    sample_id = wildcards.sample_id
    targets = [os.path.join(outdir, "bigwig", f"{sample_id}.shifted.bigWig")]
    if is_enabled("atacseqqc"):
        targets.append(os.path.join(outdir, "atacseqqc", f"{sample_id}.atacseqqc_mqc.tsv"))
    if as_bool(config.get("nfr", {}).get("fragment_counts", True)):
        targets.append(os.path.join(outdir, "nfr", f"{sample_id}.fragment_counts_mqc.tsv"))
    if is_enabled("nfr", default=True):
        targets.extend([
            os.path.join(outdir, "nfr", f"{sample_id}.nfr.bigWig"),
            os.path.join(outdir, "nfr", f"{sample_id}.mono.bigWig"),
        ])
    return targets


rule delete_tmp:
    # Clean up temporary files
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.bai"),
        fastqc = lambda wildcards: [
            os.path.join(wildcards.outdir, "trim", f"{wildcards.sample_id}_trimmed_1_fastqc.zip"),
            os.path.join(wildcards.outdir, "trim", f"{wildcards.sample_id}_trimmed_2_fastqc.zip"),
        ] if is_enabled("trimming") and TRIM_TOOL == "trim_galore" else [],
        shifted_bam_consumers = _shifted_bam_consumer_outputs
    output:
        log = os.path.join("{outdir}", "logs", "{sample_id}.deletion.log")
    params:
        fq1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2.fastq.gz"),
        fq1_fail = os.path.join("{outdir}", "trim", "{sample_id}_1.fail.fastq.gz"),
        fq2_fail = os.path.join("{outdir}", "trim", "{sample_id}_2.fail.fastq.gz"),
        raw_fq1 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_1.fastq.gz"),
        raw_fq2 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_2.fastq.gz"),
        raw_dir = os.path.join("{outdir}", "raw_merged"),
        delete_trimming = str(DELETE_TRIMMING_ON).lower(),
        shifted_bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam") if DELETE_SHIFTED_BAM_ON else "",
        shifted_bai = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam.bai") if DELETE_SHIFTED_BAM_ON else "",
    message:
        "{wildcards.sample_id}: Cleaning up temporary files"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join("{outdir}", "logs", "cleanup", "{sample_id}.cleanup.log")
    shell:
        """
        mkdir -p "$(dirname "{log}")"
        echo "[INFO] Starting cleanup for {wildcards.sample_id}" > "{log}"

        # Remove trimmed FASTQs
        if [ "{params.delete_trimming}" = "true" ]; then
            for f in "{params.fq1}" "{params.fq2}" "{params.fq1_fail}" "{params.fq2_fail}"; do
                [ -f "$f" ] && rm -f "$f" && echo "[INFO] Removed $f" >> "{log}"
            done
        else
            echo "[INFO] delete_trimming=false, keeping trimmed FASTQs." >> "{log}"
        fi

        # Remove merged raw FASTQs
        for f in "{params.raw_fq1}" "{params.raw_fq2}"; do
            { [ -e "$f" ] || [ -L "$f" ]; } && rm -f "$f" && echo "[INFO] Removed $f" >> "{log}"
        done

        # Shifted BAM (only when shift_bam.delete_shifted_bam=true; consumers are
        # already done via the shifted_bam_consumers input)
        if [ -n "{params.shifted_bam}" ]; then
            for f in "{params.shifted_bam}" "{params.shifted_bai}"; do
                [ -f "$f" ] && rm -f "$f" && echo "[INFO] Removed $f" >> "{log}"
            done
        else
            echo "[INFO] delete_shifted_bam=false, keeping shifted BAM." >> "{log}"
        fi

        # Drop raw_merged dir if now empty
        if [ -d "{params.raw_dir}" ]; then
            rmdir "{params.raw_dir}" 2>/dev/null \
                && echo "[INFO] Removed empty directory: {params.raw_dir}" >> "{log}" \
                || echo "[INFO] Directory not empty, keeping: {params.raw_dir}" >> "{log}"
        fi

        echo "[INFO] Cleanup completed for {wildcards.sample_id}" >> "{log}"
        echo "[INFO] Deletion completed for {wildcards.sample_id}. See {log}" > "{output.log}"
        """
