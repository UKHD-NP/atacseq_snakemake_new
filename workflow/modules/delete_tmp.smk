def _final_bam(wildcards, ext=""):
    """Return the path to the last permanent BAM: shifted (narrow peaks) or filtered."""
    if is_enabled("call_peaks") and CALL_PEAKS_PEAK_TYPE == "narrow" and is_enabled("shift_bam", default=True):
        return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.shifted.bam{ext}")
    return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.filtered.bam{ext}")


def _shifted_bam_is_final(wildcards=None):
    """True when _final_bam actually resolves to shifted.bam (narrow peaks + shift_bam on)."""
    return is_enabled("call_peaks") and CALL_PEAKS_PEAK_TYPE == "narrow" and is_enabled("shift_bam", default=True)


DELETE_SHIFTED_BAM_ON = as_bool(config.get("shift_bam", {}).get("delete_shifted_bam", False), default=False)


def _shifted_bam_consumer_outputs(wildcards):
    """Outputs that must exist before shifted.bam can be safely deleted - every rule
    that reads shifted.bam directly (shifted_bam_to_bigwig, NFR, ATACseqQC). Only
    populated when delete_shifted_bam is actually on, to avoid adding unnecessary
    DAG edges when this cleanup step is disabled (the default)."""
    if not (DELETE_SHIFTED_BAM_ON and _shifted_bam_is_final()):
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
        bam = lambda wildcards: _final_bam(wildcards),
        bai = lambda wildcards: _final_bam(wildcards, ".bai"),
        fastqc = lambda wildcards: [
            os.path.join(wildcards.outdir, "trim", f"{wildcards.sample_id}_trimmed_1_fastqc.zip"),
            os.path.join(wildcards.outdir, "trim", f"{wildcards.sample_id}_trimmed_2_fastqc.zip"),
        ] if is_enabled("trimming") and TRIM_TOOL == "trim_galore" else [],
        shifted_bam_consumers = lambda wildcards: _shifted_bam_consumer_outputs(wildcards)
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
        delete_trimming = lambda wildcards: str(
            is_enabled("trimming") and
            as_bool(config.get("trimming", {}).get("delete_trimming", True))
        ).lower(),
        delete_shifted_bam = lambda wildcards: str(DELETE_SHIFTED_BAM_ON and _shifted_bam_is_final()).lower()
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

        # Remove trimmed FASTQ files (conditional)
        if [ "{params.delete_trimming}" = "true" ]; then
            for f in "{params.fq1}" "{params.fq2}" "{params.fq1_fail}" "{params.fq2_fail}"; do
                if [ -f "$f" ]; then
                    rm -f "$f" && echo "[INFO] Removed $f" >> "{log}"
                fi
            done
        else
            echo "[INFO] delete_trimming=false, skipping trimmed FASTQ deletion." >> "{log}"
        fi

        # Remove merged raw FASTQs
        for f in "{params.raw_fq1}" "{params.raw_fq2}"; do
            if [ -e "$f" ] || [ -L "$f" ]; then
                rm -f "$f" && echo "[INFO] Removed $f" >> "{log}"
            fi
        done

        # Remove shifted.bam (conditional; only when it's actually the final BAM here -
        # narrow peaks + shift_bam enabled - and shift_bam.delete_shifted_bam=true).
        # All consumers (shifted.bigWig, NFR, ATACseqQC) are guaranteed done via
        # {input.shifted_bam_consumers}, so it's safe to delete at this point; it can
        # always be regenerated from filtered.bam via alignmentSieve if needed again.
        if [ "{params.delete_shifted_bam}" = "true" ]; then
            for f in "{input.bam}" "{input.bai}"; do
                if [ -f "$f" ]; then
                    rm -f "$f" && echo "[INFO] Removed $f (shift_bam.delete_shifted_bam=true)" >> "{log}"
                fi
            done
        else
            echo "[INFO] shift_bam.delete_shifted_bam=false (or not applicable), keeping shifted BAM." >> "{log}"
        fi

        # Remove raw_merged directory if empty
        if [ -d "{params.raw_dir}" ]; then
            if rmdir "{params.raw_dir}" 2>/dev/null; then
                echo "[INFO] Removed empty directory: {params.raw_dir}" >> "{log}"
            else
                echo "[INFO] Directory not empty, keeping: {params.raw_dir}" >> "{log}"
            fi
        fi

        echo "[INFO] Cleanup completed for {wildcards.sample_id}" >> "{log}"
        echo "[INFO] Deletion completed for {wildcards.sample_id}. See {log}" > "{output.log}"
        """
