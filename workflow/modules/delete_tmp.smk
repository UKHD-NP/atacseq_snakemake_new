def _final_bam(wildcards, ext=""):
    """Return the path to the last permanent BAM: shifted (narrow peaks) or filtered."""
    if is_enabled("call_peaks") and CALL_PEAKS_PEAK_TYPE == "narrow":
        return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.shifted.bam{ext}")
    return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.filtered.bam{ext}")


rule delete_tmp:
    # Clean up temporary files
    input:
        bam = lambda wildcards: _final_bam(wildcards),
        bai = lambda wildcards: _final_bam(wildcards, ".bai"),
        fastqc = lambda wildcards: [
            os.path.join(wildcards.outdir, "trim", f"{wildcards.sample_id}_trimmed_1_fastqc.zip"),
            os.path.join(wildcards.outdir, "trim", f"{wildcards.sample_id}_trimmed_2_fastqc.zip"),
        ] if is_enabled("trimming") and TRIM_TOOL == "trim_galore" else []
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
        ).lower()
    message:
        "{wildcards.sample_id}: Cleaning up temporary files"
    log:
        os.path.join("{outdir}", "logs", "cleanup", "{sample_id}.cleanup.log")
    shell:
        """
        mkdir -p $(dirname {log})
        echo "[INFO] Starting cleanup for {wildcards.sample_id}" > {log}

        # Remove trimmed FASTQ files (conditional)
        if [ "{params.delete_trimming}" = "true" ]; then
            for f in "{params.fq1}" "{params.fq2}" "{params.fq1_fail}" "{params.fq2_fail}"; do
                if [ -f "$f" ]; then
                    rm -f "$f" && echo "[INFO] Removed $f" >> {log}
                fi
            done
        else
            echo "[INFO] delete_trimming=false, skipping trimmed FASTQ deletion." >> {log}
        fi

        # Remove merged raw FASTQs
        for f in "{params.raw_fq1}" "{params.raw_fq2}"; do
            if [ -e "$f" ] || [ -L "$f" ]; then
                rm -f "$f" && echo "[INFO] Removed $f" >> {log}
            fi
        done

        # Remove raw_merged directory if empty
        if [ -d "{params.raw_dir}" ]; then
            if rmdir "{params.raw_dir}" 2>/dev/null; then
                echo "[INFO] Removed empty directory: {params.raw_dir}" >> {log}
            else
                echo "[INFO] Directory not empty, keeping: {params.raw_dir}" >> {log}
            fi
        fi

        echo "[INFO] Cleanup completed for {wildcards.sample_id}" >> {log}
        echo "[INFO] Deletion completed for {wildcards.sample_id}. See {log}" > {output.log}
        """
