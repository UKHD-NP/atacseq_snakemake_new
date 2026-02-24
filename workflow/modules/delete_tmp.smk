rule delete_tmp:
    # Clean up temporary files after analysis is complete
    input:
        done = os.path.join("{outdir}", "multiqc", "{sample_id}.multiqc.html")
    output:
        log = os.path.join("{outdir}", "logs", "{sample_id}.deletion.log")
    params:
        bam_dir = os.path.join("{outdir}", "bam"),
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

        # Start logging
        echo "[INFO] Starting cleanup for {wildcards.sample_id}" > {log}

        # Delete FASTQ files only if trimming is enabled
        if [ "{params.delete_trimming}" = "true" ]; then
            echo "[INFO] Trimming is enabled (delete_trimming={params.delete_trimming})." >> {log}

            # Remove trimmed FASTQ files (same for both fastp and trim_galore)
            if [ -f "{params.fq1}" ]; then
                rm -f "{params.fq1}" && echo "[INFO] Removed {params.fq1}" >> {log}
            fi
            if [ -f "{params.fq2}" ]; then
                rm -f "{params.fq2}" && echo "[INFO] Removed {params.fq2}" >> {log}
            fi

            # Remove failed FASTQ files (only created by fastp)
            if [ -f "{params.fq1_fail}" ]; then
                rm -f "{params.fq1_fail}" && echo "[INFO] Removed {params.fq1_fail}" >> {log}
            fi
            if [ -f "{params.fq2_fail}" ]; then
                rm -f "{params.fq2_fail}" && echo "[INFO] Removed {params.fq2_fail}" >> {log}
            fi

            echo "[INFO] Processed trimmed FASTQ files." >> {log}
        else
            echo "[INFO] Trimming is disabled (delete_trimming={params.delete_trimming})." >> {log}
            echo "[INFO] Skipping trimmed FASTQ deletion." >> {log}
        fi

        # Remove merged raw FASTQs for this sample
        if [ -e "{params.raw_fq1}" ] || [ -L "{params.raw_fq1}" ]; then
            rm -f "{params.raw_fq1}" && echo "[INFO] Removed {params.raw_fq1}" >> {log}
        fi
        if [ -e "{params.raw_fq2}" ] || [ -L "{params.raw_fq2}" ]; then
            rm -f "{params.raw_fq2}" && echo "[INFO] Removed {params.raw_fq2}" >> {log}
        fi

        # Remove raw_merged directory only when it becomes empty
        if [ -d "{params.raw_dir}" ]; then
            if rmdir "{params.raw_dir}" 2>/dev/null; then
                echo "[INFO] Removed empty raw_merged directory: {params.raw_dir}" >> {log}
            else
                echo "[INFO] raw_merged is not empty yet, keeping directory: {params.raw_dir}" >> {log}
            fi
        fi

        # Create final deletion log file
        echo "[INFO] Deletion completed for {wildcards.sample_id}" > {output.log}
        echo "[INFO] See detailed log at {log}" >> {output.log}
        echo "[INFO] Cleanup completed successfully." >> {log}
        """
