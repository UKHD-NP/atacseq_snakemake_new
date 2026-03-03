def get_input_multiqc(wildcards):
    """Collect input files for MultiQC based on enabled modules."""
    sample_id = wildcards.sample_id
    outdir = get_outdir(sample_id)

    def _path(*parts):
        return os.path.join(outdir, *parts)

    targets = []
    call_peaks_on = is_enabled("call_peaks")
    call_peaks_qc_on = call_peaks_on and as_bool(
        config.get("call_peaks", {}).get("macs3_peak_qc_plot", True),
        default=True,
    )
    feature_counts_on = call_peaks_on and is_enabled("feature_counts")
    deeptools_on = is_enabled("deeptools")
    annotation_on = is_enabled("annotate_peaks") and call_peaks_on
    ataqv_on = is_enabled("ataqv") and call_peaks_on

    # Conditionally add trimming quality control reports
    if is_enabled("trimming"):
        trimming_tool = config.get('trimming', {}).get('tool', 'fastp')

        # Add raw FastQC reports (before trimming) - for both tools
        targets.extend([
            _path("fastqc_raw", f"{sample_id}_raw_1_fastqc.zip"),
            _path("fastqc_raw", f"{sample_id}_raw_2_fastqc.zip"),
        ])

        if trimming_tool == "trim_galore":
            # Add trim_galore trimming reports
            targets.extend([
                _path("trim", f"{sample_id}_1.fastq.gz_trimming_report.txt"),
                _path("trim", f"{sample_id}_2.fastq.gz_trimming_report.txt"),
            ])
            # Add trim_galore FastQC reports on trimmed reads (ZIP files contain fastqc_data.txt that MultiQC parses)
            targets.extend([
                _path("trim", f"{sample_id}_trimmed_1_fastqc.zip"),
                _path("trim", f"{sample_id}_trimmed_2_fastqc.zip"),
            ])
        else:  # fastp
            # Add FastP quality control JSON
            targets.append(_path("trim", f"{sample_id}.fastp.json"))

    # Add Picard MarkDuplicates metrics if enabled
    if is_enabled("markduplicates"):
        targets.append(_path("bam", f"{sample_id}.markdup.sorted.MarkDuplicates.metrics.txt"))
    
    # Add pre-filter BAM stats outputs
    targets.extend(
        _path("bam", f"{sample_id}.pre_filter.bam.{ext}")
        for ext in ("stats", "flagstat", "idxstats")
    )

    # Add post-filter BAM stats outputs
    targets.extend(
        _path("bam", f"{sample_id}.filtered.bam.{ext}")
        for ext in ("stats", "flagstat", "idxstats")
    )

    # Add Picard CollectMultipleMetrics
    targets.extend([
        _path("bam", f"{sample_id}.CollectMultipleMetrics.alignment_summary_metrics"),
        _path("bam", f"{sample_id}.CollectMultipleMetrics.base_distribution_by_cycle_metrics"),
        _path("bam", f"{sample_id}.CollectMultipleMetrics.insert_size_metrics"),
        _path("bam", f"{sample_id}.CollectMultipleMetrics.quality_by_cycle_metrics"),
        _path("bam", f"{sample_id}.CollectMultipleMetrics.quality_distribution_metrics"),
    ])


    # Add peak-calling summaries.
    if call_peaks_on:
        targets.append(_path("peaks", f"{sample_id}_peaks.FRiP_mqc.tsv"))
        targets.append(_path("peaks", f"{sample_id}_peaks.count_mqc.tsv"))
    if call_peaks_qc_on:
        targets.append(_path("peaks", f"{sample_id}.macs_peakqc.summary.txt"))

    # Add annotation summaries.
    if annotation_on:
        targets.append(_path("annotation", f"{sample_id}.macs_annotatePeaks.summary.txt"))

    # Add featureCounts quantification summary.
    if feature_counts_on:
        targets.append(_path("featurecounts", f"{sample_id}.readCountInPeaks.txt.summary"))

    # Add deepTools fingerprint QC metrics.
    if deeptools_on:
        targets.append(_path("deeptools", f"{sample_id}.plotFingerprint.qcmetrics.txt"))
        targets.append(_path("deeptools", f"{sample_id}.plotFingerprint.raw_counts.txt"))
        targets.append(_path("deeptools", f"{sample_id}.scale_regions.plotProfile.tab"))
        targets.append(_path("deeptools", f"{sample_id}.fragment_size.qcmetrics.txt"))
        targets.append(_path("deeptools", f"{sample_id}.fragment_size.raw_lengths.txt"))

    # Add ataqv JSON metrics for ATAC QC.
    if ataqv_on:
        targets.append(_path("ataqv", f"{sample_id}.ataqv.json"))

    return list(dict.fromkeys(targets))

rule multiqc:
    # Aggregate quality control reports with MultiQC
    input:
        reports = get_input_multiqc
    output:
        os.path.join("{outdir}", "multiqc", "{sample_id}.multiqc.html")
    params:
        outdir = lambda wildcards: os.path.join(wildcards.outdir, "multiqc"),
        # Add extra parameters for optimized MultiQC performance
        extra_params = "--no-megaqc-upload --interactive",
        # Path to MultiQC config file
        config_file = os.path.join(workflow.basedir, "scripts", "multiqc_config.yml")
    conda:
        os.path.join(workflow.basedir, "envs", "multiqc.yml")
    message:
        "{wildcards.sample_id}: Running MultiQC"
    threads: 1
    resources:
        mem_mb = 2048
    log:
        os.path.join("{outdir}", "logs", "multiqc", "{sample_id}.multiqc.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.multiqc.benchmark.txt")
    shell:
        """
        # Set max file size limit for Python processes
        ulimit -n 2048 || true

        # Set Python environment variables to optimize memory usage
        export PYTHONHASHSEED=0
        export PYTHONIOENCODING=utf8

        # Create MultiQC output directory
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        # Run MultiQC
        multiqc {input.reports} \
            --outdir {params.outdir} \
            --filename {wildcards.sample_id}.multiqc.html \
            --config {params.config_file} \
            --force \
            {params.extra_params} \
            > {log} 2>&1 || {{ echo "[ERROR] MultiQC failed." >> {log}; exit 1; }}
        """
