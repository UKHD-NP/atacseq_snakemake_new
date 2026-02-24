TRUE_VALUES = {"true", "yes", "1", "t", "y"}


def as_bool(value, default=False):
    """
    Normalize config values to boolean.
    
    Accepts various representations of 'True':
    - True (boolean)
    - 'true' (lowercase string)
    - 'True' (capitalized string)
    - 'TRUE' (uppercase string)
    - 1 (integer)
    Returns False for other values.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in TRUE_VALUES
    if isinstance(value, int):
        return value == 1
    return default


def is_enabled(module_name, default=False):
    """
    Check if a module is enabled in the config.

    - Section missing or empty → returns *default*.
    - Section present but no 'enabled' key → True.
    - Section present with 'enabled' key → as_bool(value).

    Examples:
        is_enabled("trimming")                   # off when section missing
        is_enabled("bam_filter", default=True)    # on unless explicitly disabled
    """
    module_cfg = config.get(module_name, {})
    if not module_cfg:
        return default
    if isinstance(module_cfg, dict):
        if "enabled" not in module_cfg:
            return True
        return as_bool(module_cfg["enabled"], default=default)
    return as_bool(module_cfg, default=default)


def cfg_str(cfg_dict, key, default):
    """Read string from a config dict with fallback."""
    value = cfg_dict.get(key, default)
    if value is None:
        return str(default)
    return str(value).strip()


def get_effective_genome_size(chromsizes_path, configured_value=""):
    """Use configured genome size if set; otherwise sum chromosome sizes."""
    configured = "" if configured_value is None else str(configured_value).strip()
    if configured:
        return configured

    total = 0
    with open(chromsizes_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 2:
                total += int(fields[1])

    if total <= 0:
        fatal(f"Failed to compute effective genome size from chromsizes: {chromsizes_path}")
    return str(total)


def get_sample_rows(sample_id):
    """Return all samplesheet rows matching a sample ID."""
    sample_rows = samplesheet[samplesheet['sample_id'] == sample_id]
    if sample_rows.empty:
        fatal(f"Sample ID '{sample_id}' not found in samplesheet")
    return sample_rows


def get_outdir(sample_id):
    """Get the output directory for a given sample ID."""
    sample_rows = get_sample_rows(sample_id)
    outdirs = sample_rows['outdir'].dropna().unique().tolist()
    if len(outdirs) != 1:
        fatal(
            f"Sample ID '{sample_id}' maps to multiple outdir values in samplesheet: {outdirs}. "
            "Please keep outdir consistent for duplicated sample_id rows."
        )
    return outdirs[0]


def get_raw_lane_fastqs(wildcards):
    """
    Return all raw FASTQ lane files for a sample ID.
    Supports multiple rows per sample_id in the samplesheet.
    """
    sample_rows = get_sample_rows(wildcards.sample_id)
    fq1s = sample_rows['fq1'].tolist()
    fq2s = sample_rows['fq2'].tolist()
    if len(fq1s) != len(fq2s):
        fatal(f"Sample ID '{wildcards.sample_id}' has unequal fq1/fq2 counts: {len(fq1s)} vs {len(fq2s)}")

    invalid_fq1 = [str(path) for path in fq1s if not str(path).strip().lower().endswith(".gz")]
    invalid_fq2 = [str(path) for path in fq2s if not str(path).strip().lower().endswith(".gz")]
    if invalid_fq1 or invalid_fq2:
        invalid_parts = []
        if invalid_fq1:
            invalid_parts.append(f"fq1={invalid_fq1}")
        if invalid_fq2:
            invalid_parts.append(f"fq2={invalid_fq2}")
        fatal(
            f"Sample ID '{wildcards.sample_id}' requires gzipped FASTQ inputs (*.gz); "
            f"found non-gz paths: {'; '.join(invalid_parts)}"
        )

    return fq1s, fq2s


def get_raw_lane_fq1(wildcards):
    return get_raw_lane_fastqs(wildcards)[0]


def get_raw_lane_fq2(wildcards):
    return get_raw_lane_fastqs(wildcards)[1]


# Merge raw FASTQ lanes for duplicated sample IDs.
# For samples with one lane, this is a simple copy-through via cat.
rule merge_raw_fastqs:
    input:
        fq1 = get_raw_lane_fq1,
        fq2 = get_raw_lane_fq2
    output:
        fq1 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_2.fastq.gz")
    log:
        os.path.join("{outdir}", "logs", "raw_merge", "{sample_id}.merge_raw_fastqs.log")
    message:
        "{wildcards.sample_id}: Merging raw FASTQ lanes"
    shell:
        """
        mkdir -p $(dirname {output.fq1})
        mkdir -p $(dirname {log})
        rm -f {output.fq1} {output.fq2}

        N_R1=$(echo {input.fq1} | wc -w)
        N_R2=$(echo {input.fq2} | wc -w)

        if [ "$N_R1" -eq 1 ] && [ "$N_R2" -eq 1 ]; then
            # Single-lane sample: use symlink to avoid duplicate storage
            ln -sf "$(readlink -f {input.fq1})" {output.fq1}
            ln -sf "$(readlink -f {input.fq2})" {output.fq2}
            echo "[INFO] Mode: symlink (single lane)" > {log}
        else
            # Multi-lane sample: concatenate lanes in listed order
            cat {input.fq1} > {output.fq1}
            cat {input.fq2} > {output.fq2}
            echo "[INFO] Mode: merge (multi-lane)" > {log}
        fi

        if [ ! -s "{output.fq1}" ] || [ ! -s "{output.fq2}" ]; then
            echo "[ERROR] Merged FASTQ output is empty." >> {log}
            exit 1
        fi

        echo "[INFO] Merged R1 inputs: {input.fq1}" >> {log}
        echo "[INFO] Merged R2 inputs: {input.fq2}" >> {log}
        """


def get_paired_fq(wildcards):
    """Get merged raw FASTQ paths for a sample."""
    outdir = get_outdir(wildcards.sample_id)
    return [
        os.path.join(outdir, "raw_merged", f"{wildcards.sample_id}_merged_{read}.fastq.gz")
        for read in ("1", "2")
    ]


def get_paired_trimmed_fq(wildcards):
    """Return trimmed FASTQ paths if trimming is enabled, otherwise raw FASTQs."""
    if is_enabled("trimming"):
        # Check if outdir is available in wildcards, otherwise resolve from samplesheet
        outdir = getattr(wildcards, 'outdir', get_outdir(wildcards.sample_id))

        # Both fastp and trim_galore output to the same "trim" directory
        trim_dir = os.path.join(outdir, "trim")
        return [
            os.path.join(trim_dir, f"{wildcards.sample_id}_trimmed_{read}.fastq.gz")
            for read in ("1", "2")
        ]
    else:
        return get_paired_fq(wildcards)
        


def get_target_files(sample_ids):
    """Determine target files for the workflow based on enabled modules."""
    targets = []

    multiqc_on = is_enabled("multiqc", default=True)
    bam_filter_on = is_enabled("bam_filter", default=True)
    align_stats_on = is_enabled("align_stats", default=True) and bam_filter_on
    markdup_on = is_enabled("markduplicates")
    call_peaks_on = is_enabled("call_peaks") and bam_filter_on
    annotate_peaks_on = call_peaks_on and is_enabled("annotate_peaks")
    feature_counts_on = call_peaks_on and is_enabled("feature_counts")
    deeptools_on = is_enabled("deeptools") and bam_filter_on
    ataqv_on = is_enabled("ataqv") and call_peaks_on
    
    # Process each sample
    for sample_id in sample_ids:
        outdir = get_outdir(sample_id)

        def _path(*parts):
            return os.path.join(outdir, *parts)
        
        # MultiQC report
        if multiqc_on:
            targets.append(_path("multiqc", f"{sample_id}.multiqc.html"))

        # Core BAM-derived stats outputs.
        if align_stats_on:
            targets.extend([
                _path("bam", f"{sample_id}.filtered.bam.stats"),
                _path("bam", f"{sample_id}.filtered.bam.flagstat"),
                _path("bam", f"{sample_id}.filtered.bam.idxstats"),
            ])

        # Duplicate marking outputs (optional).
        if markdup_on:
            targets.extend([
                _path("bam", f"{sample_id}.markdup.sorted.MarkDuplicates.metrics.txt"),
            ])

        # Filtered BAM outputs (default on).
        if bam_filter_on:
            targets.extend([
                _path("bam", f"{sample_id}.filtered.bam"),
                _path("bam", f"{sample_id}.filtered.bam.bai"),
            ])

        # Peak calling and annotation outputs.
        if call_peaks_on:
            targets.extend([
                _path("peaks", f"{sample_id}_peaks.peak"),
                _path("peaks", f"{sample_id}_peaks.xls"),
                _path("peaks", f"{sample_id}.FRiP.txt"),
                _path("peaks", f"{sample_id}.macs_peakqc.summary.txt"),
                _path("peaks", f"{sample_id}.macs_peakqc.plots.pdf"),
            ])

        # featureCounts outputs (requires called peaks).
        if feature_counts_on:
            targets.extend([
                _path("featurecounts", f"{sample_id}_peaks.saf"),
                _path("featurecounts", f"{sample_id}.readCountInPeaks.txt"),
                _path("featurecounts", f"{sample_id}.readCountInPeaks.txt.summary"),
            ])

        # Peak annotation outputs (optional, requires called peaks).
        if annotate_peaks_on:
            targets.append(_path("annotation", f"{sample_id}_peaks.annotatePeaks.txt"))

        # bigWig and deepTools outputs.
        if deeptools_on:
            # Non-shifted bigWig requires flagstat from align_stats.
            if align_stats_on:
                targets.append(_path("bigwig", f"{sample_id}.bigWig"))

            targets.extend([
                _path("bam", f"{sample_id}.shifted.bam"),
                _path("bam", f"{sample_id}.shifted.bam.bai"),
                _path("bigwig", f"{sample_id}.shifted.bigWig"),
                _path("deeptools", f"{sample_id}.scale_regions.computeMatrix.mat.gz"),
                _path("deeptools", f"{sample_id}.scale_regions.computeMatrix.vals.mat.tab"),
                _path("deeptools", f"{sample_id}.reference_point.computeMatrix.mat.gz"),
                _path("deeptools", f"{sample_id}.reference_point.computeMatrix.vals.mat.tab"),
                _path("deeptools", f"{sample_id}.scale_regions.plotProfile.pdf"),
                _path("deeptools", f"{sample_id}.scale_regions.plotProfile.tab"),
                _path("deeptools", f"{sample_id}.reference_point.plotHeatmap.pdf"),
                _path("deeptools", f"{sample_id}.reference_point.plotHeatmap.mat.tab"),
                _path("deeptools", f"{sample_id}.plotFingerprint.pdf"),
                _path("deeptools", f"{sample_id}.plotFingerprint.raw.txt"),
                _path("deeptools", f"{sample_id}.plotFingerprint.qcmetrics.txt"),
            ])

        # ataqv outputs.
        if ataqv_on:
            targets.append(_path("ataqv", f"{sample_id}.ataqv.json"))
            targets.append(_path("ataqv", f"{sample_id}.mkarv_html", "index.html"))

        # Include cleanup log only when MultiQC target is enabled
        # (delete_tmp currently depends on MultiQC output)
        if multiqc_on:
            targets.append(_path("logs", f"{sample_id}.deletion.log"))

    return list(dict.fromkeys(targets))
