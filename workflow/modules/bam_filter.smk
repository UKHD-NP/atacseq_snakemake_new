def get_bam_for_filter(wildcards):
    """BAM source for filtering depends on markduplicates setting."""
    if is_enabled("markduplicates"):
        return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.markdup.sorted.bam")
    return os.path.join(wildcards.outdir, "bam", f"{wildcards.sample_id}.bam")


bam_filter_cfg = config.get("bam_filter", {})
default_bam_filter_params = "-F 0x004 -F 0x0008 -f 0x001 -F 0x0100 -q 30"

BAM_FILTER_PREDUP_PARAMS = bam_filter_cfg.get("params", default_bam_filter_params)
BAM_FILTER_MITO_NAME = str(config.get("ref", {}).get("mito_name", "chrM")).strip() or "chrM"


rule bam_filter:
    # Filters the pre-dedup BAM (MAPQ, proper pair, NM<=4, no soft-clip,
    # include_regions, orphan cleanup) to produce a "quality-filtered,
    # duplicates-retained" population, computes ENCODE library complexity QC
    # (NRF, PBC1, PBC2) on that exact population, THEN removes duplicates as
    # the final step to produce filtered.bam.
    input:
        bam = get_bam_for_filter,
        include_regions = config["ref"]["include_regions"],
        pre_filter_bam_stats = os.path.join("{outdir}", "bam", "{sample_id}.pre_filter.bam.stats")
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.bai"),
        pbc_qc = os.path.join("{outdir}", "library_complexity", "{sample_id}.pbc_qc.tsv")
    params:
        predup_params = BAM_FILTER_PREDUP_PARAMS,
        dedup_flag = "-F 0x0400",
        mito_name = BAM_FILTER_MITO_NAME,
        bamtools_script = os.path.join(workflow.basedir, "scripts", "bamtools_filter_pe.json"),
        bampe_rm_orphan_script = os.path.join(workflow.basedir, "scripts", "bampe_rm_orphan.py"),
        tmp_dir = os.path.join("{outdir}", "bam", "tmp"),
        keep_input_bam = "true" if as_bool(bam_filter_cfg.get("keep_input_bam", False), default=False) else "false",
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "{wildcards.sample_id}: Filtering BAM + library complexity QC (NRF/PBC1/PBC2)"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 57344 + (attempt - 1) * 16384,
        runtime = lambda wildcards, attempt: attempt * 1200
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.bam_filter.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bam_filter.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.bam}")"
        mkdir -p "$(dirname "{output.pbc_qc}")"
        mkdir -p "$(dirname "{log}")"
        mkdir -p "{params.tmp_dir}"

        TMP_PREDUP="{params.tmp_dir}/{wildcards.sample_id}.predup.tmp.bam"
        TMP_NAME_SORTED="{params.tmp_dir}/{wildcards.sample_id}.predup.name_sorted.tmp.bam"
        TMP_CLEANED="{params.tmp_dir}/{wildcards.sample_id}.predup.cleaned.name_sorted.tmp.bam"
        TMP_BEDPE="{params.tmp_dir}/{wildcards.sample_id}.pbc_bedpe.tmp.tsv"

        mem_per_thread=$(( {resources.mem_mb} * 7 / 10 / ({threads} + 1) / 1024 ))
        mem_per_thread=${{mem_per_thread:-1}}G
        [ "${{mem_per_thread%%G}}" -lt 1 ] 2>/dev/null && mem_per_thread=1G

        set -euo pipefail

        # Step 1: quality filter (MAPQ, proper pair, NM<=4, no soft-clip,
        # include_regions) - duplicate-exclusion flag deliberately NOT applied here.
        if command -v bamtools >/dev/null 2>&1 && [ -f "{params.bamtools_script}" ]; then
            {{
                samtools view \
                    {params.predup_params} \
                    -L "{input.include_regions}" \
                    -b "{input.bam}" \
                    | bamtools filter \
                        -out "$TMP_PREDUP" \
                        -script "{params.bamtools_script}"
            }} > "{log}" 2>&1 || {{
                echo "[ERROR] BAM filtering (samtools+bamtools) failed." >> "{log}"
                exit 1
            }}
        else
            echo "[WARNING] bamtools not available or script not found, running samtools-only filter." > "{log}"
            samtools view \
                {params.predup_params} \
                -L "{input.include_regions}" \
                -b "{input.bam}" \
                > "$TMP_PREDUP" 2>> "{log}" || {{
                echo "[ERROR] BAM filtering (samtools-only) failed." >> "{log}"
                exit 1
            }}
        fi

        if [ ! -s "$TMP_PREDUP" ]; then
            echo "[ERROR] Intermediate BAM is missing or empty: $TMP_PREDUP" >> "{log}"
            exit 1
        fi

        # Step 2: name-sort, then orphan/pair cleanup (optional; requires pysam).
        samtools sort \
            -n \
            -m $mem_per_thread \
            -T "{params.tmp_dir}/{wildcards.sample_id}.name_sort_tmp" \
            -@ {threads} \
            -o "$TMP_NAME_SORTED" \
            "$TMP_PREDUP" >> "{log}" 2>&1 || {{
            echo "[ERROR] Name sort before orphan removal failed." >> "{log}"
            exit 1
        }}

        if [ -f "{params.bampe_rm_orphan_script}" ] && command -v python3 >/dev/null 2>&1 && python3 -c "import pysam" >/dev/null 2>&1; then
            python3 "{params.bampe_rm_orphan_script}" \
                "$TMP_NAME_SORTED" \
                "$TMP_CLEANED" \
                --only_fr_pairs >> "{log}" 2>&1 || {{
                echo "[ERROR] bampe_rm_orphan.py failed." >> "{log}"
                exit 1
            }}
            PREDUP_SORTED="$TMP_CLEANED"
        else
            echo "[WARNING] Skip orphan cleanup (missing script/python3/pysam)." >> "{log}"
            PREDUP_SORTED="$TMP_NAME_SORTED"
        fi

        if [ ! -s "$PREDUP_SORTED" ]; then
            echo "[ERROR] Pre-dedup name-sorted BAM is missing or empty: $PREDUP_SORTED" >> "{log}"
            exit 1
        fi

        # Step 3: NRF, PBC1, PBC2 (ENCODE library complexity QC) - computed here,
        # before duplicates are removed, since the metric identifies duplicates
        # from genomic position (bedtools bamtobed -bedpe), not the dup flag.
        # Split into 2 sub-steps (BEDPE position file, then count) - no separate
        # script files, both stay inline in this rule.
        export LC_ALL=C

        # 3a: BAM -> BEDPE position-key TSV (mito excluded), written to disk.
        bedtools bamtobed -bedpe -i "$PREDUP_SORTED" 2>> "{log}" \
            | awk -v mito="{params.mito_name}" 'BEGIN{{OFS="\\t"}} $1 != mito {{print $1,$2,$4,$6,$9,$10}}' \
            > "$TMP_BEDPE" 2>> "{log}" || {{
            echo "[ERROR] bedtools bamtobed -bedpe / mito filter failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "$TMP_BEDPE" ]; then
            echo "[ERROR] BEDPE position file is missing or empty: $TMP_BEDPE" >> "{log}"
            exit 1
        fi

        # 3b: count duplicate positions in that TSV -> NRF, PBC1, PBC2.
        sort "$TMP_BEDPE" \
            | uniq -c \
            | awk -v sample="{wildcards.sample_id}" '
                BEGIN {{ mt=0; m0=0; m1=0; m2=0 }}
                ($1==1) {{ m1++ }}
                ($1==2) {{ m2++ }}
                {{ m0++; mt+=$1 }}
                END {{
                    nrf  = (mt>0) ? m0/mt : 0
                    pbc1 = (m0>0) ? m1/m0 : 0
                    pbc2 = (m2>0) ? m1/m2 : 0
                    print "Sample\\tTotalReadPairs\\tDistinctReadPairs\\tOneReadPair\\tTwoReadPairs\\tNRF\\tPBC1\\tPBC2"
                    printf "%s\\t%d\\t%d\\t%d\\t%d\\t%.4f\\t%.4f\\t%.4f\\n", sample, mt, m0, m1, m2, nrf, pbc1, pbc2
                }}
            ' > "{output.pbc_qc}" 2>> "{log}" || {{
            echo "[ERROR] library complexity (NRF/PBC1/PBC2) calculation failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.pbc_qc}" ]; then
            echo "[ERROR] library complexity output is missing or empty: {output.pbc_qc}" >> "{log}"
            exit 1
        fi
        echo "[INFO] $(tail -n1 "{output.pbc_qc}")" >> "{log}"

        # Step 4: remove duplicates (fixed, always-on filter) + coordinate-sort + index -> filtered.bam
        samtools view {params.dedup_flag} -b "$PREDUP_SORTED" \
            | samtools sort \
                --write-index \
                -m $mem_per_thread \
                -T "{params.tmp_dir}/{wildcards.sample_id}.filtered" \
                -@ {threads} \
                -o "{output.bam}##idx##{output.bai}" \
                - >> "{log}" 2>&1 || {{
            echo "[ERROR] Final duplicate filter/sort/index failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.bam}" ] || [ ! -s "{output.bai}" ]; then
            echo "[ERROR] Filtered BAM or BAI is missing." >> "{log}"
            exit 1
        fi

        echo "[INFO] Removing intermediate BAM files." >> "{log}"
        rm -f "$TMP_PREDUP" "$TMP_NAME_SORTED" "$TMP_CLEANED" "$TMP_BEDPE"
        if [ "{params.keep_input_bam}" = "true" ]; then
            echo "[INFO] keep_input_bam=true, preserving source BAM: {input.bam}" >> "{log}"
        else
            echo "[INFO] Removing source BAM to save space: {input.bam}" >> "{log}"
            rm -f "{input.bam}" "{input.bam}.bai"
        fi
        rm -rf "{params.tmp_dir}"
        """
