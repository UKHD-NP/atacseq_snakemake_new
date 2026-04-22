rule samtools_faidx:
    # Generate FASTA index
    input:
        config['ref']['fasta']
    output:
        config['ref']['fasta'] + ".fai"
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "Generating FASTA index"
    threads: 1
    resources:
        mem_mb = 1024
    shell:
        """
        samtools faidx {input}
        """


rule gtf2bed:
    # Generate BED file from GTF
    input:
        config['ref']['gtf']
    output:
        config['ref']['bed']
    params:
        gtf2bed_script = os.path.join(workflow.basedir, "scripts", "gtf2bed")
    conda:
        os.path.join(workflow.basedir, "envs", "gtf2bed.yml")
    message:
        "Converting GTF to BED format"
    threads: 1
    log:
        os.path.join(ref_dir, "gtf2bed", f"{config['ref']['assembly']}.log")
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})

        if [ ! -f "{params.gtf2bed_script}" ]; then
            echo "[ERROR] gtf2bed script not found: {params.gtf2bed_script}" > {log}
            exit 1
        fi

        echo "[INFO] Using gtf2bed script: {params.gtf2bed_script}" > {log}
        perl "{params.gtf2bed_script}" "{input}" > "{output}" 2>> "{log}"

        # Check if output was created successfully
        if [ ! -s "{output}" ]; then
            echo "[ERROR] Failed to create BED file: {output}" >> {log}
            exit 1
        fi

        echo "[INFO] BED file created successfully: {output}" >> {log}
        """


rule get_chromsizes:
    # Convert FASTA FAI index to UCSC-style chrom.sizes.
    input:
        fai=config["ref"]["fasta"] + ".fai"
    output:
        config["ref"]["chromsizes"]
    message:
        "Generating chromosome sizes file"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join(ref_dir, "chromsizes", f"{config['ref']['assembly']}.log")
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        cut -f 1,2 "{input.fai}" > "{output}" 2> "{log}"
        if [ ! -s "{output}" ]; then
            echo "[ERROR] chromsizes output is empty: {output}" >> "{log}"
            exit 1
        fi
        """


_BAM_FILTER_CFG = config.get("bam_filter", {}) if isinstance(config.get("bam_filter", {}), dict) else {}
_REF_CFG = config.get("ref", {}) if isinstance(config.get("ref", {}), dict) else {}

# Pattern for ALL canonical chromosomes (used by bam_filter include_regions).
# Default: human hg19/hg38 (chr1-22 + X + Y + M).
# Override via bam_filter.canonical_chroms_pattern for other genomes:
#   mouse/rat/zebrafish : "^chr([0-9]+|X|Y|M)$"
#   ENSEMBL (no chr)    : "^([0-9]+|X|Y|MT)$"
_DEFAULT_CANONICAL_PATTERN = r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$"
_CANONICAL_PATTERN = (
    str(_BAM_FILTER_CFG.get("canonical_chroms_pattern", "")).strip()
    or _DEFAULT_CANONICAL_PATTERN
)

# Pattern for TRUE autosomes only (used by ataqv --autosomal-reference-file).
# Default: human hg19/hg38 autosomes (chr1-22).
# Override via ref.autosome_pattern for other genomes:
#   mouse/rat/zebrafish : "^chr([0-9]+)$"
#   ENSEMBL (no chr)    : "^([0-9]+)$"
#   C. elegans ce11     : "^(I|II|III|IV|V)$"
_DEFAULT_AUTOSOME_PATTERN = r"^chr([1-9]|1[0-9]|2[0-2])$"
_AUTOSOME_PATTERN = (
    str(_REF_CFG.get("autosome_pattern", "")).strip()
    or _DEFAULT_AUTOSOME_PATTERN
)


rule get_canonical_chroms:
    # All canonical chromosomes (autosomes + sex + MT) for bam_filter include_regions.
    # If apply_canonical_chromosomes=true, applies canonical_chroms_pattern (default: human chr1-22/X/Y/M).
    # If apply_canonical_chromosomes=false, outputs all chromosomes from the FASTA (any reference).
    input:
        fai = config["ref"]["fasta"] + ".fai"
    output:
        config["ref"]["canonical_chroms"]
    params:
        apply_canonical = "true" if as_bool(_BAM_FILTER_CFG.get("apply_canonical_chromosomes", False), default=False) else "false",
        pattern = _CANONICAL_PATTERN
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "Extracting canonical chromosome list"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join(ref_dir, "canonical_chroms", f"{config['ref']['assembly']}.log")
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})

        if [ "{params.apply_canonical}" = "true" ]; then
            awk 'BEGIN{{FS=OFS="\\t"}} $1 ~ /{params.pattern}/ {{print $1}}' "{input.fai}" > "{output}" 2> "{log}"
        else
            awk 'BEGIN{{FS=OFS="\\t"}} {{print $1}}' "{input.fai}" > "{output}" 2> "{log}"
        fi

        if [ ! -s "{output}" ]; then
            echo "[ERROR] canonical chromosome output is empty: {output}" >> "{log}"
            exit 1
        fi
        """


rule get_autosomes:
    # True autosomes only (no sex chromosomes, no MT) for ataqv --autosomal-reference-file.
    # Derived directly from FASTA index via ref.autosome_pattern, independent of bam_filter canonical settings.
    input:
        fai = config["ref"]["fasta"] + ".fai"
    output:
        config["ref"]["autosomes"]
    params:
        pattern = _AUTOSOME_PATTERN
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "Extracting autosome chromosome list"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join(ref_dir, "autosomes", f"{config['ref']['assembly']}.log")
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})

        awk -v pattern='{params.pattern}' 'BEGIN{{FS=OFS="\\t"}} $1 ~ pattern {{print $1}}' \
            "{input.fai}" > "{output}" 2> "{log}"

        if [ ! -s "{output}" ]; then
            echo "[ERROR] autosome list is empty: {output}" >> "{log}"
            echo "[HINT] Check ref.autosome_pattern in config." >> "{log}"
            exit 1
        fi
        """


rule tss_extract:
    # Extract 1bp TSS intervals from BED annotation.
    input:
        bed=config["ref"]["bed"]
    output:
        config["ref"]["tss"]
    message:
        "Extracting TSS BED intervals"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join(ref_dir, "tss", f"{config['ref']['assembly']}.log")
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        awk 'BEGIN{{FS=OFS="\\t"}} NF>=6 {{
            start=$2; end=$3;
            if ($6 == "+") {{ end = start + 1; }}
            else if ($6 == "-") {{ start = end - 1; }}
            if (start < 0) {{ start = 0; }}
            print $1, start, end, $4, $5, $6;
        }}' "{input.bed}" > "{output}" 2> "{log}"

        if [ ! -s "{output}" ]; then
            echo "[ERROR] Failed to create TSS BED: {output}" >> "{log}"
            exit 1
        fi
        """


rule genome_blacklist_regions:
    # Build mappable genome intervals with optional canonical chromosome, blacklist, and mito filtering.
    input:
        sizes=config["ref"]["chromsizes"],
        canonical=config["ref"]["canonical_chroms"]
    output:
        config["ref"]["include_regions"]
    params:
        blacklist=config["ref"].get("blacklist", ""),
        apply_blacklist="true" if as_bool(config.get("bam_filter", {}).get("apply_blacklist", True), default=True) else "false",
        apply_canonical="true" if as_bool(config.get("bam_filter", {}).get("apply_canonical_chromosomes", False), default=False) else "false",
        mito_name=config["ref"].get("mito_name", "chrM"),
        keep_mito="true" if as_bool(config["ref"].get("keep_mito", False), default=False) else "false",
        tmp=lambda wildcards, output: f"{output}.tmp"
    conda:
        os.path.join(workflow.basedir, "envs", "bedtools.yml")
    message:
        "Generating include regions with optional canonical chromosome, blacklist, and mito filtering"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join(ref_dir, "include_regions", f"{config['ref']['assembly']}.log")
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})

        # Sort sizes file to match expected chromosome order
        LC_ALL=C sort -k1,1 "{input.sizes}" > "{params.tmp}.sizes"

        if [ "{params.apply_canonical}" = "true" ]; then
            awk 'NR==FNR {{ keep[$1]=1; next }} ($1 in keep)' "{input.canonical}" "{params.tmp}.sizes" > "{params.tmp}.canonical.sizes"
            mv "{params.tmp}.canonical.sizes" "{params.tmp}.sizes"
        fi

        if [ "{params.apply_blacklist}" = "true" ] && [ -n "{params.blacklist}" ] && [ -s "{params.blacklist}" ]; then
            LC_ALL=C sort -k1,1 -k2,2n "{params.blacklist}" > "{params.tmp}.blacklist"
            bedtools complement -i "{params.tmp}.blacklist" -g "{params.tmp}.sizes" > "{params.tmp}"
            rm -f "{params.tmp}.blacklist"
        else
            awk 'BEGIN{{OFS="\\t"}} {{print $1,0,$2}}' "{params.tmp}.sizes" > "{params.tmp}"
        fi
        rm -f "{params.tmp}.sizes"

        if [ "{params.keep_mito}" = "true" ] || [ -z "{params.mito_name}" ]; then
            mv "{params.tmp}" "{output}"
        else
            awk -v mito="{params.mito_name}" '$1 != mito' "{params.tmp}" > "{output}"
            rm -f "{params.tmp}"
        fi
        if [ ! -s "{output}" ]; then
            echo "[ERROR] include regions output is empty: {output}" >> "{log}"
            exit 1
        fi
        """
