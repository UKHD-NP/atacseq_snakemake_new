rule bedtools_genomecov:
    # Build scaled bedGraph from filtered BAM using mapped-read normalization.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam"),
        flagstat = os.path.join("{outdir}", "bam", "{sample_id}.filtered.bam.flagstat")
    output:
        bedgraph = os.path.join("{outdir}", "bigwig", "{sample_id}.bedGraph"),
        scale_factor=os.path.join("{outdir}", "bigwig", "{sample_id}.scale_factor.txt")
    conda:
        os.path.join(workflow.basedir, "envs", "bedtools.yml")
    message:
        "{wildcards.sample_id}: Generating scaled bedGraph - [Source: Filtered BAM, NOT SHIFTED]"
    threads: 2
    resources:
        mem_mb = 40960
    log:
        os.path.join("{outdir}", "logs", "bedtools", "{sample_id}.genomecov.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bedtools_genomecov.benchmark.txt")
    shell:
        """
        ulimit -n 65536
        
        mkdir -p "$(dirname "{output.bedgraph}")"
        mkdir -p "$(dirname "{log}")"

        TMP_BG="{output.bedgraph}.tmp.bg"
        MAPPED_READS=$(awk '/ mapped \\(/ && $0 !~ /primary mapped/ {{print $1; exit}}' "{input.flagstat}")
        if [ -z "$MAPPED_READS" ] || ! [[ "$MAPPED_READS" =~ ^[0-9]+$ ]] || [ "$MAPPED_READS" -le 0 ]; then
            echo "[ERROR] Could not compute scale factor from flagstat: {input.flagstat}" >> "{log}"
            exit 1
        fi
        SCALE_FACTOR=$(awk -v m="$MAPPED_READS" 'BEGIN{{printf "%.10f", 1000000/m}}')
        echo "$SCALE_FACTOR" > "{output.scale_factor}"

        bedtools genomecov \
            -ibam "{input.bam}" \
            -bg \
            -scale "$SCALE_FACTOR" \
            -pc \
            > "$TMP_BG" 2>> "{log}" || {{
            echo "[ERROR] bedtools genomecov failed." >> "{log}"
            exit 1
        }}

        sort -k1,1 -k2,2n --parallel={threads} -S 2G "$TMP_BG" > "{output.bedgraph}"

        rm -f "$TMP_BG"

        if [ ! -s "{output.bedgraph}" ]; then
            echo "[ERROR] bedGraph output is missing or empty: {output.bedgraph}" >> "{log}"
            exit 1
        fi
        """


rule ucsc_bedgraphtobigwig:
    # Convert bedGraph to bigWig.
    input:
        bedgraph=os.path.join("{outdir}", "bigwig", "{sample_id}.bedGraph"),
        chromsizes=config["ref"]["chromsizes"]
    output:
        bigwig=os.path.join("{outdir}", "bigwig", "{sample_id}.bigWig")
    conda:
        os.path.join(workflow.basedir, "envs", "ucsc_tools.yml")
    message:
        "{wildcards.sample_id}: Converting bedGraph to bigWig (Source: Filtered BAM, NOT SHIFTED)"
    threads: 1
    resources:
        mem_mb = 1024
    log:
        os.path.join("{outdir}", "logs", "ucsc", "{sample_id}.bedgraphtobigwig.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bedgraphtobigwig.benchmark.txt")
    shell:
        """
        mkdir -p "$(dirname "{output.bigwig}")"
        mkdir -p "$(dirname "{log}")"

        bedGraphToBigWig \
            "{input.bedgraph}" \
            "{input.chromsizes}" \
            "{output.bigwig}" > "{log}" 2>&1 || {{
            echo "[ERROR] bedGraphToBigWig failed." >> "{log}"
            exit 1
        }}

        if [ ! -s "{output.bigwig}" ]; then
            echo "[ERROR] bigWig output is missing or empty: {output.bigwig}" >> "{log}"
            exit 1
        fi
        """
