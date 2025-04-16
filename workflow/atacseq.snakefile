##########################################
## snakeATAC: Snakefile for DNA mapping ##
##########################################

import os
#conda_prefix = str(os.environ["CONDA_PREFIX"])

import sys
#sys.path.insert(1, conda_prefix+"/lib/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+"/site-packages")

from typing import List
import pathlib
import re
import numpy
import pandas as pd
import math
from itertools import combinations
from itertools import chain

    
# Define general variables
genome_fasta = str(config["genome_fasta"])


### working directory
home_dir = os.path.join(config["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir



# Load the sample sheet from the config file
samplesheet = pd.read_csv(config["samplesheet"], sep=",")

# Create a new column that combines 'sample' and 'replicate'
samplesheet['sample_name'] = samplesheet['sample'] + "_" + samplesheet['replicate']

# Ensure that 'sample_name' column is properly recognized
sample_list = samplesheet['sample_name'].tolist()

# Generate FILENAMES (list of all FASTQ files)
FILENAMES = samplesheet['fq1'].tolist() + samplesheet['fq2'].dropna().tolist()

# Generate RUNNAMES (filenames without suffixes, derived from fq1 and fq2)
RUNNAMES = [
    re.sub(rf"{config['fastq_suffix']}$", "", os.path.basename(fq))
    for fq in FILENAMES
]

# Generate SAMPLENAMES (unique sample names derived from 'sample_name')
SAMPLENAMES = samplesheet['sample_name'].tolist()


### Other variables
bwa_version = "bwa-mem2"
bwa_idx = ".bwt.2bit.64"
PEAKCALLER = "MACS3"
PEAKSDIR = "04_"+PEAKCALLER+"_peaks/"
SUMMARYDIR = "06_Overall_quality_and_info/"


MAPQ = str(config["MAPQ_threshold"])
genome_id = config["genomic_annotations"]["genome_id"]
BLACKLIST = str(config["blacklist"])
GROUPNAMES = SAMPLENAMES
norm_bw_average = []


if (eval(str(config["remove_duplicates"])) == True):
    DUP = "dedup"
else:
    DUP = "mdup"


# Function to check if a sample is paired-end
def is_paired_end(sample_name):
    # Select the relevant row(s) for the given sample_name
    sample_units = samplesheet[samplesheet['sample_name'] == sample_name]
    
    # Check if 'fq2' is null (i.e., missing) for all the rows corresponding to this sample
    fq2_null = sample_units["fq2"].isnull()
    
    # Determine if this sample is paired-end
    paired = ~fq2_null
    
    # Check if all associated rows are either paired-end or single-end
    all_paired = paired.all()
    all_single = (~paired).all()
    
    # Assert that the sample is either completely paired-end or single-end
    assert (
        all_single or all_paired
    ), f"Invalid units for sample {sample_name}, must be all paired-end or all single-end."
    
    return all_paired

# Function to get paired FASTQ files for a given sample
def get_paired_fq(wildcards):
    #print(f"DataFrame shape: {samplesheet.shape}")
    #print(samplesheet.head())  # Print first few rows
    #print(f"Attempting to access index for SAMPLE={wildcards.SAMPLE}")
    #if samplesheet.empty:
    #  raise ValueError("Error: The DataFrame is empty. Check input file paths.")
    #if wildcards.SAMPLE not in samplesheet.index:
    #  raise KeyError(f"Sample {wildcards.SAMPLE} not found in DataFrame.")
    sample_data = samplesheet[samplesheet['sample_name'] == wildcards.SAMPLE].iloc[0]
    if is_paired_end(wildcards.SAMPLE):
        return [sample_data['fq1'], sample_data['fq2']]  # Return a list of file paths
    else:
        raise ValueError(f"Sample {wildcards.SAMPLE} is not paired-end.")

# Function to get trimmed paired FASTQ files for a given sample
def get_paired_trimmed_fq(wildcards):
    sample_data = samplesheet[samplesheet['sample_name'] == wildcards.sample].iloc[0]
    if is_paired_end(wildcards.sample):
        if config["trimming"]["enabled"]:
            return [
                f"results/trimmed/{wildcards.sample}/{wildcards.sample}_1.fastp.fastq.gz",
                f"results/trimmed/{wildcards.sample}/{wildcards.sample}_2.fastp.fastq.gz"
            ]
        else:
            return [sample_data['fq1'], sample_data['fq2']]
    else:
        raise ValueError(f"Sample {wildcards.sample} is not paired-end.")



# Chromosome remove chr_remove_pattern
if (len(config["remove_other_chromosomes_pattern"]) > 0):
    chr_remove_pattern = '^chrM|^M|'+config["remove_other_chromosomes_pattern"]
else:
    chr_remove_pattern = '^chrM|^M'


### Generation of global wildcard_constraints
# Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
            raise ValueError("constraint_to(): Expected a list, got str instead")
    return "({})".format("|".join(values))


rm_duplicates = str(config["bam_features"]["remove_duplicates"]).lower()
if (eval(str(rm_duplicates.capitalize())) == True):
    DUP="dedup"
else:
    DUP="mdup"
    
if (eval(str(config["run_fastq_qc"])) == True):
    multiqc_fastq = "05_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.html"
else:
    multiqc_fastq = []


if (eval(str(config["peak_calling"]["call_summits"])) == True):
    SUMMITS="--call-summits"
else:
    SUMMITS=""


wildcard_constraints:
    SAMPLE = constraint_to(SAMPLENAMES),
    GROUPS = constraint_to(GROUPNAMES)

    
# Correlations outputs
peaks_comparison = []
if (len(SAMPLENAMES) > 1):
    peaks_comparison.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed"))
    peaks_comparison.append(os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_", PEAKCALLER, ".npz"])))
    peaks_comparison.append(os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_table_", PEAKCALLER, ".tsv"])))

# ruleorder: fastQC_filtered_BAM > normalized_bigWig > raw_bigWig

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all funtions
rule AAA_initialization:
    input:
        filtBAM_sorted_woMT = expand(os.path.join("02_BAM", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])), sample=SAMPLENAMES),
        filtBAM_sorted_woMT_index = expand(os.path.join("02_BAM", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])), sample=SAMPLENAMES),
        multiqc_fastq = multiqc_fastq,
        multiQC_BAM_html = os.path.join(SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/multiQC_report_BAMs_", DUP, ".html"])),
        report_fragSize_ggplot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots_ggplot.version.pdf"),
        normalized_bigWig = expand(''.join(["03_Normalization/RPM_normalized/{sample}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized.bw"]), sample=SAMPLENAMES),
        narrowPeaks_peaks = expand(os.path.join(PEAKSDIR, "{sample}_mapq{mapq}_woMT_{dup}_qValue{qValue}_peaks.narrowPeak"), sample=SAMPLENAMES, mapq=MAPQ, dup=DUP, qValue=str(config["peak_calling"]["qValue_cutoff"])),
        narrowPeaks_peaks_chr = expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks_chr.narrowPeak"])), sample=SAMPLENAMES),
        summary_file = os.path.join(SUMMARYDIR, "Counts/counts_summary.txt"),
        lorenz_plot_ggplot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf"),
        norm_bw_average = norm_bw_average,
        frip = expand(os.path.join(SUMMARYDIR, "frip/{sample}_frip_score.txt"), sample = SAMPLENAMES),
        peaks_comparison = peaks_comparison
    shell:
        """
        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================

# Generate bwa index if required [OPTIONAL]  ---------------------------------------------
rule generate_genome_index:
    input:
        genome = ancient(genome_fasta)
    output:
        genome_bwt_2bit_64 = ''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),bwa_idx]),
        genome_fai = ''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"])
    params:
        bwa = bwa_version
    threads:
        workflow.cores
    benchmark:
        "benchmarks/generate_genome_index/generate_genome_index---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating the genome index...\\n\033[0m'
        ${{CONDA_PREFIX}}/bin/{params.bwa} index {input.genome}
        samtools faidx {input.genome}
        printf '\033[1;36mGenome index done.\\n\033[0m'
        """



# cutdapat -------------------------------------------------------------------------------
rule cutadapt_PE:
    input:
        get_paired_fq
    output:
        R1_trimm = "01_trimmed_fastq/{SAMPLE}_R1_trimmed.fastq.gz",
        R2_trimm = "01_trimmed_fastq/{SAMPLE}_R2_trimmed.fastq.gz"
    params:
        sample = "{SAMPLE}",
        opts = str(config["cutadapt_trimm_options"]),
        fw_adapter_sequence = str(config["fw_adapter_sequence"]),
        rv_adapter_sequence = str(config["rv_adapter_sequence"])
    log:
        out = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.out",
        err = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.err"
    threads:
        max((workflow.cores - 1), 1)
    benchmark:
        "benchmarks/cutadapt_PE/cutadapt_PE---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: reads trimming...\\n\033[0m'
        mkdir -p 01_trimmed_fastq/logs/

        ${{CONDA_PREFIX}}/bin/cutadapt \
        --nextseq-trim=20 -j {threads} --minimum-length 20 \
        -a {params.fw_adapter_sequence} -A {params.rv_adapter_sequence} {params.opts} \
        -o {output.R1_trimm} -p {output.R2_trimm} {input[0]} {input[1]} > {log.out} 2> {log.err}
        """


# BWA mapping -----------------------------------------------------------------------------
rule BWA_PE:
    input:
        R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][0], "_trimmed.fastq.gz"])),
        R2_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][1], "_trimmed.fastq.gz"])),
        genome_index = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),bwa_idx]))
    output:
        align_summary = "02_BAM/bwa_summary/{SAMPLE}.BWA_summary.txt",
        bam = temp("02_BAM/{SAMPLE}.sorted.bam")
    params:
        bwa_opts = str(config["bwa_options"]),
        sample = "{SAMPLE}",
        bwa = bwa_version,
        genome_fasta = genome_fasta
    threads:
        max((workflow.cores - 1), 1)
    log:
        out = "02_BAM/logs/{SAMPLE}.sort.log"
    benchmark:
        "benchmarks/BWA_PE/BWA_PE---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: bwa mapping...\\n\033[0m'
        mkdir -p 02_BAM/bwa_summary

        ${{CONDA_PREFIX}}/bin/{params.bwa} mem \
        -t {threads} \
        -M {params.genome_fasta} {input.R1_trimm} {input.R2_trimm} |
            $CONDA_PREFIX/bin/samtools fixmate -m - - |
        $CONDA_PREFIX/bin/samtools sort -m 2G -T 02_BAM/{params.sample} -@ {threads} -O bam - > {output.bam} 2> {log.out};
        $CONDA_PREFIX/bin/samtools flagstat {output.bam} > {output.align_summary}
        """


# samtools mapq filter -----------------------------------------------------------------------------
rule MAPQ_MT_filter:
    input:
        source_bam = "02_BAM/{SAMPLE}.sorted.bam"
    output:
        bam_mapq_only = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), ".bam"]))),
        bam_mapq_only_sorted = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam"]))),
        bam_mapq_only_sorted_toFix = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_toFix.bam"]))),
        bam_mapq_only_sorted_index_toFix = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_toFix.bai"]))),
        bam_mapq_only_sorted_index = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bai"]))),
        idxstats_file = "02_BAM/reads_per_chromosome/{SAMPLE}_idxstats_read_per_chromosome.txt",
        bam_mapq_only_sorted_woMT = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam"]))),
        bam_mapq_only_sorted_woMT_index = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam.bai"])))
    params:
        sample = "{SAMPLE}",
        MAPQ_threshold = config["MAPQ_threshold"],
        chr_remove_pattern = chr_remove_pattern
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        fixmate_log = "02_BAM/FixMateInformation_logs/{SAMPLE}_FixMateInformation.log"
    benchmark:
        "benchmarks/MAPQ_MT_filter/MAPQ_MT_filter---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: filtering MAPQ...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools view -@ {threads} -h -q {params.MAPQ_threshold} -F 0x100 -e '([NM] <= 4) && sclen < 15' -f 3 -F 0x0008 {input.source_bam} -o {output.bam_mapq_only}

        $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted_toFix}
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_toFix} {output.bam_mapq_only_sorted_index_toFix}

        $CONDA_PREFIX/bin/gatk FixMateInformation \
        --INPUT {output.bam_mapq_only_sorted_toFix} \
        --OUTPUT {output.bam_mapq_only_sorted} \
        --ASSUME_SORTED false \
        --ADD_MATE_CIGAR true \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT &> {log.fixmate_log}

        $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} > {output.idxstats_file}

        printf '\033[1;36m{params.sample}: Removing MT from BAM...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} | cut -f 1 | grep -v -E '{params.chr_remove_pattern}' | xargs ${{CONDA_PREFIX}}/bin/samtools view -@ {threads} -b {output.bam_mapq_only_sorted} > {output.bam_mapq_only_sorted_woMT}
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_woMT} {output.bam_mapq_only_sorted_woMT_index}
        """



# gatk4 mark duplicates -----------------------------------------------------------------------------
rule gatk4_markdups:
    input:
        bam_mapq_only_sorted = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam"])),
        bam_mapq_only_sorted_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam.bai"]))
    output:
        bam_mdup = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, ".bam"])),
        bai_mdup = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, ".bai"])),
        dup_metrics = "02_BAM/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
        flagstat_filtered = os.path.join("02_BAM/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_flagstat.txt"]))
    params:
        remove_duplicates = str(config["remove_duplicates"]).strip().lower() == "true",
        sample = "{SAMPLE}"
    log:
        out = "02_BAM/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
        err = "02_BAM/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
    threads:
        workflow.cores
    benchmark:
        "benchmarks/gatk4_markdups/gatk4_markdups---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: 'standard' gatk MarkDuplicates...\\n\033[0m'

        mkdir -p 02_BAM/MarkDuplicates_metrics
        mkdir -p 02_BAM/MarkDuplicates_logs
        mkdir -p 02_BAM/flagstat

        $CONDA_PREFIX/bin/gatk MarkDuplicatesWithMateCigar \
        --INPUT {input.bam_mapq_only_sorted} \
        --OUTPUT {output.bam_mdup} \
        --REMOVE_DUPLICATES {params.remove_duplicates} \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT \
        --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

        $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
        """


# fastQC on fastq raw data ----------------------------------------------------------------------------------------
rule fastQC_trimmed_fastq:
    input:
        R1_trimm = "01_trimmed_fastq/{SAMPLE}_R1_trimmed.fastq.gz",
        R2_trimm = "01_trimmed_fastq/{SAMPLE}_R2_trimmed.fastq.gz"
    output:
        html_R1 = "05_quality_controls/trimmed_fastq_fastqc/{SAMPLE}_R1_trimmed_fastqc.html",
        zip_R1 = "05_quality_controls/trimmed_fastq_fastqc/{SAMPLE}_R1_trimmed_fastqc.zip",
        html_R2 = "05_quality_controls/trimmed_fastq_fastqc/{SAMPLE}_R2_trimmed_fastqc.html",
        zip_R2 = "05_quality_controls/trimmed_fastq_fastqc/{SAMPLE}_R2_trimmed_fastqc.zip"
    threads:
        workflow.cores
    benchmark:
        "benchmarks/fastQC_trimmed_fastq/{SAMPLE}fastQC_trimmed_fastq---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPerforming fastQC on trimmed fastq...\\n\033[0m'

        mkdir -p 05_quality_controls/trimmed_fastq_fastqc
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir 05_quality_controls/trimmed_fastq_fastqc {input.R1_trimm} {input.R2_trimm}
        """

rule multiQC_trimmed_fastq:
    input:
        fastqc_zip_trimmed_R1 = expand(os.path.join("05_quality_controls/trimmed_fastq_fastqc","{runs}_R1_trimmed_fastqc.zip"), runs = SAMPLENAMES),
        fastqc_zip_trimmed_R2 = expand(os.path.join("05_quality_controls/trimmed_fastq_fastqc","{runs}_R2_trimmed_fastqc.zip"), runs = SAMPLENAMES)
    output:
        multiqc_fastqc_report = "05_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.html"
    params:
        fastqc_zip_dir = os.path.join(home_dir, "05_quality_controls/trimmed_fastq_fastqc/"),
        out_directory = os.path.join(home_dir, "05_quality_controls/trimmed_fastq_multiQC/"),
        home_dir = home_dir,
        cutadapt_logs = os.path.join(home_dir, "01_trimmed_fastq/logs/"),
        multiqc_fastqc_report_name = "multiQC_report_trimmed_fastq.html"
    log:
        out = os.path.join(home_dir, "05_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.out"),
        err = os.path.join(home_dir, "05_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.err")
    threads: 1
    benchmark:
        "benchmarks/multiQC_trimmed_fastq/multiQC_trimmed_fastq---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPerforming multiQC of trimmed fastq...\\n\033[0m'
        printf '\033[1;36mGenerating multiQC report for trimmed fastq...\\n\033[0m'
 
        cd {params.fastqc_zip_dir}
 
        $CONDA_PREFIX/bin/multiqc -f \
        -o {params.out_directory} \
        -n {params.multiqc_fastqc_report_name} \
        --dirs-depth 2 \
        --dirs ./ {params.cutadapt_logs} > {log.err} 2> {log.out}
 
        cd {params.home_dir}
        """
        
 
# ------------------------------------------------------------------------------
#                                 END pipeline part 1
# ------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# BAM reads shifting and norm BIgWigs
rule bam_shifting_and_RPM_normalization:
    input:
        dedup_BAM = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        dedup_BAM_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
        genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
    output:
        #temp
        dedup_BAM_sortedByName = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLE}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByName.bam"]))),
        dedup_BEDPE_sortedByName = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLE}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByName.bedpe"]))),
        dedup_BED_sortedByName_shifted = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLE}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByName_shifted.bed"]))),
        dedup_BED_sortedByPos_shifted = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLE}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted.bed"]))),
        dedup_BED_sortedByPos_shifted_noBlack = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLE}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted_noBlack.bed"]))),
        dedup_BED_sortedByPos_shifted_noBlack_noIgnoreChr = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLE}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted_noBlack_noIgnoreChr.bed"]))),
        dedup_BedGraph_sortedByPos_shifted_noBlack_RPM = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLE}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted_noBlack_RPM.bedGraph"]))),
        chrSizes = temp("03_Normalization/RPM_normalized/temp_{SAMPLE}_chrSizes_from_genome.txt"),
        # keep
        norm_bw = ''.join(["03_Normalization/RPM_normalized/{SAMPLE}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized.bw"])
    params:
        sample = "{SAMPLE}",
        build_normalization = "03_Normalization/RPM_normalized/bamToBed_log",
        blacklist = BLACKLIST,
        mapq_cutoff = MAPQ,
        minFragmentLength = str(config["bam_features"]["minFragmentLength"]),
        maxFragmentLength = str(config["bam_features"]["maxFragmentLength"]),
        ignore_chr = '|'.join([re.sub('\..*$', '', i) for i in str(config["genomic_annotations"]["ignore_for_normalization"]).split(" ")])
    threads:
        max(math.floor(workflow.cores/2), 1)
    log:
        out = "03_Normalization/RPM_normalized/bamToBed_log/{SAMPLE}_bamToBed.log"
    benchmark:
        "benchmarks/bam_shifting_and_RPM_normalization/bam_shifting_and_RPM_normalization---{SAMPLE}_benchmark.txt"
    priority: 90
    shell:
        """
        printf '\033[1;36m{params.sample}: shifting reads and normalization of the data...\\n\033[0m'
        mkdir -p {params.build_normalization}

        printf '\033[1;36m{params.sample}: re-sorting by read name...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools sort -@ {threads} -n -o {output.dedup_BAM_sortedByName} {input.dedup_BAM}

        printf '\033[1;36m{params.sample}: Bam filtering and conversion to bedPE...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools view -@ {threads} -b -f 3 {output.dedup_BAM_sortedByName} | bedtools bamtobed -i stdin -cigar -bedpe > {output.dedup_BEDPE_sortedByName} 2> {log.out}

        printf '\033[1;36m{params.sample}: Shifting read fragments...\\n\033[0m'
        awk -v OFS='\\t' '{{if($9=="+"){{print $1,$2+4,$6-5,$7,1,$9}}else if($9=="-"){{print $1,$2-5,$6+4,$7,1,$9}}}}' {output.dedup_BEDPE_sortedByName} | awk '(($3 >= $2))' > {output.dedup_BED_sortedByName_shifted}
        sort -k1,1 -k2,2n {output.dedup_BED_sortedByName_shifted} | grep '+' > {output.dedup_BED_sortedByPos_shifted}

        printf '\033[1;36m{params.sample}: Filter blacklist and fragmentSize...\\n\033[0m'
        $CONDA_PREFIX/bin/bedtools intersect -a {output.dedup_BED_sortedByPos_shifted} -b {params.blacklist} -v | awk '(($3-$2 >= {params.minFragmentLength}) && ($3-$2 <= {params.maxFragmentLength}))' > {output.dedup_BED_sortedByPos_shifted_noBlack}

        printf '\033[1;36m{params.sample}: Compute and RPM-Normalize coverage...\\n\033[0m'
        cut -f 1,2,3 {output.dedup_BED_sortedByPos_shifted_noBlack} | grep -v -E '{params.ignore_chr}' > {output.dedup_BED_sortedByPos_shifted_noBlack_noIgnoreChr}
        TOTREADS=$(printf %d\\\\n $(wc -l < {output.dedup_BED_sortedByPos_shifted_noBlack_noIgnoreChr}))
        cut -f1,2 {input.genome_fai} > {output.chrSizes}
        $CONDA_PREFIX/bin/bedtools genomecov -bg -i {output.dedup_BED_sortedByPos_shifted_noBlack} -g {output.chrSizes} | awk -F "\\t" -v tot="$TOTREADS" '{{print $1,$2,$3,$4/((tot*2)/1000000)}}' OFMT="%.20f" > {output.dedup_BedGraph_sortedByPos_shifted_noBlack_RPM}

        printf '\033[1;36m{params.sample}: Computing RPM-normalized bigWig...\\n\033[0m'
        $CONDA_PREFIX/bin/bedGraphToBigWig {output.dedup_BedGraph_sortedByPos_shifted_noBlack_RPM} {output.chrSizes} {output.norm_bw}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Average bigWigs by group
rule compute_bigwigAverage:
    input:
        norm_bw = expand(''.join(["03_Normalization/RPM_normalized/{sample}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized.bw"]), sample = SAMPLENAMES)
    output:
        norm_bw_average = expand(''.join(["03_Normalization/RPM_normalized_merged/{group}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized_merged.bs", str(config["differential_TF_binding"]["merged_bigwig_binSize"]), ".bw"]), group = GROUPNAMES)
    params:
        groups_tb = str(config["differential_TF_binding"]["sample_groups_table"]),
        binSize = str(config["differential_TF_binding"]["merged_bigwig_binSize"]),
        groups = ' '.join(GROUPNAMES),
        bw_input_suffix = "_mapq"+MAPQ+"_woMT_"+DUP+"_shifted_RPM.normalized.bw",
        bw_output_suffix = "_mapq"+MAPQ+"_woMT_"+DUP+"_shifted_RPM.normalized_merged.bs"+str(config["differential_TF_binding"]["merged_bigwig_binSize"])+".bw",
        blacklist = BLACKLIST,
    threads:
        workflow.cores
    log:
        out = expand("03_Normalization/RPM_normalized_merged/log/{group}_bigwigAverage.log", group = GROUPNAMES)
    benchmark:
        "benchmarks/compute_bigwigAverage/compute_bigwigAverage---allGroups_benchmark.txt"
    priority: 50
    shell:
        """
        printf '\033[1;36mMerging bigwigs by group...\\n\033[0m'

        for i in {params.groups}
        do
            SAMPLE=$(grep $i {params.groups_tb} | cut -f 1 | uniq)
            BIGWGS=$(echo $(for s in $SAMPLE; do echo 03_Normalization/RPM_normalized/${{s}}{params.bw_input_suffix}; done))

            echo '  - '${{i}}': '$SAMPLE

            $CONDA_PREFIX/bin/bigwigAverage \
            -b $BIGWGS \
            --binSize {params.binSize} \
            --outFileName 03_Normalization/RPM_normalized_merged/${{i}}{params.bw_output_suffix} \
            --outFileFormat bigwig \
            --blackListFileName {params.blacklist} \
            -p {threads} &> 03_Normalization/RPM_normalized_merged/log/${{i}}_bigwigAverage.log
        done
        """
# ----------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------
# FastQC on BAMs
rule fastQC_BAMs:
    input:
        dedup_BAM = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        dedup_BAM_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"]))
    output:
        html = os.path.join("02_BAM_fastQC/", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, "_fastqc.html"])),
        zip = os.path.join("02_BAM_fastQC/", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, "_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = os.path.join(config["output_directory"], "02_BAM_fastQC/"),
        sample = "{SAMPLE}"
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/fastQC_BAMs/fastQC_BAMs---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Performing fastQC on deduplicated bam...\\n\033[0m'
        mkdir -p 02_BAM_fastQC
        
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir {params.fastQC_BAMs_outdir} {input.dedup_BAM}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule fragment_size_distribution:
    input:
        BAM = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        BAM_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
    output:
        plot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/{SAMPLE}_fragment_size_distribution.pdf"),
        fragmentSize_RawFragmentLengths = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{SAMPLE}_fragmentSize_RawFragmentLengths.txt"),
        fragmentSize_metrics = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{SAMPLE}_fragmentSize_metrics.txt")
    params:
        build_summary_directory = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log"),
        build_summary_directory_table = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize"),
        sample = "{SAMPLE}",
        plotFormat = "pdf",
        binSize = str(config["quality_controls"]["fragmentSize_window_length"]),
        blacklist = BLACKLIST,
        maxFragmentLength = config["bam_features"]["maxFragmentLength"]
    log:
        out = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/{SAMPLE}_fragmentSize_log.out"),
        err = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/{SAMPLE}_fragmentSize_log.err")
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/fragment_size_distribution/fragment_size_distribution---{SAMPLE}_benchmark.txt"
    priority: -1
    shell:
        """
        printf '\033[1;36m{params.sample}: Plotting the fragment size distribution...\\n\033[0m'

        mkdir -p {params.build_summary_directory}
        mkdir -p {params.build_summary_directory_table}

        $CONDA_PREFIX/bin/bamPEFragmentSize \
        -p {threads} \
        -b {input.BAM} \
        --plotFileFormat {params.plotFormat} \
        --plotTitle {params.sample} \
        --samplesLabel {params.sample} \
        --binSize {params.binSize} \
        --maxFragmentLength {params.maxFragmentLength} \
        --blackListFileName {params.blacklist} \
        --outRawFragmentLengths {output.fragmentSize_RawFragmentLengths} \
        --table {output.fragmentSize_metrics} \
        --histogram {output.plot} > {log.out} 2> {log.err}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule fragment_size_distribution_report:
    input:
        plots = expand(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/{sample}_fragment_size_distribution.pdf"), sample = SAMPLENAMES)
    output:
        replot_script = temp(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/replot_script.R")),
        report_pdf = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"),
        report_ggplot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots_ggplot.version.pdf")
    params:
        distribution_plots_pattern = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/*_fragment_size_distribution.pdf"),
        dir = os.path.join(home_dir,""),
        summary_dir = SUMMARYDIR,
        maxFragmentLength = config["bam_features"]["maxFragmentLength"]
    threads: 1
    log:
      ggplot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/ggplot_replotting.log"),
      pdfcombine = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/pdfcombine.log")
    benchmark:
        "benchmarks/fragment_size_distribution_report/fragment_size_distribution_report---benchmark.txt"
    shell:
        """
        printf '\033[1;36mMerging fragmentSizeDistribution reports in a unique PDF...\\n\033[0m'
        $CONDA_PREFIX/bin/TOBIAS MergePDF --input {params.distribution_plots_pattern} --output {output.report_pdf} &> {log.pdfcombine}

        printf '\033[1;36mReplotting fragmentSizeDistribution reports in R (ggplot version)...\\n\033[0m'
        echo "tb = do.call(rbind, lapply(list.files('{params.dir}{params.summary_dir}fragmentSizeDistribution_plots/table_and_fragmentSize', pattern = 'RawFragmentLengths', full.names = T), function(x)(read.delim(x, h=T, skip=1))))" > {output.replot_script}
        echo "n.samples = length(unique(tb[,3]))" >> {output.replot_script}
        echo "plot = ggplot2::ggplot(data = tb, ggplot2::aes(x = Size, y = Occurrences, color = Sample)) + ggplot2::geom_smooth(method = 'loess', formula = y ~ x, span = 0.05, show.legend = F, se = F, color = 'navyblue', linewidth = 0.5) + ggplot2::xlim(c(1,{params.maxFragmentLength})) + ggplot2::theme_classic() + ggplot2::facet_wrap(~Sample, scale='free', ncol = floor(sqrt(n.samples))) + ggplot2::theme(axis.ticks = ggplot2::element_line(color ='black'), axis.text = ggplot2::element_text(color = 'black'), strip.background = ggplot2::element_blank())" >> {output.replot_script}
        echo "pdf(file = '{params.dir}{output.report_ggplot}', width = floor(sqrt(n.samples)) * 2.7, height = ceiling(n.samples / floor(sqrt(n.samples))) * 1.5)" >> {output.replot_script}
        echo "print(plot)" >> {output.replot_script}
        echo "invisible(dev.off())" >> {output.replot_script}

        $CONDA_PREFIX/bin/Rscript {output.replot_script} &> {log.ggplot}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# PeakCalling on bams (MACS3 callpeak)
rule peakCalling_MACS3:
    input:
        dedup_BAM_sorted = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        dedup_BAM_sorted_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"]))
    output:
        peaks_xls = os.path.join(PEAKSDIR, ''.join(["{SAMPLE}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.xls"])),
        narrowPeaks = os.path.join(PEAKSDIR, ''.join(["{SAMPLE}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])),
        narrowPeaks_chr = os.path.join(PEAKSDIR, ''.join(["{SAMPLE}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks_chr.narrowPeak"])),
        summits = os.path.join(PEAKSDIR, ''.join(["{SAMPLE}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_summits.bed"]))
    params:
        genomeSize = str(config["genomic_annotations"]["effective_genomeSize"]),
        blacklist = BLACKLIST,
        peak_caller = PEAKCALLER.lower(),
        peaks_dir = PEAKSDIR,
        qValue = str(config["peak_calling"]["qValue_cutoff"]),
        basename = ''.join(["{SAMPLE}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"])]),
        summits = SUMMITS,
        sample = "{SAMPLE}"
    log:
        out = os.path.join(''.join([PEAKSDIR,"log/"]), ''.join(["{SAMPLE}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), ".log"]))
    benchmark:
        "benchmarks/peakCalling_MACS3/peakCalling_MACS3---{SAMPLE}_benchmark.txt"
    priority: 100
    shell:
        """
        printf '\033[1;36m{params.sample}: Calling peaks by {params.peak_caller}...\\n\033[0m'

      $CONDA_PREFIX/bin/{params.peak_caller} callpeak \
        -t {input.dedup_BAM_sorted} \
        -g {params.genomeSize} \
        -n {params.basename} \
        -q {params.qValue} \
        -f BAMPE -B \
        --nomodel --shift -75 --extsize 150 --call-summits \
        --outdir {params.peaks_dir} \
        --keep-dup all \
        {params.summits} 2> {log.out}
 
        # add chr to peak files
        $CONDA_PREFIX/bin/bedtools subtract -nonamecheck -a {output.narrowPeaks} -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > {output.narrowPeaks_chr}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Computation of the counts summary table
rule counts_summary:
    input:
        flagstat_filtered = expand(os.path.join("02_BAM/flagstat/", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, "_flagstat.txt"])), sample = SAMPLENAMES),
        dedup_BED_sortedByPos_shifted_noBlack = expand(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{sample}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted_noBlack.bed"])), sample=SAMPLENAMES),
        peaks_file = expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])), sample = SAMPLENAMES),
        norm_bw = expand("03_Normalization/RPM_normalized/{sample}_mapq{mapq}_woMT_{dup}_shifted_RPM.normalized.bw", sample=SAMPLENAMES, dup=DUP, mapq=MAPQ)
    output:
        summary_file = os.path.join(SUMMARYDIR, "Counts/counts_summary.txt"),
        summary_file_temp = temp(os.path.join(SUMMARYDIR, "Counts/summary_file.temp"))
    params:
        build_summary_directory = os.path.dirname(SUMMARYDIR),
        sample_list = str(' '.join(SAMPLENAMES)),
        peaks_dir = PEAKSDIR,
        FRiP_threshold = config["peak_calling"]["FRiP_threshold"],
        MAPQ = MAPQ,
        DUP = DUP
    threads: 1
    benchmark:
        "benchmarks/counts_summary/counts_summary---benchmark.txt"
    priority: 80
    shell:
        """
        mkdir -p {params.build_summary_directory}/Counts/subread_featureCounts_output/

        printf '\033[1;36mGeneration of a general counts summary table...\\n\033[0m'
        printf Sample'\\t'dedup_BAM'\\t'shifted_readsM'\\t'loss_post_shifting'\\t'n.peaks'\\t'FRiP.perc'\\t'FRiP.quality'\\n' > {output.summary_file}

        for NAME in {params.sample_list}
        do
            printf '\033[1;36m     - %s: adding stats to summary table...\\n\033[0m' $NAME

            mkdir -p {params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/

            woMT_BAM=$(grep mapped 02_BAM/flagstat/${{NAME}}_mapq{params.MAPQ}_sorted_woMT_{params.DUP}_flagstat.txt | head -n 1 | cut -f 1 -d ' ')

            dedupBAM=$(grep mapped 02_BAM/flagstat/${{NAME}}_mapq{params.MAPQ}_sorted_woMT_{params.DUP}_flagstat.txt | head -n 1 | cut -f 1 -d ' ')

            halfShiftedBEDPE=$(printf %d\\\\n $(wc -l < 03_Normalization/RPM_normalized/temp_${{NAME}}_mapq{params.MAPQ}_sorted_woMT_{params.DUP}_sortedByPos_shifted_noBlack.bed))
            shiftedBEDPE=$(echo "$halfShiftedBEDPE * 2" | bc)
            lossReads=$(echo "$dedupBAM - $shiftedBEDPE" | bc)

            peaks=$(wc -l {params.peaks_dir}${{NAME}}*_peaks.narrowPeak | cut -f 1 -d ' ')

            awk 'BEGIN{{FS=OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' {params.peaks_dir}${{NAME}}*.*Peak > {params.peaks_dir}${{NAME}}.saf
            FEATURECOUNTSLOG={params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/${{NAME}}.readCountInPeaks.log
            featureCounts -p -a {params.peaks_dir}${{NAME}}.saf -F SAF -o {params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/${{NAME}}.readCountInPeaks 02_BAM/${{NAME}}_mapq*_sorted_woMT_{DUP}.bam 2> ${{FEATURECOUNTSLOG}}
            rm {params.peaks_dir}${{NAME}}.saf
            frip=$(grep 'Successfully assigned alignments' ${{FEATURECOUNTSLOG}} | sed -e 's/.*(//' | sed 's/%.*$//')
            fripScore=$(echo $frip | sed 's/\\..*$//')
            fripLabel=$(if [ $fripScore -ge {params.FRiP_threshold} ]; then echo 'good'; else echo 'bad'; fi)

            printf ${{NAME}}'\\t'$dedupBAM'\\t'$shiftedBEDPE'\\t'$lossReads'\\t'$peaks'\\t'$frip'\\t'$fripLabel'\\n' >> {output.summary_file}
        done

        uniq -u {output.summary_file} > {output.summary_file_temp}
        (head -n 1 {output.summary_file_temp} && tail -n +2 {output.summary_file_temp} | sort -k 1) > {output.summary_file}
        """



# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
rule frip_score:
    input:
        bam = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        peaks = expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])), sample=SAMPLENAMES)
    output:
        frip_score = os.path.join(SUMMARYDIR, "frip/{SAMPLE}_frip_score.txt")
    shell:
        """
        reads_in_peaks=$(bedtools intersect -a {input.bam} -b {input.peaks} -bed -c -f 0.20 | awk '{{sum += $NF}} END {{print sum}}')
        total_mapped=$(samtools view -c {input.bam})
        
        frip=$(echo "scale=4; $reads_in_peaks / $total_mapped" | bc)
        
        echo "{wildcards.SAMPLE}\t$frip" > {output.frip_score}
        """
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Perform multiQC for BAMs
rule multiQC_BAMs:
    input:
        flagstat_filtered = expand(os.path.join("02_BAM/flagstat/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_flagstat.txt"])), sample = SAMPLENAMES),
        BAM_fastqc_zip = expand(os.path.join("02_BAM_fastQC/", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, "_fastqc.zip"])), sample=SAMPLENAMES),
        narrowPeaks = expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])), sample=SAMPLENAMES)
    output:
        multiqcReportBAM = os.path.join(SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/multiQC_report_BAMs_", DUP, ".html"]))
    params:
        fastQC_BAM_reports_dir = "02_BAM_fastQC/",
        picard_metrics_dir = "02_BAM/MarkDuplicates_metrics/",
        dedup_BAM_flagstat_dir = "02_BAM/flagstat/",
        macs_dir = PEAKSDIR,
        multiQC_BAM_outdir = os.path.join(config["output_directory"], SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/"]))
    log:
        out = os.path.join(SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/multiQC_report_BAMs_", DUP, ".out"])),
        err = os.path.join(SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/multiQC_report_BAMs_", DUP, ".err"]))
    benchmark:
        "benchmarks/fmultiQC_BAMs/multiQC_BAMs---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating multiQC report from deduplicated bam quality test...\\n\033[0m'

        $CONDA_PREFIX/bin/multiqc -f \
        --outdir {params.multiQC_BAM_outdir} \
        -n multiQC_report_BAMs_{DUP}.html \
        --dirs \
        02_BAM/MarkDuplicates_metrics 02_BAM/flagstat \
        {params.fastQC_BAM_reports_dir} \
        {params.picard_metrics_dir} \
        {params.dedup_BAM_flagstat_dir} \
        {params.macs_dir} > {log.out} 2> {log.err}
        """


# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of Lorenz curves
rule Lorenz_curve:
    input:
        dedup_bam_sorted = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        dedup_bam_sorted_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"]))
    output:
        lorenz_plot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_plots/{SAMPLE}_Lorenz_curve_deeptools.plotFingreprint.pdf"),
        lorenz_metrics = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_metrics/{SAMPLE}_Lorenz_quality.metrics_deeptools.plotFingreprint.txt"),
        lorenz_counts = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_counts/{SAMPLE}_Lorenz_raw.counts_deeptools.plotFingreprint.txt")
    params:
        #all_bams = ' '.join(expand(os.path.join(''.join(["02_BAM/{sample}_mapq", MAPQ ,"_sorted_woMT_", DUP, ".bam"])), sample=SAMPLENAMES)),
        #labels = ' '.join(SAMPLENAMES),
        labels = str("{SAMPLE}"),
        blacklist = BLACKLIST,
        binSize = config["quality_controls"]["plotFingerprint"]["binSize"],
        sampledRegions = config["quality_controls"]["plotFingerprint"]["sampledRegions"],
        extra_params = config["quality_controls"]["plotFingerprint"]["extra_parameters"]
    threads:
        max(math.floor((workflow.cores-1)/2), 1)
    log:
        out = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/log/{SAMPLE}_deeptools_plotFingreprint.log")
    benchmark:
        "benchmarks/Lorenz_curve/Lorenz_curve---{SAMPLE}_benchmark.txt"
    priority: -10
    shell:
        """
        printf '\033[1;36m{params.labels}: plotting Lorenz curves-Fingerprint...\\n\033[0m'

        $CONDA_PREFIX/bin/plotFingerprint \
        --bamfiles {input.dedup_bam_sorted} \
        --plotFile {output.lorenz_plot} \
        --labels {params.labels} \
        --blackListFileName {params.blacklist} \
        --binSize {params.binSize} \
        --numberOfSamples {params.sampledRegions} \
        --outQualityMetrics {output.lorenz_metrics} \
        --outRawCounts {output.lorenz_counts} \
        -p {threads} {params.extra_params} &> {log.out}
        """


rule Lorenz_curve_merge_plots:
    input:
        lorenz_plots = expand(os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_plots/{sample}_Lorenz_curve_deeptools.plotFingreprint.pdf"), sample=SAMPLENAMES)
    output:
        lorenz_plot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples_combined.pdf"),
        lorenz_plot_ggplot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf")
    params:
        lorenz_plots_pattern = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_plots/*_Lorenz_curve_deeptools.plotFingreprint.pdf"),
        dir = os.path.join(home_dir,""),
        summary_dir = SUMMARYDIR,
    log:
        pdfcombine = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/log/LorenzCurve_plotFingreprint_mergePdf.log"),
        ggplotcombine = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/log/LorenzCurve_plotFingreprint_ggplot.merge.plots.log")
    threads: 1
    benchmark:
        "benchmarks/Lorenz_curve/Lorenz_curve_merge_plots---benchmark.txt"
    priority: -5
    shell:
        """
        printf '\033[1;36mCombine Lorenz curves-Fingerprint for all samples...\\n\033[0m'
        $CONDA_PREFIX/bin/TOBIAS MergePDF --input {params.lorenz_plots_pattern} --output {output.lorenz_plot} &> {log.pdfcombine}

        printf '\033[1;36mMake combined Lorenz curves-Fingerprint plot...\\n\033[0m'
        $CONDA_PREFIX/bin/Rscript \
        -e "require(dplyr)" \
        -e "tables = list.files(path = '{params.dir}06_Overall_quality_and_info/LorenzCurve_plotFingreprint/lorenz_counts', pattern = '.plotFingreprint.txt', full.names = T)" \
        -e "combined_table = data.frame()" \
        -e "for (i in 1:length(tables)) (combined_table = rbind(combined_table, dplyr::mutate(read.delim(tables[i], skip = 2, h=F), sample = gsub('_Lorenz_raw[.]counts_deeptools[.]plotFingreprint[.]txt','',basename(tables[i]))) %>% dplyr::rename(counts = V1) %>% dplyr::arrange(counts) %>% dplyr::mutate(cumulative_sum = cumsum(counts), rank = (1:nrow(.))/nrow(.)) %>% dplyr::mutate(cumulative_sum = cumulative_sum/max(cumulative_sum))))" \
        -e "pdf('{params.dir}{output.lorenz_plot_ggplot}', width = 8, height = 6.5)" \
        -e "ggplot2::ggplot(data = combined_table, ggplot2::aes(x = rank, y = cumulative_sum, color = sample)) + ggplot2::geom_line() + ggplot2::ggtitle('Fingerprints (Lorenz curves) all samples') + ggplot2::xlim(c(0,1)) + ggplot2::xlab('Normalized rank') + ggplot2::ylab('Fraction with reference to the bin with highest coverage') + ggplot2::theme_classic() + ggplot2::theme(axis.text = ggplot2::element_text(color = 'black'), axis.ticks = ggplot2::element_line(color = 'black'))" \
        -e "invisible(dev.off())" &> {log.ggplotcombine}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Absolute peaks file and relative matrix score generation for called peaks
rule all_peaks_file_and_score_matrix:
    input:
        norm_bw = expand("03_Normalization/RPM_normalized/{sample}_mapq{mapq}_woMT_{dup}_shifted_RPM.normalized.bw", sample=SAMPLENAMES,  dup=DUP, mapq=MAPQ),
        peaks_file = expand(os.path.join(str(PEAKSDIR), "{sample}_mapq{mapq}_woMT_{dup}_qValue{qValue}_peaks.narrowPeak"), sample=SAMPLENAMES, mapq=MAPQ, dup=DUP, qValue=str(config["peak_calling"]["qValue_cutoff"]))
    output:
        concatenation_bed = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation.bed")),
        concatenation_bed_sorted = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation_sorted.bed")),
        concatenation_bed_collapsed = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation_collapsed.bed")),
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed"),
        score_matrix_peaks = os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_", PEAKCALLER, ".npz"])),
        score_matrix_peaks_table = os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_table_", PEAKCALLER, ".tsv"]))
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "Sample_comparisons/Peak_comparison/"])),
        peak_caller = PEAKCALLER,
        peaks_dir = PEAKSDIR,
        labels = ' '.join(SAMPLENAMES),
        blacklist = BLACKLIST
    threads:
        max((workflow.cores-1), 1)
    benchmark:
        "benchmarks/all_peaks_file_and_score_matrix/all_peaks_file_and_score_matrix---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating a file result of the merge of all the {params.peak_caller} peaks...\\n\033[0m'

        mkdir -p {params.make_directory}

        cat {params.peaks_dir}*.*Peak >> {output.concatenation_bed}
        sort -V -k1,1 -k2,2 -k5,5 {output.concatenation_bed} > {output.concatenation_bed_sorted}

        $CONDA_PREFIX/bin/bedtools merge -i {output.concatenation_bed_sorted} | uniq > {output.concatenation_bed_collapsed}
        sort -V -k1,1 -k2,2 -k5,5 {output.concatenation_bed_collapsed} > {output.concatenation_bed_collapsed_sorted}


        printf '\033[1;36mComputing the score matrix for all the {params.peak_caller} peaks per each sample...\\n\033[0m'

        $CONDA_PREFIX/bin/multiBigwigSummary BED-file \
        -p {threads} \
        -b {input.norm_bw} \
        -o {output.score_matrix_peaks} \
        --BED {output.concatenation_bed_collapsed_sorted} \
        --blackListFileName {params.blacklist} \
        --outRawCounts {output.score_matrix_peaks_table} \
        --labels {params.labels}
        """
# ----------------------------------------------------------------------------------------



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ------------------------------------------------------------------------------
#                                 END PIPELINE                                 #
# ------------------------------------------------------------------------------
