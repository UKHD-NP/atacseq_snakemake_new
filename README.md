# ATAC-seq Snakemake Pipeline 

 ```ascii
           /^\/^\
         _|__|  O|     
\/     /~     \_/ \     
 \____|__________/  \    
        \_______      \       
                `\     \         
                  |     |      
                 /      /    
                /     /     
              /      /     
             /     /
           /     /
          /     /
         (      (
          \      ~-____-~ -. .-.   .-. .-.   .-. .-.   .        
            ~-_           ||\|||\ /|||\|||\ /|||\|||\ /|
               ~-_        |/ \|||\|||/ \|||\|||/ \|||\|| 
                  ______  ~   `-~ `-`   `-~ `-`   `-~ `-  
```

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a widely used method to profile genome-wide chromatin accessibility. By leveraging the Tn5 transposase, which preferentially inserts sequencing adapters into open chromatin regions, ATAC-seq provides insights into regulatory elements such as promoters, enhancers, and transcription factor binding sites.

This is a modular Snakemake workflow for ATAC-seq data analysis, inspired by [snakeATAC](https://sebastian-gregoricchio.github.io/snakeATAC/) and best practices from ENCODE and nf-core ATAC-seq. It automates the analysis of ATAC-seq data from raw sequencing reads to processed results, including trimmed FASTQ files, aligned BAM files, peak calls, normalized coverage tracks, and quality control metrics. 

## Installation

1. **Place yourself in the directory where the repository should be downloaded**:
    ```bash
    cd </target/repository/folder>
    ```

2. **Clone the repository**:
   ```bash
   git clone https://github.com/UKHD-NP/atacseq_snakemake.git
   cd atacseq_snakemake
   ```
   or click on *Code > Download ZIP* on the GitHub page

3. **Install dependencies**:
   This pipeline uses Conda/Mamba for environment management. Install the required environment as follows:
   ```bash
   mamba env create -f workflow/envs/initial.yaml
   ```
4. **Activate the environment**:
    ```bash
   mamba activate atacseq
   ```

## How to Run the Pipeline

Once the environment is activated, you can run the pipeline using the following command. It is important to define all the paths to these files correctly:

1. Perform a dry run to check if everything is set up correctly:
    ```bash
    snakemake \
    --cores  \
    --use-conda \
    -s workflow/atacseq.snakefile \
    --configfile config/atacseq_config.yaml \
    --config \
    samplesheet="/path/to/samplesheet.csv" \
    output_directory="/path/to/output_directory" \
    genome_fasta="/path/to/genome.fa" \
    blacklist="/path/to/blacklist.bed" \
    -n
    ```
2. If no errors arise during the dry run *(it does not mean that the pipeline will run without any issues)*, execute the full pipeline:
    ```bash
    snakemake \
    --cores  \
    --use-conda \
    -s workflow/atacseq.snakefile \
    --configfile config/atacseq_config.yaml \
    --config \
    samplesheet="/path/to/samplesheet.csv" \
    output_directory="/path/to/output_directory" \
    genome_fasta="/path/to/genome.fa" \
    blacklist="/path/to/blacklist.bed"
    ```

### Sample Sheet Structure
The sample sheet should be a CSV file with the following columns:
| sample | replicate | fq1                       | fq2                       |
|--------|-----------|---------------------------|---------------------------|
| name   | REP1      | path/to/sample_R1.fastq.gz| path/to/sample_R2.fastq.gz|


### Parameters in configfile (atacseq_config.yaml)

| Parameter                              | Description & Usage                                                                                                                | Default Value         |
|----------------------------------------|------------------------------------------------------------------------------------------------------------------------------------|-----------------------|
| General Workflow                       |                                                                                                                                    |                       |
| samplesheet                            | Path to samplesheet CSV (sample info and FASTQ paths). Required.                                                                    | -                     |
| output_directory                       | Main output directory for results. Required.                                                                                        | -                     |
| genome_fasta                           | Reference genome FASTA file. Required.                                                                                              | -                     |
| blacklist                              | ENCODE blacklist BED file to exclude problematic genomic regions. Required.                                                         | -                     |
| fastq_suffix                           | FASTQ file extension (e.g.,.fastq.gz). Must match input files.                                                                     | ".fastq.gz"           |
| read_suffix                            | Read pair identifiers (e.g.,_R1,_R2). Match your FASTQ naming convention.                                                          | ['_R1', '_R2']        |
| Trimming                               |                                                                                                                                    |                       |
| cutadapt_trimm_options                 | Add custom trimming parameters                                                                                                     | ''                    |
| fw_adapter_sequence                    | Forward adapter sequence for Tn5 transposase trimming.                                                                             | "CTGTCTCTTATACACATCT" |
| rv_adapter_sequence                    | Reverse adapter sequence for Tn5 transposase trimming.                                                                             | "CTGTCTCTTATACACATCT" |
| BWA Mapping                            |                                                                                                                                    |                       |
| bwa_options                            | Custom alignment options for BWA-mem2.                                                                                             | ''                    |
| BAM Filtering                          |                                                                                                                                    |                       |
| remove_duplicates                      | Remove PCR duplicates using GATK MarkDuplicates. Set to False to keep duplicates.                                                  | True                  |
| MAPQ_threshold                         | Minimum mapping quality score (MAPQ). Filters ambiguous alignments (recommended: 20). Increase for stricter alignment (range 0-60) | 20                    |
| remove_other_chromosomes_pattern       | Regex pattern to exclude unwanted contigs (e.g., mitochondrial DNA, viral sequences).                                              | "CMV\|HBV..."         |
| bam_features.minFragmentLength         | Minimum fragment size for paired-end reads (filters short fragments).                                                              | 0                     |
| bam_features.maxFragmentLength         | Maximum fragment size for paired-end reads (filters overly long fragments).                                                        | 2000                  |
| Genome Annotations                     |                                                                                                                                    |                       |
| genome_id                              | Genome identifier (informational only).                                                                                            | "hg38"                |
| effective_genomeSize                   | Effective genome size for MACS3 peak calling. Must match reference genome.                                                         | 2913022398            |
| ignore_for_normalization               | Chromosomes excluded from coverage normalization (e.g., sex chromosomes, contigs).                                                 | "X Y MT..."           |
| Peak Calling                           |                                                                                                                                    |                       |
| qValue_cutoff                          | MACS3 q-value significance threshold for peak calling (lower = more stringent).                                                    | 0.05                  |
| call_summits                           | Whether to report precise peak summits (useful for motif analysis).                                                                | True                  |
| FRiP_threshold                         | Minimum Fraction of Reads in Peaks (FRiP) score (% reads in peaks).                                                                | 20                    |
| Quality Control                        |                                                                                                                                    |                       |
| fragmentSize_window_length             | Bin size for fragment size distribution plots (larger bins = smoother distributions).                                              | 1000bp                |
| multiBigwigSummary_binning_window_size | Window size for correlation heatmaps and signal comparison across samples.                                                         | 10000bp               |
| plotFingerprint.binSize                | Resolution of fingerprint plots (higher = more detail, slower processing).                                                         | 500bp                 |

## Pipeline Overview

#### **1. Initialization**
The pipeline begins with an initialization rule (`AAA_initialization`), which ensures that all required outputs are properly organized.

#### **2. Genome Index Generation**
If a genome index is not already available, the rule `generate_genome_index` creates it using BWA-mem2 and with samtools. This step ensures that the reference genome is prepared for efficient read alignment. The output includes `.bwt`, `.fai`, and other index files necessary for mapping.

#### **3. Adapter Trimming**
The `cutadapt_PE` rule trims sequencing adapters from paired-end reads using Cutadapt. 
The Cutadapt parameters in this ATAC-seq pipeline are optimized for typical Nextera-based library preparation and sequencing on NextSeq/NovaSeq platforms. The `--nextseq-trim=20` parameter addresses the two-color chemistry’s tendency to produce poly-G artifacts by trimming 3’ bases with quality scores below 20. The `--minimum-length 20` setting filters out reads shorter than 20 bp, removing adapter dimers and fragments too short for meaningful peak calling. On the basis of [this](https://pmc.ncbi.nlm.nih.gov/articles/PMC10035359/) and [this](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html) practical example was the value chosen for this parameter. The adapter sequences (set in the configfile) match the Illumina adapter sequences, ensuring proper removal of ligated adapters.

For Illumina HiSeq/MiSeq (non-NextSeq), replace `--nextseq-trim` with generic quality trimming (e.g., `-q 20`), as poly-G artifacts are absent. If using TruSeq adapters, update `-a`/`-A` to the correct sequences in the configfile. For lower-quality data (e.g., degraded samples), relax `-q` to 20–25 to avoid excessive read loss. Always validate parameters using FastQC or MultiQC to ensure adapter removal and fragment size distribution align with expectations.

![Base Content in the sequences](resources/fastqc_base_content.png)
***Figure 1: Expected [FastQC base content profiles](https://bioinformaticamente.com/2024/12/05/comprehensive-guide-to-atac-seq-data-quality-control/) for (A) raw and (B) trimmed reads. Post-trimming (B) shows reduced 3’-end G-content spikes caused by NextSeq poly-G artifacts.***

![Adapter content in the sequences](resources/fastqc_adapter_content.png)
***Figure2: [Adapter content profiles](https://bioinformaticamente.com/2024/12/05/comprehensive-guide-to-atac-seq-data-quality-control/) for (A) raw and (B) trimmed reads. Successful adapter removal eliminates >99% of Nextera sequences.***

#### **4. Read Alignment**
The `BWA_PE` rule aligns the trimmed reads to the reference genome using BWA-mem2. This aligner was chosen for its widespread use in the literature and its optimal balance of accuracy and speed showed in different [benchmarks](https://www.nature.com/articles/s41467-021-26865-w).

#### **5. BAM Filtering**
The `MAPQ_MT_filter` rule filters aligned reads based on mapping quality (MAPQ) and removes unwanted chromosomes (e.g., mitochondrial DNA). Based on [ENCODE guidelines](https://www.encodeproject.org/atac-seq/), each replicate should retain ***≥50 million non-duplicate, non-mitochondrial aligned reads*** for paired-end analysis. Additionally, duplicate reads are marked or removed using GATK's `MarkDuplicatesWithMateCigar` in the subsequent `gatk4_markdups` rule.

How BAM files and reads are filtered throughout the pipeline:
- **Remove mitochondrial DNA reads**  
    Reads mapping to `chrM` or mitochondrial contigs are excluded using SAMtools.
- **Exclude blacklisted regions**  
    Reads overlapping ENCODE blacklist regions are removed using BEDtools.
- **Remove duplicate reads**  
    PCR duplicates are marked and optionally removed using GATK `MarkDuplicatesWithMateCigar`.
- **Filter by mapping quality (MAPQ)**  
    Reads with MAPQ scores below the threshold (default: 20, set in the configfile) are discarded using SAMtools.
- **Exclude secondary alignments**  
    Only primary alignments (`-F 0x100`) are retained using SAMtools.
- **Remove unmapped reads**  
    Reads flagged as unmapped (`-F 0x4`) are filtered out using SAMtools.
- **Exclude multimapped reads**  
    Reads mapping to multiple locations are removed based on MAPQ scoring.
- **Filter by mismatches**  
    Reads with more than 4 mismatches (using the NM tag) are excluded during filtering using SAMtools.
- **Filter by fragment size**  
    Fragments outside the range of 0–2000 bp are removed using BEDtools and custom scripts.
- **Exclude improperly paired reads**  
    Only properly paired reads (`-f 3`) are retained for downstream analysis.
- **Handle paired-end read inconsistencies**  
    Reads where only one mate fails any of the above criteria are excluded during filtering.


#### **6. Quality Control on Trimmed Reads**
The `fastQC_trimmed_fastq` rule performs quality control checks on trimmed FASTQ files using FastQC. A MultiQC report (`multiQC_trimmed_fastq`) aggregates these QC metrics into a single HTML file for easy review.

#### **7. Quality Control on BAM Files**
The `fastQC_BAMs` rule assesses the quality of BAM files. A MultiQC report (`multiQC_BAMs`) consolidates these metrics alongside alignment statistics and peak-calling results.

Key metrics from [ENCODE](https://www.encodeproject.org/atac-seq/) include:  
- **Alignment rate** >95% (acceptable >80%) 
- **FRiP score** >0.3 (acceptable >0.2)  
- **Nucleosome-free regions** must be detectable in called peaks 

#### **8. Fragment Size Distribution**
The `fragment_size_distribution` rule calculates fragment size distributions for each sample, which provides insights into nucleosome positioning and library complexity. A combined plot (`fragment_size_distribution_report`) summarizes these distributions across all samples.
High-quality ATAC-seq data based on [ENCODE](https://www.encodeproject.org/atac-seq/) must show:
-**Nucleosome-free region (NFR) peak** at ~50 bp
-**Mononucleosome peak** between 147–294 bp
-Clear separation of di-/tri-nucleosome peaks (optional)


![Fragment-distribution](resources/fragmentSize_distribution_examples.svg)
***Figure 3: [Fragment distributions](https://sebastian-gregoricchio.github.io/snakeATAC/) showing (left) ideal profile with NFR and mono-/di-nucleosome peaks vs. (right) noisy data lacking clear nucleosomal patterning.***

#### **9. Read Shifting and Normalization**
The `bam_shifting_and_RPM_normalization` rule shifts paired-end reads to account for Tn5 transposase binding offset and generates RPM-normalized bigWig files for visualization in genome browsers.

#### **10. Peak Calling**
The `peakCalling_MACS3` rule identifies regions of open chromatin by calling peaks using MACS3. 
[ENCODE standards](https://www.encodeproject.org/atac-seq/) require:
-**Replicated peak files** with >150,000 peaks (acceptable >100,000)
-**Nucleosome-free regions** must be detectable in called peaks

The output includes narrowPeak files, summit locations, and associated metrics (e.g., q-values).

#### **11. Lorenz Curves**
The `Lorenz_curve` rule generates Lorenz curves (or fingerprint plots) to assess library complexity and sequencing biases across samples.

![Lorenz-curve](resources/lorenz_curve_examples.svg)
***Figure 4: [Lorenz curves](https://sebastian-gregoricchio.github.io/snakeATAC/) comparing (left) oversequenced/low-complexity sample vs. (right) ideal complex library.***

#### **12. Merging Peaks Across Samples**
The `all_peaks_file_and_score_matrix` rule merges peak files across all samples to create a unified peak set. It also generates a score matrix summarizing signal intensity at each peak for all samples.

#### **13. Summary Table Generation**
Finally, the `counts_summary` rule compiles key metrics (e.g., number of mapped reads, peaks detected, FRiP scores) into a summary table for easy interpretation.

![DAG](resources/dag.svg)


## Output Directory Structure

The pipeline generates a well-organized output directory structure. 

### **Key Directories**
- **01_trimmed_fastq/**: Contains trimmed FASTQ files and logs from Cutadapt.
- **02_BAM/**: Stores filtered and deduplicated BAM files, and alignment statistics.
- **02_BAM_fastQC/**: Stores FastQC and MultiQC files for filtered and deduplicated BAM files.
- **03_Normalization/**: Contains RPM-normalized bigWig files for visualization in genome browsers.
- **04_MACS3_peaks/**: Includes peak-calling results such as narrowPeak files (chromatin accessibility regions) and summit files (precise peak summits).
- **05_quality_controls/**: Contains FastQC and MultiQC files for trimmed FASTQ files.
- **06_Overall_quality_and_info/**: Aggregates overall quality metrics like FRiP scores, MultiQC reports summarizing QC metrics across all samples, and Lorenz curves for library complexity.

Below is a detailed breakdown:

```
output_directory/
├── 01_trimmed_fastq/
│   ├── logs/
│   │   ├── cutadapt.sample1.err
│   │   ├── cutadapt.sample1.out
│   │   └── ...
│   ├── sample1_R1_trimmed.fastq.gz
│   ├── sample1_R2_trimmed.fastq.gz
│   └── ...
├── 02_BAM/
│   ├── sample1_mapq20_sorted_woMT_dedup.bai
│   ├── sample1_mapq20_sorted_woMT_dedup.bam
│   │   ├── sample1.bam
│   │   ├── sample2.bam
│   │   └── ...
│   ├── flagstat
│   │   ├── sample1_flagstat_filtered_bam_woMT.txt
│   │   └── sample1_flagstat_UNfiltered_bam.txt
│   ├── bwa_summary/
│   │   ├── sample1.BWA_summary.txt
│   │   ├── sample2.BWA_summary.txt
│   │   └── ...
│   ├── logs/
│   │   ├── sample1.sort.log
│   │   ├── sample2.sort.log
│   │   └── ...
│   ├── MarkDuplicates_metrics/
│   │   ├── sample1_MarkDuplicates_metrics.txt
│   │   ├── sample2_MarkDuplicates_metrics.txt
│   │   └── ...
│   ├── MarkDuplicates_logs/
│   │   ├── sample1_MarkDuplicates.out
│   │   ├── sample1_MarkDuplicates.err
│   │   ├── sample2_MarkDuplicates.out
│   │   ├── sample2_MarkDuplicates.err
│   │   └── ...
│   ├── reads_per_chromosome/
│   │   ├── sample1_idxstats_read_per_chromosome.txt
│   │   ├── sample2_idxstats_read_per_chromosome.txt
│   │   └── ...
│   └── FixMateInformation_logs/
│       ├── sample1_FixMateInformation.log
│       ├── sample2_FixMateInformation.log
│       └── ...
├── 02_BAM_fastQC/
│   ├── sample1_mapq20_sorted_woMT_dedup_fastqc.zip
│   └── sample1_mapq20_sorted_woMT_dedup_fastqc.html
├── 03_Normalization/RPM_normalized
│   ├── bamToBed_log/
│   │   ├── sample1_bamToBed.log
│   │   ├── sample2_bamToBed.log
│   │   └── ...
│   └── sample1_mapq20_woMT_dedup_shifted_RPM.normalized.bw
├── 04_MACS3_peaks/
│   ├── logs/
│   │   ├── sample1_mapq20_woMT_dedup_qValue0.05.log
│   │   ├── sample2_mapq20_woMT_dedup_qValue0.05.log
│   │   └── ...
│   ├── sample1_mapq20_woMT_dedup_qValue0.05_treat_pileup.bdg
│   ├── sample1_mapq20_woMT_dedup_qValue0.05_control_lambda.bdg
│   ├── sample1_mapq20_woMT_dedup_qValue0.05_peaks.xls
│   ├── sample1_mapq20_woMT_dedup_qValue0.05_peaks.narrowPeak
│   ├── sample1_mapq20_woMT_dedup_qValue0.05_summits.bed
│   └── sample1_mapq20_woMT_dedup_qValue0.05_peaks_chr.narrowPeak
├── 05_quality_controls/
│   ├── trimmed_fastq_multiQC/
│   │   ├── multiQC_report_trimmed_fastq.out
│   │   └── multiQC_report_trimmed_fastq.err
│   └── trimmed_fastq_fastqc/
│       ├── sample1_R1_trimmed_fastqc.zip
│       └── sample1_R1_trimmed_fastqc.html
└── 06_Overall_quality_and_info/
    ├── LorenzCurve_plotFingreprint/lorenz_plots/
    │   ├── sample1_Lorenz_curve_deeptools.plotFingreprint.pdf
    │   └── ...
    ├── Counts
    │   ├── counts_summary.txt
    │   └── subread_featureCounts_output
    │       └── sample1
    │           ├── sample1.readCountInPeaks
    │           ├── sample1.readCountInPeaks.log
    │           └── sample1.readCountInPeaks.summary
    └── Sample_comparisons
        ├── multiBigWigSummary_matrix_allSamples.npz
        ├── PCA_on_BigWigs_wholeGenome.pdf
        └── Peak_comparison
            ├── all_samples_peaks_concatenation_collapsed_sorted.bed
            ├── peaks_score_matrix_all_samples_MACS3.npz
            └── peaks_score_matrix_all_samples_table_MACS3.tsv
```

## Acknowledgments

A huge thank you to **Dr. Isabell Bludau**, **Dr. Paul Kerbs**, and **Quynh Nhu Nguyen** from Heidelberg University Hospital and the German Cancer Research Center (DKFZ) for their support, feedback, and contributions to this pipeline.

## References
1. snakeATAC: [https://sebastian-gregoricchio.github.io/snakeATAC/](https://sebastian-gregoricchio.github.io/snakeATAC/)
2. ENCODE ATAC-seq Guidelines: [https://www.encodeproject.org/atac-seq/](https://www.encodeproject.org/atac-seq/)
3. nf-core ATAC-seq Workflow: [https://nf-co.re/atacseq](https://nf-co.re/atacseq)
4. Best Practices for ATAC-seq Analysis: [Genome Biology](https://doi.org/10.1186/s13059-020-1929-3)
5. [Galaxy ATAC-seq Tutorial](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html)  
6. [Cutadapt Parameter Optimization (PMC10035359)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10035359/)  
7. [FastQC Figure Source](https://bioinformaticamente.com/2024/12/05/comprehensive-guide-to-atac-seq-data-quality-control/)  

## License

This project is licensed under the MIT License—see the [`LICENSE`](MIT_License.md) file for details.