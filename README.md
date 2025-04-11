# ATAC-seq Snakemake workflow 

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

# ATAC-seq Snakemake Pipeline

A modular Snakemake workflow for ATAC-seq data analysis, inspired by [snakeATAC](https://sebastian-gregoricchio.github.io/snakeATAC/) and best practices from ENCODE and nf-core ATAC-seq.

### Key Features
1. **Adapter Trimming**: Removes sequencing adapters using Cutadapt.
2. **Alignment**: Maps reads to a reference genome using BWA-MEM2.
3. **BAM Processing**: Filters reads based on quality, removes duplicates, and excludes unwanted chromosomes.
4. **Peak Calling**: Identifies accessible chromatin regions using MACS3.
5. **Normalization**: Generates RPM-normalized coverage tracks in bigWig format.
6. **Quality Control**: Provides FRiP scores, fragment size distributions, and alignment statistics.
7. **Comprehensive Reporting**: MultiQC consolidates all QC metrics into a single report.

## Installation

1. **Place yourself in the directory where the repository should be downloaded**:
    ```bash
    cd </target/folder>
    ```

2. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/atacseq-pipeline.git
   cd atacseq_snakemake
   ```
   or click on *Code > Download ZIP* on the GitHub page

3. **Install dependencies**:
   This pipeline uses Conda/Mamba for environment management. Install the required environment as follows:
   ```bash
   mamba env create -f workflow/envs/mamba_env.yaml
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


### Parameters in configfile

| Section            | Parameter                              | Description                  | Default Value         | Usage Notes                                              |
|--------------------|----------------------------------------|------------------------------|-----------------------|----------------------------------------------------------|
| General Workflow   | samplesheet                            | Path to samplesheet CSV      | -                     | Required. Contains sample/replicate info and FASTQ paths |
|                    | output_directory                       | Main output directory        | -                     | Required. All results will be written here               |
|                    | genome_fasta                           | Reference genome FASTA       | -                     | Required. Should include .fai index                      |
|                    | blacklist                              | ENCODE blacklist BED         | -                     | Required. Removes problematic genomic regions            |
|                    | fastq_suffix                           | FASTQ file extension         | ".fastq.gz"           | Match your input file format                             |
|                    | read_suffix                            | Read pair identifiers        | ['_R1', '_R2']        | Must match FASTQ naming convention                       |
| Trimming           | cutadapt_trimm_options                 | Additional Cutadapt options  | ''                    | Add custom trimming parameters                           |
|                    | fw_adapter_sequence                    | Forward adapter              | "CTGTCTCTTATACACATCT" | Tn5 transposase adapter sequence                         |
|                    | rv_adapter_sequence                    | Reverse adapter              | "CTGTCTCTTATACACATCT" | Tn5 transposase adapter sequence                         |
| BWA Mapping        | bwa_options                            | BWA-MEM2 parameters          | ''                    | Add custom alignment options                             |
| BAM Filtering      | remove_duplicates                      | Remove PCR duplicates        | True                  | Set False to keep duplicates                             |
|                    | MAPQ_threshold                         | Minimum mapping quality      | 20                    | Increase for stricter alignment (range 0-60)             |
|                    | remove_other_chromosomes_pattern       | Excluded contigs             | "CMV\|HBV\|..."       | Regex pattern for unwanted sequences                     |
|                    | bam_features.minFragmentLength         | Minimum fragment size        | 0                     | Paired-end only                                          |
|                    | bam_features.maxFragmentLength         | Maximum fragment size        | 2000                  | Typical ATAC-seq cutoff                                  |
| Genome Annotations | genome_id                              | Genome identifier            | "hg38"                | Informational only                                       |
|                    | effective_genomeSize                   | Effective genome size        | 2913022398            | Critical for MACS3. Must match reference genome          |
|                    | ignore_for_normalization               | Excluded chromosomes         | "X Y MT M..."         | Contigs excluded from coverage normalization             |
| Peak Calling       | qValue_cutoff                          | MACS3 significance threshold | 0.05                  | Lower = more stringent                                   |
|                    | call_summits                           | Report peak summits          | True                  | Required for precise motif analysis                      |
|                    | FRiP_threshold                         | Minimum FRiP score           | 20                    | % reads in peaks (QC threshold)                          |
| Quality Control    | fragmentSize_window_length             | Fragment size bin            | 1000bp                | Larger bins smooth distributions                         |
|                    | multiBigwigSummary_binning_window_size | Correlation bin              | 10000bp               | Affects sample similarity metrics                        |
|                    | plotFingerprint.binSize                | Fingerprint resolution       | 500bp                 | Balance detail vs compute time                           |

## Pipeline Overview

#### **1. Initialization**
The pipeline begins with an initialization rule (`AAA_initialization`), which ensures that all required outputs are properly organized.

#### **2. Genome Index Generation**
If a genome index is not already available, the rule `generate_genome_index` creates it using BWA-MEM2 and with samtools. This step ensures that the reference genome is prepared for efficient read alignment. The output includes `.bwt`, `.fai`, and other index files necessary for mapping.

#### **3. Adapter Trimming**
The `cutadapt_PE` rule trims sequencing adapters from paired-end reads using Cutadapt. 
The Cutadapt parameters in this ATAC-seq pipeline are optimized for typical Nextera-based library preparation and sequencing on NextSeq/NovaSeq platforms. The `--nextseq-trim=20` parameter addresses the two-color chemistry’s tendency to produce poly-G artifacts by trimming 3’ bases with quality scores below 20. The `-q 27` threshold ensures only high-quality bases remain. The `--minimum-length 50` filters out reads shorter than 50 bp, removing adapter dimers and fragments too short for meaningful peak calling. The adapter sequences (set in the configfile) match the Illumina adapter sequences, ensuring proper removal of ligated adapters.

For **shorter read lengths**, `--minimum-length` should be reduced (e.g., 25–30 bp) to retain usable fragments. If switching to **single-end sequencing**, the `-A` parameter becomes obsolete. For **Illumina HiSeq/MiSeq** (non-NextSeq), replace `--nextseq-trim` with generic quality trimming (e.g., `-q 20`), as poly-G artifacts are absent. If using **TruSeq adapters**, update `-a`/`-A` to the correct sequences in the configfile. For **lower-quality data** (e.g., degraded samples), relax `-q` to 20–25 to avoid excessive read loss. Always validate parameters using FastQC or MultiQC to ensure adapter removal and fragment size distribution align with expectations.

#### **4. Read Alignment**
The `BWA_PE` rule aligns the trimmed reads to the reference genome using BWA-MEM2. This aligner was chosen for its widespread use in the literature and its optimal balance of accuracy and speed showed in different [benchmarks](https://www.nature.com/articles/s41467-021-26865-w).

#### **5. BAM Filtering**
The `MAPQ_MT_filter` rule filters aligned reads based on mapping quality (MAPQ) and removes unwanted chromosomes (e.g., mitochondrial DNA). Additionally, duplicate reads are marked or removed using GATK's `MarkDuplicatesWithMateCigar` in the subsequent `gatk4_markdups` rule.

How BAM files and reads are filtered throughout the pipeline:
1. **Remove mitochondrial DNA reads**  
   - Reads mapping to `chrM` or mitochondrial contigs are excluded using SAMtools.
2. **Exclude blacklisted regions**  
   - Reads overlapping ENCODE blacklist regions are removed using BEDtools.
3. **Remove duplicate reads**  
   - PCR duplicates are marked and optionally removed using GATK `MarkDuplicatesWithMateCigar`.
4. **Filter by mapping quality (MAPQ)**  
   - Reads with MAPQ scores below the threshold (default: 20, set in the configfile) are discarded using SAMtools.
5. **Exclude secondary alignments**  
   - Only primary alignments (`-F 0x100`) are retained using SAMtools.
6. **Remove unmapped reads**  
   - Reads flagged as unmapped (`-F 0x4`) are filtered out using SAMtools.
7. **Exclude multimapped reads**  
   - Reads mapping to multiple locations are removed based on MAPQ scoring.
8. **Filter by mismatches**  
   - Reads with more than 4 mismatches (using the NM tag) are excluded during filtering.
9. **Remove soft-clipped reads**  
   - Soft-clipped reads are implicitly excluded during alignment and filtering steps.
10. **Filter by fragment size**  
    - Fragments outside the range of 0–2000 bp are removed using BEDtools and custom scripts.
11. **Exclude improperly paired reads**  
    - Only properly paired reads (`-f 3`) are retained for downstream analysis.
12. **Handle paired-end read inconsistencies**  
    - Reads where only one mate fails any of the above criteria are excluded during filtering.

![DAG](https://github.com/UKHD-NP/atacseq_snakemake/dag.png)

#### **6. Quality Control on Trimmed Reads**
The `fastQC_trimmed_fastq` rule performs quality control checks on trimmed FASTQ files using FastQC. A MultiQC report (`multiQC_trimmed_fastq`) aggregates these QC metrics into a single HTML file for easy review.

#### **7. Quality Control on BAM Files**
The `fastQC_BAMs` rule assesses the quality of BAM files (e.g., duplication rates, coverage). A MultiQC report (`multiQC_BAMs`) consolidates these metrics alongside alignment statistics and peak-calling results.

#### **8. Fragment Size Distribution**
The `fragment_size_distribution` rule calculates fragment size distributions for each sample, which provides insights into nucleosome positioning and library complexity. A combined plot (`fragment_size_distribution_report`) summarizes these distributions across all samples.

#### **9. Read Shifting and Normalization**
The `bam_shifting_and_RPM_normalization` rule shifts paired-end reads to account for Tn5 transposase binding offset and generates RPM-normalized bigWig files for visualization in genome browsers.

#### **10. Peak Calling**
The `peakCalling_MACS3` rule identifies regions of open chromatin by calling peaks using MACS3. The output includes narrowPeak files, summit locations, and associated metrics (e.g., q-values).

#### **11. Lorenz Curves**
The `Lorenz_curve` rule generates Lorenz curves (or fingerprint plots) to assess library complexity and sequencing biases across samples.

#### **12. Merging Peaks Across Samples**
The `all_peaks_file_and_score_matrix` rule merges peak files across all samples to create a unified peak set. It also generates a score matrix summarizing signal intensity at each peak for all samples.

#### **13. Summary Table Generation**
Finally, the `counts_summary` rule compiles key metrics (e.g., number of mapped reads, peaks detected, FRiP scores) into a summary table for easy interpretation.

## References
1. snakeATAC: [https://sebastian-gregoricchio.github.io/snakeATAC/](https://sebastian-gregoricchio.github.io/snakeATAC/)
2. ENCODE ATAC-seq Guidelines: [https://www.encodeproject.org/atac-seq/](https://www.encodeproject.org/atac-seq/)
3. nf-core ATAC-seq Workflow: [https://nf-co.re/atacseq](https://nf-co.re/atacseq)
4. Best Practices for ATAC-seq Analysis: [Genome Biology](https://doi.org/10.1186/s13059-020-1929-3)

## License

This project is licensed under the MIT License—see the [`LICENSE`](https://github.com/UKHD-NP/atacseq_snakemake/MIT_Licence.md) file for details.