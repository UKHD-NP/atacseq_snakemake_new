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


This is a modular Snakemake workflow for ATAC-seq data analysis, inspired by [snakeATAC](https://sebastian-gregoricchio.github.io/snakeATAC/) and best practices from ENCODE and nf-core ATAC-seq. It automates the analysis of ATAC-seq data from raw sequencing reads to processed results, including trimmed FASTQ files, aligned BAM files, peak calls, normalized coverage tracks, and comprehensive quality control metrics. 


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
    --cores <number_of_cores> \
    --use-conda \
    -s workflow/atacseq.snakefile \
    --configfile config/atacseq_config.yaml \
    --config \
    samplesheet="/path/to/samplesheet.csv" \
    output_directory="/path/to/output_directory" \
    genome_fasta="/path/to/genome.fa" \
    tss_bed="/path/to/tss.bed" \
    chromsizes="/path/to/chromsizes.txt" \
    blacklist="/path/to/blacklist.bed" \
    -n
    ```
2. If no errors arise during the dry run *(it does not mean that the pipeline will run without any issues)*, execute the full pipeline:
    ```bash
    snakemake \
    --cores <number_of_cores> \
    --use-conda \
    -s workflow/atacseq.snakefile \
    --configfile config/atacseq_config.yaml \
    --config \
    samplesheet="/path/to/samplesheet.csv" \
    output_directory="/path/to/output_directory" \
    genome_fasta="/path/to/genome.fa" \
    tss_bed="/path/to/tss.bed" \
    chromsizes="/path/to/chromsizes.txt" \
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
| **General Workflow**                   |                                                                                                                                    |                       |
| samplesheet                            | Path to samplesheet CSV (sample info and FASTQ paths). **Required**.                                                               | -                     |
| output_directory                       | Main output directory for results. **Required**.                                                                                   | -                     |
| genome_fasta                           | Reference genome FASTA file. **Required**.                                                                                         | -                     |
| tss_bed                                | BED file containing transcription start site (TSS) coordinates for TSS enrichment calculation. **Required**.                       | -                     |
| chromsizes                             | Tab-delimited file with chromosome names and sizes (chr\tsize). **Required**.                                                      | -                     |
| blacklist                              | ENCODE blacklist BED file to exclude problematic genomic regions. **Required**.                                                    | -                     |
| runs_directory                         | Directory where BAM files are stored. Auto-generated from output_directory.                                                        | "{{output_directory}}/02_BAM" |
| fastq_suffix                           | FASTQ file extension (e.g., .fastq.gz). Must match input files.                                                                    | ".fastq.gz"           |
| read_suffix                            | Read pair identifiers (e.g., _R1, _R2). Match your FASTQ naming convention.                                                       | ['_R1', '_R2']        |
| **Trimming**                           |                                                                                                                                    |                       |
| cutadapt_trimm_options                 | Additional custom trimming parameters for Cutadapt.                                                                                | ''                    |
| fw_adapter_sequence                    | Forward adapter sequence for Tn5 transposase trimming (Nextera adapters).                                                          | "CTGTCTCTTATACACATCT" |
| rv_adapter_sequence                    | Reverse adapter sequence for Tn5 transposase trimming (Nextera adapters).                                                          | "CTGTCTCTTATACACATCT" |
| run_fastq_qc                           | Whether to perform FastQC and MultiQC on trimmed FASTQ files.                                                                      | True                  |
| **BWA Mapping**                        |                                                                                                                                    |                       |
| bwa_options                            | Custom alignment options for BWA-mem2.                                                                                             | ''                    |
| **BAM Filtering**                      |                                                                                                                                    |                       |
| remove_duplicates                      | Remove PCR duplicates using GATK MarkDuplicates. Set to False to keep duplicates.                                                  | True                  |
| MAPQ_threshold                         | Minimum mapping quality score (MAPQ). Filters ambiguous alignments (recommended: 20). Increase for stricter alignment (range 0-60).| 20                    |
| remove_other_chromosomes_pattern       | Regex pattern to exclude unwanted contigs (e.g., mitochondrial DNA, viral sequences, random contigs).                              | "CMV\|HBV\|HTLV..."   |
| bam_features.bam_suffix                | Suffix for processed BAM files.                                                                                                    | "_mapq20_sorted_woMT_dedup.bam" |
| bam_features.minFragmentLength         | Minimum fragment size for paired-end reads (filters short fragments).                                                              | 0                     |
| bam_features.maxFragmentLength         | Maximum fragment size for paired-end reads (filters overly long fragments).                                                        | 2000                  |
| **Genome Annotations**                 |                                                                                                                                    |                       |
| genomic_annotations.genome_id          | Genome identifier (informational only, e.g., hg38, mm10).                                                                          | "hg38"                |
| genomic_annotations.effective_genomeSize | Effective genome size for MACS3 peak calling. Must match reference genome.                                                       | 2913022398            |
| genomic_annotations.ignore_for_normalization | Chromosomes excluded from coverage normalization (e.g., sex chromosomes, contigs, mitochondrial DNA).                        | "X Y MT M KI270..."   |
| **Peak Calling**                       |                                                                                                                                    |                       |
| peak_calling.qValue_cutoff             | MACS3 q-value significance threshold for peak calling (lower = more stringent).                                                    | 0.05                  |
| peak_calling.call_summits              | Whether to report precise peak summits (useful for motif analysis and TF binding site identification).                             | True                  |
| peak_calling.FRiP_threshold            | Minimum Fraction of Reads in Peaks (FRiP) score threshold in % (ENCODE: >30% ideal, >20% acceptable).                             | 20                    |
| **Quality Control**                    |                                                                                                                                    |                       |
| quality_controls.fragmentSize_window_length | Bin size (bp) for fragment size distribution plots (larger bins = smoother distributions).                                   | 1000                  |
| quality_controls.plotFingerprint.binSize | Resolution (bp) of fingerprint/Lorenz plots (higher = more detail, slower processing).                                         | 500                   |
| quality_controls.plotFingerprint.sampledRegions | Number of genomic regions randomly sampled for Lorenz curve calculation.                                           | 500000                |
| quality_controls.plotFingerprint.extra_parameters | Additional parameters for deepTools plotFingerprint.                                                            | ""                    |
| quality_controls.calculate_pbc         | Calculate PCR Bottlenecking Coefficient (PBC) and Non-Redundant Fraction (NRF) for library complexity assessment.                  | True                  |
| quality_controls.calculate_frip        | Calculate Fraction of Reads in Peaks (FRiP) score.                                                                                | True                  |
| quality_controls.calculate_tss         | Calculate TSS enrichment score for signal-to-noise assessment.                                                                     | True                  |


## Pipeline Overview


#### **1. Initialization**
The pipeline begins with an initialization rule (`AAA_initialization`), which ensures that all required outputs are properly organized and triggers the execution of all downstream rules.


#### **2. Genome Index Generation**
If a genome index is not already available, the rule `generate_genome_index` creates it using BWA-mem2 and samtools. This step ensures that the reference genome is prepared for efficient read alignment. The output includes `.bwt.2bit.64`, `.fai`, and other index files necessary for mapping.


#### **3. Adapter Trimming**
The `cutadapt_PE` rule trims sequencing adapters from paired-end reads using Cutadapt. 
The Cutadapt parameters in this ATAC-seq pipeline are optimized for typical Nextera-based library preparation and sequencing on NextSeq/NovaSeq platforms. The `--nextseq-trim=20` parameter addresses the two-color chemistry's tendency to produce poly-G artifacts by trimming 3' bases with quality scores below 20. The `--minimum-length 20` setting filters out reads shorter than 20 bp, removing adapter dimers and fragments too short for meaningful peak calling. Based on [this](https://pmc.ncbi.nlm.nih.gov/articles/PMC10035359/) and [this](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html) practical example, these values were chosen for optimal performance.


For Illumina HiSeq/MiSeq (non-NextSeq platforms), replace `--nextseq-trim` with generic quality trimming (e.g., `-q 20`), as poly-G artifacts are absent. If using TruSeq adapters, update `-a`/`-A` to the correct sequences in the configfile. For lower-quality data (e.g., degraded samples), relax `-q` to 15-20 to avoid excessive read loss. Always validate parameters using FastQC or MultiQC to ensure adapter removal and fragment size distribution align with expectations.


![Base Content in the sequences](resources/fastqc_base_content.png)
***Figure 1: Expected [FastQC base content profiles](https://bioinformaticamente.com/2024/12/05/comprehensive-guide-to-atac-seq-data-quality-control/) for (A) raw and (B) trimmed reads. Post-trimming (B) shows reduced 3'-end G-content spikes caused by NextSeq poly-G artifacts.***


![Adapter content in the sequences](resources/fastqc_adapter_content.png)
***Figure 2: [Adapter content profiles](https://bioinformaticamente.com/2024/12/05/comprehensive-guide-to-atac-seq-data-quality-control/) for (A) raw and (B) trimmed reads. Successful adapter removal eliminates >99% of Nextera sequences.***


#### **4. Read Alignment**
The `BWA_PE` rule aligns the trimmed reads to the reference genome using BWA-mem2. This aligner was chosen for its widespread use in the literature and its optimal balance of accuracy and speed demonstrated in different [benchmarks](https://www.nature.com/articles/s41467-021-26865-w). The alignment includes proper mate-fixing with `samtools fixmate` to ensure correct mate pair information.


#### **5. BAM Filtering**
The `MAPQ_MT_filter` rule filters aligned reads based on mapping quality (MAPQ) and removes unwanted chromosomes (e.g., mitochondrial DNA). Based on [ENCODE guidelines](https://www.encodeproject.org/atac-seq/), each replicate should retain **≥50 million non-duplicate, non-mitochondrial aligned reads** for paired-end analysis. Additionally, duplicate reads are marked or removed using GATK's `MarkDuplicatesWithMateCigar` in the subsequent `gatk4_markdups` rule.


**How BAM files and reads are filtered throughout the pipeline:**
- **Remove mitochondrial DNA reads**: Reads mapping to `chrM` or mitochondrial contigs are excluded using samtools.
- **Exclude blacklisted regions**: Reads overlapping ENCODE blacklist regions are removed using BEDtools.
- **Remove duplicate reads**: PCR duplicates are marked and optionally removed using GATK `MarkDuplicatesWithMateCigar`.
- **Filter by mapping quality (MAPQ)**: Reads with MAPQ scores below the threshold (default: 20, set in the configfile) are discarded using samtools.
- **Exclude secondary alignments**: Only primary alignments (`-F 0x100`) are retained using samtools.
- **Remove unmapped reads**: Reads flagged as unmapped (`-F 0x4`) are filtered out using samtools.
- **Exclude multimapped reads**: Reads mapping to multiple locations are removed based on MAPQ scoring.
- **Filter by mismatches**: Reads with more than 4 mismatches (using the NM tag) are excluded during filtering using samtools.
- **Filter by fragment size**: Fragments outside the range of 0–2000 bp are removed using BEDtools and custom scripts.
- **Exclude improperly paired reads**: Only properly paired reads (`-f 3`) are retained for downstream analysis.
- **Handle paired-end read inconsistencies**: Reads where only one mate fails any of the above criteria are excluded during filtering.


#### **6. Quality Control on Trimmed Reads**
The `fastQC_trimmed_fastq` rule performs quality control checks on trimmed FASTQ files using FastQC. A MultiQC report (`multiQC_trimmed_fastq`) aggregates these QC metrics into a single HTML file for easy review.


#### **7. Quality Control on BAM Files**
The `fastQC_BAMs` rule assesses the quality of BAM files. A MultiQC report (`multiQC_BAMs`) consolidates these metrics alongside alignment statistics and peak-calling results.


**Key metrics from [ENCODE](https://www.encodeproject.org/atac-seq/) include:**  
- **Alignment rate**: >95% (acceptable >80%) 
- **FRiP score**: >0.3 (acceptable >0.2)  
- **TSS enrichment**: >10 (acceptable >5)
- **Library complexity** (NRF): >0.9
- **PBC1**: >0.9
- **PBC2**: >3
- **Nucleosome-free regions**: Must be detectable in called peaks 


#### **8. Fragment Size Distribution**
The `fragment_size_distribution` rule calculates fragment size distributions for each sample, which provides insights into nucleosome positioning and library complexity. A combined plot (`fragment_size_distribution_report`) summarizes these distributions across all samples using both the original deepTools plots and ggplot2 visualizations.


**High-quality ATAC-seq data based on [ENCODE](https://www.encodeproject.org/atac-seq/) must show:**
- **Nucleosome-free region (NFR) peak** at ~50 bp 
- **Mononucleosome peak** between 147–294 bp
- Clear separation of di-/tri-nucleosome peaks (desirable but not always required)


![Fragment-distribution](resources/fragmentSize_distribution_examples.svg)
***Figure 3: [Fragment distributions](https://sebastian-gregoricchio.github.io/snakeATAC/) showing (left) ideal profile with NFR and mono-/di-nucleosome peaks vs. (right) noisy data lacking clear nucleosomal patterning.***


#### **9. TSS Enrichment Score Calculation**
The `tss_plot` rule calculates the **Transcription Start Site (TSS) enrichment score**, which serves as a critical signal-to-noise quality control metric for ATAC-seq data.


**TSS Enrichment Score Definition (from [ENCODE](https://www.encodeproject.org/atac-seq/)):**

The TSS enrichment calculation is a signal-to-noise metric (based on [ATACseqQC](https://doi.org/10.1186/s12864-018-4559-3)). Reads around a reference set of TSSs are aggregated to form a distribution centered on the TSSs and extending to ±2000 bp (for a total of 4000 bp). This distribution is then normalized by taking the average read depth in the 100 bp at each of the end flanks (for a total of 200 bp of averaged data) and calculating a fold change at each position over that average read depth. 

This means the flanks should start at 1, and if there is high read signal at transcription start sites (highly accessible chromatin), there should be an increase in signal up to a peak in the middle. The signal value at the center of the distribution after normalization is taken as the TSS enrichment score.

**TSS Enrichment Score = (Average reads at TSS ± 50 bp) / (Average reads at TSS flanks ± 1900-2000 bp)**

**Quality thresholds:**
- **Excellent**: TSS enrichment >10
- **Acceptable**: TSS enrichment >5
- **Poor**: TSS enrichment <5

The pipeline outputs both a standard plot and a large version, along with the numerical TSS enrichment score for each sample.


#### **10. Library Complexity Assessment (PBC/NRF Metrics)**
The `calculate_pbc` rule computes **PCR Bottlenecking Coefficients (PBC)** and the **Non-Redundant Fraction (NRF)** to assess library complexity and determine whether additional sequencing would be beneficial.


**Metrics calculated:**
- **M_TOTAL**: Total number of mapped read pairs
- **M_DISTINCT**: Number of distinct uniquely mapped read pairs (genomically unique locations)
- **M1**: Number of genomically unique locations with exactly 1 read pair
- **M2**: Number of genomically unique locations with exactly 2 read pairs

**Formulas:**
- **NRF (Non-Redundant Fraction)** = M_DISTINCT / M_TOTAL
- **PBC1** = M1 / M_DISTINCT  
- **PBC2** = M1 / M2

**Quality thresholds from [ENCODE](https://www.encodeproject.org/atac-seq/):**
- **NRF**: >0.9 (ideal), >0.8 (acceptable), <0.7 (concerning - may indicate PCR over-amplification)
- **PBC1**: >0.9 (ideal), >0.8 (acceptable), <0.7 (concerning)
- **PBC2**: >3 (ideal), >1 (acceptable), <1 (concerning)

Low PBC/NRF values indicate library complexity issues, suggesting that most reads are PCR duplicates rather than unique genomic fragments. This can occur due to low input material, excessive PCR cycles, or poor transposition efficiency.

The `merge_pbc_results` rule combines PBC metrics from all samples into a single summary table for easy comparison.


#### **11. Read Shifting and Normalization**
The `bam_shifting_and_RPM_normalization` rule shifts paired-end reads to account for Tn5 transposase binding offset (+4 bp for forward strand, -5 bp for reverse strand) and generates RPM-normalized bigWig files for visualization in genome browsers. This shifting adjusts read positions to represent the exact Tn5 cut sites rather than the sequenced fragment ends.


#### **12. Peak Calling**
The `peakCalling_MACS3` rule identifies regions of open chromatin by calling peaks using MACS3. 


**[ENCODE standards](https://www.encodeproject.org/atac-seq/) require:**
- **Replicated peak files** with >150,000 peaks (acceptable >100,000)
- **IDR peak files** with >70,000 peaks (acceptable >50,000)
- **Nucleosome-free regions** must be detectable in called peaks

The peak calling uses the `--nomodel --shift -75 --extsize 150` parameters optimized for ATAC-seq cut-site analysis, focusing on the nucleosome-free regions. The output includes narrowPeak files, summit locations, and associated metrics (e.g., q-values). Peaks overlapping blacklist regions are automatically removed.


#### **13. Lorenz Curves (Fingerprint Plots)**
The `Lorenz_curve` rule generates Lorenz curves (or fingerprint plots) to assess library complexity and sequencing biases across samples. These plots show the cumulative fraction of reads as a function of the cumulative fraction of the genome, ranked by coverage depth.


**Interpretation:**
- **High-quality libraries**: Show a steep curve indicating that most reads are concentrated in a small fraction of the genome (open chromatin regions)
- **Low-complexity/over-sequenced libraries**: Show a diagonal curve approaching the theoretical maximum, indicating saturation
- **Poor libraries**: Show curves close to the diagonal, indicating random genomic distribution rather than enrichment

![Lorenz-curve](resources/lorenz_curve_examples.svg)
***Figure 4: [Lorenz curves](https://sebastian-gregoricchio.github.io/snakeATAC/) comparing (left) over-sequenced/low-complexity sample vs. (right) ideal complex library.***

The `Lorenz_curve_merge_plots` rule combines all individual Lorenz plots into a single PDF and creates a combined ggplot2 visualization for easy multi-sample comparison.


#### **14. Summary Table Generation**
Finally, the `counts_summary` rule compiles key metrics (e.g., number of mapped reads, peaks detected, FRiP scores) into a summary table for easy interpretation. This table includes:
- Number of deduplicated reads
- Number of shifted reads used for peak calling
- Number of peaks called
- FRiP percentage
- Quality label (good/bad) based on FRiP threshold

This provides an at-a-glance overview of the quality and results for all samples in the experiment.


![DAG](resources/dag.svg)


## Output Directory Structure


The pipeline generates a well-organized output directory structure. 


### **Key Directories**
- **01_trimmed_fastq/**: Contains trimmed FASTQ files and logs from Cutadapt.
- **02_BAM/**: Stores filtered and deduplicated BAM files, alignment statistics, and MarkDuplicates metrics.
- **02_BAM_fastQC/**: Stores FastQC and MultiQC files for filtered and deduplicated BAM files.
- **03_Normalization/**: Contains RPM-normalized bigWig files for visualization in genome browsers and shifted BED files.
- **04_MACS3_peaks/**: Includes peak-calling results such as narrowPeak files (chromatin accessibility regions), summit files (precise peak summits), and bedGraph files.
- **05_quality_controls/**: Contains FastQC and MultiQC reports for trimmed FASTQ files.
- **06_Overall_quality_and_info/**: Aggregates overall quality metrics including:
  - Fragment size distribution plots and metrics
  - TSS enrichment plots and scores
  - PBC/NRF library complexity metrics
  - Lorenz curves for library complexity assessment
  - MultiQC reports summarizing QC metrics across all samples
  - FRiP scores and counts summary table
  - Sample comparison files (merged peaks, score matrices)


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
│   ├── flagstat/
│   │   ├── sample1_mapq20_sorted_woMT_dedup_flagstat.txt
│   │   └── ...
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
│   ├── sample1_mapq20_sorted_woMT_dedup_fastqc.html
│   └── ...
├── 03_Normalization/RPM_normalized/
│   ├── bamToBed_log/
│   │   ├── sample1_bamToBed.log
│   │   ├── sample2_bamToBed.log
│   │   └── ...
│   ├── sample1_mapq20_sorted_woMT_dedup_read_sortedByPos_shifted.bed
│   ├── sample1_mapq20_woMT_dedup_shifted_RPM.normalized.bw
│   └── ...
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
│   ├── sample1_mapq20_woMT_dedup_qValue0.05_peaks_chr.narrowPeak
│   └── ...
├── 05_quality_controls/
│   ├── trimmed_fastq_multiQC/
│   │   ├── multiQC_report_trimmed_fastq.html
│   │   ├── multiQC_report_trimmed_fastq.out
│   │   └── multiQC_report_trimmed_fastq.err
│   └── trimmed_fastq_fastqc/
│       ├── sample1_R1_trimmed_fastqc.zip
│       ├── sample1_R1_trimmed_fastqc.html
│       ├── sample1_R2_trimmed_fastqc.zip
│       ├── sample1_R2_trimmed_fastqc.html
│       └── ...
└── 06_Overall_quality_and_info/
    ├── fragmentSizeDistribution_plots/
    │   ├── sample1_fragment_size_distribution.pdf
    │   ├── ALL.SAMPLES_fragmentSizeDistribution_plots.pdf
    │   ├── ALL.SAMPLES_fragmentSizeDistribution_plots_ggplot.version.pdf
    │   ├── log/
    │   │   ├── sample1_fragmentSize_log.out
    │   │   ├── sample1_fragmentSize_log.err
    │   │   ├── ggplot_replotting.log
    │   │   └── pdfcombine.log
    │   └── table_and_fragmentSize/
    │       ├── sample1_fragmentSize_RawFragmentLengths.txt
    │       ├── sample1_fragmentSize_metrics.txt
    │       └── ...
    ├── TSS_plots/
    │   ├── sample1_sample_tss-enrich.pdf
    │   ├── sample1_sample_large_tss-enrich.pdf
    │   ├── sample1_sample_tss-enrich.score.txt
    │   ├── log/
    │   │   ├── tss_plot_sample1.out
    │   │   ├── tss_plot_sample1.err
    │   │   └── ...
    │   └── ...
    ├── LorenzCurve_plotFingreprint/
    │   ├── lorenz_plots/
    │   │   ├── sample1_Lorenz_curve_deeptools.plotFingreprint.pdf
    │   │   └── ...
    │   ├── lorenz_metrics/
    │   │   ├── sample1_Lorenz_quality.metrics_deeptools.plotFingreprint.txt
    │   │   └── ...
    │   ├── lorenz_counts/
    │   │   ├── sample1_Lorenz_raw.counts_deeptools.plotFingreprint.txt
    │   │   └── ...
    │   ├── log/
    │   │   ├── sample1_deeptools_plotFingreprint.log
    │   │   ├── LorenzCurve_plotFingreprint_mergePdf.log
    │   │   └── LorenzCurve_plotFingreprint_ggplot.merge.plots.log
    │   ├── Lorenz_curve_deeptools.plotFingreprint_allSamples_combined.pdf
    │   └── Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf
    ├── Counts/
    │   ├── counts_summary.txt
    │   └── subread_featureCounts_output/
    │       └── sample1/
    │           ├── sample1.readCountInPeaks
    │           ├── sample1.readCountInPeaks.log
    │           └── sample1.readCountInPeaks.summary
    ├── PBC/
    │   ├── sample1_pbc.qc
    │   ├── sample2_pbc.qc
    │   ├── merged_pbc_metrics.tsv
    │   └── ...
    ├── multiQC_dedup_bams/
    │   ├── multiQC_report_BAMs_dedup.html
    │   ├── multiQC_report_BAMs_dedup.out
    │   └── multiQC_report_BAMs_dedup.err
    └── Sample_comparisons/
        └── Peak_comparison/
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
7. [FastQC Figure Source](https://bioinformaticamente.com/2024/12/05/comprehensive-guide-to-atac-seq-data-quality -control/)
8. [PEPATAC: Optimized Pipeline for ATAC-seq](https://doi.org/10.1093/nargab/lqab101)
9. [AIAP: ATAC-seq Integrative Analysis Package](https://doi.org/10.1016/j.gpb.2021.06.001)
10. [ataqv: Quantification and Visualization of ATAC-seq Bias](https://doi.org/10.1016/j.cels.2020.02.009)
11. [ATACseqQC: a Bioconductor package for post-alignment quality assessment of ATAC-seq data.](https://doi.org/10.1186/s12864-018-4559-3)



## License


This project is licensed under the MIT License—see the [`LICENSE`](MIT_License.md) file for details.
