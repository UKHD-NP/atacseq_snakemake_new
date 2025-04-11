# atacseq_snakemake
 ATAC-seq Snakemake workflow

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

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/atacseq-pipeline.git
   cd atacseq-pipeline
   ```
   or click on *Code > Download ZIP* on the GitHub page

2. **Install dependencies**:
   This pipeline uses Conda/Mamba for environment management. Install the required environment as follows:
   ```bash
   mamba env create -f workflow/envs/env.yaml
   ```
3. **Activate the environment**:
    ```bash
   mamba activate atacseq_pipeline
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
2. If no issues arise during the dry run, execute the full pipeline:
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

## Pipeline Overview

### Processing Steps
1. **Read Trimming**  
   Cutadapt removes adapters and low-quality bases (Q20%
- **Fragment distribution**: Peak ~200bp (nucleosomal pattern)
- **Alignment rate**: >80% mapped reads
- **Duplicate rate**: <50% (post-filtering)

## References
1. snakeATAC: [https://sebastian-gregoricchio.github.io/snakeATAC/](https://sebastian-gregoricchio.github.io/snakeATAC/)
2. ENCODE ATAC-seq Guidelines: [https://www.encodeproject.org/atac-seq/](https://www.encodeproject.org/atac-seq/)
3. nf-core ATAC-seq Workflow: [https://nf-co.re/atacseq](https://nf-co.re/atacseq)
4. Best Practices for ATAC-seq Analysis: [Genome Biology](https://doi.org/10.1186/s13059-020-1929-3)

## License

This project is licensed under the MIT License—see the `LICENSE` file for details.