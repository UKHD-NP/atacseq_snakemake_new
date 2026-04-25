suppressPackageStartupMessages({
    library(ATACseqQC)
    library(rtracklayer)
    library(GenomicAlignments)
})

bam_file    <- snakemake@input[["bam"]]
bed_file    <- snakemake@input[["bed"]]
output_file <- snakemake@output[["pt_score"]]
log_file    <- snakemake@log[[1]]
sample_id   <- snakemake@wildcards[["sample_id"]]

log_con <- file(log_file, "w")

tryCatch({
    # Read transcript annotations
    txs <- import(bed_file, format = "BED")

    # Read BAM file
    gal <- readBamFile(bamFile = bam_file, asMates = FALSE)

    # Calculate PT score (ATACseqQC implementation)
    pt <- PTscore(gal, txs)

    mean_pt   <- mean(pt$PT_score,       na.rm = TRUE)
    median_pt <- median(pt$PT_score,     na.rm = TRUE)
    mean_pro  <- mean(pt$promoter,       na.rm = TRUE)
    mean_body <- mean(pt$transcriptBody, na.rm = TRUE)

    df <- data.frame(
        Sample          = sample_id,
        PT_score_mean   = round(mean_pt,   4),
        PT_score_median = round(median_pt, 4),
        Mean_promoter   = round(mean_pro,  6),
        Mean_gene_body  = round(mean_body, 6)
    )
    write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

    qc_msg <- if (!is.nan(mean_pt) && 2^mean_pt >= 5)
        "[INFO] QC: PASS (promoter/body ratio >= 5)"
    else
        "[WARNING] QC: FAIL (promoter/body ratio < 5, expected >= 5-10)"

    writeLines(c(
        sprintf("[INFO] Transcripts processed: %d",            length(pt)),
        sprintf("[INFO] Mean PT score (log2 ratio): %.4f",     mean_pt),
        sprintf("[INFO] Median PT score (log2 ratio): %.4f",   median_pt),
        sprintf("[INFO] Equivalent mean ratio (2^PT): %.2f",   2^mean_pt),
        qc_msg
    ), log_con)

}, error = function(e) {
    writeLines(sprintf("[ERROR] %s", conditionMessage(e)), log_con)
    stop(e)
}, finally = {
    close(log_con)
})
