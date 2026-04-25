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

    # Read BAM file (shared across all three metrics)
    gal <- readBamFile(bamFile = bam_file, asMates = FALSE)

    # PT score: log2(promoter density) - log2(body density), per transcript
    # promoter_window = [TSS-2000, TSS+500]; body_window = next 2500 bp downstream
    pt <- PTscore(gal, txs)

    mean_pt   <- mean(pt$PT_score,       na.rm = TRUE)
    median_pt <- median(pt$PT_score,     na.rm = TRUE)
    mean_pro  <- mean(pt$promoter,       na.rm = TRUE)
    mean_body <- mean(pt$transcriptBody, na.rm = TRUE)

    # NFR score: log2(nf) + 1 - log2(n1 + n2), per TSS (400 bp window)
    # nf = middle 100 bp; n1/n2 = flanking 150 bp nucleosome positions
    nfr <- NFRscore(gal, txs)

    mean_nfr   <- mean(nfr$NFR_score,   na.rm = TRUE)
    median_nfr <- median(nfr$NFR_score, na.rm = TRUE)

    # TSSE score: max(LOESS-smoothed mean enrichment in sliding windows ±1000 bp of TSS)
    # normalised to depth of end flanks (100 bp each side)
    # equivalent to the ENCODE TSS enrichment score definition
    tsse <- TSSEscore(gal, txs)
    tsse_score <- tsse$TSSEscore

    df <- data.frame(
        Sample           = sample_id,
        PT_score_mean    = round(mean_pt,    4),
        PT_score_median  = round(median_pt,  4),
        Mean_promoter    = round(mean_pro,   6),
        Mean_gene_body   = round(mean_body,  6),
        NFR_score_mean   = round(mean_nfr,   4),
        NFR_score_median = round(median_nfr, 4),
        TSSE_score       = round(tsse_score, 4)
    )
    write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

    qc_msg <- if (!is.nan(mean_pt) && 2^mean_pt >= 5)
        "[INFO] QC PT: PASS (promoter/body ratio >= 5)"
    else
        "[WARNING] QC PT: FAIL (promoter/body ratio < 5, expected >= 5-10)"

    writeLines(c(
        sprintf("[INFO] Transcripts processed (PT): %d",          length(pt)),
        sprintf("[INFO] Mean PT score (log2 ratio): %.4f",        mean_pt),
        sprintf("[INFO] Median PT score (log2 ratio): %.4f",      median_pt),
        sprintf("[INFO] Equivalent mean ratio (2^PT): %.2f",      2^mean_pt),
        qc_msg,
        sprintf("[INFO] TSS windows processed (NFR): %d",         length(nfr)),
        sprintf("[INFO] Mean NFR score (log2 scale): %.4f",       mean_nfr),
        sprintf("[INFO] Median NFR score (log2 scale): %.4f",     median_nfr),
        sprintf("[INFO] TSSE score: %.4f",                        tsse_score)
    ), log_con)

}, error = function(e) {
    writeLines(sprintf("[ERROR] %s", conditionMessage(e)), log_con)
    stop(e)
}, finally = {
    close(log_con)
})
