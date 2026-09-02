suppressPackageStartupMessages({
    library(ATACseqQC)
    library(rtracklayer)
    library(GenomicAlignments)
    library(Rsamtools)
    library(GenomeInfoDb)
})

std_chr <- paste0("chr", c(1:22, "X", "Y"))

bam_file      <- snakemake@input[["bam"]]
bed_file      <- snakemake@input[["bed"]]
output_file   <- snakemake@output[["score_tsv"]]
plot_fragsize <- snakemake@output[["plot_fragsize"]]
plot_pt       <- snakemake@output[["plot_pt"]]
plot_nfr      <- snakemake@output[["plot_nfr"]]
plot_tsse     <- snakemake@output[["plot_tsse"]]
log_file      <- snakemake@log[[1]]
sample_id     <- snakemake@wildcards[["sample_id"]]

log_con <- file(log_file, "w")

tryCatch({
    # Fragment size distribution
    png(plot_fragsize, width = 900, height = 600, res = 120)
    fragSizeDist(bam_file, sample_id)
    dev.off()

    # Read transcript annotations, restricted to standard chromosomes
    # (BED spans 505 contigs incl. alt/patch/scaffold; those carry no reads
    # and just bloat every downstream seqinfo/coverage object)
    txs <- import(bed_file, format = "BED")
    txs <- keepSeqlevels(txs, std_chr, pruning.mode = "coarse")

    # Read BAM file (shared across all three metrics), scoped to standard
    # chromosomes only. The BAM header carries ~129k decoy/sponge contigs
    # from the reference used for alignment; scanning restricted to the
    # 24 real chromosomes avoids inheriting that bloat into gal's seqinfo.
    bam_targets <- scanBamHeader(bam_file)[[1]]$targets[std_chr]
    which_gr <- GRanges(names(bam_targets), IRanges(1, bam_targets))

    # PTscore/NFRscore/TSSEscore only use read positions to build coverage,
    # so read a bare GAlignments (no seq/qual/tags). readBamFile() would pull
    # "seq" and "qual" for every read across all 24 chromosomes, which for a
    # deep ATAC BAM is >100 GB of base/quality strings and OOM-kills the job.
    # BAM is already filtered upstream (no secondary/supplementary/unmapped),
    # so no flag filter here.
    bam_param <- ScanBamParam(which = which_gr, what = character(0))
    gal <- readGAlignments(bam_file, param = bam_param)
    gal <- keepSeqlevels(gal, std_chr, pruning.mode = "coarse")

    # TSSE score: max(LOESS-smoothed mean enrichment in sliding windows ±1000 bp of TSS)
    # normalised to depth of end flanks (100 bp each side)
    # equivalent to the ENCODE TSS enrichment score definition
    tsse <- TSSEscore(gal, txs)
    tsse_score <- tsse$TSSEscore

    png(plot_tsse, width = 800, height = 600, res = 120)
    plot(100 * (-9:10 - 0.5), tsse$values, type = "b",
         xlab = "Distance to TSS (bp)",
         ylab = "Aggregate TSS score",
         main = paste(sample_id, "- TSS Enrichment Score:", round(tsse_score, 2)))
    abline(v = 0,  lty = 2, col = "gray40")
    abline(h = 1,  lty = 2, col = "gray40")
    abline(h = 5,  lty = 2, col = "orange",    lwd = 1.5)
    abline(h = 7,  lty = 2, col = "darkgreen", lwd = 1.5)
    legend("topright", bty = "n", lty = 2, lwd = c(1, 1.5, 1.5),
           col = c("gray40", "orange", "darkgreen"),
           legend = c("TSS / baseline", "Minimum (5)", "ENCODE target (7)"))
    dev.off()

    rm(tsse)
    gc()

    # NFR score: log2(nf) + 1 - log2(n1 + n2), per TSS (400 bp window)
    # nf = middle 100 bp; n1/n2 = flanking 150 bp nucleosome positions
    nfr <- NFRscore(gal, txs)

    png(plot_nfr, width = 800, height = 700, res = 120)
    plot(nfr$log2meanCoverage, nfr$NFR_score,
         xlab = "log2 mean coverage",
         ylab = "Nucleosome Free Regions score",
         main = paste(sample_id, "- NFR Score"),
         xlim = c(-10, 0), ylim = c(-5, 5))
    dev.off()

    mean_nfr   <- mean(nfr$NFR_score,   na.rm = TRUE)
    median_nfr <- median(nfr$NFR_score, na.rm = TRUE)
    n_nfr      <- length(nfr)

    rm(nfr)
    gc()

    # PT score: log2(promoter density) - log2(body density), per transcript
    # promoter_window = [TSS-2000, TSS+500]; body_window = next 2500 bp downstream
    pt <- PTscore(gal, txs)

    png(plot_pt, width = 800, height = 700, res = 120)
    plot(pt$log2meanCoverage, pt$PT_score,
         xlab = "log2 mean coverage",
         ylab = "Promoter vs Transcript",
         main = paste(sample_id, "- PT Score"))
    dev.off()

    mean_pt   <- mean(pt$PT_score,       na.rm = TRUE)
    median_pt <- median(pt$PT_score,     na.rm = TRUE)
    mean_pro  <- mean(pt$promoter,       na.rm = TRUE)
    mean_body <- mean(pt$transcriptBody, na.rm = TRUE)
    n_pt      <- length(pt)

    rm(pt, gal, txs)
    gc()

    df <- data.frame(
        Sample           = sample_id,
        TSSE_score       = round(tsse_score, 4),
        NFR_score_mean   = round(mean_nfr,   4),
        NFR_score_median = round(median_nfr, 4),
        PT_score_mean    = round(mean_pt,    4),
        PT_score_median  = round(median_pt,  4),
        Mean_promoter    = round(mean_pro,   6),
        Mean_gene_body   = round(mean_body,  6)
    )
    write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

    qc_msg <- if (!is.nan(mean_pt) && 2^mean_pt >= 5)
        "[INFO] QC PT: PASS (promoter/body ratio >= 5)"
    else
        "[WARNING] QC PT: FAIL (promoter/body ratio < 5, expected >= 5-10)"

    writeLines(c(
        sprintf("[INFO] Transcripts processed (PT): %d",          n_pt),
        sprintf("[INFO] TSSE score: %.4f",                        tsse_score),
        sprintf("[INFO] TSS windows processed (NFR): %d",         n_nfr),
        sprintf("[INFO] Mean NFR score (log2 scale): %.4f",       mean_nfr),
        sprintf("[INFO] Median NFR score (log2 scale): %.4f",     median_nfr),
        sprintf("[INFO] Mean PT score (log2 ratio): %.4f",        mean_pt),
        sprintf("[INFO] Median PT score (log2 ratio): %.4f",      median_pt),
        sprintf("[INFO] Equivalent mean ratio (2^PT): %.2f",      2^mean_pt),
        qc_msg
    ), log_con)

}, error = function(e) {
    writeLines(sprintf("[ERROR] %s", conditionMessage(e)), log_con)
    stop(e)
}, finally = {
    close(log_con)
})
