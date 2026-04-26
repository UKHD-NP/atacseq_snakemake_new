rule atacseqqc_score:
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.shifted.bam.bai"),
        bed = config["ref"]["bed"],
    output:
        score_tsv    = os.path.join("{outdir}", "atacseqqc", "{sample_id}.atacseqqc_score.tsv"),
        plot_fragsize = os.path.join("{outdir}", "atacseqqc", "{sample_id}.fragsize_dist.png"),
        plot_pt       = os.path.join("{outdir}", "atacseqqc", "{sample_id}.pt_score.png"),
        plot_nfr      = os.path.join("{outdir}", "atacseqqc", "{sample_id}.nfr_score.png"),
        plot_tsse     = os.path.join("{outdir}", "atacseqqc", "{sample_id}.tsse.png"),
    conda:
        os.path.join(workflow.basedir, "envs", "atacseqqc.yml")
    message:
        "{wildcards.sample_id}: Calculating PT score, NFR score and TSSE score via ATACseqQC"
    threads: 1
    resources:
        mem_mb = 16384
    log:
        os.path.join("{outdir}", "logs", "atacseqqc", "{sample_id}.pt_score.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.pt_score.benchmark.txt")
    script:
        os.path.join(workflow.basedir, "scripts", "calc_pt_score.R")
