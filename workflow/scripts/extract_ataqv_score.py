"""Extract NFR score and TSS enrichment from ataqv JSON for MultiQC."""
import json

with open(snakemake.input.ataqv_json) as fh:
    data = json.load(fh)

if not data:
    raise RuntimeError(f"ataqv JSON is empty: {snakemake.input.ataqv_json}")

metrics = data[0].get("metrics", {})
sample_id = snakemake.wildcards.sample_id


def _fmt(val):
    if val is None:
        return "NA"
    if isinstance(val, float):
        return f"{val:.4f}"
    return str(val)


nfr_ratio      = _fmt(metrics.get("short_mononucleosomal_ratio"))
hqaa_tf        = _fmt(metrics.get("hqaa_tf_count"))
hqaa_mono      = _fmt(metrics.get("hqaa_mononucleosomal_count"))
tss_enrichment = _fmt(metrics.get("tss_enrichment"))

with open(snakemake.output.mqc_tsv, "w") as out:
    out.write("Sample\tnfr_ratio\thqaa_tf_count\thqaa_mono_count\ttss_enrichment\n")
    out.write(f"{sample_id}\t{nfr_ratio}\t{hqaa_tf}\t{hqaa_mono}\t{tss_enrichment}\n")
