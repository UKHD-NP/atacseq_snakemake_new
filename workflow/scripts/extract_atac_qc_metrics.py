"""Extract TSS enrichment and NFR metrics from ataqv JSON for MultiQC."""
import json

with open(snakemake.input.ataqv_json) as fh:
    data = json.load(fh)

metrics = data[0]["metrics"]
sample_id = snakemake.wildcards.sample_id


def _fmt(val):
    if val is None:
        return "NA"
    if isinstance(val, float):
        return f"{val:.4f}"
    return str(val)


tss_enrichment = _fmt(metrics.get("tss_enrichment"))
nfr_ratio      = _fmt(metrics.get("short_mononucleosomal_ratio"))
hqaa_tf        = _fmt(metrics.get("hqaa_tf_count"))
hqaa_mono      = _fmt(metrics.get("hqaa_mononucleosomal_count"))

with open(snakemake.output.mqc_tsv, "w") as out:
    out.write("sample\ttss_enrichment\tnfr_ratio\thqaa_tf_count\thqaa_mono_count\n")
    out.write(f"{sample_id}\t{tss_enrichment}\t{nfr_ratio}\t{hqaa_tf}\t{hqaa_mono}\n")
