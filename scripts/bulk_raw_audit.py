"""Audit raw GEO probe-to-gene processing without modifying source inputs."""
from pathlib import Path
import os, csv, gzip, io, json, re
import numpy as np
import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("FGFR3_OUT", PACKAGE_ROOT/"data/derived")).expanduser().resolve()
ROOT = Path(os.environ.get("FGFR3_PROJECT_ROOT", PACKAGE_ROOT/"data/source")).expanduser().resolve()
RAW = ROOT / "geo_validation" / "raw"
PROC = ROOT / "geo_validation" / "processed"


def clean_symbol(value):
    value = "" if pd.isna(value) else str(value).strip()
    if not value or value == "---":
        return ""
    return re.split(r"\s*///\s*|\s*//\s*|\s*;\s*", value)[0].strip()


def parse_matrix(path):
    meta, rows, samples = {}, [], None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        in_table = False
        for line in handle:
            if line.startswith("!Sample_"):
                row = next(csv.reader([line], delimiter="\t"))
                meta.setdefault(row[0][8:], []).append(row[1:])
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                row = next(csv.reader([line], delimiter="\t"))
                if samples is None:
                    samples = [x.strip('"') for x in row[1:]]
                else:
                    rows.append([row[0].strip('"')] + [float(x.strip('"')) if x.strip('"') not in ("", "null", "NA") else np.nan for x in row[1:]])
    return pd.DataFrame(rows, columns=["probe_id"] + samples).set_index("probe_id"), meta


def read_annotation(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        started = False
        lines = []
        for line in handle:
            if line.startswith("!platform_table_begin"):
                started = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if started:
                lines.append(line)
        text = "".join(lines)
        ann = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, low_memory=False)
    probe = next(c for c in ann.columns if c.lower() in {"id", "id_ref", "probe set id", "probe_id"})
    symbol = next(c for c in ann.columns if c.lower().replace("_", " ") in {"gene symbol", "gene assignment"})
    out = pd.DataFrame({"probe_id": ann[probe].astype(str), "symbol_raw": ann[symbol]})
    out["gene_symbol"] = out.symbol_raw.map(clean_symbol)
    out["ambiguous"] = out.symbol_raw.fillna("").str.contains(r"///|//|;", regex=True)
    return out[out.gene_symbol.ne("")].drop_duplicates("probe_id")


def compare_rules(accession, annotation_name):
    expression, meta = parse_matrix(RAW / f"{accession}_series_matrix.txt.gz")
    ann = read_annotation(RAW / annotation_name).set_index("probe_id")
    merged = expression.join(ann, how="inner")
    cached = pd.read_csv(PROC / f"{accession}_gene_expression_raw.csv", index_col=0)
    common_samples = cached.columns.intersection(expression.columns)
    candidates = {
        "mean_all_mapped_first_symbol": merged.groupby("gene_symbol")[common_samples].mean(),
        "median_all_mapped_first_symbol": merged.groupby("gene_symbol")[common_samples].median(),
        "mean_unambiguous_only": merged.loc[~merged.ambiguous].groupby("gene_symbol")[common_samples].mean(),
        "median_unambiguous_only": merged.loc[~merged.ambiguous].groupby("gene_symbol")[common_samples].median(),
    }
    rows = []
    for name, cand in candidates.items():
        genes = cached.index.intersection(cand.index)
        delta = cached.loc[genes, common_samples].to_numpy() - cand.loc[genes, common_samples].to_numpy()
        finite = np.isfinite(delta)
        rows.append({"cohort": accession, "candidate_rule": name, "cached_genes": len(cached), "candidate_genes": len(cand),
                     "common_genes": len(genes), "common_samples": len(common_samples), "compared_values": int(finite.sum()),
                     "max_abs_difference": float(np.nanmax(np.abs(delta))), "mean_abs_difference": float(np.nanmean(np.abs(delta))),
                     "exact_within_1e-10_fraction": float(np.mean(np.abs(delta[finite]) <= 1e-10))})
    return rows, meta


def main():
    rows, metas = [], {}
    for accession, annotation in [("GSE13507", "GSE13507_platform.annot.gz"), ("GSE32894", "GSE32894_platform.annot.gz")]:
        r, m = compare_rules(accession, annotation)
        rows.extend(r); metas[accession] = m
    results = pd.DataFrame(rows)
    results.to_csv(OUT / "public_evidence" / "bulk_probe_collapse_rule_audit.csv", index=False)
    best = results.sort_values(["cohort", "max_abs_difference", "mean_abs_difference"]).groupby("cohort").head(1)
    m = pd.read_csv(PROC / "GSE32894_sample_metadata.csv", index_col=0)
    combined = m.title.fillna("").astype(str) + " " + m.source_name_ch1.fillna("").astype(str)
    patient = combined.str.extract(r"(UC_[0-9]+)", expand=False)
    ledger = pd.DataFrame({"gsm": m.index, "title": m.title, "source_name": m.source_name_ch1, "patient_prefix": patient})
    ledger["prefix_sample_count"] = ledger.patient_prefix.map(ledger.patient_prefix.value_counts())
    ledger.to_csv(OUT / "cohort_ledgers" / "gse32894_patient_prefix_audit.csv", index=False)
    summary = {
        "best_matching_rules": best.to_dict("records"),
        "gse32894_samples": int(len(ledger)), "gse32894_mapped_prefixes": int(patient.notna().sum()),
        "gse32894_unique_prefixes": int(patient.nunique()), "gse32894_repeated_prefix_samples": int((ledger.prefix_sample_count > 1).sum()),
        "interpretation": "A UC_#### prefix is treated as a public specimen/patient identifier only when present; uniqueness is audited rather than assumed.",
    }
    (OUT / "public_evidence" / "bulk_raw_processing_audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
