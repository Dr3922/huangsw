import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("FGFR3_OUT", PACKAGE_ROOT/"data/derived")).expanduser().resolve()
CBIO = OUT / "public_evidence" / "cbio"
API = "https://www.cbioportal.org/api"
GENE = 2261


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run():
    rows = []
    session = requests.Session()
    for study_dir in sorted(p for p in CBIO.iterdir() if p.is_dir()):
        profiles_file = study_dir / "profiles.json"
        lists_file = study_dir / "lists.json"
        if not profiles_file.exists() or not lists_file.exists():
            continue
        profiles = load(profiles_file)
        lists = load(lists_file)
        discrete = [p for p in profiles if p.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION" and p.get("datatype") == "DISCRETE"]
        cna_lists = [x for x in lists if x.get("category") == "all_cases_with_cna_data"]
        status = "ok"
        records = []
        url = ""
        message = ""
        try:
            if len(discrete) != 1 or len(cna_lists) != 1:
                status = "not_available"
                message = f"discrete_profiles={len(discrete)}; cna_lists={len(cna_lists)}"
            else:
                profile_id = discrete[0]["molecularProfileId"]
                sample_list_id = cna_lists[0]["sampleListId"]
                endpoint = f"{API}/molecular-profiles/{profile_id}/molecular-data"
                params = {"sampleListId": sample_list_id, "entrezGeneId": GENE, "projection": "DETAILED"}
                response = session.get(endpoint, params=params, timeout=180)
                response.raise_for_status()
                url = response.url
                records = response.json()
                returned = {r.get("entrezGeneId", r.get("gene", {}).get("entrezGeneId")) for r in records}
                if returned - {GENE}:
                    raise RuntimeError(f"unexpected returned genes: {sorted(returned)}")
                eligible_file = study_dir / f"{sample_list_id}.json"
                eligible = set(load(eligible_file))
                samples = [r.get("sampleId") for r in records]
                if len(samples) != len(set(samples)):
                    raise RuntimeError("duplicate sample records returned")
                if not set(samples) <= eligible:
                    raise RuntimeError("records returned outside CNA sample list")
                dump(study_dir / "FGFR3_cna_gene_specific.json", records)
        except Exception as exc:
            status = "error"
            message = repr(exc)
        output = study_dir / "FGFR3_cna_gene_specific.json"
        rows.append({
            "study_id": study_dir.name,
            "status": status,
            "record_n": len(records),
            "deletion_n": sum(r.get("value") == -2 for r in records),
            "loss_n": sum(r.get("value") == -1 for r in records),
            "diploid_n": sum(r.get("value") == 0 for r in records),
            "gain_n": sum(r.get("value") == 1 for r in records),
            "amplification_n": sum(r.get("value") == 2 for r in records),
            "endpoint_url": url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "output_sha256": sha256(output) if output.exists() and status == "ok" else "",
            "message": message,
        })
        print(study_dir.name, status, len(records), flush=True)
    pd.DataFrame(rows).to_csv(OUT / "public_evidence" / "cbio_fgfr3_cna_gene_specific_acquisition.csv", index=False)


if __name__ == "__main__":
    run()
