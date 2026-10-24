from __future__ import annotations
import argparse, hashlib, json, urllib.request
from pathlib import Path

SOURCES={
"GSE169379_MIBC_snSeq.h5ad.gz": {"url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169379/suppl/GSE169379_MIBC_snSeq.h5ad.gz","sha256_uncompressed":"3944f576cdf832aea0461dabf7ededd58bb3ddf70b6effb76372fa8a61166934"},
"GSE171351_combined_visium.h5ad.gz": {"url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171351/suppl/GSE171351_combined_visium.h5ad.gz","sha256_uncompressed":"1d59b90b1da1255c62ca7c8c40eca0dcf1d5b945ec210f80e03d8476821f17d6"},
}
p=argparse.ArgumentParser(); p.add_argument("--out",required=True); p.add_argument("--download",action="store_true"); a=p.parse_args()
out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
print(json.dumps(SOURCES,indent=2))
if a.download:
    for name,x in SOURCES.items():
        target=out/name
        urllib.request.urlretrieve(x["url"],target)
        print(f"downloaded {target}; verify the uncompressed h5ad against sha256_uncompressed")
