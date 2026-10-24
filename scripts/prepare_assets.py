from pathlib import Path
import os
import anndata as ad
import numpy as np
import json, requests, hashlib, datetime
PACKAGE_ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get("FGFR3_OUT", PACKAGE_ROOT/"data/derived")).expanduser().resolve()
p=OUT/'figure_data/spatial_assets';p.mkdir(exist_ok=True)
spatial_h5ad = os.environ.get("FGFR3_SPATIAL_H5AD")
if not spatial_h5ad:
    raise SystemExit("Set FGFR3_SPATIAL_H5AD to the public GSE171351 h5ad before preparing figure assets")
a=ad.read_h5ad(spatial_h5ad,backed='r')
rows=[]
for sec in ['A1','B1','C1','D1']:
 info=a.uns['spatial'][sec]
 np.save(p/f'{sec}_hires.npy',info['images']['hires'])
 scal={k:float(v) for k,v in info['scalefactors'].items()}
 (p/f'{sec}_scalefactors.json').write_text(json.dumps(scal,indent=2))
 rows.append(dict(section=sec,bar_microns=500,nominal_spot_microns=55,bar_hires_pixels=500/55*scal['spot_diameter_fullres']*scal['tissue_hires_scalef'],**scal))
a.file.close()
import pandas as pd
pd.DataFrame(rows).to_csv(OUT/'figure_data/Figure4_scale_bar_audit.csv',index=False)
url='https://cdn.10xgenomics.com/image/upload/v1758134784/support-documents/CG000686_VisiumHDSpatialApplicationsPerformance_TechnicalNote_RevB.pdf'
r=requests.get(url,timeout=60);r.raise_for_status()
q=OUT/'public_evidence/literature/10x_Visium_spot_geometry_RevB.pdf';q.write_bytes(r.content)
(OUT/'public_evidence/literature/10x_scale_source.json').write_text(json.dumps(dict(url=url,sha256=hashlib.sha256(r.content).hexdigest(),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),evidence='Page 3 specifies Visium v1/v2 55 micrometer diameter and 100 micrometer center spacing; application uses object spot_diameter_fullres and tissue_hires_scalef',figure='500/55 * spot_diameter_fullres * tissue_hires_scalef'),indent=2))
print('spatial assets and official scale reference saved')
