from pathlib import Path
import os
import json, gzip, zipfile, xml.etree.ElementTree as ET
import anndata as ad
PACKAGE_ROOT=Path(__file__).resolve().parents[1]
ROOT=Path(os.environ.get('FGFR3_PROJECT_ROOT', PACKAGE_ROOT/'data/source')).expanduser().resolve()
out=Path(os.environ.get('FGFR3_OUT', PACKAGE_ROOT/'data/derived')).expanduser().resolve()
for acc,name in [('GSE169379','MIBC_snSeq'),('GSE171351','combined_visium')]:
 p=ROOT/f'{acc}_{name}.h5ad'
 a=ad.read_h5ad(p,backed='r')
 print(acc,a, 'obs',a.obs.columns.tolist(),'var',a.var.columns.tolist(), 'layers',list(a.layers),'uns',list(a.uns),'raw',a.raw is not None,flush=True)
 a.obs.to_csv(out/f'{acc}_source_obs.csv')
 (out/f'{acc}_object_metadata.txt').write_text(str(a)+'\n'+repr(a.uns)+'\n'+str(a.obs.describe(include='all')),encoding='utf-8')
 print(a.obs.head().to_string(),flush=True)
 for c in ['Patient','subtype','sampleID','celltype']:
  if c in a.obs: print(c,a.obs[c].value_counts().to_dict(),flush=True)
 print('X sample',a.X[:10,:10].toarray() if hasattr(a.X[:10,:10],'toarray') else a.X[:10,:10],flush=True)
 if acc=='GSE169379':print(a.obs.loc[a.obs.celltype=='Epithelial',['Patient','subtype']].value_counts().to_string(),flush=True)
 a.file.close()
base=Path(os.environ.get('FGFR3_REVISION_ROOT', PACKAGE_ROOT.parent))
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
for p in base.glob('*.docx'):
 with zipfile.ZipFile(p) as z:
  r=ET.fromstring(z.read('word/document.xml'))
  text='\n'.join(''.join(t.itertext()) for t in r.findall('.//w:p',ns))
  (out/(p.stem+'.txt')).write_text(text,encoding='utf-8')
