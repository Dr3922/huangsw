"""Exploratory CNV-window-selected epithelial sensitivity analysis.

Uses a pre-existing, source-derived heuristic label. It does not validate malignancy.
"""
from pathlib import Path
import os, hashlib, json
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

PACKAGE_ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get("FGFR3_OUT", PACKAGE_ROOT/"data/derived")).expanduser().resolve()
H5=Path(os.environ.get("FGFR3_SCRNA_H5AD", PACKAGE_ROOT/"data/source/GSE169379_MIBC_snSeq.h5ad")).expanduser().resolve()
CNV=Path(os.environ.get("FGFR3_CNV_SCORE_FILE", PACKAGE_ROOT/"data/source/gse169379_cnv_windowed_cell_scores.csv")).expanduser().resolve()
SCRIPT=Path(os.environ.get("FGFR3_CNV_CLASSIFIER_SCRIPT", PACKAGE_ROOT/"data/source/fgfr3_gse169379_cnv_support.py")).expanduser().resolve()
SEED=20261021
APM=['HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','NLRC5']
IFNG=['STAT1','IRF1','CXCL9','CXCL10','CXCL11']
FEATURES=APM+IFNG+['apm_module','ifng_execution_module','antigen_visibility_module']

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def main():
 c=pd.read_csv(CNV)
 a=ad.read_h5ad(H5)
 assert c.cell_id.is_unique and a.obs_names.is_unique
 assert set(c.cell_id)==set(a.obs_names.astype(str))
 c=c.set_index('cell_id').loc[a.obs_names]
 for col in ['Patient','celltype','subtype']:
  assert (c[col].astype(str).to_numpy()==a.obs[col].astype(str).to_numpy()).all()
 epi=np.asarray(a.obs.celltype.astype(str)=='Epithelial')
 sel=epi & c.malignant_like.astype(bool).to_numpy()
 assert sel.sum()==3013
 genes=['FGFR3']+APM+IFNG
 idx=[a.var_names.get_loc(g) for g in genes]
 x=a.X[sel]
 obs=a.obs.loc[sel,['Patient','subtype']].copy(); obs['Patient']=obs.Patient.astype(str); obs['subtype']=obs.subtype.astype(str)
 patients=sorted(obs.Patient.unique())
 rows=[]
 for (p,s),ii in obs.groupby(['Patient','subtype'],observed=True).indices.items():
  sums=np.asarray(x[ii].sum(axis=0)).ravel(); total=float(sums.sum())
  v=np.log1p(sums[idx]*10000/total)
  r={'Patient':p,'subtype':s,'n_cells':len(ii),'library_size':total,**dict(zip(genes,v))}
  r['apm_module']=np.mean([r[g] for g in APM if g!='TAP2'])
  r['ifng_execution_module']=np.mean([r[g] for g in IFNG])
  r['antigen_visibility_module']=np.mean([r[g] for g in APM+IFNG if g!='TAP2'])
  rows.append(r)
 ps=pd.DataFrame(rows)
 prow=[]
 for p,ii in obs.groupby('Patient',observed=True).indices.items():
  sums=np.asarray(x[ii].sum(axis=0)).ravel(); total=float(sums.sum()); v=np.log1p(sums[idx]*10000/total)
  r={'Patient':p,'n_cells':len(ii),'library_size':total,**dict(zip(genes,v))}
  r['apm_module']=np.mean([r[g] for g in APM if g!='TAP2']); r['ifng_execution_module']=np.mean([r[g] for g in IFNG]); r['antigen_visibility_module']=np.mean([r[g] for g in APM+IFNG if g!='TAP2']); prow.append(r)
 pp=pd.DataFrame(prow); eligible=pp[pp.n_cells>=20].copy()
 ps.to_csv(OUT/'snrna/cnv_window_selected_patient_state_pseudobulk.csv',index=False)
 pp.to_csv(OUT/'snrna/cnv_window_selected_patient_pseudobulk.csv',index=False)
 rng=np.random.default_rng(SEED); boot_idx=rng.integers(0,len(eligible),size=(5000,len(eligible)))
 results=[]
 for f in FEATURES:
  if f=='TAP2':
   results.append(dict(feature=f,n=len(eligible),estimate=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,bootstrap_effective=0,bootstrap_failed=5000,reason='source h5ad TAP2 column is all zero'))
   continue
  rho,p=stats.spearmanr(eligible.FGFR3,eligible[f]); vals=[]
  for ii in boot_idx:
   r,_=stats.spearmanr(eligible.FGFR3.to_numpy()[ii],eligible[f].to_numpy()[ii]); vals.append(r)
  vals=np.asarray(vals); valid=vals[np.isfinite(vals)]
  lo,hi=np.quantile(valid,[.025,.975])
  results.append(dict(feature=f,n=len(eligible),estimate=rho,ci_low=lo,ci_high=hi,p=p,bootstrap_effective=len(valid),bootstrap_failed=5000-len(valid),reason=''))
 rr=pd.DataFrame(results); rr['family']=np.where(rr.feature.isin(APM+IFNG),'genes','modules'); rr['q']=np.nan
 for _,ii in rr.groupby('family').groups.items():
  jj=rr.loc[ii,'p'].dropna().index; rr.loc[jj,'q']=multipletests(rr.loc[jj,'p'],method='fdr_bh')[1]
 rr.to_csv(OUT/'snrna/cnv_window_selected_patient_correlations.csv',index=False)
 # State-adjusted sensitivity among selected-cell patient-state groups with >=10 nuclei.
 d=ps[ps.n_cells>=10].copy(); models=[]
 zfg=(d.FGFR3-d.FGFR3.mean())/d.FGFR3.std(ddof=1); dum=pd.get_dummies(d.subtype,drop_first=True,dtype=float)
 for f in FEATURES:
  if f=='TAP2' or d[f].std()==0: continue
  y=(d[f]-d[f].mean())/d[f].std(ddof=1); X=sm.add_constant(pd.concat([zfg.rename('FGFR3'),dum],axis=1))
  fit=sm.OLS(y,X).fit(cov_type='cluster',cov_kwds={'groups':d.Patient,'use_correction':True,'df_correction':True},use_t=True)
  ci=fit.conf_int().loc['FGFR3']; models.append(dict(feature=f,n_groups=len(d),n_patients=d.Patient.nunique(),estimate=fit.params.FGFR3,se=fit.bse.FGFR3,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=fit.pvalues.FGFR3,inference_df=getattr(fit,'df_resid_inference',d.Patient.nunique()-1)))
 mm=pd.DataFrame(models); mm['family']=np.where(mm.feature.isin(APM+IFNG),'genes','modules'); mm['q']=np.nan
 for _,ii in mm.groupby('family').groups.items(): mm.loc[ii,'q']=multipletests(mm.loc[ii,'p'],method='fdr_bh')[1]
 mm.to_csv(OUT/'snrna/cnv_window_selected_state_adjusted_models.csv',index=False)
 counts=c[epi].groupby('Patient').agg(epithelial_n=('malignant_like','size'),selected_n=('malignant_like','sum'),mean_cnv_burden=('cnv_burden','mean')).reset_index(); counts['selected_fraction']=counts.selected_n/counts.epithelial_n
 counts.to_csv(OUT/'snrna/cnv_window_selected_patient_qc.csv',index=False)
 audit={'source_h5ad':str(H5),'source_h5ad_sha256':sha(H5),'historical_cell_score_file':str(CNV),'historical_cell_score_sha256':sha(CNV),'historical_classifier_script':str(SCRIPT),'historical_classifier_script_sha256':sha(SCRIPT),'cell_ids_exactly_match':True,'metadata_exactly_match':['Patient','celltype','subtype'],'epithelial_cells':int(epi.sum()),'heuristic_selected_cells':int(sel.sum()),'all_selected_patients':len(pp),'primary_patients_min20_selected_cells':len(eligible),'threshold_definition':'cnv_burden above 99th percentile of windowed scores in pooled Lymphocyte/Fibroblast/Myeloid/Endothelial reference cells (historical script)','interpretation':'Exploratory CNV-window-deviation heuristic only; neither CopyKAT nor inferCNV ground truth and not a definitive malignant-cell annotation.','TAP2':'not estimable because source column is all zero','seed':SEED}
 (OUT/'snrna/cnv_window_selected_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
 print(json.dumps(audit,indent=2)); print(rr.to_string(index=False)); print(mm.to_string(index=False))

if __name__=='__main__': main()
