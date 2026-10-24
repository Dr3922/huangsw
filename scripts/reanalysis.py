"""FGFR3 Discover Oncology R1. Outputs only; original inputs remain read-only."""
from pathlib import Path
import os, sys, json, hashlib, re, time, warnings
import numpy as np
import pandas as pd
from scipy import stats, sparse
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import anndata as ad

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('FGFR3_OUT', PACKAGE_ROOT/'data/derived')).expanduser().resolve()
ROOT=Path(os.environ.get('FGFR3_PROJECT_ROOT', PACKAGE_ROOT/'data/source')).expanduser().resolve()
SC=Path(os.environ.get('FGFR3_SCRNA_H5AD', ROOT/'GSE169379_MIBC_snSeq.h5ad')).expanduser().resolve()
SP=Path(os.environ.get('FGFR3_SPATIAL_H5AD', ROOT/'GSE171351_combined_visium.h5ad')).expanduser().resolve()
SEED=20261021
APM=['HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','NLRC5']
IFNG=['STAT1','IRF1','CXCL9','CXCL10','CXCL11']
LUM=['GATA3','FOXA1','KRT20','PPARG','KRT8','KRT18','MUC1','TACSTD2','UPK1A','UPK1B','UPK2','UPK3A','EPCAM']
MODULES={'apm_module':APM,'ifng_execution_module':IFNG,'antigen_visibility_module':APM+IFNG,'luminal_score':LUM}
FEATURES=APM+IFNG+list(MODULES)[:3]
for d in ['snrna','spatial','bulk','cohort_ledgers','logs','figure_data','figures','public_evidence']:(OUT/d).mkdir(exist_ok=True)

def save(df,name):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,index=False)
def note(name,text):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text,encoding='utf-8')
def z(x):
 return (x-x.mean(axis=0))/x.std(axis=0,ddof=1).replace(0,np.nan) if isinstance(x,pd.DataFrame) else (x-np.nanmean(x,axis=0))/np.nanstd(x,axis=0,ddof=1)
def bh(df,groups,p='p'):
 df=df.copy();df['q']=np.nan
 for _,ii in df.groupby(groups,dropna=False).groups.items():
  jj=df.loc[ii,p].dropna().index
  if len(jj):df.loc[jj,'q']=multipletests(df.loc[jj,p],method='fdr_bh')[1]
 return df
def fit(y,x,groups=None,coef='FGFR3'):
 y=pd.Series(y).astype(float);x=pd.DataFrame(x).astype(float)
 ok=y.notna() & x.notna().all(axis=1)
 y=y.loc[ok];x=sm.add_constant(x.loc[ok],has_constant='add')
 if len(y)<=x.shape[1]+2 or np.linalg.matrix_rank(x)<x.shape[1]:return {'n':len(y),'estimate':np.nan,'p':np.nan,'reason':'insufficient n or rank deficient'}
 m=sm.OLS(y,x)
 if groups is None:r=m.fit(cov_type='HC3',use_t=True)
 else:r=m.fit(cov_type='cluster',cov_kwds={'groups':np.asarray(groups)[ok],'use_correction':True,'df_correction':True},use_t=True)
 ci=r.conf_int().loc[coef]
 return {'n':len(y),'patients':len(np.unique(np.asarray(groups)[ok])) if groups is not None else len(y),'estimate':r.params[coef],'se':r.bse[coef],'ci_low':ci.iloc[0],'ci_high':ci.iloc[1],'p':r.pvalues[coef],'r2':r.rsquared,'df_resid':r.df_resid,'inference_df':getattr(r,'df_resid_inference',r.df_resid),'condition_number':r.condition_number}

def rankrho(x,y):
 """Rows are bootstrap replicates; pairwise finite ranks, ties averaged."""
 valid=np.isfinite(x)&np.isfinite(y)
 a=stats.rankdata(np.where(valid,x,np.nan),axis=1,nan_policy='omit')
 b=stats.rankdata(np.where(valid,y,np.nan),axis=1,nan_policy='omit')
 n=valid.sum(axis=1)
 a=a-np.nanmean(a,axis=1,keepdims=True);b=b-np.nanmean(b,axis=1,keepdims=True)
 den=np.sqrt(np.nansum(a*a,axis=1)*np.nansum(b*b,axis=1))
 rho=np.divide(np.nansum(a*b,axis=1),den,out=np.full(len(n),np.nan),where=den>0)
 rho[n<4]=np.nan
 return rho,n
def pooled_array(arr,idx):
 # arr[patient,state,FGFR3+features]. Patient multiplicity is preserved by indexing.
 samples=arr[idx]
 rhos=[];weights=[]
 for s in range(arr.shape[1]):
  rr=[];ww=[]
  for f in range(1,arr.shape[2]):
   r,n=rankrho(samples[:,:,s,0],samples[:,:,s,f]);valid=np.isfinite(r)&(n>3)&(abs(r)<1)
   rr.append(np.where(valid,np.arctanh(np.clip(r,-.999999999,.999999999)),0));ww.append(np.where(valid,n-3,0))
  rhos.append(np.array(rr).T);weights.append(np.array(ww).T)
 w=np.array(weights);rs=np.array(rhos);den=w.sum(axis=0)
 vals=np.tanh(np.divide((w*rs).sum(axis=0),den,out=np.full_like(den,np.nan,dtype=float),where=den>0))
 # Fixed requirement: all five predefined states estimable, otherwise iteration invalid.
 vals[(w>0).sum(axis=0)!=arr.shape[1]]=np.nan
 return vals
def bootstrap_summary(draws,point):
 valid=draws[np.isfinite(draws)]
 lo,hi=np.quantile(valid,[.025,.975]) if len(valid) else (np.nan,np.nan)
 return {'estimate':point,'ci_low':lo,'ci_high':hi,'bootstrap_effective':len(valid),'bootstrap_failed':len(draws)-len(valid),'negative_fraction':np.mean(valid<0) if len(valid) else np.nan,'p':(1+np.sum(abs(valid-point)>=abs(point)))/(len(valid)+1) if len(valid) else np.nan}

def snrna():
 t=time.time();a=ad.read_h5ad(SC)
 integer=np.all(a.X.data>=0) and np.all(a.X.data==np.floor(a.X.data))
 print('source',a.shape,'integer_nonnegative',integer,flush=True)
 mask=np.asarray(a.obs.celltype=='Epithelial');e=a[mask].copy();del a
 obs=e.obs.copy();obs['Patient']=obs.Patient.astype(str);obs['subtype']=obs.subtype.astype(str)
 pats=sorted(obs.Patient.unique());states=sorted(obs.subtype.unique())
 genes=['FGFR3']+sorted(set(APM+IFNG+LUM));present=[g for g in genes if g in e.var_names]
 expr=pd.DataFrame(e[:,present].X.toarray(),columns=present,index=obs.index)
 for g in set(genes)-set(present):expr[g]=np.nan
 for k,gs in MODULES.items():expr[k]=expr[[g for g in gs if g in present]].mean(axis=1)
 cell=obs[['Patient','subtype','batch']].join(expr)
 qc=obs.groupby('Patient',observed=True).agg(epithelial_cells=('subtype','size'),n_batches=('batch','nunique'))
 qc['total_UMI']=pd.Series(np.asarray(e.X.sum(axis=1)).ravel(),index=obs.index).groupby(obs.Patient).sum()
 save(qc.reset_index(),'snrna/epithelial_patient_qc.csv')
 count=obs.groupby(['Patient','subtype'],observed=True).size().reindex(pd.MultiIndex.from_product([pats,states],names=['Patient','subtype']),fill_value=0).rename('n_cells').reset_index()
 count['status']=np.where(count.n_cells==0,'observed_zero','observed')
 count['eligible_20']=count.n_cells>=20;count['eligible_50']=count.n_cells>=50
 save(count,'snrna/epithelial_patient_by_state_counts.csv')
 obs[['Patient','batch','hash_ID']].drop_duplicates().to_csv(OUT/'snrna/patient_batch_hash_ledger.csv',index=False)
 old=cell.groupby(['Patient','subtype'],observed=True)[['FGFR3']+FEATURES+['luminal_score']].mean().reset_index().merge(count,on=['Patient','subtype'])
 oldpat=cell.groupby('Patient',observed=True)[['FGFR3']+FEATURES+['luminal_score']].mean().reset_index()
 save(old,'snrna/submitted_raw_count_mean_patient_state.csv');save(oldpat,'snrna/submitted_raw_count_mean_patient.csv')
 # Genuine count pseudobulk: sum all genes per patient-state, library normalize then log1p.
 assert integer,'Count normalization requires verified nonnegative integer source'
 rows=[];prows=[]
 for p in pats:
  for s in [None]+states:
   m=np.array(obs.Patient==p,copy=True)
   if s is not None:m &=np.asarray(obs.subtype==s)
   if not m.any():continue
   sums=np.asarray(e.X[m].sum(axis=0)).ravel();total=sums.sum();v=np.log1p(sums*10000/total)
   row={'Patient':p,'subtype':s,'n_cells':m.sum(),'library_size':total}
   row.update({g:v[e.var_names.get_loc(g)] if g in present else np.nan for g in genes})
   for k,gs in MODULES.items():row[k]=np.mean([row[g] for g in gs if g in present])
   (prows if s is None else rows).append(row)
 pb=pd.DataFrame(rows);pp=pd.DataFrame(prows).drop(columns='subtype')
 save(pb,'snrna/pseudobulk_expression_patient_state.csv');save(pp,'snrna/pseudobulk_expression_patient.csv')
 note('snrna/pseudobulk_transformation_audit.md',f'''# Source and transformation audit
Source: {SC}; shape {e.shape}; {len(pats)} epithelial donors, {len(states)} states, {len(e)} nuclei.
All stored nonzero X entries tested nonnegative and integer: {integer}. No raw slot, no separate normalized layers. X has count-compatible values, not log expression. Published processing arguments: {repr(e.uns)}.
Submitted scripts used means of unnormalized X values. This preserves the submitted analysis for exact reproduction, but introduces library-depth/composition concerns. Revised sensitivity uses genuine all-gene UMI sums per donor or donor-state, log1p(10,000 * gene counts / all-gene total). Module means are calculated after normalization. No summing of log values is described as raw-count pseudobulk.
Primary minimum 20 nuclei per donor-state (sampling precision safeguard, not a universal validated threshold); sensitivities 1 and 50. All five states retained; missing groups remain missing, never zeros in expression. Source annotations are epithelial, not definitive malignant labels. Molecular CNV analyses are exploratory and do not establish malignant ground truth. Donors across batches are grouped by source Patient before inference.
''')
 rng=np.random.default_rng(SEED);idx=rng.integers(0,len(pats),size=(5000,len(pats)))
 np.save(OUT/'snrna/bootstrap_patient_indices.npy',idx)
 summary=[];detail=[];drawrows=[];influence=[];models=[];increment=[];pcorr=[]
 for scale,df,pdf in [('submitted_raw_mean',old,oldpat),('logCP10k_count_pseudobulk',pb,pp)]:
  for minimum in [1,20,50]:
   sub=df[df.n_cells>=minimum].copy()
   arr=np.full((len(pats),len(states),1+len(FEATURES)),np.nan)
   for r in sub.to_dict('records'):arr[pats.index(r['Patient']),states.index(r['subtype'])]=[r[g] for g in ['FGFR3']+FEATURES]
   point=pooled_array(arr,np.arange(len(pats))[None,:])[0];draw=pooled_array(arr,idx)
   for j,f in enumerate(FEATURES):
    meta=bootstrap_summary(draw[:,j],point[j]);meta.update(scale=scale,min_cells=minimum,feature=f,patients=sub.Patient.nunique(),patient_states=len(sub),family='genes' if f in APM+IFNG else 'modules');summary.append(meta)
    for b,val in enumerate(draw[:,j]):drawrows.append({'scale':scale,'min_cells':minimum,'feature':f,'iteration':b+1,'rho':val})
    for s in states:
     dd=sub[sub.subtype==s];r,p=stats.spearmanr(dd.FGFR3,dd[f]);n=len(dd);se=1/np.sqrt(n-3) if n>3 else np.nan
     lo,hi=np.tanh(np.arctanh(np.clip(r,-.99999999,.99999999))+np.array([-1,1])*1.96*se)
     detail.append(dict(scale=scale,min_cells=minimum,subtype=s,feature=f,n=n,rho=r,p=p,se_fisher=se,weight=n-3,ci_low=lo,ci_high=hi))
   for k,p in enumerate(pats):
    vals=pooled_array(arr,np.delete(np.arange(len(pats)),k)[None,:])[0]
    for f,v,ref in zip(FEATURES,vals,point):influence.append(dict(scale=scale,min_cells=minimum,excluded_patient=p,feature=f,rho=v,change=v-ref))
   x=z(sub[['FGFR3','luminal_score']]);dummy=pd.get_dummies(sub.subtype,drop_first=True,dtype=float)
   for f in FEATURES:
    y=z(sub[f].to_numpy());base=pd.concat([x[['luminal_score']],dummy],axis=1)
    for mod,xx in [('state_adjusted',pd.concat([x[['FGFR3']],dummy],axis=1)),('state_luminal_adjusted',pd.concat([x,dummy],axis=1))]:
     res=fit(pd.Series(y,index=sub.index),xx,sub.Patient)
     res.update(scale=scale,min_cells=minimum,feature=f,model=mod,family='genes' if f in APM+IFNG else 'modules');models.append(res)
    full=fit(pd.Series(y,index=sub.index),pd.concat([x,dummy],axis=1),sub.Patient)
    yy=pd.Series(y,index=sub.index);mb=sm.OLS(yy,sm.add_constant(base)).fit();mf=sm.OLS(yy,sm.add_constant(pd.concat([x,dummy],axis=1))).fit()
    full.update(scale=scale,min_cells=minimum,feature=f,partial_R2=(mb.ssr-mf.ssr)/mb.ssr,baseline_R2=mb.rsquared,full_R2=mf.rsquared)
    increment.append(full)
  # Patient-only correlations and patient bootstrap. All 25 independent donors.
  for f in FEATURES:
   r,p=stats.spearmanr(pdf.FGFR3,pdf[f]);dd,_=rankrho(pdf.FGFR3.to_numpy()[idx],pdf[f].to_numpy()[idx]);ss=bootstrap_summary(dd,r);ss['p']=p
   ss.update(scale=scale,feature=f,n=len(pdf),family='genes' if f in APM+IFNG else 'modules');pcorr.append(ss)
   for patient in pats:
    d=pdf[pdf.Patient!=patient];rr,ppv=stats.spearmanr(d.FGFR3,d[f]);influence.append(dict(scale=scale,min_cells='patient_only',excluded_patient=patient,feature=f,rho=rr,change=rr-r))
 smry=bh(pd.DataFrame(summary),['scale','min_cells','family']);save(smry,'snrna/within_state_bootstrap_summary.csv')
 save(pd.DataFrame(drawrows),'snrna/within_state_patient_cluster_bootstrap.csv')
 detail=pd.DataFrame(detail)
 naive=[]
 for key,d in detail.groupby(['scale','min_cells','feature']):
  d=d[np.isfinite(d.rho)&(d.n>3)&(d.rho.abs()<1)]
  zz=np.arctanh(d.rho);w=d.weight
  if len(d)==0 or w.sum()<=0:
   naive.append(dict(zip(['scale','min_cells','feature'],key),rho=np.nan,p=np.nan,ci_low=np.nan,ci_high=np.nan,invalid_independence_assumption=True));continue
  zm=np.average(zz,weights=w);se=1/np.sqrt(w.sum());lo,hi=np.tanh(zm+np.array([-1,1])*1.96*se)
  naive.append(dict(zip(['scale','min_cells','feature'],key),rho=np.tanh(zm),p=2*stats.norm.sf(abs(zm/se)),ci_low=lo,ci_high=hi,invalid_independence_assumption=True))
 save(detail,'snrna/within_state_per_state_estimates.csv');save(pd.DataFrame(naive),'snrna/submitted_naive_fisher_pooling.csv')
 save(pd.DataFrame(influence),'snrna/within_state_influence.csv');save(pd.DataFrame(influence),'snrna/epithelial_patient_influence.csv')
 mm=bh(pd.DataFrame(models),['scale','min_cells','model','family']);save(mm,'snrna/state_adjusted_robust_models.csv')
 cc=bh(pd.DataFrame(pcorr),['scale','family']);save(cc,'snrna/patient_pseudobulk_correlations.csv')
 save(pd.DataFrame(increment),'snrna/atlas_incremental_nested_models.csv')
 full=pd.concat([smry.assign(analysis='cluster_bootstrap_pooling'),mm.assign(analysis='regression'),cc.assign(analysis='patient_correlation')],ignore_index=True)
 full['evidence_level']=np.select([full.q<.05,full.p<.05,full.estimate.notna()],['BH q < .05','nominal P < .05','direction only'],default='not evaluable')
 save(full,'snrna/epithelial_full_prespecified_results.csv');save(full,'snrna/epithelial_effect_ci_plotsource.csv');save(smry,'figure_data/figure5_analysis_data.csv')
 note('snrna/within_state_methods.md','''# Within-state methods
Prespecified genes: seven APM and five IFNG; modules: APM, IFNG execution, combined antigen visibility. Correlation statistic: state Spearman rho, Fisher-z weighted by n-3, inverse transformed. Naive independent-state CI/P are reproduced only for audit and invalid for confirmatory inference because patients recur across states.
Primary revised inference: donor-cluster bootstrap, seed 20261021, 5000 draws of 25 donor IDs with replacement, preserving all eligible donor-state rows and duplicate multiplicities. All five predefined states must be evaluable; failures are reported without selective retries. Percentile 95% CI. Approximate two-sided centered bootstrap null P=(1+sum(abs(T*-T)>=abs(T)))/(B_valid+1); negative fraction is stability, not a P value. CI stability and donor deletion results supplied. Thresholds 1/20/50, submitted raw means and logCP10k count pseudobulk reported separately. Prespecified primary revised scale is normalized count pseudobulk with >=20 cells. Historical reproduction uses raw means and >=1.
Models standardize FGFR3 and outcome across eligible donor-state rows, use state fixed effects, optionally luminal mean without FGFR3, donor-cluster sandwich finite-sample correction and t inference (clusters-1 degrees of freedom). These are sensitivity models. No causal or protein-level interpretation.
''')
 note('snrna/epithelial_fdr_families.md','BH families separate 12 genes from three modules for each analysis (patient correlation, within-state bootstrap, state-adjusted model, state+luminal-adjusted model), scale and minimum-cell threshold. Sensitivities are not independent replications. Bootstrap P is approximate centered-resampling inference; patient Spearman P is asymptotic. Missing tests remain NA. No exploratory marker scanning enters these families.')
 print('snrna complete',round(time.time()-t,1),'seconds',smry.query("scale=='logCP10k_count_pseudobulk' and min_cells==20 and family=='modules'").to_dict('records'),flush=True)

if __name__=='__main__':
 globals()[sys.argv[1]]()
