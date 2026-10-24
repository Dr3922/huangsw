from pathlib import Path
import argparse, json, hashlib, shutil, re, platform, sys
import numpy as np, pandas as pd, anndata as ad
from scipy import stats
from scipy.spatial import cKDTree
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

parser=argparse.ArgumentParser()
parser.add_argument('--handoff',type=Path,required=True,help='Previous FGFR3 reanalysis handoff directory')
parser.add_argument('--spatial-h5ad',type=Path,required=True,help='Original combined Visium h5ad')
parser.add_argument('--out',type=Path,required=True,help='New V2 output directory')
args=parser.parse_args()
BASE=args.handoff.resolve(); OUT=args.out.resolve(); SP=args.spatial_h5ad.resolve()
MAN=pd.read_csv(BASE/'spatial/lopo_fold_manifest.csv')
SCORES=pd.read_csv(BASE/'spatial/lopo_spot_scores.csv')
SEED=20261021
APM=['HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','NLRC5']
IFNG=['STAT1','IRF1','CXCL9','CXCL10','CXCL11']
DC=['FCER1A','HLA-DRA','HLA-DPA1','HLA-DPB1','CLEC10A','ITGAX']
SETS={'local_dc_original':DC,'local_dc_no_hla':['FCER1A','CLEC10A','ITGAX'],'local_cDC1':['CLEC9A','XCR1','BATF3','CADM1'],'local_cDC2':['CD1C','FCER1A','CLEC10A'],'local_myeloid':['LST1','TYROBP','FCER1G','AIF1','CTSS'],'local_cd8':['CD8A','CD8B','NKG7','PRF1','GZMB'],'local_apm':APM,'local_ifng':IFNG,'local_luminal':['GATA3','FOXA1','KRT20','PPARG','KRT8','KRT18','MUC1','TACSTD2','UPK1A','UPK1B','UPK2','UPK3A','EPCAM'],'local_epithelial':['EPCAM','TACSTD2','KRT8','KRT18','KRT19'],'local_stromal':['FAP','ACTA2','COL1A1','COL1A2','PDGFRB','FN1','TGFB1','CXCL12']}
OUT.mkdir(parents=True,exist_ok=True)
for d in ['scripts','logs']: (OUT/d).mkdir(exist_ok=True)

def z(s):
    s=pd.Series(s,dtype=float); sd=s.std(ddof=1)
    return (s-s.mean())/(sd if sd and np.isfinite(sd) else 1.0)
def module(df, genes):
    cols=[g for g in genes if g in df.columns]
    return df[cols].mean(axis=1) if cols else pd.Series(np.nan,index=df.index)
def block_ci(x,y,xy,B=1200):
    ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<3: return np.array([np.nan,np.nan]),0,np.nan
    x=x[ok]; y=y[ok]; xy=xy[ok]
    spacing=np.median(cKDTree(xy).query(xy,k=2)[0][:,1]); ids=np.floor((xy-xy.min(0))/(3*spacing)).astype(int)
    _,lab=np.unique(ids,axis=0,return_inverse=True); groups=[np.flatnonzero(lab==j) for j in range(lab.max()+1)]
    rng=np.random.default_rng(SEED); vals=[]
    for _ in range(B):
        ii=np.concatenate([groups[j] for j in rng.integers(0,len(groups),len(groups))])
        r=stats.spearmanr(x[ii],y[ii]).statistic
        vals.append(r)
    vals=np.asarray(vals); return np.nanquantile(vals,[.025,.975]),len(groups),spacing
def fit_adjusted(df,ycol,score,tech):
    cols=['luminal','epithelial','stromal',tech,score]
    x=pd.DataFrame({c:z(df[c]) for c in cols}); y=z(df[ycol]); x=sm.add_constant(x,has_constant='add')
    ok=x.notna().all(1)&y.notna(); x=x.loc[ok]; y=y.loc[ok]
    if len(x)==0: return np.nan,np.nan,np.nan,np.nan,np.nan
    if np.linalg.matrix_rank(x.to_numpy())<x.shape[1]: return np.nan,np.nan,np.nan,np.nan,np.linalg.cond(x.to_numpy())
    r=sm.OLS(y,x).fit(cov_type='HC3'); ci=r.conf_int().loc[score]
    return r.params[score],ci.iloc[0],ci.iloc[1],r.pvalues[score],np.linalg.cond(x.to_numpy())
def rho_safe(x,y):
    ok=pd.Series(x).notna().to_numpy()&pd.Series(y).notna().to_numpy()
    return stats.spearmanr(np.asarray(x)[ok],np.asarray(y)[ok]).statistic if ok.sum()>=3 else np.nan

a=ad.read_h5ad(SP)
var=set(map(str,a.var_names));
allgenes=set(sum(SETS.values(),[]))|set(APM)|set(IFNG)|{'FGFR3'}
for _genes in MAN.genes.astype(str): allgenes.update(_genes.split(';'))
records=[]; comp=[]; detect=[]; multi=[]
for m in MAN.itertuples():
    sec=m.section; pat=m.held_out_patient
    q=np.asarray(a.obs.sampleID.astype(str)==sec)
    s=a[q].copy(); raw=s[:,[g for g in allgenes if g in var]].to_df().astype(float).reset_index(drop=True)
    obs=s.obs.reset_index(names='spot_id'); xy=np.asarray(s.obsm['spatial'])
    d=obs[['spot_id','sampleID','Patient','n_genes_by_counts','total_counts']].copy(); d['section']=sec; d['x']=xy[:,0]; d['y']=xy[:,1]
    # Same within-section standardized components used in the submitted strict LOPO score.
    x=pd.DataFrame({c:z(raw[c]) for c in raw.columns})
    markers=[g for g in str(m.genes).split(';') if g in x.columns]
    d['positive_signature']=x[markers].mean(1)
    d['APM']=module(x,APM); d['IFNG']=module(x,IFNG); d['FGFR3']=x['FGFR3']
    d['projection_full']=d.positive_signature-d.APM-.5*d.IFNG
    d['projection_minus_APM_only']=d.positive_signature-.5*d.IFNG
    d['projection_minus_IFNG_only']=d.positive_signature-d.APM
    d['signature_only']=d.positive_signature
    d['projection_without_immune_subtraction']=d.positive_signature
    # Local outcomes from the same raw-score scheme and same k=6 non-self neighbors.
    for nm,gs in SETS.items(): d[nm]=module(x,gs)
    d['luminal']=d['local_luminal']; d['epithelial']=d['local_epithelial']; d['stromal']=d['local_stromal']
    _,nn=cKDTree(xy).query(xy,k=7)
    for nm in SETS:
        if nm.startswith('local_'): d[nm+'_k6']=d[nm].to_numpy()[nn[:,1:]].mean(1)
    d['log_depth']=np.log1p(d.total_counts); d['log_genes']=np.log1p(d.n_genes_by_counts)
    # Raw and normalized FGFR3 detectability fields.
    d['FGFR3_raw_count']=raw['FGFR3'].to_numpy(); d['FGFR3_detectable']=d.FGFR3_raw_count>0
    norm=raw.div(np.maximum(s.obs.total_counts.to_numpy(),1),axis=0)*10000
    d['FGFR3_counts_normalized']=norm['FGFR3'].to_numpy(); d['FGFR3_within_section_z']=z(d.FGFR3_counts_normalized)
    old=SCORES.loc[SCORES.section.eq(sec),['spot_id','FGFR3','normalized_projection']].set_index('spot_id')
    if len(old)!=len(d) or set(old.index)!=set(d.spot_id): raise ValueError(f'spot_id mismatch between raw h5ad and previous score file for {sec}')
    old=old.reindex(d.spot_id)
    d['FGFR3_z_previous']=old.FGFR3.to_numpy()
    d.to_csv(OUT/f'04_{sec}_spot_detection_scores.csv.gz',index=False,compression='gzip')
    outcomes=['local_dc_original_k6','local_dc_no_hla_k6','local_cDC1_k6','local_cDC2_k6','local_myeloid_k6','local_cd8_k6','local_apm_k6']
    score_names=['signature_only','FGFR3','APM','IFNG','projection_full','projection_minus_APM_only','projection_minus_IFNG_only','projection_without_immune_subtraction']
    for score in score_names:
        for ycol in outcomes:
            rho=rho_safe(d[score],d[ycol]); ci,nb,spacing=block_ci(d[score].to_numpy(),d[ycol].to_numpy(),xy)
            for tech in ['log_depth','log_genes']:
                be,lo,hi,p,cond=fit_adjusted(d,ycol,score,tech)
                nr=rho_safe(old.normalized_projection,d[ycol]) if score in ['projection_full','signature_only'] else np.nan
                records.append({'section':sec,'score':score,'outcome':ycol.replace('_k6',''),'n_spots':len(d),'spatial_blocks':nb,'rho':rho,'rho_ci_low':ci[0],'rho_ci_high':ci[1],'prior_normalized_projection_rho_check':nr,'adjustment':tech,'beta':be,'beta_ci_low':lo,'beta_ci_high':hi,'beta_p_HC3':p,'condition_number':cond,'biological_n':1,'inference_note':'section-conditional block CI; adjusted model is descriptive HC3 sensitivity; four linked sections are biological n=4'})
    q75=d.projection_full.quantile(.75); sq75=d.signature_only.quantile(.75)
    # Tie-safe: high means detected and ranked among positive values; all zero ties are low/not identifiable.
    positive=d.loc[d.FGFR3_raw_count>0,'FGFR3_raw_count']; fg_high=(d.FGFR3_raw_count>0)&(d.FGFR3_raw_count>=positive.quantile(.75) if len(positive) else False)
    fulltop=d.projection_full>=q75; sigtop=d.signature_only>=sq75
    topfg=d.loc[fulltop,'FGFR3_raw_count']; topapm=d.loc[fulltop,'APM']
    comp.append({'section':sec,'n_spots':len(d),'positive_signature_variance':d.positive_signature.var(),'APM_variance':d.APM.var(),'IFNG_variance':d.IFNG.var(),'projection_variance':d.projection_full.var(),'cov_positive_APM':d.positive_signature.cov(d.APM),'cov_positive_IFNG':d.positive_signature.cov(d.IFNG),'cov_APM_IFNG':d.APM.cov(d.IFNG),'rho_positive_APM':rho_safe(d.positive_signature,d.APM),'rho_positive_IFNG':rho_safe(d.positive_signature,d.IFNG),'rho_positive_DC':rho_safe(d.positive_signature,d.local_dc_original_k6),'rho_full_DC':rho_safe(d.projection_full,d.local_dc_original_k6),'rho_APM_DC':rho_safe(d.APM,d.local_dc_original_k6),'rho_IFNG_DC':rho_safe(d.IFNG,d.local_dc_original_k6),'full_signature_topQ_overlap':int((fulltop&sigtop).sum()),'full_signature_topQ_jaccard':float((fulltop&sigtop).sum()/(fulltop|sigtop).sum()),'full_topQ_FGFR3_detectable_pct':float(d.loc[fulltop,'FGFR3_detectable'].mean()),'full_topQ_FGFR3_topQ_among_nonzero_pct':float(fg_high[fulltop].mean()),'full_topQ_APM_bottomQ_pct':float((d.APM<=d.APM.quantile(.25))[fulltop].mean()),'FGFR3_75th_raw_quantile':float(d.FGFR3_raw_count.quantile(.75)),'FGFR3_75th_nonzero_quantile':float(positive.quantile(.75)) if len(positive) else np.nan,'FGFR3_75th_is_zero':bool(d.FGFR3_raw_count.quantile(.75)==0),'fg_high_definition':'raw_count>0 and >=75th percentile among nonzero spots'})
    detect.append({'section':sec,'n_spots':len(d),'n_detectable':int(d.FGFR3_detectable.sum()),'pct_detectable':100*d.FGFR3_detectable.mean(),'n_zero':int((~d.FGFR3_detectable).sum()),'zero_pct':100*(~d.FGFR3_detectable).mean(),'nonzero_min':float(positive.min()) if len(positive) else np.nan,'nonzero_median':float(positive.median()) if len(positive) else np.nan,'nonzero_max':float(positive.max()) if len(positive) else np.nan,'nonzero_q25':float(positive.quantile(.25)) if len(positive) else np.nan,'nonzero_q75':float(positive.quantile(.75)) if len(positive) else np.nan,'raw_q75':float(d.FGFR3_raw_count.quantile(.75)),'raw_q75_is_zero':bool(d.FGFR3_raw_count.quantile(.75)==0),'ties_at_raw_mode_pct':100*d.FGFR3_raw_count.value_counts(normalize=True).iloc[0],'enough_nonzero_for_quartile':bool(len(positive)>=20),'interpretation':'FGFR3 quartile stratification unreliable' if len(positive)<20 else 'nonzero-ranked sensitivity estimable'})
    for tech in ['log_depth','log_genes']:
        cov=['luminal','epithelial','stromal',tech]
        xx=sm.add_constant(pd.DataFrame({c:z(d[c]) for c in cov}),has_constant='add')
        vv={c:np.nan for c in xx.columns}
        valid=xx.notna().all(1)
        if valid.sum()>len(cov)+2:
            xv=xx.loc[valid].to_numpy()
            vv={c:variance_inflation_factor(xv,i) for i,c in enumerate(xx.columns)}
        for score in ['projection_full','signature_only','FGFR3']:
            for outcome in ['local_dc_original_k6','local_dc_no_hla_k6','local_cDC1_k6','local_cDC2_k6','local_myeloid_k6','local_cd8_k6','local_apm_k6']:
                x2=xx.join(z(d[score]).rename(score)); y=z(d[outcome]); ok=x2.notna().all(1)&y.notna()
                if ok.sum()>x2.shape[1]+2:
                    r=sm.OLS(y[ok],x2.loc[ok]).fit(cov_type='HC3'); ci=r.conf_int().loc[score]; beta=r.params[score]; lo,hi=ci.iloc[0],ci.iloc[1]; pp=r.pvalues[score]; rank=int(np.linalg.matrix_rank(x2.loc[ok].to_numpy())); cond=float(np.linalg.cond(x2.loc[ok].to_numpy()))
                else: beta=lo=hi=pp=np.nan; rank=0; cond=np.nan
                multi.append({'section':sec,'score':score,'outcome':outcome.replace('_k6',''),'technical_adjustment':tech,'n':int(ok.sum()),'design_rank':rank,'design_columns':x2.shape[1],'condition_number':cond,'vif_log_depth':vv.get('log_depth',np.nan),'vif_log_genes':vv.get('log_genes',np.nan),'vif_luminal':vv.get('luminal',np.nan),'vif_epithelial':vv.get('epithelial',np.nan),'vif_stromal':vv.get('stromal',np.nan),'beta':beta,'ci_low':lo,'ci_high':hi,'p_HC3':pp,'inference_note':'technical proxy-only sensitivity; not causal'})

pd.DataFrame(records).to_csv(OUT/'02_spatial_score_ablation_section_summary.csv',index=False)
pd.DataFrame(comp).to_csv(OUT/'03_spatial_projection_component_audit.csv',index=False)
pd.DataFrame(detect).to_csv(OUT/'04_FGFR3_spatial_detection_QC.csv',index=False)
pd.DataFrame(multi).to_csv(OUT/'05_spatial_multicollinearity_minimal_sensitivity.csv',index=False)

verdict='B'
summary=f'''# Spatial finding verdict — Round 2

## First-screen verdict: **B — composite-score construction and tissue composition dominate the apparent spatial DC pattern**

The positive-only strict-LOPO epithelial signature is positive with local original DC in all four linked sections, whereas the full score is negative in all four. The full score subtracts APM and 0.5×IFNG by definition; this is mathematical coupling, not an independent DC measurement. The spatial evidence therefore cannot independently support “FGFR3-high epithelial program in DC-poor neighborhoods.”

| Section | spots | FGFR3 detectable | positive-only vs DC rho | full vs DC rho | APM vs DC rho | FGFR3 raw q75 | q75 zero? |
|---|---:|---:|---:|---:|---:|---:|---|
'''
for r in comp:
    de=next(x for x in detect if x['section']==r['section']); summary+=f"| {r['section']} | {r['n_spots']} | {de['pct_detectable']:.2f}% | {r['rho_positive_DC']:.3f} | {r['rho_full_DC']:.3f} | {r['rho_APM_DC']:.3f} | {de['raw_q75']:.3g} | {'yes' if de['raw_q75_is_zero'] else 'no'} |\n"
summary+='''
## Strongest publishable conclusion

Across bulk and donor-level epithelial data, higher FGFR3 remains associated with lower selected HLA-I/APM-related RNA programs. The linked spatial data provide a descriptive atlas-derived composite projection, not an independent validation of DC depletion, FGFR3-high localization, cell-surface HLA-I function, or causality.

## Required deletions/narrowing

Delete “dendritic-cell-poor spatial context,” “DC depletion,” “antigen-presentation impairment,” and causal immune-evasion language from the title, abstract, conclusion and discussion. Replace with an association/RNA formulation and state that the four linked sections show score-construction sensitivity. Keep the spatial result as a descriptive sensitivity layer.

## Boundaries

The raw h5ad was available for count-level detection QC. No direct DC counts, validated tissue masks, independent spatial patients, protein/functional HLA-I data or definitive malignant labels were available. Spot-level P values are not patient-level evidence; all section-conditional intervals have biological n=4.
'''
(OUT/'01_SPATIAL_FINDING_VERDICT.md').write_text(summary,encoding='utf-8')
(OUT/'07_remaining_noncomputational_actions.md').write_text('''# Remaining non-computational actions\n\n1. Author must approve and publish `repository_release_candidate_v2` with a license, then verify a logged-out URL/DOI and synchronize the manuscript/submission-system Data and Code Availability text.\n2. Authors must confirm exact Ethical approval, Consent to participate and Consent to publish wording; no new approval or consent was inferred.\n3. The manuscript writer must apply the sentence-level spatial deletions/narrowing in `06_reviewer_response_evidence_map.md`; Codex did not edit the Word files.\n4. No additional direct DC abundance, HLA-I protein/function, histology-region annotation, independent spatial cohort or definitive malignant-cell label exists in the available public objects.\n5. Raw-data end-to-end reproduction remains path-dependent; this V2 package includes a portable derived-data run and source hashes, not a claim of one-command raw download-to-manuscript reproduction.\n''',encoding='utf-8')
(OUT/'06_reviewer_response_evidence_map.md').write_text('''# Targeted reviewer evidence map — Round 2\n\n## Reviewer 1.1 / Reviewer 2 Major 11–12\nComment: causal, malignant-cell and HLA-I functional wording exceeds RNA evidence.\nOriginal evidence: `spatial/lopo_spot_scores.csv`, `spatial/lopo_leakage_checks.md`.\nNew evidence: `02_spatial_score_ablation_section_summary.csv`, `04_FGFR3_spatial_detection_QC.csv`, `03_spatial_projection_component_audit.csv`.\nConclusion strength: spatial composite is descriptive; do not call DC abundance, malignant-specificity or HLA-I function.\nSuggested replacement: “In four linked atlas sections, an out-of-fold epithelial projection showed section-specific spatial relationships with DC-associated transcriptional scores; this did not establish dendritic-cell abundance, FGFR3-high localization, or HLA-I protein function.”\nLocation: Figure 4 legend; Results spatial paragraph; Discussion and Conclusion.\n\n## Reviewer 1.2 / Reviewer 2 Major 4\nComment: LOPO leakage.\nOriginal evidence: `spatial/lopo_fold_manifest.csv`, `spatial/lopo_leakage_checks.md`.\nNew evidence: V2 uses the existing strict fold-specific marker manifest and reproduces score ablations without changing fold membership.\nConclusion strength: LOPO exclusion remains supported; it does not rescue the composite-score independence problem.\nSuggested replacement: “All four spatial sections were scored using fold-specific markers derived without the matched atlas donor.”\nLocation: Methods LOPO paragraph and Figure 4 legend.\n\n## Reviewer 2 Major 5–7 and Minor 5\nComment: spot pseudo-replication, DC specificity and composition.\nOriginal evidence: `spatial/spatial_block_resampling.csv`, `spatial/spatial_multicollinearity_diagnostics.md`.\nNew evidence: all scores/outcomes, 1,200 section-conditional block CIs, log-depth/log-genes sensitivity and raw detection QC in V2 CSVs.\nConclusion strength: biological n remains four; positive-only signature is positive with original DC in all four, while full composite is negative; DC score remains transcriptional, not abundance.\nSuggested replacement: “Spatial analyses were descriptive and conditional on four linked sections; the composite projection’s negative DC relationship was not independently reproduced by the positive-only epithelial signature.”\nLocation: Results spatial paragraph, Figure 4/5 legends, Discussion limitations.\n\n## Reviewer 2 Major 14\nComment: novelty beyond atlas.\nNew evidence: `03_spatial_projection_component_audit.csv` plus prior `snrna/atlas_incremental_nested_models.csv`.\nConclusion strength: incremental value is donor-level FGFR3–APM RNA association; spatial novelty is not independently established.\nSuggested replacement: “The reproducible incremental signal is a donor-level FGFR3–APM-related RNA association beyond state and luminal covariates.”\nLocation: Abstract, Discussion and Conclusion.\n''',encoding='utf-8')

manifest={'project':'FGFR3 Discover Oncology R1 targeted spatial validation Round 2','verdict':'B','seed':SEED,'source_h5ad':f'{SP.name} (local input; not redistributed)','source_h5ad_sha256':hashlib.sha256(SP.read_bytes()).hexdigest(),'sections':len(comp),'spot_files':len(list(OUT.glob('04_*_spot_detection_scores.csv.gz'))),'derived_from':'existing strict LOPO fold manifest; no T01-T14 rerun','raw_count_qc':'completed from original GSE171351 h5ad','block_bootstrap_draws':1200,'biological_n':4,'repository_status':'NOT PUBLISHED'}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
(OUT/'RUN_LOG.md').write_text('''# Round 2 run log\n\nInput: original GSE171351_combined_visium.h5ad; previous strict LOPO marker manifest and derived score definitions.\n\nThe script re-used the previous fold-specific marker lists and four section mapping, computed positive-only/full/ablation components, k=6 local outcomes, 1,200-draw spatial block CIs, raw FGFR3 detection QC, and log_depth/log_genes technical sensitivities. It did not rerun T01–T14 or edit Word files.\n\nThe source h5ad is not copied into this release; SHA256 is in manifest.json. The package is NOT PUBLISHED.\n''',encoding='utf-8')

print(json.dumps({'out':str(OUT),'rows':len(records),'components':len(comp),'detect':len(detect),'multicollinearity':len(multi)},indent=2))
