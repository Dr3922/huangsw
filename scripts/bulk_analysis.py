from reanalysis import *
from scipy.stats import mannwhitneyu
import csv, gzip
GEO=ROOT/'geo_validation/processed'
SETS={'Luminal-papillary':LUM,'Immune':['CD274','PDCD1','CTLA4','CD8A','CXCL9','CXCL10','STAT1','IFNG','HLA-DRA','LAG3','TIGIT'],'Basal':['KRT5','KRT6A','KRT14','TP63','EGFR','CD44','S100A2','VIM','KRT17'],'Stromal/TGFb':['TGFB1','COL1A1','COL1A2','COL3A1','ACTA2','TAGLN','FAP'],'Antigen presentation':['B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2'],'Cytolytic':['GZMA','PRF1','NKG7'],'epithelial_proxy':['EPCAM','TACSTD2','KRT8','KRT18','KRT19']}
OUTCOMES=list(SETS)[:-1]
def effect(a,b):
 a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float);a=a[np.isfinite(a)];b=b[np.isfinite(b)]
 n,m=len(a),len(b)
 if min(n,m)<2:return dict(estimate=np.nan,se=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan)
 sd=np.sqrt(((n-1)*a.var(ddof=1)+(m-1)*b.var(ddof=1))/(n+m-2));d=(a.mean()-b.mean())/sd;se=np.sqrt((n+m)/(n*m)+d*d/(2*(n+m-2)))
 return dict(estimate=d,se=se,ci_low=d-1.96*se,ci_high=d+1.96*se,p=mannwhitneyu(a,b,alternative='two-sided').pvalue)
def cohend_only(a,b):
 a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float);n=len(a);m=len(b)
 if min(n,m)<2:return np.nan
 sd=np.sqrt(((n-1)*a.var(ddof=1)+(m-1)*b.var(ddof=1))/(n+m-2))
 return (a.mean()-b.mean())/sd if sd>0 else np.nan
def meta(d):
 d=d.dropna(subset=['estimate','se']);d=d[d.se>0];k=len(d)
 if not k:return dict(k=0)
 v=d.se.to_numpy()**2;w=1/v;y=d.estimate.to_numpy();mu=np.sum(w*y)/sum(w);Q=np.sum(w*(y-mu)**2);c=sum(w)-sum(w*w)/sum(w);tau=max(0,(Q-k+1)/c) if k>1 else 0;ww=1/(v+tau);est=sum(ww*y)/sum(ww);se=np.sqrt(1/sum(ww))
 return dict(k=k,estimate=est,se=se,ci_low=est-1.96*se,ci_high=est+1.96*se,p=2*stats.norm.sf(abs(est/se)),tau2=tau,I2=100*max(0,(Q-k+1)/Q) if Q>0 else 0,Q=Q,Q_df=k-1,Q_p=stats.chi2.sf(Q,k-1) if k>1 else np.nan,estimator='DerSimonian-Laird',interval='normal Wald; k=3 unstable')
def parse_matrix_metadata(path):
 vals={}
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt',encoding='utf-8') as f:
  for line in f:
   if line.startswith('!series_matrix_table_begin'):break
   if line.startswith('!Sample_'):
    row=next(csv.reader([line],delimiter='\t'));vals.setdefault(row[0][8:],[]).append(row[1:])
 return vals
def ledger135(m):
 rows=[]
 fresh=OUT/'public_evidence/GSE13507_series_matrix.txt.gz'
 vals=parse_matrix_metadata(fresh if fresh.exists() else ROOT/'geo_validation/raw/GSE13507_series_matrix.txt.gz')
 gsm=vals['geo_accession'][0];fresh_titles=dict(zip(gsm,vals['title'][0]));fresh_source=dict(zip(gsm,vals['source_name_ch1'][0]))
 for g,r in m.iterrows():
  title=fresh_titles.get(g,str(r.title));source=fresh_source.get(g,str(r.source_name_ch1));bt=re.search(r'BT(\d+)(?:\.(\d+))?',title);pid='BT'+bt.group(1) if bt else 'unknown';spec='primary' if source=='Primary bladder cancer' else ('recurrent' if 'Recurrent' in source else 'non_tumor_or_conflicting')
  rows.append(dict(gsm=g,source_id=title,patient_id=pid,specimen_type=source,primary_vs_recurrent=spec,tumor_vs_normal='tumor' if spec in ['primary','recurrent'] else ('unknown' if 'MPP' in source else 'non_tumor'),stage=r.get('char_stage','unknown'),grade=r.get('char_grade','unknown'),molecular_context=r.get('char_invasiveness','unknown'),annotation_source='fresh GEO series matrix title/source; GSE13507 overall_design states 23 recurrent tissues from 14 patients',replicate_status='repeated patient' if spec=='recurrent' else ('index primary' if spec=='primary' else 'unknown'),primary_eligibility=spec=='primary',exclusion_reason='' if spec=='primary' else ('recurrent specimen from indexed patient' if spec=='recurrent' else 'non-tumor or conflicting control annotation'),ambiguity='BT prefix is published specimen identifier; same prefix establishes recurrence linkage. No cross-link assumed for BS/C controls.' if bt else 'patient relationship unknown; excluded from tumor analysis',cached_title=r.title,fresh_title=title))
 l=pd.DataFrame(rows).fillna('unknown');save(l,'cohort_ledgers/gse13507_sample_ledger.csv')
 primary=set(l.loc[l.primary_eligibility,'patient_id']);rec=l[l.primary_vs_recurrent=='recurrent']
 assert len(primary)==165 and rec.patient_id.nunique()==15 and set(rec.patient_id)<=primary
 note('cohort_ledgers/GSE13507_patient_identifier_discrepancy.md','GEO overall_design states 14 recurrent patients, but the 23 published recurrent titles contain 15 unique BT prefixes, all present among 165 primary titles. BT-prefix clusters are reproducible specimen-identifier links, not verified clinical patient identities. Primary165 is the main analysis; repeated-specimen sensitivity clusters by the 15 BT prefixes. Public metadata inconsistency is unresolved and must be disclosed.')
 save(rec[['gsm','source_id','patient_id']],'cohort_ledgers/gse13507_recurrent_patient_links.csv')
 return l
def run():
 rng=np.random.default_rng(SEED);results=[];cutrows=[];ctx=[];adj=[];scores=[];counts=[];ledgers=[];coverage=[];diag=[]
 for acc in ['GSE32894','GSE13507','GSE31684']:
  raw=pd.read_csv(GEO/(acc+'_gene_expression_raw.csv'),index_col=0);m=pd.read_csv(GEO/(acc+'_sample_metadata.csv'),index_col=0)
  if acc=='GSE13507':
   l=ledger135(m).set_index('gsm');m['patient_id']=l.patient_id;m['context']=m.char_invasiveness.fillna('unknown');subsets={'submitted_all256':m.index,'primary_tumor_only':l.index[l.primary_eligibility],'all_tumors_clustered':l.index[l.tumor_vs_normal=='tumor']}
  else:
   m['patient_id']=m.index if acc=='GSE31684' else m.source_name_ch1.str.extract(r'(UC_\d+)',expand=False)
   m['context']=m['char_cluster' if acc=='GSE31684' else 'char_molecular_subtype'].fillna('unknown').astype(str)
   subsets={'submitted_all':m.index,'primary_analysis':m.index}
   # UC prefixes cannot be assumed independent if repeated. Audit explicitly.
  ledgers.append(pd.DataFrame(dict(cohort=acc,sample_id=m.index,source_id=m.title,patient_id=m.patient_id,context=m.context,specimen_source=m.source_name_ch1,eligibility='see subset counts',patient_id_evidence='GEO title/source identifier; GSE31684 one cystectomy sample per published patient')))
  for subset,ids in subsets.items():
   mt=m.loc[ids].copy();n=len(mt);x=z(raw.loc[:,ids].T);d=pd.DataFrame(index=mt.index);d['FGFR3']=x.FGFR3
   for name,gs in SETS.items():
    gg=[g for g in gs if g in x];d[name]=x[gg].mean(axis=1);coverage.append(dict(cohort=acc,subset=subset,module=name,requested=len(gs),present=len(gg),genes=';'.join(gg),missing=';'.join(set(gs)-set(gg))))
   d['context']=mt.context;d['patient_id']=mt.patient_id;save(d.reset_index(names='sample_id'),f'bulk/{acc}_{subset}_restandardized_scores.csv');scores.append(d.reset_index(names='sample_id').assign(cohort=acc,subset=subset))
   grouped=subset=='all_tumors_clustered' or mt.patient_id.duplicated().any()
   cluster_ix=None
   if grouped:
    pats=sorted(d.patient_id.unique());patient_ix=[np.flatnonzero(d.patient_id.to_numpy()==p) for p in pats]
    cluster_ix=[np.concatenate([patient_ix[j] for j in rng.integers(0,len(pats),len(pats))]) for _ in range(2000)]
   known_patient=mt.patient_id.notna() & ~mt.patient_id.astype(str).isin(['unknown','NA',''])
   counts.append(dict(cohort=acc,subset=subset,samples=n,unique_patient_ids=mt.loc[known_patient,'patient_id'].nunique(),patient_mapping_unknown=int((~known_patient).sum()),excluded=len(m)-n,normalization='gene z scores recomputed on this subset',inference='patient cluster' if grouped else 'sample; one published index sample per patient'))
   for out in OUTCOMES:
    rr,pp=stats.spearmanr(d.FGFR3,d[out]);idx=rng.integers(0,n,(5000,n));boot,_=rankrho(d.FGFR3.to_numpy()[idx],d[out].to_numpy()[idx]);ss=bootstrap_summary(boot,rr);ss['p']=pp
    if grouped:
     # Patient-level resampling respects multiplicity and preserves repeated specimens.
     draw=[]
     for ix in cluster_ix:draw.append(stats.spearmanr(d.FGFR3.iloc[ix],d[out].iloc[ix]).statistic)
     ss=bootstrap_summary(np.array(draw),rr)
    results.append(dict(cohort=acc,subset=subset,outcome=out,model='continuous_spearman',n=n,**ss))
    for q in np.round(np.arange(.5,.901,.05),2):
     threshold=d.FGFR3.quantile(q);high=d.FGFR3>=threshold;ef=effect(d.loc[high,out],d.loc[~high,out]);row=dict(cohort=acc,subset=subset,outcome=out,cutoff=q,threshold=threshold,high_n=int(high.sum()),low_n=int((~high).sum()),n=n,**ef)
     if grouped:
      boots=[]
      xv=d.FGFR3.to_numpy();yv=d[out].to_numpy()
      for ix in cluster_ix:
       xx=xv[ix];yy=yv[ix];hh=xx>=np.quantile(xx,q);boots.append(cohend_only(yy[hh],yy[~hh]))
      if boots:
       bs=bootstrap_summary(np.array(boots),ef['estimate']);row.update(bs);row['se']=np.nanstd(boots,ddof=1)
     cutrows.append(row)
     if q==.75:results.append(dict(row,model='Cohens_d_q75'))
    xx=d[['FGFR3','epithelial_proxy']].join(d[['Stromal/TGFb']]);dummy=pd.get_dummies(d.context.replace('unknown',np.nan),drop_first=True,dtype=float)
    ok=d.context!='unknown';xx=xx.join(dummy).loc[ok];res=fit(d.loc[ok,out],xx,d.loc[ok,'patient_id'] if grouped else None)
    adj.append(dict(cohort=acc,subset=subset,outcome=out,model='FGFR3 + epithelial proxy + stromal + context',**res))
    for stratum,dd in d[d.context!='unknown'].groupby('context'):
     res=fit(dd[out],dd[['FGFR3']],dd.patient_id if grouped else None);ctx.append(dict(cohort=acc,subset=subset,outcome=out,stratum=stratum,**res))
    mm=sm.OLS(d[out],sm.add_constant(d[['FGFR3']])).fit();diag.append(dict(cohort=acc,subset=subset,outcome=out,shapiro_residual_p=stats.shapiro(mm.resid).pvalue,breusch_pagan_p=sm.stats.diagnostic.het_breuschpagan(mm.resid,mm.model.exog)[1],max_cooks_d=float(mm.get_influence().cooks_distance[0].max()),uncertainty='HC3 or donor-cluster sandwich; correlation bootstrap; non-normality does not invalidate rank correlation'))
  print('bulk complete',acc,flush=True)
 results=bh(pd.DataFrame(results),['cohort','subset','model']);cut=pd.DataFrame(cutrows)
 save(results,'bulk/bulk_cohort_effects.csv');save(bh(pd.DataFrame(ctx),['cohort','subset','stratum']),'bulk/bulk_within_context.csv');save(bh(pd.DataFrame(adj),['cohort','subset']),'bulk/bulk_proxy_adjusted.csv');save(cut,'bulk/bulk_cutoff_sensitivity.csv');save(pd.concat(ledgers),'cohort_ledgers/bulk_patient_sample_ledger.csv');save(pd.DataFrame(counts),'bulk/bulk_inclusion_counts.csv');save(pd.concat(scores),'figure_data/figure2_sample_scores.csv');save(pd.DataFrame(coverage),'bulk/bulk_gene_coverage.csv');save(pd.DataFrame(diag),'bulk/bulk_model_diagnostics.csv')
 save(results[results.cohort=='GSE13507'],'bulk/gse13507_recomputed_results.csv');save(cut[cut.cohort=='GSE13507'],'bulk/gse13507_sensitivity.csv')
 metas=[];loco=[]
 for family,selector in [('submitted',{'GSE32894':'submitted_all','GSE13507':'submitted_all256','GSE31684':'submitted_all'}),('revised',{'GSE32894':'primary_analysis','GSE13507':'primary_tumor_only','GSE31684':'primary_analysis'}),('all_tumors_sensitivity',{'GSE32894':'primary_analysis','GSE13507':'all_tumors_clustered','GSE31684':'primary_analysis'})]:
  eligible=cut[[s==selector[a] for a,s in zip(cut.cohort,cut.subset)]]
  for (out,q),dd in eligible.groupby(['outcome','cutoff']):
   metas.append(dict(analysis=family,outcome=out,cutoff=q,**meta(dd)))
   if q==.75:
    for acc in selector:loco.append(dict(analysis=family,outcome=out,omitted=acc,**meta(dd[dd.cohort!=acc])))
 save(bh(pd.DataFrame(metas),['analysis','cutoff']),'bulk/bulk_meta_results.csv');save(pd.DataFrame(loco),'bulk/bulk_meta_leave_one_cohort_out.csv')
 old=results[results.subset.str.startswith('submitted')];new=results[results.subset.isin(['primary_analysis','primary_tumor_only'])];compare=old.merge(new,on=['cohort','outcome','model'],suffixes=('_submitted','_revised'));save(compare,'bulk/bulk_original_vs_revised.csv')
 note('bulk/bulk_models_methods.md','''# Bulk statistical methods
Source processed log-expression matrices are read without using historical gene z scores. Each analysis subset is independently gene-standardized (sample SD, ddof=1); modules average available S1 genes. GSE13507 primary165 is main, all188 tumors sensitivity clusters by BT prefix. GEO's overall design states that the23 recurrent tissues came from14 patients, whereas the published recurrent titles contain15 unique BT prefixes; this unresolved metadata discrepancy is disclosed and the BT prefixes are treated only as reproducible specimen-linkage clusters. All256 is historical audit only, includes68 non-tumor/discordant-control samples. No false claim that a mouse expression sample was included: GSM340542 annotation conflicts are quarantined and source logged. Recurrent BT034.1 links to primary BT034 by published sample identifiers. The one-index-primary rule is deterministic.
Spearman rho with5000 percentile bootstrap intervals (2000 patient-cluster draws in repeated-sample sensitivity). Original Mann-Whitney P accompanies pooled-SD Cohen d; d CI uses submitted large-sample SE except clustered sensitivity uses2000 donor draws, cutoff recomputed each draw. Context models use HC3 sandwich/t intervals; full models add epithelial proxy, stromal/TGFb and measured context dummy variables, excluding unknown context. No claim of complete confounding removal. Missing/rank-deficient strata are NA, not zero. Meta uses submitted DerSimonian-Laird moment tau2 and normal intervals, Q,I2; three cohorts cannot reliably estimate heterogeneity. LOCO and cutoff.50-.90 are sensitivities in the same samples. BH across six module outcomes per cohort/subset/model (or context stratum; meta per cutoff), not across exploratory cutoffs. Gene coverage and assumptions diagnostics supplied. All-tumor clustered sensitivity percentile P is approximate centered bootstrap, not independent-sample Mann-Whitney inference.
''')
if __name__=='__main__':run()
