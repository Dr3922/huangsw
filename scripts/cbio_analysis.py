from bulk_analysis import effect
from reanalysis import *
from statsmodels.stats.proportion import proportion_confint
def read(p):return json.loads(p.read_text(encoding='utf-8')) if p.exists() else []
def run():
 inv=pd.read_csv(ROOT/'datasets/cbioportal_fgfr3_inventory.csv');rows=[];ledger=[];sets={};recon=[]
 excluded=['Targeted profiling','Preclinical','Variant context']
 print(inv.groupby('role').sampleCount.agg(['size','sum']).to_string(),flush=True)
 for _,r in inv.iterrows():
  sid=r.studyId;p=OUT/'public_evidence/cbio'/sid;samples=read(p/'samples.json');lists=read(p/'lists.json');muts=read(p/'FGFR3_mutations.json');cna=read(p/'FGFR3_cna_gene_specific.json');info=read(p/'study.json')
  if not isinstance(samples,list) or not samples:continue
  eligible=set();elig_cna=set()
  for li in lists:
   ss=set(read(p/(li['sampleListId']+'.json')))
   if li.get('category')=='all_cases_with_mutation_data':eligible|=ss
   if li.get('category')=='all_cases_with_cna_data':elig_cna|=ss
  # No mutation record is WT only inside a source-confirmed sequenced sample list.
  # The discrete-CNA endpoint can ignore a singular gene parameter and return all genes.
  # Always enforce the requested gene in the returned records, not only in the URL.
  muts=[m for m in muts if m.get('entrezGeneId',m.get('gene',{}).get('entrezGeneId'))==2261]
  cna=[m for m in cna if m.get('entrezGeneId',m.get('gene',{}).get('entrezGeneId'))==2261]
  # A study-level CNA sample list can span profiles or contain samples without a
  # returned value for this molecular profile. Only a returned FGFR3 value is
  # an observed denominator; absent records are unknown, not diploid.
  cna_observed={m['sampleId'] for m in cna if m.get('value') is not None}
  elig_cna &= cna_observed
  allmut={m['sampleId'] for m in muts};coding={m['sampleId'] for m in muts if m.get('mutationType') not in ['Silent','Intron','IGR',"3'UTR","5'UTR","3'Flank","5'Flank"]};protein={m['sampleId'] for m in muts if m.get('proteinChange') not in [None,'','p.?','NA']};amp={m['sampleId'] for m in cna if m.get('value')==2};gain={m['sampleId'] for m in cna if m.get('value')==1}
  patientids=set()
  for s in samples:
   smp=s['sampleId'];pid=s['patientId'];canonical=pid[:12] if pid.startswith('TCGA-') else pid;patientids.add(canonical)
   ledger.append(dict(study_id=sid,sample_id=smp,patient_id=pid,canonical_identifier=canonical,identifier_scope='TCGA global barcode' if pid.startswith('TCGA-') else 'study-specific; cross-study identity not validated',mutation_profiled=smp in eligible,FGFR3_any_mutation=smp in allmut,FGFR3_coding_mutation=smp in coding,FGFR3_protein_change=smp in protein,cna_profiled=smp in elig_cna,amplified=smp in amp,gain=smp in gain,role=r.role))
  sets[sid]=patientids
  n=len(eligible);mut=len(coding&eligible);ci=proportion_confint(mut,n,method='wilson') if n else (np.nan,np.nan)
  rows.append(dict(study_id=sid,study_name=r.studyName,role=r.role,acquisition_date='2026-10-21',source_version=info.get('citation','') if isinstance(info,dict) else '',submitted_sample_count=r.sampleCount,current_total_samples=len(samples),unique_study_patient_ids=len(patientids),mutation_eligible_n=n,coding_mutant_n=mut,any_variant_n=len(allmut&eligible),protein_change_n=len(protein&eligible),unknown_mutation_status_n=len(samples)-n,frequency=mut/n if n else np.nan,ci_low=ci[0],ci_high=ci[1],cna_eligible_n=len(elig_cna),amplified_n=len(amp&elig_cna),gain_n=len(gain&elig_cna),source='https://www.cbioportal.org/study/summary?id='+sid,independence='overlapping TCGA releases' if 'tcga' in sid else 'cross-study overlap incompletely identifiable'))
  # Submitted plot roles excludes preclinical/variant but includes all five targeted studies.
  plotted=not ('Preclinical' in r.role or 'Variant' in r.role)
  tab=plotted and 'Targeted' not in r.role
  recon.append(dict(study_id=sid,role=r.role,submitted_n=r.sampleCount,figure1_contribution=int(r.sampleCount) if plotted else 0,table1_contribution=int(r.sampleCount) if tab else 0,difference=int(r.sampleCount) if plotted and not tab else 0,include_in_main_recorded=r.include_in_main))
 df=pd.DataFrame(rows);ll=pd.DataFrame(ledger);rr=pd.DataFrame(recon);over=[]
 for a,aa in sets.items():
  for b,bb in sets.items():over.append(dict(study_a=a,study_b=b,matching_identifier_n=len(aa&bb),interpretation='same cohort' if a==b else ('confirmed TCGA barcode overlap' if a.startswith('blca_tcga') and b.startswith('blca_tcga') else 'identifier intersection only; zero cannot prove independent patients')))
 save(df,'cohort_ledgers/cbioportal_study_ledger.csv');save(ll,'cohort_ledgers/cbioportal_patient_sample_ledger.csv');save(pd.DataFrame(over),'cohort_ledgers/cbioportal_overlap_matrix.csv');save(rr,'cohort_ledgers/figure1_table1_row_reconciliation.csv');save(df,'figure_data/figure1_plot_data.csv')
 note('cohort_ledgers/figure1_table1_reconciliation.md',f'''# Submitted cohort count reconciliation
Historical row arithmetic: Figure1 sum {rr.figure1_contribution.sum()}, Table1 sum {rr.table1_contribution.sum()}, difference {rr.difference.sum()}. This discrepancy is analytical-role inclusion: five targeted profiling studies contribute912 to Figure1 but not Table1. It is NOT the number of duplicate patients. Table1 counts are study sample totals across multiple roles/releases, not deduplicated independent tumors. See each row in figure1_table1_row_reconciliation.csv. Submitted include_in_main flag was not consistently used in figure role aggregation. TCGA releases have verified overlapping patient barcodes; the overlap matrix enumerates them. Non-TCGA study-specific identifiers cannot resolve all shared participants; no global unique biological denominator is defensible. Revised forest displays each study separately, explicitly descriptive and potentially overlapping.
Current API acquisition 2026-10-21 may differ from historical portal versions. Mutation frequency uses non-silent/coding-annotation records (excluding Silent, Intron, IGR, UTR/flank) among source sequenced sample-list members; VUS included, no hotspot/pathogenicity restriction. Splice/RNA annotations, any-variant and any-protein-change counts are separately auditable from cached records. No FGFR3 mutation record outside the sequenced list is classified WT. Amplification is discrete CNA=2, gains=1, separate denominator. Wilson95% CI per study. Denominator is samples, not legacy mut_patients divided by all cases. No pooling across studies or independent-replication claim. Eligible sample denominator, study-specific unique IDs and descriptive role totals are distinct.
''')
 print('cbio complete',len(df),flush=True)
def tcga():
 sid='blca_tcga_pan_can_atlas_2018';p=OUT/'public_evidence/cbio'/sid;ll=pd.read_csv(OUT/'cohort_ledgers/cbioportal_patient_sample_ledger.csv');ll=ll[ll.study_id==sid].copy();ids=set(read(p/(sid+'_3way_complete.json')));expr={}
 for path in (OUT/'public_evidence/tcga_expression_all_sample').glob('*.json'):
  data=read(path);expr[path.stem]={r['sampleId']:r['value'] for r in data}
 x=pd.DataFrame(expr);ll=ll.set_index('sample_id');ll['complete_3way']=ll.index.isin(ids);ll['has_FGFR3_expression']=ll.index.isin(x.FGFR3.dropna().index);ll['tumor_barcode']=ll.index.str[13:15].isin(['01','02','06']);ll['primary_eligible']=ll.complete_3way&ll.mutation_profiled&ll.has_FGFR3_expression&ll.tumor_barcode
 ll['exclusion_reason']=np.where(ll.primary_eligible,'included','not complete3way / no mutation profiling / missing expression / non-tumor');save(ll.reset_index(),'cohort_ledgers/tcga_bridge_sample_ledger.csv');rows=[]
 d=x.loc[ll.index[ll.primary_eligible]].copy();d['mutant']=ll.loc[d.index,'FGFR3_coding_mutation'];d['APM']=d[APM].mean(axis=1);d['IFNG']=d[IFNG].mean(axis=1);d['luminal_corrected']=d[LUM].mean(axis=1)
 # Original TACSTD2 was erroneously queried using EPCAM Entrez4072; expose exact duplicate contribution.
 gs=[g for g in LUM if g!='TACSTD2'];d['luminal_submitted_id_error']=(d[gs].sum(axis=1)+d.EPCAM)/len(LUM)
 for sub,ok in [('all_complete',np.ones(len(d),bool)),('luminal_upper_half_corrected',d.luminal_corrected>=d.luminal_corrected.median()),('luminal_upper_half_submitted_id_error',d.luminal_submitted_id_error>=d.luminal_submitted_id_error.median())]:
  dd=d.loc[ok]
  for g in ['FGFR3']+APM+IFNG+['APM','IFNG','luminal_corrected']:
   a=dd.loc[dd.mutant,g];b=dd.loc[~dd.mutant,g];ef=effect(a,b);rows.append(dict(subset=sub,outcome=g,n=len(dd),patients=ll.loc[dd.index,'patient_id'].nunique(),mutant_n=len(a),wildtype_n=len(b),mutant_mean=a.mean(),wt_mean=b.mean(),mutant_median=a.median(),wt_median=b.median(),scale='cBioPortal all-sample reference RNASeqV2 z-score',mutation='non-silent coding; VUS included; amplification separate',**ef))
 save(bh(pd.DataFrame(rows),['subset']),'bulk/tcga_bridge_all_results.csv');save(d.reset_index(names='sample_id'),'figure_data/tcga_bridge_expression.csv')
 assert len(d)==402 and ll.loc[d.index,'patient_id'].nunique()==402
 note('bulk/tcga_bridge_original_vs_revised.md','Verified402 unique TCGA primary tumor samples/patients in the original three-way complete list, mutation-profiled with expression. Exact sample barcodes used; no aliquot expansion. Correct original profile is rna_seq_v2_mrna_median_all_sample_Zscores; alternate diploid-reference profile retained separately and not silently substituted. TACSTD2 Entrez ID corrected4070; submitted script used4072 (EPCAM), causing duplicated EPCAM in luminal score. Both subgroup definitions are supplied. Variant status is not interchangeable with expression, amplification or pathway activity. All results including nonsignificant APM contrasts retained; inspect full estimates/CI/q before interpreting mutation as recapitulating expression-defined effects.')
if __name__=='__main__':run();tcga()
