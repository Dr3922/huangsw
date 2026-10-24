from reanalysis import *
from scipy.spatial import cKDTree
from statsmodels.stats.outliers_influence import variance_inflation_factor
class old:
 @staticmethod
 def zscore_frame(df):return ((df-df.mean())/df.std().replace(0,1)).fillna(0)
 @staticmethod
 def keep_marker_gene(g):return not any(re.match(p,g) for p in [r'^RP[SL]',r'^MT-',r'^MIR',r'^SNOR',r'^LINC',r'^AL\d',r'^AC\d',r'^AP\d',r'^CTD-',r'^RNA',r'^OR\d',r'^NEAT1$'])
DC=['FCER1A','HLA-DRA','HLA-DPA1','HLA-DPB1','CLEC10A','ITGAX'];CD8=['CD8A','CD8B','NKG7','PRF1','GZMB'];STROMA=['FAP','ACTA2','COL1A1','COL1A2','PDGFRB','FN1','TGFB1','CXCL12']
CONS=['DHRS2','ERBB2','GRB7','GRHL2','KRT13','NCALD','NDRG1','PTPRQ','RNF19A','SCHLAP1','SLC2A1','SNX31','TCIRG1','UNC93B1','PTPN13']
SETS={'apm':APM,'ifng':IFNG,'dc_original':DC,'cd8':CD8,'stromal':STROMA,'luminal':LUM,'epithelial':['EPCAM','TACSTD2','KRT8','KRT18','KRT19'],'dc_no_hla':['FCER1A','CLEC10A','ITGAX'],'cDC1':['CLEC9A','XCR1','BATF3','CADM1'],'cDC2':['CD1C','FCER1A','CLEC10A'],'myeloid':['LST1','TYROBP','FCER1G','AIF1','CTSS'],'mhcii':['HLA-DRA','HLA-DPA1','HLA-DPB1','CD74'],'consensus':CONS}

def module(x,genes):return x[[g for g in genes if g in x]].mean(axis=1)
def derive(e,base,patient,global_scale=False,shuffle=False):
 mask=np.asarray(e.obs.Patient.astype(str)!=patient);tr=base.loc[mask].copy();sc=old.zscore_frame(base if global_scale else tr)
 score=sc.FGFR3-sc.apm_module-.5*sc.ifng_execution_module
 if global_scale:score=score.loc[mask]
 if shuffle:score=pd.Series(np.random.default_rng(SEED).permutation(score),index=score.index)
 cut=score.quantile(.8);hi=np.zeros(len(base),bool);hi[np.flatnonzero(mask)]=np.asarray(score>=cut);lo=mask&~hi
 mh=np.asarray(e.X[hi].mean(axis=0)).ravel();ml=np.asarray(e.X[lo].mean(axis=0)).ravel();dh=np.asarray((e.X[hi]>0).mean(axis=0)).ravel()
 t=pd.DataFrame({'gene_symbol':e.var_names.astype(str),'mean_high':mh,'mean_low':ml,'logfc_proxy':np.log1p(mh)-np.log1p(ml),'det_high':dh})
 blocked=set(APM+IFNG+['FGFR3','VEGFA','FOS','FOSB']);keep=t.gene_symbol.map(old.keep_marker_gene)&~t.gene_symbol.isin(blocked)&(t.det_high>=.1)&(t.logfc_proxy>0)
 # All genes shared in both deposited whole-transcriptome feature spaces; no spatial-expression screening.
 t=t[keep].sort_values(['logfc_proxy','det_high'],ascending=False).reset_index(drop=True);t['rank']=np.arange(1,len(t)+1);t['selected']=t['rank']<=15
 params=pd.DataFrame({'feature':tr.columns,'mean':tr.mean().to_numpy(),'sd':tr.std().to_numpy()})
 return t,cut,params,score,mask

def blockboot(x,y,xy,width,rng,B=1200):
 ids=np.floor((xy-xy.min(axis=0))/width).astype(int);_,labels=np.unique(ids,axis=0,return_inverse=True);groups=[np.flatnonzero(labels==j) for j in range(labels.max()+1)]
 draws=[];diff=[]
 for _ in range(B):
  ii=np.concatenate([groups[j] for j in rng.integers(0,len(groups),len(groups))]);xx=x[ii];yy=y[ii]
  draws.append(stats.spearmanr(xx,yy).statistic);hi=xx>=np.quantile(xx,.75);lo=xx<=np.quantile(xx,.25);diff.append(np.mean(yy[hi])-np.mean(yy[lo]))
 return np.array(draws),np.array(diff),labels,len(groups)

def run():
 rng=np.random.default_rng(SEED);a=ad.read_h5ad(SC);e=a[np.asarray(a.obs.celltype=='Epithelial')].copy();del a
 base=e[:,['FGFR3']+APM+IFNG].to_df().reset_index(drop=True);base['apm_module']=base[APM].mean(axis=1);base['ifng_execution_module']=base[IFNG].mean(axis=1);base=base[['FGFR3','apm_module','ifng_execution_module']]
 sp=ad.read_h5ad(SP)
 assert set(e.var_names)==set(sp.var_names),'Frozen feature universe differs; specify technical intersection before deriving markers'
 assert np.all(sp.X.data>=0) and np.all(sp.X.data==np.floor(sp.X.data))
 mapping=[];manifest=[];preds=[];effects=[];models=[];vifs=[];overlap=[];increment=[];controls=[]
 for sid,pat in sp.obs[['sampleID','Patient']].drop_duplicates().itertuples(index=False):
  sid=str(sid);pat=str(pat);gsm={'A1':'GSM5224027','B1':'GSM5224028','C1':'GSM5224029','D1':'GSM5224030'}[sid];soft=(OUT/'public_evidence'/f'{gsm}.soft').read_text();num=pat.replace('Bladder','')
  assert f'MIBC_B{num}' in soft and f'{sid}.jpg' in soft
  assert pat in set(e.obs.Patient.astype(str))
  mapping.append(dict(section=sid,specimen=gsm,patient=pat,site='bladder',batch=sid,mapping_evidence=f'GEO {gsm}: title MIBC_B{num}; supplement {sid}.jpg; deposited h5ad obs same identity; original gzip checksum verified',biological_independence='same atlas patient present in snRNA; excluded from this training fold',confidence='source-verified deposited pseudonymous subject mapping',raw_metadata_filename=f'public_evidence/{gsm}.soft'))
  markers,cut,params,score,mask=derive(e,base,pat);orig,_,_,_,_=derive(e,base,pat,True);shuf,_,_,_,_=derive(e,base,pat,shuffle=True)
  # Perturb held-out expression only. Strict training score/markers must remain exactly unchanged.
  pert=base.copy();pert.loc[~mask]=pert.loc[~mask]*100000+12345
  check,cc,_,_,_=derive(e,pert,pat);assert markers.equals(check) and cc==cut
  controls.append(dict(section=sid,heldout_perturbation_marker_identical=True,heldout_perturbation_cut_identical=True,shuffle_selected_overlap=len(set(shuf.head(15).gene_symbol)&set(markers.head(15).gene_symbol)),shuffle_purpose='pipeline negative control, not a patient-level null test'))
  save(markers,f'spatial/lopo_marker_lists/{sid}_all_ranked.csv');save(orig,f'spatial/lopo_marker_lists/{sid}_original_global_z.csv');save(shuf.head(15),f'spatial/lopo_marker_lists/{sid}_shuffled_control.csv');save(params,f'spatial/lopo_marker_lists/{sid}_training_standardization.csv');save(pd.DataFrame({'training_cell_index':np.flatnonzero(mask),'state_score':score}),f'spatial/lopo_marker_lists/{sid}_training_score_distribution.csv')
  selected=markers.head(15).gene_symbol.tolist();oldgenes=orig.head(15).gene_symbol.tolist();assert len(selected)==15 and pat not in set(e.obs.Patient.astype(str)[mask])
  s=sp[np.asarray(sp.obs.sampleID.astype(str)==sid)].copy();gg=sorted(set(sum(SETS.values(),[])+selected+oldgenes+['FGFR3']));gg=[g for g in gg if g in s.var_names];raw=s[:,gg].to_df().reset_index(drop=True);x=old.zscore_frame(raw)
  d=s.obs.reset_index(names='spot_id');d['section']=sid;d['x']=s.obsm['spatial'][:,0];d['y']=s.obsm['spatial'][:,1]
  for name,gs in SETS.items():d[name]=module(x,gs)
  d['signature_only']=module(x,selected);d['original_projection']=module(x,oldgenes)-d.apm-.5*d.ifng;d['projection']=d.signature_only-d.apm-.5*d.ifng;d['consensus_projection']=d.consensus-d.apm-.5*d.ifng;d['FGFR3']=x.FGFR3
  # Normalized spatial expression sensitivity, same frozen markers; original count scale retained above.
  xx=raw.div(np.asarray(s.X.sum(axis=1)).ravel(),axis=0)*10000;xx=old.zscore_frame(np.log1p(xx));d['normalized_projection']=module(xx,selected)-module(xx,APM)-.5*module(xx,IFNG)
  xy=d[['x','y']].to_numpy();dist,nn=cKDTree(xy).query(xy,k=9);assert all(i not in nn[i,1:] for i in range(len(d)))
  spacing=float(np.median(dist[:,1]));d['log_depth']=np.log1p(d.total_counts);d['log_genes']=np.log1p(d.n_genes_by_counts)
  manifest.append(dict(section=sid,held_out_patient=pat,training_patients=';'.join(sorted(set(e.obs.Patient.astype(str)[mask]))),n_training_patients=e.obs.Patient[mask].nunique(),n_training_cells=int(mask.sum()),heldout_sn_cells=int((~mask).sum()),n_spots=len(s),cutoff_q80=cut,genes=';'.join(selected),weights='equal 1/15',score_formula='mean(section gene-z selected15) - mean(APM gene-z) - 0.5*mean(IFNG gene-z)',train_normalization='FGFR3, raw-count-mean APM and IFNG modules: z parameters computed using training cells only',test_normalization='gene-wise within-section z of deposited raw counts; no outcome-based tuning',seed=SEED,marker_hyperparameters='q80; detection>=0.1; positive log1p mean difference; fixed original exclusions; top15',spatial_coordinate_unit='full-resolution image pixel',median_nn_distance_pixels=spacing,qc='all deposited in-tissue QC-passed spots; no new exclusion',original_marker_overlap=len(set(selected)&set(oldgenes))))
  for key,gs in SETS.items():
   for other,hs in SETS.items():overlap.append(dict(section=sid,signature=key,other=other,shared=';'.join(sorted(set(gs)&set(hs))),n_overlap=len(set(gs)&set(hs)),present=';'.join(g for g in gs if g in raw),missing=';'.join(g for g in gs if g not in raw)))
  for k in [4,6,8]:
   targets=['apm','cd8','dc_original','dc_no_hla','cDC1','cDC2','myeloid','mhcii','stromal']
   for target in targets:d[f'local_{target}_k{k}']=d[target].to_numpy()[nn[:,1:k+1]].mean(axis=1)
   for target in targets:
    y=d[f'local_{target}_k{k}'].to_numpy();xx=d.projection.to_numpy();rr,pp=stats.spearmanr(xx,y)
    for widthn in [3,6]:
     draws,diffs,labels,nb=blockboot(xx,y,xy,spacing*widthn,rng);bs=bootstrap_summary(draws,rr);hi=xx>=np.quantile(xx,.75);lo=xx<=np.quantile(xx,.25);delta=np.mean(y[hi])-np.mean(y[lo]);di=np.nanquantile(diffs,[.025,.975]);d[f'block_{widthn}']=labels
     effects.append(dict(section=sid,patient=pat,k=k,outcome=target,n_spots=len(d),biological_n=1,block_width_nn=widthn,block_width_pixels=spacing*widthn,n_blocks=nb,high_n=hi.sum(),low_n=lo.sum(),quartile_mean_difference=delta,difference_ci_low=di[0],difference_ci_high=di[1],descriptive_spot_p=pp,**bs))
   if k==6:
    for target in targets:
     y=d[f'local_{target}_k{k}'];covs=['luminal','epithelial','stromal','log_depth','log_genes'];cov=z(d[covs]);cov['FGFR3']=z(d.projection)
     for subset,ok in [('all',np.ones(len(d),bool)),('epithelial_top50',d.epithelial>=d.epithelial.quantile(.5)),('epithelial_top25',d.epithelial>=d.epithelial.quantile(.75))]:
      # Blocks are descriptive spatial clusters; four sections remain biological n=4.
      r=fit(z(y.loc[ok]),cov.loc[ok],d.loc[ok,'block_3']);models.append(dict(section=sid,outcome=target,subset=subset,model='projection + luminal + epithelial + stromal + log depth + log genes; spatial block cluster',**r))
     for variable in ['projection','signature_only','original_projection','consensus_projection','FGFR3','normalized_projection','luminal','apm']:
      r,p=stats.spearmanr(d[variable],y);increment.append(dict(section=sid,outcome=target,predictor=variable,model='descriptive Spearman',estimate=r,p=p,n=len(d)))
     for predictor in ['projection','signature_only','FGFR3']:
      xx0=sm.add_constant(z(d[['luminal','epithelial','stromal','log_depth','log_genes']]));xx1=xx0.join(z(d[[predictor]]));m0=sm.OLS(z(y),xx0).fit();m1=sm.OLS(z(y),xx1).fit();res=fit(z(y),xx1.drop(columns='const'),d.block_3,coef=predictor)
      increment.append(dict(section=sid,outcome=target,predictor=predictor,model='increment beyond luminal/composition/technical',partial_r2=(m1.rsquared-m0.rsquared)/(1-m0.rsquared),**res))
    vv=sm.add_constant(cov)
    for j,col in enumerate(vv):vifs.append(dict(section=sid,variable=col,vif=variance_inflation_factor(vv.to_numpy(),j),condition_number=np.linalg.cond(vv)))
  preds.append(d);print('spatial complete',sid,flush=True)
 save(pd.DataFrame(mapping),'spatial/spatial_section_patient_map.csv');save(pd.DataFrame(mapping),'spatial/lopo_patient_section_mapping.csv');save(pd.DataFrame(manifest),'spatial/lopo_fold_manifest.csv');save(pd.concat(preds),'spatial/lopo_spot_scores.csv');save(pd.DataFrame(controls),'spatial/lopo_negative_controls.csv');save(pd.DataFrame(effects),'spatial/spatial_section_effects.csv');save(pd.DataFrame(models),'spatial/spatial_composition_adjusted.csv');save(pd.DataFrame(vifs),'spatial/spatial_vif.csv');save(pd.DataFrame(overlap),'spatial/spatial_signature_overlap_coverage.csv');save(pd.DataFrame(increment),'spatial/atlas_incremental_spatial_comparison.csv')
 save(pd.DataFrame(increment).query("model == 'descriptive Spearman'"),'spatial/lopo_original_vs_corrected.csv')
 note('spatial/lopo_leakage_checks.md','''# LOPO audit
Original algorithm excludes held-out patients from differential-marker means, but computes FGFR3/APM/IFNG cell-score standardization on all epithelial cells before exclusion: normalization leakage. Revised derivation fits all three standardization parameters on training cells only, recomputes q80 cutoff and marker rankings with fixed original filters/top15, then applies frozen fold markers. Test-section gene z standardization is unsupervised and explicitly retained from submission; expression-normalized sensitivity is supplied. Whole-transcriptome gene identifiers shared between deposited objects are a technical compatibility check, not expression-driven feature selection. Consensus15 was exported AFTER original fold-specific projection; the actual original script uses fold-specific genes. Therefore evidence does not support claiming consensus leakage in the submitted implementation; consensus scores are descriptive only. Held-out expression perturbation leaves all strict marker rankings and cutoffs identical. Shuffled-training-label markers retained as pipeline negative controls, not evidence of absence of biological association.
Patient mapping is triangulated using live GSM titles MIBC_B[number], per-GSM A1/B1/C1/D1 image supplement names and checksum-matched original public h5ad patient/sampleID metadata. These are public pseudonymous atlas IDs, not clinical identity verification. All four spatial patients occur in snRNA and are excluded from their own fold; atlas is not independent external validation.
Spatial effects are reported per section. The biological n is four sections/patients, not4086. Nonself nearest neighbors k4/6/8 use full-resolution image pixels within section. All deposited spots are in-tissue. No disconnected-component annotation exists; long-distance nearest-neighbor links and tissue gaps are a limitation. Grid-block cluster bootstrap1200draws, two widths3/6 median nearest-neighbor distances, preserves within-block paired spatial structure for uncertainty conditional on this section; it is NOT a four-patient null/randomization test. P from centered bootstrap is approximate and conditional; spot P strictly descriptive. Different fixed grid widths assess sensitivity; no claim of eliminating all spatial dependence. Composition markers are proxies, not inferred cell fractions; model block-cluster intervals are conditional descriptive sensitivity.
Projection subtracts APM and IFNG directly, mathematically inducing negative associations with these components and potentially correlated immune scores. Signature-only, FGFR3-only, APM-only, generic luminal, consensus and normalized-scale comparisons are provided. No score establishes DC abundance, malignant identity, antigen presentation function or FGFR3 causality.
''')
if __name__=='__main__':run()
