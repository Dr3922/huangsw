from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd

p=argparse.ArgumentParser()
p.add_argument("--handoff", required=True)
p.add_argument("--render", action="store_true")
a=p.parse_args()
root=Path(a.handoff).resolve()
checks=[]
def check(name, ok, detail):
    checks.append({"check":name,"pass":bool(ok),"detail":str(detail)})
    if not ok: raise AssertionError(f"{name}: {detail}")

inc=pd.read_csv(root/"bulk/bulk_inclusion_counts.csv")
check("GSE13507 primary n", int(inc.query("cohort=='GSE13507' and subset=='primary_tumor_only'").samples.iloc[0])==165, 165)
check("GSE13507 eligible-tumor sensitivity n", int(inc.query("cohort=='GSE13507' and subset=='all_tumors_clustered'").samples.iloc[0])==188, 188)
meta=pd.read_csv(root/"bulk/bulk_meta_results.csv")
m=meta.query("analysis=='revised' and outcome=='Immune' and cutoff==0.75").iloc[0]
check("bulk immune q75 estimate", np.isclose(m.estimate,-0.647662,atol=5e-7),m.estimate)
sn=pd.read_csv(root/"snrna/within_state_bootstrap_summary.csv")
s=sn.query("scale=='logCP10k_count_pseudobulk' and min_cells==20 and feature=='apm_module'").iloc[0]
check("single-nucleus APM cluster bootstrap", np.isclose(s.estimate,-0.646703,atol=5e-7) and int(s.bootstrap_effective)>=4930,s.to_dict())
fold=pd.read_csv(root/"spatial/lopo_fold_manifest.csv")
check("four strict LOPO folds", len(fold)==4 and set(fold.n_training_patients)=={24},fold[["section","n_training_patients"]].to_dict("records"))
check("held-out patient excluded", all(str(r.held_out_patient) not in str(r.training_patients).split(";") for _,r in fold.iterrows()),"all folds")
tc=pd.read_csv(root/"bulk/tcga_bridge_all_results.csv")
t=tc.query("subset=='all_complete' and outcome=='FGFR3'").iloc[0]
check("TCGA bridge denominator", int(t.n)==402 and int(t.mutant_n)==58 and int(t.wildtype_n)==344,t.to_dict())

# Independently recalculate effects from included patient/sample and draw-level data.
from scipy import stats
cut=pd.read_csv(root/'bulk/bulk_cutoff_sensitivity.csv')
selector={'GSE32894':'primary_analysis','GSE13507':'primary_tumor_only','GSE31684':'primary_analysis'}
d=cut[(cut.outcome=='Immune') & (cut.cutoff==.75)]
d=d[[r.subset==selector[r.cohort] for r in d.itertuples()]]
y=d.estimate.to_numpy(); v=d.se.to_numpy()**2; w=1/v
mu=np.sum(w*y)/w.sum(); Q=np.sum(w*(y-mu)**2)
tau=max(0,(Q-len(y)+1)/(w.sum()-(w*w).sum()/w.sum()))
wr=1/(v+tau); est=np.sum(wr*y)/wr.sum(); se=np.sqrt(1/wr.sum())
check('recompute DL meta from cohort effects',len(d)==3 and np.allclose([est,se,tau],[m.estimate,m.se,m.tau2],atol=1e-12),[est,se,tau])
for cohort,subset in selector.items():
    x=pd.read_csv(root/f'bulk/{cohort}_{subset}_restandardized_scores.csv')
    high=x.FGFR3>=x.FGFR3.quantile(.75); a1=x.loc[high,'Immune']; a0=x.loc[~high,'Immune']
    sd=np.sqrt(((len(a1)-1)*a1.var(ddof=1)+(len(a0)-1)*a0.var(ddof=1))/(len(x)-2))
    cd=(a1.mean()-a0.mean())/sd
    expected=d[d.cohort==cohort].iloc[0]
    check(f'recompute {cohort} Cohen d',np.isclose(cd,expected.estimate,atol=1e-10),cd)
pat=pd.read_csv(root/'snrna/pseudobulk_expression_patient.csv')
pc=pd.read_csv(root/'snrna/patient_pseudobulk_correlations.csv')
pr=pc.query("scale=='logCP10k_count_pseudobulk' and feature=='apm_module'").iloc[0]
rho=stats.spearmanr(pat.FGFR3,pat.apm_module).statistic
check('recompute independent donor APM correlation',pat.Patient.nunique()==25 and np.isclose(rho,pr.estimate,atol=1e-12),rho)
draw=pd.read_csv(root/'snrna/within_state_patient_cluster_bootstrap.csv.gz')
v=draw.query("scale=='logCP10k_count_pseudobulk' and min_cells==20 and feature=='apm_module'").rho.dropna().to_numpy()
ci=np.quantile(v,[.025,.975]); centered=v-s.estimate
pv=(np.sum(np.abs(centered)>=abs(s.estimate))+1)/(len(v)+1)
check('recompute donor-bootstrap CI and P',len(v)==int(s.bootstrap_effective) and np.allclose(ci,[s.ci_low,s.ci_high],atol=1e-12) and np.isclose(pv,s.p),[len(v),*ci,pv])
family=sn.query("scale=='logCP10k_count_pseudobulk' and min_cells==20 and family=='modules'").dropna(subset=['p'])
if family.empty:
    family=sn.query("scale=='logCP10k_count_pseudobulk' and min_cells==20").query("feature in ['apm_module','ifng_execution_module','antigen_visibility_module']").dropna(subset=['p'])
order=np.argsort(family.p.to_numpy()); ordered=family.p.to_numpy()[order]
qv=np.minimum(1,np.minimum.accumulate((ordered*len(ordered)/np.arange(1,len(ordered)+1))[::-1])[::-1])
check('recompute module BH correction',len(family)==3 and np.allclose(qv,family.q.to_numpy()[order]),qv.tolist())
spots=pd.read_csv(root/'spatial/lopo_spot_scores.csv')
effects=pd.read_csv(root/'spatial/spatial_block_resampling.csv')
for sec,frame in spots.groupby('section'):
    rho=stats.spearmanr(frame.projection,frame.local_dc_original_k6).statistic
    expected=effects.query("section==@sec and k==6 and block_width_nn==3 and outcome=='dc_original'").iloc[0]
    check(f'recompute spatial DC rho {sec}',np.isclose(rho,expected.estimate,atol=1e-12),rho)

if a.render:
    env=os.environ.copy(); env["FGFR3_OUT"]=str(root)
    cp=subprocess.run([sys.executable,str(Path(__file__).with_name("build_figures.py"))],env=env,text=True,capture_output=True)
    check("render all figures",cp.returncode==0,cp.stdout+cp.stderr)
    expected=[root/f"figures/Figure{i}_revised.pdf" for i in range(1,6)] + list((root/"figures/supplementary").glob("*.pdf"))
    check("eleven PDFs exist",len(expected)==11 and all(x.exists() and x.stat().st_size>10000 for x in expected),len(expected))

versions={k:getattr(__import__(k),"__version__","unknown") for k in ["numpy","pandas","scipy","matplotlib","h5py"]}
out={"python":sys.version,"versions":versions,"checks":checks,"status":"PASS"}
(Path(__file__).parents[1]/"SMOKE_RESULT.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
print(json.dumps(out,indent=2))
