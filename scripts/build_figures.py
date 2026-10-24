from pathlib import Path
import json, os

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get('FGFR3_OUT', PACKAGE_ROOT/'data/derived')).expanduser().resolve()
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)
SUPP = FIG / "supplementary"
SUPP.mkdir(exist_ok=True)
BLUE, ORANGE, GREY, RED = "#2F6DAE", "#D97732", "#667085", "#B42318"
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 7, "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .6})


def panel(ax, letter):
    ax.text(-.12, 1.13, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def save(fig, name, directory=FIG):
    fig.savefig(directory / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(directory / f"{name}.svg", bbox_inches="tight")
    fig.savefig(directory / f"{name}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(directory / f"{name}_preview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def forest(ax, d, label, color, title, xlabel):
    d = d.reset_index(drop=True)
    y = np.arange(len(d))[::-1]
    ax.errorbar(d.estimate, y, xerr=[d.estimate-d.ci_low, d.ci_high-d.estimate], fmt="o", color=color, ecolor=color, ms=3, lw=.8, capsize=1.5)
    ax.axvline(0, color="#98A2B3", lw=.7, ls="--")
    ax.set_yticks(y, d[label])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.spines[["top", "right"]].set_visible(False)


def figure1():
    d = pd.read_csv(OUT / "figure_data/figure1_plot_data.csv")
    d = d[~d.role.str.contains("Preclinical|Variant", regex=True)].copy()
    d["label"] = d.study_id
    d = d.sort_values(["role", "frequency"])
    d = d.rename(columns={"frequency": "estimate"})
    d = d.reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 7.7), gridspec_kw={"width_ratios": [2.0, 1]})
    y = np.arange(len(d))[::-1]
    axes[0].set_xlim(0,1); axes[0].set_ylim(-.6,len(d)+.6); axes[0].axis("off")
    for xx,head in [(0,"Study / source ID"),(.57,"Setting / role"),(.87,"n"),(.98,"Mut")]:
        axes[0].text(xx,len(d),head,fontsize=6,fontweight="bold",ha="right" if xx>.8 else "left")
    for yy,(_,r) in zip(y,d.iterrows()):
        axes[0].text(0,yy,r.study_id,fontsize=5.3,va="center")
        axes[0].text(.57,yy,r.role.replace("Discovery / ","D / ").replace("Targeted / enriched","Targeted").replace("Clinical context","Clinical").replace("UTUC comparator","UTUC"),fontsize=5.2,va="center")
        axes[0].text(.87,yy,str(int(r.mutation_eligible_n)),fontsize=5.3,va="center",ha="right")
        axes[0].text(.98,yy,str(int(r.coding_mutant_n)),fontsize=5.3,va="center",ha="right")
    axes[1].errorbar(d.estimate,y,xerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt="o",ms=3,lw=.8,capsize=1.5,color=BLUE)
    axes[1].set_ylim(-.6,len(d)+.6); axes[1].set_yticks([])
    axes[1].set(xlim=(0,1),xlabel="Mutation proportion (Wilson 95% CI)",title="FGFR3 coding mutations")
    axes[1].spines[["top","right","left"]].set_visible(False)
    for yy,(_,r) in zip(y,d.iterrows()):
        axes[1].text(.99,yy,f"{r.estimate:.2f}" if pd.notna(r.estimate) else "not estimable",ha="right",va="center",fontsize=5.2)
    fig.text(.05,.035,"Each row is a study release; TCGA releases overlap and other cross-study overlap is unresolved.",fontsize=5.5,color=GREY)
    fig.text(.05,.018,"n = source-confirmed mutation-profiled samples; Mut = coding-mutant samples (VUS included).",fontsize=5.5,color=GREY)
    fig.subplots_adjust(wspace=.06,bottom=.09)
    save(fig, "Figure1_revised")
    c = d[d.cna_eligible_n > 0].copy()
    c["amp"] = c.amplified_n / c.cna_eligible_n
    c["gain_prop"] = c.gain_n / c.cna_eligible_n
    c = c.sort_values("amp")
    y = np.arange(len(c))[::-1]
    fig,ax=plt.subplots(figsize=(7.15,5.4)); axes=[None,ax]
    axes[1].scatter(c.amp, y, s=11, color=ORANGE, label="Amplification (+2)")
    axes[1].scatter(c.gain_prop, y, s=10, facecolors="none", edgecolors=GREY, label="Gain (+1)")
    axes[1].set_yticks(y, c.label)
    axes[1].set_xlabel("Proportion among returned FGFR3 CNA values")
    axes[1].set_title("FGFR3 copy-number endpoints (descriptive)")
    axes[1].legend(frameon=False, loc="lower right")
    axes[1].spines[["top", "right"]].set_visible(False)
    save(fig, "Supplementary_Figure_S6_CNA", SUPP)


def figure2():
    meta = pd.read_csv(OUT / "bulk/bulk_meta_results.csv")
    q75 = meta[(meta.cutoff == .75) & meta.analysis.isin(["submitted", "revised", "all_tumors_sensitivity"])].copy()
    q75["label"] = q75.analysis.str.replace("_", " ") + " | " + q75.outcome
    imm = q75[q75.outcome == "Immune"]
    main = q75[(q75.analysis == "revised") & q75.outcome.isin(["Luminal-papillary", "Immune", "Antigen presentation", "Cytolytic"])]
    adj = pd.read_csv(OUT / "bulk/bulk_proxy_adjusted.csv")
    wanted = {"GSE32894": "primary_analysis", "GSE13507": "primary_tumor_only", "GSE31684": "primary_analysis"}
    adj = adj[adj.apply(lambda r: wanted.get(r.cohort) == r.subset, axis=1) & adj.outcome.isin(["Immune", "Antigen presentation"])].copy()
    adj["label"] = adj.cohort + "\n" + adj.outcome.replace({"Antigen presentation":"APM"})
    cut = meta[(meta.analysis == "revised") & (meta.outcome == "Immune")].sort_values("cutoff")
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.4), gridspec_kw={"height_ratios": [1, .78]})
    axes = axes.ravel()
    forest(axes[0], imm, "label", ORANGE, "Immune q75:\noriginal and revised", "Cohen d")
    forest(axes[1], main, "outcome", BLUE, "Revised three-cohort\nsynthesis", "Cohen d")
    forest(axes[2], adj, "label", GREY, "Proxy/context-adjusted\nFGFR3 term", "Regression coefficient")
    axes[3].plot(cut.cutoff, cut.estimate, marker="o", ms=3, color=BLUE)
    axes[3].fill_between(cut.cutoff, cut.ci_low, cut.ci_high, color=BLUE, alpha=.18, linewidth=0)
    axes[3].axhline(0, color="#98A2B3", lw=.7, ls="--")
    axes[3].set(xlabel="FGFR3-high quantile cutoff", ylabel="Meta-analytic Cohen d", title="Revised Immune cutoff robustness")
    axes[3].spines[["top", "right"]].set_visible(False)
    axes[0].set_xlim(-.9,.05); axes[1].set_xlim(-.9,.8); axes[2].set_xlim(-.48,.08)
    for i, ax in enumerate(axes): panel(ax, chr(97+i))
    fig.subplots_adjust(wspace=.72, hspace=.72)
    save(fig, "Figure2_revised")


def figure3():
    pat = pd.read_csv(OUT / "snrna/pseudobulk_expression_patient.csv")
    boot = pd.read_csv(OUT / "snrna/within_state_bootstrap_summary.csv")
    boot = boot[(boot.scale == "logCP10k_count_pseudobulk") & (boot.min_cells == 20)].copy()
    inc = pd.read_csv(OUT / "snrna/atlas_incremental_nested_models.csv")
    inc = inc[(inc.scale == "logCP10k_count_pseudobulk") & (inc.min_cells == 20) & inc.feature.isin(["apm_module", "ifng_execution_module", "antigen_visibility_module"])].copy()
    display = {"apm_module":"APM", "ifng_execution_module":"IFN-gamma", "antigen_visibility_module":"Combined"}
    boot["display"] = boot.feature.replace(display)
    inc["display"] = inc.feature.replace(display)
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 4.0), gridspec_kw={"width_ratios": [.95, 1.15, 1]})
    ax = axes[0]
    ax.scatter(pat.FGFR3, pat.apm_module, s=13, color=BLUE, alpha=.8)
    rho, p = stats.spearmanr(pat.FGFR3, pat.apm_module)
    ax.set(xlabel="FGFR3 logCP10k pseudobulk", ylabel="APM module", title=f"Patient pseudobulk (n=25)\nSpearman ρ={rho:.2f}, p={p:.3g}")
    ax.spines[["top", "right"]].set_visible(False)
    forest(axes[1], boot, "display", ORANGE, "Within-state patient-cluster\nbootstrap", "Pooled Spearman ρ")
    for j,r in boot.reset_index(drop=True).iterrows():
        if pd.isna(r.estimate): axes[1].text(0,len(boot)-1-j,"not estimable",fontsize=5,va="center")
    forest(axes[2], inc, "display", GREY, "FGFR3 beyond state\nand luminal score", "FGFR3 coefficient")
    for i, ax in enumerate(axes): panel(ax, chr(97+i))
    fig.subplots_adjust(wspace=.65)
    save(fig, "Figure3_revised")


def spatial_axis(ax, image, d, value=None, title="", cmap="coolwarm"):
    ax.imshow(image)
    if value is not None:
        sc = ax.scatter(d.x_img, d.y_img, c=d[value], s=3.2, cmap=cmap, alpha=.82, linewidths=0)
        plt.colorbar(sc, ax=ax, fraction=.035, pad=.01)
    ax.set_title(title)
    ax.set_axis_off()


def figure4():
    asset_dir = OUT / 'figure_data/spatial_assets'
    if asset_dir.exists():
        a = None
    else:
        spatial_h5ad = os.environ.get("FGFR3_SPATIAL_H5AD")
        if not spatial_h5ad:
            raise FileNotFoundError("Figure 4 requires figure_data/spatial_assets or FGFR3_SPATIAL_H5AD")
        a = ad.read_h5ad(spatial_h5ad, backed="r")
    scores = pd.read_csv(OUT / "spatial/lopo_spot_scores.csv")
    fig, axes = plt.subplots(4, 5, figsize=(7.15, 6.6))
    for r, sec in enumerate(["A1", "B1", "C1", "D1"]):
        if a is None:
            info = {'scalefactors':json.loads((asset_dir / f'{sec}_scalefactors.json').read_text())}
            image = np.load(asset_dir / f'{sec}_hires.npy')
        else:
            info = a.uns["spatial"][sec]; image = info["images"]["hires"]
        scale = float(info["scalefactors"]["tissue_hires_scalef"])
        d = scores[scores.section == sec].copy(); d["x_img"] = d.x * scale; d["y_img"] = d.y * scale
        spatial_axis(axes[r,0], image, d, title=f"{sec} H&E")
        # Visium nominal spot diameter is 55 µm; convert through the object-provided
        # full-resolution spot diameter and hires scale factor. Bar is 500 µm.
        bar_px = 500/55 * float(info["scalefactors"]["spot_diameter_fullres"]) * scale
        x1, x2, yy = image.shape[1] - 45 - bar_px, image.shape[1] - 45, image.shape[0] - 60
        axes[r,0].plot([x1,x2],[yy,yy],color="black",lw=2)
        axes[r,0].text((x1+x2)/2,yy-18,"500 µm",ha="center",va="bottom",fontsize=5,color="black")
        spatial_axis(axes[r,1], image, d, "projection", "strict-LOPO projection")
        spatial_axis(axes[r,2], image, d, "local_dc_original_k6", "local DC-module signal")
        spatial_axis(axes[r,3], image, d, "local_cd8_k6", "local CD8 signal")
        spatial_axis(axes[r,4], image, d, "local_apm_k6", "local APM signal")
        axes[r,0].text(.02, .96, f"patient {d.Patient.iloc[0]}", transform=axes[r,0].transAxes, va="top", color="white", fontsize=5.5, bbox={"facecolor":"black","alpha":.45,"pad":1,"edgecolor":"none"})
    for i, ax in enumerate(axes[0]):
        ax.text(-.08, 1.10, chr(97+i), transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
    fig.subplots_adjust(wspace=.06, hspace=.10)
    if a is not None: a.file.close()
    save(fig, "Figure4_revised")


def figure5():
    d = pd.read_csv(OUT / "spatial/spatial_section_effects.csv")
    main = d[(d.k == 6) & (d.block_width_nn == 3) & d.outcome.isin(["dc_original", "cd8", "apm"])].copy()
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 3.1))
    for ax, out, color, letter in zip(axes, ["dc_original", "cd8", "apm"], [ORANGE, GREY, BLUE], ["a","b","c"]):
        x = main[main.outcome == out].copy(); x["label"] = x.section + " (" + x.n_spots.astype(str) + "/" + x.n_blocks.astype(str) + ")"
        forest(ax, x, "label", color, {"dc_original":"Local DC-module signal","cd8":"Local CD8 signal","apm":"Local APM signal"}[out], "Section-specific block-bootstrap ρ")
        panel(ax, letter)
    fig.text(.5, .015, "Labels show section (spots/spatial blocks); each row is one biological section.", ha="center", color=GREY, fontsize=5.5)
    fig.subplots_adjust(wspace=.9, bottom=.19)
    save(fig, "Figure5_revised")


def supplementary_figures():
    meta = pd.read_csv(OUT / "bulk/bulk_meta_results.csv")
    x = meta[(meta.analysis == "revised") & meta.outcome.isin(["Immune", "Antigen presentation", "Cytolytic"])].copy()
    fig, ax = plt.subplots(figsize=(7.15, 3.4))
    for out, color in zip(["Immune", "Antigen presentation", "Cytolytic"], [BLUE, ORANGE, GREY]):
        z = x[x.outcome == out].sort_values("cutoff")
        ax.plot(z.cutoff, z.estimate, marker="o", ms=3, label=out.replace("Antigen presentation", "APM"), color=color)
        ax.fill_between(z.cutoff, z.ci_low, z.ci_high, color=color, alpha=.12, linewidth=0)
    ax.axhline(0, color="#98A2B3", lw=.7, ls="--"); ax.set(xlabel="FGFR3-high quantile cutoff", ylabel="Meta-analytic Cohen d", title="Bulk cutoff sensitivity")
    ax.legend(frameon=False, ncol=3); ax.spines[["top", "right"]].set_visible(False); panel(ax, "a")
    save(fig, "Supplementary_Figure_S1_bulk_cutoffs", SUPP)

    d = pd.read_csv(OUT / "spatial/spatial_section_effects.csv")
    d = d[d.outcome == "dc_original"].copy(); d["setting"] = "k=" + d.k.astype(str) + ", width=" + d.block_width_nn.astype(str)
    fig, axes = plt.subplots(1, 4, figsize=(7.15, 2.6), sharex=True)
    for ax, (sec, z) in zip(axes, d.groupby("section", sort=True)):
        forest(ax, z.sort_values(["k","block_width_nn"]), "setting", ORANGE, sec, "Block-bootstrap rho")
    for i, ax in enumerate(axes): panel(ax, chr(97+i))
    fig.subplots_adjust(wspace=.72); save(fig, "Supplementary_Figure_S2_spatial_parameter_sensitivity", SUPP)

    d = pd.read_csv(OUT / "spatial/spatial_section_effects.csv")
    d = d[(d.k == 6) & (d.block_width_nn == 3) & d.outcome.isin(["dc_original","dc_no_hla","cDC1","cDC2","myeloid","mhcii"])].copy()
    labels={"dc_original":"Original DC", "dc_no_hla":"DC without HLA", "cDC1":"cDC1", "cDC2":"cDC2", "myeloid":"Myeloid", "mhcii":"MHC-II"}; d["display"] = d.outcome.replace(labels)
    fig, axes = plt.subplots(1, 4, figsize=(7.15, 3.2), sharex=True)
    for ax, (sec, z) in zip(axes, d.groupby("section", sort=True)):
        forest(ax, z, "display", BLUE, sec, "Block-bootstrap rho")
    for i, ax in enumerate(axes): panel(ax, chr(97+i))
    fig.subplots_adjust(wspace=.72); save(fig, "Supplementary_Figure_S3_alternative_DC_signatures", SUPP)

    d = pd.read_csv(OUT / "spatial/spatial_composition_adjusted.csv")
    d = d[(d.outcome.isin(["dc_original","apm","cd8"])) & (d.subset == "all")].copy(); d["display"] = d.section + " | " + d.outcome.replace({"dc_original":"DC"})
    fig, ax = plt.subplots(figsize=(7.15, 4.0)); forest(ax, d, "display", GREY, "Composition- and technical-variable-adjusted spatial models", "FGFR3 projection coefficient")
    panel(ax, "a"); save(fig, "Supplementary_Figure_S4_spatial_adjusted_models", SUPP)

    d = pd.read_csv(OUT / "bulk/bulk_meta_leave_one_cohort_out.csv")
    d = d[(d.analysis == "revised") & d.outcome.isin(["Immune","Antigen presentation"])].copy()
    d["display"] = d.outcome.replace({"Antigen presentation":"APM"}) + " | omit " + d.omitted
    fig, ax = plt.subplots(figsize=(7.15, 3.4)); forest(ax, d, "display", ORANGE, "Leave-one-cohort-out synthesis", "Cohen d")
    panel(ax, "a"); save(fig, "Supplementary_Figure_S5_bulk_influence", SUPP)


if __name__ == "__main__":
    figure1(); figure2(); figure3(); figure4(); figure5(); supplementary_figures()
    print("figures complete", flush=True)
