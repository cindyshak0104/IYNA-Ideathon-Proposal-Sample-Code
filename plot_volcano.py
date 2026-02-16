import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text  # pip install adjustText

def volcano_plot_nice(
    csv_path: str,
    out_png: str,
    title: str,
    fdr_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
    label_top_n_each_side: int = 8,
    prefer_gene: bool = True,
    use_pval: bool = False,
    legend_loc: str = 'upper right',
):
    """
    Create a nicely formatted volcano plot.
    
    Args:
        use_pval: If True, use raw p-value instead of FDR (padj).
                 Use when FDR is too conservative (all zeros).
    """
    df = pd.read_csv(csv_path)

    # Ensure y axis exists
    if "neglog10_padj" not in df.columns:
        df["neglog10_padj"] = -np.log10(df["padj"].replace(0, np.nan))

    # Clean
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FC", "neglog10_padj"])
    if "padj" not in df.columns:
        df["padj"] = 10 ** (-df["neglog10_padj"])

    # Determine significance based on use_pval flag
    if use_pval:
        # Use raw p-value instead of FDR
        sig = (df["pval"] <= fdr_thresh) & (df["log2FC"].abs() >= lfc_thresh)
        thresh_label = f"pval<{fdr_thresh}"
    else:
        # Use FDR (padj)
        sig = (df["padj"] <= fdr_thresh) & (df["log2FC"].abs() >= lfc_thresh)
        thresh_label = f"FDR<{fdr_thresh}"
    
    up = sig & (df["log2FC"] >= lfc_thresh)
    down = sig & (df["log2FC"] <= -lfc_thresh)

    # Better styling
    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter plots with better colors and sizes
    ax.scatter(
        df.loc[~sig, "log2FC"], 
        df.loc[~sig, "neglog10_padj"], 
        s=40, alpha=0.35, color='gray', label='Not significant'
    )
    ax.scatter(
        df.loc[up, "log2FC"], 
        df.loc[up, "neglog10_padj"], 
        s=80, alpha=0.8, color='#d62728', edgecolors='darkred', linewidth=0.5, label=f'Up ({thresh_label}, |log2FC|>{lfc_thresh})'
    )
    ax.scatter(
        df.loc[down, "log2FC"], 
        df.loc[down, "neglog10_padj"], 
        s=80, alpha=0.8, color='#1f77b4', edgecolors='darkblue', linewidth=0.5, label=f'Down ({thresh_label}, |log2FC|>{lfc_thresh})'
    )

    # Threshold lines with better styling
    yline = -np.log10(fdr_thresh)
    ax.axhline(yline, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.axvline(lfc_thresh, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.axvline(-lfc_thresh, color='black', linestyle='--', linewidth=1.5, alpha=0.6)

    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel("log2 fold change", fontsize=12, fontweight='bold')
    ax.set_ylabel("-log10(padj)", fontsize=12, fontweight='bold')
    ax.legend(loc=legend_loc, fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.3)

    # Choose labels: top N up + top N down (by smallest padj, tie-break by |log2FC|)
    def pick_label(r):
        if prefer_gene:
            g = str(r.get("gene", "")).strip()
            if g and g.lower() != "nan" and g != "":
                return g
        return str(r.get("id", "")).strip()

    def top_hits(mask, n):
        sub = df.loc[mask].copy()
        if sub.empty:
            return sub
        sub["absLFC"] = sub["log2FC"].abs()
        # Sort by p-value (smallest first), then by |log2FC| (largest first)
        return sub.sort_values(["padj", "absLFC"], ascending=[True, False]).head(n)

    to_label = pd.concat([
        top_hits(up, label_top_n_each_side),
        top_hits(down, label_top_n_each_side)
    ], ignore_index=True)

    # Annotate then repel - only label significant hits
    texts = []
    for _, r in to_label.iterrows():
        label = pick_label(r)
        if label:  # Only add if we have a valid label
            texts.append(
                ax.text(r["log2FC"], r["neglog10_padj"], label, fontsize=9, fontweight='bold')
            )

    # Repel labels (prevents overlap) - no arrows
    if texts:
        adjust_text(
            texts,
            arrowprops=None,  # No arrows - cleaner look
            expand_points=(1.8, 1.8),
            expand_text=(1.5, 1.5),
            force_text=(0.8, 0.8),
            force_points=(0.5, 0.5),
            ax=ax,
        )

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {out_png}")

if __name__ == "__main__":
    # Lung has limited FDR signal; relax thresholds to FDR 0.20
    # Use fewer labels to avoid overlap, legend in lower right
    volcano_plot_nice(
        "lung_TH_vs_TN_volcano.csv",
        "lung_volcano_nice.png",
        "Airway compartment (ASM): HDM vs control (asthma-associated signature)",
        fdr_thresh=0.20,
        lfc_thresh=1.0,
        label_top_n_each_side=5,  # Reduced to prevent overlap
        use_pval=False,
        legend_loc='lower right',  # Move legend to lower right blank space
    )

    # Microglia has no FDR signal; use raw p-value threshold (0.01)
    # with moderate log2FC threshold for exploratory analysis
    volcano_plot_nice(
        "microglia_5xFAD_vs_WT_volcano.csv",
        "microglia_volcano_nice.png",
        "Microglia: 5xFAD vs WT (AD-associated signature, exploratory)",
        fdr_thresh=0.01,  # Use raw p-value (pval < 0.01)
        lfc_thresh=0.5,
        label_top_n_each_side=6,
        use_pval=True,  # Use raw p-value instead of FDR
    )
