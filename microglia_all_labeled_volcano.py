import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

CSV = "microglia_5xFAD_vs_WT_volcano.csv"
OUT = "microglia_volcano_match.png"
TITLE = "Microglia: 5xFAD vs WT (p≤0.05, |log2FC|≥0.3)"

# Parameters (match Venn selection)
TOP_N = 300
LFC_MIN = 0.3

# Read full data
df_all = pd.read_csv(CSV)

# Clean and uppercase gene names
if "gene" in df_all.columns:
    df_all["gene"] = df_all["gene"].astype(str).str.strip()
    df_all = df_all[(df_all["gene"] != "") & (df_all["gene"].str.lower() != "nan")]
    df_all["gene"] = df_all["gene"].str.upper()
else:
    df_all["gene"] = df_all["id"].astype(str)

# Compute -log10(pval) for plotting (use raw p-value for microglia)
if "pval" in df_all.columns:
    df_all["neglog10_pval"] = -np.log10(df_all["pval"].replace(0, np.nan))
    ycol = "neglog10_pval"
else:
    # fallback to padj
    df_all["padj"] = df_all.get("padj", np.nan)
    df_all["neglog10_pval"] = -np.log10(df_all["padj"].replace(0, np.nan))
    ycol = "neglog10_pval"

# Drop rows missing essential values
df_all = df_all.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FC", ycol])

# Select top N by p-value among those with |log2FC| >= LFC_MIN
df_candidates = df_all[df_all["log2FC"].abs() >= LFC_MIN].copy()
if "pval" in df_candidates.columns:
    df_top = df_candidates.nsmallest(TOP_N, "pval")
else:
    df_top = df_candidates.nsmallest(TOP_N, ycol)

# split top set by direction
df_top_up = df_top[df_top["log2FC"] >= LFC_MIN].copy()
df_top_down = df_top[df_top["log2FC"] <= -LFC_MIN].copy()
up_count = len(df_top_up)
down_count = len(df_top_down)

# Plot styling and sizes to match provided example
plt.style.use('seaborn-v0_8-darkgrid')
fig, ax = plt.subplots(figsize=(10, 7))

# Not-significant points
sig_mask = (df_all["pval"] <= 0.05) & (df_all["log2FC"].abs() >= LFC_MIN) if "pval" in df_all.columns else (df_all["padj"] <= 0.05) & (df_all["log2FC"].abs() >= LFC_MIN)
ax.scatter(df_all.loc[~sig_mask, "log2FC"], df_all.loc[~sig_mask, ycol], s=40, alpha=0.35, color='gray')

# Colors requested by user
BLUE = '#2F2FFC'   # rgba(47,47,252)
RED = '#FA2F31'    # rgba(250,47,49)
DASH = '#B0B0B4'   # rgba(176,176,180)

# Significant top points colored (fill with bright colors, subtle dark edges)
ax.scatter(df_top_up["log2FC"], df_top_up[ycol], s=80, alpha=0.95, color=RED, edgecolors='#A00000', linewidth=0.5)
ax.scatter(df_top_down["log2FC"], df_top_down[ycol], s=80, alpha=0.95, color=BLUE, edgecolors='#0000A0', linewidth=0.5)

# Threshold lines (use requested dashed color)
yline = -np.log10(0.05)
ax.axhline(yline, color=DASH, linestyle='--', linewidth=1.2, alpha=0.9)
ax.axvline(LFC_MIN, color=DASH, linestyle='--', linewidth=1.2, alpha=0.9)
ax.axvline(-LFC_MIN, color=DASH, linestyle='--', linewidth=1.2, alpha=0.9)

ax.set_title(TITLE, fontsize=14, fontweight='bold', pad=10)
ax.set_xlabel('log2(Fold Change)', fontsize=12, fontweight='bold')
ax.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')

# Label only top 50 by p-value (not all 300) to reduce clutter
# But keep all 300 plotted
LABEL_TOP_N = 50
df_to_label = df_top.nsmallest(LABEL_TOP_N, "pval") if "pval" in df_top.columns else df_top.head(LABEL_TOP_N)

texts = []
for _, r in df_to_label.iterrows():
    texts.append(ax.text(r["log2FC"], r[ycol], r["gene"], fontsize=7, fontweight='bold', color='black'))

if texts:
    adjust_text(
        texts,
        arrowprops=None,
        expand_points=(1.2, 1.2),
        expand_text=(1.0, 1.0),
        force_text=(0.5, 0.5),
        force_points=(0.8, 0.8),
        ax=ax,
    )

# Legend similar to example (UP red, DOWN blue)
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='none', markeredgecolor='none', label=f'UP (n={up_count})', markerfacecolor=RED, markersize=8),
    Line2D([0], [0], marker='o', color='none', markeredgecolor='none', label=f'DOWN (n={down_count})', markerfacecolor=BLUE, markersize=8),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.95)

plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches='tight')
plt.close()
print('Saved', OUT)
