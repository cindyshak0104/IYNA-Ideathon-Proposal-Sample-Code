import pandas as pd

def clean_genes(df):
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df[(df["gene"] != "") & (df["gene"].str.lower() != "nan")]
    df["gene"] = df["gene"].str.upper()
    return df

def write_list(df, out_txt):
    genes = df["gene"].drop_duplicates().sort_values()
    genes.to_csv(out_txt, index=False, header=False)
    print(out_txt, "N =", len(genes))

# ----------- FILES -----------
LUNG_CSV = "lung_TH_vs_TN_volcano.csv"
MICRO_CSV = "microglia_5xFAD_vs_WT_volcano.csv"

# ----------- THRESHOLDS -----------
LUNG_PVAL = 0.05  # Use raw p-value (FDR signal is limited)
LUNG_LFC = 1.0

MICRO_TOP_N = 300
MICRO_LFC_MIN = 0.3
# ----------------------------------

# Lung: use raw p-values
lung = pd.read_csv(LUNG_CSV)
lung = clean_genes(lung)
lung = lung[(lung["pval"] <= LUNG_PVAL) & (lung["log2FC"].abs() >= LUNG_LFC)]
write_list(lung, "lung_list.txt")

# Microglia inclusive (ranked)
micro = pd.read_csv(MICRO_CSV)
micro = clean_genes(micro)
micro = micro[micro["log2FC"].abs() >= MICRO_LFC_MIN].copy()
micro["absLFC"] = micro["log2FC"].abs()
micro = micro.sort_values(["pval", "absLFC"], ascending=[True, False]).head(MICRO_TOP_N)
write_list(micro, "microglia_list.txt")
