# proteomics_volcano.py
# Run in IDLE: Run > Run Module (F5)
#
# Outputs volcano-ready CSVs with columns:
# id, gene, log2FC, pval, padj, neglog10_padj, n_group1, n_group2

import re
import math
import numpy as np
import pandas as pd
from scipy import stats
from pyxlsb import open_workbook



# -----------------------------
# Utilities
# -----------------------------
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR. Returns adjusted p-values (q-values)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    q = np.full(n, np.nan)
    ok = np.isfinite(p)
    p_ok = p[ok]
    m = p_ok.size
    if m == 0:
        return q
    order = np.argsort(p_ok)
    ranks = np.arange(1, m + 1)
    q_ok = (p_ok[order] * m) / ranks
    # enforce monotonicity
    q_ok = np.minimum.accumulate(q_ok[::-1])[::-1]
    q_ok = np.clip(q_ok, 0, 1)
    q[ok] = q_ok[np.argsort(order)]
    return q

def neglog10(x: float) -> float:
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return -math.log10(x)

def filter_maxquant_flags(df: pd.DataFrame) -> pd.DataFrame:
    # MaxQuant uses "+" in these columns for rows to remove.
    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in df.columns:
            df = df[df[col].fillna("") != "+"]
    return df

def filter_common_flags(df: pd.DataFrame) -> pd.DataFrame:
    # If your file is from MaxQuant, these flags exist. If not, harmless.
    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in df.columns:
            df = df[df[col].fillna("") != "+"]
    return df

def to_numeric_series(s: pd.Series) -> pd.Series:
    """Convert to numeric, handling commas and scientific notation."""
    s = s.astype(str).str.replace(",", "", regex=False)
    return pd.to_numeric(s, errors="coerce")

def compute_two_group_volcano(
    df: pd.DataFrame,
    id_col: str,
    gene_col: str | None,
    group1_cols: list[str],
    group2_cols: list[str],
    already_log2: bool,
    zeros_as_missing: bool = True,
    min_non_missing: int = 2,
    group1_name: str = "group1",
    group2_name: str = "group2",
) -> pd.DataFrame:
    """
    Computes log2FC = mean(group1) - mean(group2) on log2 scale,
    Welch t-test per protein, BH-FDR.
    """
    out = []
    pvals = []

    for _, row in df.iterrows():
        pid = str(row.get(id_col, "")).strip()
        if not pid or pid.lower() == "nan":
            continue

        gene = ""
        if gene_col and gene_col in df.columns:
            gene = str(row.get(gene_col, "")).strip()
            if gene.lower() == "nan":
                gene = ""

        x1 = to_numeric_series(row[group1_cols])
        x2 = to_numeric_series(row[group2_cols])

        if zeros_as_missing:
            x1 = x1.replace(0, np.nan)
            x2 = x2.replace(0, np.nan)

        if not already_log2:
            # intensities are not log2; transform positive values only
            x1 = x1.where(x1 > 0, np.nan)
            x2 = x2.where(x2 > 0, np.nan)
            x1 = np.log2(x1)
            x2 = np.log2(x2)

        x1 = x1[np.isfinite(x1)].to_numpy(dtype=float)
        x2 = x2[np.isfinite(x2)].to_numpy(dtype=float)

        if x1.size < min_non_missing or x2.size < min_non_missing:
            continue

        log2fc = float(np.mean(x1) - np.mean(x2))
        _, pval = stats.ttest_ind(x1, x2, equal_var=False, nan_policy="omit")

        out.append({
            "id": pid,
            "gene": gene,
            "log2FC": log2fc,
            "pval": float(pval) if np.isfinite(pval) else np.nan,
            "n_group1": int(x1.size),
            "n_group2": int(x2.size),
            "group1": group1_name,
            "group2": group2_name,
        })
        pvals.append(pval)

    res = pd.DataFrame(out)
    if res.empty:
        return res

    res["padj"] = bh_fdr(res["pval"].to_numpy())
    res["neglog10_padj"] = res["padj"].apply(neglog10)
    return res.sort_values("padj", na_position="last")

# ========== MICROGLIA: PerseusExport_proteomicRulerTransformation_v01.xlsb ==========
# 
# Uses the "Original Normalized TMT Reporter Abundances (with Imputation)" block
# with columns: AD1 AD2 AD3 B6.1 B6.2 B6.3 (ignore LPS/GIS)
#

def read_xlsb_sheet_as_rows(xlsb_path: str, sheet_index: int = 0):
    with open_workbook(xlsb_path) as wb:
        sh = wb.get_sheet(wb.sheets[sheet_index])
        for row in sh.rows():
            yield [c.v for c in row]

def run_microglia_xlsb(microglia_xlsb: str, out_csv: str):
    # We treat AD = 5xFAD, B6 = WT
    FAD_COLS = ["AD1", "AD2", "AD3"]
    WT_COLS  = ["B6.1", "B6.2", "B6.3"]

    rows = list(read_xlsb_sheet_as_rows(microglia_xlsb, sheet_index=0))

    # Find the header row that contains AD1 and B6.1
    header_idx = None
    for i in range(min(200, len(rows))):
        r = rows[i]
        if any(x == "AD1" for x in r) and any(x == "B6.1" for x in r):
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError("Microglia: could not find header row containing AD1 and B6.1.")

    header = rows[header_idx]
    col_index = {name: j for j, name in enumerate(header) if isinstance(name, str)}

    missing = [c for c in (FAD_COLS + WT_COLS) if c not in col_index]
    if missing:
        raise RuntimeError(f"Microglia: missing expected columns in XLSB header: {missing}")

    # Heuristic: UniProt-like ID somewhere in row (e.g., P46467 / Q9D7Z6)
    def find_uniprot_token(r):
        for x in r:
            if isinstance(x, str) and 5 <= len(x) <= 10 and x[0].isalpha() and any(ch.isdigit() for ch in x):
                # avoid sample labels and short words
                if x in FAD_COLS or x in WT_COLS:
                    continue
                return x
        return None

    data = []
    for r in rows[header_idx + 1:]:
        if not r or len(r) < max(col_index.values()) + 1:
            continue
        uid = find_uniprot_token(r)
        if uid is None:
            continue
        row_dict = {"id": uid, "gene": ""}

        for c in FAD_COLS + WT_COLS:
            row_dict[c] = r[col_index[c]]
        data.append(row_dict)

    df = pd.DataFrame(data)
    for c in FAD_COLS + WT_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Compute volcano: log2FC = mean(AD) - mean(B6)
    res = compute_two_group_volcano(
        df, id_col="id", gene_col="gene",
        group1_cols=FAD_COLS, group2_cols=WT_COLS,
        already_log2=False,  # values are abundances, not log2
        zeros_as_missing=True,
        group1_name="Microglia_5xFAD(AD)", group2_name="Microglia_WT(B6)"
    )

    res.to_csv(out_csv, index=False)
    print(f"[MICROGLIA] wrote {out_csv} ({len(res)} proteins)")

if __name__ == "__main__":
    run_microglia_xlsb(
        microglia_xlsb="microgliareporterabundancedata.xlsb",
        out_csv="microglia_5xFAD_vs_WT_volcano.csv"
    )
