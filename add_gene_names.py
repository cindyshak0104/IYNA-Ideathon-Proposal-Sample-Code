import pandas as pd
from pyxlsb import open_workbook

# paths (adjust if needed)
XLSB = "microgliareporterabundancedata.xlsb"
VOLC = "microglia_5xFAD_vs_WT_volcano.csv"
OUT  = "microglia_5xFAD_vs_WT_volcano.csv"  # Overwrite original file

# read first sheet rows from the XLSB and build a Uniprot->Gene map
rows = []
with open_workbook(XLSB) as wb:
    sh = wb.get_sheet(wb.sheets[0])
    for i, row in enumerate(sh.rows()):
        rows.append([c.v for c in row])
        if i > 5000:  # enough for the mapping section
            break

# header row containing "Uniprot ID" and "Gene name"
hdr_i = None
for i, r in enumerate(rows[:200]):
    if "Uniprot ID" in r and "Gene name" in r:
        hdr_i = i
        break
if hdr_i is None:
    raise RuntimeError("Couldn't find 'Uniprot ID' / 'Gene name' header row in XLSB.")

hdr = rows[hdr_i]
u_ix = hdr.index("Uniprot ID")
g_ix = hdr.index("Gene name")

mapping = {}
for r in rows[hdr_i+1:]:
    if len(r) <= max(u_ix, g_ix):
        continue
    u = r[u_ix]
    g = r[g_ix]
    if isinstance(u, str) and u and isinstance(g, str) and g:
        mapping[u.strip()] = g.strip()

print(f"Found {len(mapping)} Uniprot ID -> Gene name mappings")

# merge into volcano
df = pd.read_csv(VOLC)
df["gene"] = df["id"].map(mapping).fillna(df.get("gene"))
df.to_csv(OUT, index=False)
print("Wrote:", OUT)
print(f"Proteins with gene names: {(df['gene'].fillna('').str.strip() != '').sum()} / {len(df)}")
