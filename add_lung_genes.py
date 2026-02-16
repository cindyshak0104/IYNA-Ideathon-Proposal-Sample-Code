# add_lung_genes.py
# Map gene names from the original THMPC-Log2Table.txt to the lung volcano CSV

import pandas as pd

# Read original lung data with gene names
df_orig = pd.read_csv("THMPC-Log2Table.txt", sep="\t")
# Create a mapping: protein ID -> gene name
protein_to_gene = dict(zip(df_orig['protein'], df_orig['gene']))

print(f"Loaded {len(protein_to_gene)} protein->gene mappings from THMPC-Log2Table.txt")

# Read volcano file
df_volc = pd.read_csv("lung_TH_vs_TN_volcano.csv")

# Map genes using the protein IDs from the volcano file
df_volc['gene'] = df_volc['id'].map(protein_to_gene)

# Count coverage
has_gene = df_volc['gene'].fillna('').str.strip() != ''
print(f"Proteins with gene names: {has_gene.sum()} / {len(df_volc)} ({has_gene.sum()/len(df_volc)*100:.1f}%)")

# Save
df_volc.to_csv("lung_TH_vs_TN_volcano.csv", index=False)
print("Updated lung_TH_vs_TN_volcano.csv with gene names")

# Show examples
print("\nExample proteins with gene names:")
print(df_volc[has_gene][['id', 'gene', 'log2FC', 'pval']].head(10).to_string())
