# add_gene_symbols_uniprot.py
# Usage:
#   python add_gene_symbols_uniprot.py lung_TH_vs_TN_volcano.csv lung_with_genes.csv
#   python add_gene_symbols_uniprot.py microglia_5xFAD_vs_WT_volcano.csv microglia_with_genes.csv

import sys
import time
import re
import pandas as pd
import requests

UNIPROT_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_STATUS = "https://rest.uniprot.org/idmapping/status/"
UNIPROT_RESULTS = "https://rest.uniprot.org/idmapping/results/"

def extract_uniprot_ids(id_value: str) -> list[str]:
    """Extract UniProt accessions from an id field (handles 'P12345;Q9XXXX')."""
    if pd.isna(id_value):
        return []
    s = str(id_value)
    parts = re.split(r"[;,\s]+", s.strip())
    # UniProt accessions are typically 6-10 chars; quick filter:
    out = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        # Accept common UniProt accession patterns
        if re.fullmatch(r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]", p) or re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]", p) \
           or re.fullmatch(r"[A-Z0-9]{6,10}", p):
            out.append(p)
    return out

def submit_idmapping(uniprot_ids: list[str]) -> str:
    r = requests.post(
        UNIPROT_RUN,
        data={"from": "UniProtKB_AC-ID", "to": "Gene_Name", "ids": ",".join(uniprot_ids)},
        timeout=60,
    )
    r.raise_for_status()
    return r.json()["jobId"]

def wait_for_job(job_id: str) -> None:
    while True:
        r = requests.get(UNIPROT_STATUS + job_id, timeout=60)
        r.raise_for_status()
        j = r.json()
        if j.get("jobStatus") in ("RUNNING", "NEW"):
            time.sleep(1.0)
            continue
        # When complete, UniProt returns results availability (no RUNNING status)
        return

def fetch_results(job_id: str) -> dict[str, str]:
    """
    Returns mapping dict: UniProt accession -> gene name.
    Note: Some accessions may map to multiple gene names; UniProt returns the gene_name field.
    """
    # Request TSV for simpler parsing
    params = {"format": "tsv"}
    r = requests.get(UNIPROT_RESULTS + job_id, params=params, timeout=60)
    r.raise_for_status()
    text = r.text.strip().splitlines()
    # Expected TSV columns include: From, To
    mapping = {}
    if len(text) < 2:
        return mapping
    header = text[0].split("\t")
    from_i = header.index("From")
    to_i = header.index("To")
    for line in text[1:]:
        cols = line.split("\t")
        if len(cols) <= max(from_i, to_i):
            continue
        mapping[cols[from_i]] = cols[to_i]
    return mapping

def add_gene_column(input_csv: str, output_csv: str):
    df = pd.read_csv(input_csv)

    # Collect unique UniProt IDs
    all_ids = []
    for v in df.get("id", []):
        all_ids.extend(extract_uniprot_ids(v))
    uniq = sorted(set(all_ids))

    if not uniq:
        raise RuntimeError("No UniProt-like accessions found in the 'id' column.")

    print(f"Found {len(uniq)} unique UniProt IDs. Submitting to UniProt API...")

    # UniProt API has practical limits; batch in chunks
    mapping = {}
    chunk_size = 800
    for i in range(0, len(uniq), chunk_size):
        chunk = uniq[i:i+chunk_size]
        print(f"  Submitting chunk {i//chunk_size + 1} ({len(chunk)} IDs)...")
        job = submit_idmapping(chunk)
        print(f"    Job ID: {job}. Waiting for results...")
        wait_for_job(job)
        chunk_mapping = fetch_results(job)
        mapping.update(chunk_mapping)
        print(f"    Got {len(chunk_mapping)} mappings from this chunk.")

    print(f"Total mappings obtained: {len(mapping)}")

    # Map: if multiple IDs per row, take the first that maps
    def map_row(val):
        ids = extract_uniprot_ids(val)
        for u in ids:
            if u in mapping:
                return mapping[u]
        return ""

    # If gene already exists, only fill missing
    if "gene" not in df.columns:
        df["gene"] = ""

    df["gene"] = df["gene"].fillna("")
    fill = df["gene"].astype(str).str.strip().eq("")
    df.loc[fill, "gene"] = df.loc[fill, "id"].apply(map_row)

    df.to_csv(output_csv, index=False)
    print(f"Wrote {output_csv} (filled gene for {fill.sum()} rows where missing)")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_gene_symbols_uniprot.py <input.csv> <output.csv>")
        sys.exit(1)
    add_gene_column(sys.argv[1], sys.argv[2])
