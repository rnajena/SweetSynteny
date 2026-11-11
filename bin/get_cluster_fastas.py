# (C) 2024, Maria Schreiber, MIT license
"""Get clustered fasta files for gene of interest"""
#!/usr/bin/env python3
import os
import re
import argparse
import pandas as pd
from Bio import SeqIO

def extract_gcf(organism_name):
    match = re.search(r'(GCF_\d+\.\d+)', str(organism_name))
    return match.group(1) if match else None

def parse_fasta_records(fasta_path):
    records_index = {}
    if not os.path.exists(fasta_path):
        return records_index
    for record in SeqIO.parse(fasta_path, "fasta"):
        parts = record.id.split(":")
        if len(parts) >= 4:
            contig = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
                strand = parts[3]
                records_index[(contig, start, end, strand)] = record
            except ValueError:
                continue
    return records_index

def extract_cluster_fastas(input_tsv, neighbour_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(input_tsv, sep="\t")
    
    # Flexible column search for cluster representation
    cluster_col = None
    for col in ["cluster_label", "cluster_id", "cluster"]:
        if col in df.columns:
            cluster_col = col
            break
            
    if cluster_col is None:
        raise ValueError(f"Input TSV must contain a cluster column. Found columns: {list(df.columns)}")

    # Iterate over unique clusters (e.g. 1, 2, 3, 9999)
    for cluster_id, group in df.groupby(cluster_col):
        ncrna_out = os.path.join(output_dir, f"cluster_{cluster_id}_ncrna.mfna")
        protein_out = os.path.join(output_dir, f"cluster_{cluster_id}_protein.mfaa")
        
        ncrna_records = []
        protein_records = []
        loaded_genomes = {}
        
        for row in group.itertuples(index=False):
            org_field = getattr(row, "organism_name", getattr(row, "genome", ""))
            gcf = extract_gcf(org_field)
            if not gcf:
                continue
                
            if gcf not in loaded_genomes:
                ncrna_path = os.path.join(neighbour_dir, f"{gcf}.nb.ncrna.mfna")
                protein_path = os.path.join(neighbour_dir, f"{gcf}.nb.protein.mfaa")
                loaded_genomes[gcf] = {
                    "ncrna": parse_fasta_records(ncrna_path),
                    "protein": parse_fasta_records(protein_path)
                }
            
            contig = str(getattr(row, "contig", ""))
            start = int(getattr(row, "original_start", getattr(row, "start", 0)))
            end = int(getattr(row, "original_end", getattr(row, "end", 0)))
            orientation = str(getattr(row, "orientation", "+"))
            
            key = (contig, start, end, orientation)
            
            if key in loaded_genomes[gcf]["ncrna"]:
                rec = loaded_genomes[gcf]["ncrna"][key]
                rec.id = f"Cluster_{cluster_id}".replace(" ", "")
                ncrna_records.append(rec)
                
            if key in loaded_genomes[gcf]["protein"]:
                rec = loaded_genomes[gcf]["protein"][key]
                rec.id = f"Cluster_{cluster_id}".replace(" ", "")
                protein_records.append(rec)
                
        if ncrna_records:
            SeqIO.write(ncrna_records, ncrna_out, "fasta")
        if protein_records:
            SeqIO.write(protein_records, protein_out, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Extract multi-FASTA per cluster.")
    parser.add_argument("--input_tsv", required=True, help="Path to gene_of_interest.tsv")
    parser.add_argument("--neighbour_dir", required=True, help="Directory containing .nb.ncrna.mfna and .nb.protein.mfaa files")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    args = parser.parse_args()

    extract_cluster_fastas(args.input_tsv, args.neighbour_dir, args.output_dir)

if __name__ == "__main__":
    main()
