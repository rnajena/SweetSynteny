# (C) 2024, Maria Schreiber, MIT license
"""
This script processes mmseqEasyLinclust output (TSV file), generate a simple mapping file 
that links each input protein gene to the most statistically significant HMM domain.

postprocess_mmseqEasyLinclust.py \
    -hf mmseqEasyLinclust_output.tsv \
    -o output_path
"""
import csv
import argparse

def setup_parser():
    parser = argparse.ArgumentParser(description='Postprocessing MMseqs2 clustering for SweetSynteny')
    parser.add_argument('-i', '--input', required=True, help='MMseqs2 clusters.tsv file (2 columns)')
    parser.add_argument('-o', '--output', required=True, help='Path for output mapping file')
    parser.add_argument('-gof', '--gene_of_interest', required=True, help='Name of gene of interest (e.g., CrfA)')
    return parser

def extract_gene_name(identifier, gene_of_interest):
    """
    Extracts the gene name from complex strings like:
    NZ_CP026100.1:0_CrfA_1428029_1428165_sense_infernal
    Returns 'CrfA' if GOF is found, otherwise tries to find the name segment.
    """
    if gene_of_interest in identifier:
        return gene_of_interest
    
    # Typical pattern: Split by underscore and look for the segment after the contig:index
    # Format: CONTIG:INDEX_GENENAME_START_END_...
    parts = identifier.split('_')
    if len(parts) > 1:
        # If the first part contains a colon (the contig), the second part is usually the gene name
        return parts[1]
    
    return identifier # Fallback to full string if format is unexpected

def process_mmseqs(input_file, output_file, gene_of_interest):
    with open(input_file, 'r') as f_in, open(output_file, 'w', newline='') as f_out:
        reader = csv.reader(f_in, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')
        
        for row in reader:
            if not row or len(row) < 2:
                continue
            
            cluster_rep = row[0]
            member = row[1]
            
            # Generate the simplified label for the 3rd column
            label = extract_gene_name(cluster_rep, gene_of_interest)
            
            # Output: [Cluster, Member, Label]
            writer.writerow([cluster_rep, member, label])

def main():
    parser = setup_parser()
    args = parser.parse_args()
    
    process_mmseqs(args.input, args.output, args.gene_of_interest)

if __name__ == '__main__':
    main()
