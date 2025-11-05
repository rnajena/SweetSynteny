# (C) 2024, Maria Schreiber, MIT license
"""
This script processes hmmscan output (TSV file), generate a simple mapping file 
that links each input protein gene to the most statistically significant HMM domain.

postprocess_hmmscan.py \
    -hf hmmscan_output.tsv \
    -o output_path
"""
import csv
import argparse

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Postprocessing sRNA clustering')
    parser.add_argument('-hf', '--hmmscan_file', required=True, help='hmmscan file')
    parser.add_argument('-o', '--output_file', required=True, help='Path for output file')
    parser.add_argument('-gof', '--gene_of_interest', required=True, help='Name of gene of interest.')
    return parser

def read_file(hmmscan_file):
    results = []

    with open(hmmscan_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue  # skip header and empty lines
            
            # Split on any whitespace
            row = line.strip().split()
            # Join all columns from index 18 to the end as description
            description = ' '.join(row[18:])
            
            result = {
                'target_name': row[0],
                'target_accession': row[1],
                'query_name': row[2],
                'query_accession': row[3],
                'Evalue_fullSeq': row[4],
                'score_fullSeq': row[5],
                'bias_fullSeq': row[6],
                'Evalue_bestDom': row[7],
                'score_bestDom': row[8],
                'bias_bestDom': row[9],
                'exp': row[10],
                'reg': row[11],
                'clu': row[12],
                'ov': row[13],
                'env': row[14],
                'dom': row[15],
                'rep': row[16],
                'inc': row[17],
                'description': description
            }
            results.append(result)
    
    return results  
    
def filter_hits(results):
    """Keep only the best (lowest E-value) hit for each query_name."""

    best_hits = {}
    for entry in results:
        qname = entry['query_name']
        evalue = float(entry['Evalue_fullSeq'].replace('e', 'E'))  # handle scientific notation
        if qname not in best_hits or evalue < float(best_hits[qname]['Evalue_fullSeq'].replace('e', 'E')):
            best_hits[qname] = entry

    # To get the filtered list:
    filtered_results = list(best_hits.values())

    return filtered_results

def write_output(filtered_results, output_file, gene_of_interest):
    """Write the filtered results to a TSV file with a header."""
    
    with open(output_file, 'a', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for entry in filtered_results:
            if gene_of_interest in entry['query_name']:
                target_name = gene_of_interest
            else:
                target_name = entry['target_name']
            writer.writerow([entry['target_accession'], entry['query_name'], entry['target_name']])

def main():
    parser = setup_parser()
    args = parser.parse_args()
    
    hmmscan_lst = read_file(args.hmmscan_file)
    filtered_hmm_lst = filter_hits(hmmscan_lst)
    write_output(filtered_hmm_lst, args.output_file, args.gene_of_interest)

if __name__ == '__main__':
    main()
