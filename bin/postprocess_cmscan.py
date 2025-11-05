# (C) 2024, Maria Schreiber, MIT license
"""
This script processes cmscan output (TSV file), generate a simple mapping file 
that links each input srna gene to the most statistically significant CM from RFAM.

postprocess_cmscan.py \
    -cf cmscan_file.tsv \
    -o output_path
"""
import csv
import argparse

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Postprocessing sRNA clustering')
    parser.add_argument('-cf', '--cmscan_file', required=True, help='Cmscan file')
    parser.add_argument('-o', '--output_file', required=True, help='Path for output file')
    parser.add_argument('-gof', '--gene_of_interest', required=True, help='Name of gene of interest.')
    return parser

def read_file(cmscan_file):
    results = []

    with open(cmscan_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue  # skip header and empty lines
            
            # Split on any whitespace
            row = line.strip().split()
            
            # The fixed columns count is 17 before description (index 0 to 16)
            if len(row) < 18:
                continue  # skip malformed lines
            
            # Join all columns from index 17 to the end as description
            description = ' '.join(row[17:])
            
            result = {
                'target_name': row[0],
                'target_accession': row[1],
                'query_name': row[2],
                'query_accession': row[3],
                'mdl': row[4],
                'mdl_from': row[5],
                'mdl_to': row[6],
                'seq_from': row[7],
                'seq_to': row[8],
                'strand': row[9],
                'trunc': row[10],
                'pass': row[11],
                'gc': row[12],
                'bias': row[13],
                'score': row[14],
                'e_value': row[15],
                'inc': row[16],
                'description': description
            }
            results.append(result)
    
    return results  
    
def filter_hits(results):
    """Keep only the best (lowest E-value) hit for each query_name."""

    best_hits = {}
    for entry in results:
        qname = entry['query_name']
        evalue = float(entry['e_value'].replace('e', 'E'))  # handle scientific notation
        if qname not in best_hits:
            best_hits[qname] = entry
        elif evalue < float(best_hits[qname]['e_value'].replace('e', 'E')):
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
            writer.writerow([entry['target_accession'], entry['query_name'], target_name])

def main():
    parser = setup_parser()
    args = parser.parse_args()
    
    cmscan_lst = read_file(args.cmscan_file)
    filtered_cm_lst = filter_hits(cmscan_lst)
    write_output(filtered_cm_lst, args.output_file, args.gene_of_interest)

if __name__ == '__main__':
    main()
