# (C) 2024, Maria Schreiber, MIT license
"""
Modified postprocess_cmscan.py
Links Rfam domains from representative ncRNA sequences back to all MMseqs cluster members.
"""
import csv
import argparse
import os

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Postprocessing sRNA clustering with MMseqs mapping')
    parser.add_argument('-cf', '--cmscan_file', required=True, help='Cmscan tblout file')
    parser.add_argument('--mmseqs_map', required=True, help="ncrna_clust_cluster.tsv from MMseqs")
    parser.add_argument('-o', '--output_file', required=True, help='Path for output file')
    parser.add_argument('-gof', '--gene_of_interest', required=True, help='Name of gene of interest.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite the output file instead of appending.')
    return parser


def ensure_input_file(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f'Input file not found: {path}')
    if not os.path.isfile(path):
        raise ValueError(f'Expected a file but found: {path}')


def ensure_output_dir(path):
    output_dir = os.path.dirname(path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

def read_mmseqs_mapping(map_file):
    """
    Reads the MMseqs cluster file for ncRNAs.
    Returns a dict: {representative_id: [member1, member2, ...]}
    """
    mapping = {}
    if not os.path.exists(map_file):
        return mapping
        
    with open(map_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 2: continue
            rep, member = row[0], row[1]
            if rep not in mapping:
                mapping[rep] = []
            mapping[rep].append(member)
    return mapping

def read_file(cmscan_file):
    """Parses the cmscan --tblout output."""
    results = []
    if not os.path.exists(cmscan_file):
        return results

    with open(cmscan_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            row = line.strip().split()
            if len(row) < 18:
                continue
            
            # Index 0: target_name (Rfam ID), Index 2: query_name (Representative ID)
            result = {
                'target_name': row[0],
                'target_accession': row[1],
                'query_name': row[2],
                'e_value': row[15]
            }
            results.append(result)
    
    return results  
    
def filter_hits(results):
    """Keep only the best (lowest E-value) hit for each representative query."""
    best_hits = {}
    for entry in results:
        qname = entry['query_name']
        try:
            evalue = float(entry['e_value'])
        except ValueError:
            continue

        if qname not in best_hits or evalue < float(best_hits[qname]['e_value']):
            best_hits[qname] = entry

    return best_hits

def write_mapped_output(best_hits, mmseqs_mapping, output_file, gene_of_interest, overwrite=False):
    """Writes the results, mapping hits from reps to all cluster members."""
    mode = 'w' if overwrite else 'a'
    with open(output_file, mode, newline='') as f:
        writer = csv.writer(f, delimiter='\t')

        for rep, members in mmseqs_mapping.items():
            if rep not in best_hits:
                continue
            hit = best_hits[rep]
            for member in members:
                if gene_of_interest in member:
                    target_label = gene_of_interest
                else:
                    target_label = hit['target_name']

                writer.writerow([hit['target_accession'], member, target_label])

def main():
    parser = setup_parser()
    args = parser.parse_args()

    ensure_input_file(args.cmscan_file)
    ensure_input_file(args.mmseqs_map)
    ensure_output_dir(args.output_file)

    mmseqs_mapping = read_mmseqs_mapping(args.mmseqs_map)
    cmscan_lst = read_file(args.cmscan_file)
    best_hits_dict = filter_hits(cmscan_lst)
    write_mapped_output(best_hits_dict, mmseqs_mapping, args.output_file, args.gene_of_interest, overwrite=args.overwrite)

if __name__ == '__main__':
    main()
