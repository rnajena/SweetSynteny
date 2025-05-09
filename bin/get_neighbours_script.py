import csv
import os
from sugar import read_fts, read, Feature
import argparse

# Constants
EVALUE_THRESHOLDS = {'blast': 0.01, 'infernal': 0.05}
LEN_THRESHOLDS = {'blast': 40, 'infernal': 20}

FILE_MODES = {
    'protein': ('.protein.mfaa', ['protein_coding', 'pseudogene', 'blast']),
    'srna': ('.srna.mfna', ['tRNA', 'infernal', 'ncRNA'])
}

def validate_neighbors(neighbors):
    """Validate neighbors parameter format."""
    if ',' not in neighbors and ':' not in neighbors:
        raise ValueError("Invalid neighbors format. Use x,y or x:y")

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Process genomic features and extract neighbors')
    parser.add_argument('--hit_file', required=True, help='Path to hit file (blast/infernal)')
    parser.add_argument('--input_type', required=True, 
                      choices=['blastn', 'blastp', 'infernal', 'tblastn'],
                      help='Input file type')
    parser.add_argument('--gff_file', required=True, help='Path to GFF or GBK file')
    parser.add_argument('--fna_file', required=True, help='Path to FNA file')
    parser.add_argument('--gene_of_interest', required=True, help='Target gene identifier')
    parser.add_argument('--neighbours', type=str, required=True,
                      help='Neighbor range: x,y (genes) or x:y (nucleotides)')
    parser.add_argument('--output_path', required=True, help='Base path for output files')
    return parser

def process_hits(hit_file, input_type, gene_of_interest):
    """Process input hits file with appropriate parameters."""
    input_type = 'blast' if input_type in {'blastn', 'blastp', 'tblastn'} else input_type
    blast_fmt = 'qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq'
    
    params_map = {
        'blast': {
            'fmt': blast_fmt,
            'sort_key': lambda ft: ft.meta._blast.pident,
            'reader': lambda: read_fts(hit_file, 'blast', outfmt=blast_fmt, ftype=gene_of_interest)
        },
        'infernal': {
            'fmt': None,
            'sort_key': lambda ft: ft.meta._infernal.evalue,
            'reader': lambda: read_fts(hit_file, 'infernal', ftype=gene_of_interest)
        }
    }

    if input_type not in params_map:
        raise ValueError(f"Unsupported input type: {input_type}")

    params = params_map[input_type]
    features = params['reader']()
    features = features.select(evalue_lt=EVALUE_THRESHOLDS[input_type])
    features = features.select(len_gt=LEN_THRESHOLDS[input_type])
    features.data = sorted(features, key=params['sort_key'], reverse=True)
    return features

def remove_overlapping_features(features):
    """Remove overlapping features using efficient algorithm."""
    features.data = sorted(features, key=lambda x: (x.seqid, x.loc.start))
    features = features.remove_overlapping()
    return features

def get_orientation(strand):
    """Determine sequence orientation from strand."""
    return 'antisense' if strand == '-' else 'sense'

def generate_entry(feature, goi, counter):
    """Generate TSV entry for a feature."""
    if feature.meta._fmt == "gff":
        meta = feature.meta._gff
        try:
            name = meta.Name
        except:
            name = meta.ID

        fields = [
            f"{feature.seqid}:{counter}",
            name,
            str(feature.loc.start),
            str(feature.loc.stop),
            get_orientation(feature.loc.strand),
            meta.get('gene_biotype', feature.meta._fmt)
        ]

    elif feature.meta._fmt == "genbank":
        meta = feature.meta._genbank

        fields = [
            f"{feature.seqid}:{counter}",
            meta.get('product', meta.get('ID', goi)),
            str(feature.loc.start),
            str(feature.loc.stop),
            get_orientation(feature.loc.strand),
            meta.get('gene_biotype', feature.meta._fmt)
        ]

    elif feature.meta._fmt == "blast" or feature.meta._fmt == "infernal": #if feature.meta._fmt == "blast":
        fields = [
            f"{feature.seqid}:{counter}",
            goi,
            str(feature.loc.start),
            str(feature.loc.stop),
            get_orientation(feature.loc.strand),
            feature.meta._fmt
        ]
    return '\t'.join(fields) + '\n'

def process_neighborhood(args, merged_features, seqs):
    """Main processing pipeline for neighborhood extraction."""
    counter = 0 
    for idx, feature in enumerate(merged_features):
        if feature.type != args.gene_of_interest:
            continue
            
        neighbors = get_neighbor_features(merged_features, idx, args.neighbours)
        
        for neighbor in neighbors:
            entry = generate_entry(neighbor, args.gene_of_interest, counter)
            write_neighbour_data(neighbor, entry, args.output_path, seqs, args.input_type)
        counter = counter + 1

def get_neighbor_features(features, index, neighbors):
    """Retrieve neighboring features based on parameters."""
    if ',' in neighbors:
        center = features[index]
        up, down = map(int, neighbors.split(','))
        overlaps = []
        upstream = []
        downstream = []
        for f in features:
            # Skip the blast hit itself
            if f == center:
                continue
            # Overlap check
            if not (f.loc.stop < center.loc.start or f.loc.start > center.loc.stop):
                overlaps.append(f)
            elif f.loc.stop < center.loc.start:
                upstream.append(f)
            elif f.loc.start > center.loc.stop:
                downstream.append(f)

        # Sort upstream by end (closest first), downstream by start (closest first)
        upstream = sorted(upstream, key=lambda x: -x.loc.stop)
        downstream = sorted(downstream, key=lambda x: x.loc.start)
        
        # Take the closest up/down
        result = []
        if up > 0 and len(upstream) > 0:
            result.extend(upstream[:up])
        result.extend(overlaps)

        if down > 0 and len(downstream) > 0:
            result.extend(downstream[:down])
        
        result.append(center) # Always include the blast hit itself, in proper order
        return sorted(result, key=lambda x: x.loc.start)

    else:
        up, down = map(int, neighbors.split(':'))
        center = features[index]
        start_pos = center.loc.start - up
        stop_pos = center.loc.stop + down
        return [ft for ft in features if ft.loc.start >= start_pos and ft.loc.stop <= stop_pos]    

def write_neighbour_data(feature, entry, output_path, seqs, input_type):
    """Write feature data to appropriate output files."""

    try:
        seq = seqs.d[feature.seqid]
    except KeyError:
        for seq in seqs:
            seq.id = seq.meta._fasta.header 

    # Get sequence
    if input_type == 'tblastn' and feature.meta._fmt == 'blast':
        seq_slice = get_tblastn_sequence(seq, feature)
    else:
        seq_slice = seq[feature.loc.start:feature.loc.stop]
        if feature.loc.strand == '-':
            seq_slice = seq_slice.reverse().complement()
        
    content = seq_slice

    # Write to files
    header = f">{entry.replace('\t', '_')}"
    bio_type = entry.strip().split('\t')[-1]

    with open(f"{output_path}.tsv", 'a') as tsv_file:
        tsv_file.write(entry)

    for file_type, (suffix, categories) in FILE_MODES.items():
        if bio_type in categories:
            with open(f"{output_path}{suffix}", 'a') as out_file:
                out_file.write(header)
                content = seq_slice.translate(check_start=False) if file_type == 'protein' else seq_slice
                out_file.write(f"{content}\n")

def get_tblastn_sequence(seq, feature):
    """Special handling for tBLASTn sequences."""
    qstart = (feature.meta._blast.qstart - 1) * 3 + 1
    qend = feature.meta._blast.qend * 3
    qlen = feature.meta._blast.qlen * 3
    
    start = max(0, feature.loc.start - qstart)
    end = min(len(seq), feature.loc.stop + (qlen - qend))
    
    sequence = seq[start:end]
    return sequence.reverse().complement() if feature.loc.strand == '-' else sequence

def main():
    parser = setup_parser()
    args = parser.parse_args()
    validate_neighbors(args.neighbours)

    # Data processing pipeline
    seqs = read(args.fna_file)
    features = process_hits(args.hit_file, args.input_type, args.gene_of_interest)
    features = remove_overlapping_features(features)

    # Process genome annotations
    gff_features = read_fts(args.gff_file).select(['gene', 'pseudogene']) #'CDS'
    merged = sorted(features + gff_features, key=lambda x: x.loc.start)
    process_neighborhood(args, merged, seqs)

if __name__ == '__main__':
    main()
