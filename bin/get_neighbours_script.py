import csv
import os
import argparse
from sugar import read_fts, read, Feature

# Output file modes and feature categories
FILE_MODES = {
    'protein': ('.protein.mfaa', {'protein_coding', 'pseudogene', 'blast'}),
    'srna': ('.srna.mfna', {'tRNA', 'infernal', 'ncRNA'})
}

def validate_neighbors(neighbors):
    """Ensure the neighbors parameter is formatted correctly."""
    if ',' not in neighbors and ':' not in neighbors:
        raise ValueError('Invalid neighbors format. Use x,y or x:y')

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
    parser.add_argument('--neighbours', type=str, required=True, default='2:2',
                        help='Neighbor range: x,y (genes) or x:y (nucleotides)')
    parser.add_argument('--output_path', required=True, help='Base path for output files')
    parser.add_argument('--promoter', type=str, default='no', choices=['yes', 'no'], help='Include promoter region?')
    parser.add_argument('--promoter_len', type=int, default=100, help='Promoter region length (default: 100 nt)')
    parser.add_argument('--including_features', type=list, default=['genes', 'pseudogene', 'ncRNA'], 
                        help='Include GFF features like genes, pseudogene, ncRNA. Maybe add CDS or tRNA.')
    parser.add_argument('--evalue_threshold_blast', type=float, default=0.01, help='E-value threshold for BLAST hits')
    parser.add_argument('--evalue_threshold_infernal', type=float, default=0.05, help='E-value threshold for Infernal hits')
    parser.add_argument('--len_threshold_blast', type=int, default=40, help='Length threshold for BLAST hits')
    parser.add_argument('--len_threshold_infernal', type=int, default=20, help='Length threshold for Infernal hits')

    return parser

def process_hits(hit_file, input_type, gene_of_interest, args):
    """Read and filter hit file (BLAST/Infernal) for the gene of interest."""
    input_type_key = 'blast' if input_type in {'blastn', 'blastp', 'tblastn'} else 'infernal'
    blast_fmt = 'qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq'
    
    params_map = {
        'blast': {
            'fmt': blast_fmt,
            'sort_key': lambda ft: getattr(ft.meta._blast, 'pident', 0),
            'reader': lambda: read_fts(hit_file, 'blast', outfmt=blast_fmt, ftype=gene_of_interest)
        },
        'infernal': {
            'fmt': None,
            'sort_key': lambda ft: getattr(ft.meta._infernal, 'evalue', float('inf')),
            'reader': lambda: read_fts(hit_file, 'infernal', ftype=gene_of_interest)
        }
    }

    if input_type_key not in params_map:
        raise ValueError(f'Unsupported input type: {input_type_key}')

    params = params_map[input_type_key]
    features = params['reader']()

    # Use thresholds from args
    if input_type_key == 'blast':
        features = features.select(evalue_lt=args.evalue_threshold_blast)
        features = features.select(len_gt=args.len_threshold_blast)
    else:
        features = features.select(evalue_lt=args.evalue_threshold_infernal)
        features = features.select(len_gt=args.len_threshold_infernal)

    features.data = sorted(features, key=params['sort_key'], reverse=True)
    return features

def remove_overlapping_features(features):
    """Sort features and remove any that overlap."""
    features.data = sorted(features, key=lambda x: (x.seqid, x.loc.start))
    return features.remove_overlapping()

def get_orientation(strand):
    """Return 'sense' or 'antisense' based on strand symbol."""
    return 'antisense' if strand == '-' else 'sense'

def generate_entry(feature, goi, counter):
    """Create a TSV line for a feature, including seqid, name, coordinates, orientation, and biotype."""
    meta = getattr(feature.meta, f"_{feature.meta._fmt}", {})
    name = meta.get('Name', meta.get('ID', goi))
    start, stop = str(feature.loc.start), str(feature.loc.stop)
    orientation = get_orientation(feature.loc.strand)
    gene_biotype = meta.get('gene_biotype', feature.meta._fmt)

    fields = [
        f'{feature.seqid}:{counter}',
        name,
        start,
        stop,
        orientation,
        gene_biotype
    ]
    return '\t'.join(fields) + '\n'

def process_neighborhood(args, merged_features, seqs):
    """
    For each gene of interest, find its neighbors and write their data to output files.
    """
    counter = 0 
    for idx, feature in enumerate(merged_features):
        if feature.type != args.gene_of_interest:
            continue  # Only process features of interest
        neighbors = get_neighbor_features(merged_features, idx, args.neighbours)
        if neighbors:
            for neighbor in neighbors:
                # Generate TSV entry and write data for each neighbor
                entry = generate_entry(neighbor, args.gene_of_interest, counter)
                write_neighbour_data(
                    neighbor, entry, args.output_path, seqs, args.input_type,
                    args.promoter_len, args.promoter, args.gene_of_interest
                )
            counter += 1

def check_overlap(center, overlaps):
    """
    Only keep overlaps if the overlap covers >50% of both the center and the neighbor,
    or if they are on opposite strands.
    """
    new_list = []
    c_start, c_stop = center.loc.start, center.loc.stop
    c_length = c_stop - c_start + 1

    for gene in overlaps:
        g_start, g_end = gene.loc.start, gene.loc.stop
        g_length = g_end - g_start + 1
        overlap_start = max(c_start, g_start)
        overlap_stop = min(c_stop, g_end)
        if overlap_start > overlap_stop:
            continue  # No overlap
        overlap_length = overlap_stop - overlap_start + 1
        center_percent = overlap_length / c_length
        gene_percent = overlap_length / g_length
        # Keep if >50% overlap or on opposite strands
        if center_percent > 0.5 and gene_percent > 0.5:
            new_list.append(gene)
        elif gene.loc.strand != center.loc.strand:
            new_list.append(gene)
    return new_list

def extract_and_save_promoter_regions(args, merged_features, seqs):
    """
    Extract promoter regions for the gene of interest and write them to a new mfna file.
    """
    output_file = f"{args.output_path}.promoter.mfna"
    count = 0
    for feature in merged_features:
        if feature.type != args.gene_of_interest:
            continue

        seqid = feature.seqid
        try:
            seq = seqs.d[seqid]  # Try dict-style access
        except Exception:
            # Fallback for sugar read object (list-style)
            seq = None
            for s in seqs:
                if hasattr(s, "meta") and hasattr(s.meta, "_fasta") and getattr(s.meta._fasta, "header", None) == seqid:
                    seq = s
                    break
            if seq is None:
                continue

        # Get promoter region coordinates
        if feature.loc.strand == "+":
            promoter_start = max(0, feature.loc.start - args.promoter_len)
            promoter_seq = seq[promoter_start:feature.loc.start]
        else:  # negative strand
            promoter_start = feature.loc.stop
            promoter_end = min(len(seq), feature.loc.stop + args.promoter_len)
            promoter_seq = seq[promoter_start:promoter_end].reverse().complement()

        # Write promoter sequence to file
        header = f">{seqid}_promoter_{count} [{promoter_start+1}-{promoter_end}] ({'+' if feature.loc.strand == '+' else '-'})"
        with open(output_file, "a") as out:
            out.write(header + "\n")
            out.write(str(promoter_seq) + "\n")
        count += 1

def get_neighbor_features(features, index, neighbors):
    """
    Find neighboring features for the feature at 'index'.
    If neighbors is 'x,y', return x upstream and y downstream genes.
    If neighbors is 'x:y', return all features within x upstream and y downstream nucleotides.
    """
    center = features[index]
    if ',' in neighbors:
        up, down = map(int, neighbors.split(','))
        overlaps, upstream, downstream = [], [], []
        for f in features:
            if center.seqid != f.seqid or f == center:
                continue
            # Check for overlap
            if not (f.loc.stop < center.loc.start or f.loc.start > center.loc.stop):
                overlaps.append(f)
            elif f.loc.stop < center.loc.start:
                upstream.append(f)
            elif f.loc.start > center.loc.stop:
                downstream.append(f)
        # Filter overlaps
        overlaps_filtered = check_overlap(center, overlaps) if overlaps else []
        if overlaps and not overlaps_filtered:
            return []
        result = overlaps_filtered[:]
        # Sort and select closest upstream and downstream
        upstream = sorted(upstream, key=lambda x: -x.loc.stop)
        downstream = sorted(downstream, key=lambda x: x.loc.start)
        result.extend(upstream[:up])
        result.extend(downstream[:down])
        result.append(center)
        # Return all neighbors sorted by start position
        return sorted(result, key=lambda x: x.loc.start)
    else:
        # Nucleotide-based neighborhood
        up, down = map(int, neighbors.split(':'))
        start_pos = center.loc.start - up
        stop_pos = center.loc.stop + down
        return [ft for ft in features if ft.loc.start >= start_pos and ft.loc.stop <= stop_pos]

def write_neighbour_data(
        feature, entry, output_path, seqs,
        input_type, gene_of_interest
    ):
    """
    Write the neighbor's sequence and annotation to the appropriate output files.
    """
    # Try to get the sequence by dict or fallback to list
    seq = seqs.d.get(feature.seqid)
    if not seq:
        for seq_item in seqs:
            if hasattr(seq_item, 'meta') and hasattr(seq_item.meta, '_fasta'):
                seq = seq_item
                break
        if not seq:
            raise KeyError(f"Sequence for {feature.seqid} not found.")

    # Extract the correct sequence slice
    if input_type == 'tblastn':
        seq_slice = get_tblastn_sequence(seq, feature)
    else:
        seq_slice = seq[feature.loc.start:feature.loc.stop]
        if feature.loc.strand == '-':
            seq_slice = seq_slice.reverse().complement()

    # Prepare header and determine biotype
    header = f'>{entry.replace("\t", "_")}'
    bio_type = entry.strip().split('\t')[-1]

    # Write annotation entry to TSV
    with open(f'{output_path}.tsv', 'a') as tsv_file:
        tsv_file.write(entry)

    # Write sequence to correct output file(s)
    for file_type, (suffix, categories) in FILE_MODES.items():
        if bio_type in categories:
            with open(f'{output_path}{suffix}', 'a') as out_file:
                out_file.write(header)
                # If protein, translate; else, write as is
                content = seq_slice.translate(check_start=False) if file_type == 'protein' else seq_slice
                out_file.write(f'{content}\n')

def get_tblastn_sequence(seq, feature):
    """
    Special handling for tBLASTn hits: extract the correct nucleotide region.
    """
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

    # Read genome sequence
    seqs = read(args.fna_file)
    # Read and filter hit file using user-specified thresholds
    features = process_hits(args.hit_file, args.input_type, args.gene_of_interest, args)
    features = remove_overlapping_features(features)

    # Read genome annotation (GFF/GBK)
    gff_features = read_fts(args.gff_file).select(args.including_features)
    # Merge BLAST/infernal hits and annotation features
    merged = sorted(list(features) + list(gff_features), key=lambda x: x.loc.start)
    # Process neighborhoods and write output
    process_neighborhood(args, merged, seqs)

if __name__ == '__main__':
    main()
