import argparse
from pathlib import Path
from sugar import read_fts, read, Feature

OUTPUT_FILES = {
    'protein': ('.protein.mfaa', True),
    'srna': ('.srna.mfna', False),
}

FEATURE_TO_OUTPUT_TYPE = {
    # 'protein'
    'protein_coding': 'protein',
    'pseudogene': 'protein',
    'tblastn': 'protein',
    'CDS': 'protein',
    # 'srna' 
    'tRNA': 'srna',
    'rRNA': 'srna',
    'ncRNA': 'srna',
    'infernal': 'srna',
    'misc_RNA': 'srna',
}

BLAST_FMT = 'qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq'

def validate_neighbors(neighbors):
    """Ensure the neighbors parameter is formatted correctly."""
    if ',' not in neighbors and ':' not in neighbors:
        raise ValueError('Invalid neighbors format. Use x,y or x:y')

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Process genomic features and extract neighbors')
    parser.add_argument('--hit_file', help='Path to hit file (blast/infernal). Required if input_type is not "from_gff"')
    parser.add_argument('--input_type', required=True,
                        choices=['blastn', 'blastp', 'infernal', 'tblastn', 'from_gff'],
                        help='Input file type')
    parser.add_argument('--gff_file', required=True, help='Path to GFF/GBK file')
    parser.add_argument('--fna_file', required=True, help='Path to FNA sequence file')
    parser.add_argument('--gene_of_interest', required=True, help='Target gene identifier (e.g., ID, product name or query ID)')
    parser.add_argument('--neighbours', type=str, default='2:2',
                        help='Neighbor range: x,y (genes) or x:y (nucleotides). Default: 2:2')
    parser.add_argument('--output_path', required=True, help='Base path for output files (e.g., /path/to/results/my_run)')
    parser.add_argument('--promoter', type=str, default='no', choices=['yes', 'no'], help='Include promoter region?')
    parser.add_argument('--promoter_len', type=int, default=100, help='Promoter region length (default: 100 nt)')
    parser.add_argument('--including_features', nargs='+', default=['gene', 'pseudogene', 'ncRNA'], 
                        help='Include GFF features like genes, pseudogene, ncRNA. Maybe add CDS or tRNA.')
    parser.add_argument('--evalue_threshold_blast', type=float, default=1, help='E-value threshold for BLAST hits')
    parser.add_argument('--evalue_threshold_infernal', type=float, default=0.1, help='E-value threshold for Infernal hits')
    parser.add_argument('--len_threshold_blast', type=int, default=40, help='Length threshold for BLAST hits')
    parser.add_argument('--len_threshold_infernal', type=int, default=20, help='Length threshold for Infernal hits')
    parser.add_argument('--from_gff_feature', type=str, default="ID", help='Select which feature show match the string. Default: ID, could also be e.g. product or name')

    return parser

def process_hits(hit_file, input_type, gene_of_interest, args):
    """Read and filter hit file (BLAST/Infernal) for the gene of interest."""
    # 1. Determine the general hit typ
    input_type_key = 'blast' if input_type in {'blastn', 'blastp', 'tblastn'} else 'infernal'
    
    # 2. Define the parameter map
    params_map = {
        'blast': {
            'fmt': BLAST_FMT,
            'sort_key': lambda ft: getattr(ft.meta._blast, 'pident', 0),
            'reader': lambda: read_fts(hit_file, 'blast', outfmt=BLAST_FMT, ftype=gene_of_interest)
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

    # 3. Read features
    features = params['reader']()

    # 4. Dynamic filtering based on type
    if input_type_key == 'blast':
        features = features.select(evalue_lt=args.evalue_threshold_blast)
        features = features.select(len_gt=args.len_threshold_blast)
    else:
        features = features.select(evalue_lt=args.evalue_threshold_infernal)
        features = features.select(len_gt=args.len_threshold_infernal)

    # 5. Sort the filtered features
    features.data = sorted(features, key=params['sort_key'], reverse=True)

    return features

def remove_overlapping_features(features):
    """Sort features and remove any that overlap."""
    features.data = sorted(features, key=lambda x: (x.seqid, x.loc.start))
    return features.remove_overlapping()

def get_orientation(strand):
    """Return 'sense' or 'antisense' based on strand symbol."""
    return 'antisense' if strand == '-' else 'sense'

def _find_feature_name(meta: dict, priority_attrs: tuple, default_name: str = "") -> str:
    """Attempts to find the feature name from the metadata using a priority list."""
    for attr in priority_attrs:
        name = meta.get(attr)
        if name:
            return name
    return default_name

def _get_biotype(feature, meta: dict, args) -> str:
    """Determines the biotype based on input type and feature metadata."""
    
    # 1. Special types for hits (tblastn, infernal)
    if args.input_type == 'tblastn' and feature.meta._fmt == 'blast':
        return 'tblastn'
    if args.input_type == 'infernal' and feature.meta._fmt == 'infernal':
        return 'infernal'
        
    # 2. Preferred biotype fields (gene_biotype, feature type, format)
    return (
        meta.get('gene_biotype') or
        getattr(feature, 'type', None) or
        meta.get('_fmt') or
        'unknown' # Fallback
    )

def generate_entry(feature, counter, args):
    """Create a TSV line for a feature, including seqid, name, coordinates, orientation, and biotype."""

    # Retrieving feature.meta._fmt is the original format type (e.g. “gff”, “blast”).
    meta_key = f"_{feature.meta._fmt}"
    meta = getattr(feature.meta, meta_key, {})
    
    # Standard priority list for features outside the GFF product special case
    standard_priority = ('query', 'Name', 'ID', 'qseqid')

    # 1. Decide on a name
    if args.input_type == "from_gff" and args.from_gff_feature == "product":
        # Special case: Search for product name in GFF metadata
        target_product = args.gene_of_interest
        
        if meta.get('product') and target_product in meta['product']:
            name = target_product.replace(" ", "-")
        else:
            # Fallback search within the GFF product special case
            gff_product_priority = ('product', 'Name', 'ID', 'qseqid')
            name = _find_feature_name(meta, gff_product_priority, default_name=meta.get('ID', ''))
            
    else:
        # Standard case: BLAST/Infernal Hits or GFF search by ID
        name = _find_feature_name(meta, standard_priority)
        # Ensure that the name is not empty if ID is available.
        if not name:
             name = meta.get('ID', "")
             
    # 2. Determine your biotype
    gene_biotype = _get_biotype(feature, meta, args)

    # 3. Collect remaining fields
    start, stop = str(feature.loc.start), str(feature.loc.stop)
    orientation = get_orientation(feature.loc.strand)

    # 4. Create TSV line
    fields = [
        f'{feature.seqid}:{counter}', # Unique ID
        name,
        start,
        stop,
        orientation,
        gene_biotype
    ]

    return '\t'.join(fields)

def process_neighborhood(args, merged_features, seqs):
    """ For each gene of interest, find its neighbors and write their data to output files. """
    counter = 0 

    # Use a flag to track if the current feature is a gene of interest
    is_gene_of_interest = False
    for idx, feature in enumerate(merged_features):
        if args.input_type == 'from_gff':
            if args.from_gff_feature == "ID":
                if feature.meta._gff.ID == args.gene_of_interest:
                    is_gene_of_interest = True
            elif args.from_gff_feature == "product":
                try:
                    # Use a more specific exception for missing attributes
                    if args.gene_of_interest in feature.meta._gff.product:
                        is_gene_of_interest = True
                except AttributeError:
                    # Pass silently for features like pseudogenes that lack a 'product'
                    pass
        else: # blast and infernal
            if args.input_type.replace("tblastn","blast") == feature.meta._fmt:
                is_gene_of_interest = True

        # ---- If it's a gene of interest, process its neighbors ----
        if is_gene_of_interest:
            neighbors = get_neighbor_features(merged_features, idx, args.neighbours)
            if neighbors:
                # Generate TSV entry and write data for each neighbor
                for neighbor in neighbors:
                    entry = generate_entry(neighbor, counter, args)
                    write_neighbour_data(
                        neighbor, entry, args.output_path, seqs, args.input_type, args.gene_of_interest
                    )

            # ---- Extract and save promoter regions if requested ----
            if args.promoter == 'yes':
                entry = generate_entry(feature, counter, args)
                # Check for the sequence in the dictionary before passing it
                seq = seqs.d.get(feature.seqid)
                if seq:
                    extract_and_save_promoter_regions(args, feature, seq, entry)

            # Reset the flag and increment the counter for the next iteration
            is_gene_of_interest = False
            counter += 1

def check_overlap(center, overlaps):
    """ Only keep overlaps if the overlap covers >50% of both the center and the neighbor, or if they are on opposite strands. """
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

def get_neighbor_features(features, index, neighbors):
    """ Find neighboring features for the feature at 'index'. If neighbors is 'x,y', return x upstream and y downstream genes.
    If neighbors is 'x:y', return all features within x upstream and y downstream nucleotides. """
    
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
        if overlaps:
            overlaps_filtered = check_overlap(center, overlaps)
            if not overlaps_filtered:  # empty → gene did not pass thresholds
                return []
            result = overlaps_filtered[:]
        else:
            result = [] # no overlapping gene

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

def write_neighbour_data(feature, entry, output_path, seqs, input_type, gene_of_interest):
    """ Write the neighbor's sequence and annotation to the appropriate output files. """
    
    base_path = Path(output_path)

    # Try to get the sequence by dict or fallback to list
    seq = seqs.d.get(feature.seqid)
    if not seq:
        for seq_item in seqs:
            if hasattr(seq_item, 'meta') and hasattr(seq_item.meta, '_fasta'):
                seq = seq_item
                break
        if not seq:
            raise KeyError(f"Sequence for {feature.seqid} not found.")

    translate = True
    # Extract the correct sequence slice
    if input_type == 'tblastn' and feature.meta._fmt == 'blast':
        seq_slice = feature.meta._blast.sseq.replace("_","").replace("-","")
        translate = False
    else:
        seq_slice = seq[feature.loc.start:feature.loc.stop]
        if feature.loc.strand == '-':
            seq_slice = seq_slice.reverse().complement()

    # Prepare header and determine biotype
    header = f'>{entry.replace("\t", "_")}\n'
    bio_type = entry.strip().split('\t')[-1]

    # Write annotation entry to TSV
    with open(base_path.with_suffix('.tsv'), 'a') as tsv_file:
        tsv_file.write(entry + "\n")

    output_type = FEATURE_TO_OUTPUT_TYPE.get(bio_type)

    if output_type and output_type in OUTPUT_FILES:
        suffix, should_translate = OUTPUT_FILES[output_type]
        # Creat outputpath
        output_file_path = base_path.with_suffix(suffix)
        
        with open(output_file_path, 'a') as out_file:
            out_file.write(header)
            
            content = seq_slice
            if output_type == 'protein' and should_translate:
                # Translate Biotyp 
                content = seq_slice.translate(check_start=False, complete=False)

            out_file.write(f'{content}\n')

def extract_and_save_promoter_regions(args, feature, seq, entry):
    """ Extract promoter regions for the gene of interest and write them to a new mfna file. """
    # 1. Define file path
    output_path = Path(args.output_path)
    output_file = output_path.with_suffix('.promoter.mfna')

    # 2. Determine promoter coordinates
    promoter_start = max(0, feature.loc.start - args.promoter_len)
    promoter_end = feature.loc.start - 1

    # 3. Extract promoter sequence and header
    promoter_entry = f'{entry.split()[0]}\t{entry.split()[1]}\t{promoter_start}\t{feature.loc.start-1}\t{entry.split()[4]}\tpromoter\n'
    promoter_header = f'>{promoter_entry.replace("\t", "_").strip()}\n'
    promoter_seq = seq[promoter_start : promoter_end]
    if feature.loc.strand == '-':
        promoter_seq = promoter_seq.reverse().complement()

    # 4. Write data to file
    with open(output_file, 'a') as mfna_file:
        mfna_file.write(promoter_header)
        mfna_file.write(str(promoter_seq) + '\n')

def main():
    # 1. Initialise argument parser and read command line arguments
    parser = setup_parser()
    args = parser.parse_args()

    # 2. Validation of the neighbourhood specification (x,y or x:y)
    validate_neighbors(args.neighbours)

    # 3. Checking the required input file
    if args.input_type != 'from_gff' and not args.hit_file:
        parser.error(f'--hit_file is required for input_type={args.input_type}')

    # 4. Read genome annotation
    gff_features = read_fts(args.gff_file, fmt='gff').select(args.including_features)

    # 5. Read genome sequence
    seqs = read(args.fna_file)

    # 6. Case differentiation based on input type
    # Case A: Input is a hit report (BLAST/Infernal)
    if args.input_type in ['blastn', 'blastp', 'infernal', 'tblastn']:
        # Read and filter hit file using user-specified thresholds
        features_of_gene_of_interest = process_hits(args.hit_file, args.input_type, args.gene_of_interest, args)
        # Remove overlapping hits
        features_of_gene_of_interest = remove_overlapping_features(features_of_gene_of_interest)        
        # Merge BLAST/infernal hits and annotation features
        merged = sorted(list(features_of_gene_of_interest) + list(gff_features), key=lambda x: x.loc.start)
        # Start finding the neighbours for each gene of interest
        process_neighborhood(args, merged, seqs)
    # Case B: Input is directly a GFF feature (e.g. by ID or product name)
    elif args.input_type == 'from_gff':
        # Since no hits need to be processed, the analysis is started directly with the GFF features.
        process_neighborhood(args, gff_features, seqs)
    
    print("\n✅ Done!")

if __name__ == '__main__':
    main()

"""
TODO:
The current logic in get_neighbor_features for “x:y” (nucleotides) is incomplete. It should take into account the boundaries of the genome or not?
"""
