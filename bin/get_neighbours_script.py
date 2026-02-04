# (C) 2024, Maria Schreiber, MIT license
import argparse
from pathlib import Path
from sugar import read_fts, read, Feature

OUTPUT_FILES = {
    'protein': ('.protein.mfaa', True),
    'ncrna': ('.ncrna.mfna', False),
}

FEATURE_TO_OUTPUT_TYPE = {
    # 'protein'
    'protein_coding': 'protein',
    'pseudogene': 'protein',
    'tblastn': 'protein',
    'CDS': 'protein',
    # 'ncrna' 
    'tRNA': 'ncrna',
    'rRNA': 'ncrna',
    'ncRNA': 'ncrna',
    'infernal': 'ncrna',
    'misc_RNA': 'ncrna',
    'RNase_P_RNA': 'ncrna',
}

BLAST_FMT = 'qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq'

def validate_neighbours(neighbours):
    """Ensure the neighbours parameter is formatted correctly."""
    if ',' not in neighbours and ':' not in neighbours:
        raise ValueError('Invalid neighbours format. Use x,y or x:y')

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Process genomic features and extract neighbours\npython get_neighbours_script.py --hit_file hits.txt --input_type {blastn,blastp,infernal,tblastn,from_gff}\n--bio_type {ncRNA,protein} --gff_file GFF_FILE --fna_file FNA_FILE\n--gene_of_interest GENE_OF_INTEREST -output_path OUTPUT_PATH')
    parser.add_argument('--hit_file', help='Path to hit file (blast/infernal). Required if input_type is not "from_gff"')
    parser.add_argument('--input_type', required=True,
                        choices=['blastn', 'blastp', 'infernal', 'tblastn', 'from_gff'],
                        help='Input file type')
    parser.add_argument('--bio_type', required=True,
                        choices=['ncRNA', 'protein'], help='Bio type')
    parser.add_argument('--gff_file', required=True, help='Path to GFF/GBK file')
    parser.add_argument('--fna_file', required=True, help='Path to FNA sequence file')
    parser.add_argument('--gene_of_interest', required=True, help='Target gene identifier (e.g., ID, product name or query ID)')
    parser.add_argument('--neighbours', type=str, default='4,4',
                        help='Neighbor range: x,y (genes) or x:y (nucleotides). Default: 4,4')
    parser.add_argument('--output_path', required=True, help='Base path for output files (e.g., /path/to/results/my_run)')
    parser.add_argument('--promoter', type=str, default='no', choices=['yes', 'no'], help='Include promoter region?')
    parser.add_argument('--promoter_len', type=int, default=100, help='Promoter region length (default: 100 nt)')
    parser.add_argument('--promoter_mode', choices=['fixed','to_next_gene'],
                        default='fixed', help='How to define promoter: fixed length in nt or to next gene boundary')
    parser.add_argument('--including_features', nargs='+', default=['CDS', 'pseudogene', 'ncRNA', 'rRNA', 'tRNA'], 
                        help='Include GFF features like genes, pseudogene, ncRNA. Maybe add CDS or tRNA.')
    parser.add_argument('--evalue_threshold_blast', type=float, default=1, help='E-value threshold for BLAST hits')
    parser.add_argument('--evalue_threshold_infernal', type=float, default=0.1, help='E-value threshold for Infernal hits')
    parser.add_argument('--len_threshold_blast', type=int, default=40, help='Length threshold for BLAST hits')
    parser.add_argument('--len_threshold_infernal', type=int, default=20, help='Length threshold for Infernal hits')
    parser.add_argument('--from_gff_feature', type=str, default="ID", help='Select which feature show match the string. Default: ID, could also be e.g. product or name')
    parser.add_argument('--ignore_overlaps', action='store_true', help='If set, do not filter or exit based on overlapping features.')
    parser.add_argument('--overlap_threshold', type=float, default=0.75, help='Overlap threshold for filtering overlapping features with the gene of interest (default: 0.75)')
    parser.add_argument('--substring_search', action='store_true', help='If set, search for gene_of_interest as a substring (e.g., "SRP" finds "small_SRP").')
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
    if args.input_type == 'blastn' and feature.meta._fmt == 'blast':
        return 'blastn'
    if args.input_type == 'blastp' and feature.meta._fmt == 'blast':
        return 'blastp'
    if args.input_type == 'infernal' and feature.meta._fmt == 'infernal':
        return 'infernal'
        
    # 2. Preferred biotype fields (gene_biotype, feature type, format)
    else:
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
    """ For each gene of interest, find its neighbours and write their data to output files. """
    counter = 0 
    # Track processed locations to avoid duplicates in overlapping SRPs
    processed_locations = set()

    for idx, feature in enumerate(merged_features):
        is_gene_of_interest = False
        
        # 1. Identification logic (remains the same)
        if args.input_type == 'from_gff':
            meta = getattr(feature.meta, '_gff', None)
            if meta:
                target_value = getattr(meta, args.from_gff_feature, "")
                if isinstance(target_value, list):
                    target_value = " ".join(target_value)
                elif target_value is None:
                    target_value = ""

                if args.substring_search:
                    if args.gene_of_interest in target_value:
                        is_gene_of_interest = True
                else:
                    if args.gene_of_interest == target_value:
                        is_gene_of_interest = True
        else:
            if args.input_type.replace("tblastn","blast") == feature.meta._fmt:
                is_gene_of_interest = True
            if args.input_type.replace("blastn","blast") == feature.meta._fmt:
                is_gene_of_interest = True
            if args.input_type.replace("blastp","blast") == feature.meta._fmt:
                is_gene_of_interest = True

        # 2. Check if we have already processed this specific genomic spot
        # We use a tuple of (seqid, start, stop) as a unique key
        location_key = (feature.seqid, feature.loc.start, feature.loc.stop)
        
        # Additional check: If ignore_overlaps is TRUE, we might want to skip 
        # features that are nested inside one another to avoid duplicates.
        is_duplicate_location = False
        for loc in processed_locations:
            # Check if current feature is inside an already processed range
            # or vice versa (Overlap check)
            if feature.seqid == loc[0]:
                if not (feature.loc.stop < loc[1] or feature.loc.start > loc[2]):
                    is_duplicate_location = True
                    break

        # ---- Process if it's a match AND we haven't seen this area ----
        if is_gene_of_interest and not is_duplicate_location:

            neighbours = get_neighbor_features(merged_features, idx, args.neighbours, args.overlap_threshold, args.ignore_overlaps)
            if neighbours:
                for neighbor in neighbours:
                    entry = generate_entry(neighbor, counter, args)
                    write_neighbour_data(
                        neighbor, entry, args.output_path, seqs, args.input_type, args.bio_type, args.gene_of_interest
                    )
                
                # Mark this location as processed
                processed_locations.add(location_key)
                counter += 1
            # inside the loop where you call extract_and_save_promoter_regions
            if args.promoter == 'yes':
                coords = get_promoter_coords(feature, merged_features, args.promoter_len, args.promoter_mode)
                if coords:
                    start_pos, end_pos = coords
                    seq = seqs.d.get(feature.seqid)
                    if seq:
                        extract_and_save_promoter_regions(args, feature, seq, entry, start_pos, end_pos, Path(args.output_path))

def check_overlap(center, overlaps, overlap_threshold):
    """ Only keep overlaps if the overlap covers >75% of both the center and the neighbor, or if they are on opposite strands. """
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
        if center_percent > overlap_threshold and gene_percent > overlap_threshold:
            new_list.append(gene)
        
        elif gene.loc.strand != center.loc.strand:
            new_list.append(gene)
    
    return new_list

def get_neighbor_features(features, index, neighbours, overlap_threshold, ignore_overlaps=False):
    """ Find neighboring features for the feature at 'index'. If neighbours is 'x,y', return x upstream and y downstream genes.
    If neighbours is 'x:y', return all features within x upstream and y downstream nucleotides. """
    
    center = features[index]   
    
    if ',' in neighbours:
        up, down = map(int, neighbours.split(','))

        if center.loc.strand == '-':
            # Swap the counts so 'up' refers to higher coordinates 
            # and 'down' refers to lower coordinates.
            up, down = down, up
            
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

        result = []
        
        # Logic for overlaps
        if overlaps:
            if ignore_overlaps:
                # Treat overlaps as valid neighbours to include in output
                result.extend(overlaps)
            else:
                overlaps_filtered = check_overlap(center, overlaps, overlap_threshold=0.75)
                # If they pass the 50% rule, include them
                if overlaps_filtered:
                    result.extend(overlaps_filtered)
                # Note: We no longer return [] here if overlaps_filtered is empty

        # Sort and select closest upstream and downstream
        upstream = sorted(upstream, key=lambda x: -x.loc.stop)
        downstream = sorted(downstream, key=lambda x: x.loc.start)

        result.extend(upstream[:up])
        result.extend(downstream[:down])
        result.append(center)

        return sorted(result, key=lambda x: x.loc.start)

    else:
        # Nucleotide-based neighborhood
        up, down = map(int, neighbours.split(':'))
        if center.loc.strand == '-':
            # Swap the counts so 'up' refers to higher coordinates 
            # and 'down' refers to lower coordinates.
            up, down = down, up
        start_pos = center.loc.start - up
        stop_pos = center.loc.stop + down
        return [ft for ft in features if ft.loc.start >= start_pos and ft.loc.stop <= stop_pos]

def write_neighbour_data(feature, entry, output_path, seqs, input_type, bio_type, gene_of_interest):
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
    tsv_path = str(base_path) + '.tsv'
    with open(tsv_path, 'a', newline='') as tsv_file: 
        tsv_file.write(entry + "\n")

    output_type = FEATURE_TO_OUTPUT_TYPE.get(bio_type)

    if output_type and output_type in OUTPUT_FILES:
        suffix, should_translate = OUTPUT_FILES[output_type]
        # Creat outputpath
        output_file_path = str(base_path) + suffix
        
        with open(output_file_path, 'a') as out_file:
            out_file.write(header)
            
            content = seq_slice
            if output_type == 'protein' and should_translate:
                # Translate Biotyp 
                content = seq_slice.translate(check_start=False, complete=False)

            out_file.write(f'{content}\n')

def get_promoter_coords(center, features, promoter_len, mode='fixed'):
    """Return (start, end) inclusive promoter coords or None if invalid."""
    seqid = center.seqid
    # gather features on same seq sorted by start
    same = [f for f in features if f.seqid == seqid and f != center]
    same = sorted(same, key=lambda x: x.loc.start)

    if center.loc.strand == '+':
        promoter_end = center.loc.start - 1
        if promoter_end < 0:
            return None
        if mode == 'fixed':
            promoter_start = max(0, center.loc.start - promoter_len)
        else:  # to_next_gene
            # find nearest upstream gene ending before center.start
            upstream = [f for f in same if f.loc.stop < center.loc.start]
            promoter_start = upstream[-1].loc.stop + 1 if upstream else 0

    else:  # '-' strand: promoter is downstream of feature.stop
        promoter_start = center.loc.stop + 1
        if mode == 'fixed':
            promoter_end = promoter_start + promoter_len - 1
        else:  # to_next_gene
            # find nearest downstream gene starting after center.stop
            downstream = [f for f in same if f.loc.start > center.loc.stop]
            promoter_end = downstream[0].loc.start - 1 if downstream else None

    if promoter_end is None or promoter_start > promoter_end:
        return None
    return (promoter_start, promoter_end)

def extract_and_save_promoter_regions(args, feature, seq, entry, promoter_start, promoter_end, output_base: Path):
    """Write promoter seq defined by promoter_start..promoter_end inclusive."""
    promoter_entry = f'{entry.split()[0]}\t{entry.split()[1]}\t{promoter_start}\t{promoter_end}\t{entry.split()[4]}\tpromoter\n'
    promoter_header = f'>{promoter_entry.replace("\t", "_").strip()}\n'
    promoter_seq = seq[promoter_start : promoter_end]
    if feature.loc.strand == '-':
        promoter_seq = promoter_seq.reverse().complement()

    promoter_path = output_base.with_suffix(output_base.suffix + '.promoter.mfna')
    with promoter_path.open('a') as mfna_file:
        mfna_file.write(promoter_header)
        mfna_file.write(str(promoter_seq) + '\n')

def main():
    # 1. Initialise argument parser and read command line arguments
    parser = setup_parser()
    args = parser.parse_args()
        
    # 2. Validation of the neighbourhood specification (x,y or x:y)
    validate_neighbours(args.neighbours)

    # 3. Checking the required input file
    if args.input_type != 'from_gff' and not args.hit_file:
        parser.error(f'--hit_file is required for input_type={args.input_type}')

    # 4. Read genome annotation
    raw_gff = read_fts(args.gff_file, fmt='gff')
    gff_features = raw_gff.select(args.including_features)
    from collections import Counter
    print("All types in GFF:", Counter(getattr(f, 'type', None) for f in raw_gff))

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
    
    print("\n Done!")

if __name__ == '__main__':
<<<<<<< HEAD
    main()
=======
    main()
>>>>>>> 670fee0858c4ff7c52011902a271b3ba90d047dc
