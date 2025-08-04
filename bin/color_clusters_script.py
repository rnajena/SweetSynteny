import pandas as pd
import argparse
import colorsys
import csv
import os

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Annotate genomic features with cluster colors')
    parser.add_argument('--cluster_file', required=True, help='MMseqs2 cluster TSV file')
    parser.add_argument('--tsv_file', required=True, help='Input features TSV file')
    parser.add_argument('--output_file', required=True, help='Path for colored output file')
    parser.add_argument('--size_for_cluster', type=int, default=2, help='Defines at which cluster size, a cluster gets a color.')
    return parser

def generate_distinct_colors(num_colors):
    """Generate visually distinct colors using HSL color space."""
    colors = []
    for i in range(num_colors):
        hue = i / num_colors
        lightness = 0.5 + (i % 3) * 0.1
        saturation = 0.7 + (i % 2) * 0.2
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
        colors.append('#{:02x}{:02x}{:02x}'.format(
            int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)))
    return colors

def process_clusters(cluster_file, size_for_cluster):
    """Process cluster file and create color mappings."""
    cluster_map = {}
    with open(cluster_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for cluster_id, member_id in reader:
            cluster_map.setdefault(cluster_id, []).append(member_id)

    # Filter clusters based on the size_for_cluster parameter
    filtered_clusters = {k: v for k, v in cluster_map.items() if len(v) >= size_for_cluster}
    colors = generate_distinct_colors(len(filtered_clusters))

    # Create direct member-color mapping
    color_map = {}
    for (cluster_id, members), color in zip(filtered_clusters.items(), colors):
        for member in members:
            color_map[member] = color
    return color_map

def validate_files(*files):
    """Ensure all input files exist."""
    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input file not found: {f}")

def main():
    parser = setup_parser()
    args = parser.parse_args()
    
    # Validate input files
    validate_files(args.cluster_file, args.tsv_file)
    
    # Process cluster data and give clusters >2 a color
    color_map = process_clusters(args.cluster_file, args.size_for_cluster)

    # Process neigbours and gene_of_interest data
    df = pd.read_csv(args.tsv_file, sep='\t', header=None,
                     names=['genome', 'gene', 'start', 'end', 
                            'orientation', 'type'])
    df.columns = ['genome', 'gene', 'start', 'end', 'orientation', 'type']
    df['unique_id'] = df.apply(lambda row: '_'.join(row.astype(str)), axis=1)

    # Add color column
    df['cluster_color'] = df['unique_id'].map(color_map).fillna('#FFFFFF')
    df.drop('unique_id', axis=1, inplace=True)    # Drop the unnecessary column

    # Save results
    df.to_csv(args.output_file, sep='\t', index=False, header=False)

if __name__ == '__main__':
    main()
