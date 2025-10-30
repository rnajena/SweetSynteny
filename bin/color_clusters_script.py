# (C) 2024, Maria Schreiber, MIT license
"""
This script annotates genomic features (like proteins or sRNAs) with a distinct color derived 
from their cluster membership and generates a statistical summary visualization of the results.

color_clusters_script.py \
    --cluster_file \ 
    --tsv_file \
    --output_file \
    --gene_of_interest
    ( --size_for_cluster )
"""
import argparse
import colorsys
import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Annotate genomic features with cluster colors')
    parser.add_argument('--cluster_file', required=True, help='Input TSV file defining gene cluster membership (based on MMseqs2, cmscan, or hmmscan result).')
    parser.add_argument('--tsv_file', required=True, help='Neighbour TSV file containing features name, locations, etc.')
    parser.add_argument('--output_file', required=True, help='Path for the final output TSV file containing the original features plus the assigned cluster color.')
    parser.add_argument('--gene_of_interest', required=True, help='Identifier (name) of the gene or feature of interest.')
    parser.add_argument('--size_for_cluster', type=int, default=2, help='Minimum number of members a cluster must have to be assigned a unique color.\nSmaller clusters will be assigned the default white color.')
    return parser

def generate_distinct_colors(num_colors):
    """Generate visually distinct colors using HSL color space."""
    colors = ['#FFFFFF']
    hatch_patterns = ['--', '+', 'x', '\\', '///', '\\/...', '*']
    for i in range(num_colors):
        hue = i / num_colors
        lightness = 0.5 + (i % 3) * 0.1
        saturation = 0.7 + (i % 2) * 0.2
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
        base_color = '#{:02x}{:02x}{:02x}'.format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
        
        if base_color not in colors:
            colors.append(base_color)
        else:
            # Base color exists, try hatch patterns to find unique combo
            for hatch in hatch_patterns:
                color_hatch_combo = base_color + '_h_' + hatch
                if color_hatch_combo not in colors:
                    colors.append(color_hatch_combo)
                    break
            else:
                # If all hatch patterns used, fallback: just add base color again (or handle differently)
                colors.append(base_color)

    return colors

def process_clusters(cluster_file, size_for_cluster=2):
    """Process cluster file and create color mappings."""
    cluster_map = {}
    with open(cluster_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for cluster_id, member_id in reader:
            cluster_map.setdefault(cluster_id, []).append(member_id)

    # Filter clusters based on the size_for_cluster parameter
    filtered_clusters = {k: v for k, v in cluster_map.items() if len(v) >= size_for_cluster}
    colors = generate_distinct_colors(len(filtered_clusters))

    # Remove all '#FFFFFF' entries
    colors = [color for color in colors if color != '#FFFFFF']

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

def summary(gene_of_interest, df, output_file):
    # Summary stats
    type_counts = df['type'].value_counts()
    gene_counts = df['gene'].value_counts().get(gene_of_interest, 0)

    # Add gene_of_interest count to the type plot with a custom label
    type_labels = list(type_counts.index) + [f"{gene_of_interest}"]
    type_values = list(type_counts.values) + [gene_counts]

    # Number of unique cluster colors (total)
    color_counts = df['cluster_color'].nunique()

    # Count unique colors per genome and sort by count
    colors_per_genome = df.groupby('genome')['cluster_color'].nunique().sort_values(ascending=False)
    # Calculate mean
    mean_colors_per_genome = colors_per_genome.mean()
    
    # Indices where gene == gene_of_interest
    gene_indices = df.index[df['gene'] == gene_of_interest].tolist()
    distances = []
    for idx in gene_indices:
        # distance to upstream (previous row)
        if idx > 0 and df.loc[idx - 1, 'genome'] == df.loc[idx, 'genome']:
            dist_upstream = df.loc[idx, 'start'] - df.loc[idx - 1, 'end']
            if dist_upstream >= 0:  # distance or gap
                distances.append(dist_upstream)
        # distance to downstream (next row)
        if idx < len(df) - 1 and df.loc[idx + 1, 'genome'] == df.loc[idx, 'genome']:
            dist_downstream = df.loc[idx + 1, 'start'] - df.loc[idx, 'end']
            if dist_downstream >= 0:
                distances.append(dist_downstream)

    # Prepare figure and axes (2x2 grid, 3 plots + 1 empty)
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    axs = axs.flatten()

    # 1) Number of different types + gene of interest count (combined bar plot)
    x = np.arange(len(type_labels))
    bars = axs[0].bar(x, type_values, color='skyblue')
    axs[0].set_xticks(x)
    axs[0].set_xticklabels(type_labels, rotation=45, ha='right')
    axs[0].set_title("Number of different types + gene of interest count")
    axs[0].set_ylabel("Count")
    # add counts above bars
    for bar, val in zip(bars, type_values):
        axs[0].text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + 0.5,
                    str(val),
                    ha='center', va='bottom',
                    fontsize=10, fontweight='bold')

    # 2) Pie chart: distribution of cluster colors with gene names and colored hex slices
    num_unique_colors = df['cluster_color'].nunique()
    df['cluster_color'] = df['cluster_color'].str.split('_h_').str[0]
    color_counts_series = df['cluster_color'].value_counts()
    # Map color to first gene associated (unique)
    color_to_gene = df.drop_duplicates('cluster_color').set_index('cluster_color')['gene']
            
    def get_gene_label(c):
        if c.lower() in ["white", "#ffffff", "#fff"]:
            return "no-cluster"
        return color_to_gene.get(c, '')
    labels = [get_gene_label(c) for c in color_counts_series.index]
    colors = list(color_counts_series.index)  # hex color codes

    wedges, texts, autotexts = axs[1].pie(
        color_counts_series.values,
        labels=labels,
        autopct='%1.1f%%',
        startangle=140,
        colors=colors,
        textprops={'fontsize': 10, 'color': 'black'}
    )

    # Threshold to hide labels less than 5%
    threshold = 3.0
    for wedge, label, pct_label in zip(wedges, texts, autotexts):
        pct_text = pct_label.get_text().rstrip('%')
        if pct_text:  # check if not empty
            pct_value = float(pct_text)
            if pct_value < threshold:
                label.set_text('')       # remove label text
                pct_label.set_text('')   # remove percentage text
    axs[1].set_title(f"Distribution of cluster colors by gene (Total unique colors: {num_unique_colors})", fontsize=14)

    # 3) Number of unique colors per genome (sorted bar plot)
    axs[2].bar(colors_per_genome.index, colors_per_genome.values, color='teal')
    axs[2].set_title(f"Number of unique colors per genome locations\n(Total number of genome locations: {len(colors_per_genome)})")
    axs[2].set_ylabel("Count of unique colors")
    axs[2].set_xlabel("Genome locations")
    axs[2].tick_params(axis='x', rotation=90)
    axs[2].axhline(mean_colors_per_genome, color='red', linestyle='--', linewidth=2, label=f"Mean: {mean_colors_per_genome:.2f}")

    # Add text label on the mean line (position it slightly left and a bit above the line)
    axs[2].text(
            0,  # x position: start of x axis (first bar)
            mean_colors_per_genome + 0.1,  # y position just above the line
            f"Mean: {mean_colors_per_genome:.2f}",
            color='red',
            fontsize=10,
            fontweight='bold'
    )

    if len(colors_per_genome) > 50:
        axs[2].set_xticks([])   # hide if too many
    else:
        axs[2].tick_params(axis='x', rotation=90)

    # 4th subplot
    if distances:
        # Define bins at every 10 units, spanning from min to max distance
        bin_edges = np.arange(0, max(distances) + 10, 10)
        n, bins, patches = axs[3].hist(distances, bins=bin_edges, color='purple', alpha=0.7)
        
        # Bin centers for tick positions
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        axs[3].set_xticks(bin_centers)
        axs[3].set_xticklabels([f"{int(x)}" for x in bin_centers], rotation=90, fontsize=9)

        axs[3].set_title(f"Distance distribution around gene '{gene_of_interest}'")
        axs[3].set_xlabel("Distance (nt)")
        axs[3].set_ylabel("Frequency")
    else:
        axs[3].text(0.5, 0.5, "No adjacent gene distances found", ha='center', va='center', fontsize=12)
        axs[3].set_axis_off()

    if len(colors_per_genome) > 50:
        axs[3].set_xticks([])
    else:
        axs[3].tick_params(axis='x', rotation=90)

    plt.tight_layout()
    temp_path = os.path.dirname(output_file)
    print(temp_path)
    plt.savefig(temp_path + "/summary.png")
    plt.savefig(temp_path + "/summary.svg")

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
    
    # Assign colors # Drop the unnecessary column
    df['cluster_color'] = df['unique_id'].map(color_map).fillna('#FFFFFF')  
    df.drop('unique_id', axis=1, inplace=True)    

    df.to_csv(args.output_file, sep='\t', index=False, header=False)
    
    summary(args.gene_of_interest, df, args.output_file)

if __name__ == '__main__':
    main()
