# (C) 2024, Maria Schreiber, MIT license
'''
This script annotates genomic features (like proteins or sRNAs) with a distinct color derived 
from their cluster membership and generates a statistical summary visualization of the results.

color_clusters_script.py \
    --cluster_file \
    --tsv_file \
    --output_file \
    --gene_of_interest
    ( --size_for_cluster )
'''

import argparse
import colorsys
import csv
import json
import matplotlib
matplotlib.use('Agg')  # This tells Matplotlib NOT to use Qt/X11
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import os
import pandas as pd

def setup_parser():
    '''Configure command-line argument parser.'''
    parser = argparse.ArgumentParser(description='Annotate genomic features with cluster colors')
    parser.add_argument('--cluster_file', required=True, help='Input TSV file defining gene cluster membership (based on MMseqs2, cmscan, or hmmscan result).')
    parser.add_argument('--tsv_file', required=True, help='Neighbour TSV file containing features name, locations, etc.')
    parser.add_argument('--output_file', required=True, help='Path for the final output TSV file containing the original features plus the assigned cluster color.')
    parser.add_argument('--gene_of_interest', required=True, help='Identifier (name) of the gene or feature of interest.')
    parser.add_argument('--goi_type', required=True, help='Type of the gene or feature of interest.')
    parser.add_argument('--size_for_cluster', type=int, default=2, help='Minimum number of members a cluster must have to be assigned a unique color.\nSmaller clusters will be assigned the default white color.')
    return parser

def generate_distinct_colors(num_colors):
    '''Generate visually distinct colors using HSL color space and hatch patterns.'''
    colors = ['#FFFFFF', '#FF0000']
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
            for hatch in hatch_patterns:
                color_hatch_combo = base_color + '_h_' + hatch
                if color_hatch_combo not in colors:
                    colors.append(color_hatch_combo)
                    break
            else:
                colors.append(base_color)

    return colors

def process_clusters(cluster_file, gene_of_interest, size_for_cluster=2):
    '''Process cluster file and create color mappings.'''
    cluster_map = {}
    gene_name_map = {}
    goi_members = []

    with open(cluster_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row_index, row in enumerate(reader):
            if len(row) < 3:
                print(f'Warning: Skipping row {row_index} due to insufficient columns: {row}')
                continue

            cluster_id, member_id, gene_name = row[0], row[1], row[2]
            if gene_name == gene_of_interest:
                goi_members.append(member_id)

            cluster_map.setdefault(cluster_id, []).append(member_id)
            gene_name_map[member_id] = gene_name

    filtered_clusters = {k: v for k, v in cluster_map.items() if len(v) >= size_for_cluster}
    colors = generate_distinct_colors(len(filtered_clusters))
    colors = [color for color in colors if color != '#FFFFFF']

    color_map = {}
    for (cluster_id, members), color in zip(filtered_clusters.items(), colors):
        for member in members:
            if color != '#FF0000':
                color_map[member] = color

    for member_id in goi_members:
        color_map[member_id] = '#FF0000'

    return color_map, gene_name_map

def validate_files(*files):
    '''Ensure all input files exist.'''
    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f'Input file not found: {f}')

def summary(gene_of_interest, goi_typ, df, output_file):
    type_counts = df['type'].value_counts().to_dict() # Get counts for all types
    
    # Create the labels and values lists
    type_labels = []
    type_values = []
    
    for t, count in type_counts.items():
        if t == goi_typ:
            # Add the gene name to the category label (e.g., 'infernal (crfA)')
            type_labels.append(f'{t} ({gene_of_interest})')
        else:
            type_labels.append(t)
        type_values.append(count)

    # Count unique colors per genome and sort by count
    colors_per_genome = df.groupby('genome')['cluster_color'].nunique().sort_values(ascending=False)
    # Calculate mean
    mean_colors_per_genome = colors_per_genome.mean()
    
    # Indices where gene == gene_of_interest
    gene_indices = df.index[df['type'] == goi_typ].tolist()
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
    
    # --- Define Color Scheme ---
    cmap_qualitative = cm.get_cmap('Pastel1') # For types
    accent_color = '#E63946' # A bold red for highlights like the mean and gene_of_interest-

    # Prepare figure and axes (2x2 grid, 3 plots + 1 empty)
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    axs = axs.flatten()

    # 1) Number of different types + gene of interest count (combined bar plot)
    x = np.arange(len(type_labels))
    
    # Generate colors from the qualitative map for the types
    base_colors = cmap_qualitative(np.linspace(0, 1, len(type_labels) - 1))
    plot1_colors = list(base_colors) + [accent_color] # Use accent_color for the gene_of_interest bar
    
    bars = axs[0].bar(x, type_values, color=plot1_colors)
    
    # Highlight the gene_of_interest bar with an accent color
    bars[-1].set_color(accent_color) 
    
    axs[0].set_xticks(x)
    axs[0].set_xticklabels(type_labels, rotation=45, ha='right')
    axs[0].set_title('A) Feature Type Counts (including gene of interest)', fontsize=14, fontweight='bold')
    axs[0].set_ylabel('Count', fontsize=12)
    
    # add counts above bars
    for i, (bar, val) in enumerate(zip(bars, type_values)):
        color = 'black' if i < len(bars) - 1 else accent_color
        axs[0].text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + 0.5,
                    str(val),
                    ha='center', va='bottom',
                    fontsize=10, fontweight='bold',
                    color=color)

    # 2) Pie chart: distribution of cluster colors with gene names and colored hex slices
    df_for_pie = df.copy()
    num_unique_colors = df_for_pie['cluster_color'].nunique()
    # Ensure the part before '_h_' is the color code
    df_for_pie['cluster_color'] = df_for_pie['cluster_color'].str.split('_h_').str[0] 
    color_counts_series = df_for_pie['cluster_color'].value_counts()
    
    # Map color to first gene associated (unique) - using the original column name
    color_to_gene = df_for_pie.drop_duplicates('cluster_color').set_index('cluster_color')['gene_name_based_on_DB']

    def get_gene_label(c):
        # 1. Handle 'no-cluster' case
        if c.lower() in ['white', '#ffffff', '#fff']:
            return 'no-cluster'
            
        subset = df_for_pie[df_for_pie['cluster_color'] == c]
        if (subset['gene'] == gene_of_interest).any():
            return gene_of_interest
            
        return color_to_gene.get(c, 'Unknown')
        
    labels = [get_gene_label(c) for c in color_counts_series.index]
            
    colors = list(color_counts_series.index)  # Use the hex color codes directly for true representation

    wedges, texts, autotexts = axs[1].pie(
        color_counts_series.values,
        labels=labels,
        autopct='%1.1f%%',
        startangle=90, # Start at the top
        colors=colors,
        textprops={'fontsize': 10, 'color': 'black'},
        wedgeprops={'edgecolor': 'black', 'linewidth': 0.5, 'antialiased': True} # Better separation
    )

    # Threshold to hide labels less than 3% (kept original threshold)
    threshold = 3.0
    for _, label, pct_label in zip(wedges, texts, autotexts):
        pct_text = pct_label.get_text().rstrip('%')
        if pct_text:  # check if not empty
            pct_value = float(pct_text)
            if pct_value < threshold:
                label.set_text('')       # remove label text
                pct_label.set_text('')   # remove percentage text

    axs[1].set_title(f'B) Cluster Color Distribution (Total unique colors: {num_unique_colors})', fontsize=14, fontweight='bold')
    axs[1].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    # 3) Number of unique colors per genome (sorted bar plot)
    # Generate colors from the quantitative map (Viridis) - can represent magnitude
    num_genomes = len(colors_per_genome)
    bar_color_static = '#CCCCCC' # Define light grey
    axs[2].bar(colors_per_genome.index, colors_per_genome.values, color=bar_color_static)
    
    axs[2].set_title(f'C) Unique Colors per Genome Location\n(Total Locations: {num_genomes})', fontsize=14, fontweight='bold')
    axs[2].set_ylabel('Count of unique colors', fontsize=12)
    axs[2].set_xlabel('Genome Locations', fontsize=12)
    
    # Mean line using the accent color
    axs[2].axhline(mean_colors_per_genome, color=accent_color, linestyle='--', linewidth=2, label=f'Mean: {mean_colors_per_genome:.2f}')

    # Add text label on the mean line (position it slightly left and a bit above the line)
    right_x_pos = num_genomes - 1 
    axs[2].text(
            right_x_pos,  # x position: End of the x axis (last bar position)
            mean_colors_per_genome + 0.1,  # y position just above the line
            f'Mean: {mean_colors_per_genome:.2f}',
            color='black', # Changed color to black
            fontsize=10,
            fontweight='bold',
            ha='right' # Aligned to the right
    )

    if len(colors_per_genome) > 50:
        axs[2].set_xticks([])   # hide if too many
    else:
        axs[2].tick_params(axis='x', rotation=90)

    # 4th subplot: Distance distribution
    if distances:
        # Define bins at every 10 units
        bin_edges = np.arange(0, max(distances) + 10, 10)
        
        plot4_color = '#CCCCCC' # Light grey
        n, bins, patches = axs[3].hist(distances, bins=bin_edges, color=plot4_color, alpha=0.9, edgecolor='black')
        
        # Bin centers for tick positions
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        
        # Reduce x axis ticks by sampling every 5th bin center
        # Select every 5th tick and label
        sampled_tick_indices = np.arange(0, len(bin_centers), 5) 
        sampled_bin_centers = bin_centers[sampled_tick_indices]
        sampled_bin_labels = [f'{int(x)}' for x in sampled_bin_centers]
        
        axs[3].set_xticks(sampled_bin_centers)
        axs[3].set_xticklabels(sampled_bin_labels, rotation=90, fontsize=9)

        axs[3].set_title(f'D) Distance Distribution around gene of interest', fontsize=14, fontweight='bold')
        axs[3].set_xlabel('Intergenic Distance (nt)', fontsize=12)
        axs[3].set_ylabel('Frequency', fontsize=12)
    else:
        axs[3].text(0.5, 0.5, 'No adjacent gene distances found', ha='center', va='center', fontsize=12, color='gray')
        axs[3].set_title('D) Distance Distribution', fontsize=14, fontweight='bold')
        axs[3].set_axis_off()

    # Final layout adjustments
    plt.tight_layout(pad=3.0)
    
    # Saving (kept original saving logic)
    temp_path = os.path.dirname(output_file)
    if not temp_path:
         temp_path = '.' # ensure path is valid if output_file is just a name

    plt.savefig(os.path.join(temp_path, 'summary.png'))
    plt.savefig(os.path.join(temp_path, 'summary.svg'))

    pie_labels = [get_gene_label(c) for c in color_counts_series.index]
    pie_values = [int(v) for v in color_counts_series.values]
    pie_colors = list(color_counts_series.index) # Hex codes
    
    # Export stats to a JSON file that the HTML can read
    stats_data = {
        'plot1_types': {
            'labels': [str(l) for l in type_labels],
            'values': [int(v) for v in type_values]
        },
        'plot2_pie': {
            'labels': pie_labels,
            'values': pie_values,
            'colors': pie_colors
        },
        'plot3_genomes': {
            'labels': [str(l) for l in colors_per_genome.index],
            'values': [int(v) for v in colors_per_genome.values],
            'mean': float(colors_per_genome.mean())
        },
        'plot4_distances': {
            'values': [float(d) for d in distances]
        }
    }
    
    # Deriving directory from output_file path
    json_path = os.path.join(temp_path, 'summary.json')
    
    with open(json_path, 'w') as f:
        json.dump(stats_data, f, indent=4)

def main():
    parser = setup_parser()
    args = parser.parse_args()
    
    # Validate input files
    validate_files(args.cluster_file, args.tsv_file)
    
    # Process cluster data and give clusters >2 a color
    color_map, gene_name_map = process_clusters(args.cluster_file, args.gene_of_interest, args.size_for_cluster)

    # Process neigbours and gene_of_interest data
    df = pd.read_csv(args.tsv_file, sep='\t', header=None,
                     names=['genome', 'gene', 'start', 'end', 
                            'orientation', 'type'], dtype=str)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)

    df['unique_id'] = df.apply(lambda row: '_'.join(row.astype(str)), axis=1)
    df['cluster_color'] = df['unique_id'].map(color_map).fillna('#FFFFFF')
    df.loc[df['gene'] == args.gene_of_interest, 'cluster_color'] = '#FF0000'
    df['gene_name_based_on_DB'] = df['unique_id'].map(gene_name_map)
    df['gene_name_based_on_DB'] = df['gene_name_based_on_DB'].fillna(df['gene'])
    df.drop('unique_id', axis=1, inplace=True)

    df.to_csv(args.output_file, sep='\t', index=False, header=False)
    
    summary(args.gene_of_interest, args.goi_type, df, args.output_file)

if __name__ == '__main__':
    main()
