import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sugar import Feature, FeatureList
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

def setup_parser_func():
    """Set up and return the argument parser."""
    parser = argparse.ArgumentParser(description='Process genomic data and generate plots.')
    parser.add_argument('--scale', required=True, choices=['yes', 'no'], help='Nucleotide Scale')
    parser.add_argument('--input_file', required=True, help='Path to input file')
    parser.add_argument('--output_path', required=True, help='Path to final png or svg')
    parser.add_argument('--output_ending', required=True, choices=['png', 'svg'], help='Final plot format')
    parser.add_argument('--cluster', type=int, default=2, help='Minimal size for a cluster')
    parser.add_argument('--threshold', type=float, default=0.3, help='Similarity threshold for clustering')
    parser.add_argument('--gene_of_interest', required=True, help='gene_of_interest') 
    parser.add_argument('--name_file', required=False, help='Renaming contigs to genome names')
    parser.add_argument('--gene_lable', required=False, choices=['yes', 'no'], help='Add lables to gene [no as default]') #cut_height_args
    parser.add_argument('--cut_height_args', type=float, default=0.5, help='Cut threshold for dendogram clustering.')
    return parser

def map_ori_func(value):
    """Map orientation values."""
    return {
        'antisense': '-',
        'sense': '+'
    }.get(value, value)

def process_group_func(group, gene_of_interest):
    """Process a group of rows, adjusting orientations and positions."""
    center_row = group[group['gene'] == gene_of_interest]
    
    if not center_row.empty and center_row.iloc[0]['orientation'] == '-':
        group.loc[group['gene'] != gene_of_interest, 'orientation'] = group.loc[group['gene'] != gene_of_interest, 'orientation'].map({'+': '-', '-': '+'})
        group.loc[group['gene'] == gene_of_interest, 'orientation'] = '+'
        group = group.iloc[::-1].reset_index(drop=True)
        total_length = group['adjusted_end'].max() - group['adjusted_start'].min()
        group['adjusted_start'], group['adjusted_end'] = total_length - group['adjusted_end'] + group['adjusted_start'].min(), total_length - group['adjusted_start'] + group['adjusted_start'].min()

    return group

def prepare_dataframe_func(input_file, genome_name_tsv=""):
    """Read and prepare the input dataframe."""
    header = ['genome', 'gene', 'start', 'end', 'orientation', 'type', 'cluster_color']
    df = pd.read_csv(input_file, sep='\t', header=None, names=header)
    df['length'] = abs(df['start'] - df['end'])
    df['contig'] = df['genome'].str.split(':', n=1).str[0]
    df['orientation'] = df['orientation'].apply(map_ori_func)
    
    if genome_name_tsv:
        rename_df = pd.read_csv(genome_name_tsv, sep='\t')
        merged = df.merge(rename_df[['contig', 'organism_name']], 
                        on='contig', 
                        how='left')
        df = merged
    else:
        df['organism_name'] = ""
    
    return df

def cluster_genomes_func(df, output_path, output_ending, cut_height=0.5, cluster=2, threshold=0.3):
    """Cluster genomes based on feature matrix."""

    unique_colors = df['cluster_color'].unique()
    unique_genomes = df["genome"].unique()
    feature_matrix_color = pd.DataFrame(0, index=unique_genomes, columns=unique_colors)
    
    for _, row in df.iterrows():
        feature_matrix_color.at[row["genome"], row["cluster_color"]] = 1

    X = feature_matrix_color.values

    # Compute pairwise Jaccard distances
    distance_matrix = pdist(X, metric='jaccard')

    # Perform hierarchical clustering
    Z = linkage(distance_matrix, method='ward')    
    complete_clustering = linkage(distance_matrix, method="complete")
    average_clustering = linkage(distance_matrix, method="average")
    single_clustering = linkage(distance_matrix, method="single")

    dendrogram(complete_clustering)
    plt.tight_layout()
    plt.savefig(f"{output_path}.Tree.complete_clustering.{output_ending}")
    plt.close()
    dendrogram(average_clustering)
    plt.tight_layout()
    plt.savefig(f"{output_path}.Tree.average_clustering.{output_ending}")
    plt.close()    
    dendrogram(single_clustering)
    plt.tight_layout()
    plt.savefig(f"{output_path}.Tree.single_clustering.{output_ending}")
    plt.close()

    # Create a mapping from genome to organism_name and Get the organism names in the same order as feature_matrix_color.index
    genome_to_organism = df.set_index("genome")["organism_name"].to_dict()
    dendro_labels = [genome_to_organism.get(genome, genome) for genome in feature_matrix_color.index]
    
    # Plot dendrogram with a horizontal line at your chosen height
    plt.figure(figsize=(10, 6), dpi=750)
    dendrogram(Z, color_threshold=cut_height, labels=dendro_labels, leaf_rotation=90)
    plt.title('Clustering Dendrogram')
    plt.xlabel('Genomes')
    plt.ylabel('Jaccard Distance')

    # Choose a height (distance) to cut the tree
    plt.axhline(y=cut_height, color='r', linestyle='--', label=f'Cut at {cut_height}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_path}.Tree.{output_ending}")
    plt.close()

    # Assign clusters based on the chosen height
    l = fcluster(Z, t=cut_height, criterion='distance')
    
    # Assign -1 to genomes with only one color
    cluster_labels = pd.Series(l, index=feature_matrix_color.index)
    for genome in feature_matrix_color.index:
        if feature_matrix_color.loc[genome].sum() == 1:
            cluster_labels[genome] = -1
        if feature_matrix_color.loc[genome].sum() == 0:
            cluster_labels[genome] = -2
    
    return cluster_labels.values, feature_matrix_color

def plot_cluster_func(cluster_df, output_path, output_ending, scale, gene_of_interest, gene_lable, max_subplots=20):
    """Generate plots for a cluster."""
    genome_groups = list(cluster_df.groupby("genome"))
    num_subplots = int(np.ceil(len(genome_groups) / max_subplots))

    for subplot_idx in range(num_subplots):
        merged = {} # Merge all color dictionaries
        subplot_genomes = genome_groups[
            subplot_idx * max_subplots: (subplot_idx + 1) * max_subplots
        ]
        
        if not subplot_genomes:  # Added validation
            continue
            
        fig, axes = plt.subplots(len(subplot_genomes), 1, figsize=(10, len(subplot_genomes) * 0.9))
        axes = [axes] if len(subplot_genomes) == 1 else axes

        for i, (genome_name, genome_df) in enumerate(subplot_genomes):
            dicts_subplot = plot_genome_func(genome_df, axes[i], scale, gene_of_interest, gene_lable)
            merged.update(dicts_subplot)
            unique_organisms_name = genome_df['organism_name'].unique().tolist()[0]
            axes[i].set_title(unique_organisms_name, fontstyle='italic', loc='left')
        
        if gene_lable=='no': # Add legend to the subplot
            legend_handles = [
                Patch(facecolor=color, label=label)
                for label, color in merged.items()
                if isinstance(color, tuple) and len(color) == 3 and all(0 <= c <= 1 for c in color)
            ]
            fig.legend(handles=legend_handles, title="Legend", bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout()
        count = f"{cluster_df['cluster_label'].iloc[0]}.{subplot_idx}"
        fig.savefig(f"{output_path}_{count}.{output_ending}", bbox_inches='tight')

        plt.close(fig)


def plot_genome_func(genome_df, ax, scale, gene_of_interest, gene_lable='no'):
    """Plot a single genome."""
    df = genome_df.sort_values(by=['start'])
    df['organism_name'] = genome_df['organism_name']
    min_start = df['start'].min()
    df['start'] = df['start'].astype(int)
    df['length'] = df['length'].astype(int)
    min_start = int(df['start'].min())

    if scale.lower() == 'yes':
        df['adjusted_start'] = (df['start'] - min_start).astype(int)
        df['adjusted_end'] = (df['adjusted_start'] + df['length']).astype(int)
        seqlen_distance = (df['adjusted_end'].max() - df['adjusted_start'].min()).astype(int)
    else:
        max_distance = 0
        df['adjusted_start'] = (df['start'] - min_start) / 1000
        df['adjusted_end'] = (df['end'] - min_start) / 1000
        min_start = genome_df['start'].min()
        genome_locus_length = genome_df['end'].max() - min_start
        max_distance = max(max_distance, genome_locus_length)
        seqlen_distance = max_distance / 1000

    df = process_group_func(df, gene_of_interest)

    features = [Feature(row.type, start=row.adjusted_start, stop=row.adjusted_end, strand=row.orientation,
                        meta={'seqid': row.contig, 'organism_name': row.organism_name, 'name': row.gene, 'cluster_color': row.cluster_color})
                for _, row in df.iterrows()]

    color_dic = dict(zip(df['gene'], df['cluster_color']))

    if gene_lable=='no':
        FeatureList(features).plot_ftsviewer(ax=ax, label=None, 
                                        colorby='name', color=color_dic,
                                        seqlen=seqlen_distance, figsize=(7, 5),
                                        with_ruler=False, show=False)
    else:
        FeatureList(features).plot_ftsviewer(ax=ax, label='name', 
                                        colorby='name', color=color_dic,
                                        seqlen=seqlen_distance, figsize=(7, 5),
                                        with_ruler=False, show=False, 
                                        labels_spacing=60, fontdict={'fontsize': 7})
    
    return color_dic

def cosine_similarity(matrix):
    # Calculate norms of vectors & Scalar product matrix & External product norms
    norms = np.linalg.norm(matrix, axis=1)
    dot_products = np.dot(matrix, matrix.T)
    norm_matrix = np.outer(norms, norms)

    # Elementwise division -> Cosine similarity matrix
    similarity_matrix = np.divide(dot_products, norm_matrix, where=norm_matrix != 0)
    similarity_matrix[norm_matrix == 0] = 0.0

    return similarity_matrix

def write_labels(df, gene_of_interest, output_path):
    center_row = df[df['gene'] == gene_of_interest]
    new_df = center_row[['cluster_label', 'organism_name', 'contig', 'gene', 'start', 'end', 'orientation']]
    write_header = not os.path.exists(f'{output_path}_gene_of_interest.tsv')     # Check if the file exists to decide whether to write the header
    new_df.to_csv(f'{output_path}_gene_of_interest.tsv', mode='a', header=write_header, sep='\t', index=False)

def main():
        parser = setup_parser_func()
        args = parser.parse_args()

        #try:
        df = prepare_dataframe_func(args.input_file, args.name_file)

        labels, feature_matrix_color = cluster_genomes_func(df, args.output_path, args.output_ending, args.cut_height_args, args.cluster, args.threshold)
        
        # correlation of cluster labels to each genome
        genome_to_cluster = pd.Series(labels, index=feature_matrix_color.index).groupby(level=0).first()

        df["cluster_label"] = df["genome"].map(genome_to_cluster)

        for cluster_label, cluster_df in df.groupby("cluster_label"):
            write_labels(cluster_df, args.gene_of_interest, args.output_path)
            plot_cluster_func(cluster_df, args.output_path, args.output_ending, args.scale, args.gene_of_interest, args.gene_lable)

        cluster_stats = []

        for cluster_label, cluster_df in df.groupby("cluster_label"):
            cluster_size = cluster_df['genome'].nunique()
            color_count = cluster_df['cluster_color'].nunique()

            if cluster_label == -1:
                cluster_stats.append({
                    'cluster_label': cluster_label,
                    'genomes': cluster_size})
                continue  # Skip Noise Cluster

            # Compute cosine Similarity for Clusters
            cluster_genomes_for_cosine = feature_matrix_color.loc[genome_to_cluster[genome_to_cluster == cluster_label].index]

            if cluster_genomes_for_cosine.shape[0] > 1:
                similarity_matrix = cosine_similarity(cluster_genomes_for_cosine.values)
                avg_similarity = np.mean(similarity_matrix[np.triu_indices_from(similarity_matrix, k=1)])
            else:
                avg_similarity = 1.0

            #calculation of common genes
            cluster_genome_vectors = feature_matrix_color.loc[genome_to_cluster[genome_to_cluster == cluster_label].index]

            # extract common genes correlating to value 1 in the color/genome vector
            common_colors = np.all(cluster_genome_vectors.values == 1, axis=0)
            common_colors_count = np.sum(common_colors)  # count of common genes

            cluster_stats.append({
                'cluster_label': cluster_label,
                'genomes': cluster_size,
                'unique_genes': color_count,
                'common_genes': common_colors_count,
                'avg_cosine_similarity': round(avg_similarity, 4)
            })

        # output as csv
        stats_df = pd.DataFrame(cluster_stats)
        stats_df.to_csv(f"{args.output_path}_cluster_stats.csv", index=False)

    #except Exception as e:
    #    print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()
