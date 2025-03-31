import argparse
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from sugar import Feature, FeatureList

def setup_parser_func():
    """Set up and return the argument parser."""
    parser = argparse.ArgumentParser(description='Process genomic data and generate plots.')
    parser.add_argument('--scale', required=True, choices=['yes', 'no'], help='Nucleotide Scale')
    parser.add_argument('--input_file', required=True, help='Path to input file')
    parser.add_argument('--output_path', required=True, help='Path to final png or svg')
    parser.add_argument('--output_ending', required=True, choices=['png', 'svg'], help='Final plot format')
    parser.add_argument('--cluster', type=int, default=2, help='Minimal size for a cluster')
    parser.add_argument('--threshold', type=float, default=0.3, help='Similarity threshold for clustering')
    return parser

def map_ori_func(value):
    """Map orientation values."""
    return {
        'antisense': '-',
        'sense': '+'
    }.get(value, value)

def process_group_func(group):
    """Process a group of rows, adjusting orientations and positions."""
    infernal_row = group[group['gene'] == 'sibling']
    if not infernal_row.empty and infernal_row.iloc[0]['orientation'] == '-':
        group.loc[group['gene'] != 'sibling', 'orientation'] = group.loc[group['gene'] != 'sibling', 'orientation'].map({'+': '-', '-': '+'})
        group.loc[group['gene'] == 'sibling', 'orientation'] = '+'
        group = group.iloc[::-1].reset_index(drop=True)
        total_length = group['adjusted_end'].max() - group['adjusted_start'].min()
        group['adjusted_start'], group['adjusted_end'] = total_length - group['adjusted_end'] + group['adjusted_start'].min(), total_length - group['adjusted_start'] + group['adjusted_start'].min()
    return group

def prepare_dataframe_func(input_file):
    """Read and prepare the input dataframe."""
    header = ['genome', 'gene', 'start', 'end', 'orientation', 'type', 'cluster_color']
    df = pd.read_csv(input_file, sep='\t', header=None, names=header)
    df['length'] = abs(df['start'] - df['end'])
    df['genome_name'] = df['genome'].str.split(':', n=1).str[0]
    df['orientation'] = df['orientation'].apply(map_ori_func)
    return df

def cluster_genomes_func(df, cluster=2, threshold=0.3):
    """Cluster genomes based on feature matrix."""
    unique_colors = df['cluster_color'].unique()
    unique_genomes = df["genome"].unique()
    feature_matrix_color = pd.DataFrame(0, index=unique_genomes, columns=unique_colors)
    
    for _, row in df.iterrows():
        feature_matrix_color.at[row["genome"], row["cluster_color"]] = 1

    db = DBSCAN(eps=threshold, min_samples=cluster, metric='cosine').fit(feature_matrix_color)

    return db.labels_, feature_matrix_color

def plot_cluster_func(cluster_df, output_path, output_ending, scale, max_subplots=20):
    """Generate plots for a cluster."""
    genome_groups = list(cluster_df.groupby("genome"))
    num_subplots = int(np.ceil(len(genome_groups) / max_subplots))

    for subplot_idx in range(num_subplots):
        subplot_genomes = genome_groups[subplot_idx * max_subplots: (subplot_idx + 1) * max_subplots]
        max_distance = compute_max_distance_func(subplot_genomes, scale)
        
        fig, axes = plt.subplots(len(subplot_genomes), 1, figsize=(10, len(subplot_genomes) * 0.9))
        axes = [axes] if len(subplot_genomes) == 1 else axes

        for i, (genome_name, genome_df) in enumerate(subplot_genomes):
            plot_genome_func(genome_df, axes[i], scale, max_distance)

        plt.tight_layout()
        count = f"{cluster_df['cluster_label'].iloc[0]+1}.{subplot_idx + 1}"
        plt.savefig(f"{output_path}_{count}.{output_ending}")
        plt.close(fig)

def compute_max_distance_func(genome_groups, scale):
    """Compute the maximum distance for scaling."""
    if scale.lower() == 'no':
        return max(genome_df['end'].max() - genome_df['start'].min() for _, genome_df in genome_groups) / 1000
    return 0

def plot_genome_func(genome_df, ax, scale, max_distance):
    """Plot a single genome."""
    df = genome_df.sort_values(by=['start'])
    min_start = df['start'].min()
    
    if scale.lower() == 'yes':
        df['adjusted_start'] = df['start'] - min_start
        df['adjusted_end'] = df['adjusted_start'] + df['length']
        init_distance = df['adjusted_end'].max()
    else:
        df['adjusted_start'] = (df['start'] - min_start) / 1000
        df['adjusted_end'] = (df['end'] - min_start) / 1000
        init_distance = max_distance

    df = process_group_func(df)

    features = [Feature(row.type, start=row.adjusted_start, stop=row.adjusted_end, strand=row.orientation,
                        meta={'seqid': row.genome_name, 'name': row.gene, 'cluster_color': row.cluster_color})
                for _, row in df.iterrows()]
    color_dic = dict(zip(df['gene'], df['cluster_color']))

    # TODO if only white gene symboles, than do not plot it
    FeatureList(features).plot_ftsviewer(ax=ax, colorby='name', color=color_dic,
                                         seqlen=init_distance, figsize=(7, 5),
                                         with_ruler=False, show=False)

def cosine_similarity(matrix):
    # Calculate norms of vectors
    norms = np.linalg.norm(matrix, axis=1)

    # Scalar product matrix
    dot_products = np.dot(matrix, matrix.T)

    # External product norms
    norm_matrix = np.outer(norms, norms)

    # Elementwise division -> Cosine similarity matrix
    similarity_matrix = np.divide(dot_products, norm_matrix, where=norm_matrix != 0)

    # avoid division by 0
    similarity_matrix[norm_matrix == 0] = 0.0

    return similarity_matrix


def main():
    parser = setup_parser_func()
    args = parser.parse_args()

    try:
        df = prepare_dataframe_func(args.input_file)
        labels, feature_matrix_color = cluster_genomes_func(df, args.cluster, args.threshold)
        # correlation of cluster labels to each genome
        genome_to_cluster = pd.Series(labels, index=feature_matrix_color.index).groupby(level=0).first()
        df["cluster_label"] = df["genome"].map(genome_to_cluster)

        # old df["cluster_label"] = df["genome"].map(pd.Series(labels, index=feature_matrix.index))

        for cluster_label, cluster_df in df.groupby("cluster_label"):
            plot_cluster_func(cluster_df, args.output_path, args.output_ending, args.scale)

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

    except Exception as e:
        print(f"An error occurred: {str(e)}")
   
if __name__ == "__main__":
    main()
