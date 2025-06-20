import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sugar import Feature, FeatureList
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.cluster import KMeans, DBSCAN, HDBSCAN, Birch
from sklearn.decomposition import PCA
import sys
sys.setrecursionlimit(100000)

def setup_parser_func():
    '''Set up and return the argument parser.'''
    parser = argparse.ArgumentParser(description='Process genomic data and generate plots.')
    parser.add_argument('--scale', required=True, choices=['yes', 'no'], help='Nucleotide Scale')
    parser.add_argument('--input_file', required=True, help='Path to input file')
    parser.add_argument('--output_path', required=True, help='Path to final png or svg')
    parser.add_argument('--output_ending', required=True, choices=['png', 'svg'], help='Final plot format')
    parser.add_argument('--cluster', type=int, default=2, help='Minimal size for a cluster')
    parser.add_argument('--threshold', type=float, default=0.9, help='Similarity threshold for clustering')
    parser.add_argument('--gene_of_interest', required=True, help='gene_of_interest') 
    parser.add_argument('--name_file', required=False, help='Renaming contigs to genome names')
    parser.add_argument('--gene_lable', required=False, choices=['yes', 'no'], help='Add lables to gene [no as default]') 
    parser.add_argument('--cut_height_args', type=float, default=0, help='Cut threshold for dendogram clustering.')
    return parser

def map_ori_func(value):
    '''Map orientation values.'''
    return {
        'antisense': '-',
        'sense': '+'
    }.get(value, value)

def process_group_func(group, gene_of_interest):
    '''Process a group of rows, adjusting orientations and positions.'''
    center_row = group[group['gene'] == gene_of_interest]
    
    if not center_row.empty and center_row.iloc[0]['orientation'] == '-':
        group.loc[group['gene'] != gene_of_interest, 'orientation'] = group.loc[group['gene'] != gene_of_interest, 'orientation'].map({'+': '-', '-': '+'})
        group.loc[group['gene'] == gene_of_interest, 'orientation'] = '+'
        group = group.iloc[::-1].reset_index(drop=True)
        total_length = group['adjusted_end'].max() - group['adjusted_start'].min()
        group['adjusted_start'], group['adjusted_end'] = total_length - group['adjusted_end'] + group['adjusted_start'].min(), total_length - group['adjusted_start'] + group['adjusted_start'].min()

    return group

def prepare_dataframe_func(input_file, genome_name_tsv=''):
    '''Read and prepare the input dataframe.'''
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
        df['organism_name'] = ''
    
    return df

def cluster_genomes_func(df, output_path, output_ending, cut_height=0, cluster=2, threshold=0.1):
    '''Cluster genomes based on feature matrix.'''

    unique_colors = df['cluster_color'].unique()
    unique_genomes = df['genome'].unique()
    feature_matrix_color = pd.DataFrame(0, index=unique_genomes, columns=unique_colors)
    
    for _, row in df.iterrows():
        feature_matrix_color.at[row['genome'], row['cluster_color']] = 1

    # Create a mapping from genome to organism_name and get the organism names in the same order as feature_matrix_color.index
    genome_to_organism = df.set_index('genome')['organism_name'].to_dict()
    dendro_labels = [genome_to_organism.get(genome, genome) for genome in feature_matrix_color.index]
    
    X = feature_matrix_color.values

    # Compute pairwise Jaccard distances
    distance_matrix = pdist(X, metric='cosine') # jaccard

    def create_dendogram(distance_matrix, method, name, cut_height=0):
        # Perform hierarchical clustering
        clustering = linkage(distance_matrix, method=method)
        if cut_height==0:
            cut_height = 0.5*max(clustering[:,2])
        dendrogram(clustering, color_threshold=cut_height, labels=dendro_labels, leaf_rotation=90)
        plt.title(name)
        plt.xlabel('Genomes')
        plt.ylabel('Cosine Distance')
        plt.tight_layout()
        plt.savefig(f'{output_path}/{method}/Tree.{name}.{output_ending}')
        plt.close()
        return clustering, cut_height

    ward_clustering, ward_cut_height = create_dendogram(distance_matrix, 'ward', 'ward_clustering') # Ward variance minimization algorithm
    complete_clustering, complete_cut_height = create_dendogram(distance_matrix, 'complete', 'complete_clustering') # Farthest Point Algorithm or Voor Hees Algorithm 
    average_clustering, average_cut_height = create_dendogram(distance_matrix, 'average', 'average_clustering') # UPGMA algorithm
    single_clustering, single_cut_height = create_dendogram(distance_matrix, 'single', 'single_clustering') # Nearest Point Algorithm.
    centroid_clustering, centroid_cut_height = create_dendogram(distance_matrix, 'centroid', 'centroid_clustering') # UPGMC
    weighted_clustering, weighted_cut_height = create_dendogram(distance_matrix, 'weighted', 'weighted_clustering') # WPGMA

    def cluster_label(clustering, cut_height):
        # Assign clusters based on the chosen height
        fcluster_result = fcluster(clustering, t=cut_height, criterion='distance')
        
        # Assign -1 to genomes with only one color
        cluster_labels = pd.Series(fcluster_result, index=feature_matrix_color.index)
        for genome in feature_matrix_color.index:
            if feature_matrix_color.loc[genome].sum() == 1:
                cluster_labels[genome] = -1
            if feature_matrix_color.loc[genome].sum() == 0:
                cluster_labels[genome] = -2
        
        return cluster_labels.values

    ward = cluster_label(ward_clustering, ward_cut_height)
    complete = cluster_label(complete_clustering, complete_cut_height)
    average = cluster_label(average_clustering, average_cut_height)
    single = cluster_label(single_clustering, single_cut_height)
    centroid = cluster_label(centroid_clustering, centroid_cut_height)
    weighted = cluster_label(weighted_clustering, weighted_cut_height)

    # --- KMeans clustering ---
    print('KMeans')
    kmeans_cluster = KMeans(n_clusters=4, random_state=0)
    kmeans_labels = kmeans_cluster.fit_predict(X)
    centroids = kmeans_cluster.cluster_centers_  # Optional: for plotting centroids
    visualize_kmeans(X, kmeans_labels, output_path, centroids, title="KMeans Clustering of Genomes")
    # Optionally, assign -1 or -2 for single/no color genomes as above
    kmeans_labels_series = pd.Series(kmeans_labels, index=feature_matrix_color.index)
    for genome in feature_matrix_color.index:
        if feature_matrix_color.loc[genome].sum() == 1:
            kmeans_labels_series[genome] = -1
        if feature_matrix_color.loc[genome].sum() == 0:
            kmeans_labels_series[genome] = -2
    kmeans = kmeans_labels_series.values
    # -------------------------

    # --- DBSCAN clustering ---
    print('dbscan')
    dbscan_cluster = DBSCAN(eps=0.1, min_samples=2)  # Adjust eps and min_samples as needed
    dbscan_labels = dbscan_cluster.fit_predict(X)
    dbscan_labels_reduced, _ = visualize_pca_umap(X, dbscan_labels, output_path, 'dbscan')
    # Optionally, assign -1 or -2 for single/no color genomes as above
    dbscan_labels_series = pd.Series(dbscan_labels_reduced, index=feature_matrix_color.index)
    for genome in feature_matrix_color.index:
        if feature_matrix_color.loc[genome].sum() == 1:
            dbscan_labels_series[genome] = -1  # Or another special value
        if feature_matrix_color.loc[genome].sum() == 0:
            dbscan_labels_series[genome] = -2  # Or another special value
    dbscan = dbscan_labels_series.values
    # -------------------------

    # --- HDBSCAN clustering ---
    print('hdbscan')
    hdbscan_cluster = HDBSCAN(min_cluster_size=2)
    hdbscan_labels = hdbscan_cluster.fit_predict(X)
    hdbscan_labels_reduced, _ = visualize_pca_umap(X, hdbscan_labels, output_path, 'hdbscan')
    # Optionally, assign -1 or -2 for single/no color genomes as above
    hdbscan_labels_series = pd.Series(hdbscan_labels_reduced, index=feature_matrix_color.index)
    for genome in feature_matrix_color.index:
        if feature_matrix_color.loc[genome].sum() == 1:
            hdbscan_labels_series[genome] = -1  # Or another special value
        if feature_matrix_color.loc[genome].sum() == 0:
            hdbscan_labels_series[genome] = -2  # Or another special value
    hdbscan = hdbscan_labels_series.values
    # -------------------------

    # --- Birch clustering --- (Balanced Iterative Reducing and Clustering using Hierarchies) 
    print('birch')
    birch_cluster = Birch(n_clusters=None)
    birch_labels = birch_cluster.fit_predict(X)
    birch_labels_reduced, cluster_labels = visualize_pca_umap(X, birch_labels, output_path, 'birch')
    # Optionally, assign -1 or -2 for single/no color genomes as above
    birch_labels_series = pd.Series(cluster_labels, index=feature_matrix_color.index)
    for genome in feature_matrix_color.index:
        if feature_matrix_color.loc[genome].sum() == 1:
            birch_labels_series[genome] = -1  # Or another special value
        if feature_matrix_color.loc[genome].sum() == 0:
            birch_labels_series[genome] = -2  # Or another special value
    birch = birch_labels_series.values
    # -------------------------

    return ward, complete, average, single, centroid, weighted, kmeans, dbscan, hdbscan, cluster_labels, feature_matrix_color

def visualize_pca_umap(X, labels, output_path, folder):
    # Reduce data to 2D for plotting # Cluster the Reduced Data

    pca = PCA(n_components=2)
    reduced_PCA = pca.fit_transform(X)

    from sklearn.cluster import AgglomerativeClustering
    clustering = AgglomerativeClustering(n_clusters=2, metric='euclidean', linkage='ward')
    cluster_labels = clustering.fit_predict(reduced_PCA)

    import umap.umap_ as umap
    reducer = umap.UMAP()
    X_umap = reducer.fit_transform(X) 

    if folder == 'dbscan':
        dbscan = DBSCAN(eps=0.1, min_samples=2) 
        clusters = dbscan.fit_predict(reduced_PCA)
    
    if folder == 'hdbscan':
        hdbscan = HDBSCAN(min_cluster_size=2)
        clusters = hdbscan.fit_predict(reduced_PCA)

    if folder == 'birch':
        birch = Birch(n_clusters=None)
        cluster_labels = birch.fit_predict(reduced_PCA)
        clusters = birch.fit_predict(reduced_PCA)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(reduced_PCA[:, 0], reduced_PCA[:, 1], c=labels, cmap='viridis', alpha=0.7)
    plt.colorbar(scatter, label='Cluster Label')
    plt.title(folder)
    plt.xlabel('PCA 1')
    plt.ylabel('PCA 2')
    plt.grid()
    plt.savefig(f'{output_path}/{folder}/PCA.{folder}.png')
    plt.close()

    plt.scatter(reduced_PCA[:, 0], reduced_PCA[:, 1], c=clusters, cmap='viridis')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('Clustering on PCA-reduced data')
    plt.savefig(f'{output_path}/{folder}/PCA.REDUCED.{folder}.png')
    plt.close()

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(X_umap[:, 0], X_umap[:, 1], c=labels, cmap='viridis', alpha=0.7)
    plt.colorbar(scatter, label='Cluster Label')
    plt.title('UMAP Projection')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.grid()
    plt.savefig(f'{output_path}/{folder}/UMAP.{folder}.png')
    plt.close()

    return clusters, cluster_labels

def visualize_kmeans(X, labels, output_path, centroids=None, title="KMeans Clusters"):
    # Reduce data to 2D for plotting
    pca = PCA(n_components=2)
    reduced_PCA = pca.fit_transform(X)
    
    plt.figure(figsize=(8, 6))
    # Plot data points colored by cluster
    scatter = plt.scatter(reduced_PCA[:, 0], reduced_PCA[:, 1], c=labels, cmap='viridis', alpha=0.7, label='Genomes')
    
    # Plot centroids if available (optional)
    if centroids is not None:
        centroids_reduced = pca.transform(centroids)
        plt.scatter(centroids_reduced[:, 0], centroids_reduced[:, 1],
                    c='red', marker='X', s=200, label='Centroids')
    
    plt.title(title)
    plt.xlabel('PCA 1')
    plt.ylabel('PCA 2')
    plt.colorbar(scatter, label='Cluster Label')
    plt.legend()
    plt.grid()
    plt.savefig(f'{output_path}/kmeans/PCA.kmeans.png')
    plt.close()

def plot_cluster_func(cluster_df, output_path, output_ending, scale, gene_of_interest, gene_lable, folder, max_subplots=20):
    '''Generate plots for a cluster.'''
    genome_groups = list(cluster_df.groupby('genome'))
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
            fig.legend(handles=legend_handles, title='Legend', bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout()
        count = f'{cluster_df['cluster_label'].iloc[0]}.{subplot_idx}'
        fig.savefig(f'{output_path}/{folder}/{count}.{output_ending}', bbox_inches='tight')

        plt.close(fig)

def plot_genome_func(genome_df, ax, scale, gene_of_interest, gene_lable='no'):
    '''Plot a single genome.'''
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
    ''' Calculate norms of vectors & Scalar product matrix & External product norms and elementwise division -> Cosine similarity matrix '''

    norms = np.linalg.norm(matrix, axis=1)
    dot_products = np.dot(matrix, matrix.T)
    norm_matrix = np.outer(norms, norms)

    similarity_matrix = np.divide(dot_products, norm_matrix, where=norm_matrix != 0)
    similarity_matrix[norm_matrix == 0] = 0.0

    return similarity_matrix

def write_labels(df, gene_of_interest, output_path, folder):
    ''' Writing some stats '''
    center_row = df[df['gene'] == gene_of_interest]
    new_df = center_row[['cluster_label', 'organism_name', 'contig', 'gene', 'start', 'end', 'orientation']]
    write_header = not os.path.exists(f'{output_path}/{folder}/gene_of_interest.tsv')     # Check if the file exists to decide whether to write the header
    new_df.to_csv(f'{output_path}/{folder}/gene_of_interest.tsv', mode='a', header=write_header, sep='\t', index=False)

def main():
        parser = setup_parser_func()
        args = parser.parse_args()

        #try:
        df = prepare_dataframe_func(args.input_file, args.name_file)
        ward_labels, complete_labels, average_labels, single_labels, centroid_labels, weighted_lables, kmeans_lables, dbscan_lables, hdbscan_labels, birch_labels, feature_matrix_color = cluster_genomes_func(df, args.output_path, args.output_ending, args.cut_height_args, args.cluster, args.threshold)

        lst = ['ward', 'complete', 'average', 'single', 'centroid', 'weighted', 'kmeans', 'dbscan', 'hdbscan', 'birch']
        i = 0

        for labels in [ward_labels, complete_labels, average_labels, single_labels, centroid_labels, weighted_lables, kmeans_lables, dbscan_lables, hdbscan_labels, birch_labels]:
            # correlation of cluster labels to each genome
            genome_to_cluster = pd.Series(labels, index=feature_matrix_color.index).groupby(level=0).first()

            df['cluster_label'] = df['genome'].map(genome_to_cluster)

            for cluster_label, cluster_df in df.groupby('cluster_label'):
                write_labels(cluster_df, args.gene_of_interest, args.output_path, lst[i])
                plot_cluster_func(cluster_df, args.output_path, args.output_ending, args.scale, args.gene_of_interest, args.gene_lable, lst[i])

            cluster_stats = []

            for cluster_label, cluster_df in df.groupby('cluster_label'):
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
            stats_df.to_csv(f'{args.output_path}{lst[i]}/cluster_stats.csv', index=False)
            i = i + 1

    #except Exception as e:
    #    print(f'An error occurred: {str(e)}')

if __name__ == '__main__':
    main()
