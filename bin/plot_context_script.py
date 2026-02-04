# (C) 2024, Maria Schreiber, MIT license
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sugar import Feature, FeatureList
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, cut_tree
from scipy.stats import pearsonr
from sklearn.cluster import DBSCAN
import sys
sys.setrecursionlimit(100000)
from collections import Counter
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from sklearn.decomposition import PCA, TruncatedSVD
import seaborn as sns

def setup_parser_func():
    '''Set up and return the argument parser.'''
    parser = argparse.ArgumentParser(description='Process genomic data and generate plots.')
    parser.add_argument('--scale', required=True, choices=['yes', 'no'], help='Nucleotide Scale')
    parser.add_argument('--input_file', required=True, help='Path to input file')
    parser.add_argument('--output_path', required=True, help='Path to final png or svg')
    parser.add_argument('--gene_of_interest', required=True, help='Name of gene_of_interest')
    parser.add_argument('--gene_name', required=False, default='db', help='Use db-name [db] or gene-ids [id] for plotting.')
    parser.add_argument('--name_file', required=False, help='Renaming contigs to genome names.')
    parser.add_argument('--gene_lable', required=False, choices=['yes', 'no'], default='yes', help='Add lables to gene [yes as default]')
    parser.add_argument('--contig_cluster_file', required=False, help='Filter for Contigs your are interested in.')
    parser.add_argument('--cut_height_args', type=float, default=0, help='Cut threshold for dendogram clustering.')
    parser.add_argument('--cluster', type=int, default=2, help='Minimal size for a cluster')
    parser.add_argument('--threshold', type=float, default=0.9, help='Similarity threshold for clustering')
    parser.add_argument('--microsynteny_logo', required=False, choices=['yes', 'no'], default='no', help='Filter for Contigs your are interested in.')
    parser.add_argument('--method', required=False, default='all', 
                        choices=['dbscan', 'ward_cosinus', 'ward_euclidean', 'all'],
                        help='Choose clustering method to run (default: all)')
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
        total_length = group['end'].max() - group['start'].min()
        group['start'], group['end'] = total_length - group['end'] + group['start'].min(), total_length - group['start'] + group['start'].min()

    return group

def prepare_dataframe_func(input_file, gene_of_interest, gene_name, genome_name_tsv=''):
    '''Read and prepare the input dataframe.'''
    header = ['genome', 'gene_id', 'start', 'end', 'orientation', 'type', 'cluster_color', 'gene_db']
    df = pd.read_csv(input_file, sep='\t', header=None, names=header)

    df['length'] = abs(df['start'] - df['end'])
    df['contig'] = df['genome'].str.split(':', n=1).str[0]
    df['orientation'] = df['orientation'].apply(map_ori_func)
    df['original_start'] = df['start']
    df['original_end'] = df['end']

    # Use the DB name for logic/clustering (this keeps colors stable)
    df['gene'] = df['gene_db']

    # Create the display column explicitly
    if gene_name == 'id':
        df['display_name'] = df['gene_id']
    else:
        df['display_name'] = df['gene_db']
    
    # Ensure it's a string to avoid downstream issues
    df['display_name'] = df['display_name'].astype(str)

    if genome_name_tsv and genome_name_tsv.lower() != 'null' and os.path.exists(genome_name_tsv):
        rename_df = pd.read_csv(genome_name_tsv, sep='\t')
        merged = df.merge(rename_df[['contig', 'organism_name']],
                        on='contig',
                        how='left')
        merged['organism_name'] = merged['organism_name'].astype(str) + ' (' + merged['contig'].astype(str) + ')'
        df = merged
        organism_name_from_file = 'yes'
    else:
        df['organism_name'] = ''
        organism_name_from_file = 'no'

    groups = df.groupby('genome')
    dfs = []

    for group_name, group_df in groups: # make gene of interest sense
        regroup = process_group_func(group_df, gene_of_interest)
        dfs.append(regroup)  # add each result DataFrame to the list

    # Concatenate all the DataFrames into one
    big_df = pd.concat(dfs, ignore_index=True)
    return big_df, organism_name_from_file

def cluster_genomes_func(df, output_path, cut_height_para, dbscan_size, dbscan_threshold, selected_method):
    '''Cluster genomes based on feature matrix.'''

    unique_colors = df['cluster_color'].unique()
    unique_genomes = df['genome'].unique()
    feature_matrix_color = pd.DataFrame(0, index=unique_genomes, columns=unique_colors)

    for _, row in df.iterrows():
        feature_matrix_color.at[row['genome'], row['cluster_color']] = 1

    # Create a mapping from genome to organism_name and get the organism names in the same order as feature_matrix_color.index
    genome_to_organism = df.set_index('genome')['organism_name'].to_dict()

    X = feature_matrix_color.values.copy()
    # 1 Pre-assign cluster labels for simple cases AND dimension reduction
    preassigned = {}
    clustering_mask = []

    for genome in feature_matrix_color.index:
        row_sum = feature_matrix_color.loc[genome].sum()
        if row_sum <= 2: 
            preassigned[genome] = -1 if row_sum == 2 else -2
            clustering_mask.append(False)
        else:
            clustering_mask.append(True)
    
    clustering_mask = np.array(clustering_mask)
    X_cluster = X[clustering_mask]
    
    if X_cluster.shape[0] == 0:
        print("Warning: No genomes met the criteria for clustering (all had <= 2 genes). Skipping PCA and clustering.")
        # Return empty results or placeholders to prevent the script from crashing
        return np.array([]), pd.Series([]), pd.Series([]), feature_matrix_color
    

    def test_suitability_func(X, output_path, method='cosine'):
        '''Evaluates how "clusterable" the data is using Similarity instead of Pearson. Works for both PCA and TruncatedSVD.'''
        if X.size == 0 or X.shape[0] < 2:
            print("Warning: Not enough data for suitability test.")
            return

        # Calculate Similarity Matrix
        # Cosine is robust for sparse genomic data (presence/absence)
        if method == 'cosine':
            sim_matrix = cosine_similarity(X)
            label_str = 'Cosine Similarity'
            cmap = 'YlGnBu'
        if method == 'euclidean':
            sim_matrix = euclidean_distances(X)
            label_str = 'Euclidean Similarity'
            cmap = 'YlGnBu'

        # 2. Clustered Heatmap
        g = sns.clustermap(sim_matrix, cmap='YlGnBu', method='average', figsize=(10, 10))
        g.fig.suptitle(f'Clustered Genomic {label_str}', y=1.02)
        
        os.makedirs(output_path, exist_ok=True)
        g.savefig(f'{output_path}/{label_str}_korrelations_heatmap.png')
        plt.close()

        # 3. Distribution Plot (Without Silhouette for now)
        upper_tri = sim_matrix[np.triu_indices_from(sim_matrix, k=1)]
        plt.figure(figsize=(8, 6))
        plt.hist(upper_tri, bins=30, color='skyblue', edgecolor='black')
        plt.title(f'Distribution of Pairwise Similarities (Mean: {np.mean(upper_tri):.3f})')
        plt.xlabel(label_str)
        plt.ylabel('Frequency')
        plt.savefig(f'{output_path}/{label_str}_distribution.png')

        plt.close()

    if selected_method in ['ward_cosinus', 'all']:
        test_suitability_func(X_cluster, output_path, 'cosine')
    elif selected_method in ['ward_euclidean', 'dbscan', 'all']:
        test_suitability_func(X_cluster, output_path, 'euclidean')

    def apply_adaptive_dim_reduction(X_cluster, n_components=100, sparsity_threshold=0.9):
        '''
        Checks matrix sparsity and applies either PCA or TruncatedSVD.
        Sparsity = (Number of Zeros) / (Total Elements)
        '''
        # Calculate sparsity
        n_total = X_cluster.size
        n_zeros = np.count_nonzero(X_cluster == 0)
        sparsity = n_zeros / n_total
        
        # Ensure n_components isn't larger than our data dimensions
        n_samples, n_features = X_cluster.shape
        actual_components = min(n_components, n_samples - 1, n_features - 1)

        print(f"--- Dimensionality Reduction Report ---")
        print(f"Matrix Sparsity: {sparsity:.2%}")

        if sparsity > sparsity_threshold:
            print(f"Decision: High sparsity detected. Using TruncatedSVD.")
            # TruncatedSVD is better for sparse data as it doesn't center the mean
            model = TruncatedSVD(n_components=actual_components, random_state=42)
            reduced_data = model.fit_transform(X_cluster)
            explained_var = model.explained_variance_ratio_.sum()
        else:
            print(f"Decision: Dense data detected. Using PCA.")
            # PCA is standard for dense data
            model = PCA(n_components=actual_components, random_state=42)
            reduced_data = model.fit_transform(X_cluster)
            explained_var = model.explained_variance_ratio_.sum()

        print(f"Total Explained Variance: {explained_var:.2%}")
        print(f"---------------------------------------")
        
        return reduced_data, model

    # Dimension reduction
    reduced_data, reduction_model = apply_adaptive_dim_reduction(X_cluster)

    # Fill preassigned
    full_clusters = np.empty(len(feature_matrix_color), dtype=int)
    full_clusters[:] = -99  # default placeholder
    preassigned_indices = np.where(~clustering_mask)[0]
    for idx in preassigned_indices:
        genome = feature_matrix_color.index[idx]
        full_clusters[idx] = preassigned[genome]

    clustered_indices = np.where(clustering_mask)
    """
    Mustakim, E. Rahmi, M. R. Mundzir, S. T. Rizaldi, Okfalisa and I. Maita, "Comparison of DBSCAN and PCA-DBSCAN Algorithm for Grouping Earthquake Area," 2021 International Congress of Advanced Technology and Engineering (ICOTEN), Taiz, Yemen, 2021, pp. 1-5, doi: 10.1109/ICOTEN52080.2021.9493497. keywords: {Dimensionality reduction;Geology;Earthquakes;Clustering algorithms;Indexes;Principal component analysis;Meteorology;Climatology and Geophysics Agency;
    DBSCAN;DBSCAN-PCA;Earthquake Area;PCA},
    """
    def pca_vis(full_clusters, clustered_indices, clusters, output_path, folder, title, reduction_model, reduced_data):

        df = pd.DataFrame(reduced_data)
        sns.pairplot(df)
        plt.savefig(f'{output_path}/{folder}/PCA.pairplot.{folder}.png')
        full_clusters[clustered_indices] = clusters
        n_components_found = reduced_data.shape[1]

        plt.figure(figsize=(8, 6))
        # --- Change: Visualization only uses first two components for 2D plot ---
        plt.scatter(reduced_data[:, 0], reduced_data[:, 1], c=clusters, cmap='tab20')

        unique_clusters, counts = np.unique(clusters, return_counts=True)
        cmap = plt.get_cmap('tab20', len(unique_clusters))
        legend_handles = [
            Patch(color=cmap(i), label=f'Cluster {cluster_label} ({count} points)')
            for i, (cluster_label, count) in enumerate(zip(unique_clusters, counts))
        ]
        plt.legend(handles=legend_handles, title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left')

        # Get variance from the reduction_model (PCA or SVD)
        var1 = reduction_model.explained_variance_ratio_[0] * 100
        var2 = reduction_model.explained_variance_ratio_[1] * 100

        label_type = "PC" if isinstance(reduction_model, PCA) else "SVD"
        plt.xlabel(f'{label_type} 1 ({var1:.2f}% variance)')
        plt.ylabel(f'{label_type} 2 ({var2:.2f}% variance)')
        plt.title(f"{title}\n(Showing top 2 of {n_components_found} PCs)")
        plt.tight_layout()
        os.makedirs(f'{output_path}/{folder}', exist_ok=True)
        plt.savefig(f'{output_path}/{folder}/{label_type}.{folder}.png')
        plt.savefig(f'{output_path}/{folder}/{label_type}.{folder}.svg')
        plt.close()

    # 2 Clustering
    def density_clustering(eps, min_samples, reduced_data, pca_clusters, folder, reduction_model, output_path):
        #dbscan = HDBSCAN(min_cluster_size=5, min_samples=1, cluster_selection_epsilon=0.02
        #cluster_selection_method='eom' # Ensures stability for big clusters) #eps=eps, 
        # NOTE: with more PC followes a Curse of Dimensionality.
        # 1. Rule of Thumb: A common starting point is MinPts=2xdimensions
        #       As you add PCA components to reach 75% variance, the average Euclidean distance between points naturally increases.
        #       BUT large dataset need more than just 2x...
        #       BUT if you want small cluster, stay below 3
        # Calculate min_samples based on the dimensionality of the reduced data
        # Common rule: 2 * dimensions
        auto_min_samples = 2 * reduced_data.shape[1]
        
        dbscan = DBSCAN(eps=eps, min_samples=auto_min_samples)
        dbscan_clusters = dbscan.fit_predict(reduced_data)

        clustered_indices = np.where(clustering_mask)
        pca_clusters[clustered_indices] = dbscan_clusters

        title = f'{folder} (eps={eps}, min_samples={auto_min_samples}, n={len(X_cluster)})'
        
        # Call pca_vis with all required arguments
        pca_vis(pca_clusters, clustered_indices, dbscan_clusters, output_path, folder, title, reduction_model, reduced_data)

        dbscan_labels_series = pd.Series(pca_clusters, index=feature_matrix_color.index)
        return dbscan_labels_series.values.copy()

    def hierarchical_clustering(reduced_data, full_clusters, reduction_model, metric, method, folder, output_path, cut_height_para=0.25):
        import matplotlib.colors as mcolors
        import matplotlib.patches as mpatches
        
        # Nested helper for dendrogram coloring
        def dendogram_vis(linkage_result, all_labels, title, cut_height, clustered_labels_series, cluster_to_color):
            def link_color_func(node_id):
                if node_id < len(X_cluster): # Uses X_cluster from outer scope
                    genome_id = feature_matrix_color.index[clustering_mask][node_id]
                    cluster_id = clustered_labels_series[genome_id]
                    return cluster_to_color.get(cluster_id, "grey")
                else:
                    dist = linkage_result[node_id - len(X_cluster), 2]
                    return link_color_func(int(linkage_result[node_id - len(X_cluster), 0])) if dist < cut_height else "grey"
                        
            plt.figure(figsize=(10, 5))
            dendrogram(linkage_result, color_threshold=None, link_color_func=link_color_func, leaf_rotation=90, no_labels=True)
            
            # Legend Logic
            cluster_counts = Counter(all_labels)
            legend_handles = [mpatches.Patch(color=color, label=f"Cluster {cl} ({cluster_counts.get(cl, 0)})") 
                            for cl, color in cluster_to_color.items()]
            
            plt.legend(handles=legend_handles, title="Clusters", bbox_to_anchor=(1.05, 1), loc="upper left")
            plt.title(title)
            plt.tight_layout(rect=[0, 0, 0.85, 1])
            os.makedirs(f'{output_path}/{folder}', exist_ok=True)
            plt.savefig(f'{output_path}/{folder}/Tree.{folder}.png')
            plt.close()

        # --- Clustering Execution ---
        dist_matrix = pdist(reduced_data, metric=metric)
        linkage_result = linkage(dist_matrix, method=method)
        cut_height = cut_height_para * max(linkage_result[:, 2])
        clustered_labels = fcluster(linkage_result, t=cut_height, criterion="distance")

        # Map labels to the index
        clustered_labels_series = pd.Series(clustered_labels, index=feature_matrix_color.index[clustering_mask])

        # Re-integrate with preassigned noise/small clusters
        all_labels = pd.Series(index=feature_matrix_color.index, dtype=int)
        for genome in feature_matrix_color.index:
            all_labels[genome] = preassigned[genome] if genome in preassigned else clustered_labels_series[genome]
        
        # Assign Colors
        unique_clusters = sorted(set(clustered_labels))
        cmap = plt.get_cmap("tab20", len(unique_clusters))
        cluster_to_color = {cl: mcolors.to_hex(cmap(i)) for i, cl in enumerate(unique_clusters)}

        # Visualizations
        title = f'{folder} (n={len(X_cluster)})'
        dendogram_vis(linkage_result, all_labels, title, cut_height, clustered_labels_series, cluster_to_color)
        
        # Corrected pca_vis call
        pca_vis(full_clusters, clustered_indices, clustered_labels, output_path, folder, title, reduction_model, reduced_data)

        return all_labels

    results = {}

    if selected_method in ['dbscan', 'all']:
        print("Running Ward Cosine...")
        results['dbscan'] = density_clustering(
            dbscan_threshold, 
            dbscan_size, 
            reduced_data, 
            full_clusters.copy(), 
            'dbscan', 
            reduction_model,   # Model passed here
            output_path        # Path passed here
        )

    if selected_method in ['ward_cosinus', 'all']:
        print("Running Ward Cosine...")
        results['ward_cosinus'] = hierarchical_clustering(
            reduced_data, 
            full_clusters.copy(), 
            reduction_model,   # Positioned correctly here
            metric="cosine", 
            method="ward", 
            folder="ward_cosinus",
            output_path=output_path, 
            cut_height_para=cut_height_para
        )

    if selected_method in ['ward_euclidean', 'all']:
        print("Running Ward Euclidean...")
        results['ward_euclidean'] = hierarchical_clustering(
            reduced_data, 
            full_clusters.copy(), 
            reduction_model,   # Positioned correctly here
            metric="euclidean", 
            method="ward", 
            folder="ward_euclidean",
            output_path=output_path, 
            cut_height_para=cut_height_para
        )

    return results, feature_matrix_color

def plot_cluster_func(cluster_df, output_path, scale, gene_of_interest, gene_lable, folder, organism_name_from_file, max_subplots=15):

    '''Generate plots for a cluster.'''
    # Hatches are working now. The best way to use them is by setting ft.meta._ftsviewer_hatch attribute.
    # hatch_linewidth can be set in the same way starting with matplotlib 3.10
    # color can be set in the same way or via colorby and color arguments of plot_ftsviewer
    # needs sugar v1.0
    name2color = {gene: color for gene, color in zip(cluster_df['gene'], cluster_df['cluster_color'])}
    
    all_fts = FeatureList([
        Feature(row.type, start=row.start, stop=row.end, strand=row.orientation,
                meta={
                    'seqid': row.contig, 
                    'hit': row.genome, 
                    'organism_name': row.organism_name, 
                    'logic_name': row.gene,        # Use for alignment/color
                    'label_name': getattr(row, 'display_name', row.gene) # Fallback to gene if display_name fails
                })
                for row in cluster_df.itertuples(index=False)
        ])

    scale_bool = scale.lower() == 'yes'
    label_bool = gene_lable.lower() == 'yes'
    groups = list(all_fts.groupby('hit').values())

    for fignum in range(0, 1 + (len(groups)-1)//max_subplots):
        fts = FeatureList()
        for g in groups[fignum * max_subplots:(fignum + 1) * max_subplots]:
            fts.extend(g)
        align_features = FeatureList([f for f in fts if f.meta.get('logic_name') == gene_of_interest]) if scale_bool else None
        #align = fts.select(name=gene_of_interest) if scale else None
        fig = fts.plot_ftsviewer(
            groupby='hit',
            axlabel='organism_name' if organism_name_from_file.lower() == 'yes' else 'seqid',
            # Display the unique ID
            label=lambda ft: ft.meta['label_name'] if label_bool else None,
            align=align_features,
            crop=500,
            sharex=scale_bool,
            colorby=lambda ft: ft.meta['logic_name'],
            color=name2color,
            figsize=(10, 1.5*len(fts.groupby('hit'))), ncols=1,
            with_ruler=False,
            labels_spacing=8, fontdict={'fontsize': 7}
        )

        if not label_bool: # Add legend to the subplot
            # Track which (color, name) combinations we've already seen
            seen_combinations = set()
            legend_handles = []

            for label, color in name2color.items():
                # Create the unique combination as a tuple
                combo = (color, label)

                # Check if color is valid RGB and this exact (color+name) is new
                if (isinstance(color, tuple) and
                    len(color) == 3 and
                    all(0 <= c <= 1 for c in color) and
                    combo not in seen_combinations):

                    # Mark this combination as seen
                    seen_combinations.add(combo)

                    # Add to legend
                    patch = Patch(facecolor=color, label=label)
                    legend_handles.append(patch)

            # Add legend outside the plot
            fig.legend(handles=legend_handles, title='Legend', bbox_to_anchor=(1.05, 1), loc='upper left')

        count = f"{cluster_df['cluster_label'].iloc[0]}"
        fig.savefig(f'{output_path}/{folder}/synteny_cluster_{count}_{fignum:03d}.png', bbox_inches='tight')
        fig.savefig(f'{output_path}/{folder}/synteny_cluster_{count}_{fignum:03d}.svg', bbox_inches='tight')
        plt.close(fig)

def cosine_similarity(matrix):
    ''' Calculate Cosine similarity matrix safely with float initialization. '''
    # Convert input to float to ensure all downstream math uses floats
    matrix = matrix.astype(float) 
    
    norms = np.linalg.norm(matrix, axis=1)
    dot_products = np.dot(matrix, matrix.T)
    norm_matrix = np.outer(norms, norms)
    
    # Force the output array to be float64 to accept division results
    similarity_matrix = np.zeros_like(dot_products, dtype=float)
    
    mask = norm_matrix != 0
    np.divide(dot_products, norm_matrix, out=similarity_matrix, where=mask)
    
    return similarity_matrix

def write_labels(df, gene_of_interest, output_path, folder):
    ''' Writing some stats '''
    center_row = df[df['gene'] == gene_of_interest]
    new_df = center_row[['cluster_label', 'organism_name', 'contig', 'gene', 'original_start', 'original_end', 'orientation']]

    write_header = not os.path.exists(f'{output_path}/{folder}/gene_of_interest.tsv')     # Check if the file exists to decide whether to write the header
    new_df.to_csv(f'{output_path}/{folder}/gene_of_interest.tsv', mode='a', header=write_header, sep='\t', index=False)

def for_aln(cluster_label, clustering_algo, df, output_dir):
    filename = output_dir + "/" + clustering_algo + "/cluster_" + str(cluster_label) + ".tsv"
    colors_by_genome = df.groupby('genome')['cluster_color'].apply(lambda x: ''.join(x)).reset_index()
    colors_by_genome.to_csv(filename, sep='\t', index=False)

def main():
        parser = setup_parser_func()
        args = parser.parse_args()

        #try:
        print("Read data")
        df, organism_name_from_file = prepare_dataframe_func(args.input_file, args.gene_of_interest, args.gene_name, args.name_file)

        print("Process for clustering")
        clustering_results, feature_matrix_color = cluster_genomes_func(
            df, args.output_path, args.cut_height_args, args.cluster, args.threshold, args.method
        )

        # Now iterate only over the methods that were actually run
        for method_name, labels in clustering_results.items():
            print(f"Generating Plots and Stats for: {method_name}")
            # Reset cluster_label for each method
            genome_to_cluster = pd.Series(labels, index=feature_matrix_color.index).groupby(level=0).first()
            df['cluster_label'] = df['genome'].map(genome_to_cluster)

            unique_contigs = []  # Default: no filter
            if args.contig_cluster_file:
                filter_df = pd.read_csv(args.contig_cluster_file, header=None, names=["Contigs"])
                unique_contigs = filter_df['Contigs'].unique().tolist()

            cluster_stats = []
            for cluster_label, cluster_df in df.groupby('cluster_label'):
                if args.microsynteny_logo == 'yes':
                    for_aln(cluster_label, method_name, cluster_df, args.output_path)

                # If the list is not empty, check if any contigs in the current cluster_df are in the unique_contigs list
                if unique_contigs and not cluster_df['contig'].isin(unique_contigs).any():
                    continue # Skip to the next cluster if no contigs from the list are found

                write_labels(cluster_df, args.gene_of_interest, args.output_path, method_name)
                plot_cluster_func(cluster_df, args.output_path, args.scale, args.gene_of_interest, args.gene_lable, method_name, organism_name_from_file)
                cluster_size = cluster_df['genome'].nunique()
                color_count = cluster_df['cluster_color'].nunique()

                if cluster_label == -1:
                    cluster_stats.append({
                        'cluster_label': cluster_label,
                        'genomes': cluster_size})
                    continue  # Skip Noise Cluster

                #---------------------------------------------
                # Compute cosine Similarity for Clusters
                cluster_genomes_for_cosine = feature_matrix_color.loc[genome_to_cluster[genome_to_cluster == cluster_label].index]

                if cluster_genomes_for_cosine.shape[0] > 1:
                    similarity_matrix = cosine_similarity(cluster_genomes_for_cosine.values)
                    avg_similarity = np.mean(similarity_matrix[np.triu_indices_from(similarity_matrix, k=1)])
                else:
                    avg_similarity = 1.0

                # calculation of common genes
                cluster_genome_vectors = feature_matrix_color.loc[genome_to_cluster[genome_to_cluster == cluster_label].index]

                if cluster_genome_vectors.shape[0] > 1:
                    # 1. Cosine Similarity
                    sim_matrix = cosine_similarity(cluster_genome_vectors.values)
                    avg_cosine = np.mean(sim_matrix[np.triu_indices_from(sim_matrix, k=1)])
                    
                    # 2. Euclidean Distance
                    # pdist returns pairwise distances in a condensed flat array
                    euk_distances = pdist(cluster_genome_vectors.values, metric='euclidean')
                    avg_euclidean = np.mean(euk_distances)
                else:
                    avg_cosine = 1.0
                    avg_euclidean = 0.0

                # Extract common genes
                common_colors = np.all(cluster_genome_vectors.values == 1, axis=0)
                common_colors_count = np.sum(common_colors)

                cluster_stats.append({
                    'cluster_label': cluster_label,
                    'genomes': cluster_size,
                    'unique_genes': color_count,
                    'common_genes': common_colors_count,
                    'avg_cosine_similarity': round(avg_cosine, 4),
                    'avg_euclidean_dist': round(avg_euclidean, 4)
                })

            # output as csv
            stats_df = pd.DataFrame(cluster_stats)
            stats_df.to_csv(f'{args.output_path}{method_name}/cluster_stats.csv', index=False)

    #except Exception as e:
    #    print(f'An error occurred: {str(e)}')

if __name__ == '__main__':
    main()
