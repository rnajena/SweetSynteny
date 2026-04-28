# (C) 2024, Maria Schreiber, MIT license
import os
import argparse
import json
import logging
import pandas as pd
import numpy as np

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import matplotlib.patches as mpatches


import umap
import seaborn as sns
from sugar import Feature, FeatureList
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import HDBSCAN
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import silhouette_score
from sklearn.feature_extraction.text import TfidfTransformer
import fastcluster

CORE_GENE_THRESHOLD = 0.99      # Present in >99% of genomes
RARE_GENE_MIN_COUNT = 2         # Must appear in at least 2 genomes
DEFAULT_SIMILARITY_THRESHOLD = 0.5
MAX_LEGEND_ITEMS = 40

def setup_parser_func():
    '''Set up and return the argument parser.'''
    parser = argparse.ArgumentParser(description='Process genomic data and generate plots.')
    parser.add_argument('--scale', required=False, default='yes', choices=['yes', 'no'], help='Nucleotide Scale')
    parser.add_argument('--input_file', required=True, help='Path to input file')
    parser.add_argument('--output_path', required=True, help='Path to final png or svg')
    parser.add_argument('--svg', required=False, default='no', help='Creating svg for large data takes a lot of time.')
    parser.add_argument('--gene_of_interest', required=True, help='Name of gene_of_interest')
    parser.add_argument('--goi_type', required=True, help='Type of gene_of_interest for oriantation')
    parser.add_argument('--gene_name', required=False, default='db', help='Use db-name [db] or gene-ids [id] for plotting.')
    parser.add_argument('--name_file', required=False, help='Renaming contigs to genome names.')
    parser.add_argument('--gene_lable', required=False, choices=['yes', 'no'], default='yes', help='Add lables to gene [yes as default]')
    parser.add_argument('--contig_cluster_file', required=False, help='Filter for Contigs your are interested in.')
    parser.add_argument('--cut_height_args', type=float, default=0, help='Cut threshold for dendogram clustering.')
    parser.add_argument('--cluster', type=int, default=2, help='Minimal size for a cluster')
    parser.add_argument('--threshold', type=float, default=0.9, help='Similarity threshold for clustering')
    parser.add_argument('--microsynteny_logo', required=False, choices=['yes', 'no'], default='no', help='Filter for Contigs your are interested in.')
    parser.add_argument('--metric', required=False, default='cosine', 
                        choices=['cosine', 'euclidean'],
                        help='Choose clustering metric to run (default: cosine)')
    return parser

def map_ori_func(value):
    '''Map orientation values.'''
    return {
        'antisense': '-',
        'sense': '+'
    }.get(value, value)

def process_group_func(group, goi_type):
    '''Process a group of rows, adjusting orientations and positions.'''
    center_row = group[group['type'] == goi_type]

    if not center_row.empty and center_row.iloc[0]['orientation'] == '-':
        group.loc[group['type'] != goi_type, 'orientation'] = group.loc[group['type'] != goi_type, 'orientation'].map({'+': '-', '-': '+'})
        group.loc[group['type'] == goi_type, 'orientation'] = '+'
        group = group.iloc[:-1].reset_index(drop=True)
        total_length = group['end'].max() - group['start'].min()
        group['start'], group['end'] = total_length - group['end'] + group['start'].min(), total_length - group['start'] + group['start'].min()

    return group

def prepare_dataframe_func(input_file, goi_type, gene_name, genome_name_tsv=''):
    '''Read and prepare the input dataframe.'''
    if not os.path.exists(input_file):        
        raise FileNotFoundError(f"Input file not found: {input_file}")

    header = ['genome', 'gene_id', 'start', 'end', 'orientation', 'type', 'cluster_color', 'gene_db']
    df = pd.read_csv(input_file, sep='\t', header=None, names=header)

    if df.empty:
        raise ValueError("Input file is empty")

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

    for _, group_df in groups: # make gene of interest sense
        regroup = process_group_func(group_df, goi_type)
        dfs.append(regroup)  # add each result DataFrame to the list

    # Concatenate all the DataFrames into one
    big_df = pd.concat(dfs, ignore_index=True)
    return big_df, organism_name_from_file

def cluster_genomes_func(df, output_path, metric='cosine'):
    ''' Step by step:
    > Reduces dimensions via PCA.
    > Performs Hierarchical and Density clustering.
    > Generates Similarity Distribution, Clustermap, and PCA visualizations.
    '''

    def umap_vis(clusters, reduced_data, reduction_model, output_path, folder, title, variance_model=None, svg='no'):
        
        os.makedirs(os.path.join(output_path, folder), exist_ok=True)
        
        if hasattr(reduced_data, 'toarray'):
            data_to_plot = reduced_data.toarray()
        else:
            data_to_plot = reduced_data
    
        n_comp = data_to_plot.shape[1]
        
        # Determine Prefix and Axis Labels
        if 'UMAP' in str(type(reduction_model)):
            prefix = 'UMAP'
            cols = [f'UMAP{i+1}' for i in range(min(n_comp, 6))]
        else:
            prefix = 'PC' if 'PCA' in str(type(reduction_model)) else 'SVD'
            var = reduction_model.explained_variance_ratio_ * 100
            cols = [f'{prefix}{i+1} ({var[i]:.1f}%)' for i in range(min(n_comp, 6))]

        df_plot = pd.DataFrame(data_to_plot[:, :len(cols)], columns=cols)
        df_plot['Cluster'] = clusters.astype(str)

        # 2. PAIRPLOT (Multi-dimensional view)
        try:
            sample_size = min(len(df_plot), 2000)
            g = sns.pairplot(df_plot.sample(sample_size, random_state=42), 
                            hue='Cluster', palette='tab20', corner=True, 
                            plot_kws={'alpha': 0.5, 's': 20})
            g.fig.suptitle(f'{title} (Top {len(cols)} Comp)', y=1.02)
            save_plot(g, os.path.join(output_path, folder, f'{prefix}.pairplot'), svg)
        except Exception as e:
            logger.warning(f'Pairplot failed: {e}')

        # 3. 2D SCATTER PLOT (Primary view)
        plt.figure(figsize=(10, 6))
        ax = sns.scatterplot(data=df_plot, x=cols[0], y=cols[1], hue='Cluster', 
                            palette='tab20', alpha=0.7, s=50, edgecolor='w')
        
        # Clean up legend (limit to top 40 items to avoid overflow)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:MAX_LEGEND_ITEMS], labels[:MAX_LEGEND_ITEMS], title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2 if len(labels) > 20 else 1)
        
        plt.title(f'{title}\nTop 2 Components')
        plt.tight_layout()
        save_plot(plt, os.path.join(output_path, folder, f'{prefix}_2D'), svg)
        plt.close()

        if variance_model is not None and hasattr(variance_model, 'explained_variance_ratio_'):
            plot_variance(variance_model, output_path, folder)

    def save_plot(plt_obj, base_path, svg):
        '''Helper to save plot in multiple formats.'''
        plt_obj.savefig(f'{base_path}.png', dpi=150, bbox_inches='tight')
        if svg == 'yes':
            plt_obj.savefig(f'{base_path}.svg', bbox_inches='tight')
    
    def plot_variance(reducer, output_path, label):
        plt.figure(figsize=(8, 4))
        exp_var = reducer.explained_variance_ratio_
        plt.bar(range(1, len(exp_var)+1), exp_var, alpha=0.5, align='center', label='Individual')
        plt.step(range(1, len(exp_var)+1), np.cumsum(exp_var), where='mid', label='Cumulative')
        plt.axhline(y=0.75, color='r', linestyle='--', label='75% Threshold')
        plt.ylabel('Explained Variance Ratio')
        plt.xlabel('Principal Components')
        plt.legend(loc='best')
        plt.title(f'Variance Distribution')
        plt.savefig(f'{output_path}/{label}_variance_screen_plot.png')
        plt.savefig(f'{output_path}/{label}_variance_screen_plot.svg')
        plt.close()

    def optimize_hdbscan(X, min_cluster_range=range(2, 10), metric='jaccard'):
        ''' Sweeps through HDBSCAN parameters using a precomputed distance matrix. '''
        if X.shape[0] < 3:        
            # Not enough samples to cluster meaningfully        
            return np.zeros(X.shape[0], dtype=int)    

        # Also cap min_cluster_range based on data size    
        max_min_size = max(2, X.shape[0] // 3)    
        min_cluster_range = range(2, min(10, max_min_size))

        # 1. Precompute the distance matrix. 
        dist_matrix = pairwise_distances(X, metric=metric)
        dist_matrix = (dist_matrix + dist_matrix.T) / 2
        # Ensure the diagonal is exactly zero
        np.fill_diagonal(dist_matrix, 0)
        
        num_samples = X.shape[0]
        best_score = -1.0
        best_labels = np.zeros(num_samples) - 1 # Default to all noise
        
        for min_size in min_cluster_range:
            # Use 'precomputed' because we already did the hard work
            hdb = HDBSCAN(min_cluster_size=min_size, min_samples=1, metric='precomputed', copy=True)
            labels = hdb.fit_predict(dist_matrix)
            
            core_mask = labels != -1
            unique_labels = np.unique(labels[core_mask])
            n_clusters = len(unique_labels)
            
            if n_clusters > 1:
                # silhouette_score also needs to know we are using a distance matrix
                # We slice the distance matrix to include only non-noise points
                sub_dist = dist_matrix[core_mask][:, core_mask]
                score = silhouette_score(sub_dist, labels[core_mask], metric='precomputed')
                
                noise_ratio = 1.0 - (core_mask.sum() / num_samples)
                adjusted_score = score * (1.0 - noise_ratio) 
                
                if adjusted_score > best_score:
                    best_score = adjusted_score
                    best_labels = labels
                    best_min_size = min_size
                        
        if best_score > -1.0:
            logger.info(f'Best HDBSCAN min_cluster_size: {best_min_size} (Score: {best_score:.2f})')
            return best_labels
        else:
            # Fallback: simple fit on the precomputed distances
            return HDBSCAN(min_cluster_size=2, metric='precomputed', copy=True).fit_predict(dist_matrix)

    def get_optimal_hierarchical_clusters(X, linkage_matrix, method, metric, min_avg_sim=DEFAULT_SIMILARITY_THRESHOLD):
        # Step 1: Initial cut with stricter threshold to avoid overly broad clusters.
        # Use 'distance' criterion at t=0.5.
        # This forces tighter groupings that truly share conserved gene patterns.
        initial_labels = fcluster(linkage_matrix, t=0.5, criterion='distance')
        refined_labels = np.array(initial_labels).flatten().astype(int)
        
        # Safety Check for the IndexError you saw earlier
        if X.shape[0] != len(refined_labels):
            # Fallback: slice X to match labels if there is a mismatch
            X = X[:len(refined_labels), :]

        unique_labels = np.unique(refined_labels)
        max_label = refined_labels.max()

        for cluster_id in unique_labels:
            indices = np.where(refined_labels == cluster_id)[0]
            if len(indices) <= 1: continue
            
            cluster_data = X.iloc[indices, :] if hasattr(X, 'iloc') else X[indices, :]
            dense_cluster = cluster_data.toarray() if hasattr(cluster_data, 'toarray') else cluster_data
            
            sims = cosine_similarity(dense_cluster)
            avg_sim = sims[np.triu_indices(sims.shape[0], k=1)].mean()
            
            # If the cluster is too loose (not like your figures), split it.
            # If it's tight (like your figures), keep it!
            if avg_sim < min_avg_sim:
                logger.info(f'Refining Cluster {cluster_id} (Sim: {avg_sim:.2f})')
                # LOCAL complete linkage is the best at finding conserved blocks
                local_link = linkage(dense_cluster, method=method, metric=metric)
                local_labels = fcluster(local_link, t=0.5, criterion='distance')
                
                if len(np.unique(local_labels)) > 1:
                    for i, l_id in enumerate(local_labels):
                        refined_labels[indices[i]] = max_label + l_id
                    max_label = refined_labels.max()
                    
        return refined_labels
    
    def merge_singleton_clusters(X, labels, metric='cosine', similarity_threshold=0.3, min_cluster_size=2):
        """Merge singletons and small clusters into real clusters; enforce minimum cluster size."""
        refined_labels = labels.copy()
        unique_labels, counts = np.unique(labels, return_counts=True)
        
        # Find small clusters (singletons or below min_cluster_size)
        small_mask = np.isin(labels, unique_labels[counts <= min_cluster_size])
        real_clusters = unique_labels[counts > min_cluster_size]
        
        if len(real_clusters) == 0:
            logger.warning('No clusters meet minimum size; all will be kept but may be small.')
            return refined_labels
        
        if not any(small_mask):
            logger.info('No small clusters to merge.')
            return refined_labels

        # 1. Vectorized Centroid Calculation from real clusters
        centroids = np.array([
            X[labels == c].mean(axis=0).A1 if hasattr(X, 'toarray') else X[labels == c].mean(axis=0)
            for c in real_clusters
        ])
        
        # 2. Extract all small cluster vectors at once
        small_vectors = X[small_mask].toarray() if hasattr(X, 'toarray') else X[small_mask]

        # 3. Batch Similarity Calculation
        if metric == 'cosine':
            norm_s = small_vectors / (np.linalg.norm(small_vectors, axis=1, keepdims=True) + 1e-9)
            norm_c = centroids / (np.linalg.norm(centroids, axis=1, keepdims=True) + 1e-9)
            sims = np.dot(norm_s, norm_c.T)
        else:
            dist = euclidean_distances(small_vectors, centroids)
            sims = 1 / (1 + dist)

        # 4. Global Update: always merge small clusters into their nearest real cluster
        best_cluster_indices = np.argmax(sims, axis=1)
        
        # Force assignment to the best matching real cluster (no orphans allowed)
        refined_labels[small_mask] = real_clusters[best_cluster_indices]
        logger.info(f'Merged {np.sum(small_mask)} small cluster members into real clusters.')
                    
        return refined_labels

    os.makedirs(output_path, exist_ok=True)
    
    # --- 1. PREPARE FEATURE MATRIX ---
    unique_colors = df['cluster_color'].unique()
    unique_genomes = df['genome'].unique()
    feature_matrix = pd.DataFrame(0, index=unique_genomes, columns=unique_colors)
    feature_matrix = df.pivot_table(index='genome',     
                                    columns='cluster_color',     
                                    aggfunc='size',     
                                    fill_value=0).clip(upper=1)

    # 1. Remove 'Core' features (present in > 99% of genomes)
    core_mask = feature_matrix.mean(axis=0) < CORE_GENE_THRESHOLD

    # 2. Remove 'Ultra-rare' noise (present in only 1 or 2 genomes)
    rare_mask = feature_matrix.sum(axis=0) > RARE_GENE_MIN_COUNT
    X_cluster = feature_matrix.loc[:, core_mask & rare_mask]

    if X_cluster.shape[1] == 0:
        raise ValueError("No gene clusters remain after applying core/rare filtering. "
                         "Check your input data and threshold settings.")

    # Remove rows (genomes) that have no features after filtering
    X_cluster = X_cluster.loc[X_cluster.sum(axis=1) > 0]
    if X_cluster.shape[0] == 0:
        raise ValueError("No genomes remain after filtering for non-empty feature profiles.")

    # --- 2. TF-IDF Transformer ---
    tfidf = TfidfTransformer()
    X_weighted = tfidf.fit_transform(X_cluster)

    # --- 3. DIMENSIONALITY REDUCTION (PCA vs SVD vs UMAP) ---
    # Initial fit to determine variance distribution
    actual_components = min(X_cluster.shape[0], X_cluster.shape[1], 100)
    if actual_components < 2:
        raise ValueError("Insufficient genomes or gene clusters for dimensionality reduction.")
    '''It operates directly on sparse matrices without dense transformation. It will compress your highly 
    correlated gene clusters into a smaller set of 'topics' or 'metagenes.'
    '''
    svd_reducer = TruncatedSVD(n_components=max(1, actual_components - 1))
    X_svd_reduced = svd_reducer.fit_transform(X_weighted)
    '''Use this after Truncated SVD. UMAP is fantastic for learning the non-linear manifold of your data 
    and is highly effective at preparing data for density-based clustering.
    '''
    # For very small datasets, skip UMAP and use SVD directly
    if X_svd_reduced.shape[0] < 5 or X_svd_reduced.shape[1] < 2:
        X_reduced = X_svd_reduced
        reduction_label = 'svd'
        reduction_model = svd_reducer
    else:
        n_neighbors = min(15, max(2, X_svd_reduced.shape[0] - 1))
        n_components = min(actual_components, X_svd_reduced.shape[0] - 1, X_svd_reduced.shape[1])
        n_jobs = min(4, os.cpu_count() or 1)

        umap_reducer = umap.UMAP(
                n_neighbors=n_neighbors,
                min_dist=0.1,
                metric=metric,
                n_components=n_components,
                n_jobs=n_jobs,
                random_state=42
        )
        X_reduced = umap_reducer.fit_transform(X_svd_reduced)
        reduction_label = 'svd+UMAP'
        reduction_model = umap_reducer
    # while reducing dimensionality.
    '''TF-IDF vectors are continuous and magnitude-sensitive; Cosine distance measures the angle between these vectors, 
    which perfectly captures the similarity of their gene profiles regardless of the absolute number of genes in the location.
    '''
    if metric == 'cosine':
        # pdist expects a dense array; SVD output is already dense
        # 'cosine' distance in scipy is 1 - cosine_similarity
        dist_array = pdist(X_svd_reduced, metric='cosine')
        clustering_method = 'average' # Good for genomic conservation
    elif metric == 'euclidean':
        dist_array = pdist(X_svd_reduced, metric='euclidean')
        clustering_method = 'ward' # Minimizes variance within clusters
    else:
        raise ValueError("Metric not supported for hierarchical linkage.")

    # This matrix tells Seaborn exactly how to draw the tree
    row_linkage = fastcluster.linkage(dist_array, method=clustering_method)
    plt.figure(figsize=(6, 4))
    sns.histplot(dist_array, bins=30, kde=True, color='skyblue') 
    plt.title(f'Pairwise {metric.capitalize()} Distance Distribution')
    plt.savefig(f'{output_path}/{metric}_distance_distribution.png')
    plt.savefig(f'{output_path}/{metric}_distance_distribution.svg')
    plt.close()

    # --- 4. CLUSTERING ---
    # Option A: Hierarchical via Clustermap
    h_labels = get_optimal_hierarchical_clusters(X_reduced, row_linkage, clustering_method, metric, min_avg_sim=DEFAULT_SIMILARITY_THRESHOLD)
    # Merge singletons and enforce minimum cluster size of 2
    h_labels = merge_singleton_clusters(X_weighted, h_labels, metric=metric, similarity_threshold=0.3, min_cluster_size=2)

    # Option B: Density via HDBSCAN
    '''Because UMAP projects the data into a geometric space, you should generally use euclidean distance for HDBSCAN on UMAP embeddings, 
    even if the UMAP was built using cosine distance).'''
    d_labels = optimize_hdbscan(X_reduced, metric='euclidean')

    # --- 5. VISUALIZATION: UMAP SCATTER PLOTS ---
    # Call for Hierarchical Results
    umap_vis(
        clusters=h_labels, 
        reduced_data=X_reduced, 
        reduction_model=reduction_model, 
        output_path=output_path, 
        folder='hierarchical', 
        title=f'Hierarchical Clustering ({reduction_label} n={actual_components})',
        variance_model=svd_reducer
    )

    # Call for Density Results
    umap_vis(
        clusters=d_labels, 
        reduced_data=X_reduced, 
        reduction_model=reduction_model, 
        output_path=output_path, 
        folder='density', 
        title=f'HDBSCAN Density Clustering ({reduction_label} n={actual_components})',
        variance_model=svd_reducer
    )

    # --- 5. FINAL HEATMAP, UMAP & RETURN ---   
    binary_grey_cmap = mcolors.ListedColormap(['#FFFFFF', '#333333'])
    unique_h_labels = np.unique(h_labels)
    palette = sns.color_palette('turbo', len(unique_h_labels))
    lut = dict(zip(unique_h_labels, palette))
    row_colors = pd.Series(h_labels, index=X_cluster.index).map(lut) 
    data_to_plot = feature_matrix.loc[X_cluster.index]
    
    # 1. MANUALLY CALCULATE REORDERING TO BYPASS RECURSION CRASH
    # This computes the leaf order without drawing the full recursive tree
    # 1. Generate standalone dendrogram
    plt.figure(figsize=(12, 8))
    try:
        dendro_result = sch.dendrogram(row_linkage, labels=data_to_plot.index.astype(str), leaf_font_size=8, leaf_rotation=90)
        plt.title(f'Hierarchical Clustering Dendrogram ({metric.capitalize()} metric)', fontsize=14, fontweight='bold')
        plt.ylabel('Distance', fontsize=12)
        plt.xlabel('Sample', fontsize=12)
        
        # Add cluster legend
        cluster_patches = [mpatches.Patch(color=color, label=f'Cluster {label}') 
                         for label, color in lut.items()]
        plt.legend(handles=cluster_patches, 
                  title='Hierarchical Clusters', 
                  loc='upper right', 
                  bbox_to_anchor=(1.15, 1.0), 
                  fontsize='small')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, 'dendrogram_hierarchical.png'), dpi=150, bbox_inches='tight')
        plt.savefig(os.path.join(output_path, 'dendrogram_hierarchical.svg'), bbox_inches='tight')
        plt.close()
        logger.info('Standalone dendrogram saved successfully.')
        reordered_idx = dendro_result['leaves']
    except Exception as e:
        logger.warning(f'Dendrogram generation failed: {e}. Using default row order.')
        reordered_idx = np.arange(len(data_to_plot))
    
    # Reorder heatmap data based on dendrogram
    data_to_plot = data_to_plot.iloc[reordered_idx]
    current_row_colors = row_colors.iloc[reordered_idx]

    # 2. PLOT CLUSTERMAP WITH DENDROGRAMS (row order pre-computed to avoid recursion)
    try:
        g = sns.clustermap(data_to_plot, 
                    row_cluster=False,            # Use pre-computed dendro order
                    col_cluster=True,             # Cluster gene columns
                    row_colors=current_row_colors,
                    cmap=binary_grey_cmap, 
                    metric=metric,
                    figsize=(14, 12),
                    cbar_pos=(0.02, 0.8, 0.03, 0.15),
                    yticklabels=False,
                    xticklabels=False)
        logger.info('Heatmap with column dendrogram generated successfully.')
    except Exception as e:
        logger.warning(f'Clustermap with dendrograms failed: {e}. Generating basic heatmap.')
        g = sns.clustermap(data_to_plot,
                    row_cluster=False,
                    col_cluster=False,
                    yticklabels=False,
                    xticklabels=False,
                    row_colors=current_row_colors,
                    cmap=binary_grey_cmap,
                    metric=metric,
                    figsize=(12, 12))
    
    g.ax_cbar.set_visible(False)

    cluster_patches = [mpatches.Patch(color=color, label=f'Cluster {label}') 
                   for label, color in lut.items()]

    binary_patches = [
        mpatches.Patch(facecolor='#FFFFFF', label='Absent', edgecolor='grey', linewidth=0.5),
        mpatches.Patch(facecolor='#333333', label='Present')
    ]

    g.fig.legend(handles=cluster_patches, 
             title='Hierarchical Clusters', 
             loc='lower center', 
             bbox_to_anchor=(0.5, -0.05), 
             ncol=6, 
             fontsize='small')

    g.fig.legend(handles=binary_patches, 
             title='Genomic Context', 
             loc='upper left', 
             bbox_to_anchor=(0.05, 0.95))    
    
    plt.subplots_adjust(bottom=0.1)

    g.savefig(f'{output_path}/hierarchical_clustermap.png', bbox_inches='tight', dpi=150)
    plt.close()

    g.savefig(f'{output_path}/hierarchical_clustermap.svg', bbox_inches='tight', dpi=150)
    plt.close()
    
    active_genomes = X_cluster.index 
    
    hierarchical_series = pd.Series(h_labels, index=active_genomes)
    density_series = pd.Series(d_labels, index=active_genomes)

    return {
            'hierarchical': hierarchical_series, 
            'density': density_series,
        }, feature_matrix

def plot_cluster_func(cluster_df, goi_type, output_path, scale, gene_lable, folder, organism_name_from_file, svg, max_subplots=15):
    '''Generate plots for a cluster.'''    
    # Ensure the directory exists immediately
    final_output_dir = os.path.join(output_path, folder)
    os.makedirs(final_output_dir, exist_ok=True)

    split_data = cluster_df['cluster_color'].str.split('_h_', expand=True)
    cluster_df['cluster_color'] = split_data[0]
    if split_data.shape[1] > 1:
        cluster_df['hatch'] = split_data[1].fillna('').astype(str)
    else:
        cluster_df['hatch'] = ''

    cluster_df['hatch'] = cluster_df['hatch'].replace('nan', '')
    
    name2color = {gene: color for gene, color in zip(cluster_df['gene'], cluster_df['cluster_color'])}
    report_summary_json = []

    all_fts = FeatureList([
        Feature(row.type, start=row.start, stop=row.end, strand=row.orientation,
                meta={
                    'seqid': row.contig, 
                    'hit': row.genome, 
                    'organism_name': row.organism_name, 
                    'logic_name': row.gene,        
                    'label_name': getattr(row, 'display_name', row.gene), 
                    '_ftsviewer_hatch' : row.hatch
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
        
        align_features = FeatureList([f for f in fts if f.meta.get('type') == goi_type]) if scale_bool else None

        fig = fts.plot_ftsviewer(
            groupby='hit',
            axlabel='organism_name' if organism_name_from_file.lower() == 'yes' else 'seqid',
            label=lambda ft: ft.meta['label_name'] if label_bool else None,
            align=align_features,
            crop=500,
            sharex=scale_bool,
            colorby=lambda ft: ft.meta['logic_name'],
            color=name2color,
            figsize=(10, 1.5 * len(fts.groupby('hit'))), 
            ncols=1,
            with_ruler=False,
            labels_spacing=8, 
            fontdict={'fontsize': 7}
        )

        # Logic for Legend (only if labels are hidden)
        if not label_bool:
            seen_combinations = set()
            legend_handles = []
            for label, color in name2color.items():
                combo = (color, label)
                if isinstance(color, (str, tuple)) and combo not in seen_combinations:
                    seen_combinations.add(combo)
                    patch = Patch(facecolor=color, label=label)
                    legend_handles.append(patch)
            fig.legend(handles=legend_handles, title='Legend', bbox_to_anchor=(1.05, 1), loc='upper left')

        # Now these will run regardless of label_bool
        raw_cluster_id = cluster_df.get('cluster_label', pd.Series(['unknown'])).iloc[0]
        cluster_id = int(raw_cluster_id) if isinstance(raw_cluster_id, (np.integer,)) else raw_cluster_id
        img_filename = f'synteny_cluster_{cluster_id}_{fignum:03d}.png'
        if svg == 'yes':
            svg_filename = f'synteny_cluster_{cluster_id}_{fignum:03d}.svg'
        
        save_path_png = os.path.join(final_output_dir, img_filename)
        
        fig.savefig(save_path_png, bbox_inches='tight')
        if svg == 'yes':
            fig.savefig(os.path.join(final_output_dir, svg_filename), bbox_inches='tight')
        
        summary_json_entry = {
            'cluster_id': cluster_id,
            'image_path': img_filename,
            'genomes': cluster_df['organism_name'].unique().tolist(),
            'gene_count': int(len(cluster_df)),
            'sub_plot_index': int(fignum)
        }
        report_summary_json.append(summary_json_entry)
        plt.close(fig)

    return report_summary_json

def compute_cosine_similarity_matrix(matrix):
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

def write_labels(df, goi_type, output_path, folder):
    ''' Writing some stats '''
    output_dir = os.path.join(output_path, folder)
    os.makedirs(output_dir, exist_ok=True)
    center_row = df[df['type'] == goi_type]
    new_df = center_row[['cluster_label', 'organism_name', 'contig', 'gene', 'original_start', 'original_end', 'orientation']]

    output_file = os.path.join(output_dir, 'gene_of_interest.tsv')
    write_header = not os.path.exists(output_file)
    new_df.to_csv(output_file, mode='a', header=write_header, sep='\t', index=False)

def for_aln(cluster_label, clustering_algo, df, output_dir):
    target_dir = os.path.join(os.path.abspath(output_dir), clustering_algo)
    
    if not os.path.exists(target_dir):
        os.makedirs(target_dir, exist_ok=True)
        
    filename = os.path.join(target_dir, f'cluster_{cluster_label}.tsv')
    
    colors_by_genome = df.groupby('genome')['cluster_color'].apply(lambda x: ''.join(x)).reset_index()
    
    colors_by_genome.to_csv(filename, sep='\t', index=False)

def main():
    try:
        parser = setup_parser_func()
        args = parser.parse_args()

        logger.info('Reading data...')
        df, organism_name_from_file = prepare_dataframe_func(args.input_file, args.goi_type, args.gene_name, args.name_file)

        os.makedirs(args.output_path, exist_ok=True)
        logger.info('Process for clustering')
        clustering_results, feature_matrix_color = cluster_genomes_func(
            df, args.output_path, metric=args.metric
        )

        # Now iterate only over the metrics that were actually run
        for metric_name, labels in clustering_results.items():

            logger.info(f'Generating Plots and Stats for: {metric_name}')
            # Reset cluster_label for each metric
            genome_to_cluster = pd.Series(labels, index=feature_matrix_color.index).groupby(level=0).first()
            df['cluster_label'] = df['genome'].map(genome_to_cluster)

            unique_contigs = []  # Default: no filter
            if args.contig_cluster_file:
                filter_df = pd.read_csv(args.contig_cluster_file, header=None, names=['Contigs'])
                unique_contigs = filter_df['Contigs'].unique().tolist()

            cluster_stats = []
            master_summary = []
            for cluster_label, cluster_df in df.groupby('cluster_label'):
                if str(args.microsynteny_logo).lower() == 'yes':
                    for_aln(cluster_label, metric_name, cluster_df, args.output_path)

                # If the list is not empty, check if any contigs in the current cluster_df are in the unique_contigs list
                if unique_contigs and not cluster_df['contig'].isin(unique_contigs).any():
                    continue # Skip to the next cluster if no contigs from the list are found

                write_labels(cluster_df, args.goi_type, args.output_path, metric_name)
                cluster_entries = plot_cluster_func(cluster_df, args.goi_type, args.output_path, args.scale, args.gene_lable, metric_name, organism_name_from_file, args.svg)
                master_summary.extend(cluster_entries)

                cluster_size = cluster_df['genome'].nunique()
                color_count = cluster_df['cluster_color'].nunique()

                if cluster_label == -1:
                    cluster_stats.append({
                        'cluster_label': cluster_label,
                        'genomes': cluster_size})
                    continue  # Skip Noise Cluster

                #---------------------------------------------
                # Compute cosine Similarity for Clusters
                # calculation of common genes
                cluster_genome_vectors = feature_matrix_color.loc[genome_to_cluster[genome_to_cluster == cluster_label].index]

                if cluster_genome_vectors.shape[0] > 1:
                    # 1. Cosine Similarity
                    sim_matrix = compute_cosine_similarity_matrix(cluster_genome_vectors.values)
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
            stats_dir = os.path.join(args.output_path, metric_name)
            os.makedirs(stats_dir, exist_ok=True)
            stats_df = pd.DataFrame(cluster_stats)
            stats_df.to_csv(os.path.join(stats_dir, 'cluster_stats.csv'), index=False)
            final_json_path = os.path.join(args.output_path, f'summary_{metric_name}.json')
            with open(final_json_path, 'w') as f:
                json.dump(master_summary, f, indent=4)

        logger.info('Processing completed successfully!')
    except Exception as e:
        logger.error(f'An error occurred: {str(e)}')
        raise

if __name__ == '__main__':
    main()
