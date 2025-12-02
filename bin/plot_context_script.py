# (C) 2024, Maria Schreiber, MIT license
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sugar import Feature, FeatureList
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, cut_tree
from sklearn.cluster import DBSCAN, Birch
from sklearn.decomposition import PCA
import sys
sys.setrecursionlimit(100000)
from collections import Counter
from sklearn.metrics.pairwise import cosine_similarity


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
    parser.add_argument('--gene_name', required=False, default='db', help='gene_of_interest')
    parser.add_argument('--name_file', required=False, help='Renaming contigs to genome names')
    parser.add_argument('--gene_lable', required=False, choices=['yes', 'no'], default='yes', help='Add lables to gene [yes as default]')
    parser.add_argument('--cut_height_args', type=float, default=0, help='Cut threshold for dendogram clustering.')
    parser.add_argument('--contig_cluster_file', required=False, help='Filter for Contigs your are interested in.')
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

    if gene_name == 'db':
        df['gene'] = df['gene_db']

    if gene_name == 'id':
        df['gene'] = df['gene_id']
    print(df)
    if genome_name_tsv:
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

def cluster_genomes_func(df, output_path, output_ending, cut_height_para, cluster=2, threshold=0.1):
    '''Cluster genomes based on feature matrix.'''

    unique_colors = df['cluster_color'].unique()
    unique_genomes = df['genome'].unique()
    feature_matrix_color = pd.DataFrame(0, index=unique_genomes, columns=unique_colors)

    for _, row in df.iterrows():
        feature_matrix_color.at[row['genome'], row['cluster_color']] = 1

    # Create a mapping from genome to organism_name and get the organism names in the same order as feature_matrix_color.index
    genome_to_organism = df.set_index('genome')['organism_name'].to_dict()

    X = feature_matrix_color.values

    def refine_clusters_by_distance(X_cluster, clustered_labels_series, refinement_cut_factor=0.5):
        """Performs secondary, localized hierarchical clustering within each primary cluster
        to refine assignments."""

        # Ensure X_cluster is a DataFrame for robust index access
        if not isinstance(X_cluster, pd.DataFrame):
            X_cluster = pd.DataFrame(X_cluster, index=clustered_labels_series.index)

        unique_primary_clusters = sorted(clustered_labels_series.unique())
        current_max_label = clustered_labels_series.max()
        new_labels = clustered_labels_series.copy()

        for cluster_id in unique_primary_clusters:
            # 1. Isolate the data for the current primary cluster
            cluster_indices = clustered_labels_series[clustered_labels_series == cluster_id].index
            X_sub = X_cluster.loc[cluster_indices]

            # Only recluster if there are enough points
            if len(X_sub) <= 2:
                continue

            try:
                # 2. Perform new, localized hierarchical clustering
                # Use .values to explicitly pass a NumPy array to pdist
                matrix_sub = pdist(X_sub.values, metric='euclidean')
                clustering_sub = linkage(matrix_sub, method='average') # Assuming the same method

                # 3. Calculate a local cut height
                max_local_dist = max(clustering_sub[:, 2]) if clustering_sub.size > 0 else 0

                # Apply cut only if there is structure to cut
                if max_local_dist > 0.001:
                    local_cut_height = refinement_cut_factor * max_local_dist

                    # 4. Find new sub-cluster labels
                    local_sub_labels = fcluster(clustering_sub, t=local_cut_height, criterion="distance")
                    local_sub_labels_series = pd.Series(local_sub_labels, index=X_sub.index)

                    # Check if refinement actually created sub-clusters
                    num_sub_clusters = local_sub_labels_series.max()

                    if num_sub_clusters > 1:
                        # 5. Relabel sub-clusters with globally unique IDs

                        # The first sub-cluster keeps the original ID, others get new IDs
                        new_sub_cluster_labels = local_sub_labels_series.copy()

                        for sub_label in range(2, num_sub_clusters + 1):
                            current_max_label += 1
                            # Assign new global ID to the sub-cluster
                            new_sub_cluster_labels[local_sub_labels_series == sub_label] = current_max_label

                        # Update the global label series
                        new_labels.loc[cluster_indices] = new_sub_cluster_labels
                        print(f"  -> Split into {num_sub_clusters} sub-clusters. New global max ID: {current_max_label}")
                    else:
                        print("  -> No significant sub-structure found for refinement.")

                else:
                    print("  -> Trivial cluster structure (max internal distance is zero).")

            except ValueError as e:
                # Catch cases where linkage fails on very small or identical subsets
                print(f"  -> Could not recluster (Error: {e}). Skipping.")
                continue

        print("--- Refinement Complete ---")
        return new_labels

    def create_dendrogram_with_preassign(
        X, feature_matrix_color, metric, method, name, output_path, cut_height_para=0.25, filter_row_sum_1=True
        ):
        # 1 Pre-assign cluster labels for simple cases
        preassigned = {}
        clustering_mask = []

        for genome in feature_matrix_color.index:
            row_sum = feature_matrix_color.loc[genome].sum()
            if filter_row_sum_1 and row_sum == 1:
                preassigned[genome] = -2
                clustering_mask.append(False)
            elif filter_row_sum_1 and row_sum == 2:
                preassigned[genome] = -1
                clustering_mask.append(False)
            else:
                clustering_mask.append(True)
        clustering_mask = np.array(clustering_mask)
        X_cluster = X[clustering_mask]

        # 2 Primary hierarchical clustering
        matrix = pdist(X_cluster, metric=metric)
        clustering = linkage(matrix, method=method)

        max_dist = max(clustering[:, 2]) if clustering.size > 0 else 0
        cut_height = cut_height_para * max(clustering[:, 2])

        # 3 Cluster assignments
        clustered_labels = fcluster(clustering, t=cut_height, criterion="distance")

        if filter_row_sum_1:
            clustered_labels_series = pd.Series(
                clustered_labels, index=feature_matrix_color.index[clustering_mask]
            )

        if filter_row_sum_1 == False:
            clustered_labels_series = pd.Series(
                clustered_labels, index=feature_matrix_color.index
            )

        # Apply secondary clustering within primary cluster
        refined_labels_series = refine_clusters_by_distance(
            X_cluster, clustered_labels_series, refinement_cut_factor=0.5
        )

        # Re-integrate pre-assigned and clustered labels
        all_labels = pd.Series(index=feature_matrix_color.index, dtype=int)
        for genome in feature_matrix_color.index:
            if genome in preassigned:
                all_labels[genome] = preassigned[genome]
            else:
                all_labels[genome] = clustered_labels_series[genome]
        all_labels = all_labels.astype(int)

        # 4 Assign unique colors to clusters
        import matplotlib.colors as mcolors
        import matplotlib.patches as mpatches
        unique_clusters = sorted(set(clustered_labels))
        cmap = plt.get_cmap("tab20", len(unique_clusters))
        cluster_to_color = {
            cl: mcolors.to_hex(cmap(i)) for i, cl in enumerate(unique_clusters)
        }

        # 5 Node color
        def link_color_func(node_id):
            # leaf node
            if node_id < len(X_cluster):
                genome_id = feature_matrix_color.index[clustering_mask][node_id]
                cluster_id = clustered_labels_series[genome_id]
                return cluster_to_color[cluster_id]
            else:
                # Internal node
                dist = clustering[node_id - len(X_cluster), 2]
                if dist < cut_height:
                    left_child = int(clustering[node_id - len(X_cluster), 0])
                    return link_color_func(left_child)
                else:
                    return "grey"

        # 6 Plot dendrogram
        plt.figure(figsize=(10, 5))
        dendrogram(
            clustering,
            color_threshold=None,   # disable default coloring
            link_color_func=link_color_func,
            leaf_rotation=90,
            no_labels=True
        )

        cluster_counts = Counter(all_labels)
        # Create legend handles with cluster IDs and colors
        legend_handles = []
        for cluster_id, color in cluster_to_color.items():
            count = cluster_counts.get(cluster_id, 0)
            if cluster_id < 0:
                label = f"Preassigned {cluster_id}"
            else:
                label = f"Cluster {cluster_id} ({count} leaves)"
            patch = mpatches.Patch(color=color, label=label)
            legend_handles.append(patch)

        # Add legend to the plot (place it outside right of the figure)
        plt.legend(handles=legend_handles, title="Clusters", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.title(name)
        plt.xlabel("Genomes")
        plt.ylabel("Distance")
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        os.makedirs(f'{output_path}/{name.replace("_clustering", "")}', exist_ok=True)
        plt.savefig(f'{output_path}/{name.replace("_clustering", "")}/Tree.{name}.svg')
        plt.savefig(f'{output_path}/{name.replace("_clustering", "")}/Tree.{name}.png')
        plt.close()

        return clustering, cut_height, all_labels

    def non_or_one_color(cluster_labels, feature_matrix_color):
        for genome in feature_matrix_color.index:
            if feature_matrix_color.loc[genome].sum() == 1:
                cluster_labels[genome] = -1
            if feature_matrix_color.loc[genome].sum() == 0:
                cluster_labels[genome] = -2
        return cluster_labels.values

    def cluster_label(clustering, cut_height):
        # Assign clusters based on the chosen height
        fcluster_result = fcluster(clustering, t=cut_height, criterion='distance')
        cluster_labels_series = pd.Series(fcluster_result, index=feature_matrix_color.index)
        cluster_labels_values = non_or_one_color(cluster_labels_series, feature_matrix_color)

        return cluster_labels_values

    dbscan_labels_reduced, _ = dimensionreduction_visualize_pca(feature_matrix_color, output_path, 'dbscan', cluster, threshold)
    dbscan_labels_series = pd.Series(dbscan_labels_reduced, index=feature_matrix_color.index)
    dbscan = dbscan_labels_series.values
    #----------------------

    # Birch clustering (Balanced Iterative Reducing and Clustering using Hierarchies)
    #birch_labels_reduced, _ = dimensionreduction_visualize_pca(feature_matrix_color, output_path, 'birch', cluster, threshold)
    #birch_labels_series = pd.Series(birch_labels_reduced, index=feature_matrix_color.index)
    #birch = birch_labels_series.values
    #----------------------

    ward_clustering_cosinus, ward_cut_height_cosinus, ward_cosinus = create_dendrogram_with_preassign(
        X,                       # full feature matrix used for clustering
        feature_matrix_color,    # dataframe with your "colors"
        metric="cosine",         # cosine similarity distance
        method="ward",
        name="ward_cosinus_clustering",
        output_path=output_path,
        cut_height_para=cut_height_para
    )

    ward_clustering_jaccard, ward_cut_height_jaccard, ward_jaccard = create_dendrogram_with_preassign(
        X,
        feature_matrix_color,
        metric="jaccard",
        method="ward",
        name="ward_jaccard_clustering",
        output_path=output_path,
        cut_height_para=cut_height_para
    )

    return dbscan, ward_cosinus, ward_jaccard, feature_matrix_color # ward_dice, birch,

def dimensionreduction_visualize_pca(feature_matrix_color, output_path, folder, cluster=2, threshold=0.2):
    # Reduce data to 2D for plotting # Cluster the Reduced Data
    X = feature_matrix_color.values

    preassigned = {}  # genome -> cluster id
    clustering_mask = []  # keep track of which rows go into clustering
    for genome in feature_matrix_color.index:
        row_sum = feature_matrix_color.loc[genome].sum()
        if row_sum == 1: # White is a color = no color
            preassigned[genome] = -2
            clustering_mask.append(False)
        elif row_sum == 2:
            preassigned[genome] = -1   # only one color
            clustering_mask.append(False)
        else:
            clustering_mask.append(True)

    clustering_mask = np.array(clustering_mask)
    X_cluster = X[clustering_mask]

    pca = PCA(n_components=2)
    reduced_PCA = pca.fit_transform(X_cluster)

    if folder == 'dbscan':
        dbscan = DBSCAN(eps=threshold, min_samples=cluster)
        clusters = dbscan.fit_predict(reduced_PCA)

    #if folder == 'birch':
    #    birch = Birch(n_clusters=None)
    #    clusters = birch.fit_predict(reduced_PCA)

    # Now merge preassigned clusters into full cluster array
    full_clusters = np.empty(len(feature_matrix_color), dtype=int)
    full_clusters[:] = -99  # default placeholder

    # Fill preassigned
    preassigned_indices = np.where(~clustering_mask)[0]
    for idx in preassigned_indices:
        genome = feature_matrix_color.index[idx]
        full_clusters[idx] = preassigned[genome]

    # Fill newly assigned clusters
    clustered_indices = np.where(clustering_mask)
    full_clusters[clustered_indices] = clusters

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(reduced_PCA[:, 0], reduced_PCA[:, 1], c=clusters, cmap='tab20')
    unique_clusters, counts = np.unique(clusters, return_counts=True)
    cmap = plt.get_cmap('tab20', len(unique_clusters))
    legend_handles = [
        Patch(color=cmap(i), label=f'Cluster {cluster_label} ({count} points)')
        for i, (cluster_label, count) in enumerate(zip(unique_clusters, counts))
    ]
    plt.legend(handles=legend_handles, title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left')
    variance_1 = pca.explained_variance_ratio_[0] * 100
    variance_2 = pca.explained_variance_ratio_[1] * 100

    plt.xlabel(f'PC 1 ({variance_1:.2f}% variance)')
    plt.ylabel(f'PC 2 ({variance_2:.2f}% variance)')
    plt.title('Clustering on PCA-reduced data')
    plt.tight_layout()
    os.makedirs(f'{output_path}/{folder}', exist_ok=True)
    plt.savefig(f'{output_path}/{folder}/PCA.REDUCED.{folder}.png')
    plt.savefig(f'{output_path}/{folder}/PCA.REDUCED.{folder}.svg')
    plt.close()

    return full_clusters, clusters

def plot_cluster_func(cluster_df, output_path, output_ending, scale, gene_of_interest, gene_lable, folder, organism_name_from_file, max_subplots=15):

    '''Generate plots for a cluster.'''
    # Hatches are working now. The best way to use them is by setting ft.meta._ftsviewer_hatch attribute.
    # hatch_linewidth can be set in the same way starting with matplotlib 3.10
    # color can be set in the same way or via colorby and color arguments of plot_ftsviewer
    # needs sugar v1.0
    name2color = {gene: color for gene, color in zip(cluster_df['gene'], cluster_df['cluster_color'])}
    fts = FeatureList([
        Feature(row.type, start=row.start, stop=row.end, strand=row.orientation,
                meta={'seqid': row.contig, 'hit': row.genome, 'organism_name': row.organism_name, 'name': row.gene})
                for row in cluster_df.itertuples(index=False)
        ])

    scale = scale.lower() == 'yes'
    label = gene_lable.lower() == 'yes'
    align = fts.select(name=gene_of_interest) if scale else None
    if organism_name_from_file == 'yes':
        fig = fts.plot_ftsviewer(
            groupby='hit',
            axlabel='organism_name',
            label='name' if label else None,
            align=align, crop=500,
            sharex=scale, xticks=False,
            colorby='name', color=name2color,
            figsize=(10, 1.5*len(fts.groupby('hit'))), ncols=1,
            with_ruler=False,
            labels_spacing=8, fontdict={'fontsize': 7}
        )
    if organism_name_from_file == 'no':
        fig = fts.plot_ftsviewer(
            groupby='hit',
            axlabel='seqid',
            label='name' if label else None,
            align=align, crop=500,
            sharex=scale, xticks=False,
            colorby='name', color=name2color,
            figsize=(10, 1.5*len(fts.groupby('hit'))), ncols=1,
            with_ruler=False,
            labels_spacing=8, fontdict={'fontsize': 7}
        )

    if not label: # Add legend to the subplot
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

        #legend_handles = [
        #    Patch(facecolor=color, label=label)
        #    for label, color in name2color.items()
        #    if isinstance(color, tuple) and len(color) == 3 and all(0 <= c <= 1 for c in color)
        #]
        #fig.legend(handles=legend_handles, title='Legend', bbox_to_anchor=(1.05, 1), loc='upper left')

    count = f"{cluster_df['cluster_label'].iloc[0]}"
    fig.savefig(f'{output_path}/{folder}/synteny_cluster_{count}.png', bbox_inches='tight')
    fig.savefig(f'{output_path}/{folder}/synteny_cluster_{count}.svg', bbox_inches='tight')
    plt.close(fig)

# def old_plot_cluster_func(cluster_df, output_path, output_ending, scale, gene_of_interest, gene_lable, folder, max_subplots=15):
#     genome_groups = list(cluster_df.groupby('genome'))
#     num_subplots = int(np.ceil(len(genome_groups) / max_subplots))
#     for subplot_idx in range(num_subplots):
#         merged = {} # Merge all color dictionaries
#         subplot_genomes = genome_groups[
#             subplot_idx * max_subplots: (subplot_idx + 1) * max_subplots
#         ]
#         if not subplot_genomes:  # Added validation
#             continue
#         fig, axes = plt.subplots(len(subplot_genomes), 1, figsize=(12, len(subplot_genomes) * 1.2))
#         axes = [axes] if len(subplot_genomes) == 1 else axes
#         for i, (genome_name, genome_df) in enumerate(subplot_genomes):
#             dicts_subplot, df_subplot = plot_genome_func(genome_df, axes[i], scale, gene_of_interest, gene_lable)
#             merged.update(dicts_subplot)
#             unique_organisms_name = genome_df['organism_name'].unique().tolist()[0]
#             axes[i].set_title(unique_organisms_name, fontstyle='italic', loc='left')
#         if gene_lable=='no': # Add legend to the subplot
#             legend_handles = [
#                 Patch(facecolor=color, label=label)
#                 for label, color in merged.items()
#                 if isinstance(color, tuple) and len(color) == 3 and all(0 <= c <= 1 for c in color)
#             ]
#             fig.legend(handles=legend_handles, title='Legend', bbox_to_anchor=(1.05, 1), loc='upper left')
#         # If further spacing needed:
#         plt.subplots_adjust(left=0.1, right=0.8, hspace=0.75)  # shrink plot area for
#         count = f"{cluster_df['cluster_label'].iloc[0]}.{subplot_idx}"
#         fig.savefig(f'{output_path}/{folder}/{count}.png', bbox_inches='tight')
#         fig.savefig(f'{output_path}/{folder}/{count}.svg', bbox_inches='tight')
#         plt.close(fig)
#
# def scale_genes(df, gene_name, target_length=50):
#     # Find the row with the specified gene
#     goi_row = df[df['gene'] == gene_name].iloc[0]
#     # Calculate the current length of the CrfA gene
#     crfa_length = goi_row['end'] - goi_row['start'] + 1
#     # Calculate scale factor
#     scale_factor = target_length / crfa_length
#     # Reference position is the minimum start to maintain relative distances
#     ref_start = df['start'].min()
#     # Function to calculate adjusted positions
#     def adjusted_positions(row):
#         original_length = row['end'] - row['start'] + 1
#         scaled_length = int(original_length * scale_factor)
#         scaled_start = int((row['start'] - ref_start) * scale_factor)
#         scaled_end = scaled_start + scaled_length - 1
#         return pd.Series({'adjusted_start': scaled_start, 'adjusted_end': scaled_end})
#     # Apply the function and assign new columns
#     df[['adjusted_start', 'adjusted_end']] = df.apply(adjusted_positions, axis=1)
#     return df
#
# def plot_genome_func(genome_df, ax, scale, gene_of_interest, gene_lable='no'):
#     '''Plot a single genome.'''
#     df = genome_df.sort_values(by=['start'])
#     df['organism_name'] = genome_df['organism_name']
#     min_start = df['start'].min()
#     df['start'] = df['start'].astype(int)
#     df['length'] = df['length'].astype(int)
#     min_start = int(df['start'].min())
#     if scale.lower() == 'yes':
#         df = scale_genes(df, gene_of_interest)
#         seqlen_distance = (df['adjusted_end'].max() - df['adjusted_start'].min()).astype(int)# Drop the original 'start' and 'end' columns
#     else:
#         max_distance = 0
#         df['adjusted_start'] = (df['start'] - min_start) / 1000
#         df['adjusted_end'] = (df['end'] - min_start) / 1000
#         min_start = genome_df['start'].min()
#         genome_locus_length = genome_df['end'].max() - min_start
#         max_distance = max(max_distance, genome_locus_length)
#         seqlen_distance = max_distance / 1000
#     features = [Feature(row.type, start=row.adjusted_start, stop=row.adjusted_end, strand=row.orientation,
#                         meta={'seqid': row.contig, 'organism_name': row.organism_name, 'name': row.gene, 'cluster_color': row.cluster_color})
#                 for _, row in df.iterrows()]
#     color_dic = dict(zip(df['gene'], df['cluster_color']))
#     # New dictionaries to store cleaned colors and hatches
#     clean_color_dic = {}
#     hatch_dic = {}
#     for gene, color_hatch in color_dic.items():
#         if '_h_' in color_hatch:
#             # Split into color and hatch parts
#             color, hatch = color_hatch.split('_h_')
#             clean_color_dic[gene] = color
#             hatch_dic[gene] = hatch
#         else:
#             clean_color_dic[gene] = color_hatch  # just hex color, no hatch
#     FeatureList(features).plot_ftsviewer(
#         ax=ax, label=None if gene_lable=='no' else 'name',
#         colorby='name', color=clean_color_dic,
#         seqlen=seqlen_distance, figsize=(7, 5),
#         with_ruler=False, show=False,
#         labels_spacing=60, fontdict={'fontsize': 7})
#     # color dic for legend
#     return clean_color_dic, df

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
        print(df)
        print("Process for clustering")
        dbscan_lables, ward_labels_cosinus, ward_labels_jaccard, feature_matrix_color = cluster_genomes_func(df, args.output_path, args.output_ending, args.cut_height_args, args.cluster, args.threshold)
        print("Finished clustering")
        lst = ['ward_cosinus', 'ward_jaccard', 'dbscan'] # 'dice' 'birch' 'complete', 'average', 'single', 'centroid', 'weighted', 'kmeans', 'hdbscan',
        i = 0

        for labels in [ward_labels_cosinus, ward_labels_jaccard, dbscan_lables]: # , ward_labels_dice, birch_labels
            print("Plots " + str(lst[i]))
            # correlation of cluster labels to each genome
            genome_to_cluster = pd.Series(labels, index=feature_matrix_color.index).groupby(level=0).first()

            df['cluster_label'] = df['genome'].map(genome_to_cluster)
            unique_contigs = []  # Default: no filter
            if args.contig_cluster_file:
                filter_df = pd.read_csv(args.contig_cluster_file, header=None, names=["Contigs"])
                unique_contigs = filter_df['Contigs'].unique().tolist()

            cluster_stats = []
            for cluster_label, cluster_df in df.groupby('cluster_label'):

                for_aln(cluster_label, lst[i], cluster_df, args.output_path)

                # If the list is not empty, check if any contigs in the current cluster_df are in the unique_contigs list
                if unique_contigs and not cluster_df['contig'].isin(unique_contigs).any():
                    continue # Skip to the next cluster if no contigs from the list are found

                write_labels(cluster_df, args.gene_of_interest, args.output_path, lst[i])
                plot_cluster_func(cluster_df, args.output_path, args.output_ending, args.scale, args.gene_of_interest, args.gene_lable, lst[i], organism_name_from_file)
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
