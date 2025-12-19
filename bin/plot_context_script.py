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
from scipy.stats import pearsonr
from sklearn.cluster import DBSCAN, Birch
from sklearn.decomposition import PCA
import sys
sys.setrecursionlimit(100000)
from collections import Counter
from sklearn.metrics.pairwise import cosine_similarity
import seaborn as sns

def setup_parser_func():
    '''Set up and return the argument parser.'''
    parser = argparse.ArgumentParser(description='Process genomic data and generate plots.')
    parser.add_argument('--scale', required=True, choices=['yes', 'no'], help='Nucleotide Scale')
    parser.add_argument('--input_file', required=True, help='Path to input file')
    parser.add_argument('--output_path', required=True, help='Path to final png or svg')
    parser.add_argument('--cluster', type=int, default=2, help='Minimal size for a cluster')
    parser.add_argument('--threshold', type=float, default=0.9, help='Similarity threshold for clustering')
    parser.add_argument('--gene_of_interest', required=True, help='gene_of_interest')
    parser.add_argument('--gene_name', required=False, default='db', help='gene_of_interest')
    parser.add_argument('--name_file', required=False, help='Renaming contigs to genome names')
    parser.add_argument('--gene_lable', required=False, choices=['yes', 'no'], default='yes', help='Add lables to gene [yes as default]')
    parser.add_argument('--cut_height_args', type=float, default=0, help='Cut threshold for dendogram clustering.')
    parser.add_argument('--contig_cluster_file', required=False, help='Filter for Contigs your are interested in.')
    parser.add_argument('--microsynteny_logo', required=False, choices=['yes', 'no'], default='no', help='Filter for Contigs your are interested in.')
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
        df['gene'] = df['gene_db'] # Name from Pfam or Rfam

    if gene_name == 'id':
        df['gene'] = df['gene_id']

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

def cluster_genomes_func(df, output_path, cut_height_para, dbscan_size=2, dbscan_threshold=0.1):
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
    
    def test_suitablty_for_pca(X):

        first_df = pd.DataFrame(X)
        # 1. Finde alle Spalten, die nur einen einzigen eindeutigen Wert enthalten
        constant_cols = [col for col in first_df.columns if first_df[col].nunique() <= 1]
        df = first_df.drop(columns=constant_cols)

        # Pearson-Korrelationsmatrix
        correlation_matrix = df.corr(method='pearson')

        plt.figure(figsize=(8, 6))
        sns.heatmap(
            correlation_matrix,
            annot=False,
            cmap='coolwarm',
            linewidths=.5,
            cbar_kws={'label': 'Pearson Korrelationskoeffizient'}
        )
        plt.title('Korrelations-Heatmap')
        plt.savefig(f'{output_path}/korrelations_heatmap.png')
        plt.savefig(f'{output_path}/korrelations_heatmap.svg')
        plt.close()

        df_cols = df.columns
        p_values_matrix = np.zeros((len(df_cols), len(df_cols)))

        for i in range(len(df_cols)):
            for j in range(len(df_cols)):
                # Die Korrelation einer Variable mit sich selbst ist 1, P-Wert ist nicht anwendbar
                if i == j:
                    p_values_matrix[i, j] = 0.0
                    continue

                # pearsonr gibt (r, p-value) 
                r, p_value = pearsonr(df[df_cols[i]], df[df_cols[j]])
                p_values_matrix[i, j] = p_value

        p_values_df = pd.DataFrame(p_values_matrix, index=df_cols, columns=df_cols)

        plt.figure(figsize=(8, 6))
        sns.heatmap(
            p_values_df,
            annot=False,
            cmap='viridis_r',        # Dark color = low p (signifikant)
            vmax=0.05,               # Alle Werte > 0.05 sind blass
            linewidths=.5,
            cbar_kws={'label': 'P-Wert (Signifikanzniveau)', 'ticks': [0.0, 0.01, 0.02, 0.03, 0.04, 0.05]}
        )
        plt.title('P-value')
        plt.savefig(f'{output_path}/p_value.png')
        plt.savefig(f'{output_path}/p_value.svg')
        plt.close()

    clustering_mask = np.array(clustering_mask)
    X_cluster = X[clustering_mask]
    test_suitablty_for_pca(X_cluster)
    # Dimension reduction
    pca = PCA(n_components=2)
    reduced_PCA = pca.fit_transform(X_cluster)
    full_clusters = np.empty(len(feature_matrix_color), dtype=int)
    full_clusters[:] = -99  # default placeholder

    # Fill preassigned
    preassigned_indices = np.where(~clustering_mask)[0]
    for idx in preassigned_indices:
        genome = feature_matrix_color.index[idx]
        full_clusters[idx] = preassigned[genome]

    clustered_indices = np.where(clustering_mask)
    
    def pca_vis(full_clusters, clustered_indices, clusters, output_path, folder, title):
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
        plt.title(title)
        plt.tight_layout()
        os.makedirs(f'{output_path}/{folder}', exist_ok=True)
        plt.savefig(f'{output_path}/{folder}/PCA.REDUCED.{folder}.png')
        plt.savefig(f'{output_path}/{folder}/PCA.REDUCED.{folder}.svg')
        plt.close()

    # 2 Clustering
    def density_clustering(eps, min_samples, reduced_PCA, pca_clusters, folder):

        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        dbscan_clusters = dbscan.fit_predict(reduced_PCA)

        clustered_indices = np.where(clustering_mask)
        pca_clusters[clustered_indices] = dbscan_clusters

        title = f'{folder} with eps={eps}, min_samples={min_samples} on PCA-reduced Data (n={len(X_cluster)})'
        pca_vis(pca_clusters, clustered_indices, dbscan_clusters, output_path, folder, title)

        dbscan_labels_series = pd.Series(pca_clusters, index=feature_matrix_color.index)
        dbscan = dbscan_labels_series.values.copy()
        return dbscan

    dbscan_result = density_clustering(dbscan_threshold, dbscan_size, reduced_PCA, full_clusters, 'dbscan')

    def hierarchical_clustering(reduced_PCA, feature_matrix_color, pca_clusters, metric, method, folder, output_path, cut_height_para=0.25, filter_row_sum_1=True):
        import matplotlib.colors as mcolors
        import matplotlib.patches as mpatches
        
        def dendogram_vis(linkage_result, all_labels, title, cut_height, clustered_labels_series, cluster_to_color):

            # Node color
            def link_color_func(node_id):
                # leaf node
                if node_id < len(X_cluster):
                    genome_id = feature_matrix_color.index[clustering_mask][node_id]
                    cluster_id = clustered_labels_series[genome_id]
                    return cluster_to_color[cluster_id]
                else:
                    # Internal node
                    dist = linkage_result[node_id - len(X_cluster), 2]
                    if dist < cut_height:
                        left_child = int(linkage_result[node_id - len(X_cluster), 0])
                        return link_color_func(left_child)
                    else:
                        return "grey"
                        
            plt.figure(figsize=(10, 5))
            dendrogram(
                linkage_result,
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
            plt.title(title)
            plt.xlabel("Genomes")
            plt.ylabel("Distance")
            plt.tight_layout(rect=[0, 0, 0.85, 1])
            os.makedirs(f'{output_path}/{folder.replace("_clustering", "")}', exist_ok=True)
            plt.savefig(f'{output_path}/{folder.replace("_clustering", "")}/Tree.{folder}.svg')
            plt.savefig(f'{output_path}/{folder.replace("_clustering", "")}/Tree.{folder}.png')
            plt.close()

        # Primary hierarchical clustering
        dist_matrix = pdist(reduced_PCA, metric=metric)
        linkage_result = linkage(dist_matrix, method=method)

        max_dist = max(linkage_result[:, 2]) if linkage_result.size > 0 else 0
        cut_height = cut_height_para * max(linkage_result[:, 2])

        # 3 Cluster assignments
        clustered_labels = fcluster(linkage_result, t=cut_height, criterion="distance")

        if filter_row_sum_1:
            clustered_labels_series = pd.Series(
                clustered_labels, index=feature_matrix_color.index[clustering_mask]
            )

        if filter_row_sum_1 == False:
            clustered_labels_series = pd.Series(
                clustered_labels, index=feature_matrix_color.index
            )

        # Re-integrate pre-assigned and clustered labels
        all_labels = pd.Series(index=feature_matrix_color.index, dtype=int)
        for genome in feature_matrix_color.index:
            if genome in preassigned:
                all_labels[genome] = preassigned[genome]
            else:
                all_labels[genome] = clustered_labels_series[genome]
        all_labels = all_labels.astype(int).copy()

        # 4 Assign unique colors to clusters
        unique_clusters = sorted(set(clustered_labels))
        cmap = plt.get_cmap("tab20", len(unique_clusters))
        cluster_to_color = {
            cl: mcolors.to_hex(cmap(i)) for i, cl in enumerate(unique_clusters)
        }

        # 5 Plot dendrogram and pca
        title = f'{folder} with cut_height={cut_height_para} on PCA-reduced Data (n={len(X_cluster)})'
        dendogram_vis(linkage_result, all_labels, title, cut_height, clustered_labels_series, cluster_to_color)
        pca_vis(pca_clusters, clustered_indices, clustered_labels, output_path, folder, title)

        return all_labels

    ward_cosinus_result = hierarchical_clustering(
        reduced_PCA,             # dimension reduced matrix used for clustering
        feature_matrix_color,    # dataframe with your "colors"
        full_clusters,
        metric="cosine",         # cosine similarity distance
        method="ward",
        folder="ward_cosinus",
        output_path=output_path,
        cut_height_para=cut_height_para
    )

    ward_euclidean_result = hierarchical_clustering(
        reduced_PCA,
        feature_matrix_color,
        full_clusters,
        metric="euclidean",
        method="ward",
        folder="ward_euclidean",
        output_path=output_path,
        cut_height_para=cut_height_para
    )

    return dbscan_result, ward_cosinus_result, ward_euclidean_result, feature_matrix_color

def plot_cluster_func(cluster_df, output_path, scale, gene_of_interest, gene_lable, folder, organism_name_from_file, max_subplots=15):

    '''Generate plots for a cluster.'''
    # Hatches are working now. The best way to use them is by setting ft.meta._ftsviewer_hatch attribute.
    # hatch_linewidth can be set in the same way starting with matplotlib 3.10
    # color can be set in the same way or via colorby and color arguments of plot_ftsviewer
    # needs sugar v1.0
    name2color = {gene: color for gene, color in zip(cluster_df['gene'], cluster_df['cluster_color'])}
    all_fts = FeatureList([
        Feature(row.type, start=row.start, stop=row.end, strand=row.orientation,
                meta={'seqid': row.contig, 'hit': row.genome, 'organism_name': row.organism_name, 'name': row.gene})
                for row in cluster_df.itertuples(index=False)
        ])

    scale = scale.lower() == 'yes'
    label = gene_lable.lower() == 'yes'
    groups = list(all_fts.groupby('hit').values())
    for fignum in range(0, 1 + (len(groups)-1)//max_subplots):
        fts = FeatureList()
        for g in groups[fignum * max_subplots:(fignum + 1) * max_subplots]:
            fts.extend(g)
        align = fts.select(name=gene_of_interest) if scale else None
        fig = fts.plot_ftsviewer(
            groupby='hit',
            axlabel='organism_name' if organism_name_from_file.lower() == 'yes' else 'seqid',
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

        count = f"{cluster_df['cluster_label'].iloc[0]}"
        fig.savefig(f'{output_path}/{folder}/synteny_cluster_{count}_{fignum:03d}.png', bbox_inches='tight')
        fig.savefig(f'{output_path}/{folder}/synteny_cluster_{count}_{fignum:03d}.svg', bbox_inches='tight')
        plt.close(fig)

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

        print("Process for clustering")
        dbscan_lables, ward_labels_cosinus, ward_labels_euclidean, feature_matrix_color = cluster_genomes_func(df, args.output_path, args.cut_height_args, args.cluster, args.threshold)

        print("Finished clustering")
        lst = ['ward_cosinus', 'ward_euclidean', 'dbscan'] # 'dice' 'birch' 'complete', 'average', 'single', 'centroid', 'weighted', 'kmeans', 'hdbscan',
        i = 0

        for labels in [ward_labels_cosinus, ward_labels_euclidean, dbscan_lables]: # , ward_labels_dice, birch_labels
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
                if args.microsynteny_logo == 'yes':
                    for_aln(cluster_label, lst[i], cluster_df, args.output_path)

                # If the list is not empty, check if any contigs in the current cluster_df are in the unique_contigs list
                if unique_contigs and not cluster_df['contig'].isin(unique_contigs).any():
                    continue # Skip to the next cluster if no contigs from the list are found

                write_labels(cluster_df, args.gene_of_interest, args.output_path, lst[i])
                plot_cluster_func(cluster_df, args.output_path, args.scale, args.gene_of_interest, args.gene_lable, lst[i], organism_name_from_file)
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
