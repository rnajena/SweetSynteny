#!/usr/bin/env python3
# (C) 2024, Maria Schreiber, MIT license
"""Cluster gene neighbourhoods and generate synteny/cluster summaries."""
import sys
import argparse
import colorsys
import json
import logging
import os
import warnings
from pathlib import Path
import subprocess
from io import StringIO
from Bio import SeqIO
try:
    import fastcluster
except ImportError:
    fastcluster = None

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

try:
    import umap
except ImportError:
    umap = None

try:
    import networkx as nx
except ImportError:
    nx = None

from matplotlib.patches import Patch
from scipy.cluster import hierarchy as sch
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import pdist
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.metrics.pairwise import cosine_similarity, pairwise_distances

try:
    from sugar import Feature, FeatureList
except ImportError:
    Feature = None
    FeatureList = None

from plot_colors import load_colors
from plot_microsynteny_logo import (
    align_microsynteny_neighborhoods,
    plot_microsynteny_logo,
)
from microsynteny_stats import plot_cluster_statistics
try:
    from get_cluster_fastas import extract_cluster_fastas
except ImportError:
    extract_cluster_fastas = None

CORE_GENE_THRESHOLD = 0.99
RARE_GENE_MIN_COUNT = 2
DEFAULT_CUT_HEIGHT = 0.70   # only used as a fixed cut if the caller passes 0; the default is now an automatic cut
DEFAULT_MIN_CLUSTER_SIZE = 2
DEFAULT_MIN_COVERAGE = 0.50
DEFAULT_MIN_CLUSTER_COSINE = 0.5
DEFAULT_ASSIGN_MIN_SIMILARITY = 0.5
DEFAULT_MIN_SHARED_FAMILIES = 2      # a cluster must share >= this many families ...
DEFAULT_SIGNATURE_FRACTION = 0.6     # ... each present in >= this fraction of its members
DEFAULT_MIN_FEATURES = 3             # neighbourhoods with fewer eligible families are not clustered
DEFAULT_SUPERCLUSTER_SIM = 0.30      # loose second level (0 = off)
MAX_LEGEND_ITEMS = 40
DEFAULT_MAX_HIERARCHICAL_SAMPLES = 10000
DEFAULT_MAX_HEATMAP_FEATURES = 500

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

def setup_parser_func():
    parser = argparse.ArgumentParser(description="Cluster genomic neighbourhoods and generate plots.")
    parser.add_argument("--scale", choices=["yes", "no"], default="yes")
    parser.add_argument("--input_file", required=True, help="Input TSV file.")
    parser.add_argument("--output_path", required=True, help="Output directory.")
    parser.add_argument("--neighbour_dir", default=None, help="Directory containing neighbourhood FASTA files for cluster FASTA extraction.")
    parser.add_argument("--svg", choices=["yes", "no"], default="no")
    parser.add_argument("--gene_of_interest", required=True)
    parser.add_argument("--goi_type", required=True)
    parser.add_argument("--gene_name", choices=["db", "id"], default="db")
    parser.add_argument("--name_file", default=None)
    parser.add_argument("--gene_lable", "--gene_label", dest="gene_lable", choices=["yes", "no"], default="yes")
    parser.add_argument("--contig_cluster_file", default=None)
    parser.add_argument("--cut_height_args", type=float, default=None, help="Fixed cosine-distance cut for the dendrogram. Omit (or pass 0) to choose the cut automatically (see --min_coverage / --min_cluster_cosine).")
    parser.add_argument("--min_cluster_size", type=int, default=DEFAULT_MIN_CLUSTER_SIZE, help="Minimum number of neighbourhoods required to form a cluster (2 gives the highest coverage; 5 leaves many neighbourhoods unclustered).")
    parser.add_argument("--min_coverage", type=float, default=DEFAULT_MIN_COVERAGE, help="Automatic cut: fraction of neighbourhoods that should fall in clusters of at least --min_cluster_size.")
    parser.add_argument("--min_cluster_cosine", type=float, default=DEFAULT_MIN_CLUSTER_COSINE, help="Automatic cut: minimum mean pairwise cosine similarity allowed in every cluster.")
    parser.add_argument("--assign_leftovers", choices=["yes", "no"], default="yes", help="Attach unclustered neighbourhoods to the most similar cluster when the mean cosine similarity to it is >= --assign_min_similarity.")
    parser.add_argument("--assign_min_similarity", type=float, default=DEFAULT_ASSIGN_MIN_SIMILARITY, help="Minimum mean cosine similarity to a cluster for a leftover neighbourhood to be attached to it.")
    parser.add_argument("--density_clustering", choices=["yes", "no", "only"], default="yes", help="HDBSCAN (written to the 'density' sub-folder). yes (default) = run it next to the hierarchical clustering and write clustering_comparison.tsv (ARI + per-neighbourhood agreement); no = hierarchical only (HDBSCAN still runs if hierarchical clustering is skipped for size); only = HDBSCAN only.")
    parser.add_argument("--hdbscan_selection", choices=["eom", "leaf"], default="eom", help="HDBSCAN cluster_selection_method: eom (fewer, larger clusters) or leaf (finer, more homogeneous clusters).")
    parser.add_argument("--logo_for_noise", choices=["yes", "no"], default="no", help="Also draw logos / statistics for the noise group (-1) and the paralog overflow group (9999). Off by default: they are not clusters, and aligning thousands of unrelated neighbourhoods is slow.")
    parser.add_argument("--adjacency_across_gaps", choices=["yes", "no"], default="yes", help="yes: adjacent-gene pairs are formed between consecutive coloured genes, ignoring uncoloured (white) genes in between. no: pairs require physical adjacency (an uncoloured gene breaks the pair).")
    parser.add_argument("--microsynteny_logo", choices=["yes", "no"], default="no")
    parser.add_argument("--metric", choices=["cosine", "euclidean"], default="cosine")
    parser.add_argument("--remove_core_genes", choices=["yes", "no"], default="yes")
    parser.add_argument("--remove_rare_genes", choices=["yes", "no"], default="yes")
    parser.add_argument("--core_threshold", type=float, default=CORE_GENE_THRESHOLD)
    parser.add_argument("--min_gene_prevalence", type=_parse_float_or_none, default=None, metavar="0-1", help="Fraction representing the top X%% of most frequent genes to include (e.g., 0.25 for top 25%%).")
    parser.add_argument("--rare_min_count", type=int, default=RARE_GENE_MIN_COUNT)
    parser.add_argument("--max_hierarchical_samples", type=int, default=DEFAULT_MAX_HIERARCHICAL_SAMPLES)
    parser.add_argument("--max_heatmap_features", type=int, default=DEFAULT_MAX_HEATMAP_FEATURES)
    parser.add_argument("--max_hdbscan_candidates", type=int, default=10)
    parser.add_argument("--min_samples", type=int, default=1, help="HDBSCAN min_samples (higher = more noise). Was hard-coded to 1; the min_samples in your JSON was never read.")
    parser.add_argument("--noise_penalty", type=float, default=0.25, help="Weight of the noise fraction when HDBSCAN picks min_cluster_size (higher = prefers fewer noise points).")
    parser.add_argument("--keep_small_hier_clusters", choices=["yes", "no"], default="no", help="no = hierarchical clusters smaller than --min_cluster_size become noise (-1).")
    parser.add_argument("--max_silhouette_samples", type=int, default=1000)
    parser.add_argument("--color_scheme_file", default=None)
    parser.add_argument("--paralog_resolution", choices=["none", "reassign", "overflow_only"], default="none")
    parser.add_argument("--feature_mode", choices=["content", "combined", "synteny"], default="combined", help="Gene-neighbourhood representation: family content, or content plus GOI-relative position/orientation and adjacent-gene synteny.")
    parser.add_argument("--position_weight", type=float, default=0.1, help="Weight of GOI-relative position/orientation features in combined mode.")
    parser.add_argument("--synteny_weight", type=float, default=0.3, help="Weight of adjacent-gene pair features in combined/synteny mode.")
    parser.add_argument("--collapse_overlaps", choices=["yes", "no"], default="yes", help="Merge features with the same family name that overlap (e.g. infernal hit + GFF annotation of the same gene).")
    parser.add_argument("--overlap_fraction", type=float, default=0.5, help="Overlap (fraction of the shorter feature) needed to call two same-family features duplicates.")
    parser.add_argument("--paralog_sim_threshold", type=float, default=0.5)
    parser.add_argument("--linkage_method", choices=["average", "ward"], default="average", help="ward requires --metric euclidean.")
    parser.add_argument("--cohesion_stat", choices=["worst", "median", "mean"], default="worst", help="Statistic over clusters that --min_cluster_cosine is applied to in the automatic cut. worst = every cluster must reach it (strict; one weak pair holds the cut low); median/mean = typical cluster must reach it (allows a higher cut / more coverage).")
    parser.add_argument("--min_shared_families", type=int, default=DEFAULT_MIN_SHARED_FAMILIES, help="A cluster is kept only if at least this many families are present in >= --signature_fraction of its members (0 = off). Removes clusters glued together by a single shared gene.")
    parser.add_argument("--signature_fraction", type=float, default=DEFAULT_SIGNATURE_FRACTION)
    parser.add_argument("--min_features_per_neighbourhood", type=int, default=DEFAULT_MIN_FEATURES, help="Neighbourhoods with fewer eligible families than this (e.g. very short contigs) are not clustered; they are listed in excluded_neighbourhoods.tsv.")
    parser.add_argument("--rare_unit", choices=["genome", "neighbourhood"], default="genome", help="Count --rare_min_count in distinct genomes (default) or in neighbourhoods. 'genome' stops multi-contig/near-identical genomes from making a family look recurrent.")
    parser.add_argument("--supercluster_similarity", type=float, default=DEFAULT_SUPERCLUSTER_SIM, help="Mean-similarity level of a second, looser grouping written as the 'supercluster' column of cluster_membership.tsv (0 = off).")
    parser.add_argument("--network_similarity_threshold", type=float, default=0.5, help="Minimum cosine similarity threshold to draw an edge in the network graph.")
    return parser

def _parse_float_or_none(value):
    if value is None or str(value).lower() in ("none", "null", ""):
        return None
    return float(value)

def map_ori_func(value):
    return {"antisense": "-", "sense": "+"}.get(value, value)

def collapse_overlapping_duplicates(group, goi_type, min_overlap=0.5):
    """Drop features that are a second annotation of the same locus.

    Two features count as duplicates if they have the same family name ('gene') and overlap by
    >= min_overlap of the shorter one (e.g. an infernal hit and the GFF gene it falls in).
    The GOI-type feature is kept; otherwise the longer one is kept.
    """
    group = group.reset_index(drop=True)
    drop = set()
    order = group.sort_values("start", kind="stable").index.tolist()
    for ai, i in enumerate(order):
        if i in drop:
            continue
        for j in order[ai + 1:]:
            if j in drop:
                continue
            if group.at[j, "start"] >= group.at[i, "end"]:
                break                                   # sorted by start: nothing further overlaps i
            if group.at[i, "gene"] != group.at[j, "gene"]:
                continue
            ov = min(group.at[i, "end"], group.at[j, "end"]) - max(group.at[i, "start"], group.at[j, "start"])
            shorter = min(group.at[i, "end"] - group.at[i, "start"], group.at[j, "end"] - group.at[j, "start"])
            if shorter <= 0 or ov / shorter < min_overlap:
                continue
            i_goi = str(group.at[i, "type"]) == str(goi_type)
            j_goi = str(group.at[j, "type"]) == str(goi_type)
            if i_goi == j_goi:
                loser = j if (group.at[i, "end"] - group.at[i, "start"]) >= (group.at[j, "end"] - group.at[j, "start"]) else i
            else:
                loser = j if i_goi else i
            drop.add(loser)
            if loser == i:
                break
    if drop:
        logger.info("Collapsed %d overlapping duplicate feature(s) in %s", len(drop), group["genome"].iloc[0])
    return group.drop(index=list(drop)).reset_index(drop=True)


def process_group_func(group, goi_type, collapse_overlaps=True, min_overlap=0.5):
    group = group.copy()
    if collapse_overlaps:
        group = collapse_overlapping_duplicates(group, goi_type, min_overlap)
    center_rows = group[group["type"].astype(str) == str(goi_type)]
    if center_rows.empty:
        group.sort_values("start", kind="stable", inplace=True)
        group = group.reset_index(drop=True)
        group["rel_pos"] = np.arange(len(group))
        return group

    center_index = center_rows.index[0]
    center_orientation = group.loc[center_index, "orientation"]
    if center_orientation == "-":
        flipped = group["orientation"].map({"+": "-", "-": "+"})
        group["orientation"] = flipped.where(flipped.notna(), group["orientation"])
        group.loc[center_index, "orientation"] = "+"
        left = float(group["start"].min())
        right = float(group["end"].max())
        old_start = group["start"].copy()
        old_end = group["end"].copy()
        group["start"] = right - old_end
        group["end"] = right - old_start

    group = group.sort_values("start", kind="stable").reset_index(drop=True)
    goi_matches = group.index[group["type"].astype(str) == str(goi_type)]
    goi_idx = goi_matches[0] if not goi_matches.empty else 0
    group["rel_pos"] = group.index - goi_idx
    return group

def _read_name_mapping(genome_name_tsv):
    if not genome_name_tsv or str(genome_name_tsv).lower() == "null" or not os.path.isfile(genome_name_tsv):
        return None
    rename_df = pd.read_csv(genome_name_tsv, sep="\t", dtype=str, keep_default_na=False)
    if not {"contig", "organism_name"}.issubset(rename_df.columns):
        return None
    rename_df = rename_df[["contig", "organism_name"]].drop_duplicates(subset=["contig"])
    return rename_df.set_index("contig")["organism_name"]

def prepare_dataframe_func(input_file, goi_type, gene_name, genome_name_tsv=None,
                           collapse_overlaps=True, min_overlap=0.5):
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    raw_df = pd.read_csv(input_file, sep="\t", header=None, dtype=str, keep_default_na=False)
    if raw_df.empty:
        raise ValueError("Input file is empty.")

    compact_columns = [
        "neighbourhood_id", "feature_id", "gene", "start", "end",
        "orientation", "type", "genome_id", "cluster_id",
        "cluster_member_id", "cluster_gene_name", "cluster_color",
    ]
    is_compact_format = False
    if raw_df.shape[1] >= 12:
        start_values = pd.to_numeric(raw_df.iloc[:, 3], errors="coerce")
        end_values = pd.to_numeric(raw_df.iloc[:, 4], errors="coerce")
        color_values = raw_df.iloc[:, 11].astype(str).str.strip()
        valid_color = color_values.eq("") | color_values.str.fullmatch(r"#[0-9A-Fa-f]{6}") | color_values.str.contains(r"_h_", na=False)
        is_compact_format = bool(start_values.notna().all() and end_values.notna().all() and valid_color.all())

    name_map = _read_name_mapping(genome_name_tsv)
    if is_compact_format:
        source = raw_df.iloc[:, :12].copy()
        source.columns = compact_columns
        df = pd.DataFrame(index=source.index)
        df["genome"] = source["neighbourhood_id"].astype(str).str.strip()
        df["gene_id"] = source["gene"].astype(str).str.strip()
        df["start"] = pd.to_numeric(source["start"], errors="raise")
        df["end"] = pd.to_numeric(source["end"], errors="raise")
        df["orientation"] = source["orientation"].astype(str).str.strip()
        df["type"] = source["type"].astype(str).str.strip()
        df["cluster_color"] = source["cluster_color"].astype(str).str.strip().replace("", "#FFFFFF")
        family_name = source["cluster_gene_name"].astype(str).str.strip()
        original_gene = source["gene"].astype(str).str.strip()
        df["gene_db"] = family_name.where(family_name.ne(""), original_gene)
        df["neighbourhood_id"] = source["neighbourhood_id"].astype(str)
        df["feature_id"] = source["feature_id"].astype(str)
        df["genome_id"] = source["genome_id"].astype(str).str.strip()
        df["cluster_id"] = source["cluster_id"].astype(str).str.strip()
        df["cluster_member_id"] = source["cluster_member_id"].astype(str)
        df["cluster_gene_name"] = source["cluster_gene_name"].astype(str)
        df["contig"] = source["feature_id"].astype(str).str.split(":").str[0]

        family_id = source["cluster_id"].astype(str).str.strip()
        gene_accession = source["gene"].astype(str).str.strip()
        feature_id = source["feature_id"].astype(str).str.strip()
        fallback_id = gene_accession.where(gene_accession.ne(""), feature_id)
        df["feature_key"] = np.where(family_id.ne(""), "family:" + family_id, "gene:" + fallback_id)

        fallback_label = df["genome_id"].astype(str) + " | " + df["neighbourhood_id"].astype(str)
        df["organism_name"] = fallback_label
        if name_map is not None:
            mapped_name = df["contig"].map(name_map).fillna("").astype(str).str.strip()
            has_name = mapped_name.ne("")
            df.loc[has_name, "organism_name"] = mapped_name[has_name] + " | " + df.loc[has_name, "genome_id"].astype(str) + " | " + df.loc[has_name, "neighbourhood_id"].astype(str)
        organism_name_from_file = "yes"
    else:
        expected_cols = ["genome", "gene_id", "start", "end", "orientation", "type", "cluster_color", "gene_db"]
        extra_cols = [f"extra_{i + 1}" for i in range(raw_df.shape[1] - len(expected_cols))]
        raw_df.columns = expected_cols + extra_cols
        df = raw_df.copy()
        df["start"] = pd.to_numeric(df["start"], errors="raise")
        df["end"] = pd.to_numeric(df["end"], errors="raise")
        df["genome"] = df["genome"].astype(str).str.strip()
        df["gene_id"] = df["gene_id"].astype(str).str.strip()
        df["gene_db"] = df["gene_db"].astype(str).str.strip()
        df["cluster_color"] = df["cluster_color"].astype(str).str.strip().replace("", "#FFFFFF")
        df["contig"] = df["genome"].astype(str).str.split(":", n=1).str[0]
        annotation = df["gene_db"].astype(str).str.strip()
        accession = df["gene_id"].astype(str).str.strip()
        df["feature_key"] = np.where(annotation.ne(""), "annotation:" + annotation, "gene:" + accession)
        if name_map is not None:
            mapped_name = df["contig"].map(name_map).fillna("").astype(str).str.strip()
            has_name = mapped_name.ne("")
            df["organism_name"] = df["contig"].astype(str)
            df.loc[has_name, "organism_name"] = mapped_name[has_name] + " (" + df.loc[has_name, "contig"] + ")"
            organism_name_from_file = "yes"
        else:
            df["organism_name"] = df["contig"].astype(str)
            organism_name_from_file = "no"
        if "genome_id" not in df.columns:
            df["genome_id"] = ""

    if "original_start" not in df.columns:
        df["original_start"] = df["start"].copy()
    if "original_end" not in df.columns:
        df["original_end"] = df["end"].copy()

    df["length"] = (df["end"] - df["start"]).abs()
    df["bit_score"] = np.nan if "bit_score" not in df.columns else pd.to_numeric(df["bit_score"], errors="coerce")
    df["orientation"] = df["orientation"].apply(map_ori_func)
    df["orientation_orig"] = df["orientation"]
    df["gene"] = df["gene_db"].astype(str).replace("", np.nan).fillna(df["gene_id"]).astype(str)
    df["display_name"] = df["gene_id"].astype(str) if gene_name == "id" else df["gene_db"].astype(str).replace("", np.nan).fillna(df["gene_id"]).astype(str)

    processed = [process_group_func(group_df, goi_type, collapse_overlaps, min_overlap) for _, group_df in df.groupby("genome", sort=False, dropna=False)]
    big_df = pd.concat(processed, ignore_index=True)
    return big_df, organism_name_from_file

def save_figure(fig, base_path, svg="no", dpi=150):
    os.makedirs(os.path.dirname(base_path) or ".", exist_ok=True)
    fig.savefig(f"{base_path}.png", dpi=dpi, bbox_inches="tight")
    if str(svg).strip().lower() == "yes":
        fig.savefig(f"{base_path}.svg", bbox_inches="tight")

def make_cluster_palette(labels, noise_color="#BDBDBD"):
    """One distinct colour per cluster (golden-angle hues, alternating lightness) so that
    ~60 clusters do not recycle tab20; noise (-1) is always light grey."""
    ids = sorted({int(x) for x in np.asarray(labels).astype(int) if int(x) != -1})
    palette = {}
    for i, c in enumerate(ids):
        rgb = colorsys.hls_to_rgb((i * 0.618033988749895) % 1.0, 0.38 + 0.18 * (i % 2), 0.70)
        palette[str(c)] = mcolors.to_hex(rgb)
    palette["-1"] = noise_color
    return palette


def add_nt_scale_bar(fig, length_nt=10_000, margin_in=0.75):
    if not getattr(fig, "axes", None) or length_nt <= 0:
        return
    fig_h = fig.get_figheight()
    frac = margin_in / fig_h                       # same absolute margin for 3 or 15 rows
    try:
        fig.subplots_adjust(bottom=frac)
    except Exception:
        pass
    ruler_ax = fig.add_axes([0.30, 0.1 * frac, 0.40, 0.55 * frac])
    ruler_ax.set_xlim(0, length_nt)
    ruler_ax.set_ylim(0, 1)
    ruler_ax.set_yticks([])
    tick_values = np.linspace(0, length_nt, 5)     # 0, 2.5k, 5k, 7.5k, 10k: no crowded 0/1,000
    ruler_ax.set_xticks(tick_values)
    ruler_ax.set_xticklabels([f"{int(v):,}" for v in tick_values], fontsize=8)
    ruler_ax.set_xlabel("Nucleotide scale (nt)", fontsize=8, labelpad=1)
    ruler_ax.tick_params(axis="x", direction="out", length=4, pad=1)
    ruler_ax.tick_params(axis="y", left=False, labelleft=False)
    for spine in ("left", "right", "top"):
        ruler_ax.spines[spine].set_visible(False)
    ruler_ax.spines["bottom"].set_position(("data", 0.35))
    ruler_ax.patch.set_alpha(0)

def network_graph_vis(clusters, X_dense, output_path, folder, title, svg="no", cluster_palette=None, similarity_threshold=0.5, max_nodes=500):
    if nx is None:
        logger.warning("networkx is not installed; skipping network graph visualization for %s.", folder)
        return
    out_dir = os.path.join(output_path, folder)
    os.makedirs(out_dir, exist_ok=True)
    
    n_samples = X_dense.shape[0]
    if n_samples > max_nodes:
        logger.warning("Too many samples (%d) for network graph (max %d); skipping to avoid visual clutter.", n_samples, max_nodes)
        return
        
    # Compute pair-wise similarities
    S = cosine_similarity(X_dense)
    
    # Build graph
    G = nx.Graph()
    for i in range(n_samples):
        G.add_node(i, cluster=str(clusters[i]))
        
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            if S[i, j] >= similarity_threshold:
                G.add_edge(i, j, weight=float(S[i, j]))
                
    if G.number_of_edges() == 0:
        logger.warning("No edges passed the similarity threshold (%.2f) for network graph in %s.", similarity_threshold, folder)
        return
        
    fig, ax = plt.subplots(figsize=(10, 10))
    try:
        # Compute force-directed spring layout
        pos = nx.spring_layout(G, weight='weight', seed=42)
        if cluster_palette is None:
            cluster_palette = make_cluster_palette(clusters)
            
        node_colors = [cluster_palette.get(str(clusters[i]), "#BDBDBD") for i in G.nodes()]
        
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=40, alpha=0.8, ax=ax, linewidths=0.2, edgecolors="white")
        nx.draw_networkx_edges(G, pos, alpha=0.15, width=0.5, edge_color="gray", ax=ax)
        
        ax.set_title(f"{title}\nForce-Directed Network Graph (similarity \\geq {similarity_threshold})")
        ax.axis("off")
        fig.tight_layout()
        save_figure(fig, os.path.join(out_dir, "network_graph"), svg)
    except Exception as exc:
        logger.warning("Network graph generation failed for %s: %s", folder, exc)
    finally:
        plt.close(fig)

def umap_vis(clusters, reduced_data, reduction_label, output_path, folder, title, svg="no", variance_model=None, colors=None, cluster_palette=None):
    colors = colors or load_colors()  
    if cluster_palette is None:
        cluster_palette = make_cluster_palette(clusters)
    out_dir = os.path.join(output_path, folder)
    os.makedirs(out_dir, exist_ok=True)
    data = np.asarray(reduced_data)
    if data.ndim != 2 or data.shape[0] != len(clusters) or data.shape[1] < 2:
        return

    columns = [f"UMAP{i + 1}" for i in range(data.shape[1])] if reduction_label.lower() == "umap" else [f"SVD{i + 1}" for i in range(data.shape[1])]
    n_pair_components = min(4, data.shape[1])
    plot_df = pd.DataFrame(data[:, :n_pair_components], columns=columns[:n_pair_components])
    plot_df["Cluster"] = np.asarray(clusters).astype(str)

    try:
        scatter_df = pd.DataFrame(data[:, :2], columns=columns[:2])
        scatter_df["Cluster"] = np.asarray(clusters).astype(str)
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.scatterplot(data=scatter_df, x=columns[0], y=columns[1], hue="Cluster", palette=cluster_palette, alpha=0.75, s=36, edgecolor="white", linewidth=0.3, ax=ax)
        handles, labels = ax.get_legend_handles_labels()
        if labels and labels[0] == "Cluster":
            handles, labels = handles[1:], labels[1:]
        if handles:
            ax.legend(handles[:MAX_LEGEND_ITEMS], labels[:MAX_LEGEND_ITEMS], title="Clusters", bbox_to_anchor=(1.02, 1), loc="upper left")
        ax.set_title(f"{title}\nFirst two components")
        fig.tight_layout()
        save_figure(fig, os.path.join(out_dir, f"{reduction_label}_2D"), svg)
        plt.close(fig)
    except Exception as exc:
        logger.warning("2D scatter plot failed for %s: %s", folder, exc)

def sim_to_height(sim, metric="cosine"):
    """Mean similarity -> dendrogram height (cosine distance, or euclidean on L2-normalised vectors)."""
    d = 1.0 - float(sim)
    return float(np.sqrt(2.0 * d)) if metric == "euclidean" else d

def height_grid(metric="cosine"):
    """Legacy fixed grid (cosine distance 0.30-0.95), kept for compatibility."""
    grid = np.arange(0.30, 0.9501, 0.025)
    return np.sqrt(2.0 * grid) if metric == "euclidean" else grid

def candidate_heights(linkage_matrix, metric="cosine", linkage_method="average", max_candidates=250):
    """Cut candidates = the actual merge heights of the dendrogram (no coarse grid)."""
    h = np.unique(np.round(linkage_matrix[:, 2], 3))
    h = h[h > 0]
    if linkage_method != "ward":
        lo, hi = sim_to_height(0.70, metric), sim_to_height(0.05, metric)
        sel = h[(h >= lo) & (h <= hi)]
        h = sel if len(sel) else h
    if len(h) > max_candidates:
        h = np.unique(np.round(np.quantile(h, np.linspace(0, 1, max_candidates)), 3))
    return h

def require_signature(labels, M, min_shared=DEFAULT_MIN_SHARED_FAMILIES,
                      frac=DEFAULT_SIGNATURE_FRACTION, noise=-1):
    """Demote clusters that do not share a gene signature to noise.

    M is the binary neighbourhood x family matrix (rows aligned with labels). A cluster of >= 2
    members is kept only if >= min_shared families are present in >= frac of its members.
    """
    labels = np.asarray(labels).copy()
    if M is None or min_shared is None or min_shared <= 0:
        return labels
    Mb = np.asarray(M) > 0
    for c in np.unique(labels[labels != noise]):
        idx = labels == c
        if idx.sum() < 2:
            continue
        if int((Mb[idx].mean(axis=0) >= frac).sum()) < min_shared:
            labels[idx] = noise
    return labels

def cut_labels(linkage_matrix, height, min_cluster_size, M=None,
               min_shared=0, signature_fraction=DEFAULT_SIGNATURE_FRACTION):
    """Cut the tree, send too-small clusters to -1, then clusters without a gene signature to -1."""
    lab = fcluster(linkage_matrix, t=float(height), criterion="distance").astype(int)
    sizes = pd.Series(lab).value_counts()
    lab = np.where(np.isin(lab, sizes[sizes < min_cluster_size].index), -1, lab)
    return require_signature(lab, M, min_shared, signature_fraction)

def auto_cut_height(linkage_matrix, X_dense, min_cluster_size=DEFAULT_MIN_CLUSTER_SIZE,
                    min_coverage=DEFAULT_MIN_COVERAGE, min_cluster_cosine=DEFAULT_MIN_CLUSTER_COSINE,
                    heights=None, log=True, M=None, min_shared=0,
                    signature_fraction=DEFAULT_SIGNATURE_FRACTION,
                    metric="cosine", linkage_method="average", cohesion_stat="worst"):
    """Choose the dendrogram cut automatically.

    Every candidate height (all merge heights) is evaluated AFTER the same filters that the final
    labels get (minimum size, gene signature). coverage = fraction of neighbourhoods in kept
    clusters; coherence = mean pairwise cosine similarity inside each kept cluster.

    1. coverage >= min_coverage AND worst cluster >= min_cluster_cosine -> most coherent such cut.
    2. otherwise coherence wins: among cuts whose worst cluster >= min_cluster_cosine take the
       one with the highest coverage (previously coverage x coherence was maximised, which
       silently chose loose cuts when the two targets conflicted).
    3. otherwise the best coverage x coherence.
    Returns (height, scan DataFrame with a 'selected' column).
    """
    if heights is None:
        heights = candidate_heights(linkage_matrix, metric, linkage_method)
    S = cosine_similarity(X_dense)
    n = X_dense.shape[0]
    rows = []
    for h in heights:
        lab = cut_labels(linkage_matrix, h, min_cluster_size, M, min_shared, signature_fraction)
        ids = [c for c in np.unique(lab) if c != -1]
        row = dict(cut_height=round(float(h), 3), n_clusters_ge_min=len(ids),
                   fraction_in_clusters=0.0, largest_cluster=0, weighted_mean_cosine=np.nan,
                   median_cluster_cosine=np.nan, worst_cluster_cosine=np.nan)
        if ids:
            sizes = np.array([(lab == c).sum() for c in ids])
            coh = []
            for c in ids:
                idx = np.where(lab == c)[0]
                coh.append(float(S[np.ix_(idx, idx)][np.triu_indices(len(idx), 1)].mean()) if len(idx) > 1 else 1.0)
            row.update(fraction_in_clusters=float(sizes.sum() / n), largest_cluster=int(sizes.max()),
                       weighted_mean_cosine=float(np.average(coh, weights=sizes)),
                       median_cluster_cosine=float(np.median(coh)), worst_cluster_cosine=float(min(coh)))
        rows.append(row)
    scan = pd.DataFrame(rows)
    valid = scan.dropna(subset=["weighted_mean_cosine"])
    if valid.empty:
        scan["selected"] = False
        logger.warning("No cut produced any cluster; using the median candidate height.")
        return float(np.median(heights)), scan
    ok_cov = valid.fraction_in_clusters >= min_coverage
    _col = {"worst": "worst_cluster_cosine", "median": "median_cluster_cosine", "mean": "weighted_mean_cosine"}[cohesion_stat]
    ok_coh = valid[_col] >= min_cluster_cosine
    if (ok_cov & ok_coh).any():
        best = valid[ok_cov & ok_coh].sort_values(["weighted_mean_cosine", "cut_height"], ascending=[False, True]).iloc[0]
        reason, warn = "meets coverage and coherence targets", False
    elif ok_coh.any():
        best = valid[ok_coh].sort_values(["fraction_in_clusters", "cut_height"], ascending=[False, True]).iloc[0]
        reason, warn = "coverage target not reachable at the requested coherence (%s cluster cosine >= %.2f); maximised coverage subject to coherence" % (cohesion_stat, min_cluster_cosine), True
    else:
        best = valid.assign(_score=valid.fraction_in_clusters * valid.weighted_mean_cosine).sort_values("_score", ascending=False).iloc[0]
        reason, warn = "no cut met the coherence target; using best coverage x coherence", True
    if log:
        (logger.warning if warn else logger.info)(
            "Automatic cut height %.3f (%s): %.0f%% of neighbourhoods in %d clusters, mean cosine %.2f, median %.2f, worst cluster %.2f.",
            best.cut_height, reason, 100 * best.fraction_in_clusters, int(best.n_clusters_ge_min),
            best.weighted_mean_cosine, best.median_cluster_cosine, best.worst_cluster_cosine)
    scan["selected"] = np.isclose(scan["cut_height"], best.cut_height)
    return float(best.cut_height), scan

def recluster_overflow(labels, X_dense, active_samples, df,
                       overflow_id=9999, min_size=5, min_sim=0.5, start_id=10000):
    labels = np.asarray(labels).copy()
    idx = np.where(labels == overflow_id)[0]
    if len(idx) < max(2, min_size):
        return labels
    sub = X_dense[idx]
    Z = sch.linkage(pdist(sub, metric="cosine"), method="average")
    sub_lab = fcluster(Z, t=1 - min_sim, criterion="distance")
    g = df.groupby("genome")["genome_id"].first()
    genome_of = [g.get(active_samples[i], active_samples[i]) for i in idx]
    next_id = start_id
    for c in np.unique(sub_lab):
        members = np.where(sub_lab == c)[0]
        # keep one neighbourhood per genome (the one closest to the cluster centroid)
        cen = sub[members].mean(axis=0, keepdims=True)
        sims = cosine_similarity(sub[members], cen).ravel()
        best = {}
        for m, s in zip(members, sims):
            gg = genome_of[m]
            if gg not in best or s > sims[list(members).index(best[gg])]:
                best[gg] = m
        keep = list(best.values())
        if len(keep) >= min_size:
            labels[idx[keep]] = next_id
            next_id += 1
    return labels

def assign_leftovers(labels, X_dense, min_similarity=DEFAULT_ASSIGN_MIN_SIMILARITY, noise_labels=(-1,),
                     M=None, min_shared=0, signature_fraction=DEFAULT_SIGNATURE_FRACTION):
    """Attach noise points to the most similar cluster (mean cosine similarity to its core members).

    Decisions use only the original core members, so the result does not depend on point order.
    If M (binary neighbourhood x family matrix) is given, a leftover may only join a cluster whose
    gene signature (families in >= signature_fraction of the core members) it shares by at least
    min(min_shared, size of that signature) families.
    Returns (new_labels, status array of 'core'/'attached'/'unassigned').
    """
    labels = np.asarray(labels).copy()
    status = np.where(np.isin(labels, noise_labels), "unassigned", "core").astype(object)
    clusters = [c for c in np.unique(labels) if c not in noise_labels]
    if not clusters:
        return labels, status
    S = cosine_similarity(X_dense)
    core_labels = labels.copy()
    Mb = None if (M is None or not min_shared or min_shared <= 0) else (np.asarray(M) > 0)
    signature = {}
    if Mb is not None:
        for c in clusters:
            signature[c] = np.where(Mb[core_labels == c].mean(axis=0) >= signature_fraction)[0]
    for i in np.where(np.isin(core_labels, noise_labels))[0]:
        sims = {}
        for c in clusters:
            sig = signature.get(c)
            if sig is not None and len(sig) > 0 and Mb[i, sig].sum() < min(min_shared, len(sig)):
                continue
            sims[c] = S[i, core_labels == c].mean()
        if not sims:
            continue
        best = max(sims, key=sims.get)
        if sims[best] >= min_similarity:
            labels[i] = best
            status[i] = "attached"
    return labels, status


def rare_support(M, df, unit="genome"):
    """Per-family support used for the rare-gene filter: number of distinct genomes (or neighbourhoods)."""
    if unit != "genome" or "genome_id" not in df.columns:
        return M.sum(axis=0)
    gmap = df.drop_duplicates("genome").set_index("genome")["genome_id"].astype(str).str.strip()
    nb = pd.Series(M.index.astype(str), index=M.index)
    g = nb.map(lambda x: gmap.get(x, "")).where(lambda v: v != "", nb)
    return M.groupby(g.values).max().sum(axis=0)


def assignment_table(labels, status, X_dense, sample_ids, method):
    """Per-neighbourhood confidence: mean cosine to its own cluster (excluding itself) and to the best other cluster."""
    labels = np.asarray(labels)
    S = cosine_similarity(X_dense)
    rows = []
    clusters = [c for c in np.unique(labels) if c != -1]
    for i, sid in enumerate(sample_ids):
        own = np.nan
        if labels[i] != -1:
            members = np.where((labels == labels[i]) & (np.arange(len(labels)) != i))[0]
            own = float(S[i, members].mean()) if len(members) else np.nan
        others = [S[i, labels == c].mean() for c in clusters if c != labels[i]]
        rows.append(dict(neighbourhood_id=str(sid), method=method, cluster=int(labels[i]), status=status[i],
                         mean_cosine_to_own_cluster=own,
                         best_other_cluster_cosine=float(max(others)) if others else np.nan))
    return pd.DataFrame(rows)


def optimize_hdbscan(X, min_cluster_size=5, metric="cosine", max_candidate_size=10, max_silhouette_samples=1000,
                     noise_penalty=0.25, min_samples=1, selection_method="eom"):
    """HDBSCAN with min_cluster_size chosen by (silhouette - noise_penalty * noise fraction).

    The pairwise distance matrix is computed ONCE and passed as metric="precomputed" (the old code
    recomputed it for every candidate size). Returns int labels, -1 = noise; never raises.
    """
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=np.float64)
    n_samples = X.shape[0]
    noise = np.full(n_samples, -1, dtype=int)
    min_cluster_size = max(int(min_cluster_size), 2)           # sklearn requires >= 2
    if n_samples < 3 or min_cluster_size > n_samples:
        return noise
    min_samples = max(1, min(int(min_samples), n_samples))     # sklearn requires <= n_samples

    D = pairwise_distances(X, metric=metric)
    np.fill_diagonal(D, 0.0)
    np.clip(D, 0.0, None, out=D)                               # cosine can give tiny negative values

    def _fit(size, single=False):
        return HDBSCAN(min_cluster_size=size, min_samples=min_samples, metric="precomputed",
                   cluster_selection_method=selection_method,
                   allow_single_cluster=single, copy=True).fit_predict(D)

    best_labels, best_score = None, -np.inf
    top = min(max(int(max_candidate_size), min_cluster_size), n_samples)
    for size in range(min_cluster_size, top + 1):
        try:
            labels = _fit(size)
        except Exception as exc:
            logger.debug("HDBSCAN failed for min_cluster_size=%d: %s", size, exc)
            continue
        inlier = labels != -1
        if len(np.unique(labels[inlier])) < 2:
            continue                                           # silhouette needs >= 2 clusters
        try:
            idx = np.flatnonzero(inlier)
            score = silhouette_score(D[np.ix_(idx, idx)], labels[idx], metric="precomputed",
                                     sample_size=min(max_silhouette_samples, len(idx)), random_state=42)
            adjusted = float(score - noise_penalty * (1.0 - inlier.mean()))
        except Exception as exc:
            logger.debug("silhouette failed for min_cluster_size=%d: %s", size, exc)
            continue
        if adjusted > best_score:
            best_score, best_labels = adjusted, labels.copy()

    if best_labels is not None:
        return best_labels
    # no candidate gave >= 2 clusters (e.g. one homogeneous group): use the requested size as is
    try:
        labels = _fit(min_cluster_size)
        if not (labels != -1).any():                           # everything is noise: maybe ONE homogeneous group
            labels = _fit(min_cluster_size, single=True)
    except Exception as exc:
        logger.warning("HDBSCAN failed (%s); all neighbourhoods set to noise.", exc)
        return noise
    logger.info("HDBSCAN: no min_cluster_size in [%d, %d] gave >= 2 clusters; using %d (%d cluster(s)).",
                min_cluster_size, top, min_cluster_size, len(set(labels[labels != -1])))
    return labels

def resolve_paralog_conflicts(labels, active_samples, df, X_svd, strategy="reassign", sim_threshold=0.5, overflow_cluster_id=9999):
    if strategy == "none":
        return labels

    labels_series = pd.Series(labels, index=active_samples, copy=True)
    sample2genome = (
        df.groupby("genome")["genome_id"]
        .first()
        .replace("", np.nan)
        .fillna(df.groupby("genome")["contig"].first())
        .to_dict()
    )

    centroids = {}
    for cid in [
        c for c in labels_series.unique()
        if c != -1 and c != overflow_cluster_id
    ]:
        member_indices = np.where(labels_series.to_numpy() == cid)[0]
        if len(member_indices) > 0:
            centroids[cid] = X_svd[member_indices].mean(axis=0)

    for cid in centroids:
        cluster_samples = labels_series[labels_series == cid].index
        conflicting_genomes = (
            pd.Series([sample2genome.get(s, s) for s in cluster_samples])
            .value_counts()
        )
        conflicting_genomes = conflicting_genomes[
            conflicting_genomes > 1
        ].index.tolist()

        for genome in conflicting_genomes:
            genome_samples = [
                s for s in cluster_samples
                if sample2genome.get(s, s) == genome
            ]
            genome_indices = [
                active_samples.get_loc(s) for s in genome_samples
            ]
            similarities = cosine_similarity(
                X_svd[genome_indices],
                centroids[cid].reshape(1, -1),
            ).flatten()
            best_idx = int(np.argmax(similarities))

            for i, sample in enumerate(genome_samples):
                if i == best_idx:
                    continue

                sample_idx = active_samples.get_loc(sample)
                new_cluster = None

                if strategy == "reassign":
                    for alt_cid, alt_centroid in centroids.items():
                        if alt_cid == cid:
                            continue

                        alternative_genomes = {
                            sample2genome.get(s, s)
                            for s in labels_series[labels_series == alt_cid].index
                        }
                        if genome not in alternative_genomes:
                            similarity = float(
                                cosine_similarity(
                                    X_svd[sample_idx].reshape(1, -1),
                                    alt_centroid.reshape(1, -1),
                                )[0, 0]
                            )
                            if similarity >= sim_threshold:
                                new_cluster = alt_cid
                                break

                labels_series.loc[sample] = (
                    new_cluster if new_cluster is not None else overflow_cluster_id
                )

    return labels_series.to_numpy()

def save_dendrogram_safe(
    linkage_matrix,
    sample_labels,
    cluster_labels,
    cluster_colors,
    output_file,
    metric,
    svg="no",
    max_full_samples=500,
    truncated_leaves=50,
    cut_height=DEFAULT_CUT_HEIGHT,
    supercluster_height=None,
):
    n_samples = len(sample_labels)
    labels = np.asarray(cluster_labels, dtype=int)
    descendant_labels = {i: {int(labels[i])} for i in range(n_samples)}
    descendant_leaves = {i: [i] for i in range(n_samples)} if n_samples <= max_full_samples else None
    link_colors = {}
    neutral_color = "#000000"           
    cluster_top = {}                    # cluster -> (height, leaves) of its highest below-cut node
    noise_color = cluster_colors.get(-1, cluster_colors.get("-1", "#BDBDBD"))

    for merge_index, row in enumerate(linkage_matrix):
        node_id = n_samples + merge_index
        a, b = int(row[0]), int(row[1])
        child_labels = descendant_labels[a] | descendant_labels[b]
        descendant_labels[node_id] = child_labels
        if descendant_leaves is not None:
            descendant_leaves[node_id] = descendant_leaves[a] + descendant_leaves[b]
        single = row[2] <= cut_height and len(child_labels) == 1
        if single:
            lab = next(iter(child_labels))
            link_colors[node_id] = noise_color if lab == -1 else cluster_colors.get(lab, neutral_color)
            if lab != -1 and descendant_leaves is not None:
                cluster_top[lab] = (row[2], descendant_leaves[node_id])
        else:
            link_colors[node_id] = neutral_color

    fig, ax = plt.subplots(figsize=(max(12, n_samples * 0.07), 8))
    try:
        if n_samples <= max_full_samples:
            result = sch.dendrogram(
                linkage_matrix,
                labels=[str(x) for x in sample_labels],
                leaf_rotation=90,
                leaf_font_size=6,
                link_color_func=lambda n: link_colors.get(int(n), neutral_color),
                ax=ax,
            )
            for tick, leaf_index in zip(ax.get_xticklabels(), result["leaves"]):
                lab = int(labels[leaf_index])
                tick.set_color(noise_color if lab == -1 else cluster_colors.get(lab, neutral_color))
            # cluster ids above each cluster, so colours alone do not have to identify them
            pos = {leaf: k for k, leaf in enumerate(result["leaves"])}
            for lab, (height, leaves) in cluster_top.items():
                x = 5 + 10 * float(np.mean([pos[l] for l in leaves]))
                ax.text(x, height + 0.01, str(lab), ha="center", va="bottom", fontsize=6, color="black")
        else:
            sch.dendrogram(
                linkage_matrix,
                truncate_mode="lastp",
                p=truncated_leaves,
                show_leaf_counts=True,
                no_labels=True,
                link_color_func=lambda n: link_colors.get(int(n), neutral_color),
                ax=ax,
            )

        if cut_height is not None and np.isfinite(cut_height):
            ax.axhline(y=cut_height, color="red", linestyle="--", linewidth=1.5,
                       label=f"Cluster cut ({cut_height:.2f})")
        if supercluster_height is not None and np.isfinite(supercluster_height):
            ax.axhline(y=supercluster_height, color="tab:blue", linestyle=":", linewidth=1.2,
                       label=f"Supercluster cut ({supercluster_height:.2f})")
        ax.set_ylabel(f"{metric} distance (average linkage)" if metric else "distance")
        ax.legend(loc="upper right")

        fig.tight_layout()
        save_figure(fig, os.path.splitext(output_file)[0], svg)
    except Exception as exc:
        logger.warning("Dendrogram generation failed: %s", exc)
    finally:
        plt.close(fig)


def plot_hierarchical_heatmap(
    feature_matrix,
    active_samples,
    linkage_matrix,
    cluster_labels,
    palette_map,
    output_path,
    svg="no",
    max_features=DEFAULT_MAX_HEATMAP_FEATURES,
    colors=None,
):
    colors = colors or load_colors()
    data = feature_matrix.loc[active_samples].copy()
    data = data.loc[:, data.nunique(axis=0) > 1]
    if data.shape[1] == 0:
        return
    if data.shape[1] > max_features:
        data = data.loc[
            :,
            data.var(axis=0).sort_values(ascending=False)
            .head(max_features).index,
        ]
    ordered_samples = active_samples[sch.leaves_list(linkage_matrix)]
    data = data.loc[ordered_samples]
    row_colors = (
        pd.Series(cluster_labels, index=active_samples)
        .map(palette_map)
        .loc[ordered_samples]
    )
    binary_cmap = mcolors.ListedColormap([
        colors["heatmap_absent"],
        colors["heatmap_present"],
    ])
    try:
        grid = sns.clustermap(
            data,
            row_cluster=False,
            col_cluster=(data.shape[1] > 1),
            row_colors=row_colors,
            cmap=binary_cmap,
            vmin=0,
            vmax=1,
            figsize=(14, 10),
            cbar_pos=None,
            yticklabels=False,
            xticklabels=False,
        )
        grid.savefig(
            os.path.join(output_path, "hierarchical_clustermap.png"),
            dpi=150,
            bbox_inches="tight",
        )
        plt.close(grid.fig)
    except Exception:
        pass

def build_biological_feature_matrix(
    df,
    eligible_families,
    feature_mode="combined",
    position_weight=0.1,
    synteny_weight=0.3,
):
    if feature_mode not in {"content", "combined", "synteny"}:
        raise ValueError(f"Unknown feature_mode: {feature_mode}")
    eligible = set(map(str, eligible_families))
    row_features = {}
    for sample_id, group in df.groupby("genome", sort=False, dropna=False):
        group = (
            group.sort_values("rel_pos", kind="stable")
            if "rel_pos" in group.columns
            else group.sort_values("start", kind="stable")
        )
        records = []
        tokens = {}
        for row in group.itertuples(index=False):
            family = str(getattr(row, "feature_key"))
            colour = str(getattr(row, "cluster_color", "")).strip().upper()
            if family not in eligible or colour == "#FFFFFF":
                records.append(None)          # uncoloured genes carry no feature and break strict adjacency
                continue
            orientation = str(getattr(row, "orientation", "."))
            rel_pos = int(getattr(row, "rel_pos", 0))
            records.append((family, orientation, rel_pos))
            if feature_mode != "synteny":
                gene_token = f"GENE::{family}"
                tokens[gene_token] = max(tokens.get(gene_token, 0.0), 1.0)
            if feature_mode in {"combined", "synteny"}:
                pos_token = f"POSITION::{rel_pos}::{family}::{orientation}"
                tokens[pos_token] = max(
                    tokens.get(pos_token, 0.0),
                    float(position_weight),
                )
        if feature_mode in {"combined", "synteny"}:
            for left, right in zip(records[:-1], records[1:]):
                if left is None or right is None:
                    continue
                lf, lo, _ = left
                rf, ro, _ = right
                pair_token = (
                    f"ADJACENCY::{lf}::{lo}>{rf}::{ro}"
                )
                tokens[pair_token] = (
                    tokens.get(pair_token, 0.0) + float(synteny_weight)
                )
        row_features[str(sample_id)] = tokens

    matrix = pd.DataFrame.from_dict(
        row_features,
        orient="index",
    ).fillna(0.0)
    matrix.index = matrix.index.astype(str)
    matrix = matrix.sort_index(axis=1)
    return matrix.astype(np.float32)

def describe_cluster_taxonomy(cluster_df):
    sample_df = cluster_df.drop_duplicates(subset=["genome"]).copy()
    if "organism_name" in sample_df.columns:
        names = (
            sample_df["organism_name"]
            .fillna("")
            .astype(str)
            .str.split("|", n=1)
            .str[0]
            .str.strip()
        )
        plausible = names.ne("") & ~names.str.match(
            r"^(GCF|GCA|NZ_|NC_|CP_)",
            case=False,
        )
        genera = names[plausible].str.split().str[0].replace("", np.nan).dropna()
    else:
        genera = pd.Series(dtype=str)
    if genera.empty:
        return {
            "n_genera": 0,
            "dominant_genus": "",
            "dominant_genus_fraction": np.nan,
        }
    counts = genera.value_counts()
    return {
        "n_genera": int(counts.size),
        "dominant_genus": str(counts.index[0]),
        "dominant_genus_fraction": float(counts.iloc[0] / counts.sum()),
    }

def cluster_genomes_func(
    df,
    output_path,
    metric="cosine",
    remove_core="yes",
    remove_rare="yes",
    core_threshold=CORE_GENE_THRESHOLD,
    rare_min_count=RARE_GENE_MIN_COUNT,
    svg="no",
    cut_height=None,
    min_cluster_size=DEFAULT_MIN_CLUSTER_SIZE,
    max_hierarchical_samples=DEFAULT_MAX_HIERARCHICAL_SAMPLES,
    max_heatmap_features=DEFAULT_MAX_HEATMAP_FEATURES,
    max_hdbscan_candidates=10,
    max_silhouette_samples=1000,
    min_gene_prevalence=None,
    paralog_resolution="none",
    paralog_sim_threshold=0.5,
    feature_mode="combined",
    position_weight=0.1,
    synteny_weight=0.3,
    colors=None,
    min_samples=1,
    noise_penalty=0.25,
    keep_small_hier_clusters="no",
    min_coverage=DEFAULT_MIN_COVERAGE,
    min_cluster_cosine=DEFAULT_MIN_CLUSTER_COSINE,
    assign_leftover_neighbourhoods="yes",
    assign_min_similarity=DEFAULT_ASSIGN_MIN_SIMILARITY,
    density_clustering="yes",
    hdbscan_selection="eom",
    adjacency_across_gaps="yes",
    linkage_method="average",
    cohesion_stat="worst",
    min_shared_families=DEFAULT_MIN_SHARED_FAMILIES,
    signature_fraction=DEFAULT_SIGNATURE_FRACTION,
    min_features_per_neighbourhood=DEFAULT_MIN_FEATURES,
    rare_unit="genome",
    supercluster_similarity=DEFAULT_SUPERCLUSTER_SIM,
):
    colors = colors or load_colors()
    os.makedirs(output_path, exist_ok=True)
    run_info = {"methods": [], "cut_height_mode": None, "cut_height": None}
    if linkage_method == "ward" and metric != "euclidean":
        raise ValueError("--linkage_method ward is only valid with --metric euclidean.")
    filtered_df = (
        df[
            df["cluster_color"].astype(str).str.strip().str.upper()
            != "#FFFFFF"
        ].copy()
        if "cluster_color" in df.columns
        else df.copy()
    )

    feature_matrix = (
        filtered_df.pivot_table(
            index="genome",
            columns="feature_key",
            aggfunc="size",
            fill_value=0,
        )
        .clip(upper=1)
        .astype(np.uint8)
    )
    analysis_matrix = feature_matrix.copy()

    if min_gene_prevalence is not None:
        gene_means = feature_matrix.mean(axis=0)
        val = float(min_gene_prevalence)
        threshold_val = (
            gene_means.quantile(1.0 - val)
            if 0.0 < val <= 1.0
            else val
        )
        analysis_matrix = feature_matrix.loc[
            :,
            gene_means >= threshold_val,
        ].copy()
    else:
        if str(remove_core).lower() == "yes":
            analysis_matrix = analysis_matrix.loc[
                :,
                analysis_matrix.mean(axis=0) < core_threshold,
            ]
        if str(remove_rare).lower() == "yes":
            analysis_matrix = analysis_matrix.loc[
                :,
                rare_support(analysis_matrix, df, rare_unit) >= rare_min_count,
            ]

    analysis_matrix = analysis_matrix.loc[
        analysis_matrix.sum(axis=1) > 0
    ].copy()

    if min_features_per_neighbourhood and min_features_per_neighbourhood > 1:
        n_feat = (analysis_matrix > 0).sum(axis=1)
        weak = n_feat[n_feat < min_features_per_neighbourhood]
        if len(weak):
            logger.warning("%d neighbourhood(s) have < %d eligible families (short contig / little context) and are not clustered; see excluded_neighbourhoods.tsv.",
                           len(weak), min_features_per_neighbourhood)
            weak.rename("n_eligible_families").rename_axis("neighbourhood_id").reset_index().to_csv(
                os.path.join(output_path, "excluded_neighbourhoods.tsv"), sep="\t", index=False)
            analysis_matrix = analysis_matrix.loc[n_feat >= min_features_per_neighbourhood].copy()

    if analysis_matrix.empty:
        raise ValueError(
            "No neighbourhoods/features remain after feature filtering. "
            "Relax core/rare filtering or check input colours."
        )

    eligible_families = analysis_matrix.columns.astype(str)
    # filtered_df has the uncoloured genes removed, so pairs are formed across the gaps they leave;
    # the full table keeps them (as non-features) so that an uncoloured gene breaks the pair.
    feature_source = filtered_df if str(adjacency_across_gaps).lower() == "yes" else df
    clustering_matrix = build_biological_feature_matrix(
        feature_source,
        eligible_families,
        feature_mode=feature_mode,
        position_weight=position_weight,
        synteny_weight=synteny_weight,
    )
    active_samples = analysis_matrix.index.intersection(
        clustering_matrix.index,
        sort=False,
    )
    analysis_matrix = analysis_matrix.loc[active_samples].copy()
    clustering_matrix = clustering_matrix.loc[active_samples].copy()
    clustering_matrix = clustering_matrix.loc[
        :,
        clustering_matrix.sum(axis=0) > 0,
    ]

    if clustering_matrix.shape[1] == 0:
        raise ValueError(
            "The selected neighbourhood representation produced no features."
        )

    n_samples = clustering_matrix.shape[0]
    if n_samples < 2:
        raise ValueError(
            f"At least two neighbourhoods are required for clustering; found {n_samples}."
        )

    X_clustering = TfidfTransformer(
        norm="l2",
        use_idf=True,
        smooth_idf=True,
    ).fit_transform(clustering_matrix)
    X_dense = X_clustering.toarray().astype(np.float64, copy=False)

    visual_label = "UMAP"
    if n_samples >= 10 and umap is not None:
        X_visual = umap.UMAP(
            n_neighbors=min(15, n_samples - 1),
            min_dist=0.1,
            metric=metric,
            n_components=2,
        ).fit_transform(X_clustering)
    else:
        # too few samples for UMAP (or umap missing): project onto the first two SVD components
        visual_label = "SVD"
        n_comp = min(2, X_dense.shape[0] - 1, X_dense.shape[1] - 1)
        X_visual = (
            TruncatedSVD(n_components=n_comp, random_state=42).fit_transform(X_dense)
            if n_comp >= 2 else X_dense[:, :1]   # umap_vis skips plots with < 2 components
        )

    clustering_results = {}
    assignment_tables = []
    hierarchical_cluster_palette = None

    _dens = str(density_clustering).lower()
    if _dens != "only" and n_samples <= max_hierarchical_samples:
        condensed_distances = pdist(X_dense, metric=metric)
        _linkage = fastcluster.linkage if fastcluster is not None else sch.linkage
        row_linkage = _linkage(condensed_distances, method=linkage_method)

        M_bin = analysis_matrix.to_numpy() > 0          # rows aligned with X_dense / active_samples
        _sig = dict(M=M_bin, min_shared=min_shared_families, signature_fraction=signature_fraction)

        _auto_h, _scan = auto_cut_height(
            row_linkage, X_dense, min_cluster_size, min_coverage, min_cluster_cosine,
            log=(cut_height is None or float(cut_height) <= 0),
            metric=metric, linkage_method=linkage_method, cohesion_stat=cohesion_stat, **_sig,
        )
        if cut_height is None or float(cut_height) <= 0:
            cut_height = _auto_h
            run_info["cut_height_mode"] = "auto"
        else:
            cut_height = float(cut_height)
            run_info["cut_height_mode"] = "fixed"
            _scan["selected"] = False
        run_info["cut_height"] = float(cut_height)
        _scan.to_csv(os.path.join(output_path, "hierarchical", "hierarchical_cut_scan.csv"), index=False)
        _eff_min = 1 if str(keep_small_hier_clusters).lower() == "yes" else min_cluster_size
        h_labels = cut_labels(row_linkage, cut_height, _eff_min, **_sig)
        logger.info("Hierarchical: %d neighbourhoods are in no cluster after the size/signature filters.",
                    int((h_labels == -1).sum()))

        # optional looser second level
        super_height = None
        if supercluster_similarity and supercluster_similarity > 0:
            super_height = sim_to_height(supercluster_similarity, metric)
            s_lab = cut_labels(row_linkage, super_height, 2)
            run_info["supercluster"] = pd.Series(s_lab, index=active_samples.astype(str))

        h_status = np.where(h_labels == -1, "unassigned", "core").astype(object)
        if str(assign_leftover_neighbourhoods).lower() == "yes" and (h_labels == -1).any():
            h_labels, h_status = assign_leftovers(h_labels, X_dense, assign_min_similarity, **_sig)
            logger.info("Hierarchical: attached %d leftover neighbourhood(s) to existing clusters; %d remain unassigned.",
                        int((h_status == "attached").sum()), int((h_labels == -1).sum()))

        h_labels_before_paralog = h_labels.copy()
        if paralog_resolution != "none":
            logger.warning(
                "Paralog reassignment is a post-processing operation and can "
                "alter biological clusters; applying requested strategy %s.",
                paralog_resolution,
            )
            h_labels = resolve_paralog_conflicts(
                h_labels,
                active_samples,
                df,
                analysis_matrix.to_numpy(dtype=np.float32),
                paralog_resolution,
                paralog_sim_threshold,
            )
            
            h_labels = recluster_overflow(h_labels, X_dense, active_samples, df, min_size=5, min_sim=0.5)

        clustering_results["hierarchical"] = pd.Series(
            h_labels,
            index=active_samples,
            name="cluster",
        )
        hierarchical_cluster_palette = make_cluster_palette(h_labels)
        palette_map = {
            int(k): v
            for k, v in hierarchical_cluster_palette.items()
        }

        save_dendrogram_safe(
            row_linkage,
            active_samples.astype(str),
            h_labels,
            palette_map,
            os.path.join(output_path, "hierarchical", "dendrogram_hierarchical"), 
            metric,
            svg,
            cut_height=float(cut_height),
            supercluster_height=super_height,
        )
        h_status = np.where(h_labels != h_labels_before_paralog, "paralog_reassigned", h_status).astype(object)
        assignment_tables.append(assignment_table(h_labels, h_status, X_dense, active_samples, "hierarchical"))
        run_info["methods"].append("hierarchical_average_linkage")

        plot_hierarchical_heatmap(
            feature_matrix,
            active_samples,
            row_linkage,
            h_labels,
            palette_map,
            os.path.join(output_path, "hierarchical"),
            svg,
            max_heatmap_features,
            colors,
        )
    elif _dens != "only":
        logger.warning(
            "Skipping hierarchical clustering: %d samples exceeds limit %d; running HDBSCAN instead.",
            n_samples,
            max_hierarchical_samples,
        )

    run_density = _dens in ("yes", "only") or "hierarchical" not in clustering_results
    if run_density:
        d_labels = optimize_hdbscan(
            X_clustering,
            min_cluster_size,
            metric,
            max_hdbscan_candidates,
            max_silhouette_samples,
            noise_penalty=noise_penalty,
            min_samples=min_samples,
            selection_method=hdbscan_selection,
        )
        logger.info("HDBSCAN: %d cluster(s), %d neighbourhood(s) flagged as noise.",
                    len(set(d_labels[d_labels != -1])), int((d_labels == -1).sum()))
        d_status = np.where(d_labels == -1, "unassigned", "core").astype(object)
        if str(assign_leftover_neighbourhoods).lower() == "yes" and (d_labels == -1).any() and (d_labels != -1).any():
            d_labels, d_status = assign_leftovers(d_labels, X_dense, assign_min_similarity,
                                                  M=analysis_matrix.to_numpy() > 0, min_shared=min_shared_families,
                                                  signature_fraction=signature_fraction)
        d_labels_before_paralog = d_labels.copy()
        if paralog_resolution != "none":
            # same rule as for the hierarchical labels: one neighbourhood per genome per cluster
            d_labels = resolve_paralog_conflicts(
                d_labels,
                active_samples,
                df,
                analysis_matrix.to_numpy(dtype=np.float32),
                paralog_resolution,
                paralog_sim_threshold,
            )
            d_status = np.where(d_labels != d_labels_before_paralog, "paralog_reassigned", d_status).astype(object)
        clustering_results["density"] = pd.Series(
            d_labels,
            index=active_samples,
            name="cluster",
        )
        assignment_tables.append(assignment_table(d_labels, d_status, X_dense, active_samples, "density"))
        run_info["methods"].append("hdbscan")

        if "hierarchical" in clustering_results:
            h_ser = clustering_results["hierarchical"]
            ari = float(adjusted_rand_score(h_ser.to_numpy(), d_labels))
            run_info["ari_hierarchical_vs_density"] = round(ari, 3)
            logger.info("Agreement between hierarchical and HDBSCAN clusters: ARI = %.2f", ari)
            # neighbourhoods are 'consistent' when their hierarchical cluster maps onto one HDBSCAN cluster
            cross = pd.crosstab(h_ser.to_numpy(), d_labels)
            dominant = cross.idxmax(axis=1)
            agree = [(d_labels[i] == dominant.get(h_ser.iloc[i], None)) for i in range(len(h_ser))]
            pd.DataFrame({
                "neighbourhood_id": h_ser.index.astype(str),
                "hierarchical": h_ser.to_numpy(),
                "density": d_labels,
                "agree_with_majority_mapping": agree,
            }).to_csv(os.path.join(output_path, "clustering_comparison.tsv"), sep="\t", index=False)

    for method_name, labels in clustering_results.items():
        method_palette = (
            hierarchical_cluster_palette
            if method_name == "hierarchical" and hierarchical_cluster_palette is not None
            else make_cluster_palette(labels.to_numpy())
        )
        umap_vis(
            labels.to_numpy(),
            X_visual,
            visual_label,
            output_path,
            method_name,
            f"{method_name.capitalize()} clustering",
            svg,
            cluster_palette=method_palette,
        )

        network_graph_vis(
            labels.to_numpy(),
            X_dense,
            output_path,
            method_name,
            f"{method_name.capitalize()} clustering",
            svg,
            cluster_palette=method_palette,
        )

    # neighbourhoods without any eligible family were silently dropped before; report them as noise
    all_ids = pd.Index(df["genome"].astype(str).unique())
    n_dropped = len(all_ids.difference(active_samples.astype(str)))
    if n_dropped:
        logger.warning("%d neighbourhood(s) had no family passing the core/rare filters; assigned to noise (-1).", n_dropped)
    for _m in list(clustering_results):
        s_ = clustering_results[_m]; s_.index = s_.index.astype(str)
        clustering_results[_m] = s_.reindex(all_ids, fill_value=-1).astype(int)

    if assignment_tables:
        pd.concat(assignment_tables, ignore_index=True).to_csv(
            os.path.join(output_path, "neighbourhood_assignments.tsv"), sep="\t", index=False)
    if "supercluster" in run_info:
        _sup = run_info["supercluster"]
        _sup.index = _sup.index.astype(str)
        run_info["supercluster"] = _sup.reindex(all_ids, fill_value=-1).astype(int)
    run_info["n_neighbourhoods"] = int(len(all_ids))
    run_info["n_unassigned"] = {m: int((l == -1).sum()) for m, l in clustering_results.items()}
    run_info["n_clusters"] = {m: int(len(set(l[l != -1]))) for m, l in clustering_results.items()}
    return clustering_results, feature_matrix, analysis_matrix, run_info

def plot_cluster_func(cluster_df, goi_type, output_path, scale, gene_lable, folder, organism_name_from_file, svg, max_subplots=15):
    if Feature is None or FeatureList is None:
        return []
    cluster_df = cluster_df.copy()
    final_output_dir = os.path.join(output_path, folder, "synteny_cluster_plots")
    os.makedirs(final_output_dir, exist_ok=True)
    split_data = cluster_df["cluster_color"].astype(str).str.split("_h_", expand=True)
    cluster_df["cluster_color"] = split_data[0].replace("", "#FFFFFF")
    cluster_df["hatch"] = split_data[1].fillna("").astype(str) if split_data.shape[1] > 1 else ""
    name2color = {str(gene): color for gene, color in zip(cluster_df["gene"], cluster_df["cluster_color"])}

    all_fts = FeatureList([
        Feature(row.type, start=row.start, stop=row.end, strand=row.orientation, meta={"seqid": row.contig, "hit": row.genome, "organism_name": row.organism_name, "logic_name": row.gene, "label_name": getattr(row, "display_name", row.gene), "_ftsviewer_hatch": row.hatch})
        for row in cluster_df.itertuples(index=False)
    ])
    scale_bool = str(scale).lower() == "yes"
    label_bool = str(gene_lable).lower() == "yes"
    groups = list(all_fts.groupby("hit").values())
    report_summary_json = []
    if not groups:
        return report_summary_json

    for fignum in range((len(groups) + max_subplots - 1) // max_subplots):
        fts = FeatureList([f for g in groups[fignum * max_subplots:(fignum + 1) * max_subplots] for f in g])
        align_features = FeatureList([f for f in fts if f.type == goi_type]) if scale_bool else None
        fig = fts.plot_ftsviewer(groupby="hit", axlabel="organism_name" if organism_name_from_file.lower() == "yes" else "seqid", label=(lambda ft: ft.meta["label_name"]) if label_bool else (lambda ft: None), align=align_features, crop=500, sharex=scale_bool, colorby=lambda ft: ft.meta["logic_name"], color=name2color, figsize=(10, 1.5 * len(groups[fignum * max_subplots:(fignum + 1) * max_subplots])), ncols=1, with_ruler=False, labels_spacing=8, fontdict={"fontsize": 7})
        add_nt_scale_bar(fig, length_nt=10_000)
        _axes = [a for a in fig.axes if a.get_visible() and a.get_position().width > 0.2]
        if _axes:   # sit just above the first row title, so bbox_inches="tight" adds no blank band
            fig.text(min(a.get_position().x0 for a in _axes),
                     max(a.get_position().y1 for a in _axes) + 0.3 / fig.get_figheight(),
                     "A", fontsize=18, fontweight="bold", va="bottom", ha="left")
        base_name = f"synteny_cluster_{int(cluster_df['cluster_label'].iloc[0])}_{fignum:03d}_A"
        save_figure(fig, os.path.join(final_output_dir, base_name), svg)
        plt.close(fig)
    return report_summary_json

def write_labels(df, goi_type, output_path, folder):
    output_dir = os.path.join(output_path, folder)
    os.makedirs(output_dir, exist_ok=True)
    center_rows = df[df["type"].astype(str) == str(goi_type)]
    new_df = center_rows[["cluster_label", "organism_name", "contig", "gene", "original_start", "original_end"]].copy()
    new_df["orientation"] = center_rows.get("orientation_orig", center_rows["orientation"])
    output_file = os.path.join(output_dir, "gene_of_interest.tsv")
    new_df.to_csv(output_file, mode="a", header=not os.path.exists(output_file), sep="\t", index=False)

def conserved_positions(cluster_df, threshold=0.5):
    """Number of GOI-relative positions where ONE family occupies >= threshold of the neighbourhoods.
    0-1 = the cluster shares genes but not a syntenic block (diverse); >=3 = a real conserved block."""
    d = cluster_df[(cluster_df["rel_pos"] != 0) & (cluster_df["cluster_color"].astype(str).str.upper() != "#FFFFFF")]
    n = cluster_df["genome"].nunique()
    if d.empty or n < 2:
        return 0, 0.0
    cnt = d.groupby(["rel_pos", "feature_key"])["genome"].nunique() / n
    best = cnt.groupby(level=0).max()
    return int((best >= threshold).sum()), float(best.mean() if len(best) else 0.0)


def summarize_cluster(cluster_label, cluster_df, feature_matrix, analysis_matrix):
    sample_ids = cluster_df["genome"].dropna().astype(str).drop_duplicates().tolist()
    full_valid_ids = [s for s in sample_ids if s in feature_matrix.index]
    n_neighbourhoods = len(full_valid_ids)
    avg_cosine_similarity, avg_euclidean_dist = np.nan, np.nan
    analysis_ids = [s for s in sample_ids if s in analysis_matrix.index]
    metric_vectors = analysis_matrix.loc[analysis_ids].to_numpy(dtype=np.float64)
    if cluster_label not in (-1, 9999) and len(metric_vectors) > 1:
        sims = cosine_similarity(metric_vectors)
        avg_cosine_similarity = float(sims[np.triu_indices_from(sims, k=1)].mean())
        avg_euclidean_dist = float(pdist(metric_vectors, metric="euclidean").mean())
    taxonomy = describe_cluster_taxonomy(cluster_df)
    genome_count = int(cluster_df.drop_duplicates(subset=["genome"]).get("genome_id", pd.Series(dtype=str)).nunique()) if "genome_id" in cluster_df.columns else n_neighbourhoods
    return {
        "cluster_label": int(cluster_label),
        "cluster_type": "paralog_overflow" if cluster_label == 9999 else ("noise" if cluster_label == -1 else ("singleton" if n_neighbourhoods == 1 else "cluster")),
        "n_neighbourhoods": int(n_neighbourhoods),
        "n_genomes": genome_count,
        "n_gene_records": int(len(cluster_df)),
        "avg_cosine_similarity": float(avg_cosine_similarity) if np.isfinite(avg_cosine_similarity) else None,
        "avg_euclidean_dist": float(avg_euclidean_dist) if np.isfinite(avg_euclidean_dist) else None,
        "n_conserved_positions": (_cp := conserved_positions(cluster_df))[0],
        "cluster_quality": ("noise" if cluster_label in (-1, 9999) else ("syntenic" if _cp[0] >= 3 else ("partial" if _cp[0] == 2 else "content_only"))),
        "mean_top_family_fraction_per_position": round(_cp[1], 3),
        **taxonomy,
    }

def main():
    try:
        args = setup_parser_func().parse_args()
        if args.cut_height_args is not None and args.cut_height_args <= 0:
            args.cut_height_args = None      # 0 / negative -> automatic cut
        if args.linkage_method == "ward" and args.metric != "euclidean":
            raise SystemExit("--linkage_method ward requires --metric euclidean.")
        if args.metric == "euclidean" and args.cut_height_args is not None:
            logger.warning("--cut_height_args is a distance in the chosen metric; with euclidean distance on "
                           "L2-normalised vectors the range is 0-1.41 (cosine distance 0.7 ~ euclidean 1.18).")

        logger.info("Reading input data...")
        df, organism_name_from_file = prepare_dataframe_func(
            args.input_file,
            args.goi_type,
            args.gene_name,
            args.name_file,
            collapse_overlaps=(args.collapse_overlaps == "yes"),
            min_overlap=args.overlap_fraction,
        )

        clustering_results, feature_matrix, analysis_matrix, run_info = cluster_genomes_func(
            df=df,
            output_path=args.output_path,
            metric=args.metric,
            remove_core=args.remove_core_genes,
            remove_rare=args.remove_rare_genes,
            core_threshold=args.core_threshold,
            rare_min_count=args.rare_min_count,
            min_cluster_size=args.min_cluster_size,
            max_hierarchical_samples=args.max_hierarchical_samples,
            max_hdbscan_candidates=args.max_hdbscan_candidates,
            max_silhouette_samples=args.max_silhouette_samples,
            min_gene_prevalence=args.min_gene_prevalence,
            paralog_resolution=args.paralog_resolution,
            paralog_sim_threshold=args.paralog_sim_threshold,
            feature_mode=args.feature_mode,
            position_weight=args.position_weight,
            synteny_weight=args.synteny_weight,
            min_samples=args.min_samples,
            noise_penalty=args.noise_penalty,
            keep_small_hier_clusters=args.keep_small_hier_clusters,
            cut_height=args.cut_height_args,
            svg=args.svg,
            min_coverage=args.min_coverage,
            min_cluster_cosine=args.min_cluster_cosine,
            assign_leftover_neighbourhoods=args.assign_leftovers,
            assign_min_similarity=args.assign_min_similarity,
            density_clustering=args.density_clustering,
            hdbscan_selection=args.hdbscan_selection,
            adjacency_across_gaps=args.adjacency_across_gaps,
            linkage_method=args.linkage_method,
            cohesion_stat=args.cohesion_stat,
            min_shared_families=args.min_shared_families,
            signature_fraction=args.signature_fraction,
            min_features_per_neighbourhood=args.min_features_per_neighbourhood,
            rare_unit=args.rare_unit,
            supercluster_similarity=args.supercluster_similarity,
        )

        for metric_name, labels in clustering_results.items():
            metric_dir = os.path.join(args.output_path, metric_name)
            os.makedirs(metric_dir, exist_ok=True)
            goi_path = os.path.join(metric_dir, "gene_of_interest.tsv")
            if os.path.exists(goi_path):
                os.remove(goi_path)

            metric_df = df.copy()
            metric_df["cluster_label"] = metric_df["genome"].map(labels)
            center_rows = metric_df[
                metric_df["type"].astype(str) == str(args.goi_type)
            ].copy()
            membership_cols = [
                c for c in [
                    "cluster_label",
                    "genome",
                    "genome_id",
                    "organism_name",
                    "contig",
                    "gene",
                    "original_start",
                    "original_end",
                    "orientation_orig",
                ]
                if c in center_rows.columns
            ]
            membership = (
                center_rows[membership_cols]
                .drop_duplicates(subset=["genome"])
                .sort_values(
                    ["cluster_label", "organism_name"],
                    kind="stable",
                )
            )
            membership = membership.rename(columns={
                "genome": "neighbourhood_id",
                "orientation_orig": "orientation",
            })
            if metric_name == "hierarchical" and "supercluster" in run_info:
                membership["supercluster"] = membership["neighbourhood_id"].astype(str).map(run_info["supercluster"]).fillna(-1).astype(int)
            membership.to_csv(
                os.path.join(metric_dir, "cluster_membership.tsv"),
                sep="\t",
                index=False,
            )

            cluster_stats = []
            usable_df = metric_df.dropna(subset=["cluster_label"]).copy()
            for cluster_label, cluster_df in usable_df.groupby(
                "cluster_label",
                sort=True,
            ):
                cluster_label = int(cluster_label)      
                write_labels(
                    cluster_df,
                    args.goi_type,
                    args.output_path,
                    metric_name,
                )
                plot_cluster_func(
                    cluster_df,
                    args.goi_type,
                    args.output_path,
                    args.scale,
                    args.gene_lable,
                    metric_name,
                    organism_name_from_file,
                    args.svg,
                )
                cluster_stats.append(
                    summarize_cluster(
                        cluster_label,
                        cluster_df,
                        feature_matrix,
                        analysis_matrix,
                    )
                )

                # Prepare required columns and directories for stats and logos
                cluster_logo_dir = os.path.join(
                    args.output_path,
                    metric_name,
                    "microsynteny_logo",
                )
                cluster_df = cluster_df.copy()
                cluster_df["gene_with_orientation"] = (
                    cluster_df["gene"].astype(str)
                    + " ("
                    + cluster_df["orientation"].astype(str)
                    + ")"
                )

                if cluster_label in (-1, 9999) and args.logo_for_noise != "yes":
                    continue      # noise / overflow are not clusters: no logo or statistics

                if args.microsynteny_logo.lower() == "yes":
                    plot_microsynteny_logo(
                        cluster_label=cluster_label,
                        clustering_algo=metric_name,
                        df=cluster_df,
                        output_dir=cluster_logo_dir,
                        goi_type=args.goi_type,
                        svg=(args.svg.lower() == "yes"),
                        colors=None,
                        genome_col="genome",
                        rel_pos_col="rel_pos",
                        family_col=(
                            "gene_orientation"
                            if "gene_orientation" in cluster_df.columns
                            else "gene_with_orientation"
                        ),
                        color_col="cluster_color",
                        max_flanking_genes=None,
                        drop_unannotated=False,      # Keep white/unannotated genes
                        min_column_occupancy=0.0,    # Do not hide low-occupancy columns
                        max_columns_per_side=None,   # Do not cut columns (default used to be 8)
                        min_family_fraction=0.0,     # Do not drop rare families
                        small_sample_correction=False,
                        suffix="B", panel_letter="B",
                    )

                # Move plot_cluster_statistics INSIDE the loop and guard against NaN/Inf axis limits
                try:
                    plot_cluster_statistics(
                        cluster_df, 
                        cluster_logo_dir, 
                        cluster_label, 
                        metric_name,
                        family_col="gene_with_orientation",
                        drop_unannotated=False,        # same settings as the logo
                        min_family_fraction=0.0,
                        small_sample_correction=False,
                        suffix="C", panel_letter="C",
                    )
                except Exception as exc:
                    logger.warning(
                        "Could not plot cluster statistics for cluster %s in %s: %s",
                        cluster_label, metric_name, exc
                    )

            pd.DataFrame(cluster_stats).to_csv(
                os.path.join(metric_dir, "cluster_stats.csv"),
                index=False,
            )
            if (
                args.neighbour_dir
                and os.path.exists(args.neighbour_dir)
                and extract_cluster_fastas is not None
            ):
                try:
                    logger.info(
                        "Extracting cluster FASTAs for %s into %s",
                        metric_name,
                        metric_dir,
                    )
                    extract_cluster_fastas(
                        input_tsv=goi_path,
                        neighbour_dir=args.neighbour_dir,
                        output_dir=metric_dir,
                    )
                except Exception as exc:
                    logger.warning(
                        "Could not extract cluster FASTAs for %s: %s",
                        metric_name,
                        exc,
                    )

        with open(
            os.path.join(args.output_path, "clustering_parameters.json"),
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump({
                "methods": run_info.get("methods", []),
                "feature_mode": args.feature_mode,
                "cut_height_mode": run_info.get("cut_height_mode"),
                "cut_height": run_info.get("cut_height", args.cut_height_args),
                "min_cluster_size": args.min_cluster_size,
                "min_coverage": args.min_coverage,
                "min_cluster_cosine": args.min_cluster_cosine,
                "assign_leftovers": args.assign_leftovers,
                "assign_min_similarity": args.assign_min_similarity,
                "cohesion_stat": args.cohesion_stat,
                "min_shared_families": args.min_shared_families,
                "signature_fraction": args.signature_fraction,
                "min_features_per_neighbourhood": args.min_features_per_neighbourhood,
                "rare_unit": args.rare_unit,
                "supercluster_similarity": args.supercluster_similarity,
                "linkage_method": args.linkage_method,
                "density_clustering": args.density_clustering,
                "hdbscan_selection": args.hdbscan_selection,
                "min_samples": args.min_samples,
                "noise_penalty": args.noise_penalty,
                "adjacency_across_gaps": args.adjacency_across_gaps,
                "ari_hierarchical_vs_density": run_info.get("ari_hierarchical_vs_density"),
                "n_clusters": run_info.get("n_clusters"),
                "n_unassigned": run_info.get("n_unassigned"),
                "metric": args.metric,
                "remove_core_genes": args.remove_core_genes,
                "core_threshold": args.core_threshold,
                "remove_rare_genes": args.remove_rare_genes,
                "rare_min_count": args.rare_min_count,
                "paralog_resolution": args.paralog_resolution,
                "position_weight": args.position_weight,
                "synteny_weight": args.synteny_weight,
            }, handle, indent=2)

        report_path = os.path.join(args.output_path, "comand.txt")
        command_line = " ".join(sys.argv)
        with open(report_path, "w", encoding="utf-8", errors="replace") as f:
            f.write(f"Command:\n{command_line}\n")

        logger.info("Processing completed successfully.")

    except Exception:
        logger.exception("Processing failed.")
        raise

if __name__ == "__main__":
    main()
