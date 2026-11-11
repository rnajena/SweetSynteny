# (C) 2024, Maria Schreiber, MIT license
"""Color neighbourhood features by cluster membership and create summary plots."""

import argparse
import colorsys
import csv
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
# Keep SVG labels as editable text in Inkscape rather than vector outlines.
matplotlib.rcParams["svg.fonttype"] = "none"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_colors import load_colors


WHITE = "#FFFFFF"

OUTPUT_COLUMNS = [
    "neighbourhood_id",
    "feature_id",
    "gene",
    "start",
    "end",
    "orientation",
    "type",
    "genome_id",
    "cluster_id",
    "cluster_member_id",
    "cluster_gene_name",
    "cluster_color",
]


def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Annotate neighbourhood features with cluster colors and summaries."
    )
    parser.add_argument(
        "--cluster_file",
        required=True,
        help="TSV containing cluster ID, member ID, and gene/family name in its first 3 columns.",
    )
    parser.add_argument(
        "--tsv_file",
        required=True,
        help="Headerless 8-column neighbourhood TSV including genome_id.",
    )
    parser.add_argument(
        "--output_file",
        required=True,
        help="Output path for the compact, colored neighbourhood TSV.",
    )
    parser.add_argument(
        "--gene_of_interest",
        required=True,
        help="Gene/feature name to highlight.",
    )
    parser.add_argument(
        "--goi_type",
        required=True,
        help="Feature type for the gene of interest, e.g. infernal.",
    )
    parser.add_argument(
        "--size_for_cluster",
        type=int,
        default=2,
        help="Minimum distinct members for a cluster to receive a color (default: 2).",
    )
    parser.add_argument(
        "--color_scheme_file",
        default=None,
        help="Optional JSON file defining plot colors.",
    )
    return parser


def validate_files(*files: str) -> None:
    """Raise a clear error if any required input file is missing."""
    for filename in files:
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"Input file not found: {filename}")


def generate_distinct_colors(
    num_colors: int,
    special_colors: set[str],
) -> list[str]:
    """Generate deterministic colors, excluding reserved special colors.

    The reserved colors are provided by the caller. They include the GOI
    color loaded from color.json and the white no-cluster color.
    """
    colors: list[str] = []
    reserved = {color.upper() for color in special_colors}
    i = 0

    # The golden-ratio hue increment spreads colors across the HSV color wheel.
    while len(colors) < num_colors:
        hue = (i * 0.618033988749895) % 1.0
        red, green, blue = colorsys.hsv_to_rgb(hue, 0.75, 0.85)

        color = "#{:02X}{:02X}{:02X}".format(
            int(red * 255),
            int(green * 255),
            int(blue * 255),
        )

        if color.upper() not in reserved:
            colors.append(color)

        i += 1

    return colors


def process_clusters(
    cluster_file: str,
    gene_of_interest: str,
    size_for_cluster: int,
    goi_color: str,
    special_colors: set[str],
) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Return color, annotation, and cluster-ID maps keyed by original member ID."""
    if size_for_cluster < 1:
        raise ValueError("--size_for_cluster must be at least 1.")

    goi_color = goi_color.upper()
    reserved_colors = {color.upper() for color in special_colors}
    reserved_colors.add(WHITE)
    reserved_colors.add(goi_color)

    cluster_map: dict[str, list[str]] = {}
    gene_name_map: dict[str, str] = {}
    member_cluster_map: dict[str, str] = {}
    goi_members: set[str] = set()

    with open(
        cluster_file,
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row_index, row in enumerate(reader, start=1):
            if len(row) < 3:
                print(
                    f"Warning: skipping cluster row {row_index} "
                    f"with fewer than 3 columns: {row}"
                )
                continue

            cluster_id, member_id, gene_name = (
                item.strip() for item in row[:3]
            )

            if not cluster_id or not member_id:
                print(
                    f"Warning: skipping cluster row {row_index} "
                    "with an empty cluster/member ID."
                )
                continue

            cluster_map.setdefault(cluster_id, []).append(member_id)
            gene_name_map[member_id] = gene_name
            member_cluster_map[member_id] = cluster_id

            if (
                gene_name.casefold() == gene_of_interest.casefold()
                or member_id.casefold() == gene_of_interest.casefold()
            ):
                goi_members.add(member_id)

    filtered_clusters = {
        cluster_id: list(dict.fromkeys(members))
        for cluster_id, members in cluster_map.items()
        if len(set(members)) >= size_for_cluster
    }

    colors = generate_distinct_colors(
        len(filtered_clusters),
        reserved_colors,
    )

    color_map: dict[str, str] = {}

    for (cluster_id, members), color in zip(
        filtered_clusters.items(),
        colors,
    ):
        for member_id in members:
            color_map[member_id] = color

    # The configured JSON color always takes precedence for GOI members.
    for member_id in goi_members:
        color_map[member_id] = goi_color

    print(f"Total clusters: {len(cluster_map)}")
    print(f"Clusters meeting size threshold: {len(filtered_clusters)}")
    print(f"Features assigned a color: {len(color_map)}")

    return color_map, gene_name_map, member_cluster_map


def read_neighbourhood_tsv(filename: str) -> pd.DataFrame:
    """Read the current 8-column TSV and validate its schema."""
    df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        dtype=str,
        keep_default_na=False,
    )

    if df.empty:
        raise ValueError(f"No neighbourhood records found in {filename}")

    if df.shape[1] != 8:
        raise ValueError(
            f"Expected an 8-column neighbourhood TSV, found "
            f"{df.shape[1]} columns in {filename}. "
            "Expected: neighbourhood_id, feature_id, gene, start, end, "
            "orientation, type, genome_id."
        )

    df.columns = [
        "neighbourhood_id",
        "feature_id",
        "gene",
        "start",
        "end",
        "orientation",
        "type",
        "genome_id",
    ]

    return df


def match_cluster_members(
    feature_ids: pd.Series,
    member_cluster_map: dict[str, str],
) -> dict[str, str]:
    """Match stable neighbourhood IDs to longer cluster member IDs.

    The longer cluster member ID typically starts with the full feature ID
    and then appends a suffix. Prefix matching uses an underscore boundary
    to avoid accidental partial matches.
    """
    member_ids = list(member_cluster_map)
    matched: dict[str, str] = {}

    for feature_id in feature_ids.astype(str).unique():
        candidates = [
            member_id
            for member_id in member_ids
            if member_id == feature_id
            or member_id.startswith(feature_id + "_")
        ]

        if candidates:
            # Choose the most specific match in the unlikely case of ambiguity.
            matched[feature_id] = max(candidates, key=len)

    return matched


def _safe_int(value: object, label: str) -> int:
    """Convert a value to int, raising a descriptive error if invalid."""
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid integer in {label}: {value!r}"
        ) from exc


def _contig_from_feature_id(feature_id: object) -> str:
    return str(feature_id).split(":", 1)[0]


def _make_neighbourhood_label(row: pd.Series) -> str:
    return f"{row['genome']} | {row['neighbourhood_id']}"


def _get_distances_around_goi(
    df: pd.DataFrame,
    goi_type: str,
) -> list[int]:
    """Get adjacent intergenic distances within each genome/neighbourhood/contig."""
    distances: list[int] = []
    group_columns = ["genome", "neighbourhood_id", "contig"]

    for _, group in df.groupby(
        group_columns,
        sort=False,
        dropna=False,
    ):
        group = group.sort_values(
            ["start", "end"],
            kind="stable",
        ).reset_index(drop=True)

        goi_indices = group.index[
            group["type"].astype(str).eq(goi_type)
        ].tolist()

        for index in goi_indices:
            current = group.iloc[index]

            if index > 0:
                previous = group.iloc[index - 1]
                distance = int(current["start"]) - int(previous["end"])

                if distance >= 0:
                    distances.append(distance)

            if index + 1 < len(group):
                following = group.iloc[index + 1]
                distance = int(following["start"]) - int(current["end"])

                if distance >= 0:
                    distances.append(distance)

    return distances


def _feature_type_plot(
    ax,
    type_labels: list[str],
    type_values: list[int],
    goi_label: str,
    accent_color: str,
    normal_color: str,
) -> None:
    """Highlight only the GOI feature-type bar; use the normal color for others."""
    x = np.arange(len(type_labels))

    # Highlight only the GOI bar.
    bar_colors = [
        accent_color if label == goi_label else normal_color
        for label in type_labels
    ]

    bars = ax.bar(x, type_values, color=bar_colors)

    ax.set_xticks(x)
    ax.set_xticklabels(type_labels, rotation=45, ha="right")
    ax.set_title(
        "A) Feature Type Counts",
        fontsize=13,
        fontweight="bold",
    )
    ax.set_ylabel("Count")

    for bar, value in zip(bars, type_values):
        ax.annotate(
            str(value),
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
        )


def _pie_label_for_color(
    color: str,
    df_for_pie: pd.DataFrame,
    gene_of_interest: str,
    color_to_gene: dict[str, str],
    goi_color: str,
) -> str:
    """Return a meaningful pie-chart label for a cluster color."""
    normalized_color = color.split("_h_")[0].upper()

    if normalized_color in {WHITE, "#FFF"}:
        return "no-cluster"

    # The GOI color comes from color.json
    if normalized_color == goi_color.upper():
        return gene_of_interest

    subset = df_for_pie[
        df_for_pie["cluster_color"].str.split("_h_").str[0].str.upper().eq(normalized_color)
    ]

    if subset["gene"].astype(str).eq(gene_of_interest).any():
        return gene_of_interest

    return str(color_to_gene.get(normalized_color, "Unknown"))


def summary(
    gene_of_interest: str,
    goi_type: str,
    input_df: pd.DataFrame,
    output_file: str,
    colors=None,
) -> None:
    """Write summary plots (PNG/SVG), summary.txt, and summary.json."""
    df = input_df.copy()

    # Use the same color scheme that was loaded by main().
    colors = colors or load_colors()

    goi_color = colors["goi"].upper()
    special_colors = {WHITE, goi_color}

    output_dir = os.path.dirname(os.path.abspath(output_file))
    os.makedirs(output_dir, exist_ok=True)

    required_columns = (
        "genome_id",
        "neighbourhood_id",
        "feature_id",
        "gene",
        "start",
        "end",
        "type",
        "cluster_color",
    )

    for required in required_columns:
        if required not in df.columns:
            raise ValueError(
                f"Cannot create summary: required column {required!r} is missing."
            )

    df["cluster_color"] = (
        df["cluster_color"]
        .fillna(WHITE)
        .astype(str)
        .str.split("_h_")
        .str[0]
        .str.upper()
        .replace({"#FFF": WHITE})
    )

    df["genome"] = (
        df["genome_id"]
        .replace("", "unknown")
        .fillna("unknown")
        .astype(str)
    )

    df["contig"] = df["feature_id"].map(_contig_from_feature_id)
    df["start"] = pd.to_numeric(
        df["start"],
        errors="raise",
    ).astype(int)
    df["end"] = pd.to_numeric(
        df["end"],
        errors="raise",
    ).astype(int)

    df["location"] = df.apply(_make_neighbourhood_label, axis=1)

    if "gene_name_based_on_DB" not in df.columns:
        if "cluster_gene_name" in df.columns:
            df["gene_name_based_on_DB"] = (
                df["cluster_gene_name"]
                .replace("", np.nan)
                .fillna(df["gene"])
            )
        else:
            df["gene_name_based_on_DB"] = df["gene"]

    # Figure A: feature-type counts.
    type_counts = df["type"].astype(str).value_counts()
    type_labels: list[str] = []
    type_values: list[int] = []
    goi_plot_label = goi_type

    for feature_type, count in type_counts.items():
        label = (
            f"{feature_type} ({gene_of_interest})"
            if feature_type == goi_type
            else str(feature_type)
        )

        if feature_type == goi_type:
            goi_plot_label = label

        type_labels.append(label)
        type_values.append(int(count))

    # Exclude the configured GOI color and white/no-cluster color.
    cluster_assigned = ~df["cluster_color"].isin(special_colors)

    location_order = df["location"].drop_duplicates().tolist()

    colors_per_location = (
        df.loc[cluster_assigned]
        .groupby("location")["cluster_color"]
        .nunique()
        .reindex(location_order, fill_value=0)
        .sort_values(ascending=False, kind="stable")
    )

    mean_colors_per_location = (
        float(colors_per_location.mean())
        if not colors_per_location.empty
        else 0.0
    )

    genome_order = df["genome"].drop_duplicates().tolist()

    colors_per_genome = (
        df.loc[cluster_assigned]
        .groupby("genome")["cluster_color"]
        .nunique()
        .reindex(genome_order, fill_value=0)
        .sort_values(ascending=False, kind="stable")
    )

    mean_colors_per_genome = (
        float(colors_per_genome.mean())
        if not colors_per_genome.empty
        else 0.0
    )

    distances = _get_distances_around_goi(df, goi_type)
    accent_color = goi_color

    num_locations = len(colors_per_location)
    if num_locations > 50:
        fig_width = 14.0
    else:
        fig_width = max(
            14.0,
            min(0.42 * max(num_locations, 1), 100.0),
        )

    fig, axs = plt.subplots(
        2,
        2,
        figsize=(fig_width, 10),
    )
    axs = axs.flatten()

    _feature_type_plot(
        axs[0],
        type_labels,
        type_values,
        goi_plot_label,
        accent_color=goi_color,
        normal_color=colors["bars"],
    )

    # Figure B: distribution of colors in the output.
    df_for_pie = df.copy()
    color_counts = df_for_pie["cluster_color"].value_counts()

    color_to_gene = {
        str(color).upper(): str(gene)
        for color, gene in (
            df_for_pie.drop_duplicates("cluster_color")
            .set_index("cluster_color")["gene_name_based_on_DB"]
            .to_dict()
            .items()
        )
    }

    pie_labels = [
        _pie_label_for_color(
            str(color),
            df_for_pie,
            gene_of_interest,
            color_to_gene,
            goi_color,
        )
        for color in color_counts.index
    ]

    pie_colors = [str(color) for color in color_counts.index]

    if len(color_counts):
        wedges, texts, autotexts = axs[1].pie(
            color_counts.values,
            labels=pie_labels,
            autopct="%1.1f%%",
            startangle=90,
            colors=pie_colors,
            textprops={"fontsize": 8, "color": "black"},
            wedgeprops={"edgecolor": "black", "linewidth": 0.5},
        )

        for label, pct_label in zip(texts, autotexts):
            pct_text = pct_label.get_text().rstrip("%")

            if pct_text:
                try:
                    if float(pct_text) < 3.0:
                        label.set_text("")
                        pct_label.set_text("")
                except ValueError:
                    pass
    else:
        axs[1].text(
            0.5,
            0.5,
            "No data",
            ha="center",
            va="center",
        )

    axs[1].set_title(
        f"B) Cluster Color Distribution "
        f"(unique colors: {df['cluster_color'].nunique()})",
        fontsize=13,
        fontweight="bold",
    )
    axs[1].axis("equal")

    # --- Figure C: Ein Bar pro Genom + Neighbourhood mit Logik für > 50 Locations ---
    location_labels = [
        str(location) for location in colors_per_location.index
    ]
    location_values = colors_per_location.to_numpy(dtype=int)
    positions = np.arange(len(location_values))

    bars = axs[2].bar(
        positions,
        location_values,
        color=colors["bars"],
    )

    axs[2].set_title(
        f"C) Distinct Cluster Colors per Neighbourhood (n={num_locations})",
        fontsize=13,
        fontweight="bold",
    )
    axs[2].set_ylabel("Distinct cluster colors (excluding white/GOI color)")

    if num_locations > 50:
        axs[2].set_xticks([])
        axs[2].set_xlabel("Genomic Neighbourhoods (>50 locations; labels hidden)")
    else:
        axs[2].set_xlabel("Genome accession | Neighbourhood ID")
        axs[2].set_xticks(positions)
        axs[2].set_xticklabels(
            location_labels,
            rotation=90,
            ha="center",
            fontsize=7,
        )
        for bar, count in zip(bars, location_values):
            axs[2].annotate(
                str(int(count)),
                xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
            )

    if num_locations:
        axs[2].axhline(
            mean_colors_per_location,
            color=accent_color,
            linestyle="--",
            linewidth=1.8,
        )

        axs[2].text(
            max(num_locations - 1, 0),
            mean_colors_per_location + 0.1,
            f"Mean: {mean_colors_per_location:.2f}",
            fontsize=8,
            fontweight="bold",
            ha="right",
            va="bottom",
        )

    axs[2].margins(x=0.01, y=0.15)

    # --- Figure D: Intergenische Distanzen mit größeren Bins ---
    if distances:
        max_distance = max(distances)
        bin_step = 50  # Vergrößerte Bin-Größe (z.B. 50 nt statt 10 nt)
        bin_edges = np.arange(0, max_distance + bin_step, bin_step)

        if len(bin_edges) < 2:
            bin_edges = np.array([0, bin_step])

        axs[3].hist(
            distances,
            bins=bin_edges,
            color=colors["histogram"],
            alpha=0.9,
            edgecolor=colors["text"],
        )

        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        tick_step = max(1, len(bin_centers) // 10)
        tick_indices = np.arange(0, len(bin_centers), tick_step)

        axs[3].set_xticks(bin_centers[tick_indices])
        axs[3].set_xticklabels(
            [
                str(int(value))
                for value in bin_centers[tick_indices]
            ],
            rotation=45,
            ha="right",
            fontsize=8,
        )
        axs[3].set_title(
            "D) Distance Distribution around Gene of Interest",
            fontsize=13,
            fontweight="bold",
        )
        axs[3].set_xlabel("Intergenic distance (nt)")
        axs[3].set_ylabel("Frequency")
    else:
        axs[3].text(
            0.5,
            0.5,
            "No adjacent gene distances found",
            ha="center",
            va="center",
            color="gray",
        )
        axs[3].set_title(
            "D) Distance Distribution",
            fontsize=13,
            fontweight="bold",
        )
        axs[3].set_axis_off()

    fig.tight_layout(pad=3.0)
    fig.savefig(
        os.path.join(output_dir, "summary.png"),
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        os.path.join(output_dir, "summary.svg"),
        bbox_inches="tight",
    )
    plt.close(fig)

    # Summary statistics.
    white_mask = df["cluster_color"].eq(WHITE)
    goi_color_mask = df["cluster_color"].eq(goi_color)

    goi_mask = (
        df["gene"].astype(str).eq(gene_of_interest)
        | df["feature_id"].astype(str).str.contains(
            gene_of_interest,
            case=False,
            regex=False,
            na=False,
        )
    )

    feature_type_counts = {
        str(key): int(value)
        for key, value in df["type"].value_counts(dropna=False).items()
    }

    unique_features = int(
        df.drop_duplicates(
            subset=["genome", "feature_id"]
        ).shape[0]
    )
    unique_genomes = int(df["genome"].nunique())
    unique_contig_ids = int(df["contig"].nunique())
    genome_contig_combinations = int(
        df[["genome", "contig"]].drop_duplicates().shape[0]
    )
    unique_neighbourhoods = int(
        df[["genome", "neighbourhood_id"]].drop_duplicates().shape[0]
    )
    unique_cluster_colors = int(
        df.loc[cluster_assigned, "cluster_color"].nunique()
    )
    unique_colors_in_output = int(df["cluster_color"].nunique())
    neighbourhoods_with_colors = int((colors_per_location > 0).sum())
    max_colors_at_location = (
        int(colors_per_location.max())
        if not colors_per_location.empty
        else 0
    )
    location_with_max_colors = (
        str(colors_per_location.idxmax())
        if not colors_per_location.empty
        else "None"
    )

    if "cluster_gene_name" in df.columns:
        unique_database_gene_labels = int(
            df["cluster_gene_name"]
            .replace("", np.nan)
            .dropna()
            .nunique()
        )
    else:
        unique_database_gene_labels = int(
            df["gene_name_based_on_DB"]
            .replace("", np.nan)
            .dropna()
            .nunique()
        )

    summary_stats = {
        "total_feature_rows": int(len(df)),
        "unique_features": unique_features,
        "unique_genomes": unique_genomes,
        "unique_contig_ids": unique_contig_ids,
        "genome_contig_combinations": genome_contig_combinations,
        "unique_neighbourhoods": unique_neighbourhoods,
        "unique_cluster_colors_excluding_white_and_goi": unique_cluster_colors,
        "unique_colors_including_white_and_goi": unique_colors_in_output,
        "feature_rows_assigned_a_cluster_color": int(cluster_assigned.sum()),
        "feature_rows_with_white_no_cluster_color": int(white_mask.sum()),
        "feature_rows_highlighted_with_goi_color": int(goi_color_mask.sum()),
        "gene_of_interest_matching_rows": int(goi_mask.sum()),
        "neighbourhoods_with_cluster_colors": neighbourhoods_with_colors,
        "mean_cluster_colors_per_neighbourhood": round(
            mean_colors_per_location,
            4,
        ),
        "max_cluster_colors_in_neighbourhood": max_colors_at_location,
        "neighbourhood_with_max_cluster_colors": location_with_max_colors,
        "unique_database_gene_labels": unique_database_gene_labels,
        "feature_type_counts": feature_type_counts,
    }

    text_labels = {
        "total_feature_rows": "Total feature rows",
        "unique_features": "Unique features (genome + feature ID)",
        "unique_genomes": "Unique genomes",
        "unique_contig_ids": "Unique contig IDs (not genome-qualified)",
        "genome_contig_combinations": "Genome-contig combinations",
        "unique_neighbourhoods": "Unique genome-neighbourhood combinations",
        "unique_cluster_colors_excluding_white_and_goi":
            "Distinct cluster colors (excluding white and GOI color)",
        "unique_colors_including_white_and_goi":
            "Distinct colors present in the output",
        "feature_rows_assigned_a_cluster_color":
            "Feature rows assigned a cluster color",
        "feature_rows_with_white_no_cluster_color":
            "Feature rows with white/no-cluster color",
        "feature_rows_highlighted_with_goi_color":
            "Feature rows highlighted with the configured GOI color",
        "gene_of_interest_matching_rows":
            "Rows matching the gene of interest",
        "neighbourhoods_with_cluster_colors":
            "Neighbourhoods with at least one cluster color",
        "mean_cluster_colors_per_neighbourhood":
            "Mean cluster colors per neighbourhood",
        "max_cluster_colors_in_neighbourhood":
            "Maximum cluster colors in one neighbourhood",
        "neighbourhood_with_max_cluster_colors":
            "Neighbourhood with the maximum number of cluster colors",
        "unique_database_gene_labels":
            "Unique database gene/family labels",
    }

    summary_txt_path = os.path.join(output_dir, "summary.txt")

    with open(summary_txt_path, "w", encoding="utf-8") as handle:
        handle.write("GENOMIC NEIGHBOURHOOD CLUSTERING SUMMARY\n")
        handle.write("=" * 44 + "\n\nGENERAL STATISTICS\n")
        handle.write("-" * 44 + "\n")

        for key, label in text_labels.items():
            handle.write(f"{label}: {summary_stats[key]}\n")

        handle.write("\nFEATURE TYPE COUNTS\n" + "-" * 44 + "\n")

        for feature_type, count in sorted(feature_type_counts.items()):
            handle.write(f"{feature_type}: {count}\n")

        handle.write(
            "\nCLUSTER COLORS BY GENOME + NEIGHBOURHOOD\n"
            + "-" * 44
            + "\n"
        )
        handle.write("Genome | Neighbourhood ID\tDistinct cluster colors\n")

        for location, count in colors_per_location.items():
            handle.write(f"{location}\t{int(count)}\n")

    stats_data = {
        "summary": summary_stats,
        "plot1_types": {
            "labels": type_labels,
            "values": type_values,
        },
        "plot2_pie": {
            "labels": pie_labels,
            "values": [
                int(value) for value in color_counts.values
            ],
            "colors": pie_colors,
        },
        "plot3_genomes": {
            "labels": [
                str(label) for label in colors_per_genome.index
            ],
            "values": [
                int(value) for value in colors_per_genome.values
            ],
            "mean": mean_colors_per_genome,
        },
        "plot3_locations": {
            "labels": [
                str(label) for label in colors_per_location.index
            ],
            "values": [
                int(value) for value in colors_per_location.values
            ],
            "mean": mean_colors_per_genome,
        },
        "plot4_distances": {
            "values": [float(value) for value in distances],
        },
    }

    with open(
        os.path.join(output_dir, "summary.json"),
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(stats_data, handle, indent=4)


def main() -> None:
    parser = setup_parser()
    args = parser.parse_args()

    validate_files(args.cluster_file, args.tsv_file)

    colors = load_colors(args.color_scheme_file)
    goi_color = colors["goi"].upper()

    special_colors = {
        WHITE,
        goi_color,
    }

    color_map, gene_name_map, member_cluster_map = process_clusters(
        args.cluster_file,
        args.gene_of_interest,
        args.size_for_cluster,
        goi_color=goi_color,
        special_colors=special_colors,
    )

    df = read_neighbourhood_tsv(args.tsv_file)
    df.reset_index(drop=True, inplace=True)

    df["genome_id"] = (
        df["genome_id"]
        .replace("", "unknown")
        .fillna("unknown")
    )
    df["genome"] = df["genome_id"].astype(str)
    df["contig"] = df["feature_id"].map(_contig_from_feature_id)
    df["start"] = df["start"].map(
        lambda value: _safe_int(value, "start")
    )
    df["end"] = df["end"].map(
        lambda value: _safe_int(value, "end")
    )

    matched_member_by_feature = match_cluster_members(
        df["feature_id"],
        member_cluster_map,
    )

    df["cluster_member_id"] = (
        df["feature_id"]
        .map(matched_member_by_feature)
        .fillna("")
    )
    df["cluster_id"] = (
        df["cluster_member_id"]
        .map(member_cluster_map)
        .fillna("")
    )
    df["cluster_gene_name"] = (
        df["cluster_member_id"]
        .map(gene_name_map)
        .fillna("")
    )
    df["cluster_color"] = (
        df["cluster_member_id"]
        .map(color_map)
        .fillna(WHITE)
    )

    goi_mask = (
        df["gene"].astype(str).str.casefold()
        == args.gene_of_interest.casefold()
    ) | df["feature_id"].astype(str).str.contains(
        args.gene_of_interest,
        case=False,
        regex=False,
        na=False,
    )

    df.loc[goi_mask, "cluster_color"] = goi_color

    hatches = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']
    cluster_hatch = pd.Series("", index=df.index)

    colored_mask = (
        (df["cluster_color"].str.upper() != WHITE) &
        (df["cluster_color"].str.upper() != goi_color) &
        (~goi_mask)
    )
    eligible_df = df[colored_mask]

    if not eligible_df.empty:
        group_col = "cluster_id" if "cluster_id" in eligible_df.columns and eligible_df["cluster_id"].ne("").any() else "cluster_gene_name"
        if eligible_df[group_col].eq("").all():
            group_col = "gene"

        cluster_counts = eligible_df[group_col].value_counts()
        num_top_clusters = max(1, int(np.ceil(0.25 * len(cluster_counts))))
        top_clusters = cluster_counts.head(num_top_clusters).index.tolist()

        rng = np.random.default_rng(42)
        cluster_to_hatch = {
            cid: rng.choice(hatches)
            for cid in top_clusters
        }

        for idx, row in eligible_df.iterrows():
            cid = row[group_col]
            if cid in cluster_to_hatch:
                cluster_hatch.at[idx] = cluster_to_hatch[cid]

    df["cluster_color"] = np.where(
        cluster_hatch.ne(""),
        df["cluster_color"] + "_h_" + cluster_hatch,
        df["cluster_color"]
    )

    df["gene_name_based_on_DB"] = (
        df["cluster_gene_name"]
        .replace("", np.nan)
        .fillna(df["gene"])
    )

    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df[OUTPUT_COLUMNS].to_csv(
        output_path,
        sep="\t",
        index=False,
        header=False,
    )
    summary(
        args.gene_of_interest,
        args.goi_type,
        df,
        str(output_path),
        colors=colors,
    )

    print(f"Wrote colored neighbourhood TSV: {output_path}")
    print(f"Wrote summary files to: {output_path.parent.resolve()}")


if __name__ == "__main__":
    main()
