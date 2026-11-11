# (C) 2024, Maria Schreiber, MIT license
"""
microsynteny_stats.py  (v3)

Conservation statistics and a summary figure (conservation map + neighbour-similarity
histogram) for one synteny cluster. Uses exactly the same alignment settings as the
logo, so numbers in the two figures agree. Works for any clustering method
(hierarchical or HDBSCAN/"density"); noise labels (-1, 9999) are titled accordingly.

Outputs (file-name suffix letter configurable via `suffix`):
    microsynteny_stats_<algo>_cluster_<id>_C.png/.svg        2-panel figure (map | histogram)
    microsynteny_stats_<algo>_cluster_<id>_summary.csv
    microsynteny_stats_<algo>_cluster_<id>_families.csv
    microsynteny_stats_<algo>_cluster_<id>_similarity_cosine.csv / _jaccard.csv
"""
import os
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from plot_microsynteny_logo import (
    WHITE, parse_color_and_hatch, align_microsynteny_neighborhoods, column_information,
)

logger = logging.getLogger(__name__)

NOISE_LABELS = (-1, 9999)     # -1 = HDBSCAN / hierarchical noise, 9999 = paralog overflow


def _cluster_title(cluster_label):
    try:
        if int(cluster_label) == -1:
            return "unclustered (noise)"
        if int(cluster_label) == 9999:
            return "paralog overflow"
    except (TypeError, ValueError):
        pass
    return f"cluster {cluster_label}"

def _lcs(a, b):
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    for x in a:
        cur = [0]
        for j, y in enumerate(b):
            cur.append(prev[j] + 1 if x == y else max(prev[j + 1], cur[j]))
        prev = cur
    return prev[-1]


def _lcs_score(a, b):
    """a, b = (left_tokens_outward, right_tokens_outward); order-aware similarity in [0, 1]."""
    den = max(len(a[0]), len(b[0])) + max(len(a[1]), len(b[1]))
    return (_lcs(a[0], b[0]) + _lcs(a[1], b[1])) / den if den else np.nan

def _neighbourhoods(aligned, genome_col, use_strand):
    """{genome: (left_outward, right_outward)} built from the aligned table (GOI excluded).

    Gap columns (token NaN) are dropped, so the order of genes on each side is exactly
    what the logo shows.
    """
    body = aligned[(aligned["aligned_pos"] != 0) & aligned["token"].notna()].copy()
    body["_t"] = body["token"] + body["strand"].fillna("+") if use_strand else body["token"]
    hoods = {g: ([], []) for g in aligned[genome_col].unique()}
    left = body[body["aligned_pos"] < 0].sort_values("aligned_pos", ascending=False)
    right = body[body["aligned_pos"] > 0].sort_values("aligned_pos")
    for g, t in zip(left[genome_col], left["_t"]):
        hoods[g][0].append(t)
    for g, t in zip(right[genome_col], right["_t"]):
        hoods[g][1].append(t)
    return hoods

def _presence_matrix(hoods):
    keys = list(hoods)
    toks = sorted({t for l, r in hoods.values() for t in (*l, *r)})
    idx = {t: i for i, t in enumerate(toks)}
    M = np.zeros((len(keys), len(toks)))
    for i, k in enumerate(keys):
        for t in (*hoods[k][0], *hoods[k][1]):
            M[i, idx[t]] = 1.0
    return keys, M

def _set_similarities(M):
    """Vectorised Jaccard and cosine on presence/absence vectors (GOI excluded).

    Neighbourhoods without any flanking gene have no defined similarity -> NaN row/column
    (including the diagonal), so they never inflate or deflate the means.
    """
    inter = M @ M.T
    size = M.sum(axis=1)
    union = size[:, None] + size[None, :] - inter
    norm = np.sqrt(np.outer(size, size))
    with np.errstate(divide="ignore", invalid="ignore"):
        J = np.where(union > 0, inter / union, np.nan)
        C = np.where(norm > 0, inter / norm, np.nan)
    ok = size > 0
    np.fill_diagonal(J, np.where(ok, 1.0, np.nan))
    np.fill_diagonal(C, np.where(ok, 1.0, np.nan))
    return J, C

def _nanmean(x, default=np.nan):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(x.mean()) if x.size else default

def _short_label(df, genome_col):
    if "organism_name" in df.columns:
        lab = df.groupby(genome_col)["organism_name"].first().astype(str).str.split(r"\s*\|\s*").str[0]
        lab = lab.str.replace(r"^([A-Z])[a-z]+\s+", r"\1. ", regex=True)
    else:
        lab = pd.Series({g: str(g) for g in df[genome_col].unique()})
    if lab.duplicated().any():
        lab = lab + " (" + lab.index.astype(str) + ")"
    return lab.to_dict()

def _family_table(aligned, genome_col):
    """Per-family conservation, using the SAME aligned columns as the logo."""
    n = aligned[genome_col].nunique()
    body = aligned[(aligned["aligned_pos"] != 0) & aligned["token"].notna()]
    t = body.groupby("token").agg(
        n=(genome_col, "nunique"),
        median_pos=("aligned_pos", "median"),
        sd_pos=("aligned_pos", "std"),
        same_strand=("strand", lambda s: (s == "+").mean()),
    )
    t["prevalence"] = t["n"] / n
    raw = body.dropna(subset=["colour"]).groupby("token")["colour"].first()
    parsed = raw.map(parse_color_and_hatch)
    t["colour"] = parsed.map(lambda x: x[0]).reindex(t.index).fillna("#808080")
    t["hatch"] = parsed.map(lambda x: x[1] or "").reindex(t.index).fillna("")
    t["unannotated"] = t["colour"].str.lower().isin(WHITE)
    return t.sort_values(["prevalence", "median_pos"], ascending=[False, True])

def _save(fig, base, dpi, svg):
    fig.savefig(base + ".png", dpi=dpi, bbox_inches="tight")
    if svg:
        fig.savefig(base + ".svg", bbox_inches="tight")
    plt.close(fig)


def _letter(fig, letter):
    if letter:
        fig.text(0.005, 0.995, letter, fontsize=18, fontweight="bold", va="top", ha="left")


def _sample_pairs(n, max_pairs, seed):
    """Row/column indices of (at most max_pairs) distinct pairs i < j, without building a list."""
    iu = np.triu_indices(n, 1)
    total = len(iu[0])
    if total <= max_pairs:
        return iu, total
    pick = np.random.default_rng(seed).choice(total, max_pairs, replace=False)
    return (iu[0][pick], iu[1][pick]), total

def _draw_conservation_map(ax, fam, title, label_min_prevalence=None):
    """Fraction of neighbourhoods containing a family vs its median aligned position.

    Bubbles use the logo's colours AND hatches. Unannotated genes (white, mostly one-off
    accessions) are not drawn. label_min_prevalence=None -> no names next to the bubbles.
    """
    shown = fam[~fam["unannotated"]] if "unannotated" in fam else fam
    n_hidden = len(fam) - len(shown)
    if len(shown):
        with plt.rc_context({"hatch.linewidth": 1.0}):
            # one call per hatch style (hatch is a collection-level property)
            for hatch, grp in shown.groupby("hatch"):
                size = 50 + 250 * grp["prevalence"].to_numpy()
                xy = np.column_stack([grp["median_pos"].to_numpy(), grp["prevalence"].to_numpy()])
                ax.scatter(xy[:, 0], xy[:, 1], s=size, c=list(grp["colour"]),
                           edgecolor="white", linewidth=0, hatch=(hatch * 2 if hatch else None),
                           alpha=0.95, zorder=3)
                ax.scatter(xy[:, 0], xy[:, 1], s=size, facecolor="none",
                           edgecolor="black", linewidth=0.6, zorder=4)
        if label_min_prevalence is not None:
            offsets = [(6, 6), (6, -13), (-6, 6), (-6, -13)]
            for i, (tok, r) in enumerate(shown[shown.prevalence >= label_min_prevalence].iterrows()):
                dx, dy = offsets[i % len(offsets)]
                ax.annotate(tok, (r.median_pos, r.prevalence), xytext=(dx, dy), textcoords="offset points",
                            ha="left" if dx > 0 else "right", fontsize=8, zorder=5,
                            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.75))
        lo, hi = shown.median_pos.min(), shown.median_pos.max()
        ax.set_xlim(min(lo, 0) - 1, max(hi, 0) + 1)
    else:
        ax.set_xlim(-1, 1)
    ax.scatter([0], [1.0], marker="*", s=260, c="#d9a441", edgecolor="black", lw=0.7, zorder=5)
    ax.annotate("GOI", (0, 1.0), xytext=(8, -3), textcoords="offset points", fontsize=8, fontweight="bold")
    ax.axvline(0, color="gray", ls="--", lw=1.0, alpha=0.7)
    ax.set_ylim(0, 1.08)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Median aligned position relative to GOI", fontsize=10)
    ax.set_ylabel("Fraction of neighbourhoods", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.3, linestyle=":", linewidth=0.5)
    if n_hidden:
        ax.text(0.99, 0.015, f"{n_hidden} unannotated genes not shown", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7, color="0.4")


def _draw_similarity_hist(ax, cos_vals, title):
    cos_vals = np.asarray(cos_vals, dtype=float)
    cos_vals = cos_vals[np.isfinite(cos_vals)]
    bins = np.linspace(0, 1, 11)
    if cos_vals.size:
        ax.hist(cos_vals, bins=bins, color="#4c72b0", edgecolor="black", lw=0.6, alpha=0.85)
        ax.axvline(cos_vals.mean(), color="red", ls="--", lw=1.5, alpha=0.8,
                   label=f"Mean = {cos_vals.mean():.2f}")
        ax.legend(framealpha=0.9, fontsize=9)
    else:
        ax.text(0.5, 0.5, "no pairs to compare", transform=ax.transAxes, ha="center", va="center",
                color="0.4")
    ax.set_xlabel("Pairwise cosine similarity of gene content", fontsize=10)
    ax.set_ylabel("Number of neighbourhood pairs", fontsize=10)
    ax.set_xlim(0, 1)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y", linestyle=":", linewidth=0.5)

def plot_cluster_statistics(
    df, output_dir, cluster_label, clustering_algo,
    genome_col="genome", rel_pos_col="rel_pos", family_col="gene_with_orientation",
    color_col="cluster_color", use_strand=False, max_pairs=3000,
    dpi=300, svg=True, seed=1,
    # keep these identical to the plot_microsynteny_logo call:
    drop_unannotated=False, min_family_fraction=0.0, small_sample_correction=False,
    match_score=3, mismatch_score=-0.9, gap_score=-1,
    label_min_prevalence=None,      # None = no names on the bubbles (use the logo legend)
    suffix="C", panel_letter="C",
):
    os.makedirs(output_dir, exist_ok=True)

    try:
        aligned = align_microsynteny_neighborhoods(
            df, genome_col, rel_pos_col, family_col, color_col=color_col,
            drop_unannotated=drop_unannotated, min_family_fraction=min_family_fraction,
            match_score=match_score, mismatch_score=mismatch_score, gap_score=gap_score)
    except ValueError as exc:       # aligner raises when no neighbourhood contains the GOI
        logger.warning("Cluster %s: %s; statistics skipped.", cluster_label, exc)
        return None
    if aligned.empty:
        logger.warning("Cluster %s has 0 usable neighbourhoods; statistics skipped.", cluster_label)
        return None
    genomes = list(aligned[genome_col].unique())
    n = len(genomes)

    labels = _short_label(df, genome_col)
    hoods = _neighbourhoods(aligned, genome_col, use_strand)
    keys, M = _presence_matrix(hoods)

    # --- pairwise similarity: set-based ones for ALL pairs (vectorised), LCS on a sample if huge
    J, C = _set_similarities(M)
    iu = np.triu_indices(n, 1)
    cos_vals = C[iu]
    jac = _nanmean(J[iu], default=1.0 if n < 2 else np.nan)
    cos_mean = _nanmean(cos_vals, default=1.0 if n < 2 else np.nan)

    (pi, pj), n_pairs_total = _sample_pairs(n, max_pairs, seed)
    lcs_vals = [_lcs_score(hoods[keys[i]], hoods[keys[j]]) for i, j in zip(pi, pj)]
    lcs_mean = _nanmean(lcs_vals, default=1.0 if n < 2 else np.nan)
    sim = float(np.nanmean([jac, lcs_mean])) if np.isfinite([jac, lcs_mean]).any() else np.nan

    fam = _family_table(aligned, genome_col)

    # --- information content (same alignment / settings as the logo)
    info, S = column_information(aligned, genome_col, small_sample_correction=small_sample_correction)
    col_bits = info.groupby("aligned_pos")["height"].sum() if len(info) else pd.Series(dtype=float)

    summary = dict(
        cluster=cluster_label,
        n_neighbourhoods=n,
        n_families=len(fam),
        total_bits=round(float(col_bits.sum()), 3),
        bits_per_column=round(float(col_bits.mean()), 3) if len(col_bits) else 0.0,
        max_possible_bits_per_column=round(float(np.log2(S)), 3) if S > 1 else 0.0,
        core_families_ge80pct=int((fam.prevalence >= 0.8).sum()),
        mean_pairwise_jaccard=round(jac, 4),
        mean_pairwise_lcs=round(lcs_mean, 4),
        mean_pairwise_cosine=round(cos_mean, 4),
        synteny_similarity=round(sim, 4),
        n_pairs_lcs=len(pi),
        n_pairs_total=n_pairs_total,
    )

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    _draw_conservation_map(axs[0], fam, "Conservation map", label_min_prevalence)
    _draw_similarity_hist(axs[1], cos_vals if n > 1 else np.array([1.0]), "Neighbour similarity")
    fig.suptitle(
        f"Microsynteny statistics - {clustering_algo} - {_cluster_title(cluster_label)} (n = {n})\n"
        f"mean Jaccard = {jac:.2f}, mean LCS = {lcs_mean:.2f}, mean cosine = {cos_mean:.2f}; "
        f"{summary['total_bits']:.1f} bits total", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    _letter(fig, panel_letter)

    base = os.path.join(output_dir, f"microsynteny_stats_{clustering_algo}_cluster_{cluster_label}")
    _save(fig, base + (f"_{suffix}" if suffix else ""), dpi, svg)

    pd.DataFrame([summary]).to_csv(base + "_summary.csv", index=False)
    fam.drop(columns=["colour"]).rename_axis("family").to_csv(base + "_families.csv")   # keep family names!
    names = [labels.get(k, str(k)) for k in keys]
    pd.DataFrame(C, index=names, columns=names).to_csv(base + "_similarity_cosine.csv")

    return summary
