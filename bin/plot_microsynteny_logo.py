# (C) 2024, Maria Schreiber, MIT license
"""
Reference-guided alignment and information-content logo of microsyntenic
gene neighbourhoods around a gene of interest (GOI).
"""

import os
import re
import logging
import random
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

logger = logging.getLogger(__name__)

GAP = None
WHITE = {"#ffffff", "#fff", "white"}
UNANNOTATED_GREY = "#d0d0d0"
_ORI_SUFFIX = re.compile(r"\s*\(([+-])\)\s*$")


# ============================================================
# Helpers
# ============================================================

def parse_color_and_hatch(color_val):
    """Split colour string with optional _h_ hatch into base colour and hatch."""
    if not isinstance(color_val, str):
        return "#808080", None
    parts = color_val.split("_h_")
    base_color = parts[0].strip()
    if not base_color or base_color.lower() == "nan":
        base_color = "#808080"
    hatch = parts[1].strip() if len(parts) > 1 else None
    if not hatch or hatch.lower() == "nan":
        hatch = None
    return base_color, hatch


def _is_white(color_val):
    if not isinstance(color_val, str):
        return False
    return parse_color_and_hatch(color_val)[0].lower() in WHITE


def _text_color(hex_color):
    try:
        r, g, b = mcolors.to_rgb(hex_color)
    except ValueError:
        return "black"
    return "black" if (0.299 * r + 0.587 * g + 0.114 * b) > 0.55 else "white"


def _prepare(df, genome_col, rel_pos_col, family_col, color_col,
             strip_orientation, drop_unannotated):
    """Return a clean frame with columns: genome_col, rel_pos_col, token, strand, color."""
    data = df.dropna(subset=[genome_col, rel_pos_col, family_col]).copy()
    data[rel_pos_col] = pd.to_numeric(data[rel_pos_col], errors="coerce")
    data = data.dropna(subset=[rel_pos_col])
    fam = data[family_col].astype(str)

    if strip_orientation:
        parsed = fam.str.extract(_ORI_SUFFIX)[0]
        data["token"] = fam.str.replace(_ORI_SUFFIX, "", regex=True).str.strip()
    else:
        parsed = pd.Series(np.nan, index=data.index)
        data["token"] = fam
    if "orientation" in data.columns:
        strand = data["orientation"].map({"+": "+", "-": "-", "sense": "+", "antisense": "-"})
        data["strand"] = strand.fillna(parsed).fillna("+")
    else:
        data["strand"] = parsed.fillna("+")

    if color_col in data.columns:
        data["_color"] = data[color_col]
    else:
        data["_color"] = np.nan

    if drop_unannotated:
        bad = data["token"].str.lower().isin(["", "nan", "none", "no_hit", "unknown"])
        bad |= data["_color"].map(_is_white)
        n_drop = int(bad.sum())
        # never drop the GOI row
        bad &= data[rel_pos_col] != 0
        if n_drop:
            logger.info("Removed %d unannotated (white) genes before alignment.", int(bad.sum()))
        data = data[~bad].copy()
    return data


# ============================================================
# 1. Global alignment of gene-family tokens
# ============================================================

def global_token_alignment(reference, query, match_score=3, mismatch_score=-0.9, gap_score=-1,
                           free_end_gaps=True):
    """Needleman-Wunsch on lists of family tokens. Returns two lists with None for gaps.

    free_end_gaps: trailing gaps (far from the GOI) cost nothing, so genes next to the
    GOI are kept together instead of being shifted by a gap at the GOI side.
    """
    n, m = len(reference), len(query)
    S = np.zeros((n + 1, m + 1))
    T = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up (gap in query), 2 left (gap in ref)
    S[1:, 0] = np.arange(1, n + 1) * gap_score
    T[1:, 0] = 1
    S[0, 1:] = np.arange(1, m + 1) * gap_score
    T[0, 1:] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = S[i - 1, j - 1] + (match_score if reference[i - 1] == query[j - 1] else mismatch_score)
            u = S[i - 1, j] + (0 if (free_end_gaps and j == m) else gap_score)
            l = S[i, j - 1] + (0 if (free_end_gaps and i == n) else gap_score)
            b = max(d, u, l)
            S[i, j] = b
            T[i, j] = 0 if d == b else (1 if u == b else 2)
    ar, aq = [], []
    i, j = n, m
    while i > 0 or j > 0:
        t = T[i, j]
        if i > 0 and j > 0 and t == 0:
            ar.append(reference[i - 1]); aq.append(query[j - 1]); i -= 1; j -= 1
        elif i > 0 and (j == 0 or t == 1):
            ar.append(reference[i - 1]); aq.append(None); i -= 1
        else:
            ar.append(None); aq.append(query[j - 1]); j -= 1
    return ar[::-1], aq[::-1]


def _pair_score(a, b, match_score=3, mismatch_score=-0.9, gap_score=-1):
    ar, aq = global_token_alignment(a, b, match_score, mismatch_score, gap_score)
    return sum(match_score if (x is not None and x == y) else (gap_score if (x is None or y is None) else mismatch_score)
               for x, y in zip(ar, aq))


def _choose_center(sequences, match_score=3, mismatch_score=-0.9, gap_score=-1,
                   n_candidates=25, n_targets=60, seed=1):
    """Most central sequence = highest summed alignment score to a sample of the others."""
    idx = list(range(len(sequences)))
    rng = random.Random(seed)
    cands = rng.sample(idx, min(n_candidates, len(idx)))
    cands.append(max(idx, key=lambda i: len(sequences[i])))
    targets = rng.sample(idx, min(n_targets, len(idx)))
    best, best_score = cands[0], -np.inf
    for c in set(cands):
        s = sum(_pair_score(sequences[c], sequences[t], match_score, mismatch_score, gap_score)
                for t in targets if t != c)
        if s > best_score:
            best, best_score = c, s
    return best


# ============================================================
# 2. Centre-star multiple alignment (with aligned insertions)
# ============================================================

def align_neighbourhood_group(sequences, match_score=3, mismatch_score=-0.9, gap_score=-1):
    """Align token lists to a central sequence; insertions are aligned within each slot."""
    if not sequences:
        return []
    centre = sequences[_choose_center(sequences, match_score, mismatch_score, gap_score)]
    if not centre:
        # centre has no genes; fall back to the longest
        centre = max(sequences, key=len)
        if not centre:
            return [[] for _ in sequences]
    n_ref = len(centre)

    maps = []
    slot_lists = [[] for _ in range(n_ref + 1)]
    for seq in sequences:
        ar, aq = global_token_alignment(centre, seq, match_score, mismatch_score, gap_score)
        ref_map = [None] * n_ref
        ins = [[] for _ in range(n_ref + 1)]
        p = 0
        for r, q in zip(ar, aq):
            if r is None:
                ins[p].append(q)
            else:
                if q is not None:
                    ref_map[p] = q
                p += 1
        maps.append((ref_map, ins))
        for s in range(n_ref + 1):
            if ins[s]:
                slot_lists[s].append(ins[s])

    # Build consensus insertion columns per slot (most common insertion lists first).
    slot_cols = []
    for lists in slot_lists:
        cols = []
        ordered = [l for l, _ in Counter(tuple(l) for l in lists).most_common()]
        for l in ordered:
            l = list(l)
            if not cols:
                cols = l[:]
                continue
            ar, aq = global_token_alignment(cols, l, match_score, mismatch_score, gap_score)
            cols = [r if r is not None else q for r, q in zip(ar, aq)]
        slot_cols.append(cols)

    out = []
    for ref_map, ins in maps:
        row = []
        for s in range(n_ref + 1):
            cols = slot_cols[s]
            filled = [None] * len(cols)
            if ins[s]:
                ar, aq = global_token_alignment(cols, ins[s], match_score, mismatch_score, gap_score)
                k = 0
                for r, q in zip(ar, aq):
                    if r is not None:
                        if q is not None:
                            filled[k] = q
                        k += 1
            row.extend(filled)
            if s < n_ref:
                row.append(ref_map[s])
        out.append(row)
    return out


# ============================================================
# 3. Align complete neighbourhoods around the GOI
# ============================================================

def align_microsynteny_neighborhoods(
    df,
    genome_col="genome",
    rel_pos_col="rel_pos",
    family_col="gene_family",
    goi_pos=0,
    max_flanking_genes=None,
    color_col="cluster_color",
    strip_orientation=True,
    drop_unannotated=True,
    min_family_fraction=0.02,
    match_score=3,
    mismatch_score=-0.9,
    gap_score=-1,
):
    """Return a long dataframe: one row per (neighbourhood, aligned column).

    Columns: genome_col, aligned_pos, token (family without orientation),
    strand, colour, and family_col (kept = token, for backwards compatibility).
    """
    for c in (genome_col, rel_pos_col, family_col):
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
    data = _prepare(df, genome_col, rel_pos_col, family_col, color_col,
                    strip_orientation, drop_unannotated)
    if max_flanking_genes is not None:
        data = data[data[rel_pos_col].between(goi_pos - max_flanking_genes,
                                              goi_pos + max_flanking_genes)].copy()
    if min_family_fraction:
        # families seen in < min_family_fraction of neighbourhoods would only create
        # near-empty columns; remove them before alignment (GOI is always kept)
        n_h = data[genome_col].nunique()
        prev = data.groupby("token")[genome_col].nunique() / n_h
        rare = set(prev[prev < min_family_fraction].index)
        data = data[~(data["token"].isin(rare) & (data[rel_pos_col] != goi_pos))].copy()

    hoods = {}
    for g, grp in data.groupby(genome_col, sort=False):
        goi = grp[grp[rel_pos_col] == goi_pos]
        if goi.empty:
            logger.warning("No GOI found for neighbourhood %s. Skipping.", g)
            continue
        left = grp[grp[rel_pos_col] < goi_pos].sort_values(rel_pos_col, ascending=False)
        right = grp[grp[rel_pos_col] > goi_pos].sort_values(rel_pos_col)
        hoods[g] = dict(goi=goi.iloc[0], left=left, right=right)
    if not hoods:
        raise ValueError("No valid neighbourhoods containing the GOI were found.")

    genomes = list(hoods)
    al_left = align_neighbourhood_group([hoods[g]["left"]["token"].tolist() for g in genomes],
                                        match_score, mismatch_score, gap_score)
    al_right = align_neighbourhood_group([hoods[g]["right"]["token"].tolist() for g in genomes],
                                         match_score, mismatch_score, gap_score)
    lw = max((len(s) for s in al_left), default=0)

    rows = []
    for gi, g in enumerate(genomes):
        h = hoods[g]
        # left alignment runs outward from the GOI; reverse into genomic order
        left_tokens = list(reversed(al_left[gi] + [None] * (lw - len(al_left[gi]))))
        # map tokens back to strand/colour (tokens are unique per side in a neighbourhood
        # except for tandem repeats -> consume in order)
        def attach(tokens, frame):
            pool = defaultdict(list)
            for _, r in frame.iterrows():
                pool[r["token"]].append((r["strand"], r["_color"]))
            res = []
            for t in tokens:
                if t is None:
                    res.append((None, None, None))
                else:
                    s, c = pool[t].pop(0)
                    res.append((t, s, c))
            return res

        # left tokens in original (outward) order for pool consumption
        outward = al_left[gi] + [None] * (lw - len(al_left[gi]))
        left_attached = list(reversed(attach(outward, h["left"])))
        right_attached = attach(al_right[gi], h["right"])
        goi_att = [(h["goi"]["token"], h["goi"]["strand"], h["goi"]["_color"])]
        for pos, (t, s, c) in enumerate(left_attached + goi_att + right_attached):
            rows.append({genome_col: g, "aligned_pos": pos - lw, "token": t,
                         "strand": s, "colour": c, family_col: t})
    return pd.DataFrame(rows)


# ============================================================
# 4. Information content
# ============================================================

def column_information(aligned, genome_col, alphabet_size=None, mode="bits",
                       small_sample_correction=False):
    """Per-column, per-token heights. Gaps count in the entropy but are not drawn."""
    n = aligned[genome_col].nunique()
    tokens = aligned["token"].dropna().unique()
    S = alphabet_size or (len(tokens) + 1)  # +1 for the gap symbol
    recs = []
    for pos, col in aligned.groupby("aligned_pos"):
        counts = col["token"].fillna("<gap>").value_counts()
        p = counts / n
        H = float(-(p * np.log2(p)).sum())
        # Schneider small-sample correction. With n ~ 10 and ~15 families it removes ~1 bit
        # per column and clamps most columns to 0 (= invisible bars), so it is opt-in.
        e_n = (S - 1) / (2 * np.log(2) * n) if small_sample_correction else 0.0
        R = max(np.log2(S) - (H + e_n), 0.0)
        for tok, c in counts.items():
            if tok == "<gap>":
                continue
            sub = col[col["token"] == tok]
            f_plus = float((sub["strand"] == "+").mean())
            height = (c / n) * R if mode == "bits" else c / n
            recs.append(dict(aligned_pos=pos, token=tok, count=int(c), proportion=c / n,
                             column_bits=R, height=height, frac_plus=f_plus,
                             colour=sub["colour"].dropna().iloc[0] if sub["colour"].notna().any() else "#808080"))
    return pd.DataFrame(recs), S


# ============================================================
# 5. Colours
# ============================================================

def get_family_colors(df, family_col, color_col=None, colors=None):
    """Consistent family -> colour mapping (explicit dict, then colour column, then tab20)."""
    fams = df[family_col].dropna().unique()
    out = dict(colors) if colors else {}
    if color_col is not None and color_col in df.columns:
        for _, r in df[[family_col, color_col]].dropna().iterrows():
            if r[family_col] not in out:
                out[r[family_col]] = parse_color_and_hatch(r[color_col])[0]
    cmap = plt.get_cmap("tab20")
    for i, f in enumerate(fams):
        out.setdefault(f, mcolors.to_hex(cmap(i % 20)))
    return out


# ============================================================
# 6. Plot
# ============================================================

def plot_microsynteny_logo(
    cluster_label,
    clustering_algo,
    df,
    output_dir,
    goi_type=None,
    svg=False,
    colors=None,
    genome_col="genome",
    rel_pos_col="rel_pos",
    family_col="gene_family",
    color_col="cluster_color",
    max_flanking_genes=None,
    figsize=None,
    dpi=300,
    show=False,
    logo_mode="bits",            # "bits" or "proportion"
    strand_mode="split",         # "split" (same strand up / opposite down) or "ignore"
    drop_unannotated=False,      # remove white / no-family genes
    strip_orientation=True,      # align on family name only
    min_column_occupancy=0.0,    # 0 = show every column
    max_columns_per_side=None,   # None = show every column
    min_family_fraction=0.0,     # 0 = keep every family
    alphabet_size=None,          # None = number of distinct families + gap
    label_genes=False,           # write family names inside tall bars (legend is used by default)
    match_score=3,
    mismatch_score=-0.9,         # > 2*gap_score lets different families share a column
    gap_score=-1,
    label_min_height=0.06,       # fraction of the y-range needed for a label
    max_legend_items=40,
    save_table=True,
    small_sample_correction=False,  # Schneider correction (hides low-information bars)
    min_bar_frac=0.03,           # every bar is drawn at least this fraction of the tallest column
    suffix="B",                  # file-name suffix: ..._B.png (None/"" = no suffix)
    panel_letter="B",            # bold letter in the figure corner (None = off)
):
    os.makedirs(output_dir, exist_ok=True)
    aligned = align_microsynteny_neighborhoods(
        df, genome_col, rel_pos_col, family_col, max_flanking_genes=max_flanking_genes,
        color_col=color_col, strip_orientation=strip_orientation, drop_unannotated=drop_unannotated,
        min_family_fraction=min_family_fraction,
        match_score=match_score, mismatch_score=mismatch_score, gap_score=gap_score)
    if aligned.empty:
        logger.warning("No aligned neighbourhoods for cluster %s.", cluster_label)
        return None
    n = aligned[genome_col].nunique()

    # hide rare columns (always keep the GOI)
    occ = aligned.groupby("aligned_pos")["token"].apply(lambda s: s.notna().mean())
    cand = occ[occ >= min_column_occupancy]
    keep_set = {0}
    for side in (cand[cand.index < 0], cand[cand.index > 0]):
        if max_columns_per_side:
            order = sorted(side.index, key=lambda p: (-side[p], abs(p)))
            side = side.loc[order[:max_columns_per_side]]
        keep_set |= set(side.index)
    keep = sorted(keep_set)
    aligned = aligned[aligned["aligned_pos"].isin(keep)].copy()
    remap = {p: i for i, p in enumerate(keep)}
    goi_x = remap[0]
    aligned["x"] = aligned["aligned_pos"].map(remap)

    info, S = column_information(aligned, genome_col, alphabet_size, logo_mode,
                                 small_sample_correction=small_sample_correction)
    info["x"] = info["aligned_pos"].map(remap)
    fam_colors = get_family_colors(aligned.rename(columns={"colour": "_c"}), "token", "_c", colors)

    # layout
    xs = list(range(len(keep)))
    if figsize is None:
        figsize = (max(8, len(xs) * 0.75), 6)
    fig, ax = plt.subplots(figsize=figsize)

    # Smallest height we draw: bars with ~0 information (diverse / sparse columns) would
    # otherwise be invisible. Only the drawing is floored; the saved table keeps true values.
    ref = info.groupby("x")["height"].sum().max() if len(info) else 0.0
    floor = min_bar_frac * (ref if ref and ref > 0 else 1.0)

    ymax_up = ymax_dn = 0.0
    seg = []
    for x, col in info.groupby("x"):
        is_goi = (x == goi_x)
        up, down = [], []
        for _, r in col.iterrows():
            h = r["height"]
            if is_goi or strand_mode == "ignore":
                up.append((r["token"], max(h, floor), r))
            else:
                if r["frac_plus"] > 0:
                    up.append((r["token"], max(h * r["frac_plus"], floor), r))
                if r["frac_plus"] < 1:
                    down.append((r["token"], max(h * (1 - r["frac_plus"]), floor), r))
        for side, items in (("up", up), ("down", down)):
            items.sort(key=lambda t: t[1])            # smallest next to the axis
            base = 0.0
            for tok, h, r in items:
                seg.append((x, tok, base, h, side, r))
                base += h
            if side == "up":
                ymax_up = max(ymax_up, base)
            else:
                ymax_dn = max(ymax_dn, base)
    ymax_up = max(ymax_up, 1e-6) * 1.08
    ymax_dn = ymax_dn * 1.08
    two_sided = (strand_mode == "split")
    # no opposite-strand genes -> do not waste half of the plot on an empty mirror
    y_lo = -max(ymax_dn, 0.06 * ymax_up) if two_sided else 0

    for x, tok, base, h, side, r in seg:
        base_color, hatch = parse_color_and_hatch(r["colour"])
        base_color = fam_colors.get(tok, base_color) if not isinstance(r["colour"], str) else base_color
        if _is_white(base_color):               # unannotated genes: white would be invisible
            base_color, hatch = UNANNOTATED_GREY, None
        sign = 1 if side == "up" else -1
        kw = dict(x=x, height=sign * h, bottom=sign * base, width=0.85, color=base_color,
                  edgecolor="white", linewidth=0.3)
        if hatch:
            kw["hatch"] = hatch
        ax.bar(**kw)
        if label_genes and h >= label_min_height * ymax_up:
            ax.text(x, sign * (base + h / 2), tok[:11], ha="center", va="center",
                    fontsize=6, color=_text_color(base_color), clip_on=True,
                    bbox=dict(boxstyle="round,pad=0.12", fc=base_color, ec="none", alpha=0.9))

    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(goi_x, color="black", ls="--", lw=1, alpha=0.6)
    ax.set_ylim(y_lo, ymax_up)
    ax.set_xticks(xs)
    ax.set_xticklabels(["GOI" if x == goi_x else f"{x - goi_x:+d}" for x in xs], rotation=45, ha="right")
    ax.set_xlabel("Aligned gene position relative to GOI")
    unit = "Information content (bits)" if logo_mode == "bits" else "Proportion of neighbourhoods"
    ax.set_ylabel(unit + ("\nsame strand as GOI (up) / opposite (down)" if two_sided else ""))
    ticks = ax.get_yticks()
    ax.set_yticks(ticks); ax.set_yticklabels([f"{abs(t):g}" for t in ticks]); ax.set_ylim(y_lo, ymax_up)

    title = f"Microsynteny logo - {clustering_algo} - cluster {cluster_label}"
    if goi_type is not None:
        title += f" - {goi_type}"
    title += f"  (n = {n})"
    ax.set_title(title)

    # legend: families ordered by how many neighbourhoods carry them
    freq = aligned.dropna(subset=["token"]).groupby("token")[genome_col].nunique().sort_values(ascending=False)

    def _first_colour(f):
        sub = aligned[aligned["token"] == f]["colour"].dropna()
        return sub.iloc[0] if len(sub) else None

    annotated = [f for f in freq.index if not _is_white(_first_colour(f))]
    handles = []
    for fam in annotated[:max_legend_items]:
        c = _first_colour(fam)
        bc, hv = parse_color_and_hatch(c) if c is not None else (fam_colors[fam], None)
        kw = dict(facecolor=bc, edgecolor="white", linewidth=0.3, label=str(fam))
        if hv:
            kw["hatch"] = hv
        handles.append(Patch(**kw))
    if len(annotated) > max_legend_items:
        handles.append(Patch(facecolor="#CCCCCC", label=f"+{len(annotated) - max_legend_items} more families"))
    if len(annotated) < len(freq):
        handles.append(Patch(facecolor=UNANNOTATED_GREY, label="Unannotated"))
    ax.legend(handles=handles, title="Gene family", bbox_to_anchor=(1.02, 1),
              loc="upper left", frameon=False, fontsize="small")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    if panel_letter:
        fig.text(0.005, 0.995, panel_letter, fontsize=18, fontweight="bold", va="top", ha="left")
    base = os.path.join(output_dir, f"microsynteny_{clustering_algo}_cluster_{cluster_label}")
    fig_base = base + (f"_{suffix}" if suffix else "")
    fig.savefig(fig_base + ".png", dpi=dpi, bbox_inches="tight")
    if svg:
        fig.savefig(fig_base + ".svg", format="svg", bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    if save_table:
        info.drop(columns=["x"]).to_csv(base + "_columns.csv", index=False)
    logger.info("Microsynteny logo saved to %s.png", fig_base)
    return info
