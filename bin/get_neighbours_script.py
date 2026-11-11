
# (C) 2024, Maria Schreiber, MIT license
""" Revised genomic neighbourhood extraction script. """

import argparse
import logging
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from Bio.Seq import Seq
from sugar import read, read_fts


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


OUTPUT_FILES = {
    "protein": ".protein.mfaa",
    "ncrna": ".ncrna.mfna",
    "nucleotide": ".nucleotide.mfna",
}

PROMOTER_SUFFIX = ".promoter.mfna"

FEATURE_TO_OUTPUT_TYPE = {
    "protein_coding": "protein",
    "CDS": "protein",
    "pseudogene": "nucleotide",
    "tRNA": "ncrna",
    "rRNA": "ncrna",
    "ncRNA": "ncrna",
    "infernal": "ncrna",
    "blastn": "ncrna",
    "misc_RNA": "ncrna",
    "RNase_P_RNA": "ncrna",
}

BLAST_FMT = (
    "qseqid sseqid bitscore evalue pident length mismatch gapopen "
    "qstart qend qlen sstart send sstrand slen qseq sseq"
)


# ---------------------------------------------------------------------------
# Structured internal data
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FeatureKey:
    """Stable biological feature identity."""

    seqid: str
    start: int
    stop: int
    strand: str
    feature_type: str
    annotation_id: str = ""

    def as_identifier(self) -> str:
        """Create a deterministic identifier suitable for FASTA headers."""
        parts = [
            self.seqid,
            str(self.start),
            str(self.stop),
            self.strand,
            self.feature_type,
        ]

        if self.annotation_id:
            parts.append(self.annotation_id)

        return _safe_identifier(":".join(parts))


@dataclass(frozen=True)
class AnnotationEntry:
    """One annotation record belonging to a neighbourhood."""

    identifier: str
    name: str
    start: int
    stop: int
    orientation: str
    biotype: str
    member_id: str = ""
    genome_id: str = ""

    def to_tsv(self) -> str:
        return "\t".join(
            (
                self.member_id,
                self.identifier,
                self.name,
                str(self.start),
                str(self.stop),
                self.orientation,
                self.biotype,
                self.genome_id,
            )
        )

    def to_fasta_header(self) -> str:
        return _safe_identifier(
            "\t".join(
                (
                    self.identifier,
                    self.name,
                    str(self.start),
                    str(self.stop),
                    self.orientation,
                    self.biotype,
                )
            )
        )


# ---------------------------------------------------------------------------
# Argument parsing and validation
# ---------------------------------------------------------------------------

def validate_neighbours(neighbours: str) -> tuple[int, int]:
    """Validate x,y gene counts or x:y nucleotide distances."""
    if neighbours.count(",") + neighbours.count(":") != 1:
        raise ValueError(
            "Invalid neighbours format. Use x,y for gene counts "
            "or x:y for nucleotide distances."
        )

    separator = "," if "," in neighbours else ":"
    parts = neighbours.split(separator)

    try:
        up, down = map(int, parts)
    except ValueError as exc:
        raise ValueError(
            "Neighbour values must be integers."
        ) from exc

    if up < 0 or down < 0:
        raise ValueError("Neighbour values must be non-negative.")

    return up, down


def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract genomic neighbourhoods and sequences around "
            "BLAST/Infernal hits or selected GFF features."
        )
    )

    parser.add_argument("--hit_file", help="BLAST or Infernal hit file.")

    parser.add_argument(
        "--input_type",
        required=True,
        choices=["blastn", "blastp", "infernal", "tblastn", "from_gff"],
    )

    parser.add_argument(
        "--bio_type",
        required=True,
        choices=["ncRNA", "protein"],
        help="Biological type of the target/query.",
    )

    parser.add_argument("--gff_file", required=True)
    parser.add_argument("--fna_file", required=True)

    parser.add_argument(
        "--gene_of_interest",
        required=True,
        help="Target query ID, feature ID, product name, or feature name.",
    )

    parser.add_argument(
        "--neighbours",
        default="4,4",
        help="Neighbour range: x,y genes or x:y nucleotides.",
    )

    parser.add_argument("--output_path", required=True)

    parser.add_argument(
        "--promoter",
        choices=["yes", "no"],
        default="no",
    )

    parser.add_argument("--promoter_len", type=int, default=100)

    parser.add_argument(
        "--promoter_mode",
        choices=["fixed", "to_next_gene"],
        default="fixed",
        help=(
            "Use a fixed promoter length or extend to the nearest "
            "upstream feature boundary."
        ),
    )

    parser.add_argument(
        "--including_features",
        nargs="+",
        default=[
            "CDS",
            "pseudogene",
            "ncRNA",
            "rRNA",
            "tRNA",
            "infernal",
        ],
        help="GFF feature types to include.",
    )

    parser.add_argument(
        "--filter_by",
        choices=["evalue", "bitscore", "none"],
        default="evalue",
    )

    parser.add_argument(
        "--evalue_threshold_blast",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--evalue_threshold_infernal",
        type=float,
        default=0.1,
    )

    parser.add_argument(
        "--bitscore_threshold",
        type=float,
        default=50.0,
        help=(
            "Minimum BLAST bitscore or Infernal score when "
            "--filter_by bitscore is selected."
        ),
    )

    parser.add_argument("--len_threshold_blast", type=int, default=40)
    parser.add_argument("--len_threshold_infernal", type=int, default=20)

    parser.add_argument(
        "--from_gff_feature",
        default="ID",
        help="GFF metadata attribute used to identify the target.",
    )

    parser.add_argument(
        "--include_all_overlaps",
        action="store_true",
        help="Include all overlapping features as neighbours.",
    )

    parser.add_argument(
        "--overlap_threshold",
        type=float,
        default=0.75,
        help=(
            "Minimum overlap fraction required for same-strand "
            "features to be retained."
        ),
    )

    parser.add_argument(
        "--substring_search",
        action="store_true",
        help="Match gene_of_interest as a substring of the selected attribute.",
    )
    
    parser.add_argument(
        "--genome_id",
        default=None,
        help="Genome accession to write in the output TSV, e.g. GCF_000022005.1",
    )
    
    return parser


def validate_args(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
) -> None:
    try:
        validate_neighbours(args.neighbours)
    except ValueError as exc:
        parser.error(str(exc))

    if args.input_type != "from_gff" and not args.hit_file:
        parser.error(
            f"--hit_file is required for input_type={args.input_type}"
        )

    if not args.including_features:
        parser.error("--including_features must contain at least one feature type.")

    if not 0.0 <= args.overlap_threshold <= 1.0:
        parser.error("--overlap_threshold must be between 0 and 1.")

    if args.promoter_len < 0:
        parser.error("--promoter_len must be non-negative.")

    if args.len_threshold_blast < 0 or args.len_threshold_infernal < 0:
        parser.error("Length thresholds must be non-negative.")

    if (
        args.evalue_threshold_blast < 0
        or args.evalue_threshold_infernal < 0
    ):
        parser.error("E-value thresholds must be non-negative.")

    if args.bitscore_threshold < 0:
        parser.error("--bitscore_threshold must be non-negative.")

    if not math.isfinite(args.overlap_threshold):
        parser.error("--overlap_threshold must be finite.")

    for label, path in (
        ("GFF", args.gff_file),
        ("FASTA", args.fna_file),
    ):
        if not Path(path).is_file():
            parser.error(f"{label} input file does not exist: {path}")

    if args.input_type != "from_gff" and not Path(args.hit_file).is_file():
        parser.error(f"Hit file does not exist: {args.hit_file}")


# ---------------------------------------------------------------------------
# Metadata helpers
# ---------------------------------------------------------------------------

def _get_meta(feature: Any, fmt: Optional[str]) -> Any:
    if not fmt:
        return None

    return getattr(feature.meta, f"_{fmt}", None)


def _meta_value(
    meta: Any,
    *names: str,
    default: Any = None,
) -> Any:
    if meta is None:
        return default

    for name in names:
        if isinstance(meta, dict):
            value = meta.get(name)
        else:
            value = getattr(meta, name, None)

        if value is not None:
            return value

    return default


def _as_text(value: Any) -> str:
    """Normalize metadata values, including list-valued attributes."""
    if isinstance(value, (list, tuple)):
        return " ".join(str(item) for item in value)

    return str(value) if value is not None else ""


def _safe_identifier(value: Any) -> str:
    """Make a deterministic, single-token identifier."""
    return re.sub(r"[\s\t\r\n]+", "_", str(value))


def _feature_metadata(feature: Any) -> Any:
    fmt = getattr(feature.meta, "_fmt", None)
    return _get_meta(feature, fmt)


def _feature_annotation_id(feature: Any) -> str:
    """Return a source annotation ID where available."""
    for fmt in ("gff", "blast", "infernal"):
        meta = _get_meta(feature, fmt)
        value = _meta_value(meta, "ID", "id", "feature_id")

        if value:
            return _as_text(value)

    qualifiers = getattr(feature, "qualifiers", None)

    if isinstance(qualifiers, dict):
        for name in ("ID", "id", "feature_id"):
            if qualifiers.get(name):
                return _as_text(qualifiers[name])

    return ""


def _feature_key(feature: Any) -> FeatureKey:
    return FeatureKey(
        seqid=str(feature.seqid),
        start=int(feature.loc.start),
        stop=int(feature.loc.stop),
        strand=str(feature.loc.strand or "."),
        feature_type=str(getattr(feature, "type", None) or "unknown"),
        annotation_id=_feature_annotation_id(feature),
    )


def _feature_length(feature: Any) -> int:
    return max(0, int(feature.loc.stop) - int(feature.loc.start))


def _numeric_value(value: Any, label: str, feature: Any) -> float:
    try:
        number = float(value)

        if not math.isfinite(number):
            raise ValueError("Value is not finite.")

        return number

    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Could not read numeric {label} for feature "
            f"{feature.seqid}:{feature.loc.start}-{feature.loc.stop}."
        ) from exc

def _get_genome_id(args: argparse.Namespace) -> str:
    """Get an explicitly supplied genome ID or infer it from the input filename."""

    if args.genome_id and args.genome_id.strip():
        return args.genome_id.strip()

    source_path = args.hit_file if args.hit_file else args.gff_file
    filename = Path(source_path).name

    if filename.lower().endswith(".tsv"):
        return filename[:-4]

    return Path(filename).stem

# ---------------------------------------------------------------------------
# Score handling
# ---------------------------------------------------------------------------

def _feature_evalue(feature: Any, input_type_key: str) -> float:
    fmt = "blast" if input_type_key == "blast" else "infernal"
    meta = _get_meta(feature, fmt)
    value = _meta_value(meta, "evalue", "e_value")

    if value is None:
        raise ValueError(
            f"No E-value found in {fmt} metadata for feature "
            f"{feature.seqid}:{feature.loc.start}-{feature.loc.stop}."
        )

    return _numeric_value(value, "E-value", feature)


def _feature_bitscore(
    feature: Any,
    input_type_key: Optional[str] = None,
) -> Optional[float]:
    if input_type_key == "blast":
        formats = ("blast",)
    elif input_type_key == "infernal":
        formats = ("infernal",)
    else:
        formats = ("blast", "infernal")

    for fmt in formats:
        meta = _get_meta(feature, fmt)

        value = _meta_value(meta, "bitscore", "bit_score", "score")

        if value is not None:
            try:
                number = float(value)
                if math.isfinite(number):
                    return number
            except (TypeError, ValueError):
                pass

    qualifiers = getattr(feature, "qualifiers", None)

    if isinstance(qualifiers, dict):
        for key in ("bitscore", "bit_score", "score"):
            value = qualifiers.get(key)

            if isinstance(value, (list, tuple)):
                value = value[0] if value else None

            if value is not None:
                try:
                    number = float(value)
                    if math.isfinite(number):
                        return number
                except (TypeError, ValueError):
                    pass

    return None


def _feature_quality(
    feature: Any,
    input_type_key: str,
) -> Optional[float]:
    """Higher is better."""
    if input_type_key == "blast":
        return _feature_bitscore(feature, input_type_key)

    return -_feature_evalue(feature, input_type_key)


def _score_sort_key(
    feature: Any,
    input_type_key: str,
) -> tuple:
    """Sort best-first when a score exists, with deterministic tie-breaking."""
    quality = _feature_quality(feature, input_type_key)

    return (
        quality is None,
        -(quality if quality is not None else 0.0),
        str(feature.seqid),
        int(feature.loc.start),
        int(feature.loc.stop),
        str(feature.loc.strand),
    )


def process_hits(
    hit_file: str,
    input_type: str,
    gene_of_interest: str,
    args: argparse.Namespace,
) -> list:
    """Read and filter BLAST/Infernal hits."""
    input_type_key = (
        "blast"
        if input_type in {"blastn", "blastp", "tblastn"}
        else "infernal"
    )

    if input_type_key == "blast":
        features = read_fts(
            hit_file,
            "blast",
            outfmt=BLAST_FMT,
            ftype=gene_of_interest,
        )
        min_length = args.len_threshold_blast
        max_evalue = args.evalue_threshold_blast
    else:
        features = read_fts(
            hit_file,
            "infernal",
            ftype=gene_of_interest,
        )
        min_length = args.len_threshold_infernal
        max_evalue = args.evalue_threshold_infernal

    filtered = []

    for feature in features:
        if _feature_length(feature) <= min_length:
            continue

        if args.filter_by == "evalue":
            evalue = _feature_evalue(feature, input_type_key)
            if evalue > max_evalue:
                continue

        elif args.filter_by == "bitscore":
            score = _feature_bitscore(feature, input_type_key)

            if score is None:
                raise ValueError(
                    f"No numeric score found for {input_type_key} hit "
                    f"{feature.seqid}:{feature.loc.start}-{feature.loc.stop}."
                )

            if score < args.bitscore_threshold:
                continue

        filtered.append(feature)

    filtered.sort(key=lambda ft: _score_sort_key(ft, input_type_key))

    verified = []

    for feature in filtered:
        if _hit_matches_target(feature, gene_of_interest):
            verified.append(feature)
        else:
            meta = _get_meta(feature, "blast") or _get_meta(
                feature, "infernal"
            )
            candidates = _hit_query_names(feature, meta)

            if candidates:
                logger.warning(
                    "Skipping hit whose query identifier %r does not match "
                    "requested target %r.",
                    candidates,
                    gene_of_interest,
                )

    if not filtered:
        logger.warning("No hits remain after filtering.")

    logger.info(
        "Hits retained after filtering and target verification: %d "
        "(type=%s, filter_by=%s)",
        len(verified),
        input_type,
        args.filter_by,
    )

    return verified


def _hit_query_names(feature: Any, meta: Any = None) -> list[str]:
    if meta is None:
        meta = _get_meta(feature, "blast") or _get_meta(feature, "infernal")

    names = []

    for key in ("query", "qseqid", "query_id", "query_name"):
        value = _meta_value(meta, key)

        if value is not None:
            text = _as_text(value)
            if text:
                names.append(text)

    qualifiers = getattr(feature, "qualifiers", None)

    if isinstance(qualifiers, dict):
        for key in ("query", "qseqid", "query_id", "query_name"):
            value = qualifiers.get(key)

            if value is not None:
                text = _as_text(value)
                if text:
                    names.append(text)

    return list(dict.fromkeys(names))


def _hit_matches_target(feature: Any, gene_of_interest: str) -> bool:
    """Check query names when available; otherwise trust ftype filtering."""
    names = _hit_query_names(feature)
    return not names or gene_of_interest in names


def remove_overlapping_features(features: list, input_type: str) -> list:
    """Remove overlapping hits, retaining the highest-quality hit."""
    if not features:
        return []

    input_type_key = "infernal" if input_type == "infernal" else "blast"

    ranked = sorted(
        features,
        key=lambda ft: _score_sort_key(ft, input_type_key),
    )

    selected = []

    for feature in ranked:
        overlaps_existing = any(
            feature.seqid == other.seqid
            and int(feature.loc.start) < int(other.loc.stop)
            and int(other.loc.start) < int(feature.loc.stop)
            for other in selected
        )

        if not overlaps_existing:
            selected.append(feature)

    selected.sort(
        key=lambda ft: (
            str(ft.seqid),
            int(ft.loc.start),
            int(ft.loc.stop),
        )
    )

    logger.info("Non-overlapping hits retained: %d", len(selected))
    return selected


# ---------------------------------------------------------------------------
# Feature labels and annotation rows
# ---------------------------------------------------------------------------

def get_orientation(strand: str) -> str:
    return "antisense" if strand == "-" else "sense"


def _find_feature_name(
    meta: Any,
    priority_attrs: tuple[str, ...],
    default_name: str = "",
) -> str:
    for attr in priority_attrs:
        name = _meta_value(meta, attr)
        if name:
            return _as_text(name)

    return default_name


def _get_biotype(
    feature: Any,
    meta: Any,
    args: argparse.Namespace,
) -> str:
    if args.input_type == "tblastn" and _get_meta(feature, "blast") is not None:
        return "tblastn"

    if args.input_type == "blastn" and _get_meta(feature, "blast") is not None:
        return "blastn"

    if args.input_type == "blastp" and _get_meta(feature, "blast") is not None:
        return "blastp"

    if args.input_type == "infernal" and _get_meta(feature, "infernal") is not None:
        return "infernal"

    return (
        _as_text(_meta_value(meta, "gene_biotype"))
        or str(getattr(feature, "type", None) or "")
        or _as_text(_meta_value(meta, "_fmt"))
        or "unknown"
    )


def generate_entry(
    feature: Any,
    args: argparse.Namespace,
) -> AnnotationEntry:
    """Build a structured, stable annotation record."""
    fmt = getattr(feature.meta, "_fmt", None)
    meta = _get_meta(feature, fmt) if fmt else None

    if args.input_type == "from_gff":
        name_priority = (
            args.from_gff_feature,
            "Name",
            "ID",
            "product",
        )
    else:
        name_priority = ("query", "qseqid", "Name", "ID")

    name = _find_feature_name(meta, name_priority, default_name="unnamed")

    if args.input_type == "from_gff" and args.from_gff_feature == "product":
        product = _meta_value(meta, "product")
        if product and args.gene_of_interest in str(product):
            name = args.gene_of_interest

    biotype = _get_biotype(feature, meta, args)
    key = _feature_key(feature)

    return AnnotationEntry(
        identifier=key.as_identifier(),
        name=name,
        start=key.start,
        stop=key.stop,
        orientation=get_orientation(key.strand),
        biotype=biotype,
    )


# ---------------------------------------------------------------------------
# Feature overlap and neighbourhood selection
# ---------------------------------------------------------------------------

def check_overlap(
    center: Any,
    overlaps: list,
    overlap_threshold: float,
) -> list:
    """Retain sufficiently overlapping same-strand features and opposite-strand overlaps."""
    retained = []

    c_start = int(center.loc.start)
    c_stop = int(center.loc.stop)
    c_length = c_stop - c_start

    if c_length <= 0:
        return retained

    for gene in overlaps:
        g_start = int(gene.loc.start)
        g_stop = int(gene.loc.stop)
        g_length = g_stop - g_start

        if g_length <= 0:
            continue

        overlap_length = min(c_stop, g_stop) - max(c_start, g_start)

        if overlap_length <= 0:
            continue

        center_fraction = overlap_length / c_length
        gene_fraction = overlap_length / g_length

        if (
            center_fraction >= overlap_threshold
            and gene_fraction >= overlap_threshold
        ):
            retained.append(gene)
        elif gene.loc.strand != center.loc.strand:
            retained.append(gene)

    return retained



def get_neighbor_features(
    features: list,
    center: Any,
    neighbours: str,
    overlap_threshold: float,
    include_all_overlaps: bool = False,
) -> list:
    """
    Find neighbours around a centre feature.

    For count mode (e.g. --neighbours 3,3), strand orientation is ignored:
    select the N closest features by genomic coordinates on each side.

    Circular wraparound is not applied in count mode.
    """
    same_contig = [
        ft for ft in features
        if ft.seqid == center.seqid
    ]

    c_start = int(center.loc.start)
    c_stop = int(center.loc.stop)

    # Features strictly before or after the centre.
    before = [
        ft for ft in same_contig
        if ft is not center and int(ft.loc.stop) <= c_start
    ]

    after = [
        ft for ft in same_contig
        if ft is not center and int(ft.loc.start) >= c_stop
    ]

    # Nearest first, independent of strand.
    before.sort(
        key=lambda ft: (
            -int(ft.loc.stop),
            -int(ft.loc.start),
        )
    )
    after.sort(
        key=lambda ft: (
            int(ft.loc.start),
            int(ft.loc.stop),
        )
    )

    if "," in neighbours:
        left_count, right_count = map(int, neighbours.split(","))

        # Find overlapping features separately.
        overlaps = [
            ft for ft in same_contig
            if ft is not center
            and int(ft.loc.start) < c_stop
            and c_start < int(ft.loc.stop)
        ]

        retained_overlaps = (
            overlaps
            if include_all_overlaps
            else check_overlap(center, overlaps, overlap_threshold)
        )

        result = (
            retained_overlaps
            + before[:left_count]
            + after[:right_count]
            + [center]
        )

    else:
        # Preserve the existing base-pair window mode.
        upstream_nt, downstream_nt = map(int, neighbours.split(":"))

        start_pos = max(0, c_start - upstream_nt)
        stop_pos = c_stop + downstream_nt

        result = [
            ft for ft in same_contig
            if int(ft.loc.start) < stop_pos
            and int(ft.loc.stop) > start_pos
        ]

        if not any(ft is center for ft in result):
            result.append(center)

    # Deduplicate and return in genomic order.
    unique = {}
    for ft in result:
        unique[_feature_key(ft)] = ft

    return sorted(
        unique.values(),
        key=lambda ft: (
            str(ft.seqid),
            int(ft.loc.start),
            int(ft.loc.stop),
            str(ft.loc.strand),
            _feature_annotation_id(ft),
        ),
    )

def _passes_same_strand_overlap(
    hit: Any,
    feature: Any,
    overlap_threshold: float,
) -> bool:
    """Return True if two features have the same strand and
    overlap by more than the threshold relative to both lengths.
    """
    if str(hit.seqid) != str(feature.seqid):
        return False

    if str(hit.loc.strand) != str(feature.loc.strand):
        return False

    hit_start = int(hit.loc.start)
    hit_stop = int(hit.loc.stop)
    feature_start = int(feature.loc.start)
    feature_stop = int(feature.loc.stop)

    hit_length = hit_stop - hit_start
    feature_length = feature_stop - feature_start

    if hit_length <= 0 or feature_length <= 0:
        return False

    overlap = (
        min(hit_stop, feature_stop)
        - max(hit_start, feature_start)
    )

    if overlap <= 0:
        return False

    hit_fraction = overlap / hit_length
    feature_fraction = overlap / feature_length

    # Strictly greater than the threshold, not greater than or equal.
    return (
        hit_fraction > overlap_threshold
        and feature_fraction > overlap_threshold
    )

def _find_overlapping_annotation(
    hit: Any,
    gff_features: list,
    overlap_threshold: Optional[float] = None,
    require_same_strand: bool = False,
) -> Optional[Any]:
    """Map a hit to its best-overlapping GFF feature.

    When require_same_strand is True, only features on the same
    strand that exceed the reciprocal-overlap threshold are eligible.
    """
    candidates = []

    hit_start = int(hit.loc.start)
    hit_stop = int(hit.loc.stop)
    hit_length = hit_stop - hit_start

    if hit_length <= 0:
        return None

    for feature in gff_features:
        if feature.seqid != hit.seqid:
            continue

        start = int(feature.loc.start)
        stop = int(feature.loc.stop)
        length = stop - start

        if length <= 0:
            continue

        overlap = min(hit_stop, stop) - max(hit_start, start)

        if overlap <= 0:
            continue

        if require_same_strand:
            if overlap_threshold is None:
                raise ValueError(
                    "overlap_threshold is required for strict matching."
                )

            if not _passes_same_strand_overlap(
                hit, feature, overlap_threshold
            ):
                continue

        reciprocal_overlap = min(
            overlap / hit_length,
            overlap / length,
        )

        candidates.append(
            (
                -reciprocal_overlap,
                abs(length - hit_length),
                start,
                stop,
                _feature_annotation_id(feature),
                feature,
            )
        )

    if not candidates:
        return None

    candidates.sort(key=lambda item: item[:-1])
    return candidates[0][-1]
# ---------------------------------------------------------------------------
# Sequence validation and extraction
# ---------------------------------------------------------------------------

def _get_sequence(seqs: Any, seqid: str) -> Any:
    sequence = seqs.d.get(seqid)

    if sequence is None:
        raise KeyError(
            f"Sequence ID {seqid!r} was not found in the FASTA file. "
            "Check that GFF seqids and FASTA record IDs match."
        )

    return sequence


def _sequence_to_string(sequence: Any) -> str:
    return str(sequence)


def _reverse_complement(sequence: Any) -> Seq:
    return Seq(str(sequence)).reverse_complement()


def _feature_phase(feature: Any) -> Optional[int]:
    """Read phase/frame metadata where available."""
    for fmt in ("gff", None):
        meta = _get_meta(feature, fmt) if fmt else None

        value = _meta_value(meta, "phase", "frame", "codon_start")

        if value is not None:
            try:
                phase = int(value)

                if isinstance(meta, dict) and "codon_start" in meta:
                    phase -= 1

                if phase in (0, 1, 2):
                    return phase
            except (TypeError, ValueError):
                pass

    qualifiers = getattr(feature, "qualifiers", None)

    if isinstance(qualifiers, dict):
        for key in ("codon_start", "phase", "frame"):
            value = qualifiers.get(key)

            if isinstance(value, (list, tuple)):
                value = value[0] if value else None

            if value is not None:
                try:
                    phase = int(value)

                    if key == "codon_start":
                        phase -= 1

                    if phase in (0, 1, 2):
                        return phase
                except (TypeError, ValueError):
                    pass

    return None


def _extract_genomic_feature(
    feature: Any,
    sequence: Any,
    circular_seqids: Optional[set[str]] = None,
) -> Seq:
    """Extract a feature, including a single circular-origin wrap."""
    circular_seqids = circular_seqids or set()

    location = getattr(feature, "location", None)

    if location is not None and hasattr(location, "extract"):
        try:
            return Seq(str(location.extract(Seq(str(sequence)))))
        except (TypeError, ValueError, AttributeError):
            logger.debug(
                "Compound-location extraction unavailable; "
                "using sugar feature coordinates."
            )

    start = int(feature.loc.start)
    stop = int(feature.loc.stop)
    sequence_length = len(sequence)
    seqid = str(feature.seqid)

    if start < 0 or stop <= start:
        raise ValueError(f"Invalid interval {seqid}:{start}-{stop}.")

    if stop <= sequence_length:
        extracted = Seq(str(sequence[start:stop]))
    elif (
        seqid in circular_seqids
        and start < sequence_length
        and stop - start <= sequence_length
    ):
        extracted = Seq(
            str(sequence[start:sequence_length])
            + str(sequence[0:stop - sequence_length])
        )
    else:
        raise ValueError(
            f"Invalid interval {seqid}:{start}-{stop} "
            f"for sequence length {sequence_length}; "
            f"circular={seqid in circular_seqids}."
        )

    if feature.loc.strand == "-":
        extracted = extracted.reverse_complement()

    return extracted


def _translate_cds(feature: Any, nucleotide_sequence: Seq) -> Seq:
    """Translate a CDS after applying a known phase/codon-start offset."""
    phase = _feature_phase(feature)

    if phase is not None:
        nucleotide_sequence = nucleotide_sequence[phase:]

    usable_length = (len(nucleotide_sequence) // 3) * 3
    nucleotide_sequence = nucleotide_sequence[:usable_length]

    return nucleotide_sequence.translate()


# ---------------------------------------------------------------------------
# Transactional output writer
# ---------------------------------------------------------------------------

class OutputWriter:
    """Write output files to temporary paths and replace final files on commit."""

    def __init__(
        self,
        output_path: str,
        circular_seqids: Optional[set[str]] = None,
    ):
        self.base_path = Path(output_path)
        self.circular_seqids = set(circular_seqids or set())
        self.base_path.parent.mkdir(parents=True, exist_ok=True)

        self.final_paths = {
            "tsv": Path(f"{self.base_path}.tsv"),
            "protein": Path(f"{self.base_path}{OUTPUT_FILES['protein']}"),
            "ncrna": Path(f"{self.base_path}{OUTPUT_FILES['ncrna']}"),
            "nucleotide": Path(
                f"{self.base_path}{OUTPUT_FILES['nucleotide']}"
            ),
            "promoter": Path(f"{self.base_path}{PROMOTER_SUFFIX}"),
        }

        self.temp_paths = {
            key: path.with_name(path.name + ".tmp")
            for key, path in self.final_paths.items()
        }

        self.handles: dict[str, Any] = {}

        # Sequence records are written once per feature and biotype.
        self.written_entries: set[tuple] = set()

        # Annotation rows are tracked separately by neighbourhood.
        self.written_annotations: set[tuple] = set()

        self.written_promoters: set[tuple] = set()
        self._committed = False

    def __enter__(self):
        try:
            for key, path in self.temp_paths.items():
                path.parent.mkdir(parents=True, exist_ok=True)
                self.handles[key] = path.open("w", encoding="utf-8")
        except Exception:
            self._close_handles()
            self._remove_temporary_files()
            raise

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._close_handles()

        if exc_type is not None or not self._committed:
            self._remove_temporary_files()

        return False

    def _close_handles(self):
        for handle in self.handles.values():
            if not handle.closed:
                handle.close()
        self.handles.clear()

    def _remove_temporary_files(self):
        for path in self.temp_paths.values():
            try:
                path.unlink(missing_ok=True)
            except OSError:
                logger.warning("Could not remove temporary file %s", path)

    def commit(self):
        """Replace final outputs with completed temporary files."""
        self._close_handles()

        for key, temp_path in self.temp_paths.items():
            os.replace(temp_path, self.final_paths[key])

        self._committed = True

    def write_feature(
        self,
        feature: Any,
        entry: AnnotationEntry,
        seqs: Any,
        input_type: str,
        bio_type: str,
    ):
        """Write neighbourhood membership rows and deduplicated sequences."""
        key = _feature_key(feature)
        biotype = entry.biotype

        # Record each feature's membership in each neighbourhood.
        annotation_key = (entry.member_id, key, biotype)

        if annotation_key not in self.written_annotations:
            self.handles["tsv"].write(entry.to_tsv() + "\n")
            self.written_annotations.add(annotation_key)

        # Avoid writing the same sequence record repeatedly.
        sequence_key = (key, biotype)

        if sequence_key in self.written_entries:
            return

        sequence = _get_sequence(seqs, str(feature.seqid))

        if (
            input_type == "tblastn"
            and _get_meta(feature, "blast") is not None
        ):
            blast_meta = _get_meta(feature, "blast")
            aligned_subject = _meta_value(blast_meta, "sseq")

            if not aligned_subject:
                raise ValueError(
                    "tblastn hit has no subject alignment sequence (sseq)."
                )

            content = Seq(
                str(aligned_subject).replace("_", "").replace("-", "")
            )
            output_type = "protein"

        else:
            content = _extract_genomic_feature(
                feature,
                sequence,
                circular_seqids=self.circular_seqids,
            )

            output_type = FEATURE_TO_OUTPUT_TYPE.get(biotype)

            if output_type is None:
                if biotype in {"blastn", "infernal"}:
                    output_type = "ncrna"
                elif biotype in {"blastp", "tblastn"}:
                    output_type = "protein"
                elif biotype == "unknown":
                    output_type = (
                        "protein" if bio_type == "protein" else "ncrna"
                    )

            if output_type == "protein":
                content = _translate_cds(feature, content)

        if output_type not in OUTPUT_FILES:
            logger.warning(
                "No FASTA output mapping for biotype %r; "
                "annotation row written only.",
                biotype,
            )
            self.written_entries.add(sequence_key)
            return

        self.handles[output_type].write(
            f">{entry.to_fasta_header()}\n{content}\n"
        )
        self.written_entries.add(sequence_key)

    def write_promoter(
        self,
        feature: Any,
        entry: AnnotationEntry,
        sequence: Any,
        start: int,
        stop: int,
    ):
        """Write a promoter interval using half-open coordinates."""
        start = max(0, int(start))
        stop = min(len(sequence), int(stop))

        if start >= stop:
            logger.warning(
                "Skipping empty promoter for %s:%d-%d",
                feature.seqid,
                start,
                stop,
            )
            return

        promoter_key = (_feature_key(feature), start, stop)

        if promoter_key in self.written_promoters:
            return

        promoter_sequence = Seq(str(sequence[start:stop]))

        if feature.loc.strand == "-":
            promoter_sequence = promoter_sequence.reverse_complement()

        promoter_entry = AnnotationEntry(
            identifier=entry.identifier,
            name=entry.name,
            start=start,
            stop=stop,
            orientation=entry.orientation,
            biotype="promoter",
            member_id=entry.member_id,
        )

        self.handles["promoter"].write(
            f">{promoter_entry.to_fasta_header()}\n"
            f"{promoter_sequence}\n"
        )
        self.written_promoters.add(promoter_key)


# ---------------------------------------------------------------------------
# Promoter extraction
# ---------------------------------------------------------------------------

def get_promoter_coords(
    center: Any,
    features: list,
    promoter_len: int,
    mode: str = "fixed",
    sequence_length: Optional[int] = None,
) -> Optional[tuple[int, int]]:
    """Return promoter coordinates as a half-open interval."""
    center_start = int(center.loc.start)
    center_stop = int(center.loc.stop)

    same_contig = [
        ft for ft in features
        if ft.seqid == center.seqid and ft is not center
    ]
    same_contig.sort(key=lambda ft: int(ft.loc.start))

    if center.loc.strand == "+":
        promoter_stop = center_start

        if mode == "fixed":
            promoter_start = max(0, center_start - promoter_len)
        else:
            upstream = [
                ft for ft in same_contig
                if int(ft.loc.stop) <= center_start
            ]
            promoter_start = (
                max(int(ft.loc.stop) for ft in upstream)
                if upstream else 0
            )

        if sequence_length is not None:
            promoter_start = max(0, min(promoter_start, sequence_length))
            promoter_stop = max(0, min(promoter_stop, sequence_length))

        if promoter_start >= promoter_stop:
            return None

        return promoter_start, promoter_stop

    promoter_start = center_stop

    if mode == "fixed":
        promoter_stop = center_stop + promoter_len
    else:
        downstream = [
            ft for ft in same_contig
            if int(ft.loc.start) >= center_stop
        ]
        promoter_stop = (
            min(int(ft.loc.start) for ft in downstream)
            if downstream
            else (
                sequence_length
                if sequence_length is not None
                else center_stop + promoter_len
            )
        )

    if sequence_length is not None:
        promoter_start = max(0, min(promoter_start, sequence_length))
        promoter_stop = max(0, min(promoter_stop, sequence_length))

    if promoter_start >= promoter_stop:
        return None

    return promoter_start, promoter_stop


# ---------------------------------------------------------------------------
# Target identification
# ---------------------------------------------------------------------------

def _is_target_feature(
    feature: Any,
    args: argparse.Namespace,
) -> bool:
    if args.input_type != "from_gff":
        fmt = str(getattr(feature.meta, "_fmt", "")).lower()

        return fmt in {
            "blast",
            "blastn",
            "blastp",
            "tblastn",
            "infernal",
        } and _hit_matches_target(feature, args.gene_of_interest)

    meta = _get_meta(feature, "gff")

    if meta is None:
        return False

    target_value = _as_text(
        _meta_value(meta, args.from_gff_feature, default="")
    )

    if args.substring_search:
        return args.gene_of_interest in target_value

    return args.gene_of_interest == target_value


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def validate_features_against_fasta(
    gff_features: list,
    hit_features: list,
    seqs: Any,
    circular_seqids: Optional[set[str]] = None,
) -> None:
    """Validate contig IDs and intervals before creating output files."""
    circular_seqids = circular_seqids or set()
    fasta_ids = set(seqs.d.keys())

    gff_ids = {str(ft.seqid) for ft in gff_features}
    hit_ids = {str(ft.seqid) for ft in hit_features}

    unmatched_gff = sorted(gff_ids - fasta_ids)
    unmatched_hits = sorted(hit_ids - fasta_ids)

    if unmatched_gff or unmatched_hits:
        messages = [
            "Input contig validation failed: GFF/hit sequence IDs "
            "must match FASTA record IDs."
        ]

        if unmatched_gff:
            messages.append(
                "GFF contigs missing from FASTA: " + ", ".join(unmatched_gff)
            )

        if unmatched_hits:
            messages.append(
                "Hit contigs missing from FASTA: " + ", ".join(unmatched_hits)
            )

        raise ValueError("\n".join(messages))

    invalid_intervals = []

    for label, features in (("GFF", gff_features), ("hit", hit_features)):
        for feature in features:
            seqid = str(feature.seqid)
            start = int(feature.loc.start)
            stop = int(feature.loc.stop)
            sequence_length = len(seqs.d[seqid])
            span = stop - start

            ordinary_interval = 0 <= start < stop <= sequence_length

            valid_circular_wrap = (
                seqid in circular_seqids
                and 0 <= start < sequence_length
                and stop > sequence_length
                and 0 < span <= sequence_length
            )

            if not (ordinary_interval or valid_circular_wrap):
                invalid_intervals.append(
                    f"{label} {seqid}:{start}-{stop} "
                    f"(sequence length {sequence_length}, "
                    f"circular={seqid in circular_seqids})"
                )

    if invalid_intervals:
        raise ValueError(
            "Invalid feature coordinates detected:\n  "
            + "\n  ".join(invalid_intervals)
        )

def passes_existing_validation(feature, accepted_features):
    """Return False if the feature overlaps an already accepted feature."""
    feature_start = int(feature.loc.start)
    feature_stop = int(feature.loc.stop)
    feature_seqid = feature.seqid

    for accepted in accepted_features:
        if accepted.seqid != feature_seqid:
            continue

        accepted_start = int(accepted.loc.start)
        accepted_stop = int(accepted.loc.stop)

        # Coordinates overlap if the intervals intersect.
        if feature_start <= accepted_stop and accepted_start <= feature_stop:
            return False

    return True

def validate_selected_gff_features(raw_gff, selected_features, args):
    validated_features = []

    for feature in selected_features:
        feature_type = str(getattr(feature, "type", "")).lower()

        # Always preserve Infernal and BLAST hits, even when overlapping.
        if feature_type in {"infernal", "blast", "blast_hit"}:
            validated_features.append(feature)
            continue

        # Apply overlap filtering to other feature types.
        if passes_existing_validation(feature, validated_features):
            validated_features.append(feature)

    return validated_features


def get_circular_seqids(raw_gff: list) -> set[str]:
    """Return sequence IDs whose GFF region feature marks them circular."""
    circular = set()

    for feature in raw_gff:
        if str(getattr(feature, "type", "")).lower() == "infernal":
            logger.info(
                "Infernal feature: seqid=%s start=%s stop=%s metadata=%s",
                feature.seqid,
                feature.loc.start,
                feature.loc.stop,
                getattr(feature.meta, "_gff", None),
            )

        meta = _get_meta(feature, "gff")
        value = _meta_value(meta, "Is_circular", "is_circular")

        if _as_text(value).strip().lower() in {"true", "yes", "1"}:
            circular.add(str(feature.seqid))

    return circular


# ---------------------------------------------------------------------------
# Neighbourhood processing
# ---------------------------------------------------------------------------

def process_neighborhood(
    args: argparse.Namespace,
    gff_features: list,
    hit_features: list,
    seqs: Any,
    writer: OutputWriter,
) -> None:
    """Process GFF targets or search-hit targets."""

    genome_id = _get_genome_id(args)
    logger.info("Genome ID: %s", genome_id)

    processed_targets = set()
    target_count = 0
    neighbourhood_id = 0

    gff_features = sorted(
        gff_features,
        key=lambda ft: (
            str(ft.seqid),
            int(ft.loc.start),
            int(ft.loc.stop),
            str(ft.loc.strand),
            _feature_annotation_id(ft),
        ),
    )

    hit_features = sorted(
        hit_features,
        key=lambda ft: (
            str(ft.seqid),
            int(ft.loc.start),
            int(ft.loc.stop),
        ),
    )

    if args.input_type == "from_gff":
        targets = [
            feature for feature in gff_features
            if _is_target_feature(feature, args)
        ]
    else:
        targets = [
            feature for feature in hit_features
            if _is_target_feature(feature, args)
        ]

    logger.info("Query hits: %d", len(hit_features))
    logger.info("Annotated features available: %d", len(gff_features))
    logger.info("Targets selected: %d", len(targets))

    genome_id = _get_genome_id(args)
    logger.info("Genome ID: %s", genome_id)

    processed_targets = set()
    target_count = 0

    for target in targets:
        target_key = _feature_key(target)

        if target_key in processed_targets:
            logger.debug("Skipping duplicate target: %s", target_key)
            continue

        neighbourhood_id += 1

        mapped_annotation = None

        if args.input_type == "from_gff":
            center = target

        else:
            if args.input_type == "infernal":
                # Infernal-to-GFF mapping requires the same strand
                # and > overlap_threshold coverage of both features.
                mapped_annotation = _find_overlapping_annotation(
                    target,
                    gff_features,
                    overlap_threshold=args.overlap_threshold,
                    require_same_strand=True,
                )
            else:
                # Preserve existing behaviour for BLAST hits.
                mapped_annotation = _find_overlapping_annotation(
                    target,
                    gff_features,
                )

            center = mapped_annotation

            if center is None:
                center = target

                if args.input_type == "infernal":
                    logger.info(
                        "No same-strand GFF annotation with reciprocal "
                        "overlap > %.1f%% for Infernal hit "
                        "%s:%s-%s; using hit coordinates as the centre.",
                        args.overlap_threshold * 100,
                        target.seqid,
                        target.loc.start,
                        target.loc.stop,
                    )
                else:
                    logger.info(
                        "No overlapping GFF annotation for hit "
                        "%s:%s-%s; using hit coordinates as the centre.",
                        target.seqid,
                        target.loc.start,
                        target.loc.stop,
                    )
            else:
                logger.info(
                    "Mapped hit %s:%s-%s to GFF feature %s:%s-%s.",
                    target.seqid,
                    target.loc.start,
                    target.loc.stop,
                    center.seqid,
                    center.loc.start,
                    center.loc.stop,
                )

        neighbours = get_neighbor_features(
            gff_features,
            center,
            args.neighbours,
            args.overlap_threshold,
            include_all_overlaps=args.include_all_overlaps,
        )

        if not any(feature is center for feature in neighbours):
            neighbours.append(center)

        if args.input_type == "infernal":
            hit_start = int(target.loc.start)
            hit_stop = int(target.loc.stop)

            # If a GFF feature overlaps the Infernal hit, retain it
            # only if it passes the strict same-strand overlap filter.
            # Non-overlapping neighbours remain unaffected.
            neighbours = [
                feature
                for feature in neighbours
                if (
                    feature is target
                    or not (
                        feature.seqid == target.seqid
                        and int(feature.loc.start) < hit_stop
                        and hit_start < int(feature.loc.stop)
                    )
                    or _passes_same_strand_overlap(
                        target,
                        feature,
                        args.overlap_threshold,
                    )
                )
            ]

            # When a matching GFF annotation exists, output the
            # Infernal hit as a separate record as well.
            if mapped_annotation is not None:
                if not any(feature is target for feature in neighbours):
                    neighbours.append(target)

            # Keep deterministic genomic ordering.
            neighbours.sort(
                key=lambda ft: (
                    str(ft.seqid),
                    int(ft.loc.start),
                    int(ft.loc.stop),
                    str(ft.loc.strand),
                    _feature_annotation_id(ft),
                )
            )
        logger.info(
            "Selected %d features, including the centre, for %s:%s-%s",
            len(neighbours),
            center.seqid,
            center.loc.start,
            center.loc.stop,
        )

        for index, feature in enumerate(
            sorted(
                neighbours,
                key=lambda ft: (
                    str(ft.seqid),
                    int(ft.loc.start),
                    int(ft.loc.stop),
                ),
            ),
            start=1,
        ):
            logger.info(
                "  Neighbour %d: seqid=%s start=%s stop=%s "
                "strand=%s type=%s id=%s",
                index,
                feature.seqid,
                feature.loc.start,
                feature.loc.stop,
                feature.loc.strand,
                getattr(feature, "type", "unknown"),
                _feature_annotation_id(feature),
            )

        for neighbour in neighbours:
            entry = generate_entry(neighbour, args)

            # For Infernal hits, label the hit and its matched GFF
            # annotation with the requested query name (e.g. "sis").
            if (
                args.input_type == "infernal"
                and (
                    neighbour is target
                    or neighbour is mapped_annotation
                )
            ):
                entry = AnnotationEntry(
                    identifier=entry.identifier,
                    name=args.gene_of_interest,
                    start=entry.start,
                    stop=entry.stop,
                    orientation=entry.orientation,
                    biotype=entry.biotype,
                )

            entry = AnnotationEntry(
                identifier=entry.identifier,
                name=entry.name,
                start=entry.start,
                stop=entry.stop,
                orientation=entry.orientation,
                biotype=entry.biotype,
                member_id=f"{neighbour.seqid} ({neighbourhood_id})",
                genome_id=genome_id,
            )

            writer.write_feature(
                neighbour,
                entry,
                seqs,
                args.input_type,
                args.bio_type,
            )

        if args.promoter == "yes":
            sequence = _get_sequence(seqs, str(center.seqid))
            coords = get_promoter_coords(
                center,
                gff_features,
                args.promoter_len,
                args.promoter_mode,
                sequence_length=len(sequence),
            )

            if coords is not None:
                writer.write_promoter(
                    center,
                    generate_entry(center, args),
                    sequence,
                    coords[0],
                    coords[1],
                )

        processed_targets.add(target_key)
        target_count += 1

    if target_count == 0:
        logger.warning(
            "No target features found for %r. Check the identifier, "
            "GFF attributes, input type, and filtering thresholds.",
            args.gene_of_interest,
        )
    else:
        logger.info("Processed %d target feature(s).", target_count)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    logger.info("Executing script: %s", os.path.abspath(__file__))

    parser = setup_parser()
    args = parser.parse_args()
    validate_args(args, parser)

    logger.info("Processing target: %s", args.gene_of_interest)
    logger.info("Neighbour specification: %s", args.neighbours)
    logger.info("Score filter: %s", args.filter_by)

    logger.info("Reading GFF annotation...")
    raw_gff = read_fts(args.gff_file, fmt="gff")

    from collections import Counter

    logger.info(
        "Feature types in GFF: %s",
        Counter(getattr(feature, "type", None) for feature in raw_gff),
    )

    circular_seqids = get_circular_seqids(raw_gff)

    logger.info("Circular sequences detected in GFF: %d", len(circular_seqids))

    if circular_seqids:
        logger.info(
            "Circular sequence IDs: %s",
            ", ".join(sorted(circular_seqids)),
        )

    gff_features = validate_selected_gff_features(
        raw_gff,
        raw_gff.select(args.including_features),
        args,
    )

    logger.info("Reading genome FASTA...")
    seqs = read(args.fna_file)

    if args.input_type in {"blastn", "blastp", "infernal", "tblastn"}:
        logger.info("Reading and filtering hit file...")
        hit_features = process_hits(
            args.hit_file,
            args.input_type,
            args.gene_of_interest,
            args,
        )
        hit_features = remove_overlapping_features(
            hit_features,
            args.input_type,
        )
    else:
        hit_features = []

    validate_features_against_fasta(
        gff_features,
        hit_features,
        seqs,
        circular_seqids=circular_seqids,
    )

    logger.info("Input validation completed successfully.")

    with OutputWriter(
        args.output_path,
        circular_seqids=circular_seqids,
    ) as writer:
        process_neighborhood(
            args,
            gff_features,
            hit_features,
            seqs,
            writer,
        )
        writer.commit()

    logger.info("Processing completed successfully.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        logger.exception("An error occurred during processing.")
        raise
