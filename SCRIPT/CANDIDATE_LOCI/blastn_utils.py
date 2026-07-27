"""Strict BLASTN evidence parsing and interval matching for candidate loci.

The accepted table is headerless BLAST outfmt 6 with these first ten fields:
qseqid, sseqid, qlen, length, qstart, qend, sstart, send, nident, pident.
Every nonblank row must have the same width of at least twelve columns.

Production command contract::

    blastn ... -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore positive"
"""

import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Optional, Union

from attrs import define, field


@define(frozen=True, slots=True)
class ParametersBlastn:
    min_coverage: float = field(default=0.85)
    min_identity: float = field(default=0.95)

    def __attrs_post_init__(self) -> None:
        for name, value in (
            ("min_coverage", self.min_coverage),
            ("min_identity", self.min_identity),
        ):
            if (
                not isinstance(value, (int, float))
                or not math.isfinite(value)
                or not 0 <= value <= 1
            ):
                raise ValueError(f"{name} must be a finite fraction in [0, 1], got {value!r}")


@define(frozen=True, slots=True)
class BlastnHit:
    query_id: str
    chr_id: str
    query_length: int
    alignment_length: int
    nident: int
    query_start: int
    query_end: int
    start: int
    end: int
    strand: int
    coverage: float
    identity: float
    quality: float


@define(slots=True)
class BlastnStats:
    raw_rows: int = 0
    qualifying_rows: int = 0
    retained_hits: int = 0
    expanded_protein_loci: int = 0
    accepted_nucleotide_loci: int = 0
    discarded_overlap: int = 0
    discarded_filter: int = 0
    clipped_expansions: int = 0


@define(frozen=True, slots=True)
class PreparedBlastn:
    path: str
    parameters: ParametersBlastn
    hits_by_chr: Mapping[str, tuple[BlastnHit, ...]]
    stats: BlastnStats


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    """Return inclusive overlap length; adjacent intervals have zero overlap."""
    return max(0, min(end_a, end_b) - max(start_a, start_b) + 1)


def _parse_int(value: str, field_name: str, path: str, line_number: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{path}:{line_number}: invalid integer in {field_name}: {value!r}"
        ) from exc


def _parse_float(value: str, field_name: str, path: str, line_number: int) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{path}:{line_number}: invalid number in {field_name}: {value!r}"
        ) from exc
    if not math.isfinite(parsed):
        raise ValueError(
            f"{path}:{line_number}: {field_name} must be finite, got {value!r}"
        )
    return parsed


def _hit_sort_key(hit: BlastnHit) -> tuple:
    return (-hit.quality, hit.start, hit.end, hit.query_id, hit.chr_id, hit.strand)


def prepare_blastn(
    blastn_file: Union[str, Path],
    model_metadata: Mapping[str, object],
    params: Optional[ParametersBlastn] = None,
    chr: Optional[str] = None,
) -> PreparedBlastn:
    """Scan, validate, filter, deduplicate, and rank a BLASTN table once."""
    parameters = ParametersBlastn() if params is None else params
    parameters.__attrs_post_init__()
    path = str(blastn_file)
    stats = BlastnStats()
    expected_width: Optional[int] = None
    retained: dict[str, set[BlastnHit]] = defaultdict(set)

    try:
        file_handle = open(blastn_file, "r", newline="")
    except (OSError, TypeError) as exc:
        raise FileNotFoundError(f"BLASTN table is missing or unreadable: {path}") from exc

    with file_handle:
        for line_number, raw_line in enumerate(file_handle, start=1):
            if not raw_line.strip():
                continue
            stats.raw_rows += 1
            row = raw_line.rstrip("\r\n").split("\t")
            if expected_width is None:
                expected_width = len(row)
                if expected_width < 12:
                    raise ValueError(
                        f"{path}:{line_number}: BLASTN row has {expected_width} fields; "
                        "at least 12 are required"
                    )
            elif len(row) != expected_width:
                raise ValueError(
                    f"{path}:{line_number}: BLASTN row has {len(row)} fields; "
                    f"expected {expected_width}"
                )

            query_id, chr_id = row[0], row[1]
            qlen = _parse_int(row[2], "qlen", path, line_number)
            length = _parse_int(row[3], "length", path, line_number)
            qstart = _parse_int(row[4], "qstart", path, line_number)
            qend = _parse_int(row[5], "qend", path, line_number)
            sstart = _parse_int(row[6], "sstart", path, line_number)
            send = _parse_int(row[7], "send", path, line_number)
            nident = _parse_int(row[8], "nident", path, line_number)
            pident = _parse_float(row[9], "pident", path, line_number)

            if qlen <= 0 or length <= 0:
                raise ValueError(f"{path}:{line_number}: qlen and length must be positive")
            if min(qstart, qend, sstart, send) <= 0:
                raise ValueError(
                    f"{path}:{line_number}: BLASTN coordinates must be positive 1-based integers"
                )
            query_start, query_end = sorted((qstart, qend))
            if query_start < 1 or query_end > qlen:
                raise ValueError(
                    f"{path}:{line_number}: query coordinates {query_start}-{query_end} "
                    f"fall outside 1..{qlen}"
                )
            if not 0 <= pident <= 100:
                raise ValueError(
                    f"{path}:{line_number}: pident must be in [0, 100], got {pident}"
                )

            coverage = length / qlen
            identity = pident / 100.0
            if coverage < parameters.min_coverage or identity < parameters.min_identity:
                continue
            stats.qualifying_rows += 1
            if chr is not None and chr_id != chr:
                continue
            if query_id not in model_metadata:
                raise ValueError(
                    f"{path}:{line_number}: qualifying BLASTN query/model ID "
                    f"{query_id!r} is absent from the reference GFF"
                )
            start, end = sorted((sstart, send))
            retained[chr_id].add(
                BlastnHit(
                    query_id=query_id,
                    chr_id=chr_id,
                    query_length=qlen,
                    alignment_length=length,
                    nident=nident,
                    query_start=query_start,
                    query_end=query_end,
                    start=start,
                    end=end,
                    strand=1 if sstart <= send else -1,
                    coverage=coverage,
                    identity=identity,
                    quality=coverage * identity,
                )
            )

    hits_by_chr = {
        chr_id: tuple(sorted(hits, key=_hit_sort_key))
        for chr_id, hits in retained.items()
    }
    stats.retained_hits = sum(len(hits) for hits in hits_by_chr.values())
    return PreparedBlastn(path, parameters, hits_by_chr, stats)


def select_same_model_candidate(hit: BlastnHit, candidates: Iterable[object]):
    """Choose the overlapping same-model protein candidate deterministically."""
    matches = []
    for candidate in candidates:
        if candidate.prot_id != hit.query_id:
            continue
        overlap = interval_overlap(
            hit.start, hit.end, candidate.chr_bounds.start, candidate.chr_bounds.end
        )
        if overlap:
            matches.append((candidate, overlap))
    if not matches:
        return None
    return min(
        matches,
        key=lambda item: (
            -item[1],
            -item[0].score,
            item[0].chr_bounds.start,
            item[0].chr_bounds.end,
        ),
    )[0]


def blastn_desired_expansion(hit: BlastnHit, candidates: Iterable[object]):
    """Return candidate and genomic left/right demand for a same-model hit."""
    candidate = select_same_model_candidate(hit, candidates)
    if candidate is None:
        return None
    return (
        candidate,
        max(0, candidate.chr_bounds.start - hit.start),
        max(0, hit.end - candidate.chr_bounds.end),
    )


def blastn_new_loci_discovery(
    hit: BlastnHit,
    protein_candidates: Iterable[object],
    accepted_blastn_candidates: Iterable[object],
) -> bool:
    """Whether a ranked hit is unblocked by accepted protein/BLASTN intervals."""
    for candidate in (*tuple(protein_candidates), *tuple(accepted_blastn_candidates)):
        if interval_overlap(
            hit.start, hit.end, candidate.chr_bounds.start, candidate.chr_bounds.end
        ):
            return False
    return True


__all__ = [
    "BlastnHit",
    "BlastnStats",
    "ParametersBlastn",
    "PreparedBlastn",
    "blastn_desired_expansion",
    "blastn_new_loci_discovery",
    "interval_overlap",
    "prepare_blastn",
    "select_same_model_candidate",
]
