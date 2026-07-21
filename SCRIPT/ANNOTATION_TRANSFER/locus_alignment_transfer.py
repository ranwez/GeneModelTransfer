#!/usr/bin/env python3
"""Transfer one transcript's CDS boundaries between nucleotide loci.

The public functions in this module deliberately separate parsing, eligibility,
alignment traversal, projection, and process orchestration.  This keeps the
coordinate-sensitive code testable without requiring BLAST+ or MAFFT.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence


BLAST_OUTFMT_FIELDS = (
    "qseqid sseqid qlen length qstart qend sstart send "
    "nident pident bitscore evalue"
)
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


class TransferError(Exception):
    """Base class for expected transfer failures."""


class InputValidationError(TransferError):
    """An input FASTA, GFF, or command-line value is invalid."""


class AlignmentError(TransferError):
    """BLAST, MAFFT, or an alignment result is invalid."""


class ProjectionError(TransferError):
    """CDS boundaries cannot be projected without violating the contract."""


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    sequence: str


@dataclass(frozen=True)
class BlastHsp:
    query_id: str
    subject_id: str
    query_length: int
    alignment_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    identical: int
    percent_identity: float
    bit_score: float
    evalue: float


@dataclass(frozen=True)
class CdsFeature:
    start: int
    end: int
    strand: str
    phase: str


@dataclass(frozen=True)
class ModelAnnotation:
    gene_id: str
    transcript_id: str
    sequence_id: str
    gene_start: int
    gene_end: int
    strand: str
    cds: tuple[CdsFeature, ...]


@dataclass(frozen=True)
class EligibilityDecision:
    eligible: bool
    selected_hsp: BlastHsp | None
    coverage: float | None
    identity: float | None
    orientation: str | None
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class TransferResult:
    status: str
    exit_code: int
    reason: str
    output_gff: Path
    diagnostics_json: Path
    projected_cds: tuple[CdsFeature, ...] = ()


@dataclass(frozen=True)
class _GffRow:
    sequence_id: str
    source: str
    feature_type: str
    start: int
    end: int
    strand: str
    phase: str
    attributes: Mapping[str, str]
    line_number: int


def read_single_fasta(path: str | os.PathLike[str]) -> FastaRecord:
    """Read exactly one non-empty IUPAC DNA FASTA record."""

    fasta_path = Path(path)
    try:
        lines = fasta_path.read_text().splitlines()
    except OSError as exc:
        raise InputValidationError(f"Cannot read FASTA {fasta_path}: {exc}") from exc

    headers: list[tuple[int, str]] = []
    sequence_parts: list[str] = []
    seen_header = False
    for line_number, raw_line in enumerate(lines, 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            identifier = line[1:].split(None, 1)[0] if line[1:].strip() else ""
            headers.append((line_number, identifier))
            seen_header = True
            continue
        if not seen_header:
            raise InputValidationError(f"FASTA {fasta_path} has sequence before its header")
        sequence_parts.append("".join(line.split()))

    if len(headers) != 1:
        raise InputValidationError(f"FASTA {fasta_path} must contain exactly one record")
    if not headers[0][1]:
        raise InputValidationError(f"FASTA {fasta_path} has an empty record identifier")
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise InputValidationError(f"FASTA {fasta_path} record must be non-empty")
    invalid = sorted(set(sequence) - DNA_ALPHABET)
    if invalid:
        raise InputValidationError(
            f"FASTA {fasta_path} contains non-IUPAC DNA symbols: {''.join(invalid)}"
        )
    return FastaRecord(headers[0][1], sequence)


def _parse_attributes(text: str, path: Path, line_number: int) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise InputValidationError(
                f"Malformed GFF attributes in {path} line {line_number}: {item!r}"
            )
        key, value = item.split("=", 1)
        if not key or not value or key in attributes:
            raise InputValidationError(
                f"Malformed GFF attributes in {path} line {line_number}: {item!r}"
            )
        attributes[key] = value
    return attributes


def _read_gff_rows(path: str | os.PathLike[str]) -> list[_GffRow]:
    gff_path = Path(path)
    try:
        lines = gff_path.read_text().splitlines()
    except OSError as exc:
        raise InputValidationError(f"Cannot read GFF {gff_path}: {exc}") from exc

    rows: list[_GffRow] = []
    for line_number, line in enumerate(lines, 1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise InputValidationError(
                f"Malformed GFF row in {gff_path} line {line_number}: expected 9 columns"
            )
        try:
            start = int(fields[3])
            end = int(fields[4])
        except ValueError as exc:
            raise InputValidationError(
                f"Malformed GFF coordinates in {gff_path} line {line_number}"
            ) from exc
        if start < 1 or end < start:
            raise InputValidationError(
                f"Invalid GFF interval in {gff_path} line {line_number}: {start}-{end}"
            )
        rows.append(
            _GffRow(
                sequence_id=fields[0],
                source=fields[1],
                feature_type=fields[2],
                start=start,
                end=end,
                strand=fields[6],
                phase=fields[7],
                attributes=_parse_attributes(fields[8], gff_path, line_number),
                line_number=line_number,
            )
        )
    if not rows:
        raise InputValidationError(f"GFF {gff_path} contains no feature rows")
    return rows


def _select_and_validate_model(
    rows: Sequence[_GffRow],
    sequence_length: int | None,
    model_id: str | None,
) -> ModelAnnotation:
    if sequence_length is not None and sequence_length < 1:
        raise InputValidationError("Sequence length must be positive")

    gene_rows = [row for row in rows if row.feature_type == "gene"]
    if any("ID" not in row.attributes for row in gene_rows):
        raise InputValidationError("Every gene row must have an ID")
    genes = gene_rows
    if model_id is None:
        if len(genes) != 1:
            raise InputValidationError("Model GFF must contain exactly one gene")
        gene = genes[0]
        model_id = gene.attributes["ID"]
    else:
        selected = [row for row in genes if row.attributes["ID"] == model_id]
        if len(selected) != 1:
            raise InputValidationError(
                f"Model GFF must contain exactly one gene with ID={model_id}"
            )
        gene = selected[0]

    transcripts = [
        row
        for row in rows
        if row.feature_type == "mRNA" and row.attributes.get("Parent") == model_id
    ]
    if len(transcripts) != 1:
        raise InputValidationError(
            f"Gene {model_id} must have exactly one directly parented mRNA"
        )
    transcript = transcripts[0]
    transcript_id = transcript.attributes.get("ID")
    if not transcript_id:
        raise InputValidationError(f"mRNA for gene {model_id} has no ID")

    cds_rows = [
        row
        for row in rows
        if row.feature_type == "CDS" and row.attributes.get("Parent") == transcript_id
    ]
    if not cds_rows:
        raise InputValidationError(f"mRNA {transcript_id} must have at least one CDS")
    ambiguous_cds = [
        row
        for row in rows
        if row.feature_type == "CDS"
        and transcript_id in row.attributes.get("Parent", "").split(",")
        and row.attributes.get("Parent") != transcript_id
    ]
    if ambiguous_cds:
        raise InputValidationError(f"Every CDS for {transcript_id} must belong only to that mRNA")
    if len(genes) == 1:
        all_transcripts = [row for row in rows if row.feature_type == "mRNA"]
        all_cds = [row for row in rows if row.feature_type == "CDS"]
        if all_transcripts != transcripts or all_cds != cds_rows:
            raise InputValidationError(
                "Single-model GFF contains an orphan or unrelated mRNA/CDS feature"
            )

    relevant = [gene, transcript, *cds_rows]
    if gene.strand not in {"+", "-"}:
        raise InputValidationError(f"Gene {model_id} must have strand + or -")
    for row in relevant:
        if row.sequence_id != gene.sequence_id:
            raise InputValidationError(f"Features for gene {model_id} use mixed sequence IDs")
        if row.strand != gene.strand:
            raise InputValidationError(f"Features for gene {model_id} use mixed strands")
        if row.start < gene.start or row.end > gene.end:
            raise InputValidationError(f"Feature outside gene {model_id} interval")
    if transcript.start < gene.start or transcript.end > gene.end:
        raise InputValidationError(f"mRNA {transcript_id} lies outside gene {model_id}")
    for row in cds_rows:
        if row.start < transcript.start or row.end > transcript.end:
            raise InputValidationError(f"CDS lies outside mRNA {transcript_id}")
        if row.phase not in {"0", "1", "2", "."}:
            raise InputValidationError(f"CDS in line {row.line_number} has invalid phase")

    genomic_cds = sorted(cds_rows, key=lambda row: (row.start, row.end))
    for previous, current in zip(genomic_cds, genomic_cds[1:]):
        if current.start <= previous.end:
            raise InputValidationError(f"CDS intervals for {transcript_id} overlap")
    if sequence_length is not None and gene.end > sequence_length:
        raise InputValidationError(f"Gene {model_id} lies outside the sequence")

    return ModelAnnotation(
        gene_id=model_id,
        transcript_id=transcript_id,
        sequence_id=gene.sequence_id,
        gene_start=gene.start,
        gene_end=gene.end,
        strand=gene.strand,
        cds=tuple(
            CdsFeature(row.start, row.end, row.strand, row.phase) for row in genomic_cds
        ),
    )


def parse_model_gff(
    path: str | os.PathLike[str], sequence_length: int | None = None
) -> ModelAnnotation:
    """Parse a GFF that contains exactly one single-transcript model."""

    return _select_and_validate_model(_read_gff_rows(path), sequence_length, None)


def normalize_model_gff(
    path: str | os.PathLike[str], model_id: str, locus_length: int
) -> ModelAnnotation:
    """Select a model and convert genomic coordinates to transcription-locus coordinates."""

    source = _select_and_validate_model(_read_gff_rows(path), None, model_id)
    gene_length = source.gene_end - source.gene_start + 1
    if gene_length != locus_length:
        raise InputValidationError(
            f"Gene {model_id} length {gene_length} does not match locus length {locus_length}"
        )

    normalized: list[CdsFeature] = []
    for cds in source.cds:
        if source.strand == "+":
            start = cds.start - source.gene_start + 1
            end = cds.end - source.gene_start + 1
        else:
            start = source.gene_end - cds.end + 1
            end = source.gene_end - cds.start + 1
        normalized.append(CdsFeature(start, end, "+", cds.phase))
    normalized.sort(key=lambda item: (item.start, item.end))
    return ModelAnnotation(
        gene_id=source.gene_id,
        transcript_id=source.transcript_id,
        sequence_id=source.gene_id,
        gene_start=1,
        gene_end=locus_length,
        strand="+",
        cds=tuple(normalized),
    )


def parse_blast_tabular(text: str) -> list[BlastHsp]:
    """Parse the fixed twelve-column V3 nucleotide BLAST format."""

    hsps: list[BlastHsp] = []
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 12:
            raise AlignmentError(
                f"Malformed BLAST output line {line_number}: expected 12 columns"
            )
        try:
            hsp = BlastHsp(
                query_id=fields[0],
                subject_id=fields[1],
                query_length=int(fields[2]),
                alignment_length=int(fields[3]),
                query_start=int(fields[4]),
                query_end=int(fields[5]),
                subject_start=int(fields[6]),
                subject_end=int(fields[7]),
                identical=int(fields[8]),
                percent_identity=float(fields[9]),
                bit_score=float(fields[10]),
                evalue=float(fields[11]),
            )
        except ValueError as exc:
            raise AlignmentError(f"Malformed numeric BLAST value on line {line_number}") from exc
        if (
            hsp.query_length < 1
            or hsp.alignment_length < 1
            or hsp.query_start < 1
            or hsp.query_end < 1
            or hsp.subject_start < 1
            or hsp.subject_end < 1
            or hsp.identical < 0
            or not (0.0 <= hsp.percent_identity <= 100.0)
            or hsp.bit_score < 0.0
            or hsp.evalue < 0.0
        ):
            raise AlignmentError(f"Invalid BLAST values on line {line_number}")
        hsps.append(hsp)
    return hsps


def _hsp_rank(hsp: BlastHsp) -> tuple[object, ...]:
    return (
        -hsp.bit_score,
        hsp.evalue,
        -hsp.alignment_length,
        -hsp.percent_identity,
        hsp.query_start,
        hsp.query_end,
        hsp.subject_start,
        hsp.subject_end,
        hsp.query_id,
        hsp.subject_id,
    )


def assess_eligibility(
    hsps: Iterable[BlastHsp],
    minimum_coverage: float = 0.85,
    minimum_identity: float = 0.95,
) -> EligibilityDecision:
    """Rank individual HSPs deterministically and gate only the best one."""

    if not 0.0 <= minimum_coverage <= 1.0:
        raise InputValidationError("Minimum coverage must be between 0 and 1")
    if not 0.0 <= minimum_identity <= 1.0:
        raise InputValidationError("Minimum identity must be between 0 and 1")
    candidates = list(hsps)
    if not candidates:
        return EligibilityDecision(False, None, None, None, None, ("no_significant_hit",))
    selected = min(candidates, key=_hsp_rank)
    if selected.query_length < 1:
        raise AlignmentError("Selected BLAST HSP has an invalid query length")
    coverage = selected.alignment_length / selected.query_length
    identity = selected.percent_identity / 100.0
    orientation = "forward" if selected.subject_end >= selected.subject_start else "reverse"
    reasons: list[str] = []
    if coverage < minimum_coverage:
        reasons.append("coverage_below_threshold")
    if identity < minimum_identity:
        reasons.append("identity_below_threshold")
    return EligibilityDecision(
        not reasons,
        selected,
        coverage,
        identity,
        orientation,
        tuple(reasons),
    )


def projection_candidate_positions(cds_features: Iterable[CdsFeature]) -> set[int]:
    """Return only CDS endpoints and bounded inward 3/6/9-nt alternatives."""

    positions: set[int] = set()
    for cds in cds_features:
        if cds.start < 1 or cds.end < cds.start:
            raise ProjectionError(f"Invalid CDS interval {cds.start}-{cds.end}")
        positions.update((cds.start, cds.end))
        for offset in (3, 6, 9):
            if cds.start + offset <= cds.end:
                positions.add(cds.start + offset)
            if cds.end - offset >= cds.start:
                positions.add(cds.end - offset)
    return positions


def map_requested_model_positions(
    model_alignment: str,
    target_alignment: str,
    requested_positions: Iterable[int],
    expected_model: str | None = None,
    expected_target: str | None = None,
) -> dict[int, int | None]:
    """Traverse one alignment and retain mappings only for requested model positions."""

    if len(model_alignment) != len(target_alignment):
        raise AlignmentError("Aligned model and target must have the same length")
    requested = set(requested_positions)
    if any(not isinstance(position, int) for position in requested):
        raise ProjectionError("Every requested model position must be an integer")

    ungapped_model_length = sum(character != "-" for character in model_alignment)
    if any(position < 1 or position > ungapped_model_length for position in requested):
        raise ProjectionError("Every requested model position must lie inside the model")

    mapping: dict[int, int | None] = {}
    model_position = 0
    target_position = 0
    model_ungapped: list[str] | None = [] if expected_model is not None else None
    target_ungapped: list[str] | None = [] if expected_target is not None else None
    for model_base, target_base in zip(model_alignment, target_alignment):
        if model_base != "-":
            model_position += 1
            if model_ungapped is not None:
                model_ungapped.append(model_base)
        if target_base != "-":
            target_position += 1
            if target_ungapped is not None:
                target_ungapped.append(target_base)
        if model_base != "-" and model_position in requested:
            mapping[model_position] = target_position if target_base != "-" else None

    if model_ungapped is not None and "".join(model_ungapped).upper() != expected_model.upper():
        raise AlignmentError("Ungapped aligned model does not reproduce the model input")
    if target_ungapped is not None and "".join(target_ungapped).upper() != expected_target.upper():
        raise AlignmentError("Ungapped aligned target does not reproduce the target input")
    if mapping.keys() != requested:
        missing = sorted(requested - mapping.keys())
        raise ProjectionError(f"Missing requested model positions: {missing}")
    return mapping


def reverse_coordinate(position: int, sequence_length: int) -> int:
    if sequence_length < 1 or position < 1 or position > sequence_length:
        raise ProjectionError("Coordinate to reverse lies outside the target")
    return sequence_length - position + 1


def _mapped_boundary(
    cds: CdsFeature,
    boundary: str,
    mapping: Mapping[int, int | None],
) -> tuple[int, int]:
    endpoint = cds.start if boundary == "start" else cds.end
    candidates = [endpoint]
    for offset in (3, 6, 9):
        candidate = endpoint + offset if boundary == "start" else endpoint - offset
        if cds.start <= candidate <= cds.end:
            candidates.append(candidate)
    for candidate in candidates:
        if candidate not in mapping:
            raise ProjectionError(f"Missing requested {boundary} boundary candidate {candidate}")
        target_position = mapping[candidate]
        if target_position is not None:
            return target_position, candidate
    raise ProjectionError(
        f"CDS {boundary} boundary {endpoint} cannot be mapped within three codons"
    )


def _project_cds_with_snaps(
    cds_features: Iterable[CdsFeature],
    mapping: Mapping[int, int | None],
    target_length: int,
    orientation: str,
) -> tuple[list[CdsFeature], list[dict[str, int | str]]]:
    if orientation not in {"forward", "reverse"}:
        raise ProjectionError("Projection orientation must be forward or reverse")
    if target_length < 1:
        raise ProjectionError("Target length must be positive")

    source_features = sorted(cds_features, key=lambda item: (item.start, item.end))
    if not source_features:
        raise ProjectionError("At least one CDS is required")
    projected_in_model_order: list[CdsFeature] = []
    snaps: list[dict[str, int | str]] = []
    prepared_intervals: list[tuple[int, int]] = []
    for cds_index, cds in enumerate(source_features, 1):
        mapped_start, used_start = _mapped_boundary(cds, "start", mapping)
        mapped_end, used_end = _mapped_boundary(cds, "end", mapping)
        if used_start != cds.start:
            snaps.append(
                {
                    "cds": cds_index,
                    "boundary": "start",
                    "requested_model_position": cds.start,
                    "used_model_position": used_start,
                    "offset": used_start - cds.start,
                }
            )
        if used_end != cds.end:
            snaps.append(
                {
                    "cds": cds_index,
                    "boundary": "end",
                    "requested_model_position": cds.end,
                    "used_model_position": used_end,
                    "offset": used_end - cds.end,
                }
            )
        prepared_start, prepared_end = sorted((mapped_start, mapped_end))
        if prepared_start < 1 or prepared_end > target_length or prepared_start > prepared_end:
            raise ProjectionError("Projected CDS lies outside the target")
        prepared_intervals.append((prepared_start, prepared_end))
        if orientation == "forward":
            start, end, strand = prepared_start, prepared_end, "+"
        else:
            start = reverse_coordinate(prepared_end, target_length)
            end = reverse_coordinate(prepared_start, target_length)
            strand = "-"
        projected_in_model_order.append(CdsFeature(start, end, strand, cds.phase))

    for previous, current in zip(prepared_intervals, prepared_intervals[1:]):
        if current[0] <= previous[1]:
            raise ProjectionError("Projected CDS intervals overlap or are inconsistently ordered")
    unique_intervals = {(item.start, item.end) for item in projected_in_model_order}
    if len(unique_intervals) != len(projected_in_model_order):
        raise ProjectionError("Projected CDS intervals are duplicated")
    return sorted(projected_in_model_order, key=lambda item: (item.start, item.end)), snaps


def project_cds_features(
    cds_features: Iterable[CdsFeature],
    mapping: Mapping[int, int | None],
    target_length: int,
    orientation: str,
) -> list[CdsFeature]:
    projected, _ = _project_cds_with_snaps(
        cds_features, mapping, target_length, orientation
    )
    return projected


def render_draft_gff(
    target_id: str,
    cds_features: Iterable[CdsFeature],
    source: str = "locusAlignment",
) -> str:
    """Render deterministic locus-relative gene + CDS GFF records."""

    features = sorted(cds_features, key=lambda item: (item.start, item.end))
    if not features:
        raise ProjectionError("Cannot render a draft without CDS features")
    strands = {feature.strand for feature in features}
    if len(strands) != 1 or not strands <= {"+", "-"}:
        raise ProjectionError("Draft CDS features must have one valid strand")
    for previous, current in zip(features, features[1:]):
        if current.start <= previous.end:
            raise ProjectionError("Draft CDS features must not overlap")
    identifier = f"{target_id}"
    strand = features[0].strand
    rows = [
        f"{target_id}\t{source}\tgene\t{features[0].start}\t{features[-1].end}"
        f"\t.\t{strand}\t.\tID={identifier}"
    ]
    for index, feature in enumerate(features, 1):
        rows.append(
            f"{target_id}\t{source}\tCDS\t{feature.start}\t{feature.end}\t.\t"
            f"{feature.strand}\t{feature.phase}\tID={identifier}_cds{index};Parent={identifier}"
        )
    return "\n".join(rows) + "\n"


def render_normalized_model_gff(annotation: ModelAnnotation) -> str:
    """Render a complete one-transcript model in locus transcription orientation."""

    if annotation.gene_start != 1 or annotation.strand != "+":
        raise InputValidationError("Reference-locus annotation must be normalized to strand +")
    rows = [
        f"{annotation.gene_id}\tlocusAlignmentReference\tgene\t1\t{annotation.gene_end}"
        f"\t.\t+\t.\tID={annotation.gene_id}",
        f"{annotation.gene_id}\tlocusAlignmentReference\tmRNA\t1\t{annotation.gene_end}"
        f"\t.\t+\t.\tID={annotation.transcript_id};Parent={annotation.gene_id}",
    ]
    for index, cds in enumerate(annotation.cds, 1):
        rows.append(
            f"{annotation.gene_id}\tlocusAlignmentReference\tCDS\t{cds.start}\t{cds.end}"
            f"\t.\t+\t{cds.phase}\tID={annotation.gene_id}_cds{index};"
            f"Parent={annotation.transcript_id}"
        )
    return "\n".join(rows) + "\n"


def reverse_complement(sequence: str) -> str:
    try:
        return sequence.translate(COMPLEMENT)[::-1]
    except Exception as exc:  # pragma: no cover - input validation normally prevents this
        raise InputValidationError("Cannot reverse-complement target sequence") from exc


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    _atomic_write(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _tool_version(executable: str) -> str:
    variants = ([executable, "-version"], [executable, "--version"])
    for command in variants:
        try:
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
        except OSError:
            return "unavailable"
        text = (completed.stdout or completed.stderr).strip()
        if completed.returncode == 0 and text:
            return text.splitlines()[0]
    return "unknown"


def _default_blast_runner(
    model_fasta: Path, target_fasta: Path, executable: str
) -> Callable[[FastaRecord, FastaRecord], list[BlastHsp]]:
    command = [
        executable,
        "-query",
        str(model_fasta),
        "-subject",
        str(target_fasta),
        "-outfmt",
        f"6 {BLAST_OUTFMT_FIELDS}",
    ]

    def run(_model: FastaRecord, _target: FastaRecord) -> list[BlastHsp]:
        try:
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
        except OSError as exc:
            raise AlignmentError(f"BLAST could not be executed: {exc}") from exc
        if completed.returncode != 0:
            message = completed.stderr.strip() or "no error message"
            raise AlignmentError(f"BLAST failed with exit {completed.returncode}: {message}")
        return parse_blast_tabular(completed.stdout)

    return run


def _parse_two_record_alignment(text: str) -> tuple[str, str]:
    records: list[str] = []
    sequence: list[str] = []
    seen_header = False
    for line in text.splitlines():
        if line.startswith(">"):
            if seen_header:
                records.append("".join(sequence))
            seen_header = True
            sequence = []
        elif line.strip():
            if not seen_header:
                raise AlignmentError("MAFFT output has sequence before a FASTA header")
            sequence.append(line.strip())
    if seen_header:
        records.append("".join(sequence))
    if len(records) != 2 or any(not record for record in records):
        raise AlignmentError("MAFFT output must contain exactly two non-empty records")
    return records[0], records[1]


def _default_mafft_runner(
    executable: str,
) -> Callable[[str, str], tuple[str, str]]:
    def run(model_sequence: str, target_sequence: str) -> tuple[str, str]:
        with tempfile.TemporaryDirectory(prefix="locusAlignment-mafft-") as directory:
            input_path = Path(directory) / "input.fasta"
            input_path.write_text(
                f">model\n{model_sequence}\n>target\n{target_sequence}\n"
            )
            command = [executable, "--auto", str(input_path)]
            try:
                completed = subprocess.run(command, capture_output=True, text=True, check=False)
            except OSError as exc:
                raise AlignmentError(f"MAFFT could not be executed: {exc}") from exc
            if completed.returncode != 0:
                message = completed.stderr.strip() or "no error message"
                raise AlignmentError(f"MAFFT failed with exit {completed.returncode}: {message}")
            return _parse_two_record_alignment(completed.stdout)

    return run


def _hsp_payload(hsp: BlastHsp | None) -> dict[str, object] | None:
    return asdict(hsp) if hsp is not None else None


def run_transfer(
    model_fasta: str | os.PathLike[str],
    model_gff: str | os.PathLike[str],
    target_fasta: str | os.PathLike[str],
    output_gff: str | os.PathLike[str],
    diagnostics_json: str | os.PathLike[str],
    minimum_coverage: float = 0.85,
    minimum_identity: float = 0.95,
    blast_executable: str = "blastn",
    mafft_executable: str = "mafft",
    blast_runner: Callable[[FastaRecord, FastaRecord], Iterable[BlastHsp]] | None = None,
    mafft_runner: Callable[[str, str], tuple[str, str]] | None = None,
) -> TransferResult:
    """Validate, gate, align, project, and atomically write one transfer."""

    started = time.perf_counter()
    model_fasta_path = Path(model_fasta)
    target_fasta_path = Path(target_fasta)
    output_path = Path(output_gff)
    diagnostics_path = Path(diagnostics_json)
    diagnostics: dict[str, object] = {
        "status": "error",
        "reason": "not_started",
        "thresholds": {
            "minimum_coverage": minimum_coverage,
            "minimum_identity": minimum_identity,
        },
        "commands": {
            "blast": [
                blast_executable,
                "-query",
                str(model_fasta_path),
                "-subject",
                str(target_fasta_path),
                "-outfmt",
                f"6 {BLAST_OUTFMT_FIELDS}",
            ],
            "mafft": [mafft_executable, "--auto", "<two-sequence-fasta>"],
        },
        "versions": {
            "blast": "not_invoked" if blast_runner is None else "injected",
            "mafft": "not_invoked" if mafft_runner is None else "injected",
        },
        "tool_exit_states": {"blast": "not_invoked", "mafft": "not_invoked"},
        "boundary_snaps": [],
        "timing_seconds": {},
    }

    try:
        if not 0.0 <= minimum_coverage <= 1.0 or not 0.0 <= minimum_identity <= 1.0:
            raise InputValidationError("Coverage and identity thresholds must be between 0 and 1")
        model = read_single_fasta(model_fasta_path)
        target = read_single_fasta(target_fasta_path)
        annotation = normalize_model_gff(model_gff, model.identifier, len(model.sequence))
        diagnostics.update(
            {
                "model_id": model.identifier,
                "target_id": target.identifier,
                "model_length": len(model.sequence),
                "target_length": len(target.sequence),
            }
        )

        effective_blast_runner = blast_runner or _default_blast_runner(
            model_fasta_path, target_fasta_path, blast_executable
        )
        if blast_runner is None:
            diagnostics["versions"]["blast"] = _tool_version(blast_executable)  # type: ignore[index]
        blast_started = time.perf_counter()
        diagnostics["tool_exit_states"]["blast"] = "running"  # type: ignore[index]
        hsps = list(effective_blast_runner(model, target))
        diagnostics["tool_exit_states"]["blast"] = "success"  # type: ignore[index]
        diagnostics["timing_seconds"]["blast"] = time.perf_counter() - blast_started  # type: ignore[index]
        for hsp in hsps:
            if hsp.query_id != model.identifier or hsp.subject_id != target.identifier:
                raise AlignmentError("BLAST HSP identifiers do not match the input FASTA records")
            if hsp.query_length != len(model.sequence):
                raise AlignmentError("BLAST HSP query length does not match the model FASTA")

        decision = assess_eligibility(hsps, minimum_coverage, minimum_identity)
        diagnostics.update(
            {
                "selected_hsp": _hsp_payload(decision.selected_hsp),
                "coverage": decision.coverage,
                "identity": decision.identity,
                "relative_orientation": decision.orientation,
            }
        )
        if not decision.eligible:
            reason = decision.reasons[0]
            diagnostics.update(
                {
                    "status": "ineligible",
                    "reason": reason,
                    "reasons": list(decision.reasons),
                    "timing_seconds": {**diagnostics["timing_seconds"], "total": time.perf_counter() - started},  # type: ignore[dict-item]
                }
            )
            _atomic_write(output_path, "")
            _atomic_json(diagnostics_path, diagnostics)
            return TransferResult("ineligible", 0, reason, output_path, diagnostics_path)

        assert decision.orientation is not None
        prepared_target = (
            target.sequence
            if decision.orientation == "forward"
            else reverse_complement(target.sequence)
        )
        effective_mafft_runner = mafft_runner or _default_mafft_runner(mafft_executable)
        if mafft_runner is None:
            diagnostics["versions"]["mafft"] = _tool_version(mafft_executable)  # type: ignore[index]
        mafft_started = time.perf_counter()
        diagnostics["tool_exit_states"]["mafft"] = "running"  # type: ignore[index]
        model_alignment, target_alignment = effective_mafft_runner(
            model.sequence, prepared_target
        )
        diagnostics["tool_exit_states"]["mafft"] = "success"  # type: ignore[index]
        diagnostics["timing_seconds"]["mafft"] = time.perf_counter() - mafft_started  # type: ignore[index]

        try:
          requested = projection_candidate_positions(annotation.cds)
          mapping = map_requested_model_positions(
              model_alignment,
              target_alignment,
              requested,
              expected_model=model.sequence,
              expected_target=prepared_target,
          )
          projected, snaps = _project_cds_with_snaps(
              annotation.cds,
              mapping,
              len(target.sequence),
              decision.orientation,
          )
          diagnostics["boundary_snaps"] = snaps
          rendered = render_draft_gff(target.identifier, projected)
          diagnostics.update(
              {
                  "status": "success",
                  "reason": "projected",
                  "timing_seconds": {**diagnostics["timing_seconds"], "total": time.perf_counter() - started},  # type: ignore[dict-item]
              }
          )
        except TransferError as exc:
          print(f"WARNING: locusAlignment: {exc}", file=sys.stderr)
          rendered=""
          projected=()
          diagnostics.update(
              {
                "status": "failure",
                "reason": "projection failed"
              }
          )

        _atomic_write(output_path, rendered)
        _atomic_json(diagnostics_path, diagnostics)
        return TransferResult(
            "success",
            0,
            "projected",
            output_path,
            diagnostics_path,
            tuple(projected),
        )
    except Exception as exc:
        if output_path.exists():
            try:
                output_path.unlink()
            except OSError:
                _atomic_write(output_path, "")
        tool_states = diagnostics["tool_exit_states"]
        if isinstance(tool_states, dict):
            for tool, state in list(tool_states.items()):
                if state == "running":
                    tool_states[tool] = "error"
        diagnostics.update(
            {
                "status": "error",
                "reason": exc.__class__.__name__,
                "error": str(exc),
                "timing_seconds": {**diagnostics["timing_seconds"], "total": time.perf_counter() - started},  # type: ignore[dict-item]
            }
        )
        try:
            _atomic_json(diagnostics_path, diagnostics)
        except OSError:
            pass
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Transfer one model transcript's CDS boundaries between nucleotide loci"
    )
    parser.add_argument("--model-fasta", required=True)
    parser.add_argument("--model-gff", required=True)
    parser.add_argument("--target-fasta", required=True)
    parser.add_argument("--output-gff", required=True)
    parser.add_argument("--diagnostics-json", required=True)
    parser.add_argument("--minimum-coverage", type=float, default=0.85)
    parser.add_argument("--minimum-identity", type=float, default=0.95)
    parser.add_argument("--blast-executable", default="blastn")
    parser.add_argument("--mafft-executable", default="mafft")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _argument_parser().parse_args(argv)
    try:
        _result = run_transfer(
            model_fasta=arguments.model_fasta,
            model_gff=arguments.model_gff,
            target_fasta=arguments.target_fasta,
            output_gff=arguments.output_gff,
            diagnostics_json=arguments.diagnostics_json,
            minimum_coverage=arguments.minimum_coverage,
            minimum_identity=arguments.minimum_identity,
            blast_executable=arguments.blast_executable,
            mafft_executable=arguments.mafft_executable,
        )
    except TransferError:
        return 0
    except Exception as exc:  # unexpected failures are also unscoreable
        print(f"locusAlignment: unexpected error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
