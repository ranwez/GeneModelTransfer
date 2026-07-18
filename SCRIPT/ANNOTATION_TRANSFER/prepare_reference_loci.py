#!/usr/bin/env python3
"""Create normalized per-model GFF and checksum provenance for REF_LOCI."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from locus_alignment_transfer import (
    InputValidationError,
    _read_gff_rows,
    normalize_model_gff,
    read_single_fasta,
    render_normalized_model_gff,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare_reference_loci(reference_gff: Path, loci_dir: Path, output_dir: Path) -> None:
    rows = _read_gff_rows(reference_gff)
    genes = {
        row.attributes["ID"]: row
        for row in rows
        if row.feature_type == "gene" and "ID" in row.attributes
    }
    if not genes:
        raise InputValidationError("Reference GFF contains no identified genes")

    output_dir.mkdir(parents=True, exist_ok=True)
    source_gff_sha256 = _sha256(reference_gff)
    provenance_rows = [
        "model_id\tsource_sequence\tsource_start\tsource_end\tsource_strand\t"
        "locus_sha256\tnormalized_gff_sha256\tsource_gff_sha256"
    ]
    locus_paths = sorted(path for path in loci_dir.iterdir() if path.is_file())
    if not locus_paths:
        raise InputValidationError(f"No reference loci found in {loci_dir}")
    for locus_path in locus_paths:
        record = read_single_fasta(locus_path)
        if record.identifier != locus_path.name:
            raise InputValidationError(
                f"Reference locus filename {locus_path.name} does not match FASTA ID {record.identifier}"
            )
        if record.identifier not in genes:
            raise InputValidationError(
                f"Reference locus {record.identifier} has no matching gene in {reference_gff}"
            )
        annotation = normalize_model_gff(reference_gff, record.identifier, len(record.sequence))
        output_path = output_dir / f"{record.identifier}.gff"
        output_path.write_text(render_normalized_model_gff(annotation))
        gene = genes[record.identifier]
        provenance_rows.append(
            "\t".join(
                (
                    record.identifier,
                    gene.sequence_id,
                    str(gene.start),
                    str(gene.end),
                    gene.strand,
                    _sha256(locus_path),
                    _sha256(output_path),
                    source_gff_sha256,
                )
            )
        )
    (output_dir.parent / "REF_LOCI_PROVENANCE.tsv").write_text(
        "\n".join(provenance_rows) + "\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_gff", type=Path)
    parser.add_argument("loci_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    arguments = parser.parse_args()
    prepare_reference_loci(arguments.reference_gff, arguments.loci_dir, arguments.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
