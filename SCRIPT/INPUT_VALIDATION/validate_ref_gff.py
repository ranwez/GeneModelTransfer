#!/usr/bin/env python3
"""Validate the reference GFF and its reference-side identifiers."""

import argparse
import logging
from pathlib import Path
from typing import Optional

from attrs import define, field


class ValidationError(ValueError):
    """Raised when an input validation rule is violated."""


@define(frozen=True, slots=True)
class GffFeature:
    """Represent a GFF feature required for reference model validation."""

    line_number: int
    seqid: str
    feature_type: str
    start: int
    end: int
    strand: str
    attributes: dict[str, list[str]]


@define(slots=True)
class ReferenceGffValidator:
    """Validate reference gene models and related reference inputs."""

    gff_path: Path
    ref_fai_path: Optional[Path] = None
    ref_locus_info_path: Optional[Path] = None

    features: list[GffFeature] = field(factory=list, init=False)
    genes: dict[str, GffFeature] = field(factory=dict, init=False)
    mrnas: dict[str, GffFeature] = field(factory=dict, init=False)
    cds: list[GffFeature] = field(factory=list, init=False)
    mrnas_by_gene: dict[str, list[GffFeature]] = field(factory=dict, init=False)
    cds_by_mrna: dict[str, list[GffFeature]] = field(factory=dict, init=False)
    seqids: set[str] = field(factory=set, init=False)

    def validate(self) -> None:
        """Run all applicable reference GFF validations."""
        self._read_gff()
        self._index_models()
        self._validate_models()
        self._validate_reference_seqids()
        self._validate_info_locus()

    def write_report(self, output_path: Path) -> None:
        """Write a compact validation summary."""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        rows = [
            ("status", "OK"),
            ("gff", str(self.gff_path)),
            ("genes", str(len(self.genes))),
            ("mrnas", str(len(self.mrnas))),
            ("cds", str(len(self.cds))),
            ("seqids", str(len(self.seqids))),
            (
                "ref_genome_seqids",
                "OK" if self.ref_fai_path is not None else "SKIPPED",
            ),
            (
                "ref_locus_info_gene_ids",
                "OK" if self.ref_locus_info_path is not None else "SKIPPED",
            ),
        ]

        with output_path.open("w", encoding="utf-8") as handle:
            for key, value in rows:
                handle.write(f"{key}\t{value}\n")

    def _read_gff(self) -> None:
        """Read GFF feature rows and validate their basic structure."""
        with self.gff_path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")

                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) != 9:
                    raise ValidationError(
                        f"{self.gff_path}:{line_number}: expected 9 columns, "
                        f"found {len(fields)}"
                    )

                try:
                    start = int(fields[3])
                    end = int(fields[4])
                except ValueError as exc:
                    raise ValidationError(
                        f"{self.gff_path}:{line_number}: "
                        "start and end must be integers"
                    ) from exc

                if start < 1:
                    raise ValidationError(
                        f"{self.gff_path}:{line_number}: start must be >= 1"
                    )

                if start > end:
                    raise ValidationError(
                        f"{self.gff_path}:{line_number}: start must be <= end"
                    )

                feature = GffFeature(
                    line_number=line_number,
                    seqid=fields[0],
                    feature_type=fields[2],
                    start=start,
                    end=end,
                    strand=fields[6],
                    attributes=self._parse_attributes(fields[8]),
                )

                self.features.append(feature)
                self.seqids.add(feature.seqid)

        if not self.features:
            raise ValidationError(
                f"{self.gff_path}: no GFF feature rows found"
            )

    def _index_models(self) -> None:
        """Index genes, mRNAs, CDS features and their parent relationships."""
        for feature in self.features:
            if feature.feature_type == "gene":
                gene_id = self._single_attribute(feature, "ID")
                self._add_unique_feature(self.genes, gene_id, feature, "gene")
                self._validate_strand(feature)

            elif feature.feature_type == "mRNA":
                mrna_id = self._single_attribute(feature, "ID")
                self._single_attribute(feature, "Parent")
                self._add_unique_feature(self.mrnas, mrna_id, feature, "mRNA")

            elif feature.feature_type == "CDS":
                self._single_attribute(feature, "Parent")
                self.cds.append(feature)

        if not self.genes:
            raise ValidationError(
                f"{self.gff_path}: no gene features found"
            )

        for mrna in self.mrnas.values():
            gene_id = self._single_attribute(mrna, "Parent")

            if gene_id not in self.genes:
                raise ValidationError(
                    f"{self.gff_path}:{mrna.line_number}: "
                    f"mRNA Parent {gene_id!r} does not match any gene"
                )

            self.mrnas_by_gene.setdefault(gene_id, []).append(mrna)

        for cds in self.cds:
            mrna_id = self._single_attribute(cds, "Parent")

            if mrna_id not in self.mrnas:
                raise ValidationError(
                    f"{self.gff_path}:{cds.line_number}: "
                    f"CDS Parent {mrna_id!r} does not match any mRNA"
                )

            self.cds_by_mrna.setdefault(mrna_id, []).append(cds)

    def _validate_models(self) -> None:
        """Validate gene, mRNA and CDS model relationships."""
        for gene_id, gene in self.genes.items():
            mrnas = self.mrnas_by_gene.get(gene_id, [])

            if len(mrnas) != 1:
                raise ValidationError(
                    f"{self.gff_path}: gene {gene_id!r} must have exactly "
                    f"one mRNA; found {len(mrnas)}"
                )

            mrna = mrnas[0]
            mrna_id = self._single_attribute(mrna, "ID")

            self._validate_child(mrna, gene, "mRNA")

            cds_features = self.cds_by_mrna.get(mrna_id, [])
            if not cds_features:
                raise ValidationError(
                    f"{self.gff_path}: mRNA {mrna_id!r} has no CDS"
                )

            for cds in cds_features:
                self._validate_child(cds, mrna, "CDS")

            self._validate_cds_overlap(mrna_id, cds_features)

    def _validate_child(
        self,
        child: GffFeature,
        parent: GffFeature,
        label: str,
    ) -> None:
        """Validate sequence, strand and containment of a child feature."""
        self._validate_strand(child)

        if child.seqid != parent.seqid:
            raise ValidationError(
                f"{self.gff_path}:{child.line_number}: "
                f"{label} seqid {child.seqid!r} differs from its parent "
                f"seqid {parent.seqid!r}"
            )

        if child.strand != parent.strand:
            raise ValidationError(
                f"{self.gff_path}:{child.line_number}: "
                f"{label} strand {child.strand!r} differs from its parent "
                f"strand {parent.strand!r}"
            )

        if child.start < parent.start or child.end > parent.end:
            raise ValidationError(
                f"{self.gff_path}:{child.line_number}: "
                f"{label} coordinates {child.start}-{child.end} are not "
                f"contained in parent coordinates "
                f"{parent.start}-{parent.end}"
            )

    def _validate_cds_overlap(
        self,
        mrna_id: str,
        cds_features: list[GffFeature],
    ) -> None:
        """Ensure CDS features from the same mRNA do not overlap."""
        ordered = sorted(cds_features, key=lambda feature: feature.start)

        for previous, current in zip(ordered, ordered[1:]):
            if current.start <= previous.end:
                raise ValidationError(
                    f"{self.gff_path}: overlapping CDS features for "
                    f"mRNA {mrna_id!r} at lines "
                    f"{previous.line_number} and {current.line_number}"
                )

    def _validate_reference_seqids(self) -> None:
        """Check that every GFF seqid exists in the reference FASTA index."""
        if self.ref_fai_path is None:
            return

        fasta_seqids = self._read_fai_seqids(self.ref_fai_path)
        missing = sorted(self.seqids - fasta_seqids)

        if missing:
            raise ValidationError(
                f"{self.gff_path}: {len(missing)} seqid(s) are absent from "
                f"{self.ref_fai_path}: {self._format_examples(missing)}"
            )

    def _validate_info_locus(self) -> None:
        """Check that every GFF gene is represented in the locus information."""
        if self.ref_locus_info_path is None:
            return

        info_gene_ids = self._read_info_locus_ids(self.ref_locus_info_path)
        missing = sorted(set(self.genes) - info_gene_ids)

        if missing:
            raise ValidationError(
                f"{self.gff_path}: {len(missing)} gene ID(s) are absent from "
                f"{self.ref_locus_info_path}: "
                f"{self._format_examples(missing)}"
            )

    def _single_attribute(
        self,
        feature: GffFeature,
        attribute: str,
    ) -> str:
        """Return an attribute that must contain exactly one value."""
        values = feature.attributes.get(attribute, [])

        if len(values) != 1 or not values[0]:
            raise ValidationError(
                f"{self.gff_path}:{feature.line_number}: "
                f"{feature.feature_type} must have exactly one "
                f"{attribute} attribute"
            )

        value = values[0]

        if "," in value:
            raise ValidationError(
                f"{self.gff_path}:{feature.line_number}: "
                f"{feature.feature_type} {attribute} must reference "
                "exactly one value"
            )

        return value

    def _add_unique_feature(
        self,
        index: dict[str, GffFeature],
        feature_id: str,
        feature: GffFeature,
        label: str,
    ) -> None:
        """Add a feature to an ID index and reject duplicate IDs."""
        if feature_id in index:
            raise ValidationError(
                f"{self.gff_path}:{feature.line_number}: duplicate "
                f"{label} ID {feature_id!r}"
            )

        index[feature_id] = feature

    def _validate_strand(self, feature: GffFeature) -> None:
        """Require a forward or reverse strand for model features."""
        if feature.strand not in {"+", "-"}:
            raise ValidationError(
                f"{self.gff_path}:{feature.line_number}: "
                f"{feature.feature_type} strand must be '+' or '-'"
            )

    @staticmethod
    def _parse_attributes(raw_attributes: str) -> dict[str, list[str]]:
        """Parse GFF attributes while preserving repeated keys."""
        attributes: dict[str, list[str]] = {}

        for item in raw_attributes.split(";"):
            item = item.strip()

            if not item or "=" not in item:
                continue

            key, value = item.split("=", 1)
            attributes.setdefault(key.strip(), []).append(value.strip())

        return attributes

    @staticmethod
    def _read_fai_seqids(fai_path: Path) -> set[str]:
        """Read sequence identifiers from a FASTA index."""
        seqids = set()

        with fai_path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\r\n").split("\t")

                if len(fields) < 2:
                    raise ValidationError(
                        f"{fai_path}:{line_number}: malformed FAI row"
                    )

                seqids.add(fields[0])

        return seqids

    @staticmethod
    def _read_info_locus_ids(info_path: Path) -> set[str]:
        """Read gene identifiers from the first TSV column."""
        gene_ids = set()

        with info_path.open(encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                gene_ids.add(line.split("\t", 1)[0])

        return gene_ids

    @staticmethod
    def _format_examples(values: list[str], limit: int = 10) -> str:
        """Format a compact subset of invalid identifiers."""
        examples = ", ".join(values[:limit])

        if len(values) > limit:
            examples += ", ..."

        return examples


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate an LRRtransfer reference GFF."
    )
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--ref-fai", default="")
    parser.add_argument("--ref-locus-info", default="")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def optional_path(value: str) -> Optional[Path]:
    """Convert a non-empty argument to a Path."""
    return Path(value) if value else None


def configure_logging() -> None:
    """Configure command-line logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )


def run_validation(args: argparse.Namespace) -> int:
    """Run validation and return a process exit code."""
    validator = ReferenceGffValidator(
        gff_path=args.gff,
        ref_fai_path=optional_path(args.ref_fai),
        ref_locus_info_path=optional_path(args.ref_locus_info),
    )

    try:
        validator.validate()
        validator.write_report(args.output)
    except (ValidationError, OSError) as exc:
        logging.error("%s", exc)
        return 1

    logging.info("Reference GFF validation passed.")
    return 0


def main() -> int:
    """Run the command-line interface."""
    args = parse_args()
    configure_logging()
    return run_validation(args)


if __name__ == "__main__":
    raise SystemExit(main())
