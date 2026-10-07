#!/usr/bin/env python3
"""Validate an external LRRome and its consistency with the reference GFF."""

import argparse
import logging
from pathlib import Path
from typing import Optional

from attrs import define, field


class ValidationError(ValueError):
    """Raised when an input validation rule is violated."""


@define(slots=True)
class LrromeValidator:
    """Validate the structure and model identifiers of an LRRome."""

    lrrome_path: Path
    reference_gff_path: Path
    model_ids: set[str] = field(factory=set, init=False)
    exon_file_count: int = field(default=0, init=False)

    def validate(self) -> None:
        """Run all LRRome validation checks."""
        self._validate_required_components()

        self.model_ids = self._read_fasta_ids(
            self.lrrome_path / "REF_proteins.fasta"
        )

        if not self.model_ids:
            raise ValidationError(
                f"{self.lrrome_path}/REF_proteins.fasta contains no sequences"
            )

        self._validate_fasta_ids(
            self.lrrome_path / "REF_loci.fasta",
            self.model_ids,
            "REF_loci.fasta",
        )

        self._validate_model_directory(
            self.lrrome_path / "REF_PEP",
            self.model_ids,
        )
        self._validate_model_directory(
            self.lrrome_path / "REF_cDNA",
            self.model_ids,
        )
        self._validate_model_directory(
            self.lrrome_path / "REF_LOCI",
            self.model_ids,
        )

        self._validate_gff_directory()
        self._validate_exon_directory()
        self._validate_reference_gff()

    def write_report(self, output_path: Path) -> None:
        """Write a compact LRRome validation summary."""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        rows = [
            ("status", "OK"),
            ("lrrome", str(self.lrrome_path)),
            ("models", str(len(self.model_ids))),
            ("exon_files", str(self.exon_file_count)),
            ("ref_gff_gene_ids", "OK"),
        ]

        with output_path.open("w", encoding="utf-8") as handle:
            for key, value in rows:
                handle.write(f"{key}\t{value}\n")

    def _validate_required_components(self) -> None:
        """Check that all required LRRome files and directories exist."""
        required_files = [
            "REF_proteins.fasta",
            "REF_loci.fasta",
        ]
        required_directories = [
            "REF_PEP",
            "REF_EXONS",
            "REF_cDNA",
            "REF_LOCI",
            "REF_LOCI_GFF",
        ]

        for name in required_files:
            path = self.lrrome_path / name

            if not path.is_file():
                raise ValidationError(
                    f"Missing LRRome file: {path}"
                )

            if path.stat().st_size == 0:
                raise ValidationError(
                    f"LRRome file is empty: {path}"
                )

        for name in required_directories:
            path = self.lrrome_path / name

            if not path.is_dir():
                raise ValidationError(
                    f"Missing LRRome directory: {path}"
                )

    def _validate_fasta_ids(
        self,
        fasta_path: Path,
        expected_ids: set[str],
        label: str,
    ) -> None:
        """Check that FASTA identifiers exactly match expected model IDs."""
        observed_ids = self._read_fasta_ids(fasta_path)
        self._require_equal_ids(expected_ids, observed_ids, label)

    def _validate_model_directory(
        self,
        directory: Path,
        expected_ids: set[str],
    ) -> None:
        """Check per-model FASTA files and their internal identifiers."""
        observed_ids = self._directory_file_names(directory)
        self._require_equal_ids(
            expected_ids,
            observed_ids,
            directory.name,
        )

        for model_id in expected_ids:
            path = directory / model_id
            fasta_ids = self._read_fasta_ids(path)

            if fasta_ids != {model_id}:
                raise ValidationError(
                    f"{path}: expected a single FASTA record named "
                    f"{model_id!r}, found {sorted(fasta_ids)}"
                )

    def _validate_gff_directory(self) -> None:
        """Check that REF_LOCI_GFF contains one GFF per model."""
        directory = self.lrrome_path / "REF_LOCI_GFF"
        observed_files = self._directory_file_names(directory)
        expected_files = {
            f"{model_id}.gff"
            for model_id in self.model_ids
        }

        self._require_equal_ids(
            expected_files,
            observed_files,
            "REF_LOCI_GFF",
        )

    def _validate_exon_directory(self) -> None:
        """Check that every model has at least one reference exon file."""
        directory = self.lrrome_path / "REF_EXONS"
        exon_files = self._directory_file_names(directory)

        if not exon_files:
            raise ValidationError(
                f"{directory}: no exon files found"
            )

        model_counts = {
            model_id: 0
            for model_id in self.model_ids
        }

        for filename in exon_files:
            model_id = self._match_exon_model(filename)

            if model_id is None:
                raise ValidationError(
                    f"{directory}: exon file {filename!r} does not "
                    "match any LRRome model"
                )

            fasta_ids = self._read_fasta_ids(directory / filename)

            if fasta_ids != {filename}:
                raise ValidationError(
                    f"{directory / filename}: expected a single FASTA "
                    f"record named {filename!r}, found {sorted(fasta_ids)}"
                )

            model_counts[model_id] += 1

        missing = sorted(
            model_id
            for model_id, count in model_counts.items()
            if count == 0
        )

        if missing:
            raise ValidationError(
                f"{directory}: {len(missing)} model(s) have no exon "
                f"file: {self._format_examples(missing)}"
            )

        self.exon_file_count = len(exon_files)

    def _validate_reference_gff(self) -> None:
        """Check that LRRome model IDs exactly match GFF gene IDs."""
        gff_gene_ids = self._read_gff_gene_ids(
            self.reference_gff_path
        )

        self._require_equal_ids(
            gff_gene_ids,
            self.model_ids,
            "reference GFF vs LRRome",
        )

    def _match_exon_model(self, filename: str) -> Optional[str]:
        """Return the model corresponding to an exon filename."""
        matches = [
            model_id
            for model_id in self.model_ids
            if any(
                filename.startswith(f"{model_id}{separator}")
                for separator in ("_", ":", "-")
            )
        ]

        if not matches:
            return None

        return max(matches, key=len)

    @staticmethod
    def _read_fasta_ids(fasta_path: Path) -> set[str]:
        """Read unique identifiers from a FASTA file."""
        identifiers = []

        with fasta_path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue

                identifier = line[1:].strip().split(maxsplit=1)[0]

                if not identifier:
                    raise ValidationError(
                        f"{fasta_path}: empty FASTA identifier"
                    )

                identifiers.append(identifier)

        if len(identifiers) != len(set(identifiers)):
            raise ValidationError(
                f"{fasta_path}: duplicate FASTA identifiers"
            )

        return set(identifiers)

    @staticmethod
    def _directory_file_names(directory: Path) -> set[str]:
        """Return non-hidden regular filenames from a directory."""
        names = set()

        for path in directory.iterdir():
            if path.name.startswith("."):
                continue

            if not path.is_file():
                raise ValidationError(
                    f"{directory}: unexpected directory entry {path.name!r}"
                )

            if path.stat().st_size == 0:
                raise ValidationError(
                    f"{directory}: empty file {path.name!r}"
                )

            names.add(path.name)

        return names

    @staticmethod
    def _read_gff_gene_ids(gff_path: Path) -> set[str]:
        """Read gene identifiers from a reference GFF."""
        gene_ids = set()

        with gff_path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")

                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) != 9 or fields[2] != "gene":
                    continue

                gene_id = LrromeValidator._get_gff_attribute(
                    fields[8],
                    "ID",
                )

                if gene_id is None:
                    raise ValidationError(
                        f"{gff_path}:{line_number}: gene has no ID"
                    )

                gene_ids.add(gene_id)

        return gene_ids

    @staticmethod
    def _get_gff_attribute(
        attributes: str,
        key: str,
    ) -> Optional[str]:
        """Return a GFF attribute value."""
        prefix = f"{key}="

        for field_value in attributes.split(";"):
            field_value = field_value.strip()

            if field_value.startswith(prefix):
                return field_value[len(prefix):]

        return None

    @staticmethod
    def _require_equal_ids(
        expected: set[str],
        observed: set[str],
        label: str,
    ) -> None:
        """Require equality between two identifier sets."""
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)

        if not missing and not extra:
            return

        details = []

        if missing:
            details.append(
                "missing: "
                + LrromeValidator._format_examples(missing)
            )

        if extra:
            details.append(
                "unexpected: "
                + LrromeValidator._format_examples(extra)
            )

        raise ValidationError(
            f"{label}: identifier mismatch ({'; '.join(details)})"
        )

    @staticmethod
    def _format_examples(
        values: list[str],
        limit: int = 10,
    ) -> str:
        """Format a compact subset of identifiers."""
        examples = ", ".join(values[:limit])

        if len(values) > limit:
            examples += ", ..."

        return examples


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate an external LRRome."
    )
    parser.add_argument(
        "--lrrome",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure command-line logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )


def run_validation(args: argparse.Namespace) -> int:
    """Run validation and return a process exit code."""
    validator = LrromeValidator(
        lrrome_path=args.lrrome,
        reference_gff_path=args.gff,
    )

    try:
        validator.validate()
        validator.write_report(args.output)
    except (ValidationError, OSError) as exc:
        logging.error("%s", exc)
        return 1

    logging.info("LRRome validation passed.")
    return 0


def main() -> int:
    """Run the command-line interface."""
    args = parse_args()
    configure_logging()
    return run_validation(args)


if __name__ == "__main__":
    raise SystemExit(main())
