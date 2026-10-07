#!/usr/bin/env python3
"""Validate external BLAST result identifiers against workflow inputs."""

import argparse
import logging
from pathlib import Path
from typing import Optional

from attrs import define, field


class ValidationError(ValueError):
    """Raised when an input validation rule is violated."""


@define(slots=True)
class BlastResultsValidator:
    """Validate query and subject identifiers in external BLAST results."""

    reference_gff_path: Path
    target_fai_path: Path
    tblastn_path: Optional[Path] = None
    blastn_path: Optional[Path] = None

    summaries: dict[str, tuple[int, int, int]] = field(
        factory=dict,
        init=False,
    )

    def validate(self) -> None:
        """Run all applicable BLAST identifier checks."""
        if self.tblastn_path is None and self.blastn_path is None:
            raise ValidationError(
                "At least one BLAST result file must be provided"
            )

        gene_ids = self._read_gff_gene_ids(self.reference_gff_path)
        target_seqids = self._read_fai_seqids(self.target_fai_path)

        if self.tblastn_path is not None:
            self._validate_blast_file(
                label="TBLASTN",
                blast_path=self.tblastn_path,
                reference_gene_ids=gene_ids,
                target_seqids=target_seqids,
            )

        if self.blastn_path is not None:
            self._validate_blast_file(
                label="BLASTN",
                blast_path=self.blastn_path,
                reference_gene_ids=gene_ids,
                target_seqids=target_seqids,
            )

    def write_report(self, output_path: Path) -> None:
        """Write a compact BLAST validation summary."""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        rows = [("status", "OK")]

        for label in ("TBLASTN", "BLASTN"):
            key = label.lower()

            if label not in self.summaries:
                rows.append((key, "SKIPPED"))
                continue

            row_count, query_count, subject_count = self.summaries[label]
            rows.extend(
                [
                    (key, "OK"),
                    (f"{key}_rows", str(row_count)),
                    (f"{key}_query_ids", str(query_count)),
                    (f"{key}_subject_seqids", str(subject_count)),
                ]
            )

        with output_path.open("w", encoding="utf-8") as handle:
            for key, value in rows:
                handle.write(f"{key}\t{value}\n")

    def _validate_blast_file(
        self,
        label: str,
        blast_path: Path,
        reference_gene_ids: set[str],
        target_seqids: set[str],
    ) -> None:
        """Validate query and subject identifiers from one BLAST table."""
        row_count, query_ids, subject_ids = self._read_blast_ids(
            blast_path
        )

        missing_queries = sorted(query_ids - reference_gene_ids)
        if missing_queries:
            raise ValidationError(
                f"{blast_path}: {len(missing_queries)} {label} "
                "query ID(s) are absent from the reference GFF: "
                f"{self._format_examples(missing_queries)}"
            )

        missing_subjects = sorted(subject_ids - target_seqids)
        if missing_subjects:
            raise ValidationError(
                f"{blast_path}: {len(missing_subjects)} {label} "
                "subject seqid(s) are absent from the target genome: "
                f"{self._format_examples(missing_subjects)}"
            )

        self.summaries[label] = (
            row_count,
            len(query_ids),
            len(subject_ids),
        )

    @staticmethod
    def _read_blast_ids(
        blast_path: Path,
    ) -> tuple[int, set[str], set[str]]:
        """Read qseqid and sseqid from a BLAST tabular file."""
        query_ids = set()
        subject_ids = set()
        row_count = 0

        with blast_path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")

                if not line:
                    continue

                fields = line.split("\t")

                if len(fields) < 2:
                    raise ValidationError(
                        f"{blast_path}:{line_number}: "
                        "BLAST row must contain at least qseqid and sseqid"
                    )

                qseqid = fields[0].strip()
                sseqid = fields[1].strip()

                if not qseqid or not sseqid:
                    raise ValidationError(
                        f"{blast_path}:{line_number}: "
                        "qseqid and sseqid must not be empty"
                    )

                query_ids.add(qseqid)
                subject_ids.add(sseqid)
                row_count += 1

        if row_count == 0:
            raise ValidationError(
                f"{blast_path}: no BLAST result rows found"
            )

        return row_count, query_ids, subject_ids

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

                gene_id = BlastResultsValidator._get_gff_attribute(
                    fields[8],
                    "ID",
                )

                if gene_id is None:
                    raise ValidationError(
                        f"{gff_path}:{line_number}: gene has no ID"
                    )

                gene_ids.add(gene_id)

        if not gene_ids:
            raise ValidationError(
                f"{gff_path}: no gene identifiers found"
            )

        return gene_ids

    @staticmethod
    def _read_fai_seqids(fai_path: Path) -> set[str]:
        """Read target sequence identifiers from a FASTA index."""
        seqids = set()

        with fai_path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                fields = raw_line.rstrip("\r\n").split("\t")

                if len(fields) < 2:
                    raise ValidationError(
                        f"{fai_path}:{line_number}: malformed FAI row"
                    )

                seqids.add(fields[0])

        if not seqids:
            raise ValidationError(
                f"{fai_path}: no sequence identifiers found"
            )

        return seqids

    @staticmethod
    def _get_gff_attribute(
        attributes: str,
        key: str,
    ) -> Optional[str]:
        """Return a GFF attribute value."""
        prefix = f"{key}="

        for item in attributes.split(";"):
            item = item.strip()

            if item.startswith(prefix):
                return item[len(prefix):]

        return None

    @staticmethod
    def _format_examples(
        values: list[str],
        limit: int = 10,
    ) -> str:
        """Format a compact subset of invalid identifiers."""
        examples = ", ".join(values[:limit])

        if len(values) > limit:
            examples += ", ..."

        return examples


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate external TBLASTN and BLASTN result identifiers."
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--target-fai",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--tblastn",
        default="",
    )
    parser.add_argument(
        "--blastn",
        default="",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )
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
    validator = BlastResultsValidator(
        reference_gff_path=args.gff,
        target_fai_path=args.target_fai,
        tblastn_path=optional_path(args.tblastn),
        blastn_path=optional_path(args.blastn),
    )

    try:
        validator.validate()
        validator.write_report(args.output)
    except (ValidationError, OSError) as exc:
        logging.error("%s", exc)
        return 1

    logging.info("BLAST result validation passed.")
    return 0


def main() -> int:
    """Run the command-line interface."""
    args = parse_args()
    configure_logging()
    return run_validation(args)


if __name__ == "__main__":
    raise SystemExit(main())
