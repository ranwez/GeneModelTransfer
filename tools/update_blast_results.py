#!/usr/bin/env python3
"""Update a BLAST table by removing and replacing query results."""

import argparse
import logging
import sys
from contextlib import nullcontext
from pathlib import Path
from typing import Optional, TextIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Update a headerless BLAST TSV using qseqid from column 1."
    )
    parser.add_argument(
        "-i", "--init-blast", required=True, type=Path,
        help="Initial BLAST tabular file."
    )
    parser.add_argument(
        "--queries-to-rm", type=Path,
        help="File containing one qseqid per line to remove."
    )
    parser.add_argument(
        "--blast-update", type=Path,
        help="BLAST rows replacing existing queries or adding new queries."
    )
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Output BLAST file. Defaults to stdout; '-' also means stdout."
    )

    args = parser.parse_args()
    if args.queries_to_rm is None and args.blast_update is None:
        parser.error(
            "at least one of --queries-to-rm or --blast-update is required"
        )
    return args


def read_query_ids(path: Path) -> set[str]:
    """Read one query identifier per non-empty line."""
    with path.open(encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def validate_blast_columns(path: Path) -> Optional[int]:
    """Return the column count and reject inconsistent BLAST rows."""
    expected: Optional[int] = None

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue

            fields = line.split("\t")
            if not fields[0].strip():
                raise ValueError(
                    f"{path}:{line_number}: empty qseqid in first column"
                )

            if expected is None:
                expected = len(fields)
            elif len(fields) != expected:
                raise ValueError(
                    f"{path}:{line_number}: found {len(fields)} columns, "
                    f"expected {expected}"
                )

    return expected


def load_blast_update(path: Path) -> tuple[list[str], set[str], int]:
    """Load BLAST update rows and their query identifiers."""
    column_count = validate_blast_columns(path)
    if column_count is None:
        raise ValueError(f"{path}: BLAST update contains no rows")

    lines = []
    query_ids = set()

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line:
                lines.append(line)
                query_ids.add(line.split("\t", 1)[0].strip())

    return lines, query_ids, column_count


def merge_blast(
    init_blast: Path,
    remove_queries: set[str],
    update_queries: set[str],
    update_lines: list[str],
    output_handle: TextIO,
) -> tuple[int, int]:
    """Write unchanged initial rows followed by the BLAST update."""
    removed_rows = 0
    replaced_rows = 0

    with init_blast.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue

            qseqid = line.split("\t", 1)[0].strip()
            if qseqid in remove_queries:
                removed_rows += 1
            elif qseqid in update_queries:
                replaced_rows += 1
            else:
                output_handle.write(f"{line}\n")

    for line in update_lines:
        output_handle.write(f"{line}\n")

    return removed_rows, replaced_rows


def run(args: argparse.Namespace) -> None:
    """Validate inputs and write the updated BLAST table."""
    remove_queries = (
        read_query_ids(args.queries_to_rm)
        if args.queries_to_rm is not None
        else set()
    )

    update_lines: list[str] = []
    update_queries: set[str] = set()
    update_columns: Optional[int] = None

    if args.blast_update is not None:
        update_lines, update_queries, update_columns = load_blast_update(
            args.blast_update
        )

    overlap = remove_queries & update_queries
    if overlap:
        raise ValueError(
            "queries cannot be both removed and updated: "
            + ", ".join(sorted(overlap)[:10])
        )

    if not remove_queries and not update_lines:
        raise ValueError("no query removal or BLAST update was provided")

    init_columns = validate_blast_columns(args.init_blast)
    if (
        init_columns is not None
        and update_columns is not None
        and init_columns != update_columns
    ):
        raise ValueError(
            "initial BLAST and BLAST update have different column counts: "
            f"{init_columns} != {update_columns}"
        )

    if (
        args.output is not None
        and str(args.output) != "-"
        and args.output.resolve() == args.init_blast.resolve()
    ):
        raise ValueError("output must differ from --init-blast")

    output_context = (
        nullcontext(sys.stdout)
        if args.output is None or str(args.output) == "-"
        else args.output.open("w", encoding="utf-8")
    )

    with output_context as output_handle:
        removed_rows, replaced_rows = merge_blast(
            args.init_blast,
            remove_queries,
            update_queries,
            update_lines,
            output_handle,
        )

    logging.info(
        "Removed %d row(s), replaced %d row(s), added %d update row(s).",
        removed_rows,
        replaced_rows,
        len(update_lines),
    )


def main() -> int:
    """Run the command-line interface."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()

    try:
        run(args)
    except (OSError, ValueError) as exc:
        logging.error("%s", exc)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
