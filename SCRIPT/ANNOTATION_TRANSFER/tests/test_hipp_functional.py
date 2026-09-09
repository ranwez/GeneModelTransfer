from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from conftest import REPOSITORY_ROOT, TEST_DIR
from locus_alignment_transfer import (
    normalize_model_gff,
    read_single_fasta,
    render_normalized_model_gff,
    run_transfer,
)


SOURCE_GFF = REPOSITORY_ROOT / "HIPP/HIPP_NIP_KIT/OsjNip_HIPPHPP_sample.gff"


def _cds_rows(path: Path):
    return [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#") and line.split("\t")[2] == "CDS"
    ]


def _translated_protein(fasta_path: Path, gff_path: Path) -> str:
    sequence = read_single_fasta(fasta_path).sequence
    rows = _cds_rows(gff_path)
    assert rows
    strand = rows[0][6]
    rows.sort(key=lambda row: int(row[3]), reverse=strand == "-")
    coding = "".join(sequence[int(row[3]) - 1 : int(row[4])] for row in rows)
    if strand == "-":
        complement = str.maketrans("ACGT", "TGCA")
        coding = coding.translate(complement)[::-1]
    assert len(coding) % 3 == 0
    bases = "TCAG"
    amino_acids = (
        "FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG"
    )
    codon_table = dict(
        zip((a + b + c for a in bases for b in bases for c in bases), amino_acids)
    )
    return "".join(codon_table.get(coding[index : index + 3], "X") for index in range(0, len(coding), 3))


def _needleman_wunsch_identity(model: str, predicted: str) -> float:
    """Global Needleman-Wunsch identity divided by model protein length."""

    gap = -2
    scores = [[0] * (len(predicted) + 1) for _ in range(len(model) + 1)]
    trace = [[0] * (len(predicted) + 1) for _ in range(len(model) + 1)]
    for row in range(1, len(model) + 1):
        scores[row][0] = row * gap
        trace[row][0] = 1
    for column in range(1, len(predicted) + 1):
        scores[0][column] = column * gap
        trace[0][column] = 2
    for row in range(1, len(model) + 1):
        for column in range(1, len(predicted) + 1):
            options = (
                scores[row - 1][column - 1]
                + (2 if model[row - 1] == predicted[column - 1] else -1),
                scores[row - 1][column] + gap,
                scores[row][column - 1] + gap,
            )
            scores[row][column] = max(options)
            trace[row][column] = options.index(scores[row][column])
    identical = 0
    row, column = len(model), len(predicted)
    while row or column:
        direction = trace[row][column]
        if row and column and direction == 0:
            identical += model[row - 1] == predicted[column - 1]
            row -= 1
            column -= 1
        elif row and (not column or direction == 1):
            row -= 1
        else:
            column -= 1
    return identical / len(model)


@pytest.mark.external_tools
def test_hipp_cases_meet_nucleotide_and_protein_contracts(tmp_path):
    assert shutil.which("blastn"), "BLAST+ is required for the external contract"
    assert shutil.which("mafft"), "MAFFT is required for the external contract"
    with (TEST_DIR / "data/hipp_cases.tsv").open(newline="") as stream:
        cases = list(csv.DictReader(stream, delimiter="\t"))

    for case in cases:
        model_fasta = REPOSITORY_ROOT / case["model_fasta"]
        target_fasta = REPOSITORY_ROOT / case["target_fasta"]
        model_record = read_single_fasta(model_fasta)
        normalized = normalize_model_gff(SOURCE_GFF, case["model_id"], len(model_record.sequence))
        normalized_gff = tmp_path / f"{case['case_id']}.model.gff"
        normalized_gff.write_text(render_normalized_model_gff(normalized))
        output_gff = tmp_path / f"{case['case_id']}.gff"
        diagnostics = tmp_path / f"{case['case_id']}.json"

        result = run_transfer(
            model_fasta,
            normalized_gff,
            target_fasta,
            output_gff,
            diagnostics,
        )
        expected_eligible = case["expected_eligible"] == "true"
        assert (result.status == "success") is expected_eligible
        if not expected_eligible:
            assert result.reason == "no_significant_hit"
            assert output_gff.stat().st_size == 0
            status = json.loads(diagnostics.read_text())
            assert status["tool_exit_states"]["mafft"] == "not_invoked"
            assert status["versions"]["mafft"] == "not_invoked"
            continue

        assert case["expected_orientation"] == "forward"
        model_protein = _translated_protein(model_fasta, normalized_gff)
        predicted_protein = _translated_protein(target_fasta, output_gff)
        assert "*" not in predicted_protein[:-1]
        identity = _needleman_wunsch_identity(model_protein, predicted_protein)
        assert identity > 0.95
        if case["case_id"] == "FT-01":
            assert [row[3:8] for row in _cds_rows(output_gff)] == [
                row[3:8] for row in _cds_rows(normalized_gff)
            ]
            assert predicted_protein == model_protein
