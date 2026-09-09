from __future__ import annotations

import csv
import hashlib
from pathlib import Path

from conftest import REPOSITORY_ROOT, TEST_DIR


DATA_DIR = TEST_DIR / "data"
SOURCE_GFF = REPOSITORY_ROOT / "HIPP/HIPP_NIP_KIT/OsjNip_HIPPHPP_sample.gff"


def _fasta_length(path: Path) -> int:
    lines = path.read_text().splitlines()
    assert len([line for line in lines if line.startswith(">")]) == 1
    return sum(len(line.strip()) for line in lines if not line.startswith(">"))


def test_pinned_hipp_sources_have_expected_checksums():
    for line in (DATA_DIR / "hipp_sources.sha256").read_text().splitlines():
        expected, relative_path = line.split(maxsplit=1)
        content = (REPOSITORY_ROOT / relative_path).read_bytes()
        assert hashlib.sha256(content).hexdigest() == expected


def test_functional_case_manifest_references_single_record_fastas():
    with (DATA_DIR / "hipp_cases.tsv").open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    assert [row["case_id"] for row in rows] == ["FT-01", "FT-02", "FT-03", "FT-04"]
    for row in rows:
        assert _fasta_length(REPOSITORY_ROOT / row["model_fasta"]) == int(row["verified_best_qlen"])
        assert _fasta_length(REPOSITORY_ROOT / row["target_fasta"]) > 0


def test_authoritative_gff_contains_complete_selected_models():
    expected_cds_counts = {
        "OsjNip_Chr04_02891360": 5,
        "OsjNip_Chr01_23324465": 3,
        "OsjNip_Chr04_19814005": 2,
    }
    records = [
        line.split("\t")
        for line in SOURCE_GFF.read_text().splitlines()
        if line and not line.startswith("#")
    ]

    for model_id, expected_cds_count in expected_cds_counts.items():
        gene = [row for row in records if row[2] == "gene" and f"ID={model_id}" in row[8]]
        mrna = [row for row in records if row[2] == "mRNA" and f"ID={model_id}_rna" in row[8]]
        cds = [row for row in records if row[2] == "CDS" and f"Parent={model_id}_rna" in row[8]]
        assert len(gene) == 1
        assert len(mrna) == 1
        assert len(cds) == expected_cds_count
