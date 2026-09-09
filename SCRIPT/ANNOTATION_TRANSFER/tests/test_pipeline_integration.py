from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

from conftest import REPOSITORY_ROOT


GENE_PREDICTION = REPOSITORY_ROOT / "bin/genePrediction.sh"
CREATE_LRROME = REPOSITORY_ROOT / "bin/create_LRRome.sh"
SNAKEFILE = REPOSITORY_ROOT / "SNAKEMAKE/LRR_transfer.smk"


@pytest.mark.integration
def test_pipeline_declares_and_converts_locus_alignment_exactly_once():
    script = GENE_PREDICTION.read_text()
    match = re.search(r'methods="([^"]+)"', script)
    assert match is not None, 'No methods="..." declaration found'
    methods = match.group(1).split()
    assert "locusAlignment" in methods

    snakefile = SNAKEFILE.read_text()
    assert "annotate_one_{split_id}_locusAlignment.gff" in snakefile
    assert "annot_locusAlignment.gff" in snakefile


@pytest.mark.integration
@pytest.mark.parametrize("script", [GENE_PREDICTION, CREATE_LRROME])
def test_modified_shell_entrypoints_remain_syntactically_valid(script: Path):
    completed = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr


@pytest.mark.integration
def test_lrrome_builder_creates_exact_locus_resources_and_provenance():
    script = CREATE_LRROME.read_text()
    assert "CANDIDATE_LOCI/extract_loci.sh" in script
    assert "REF_LOCI_GFF" in script
    assert "prepare_reference_loci.py" in script
