#!/usr/bin/env bash
set -euo pipefail
python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests
python -m pytest -q SCRIPT/CANDIDATE_LOCI/tests