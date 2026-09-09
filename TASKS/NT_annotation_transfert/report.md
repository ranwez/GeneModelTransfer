# Task Completion Report: Nucleotide locus-alignment annotation transfer

**COMPLETED**: 2026-07-15  
**SPECIFICATION**: [`planning_NT_transfert_V3.md`](planning_NT_transfert_V3.md)

---

## Summary

Implemented the standalone `locusAlignment` Python method with strict input validation,
deterministic single-HSP BLAST gating, sparse MAFFT coordinate mapping, in-frame boundary
fallbacks, reverse-orientation support, atomic GFF/JSON outputs, and explicit
success/ineligible/error semantics. Added exact-locus LRRome resource generation and a
minimal pipeline integration that converts the standalone locus-relative draft to genome
coordinates exactly once before the existing improvement/scoring loop.

## Files Created

- `SCRIPT/ANNOTATION_TRANSFER/locus_alignment_transfer.py` - standalone core, injectable
  orchestration boundaries, real BLAST+/MAFFT adapters, diagnostics, and CLI.
- `SCRIPT/ANNOTATION_TRANSFER/prepare_reference_loci.py` - normalized per-model
  `REF_LOCI_GFF` generation plus source/locus/GFF checksum provenance.
- `SCRIPT/ANNOTATION_TRANSFER/tests/test_hipp_functional.py` - real BLAST+/MAFFT HIPP
  decisions and protein-level acceptance tests.
- `SCRIPT/ANNOTATION_TRANSFER/tests/test_pipeline_integration.py` - resource, shell syntax,
  method declaration, and exactly-once coordinate-conversion integration checks.
- `SCRIPT/ANNOTATION_TRANSFER/tests/test_reverse_transfer.py` - synthetic reverse-HSP
  orchestration and negative-strand projection regression.
- `TASKS/NT_annotation_transfert/report.md` - this report.

## Files Modified

- `bin/create_LRRome.sh` - creates exact strand-aware `REF_LOCI`, normalized
  `REF_LOCI_GFF`, and provenance resources.
- `bin/genePrediction.sh` - adds `locusAlignment`, invokes the standalone CLI, performs
  one locus-to-genome conversion, and includes the method in empty-target and
  `best2rounds` method loops.
- `SNAKEMAKE/LRR_transfer.smk` - declares per-pair and merged `locusAlignment` outputs.
- `SCRIPT/ANNOTATION_TRANSFER/tests/test_locus_alignment_transfer.py` - adds coverage for
  a valid CDS accompanied by an orphan CDS.

## Verification

| Command | Result | Notes |
|---------|--------|-------|
| `python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests` | passed | 56 passed in 2.72 s; includes pure, external-tool, reverse, and integration groups with no skips. |
| `python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests/test_hipp_functional.py SCRIPT/ANNOTATION_TRANSFER/tests/test_pipeline_integration.py` | passed | 5 passed; FT-01..FT-04 and pipeline/resource assertions. |
| `sha256sum --check SCRIPT/ANNOTATION_TRANSFER/tests/data/hipp_sources.sha256` | passed | All authoritative HIPP inputs retain their pinned checksums. |
| `python -m compileall -q SCRIPT/ANNOTATION_TRANSFER` | passed | New Python files compile. |
| `bash -n bin/create_LRRome.sh bin/genePrediction.sh` | passed | Modified shell entrypoints are syntactically valid. |
| `git diff --check` | passed | No whitespace errors in tracked changes. |
| `python -m pytest -q` | failed outside task scope | Legacy collection segfaults while importing NumPy through Biopython in `SCRIPT/Canonical_gene_model_test.py`; no Python assertion is reached. |
| `python -m pytest -q SCRIPT/CANDIDATE_LOCI/tests` | failed outside task scope | Existing Polars code segfaults in `prepare_blast_for_HSP_export`; unrelated files were not changed. |
| Existing `Extract_sequences_from_genome.py` protein extraction | failed outside task scope | The same compiled-library environment issue segfaults on both model and projected HIPP inputs. Protein acceptance was independently checked in memory with a dependency-free standard-code translator and named Needleman-Wunsch global alignment. |

`pytest 9.1.1` was installed into the active `pipeline` conda environment at the user's
request. Validation used Python 3.11.9, BLAST+ 2.13.0+, and MAFFT 7.526.

## Benchmark Results

Commands used `/usr/bin/time -v` around the production CLI. BLAST, MAFFT, and total values
come from the structured diagnostics; peak RSS and wall time come from `/usr/bin/time`.

| Case | Model nt | Target nt | Decision | BLAST s | MAFFT s | Total s | Wall | Peak RSS KiB | Protein result |
|------|---------:|----------:|----------|--------:|--------:|--------:|------|-------------:|----------------|
| FT-01 | 1403 | 1403 | success, forward | 0.084565 | 0.434182 | 0.753645 | 0:00.83 | 142280 | 249/249 identical; 100%; no internal stop |
| FT-02 | 661 | 1246 | success, forward | 0.087128 | 0.405997 | 0.744935 | 0:00.81 | 132672 | 131/132 identical; 99.2424%; no internal stop |
| FT-03 | 1678 | 2279 | success, forward | 0.090995 | 0.485277 | 0.827047 | 0:00.89 | 180040 | 208/208 identical; 100%; no internal stop |
| FT-04 | 1403 | 4056 | ineligible, no hit | 0.086397 | not invoked | 0.168855 | 0:00.24 | 112264 | zero-byte GFF |

FT-01 projected CDS coordinates, strand, and phases exactly match the normalized model.
FT-04 diagnostics record both MAFFT state and version as `not_invoked`; even a MAFFT
version subprocess is avoided for ineligible pairs.

## Deviations From Specification

- A full mini-dataset `genePrediction.sh` regression and downstream `improve_annot` run
  could not be completed because the available legacy runtime segfaults in its existing
  NumPy/Biopython/Polars dependencies and no self-contained pipeline mini-fixture is
  present. The integration path is instead covered by syntax checks and focused tests
  asserting one conversion before the existing generic improvement/scoring loop.
- The authoritative HIPP GFF stores CDS phase as unknown (`.`). Unknown phases are
  preserved exactly; explicit valid phases `0`, `1`, and `2` are also preserved. No phase
  recalculation was introduced.

## Risks And Follow-Up

- Prebuilt LRRomes created before this change do not contain `REF_LOCI` and
  `REF_LOCI_GFF`; they must be rebuilt (or upgraded with the new preparation helper)
  before `locusAlignment` can run.
- Repair or recreate the compiled scientific-Python environment, then run the full legacy
  suite and a representative end-to-end pipeline mini-dataset through improvement,
  comments, extraction, scoring, merge, and two-round cleanup.
- The pipeline emits per-pair `*_locusAlignment.diagnostics.json` files in addition to the
  declared GFF output. These are intentional benchmark/debug artifacts but are not yet a
  declared Snakemake output.

## Implementation Self-Critique

The standalone and biological acceptance layers are strong, but a maintained,
self-contained pipeline integration fixture would have allowed verification of downstream
improvement/scoring and unchanged legacy method outputs without depending on the broken
host scientific stack. The workspace patch helper also failed intermittently with a
network-namespace sandbox error, requiring the same unified patches to be applied through
Git's patch engine; this increased implementation time without changing the resulting
scope.

## Tentative Git Message

- Add nucleotide locus-alignment annotation transfer
- Implement sparse BLAST/MAFFT CDS projection, LRRome resources, pipeline integration,
  HIPP biological regressions, diagnostics, and benchmarks.

## Notes

- The user-supplied `HIPP/`, V3 specification, and original generated tests were preserved.
- No existing annotation/scoring Python implementation was modified; the core and resource
  logic live in new standalone files.
- BLAST eligibility evaluates one deterministically ranked individual HSP with inclusive,
  unrounded thresholds. It never aggregates HSPs.
