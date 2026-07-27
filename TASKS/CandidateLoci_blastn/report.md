# Task Completion Report: BLASTN-guided candidate-locus discovery and expansion

**COMPLETED**: 2026-07-27
**SPECIFICATION**: [V2.md](V2.md)

---

## Summary

Implemented the opt-in BLASTN evidence path for expanding accepted protein-derived
candidate loci and discovering nucleotide-only candidates. The implementation adds
strict one-pass BLASTN parsing, deterministic conflict handling, validated reference
protein lengths, source-aware GFF output, CLI staging/rollback, diagnostics, and
focused regression coverage while preserving the no-BLASTN path. A follow-up
regression now jointly locks the `blastp_default.tsv` recovery diagnostic, origin,
exact interval, and gene+mRNA-only GFF structure.

## Files Created

- `SCRIPT/CANDIDATE_LOCI/blastn_utils.py` - immutable BLASTN configuration/hit
  records, strict parser/filter preparation, ranking, deduplication, and pure
  expansion/discovery matching helpers.
- `SCRIPT/CANDIDATE_LOCI/tests/test_blastn_utils.py` - focused parser, metadata,
  interval, integration, scoring, serialization, CLI, and failure-safety tests.
- `TASKS/CandidateLoci_blastn/report.md` - this implementation report.

## Files Modified

- `SCRIPT/CANDIDATE_LOCI/candidate_loci.py` - keyword-only BLASTN API integration,
  nucleotide-only candidate construction, source-aware output, evidence expansion,
  fixed-bound discovery, neighbor reconciliation, filtering, and diagnostics.
- `SCRIPT/CANDIDATE_LOCI/gff_utils.py` - optional strict single-mRNA/CDS validation
  and inclusive CDS-derived amino-acid protein length metadata.
- `SCRIPT/candidate_loci_VR.py` - `--blastn_table`, input preflight, pass-through,
  and staged paired-output replacement with rollback and unique cleanup.

The two finalized HIPP tests in `test_candidate_loci.py` were not modified.

## Verification

| Command | Result | Notes |
|---------|--------|-------|
| `python -m pytest -q SCRIPT/CANDIDATE_LOCI/tests/test_candidate_loci.py -k 'HIPP_blastn'` (baseline) | failed as expected | Exactly two `TypeError` failures: unsupported `blastn_file`. |
| `bash bin/run_tests.sh` (baseline) | expected partial failure | Annotation: 56 passed. Candidate loci: 25 passed, 1 skipped, exactly the 2 expected HIPP failures. |
| `python -m py_compile SCRIPT/CANDIDATE_LOCI/blastn_utils.py SCRIPT/CANDIDATE_LOCI/gff_utils.py SCRIPT/CANDIDATE_LOCI/candidate_loci.py SCRIPT/candidate_loci_VR.py` | passed | New and integrated modules compile. |
| `python -m pytest -q SCRIPT/CANDIDATE_LOCI/tests/test_candidate_loci.py -k 'HIPP_blastn'` | passed | Both immutable exact-bound regressions pass. |
| `python -m pytest -q SCRIPT/CANDIDATE_LOCI/tests/test_blastn_utils.py` | passed | 37 focused tests, including the coherent default-BLASTP recovery, CLI, and output-failure cases. |
| `bash bin/run_tests.sh` | passed | Annotation: 56 passed. Candidate loci: 64 passed, 1 pre-existing skip. |
| CLI with `blastp_softmax.tsv` and `--blastn_table` | passed | 7 GFF genes and 7 associations; all 7 protein-derived. |
| CLI with `blastp_default.tsv` and `--blastn_table` | passed | 7 GFF genes and 7 associations; 6 protein-derived plus 1 nucleotide-only. |
| `test_HIPP_default_blastp_recovers_one_true_blastn_only_locus` | passed | BLASTP-only lacks `OsjNip_Chr04_19836437`; BLASTN reports `expanded=1, new=1` and returns exactly one `candidateLoci_blastn` gene+mRNA pair at `19835501-19836437`. |
| CLI with `blastp_default.tsv` and no BLASTN | passed | Legacy 6 loci and legacy `candidateLoci` source preserved. |
| `git diff --check` | passed | No whitespace errors. |

The final candidate-locus run retains only the repository's pre-existing explicit
skip. The reported Polars `dtypes` deprecation warnings predate this feature.

## Benchmark Results

No production-representative wheat BLASTN table was present, so the authorized
deterministic synthetic fallback was used. Both 14-column tables contain the same
100 qualifying records (34/33/33 across three chromosomes); the additional rows
are valid but non-qualifying.

| Input | Bytes | Rows | Qualifying/retained | Parse wall time | Peak RSS |
|-------|------:|-----:|--------------------:|----------------:|---------:|
| Synthetic small | 6,289,298 | 100,000 | 100/100 | 0.21 s | 14,084 KB |
| Synthetic 10x | 64,889,298 | 1,000,000 | 100/100 | 1.66 s | 15,740 KB |

The 10x-row scan took about 7.9x as long while peak RSS increased by only
1,656 KB (about 11.8%), supporting the design assumption that raw rows are not
retained. The fixture end-to-end timing over 30 invocations averaged 0.01770 s
without BLASTN and 0.02837 s with BLASTN parsing/integration.

Fixture diagnostics were:

- `blastp_default.tsv` (the required discovery case): 24 raw rows, 7 qualifying
  hits, 1 expanded protein locus, **1 new nucleotide-only locus**, 0
  overlap/filter discards, and 0 clipped expansions;
- `blastp_softmax.tsv` (the expansion-only case): 24 raw rows, 7 qualifying
  hits, 3 expanded protein loci, 0 new loci, 0 overlap/filter discards, and 0
  clipped expansions.

## Deviations From Specification

None. The benchmark used the explicitly permitted synthetic fallback because no
production wheat table was available in the repository.

## Risks And Follow-Up

- Repeat the benchmark with a production wheat BLASTN table when one is available;
  the synthetic evidence validates scaling behavior but not the real qualifying-hit
  fraction.
- Paired CLI output replacement is staged and rolls back caught failures. As with
  ordinary filesystem replacement, process termination or storage failure during
  the final replacement sequence is not a true multi-file transaction.
- Existing GFF parsing emits a Polars deprecation warning for `dtypes`; this is
  unrelated and can be cleaned up separately.

## Implementation Self-Critique

The initial focused coverage asserted BLASTN-only scoring/source and CLI source
counts separately, but it should have jointly asserted the default fixture diagnostic,
origin, exact bounds, and two-feature serialization from the start. The added coherent
regression closes that observability gap. Building all strict GFF cases before
integration would also have shortened the feedback loop around error classification. A reusable CLI paired-output
transaction helper would also make injected rollback-failure testing more direct;
the current subprocess coverage verifies preflight and staging failures but not
every possible `os.replace` failure position. Finally, real wheat data would have
made the memory conclusion stronger than the required synthetic fallback.

## Tentative Git Message

- `Add BLASTN-guided candidate locus discovery`
- Parse qualifying BLASTN evidence once, expand or discover deterministic loci,
  validate GFF protein lengths, add source-aware output/CLI safety, and cover the
  behavior with focused and full-suite tests, including the exact default-BLASTP
  missing-locus recovery contract.

## Notes

- Production input format is documented in `blastn_utils.py` as:
  `blastn ... -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore positive"`.
- The nucleotide-only `OsjNip_Chr04_19836437` regression obtains
  `prot_len == 267`, `pc_similarity == 1.0`, and `score == 267.0`. It also
  requires `new=1`, source `candidateLoci_blastn`, and exactly `gene` plus `mRNA`
  features with no CDS.
- The integration counter is incremented at candidate acceptance and checked against
  the count of returned BLASTN-origin candidates before diagnostics are emitted.
