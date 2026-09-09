# Task: Nucleotide locus-alignment annotation transfer

**STATUS**: Draft V3 - reviewer answers integrated; unit-test contract created; implementation awaiting approval  
**CREATED**: 2026-07-15  
**PREVIOUS SPECIFICATION**: [`planning_NT_transfert_V2.md`](planning_NT_transfert_V2.md)  
**TEST SUITE**: [`../../SCRIPT/ANNOTATION_TRANSFER/tests/`](../../SCRIPT/ANNOTATION_TRANSFER/tests/)  
**IMPLEMENTATION REPORT**: `TASKS/NT_annotation_transfert/report.md`, to be created from `TASKS/_TEMPLATE/report.md`

## 1. Authorization and scope

The reviewer authorized creation of the pre-implementation unit tests and supplied the authoritative model annotation `HIPP/HIPP_NIP_KIT/OsjNip_HIPPHPP_202606v1.gff`. This V3 records every reviewer decision from V2, corrects the locus-orientation contract based on the supplied extraction command, and defines how the new tests drive and validate implementation.

This phase creates tests and this V3 only. It does not implement `locusAlignment`, modify `genePrediction.sh`, or build LRRome resources. Those changes start after explicit approval of V3.

The goal remains a standalone Python method that gates a model/target locus pair with nucleotide BLAST, aligns eligible pairs with MAFFT, projects one transcript's CDS boundaries, emits a raw minimal GFF, and then participates in the existing prediction-method comparison.

## 2. Decisions incorporated from V2

| Decision | V3 contract |
|---|---|
| Implementation | Python core and direct CLI at `SCRIPT/ANNOTATION_TRANSFER/locus_alignment_transfer.py`; method name `locusAlignment`. |
| BLAST eligibility | Select one best individual HSP; coverage is `length / qlen`; require coverage `>= 0.85` and `pident / 100 >= 0.95`. |
| BLAST format/ranking | Use tabular `blastn -outfmt 6` with `qseqid sseqid qlen length qstart qend sstart send nident pident bitscore evalue`. This retains the pipeline's familiar leading columns and adds the fields required for deterministic ranking/diagnostics. Rank by bit score descending, e-value ascending, aligned length descending, identity descending, then `(qstart, qend, sstart, send, qseqid, sseqid)` ascending. |
| Boundary on target gap | Search inward in the same CDS and reading frame, at offsets 3, 6, and 9 model nucleotides. A start boundary searches forward; an end boundary searches backward. Use the first mapped position. Reject when none exists. Never snap by 1 or 2 nt and never cross the CDS boundary. |
| Ineligible status | Exit `0`, write structured status `ineligible`, create an intentionally empty output GFF, and do not invoke MAFFT. |
| Unexpected failure | Non-zero CLI exit; status `error`; remove or empty any partial GFF so the pipeline aborts rather than scoring it. |
| Reference annotation | Authoritative source is `HIPP/HIPP_NIP_KIT/OsjNip_HIPPHPP_202606v1.gff`. |
| Reference flanks | No flanks: the reference locus is exactly the gene interval from the authoritative GFF. |
| Reference extraction | Use `SCRIPT/CANDIDATE_LOCI/extract_loci.sh <gff> <genome> REF_LOCI`. LRRome construction must reproduce this deterministically. |
| Draft schema | Minimal gene + CDS structure, deterministic IDs, no mRNA record. Preserve valid input phases through projection unless focused downstream compatibility tests require recalculation; document any change. |
| MAFFT | Invoke `mafft --auto` for production and external tests. Record command/version. Tests of pure coordinate logic do not depend on MAFFT output. |

## 3. Corrected orientation and coordinate contract

V2 described reference FASTA as genomic orientation. The supplied extraction command proves otherwise: `extract_loci.sh` calls `bedtools getfasta -s`, which reverse-complements negative-strand genes. Therefore:

1. `REF_LOCI/<model_id>` contains exactly the gene interval in **model transcription orientation**.
2. Genome GFF features must be converted to 1-based locus-relative coordinates before transfer.
3. For a positive-strand source gene with genomic start `G`, a coordinate `p` becomes `p - G + 1`.
4. For a negative-strand source gene with genomic end `E`, an interval `[start, end]` becomes `[E - end + 1, E - start + 1]`.
5. After normalization, the model gene and CDS features are `+` relative to the model FASTA, including models that were `-` on the reference chromosome.
6. Candidate target FASTA is also extracted with strand awareness, but the method must not trust naming or candidate strand as alignment orientation. The selected nucleotide HSP determines whether the target must be reverse-complemented for MAFFT.
7. Projection first produces coordinates relative to the original target FASTA. If BLAST orientation is reverse, prepared-target coordinate `p` becomes `target_length - p + 1`, and the projected strand becomes `-`.
8. Standalone output is locus-relative. Pipeline integration must convert that draft once with the existing `gff_target_to_genome` logic before passing it to `improve_annot`, which currently expects genome-coordinate drafts and converts them back internally. An integration test must catch accidental double conversion.

Coordinates are always 1-based inclusive and serialized with `start <= end`. CDS records are serialized in genomic/locus coordinate order; transcription order is obtained from strand.

## 4. Standalone Python contract fixed by the tests

The implementation module must expose testable core types/functions equivalent to the names imported by the dedicated test suite:

- immutable value objects `BlastHsp`, `CdsFeature`, FASTA/model annotation records, eligibility decision, and transfer result;
- explicit exceptions `InputValidationError`, `AlignmentError`, and `ProjectionError`;
- `parse_blast_tabular`, `assess_eligibility`;
- `read_single_fasta`, `parse_model_gff`, `normalize_model_gff`;
- `projection_candidate_positions`, `map_requested_model_positions`, `reverse_coordinate`, `project_cds_features`;
- `render_draft_gff`; and
- `run_transfer` with injectable `blast_runner` and `mafft_runner` process boundaries.

`map_requested_model_positions` must return exactly the requested keys and must not construct or retain a full model-position dictionary internally. `projection_candidate_positions` is responsible for bounding the requested set to CDS endpoints and the approved in-frame fallbacks.

Names may change only by updating V3 and tests together before implementation review. Equivalent hidden behavior behind an untestable CLI is not acceptable.

The CLI must accept explicit model FASTA, model GFF, target FASTA, output GFF, diagnostics JSON, minimum coverage, minimum identity, BLAST executable, and MAFFT executable. Defaults are `0.85`, `0.95`, `blastn`, and `mafft`. Output replacement must be atomic on success. Diagnostics must include status/reason, IDs, sequence lengths, thresholds, selected raw HSP, unrounded coverage/identity, relative orientation, commands, versions, tool exit states, boundary snaps, and timing.

## 5. Validation algorithm

### 5.1 Input validation

- Each FASTA contains exactly one non-empty uniquely identified record.
- Accept case-insensitive DNA and IUPAC `ACGTRYSWKMBDHVN`; normalize sequence case without changing length.
- The selected model has exactly one gene, one mRNA, and at least one CDS.
- Every CDS belongs to that mRNA, is within the gene/locus, has a valid consistent strand, and does not overlap another CDS.
- Reject malformed rows, orphan CDS, mixed strands, multiple transcripts, out-of-range coordinates, or empty annotations before invoking external tools.

### 5.2 BLAST gate

Run one model-query/target-subject nucleotide search. Parse all significant HSP rows, rank using Section 2, and evaluate only the selected HSP. Do not aggregate HSPs. Do not round before comparison. Subject coordinates increasing means forward; decreasing means reverse. No hit is ordinary ineligibility.

### 5.3 Alignment and projection

Reverse-complement target only for a reverse selected HSP, then run `mafft --auto` on exactly two sequences. Verify that removing gaps reproduces the prepared inputs case-insensitively.

Before traversing the alignment, `projection_candidate_positions` collects only each CDS start and end plus valid inward same-frame alternatives at 3, 6, and 9 model nucleotides that remain inside that CDS. Traverse alignment columns once with 1-based model and prepared-target counters. Advance a counter only for a non-gap character. `map_requested_model_positions` records a target coordinate, or `None` opposite a target gap, only when the current model coordinate is in the requested set; it must not allocate entries for other model positions. At the same time, validate that ungapped aligned sequences reproduce both prepared inputs.

Project exact CDS endpoints from this sparse result. If an endpoint is `None`, try its previously collected inward candidates in nearest-first order. Record every snap in diagnostics. Reject a missing requested coordinate, zero-length, out-of-bounds, overlapping, duplicated, or inconsistently ordered projection. Alignment traversal remains O(alignment length) time, while additional coordinate-map memory is O(number of CDS) because each CDS contributes at most eight distinct requested positions.

### 5.4 Outputs

Success creates a deterministic gene + CDS GFF and status `success`. Ineligibility creates a zero-byte GFF and status `ineligible`. Validation, process, malformed-output, and unprojectable-boundary failures are errors: the CLI returns non-zero and leaves no partial scoreable draft.

## 6. Created test suite and how to use it

The dedicated suite is `SCRIPT/ANNOTATION_TRANSFER/tests/`:

| File | Purpose | When it must pass |
|---|---|---|
| `test_hipp_source_data.py` | Pins authoritative GFF/FASTA checksums, validates all four case paths, single-record FASTA assumptions, and complete gene/mRNA/CDS records for selected models. | Now and throughout implementation. |
| `test_locus_alignment_transfer.py::TestBlastEligibility` | Threshold inclusivity, no-hit outcome, one-HSP ranking/ties, orientation, and malformed BLAST handling. | First implementation slice, before calling tools. |
| `test_locus_alignment_transfer.py::TestCoordinateProjection` | Sparse coordinate use by projection, gap traversal, one-base/terminal coordinates, inward same-frame snapping, reverse transform, strand/order, and invalid projections. | Before MAFFT subprocess integration. |
| `test_sparse_projection_mapping.py` | Candidate-set construction, CDS-boundary containment, requested-position validation, and proof that unrequested locus positions are not retained. | Before MAFFT subprocess integration. |
| `test_locus_alignment_transfer.py::TestInputAndGffContracts` | FASTA/IUPAC validation, invalid GFF rejection, positive/negative genome-to-locus normalization, and deterministic minimal GFF. | Before functional HIPP cases. |
| `test_transfer_orchestration.py` | Injected process boundaries; success, ineligibility/empty output/no MAFFT, and visible BLAST/MAFFT failures with no partial output. | Before exposing the CLI. |
| `data/hipp_cases.tsv` | Machine-readable expected decisions and observed best-HSP metrics for FT-01..FT-04. | External-tool layer and regression review. |
| `data/hipp_sources.sha256` | Detects unreviewed mutation of HIPP evidence. | Before every external-tool validation. |

### 6.1 Red-green implementation sequence

The source module does not exist yet, so importing the contract tests is intentionally red. This is not a completed validation run. Implement in these gates:

1. Run the source-data tests; they must pass before coding against the HIPP evidence.
2. Implement parsing/value objects until input and BLAST unit groups pass.
3. Implement sparse pure mapping/projection until the coordinate group passes without BLAST or MAFFT installed.
4. Implement GFF rendering and orchestration until the complete pure suite passes.
5. Add real BLAST/MAFFT adapters and the external HIPP functional tests described below.
6. Only then integrate with `genePrediction.sh` and run integration/regression tests.

Recommended commands from the repository root:

```bash
python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests/test_hipp_source_data.py
python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests/test_locus_alignment_transfer.py
python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests/test_sparse_projection_mapping.py
python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests/test_transfer_orchestration.py
python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests
python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests -m external_tools
python -m pytest -q SCRIPT/ANNOTATION_TRANSFER/tests -m integration
```

Tests must not be skipped merely because implementation is missing. External/integration markers may skip only when deliberately running the pure subset; required-tool skips do not satisfy completion criteria.

## 7. HIPP functional expectations, now verified

Local nucleotide BLAST validation on 2026-07-15 established:

| ID | Pair | Best individual HSP | Required functional assertion |
|---|---|---|---|
| FT-01 | `OsjNip_Chr04_02891360` to itself | `length=1403`, `qlen=1403`, `pident=100`, forward | Eligible; exact normalized CDS/phase; identical translated protein. |
| FT-02 | `OsjNip_Chr01_23324465` to `OSJkit_chr01_0023883423` | `length=661`, `qlen=661`, `pident=99.697`, forward | Eligible; valid projection; predicted protein global identity >95%; no internal stop. |
| FT-03 | `OsjNip_Chr04_19814005` to `OSJkit_chr04_0020834228` | `length=1679`, `qlen=1678`, `pident=99.762`, forward | Eligible; validates a model/target that are negative on their source genomes but forward after strand-aware extraction; valid projection and protein identity >95%. |
| FT-04 | `OsjNip_Chr04_02891360` to `OSJkit_chr04_0003866905` | no significant nucleotide HSP | Ineligible; zero-byte GFF; reason `no_significant_hit`; MAFFT not invoked. |

FT-03 is no longer called a reverse-alignment case: observed nucleotide coordinates are forward because both negative-genome loci were reverse-complemented during extraction. Reverse alignment remains mandatory and is covered deterministically by synthetic unit tests. Add a real reverse HIPP case only if one is identified without changing thresholds.

External tests must compute current BLAST decisions and biological criteria; the observed numeric values are regression evidence, not byte-for-byte BLAST-version goldens. A material decision change requires review. Protein identity is identical aligned amino-acid columns divided by model protein length in a named global aligner, strictly `> 0.95`; internal stops fail.

## 8. Pipeline and LRRome validation still required

After standalone tests pass:

1. Extend LRRome creation to distribute exact-gene `REF_LOCI/<model_id>` and normalized locus-relative `REF_LOCI_GFF/<model_id>.gff`, with source chromosome/start/end/strand and checksum provenance.
2. Add `locusAlignment` to `methods` and a minimal invocation block in `bin/genePrediction.sh`.
3. Convert standalone locus-relative output to genome coordinates exactly once before the existing improvement/scoring loop.
4. Verify success produces `${outfile}_locusAlignment.gff`; ineligibility produces its empty artifact without changing existing method results.
5. Verify empty-target handling and `best2rounds` cleanup include the method.
6. Run a mini-dataset regression proving existing per-method outputs are unchanged, except a legitimate new best-method selection.
7. Run standalone output through `improve_annot`, `set_gff_comments`, protein extraction, and scoring without schema/phase loss.

## 9. Completion criteria

- [ ] V3 is explicitly approved before production implementation.
- [ ] HIPP source-data tests pass and checksums are reviewed.
- [ ] All created pure contract tests pass without skips.
- [ ] BLAST uses one deterministically ranked HSP and inclusive unrounded thresholds.
- [ ] Mapping retains only CDS endpoints and valid inward 3/6/9 candidates; additional mapping memory is O(number of CDS), not O(locus length).
- [ ] Same-frame boundary snapping uses only inward offsets 3/6/9 and is diagnostic.
- [ ] FT-01..FT-04 automated functional criteria pass with BLAST+ and MAFFT installed.
- [ ] A true reverse-orientation synthetic case passes; no false claim is made about FT-03 orientation.
- [ ] Ineligible is exit 0 + empty GFF; unexpected errors are non-zero and unscoreable.
- [ ] LRRome resources use exact gene intervals, strand-aware FASTA, normalized GFF, and provenance.
- [ ] Pipeline coordinate conversion occurs exactly once and all integration/regression tests pass.
- [ ] Benchmark records BLAST, MAFFT/not-invoked, total time, memory where available, locus lengths, commands, and versions for FT-01..FT-04.
- [ ] `TASKS/NT_annotation_transfert/report.md` records commands, versions, results, skips, deviations, benchmark evidence, risks, and a tentative Git message.

## 10. Deferred/non-goals

Do not aggregate BLAST HSPs, support multiple transcripts, project non-CDS features, replace downstream annotation polishing, alter existing method scoring, or optimize general locus extraction in this task. Findings about extraction performance belong in a separately approved `TASKS/planning_speedUpLocus_extraction_byLRRome.md` task.
