# Task: Nucleotide locus-alignment annotation transfer

**STATUS**: Draft V2 - awaiting human validation  
**CREATED**: 2026-07-15  
**SCOPE**: Standalone nucleotide-level annotation transfer, focused tests, and minimal integration in `bin/genePrediction.sh`  
**PRIMARY GOAL**: Add a nucleotide locus-alignment method that transfers a single-mRNA CDS annotation between highly similar loci and participates in the existing prediction-method comparison.

**PREVIOUS SPECIFICATION**: [`planning_NT_transfert.md`](planning_NT_transfert.md)  
**REPORT REQUIRED AT IMPLEMENTATION COMPLETION**: `TASKS/NT_annotation_transfert/report.md`, created from `TASKS/_TEMPLATE/report.md`  
**APPROVAL GATE**: This V2 is a specification only. Do not implement code, create fixtures, modify the pipeline, or add LRRome data until a human explicitly approves this V2 and resolves the blocking open questions.

---

## 1. Goal

The current LRRtransfer strategies can miss short exons or introduce errors where a model locus and a target locus are highly similar at nucleotide level. Implement a new method, provisionally named `locusAlignment`, that:

1. determines whether the complete model and target loci are sufficiently similar with nucleotide BLAST;
2. aligns eligible locus pairs with MAFFT;
3. projects the model CDS boundaries through the nucleotide alignment onto the target locus;
4. emits a raw draft GFF compatible with the existing `improve_annot`, scoring, and best-method selection code in `bin/genePrediction.sh`; and
5. is independently and reproducibly tested before pipeline integration.

The implementation must remain easy to test outside the complete Snakemake workflow. It must preserve coordinate, strand, CDS ordering, and output-status information without silently guessing when input or alignment ambiguity makes a valid transfer impossible.

## 2. Non-Goals

- Do not implement anything during review of this V2.
- Do not replace, remove, or change the scoring of existing `mapping` or Exonerate methods.
- Do not change `SNAKEMAKE/LRR_transfer.smk` unless later evidence proves that the agreed interface cannot be integrated solely through `bin/genePrediction.sh` and LRRome content.
- Do not optimize general locus extraction in this task.
- Do not modify `Extract_sequences_from_genome.py` to consume pre-extracted LRRome loci. Record relevant findings in the future specification `TASKS/planning_speedUpLocus_extraction_byLRRome.md` (create or correct that specification only in a separately approved task).
- Do not support multiple transcripts, overlapping CDS features, trans-splicing, or non-CDS feature projection in the first implementation.
- Do not treat protein similarity alone as permission to run the method; eligibility is defined by nucleotide-locus similarity.
- Do not copy all of `HIPP/HIPP_NIP_KIT` into the test suite. Tests must use a documented minimal subset.
- Do not use the pre-existing files in `HIPP/HIPP_NIP_KIT/_HIPP_ALI` as golden MAFFT output: they are exploratory evidence, while the implemented command and deterministic expected results must be controlled by the tests.

## 3. Current State

- `SNAKEMAKE/LRR_transfer.smk` drives candidate-locus discovery and gene prediction.
- `bin/genePrediction.sh` receives a target/model pair, creates one working directory per method, improves each non-empty draft GFF, scores it, and selects the best result.
- The current method list is:

  ```text
  mapping cdna2genome cdna2genomeExon cds2genome cds2genomeExon prot2genome prot2genomeExon
  ```

- `bin/genePrediction.sh` currently reads `REF_PEP`, `REF_EXONS`, and `REF_cDNA` from LRRome. There is no defined `REF_LOCI`/reference-locus-GFF contract yet.
- Target loci are already supplied to `genePrediction.sh` as individual FASTA files. Existing coordinate conversion helpers move annotations between extracted-target and genome coordinates and account for strand orientation.
- `HIPP/HIPP_NIP_KIT` contains individual reference-locus FASTA files, target-locus FASTA files, exploratory pairwise alignments, a candidate-locus GFF, and protein-to-genome BLAST results.
- `HIPP/HIPP_NIP_KIT` does **not** contain the reference model GFF required to define the CDS boundaries for the named `OsjNip_*` examples. `blast_refProt.tsv` is protein-to-genome evidence and cannot serve as the nucleotide eligibility result or as the model annotation.
- Existing Python unit tests use pytest and keep reduced fixture data below the corresponding script's `tests/data` directory.
- The repository report template is `TASKS/_TEMPLATE/report.md`. The path `TASKS/report.md` mentioned in the request does not exist as of this V2.

## 4. Proposed Behavior

### 4.1 Standalone method contract

Create a standalone command whose core logic is testable as Python functions. The provisional locations are:

- `SCRIPT/ANNOTATION_TRANSFER/locus_alignment_transfer.py` - validation, eligibility evaluation, alignment-coordinate mapping, and draft-GFF generation;
- an optional thin shell wrapper only if it materially simplifies calls from `bin/genePrediction.sh`.

The final command name and CLI spelling are subject to open question OQ-1, but the interface must accept these explicit inputs:

| Input | Required contract |
|---|---|
| Model locus FASTA | One record; sequence in extracted-locus genomic orientation |
| Model locus GFF | GFF3-like, 1-based inclusive locus-relative coordinates; exactly one gene, one mRNA, and one or more non-overlapping CDS features |
| Target locus FASTA | One record; sequence in extracted-locus genomic orientation |
| Output draft GFF | Explicit path; written atomically or left absent/empty on ineligible pairs according to OQ-4 |
| Thresholds | Defaults: model-locus coverage `0.85`, nucleotide identity `0.95`; configurable for unit tests and future tuning |
| Diagnostic output | Machine-readable eligibility/status details including selected BLAST record, coverage, identity, orientation, and failure reason |

The command must return distinct outcomes for: successful transfer, pair below eligibility thresholds, invalid input, alignment failure, and unprojectable CDS boundary. Expected exit-code and empty-output semantics are resolved in OQ-4.

### 4.2 Input validation

Before BLAST or MAFFT, the method must validate:

- both FASTA files contain exactly one non-empty record;
- FASTA identifiers are unambiguous and retained in diagnostics;
- the model GFF uses 1-based inclusive coordinates within the model sequence;
- the model contains exactly one mRNA and at least one CDS;
- every CDS belongs to that mRNA, has `start <= end`, and lies inside the locus;
- CDS features do not overlap;
- all model CDS features have the same valid strand (`+` or `-`);
- CDS order can be derived both genomically and in transcription order;
- input parse errors never yield a partial GFF.

Lowercase bases must be accepted. The behavior for ambiguous IUPAC symbols must be documented and tested; BLAST/MAFFT may handle them, but the projection code must not alter sequence length or coordinates.

### 4.3 Nucleotide BLAST eligibility gate

Run nucleotide BLAST directly between the model locus (query) and target locus (subject), using an output format that includes at least query length, aligned length, query start/end, subject start/end, percent identity, bit score, and e-value.

Provisional V2 rule, matching the original specification literally:

1. rank individual HSPs by descending bit score, then descending aligned length, then descending percent identity, with a documented deterministic final tie-break;
2. select the best individual HSP;
3. compute `coverage = aligned_length / model_locus_length`;
4. declare the pair eligible only when `coverage >= 0.85` **and** `pident / 100 >= 0.95`;
5. infer relative orientation from the selected subject coordinates (`sstart <= send` means same sequence orientation; `sstart > send` means reverse orientation).

Threshold comparisons are inclusive. Do not round before comparison. The raw values, selected HSP, command/version, and decision reason must be present in diagnostic output. Whether several collinear HSPs should be aggregated instead is a key decision in OQ-2.

If the gate fails, MAFFT must not run and no draft annotation may enter the method-scoring loop.

### 4.4 MAFFT alignment and orientation

For an eligible pair:

1. retain the model sequence in its input genomic orientation;
2. reverse-complement the target sequence before alignment when BLAST indicates reverse relative orientation;
3. align exactly the two prepared locus sequences with a documented, deterministic MAFFT invocation;
4. verify that the ungapped aligned sequences reproduce the exact prepared inputs (case-insensitive); and
5. retain the transform needed to convert coordinates from aligned/prepared target orientation back to the original target-locus coordinates.

The projected target gene strand is the model gene strand when the two locus sequences have the same orientation and the opposite strand when the target was reverse-complemented. Output coordinates must always satisfy `start <= end`, irrespective of strand.

MAFFT stderr, version, and exit status must be captured for diagnostics without allowing nondeterministic temporary paths into golden outputs.

### 4.5 CDS boundary projection

Coordinate projection must be implemented independently of process execution so it can be exhaustively unit tested with small in-memory alignments.

Traverse alignment columns from left to right while maintaining 1-based model and prepared-target coordinates. Increment a sequence coordinate only when that sequence has a non-gap character in the current column. Build a mapping from each non-gap model position to its homologous prepared-target position or to an explicit `unmappable` value when the model base aligns to a target gap.

For each CDS, project both inclusive boundary positions, convert projected target coordinates back to original target-locus coordinates when reverse-complemented, normalize `start <= end`, and order output CDS features genomically. The output strand and transcription order must remain correct.

The first implementation must not silently shift an endpoint whose exact alignment column contains a target gap. Recommended behavior is to mark the transfer unprojectable and emit no draft, because choosing the preceding/following base changes exon length and can conceal a frameshift. OQ-3 offers alternatives for human decision.

After projection, reject the transfer if any CDS is zero-length, out of bounds, duplicated, overlapping, reordered inconsistently, or lacks both mapped endpoints.

### 4.6 Raw draft GFF

On success, write the minimal draft-GFF structure consumed by `improve_annot`: one gene record followed by ordered CDS records. The implementation must confirm from `improve_annot`/`set_gff_comments` tests whether an mRNA record is accepted or needed before finalizing the schema.

At minimum, the raw result must have:

- target-locus sequence ID in column 1;
- a stable source label and feature types;
- 1-based inclusive target-locus coordinates;
- projected strand;
- a documented phase policy compatible with downstream correction;
- deterministic IDs/Parent attributes; and
- an eligibility/method diagnostic outside the GFF or in comments that downstream parsing safely ignores.

Golden tests compare normalized GFF semantics (feature type, sequence ID, coordinates, strand, phase, and parentage). A byte-for-byte assertion may additionally be used once output serialization is fixed.

### 4.7 LRRome data contract

Extend LRRome with pre-extracted model-locus data, provisionally:

```text
LRRome/
  REF_LOCI/<model_id>.fasta
  REF_LOCI_GFF/<model_id>.gff
```

Each model ID used by `REF_PEP`, `REF_cDNA`, and `REF_EXONS` must resolve deterministically to one locus FASTA and one locus-relative GFF. Define and validate:

- filename extension and exact ID matching;
- extraction flank policy;
- whether the FASTA remains in genomic orientation;
- conversion of genome GFF coordinates into locus-relative 1-based inclusive coordinates;
- preservation of the original chromosome, offset, strand, and source annotation as provenance metadata;
- behavior when either file is missing or malformed; and
- whether `bin/create_LRRome.sh` builds these resources or they are distributed precomputed (OQ-5).

No full-genome extraction optimization is included in this task. The new files are inputs to the new method only.

### 4.8 Pipeline integration

After standalone tests pass:

1. add `locusAlignment` to the `methods` list in `bin/genePrediction.sh`;
2. create its temporary method directory through the existing loop;
3. invoke the standalone command with the selected model locus/GFF and existing target locus;
4. place a successful raw draft at `locusAlignment/${target}_draft.gff`;
5. leave the method with no scoreable draft for ineligible pairs;
6. allow the existing `improve_annot`, `set_gff_comments`, protein-similarity scoring, per-method output, and best-method selection loops to handle it identically to other methods; and
7. ensure empty-target and second-round logic create/remove the additional per-method files consistently.

The original V1 phrase “the sole modification ... will be to add a new method in the list” is not executable literally: adding the list item creates a directory/output but does not generate its draft. V2 confines pipeline edits to the method list plus the smallest invocation block in `bin/genePrediction.sh`; no broader workflow edit is planned.

### 4.9 Failure and reproducibility behavior

- Ineligibility is an expected method-level outcome, not a whole-pipeline failure.
- Invalid data, missing executables, or BLAST/MAFFT process failures must be distinguishable from ordinary ineligibility and must not be silently converted into it.
- All subprocess arguments must be passed safely without shell interpolation in the Python core.
- Temporary files must use an isolated temporary directory and be cleaned on normal completion; a debug option may retain them.
- Repeated runs with the same inputs and tool versions must produce semantically identical GFF and diagnostic results.
- External tool names/paths must be configurable or resolved through the execution environment consistently with existing pipeline conventions.

## 5. Affected Files

Only this V2 specification is created in the current phase. The following are **planned implementation-phase paths**, subject to approval and OQ-1/OQ-5:

- `TASKS/NT_annotation_transfert/planning_NT_transfert_V2.md` - this traceable V2 specification.
- `SCRIPT/ANNOTATION_TRANSFER/locus_alignment_transfer.py` - proposed standalone implementation.
- `SCRIPT/ANNOTATION_TRANSFER/tests/test_locus_alignment_transfer.py` - proposed unit/functional tests.
- `SCRIPT/ANNOTATION_TRANSFER/tests/data/hipp_nip_kit_subset/` - proposed minimal, provenance-documented HIPP-derived fixtures and expected outputs/criteria.
- `bin/genePrediction.sh` - add the method and minimal invocation only.
- `bin/create_LRRome.sh` and/or LRRome documentation/data layout - establish the approved model-locus FASTA/GFF contract.
- dependency/environment documentation - declare supported BLAST+, MAFFT, Python, and pytest requirements if not already declared.
- `TASKS/NT_annotation_transfert/report.md` - mandatory full completion report created from `TASKS/_TEMPLATE/report.md`.
- `TASKS/planning_speedUpLocus_extraction_byLRRome.md` - future, separate task specification; only relevant discoveries are to be noted for later work.

## 6. Design Constraints And Invariants

- GFF and locus coordinates are 1-based and inclusive; alignment counters must state their initialization and update order explicitly in code/tests.
- FASTA sequence length excludes headers and whitespace and must equal the coordinate domain used by its locus-relative GFF.
- Eligibility means both thresholds pass; `85%` coverage and `95%` identity are inclusive boundaries.
- The eligibility denominator is model-locus length, not target length, CDS length, or aligned ungapped length.
- A failed eligibility gate never invokes MAFFT and never yields a scoreable draft.
- The model annotation has exactly one mRNA with non-overlapping CDS features.
- Input CDS coordinates and strands are never mutated in place.
- Every successful output CDS boundary maps to an actual non-gap target nucleotide.
- Output CDS features remain non-empty, within target bounds, non-overlapping, and correctly ordered.
- Reverse-orientation transfer must transform both coordinates and strand exactly once.
- The projected annotation is raw input to existing improvement logic; the new method must not duplicate splice-site, missing-start, or scoring improvements already performed downstream.
- Existing method results and selection behavior must remain unchanged when the new method is ineligible.
- Test fixtures must be small, committed, license-compatible, and traceable to their source path and extraction/transformation procedure.
- No implementation begins before human approval of this V2 and resolution of blocking choices.

## 7. Risks And Edge Cases

- **Missing source annotation**: named HIPP reference loci lack their model GFF in the kit, so exact CDS projection tests cannot yet be constructed reproducibly.
- **Ambiguous “best BLAST hit”**: a single-HSP rule may reject loci whose similarity is split across indels/repeats; aggregation risks admitting rearranged or duplicated matches.
- **Repeated loci**: the highest-bit-score HSP may cover a repeated segment rather than preserve whole-locus collinearity.
- **Boundary gaps**: an insertion/deletion exactly at a CDS boundary has no unique homologous target coordinate.
- **Large internal indels**: both thresholds may pass even though a CDS contains a disruptive indel; downstream protein validation is therefore required.
- **Reverse strand/off-by-one errors**: extracted-locus and genomic orientation can be confused, especially when converting `p -> target_length - p + 1`.
- **Flank asymmetry**: different reference/target extraction flanks can lower whole-locus coverage despite a nearly identical gene.
- **Multiple FASTA records or duplicate IDs**: must fail validation rather than select one implicitly.
- **Model annotation irregularities**: multi-mRNA genes, overlapping CDS, missing parents, mixed strands, or phase inconsistencies are outside first-version support.
- **Tool-version differences**: MAFFT/BLAST output selection can vary; commands and versions must be reported and golden tests should assert biological semantics rather than exploratory alignment bytes.
- **Method-selection regression**: a structurally plausible but frameshifted transfer could score competitively; translated-protein validation and integration regression tests are mandatory.
- **Dependency availability**: standalone test layers must distinguish pure coordinate tests from external-tool tests so useful unit tests still run when BLAST+/MAFFT are unavailable.
- **Naming**: the repository uses both “transfer” and “transfert”; new code/API names should use consistent English spelling while retaining existing task paths.

## 8. Validation Plan

Selected verification tier: `custom` (pure unit tests, external-tool functional tests, and focused pipeline integration tests)

### 8.1 Test fixture policy

During implementation, create `SCRIPT/ANNOTATION_TRANSFER/tests/data/hipp_nip_kit_subset/` rather than reading mutable exploratory data directly from `HIPP/HIPP_NIP_KIT`. Include a README/manifest recording for every copied or derived file:

- original repository path and source identifier;
- checksum of the source content;
- whether the file is copied verbatim, renamed, trimmed, reverse-complemented, or newly authored;
- extraction command and coordinates for derived content;
- source and license/provenance of model GFF and model protein;
- expected outcome type (`exact`, `criterion`, or `reject`); and
- the reason the case is retained.

Keep only files needed by the cases below. Expected GFF/protein files must be committed when the exact answer is known. When it is not known, commit an explicit machine-readable criteria file and make the test compute the criterion; do not rely on a README-only manual assertion.

### 8.2 Required HIPP-derived functional cases

| ID | Model | Target | Purpose | Expected result |
|---|---|---|---|---|
| FT-01 | `OsjNip_Chr04_02891360` | the same `OsjNip_Chr04_02891360` locus | Exact identity; eligibility and transfer baseline | BLAST gate passes with 100% coverage/identity; raw projected CDS coordinates/strand/phase match the normalized model annotation exactly; translated protein is identical |
| FT-02 | `OsjNip_Chr01_23324465` | `OSJkit_chr01_0023883423` | Similar non-identical, same-orientation pair | Gate passes; output is structurally valid; predicted protein has **>95% global amino-acid identity** to the model protein under the metric defined below |
| FT-03 | `OsjNip_Chr04_19814005` | `OSJkit_chr04_0020834228` | Similar non-identical, reverse-orientation pair | Gate passes; output strand is reversed exactly once; structurally valid; predicted protein has **>95% global amino-acid identity** to model |
| FT-04 | `OsjNip_Chr04_02891360` | `OSJkit_chr04_0003866905` | Divergent negative control | Nucleotide gate fails; MAFFT is not invoked; no scoreable draft GFF is produced; diagnostic reason identifies the failed threshold(s) |

FT-02/FT-03 eligibility is a hypothesis from exploratory data, not proven by `blast_refProt.tsv`. Fixture construction must first run the specified nucleotide gate. If either pair does not meet the fixed nucleotide thresholds, do not weaken the thresholds: document the result in the report and replace/add a demonstrably eligible HIPP pair only after human approval.

For the protein criterion, “>95% identity” means a reproducible end-to-end/global alignment of translated predicted CDS versus the model protein, with identity calculated as identical aligned amino-acid columns divided by the model protein length. Internal stops fail. The implementation report must name the alignment implementation and treatment of terminal stop characters. Exact identity is `> 0.95`, not `>= 0.95`, matching the user-requested acceptance wording.

### 8.3 Required pure unit tests

- exact threshold boundaries (`0.85`, `0.95`) pass; values immediately below either threshold fail;
- deterministic best-HSP ranking and tie handling;
- no MAFFT call after eligibility rejection (mocked process boundary);
- forward and reverse subject orientation inference;
- coordinate traversal with no gaps, model-only gaps, target-only gaps, and gaps inside a CDS;
- one-base and other short exons are preserved;
- first and last locus nucleotide boundaries have no off-by-one error;
- reverse-coordinate transform uses `target_length - position + 1` and preserves `start <= end`;
- plus- and minus-strand multi-CDS ordering;
- exact boundary aligned to a target gap follows the approved OQ-3 policy;
- lowercase and supported ambiguous nucleotide input;
- deterministic normalized GFF serialization;
- rejection of empty/multi-record FASTA, out-of-range GFF, no CDS, multiple mRNAs, overlapping CDS, mixed strands, and orphan CDS;
- subprocess failure and malformed BLAST/MAFFT output produce explicit error outcomes and no partial draft.

### 8.4 Required integration tests

- successful standalone output is accepted by `improve_annot` and `set_gff_comments` without schema loss;
- an eligible case produces `${outfile}_locusAlignment.gff` and participates in best-method comparison;
- an ineligible case produces the approved empty/absent per-method artifact and leaves all existing method scores/results unchanged;
- empty target-locus handling creates the new method's expected empty output;
- `best2rounds` cleanup/re-entry includes the new method without stale files;
- at least one existing mini-dataset regression run confirms pre-existing outputs other than the new method/best selection are unchanged.

Commands expected (final paths may change under OQ-1):

```bash
pytest -q SCRIPT/ANNOTATION_TRANSFER/tests/test_locus_alignment_transfer.py
pytest -q SCRIPT/ANNOTATION_TRANSFER/tests -m external_tools
pytest -q SCRIPT/ANNOTATION_TRANSFER/tests -m integration
```

The implementation report must record exact commands, dependency versions, results, skipped tests and reasons, and paths to retained diagnostic logs. A test that is skipped because a required tool is missing does not satisfy completion criteria.

## 9. Benchmark Plan

Benchmark evidence is required but is not a performance acceptance gate for the initial feature.

Run the standalone method on FT-01 through FT-04 and report, per pair:

- BLAST wall time;
- MAFFT wall time (or “not invoked” for rejection);
- total wall time;
- peak resident memory when the environment provides a reproducible measurement;
- locus lengths and eligibility decision; and
- BLAST+, MAFFT, Python, CPU/thread, and command configuration.

Confirm specifically that FT-04 stops after BLAST. Compare one focused `genePrediction.sh` run before/after integration when feasible and report the incremental time. These measurements characterize cost and catch accidental full-genome inputs; they do not belong to the deferred locus-extraction optimization.

## 10. Open Questions

Questions marked **blocking** require human resolution before implementation.

### OQ-1 (**blocking**) - standalone implementation/API layout

**Recommended: Python core plus direct CLI in `SCRIPT/ANNOTATION_TRANSFER/`.**

- Pros: pure mapping logic is easy to unit test; robust parsing/status data; safe subprocess execution; matches existing pytest practice.
- Cons: introduces a new Python module/CLI surface and requires its dependencies to be documented.

Alternative: shell/AWK orchestration.

- Pros: stylistically close to `genePrediction.sh`; minimal Python packaging work.
- Cons: harder to validate GFF rigorously, test boundary mapping, return structured diagnostics, and handle errors safely.

Decision needed: approve the recommended layout and provisional method name `locusAlignment`, or provide preferred paths/name.

=> use python core with direct CLI in `SCRIPT/ANNOTATION_TRANSFER/locus_alignment_transfer.py` and retain the provisional method name `locusAlignment`.

### OQ-2 (**blocking**) - meaning of “best BLAST hit” and coverage

**Recommended for strict traceability: one best individual HSP, `length / qlen`.**

- Pros: directly matches the original wording; simple, deterministic, conservative.
- Cons: may reject biologically suitable loci when similarity is split into multiple collinear HSPs.

Alternative: aggregate non-overlapping, collinear HSP query coverage and compute length-weighted identity.

- Pros: tolerates long insertions/deletions and split alignments.
- Cons: requires chaining rules and may combine repeats or rearrangements; more edge cases and tests.

Decision needed: retain the provisional single-HSP definition or specify an HSP-chaining policy.

=> use one best individual HSP, `length / qlen`. Use the tabular output from `blastn -outfmt 6` (using exact same output format as elsewhere in the pipeline for readibility); test if there are signifcant hits, sort by bit score or evalue (as done elsewhere) and decide if the best hit meets the criteria using the tabular output for this hit.

### OQ-3 (**blocking**) - CDS boundary aligned to a target gap

**Recommended: reject that transfer as unprojectable.**

- Pros: never invents an exon boundary; exposes a biologically meaningful ambiguity.
- Cons: loses some otherwise useful near-identical transfers.

Alternative A: move inward to the nearest aligned target base within the projected CDS.

- Pros: avoids expanding CDS and may preserve a plausible exon.
- Cons: shortens CDS and can alter frame/protein without splice evidence.

Alternative B: use the current/preceding coordinate implied by the original counter description.

- Pros: closest to the initial algorithm sketch.
- Cons: update order is easy to misinterpret and can select a base outside the homologous boundary.

Decision needed: choose rejection or an exact, testable snapping rule.
=> find alternatives in the same reading frame, up to 3 codons away (in the same CDS):
- for a start boundary, move forward to the next aligned target base within the same reading frame, up to 3 codons away; if none, reject.
- for an end boundary, move backward to the previous aligned target base within the same reading frame; if none, reject.
(the idea is that we allow missing part of the CURRENT CDS up to 3 codons at the start and end, this should not impact the translation of the predicted protein (just removing a few amino acids))
So the move should remain in the same reading frame and same CDS of the model protein in both cases. Double check the logic. But this should remain simple, easy to read (annotation polishing is done elsewhere).
### OQ-4 (**blocking**) - ineligible versus error status contract

**Recommended: exit `0` with structured status `ineligible` and an intentionally empty output for threshold rejection; non-zero exit for invalid input/tool/alignment/projection errors.**

- Pros: expected method rejection does not abort the pipeline; real failures remain visible; current method loop already understands empty drafts.
- Cons: callers must inspect both status and output; an empty file can be mistaken for missing work without diagnostics.

Alternative: distinct non-zero code for ineligibility that `genePrediction.sh` catches explicitly.

- Pros: very explicit control flow.
- Cons: more shell error handling under `set -e`; easy to accidentally abort the full job.

Decision needed: approve the recommended contract and decide whether ineligible output is empty or absent before pipeline normalization.
=> I approve the recommended contract and decide that ineligible output is empty (not absent) before pipeline normalization and non zero is for unexpected errors that should abort the pipeline.

### OQ-5 (**blocking**) - source and construction of reference locus GFF/FASTA

**Recommended: generate both deterministically in `bin/create_LRRome.sh` from the authoritative reference genome/GFF, while distributing them in each completed LRRome.**

- Pros: reproducible provenance; all model IDs remain synchronized; no runtime whole-genome extraction.
- Cons: expands LRRome creation and storage; requires a stable flank policy and migration for existing LRRomes.

Alternative: maintain a separately precomputed/manual LRRome add-on.

- Pros: isolates this task from LRRome builder changes and allows early testing.
- Cons: risks stale/missing files and weaker reproducibility.

Decision/data needed: identify the authoritative reference GFF containing the HIPP `OsjNip_*` model annotations, approve a flank policy, and choose builder-generated versus separately distributed resources.
=> I provided the GFF file HIPP/HIPP_NIP_KIT/OsjNip_HIPPHPP_202606v1.gff of the model annotations.
=> Flanking policy: no flanks for the reference locus, just the exact coordinates of the gene in the GFF. As done with the provided REF_LOCI using :
```bash
SCRIPT/CANDIDATE_LOCI/extract_loci.sh OsjNip_HIPPHPP_202606v1.gff Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fasta REF_LOCI
```
intergenic region may be more variable and interfere with the process. For candidate loci there are flanking regions since the exact position of the gene is not known yet, this is done using a separate script (candidate_loci.sh). But for the reference locus we want to use the exact coordinates of the gene in the GFF, this also simplify the implementation reducing the coordinate shifting.

### OQ-6 - raw draft GFF schema and phase

**Recommended: mirror the minimal gene+CDS draft structure currently accepted from existing methods, use deterministic IDs, and preserve/recalculate phase only according to verified downstream expectations.**

- Pros: smallest integration surface; existing improvement/scoring path remains authoritative.
- Cons: may omit transcript metadata useful in standalone output.

Alternative: emit full gene+mRNA+CDS GFF3.

- Pros: self-contained and standards-oriented.
- Cons: may require changes in shell helpers that currently expect minimal drafts.

Decision may be made from focused compatibility tests during implementation, but any deviation must be documented in the report.

=> accept the recommended minimal gene+CDS draft structure, use deterministic IDs, and preserve/recalculate phase only according to verified downstream expectations. 
### OQ-7 - MAFFT mode and thread policy

**Recommended: fixed deterministic two-sequence parameters and one thread for tests; configurable threads in production.**

- Pros: reproducible golden tests and predictable resource use.
- Cons: single-thread production may be slower unless overridden.

Alternative: `--auto` with environment-dependent threads.

- Pros: simple and potentially faster.
- Cons: algorithm/version/thread choices can change output and complicate reproducibility.

The implementation report must record the final invocation and supported MAFFT version range.
=> use mafft with no parameters (default); default to --auto since depending on sequence length mafft will choose the appropriate method. This should no break the tests in most cases (easy alignment when this strategy is used).

## 11. Implementation Plan

Implementation starts only after V2 approval and blocking decisions.

1. Resolve OQ-1 through OQ-5 and record decisions in an approved V3 or an explicit decision addendum linked from this file.
2. Locate the authoritative model GFF/protein inputs, define reference-locus extraction/flanks, and build a provenance manifest for the four required HIPP cases.
3. Create the minimal dedicated fixture directory with copied/derived FASTA, model GFF/protein, exact expected outputs where known, and machine-readable criteria otherwise.
4. Write failing pure unit tests for input validation, BLAST decision logic, alignment-coordinate mapping, reverse orientation, boundary gaps, and GFF serialization.
5. Implement the Python core and CLI contract without pipeline changes.
6. Add external-tool functional tests for the four required cases; verify the nucleotide eligibility assumptions and resolve fixture substitutions through human review if needed.
7. Validate raw draft compatibility with existing `improve_annot` and scoring helpers; finalize GFF schema/phase policy.
8. Add the approved LRRome reference-locus FASTA/GFF layout and validation/build path.
9. Integrate `locusAlignment` into `bin/genePrediction.sh` with the minimal list, invocation, empty-target, and second-round changes.
10. Run pure, external-tool, integration, and existing regression tests; run the benchmark plan.
11. Update user/developer documentation for inputs, dependencies, statuses, thresholds, diagnostics, and limitations.
12. Create `TASKS/NT_annotation_transfert/report.md` from `TASKS/_TEMPLATE/report.md` and fully populate every section, including exact verification evidence, benchmark results, deviations, risks, self-critique, and tentative Git message.

## 12. Completion Criteria

- [ ] Human approval of this V2 (or its approved successor/addendum) is recorded before implementation begins.
- [ ] OQ-1 through OQ-5 have explicit decisions and the model annotation/protein source is identified.
- [ ] A standalone, documented command implements the approved BLAST gate, MAFFT alignment, boundary projection, orientation transform, and draft-GFF output.
- [ ] Threshold and mapping logic is independently unit tested, including inclusive boundaries and off-by-one/reverse cases.
- [ ] The minimal HIPP-derived fixture subset has checksummed provenance and no unnecessary bulk data.
- [ ] FT-01 exact identity, FT-02 same orientation, FT-03 reverse orientation, and FT-04 rejection all meet their approved automated expectations.
- [ ] Unknown exact outputs use executable criteria; FT-02/FT-03 predicted proteins are globally **>95% identical** to their model proteins with no internal stop.
- [ ] Ineligible pairs skip MAFFT and do not enter method scoring.
- [ ] Invalid inputs and tool/projection failures are visible and cannot masquerade as ordinary ineligibility.
- [ ] LRRome provides validated, deterministic reference-locus FASTA and locus-relative GFF inputs for every applicable model.
- [ ] `locusAlignment` participates in existing improvement, scoring, per-method output, best selection, empty-target handling, and second-round cleanup.
- [ ] Existing method outputs remain unchanged for the focused regression dataset except where the new method legitimately changes the selected best annotation.
- [ ] All required tests pass with no required-tool skips, and exact commands/versions/results are recorded.
- [ ] Benchmark evidence is recorded as specified.
- [ ] Deferred locus-extraction findings are captured for `TASKS/planning_speedUpLocus_extraction_byLRRome.md` without implementing that optimization here.
- [ ] A full `TASKS/NT_annotation_transfert/report.md` is created from `TASKS/_TEMPLATE/report.md` and links to the final approved specification.

## 13. Traceability Matrix

| Requirement ID | V1/request source | V2 specification | Verification/evidence |
|---|---|---|---|
| TR-01 | No coding; plan first | Approval gate; Sections 1-2 | Human-approved spec before implementation commit |
| TR-02 | Obtain model and target loci | Sections 4.1, 4.7 | Input-validation tests; LRRome contract tests |
| TR-03 | BLAST: >=85% model coverage and >=95% identity | Section 4.3; OQ-2 | Threshold unit tests; FT-01 to FT-04 diagnostics |
| TR-04 | Align eligible loci with MAFFT | Section 4.4 | External-tool tests; MAFFT skipped in FT-04 |
| TR-05 | Transfer ordered, non-overlapping single-mRNA CDS through alignment | Sections 4.2, 4.5 | Pure coordinate tests; normalized GFF assertions |
| TR-06 | Handle same and opposite strands | Sections 4.4-4.5 | FT-02, FT-03; reverse-coordinate unit tests |
| TR-07 | Add a new `genePrediction.sh` method | Section 4.8 | Integration tests for output/scoring/best selection |
| TR-08 | Implement standalone before integration | Sections 4.1, 11 | Report chronology and standalone test results |
| TR-09 | Add reference locus/GFF data to LRRome | Section 4.7; OQ-5 | LRRome validation and provenance manifest |
| TR-10 | Exact identical-locus test | FT-01 | Exact normalized CDS and protein assertions |
| TR-11 | Similar same-strand and opposite-strand tests | FT-02, FT-03 | Structural assertions and >95% protein identity |
| TR-12 | Divergent pair rejected by gate | FT-04 | Gate failure, MAFFT mock/diagnostic, no draft |
| TR-13 | Dedicated reduced fixture subset with expected output/criteria | Section 8.1 | Committed manifest, checksums, golden/criteria files |
| TR-14 | Defer faster LRRome locus extraction | Sections 2, 5, 12 | Reported follow-up only; no extraction optimization diff |
| TR-15 | Full implementation report in task folder | Metadata; Sections 11-12 | Completed `TASKS/NT_annotation_transfert/report.md` from actual repository template |

### V2 clarifications relative to the original planning file

- Corrected the completion-report template path from the requested but absent `TASKS/report.md` to the repository's actual `TASKS/_TEMPLATE/report.md`; the completed report itself belongs in this task folder.
- Made explicit that integration needs both a method-list entry and a small invocation block in `bin/genePrediction.sh`.
- Distinguished protein BLAST evidence already present in the HIPP kit from the required new nucleotide-locus eligibility BLAST.
- Added explicit coordinate, orientation, failure, reproducibility, fixture-provenance, validation, and benchmark contracts.
- Preserved uncertain algorithmic choices as blocking questions with recommended alternatives rather than silently deciding them.
