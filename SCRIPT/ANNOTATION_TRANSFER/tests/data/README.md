# NT annotation-transfer test data

`hipp_cases.tsv` records the four HIPP cases selected for external-tool validation. The BLAST values were verified on 2026-07-15 with the V3 command contract. The paths intentionally point to the repository's authoritative HIPP kit: the pure unit suite uses synthetic in-memory data, while the external functional layer validates these immutable, checksummed source files before deriving temporary locus-relative annotations.

`hipp_sources.sha256` pins every HIPP input used by those cases. Run `sha256sum --check SCRIPT/ANNOTATION_TRANSFER/tests/data/hipp_sources.sha256` from the repository root before external tests. A checksum change requires deliberate review and regeneration of the recorded expectations; it must not be accepted by weakening assertions.

The reference loci were extracted with `SCRIPT/CANDIDATE_LOCI/extract_loci.sh`. Because that script uses `bedtools getfasta -s`, negative-strand genes are reverse-complemented and all reference FASTA records are in model transcription orientation. Genome GFF coordinates must therefore be normalized to 1-based locus-relative coordinates and negative-strand features must be reverse-transformed before transfer.
