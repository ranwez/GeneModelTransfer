#!/bin/bash
# Check that LRRtransfer input paths are accessible and non-empty.

set -euo pipefail


usage() {
  cat <<EOF
Usage:
  $0 \
    --target-genome FILE \
    --ref-gff FILE \
    [--ref-genome FILE] \
    [--ref-locus-info FILE] \
    [--lrrome DIR] \
    [--tblastn-results FILE] \
    [--blastn-results FILE] \
    --output FILE
EOF
}


check_file() {
  local label="$1"
  local path="$2"

  if [[ ! -f "$path" ]]; then
    echo "ERROR: ${label} is not a file: ${path}" >&2
    return 1
  fi

  if [[ ! -r "$path" ]]; then
    echo "ERROR: ${label} is not readable: ${path}" >&2
    return 1
  fi

  if [[ ! -s "$path" ]]; then
    echo "ERROR: ${label} is empty: ${path}" >&2
    return 1
  fi
}


check_directory() {
  local label="$1"
  local path="$2"

  if [[ ! -d "$path" ]]; then
    echo "ERROR: ${label} is not a directory: ${path}" >&2
    return 1
  fi

  if [[ ! -r "$path" || ! -x "$path" ]]; then
    echo "ERROR: ${label} is not accessible: ${path}" >&2
    return 1
  fi
}


TARGET_GENOME=""
REF_GENOME=""
REF_GFF=""
REF_LOCUS_INFO=""
LRROME=""
TBLASTN_RESULTS=""
BLASTN_RESULTS=""
OUTPUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target-genome)
      TARGET_GENOME="$2"
      shift 2
      ;;
    --ref-genome)
      REF_GENOME="$2"
      shift 2
      ;;
    --ref-gff)
      REF_GFF="$2"
      shift 2
      ;;
    --ref-locus-info)
      REF_LOCUS_INFO="$2"
      shift 2
      ;;
    --lrrome)
      LRROME="$2"
      shift 2
      ;;
    --tblastn-results)
      TBLASTN_RESULTS="$2"
      shift 2
      ;;
    --blastn-results)
      BLASTN_RESULTS="$2"
      shift 2
      ;;
    --output)
      OUTPUT="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$TARGET_GENOME" || -z "$REF_GFF" || -z "$OUTPUT" ]]; then
  usage >&2
  exit 2
fi

check_file "target genome" "$TARGET_GENOME"
check_file "reference GFF" "$REF_GFF"

if [[ -n "$REF_GENOME" ]]; then
  check_file "reference genome" "$REF_GENOME"
fi

if [[ -n "$REF_LOCUS_INFO" ]]; then
  check_file "reference locus information" "$REF_LOCUS_INFO"
fi

if [[ -n "$TBLASTN_RESULTS" ]]; then
  check_file "TBLASTN results" "$TBLASTN_RESULTS"
fi

if [[ -n "$BLASTN_RESULTS" ]]; then
  check_file "BLASTN results" "$BLASTN_RESULTS"
fi

if [[ -n "$LRROME" ]]; then
  check_directory "LRRome" "$LRROME"
fi

{
  echo "LRRtransfer input files are accessible."
  echo "target_genome: $TARGET_GENOME"
  echo "ref_gff: $REF_GFF"
  echo "ref_genome: ${REF_GENOME:-not provided}"
  echo "ref_locus_info: ${REF_LOCUS_INFO:-not provided}"
  echo "lrrome: ${LRROME:-not provided}"
  echo "tblastn_results: ${TBLASTN_RESULTS:-not provided}"
  echo "blastn_results: ${BLASTN_RESULTS:-not provided}"
} > "$OUTPUT"
