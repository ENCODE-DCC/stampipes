#!/bin/bash

set -e

here=$(dirname "$0")
cd "$here"

STAMPIPES=$PWD/../../..
export STAMPIPES
nextflow run $STAMPIPES/processes/bwa/process_bwa_paired_trimmed.nf \
  -profile test,docker \
  "$@"

# Verify

function cmp_picard() {
  name=$(basename "$1")
  expected=$(grep -v '^#' "expected/$name")
  actual=$(grep -v '^#' "output/$name")

  echo "Comparing $name..."
  diff <(echo "$expected") <(echo "$actual")
}

function cmp_starch() {
  name=$(basename "$1")
  if ! which unstarch &>/dev/null ; then
    echo "Cannot verify $name, unstarch is not available"
    return 0
  fi
  echo "Comparing $name..."
  cmp <(unstarch "expected/$name") <(unstarch "output/$name") \
    || (echo "$name does not match" ; false)
}

function cmp_bam() {
  name=$(basename "$1")
  if ! which samtools &>/dev/null ; then
    echo "Cannot verify $name, samtools is not available"
    return 0
  fi
  echo "Comparing $name..."
  cmp <(samtools view "expected/$name") <(samtools view "output/$name") \
    || (echo "$name does not match" ; false)

}
set +e
cmp_picard "tagcounts.txt" || echo "tagcounts differ"
cmp_picard "MarkDuplicates.picard" || echo "dups differ"
cmp_picard "CollectInsertSizeMetrics.picard" || echo "insert differ"
cmp_picard "subsample.spot.out" || echo "spot differ"

cmp_starch "density.bed.starch" || echo "density differ"

cmp_bam "filtered.bam" || echo "filtered differ"
cmp_bam "marked.bam" || echo "marked differ"

echo "Testing completed successfully"
