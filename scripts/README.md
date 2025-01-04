# Scripts

## 1. filter_trim.py

This Python script parses read 1 fastq files (the read that should contain the transposon sequence) and returns filtered reads which contain transposon sequence, UMIs, and indices of reads which passed the filter.

## 2. run_bowtie.sh

This bash script runs bowtie to map fastq reads to reference genome. Please prepare the bowtie index for the reference genome prior to running this step