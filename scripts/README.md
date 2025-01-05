# Scripts

> **Note**: set all scripts to be executable with the command `chmod +x {script_name}`

## 0. rename_fastq.sh

This script renames fasts in a way that is compatible with the Snakemake workflow. Essentially if fastqs follow the name format `{sample_name}_R1_L00{lane_num}.fastq`, they are renamed to `{sample_name}_L00{lane_num}_R1.fastq`. 

> **Usage**: `bash rename_fastq.sh /path_to_your_directory`

## 1. filter_trim.py

This Python script parses read 1 fastq files (the read that should contain the transposon sequence) and returns filtered reads which contain transposon sequence, UMIs, and indices of reads which passed the filter.

> **Usage**: `python filter_trim.py --input /path_to_raw_data/filename.fastq --output_dir /output_dir`


## 2. run_bowtie.sh

This bash script runs bowtie to map fastq reads to reference genome. Please prepare the bowtie index for the reference genome prior to running this step

> **Usage**: `bash run_bowtie.sh --input /path_to_fastq/filename.fastq --output_dir /output_dir --reference /path_to_reference_index --threads 4`


## 3. process_sam.py

This Python script parses sam files and reads in UMIs and indices of reads passing filters and returns position counts matrix, with raw and UMI corrected counts. 

> **Usage**: `python process_sam.py \
            --sam_r1  /path_to_sam/filename_R1.sam \
            --sam_r2 /path_to_sam/filename_R2.sam \
            --umi_list /path_to_filtered_reads/UMI_filename.txt \
            --indices_pf /path_to_filtered_reads/PF_filename.index \
            --output_dir /path_to_results`