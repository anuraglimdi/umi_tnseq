# Transposon sequencing with unique molecular identifiers

Transposon sequencing (TnSeq) is a powerful method used in bacterial genomics for genetic screens to identify essential genes or genes central to a biological process of interest. Transposons, mobile genetic elements found across the tree of life, are engineered into a plasmid based cloning system. These transposons can jump out of the plasmid into the bacterial genome and disrupt genes. 

## Genetic screening using transposon sequencing
TnSeq is based on the idea that if genes can tolerate transposon insertions, they're not essential for a biological process of interest. By comparing the profile of transposon insertions using high-throughput sequencing before and after selection (let's say growth in a particular condition mimicking biofilms in the lungs), it is possible to identify genes where insertions are absent after but not prior to selection. 

The challenge with TnSeq is that the readout can be biased due to PCR amplification, either due to jackpot events where certain templates are overrepresented by getting lucky during initial rounds of amplification, or due to more consistent sequence related biases, often linked to GC content. These confounders make inferences of gene essentiality less reliable. 

## Correcting for PCR amplification biases
In my PhD work ([Changing fitness effects of mutations through long-term bacterial evolution](https://www.science.org/doi/abs/10.1126/science.add1417)), I developed a method, called UMI-TnSeq, to correct for PCR biases by tagging every read with a unique molecular identifier (see Supplementary Information at the link above). If two reads share the same unique molecular identifier and map to the same insertion site, we can be confident that this read is a PCR duplicate, and not a biological duplicate. The design of the UMI-TnSeq method contains two set of unique identifiers for a read:
- A unique DNA tag of length 10 (called the unique molecular identifer) at the 5' end of the read
- The 3' end cut site during the tagmentation step during library preparation

## How to use `umi_tnseq` codebase
This repository contains scripts to process paired-end Illumina sequencing data generated using the UMI-TnSeq method. If you'd like to run scripts individually, detailed usage instructions are in `scripts/README.md`. Alternatively, you can run the entire analysis pipeline using Snakemake. See instructions below:

- Organize the `data` directory as follows, placing all raw fastq files in the `raw` sub-directory

    ```
    ├── data
    │   ├── filtered
    │   ├── raw
    │   ├── sam
    ```
- Rename the fastqs to be compatible with the Snakemake workflow
    ```
    bash scripts/rename_fastq.sh /data/raw
    ```
- Install snakemake (instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html))
- Create environment to run the snakelike pipeline

    ```conda env create -f envs/environment.yaml```
- Run the pipeline with this command (which you may need to modify slightly depending on whether it is run locally or on SLURM)

    ```snakemake -p --snakefile SNAKEFILE --cores 4 --use-conda --conda-frontend conda```
- Final position counts matrices with raw and UMI corrected counts will be written to the `results` directory 
