# FARMED_onsite\

## Installation:\
1) Clone the repository: `git clone https://git.sciensano.be/brbloemen/FARMED_onsite`
2) Change directory into repository: `cd FARMED_onsite`
3) Create conda environment: `conda env create -f FARMED_onsite.yml`
4) Download required databases: refer to config.yaml file for instructions

## Usage:
1) Activate conda environment: `conda activate FARMED_onsite`
2) Download sequencing data, which you can find here: [https://www.ncbi.nlm.nih.gov/sra/PRJNA1011201](https://www.ncbi.nlm.nih.gov/sra/PRJNA1011201), move all the fastq.gz files to ./data
3) execute workflow: `snakemake --cores {threads} --latency-wait 30`
