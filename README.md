# FARMED_onsite

## Installation:
1) Clone the repository: `git clone https://github.com/brambloemen/FARMED_onsite`
2) Change directory into repository: `cd FARMED_onsite`
3) Create conda environment: `conda env create -f FARMED_onsite.yml`

## Usage:
1) Activate conda environment: `conda activate FARMED_onsite`
2) Download sequencing data, which you can find here: [https://www.ncbi.nlm.nih.gov/sra/PRJNA1011201](https://www.ncbi.nlm.nih.gov/sra/PRJNA1011201), move all the fastq.gz files to ./data
3) execute workflow: `snakemake --cores {threads} --latency-wait 30`
