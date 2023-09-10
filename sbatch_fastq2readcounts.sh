#!/bin/bash
#SBATCH -D /home/ajmuhich/fastq2readcounts
#SBATCH -o /home/ajmuhich/slurm-log/alignment_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/alignment_stderr-%j.txt
#SBATCH -J alignment
#SBATCH -t 90:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load conda
conda init bash

#run trimmomatic
conda activate trimmomatic
bash scripts/trim.sh
conda deactivate

#run qc
conda activate MultiQC
bash scripts/1_QC.sh
conda deactivate

#run alignments
conda activate hisat2
bash scripts/2_alignment.sh
conda deactivate

#generate alignment_summary.csv
bash scripts/2a_alignment_summary1.sh
module load R
R CMD BATCH scripts/2b_alignment_summary2.R

#generate readcounts from .bams
R CMD BATCH scripts/3_readcounts.R

