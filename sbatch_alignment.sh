#!/bin/bash
#SBATCH -D /home/ajmuhich/alignment
#SBATCH -o /home/ajmuhich/slurm-log/alignment_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/alignment_stderr-%j.txt
#SBATCH -J alignment
#SBATCH -t 90:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

#run alignments
bash /home/ajmuhich/scripts/2_alignment.sh

#generate alignment_summary.csv
bash /home/ajmuhich/scripts/2a_alignment_summary1.sh
R /home/ajmuhich/scripts/2b_alignment_summary2.R

#generate readcounts from .bams
R /home/ajmuhich/scripts/3_readcounts.R