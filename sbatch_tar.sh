#!/bin/bash
#SBATCH -D /group/kliebengrp/ajmuhich
#SBATCH -o /home/ajmuhich/slurm-log/tar_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/tar_stderr-%j.txt
#SBATCH -J tar
#SBATCH -t 90:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

bash scripts/tar.sh

