################# QC - MultiQC
################# Anna Jo Muhich
################# February 2023

### Versions
#pip 23.0.1 from /usr/local/lib/python3.9/site-packages/pip (python 3.9)
#Python 2.7.16
#multiqc, version 1.12

### Set up environment
#conda create -n MultiQC
conda activate MultiQC

# Put your fastq files to run into the /input/AtBcfastq directory
cd ~/UCDavis/Klieb_Lab/Projects/Cucurbit/Cuc_RNAseq_Pilot/input/PlantBcfastq
gunzip *.fastq.gz
mkdir ../../qc/fastqc_out
fastqc --threads 3 -o ../../qc/fastqc_out *.fastq 
cd ../../qc
multiqc fastqc_out

conda deactivate

### Proceed to alignment.sh