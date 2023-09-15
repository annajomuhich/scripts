################# QC - MultiQC
################# Anna Jo Muhich
################# September 2023

### Versions
#pip 23.0.1 from /usr/local/lib/python3.9/site-packages/pip (python 3.9)
#Python 2.7.16
#multiqc, version 1.12

### Set up environment
#conda create -n MultiQC
eval "$(conda shell.bash hook)"
conda activate MultiQC
#conda install MultiQC
#conda install fastqc

### Run QC
cd fastq
# make qc directories
mkdir ~/fastq2readcounts/qc
mkdir ~/fastq2readcounts/qc/fastqc_out
# run fastqc
fastqc --threads 8 -o ~/fastq2readcounts/qc/fastqc_out *.fastq 
cd ~/fastq2readcounts/qc
multiqc fastqc_out

#nav back to fastq2readcounts to proceed
cd ~/fastq2readcounts/

### Run QC on raw fastq
cd raw_fastq
# make qc directories
mkdir ~/fastq2readcounts/qc_raw
mkdir ~/fastq2readcounts/qc_raw/fastqc_out
# run fastqc
fastqc --threads 8 -o ~/fastq2readcounts/qc_raw/fastqc_out *.fastq 
cd ~/fastq2readcounts/qc_raw
multiqc fastqc_out

#nav back to fastq2readcounts to proceed
cd ~/fastq2readcounts/
conda deactivate

### Proceed to alignment.sh