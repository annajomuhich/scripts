################# QC - MultiQC
################# Anna Jo Muhich
################# September 2023

### Versions
#pip 23.0.1 from /usr/local/lib/python3.9/site-packages/pip (python 3.9)
#Python 2.7.16
#multiqc, version 1.12

### Set up environment
module load conda
#conda create -n MultiQC
conda activate MultiQC
#conda install MultiQC

### Generate a list of the samples you want to run and save as file_list.txt
cd fastq2readcounts
#assign your file_list.txt to files variable
readarray -t files < file_list.txt

### Get fastq files
# set up symlink to your fastq files that is hosted in shared lab directory
ln -s /group/kliebengrp/ajmuhich/fastq fastq
cd fastq
# loop to download fastq files using your file list
for file in "${files[@]}"
do
  # Download R1 and R2 for the sample. Change path as needed
  wget -nv "http://slimsdata.genomecenter.ucdavis.edu/Data/04p07y38wc/Unaligned/Project_DKAM_BOS_1/${file}_R1.fastq.gz"
  wget -nv "http://slimsdata.genomecenter.ucdavis.edu/Data/04p07y38wc/Unaligned/Project_DKAM_BOS_1/${file}_R2.fastq.gz"
done
# unzip the files
gunzip *.fastq.gz

### Run QC
# make qc directories
mkdir ~/fastq2readcounts/qc
mkdir ~/fastq2readcounts/qc/fastqc_out
# run fastqc
fastqc --threads 3 -o ~/fastq2readcounts/qc/fastqc_out *.fastq 
cd ~/fastq2readcounts/qc
multiqc fastqc_out

#nav back to fastq2readcounts to proceed
cd ~/fastq2readcounts/
conda deactivate

### Proceed to alignment.sh