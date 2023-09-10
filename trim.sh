################# Sequence Trimming - Trimmomatic
################# Anna Jo Muhich
################# September 2023

### Set up environment
#conda create -n trimmomatic
conda activate trimmomatic
#conda install -c bioconda trimmomatic

### Generate a list of the samples you want to run and save as file_list.txt
cd fastq2readcounts
#assign your file_list.txt to files variable
readarray -t files < file_list.txt

### Get fastq files
# set up symlink to your fastq files that is hosted in shared lab directory
#mkdir /group/kliebengrp/ajmuhich/raw_fastq
#ln -s /group/kliebengrp/ajmuhich/raw_fastq raw_fastq
#cd fastq
# loop to download fastq files using your file list
#for file in "${files[@]}"
#do
  # Download R1 and R2 for the sample. Change path as needed
#  wget -nv "http://slimsdata.genomecenter.ucdavis.edu/Data/04p07y38wc/Unaligned/Project_DKAM_BOS_1/${file}_R1.fastq.gz"
#  wget -nv "http://slimsdata.genomecenter.ucdavis.edu/Data/04p07y38wc/Unaligned/Project_DKAM_BOS_1/${file}_R2.fastq.gz"
#done
# unzip the files
#gunzip *.fastq.gz

### run trimmomatic
for file in "${files[@]}"
do
  echo ' '
  echo 'Trimming' $file '...'
  echo ' '
  trimmomatic PE -threads 8 ${file}_R1.fastq ${file}_R2.fastq \
  ${file}_R1_trimmed_paired.fastq ${file}_R1_trimmed_unpaired.fastq \
  ${file}_R2_trimmed_paired.fastq ${file}_R2_trimmed_unpaired.fastq \
  ILLUMINACLIP:~/fastq2readcounts/reference/adapters.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
done

### reorganize fastqs
#move trimmed paired fastq to a new directory for downstream work
mkdir /group/kliebengrp/ajmuhich/fastq
ln -s /group/kliebengrp/ajmuhich/fastq fastq
mv *trimmed_paired.fastq fastq
#move trimmed unpaired fastq to a new directory for archival
mkdir /group/kliebengrp/ajmuhich/unpaired_fastq
ln -s /group/kliebengrp/ajmuhich/unpaired_fastq unpaired_fastq
mv *trimmed_unpaired.fastq unpaired_fastq

### compress unused fastqs to conserve space
cd ~/fastq2readcounts
tar -czvf raw_fastq.tar.gz raw_fastq
tar -czvf unpaired_fastq.tar.gz unpaired_fastq

conda deactivate