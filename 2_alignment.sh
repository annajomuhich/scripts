################# Alignments - HiSat2
################# Cucurbit RNAseq
################# Anna Jo Muhich
################# August 2023

###Versions
#hisat2 version 2.2.1

###Set up environment
#conda create -n hisat2
eval "$(conda shell.bash hook)"
conda activate hisat2
#conda install hisat2
#conda install samtools

####### First, generate a list of the samples you want to run and save as file_list.txt
cd ~/fastq2readcounts
#assign your file_list.txt to files variable
readarray -t files < file_list.txt

#set up needed directories
#make symlink to bams on shared storage
mkdir /group/kliebengrp/ajmuhich/bams
ln -s /group/kliebengrp/ajmuhich/bams bams
mkdir readcounts

#navigate to bams output
cd ~/fastq2readcounts/bams

####### LOOP to align each sample 2 fastqs to host and Bcin genome, convert to .bam
for file in "${files[@]}"
do
  #make new directory for the file
  mkdir ${file}
  cd ${file}
  #perform Host alignment
  echo ' '
  echo 'aligning' $file 'to Host genome...'
  echo ' '
  hisat2  \
     -p 8 \
      -t \
      --summary-file 'Host_log.txt' \
      -I 5 \
      -5 10 \
      -3 5 \
      -x  ~/fastq2readcounts/reference/Pv_index/Pv  \
      -1 ~/fastq2readcounts/fastq/${file}_R1_trimmed_paired.fastq\
      -2 ~/fastq2readcounts/fastq/${file}_R2_trimmed_paired.fastq \
      --score-min L,0,-0.85 \
      -S ${file}_Host.sam \
      --un-conc ~/fastq2readcounts/fastq/${file}_unmapped.fastq
  #convert Host sam to bam
  echo ' '
  echo 'converting .sam to .bam...'
  echo ' '
  samtools view -b ${file}_Host.sam > ${file}_Host.bam
  rm ${file}_Host.sam
  #perform Bcin alignment
  echo ' '
  echo 'aligning Host unmapped reads from' $file 'to Bcin genome...'
  echo ' '
  hisat2 \
    -p 8 \
    -t \
    --summary-file 'Bcin_log.txt' \
    -I 5 \
    -5 10 \
    -3 5 \
    -x ~/fastq2readcounts/reference/Bcin_toplevelDNA_index/Bcin_toplevelDNA \
    -1 ~/fastq2readcounts/fastq/${file}_unmapped.1.fastq \
    -2 ~/fastq2readcounts/fastq/${file}_unmapped.2.fastq \
    --score-min L,0,-0.85 \
    -S ${file}_Bcin.sam
  #convert Bcin sam to bam
  echo ' '
  echo 'converting .sam to .bam...'
  echo ' '
  samtools view -b ${file}_Bcin.sam > ${file}_Bcin.bam
  rm ${file}_Bcin.sam
  #navigate back to bams
  cd ~/fastq2readcounts/bams
done

#navigate back to main
cd ~/fastq2readcounts/
conda deactivate
########## Proceed to readcounts.R