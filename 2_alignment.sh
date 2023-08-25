################# Alignments - HiSat2
################# Cucurbit RNAseq
################# Anna Jo Muhich
################# August 2023

###Versions
#hisat2 version 2.2.1

###Set up environment
module load conda
#conda create -n hisat2
conda activate hisat2
#conda install hisat2

####### First, generate a list of the samples you want to run and save as file_list.txt
cd fastq2readcounts
#assign your file_list.txt to files variable
readarray -t files < file_list.txt

#set up needed directories
mkdir fastq
mkdir bams
mkdir readcounts

####### LOOP to align each sample 2 fastqs to host and Bcin genome, convert to .bam, and remove fastqs
cd fastq
for file in "${files[@]}"
do
  # Download R1 and R2 for the sample. Change path as needed
  wget -nv "http://slimsdata.genomecenter.ucdavis.edu/Data/04p07y38wc/Unaligned/Project_DKAM_BOS_1/${file}_1_R1.fastq.gz"
  wget -nv "http://slimsdata.genomecenter.ucdavis.edu/Data/04p07y38wc/Unaligned/Project_DKAM_BOS_1/${file}_1_R2.fastq.gz"
  # unzip Data
  echo ' '
  echo 'unzipping' $file '...'
  echo ' '
  gunzip *.gz
  #navigate to bams output
  cd ../bams
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
      -x  ../../reference/ChineseLong_DNA_index/ChineseLong_DNA  \
      -1 ../../fastq/${file}_1_R1.fastq \
      -2 ../../fastq/${file}_1_R2.fastq \
      --score-min L,0,-0.6 \
      -S ${file}_Host.sam \
      --un-conc ../../fastq/${file}_unmapped.fastq
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
    -x ../../reference/Bcin_toplevelDNA_index/Bcin_toplevelDNA \
    -1 ../../fastq/${file}_unmapped.1.fastq \
    -2 ../../fastq/${file}_unmapped.2.fastq \
    --score-min L,0,-0.85 \
    -S ${file}_Bcin.sam
  #convert Bcin sam to bam
  echo ' '
  echo 'converting .sam to .bam...'
  echo ' '
  samtools view -b ${file}_Bcin.sam > ${file}_Bcin.bam
  rm ${file}_Bcin.sam
  #navigate back to fastq location
  cd ../../fastq
  #clean out fastqs
  rm ${file}*
done
#navigate back to main
cd ../
conda deactivate
########## Proceed to readcounts.R