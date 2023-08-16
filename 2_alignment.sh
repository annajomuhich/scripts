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
cd alignment
#assign your file_list.txt to files variable
readarray -t files < file_list.txt

####### LOOP UNDER CONSTRUCTION
cd fastq
for file in "${files[@]}"
do
  # Download R1 and R2 for the sample. Change path as needed
  wget "slimsdata.genomecenter.ucdavis.edu/Data/gyfdugxst/Unaligned/Project_DKAA_GSLOHBRASS/${file}_L003_R1_001.fastq.gz"
  wget "slimsdata.genomecenter.ucdavis.edu/Data/gyfdugxst/Unaligned/Project_DKAA_GSLOHBRASS/${file}_L003_R2_001.fastq.gz"
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
  echo 'aligning' $file 'to Chinese Long DNA genome...'
  echo ' '
  hisat2  \
     -p 8 \
      -t \
      --summary-file 'ChineseLong_DNA_log.txt' \
      -I 5 \
      -5 10 \
      -3 5 \
      -x  ../../reference/ChineseLong_DNA_index/ChineseLong_DNA  \
      -1 ../../fastq/${file}_L003_R1_001.fastq \
      -2 ../../fastq/${file}_L003_R2_001.fastq \
      --score-min L,0,-0.6 \
      -S ${file}_ChineseLong_DNA.sam \
      --un-conc ../../fastq/${file}_ChineseLong_DNA_unmapped.fastq
  #convert Host sam to bam
  echo ' '
  echo 'converting .sam to .bam...'
  echo ' '
  samtools view -b ${file}_ChineseLong_DNA.sam > ${file}_ChineseLong_DNA.bam
  rm ${file}_ChineseLong_DNA.sam
  #perform Bcin alignment
  echo ' '
  echo 'aligning ChineseLong DNA unmapped reads from' $file 'to Bc genome...'
  echo ' '
  hisat2 \
    -p 8 \
    -t \
    --summary-file 'Bcin_ChineseLong_DNA_log.txt' \
    -I 5 \
    -5 10 \
    -3 5 \
    -x ../../reference/Bcin_toplevelDNA_index/Bcin_toplevelDNA \
    -1 ../../fastq/${file}_ChineseLong_DNA_unmapped.1.fastq \
    -2 ../../fastq/${file}_ChineseLong_DNA_unmapped.2.fastq \
    --score-min L,0,-0.85 \
    -S ${file}_ChineseLong_DNA_Bcin.sam
  #convert Bcin sam to bam
  echo ' '
  echo 'converting .sam to .bam...'
  echo ' '
  samtools view -b ${file}_ChineseLong_DNA_Bcin.sam > ${file}_ChineseLong_DNA_Bcin.bam
  rm ${file}_ChineseLong_DNA_Bcin.sam
  #navigate back to fastq location
  cd ../../fastq
done

conda deactivate
########## Proceed to readcounts.R