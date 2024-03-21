cd /group/kliebengrp/ajmuhich/UCC1
echo "compressing fastq 1/3"
tar -czvf fastq.tar.gz fastq
rm -r fastq

echo "compressing unpaired_fastq"
tar -czvf unpaired_fastq.tar.gz unpaired_fastq
rm -r unpaired_fastq

echo "compressing raw_fastq"
tar -czvf raw_fastq.tar.gz raw_fastq
rm -r raw_fastq