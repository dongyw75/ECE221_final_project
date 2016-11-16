#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1:ppn=16
#PBS -l mem=72gb
#PBS -j oe
#PBS -A ged
#PBS -M ljcohen@msu.edu
#PBS -m ae
  
module load STAR/2.5.1b
READS=/mnt/home/ljcohen/flies_space/sra/
ALIGNMENTS=/mnt/home/ljcohen/flies_space/alignments/
REFERENCE=/mnt/home/ljcohen/flies_space/ref/
for fn in $READS/*_1.fastq
do
  # get just the file (without the path)
  base=`basename $fn`;
  echo $base
  # the read filename, without the _R1_001.fastq.gz suffix
  rf=${base%_1.fastq};
  echo $rf
  STAR --runThreadN 14 --genomeDir $REFERENCE/star_index \
    --readFilesIn $READS/${rf}_1.fastq $READS/${rf}_2.fastq \
    --outFileNamePrefix $ALIGNMENTS/${rf} --outSAMtype BAM Unsorted
done  
