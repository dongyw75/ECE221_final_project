#!/bin/bash
#PBS -l walltime=6:00:00,nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -j oe
#PBS -A ged
#PBS -M ljcohen@msu.edu
#PBS -m ae


READS=/mnt/home/ljcohen/flies_space/sra/
REFERENCE=/mnt/home/ljcohen/flies_space/ref/
SALMON=/mnt/home/ljcohen/flies_space/salmon/


for fn in $READS/*_1.fastq
do

  # get just the file (without the path)
  base=`basename $fn`;
  echo $base
  # the read filename, without the _R1_001.fastq.gz suffix
  rf=${base%_1.fastq};
  echo $rf
  # Do whatever we want with it
  salmon quant -i $REFERENCE/salmon_index -p 4 -l IU -1 $READS/${rf}_1.fastq -2 $READS/${rf}_2.fastq -o $SALMON/${rf}.quant

done
