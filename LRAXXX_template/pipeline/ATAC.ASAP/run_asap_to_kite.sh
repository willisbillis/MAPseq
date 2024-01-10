#!/bin/bash
fastq_path=/home/boss_lab/Projects/Scharer_sc/LRA.MAPseq/LRA001.Oct2023/data/LRA001_ATAC/outs/fastq_path
NCPU=48

python asap_to_kite_v2.py -f $fastq_path -s LRA001_ASAP -o LRA001_ASAP -j TotalSeqB -c $NCPU
