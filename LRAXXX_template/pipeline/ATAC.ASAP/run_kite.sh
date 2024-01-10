#!/bin/bash
#
# run_kite.sh - written by MEW (https://github.com/willisbillis) Dec 2023
# This script runs the kite pipeline formatted for the LRA MAPseq project 
# started in 2023.

# NOTICE: At this point, the cellranger mkfastq-atac to generate fastqs and
#	    asap_to_kite.py generation of kite formatted fastqs should have been run.
#	    The ASAP fastqs are in the PROJECT_DIR/pipeline/ATAC.ASAP/ASAP directory
#	    under their own separate directories named after each sample in the run.


# This is the highest level directory for the LRA run
PROJECT_DIR=/home/boss_lab/Projects/Scharer_sc/LRA.MAPseq/LRA001.Oct2023
# This is the list of space separated sample names of the LRA samples, which should look something like LRAXXX_ASAP_XXX
SAMPLE_LIST=("LRA001_ASAP")
# The number of cores used when running tools
NCPU=32

for sample in ${SAMPLE_LIST[@]}; do
  echo $sample
  WD=$PROJECT_DIR/pipeline/ATAC.ASAP/ASAP/$sample
  barcodes_csv=$WD/FeatureBarcodes.csv
  # cut just the HTO names and sequences from the full table
  cut -d, -f2,5 < $PROJECT_DIR/pipeline/ATAC.ASAP/ASAP/HTOB_master_feature_ref.csv > $barcodes_csv
  # generate the mismatch FASTA and t2g files (for following commands, see tutorial at https://github.com/pachterlab/kite)
  cd $WD
  python $PROJECT_DIR/pipeline/ATAC.ASAP/tools/kite/featuremap/featuremap.py $barcodes_csv --header
  # build kallisto index with mismatch fasta and a k-mer length equal to the length of the Feature Barcodes (of the HTO)
  kallisto index -i $WD/FeaturesMismatch.idx -k 15 $WD/FeaturesMismatch.fa
  # pseudoalign the reads
  kallisto bus -i $WD/FeaturesMismatch.idx -o $WD -x 10xv3 -t $NCPU $WD/*fastq.gz
  # run bustools (note we are NOT running the whitelist filtering command from the tutorial,
  # we are trusting the barcodes we have are good 10x barcodes. Removed because this filtered too many cells in the past)
  bustools sort -t $NCPU -o $WD/output_sorted.bus $WD/output.bus
  mkdir $WD/featurecounts/
  # generate the counts matrix using bustools
  bustools count -o $WD/featurecounts/featurecounts --genecounts -g $WD/FeaturesMismatch.t2g -e $WD/matrix.ec -t $WD/transcripts.txt $WD/output_sorted.bus
done
