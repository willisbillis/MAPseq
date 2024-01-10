#!/bin/bash
PROJECT_PATH=/home/boss_lab/Projects/Scharer_sc/LRA.MAPseq
PROJECT_NAME=LRA_all
REF_PATH=/home/boss_lab/Apps/genomes/cellranger-ATAC2/hg38
OUTPUT_DIR=$PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP/ATAC
OUTPUT_FILE=$OUTPUT_DIR/cellranger_mapping.log

MEM=128
NCPU=32

################################################################################
mkdir -p $OUTPUT_DIR
AGGR_CSV=$OUTPUT_DIR/aggr.csv

CR_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')
echo "$(date) Running Cell Ranger ATAC version $CR_version using binary $(which cellranger-atac)" >> $OUTPUT_FILE

cd $OUTPUT_DIR
printf '%s\n' library_id fragments cells | paste -sd ',' >> $AGGR_CSV
for sample_path in $(ls -d $PROJECT_PATH/LRA0*/pipeline/ATAC.ASAP/ATAC/*/)
do
if [[ $sample_path != *"_aggr"* ]]; then
  sample_name=$(basename $sample_path)
  fragments_file=$(realpath $sample_path)/outs/fragments.tsv.gz
  cells_file=$(realpath $sample_path)/outs/singlecell.csv
  printf '%s\n' $sample_name $fragments_file $cells_file | paste -sd ',' >> $AGGR_CSV
fi
done

echo "$(date) Running sample aggregation..." >> $OUTPUT_FILE
cellranger-atac aggr --id=${PROJECT_NAME}_aggr --csv $AGGR_CSV --reference $REF_PATH --localcores=$NCPU --localmem=$MEM
