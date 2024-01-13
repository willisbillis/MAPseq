#!/bin/bash
PROJECT_PATH=/home/boss_lab/Projects/Scharer_sc/LRA.MAPseq
PROJECT_NAME=LRA_all
OUTPUT_DIR=$PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB
OUTPUT_FILE=$OUTPUT_DIR/cellranger_aggr.log

MEM=128
NCPU=32

################################################################################
mkdir -p $OUTPUT_DIR
AGGR_CSV=$OUTPUT_DIR/aggr.csv

CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" >> $OUTPUT_FILE

cd $OUTPUT_DIR
printf '%s\n' sample_id molecule_h5 | paste -sd ',' >> $AGGR_CSV
for sample_path in $(ls -d $PROJECT_PATH/LRA0*/pipeline/RNA.FB*/*/)
do
if [[ $sample_path/.. == *".BCR"* ]]; then
  sample_name=$(basename $sample_path)
  molecule_file=$(realpath $sample_path)/outs/per_sample_outs/$sample_name/count/sample_molecule_info.h5
  printf '%s\n' $sample_name $molecule_file | paste -sd ',' >> $AGGR_CSV
else
  sample_name=$(basename $sample_path)
  molecule_file=$(realpath $sample_path)/outs/molecule_info.h5
  printf '%s\n' $sample_name $molecule_file | paste -sd ',' >> $AGGR_CSV
fi
done

echo "$(date) Running sample aggregation..." >> $OUTPUT_FILE
cellranger aggr --id=${PROJECT_NAME}_aggr \
	--csv $AGGR_CSV  --normalize none \
	--localcores=$NCPU --localmem=$MEM
