#!/bin/bash
#
# run_cellranger.aggr_rna.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger aggr pipeline formatted for the RNA + HTO + ADT + VDJ
# side of the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has run through all steps in at least two
#     RNA + HTO + ADT (+ VDJ) runs and is ready to aggregate the runs for a
#     single project.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
SAMPLES_ARRAY=($(ls -d $PROJECT_PATH/*/pipeline/RNA.FB.VDJ/*/))
OUTPUT_DIR=$PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ
OUTPUT_FILE=$OUTPUT_DIR/cellranger_rna_aggr.log
AGGR_CSV=$OUTPUT_DIR/aggr.csv
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" >> $OUTPUT_FILE

## Run test to assess if there is VDJ information for all (or no) runs
vdj_detected_array=()
for sample_path in "${SAMPLES_ARRAY[@]}"; do
  multi_dir=$sample_path/outs/per_sample_outs
  count_file=$sample_path/outs/molecule_info.h5
  if [ -d $multi_dir ]; then
    detected=1
  elif [ -f $count_file ]; then
    detected=0
  else
    echo "Detected incomplete run for sample $sample_name. Please check \
      directory ($sample_path) and try again."
    exit 1
  fi
  vdj_detected_array+=($detected)
done

# Total up VDJ runs detected
tot=0
for i in ${vdj_detected_array[@]}; do
  let tot+=$i
done

# Set multi aggregation mode - aggr.csv must be set up differently if so
if [ $tot == ${#vdj_detected_array[@]} ]; then
  multi_mode=1
elif [ $tot == 0 ]; then
  multi_mode=0
else
  echo "Detected a mix of runs with and without VDJ information. All samples \
must have same libraries. Either move desired libraries to a separate run \
directory or wait until consistent modalities are available."
  echo -e "\nSamples and if they have VDJ libraries (1) or not (0).\n"
  
  printf '%s\n' Sample VDJ_library | paste -sd ' '
  for idx in $(seq 1 ${#vdj_detected_array[@]}); do
    zero_idx=$(($idx-1))
    printf '%s\n' $SAMPLES_ARRAY[$zero_idx] $vdj_detected_array[$zero_idx] | paste -sd ' '
  done
fi

## Create aggregation csv
if [ $multi_mode ] && ! [ -f $AGGR_CSV ]; then
  printf '%s\n' sample_id molecule_h5 donor origin | paste -sd ',' >> $AGGR_CSV
  for sample_path in "${SAMPLES_ARRAY[@]}"; do
    sample_name=$(basename $sample_path)
    multi_molecule_file=$sample_path/outs/per_sample_outs/$sample_name/count/sample_molecule_info.h5
    printf '%s\n' $sample_name $multi_molecule_file D0 PBMC | paste -sd ',' >> $AGGR_CSV
  done
  echo "[WARNING] Aggregation not run. Aggregating with VDJ libraries requires a donor \
and origin specific sample sheet. One has been created with dummy values at $AGGR_CSV. \
Please EDIT FIRST, then run this aggregation again."
  exit 1
elif [ $multi_mode ] && [ -f $AGGR_CSV ]; then
  echo "Multi aggregation csv detected. Using this for aggregation and assuming it is accurate."
else
  printf '%s\n' sample_id molecule_h5 | paste -sd ',' >> $AGGR_CSV
  for sample_path in "${SAMPLES_ARRAY[@]}"; do
    sample_name=$(basename $sample_path)
    count_molecule_file=$sample_path/outs/molecule_info.h5
    printf '%s\n' $sample_name $count_molecule_file | paste -sd ',' >> $AGGR_CSV
  done
fi

## Run aggregation pipeline
echo "$(date) Running sample aggregation..." >> $OUTPUT_FILE
cellranger aggr --id=${PROJECT_NAME}_aggr \
	--csv $AGGR_CSV  --normalize none \
	--localcores=$NCPU --localmem=$MEM
