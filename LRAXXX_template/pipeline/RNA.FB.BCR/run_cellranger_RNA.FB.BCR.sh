#!/bin/bash
#
# run_cellranger_RNA.FB_multi.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger multi pipeline formatted for the RNA + HTO + ADT + BCR side of
# the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has generated the fastqs and
#       has formatted a feature reference table for the 
#       HTO/ADT demultiplexing.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path/$RNA_FLOWCELL_ID
OUTPUT_DIR=$PROJECT_PATH/pipeline/RNA.FB.BCR
OUTPUT_FILE=$OUTPUT_DIR/cellranger_rna.fb.bcr_mapping.log
################################################################################
mkdir -p $OUTPUT_DIR
sample_name_col=$(cut -d, -f2 $PROJECT_PATH/data/${PROJECT_NAME}.RNA.sampleManifest.csv)
sample_names=$(printf -- '%s\n' "${sample_name_col[@]}" | grep -v Sample | uniq)
echo $sample_names # verbose for testing
bcr_samples=$(printf -- '%s\n' "${sample_names[@]}" | grep '.*${BCR_NAMING_ID}.*')
rna_samples=$(printf -- '%s\n' "${sample_names[@]}" | grep '.*${GEX_NAMING_ID}.*')
echo $rna_samples # verbose for testing

CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" >> $OUTPUT_FILE
GEX_REF_version=$(echo $GEX_REF_PATH | grep -Po '(?<=refdata-gex-)[^;]+')
echo "$(date) Using transcriptome reference $GEX_REF_version located at $GEX_REF_PATH" >> $OUTPUT_FILE
echo "$(date) Using HTO/ADT feature reference located at $GEX_FEAT_REF_PATH" >> $OUTPUT_FILE
if [ ${#bcr_samples[@]} != 0 ]; then
    VDJ_REF_version=$(echo $VDJ_REF_PATH | grep -Po '(?<=refdata-cellranger-vdj-)[^;]+')
    echo "$(date) Using vdj reference $VDJ_REF_version located at $VDJ_REF_PATH" >> $OUTPUT_FILE
fi
cd $OUTPUT_DIR

for sample in "${rna_samples[@]}"; do
    echo "$(date) Running sample ${sample}..." >> $OUTPUT_FILE

    if [ ${#bcr_samples[@]} != 0 ]; then
        ## BCR SAMPLES DETECTED, RUN CR MULTI
        # Create the config csv for the sample being run
        SAMPLE_CONFIG_CSV=$OUTPUT_DIR/${sample}_config.csv
        echo "[gene-expression]" >> $SAMPLE_CONFIG_CSV
        printf '%s\n' reference $GEX_REF_PATH | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        echo "[feature]" >> $SAMPLE_CONFIG_CSV
        printf '%s\n' reference $GEX_FEAT_REF_PATH | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        echo "[vdj]" >> $SAMPLE_CONFIG_CSV
        printf '%s\n' reference $VDJ_REF_PATH | paste -sd ',' >> $SAMPLE_CONFIG_CSV

        echo "[libraries]" >> $SAMPLE_CONFIG_CSV
        printf '%s\n' fastq_id fastqs feature_types | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        gex_sample=$sample
        hashtag_sample=$(sed -n -e 's/GEX/CSP/p' $sample) # this sed command substitutes the 'GEX' part of the sample name with 'CSP'
        vdj_sample=$(sed -n -e 's/GEX/XP/p' $sample) # this sed command substitutes the 'GEX' part of the sample name with 'XP'
        printf '%s\n' $gex_sample $FASTQ_PATH Gene Expression | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        printf '%s\n' $hashtag_sample $FASTQ_PATH Antibody Capture | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        printf '%s\n' $vdj_sample $FASTQ_PATH VDJ | paste -sd ',' >> $SAMPLE_CONFIG_CSV

        # Run the Cell Ranger multi command for the sample
        cellranger multi --id $sample \
        --csv $SAMPLE_CONFIG_CSV \
        --localcores $NCPU --localmem $MEM
    else
        ## NO BCR SAMPLES DETECTED, RUN CR COUNT
        # Create the config csv for the sample being run
        SAMPLE_CONFIG_CSV=$OUTPUT_DIR/${sample}_config.csv
        printf '%s\n' fastqs sample library_type | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        fq_path=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path/$RNA_FLOWCELL_ID
        gex_sample=$sample
        hashtag_sample=$(sed -n -e 's/GEX/CSP/p' $sample) # this sed command substitutes the 'GEX' part of the sample name with 'CSP'
        printf '%s\n' $FASTQ_PATH $gex_sample Gene Expression | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        printf '%s\n' $FASTQ_PATH $hashtag_sample Antibody Capture | paste -sd ',' >> $SAMPLE_CONFIG_CSV

        # Run the Cell Ranger count command for the sample
        cellranger count --id $sample \
            --transcriptome $GEX_REF_PATH --feature-ref $GEX_FEAT_REF_PATH \
            --libraries $SAMPLE_CONFIG_CSV \
            --localcores $NCPU --localmem $MEM
    fi
done