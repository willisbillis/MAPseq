#!/bin/bash
#
# run_cellranger_RNA.FB.VDJ.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger multi pipeline formatted for the RNA + HTO + ADT (+ VDJ) side of
# the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has generated the fastqs and
#       has formatted a feature reference table for the 
#       HTO/ADT demultiplexing.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs
OUTPUT_DIR=$PROJECT_PATH/pipeline/RNA.FB.VDJ
OUTPUT_FILE=$OUTPUT_DIR/cellranger_rna.fb.vdj_mapping.log
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/reports

sample_name_col=$(cut -d, -f2 $PROJECT_PATH/data/${PROJECT_NAME}.RNA.sampleManifest.csv)
sample_names=$(printf -- '%s ' "${sample_name_col[@]}" | grep -v Sample | uniq)

vdj_samples=($(printf -- '%s ' "${sample_names[@]}" | grep .*${VDJ_NAMING_ID}.*))
rna_samples=($(printf -- '%s ' "${sample_names[@]}" | grep .*${GEX_NAMING_ID}.*))

CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" >> $OUTPUT_FILE
GEX_REF_version=$(echo $GEX_REF_PATH | grep -Po '(?<=refdata-gex-)[^;]+')
echo "$(date) Using transcriptome reference $GEX_REF_version located at $GEX_REF_PATH" >> $OUTPUT_FILE
echo "$(date) Using HTO/ADT feature reference located at $GEX_FEAT_REF_PATH" >> $OUTPUT_FILE
if [ ${#vdj_samples[@]} != 0 ]; then
    VDJ_REF_version=$(echo $VDJ_REF_PATH | grep -Po '(?<=refdata-cellranger-vdj-)[^;]+')
    echo "$(date) Using vdj reference $VDJ_REF_version located at $VDJ_REF_PATH" >> $OUTPUT_FILE
fi

for sample in "${rna_samples[@]}"; do
    echo "$(date) Running sample ${sample}..." >> $OUTPUT_FILE
    SAMPLE_CONFIG_CSV=$OUTPUT_DIR/${sample}_config.csv

    if [ ${#vdj_samples[@]} != 0 ]; then
        ## VDJ SAMPLES DETECTED, RUN CR MULTI
        # Create the config csv for the sample being run
        echo "[gene-expression]" >> $SAMPLE_CONFIG_CSV
        printf '%s\n' reference $GEX_REF_PATH | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        printf '%s\n' create-bam true | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        echo "[feature]" >> $SAMPLE_CONFIG_CSV
        printf '%s\n' reference $GEX_FEAT_REF_PATH | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        echo "[vdj]" >> $SAMPLE_CONFIG_CSV
        printf '%s\n' reference $VDJ_REF_PATH | paste -sd ',' >> $SAMPLE_CONFIG_CSV

        echo "[libraries]" >> $SAMPLE_CONFIG_CSV
        printf '%s\n' fastq_id fastqs feature_types | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        gex_sample=$sample
        # this sed command substitutes the $GEX_NAMING_ID part of the sample name with $GEX_FEAT_NAMING_ID
        hashtag_sample=$(echo $sample | sed -n -e "s/$GEX_NAMING_ID/$GEX_FEAT_NAMING_ID/p")
        # this sed command substitutes the $GEX_NAMING_ID part of the sample name with $VDJ_NAMING_ID
        vdj_sample=$(echo $sample | sed -n -e "s/$GEX_NAMING_ID/$VDJ_NAMING_ID/p")
        printf '%s\n' $gex_sample $FASTQ_PATH 'Gene Expression' | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        printf '%s\n' $hashtag_sample $FASTQ_PATH 'Antibody Capture' | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        printf '%s\n' $vdj_sample $FASTQ_PATH VDJ | paste -sd ',' >> $SAMPLE_CONFIG_CSV

        # Run the Cell Ranger multi command for the sample
        cellranger multi --id $(echo $sample | sed -n -e "s/$GEX_NAMING_ID//p") \
        --csv $SAMPLE_CONFIG_CSV \
        --localcores $NCPU --localmem $MEM

        cp $sample/outs/per_sample_outs/$sample/web_summary.html $OUTPUT_DIR/reports/mapping.report_${sample}.html
    else
        ## NO VDJ SAMPLES DETECTED, RUN CR COUNT
        # Create the config csv for the sample being run
        printf '%s\n' fastqs sample library_type | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        gex_sample=$sample
        # this sed command substitutes the $GEX_NAMING_ID part of the sample name with $GEX_FEAT_NAMING_ID
        hashtag_sample=$(echo $sample | sed -n -e "s/$GEX_NAMING_ID/$GEX_FEAT_NAMING_ID/p")
        printf '%s\n' $FASTQ_PATH $gex_sample 'Gene Expression' | paste -sd ',' >> $SAMPLE_CONFIG_CSV
        printf '%s\n' $FASTQ_PATH $hashtag_sample 'Antibody Capture' | paste -sd ',' >> $SAMPLE_CONFIG_CSV

        # Run the Cell Ranger count command for the sample
        run_id=$(echo $sample | sed -n -e "s/$GEX_NAMING_ID//p")
        cellranger count --id $run_id \
            --create-bam=true \
            --transcriptome $GEX_REF_PATH --feature-ref $GEX_FEAT_REF_PATH \
            --libraries $SAMPLE_CONFIG_CSV \
            --localcores $NCPU --localmem $MEM
        
        cp $run_id/outs/web_summary.html $OUTPUT_DIR/reports/mapping.report_${run_id}.html
    fi
done