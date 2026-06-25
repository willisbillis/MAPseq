#!/bin/bash
# generate_hla_vcf.sh
# written by Gemini 3 Flash (Preview) (GitHub Copilot) Jun 2026
# This script generates a multi-sample VCF from known HLA alleles for souporcell demultiplexing.
#
# INPUT: hla_alleles.csv (format: donor_id,hla_allele)
#        Example: 
#        DonorA,HLA-A*02:01
#        DonorA,HLA-B*07:02
#        DonorB,HLA-A*01:01
#        DonorB,HLA-B*08:01

################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this script
HLA_ALLELE_FILE=$PROJECT_PATH/pipeline/RNA.FB.VDJ/hla_alleles.csv
HLA_DB_URL="https://raw.githubusercontent.com/ANHIG/IMGTHLA/master/hla_gen.fasta"
HLA_DB=$PROJECT_PATH/pipeline/RNA.FB.VDJ/hla_gen.fasta
REF_FASTA=$GEX_REF_PATH/fasta/genome.fa
OUTPUT_DIR=$PROJECT_PATH/pipeline/RNA.FB.VDJ/hla_vcf_gen
MULTI_VCF=$OUTPUT_DIR/merged_hla_genotypes.vcf.gz

# Tool paths (assumed to be in PATH or conda env)
# minimap2, samtools, bcftools
################################################################################

# Preflight checks
if [ ! -f $HLA_ALLELE_FILE ]; then
  echo "ERROR: HLA allele mapping file not found at $HLA_ALLELE_FILE"
  echo "Please create this file with 'donor_id,hla_allele' format."
  exit 1
fi

if [ ! -f $REF_FASTA ]; then
  echo "ERROR: Reference genome FASTA not found at $REF_FASTA"
  exit 1
fi

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# 1. Download HLA database if not present
if [ ! -f $HLA_DB ]; then
  echo "Downloading IPD-IMGT/HLA database..."
  wget -O $HLA_DB $HLA_DB_URL
fi

# 2. Extract sequences for each donor
echo "Extracting HLA sequences and calling variants..."
DONORS=($(tail -n +2 $HLA_ALLELE_FILE | cut -d',' -f1 | sort | uniq))

for donor in "${DONORS[@]}"; do
  echo "Processing donor: $donor"
  DONOR_FASTA=$OUTPUT_DIR/${donor}_hla.fasta
  rm -f $DONOR_FASTA
  
  ALLELES=($(grep "^$donor," $HLA_ALLELE_FILE | cut -d',' -f2))
  
  for allele in "${ALLELES[@]}"; do
    # Remove 'HLA-' prefix if present for matching in hla_gen.fasta
    match_str=$(echo $allele | sed 's/^HLA-//')
    # Matching header like ">HLA:HLA00001 A*01:01:01:01"
    # We use awk to extract the sequence for the matched header
    awk -v pattern=" $match_str" '
      BEGIN { RS=">"; FS="\n" }
      $1 ~ pattern { print ">"$0 }
    ' $HLA_DB >> $DONOR_FASTA
  done

  if [ -s $DONOR_FASTA ]; then
    # 3. Virtual Alignment to Reference Genome
    minimap2 -ax asm5 -t $NCPU $REF_FASTA $DONOR_FASTA | \
      samtools view -uS - | \
      samtools sort -o ${donor}_hla_sorted.bam -
    samtools index ${donor}_hla_sorted.bam

    # 4. Variant Calling & Filtering (SNPs only)
    bcftools mpileup -f $REF_FASTA ${donor}_hla_sorted.bam | \
      bcftools call -mv -Ov | \
      bcftools view -v snps -o ${donor}_variants.vcf
    
    # Set sample name in VCF
    echo "$donor" > ${donor}_sample_name.txt
    bcftools reheader -s ${donor}_sample_name.txt ${donor}_variants.vcf -o ${donor}_final.vcf
    bgzip -f ${donor}_final.vcf
    bcftools index ${donor}_final.vcf.gz
  else
    echo "WARNING: No sequences found for donor $donor"
  fi
done

# 5. Sample Merging
echo "Merging donor VCFs..."
VCF_LIST=$(ls ${OUTPUT_DIR}/*_final.vcf.gz)
bcftools merge -Oz -o $MULTI_VCF $VCF_LIST
bcftools index $MULTI_VCF

echo "HLA Multi-Sample VCF generated at $MULTI_VCF"
