#!/bin/bash
#SBATCH --job-name=ExtractVCFSequences
#SBATCH --output=ExtractVCFSequences_%j.log
#SBATCH --error=ExtractVCFSequences_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load required modules
module load BCFtools

# Define variables
VCF_DIR="4_processing/GVCF"
OUTPUT_DIR="extracted_sequences"
OUTPUT_DIR_CONSENSUS="${OUTPUT_DIR}/consensus_sequences"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR_CONSENSUS}

# Define the region to extract
REGIONS=("3:247434-251178")

# Initialize or clear region-specific consensus FASTA file
for REGION in "${REGIONS[@]}"; do
    IFS=':' read -ra ADDR <<< "$REGION"
    CHROM=${ADDR[0]}
    IFS='-' read -ra POS <<< "${ADDR[1]}"
    START=${POS[0]}
    END=${POS[1]}
    REGION_FASTA_FILE="${OUTPUT_DIR_CONSENSUS}/vcf_sequences_${CHROM}_${START}_${END}.fasta"
    > ${REGION_FASTA_FILE} # Clear or initialize the file
    echo "Initialized VCF-based sequence FASTA file: ${REGION_FASTA_FILE}"
done

# Create an empty file to track processed VCF files
PROCESSED_VCF_FILES="${OUTPUT_DIR_CONSENSUS}/processed_vcf_files.txt"
> ${PROCESSED_VCF_FILES}

# Loop through VCF files and extract variant sequences for each strain
for VCF_FILE in $(ls ${VCF_DIR}/*.vcf.gz | grep -v 'sorted'); do
    STRAIN_NAME=$(basename ${VCF_FILE} | cut -d'.' -f1)
    echo "Processing VCF file: ${VCF_FILE} (Strain: ${STRAIN_NAME})"
    
    # Add the VCF file name to the processed file list
    echo "${VCF_FILE}" >> ${PROCESSED_VCF_FILES}

    for REGION in "${REGIONS[@]}"; do
        IFS=':' read -ra ADDR <<< "$REGION"
        CHROM=${ADDR[0]}
        IFS='-' read -ra POS <<< "${ADDR[1]}"
        START=${POS[0]}
        END=${POS[1]}
        REGION_STR="${CHROM}:${START}-${END}"
        REGION_FASTA_FILE="${OUTPUT_DIR_CONSENSUS}/vcf_sequences_${CHROM}_${START}_${END}.fasta"

        # Filter VCF for the specific region and extract sequences
        echo "Extracting variant information for region: ${REGION_STR}"

        # Extract sequence based on genotype
        bcftools query -r ${REGION_STR} -f '[%CHROM:%POS %REF %ALT %GT\n]' ${VCF_FILE} > temp_variants.txt

        # Build sequence for the strain by interpreting the genotypes
        echo ">${STRAIN_NAME}" >> ${REGION_FASTA_FILE}
        SEQUENCE=""
        while read -r LINE; do
            CHROM_POS=$(echo "$LINE" | cut -d ' ' -f 1)
            REF=$(echo "$LINE" | cut -d ' ' -f 2)
            ALT=$(echo "$LINE" | cut -d ' ' -f 3)
            GENOTYPE=$(echo "$LINE" | cut -d ' ' -f 4
