#!/bin/bash
#SBATCH --job-name=ExtractConsensusSequences
#SBATCH --output=ExtractConsensusSequences_%j.log
#SBATCH --error=ExtractConsensusSequences_%j.err
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
REFERENCE_GENOME="0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"
OUTPUT_DIR="extracted_sequences"
OUTPUT_DIR_CONSENSUS="${OUTPUT_DIR}/consensus_sequences"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR_CONSENSUS}

# Define the region to extract
REGIONS=("3:247434-251178")

# Process each VCF file and extract the consensus sequence
for VCF_FILE in $(ls ${VCF_DIR}/*.vcf.gz | grep -v 'sorted'); do
    STRAIN_NAME=$(basename ${VCF_FILE} | cut -d'.' -f1)
    echo "Processing VCF file: ${VCF_FILE} (Strain: ${STRAIN_NAME})"

    for REGION in "${REGIONS[@]}"; do
        IFS=':' read -ra ADDR <<< "$REGION"
        CHROM=${ADDR[0]}
        IFS='-' read -ra POS <<< "${ADDR[1]}"
        START=${POS[0]}
        END=${POS[1]}
        REGION_STR="${CHROM}:${START}-${END}"
        REGION_FASTA_FILE="${OUTPUT_DIR_CONSENSUS}/consensus_${STRAIN_NAME}_${CHROM}_${START}_${END}.fasta"

        # Use bcftools consensus to apply the variants from the VCF file to the reference sequence
        bcftools consensus -f ${REFERENCE_GENOME} -r ${REGION_STR} ${VCF_FILE} > ${REGION_FASTA_FILE}

        echo "Consensus sequence for strain ${STRAIN_NAME} in region ${REGION_STR} written to ${REGION_FASTA_FILE}"
    done
done

echo "Consensus sequence extraction completed."
