#!/bin/bash
#SBATCH --job-name=ExtractSequences
#SBATCH --output=ExtractSequences_%j.log
#SBATCH --error=ExtractSequences_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load required modules
module load BCFtools/1.9-GCC-8.2.0-2.31.1
module load BEDTools/2.28.0-GCC-8.2.0-2.31.1

# Define variables
REFERENCE_GENOME="0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"
VCF_DIR="4_processing/GVCF"
OUTPUT_DIR="extracted_sequences"
OUTPUT_DIR_CONSENSUS="${OUTPUT_DIR}/consensus_sequences"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR_CONSENSUS}

# Define regions to extract - 500bp flanking regions
REGIONS=("3:247434-251178")

# Extract SNP positions and create BED files
for REGION in "${REGIONS[@]}"; do
    IFS=':' read -ra ADDR <<< "$REGION"
    CHROM=${ADDR[0]}
    IFS='-' read -ra POS <<< "${ADDR[1]}"
    START=${POS[0]}
    END=${POS[1]}
    
    # Create a bed file for each region
    echo -e "${CHROM}\t${START}\t${END}" > "${OUTPUT_DIR}/${CHROM}_${START}_${END}.bed"
    
    # Use bedtools to extract sequences
    bedtools getfasta -fi "${REFERENCE_GENOME}" -bed "${OUTPUT_DIR}/${CHROM}_${START}_${END}.bed" -fo "${OUTPUT_DIR}/${CHROM}_${START}_${END}.fasta"
done

echo "Reference sequence extraction completed."

# Initialize or clear region-specific consensus FASTA files
for REGION in "${REGIONS[@]}"; do
    IFS=':' read -ra ADDR <<< "$REGION"
    CHROM=${ADDR[0]}
    IFS='-' read -ra POS <<< "${ADDR[1]}"
    START=${POS[0]}
    END=${POS[1]}
    REGION_FASTA_FILE="${OUTPUT_DIR_CONSENSUS}/consensus_${CHROM}_${START}_${END}.fasta"
    > ${REGION_FASTA_FILE} # Clear the file if it exists or create it if not
done

# Iterate over VCF files that have 'sorted_md' in their names
for VCF_FILE in ${VCF_DIR}/*sorted_md*.vcf.gz; do
    STRAIN_NAME=$(basename ${VCF_FILE} | cut -d'.' -f1)

    # Ensure the VCF file is indexed
    if [ ! -f "${VCF_FILE}.tbi" ]; then
        bcftools index ${VCF_FILE}
    fi

    for REGION in "${REGIONS[@]}"; do
        # Define variables for region-specific operations
        IFS=':' read -ra ADDR <<< "$REGION"
        CHROM=${ADDR[0]}
        IFS='-' read -ra POS <<< "${ADDR[1]}"
        START=${POS[0]}
        END=${POS[1]}
        REGION_STR="${CHROM}:${START}-${END}"
        REGION_FILE="${OUTPUT_DIR}/${CHROM}_${START}_${END}.fasta"
        TEMP_VCF="${OUTPUT_DIR}/${STRAIN_NAME}_${CHROM}_${START}_${END}.vcf"
        REGION_FASTA_FILE="${OUTPUT_DIR_CONSENSUS}/consensus_${CHROM}_${START}_${END}.fasta"

        # Filter VCF for the specific region
        bcftools view -Oz -o ${TEMP_VCF}.gz -r ${REGION_STR} ${VCF_FILE}
        bcftools index ${TEMP_VCF}.gz

        # Generate a temporary consensus sequence
        bcftools consensus -f ${REGION_FILE} -o temp_consensus.fasta ${TEMP_VCF}.gz

        # Add header with strain name to the temporary consensus sequence
        echo ">${STRAIN_NAME}" > temp_header.fasta
        cat temp_consensus.fasta >> temp_header.fasta

        # Append the consensus sequence with header to the region-specific FASTA file
        cat temp_header.fasta >> ${REGION_FASTA_FILE}

        # Clean up temporary files
        rm ${TEMP_VCF}.gz ${TEMP_VCF}.gz.csi temp_consensus.fasta temp_header.fasta
    done
done

echo "Consensus sequence compilation for individual strains completed."
