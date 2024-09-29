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
BED_DIR="${OUTPUT_DIR}/bed_files"
FASTA_DIR="${OUTPUT_DIR}/fasta_files"
VCF_OUTPUT_DIR="${OUTPUT_DIR}/vcf_files"
CONSENSUS_DIR="${OUTPUT_DIR}/consensus_sequences"
COMBINED_OUTPUT_FILE="${OUTPUT_DIR}/combined_output.tsv"

# Create output directories
mkdir -p ${BED_DIR} ${FASTA_DIR} ${VCF_OUTPUT_DIR} ${CONSENSUS_DIR}

# Define regions to extract
REGIONS=("3:247434-251178")

# Initialize the combined output file with a header
echo -e "SeqName\tChromosome\tLocation\tReferenceAllele" > ${COMBINED_OUTPUT_FILE}

# Prepare a temporary file to store combined genotype data
TEMP_COMBINED_FILE="temp_combined_output.tsv"
> ${TEMP_COMBINED_FILE}  # Clear the file before starting

# Process only files with "sorted_md" in their names
for VCF_FILE in ${VCF_DIR}/*sorted_md*.vcf.gz; do
    # Extract the strain name from the VCF filename
    STRAIN_NAME=$(basename ${VCF_FILE} | cut -d'.' -f1)
    
    # Make sure the VCF file is indexed
    if [ ! -f "${VCF_FILE}.tbi" ]; then
        bcftools index ${VCF_FILE}
    fi
    
    # Add strain name to the header of the combined output file
    sed -i "1s/$/\t${STRAIN_NAME}/" ${COMBINED_OUTPUT_FILE}

    for REGION in "${REGIONS[@]}"; do
        # Define variables for region-specific operations
        IFS=':' read -ra ADDR <<< "$REGION"
        CHROM=${ADDR[0]}
        IFS='-' read -ra POS <<< "${ADDR[1]}"
        START=${POS[0]}
        END=${POS[1]}
        REGION_STR="${CHROM}:${START}-${END}"
        
        # Define temporary VCF and FASTA files
        TEMP_VCF="${VCF_OUTPUT_DIR}/${STRAIN_NAME}_${CHROM}_${START}_${END}.vcf"
        REGION_FASTA_FILE="${CONSENSUS_DIR}/consensus_${CHROM}_${START}_${END}.fasta"
        
        # Filter the VCF for the specific region
        bcftools view -Oz -o ${TEMP_VCF}.gz -r ${REGION_STR} ${VCF_FILE}
        bcftools index ${TEMP_VCF}.gz
        
        # Extract genotype information and append to temporary combined file
        bcftools query -f '%CHROM\t%POS\t%REF[\t%GT]\n' -r ${REGION_STR} ${TEMP_VCF}.gz | \
        awk -v strain=${STRAIN_NAME} '{print $1 "\t" $2 "\t" $3 "\t" $4}' | \
        paste -d '\t' - ${TEMP_COMBINED_FILE} > temp_output.tsv && mv temp_output.tsv ${TEMP_COMBINED_FILE}
    done
done

# Combine the final genotype data into the combined output file
awk 'FNR==NR{a[$1 FS $2]=$3; next} {print $0, a[$1 FS $2]}' ${TEMP_COMBINED_FILE} >> ${COMBINED_OUTPUT_FILE}

# Clean up temporary files
rm ${TEMP_COMBINED_FILE}

echo "Combined genotype extraction and sequence alignment completed."
