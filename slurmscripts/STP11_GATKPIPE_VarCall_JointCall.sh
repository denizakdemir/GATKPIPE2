#!/bin/bash
#SBATCH --job-name=ExtractSequences
#SBATCH --output=ExtractSequences_%j.log
#SBATCH --error=ExtractSequences_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load BCFtools module
module load BCFtools

# Define variables
VCF_FILE="4_processing/GVCF/final_joint_called.vcf"
OUTPUT_DIR="extracted_sequences"
REGION="3:247434-251178"  # Define the region you are interested in
OUTPUT_FASTA="${OUTPUT_DIR}/extracted_sequences_${REGION}.fasta"

# Create the output directory if it does not exist
mkdir -p ${OUTPUT_DIR}

# Extract the variant sequences from the VCF file for the specified region
echo "Extracting sequences for region ${REGION} from VCF file: ${VCF_FILE}"
bcftools consensus -f $VCF_FILE -r $REGION -o $OUTPUT_FASTA

# Check if the sequence extraction was successful
if [[ -s ${OUTPUT_FASTA} ]]; then
    echo "Sequence extraction completed successfully. Output saved to: ${OUTPUT_FASTA}"
else
    echo "Error: No sequences were extracted. Please check if the region contains variants or if the VCF file is correct."
fi

