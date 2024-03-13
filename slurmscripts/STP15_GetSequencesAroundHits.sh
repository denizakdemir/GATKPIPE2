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
VCF_FILE="5_annotation/snpEff_snps/geno_filtered_snps.vcf"
OUTPUT_DIR="extracted_sequences"
mkdir -p ${OUTPUT_DIR}


# Define regions to extract - 500bp flanking regions
REGIONS=("3:1981024-1982024" "7:2122486-2123486" "2:3164088-3165088")

# Extract SNP positions and create BED files
for REGION in "${REGIONS[@]}"; do
    # Extract the chromosome and position
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
# Define variables
REFERENCE_GENOME="0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"
VCF_DIR="4_processing/GVCF"
OUTPUT_DIR="extracted_sequences"
mkdir -p ${OUTPUT_DIR}

# Define SNP positions
SNP_POSITIONS=("3:1981524" "7:2122986" "2:3164588")

# Iterate over individual VCF files
for VCF_FILE in ${VCF_DIR}/S*.vcf.gz; do
    # Extract the strain name from the VCF filename
    STRAIN_NAME=$(basename ${VCF_FILE} | cut -d'.' -f1)
    
    # Make sure the VCF file is indexed
    if [ ! -f "${VCF_FILE}.tbi" ]; then
        bcftools index ${VCF_FILE}
    fi
    
    # Extract genotypes for each SNP position for the current strain
    for SNP_POS in "${SNP_POSITIONS[@]}"; do
        # Define output file name based on strain and SNP position
        OUTPUT_FILE="${OUTPUT_DIR}/${STRAIN_NAME}_$(echo ${SNP_POS} | tr ':' '_')_genotypes.txt"
        
        # Use bcftools query to extract genotypes
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -r ${SNP_POS} ${VCF_FILE} > ${OUTPUT_FILE}
    done
done

echo "Genotype extraction for individual strains completed."


# Define directories and files
REFERENCE_GENOME="0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"
VCF_DIR="4_processing/GVCF"
OUTPUT_DIR_CONSENSUS="${OUTPUT_DIR}/consensus_sequences"
mkdir -p ${OUTPUT_DIR_CONSENSUS}

# Define regions (same as before)
REGIONS=("3:1981024-1982024" "7:2122486-2123486" "2:3164088-3165088")

# Iterate over individual VCF files to create consensus sequences
for VCF_FILE in ${VCF_DIR}/S*.vcf.gz; do
    STRAIN_NAME=$(basename ${VCF_FILE} | cut -d'.' -f1)

    # Index the VCF file if not already indexed
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
        CONSENSUS_FILE="${OUTPUT_DIR_CONSENSUS}/${STRAIN_NAME}_${CHROM}_${START}_${END}_consensus.fasta"

        # Use bcftools consensus to apply the VCF variants to the reference sequence
        bcftools consensus -f ${REGION_FILE} -r ${REGION_STR} ${VCF_FILE} -o ${CONSENSUS_FILE}
    done
done

echo "Consensus sequence extraction for individual strains completed."
