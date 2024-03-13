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
VCF_FILE="5_annotation/snpEff_snps/geno_filtered_snps.ann.vcf.gz"
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

echo "Sequence extraction completed."


# Define variables
OUTPUT_DIR="extracted_sequences"
mkdir -p ${OUTPUT_DIR}

# Define SNP positions
SNP_POSITIONS=("3:1981524" "7:2122986" "2:3164588")

# Extract genotypes for each SNP position
for SNP_POS in "${SNP_POSITIONS[@]}"; do
    # Use bcftools query to extract genotypes
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -r $SNP_POS $VCF_FILE > "${OUTPUT_DIR}/${SNP_POS}_genotypes.txt"
done

echo "Genotype extraction completed."
