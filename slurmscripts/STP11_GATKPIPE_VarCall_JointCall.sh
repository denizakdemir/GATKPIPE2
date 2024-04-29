#!/bin/bash
#SBATCH --job-name=GenomicAnalysis
#SBATCH --output=GenomicAnalysis_%j.log
#SBATCH --error=GenomicAnalysis_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load GATK module
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8 

# Define variables
REF_GENOME="0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"
GVCF_DIR="4_processing/GVCF"
DB_PATH="$GVCF_DIR/GenomicsDB"
FINAL_VCF="$GVCF_DIR/final_joint_called.vcf"
ORIGINAL_INTERVAL_LIST="$GVCF_DIR/interval.list"
INTERVAL_LIST="$GVCF_DIR/formatted.interval_list"
TMP_DIR="$GVCF_DIR/tmp"

# Perform joint genotyping
gatk GenotypeGVCFs \
   -R $REF_GENOME \
   -V gendb://$DB_PATH \
   -O $FINAL_VCF \
   -L $INTERVAL_LIST

echo "Joint variant calling completed."
