#!/bin/bash
#SBATCH --job-name=GenomicAnalysis
#SBATCH --output=GenomicAnalysis_%j.log
#SBATCH --error=GenomicAnalysis_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load necessary modules
echo "Loading modules..."
module load Java/1.8.0_212
module load R/4.2.0-foss-2021b
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8 

# Define reference genome
echo "Setting reference genome..."
ref="/home/deniz/GenomicData/github/GATKPIPE/0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"

# Navigate to the directory containing the BAM files
echo "Changing to BAM files directory..."
cd 3_mapping/aligned_samples 

# Initial variant calling and filtration
echo "Creating Initial_variants directory..."
mkdir -p ../../4_processing/Initial_variants

MAX_JOBS=30  # Control the maximum number of parallel jobs
running_jobs=0

echo "Starting variant calling and filtration..."
for i in *_sorted_md.bam; do
  base=$(basename ${i} _sorted_md.bam)

  echo "Processing file: $i"

  # Variant calling
  echo "Running HaplotypeCaller for $i..."
  (
  gatk HaplotypeCaller -R $ref -I $i -O "../../4_processing/Initial_variants/${base}_raw_variants.vcf"

  # SNP selection and filtration
  echo "Selecting and filtering SNPs for $i..."
  gatk SelectVariants -R $ref -V "../../4_processing/Initial_variants/${base}_raw_variants.vcf" -select-type SNP -O "../../4_processing/Initial_variants/${base}_raw_snps.vcf"
  gatk VariantFiltration -R $ref -V "../../4_processing/Initial_variants/${base}_raw_snps.vcf" -O "../../4_processing/Initial_variants/${base}_filtered_snps.vcf" \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "SOR_filter" -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
  gatk SelectVariants --exclude-filtered -V "../../4_processing/Initial_variants/${base}_filtered_snps.vcf" -O "../../4_processing/Initial_variants/${base}_bqsr_snps.vcf"

  # Indel selection and filtration
  echo "Selecting and filtering Indels for $i..."
  gatk SelectVariants -R $ref -V "../../4_processing/Initial_variants/${base}_raw_variants.vcf" -select-type INDEL -O "../../4_processing/Initial_variants/${base}_raw_indels.vcf"
  gatk VariantFiltration -R $ref -V "../../4_processing/Initial_variants/${base}_raw_indels.vcf" -O "../../4_processing/Initial_variants/${base}_filtered_indels.vcf" \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0"
  gatk SelectVariants --exclude-filtered -V "../../4_processing/Initial_variants/${base}_filtered_indels.vcf" -O "../../4_processing/Initial_variants/${base}_bqsr_indels.vcf" 
  ) &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1  # wait for 1 second before checking again
            running_jobs=$(jobs -p | wc -l)  # update the count of running jobs
        done
    fi
done
wait

# Base Quality Score Recalibration (BQSR)
echo "Creating BQSR directory and starting Base Quality Score Recalibration..."
mkdir -p ../../4_processing/BQSR

running_jobs=0
MAX_JOBS=10 # runs out of memory if you increase this number
for i in *_sorted_md.bam; do
  base=$(basename ${i} _sorted_md.bam)
  echo "Recalibrating base quality for $i..."
  (
  # Base recalibration
  gatk BaseRecalibrator -R $ref -I $i --known-sites "../../4_processing/Initial_variants/${base}_bqsr_snps.vcf" --known-sites "../../4_processing/Initial_variants/${base}_bqsr_indels.vcf" -O "../../4_processing/BQSR/${base}_recal_data.table"
  gatk ApplyBQSR -R $ref -I $i -bqsr "../../4_processing/BQSR/${base}_recal_data.table" -O "../../4_processing/BQSR/${base}_recal.bam"

  # Post-recalibration analysis
  gatk BaseRecalibrator -R $ref -I "../../4_processing/BQSR/${base}_recal.bam" --known-sites "../../4_processing/Initial_variants/${base}_bqsr_snps.vcf" --known-sites "../../4_processing/Initial_variants/${base}_bqsr_indels.vcf" -O "../../4_processing/BQSR/${base}_post_recal_data.table"
  gatk AnalyzeCovariates -before "../../4_processing/BQSR/${base}_recal_data.table" -after "../../4_processing/BQSR/${base}_post_recal_data.table" -plots "../../4_processing/BQSR/${base}_recalibration_plots.pdf" 
  ) &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1  # wait for 1 second before checking again
            running_jobs=$(jobs -p | wc -l)  # update the count of running jobs
        done
    fi
done
wait

echo "Genomic analysis completed."
