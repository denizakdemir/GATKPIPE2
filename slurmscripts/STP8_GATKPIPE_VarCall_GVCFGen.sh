#!/bin/bash
#SBATCH --job-name=GVCF_Generation
#SBATCH --output=GVCF_Generation_%j.log
#SBATCH --error=GVCF_Generation_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load necessary modules
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8 

cd 4_processing/BQSR
# Define reference genome
ref="/home/deniz/GenomicData/github/GATKPIPE/0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"

# Create GVCF directory
echo "Creating GVCF directory..."
mkdir -p ../GVCF

MAX_JOBS=30  # Control the maximum number of parallel jobs
running_jobs=0

# GVCF generation
echo "Starting GVCF generation..."
for i in *_recal.bam; do
  base=$(basename ${i} _recal.bam)
  output_file="../GVCF/${base}.g.vcf.gz"

  # Check if the output file already exists
  if [ ! -f "$output_file" ]; then
    echo "Processing file: $i"

    gatk --java-options "-Xmx4g" HaplotypeCaller \
     -R $ref \
     -I $i \
     -O "$output_file" \
     -ploidy 1 \
     --max-alternate-alleles 2 \
     -ERC GVCF &
    
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1  # wait for 1 second before checking again
            running_jobs=$(jobs -p | wc -l)  # update the count of running jobs
        done
    fi
  else
    echo "Skipping already processed file: $i"
  fi
done
wait

echo "GVCF generation completed."



echo "GVCF generation completed."
