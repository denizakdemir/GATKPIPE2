#!/bin/bash
#SBATCH --job-name=GenomicDB
#SBATCH --output=GenomicDB_%j.log
#SBATCH --error=GenomicDB_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=32:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load modules
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8 
module load SAMtools/1.16.1-GCC-11.3.0

# Relative paths from the script's current working directory
REF_DIR="./0_index/referenceIPO323"
REF_GENOME="$REF_DIR/Zymoseptoria_tritici.MG2.dna.toplevel.fa"

GVCF_DIR="./4_processing/GVCF"
DB_PATH="$GVCF_DIR/GenomicsDB"
FINAL_VCF="$GVCF_DIR/final_joint_called.vcf"
TMP_DIR="$GVCF_DIR/tmp"

# Ensure the script is executed from the correct directory
if [[ $PWD != */GATKPIPE2 ]]; then
  echo "This script must be run from the /home/deniz/GenomicData/github/GATKPIPE2 directory"
  exit 1
fi

# Index the reference genome using samtools
samtools faidx $REF_GENOME

# Create the sequence dictionary for the reference genome using GATK
gatk CreateSequenceDictionary -R $REF_GENOME

# Generate a BED file from the reference genome fasta index
awk -v FS="\t" -v OFS="\t" '{print $1 FS "0" FS ($2-1)}' "${REF_GENOME}.fai" > "$GVCF_DIR/genome.bed"

# Update INTERVAL_LIST variable to point to the newly created interval list
INTERVAL_LIST="$GVCF_DIR/genome.bed"
NEW_INTERVAL_LIST="$GVCF_DIR/formatted.interval_list"

# Reformat interval list
awk '{print $1 ":" $2 "-" $3}' $INTERVAL_LIST > $NEW_INTERVAL_LIST

# Create the temporary directory if it does not exist
if [ ! -d "$TMP_DIR" ]; then
  mkdir -p "$TMP_DIR"
  chmod 777 "$TMP_DIR"
fi

GenomicsDBImport command
gatk GenomicsDBImport \
   --genomicsdb-workspace-path $DB_PATH \
   -L $NEW_INTERVAL_LIST \
   --sample-name-map $GVCF_DIR/zymoseptoria.sample_map \
   --batch-size 16 \
   --overwrite-existing-genomicsdb-workspace \
   --tmp-dir $TMP_DIR \
   --reader-threads $SLURM_CPUS_PER_TASK

echo "Genomic DB created"
