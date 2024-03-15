#install.packages("vcfR")
#install.packages("dplyr")
library(vcfR)
library(dplyr)


# Load necessary libraries
library(vcfR)
library(dplyr)
library(tidyr)

# Read the VCF file
vcf_file <- "~/Desktop/geno_filtered_snps.ann.vcf"
vcf <- read.vcfR(vcf_file)

# Extracting map information
map_info <- data.frame(
  CHROM = vcf@fix[, 1],
  POS = vcf@fix[, 2],
  ID = vcf@fix[, 3],
  REF = vcf@fix[, 4],
  ALT = vcf@fix[, 5]
)

table(map_info$CHROM)

save(map_info, file = "map_info.RData")



# Extracting data from VCF
geno <- extract.gt(vcf, element = "GT", as.numeric = TRUE)

save(geno, file = "genoGT.RData")


load("genoGT.RData")
str(geno)

colnames(geno)[1:100]
colnames(geno)[duplicated(colnames(geno))]
geno[1:5, 1:8]
dim(geno)
# Identifying monomorphic loci
library(dplyr)
geno[1:5,1:5]
# Assuming 'geno' is your data frame
monomorphic_loci <- apply(geno, 1, function(x) length(unique(x, na.rm = TRUE)) == 1)
sum(monomorphic_loci) 

# Remove monomorphic loci
geno <- geno[!monomorphic_loci, ]
dim(geno)


# Step 2: Calculate Missingness for Each SNP
missingness <- apply(geno, 1, function(x) mean(is.na(x))) * 100
head(missingness)
sum(missingness > 30)

# Assuming your data is in a dataframe named 'geno'
# Replace 'geno' with the actual name of your dataframe

library(dplyr)
library(tidyr)
library(progress)

# Function to count alleles
count_alleles <- function(x) {
  alleles <- unlist(strsplit(na.omit(x), split="/"))
  table(factor(alleles, levels = c("0", "1")))
}

# Function to determine minor allele count for a genotype
minor_allele_count <- function(genotype, minor_allele) {
  if (is.na(genotype)) {
    return(NA)
  }
  alleles <- unlist(strsplit(genotype, split="/", fixed=TRUE))
  sum(alleles == minor_allele, na.rm=TRUE)
}

# Initialize progress bar
num_snps <- nrow(geno)


# Apply the functions to each SNP
results <- apply(geno[1:10000, ], 1, function(snp) {
  allele_counts <- count_alleles(snp)
  if(length(allele_counts) < 2) {
    return(c(MAF=NA, Minor_Allele_Counts=rep(NA, length(snp))))
  }
  minor_allele <- names(which.min(allele_counts))
  maf <- min(allele_counts) / sum(allele_counts)
  minor_allele_counts <- sapply(snp, minor_allele_count, minor_allele)
  c(MAF=maf, Minor_Allele_Counts=minor_allele_counts)
})

# Convert results to a data frame
results_df <- as.data.frame(t(results))

# Add row names as a column (if you want SNP identifiers)
results_df$SNP <- rownames(results_df)

results_df[1:10,1:10]

hist(results_df$MAF)
# Calculating the numbers (variant count)
variant_count <- vcf@fix[, 1:8] %>% nrow()

# Calculating missingness
missingness <- apply(geno, 1, function(x) sum(is.na(x)) / length(x))

# Calculating MAF (Minor Allele Frequency)
calculate_maf <- function(geno_row) {
  alleles <- table(geno_row, useNA = "ifany")
  freq <- alleles / sum(alleles)
  min(freq[names(freq) != "<NA>"])
}


maf <- apply(geno, 1, calculate_maf)

# Creating a data frame for the report
report <- data.frame(
  Variant_ID = vcf@fix[, 1],
  Missingness = missingness,
  MAF = maf
)

# Displaying a snippet of the report
head(report)

# Writing the full report to a CSV file
write.csv(report, "vcf_report.csv", row.names = FALSE)



if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
BiocManager::install("GenomicRanges")

library(VariantAnnotation)
library(GenomicRanges)