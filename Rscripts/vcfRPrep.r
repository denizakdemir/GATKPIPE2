install.packages("vcfR")
install.packages("data.table")

library(vcfR)
library(data.table)

vcf <- read.vcfR("~/Desktop/geno_filtered_snps.ann.vcf")
summary(vcf)
# Extract the genotype matrix
geno_matrix <- extract.gt(vcf, element = c('GT'))
dim(geno_matrix)
missingSNPSs<-rowSums(is.na(geno_matrix))

hist(missingSNPSs, breaks = 100, col = "lightblue", main = "Histogram of missing SNPs", xlab = "Number of missing SNPs")
# filter by missing data
geno_matrix <- geno_matrix[rowSums(is.na(geno_matrix)) < 0.2 * ncol(geno_matrix),]
dim(geno_matrix)

# filter by allele frequency > 0.05 and < 0.95
geno_matrix <- geno_matrix[rowMeans(geno_matrix)/2 > 0.05 & rowMeans(geno_matrix)/2 < 0.95, ]

# show the number of SNPs and individuals
dim(geno_matrix)
geno_matrix[1, ]
