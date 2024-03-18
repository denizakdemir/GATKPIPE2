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
geno[1:5, 1:8]
max(na.omit(geno[3, ]))
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
sum(missingness <= 30)
geno <- geno[missingness <= 30, ]
dim(geno)
# Assuming your data is in a dataframe named 'geno'
# Replace 'geno' with the actual name of your dataframe

library(dplyr)
library(tidyr)
library(progress)

# Function to count alleles
allele_counts<-apply(geno, 1, function(snp) {
  mean(snp, na.rm=TRUE)
})

hist(allele_counts, breaks=100)
# if an allele count is more than 0.5,it needs to be flipped
for (i in 1:nrow(geno)) {
  if(mean(geno[i,], na.rm=TRUE) >= 0.5) {
    geno[i,] <- 1 - geno[i,]
  }
}


# Function to count alleles
allele_counts<-apply(geno, 1, function(snp) {
  mean(snp, na.rm=TRUE)
})

hist(allele_counts, breaks=100)

# remove alleles that have a frequency of 0.05 or less
geno <- geno[allele_counts > 0.05, ]
dim(geno)

# Duplicate genotypes.
# Take the consensus genotype for each SNP (mean of the two genotypes if both are present, otherwise the one that is present)
#  the duplicate names differ by one character, so use partial matching
# use pmatch(x, table, nomatch = NA_integer_, duplicates.ok = FALSE)

# Remove duplicate Genotypes
colnamesGeno<-colnames(geno)
colnamesGeno<-gsub("X", "", colnamesGeno)

# remove the last character from the column names
colnamesGeno<-substr(colnamesGeno, 1, nchar(colnamesGeno)-13)

#see which column names are duplicated
duplicatedColnames<-duplicated(colnamesGeno)

# now for each duplicated column name, find the index of the first occurence of the column name and the sec
for (i in 1:length(duplicatedColnames)) {
  if(duplicatedColnames[i]) {
    firstOccurence<-which(colnamesGeno==colnamesGeno[i])[1]
    secondOccurence<-which(colnamesGeno==colnamesGeno[i])[2]
    # now take the mean of the two columns
    geno[,firstOccurence]<-rowMeans(geno[,c(firstOccurence, secondOccurence)], na.rm=TRUE)
    # remove the second column
    }

  }

dim(geno)
install.packages("stringdist")
# now use the following example to go over each column and see if the remaning 
# stringdist::amatch('fu', c('foo','bar'), maxDist=2)

for (colname in colnamesGeno) {
  if(length(stringdist::amatch(colname, setdiff(colnamesGeno,colname), maxDist=2))>1) {
    print(colname)
  }
}

# impute each snp with the mean of the snp (snps are in rows)
for (i in 1:nrow(geno)) {
  geno[i,]<-ifelse(is.na(geno[i,]), mean(geno[i,], na.rm=TRUE), geno[i,])
}

# calculate GRM
grm<-rrBLUP::A.mat(t(geno)*2-1)
grm[1:5,1:5]

heatmap(grm)


svdG<-svd(scale(t(geno), scale=FALSE, center=TRUE), nu=5,nv=5)

PCGeno<-scale(t(geno), scale=FALSE, center=TRUE)%*%svdG$v

PCGeno<-cbind(data.frame(Taxa=rownames(PCGeno)), PCGeno)
colnames(PCGeno)[2:6]<-paste("PC",1:5, sep="")

pcs<-PCGeno

# View the first few rows of the data
head(pcs)

plot(pcs$PC1, pcs$PC2)
library(dplyr)

# Assuming genotype_df is your dataframe
filtered_pcs <- pcs %>%
  filter(!grepl("EKDN2300427", Taxa)) # Remove rows where 'Taxa' contains 'EKDN2300427'

dim(filtered_pcs)
excludedTaxa<-c("S101_EKDN220047288-1A_HMVJWDSX5_L2")

filtered_pcs<-filtered_pcs[!filtered_pcs$Taxa%in%excludedTaxa,]
dim(filtered_pcs)

# Remove duplicate Genotypes
colnamesGeno<-rownames(filtered_pcs)

# remove the last character from the column names
colnamesGeno<-substr(colnamesGeno, 1, nchar(colnamesGeno)-13)

#see which column names are duplicated
duplicatedColnames<-duplicated(colnamesGeno)

filtered_pcs<-filtered_pcs[!duplicatedColnames,]
dim(filtered_pcs)

plot(filtered_pcs$PC1, filtered_pcs$PC2)


# read this file "/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE/GenotypeTables/isolated\ list\ 2021-2022\ for\ sequencing.xlsx"

# read the file
isolatedlist<-readxl::read_xlsx("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE/GenotypeTables/isolated\ list\ 2021-2022\ for\ sequencing.xlsx")


# have isolate and Name columns in the isolatedlist dataframe
# the name column has the S1, S2, ...S100 names
# now by matching the Name column with the TaxaShort column of the filtered_pcs_nonduplicated dataframe, get the isolate column of the isolatedlist dataframe

filtered_pcs$Isolate<-isolatedlist$Isolate[match(filtered_pcs$TaxaShort, isolatedlist$Name)]


# now plot the first two PCs, color by the PC3.

pcaplot<-ggplot(filtered_pcs, aes(x=PC1, y=PC2, color=PC3)) + geom_point() + theme_minimal()+
    geom_text(aes(label=Isolate),hjust=0, vjust=0)


# save the image as PCAPLOT.png
ggsave("PCAPLOT.png", pcaplot, dpi=600)

