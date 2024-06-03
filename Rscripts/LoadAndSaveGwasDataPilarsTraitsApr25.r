# Load necessary packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(
  tidyverse, data.table, ggplot2, dplyr, stringr, 
  genetics, scatterplot3d, snpStats, rtracklayer, 
  GenomicRanges, IRanges, AssocTests, BiocManager
)

# Install Bioconductor packages if not already installed
#BiocManager::install(c("snpStats", "rtracklayer", "GenomicRanges", "IRanges"))
pacman::p_load_gh("SFUStatgen/LDheatmap")  # For LDheatmap from GitHub
source("gapit_functions.txt")  # GAPIT functions

# Read numeric genotypes file
geno <- fread("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/6_tassel_analysis/numeric_allele_counts_filtered.txt", skip = 1, header = TRUE, sep = "\t", na.strings = c("", "NA",0.5,"0.5"))
geno <- as.data.frame(geno)
# find all genotypes that have _sorted in their names

namesInGeno<- geno[,1]
has_sorted<-grepl("_sorted", namesInGeno)

geno<-geno[!has_sorted,]
dim(geno)


geno[1:5,1:5]
dim(geno)
gc()
# Remove duplicate Genotypes
colnames(geno)[1]<-"Taxa"
namesGeno<-geno$Taxa

# remove the last characters from the column names
namesGeno<-substr(namesGeno, 1, nchar(namesGeno)-13)

#see which column names are duplicated
duplicatednames<-duplicated(namesGeno)

# now for each duplicated column name, find the index of the first occurence of the column name and the sec
for (i in 1:length(duplicatednames)) {
    print(i)
  if(duplicatednames[i]) {
    firstOccurence<-which(namesGeno==namesGeno[i])[1]
    secondOccurence<-which(namesGeno==namesGeno[i])[2]
    # now take the mean of the two columns
    geno[firstOccurence,2:ncol(geno)]<-colMeans(geno[c(firstOccurence, secondOccurence),2:ncol(geno)], na.rm=TRUE)
    # remove the second column
    }
  }

sum(duplicatednames)
geno$Taxa[duplicatednames]


geno<-geno[!duplicatednames,]
rownames(geno)<-geno$Taxa
geno<-geno[,-1]
geno<-as.matrix(geno)*2-1
rownames(geno)

genoImp<-rrBLUP::A.mat(geno, impute.method = "mean", return.imputed = TRUE, min.MAF = .1, max.missing = .3)$imputed

gc()

svdG<-svd(scale(genoImp, scale=FALSE, center=TRUE), nu=5,nv=5)
pcs<-scale(genoImp, scale=FALSE, center=TRUE)%*%svdG$v
pcs<-cbind(data.frame(Taxa=rownames(pcs)), pcs)
colnames(pcs)[2:6]<-paste("PC",1:5, sep="")

pcs$Year<-ifelse(grepl("EKDN2300427", pcs$Taxa), "2023", "2022")

filtered_pcs<- pcs[pcs$Year=="2022",]

## ----match_genotype_phenotype---------------------------------------------------------------
# one of the names have an extra "combined_" in it, remove it
filtered_pcs$TaxaShort<-gsub("combined_", "", filtered_pcs$Taxa)
filtered_pcs$TaxaShort

filtered_pcs$TaxaShort<-substr(filtered_pcs$Taxa, 1, 4)

# remove _ and anything that comes after that from the TaxaShort column
filtered_pcs$TaxaShort<-gsub("_.*", "", filtered_pcs$TaxaShort)

# sort the data as S1, S2,...
filtered_pcs<-filtered_pcs[order(as.numeric(gsub("S", "", filtered_pcs$TaxaShort))),]



# read the file
isolateslist<-readxl::read_xlsx("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/GenotypeTables/isolated\ list\ 2021-2022\ for\ sequencing.xlsx")


# have isolate and Name columns in the isolateslist dataframe
# the name column has the S1, S2, ...S100 names
# now by matching the Name column with the TaxaShort column of the filtered_pcs dataframe, get the isolate column of the isolateslist dataframe

filtered_pcs$Isolate<-isolateslist$Isolate[match(filtered_pcs$TaxaShort, isolateslist$Name)]






# Load phenotypes and ensure matching with genotypes
phenotypes<-read.csv("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/PhenoData/20240423_R3E3_R1&R2_GWAS_astringent&loose.csv", header=TRUE)

filtered_pcs$Isolate<-gsub(",", ".", filtered_pcs$Isolate)

phenotypes$Zt_Strain[1:10]
filtered_pcs$Isolate[1:10]
###
phenotypes$Zt_Strain <- gsub("-", ".", phenotypes$Zt_Strain)



# Correct isolate names to match
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_Conil_Fer_L1"]<-"22_ConilFer_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_Jerez Val_L1"]<-"22_JerezVal_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EcijaSecOrt _L1"]<-"22_EcijaSecOrt_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2" ## Check this one
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EscIca_L2"]<-"22_EscIca_L2"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EcijaSecSha_L1"]<-"22_EcijaSecSah_L1"



## -------------------------------------------------------------------------------------------
# put that name in the phenotype file in a new column
for (i in 1:nrow(phenotypes)){
  for (j in 1:nrow(filtered_pcs)){
    if (grepl(filtered_pcs$Isolate[j], phenotypes$Zt_Strain[i])){
      phenotypes$Isolate[i]<-filtered_pcs$Isolate[j]
    }
  }
}
table(phenotypes$Isolate)

head(phenotypes)
PhenoGWAS<-phenotypes[,c("Isolate", "BReplicate", "TReplicate", "Endophyte","Growth_endophyte", "Area_inhibition", "Radius_inhibition", "Ratio_inhibition", "Percentage_inhibition")]


# read snp_genotypeSummary3.txt data 
SNPSummary<-read.table("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/6_tassel_analysis/snp_genotypeSummary3.txt", header = TRUE, sep = "\t", as.is = TRUE, skip=0)
SNPSummary[1:5,]


# Now a map file 
SNPMAP<- data.frame(
  rs = SNPSummary$Site.Name,
  chrom = SNPSummary$Chromosome,
  pos = SNPSummary$Physical.Position
)

SNPMAP<-SNPMAP[SNPMAP$rs%in%colnames(genoImp),]
dim(SNPMAP)

gc()

genoImp <- data.frame(marker=SNPMAP$rs,chrom=SNPMAP$chrom,pos=SNPMAP$pos,t(genoImp))

gc()
genoImp[1:5,1:5]

colnames(PhenoGWAS)[1]<-"line"
PhenoGWAS$BReplicate<-as.factor(PhenoGWAS$BReplicate)
PhenoGWAS$TReplicate<-as.factor(PhenoGWAS$TReplicate)
PhenoGWAS$line<-as.factor(PhenoGWAS$line)
PhenoGWAS$Endophyte <-as.factor(PhenoGWAS$Endophyte)



PhenoGWASR3E3<-PhenoGWAS[PhenoGWAS$Endophyte=="R3E3",]
PhenoGWASR3E3<-na.omit(PhenoGWASR3E3)

PhenoGWASR3E5<-PhenoGWAS[PhenoGWAS$Endophyte=="R3E5",]
PhenoGWASR3E5<-na.omit(PhenoGWASR3E5)

GenoMat<-genoImp[,-c(1:3)]
GenoMat<-t(GenoMat)+1
GenoMat<-as.data.frame(GenoMat)

# rownames as the first column as TAxa
GenoMat$Taxa<-rownames(GenoMat)
rownames(GenoMat)<-NULL
# make taxa the first column
GenoMat<-GenoMat[,c(ncol(GenoMat), 1:(ncol(GenoMat)-1))]

gc()
GenoMat$Taxa<-gsub("\\.","-", GenoMat$Taxa)
GenoMat<-GenoMat[GenoMat$Taxa%in%filtered_pcs$Taxa,]
GenoMat<-GenoMat[match(filtered_pcs$Taxa,GenoMat$Taxa),]

GenoMat$Taxa<-filtered_pcs$Isolate

GenoMat[1:5,1:5]

PhenoGWASR3E3$line<-factor(PhenoGWASR3E3$line, levels=GenoMat$Taxa)
PhenoGWASR3E5$line<-factor(PhenoGWASR3E5$line, levels=GenoMat$Taxa)
PhenoGWAS$line<-factor(PhenoGWAS$line, levels=GenoMat$Taxa)

table(SNPMAP$chrom)
# filter to first 13 chromosomes
SNPMAP<-SNPMAP[SNPMAP$chrom<=13,]
GenoMat<-GenoMat[,c("Taxa", SNPMAP$rs)]


# Save all data needed for GWAS
new_dir <- "DataforPillarApr25"
if (!dir.exists(new_dir)) dir.create(new_dir, recursive = TRUE)
save(PhenoGWASR3E3, GenoMat, SNPMAP, file = paste0(new_dir, "/DataforPillar.RData"))

rm(list=ls())
gc()
###############################################

setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts")

# Load necessary libraries
if (!require(pacman)) install.packages("pacman")
pacman::p_load(
  tidyverse, data.table, ggplot2, dplyr, stringr, 
  genetics, scatterplot3d, snpStats, rtracklayer, 
  GenomicRanges, IRanges, AssocTests, BiocManager, lme4
)

# Load external GAPIT functions
source("gapit_functions.txt")

# Load data from the previous preparation stage
load("DataforPillarApr25/DataforPillar.RData")

# Ensure phenotype data and genotype data Taxa columns are aligned
GenoMat$Taxa <- as.character(GenoMat$Taxa)
PhenoGWASR3E3$line <- as.character(PhenoGWASR3E3$line)




# Assuming GenoMat and PhenoGWAS are correctly matched and prepped
# Perform SNP thinning based on correlation
library(pbapply)
library(dplyr)

GenoMat[1:5, 1:5]
GenoImp<-GenoMat[,-1]-1
rownames(GenoImp)<-GenoMat$Taxa
gc()
GenoImp<-rrBLUP::A.mat(GenoImp, impute.method = "mean", return.imputed = TRUE, min.MAF = .1, max.missing = .3)$imputed+1
# read snp_genotypeSummary3.txt data 
SNPSummary<-read.table("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/6_tassel_analysis/snp_genotypeSummary3.txt", header = TRUE, sep = "\t", as.is = TRUE, skip=0)
SNPSummary[1:5,]
SNPSummary<-SNPSummary[SNPSummary$Site.Name %in% colnames(GenoImp),]
Kmat<-rrBLUP::A.mat(GenoImp-1, return.imputed = FALSE)
Kmat[1:5,1:5]
# Now a map file 
SNPMAP<- data.frame(
  Name = SNPSummary$Site.Name,
  Chromosome = SNPSummary$Chromosome,
  Position = SNPSummary$Physical.Position
)
library(dplyr)
library(pbapply)  # For progress bars
dim(SNPMAP)

find_highest_corr_snp <- function(df, GenoMat, distance) {
  # Sort the dataframe by chromosome and position
  df <- df %>%
    arrange(Chromosome, Position)

  # Split the dataframe by chromosome
  chromosome_list <- split(df, df$Chromosome)

  # Function to process each chromosome
  process_chromosome <- function(chromo_df) {
    # Create regions within the chromosome
    chromo_df <- chromo_df %>%
      mutate(region = floor(Position / distance))

    # Split the chromosome dataframe by region
    region_list <- split(chromo_df, chromo_df$region)

    # Function to process each region
    process_region <- function(region_df) {
      # Filter SNPs that are present in the GenoMat column names
      valid_snp_names <- region_df$Name[region_df$Name %in% colnames(GenoMat)]
      if (length(valid_snp_names) < 2) {
        return(region_df[1, ])
      }
      region_data <- GenoMat[, valid_snp_names, drop = FALSE]
      corr_matrix <- cor(region_data, use = "pairwise.complete.obs")
      if (all(is.na(corr_matrix))) {
        return(region_df[1, ])
      }
      avg_corr <- colMeans(corr_matrix, na.rm = TRUE)
      highest_corr_snp <- names(which.max(avg_corr))
      return(region_df %>% filter(Name == highest_corr_snp))
    }

    # Apply the process_region to each region and combine results
    highest_corr_snps_list <- pbapply::pblapply(region_list, process_region)
    highest_corr_snps <- bind_rows(highest_corr_snps_list)
    return(highest_corr_snps)
  }

  # Apply the process_chromosome to each chromosome and combine results
  all_highest_corr_snps <- lapply(chromosome_list, process_chromosome)
  final_result <- bind_rows(all_highest_corr_snps)
  return(final_result)
}


# Assuming SNPMAP and GenoMat are already loaded
# Apply thinning function to your mapping dataframe for the first 13 chromosomes
mapping13 <- SNPMAP %>% 
  filter(Chromosome <= 13)

# Assuming GenoMat is appropriately indexed or adjusted for the data subset
GenoMat13 <- GenoImp[, mapping13$Name]

# Applying the function with a progress bar
mapping_highest_corr <- find_highest_corr_snp(mapping13, GenoMat13, 200) # Adjust the distance as neededGenoMat13_thinned <- GenoMat13[, c("Taxa", mapping13_thinned$Name)]

mapping_highest_corr<-na.omit(mapping_highest_corr)

dim(mapping_highest_corr)

# Filter the genotype matrix based on the selected SNPs
GenoMat13_thinned <- GenoMat13[, mapping_highest_corr$Name]
dim(GenoMat13_thinned)
# filter mapping13
mapping13_thinned <- mapping13[mapping13$Name %in% colnames(GenoMat13_thinned),]
dim(mapping13_thinned)
##########################save the data neede for GWAS
save(PhenoGWASR3E3, GenoMat13_thinned,mapping13_thinned,file="DataforPillarApr25/DataforPillarGWAS.RData")

# calculate the Kmat
Kmat<-rrBLUP::A.mat(GenoMat13_thinned[, -1]-1, min.MAF = .1, max.missing = .3)

# save Kmat, GenoMat13_thinned, mapping13_thinned 
save(Kmat, GenoMat13_thinned, mapping13_thinned, file = "DataforPillarApr25/GenoDataSeptoria100.RData")

