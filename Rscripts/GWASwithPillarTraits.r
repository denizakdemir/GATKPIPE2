setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts")

# Load necessary libraries
if (!require(pacman)) install.packages("pacman")
pacman::p_load(
  tidyverse, data.table, ggplot2, dplyr, stringr, 
  genetics, scatterplot3d, snpStats, rtracklayer, 
  GenomicRanges, IRanges, AssocTests, BiocManager, lme4
)

# Load external GAPIT functions
source("http://zzlab.net/GAPIT/gapit_functions.txt")

# Load data from the previous preparation stage
load("DataforPillar/DataforPillar.RData")

# Ensure phenotype data and genotype data Taxa columns are aligned
GenoMat$Taxa <- as.character(GenoMat$Taxa)
PhenoGWAS$line <- as.character(PhenoGWAS$line)
PhenoGWASR3E3$line <- as.character(PhenoGWASR3E3$line)
PhenoGWASR3E5$line <- as.character(PhenoGWASR3E5$line)

intersect(PhenoGWAS$line, GenoMat$Taxa)
setdiff(PhenoGWAS$line, GenoMat$Taxa)
setdiff(GenoMat$Taxa,PhenoGWAS$line)


GenoMat$Taxa[GenoMat$Taxa=="22_Conil_Fer_L1"]<-"22_ConilFer_L1"
GenoMat$Taxa[GenoMat$Taxa=="22_Jerez Val_L1"]<-"22_JerezVal_L1"
GenoMat$Taxa[GenoMat$Taxa=="22_EcijaSecOrt _L1"]<-"22_EcijaSecOrt_L1"
GenoMat$Taxa[GenoMat$Taxa=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2" ## Check this one
GenoMat$Taxa[GenoMat$Taxa=="22_EscIca_L2"]<-"22_EscIca_L2"
GenoMat$Taxa[GenoMat$Taxa=="22_EcijaSecSha_L1"]<-"22_EcijaSecSah_L1"

length(unique(PhenoGWAS$line))
intersect(PhenoGWAS$line, GenoMat$Taxa)
setdiff(PhenoGWAS$line, GenoMat$Taxa)
setdiff(GenoMat$Taxa,PhenoGWAS$line)


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



library(dplyr)
library(pbapply)
library(stats)

# find_region_stats <- function(df, GenoMat, distance) {
#   # Sort the dataframe by chromosome and position
#   df <- df %>%
#     arrange(Chromosome, Position)

#   # Split the dataframe by chromosome
#   chromosome_list <- split(df, df$Chromosome)

#   # Function to process each chromosome
#   process_chromosome <- function(chromo_df, chromo_key) {
#     # Create regions within the chromosome
#     chromo_df <- chromo_df %>%
#       mutate(region = floor(Position / distance))

#     # Split the chromosome dataframe by region
#     region_list <- split(chromo_df, chromo_df$region)

#     # Function to process each region
#     process_region <- function(region_df, region_key) {
#       # Calculate start and end positions
#       start_pos <- min(region_df$Position)
#       end_pos <- max(region_df$Position)

#       # Filter SNPs that are present in the GenoMat column names
#       valid_snp_names <- region_df$Name[region_df$Name %in% colnames(GenoMat)]
#       if (length(valid_snp_names) < 1) {
#         return(data.frame(Chromosome = chromo_key, Region = region_key, Start = start_pos, End = end_pos, PC1 = NA))
#       }

#       # Get genotype data for these SNPs
#       region_data <- GenoMat[, valid_snp_names, drop = FALSE]

#       # Perform PCA on the genotype data
#       pca <- prcomp(region_data, scale. = FALSE)
#       pc1 <- as.data.frame(t(pca$x[, 1]))  # First principal component across genotypes
#       colnames(pc1) <- rownames(GenoMat)
#       print(dim(pc1))
#       out<-data.frame(Chromosome = chromo_key, Region = region_key, Start = start_pos, End = end_pos)
#       out<-cbind(out, pc1)
#       # Return data frame with chromosome, region, start, end, and the first PC
#       # print(out$Chromosome)
#       # print(out$Region)
#       # print(out$Start)
#       # print(out$End)
#       # print(out[,1:10])
#       return(out)
#     }

#     # Apply the process_region to each region and combine results
#     region_stats_list <- pbapply::pblapply(region_list, process_region, region_key = names(region_list))
#     region_stats <- bind_rows(region_stats_list)
#     return(region_stats)
#   }

#   # Apply the process_chromosome to each chromosome and combine results
#   all_region_stats <- mapply(process_chromosome, chromosome_list, names(chromosome_list), SIMPLIFY = FALSE)
#   final_result <- bind_rows(all_region_stats)
#   return(final_result)
# }


# Assuming SNPMAP and GenoMat are already loaded
# Apply thinning function to your mapping dataframe for the first 13 chromosomes
mapping13 <- SNPMAP %>% 
  filter(Chromosome <= 13)

# Assuming GenoMat is appropriately indexed or adjusted for the data subset
GenoMat13 <- GenoImp[, mapping13$Name]

# Applying the function with a progress bar
mapping_highest_corr <- find_highest_corr_snp(mapping13, GenoMat13, 200) # Adjust the distance as neededGenoMat13_thinned <- GenoMat13[, c("Taxa", mapping13_thinned$Name)]
#PCAbyRegion <- find_region_stats(mapping13, GenoMat13, 500) # Adjust the distance as needed

#str(PCAbyRegion)

#PCAbyRegionMat<-as.data.frame(PCAbyRegion)


dim(mapping_highest_corr)
head(mapping_highest_corr)
mapping_highest_corr<-na.omit(mapping_highest_corr)
dim(mapping_highest_corr)

# There are a couple of snps I want to keep in add them to mapping_highest_corr
# they are here: Rscripts/GWASforPilarTraitsMarch27thThinnedSnps/GAPIT.Association.Filter_GWAS_results.csv

# read the file
GWASresults<-read.table("GWASforPilarTraitsMarch27thThinnedSnps/GAPIT.Association.Filter_GWAS_results.csv", header = TRUE, sep = ",", as.is = TRUE, skip=0)
GWASresults$SNP
GWASresults$Chr
GWASresults$Position
# get the unique snps
SNPStoKEEP<-unique(GWASresults$SNP)

MaptoAdd<-mapping13[mapping13$Name %in% SNPStoKEEP,]
mapping_highest_corr<-rbind(mapping_highest_corr, MaptoAdd)
# sort the mapping_highest_corr, chr and pos
mapping_highest_corr<-mapping_highest_corr[order(mapping_highest_corr$Chromosome, mapping_highest_corr$Position),]



# Filter the genotype matrix based on the selected SNPs
GenoMat13_thinned <- GenoMat13[, mapping_highest_corr$Name]
dim(GenoMat13_thinned)
# filter mapping13
mapping13_thinned <- mapping13[mapping13$Name %in% colnames(GenoMat13_thinned),]

##########################save the data neede for GWAS
save(PhenoGWAS,PhenoGWASR3E3,PhenoGWASR3E5,GenoMat13_thinned,mapping_highest_corr,file="DataforPillar/DataforPillarGWAS.RData")


load("DataforPillar/DataforPillarGWAS.RData")

# Calculate BLUPs for specified traits using lmer from the lme4 package
calculate_blups <- function(data, trait) {
  model <- lmer(formula(paste(trait, "~ BReplicate + TReplicate+ (1|line)")), data = data)
  return(ranef(model)$line[,1])
}
library(lme4)
traits <- c("Radius_inhibition", "Area_inhibition", "Ratio_inhibition", "Percentage_inhibition")
blups <- sapply(traits, calculate_blups, data = PhenoGWASR3E3)
dim(blups)
GenoMat13_thinned[1:5,1:5]
Kmat<-rrBLUP::A.mat(GenoMat13_thinned[, -1]-1, min.MAF = .1, max.missing = .3)
# Put these BLUPs in a dataframe
BLUPS <- data.frame(Isolate = rownames(GenoMat13_thinned), blups)

# Run GAPIT for GWAS analysis
models <- c("FarmCPU", "Blink")

# Save all data needed for GWAS
new_dir <- "GWAS_resultsForPilarsTraits"
if (!dir.exists(new_dir)) dir.create(new_dir, recursive = TRUE)
GenoMat13_thinned<-as.data.frame(GenoMat13_thinned)
GenoMat13_thinned$Taxa<-rownames(GenoMat13_thinned)
# put taxa as the first column
GenoMat13_thinned<-GenoMat13_thinned[,c(ncol(GenoMat13_thinned), 1:(ncol(GenoMat13_thinned)-1))]

setwd("GWAS_resultsForPilarsTraits")
dim(BLUPS)
dim(GenoMat13_thinned)
dim(mapping13_thinned)
dimnames(Kmat)
dimnames(GenoMat13_thinned)
GenoMat13_thinned$Taxa<-NULL
GenoMat13_thinned$Taxa<-rownames(GenoMat13_thinned)
GenoMat13_thinned<-GenoMat13_thinned[,c(ncol(GenoMat13_thinned), 1:(ncol(GenoMat13_thinned)-1))]
GenoMat13_thinned[1:5,1:5]

Kmat<-as.data.frame(Kmat)
Kmat$Taxa<-rownames(Kmat)
Kmat<-Kmat[,c(ncol(Kmat), 1:(ncol(Kmat)-1))]
dim(Kmat)
rownames(Kmat)<-NULL
colnames(Kmat)<-NULL

Kmat[1:5,1:5]
GenoMat13_thinned[1:5,1:5]

  tmp <- capture.output({
   GAPIT(
      Y = BLUPS,
      GD = GenoMat13_thinned,
      GM = mapping13_thinned,
      KI = NULL,
      CV = NULL,
      PCA.total = 0,
      model = models,
      file.output = TRUE,
      cutOff = 0.1
    )
  })




# ########################

# sampleTrain<-sample(BLUPS$Taxa, nrow(BLUPS)*0.8)

# # Clean up and display summary
# rm(tmp)
# gc()
# print(results)



# # Load the NAM package
# if (!require(NAM)) install.packages("NAM", repos = "https://cran.r-project.org")
# library(NAM)

# # https://rdrr.io/cran/NAM/man/gwas.html



# gwasNAM<-gwas(y=BLUPS$Radius_inhibition,
#               gen=as.matrix(GenoMat13_thinned[,-1]))
# )

# str(gwasNAM)
# head(gwasNAM$PolyTest)
# pvals<-data.frame(pval=gwasNAM$PolyTest$pval)
# pvals$SNP<-rownames(pvals)
# pvals$pval[pvals$pval<0]<-1

# png("GWAS_pvalues.png", width=800, height=600)
# plot(-log10(pvals$pval), type="h", xlab="SNP", ylab="p-value", main="GWAS p-values")
# abline(h=-log10(0.05), col="red", lty=2)
# dev.off()
