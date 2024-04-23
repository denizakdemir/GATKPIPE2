setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts")
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
source("http://zzlab.net/GAPIT/gapit_functions.txt")  # GAPIT functions

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
dim(mapping13_thinned)
dim(GenoMat13_thinned)
GenoMat13_thinned[,ncol(GenoMat13_thinned)]
  tmp <- capture.output({
   GAPIT(
      Y = BLUPS,
      GD = GenoMat13_thinned,
      GM = mapping13_thinned,
      KI = Kmat,
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
