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
source("gapit_functions.txt")  # GAPIT functions

load("DataforPillarApr25_filtered_loose/DataforPillarGWAS.RData")
length(unique(PhenoGWASR3E3$line))
GenoMat13_thinned[1:5,1:5]
Kmat<-rrBLUP::A.mat(GenoMat13_thinned-1)


intersectNames<-intersect(PhenoGWASR3E3$line, rownames(Kmat))
Kmat<-Kmat[intersectNames,intersectNames]
GenoMat13_thinned<-GenoMat13_thinned[rownames(GenoMat13_thinned) %in% intersectNames,]
PhenoGWASR3E3<-PhenoGWASR3E3[PhenoGWASR3E3$line %in% intersectNames,]
PhenoGWASR3E3$line<-factor(PhenoGWASR3E3$line, levels=rownames(Kmat))
# Calculate BLUPs for specified traits using lmer from the lme4 package
calculate_blups <- function(data, trait) {
  model <- lmer(formula(paste(trait, "~ BReplicate + TReplicate+ (1|line)")), data = data)
  return(ranef(model)$line[,1])
}
library(lme4)



traits <- c("Growth_endophyte","Radius_inhibition", "Area_inhibition", "Ratio_inhibition", "Percentage_inhibition")
blups <- sapply(traits, calculate_blups, data = PhenoGWASR3E3)
dim(blups)


GenoMat13_thinned[1:5,1:5]
# Put these BLUPs in a dataframe
BLUPS <- data.frame(Isolate = rownames(GenoMat13_thinned), blups)

# Run GAPIT for GWAS analysis
models <- c("FarmCPU", "Blink")

# Save all data needed for GWAS
new_dir <- "GWAS_resultsForPilarsTraitsApr25_WithFIlter_loose"
if (!dir.exists(new_dir)) dir.create(new_dir, recursive = TRUE)
GenoMat13_thinned<-as.data.frame(GenoMat13_thinned)
GenoMat13_thinned$Taxa<-rownames(GenoMat13_thinned)
# put taxa as the first column
GenoMat13_thinned<-GenoMat13_thinned[,c(ncol(GenoMat13_thinned), 1:(ncol(GenoMat13_thinned)-1))]

setwd("GWAS_resultsForPilarsTraitsApr25_WithFIlter_loose")
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

# save the data needed for GWAS
save(BLUPS, GenoMat13_thinned, mapping13_thinned, Kmat, file="GWASdataBefGapit.RData")

###############
getwd()
#setwd("Rscripts/GWAS_resultsForPilarsTraitsApr25_WithFIlter_loose")
load("GWASdataBefGapit.RData")
source("../gapit_functions.txt")  # GAPIT functions
models=c("FarmCPU", "Blink")


# Assuming you start in the root directory where you want to create subdirectories
initial_dir <- getwd()

# Array to hold the number of principal components to be used in each analysis
pc_values <- 0:5

# Model types to be used
models <- c("FarmCPU", "Blink")

# Loop through each principal component setting
for (pc in pc_values) {
    # Create a directory name based on the principal component
    dir_name <- paste("GWAS_resultsForPilarsTraitsApr25_", pc, "PC", sep="")
    
    # Check if the directory exists, if not create it
    if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE)
    }
    
    # Set the working directory to the new directory
    setwd(file.path(initial_dir, dir_name))
    
    # Run the GAPIT function
    tmp <- capture.output({
        GAPIT(
            Y = BLUPS,
            GD = GenoMat13_thinned,
            GM = mapping13_thinned,
            KI = Kmat,
            CV = NULL,
            PCA.total = pc,
            model = models,
            file.output = TRUE,
            cutOff = 0.1
        )
    })
    
    # Go back to the initial directory
    setwd(initial_dir)
}

# Optionally print a message when the loop is complete
cat("GAPIT analysis complete for all principal component settings.\n")




############################################################################################################






# # Load the NAM package
# if (!require(NAM)) install.packages("NAM", repos = "https://cran.r-project.org")
# library(NAM)

# # https://rdrr.io/cran/NAM/man/gwas.html

# MAP<-mapping13_thinned

# gwasNAM<-gwas(y=BLUPS$Radius_inhibition,
#               gen=as.matrix(GenoMat13_thinned[,-1]), MAP = mapping13_thinned)


# str(gwasNAM)
# plot(gwasNAM)
# head(gwasNAM$PolyTest)
# hist(gwasNAM$PolyTest$pval)
# pvals<-data.frame(pval=gwasNAM$PolyTest$pval)
# pvals$SNP<-rownames(pvals)
# pvals$pval[pvals$pval<0]<-1
# pvals$pval[pvals$pval>1]<-1

# png("GWAS_pvalues.png", width=800, height=600)
# plot(-log10(pvals$pval), type="h", xlab="SNP", ylab="p-value", main="GWAS p-values")
# abline(h=-log10(0.05), col="red", lty=2)
# dev.off()



