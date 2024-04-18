






# Load the pacman package for efficient package management
if (!require(pacman)) install.packages("pacman")

# Use pacman to load or install and load packages
pacman::p_load(
  tidyverse,       # For data manipulation and visualization
  data.table,      # For fast data manipulation
  ggplot2,         # For data visualization
  dplyr,           # For data manipulation
  stringr,         # For string operations
  genetics,        # For genetic data analysis
  scatterplot3d,   # For 3D scatter plots
  snpStats,        # For SNP statistics
  rtracklayer,     # For reading/writing genomic data
  GenomicRanges,   # For representing genomic ranges and operations on them
  IRanges,         # For integer ranges operations
  AssocTests       # For association tests
)

# Bioconductor package management with BiocManager
pacman::p_load(BiocManager)
pacsBIO <- c("snpStats", "rtracklayer", "GenomicRanges", "IRanges")

# Install Bioconductor packages if not already installed, using BiocManager
pacman::p_load_gh("SFUStatgen/LDheatmap") # For loading LDheatmap from GitHub directly

# Check and install (if necessary) Bioconductor packages
# for (pkg in pacsBIO) {
# 
#     BiocManager::install(pkg)
#   
# }

# If the external GAPIT functions are not available as a GitHub repository, source them directly from the web
source("http://zzlab.net/GAPIT/gapit_functions.txt")

# create a directory to save the report files if the directory does not exist
if (!file.exists("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/reportfiles")){
dir.create("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/reportfiles", showWarnings = FALSE)
}


## ----load_genomic_data, eval=FALSE----------------------------------------------------------
## # read numeric genotypes file:
geno <- fread("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/6_tassel_analysis/numeric_allele_counts_filtered.txt", skip = 1, header = TRUE, sep = "\t", na.strings = c("", "NA",0.5,"0.5"))
geno<-as.data.frame(geno)

geno[,-1][geno[,-1] >1] <- NA

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

dim(geno)
sum(duplicatednames)
geno$Taxa[duplicatednames]
geno<-geno[!duplicatednames,]
dim(geno)
rownames(geno)<-geno$Taxa
geno<-geno[,-1]
dim(geno)
geno<-as.matrix(geno)*2-1


## ----imputation, eval=FALSE-----------------------------------------------------------------
geno<-rrBLUP::A.mat(geno, impute.method = "mean", return.imputed = TRUE, min.MAF = .1, max.missing = .3)$imputed
dim(geno)
save(geno, file="/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/reportfiles/geno.RData")
gc()

KMat<-rrBLUP::A.mat(geno, impute.method = "mean", return.imputed = FALSE, min.MAF = .1, max.missing = .3)

save(KMat, file="/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/reportfiles/KMat.RData")
gc()


## ----pca------------------------------------------------------------------------------------
load("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/reportfiles/geno.RData")
svdG<-svd(scale(geno, scale=FALSE, center=TRUE), nu=5,nv=5)
pcs<-scale(geno, scale=FALSE, center=TRUE)%*%svdG$v
pcs<-cbind(data.frame(Taxa=rownames(pcs)), pcs)
colnames(pcs)[2:6]<-paste("PC",1:5, sep="")


## ----plot_pca-------------------------------------------------------------------------------

# Create a ggplot object
# Assuming you have the 'ggplot2' library loaded
library(ggplot2)
library(ggrepel)
# Create a ggplot object with adjusted x and y-axis limits
# Assuming pcs is your data frame and it's already loaded
p <- ggplot(pcs, aes(x=PC1, y=PC2, label=Taxa)) +
  geom_point() +
  geom_text_repel(aes(label=Taxa), size=.5, 
                  box.padding = .1, 
                  point.padding = 0.1,
                  segment.color = 'grey50', max.overlaps = 120) +
  xlim(min(pcs$PC1) - 5, max(pcs$PC1) + 5) + # Adjust xlim
  ylim(min(pcs$PC2) - 5, max(pcs$PC2) + 5)   # Adjust ylim





# Convert to plotly object
#p_plotly <- plotly::ggplotly(p)

# Display the plot
#p_plotly


## -------------------------------------------------------------------------------------------
# use  ('Taxa' contains 'EKDN2300427') as color for the points (grepl("EKDN2300427", Taxa)).
# Create a ggplot object with smaller font size for labels

pcs$Year<-ifelse(grepl("EKDN2300427", pcs$Taxa), "2023", "2022")
p <- ggplot(pcs, aes(x=PC1, y=PC2, color=Year)) +
  geom_point()  +
  geom_text_repel(aes(label=Taxa), size=.5, 
                  box.padding = .1, 
                  point.padding = 0.1,
                  segment.color = 'grey50', max.overlaps = 120) +
  xlim(min(pcs$PC1) - 10, max(pcs$PC1) + 10) + # Adjust xlim
  ylim(min(pcs$PC2) - 10, max(pcs$PC2) + 10)   # Adjust ylim


# save the plot
ggsave("pcaplot100withnames.png", p, device="png")



## -------------------------------------------------------------------------------------------
# Assuming genotype_df is your dataframe
filtered_pcs <- pcs %>%
  filter(!grepl("EKDN2300427", Taxa)) # Remove rows where 'Taxa' contains 'EKDN2300427'

p<-ggplot(filtered_pcs, aes(x=PC1, y=PC2, label=Taxa)) + geom_point()  +
  geom_text_repel(aes(label=Taxa), size=.5, 
                  box.padding = .1, 
                  point.padding = 0.1,
                  segment.color = 'grey50', max.overlaps = 120) +
  xlim(min(pcs$PC1) - 10, max(pcs$PC1) + 10) + # Adjust xlim
  ylim(min(pcs$PC2) - 10, max(pcs$PC2) + 10)   # Adjust ylim
p



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
isolatedlist<-readxl::read_xlsx("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/GenotypeTables/isolated\ list\ 2021-2022\ for\ sequencing.xlsx")


# have isolate and Name columns in the isolatedlist dataframe
# the name column has the S1, S2, ...S100 names
# now by matching the Name column with the TaxaShort column of the filtered_pcs dataframe, get the isolate column of the isolatedlist dataframe

filtered_pcs$Isolate<-isolatedlist$Isolate[match(filtered_pcs$TaxaShort, isolatedlist$Name)]


# now plot the first two PCs, color by the PC3.
pcaplot<-ggplot(filtered_pcs, aes(x=PC1, y=PC2, color=PC3)) + geom_point() + theme_minimal() +
  geom_text_repel(aes(label=Taxa), size=.5, 
                  box.padding = .1, 
                  point.padding = 0.1,
                  segment.color = 'grey50', max.overlaps = 120) +
  xlim(min(pcs$PC1) - 5, max(pcs$PC1) + 5) + # Adjust xlim
  ylim(min(pcs$PC2) - 5, max(pcs$PC2) + 5)   # Adjust ylim

# Convert to plotly object

ggsave("pcaplot100withnamesinPheno.png", pcaplot, device="png")



## -------------------------------------------------------------------------------------------
 
# read the file
phenotypes<-read.csv("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/PhenoData/20240323_R3E3&R3E5_R1&R2_AreaRatioPercScalesRemoval_GWAS.csv", header=TRUE)
head(phenotypes)
colnames(phenotypes)


filtered_pcs$Isolate[1:10]

# convert all commas to dots in filtered_pcs$Isolate
filtered_pcs$Isolate<-gsub(",", ".", filtered_pcs$Isolate)

phenotypes$Zt_Strain[1:10]
filtered_pcs$Isolate[1:10]
###
phenotypes$Zt_Strain <- gsub("-", ".", phenotypes$Zt_Strain)

dim(geno)
geno[1:5,1:5]


intersect(phenotypes$Zt_Strain, filtered_pcs$Isolate)
setdiff(filtered_pcs$Isolate, phenotypes$Zt_Strain)
setdiff(phenotypes$Zt_Strain,filtered_pcs$Isolate)

filtered_pcs$Isolate[filtered_pcs$Isolate=="22_Conil_Fer_L1"]<-"22_ConilFer_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_Jerez Val_L1"]<-"22_JerezVal_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EcijaSecOrt _L1"]<-"22_EcijaSecOrt_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2" ## Check this one
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EscIca_L2"]<-"22_EscIca_L2"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EcijaSecSha_L1"]<-"22_EcijaSecSah_L1"

intersect(phenotypes$Zt_Strain, filtered_pcs$Isolate)
setdiff(filtered_pcs$Isolate, phenotypes$Zt_Strain)
setdiff(phenotypes$Zt_Strain,filtered_pcs$Isolate)


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


intersect(colnames(phenotypes), c("Isolate", "BReplicate", "TReplicate", "Endophyte","Area_Endophyte", "Area_inhibition", "Radius_inhibition", "Ratio_inhibition", "Ratio_Scales1", "Ratio_Scales2", "Percentage_inhibition"))
PhenoGWAS<-phenotypes[,c("Isolate", "BReplicate", "TReplicate", "Endophyte","Area_Endophyte", "Area_inhibition", "Radius_inhibition", "Ratio_inhibition", "Ratio_Scales1", "Ratio_Scales2", "Percentage_inhibition")]


# add the first three PCs as fixed effects
#PhenoGWAS<-merge(PhenoGWAS, filtered_pcs[,c("Isolate", "PC1", "PC2", "PC3")], by="Isolate")


table(PhenoGWAS$Endophyte)
table(PhenoGWAS$Endophyte, PhenoGWAS$BReplicate)
table(PhenoGWAS$Endophyte, PhenoGWAS$TReplicate)
table(PhenoGWAS$TReplicate, PhenoGWAS$BReplicate)
PhenoGWAS[1:5,]


## -------------------------------------------------------------------------------------------

# read snp_genotypeSummary3.txt data 
SNPSummary<-read.table("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/6_tassel_analysis/snp_genotypeSummary3.txt", header = TRUE, sep = "\t", as.is = TRUE, skip=0)
SNPSummary[1:5,]


# Now a map file 
SNPMAP<- data.frame(
  rs = SNPSummary$Site.Name,
  chrom = SNPSummary$Chromosome,
  pos = SNPSummary$Physical.Position
)

SNPMAP<-SNPMAP[SNPMAP$rs%in%colnames(geno),]
dim(SNPMAP)

geno <- data.frame(marker=SNPMAP$rs,chrom=SNPMAP$chrom,pos=SNPMAP$pos,t(geno))

geno[1:5,1:5]
dim(geno)


PhenoGWAS[1:5,]
colnames(PhenoGWAS)[1]<-"line"
PhenoGWAS$BReplicate<-as.factor(PhenoGWAS$BReplicate)
PhenoGWAS$TReplicate<-as.factor(PhenoGWAS$TReplicate)
PhenoGWAS$line<-as.factor(PhenoGWAS$line)
PhenoGWAS$Endophyte <-as.factor(PhenoGWAS$Endophyte)
PhenoGWAS[1:5,]




## -------------------------------------------------------------------------------------------
# format data for GAPIT
# phenotypes
# we will get the BLUPS for each of the phenotypes for the isolates (taxa) using lmer
table(PhenoGWAS$Endophyte)
PhenoGWASR3E3<-PhenoGWAS[PhenoGWAS$Endophyte=="R3E3",]
PhenoGWASR3E3<-na.omit(PhenoGWASR3E3)



GenoMat<-geno[,-c(1:3)]
GenoMat<-t(GenoMat)+1
GenoMat<-as.data.frame(GenoMat)

# rownames as the first column as TAxa
GenoMat$Taxa<-rownames(GenoMat)
rownames(GenoMat)<-NULL
# make taxa the first column
GenoMat<-GenoMat[,c(ncol(GenoMat), 1:(ncol(GenoMat)-1))]


GenoMat[1:5,1:5]
#str(GenoMat)
sum(is.na(GenoMat))
GenoMat$Taxa[1:5]
filtered_pcs$Isolate[1:5]
filtered_pcs$Taxa[1:5]
GenoMat$Taxa<-gsub("\\.","-", GenoMat$Taxa)
GenoMat<-GenoMat[GenoMat$Taxa%in%filtered_pcs$Taxa,]
GenoMat<-GenoMat[match(filtered_pcs$Taxa,GenoMat$Taxa),]

GenoMat$Taxa<-filtered_pcs$Isolate

GenoMat[1:5,1:5]
PhenoGWASR3E3$line<-factor(PhenoGWASR3E3$line, levels=GenoMat$Taxa)

################################################




# we need the BLUPs for the following traits "Radius_inhibition", "Area_Endophyte", "Area_inhibition", "Ratio_inhibition", "Ratio_Scales1", "Ratio_Scales2", "Percentage_inhibition"
ModelLMERadius_inhibition<-lmer(Radius_inhibition~BReplicate+TReplicate+(1|line), data=PhenoGWASR3E3)
BLUPRadius_inhibitionR3E3<-ranef(ModelLMERadius_inhibition)$line[,1]
ModelLMEArea_Endophyte<-lmer(Area_Endophyte~BReplicate+TReplicate+(1|line), data=PhenoGWASR3E3)
BLUPArea_EndophyteR3E3<-ranef(ModelLMEArea_Endophyte)$line[,1]
ModelLMEArea_inhibition<-lmer(Area_inhibition~BReplicate+TReplicate+(1|line), data=PhenoGWASR3E3)
BLUPArea_inhibitionR3E3<-ranef(ModelLMEArea_inhibition)$line[,1]
ModelLMERatio_inhibition<-lmer(Ratio_inhibition~BReplicate+TReplicate+(1|line), data=PhenoGWASR3E3)
BLUPRatio_inhibitionR3E3<-ranef(ModelLMERatio_inhibition)$line[,1]
ModelLMERatio_Scales1<-lmer(Ratio_Scales1~BReplicate+TReplicate+(1|line), data=PhenoGWASR3E3)
BLUPRatio_Scales1R3E3<-ranef(ModelLMERatio_Scales1)$line[,1]
ModelLMERatio_Scales2<-lmer(Ratio_Scales2~BReplicate+TReplicate+(1|line), data=PhenoGWASR3E3)
BLUPRatio_Scales2R3E3<-ranef(ModelLMERatio_Scales2)$line[,1]
ModelLMEPercentage_inhibition<-lmer(Percentage_inhibition~BReplicate+TReplicate+(1|line), data=PhenoGWASR3E3)
BLUPPercentage_inhibitionR3E3<-ranef(ModelLMEPercentage_inhibition)$line[,1]





# put these in a dataframe
BLUPS<-data.frame(Isolate=rownames(ranef(ModelLMERadius_inhibition)$line), BLUPRadius_inhibitionR3E3, BLUPArea_EndophyteR3E3, BLUPArea_inhibitionR3E3, BLUPRatio_inhibitionR3E3, BLUPRatio_Scales1R3E3, BLUPRatio_Scales2R3E3, BLUPPercentage_inhibitionR3E3)
BLUPS[1:5,]
# first column is the isolate name, change this to Taxa
colnames(BLUPS)[1]<-"Taxa"


## ----eval=TRUE, include=TRUE----------------------------------------------------------------

source("http://zzlab.net/GAPIT/gapit_functions.txt")
pacman::p_load(AssocTests)


mapping<-data.frame(Name=SNPMAP$rs, Chromosome=SNPMAP$chrom, Position=SNPMAP$pos)
models = c("FarmCPU","Blink")
# only the first 13 chromosomes
rs13<-mapping[mapping$Chromosome<=13,]$Name
GenoMat13<-GenoMat[,c("Taxa", rs13)]
mapping<-mapping[mapping$Chromosome<=13,]
# make a folder to save the output files "GWASforPilarTraitsMarch27th"
if (!file.exists("GWASforPilarTraitsMarch27th")){
dir.create("GWASforPilarTraitsMarch27th", showWarnings = FALSE)
}



#############################################################################
### thinning again but another method
setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts")


if (!file.exists("GWASforPilarTraitsMarch27thThinnedSnpsWithCor500")){
  dir.create("GWASforPilarTraitsMarch27thThinnedSnpsWithCor500", showWarnings = FALSE)
}

setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/GWASforPilarTraitsMarch27thThinnedSnpsWithCor500")


library(dplyr)
library(tidyr)
library(purrr)
library(pbapply)

find_highest_corr_snp <- function(df, GenoMat, distance) {
  # Preparing the data frame
  df <- df %>%
    arrange(Chromosome, Position) %>%
    mutate(region = floor(Position / distance))
  
  # Splitting the dataframe into a list of dataframes for each region
  df_list <- split(df, list(df$Chromosome, df$region))
  
  # Function to process each region
  process_region <- function(region_df) {
    region_snp_names <- region_df$Name
    region_data <- GenoMat[, region_snp_names, drop = FALSE]
    
    if(ncol(region_data) < 2) {
      return(region_df[1,])
    }
    
    corr_matrix <- cor(region_data, use = "pairwise.complete.obs")
    avg_corr <- colMeans(corr_matrix, na.rm = TRUE)
    highest_corr_snp <- names(which.max(avg_corr))
    
    region_df %>% filter(Name == highest_corr_snp)
  }
  
  # Applying the function to each region with a progress bar
  highest_corr_snps_list <- pbapply::pblapply(df_list, process_region)
  
  # Combining the results into a single data frame
  highest_corr_snps <- bind_rows(highest_corr_snps_list)
  
  highest_corr_snps
}

# Apply thinning function to your mapping dataframe for the first 13 chromosomes
mapping13 <- mapping %>% 
  filter(Chromosome <= 13)

GenoMat13<-GenoMat[,c("Taxa", mapping13$Name)]  

# Applying the function with a progress bar
mapping_highest_corr <- find_highest_corr_snp(mapping13, GenoMat13, 500) # Adjust the distance as needed

mapping_highest_corr<-na.omit(mapping_highest_corr)
dim(mapping_highest_corr)
mapping13_thinned<-mapping13[mapping13$Name%in%mapping_highest_corr$Name,]
# Apply the function
# Make sure `mapping` and `genoData` are properly defined before this step
dim(mapping13_thinned)

GenoMat13_thinned <- GenoMat13[, c("Taxa", mapping13_thinned$Name)]
dim(GenoMat13_thinned)
sum(is.na(GenoMat13_thinned))
dim(mapping13_thinned)
sum(is.na(BLUPS))
tmp<-capture.output({
  modK <- GAPIT(Y=BLUPS,
                GD=GenoMat13_thinned,
                GM=mapping13_thinned,
                KI=NULL,
                CV=NULL,
                PCA.total=2,
                model=models,
                file.output=TRUE,
                cutOff = 0.1)
})

rm(tmp)
gc()
modK


