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



## ----load_genomic_data, eval=FALSE----------------------------------------------------------
## # read numeric genotypes file:
geno <- fread("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/6_tassel_analysis/numeric_allele_counts_filtered.txt", skip = 1, header = TRUE, sep = "\t", na.strings = c("", "NA",0.5,"0.5"))
geno<-as.data.frame(geno)
geno[,-1][geno[,-1] >1] <- NA

geno[1:5,1:5]
table(geno[,2])
table(c(unlist(geno[,2:1000])), useNA = "always")

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


p





## -------------------------------------------------------------------------------------------
# use  ('Taxa' contains 'EKDN2300427') as color for the points (grepl("EKDN2300427", Taxa)).
# Create a ggplot object with smaller font size for labels

pcs$Year<-ifelse(grepl("EKDN2300427", pcs$Taxa), "2022", "2023")
p <- ggplot(pcs, aes(x=PC1, y=PC2, color=Year)) +
  geom_point()  +
  geom_text_repel(aes(label=Taxa), size=.5, 
                  box.padding = .1, 
                  point.padding = 0.1,
                  segment.color = 'grey50', max.overlaps = 120) +
  xlim(min(pcs$PC1) - 5, max(pcs$PC1) + 5) + # Adjust xlim
  ylim(min(pcs$PC2) - 5, max(pcs$PC2) + 5)   # Adjust ylim
p

## -------------------------------------------------------------------------------------------
# Assuming genotype_df is your dataframe
filtered_pcs <- pcs %>%
  filter(!grepl("EKDN2300427", Taxa)) # Remove rows where 'Taxa' contains 'EKDN2300427'

p<-ggplot(filtered_pcs, aes(x=PC1, y=PC2, label=Taxa)) + geom_point()  +
  geom_text_repel(aes(label=Taxa), size=.5, 
                  box.padding = .1, 
                  point.padding = 0.1,
                  segment.color = 'grey50', max.overlaps = 120) +
  xlim(min(pcs$PC1) - 5, max(pcs$PC1) + 5) + # Adjust xlim
  ylim(min(pcs$PC2) - 5, max(pcs$PC2) + 5)   # Adjust ylim
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
pcaplot

ggsave("pcaplot100withnamesinPheno.png", pcaplot, device="png")



## -------------------------------------------------------------------------------------------
 
# read the file


# load the phenotypic data from "/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/PhenoData/raw_phenotypes.csv"

# read the file
phenotypes<-read.csv("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/PhenoData/raw_phenotypes.csv")
head(phenotypes)
phenotypes$Picture[1:10]

library(dplyr)
library(stringr)

# Assuming phenotypes is your dataframe and Picture is the column of interest
phenotypes <- phenotypes %>%
  mutate(
    WheatGeno = str_extract(Picture, "^[^_]+"), # Extracts everything before the first underscore
    Year = str_extract(Picture, "(?<=_)\\d+") # Extracts digits following the first underscore
  )



# convert all commas to dots in filtered_pcs$Isolate
filtered_pcs$Isolate<-gsub(",", ".", filtered_pcs$Isolate)
head(phenotypes)



phenotypes$Picture[1:10]
filtered_pcs$Isolate[1:10]
###

# go over each row of the phenotypes dataframe and find the Isolate name that is in the picture name
# put that name in the phenotype file in a new column
for (i in 1:nrow(phenotypes)){
  for (j in 1:nrow(filtered_pcs)){
    if (grepl(filtered_pcs$Isolate[j], phenotypes$Picture[i])){
      phenotypes$Isolate[i]<-filtered_pcs$Isolate[j]
    }
  }
}
table(phenotypes$Isolate)
sum(is.na(phenotypes$Isolate))

length(unique(phenotypes$Isolate))
dim(geno)
geno[1:5,1:5]

filtered_pcs$Isolate<-gsub(" ","", filtered_pcs$Isolate)
intersect(phenotypes$Isolate, filtered_pcs$Isolate)
setdiff(filtered_pcs$Isolate, phenotypes$Isolate)
setdiff(phenotypes$Isolate,filtered_pcs$Isolate)





intersect(colnames(phenotypes), c(c("Isolate", "REP", "leave_id","WheatGeno", "Year", "PLACL", "pycnidiaPerCm2Leaf",
"pycnidiaPerCm2Lesion", "necrosisArea", "pycnidiaCount", "leafArea", "meanPycnidiaArea","pycnidiaGreyValue")))
PhenoGWAS<-phenotypes[,c(c("Isolate", "REP", "leave_id","WheatGeno", "Year", "PLACL", "pycnidiaPerCm2Leaf",
"pycnidiaPerCm2Lesion", "necrosisArea", "pycnidiaCount", "leafArea", "meanPycnidiaArea","pycnidiaGreyValue"))]



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



## ----------------------


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
PhenoGWAS$line<-factor(PhenoGWAS$line, levels=GenoMat$Taxa)
# we need the BLUPs for the following traits PLACL"               
# [7] "pycnidiaPerCm2Leaf"   "pycnidiaPerCm2Lesion" "necrosisArea"        
#[10] "pycnidiaCount"        "leafArea"             "meanPycnidiaArea"    
#[13] "pycnidiaGreyValue"   

# These are the fixed effects:  "REP"                  "leave_id"            
# [4] "WheatGeno"            "Year" 

PhenoGWAS$REP<-as.factor(PhenoGWAS$REP)
PhenoGWAS$leave_id<-as.factor(PhenoGWAS$leave_id)
PhenoGWAS$WheatGeno<-as.factor(PhenoGWAS$WheatGeno)
PhenoGWAS$Year<-as.factor(PhenoGWAS$Year)

library(lme4)
ModelLMEPLACL<-lmer(PLACL~REP+leave_id+WheatGeno+Year+(1|line), data=PhenoGWAS)
BLUPPLACL<-ranef(ModelLMEPLACL)$line[,1]
ModelLMEpycnidiaPerCm2Leaf<-lmer(pycnidiaPerCm2Leaf~REP+leave_id+WheatGeno+Year+(1|line), data=PhenoGWAS)
BLUPpycnidiaPerCm2Leaf<-ranef(ModelLMEpycnidiaPerCm2Leaf)$line[,1]
ModelLMEpycnidiaPerCm2Lesion<-lmer(pycnidiaPerCm2Lesion~REP+leave_id+WheatGeno+Year+(1|line), data=PhenoGWAS)
BLUPpycnidiaPerCm2Lesion<-ranef(ModelLMEpycnidiaPerCm2Lesion)$line[,1]
ModelLMEnecrosisArea<-lmer(necrosisArea~REP+leave_id+WheatGeno+Year+(1|line), data=PhenoGWAS)
BLUPnecrosisArea<-ranef(ModelLMEnecrosisArea)$line[,1]
ModelLMEpycnidiaCount<-lmer(pycnidiaCount~REP+leave_id+WheatGeno+Year+(1|line), data=PhenoGWAS)
BLUPpycnidiaCount<-ranef(ModelLMEpycnidiaCount)$line[,1]
ModelLMEleafArea<-lmer(leafArea~REP+leave_id+WheatGeno+Year+(1|line), data=PhenoGWAS)
BLUPleafArea<-ranef(ModelLMEleafArea)$line[,1]
ModelLMEmeanPycnidiaArea<-lmer(meanPycnidiaArea~REP+leave_id+WheatGeno+Year+(1|line), data=PhenoGWAS)
BLUPmeanPycnidiaArea<-ranef(ModelLMEmeanPycnidiaArea)$line[,1]
ModelLMEpycnidiaGreyValue<-lmer(pycnidiaGreyValue~REP+leave_id+WheatGeno+Year+(1|line), data=PhenoGWAS)
BLUPpycnidiaGreyValue<-ranef(ModelLMEpycnidiaGreyValue)$line[,1]





# put these in a dataframe
BLUPS<-data.frame(Isolate=rownames(ranef(ModelLMEPLACL)$line),BLUPPLACL, BLUPpycnidiaPerCm2Leaf, BLUPpycnidiaPerCm2Lesion, BLUPnecrosisArea, BLUPpycnidiaCount, BLUPleafArea, BLUPmeanPycnidiaArea, BLUPpycnidiaGreyValue)
BLUPS[1:5,]
# first column is the isolate name, change this to Taxa
colnames(BLUPS)[1]<-"Taxa"


## ----eval=TRUE, include=TRUE----------------------------------------------------------------
#############################################################################
### thinning again but another method
setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts")


source("http://zzlab.net/GAPIT/gapit_functions.txt")
pacman::p_load(AssocTests)


mapping<-data.frame(Name=SNPMAP$rs, Chromosome=SNPMAP$chrom, Position=SNPMAP$pos)
models = c("FarmCPU","Blink")
# only the first 13 chromosomes
rs13<-mapping[mapping$Chromosome<=13,]$Name
GenoMat13<-GenoMat[,c("Taxa", rs13)]
mapping<-mapping[mapping$Chromosome<=13,]
# make a folder to save the output files "GWASforPilarTraitsMarch27th"
if (!file.exists("GWASforGreenHouseTraitsMarch27th")){
dir.create("GWASforGreenHouseTraitsMarch27th", showWarnings = FALSE)
}


setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/GWASforGreenHouseTraitsMarch27th")
tmp<-capture.output({
modK <- GAPIT(Y=BLUPS,
                GD=GenoMat13,
                GM=mapping,
                KI=NULL,
                CV=NULL,
                PCA.total=3,
                model=models,
                file.output=TRUE,
              cutOff = 0.1)
})

rm(tmp)
gc()


######### Another time but with thinned snps
setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts")


if (!file.exists("GWASforGreenHouseTraitsMarch27thThinnedSnps")){
  dir.create("GWASforGreenHouseTraitsMarch27thThinnedSnps", showWarnings = FALSE)
}

setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/GWASforGreenHouseTraitsMarch27thThinnedSnps")


library(dplyr)

# Function to thin SNPs within each chromosome
thin_snps <- function(df, distance) {
  df %>% 
    arrange(Chromosome, Position) %>%
    group_by(Chromosome) %>%
    mutate(kb_block = floor(Position / distance)) %>%
    group_by(Chromosome, kb_block) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::select(-kb_block)
}

# Apply thinning function to your mapping dataframe for the first 13 chromosomes
mapping_thinned <- mapping %>% 
  filter(Chromosome <= 13) %>% 
  thin_snps(distance = 100) %>%as.data.frame() # Thinning within each 100 kb
sum(is.na(mapping_thinned))
dim(mapping_thinned)
# Update SNP list and GenoMat accordingly
rs13_thinned <- mapping_thinned$Name
GenoMat13_thinned <- GenoMat[, c("Taxa", rs13_thinned)]
GenoMat13_thinned[1:5,1:5]
#str(GenoMat13_thinned)
dim(GenoMat13_thinned)
sum(is.na(GenoMat13_thinned))
dim(mapping_thinned)
sum(is.na(BLUPS))
tmp<-capture.output({
  modK <- GAPIT(Y=BLUPS,
                GD=GenoMat13_thinned,
                GM=mapping_thinned,
                KI=NULL,
                CV=NULL,
                PCA.total=3,
                model=models,
                file.output=TRUE)
})

rm(tmp)
gc()
modK

#############################################################################
### thinning again but another method
setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts")


if (!file.exists("GWASforGreenHouseTraitsMarch27thThinnedSnpsWithCor")){
  dir.create("GWASforGreenHouseTraitsMarch27thThinnedSnpsWithCor", showWarnings = FALSE)
}

setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/GWASforGreenHouseTraitsMarch27thThinnedSnpsWithCor")


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

mapping_thinned<-mapping13[mapping13$Name%in%mapping_highest_corr$Name,]
# Apply the function
# Make sure `mapping` and `genoData` are properly defined before this step
dim(mapping_thinned)
# Process to update `GenoMat` and `rs13_thinned` based on `mapping_thinned`
GenoMat13_thinned <- GenoMat13[, c("Taxa", mapping_thinned$Name)]
dim(GenoMat13_thinned)
sum(is.na(GenoMat13_thinned))
dim(mapping_thinned)
sum(is.na(BLUPS))
tmp<-capture.output({
  modK <- GAPIT(Y=BLUPS,
                GD=GenoMat13_thinned,
                GM=mapping_thinned,
                KI=NULL,
                CV=NULL,
                PCA.total=4,
                model=models,
                file.output=TRUE)
})

rm(tmp)
gc()
modK

#######
