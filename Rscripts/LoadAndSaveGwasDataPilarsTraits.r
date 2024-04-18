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

# Read numeric genotypes file
geno <- fread("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/6_tassel_analysis/numeric_allele_counts_filtered.txt", skip = 1, header = TRUE, sep = "\t", na.strings = c("", "NA",0.5,"0.5"))
geno <- as.data.frame(geno)

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
phenotypes<-read.csv("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/PhenoData/20240323_R3E3&R3E5_R1&R2_AreaRatioPercScalesRemoval_GWAS.csv", header=TRUE)

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


PhenoGWAS<-phenotypes[,c("Isolate", "BReplicate", "TReplicate", "Endophyte","Area_Endophyte", "Area_inhibition", "Radius_inhibition", "Ratio_inhibition", "Ratio_Scales1", "Ratio_Scales2", "Percentage_inhibition")]


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

colnames(PhenoGWAS)[1]<-"line"
PhenoGWAS$BReplicate<-as.factor(PhenoGWAS$BReplicate)
PhenoGWAS$TReplicate<-as.factor(PhenoGWAS$TReplicate)
PhenoGWAS$line<-as.factor(PhenoGWAS$line)
PhenoGWAS$Endophyte <-as.factor(PhenoGWAS$Endophyte)



PhenoGWASR3E3<-PhenoGWAS[PhenoGWAS$Endophyte=="R3E3",]
PhenoGWASR3E3<-na.omit(PhenoGWASR3E3)

PhenoGWASR3E5<-PhenoGWAS[PhenoGWAS$Endophyte=="R3E5",]
PhenoGWASR3E5<-na.omit(PhenoGWASR3E5)

GenoMat<-geno[,-c(1:3)]
GenoMat<-t(GenoMat)+1
GenoMat<-as.data.frame(GenoMat)

# rownames as the first column as TAxa
GenoMat$Taxa<-rownames(GenoMat)
rownames(GenoMat)<-NULL
# make taxa the first column
GenoMat<-GenoMat[,c(ncol(GenoMat), 1:(ncol(GenoMat)-1))]


GenoMat$Taxa<-gsub("\\.","-", GenoMat$Taxa)
GenoMat<-GenoMat[GenoMat$Taxa%in%filtered_pcs$Taxa,]
GenoMat<-GenoMat[match(filtered_pcs$Taxa,GenoMat$Taxa),]

GenoMat$Taxa<-filtered_pcs$Isolate

GenoMat[1:5,1:5]

PhenoGWASR3E3$line<-factor(PhenoGWASR3E3$line, levels=GenoMat$Taxa)
PhenoGWASR3E5$line<-factor(PhenoGWASR3E5$line, levels=GenoMat$Taxa)
PhenoGWAS$line<-factor(PhenoGWAS$line, levels=GenoMat$Taxa)

# Save all data needed for GWAS
new_dir <- "DataforPillar"
if (!dir.exists(new_dir)) dir.create(new_dir, recursive = TRUE)
save(PhenoGWASR3E3, PhenoGWASR3E5, PhenoGWAS, GenoMat, file = paste0(new_dir, "/DataforPillar.RData"))

