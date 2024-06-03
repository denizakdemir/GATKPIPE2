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
geno <- fread("6_tassel_analysis/numeric_allele_counts_filtered.txt", skip = 1, header = TRUE, sep = "\t", na.strings = c("", "NA",0.5,"0.5"))
geno <- as.data.frame(geno)
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
dublicatedGenotypeNames<-geno$Taxa[duplicatednames]


geno<-geno[!duplicatednames,]
rownames(geno)<-geno$Taxa
geno<-geno[,-1]
geno<-as.matrix(geno)*2-1
rownames(geno)
# remove any rows that have "sorted" in it
removeNames<-grepl("sorted", rownames(geno))
geno<-geno[!removeNames, ]

dim(geno)

genoImp<-rrBLUP::A.mat(geno, impute.method = "mean", return.imputed = TRUE, min.MAF = .1, max.missing = .5)$imputed
dim(genoImp)

save(genoImp, file="ImputedGeno119Genotypes.RData")
rm(geno)
gc()

svdG<-svd(scale(genoImp, scale=FALSE, center=TRUE), nu=5,nv=5)
pcs<-scale(genoImp, scale=FALSE, center=TRUE)%*%svdG$v
pcs<-cbind(data.frame(Taxa=rownames(pcs)), pcs)
colnames(pcs)[2:6]<-paste("PC",1:5, sep="")

pcs$Year<-ifelse(grepl("EKDN2300427", pcs$Taxa), "2023", "2022")

save(pcs, file="PCsforGeno119Genotypes.RData")

rownames(pcs)





pcs2022 <- pcs[pcs$Year=="2022",]
pcs2023 <- pcs[pcs$Year=="2023",]
## ----match_genotype_phenotype---------------------------------------------------------------
# one of the names have an extra "combined_" in it, remove it
pcs2022$TaxaShort<-gsub("combined_", "", pcs2022$Taxa)
pcs2022$TaxaShort

pcs2023$TaxaShort<-gsub("combined_", "", pcs2023$Taxa)
pcs2023$TaxaShort

pcs2022$TaxaShort<-substr(pcs2022$Taxa, 1, 4)
pcs2023$TaxaShort<-substr(pcs2023$Taxa, 1, 4)

# remove _ and anything that comes after that from the TaxaShort column
pcs2022$TaxaShort<-gsub("_.*", "", pcs2022$TaxaShort)
pcs2023$TaxaShort<-gsub("_.*", "", pcs2023$TaxaShort)

# sort the data as S1, S2,...
pcs2022<-pcs2022[order(as.numeric(gsub("S", "", pcs2022$TaxaShort))),]
pcs2022$Taxa

pcs2023<-pcs2022[order(as.numeric(gsub("S", "", pcs2023$TaxaShort))),]
pcs2023$Taxa


# read the file
isolateslist1<-readxl::read_xlsx("GenotypeTables/isolated\ list\ 2021-2022\ for\ sequencing.xlsx")
dim(isolateslist1)
isolateslist2<-readxl::read_xlsx("GenotypeTables/isolated list 2023 for sequencing.xlsx")
dim(isolateslist2)

# have isolate and Name columns in the isolateslist dataframe
# the name column has the S1, S2, ...S100 names
# now by matching the Name column with the TaxaShort column of the pcs2022 dataframe, get the isolate column of the isolateslist dataframe

pcs2022$Isolate<-isolateslist1$Isolate[match(pcs2022$TaxaShort, isolateslist1$Name)]
pcs2022$Isolate<-gsub(",", ".", pcs2022$Isolate)


dim(pcs2022)

pcs2023$Isolate<-isolateslist1$Isolate[match(pcs2023$TaxaShort, isolateslist1$Name)]
pcs2023$Isolate<-gsub(",", ".", pcs2023$Isolate)


dim(pcs2023)

pcs2023

# Load phenotypes and ensure matching with genotypes
phenotypes<-read.csv("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/PhenoData/20240423_R3E3_R1&R2_GWAS_astringent&loose.csv", header=TRUE)


phenotypes$Zt_Strain[1:10]
pcs2023$Isolate[1:10]
pcs2022$Isolate[1:10]

###
phenotypes$Zt_Strain <- gsub("-", ".", phenotypes$Zt_Strain)



# Correct isolate names to match
pcs2022$Isolate[pcs2022$Isolate=="22_Conil_Fer_L1"]<-"22_ConilFer_L1"
pcs2022$Isolate[pcs2022$Isolate=="22_Jerez Val_L1"]<-"22_JerezVal_L1"
pcs2022$Isolate[pcs2022$Isolate=="22_EcijaSecOrt _L1"]<-"22_EcijaSecOrt_L1"
pcs2022$Isolate[pcs2022$Isolate=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2" ## Check this one
pcs2022$Isolate[pcs2022$Isolate=="22_EscIca_L2"]<-"22_EscIca_L2"
pcs2022$Isolate[pcs2022$Isolate=="22_EcijaSecSha_L1"]<-"22_EcijaSecSah_L1"



## -------------------------------------------------------------------------------------------
# put that name in the phenotype file in a new column
for (i in 1:nrow(phenotypes)){
  for (j in 1:nrow(pcs2022)){
    if (grepl(pcs2022$Isolate[j], phenotypes$Zt_Strain[i])){
      phenotypes$Isolate[i]<-pcs2022$Isolate[j]
    }
  }
}
table(phenotypes$Isolate)

head(phenotypes)


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
SNPMAP<-SNPMAP[SNPMAP$chrom<=13,]
dim(genoImp)
genoImp[1:5,1:5]
genoImp<-genoImp[, SNPMAP$rs]
dim(genoImp)
save(phenotypes, genoImp, SNPMAP,pcs2022,pcs2023, file="Pheno_geno_Map_PillarTraits.RData")


#############################################
