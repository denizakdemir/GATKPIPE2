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

load("ImputedGeno119Genotypes.RData")
load("PCsforGeno119Genotypes.RData")

rownames(pcs)
dim(genoImp)


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


# Load phenotypes and ensure matching with genotypes
phenotypes<-read.csv("PhenoData/raw_phenotypes.csv" , header=TRUE)
phenotypes[1:5,]

library(magrittr)
library(tidyverse)
# Extract information from the Picture column
phenotypes <- phenotypes %>%
  mutate(
    WheatGeno = str_extract(Picture, "^[^_]+"),
    Year = str_extract(Picture, "(?<=_)\\d+")
  )

# Convert commas to dots in filtered_pcs$Isolate
pcs2022$Isolate <- gsub(",", ".", pcs2022$Isolate)
pcs2023$Isolate <- gsub(",", ".", pcs2023$Isolate)

# Match phenotypes with genotypes
for (i in 1:nrow(phenotypes)) {
  for (j in 1:nrow(pcs2022)) {
    if (grepl(pcs2022$Isolate[j], phenotypes$Picture[i])) {
      phenotypes$Isolate[i] <- pcs2022$Isolate[j]
    }
  }
}

# Ensure all required columns exist in the phenotypes dataframe
required_columns <- c("Isolate", "leave_id", "REP", "WheatGeno", "Year", "PLACL", "pycnidiaPerCm2Leaf",
                      "pycnidiaPerCm2Lesion")


missing_columns <- setdiff(required_columns, colnames(phenotypes))
if (length(missing_columns) > 0) {
  stop(paste("The following required columns are missing from the phenotypes dataframe:", paste(missing_columns, collapse = ", ")))
}

# Select relevant columns for GWAS analysis
phenotypes <- phenotypes[,required_columns]


pcs2023$Isolate[1:10]
pcs2022$Isolate[1:10]


# Correct isolate names to match
pcs2022$Isolate[pcs2022$Isolate=="22_Conil_Fer_L1"]<-"22_ConilFer_L1"
pcs2022$Isolate[pcs2022$Isolate=="22_Jerez Val_L1"]<-"22_JerezVal_L1"
pcs2022$Isolate[pcs2022$Isolate=="22_EcijaSecOrt _L1"]<-"22_EcijaSecOrt_L1"
pcs2022$Isolate[pcs2022$Isolate=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2" ## Check this one
pcs2022$Isolate[pcs2022$Isolate=="22_EscIca_L2"]<-"22_EscIca_L2"
pcs2022$Isolate[pcs2022$Isolate=="22_EcijaSecSha_L1"]<-"22_EcijaSecSah_L1"





SNPSummary<-read.table("6_tassel_analysis/snp_genotypeSummary3.txt", header = TRUE, sep = "\t", as.is = TRUE, skip=0)
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



save(phenotypes, genoImp, SNPMAP,pcs2022,pcs2023, file="Pheno_geno_Map_GreenHouseTraits.RData")


