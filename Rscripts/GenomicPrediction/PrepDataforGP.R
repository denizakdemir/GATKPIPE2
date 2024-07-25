setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/GenomicPrediction")
# Load necessary packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(
  tidyverse, data.table, ggplot2, dplyr, stringr,scatterplot3d
)


load("Pheno_geno_Map_PillarTraits.RData")


genoImp[1:5,1:5]
pcs2022[1:5,]
pcs2023[1:5,]

Pcs<-rbind(pcs2022,pcs2023)
intersect(Pcs$Taxa, rownames(genoImp))

pcs2023$Taxa[1:5]
rownames(genoImp)[1:5]
intersect(pcs2023$Taxa,rownames(genoImp))
intersect(pcs2023$Taxa, pcs2022$Taxa)
intersect(rownames(genoImp),Pcs$Taxa)
Pcs<-Pcs[match(rownames(genoImp),Pcs$Taxa),]
Pcs[1:10,]
rownames(genoImp)<-Pcs$Isolate


mix1<-c("22_EcijaSec83Ica_L2", "22_EcijaSecCris_L1", "22_EcijaRegTej_L1")
intersect(mix1, rownames(genoImp))
mix2<-c("22_CorKiko_L1", "22_CorCale_L1", "22_Cor3927_L1")
intersect(mix2, rownames(genoImp))
mix3<-c("22_ConilAmi_L1", "22_Conil3806_L1", "22_Jerez3927_L1")
intersect(mix3, rownames(genoImp))
mix4<-c("22_CorVal_L1")
intersect(mix4, rownames(genoImp))

Mixes<-list(mix1, mix2, mix3, mix4)

gentmeangeno<-function(mix){
  mix<-intersect(mix, rownames(genoImp))
  print(mix)
  matMix<-genoImp[rownames(genoImp)%in%mix,]
  if (length(mix)>1){
   meanGeno<-colMeans(matMix)
  } else {
    meanGeno<-matMix
  }
  return(meanGeno)
}

MixgenoImp<-t(sapply(Mixes, gentmeangeno))

dim(MixgenoImp)
MixgenoImp[1:4,1:5]
rownames(MixgenoImp)<-c("mix1", "mix2", "mix3", "mix4")


genoImpwithMixes<-rbind(genoImp,MixgenoImp)
dim(genoImpwithMixes)
genoImpwithMixes<-genoImpwithMixes[!duplicated(rownames(genoImpwithMixes)),]

KmatSepMix<-rrBLUP::A.mat(genoImpwithMixes-1)
heatmap(KmatSepMix, symm = TRUE, scale = "none", col = colorRampPalette(c("white", "blue"))(100), margins = c(10, 10))
round(KmatSepMix, 2)


# load wheat Data
load("WheatData/SNPmatrixFiltSubset_inbreed.RData")

ls()

GenoWheat<-snp_matrix_sub
MapWheat<-MAP
rm(snp_matrix_sub, MAP)
GenoWheat[1:5,1:5]
dim(GenoWheat)
KmatWheat<-rrBLUP::A.mat(GenoWheat-1)

heatmap(KmatWheat, symm = TRUE, scale = "none", col = colorRampPalette(c("white", "blue"))(100), margins = c(10, 10))

rownames(KmatWheat)
## phenotypic data 12cm

# data is at WheatData/Output_reps_unidas_12cm_mock_checks.xlsx
# load data
library(readxl)
Pheno12cm<-read_excel("WheatData/Output_reps_unidas_12cm_mock_checks.xlsx")

Pheno12cm[1:5,1:7]
# make sure these are turned into 3 number character strings
Pheno12cm$Plant<-as.character(Pheno12cm$Plant)
Pheno12cm$Plant<-str_pad(Pheno12cm$Plant, 3, pad = "0")
Pheno12cm$Plant<-paste("PyrSep.", Pheno12cm$Plant, sep = "")
intersect(Pheno12cm$Plant, rownames(KmatWheat))
setdiff(Pheno12cm$Plant, rownames(KmatWheat))
setdiff(rownames(KmatWheat), Pheno12cm$Plant)
length(unique(Pheno12cm$Plant))

Pheno17cm<-read_excel("WheatData/Output_reps_unidas_17cm_mock_checks.xlsx")
Pheno17cm$Plant<-as.character(Pheno17cm$Plant)
Pheno17cm$Plant<-str_pad(Pheno17cm$Plant, 3, pad = "0")
Pheno17cm$Plant<-paste("PyrSep.", Pheno17cm$Plant, sep = "")
intersect(Pheno17cm$Plant, rownames(KmatWheat))
setdiff(Pheno17cm$Plant, rownames(KmatWheat))
setdiff(rownames(KmatWheat), Pheno17cm$Plant)

# save wheat and septoria data
save(GenoWheat, MapWheat, KmatWheat, Pheno12cm, Pheno17cm,genoImp,SNPMAP,KmatSepMix, file = "DataforGenomicPrediction.RData")

rm(list = ls())
gc()
load("DataforGenomicPrediction.RData")

ls()


##########