setwd("~/Google Drive/My Drive/Akdemir_Deniz/github/GATKPIPE2/Rscripts")


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

# create a directory to save the report files if the directory does not exist
if (!file.exists("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/FiguresForReport")){
dir.create("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/FiguresForReport", showWarnings = FALSE)
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
gc()
KMat<-rrBLUP::A.mat(geno, impute.method = "mean", return.imputed = FALSE, min.MAF = .1, max.missing = .3)
gc()


## ----pca------------------------------------------------------------------------------------
svdG<-svd(scale(geno, scale=FALSE, center=TRUE), nu=5,nv=5)
pcs<-scale(geno, scale=FALSE, center=TRUE)%*%svdG$v
pcs<-cbind(data.frame(Taxa=rownames(pcs)), pcs)
colnames(pcs)[2:6]<-paste("PC",1:5, sep="")
# percentage of variance explained by each PC
percvar<-svdG$d^2/sum(svdG$d^2)*100
percvar<-round(percvar, 2)[1:5]
# I need a nice looking pcA plot with the percentage of variance explained by each PC on the axes.
# plus plot Kmat heatmap
# add the screeplot as well
# the plot should be saved in the FiguresForReport directory. One plot for all of them. With a, b, c as figure titles

# plot the PCA
p1<-ggplot(pcs, aes(x=PC1, y=PC2))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+
  labs(title="PCA", x=paste("PC1 (", percvar[1], "%)", sep=""), y=paste("PC2 (", percvar[2], "%)", sep=""))

ggsave("FiguresForReport/PCA.png",p1, width=5, height=5, dpi=300)

# plot screeeplot

p2<-ggplot(data.frame(PC=1:5, percvar=percvar), aes(x=PC, y=percvar))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+
  labs(title="Screeplot", x="PC", y="Percentage of variance explained")
ggsave("FiguresForReport/Screeplot.png",p2, width=5, height=5, dpi=300)

# Kmat as image plot
# use corrplot
library(corrplot)
diag(KMat)<-0
corrplot::corrplot(KMat^.2, method="color", type="upper", tl.col="black", tl.srt=45, is.corr=FALSE, tl.cex=0.5)  




# combine the three plots
gridExtra::grid.arrange(p1, p2, ncol=3)




## -------------------------------------------------------------------------------------------
# use  ('Taxa' contains 'EKDN2300427') as color for the points (grepl("EKDN2300427", Taxa)).
# Create a ggplot object with smaller font size for labels

pcs$Year<-ifelse(grepl("EKDN2300427", pcs$Taxa), "2023", "2022")




## -------------------------------------------------------------------------------------------
# Assuming genotype_df is your dataframe
filtered_pcs <- pcs %>%
  filter(!grepl("EKDN2300427", Taxa)) # Remove rows where 'Taxa' contains 'EKDN2300427'



## ----match_genotype_phenotype---------------------------------------------------------------
# one of the names have an extra "combined_" in it, remove it
filtered_pcs$TaxaShort<-gsub("combined_", "", filtered_pcs$Taxa)

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





## -------------------------------------------------------------------------------------------
 
# read the file
phenotypes<-read.csv("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/PhenoData/20240323_R3E3&R3E5_R1&R2_AreaRatioPercScalesRemoval_GWAS.csv", header=TRUE)

head(phenotypes)

heatmap(cor(phenotypes[,5:12], use="pairwise.complete.obs"))
filtered_pcs$Isolate[1:10]

# convert all commas to dots in filtered_pcs$Isolate
filtered_pcs$Isolate<-gsub(",", ".", filtered_pcs$Isolate)

phenotypes$Zt_Strain <- gsub("-", ".", phenotypes$Zt_Strain)


filtered_pcs$Isolate[filtered_pcs$Isolate=="22_Conil_Fer_L1"]<-"22_ConilFer_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_Jerez Val_L1"]<-"22_JerezVal_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EcijaSecOrt _L1"]<-"22_EcijaSecOrt_L1"
filtered_pcs$Isolate[filtered_pcs$Isolate=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2"
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

# now I need PCs for GenoMat
svdG<-svd(scale(GenoMat[,-1], scale=FALSE, center=TRUE), nu=5,nv=5)
pcs<-scale(GenoMat[,-1], scale=FALSE, center=TRUE)%*%svdG$v
pcs<-cbind(data.frame(Taxa=GenoMat$Taxa), pcs)
colnames(pcs)[2:6]<-paste("PC",1:5, sep="")
# percentage of variance explained by each PC
percvar<-svdG$d^2/sum(svdG$d^2)*100
plot(percvar)
# calculate the Kmat
KMat<-rrBLUP::A.mat(GenoMat[,-1]-1, impute.method = "mean", return.imputed = FALSE, min.MAF = .1, max.missing = .3)

rownames(KMat)<-colnames(KMat)<-GenoMat$Taxa
GenoMat[1:5,1:5]
# plot the PCA
p1<-ggplot(pcs, aes(x=PC1, y=PC2))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+
  labs(title="PCA", x=paste("PC1 (", percvar[1], "%)", sep=""), y=paste("PC2 (", percvar[2], "%)", sep=""))

ggsave("FiguresForReport/PCA_GenoMat.png",p1, width=5, height=5, dpi=300)

# plot screeeplot
p2<-ggplot(data.frame(PC=1:5, percvar=percvar[1:5]), aes(x=PC, y=percvar))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+
  labs(title="Screeplot", x="PC", y="Percentage of variance explained")

# save the screeplot
ggsave("FiguresForReport/Screeplot_GenoMat.png",p2, width=5, height=5, dpi=300)

# Kmat as image plot
# use corrplot
library(corrplot)
diag(KMat)<-0
max(KMat)

corrplot::corrplot(KMat, method="circle", type="upper", tl.col="black", tl.srt=45, is.corr=FALSE, tl.cex=0.5,col=colorRampPalette(c("lightblue","blue","darkblue"))(100))

# save the above plot as png
png("FiguresForReport/Kmat_GenoMat.png", width=5, height=5, units="in", res=600)
corrplot::corrplot(KMat, method="square", type="upper", tl.col="black", tl.srt=45, is.corr=FALSE, tl.cex=0.5,col=colorRampPalette(c("lightblue","blue","darkblue"))(100), cl.cex = 0.5)
dev.off()

# use cowplot to combine the three plots,use labels a, b, c
# we need to read each as png and then combine them
library(png)

p1<-readPNG("FiguresForReport/PCA_GenoMat.png")
p2<-readPNG("FiguresForReport/Screeplot_GenoMat.png")
p3<-readPNG("FiguresForReport/Kmat_GenoMat.png")

# combine the three plots, Kmat on top, spanning the whole width, PCA and screeplot below in one row 
# we also need subfigure labels
library(cowplot)
p1<-rasterGrob(p1, interpolate=TRUE)
p2<-rasterGrob(p2, interpolate=TRUE)
p3<-rasterGrob(p3, interpolate=TRUE)

# combine the three plots, Kmat on top, spanning the whole width, PCA and screeplot below in one row
# we also need subfigure labels
p<-plot_grid(p3, plot_grid(p1, p2, ncol=2, labels=c("b", "c"), rel_widths=c(1,1)), labels="a", ncol=1, rel_heights=c(1, 0.5))
p
# save the combined plot
ggsave("FiguresForReport/Combined_GenoMat.png", p, width=10, height=10, dpi=600)


#####################################

library(lme4)

# get the BLUPs for the phenotypes

PhenoGWASR3E3$line<-factor(PhenoGWASR3E3$line, levels=GenoMat$Taxa)
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



# for each signoficant SNP, I would like to plot BLUPS against the SNP
# 	SNP	Chr	Pos	P.value	MAF	traits
# 44916	S2_3164606	2	3164606	2.077823e-07	0.1236207	FarmCPU.BLUPArea_inhibitionR3E3
# 449161	S2_3164606	2	3164606	1.219242e-07	0.1236207	FarmCPU.BLUPPercentage_inhibitionR3E3
# 449162	S2_3164606	2	3164606	2.077823e-07	0.1236207	BLINK.BLUPArea_inhibitionR3E3
# 449163	S2_3164606	 2	3164606	6.131392e-11	0.1236207	BLINK.BLUPPercentage_inhibitionR3E3
# 159597	S12_1035503	12	1035503	8.168641e-09	0.4676106	BLINK.BLUPPercentage_inhibitionR3E3


SNPS2_3164606<-GenoMat[,"S2_3164606"]
SNPS12_1035503<-GenoMat[,"S12_1035503"]

# plot the BLUPS against the SNPs
p1<-ggplot(BLUPS, aes(x=SNPS2_3164606, y=BLUPArea_inhibitionR3E3))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+
  labs(title="BLUPArea_inhibitionR3E3 vs SNP S2_3164606", x="SNP S2_3164606", y="BLUPArea_inhibitionR3E3")

p2<-ggplot(BLUPS, aes(x=SNPS2_3164606, y=BLUPPercentage_inhibitionR3E3))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+
  labs(title="BLUPPercentage_inhibitionR3E3 vs SNP S2_3164606", x="SNP S2_3164606", y="BLUPPercentage_inhibitionR3E3")

p3<-ggplot(BLUPS, aes(x=SNPS12_1035503, y=BLUPPercentage_inhibitionR3E3))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+
  labs(title="BLUPPercentage_inhibitionR3E3 vs SNP S12_1035503", x="SNP S12_1035503", y="BLUPPercentage_inhibitionR3E3")

# save the plots
ggsave("FiguresForReport/BLUPArea_inhibitionR3E3_vs_SNP_S2_3164606.png",p1, width=5, height=5, dpi=300)
ggsave("FiguresForReport/BLUPPercentage_inhibitionR3E3_vs_SNP_S2_3164606.png",p2, width=5, height=5, dpi=300)
ggsave("FiguresForReport/BLUPPercentage_inhibitionR3E3_vs_SNP_S12_1035503.png",p3, width=5, height=5, dpi=300)

# now group Taxa into Groups based on the SNPs (number of alleles are diffeerent for each SNP)

# first SNP
# S2_3164606
# group Taxa based on this SNP


# put this information on a table 

SNPS2_3164606Grouping<-data.frame(Taxa=GenoMat$Taxa, SNPS2_3164606)

SNPS12_1035503Grouping<-data.frame(Taxa=GenoMat$Taxa, SNPS12_1035503)

# one table 
SNPGrouping<-merge(SNPS2_3164606Grouping, SNPS12_1035503Grouping, by="Taxa")

# save the table
write.csv(SNPGrouping, "FiguresForReport/SNPGrouping.csv", row.names=FALSE)



