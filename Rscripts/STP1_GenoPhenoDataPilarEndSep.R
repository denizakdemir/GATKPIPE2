# STP1 starts with data preporocessed by TASSEL


# Set the working directory to the location of your file, if necessary
setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/ProcessedVCFs/SNPGENOTASSELOUTPUT/")


# read numeric genotypes file:
library(data.table)
GenoMat <- fread("6_tassel_analysis/numeric_allele_counts_filtered.txt", skip = 1, header = TRUE, sep = "\t")
#GenoMat[1:5,1:5]

dim(GenoMat)

#apply(GenoMat, 1, function(x){sum(is.na(x))})

GenoMatImp<-as.data.frame(GenoMat)
GenoMatImp<-GenoMatImp[, c(1, sample(2:ncol(GenoMatImp), 10000))]

GenoMatImp[is.na(GenoMatImp)]<-1

rownames(GenoMatImp)<-GenoMatImp[,1]
GenoMatImp<-GenoMatImp[,-1]
svdG<-svd(scale(GenoMatImp, scale=FALSE, center=TRUE), nu=5,nv=5)

PCGeno<-as.matrix(GenoMatImp)%*%svdG$v

str(PCGeno)

PCGeno<-cbind(data.frame(Taxa=rownames(PCGeno)), PCGeno)
colnames(PCGeno)[2:6]<-paste("PC",1:5, sep="")
head(PCGeno)

pcs<-PCGeno

# View the first few rows of the data
head(pcs)

dim(GenoMat)

library(dplyr)


dim(filtered_pcs)
excludedTaxa<-c("S101_EKDN220047288-1A_HMVJWDSX5_L2")

filtered_pcs<-filtered_pcs[!filtered_pcs$Taxa%in%excludedTaxa,]
dim(filtered_pcs)

# taxa names start with S1, S2,....S100. But some of them are not unique.
# Get me the counts of each of S1, S2, ...S100 in the taxa column
table(substr(filtered_pcs$Taxa, 1, 4))

# show the ones that have more than one entry
sort(table(substr(filtered_pcs$Taxa, 1, 4))>1)

# Get samples names that have more than one entry
duplicatedTaxa<-names(sort(table(substr(filtered_pcs$Taxa, 1, 3))>1))[sort(table(substr(filtered_pcs$Taxa, 1, 3))>1)]


# get only those rows of the pcs dataframe that have the duplicated taxa, partial matching
# use this list:
#  S12_   S19_ S20_   S23_   S3_   S42_   S57_   S6_   S7_ 
excludeNames<-c("S12_", "S19_", "S20_", "S23_", "S3_", "S42_", "S57_", "S6_", "S7_")
# Get the rows that have the duplicated taxa, partial matching
filtered_pcs_duplicated<-filtered_pcs[grepl(paste(excludeNames, collapse="|"), filtered_pcs$Taxa),]

# get only one of the duplicated taxa and reduce the data to this nonduplicated data
filtered_pcs_nonduplicated<-filtered_pcs[!grepl(paste(excludeNames, collapse="|"), filtered_pcs$Taxa),]
dim(filtered_pcs_nonduplicated)

#Now get one of each of the rows that have the duplicated taxa
for (i in 1:length(excludeNames)){
  filtered_pcs_nonduplicated<-rbind(filtered_pcs_nonduplicated, filtered_pcs_duplicated[grepl(excludeNames[i], filtered_pcs_duplicated$Taxa),][1,])
}

dim(filtered_pcs_nonduplicated)


# plot the first two PCs, clor by the third pc
library(ggplot2)
# add the names of the samples to the plot
ggplot(filtered_pcs_nonduplicated, aes(x=PC1, y=PC2, color=PC3)) + geom_point() + theme_minimal()+
    geom_text(aes(label=Taxa),hjust=0, vjust=0)



# now add a column name in the filtered_pcs_nonduplicated dataframe that has only the S1, S2, S3, ...S100 part of the Taxa column
filtered_pcs_nonduplicated$TaxaShort<-substr(filtered_pcs_nonduplicated$Taxa, 1, 4)

# remove _ and anything that comes after that from the TaxaShort column
filtered_pcs_nonduplicated$TaxaShort<-gsub("_.*", "", filtered_pcs_nonduplicated$TaxaShort)

# sort the data as S1, S2,...
filtered_pcs_nonduplicated<-filtered_pcs_nonduplicated[order(as.numeric(gsub("S", "", filtered_pcs_nonduplicated$TaxaShort))),]

filtered_pcs_nonduplicated



# read this file "/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE/GenotypeTables/isolated\ list\ 2021-2022\ for\ sequencing.xlsx"

# read the file
isolatedlist<-readxl::read_xlsx("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE/GenotypeTables/isolated\ list\ 2021-2022\ for\ sequencing.xlsx")


# have isolate and Name columns in the isolatedlist dataframe
# the name column has the S1, S2, ...S100 names
# now by matching the Name column with the TaxaShort column of the filtered_pcs_nonduplicated dataframe, get the isolate column of the isolatedlist dataframe

filtered_pcs_nonduplicated$Isolate<-isolatedlist$Isolate[match(filtered_pcs_nonduplicated$TaxaShort, isolatedlist$Name)]


# now plot the first two PCs, color by the PC3.

pcaplot<-ggplot(filtered_pcs_nonduplicated, aes(x=PC1, y=PC2, color=PC3)) + geom_point() + theme_minimal()+
    geom_text(aes(label=Isolate),hjust=0, vjust=0)


# save the image as PCAPLOT.png
ggsave("PCAPLOT.png", pcaplot, dpi=600)


filtered_pcs_nonduplicated





# read the file 
library(readxl)

phenotypes<-read.csv("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/PhenoData/20240206_R3E3&R3E5_R1&R2_processed_data.csv", sep=";")
head(phenotypes)
table(phenotypes$Endophyte)
phenotype<-phenotypes[phenotypes$Endophyte%in%"R3E5",]
table(phenotype$Endophyte)
phenotype$Endophyte<-NULL
phenotypes$Zt_Strain[1:10]

library(dplyr)
library(stringr)


filtered_pcs_nonduplicated$Isolate[1:10]

# convert all commas to dots in filtered_pcs_nonduplicated$Isolate
filtered_pcs_nonduplicated$Isolate<-gsub(",", ".", filtered_pcs_nonduplicated$Isolate)

phenotypes$Zt_Strain[1:10]
filtered_pcs_nonduplicated$Isolate[1:10]
###
phenotypes$Zt_Strain <- gsub("-", ".", phenotypes$Zt_Strain)



# read snp_genotypeSummary3.txt data 
SNPSummary<-read.table("6_tassel_analysis/snp_genotypeSummary3.txt", header = TRUE, sep = "\t", as.is = TRUE, skip=0)
SNPSummary[1:5,]
dim(SNPSummary)
MAF<-SNPSummary$Minor.Allele.Proportion
MissingProp<-SNPSummary$Proportion.Missing


# Create a dataframe for MAF and MissingProp
df <- data.frame(MAF = MAF, MissingProp = MissingProp)

# Plot histogram of MAF
p1<-ggplot(df, aes(x = MAF)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Minor Allele Frequency (MAF)",
       x = "MAF",
       y = "Frequency") +
  theme_minimal()

# Plot histogram of MissingProp
p2<-ggplot(df, aes(x = MissingProp)) +
  geom_histogram(binwidth = 0.05, fill = "lightgreen", color = "black") +
  labs(title = "Histogram of Proportion of Missing Data (MissingProp)",
       x = "MissingProp",
       y = "Frequency") +
  theme_minimal()

# Combine the two plots
library(gridExtra)
p3<-grid.arrange(p1, p2, ncol = 2)

# save the plot
ggsave("MAF_MissingProp.png", p3, width = 10, height = 5, units = "in")




# read snp_genotypeSummary4.txt data 
TAXASummary<-read.table("6_tassel_analysis/snp_genotypeSummary4.txt", header = TRUE, sep = "\t", as.is = TRUE, skip=1)

library(data.table)
GenoMat <- fread("6_tassel_analysis/numeric_allele_counts_filtered.txt", skip = 1, header = TRUE, sep = "\t")
#GenoMat[1:5,1:5]
dim(GenoMat)


keepCols <- c(TRUE, (MAF > 0.05 & MissingProp < 0.1))

table(keepCols)

# Subset GenoMat to keep only the desired columns
GenoMat <- GenoMat[, ..keepCols]
dim(GenoMat)

#keepRows<-GenoMat[apply(GenoMat,1, function(x){sum(is.na(x))<3000000}),]
#GenoMat<-GenoMat[..keepRows, ]

GenoMat<-as.data.frame(GenoMat)
GenoMat<-GenoMat[GenoMat[,1]%in%filtered_pcs_nonduplicated$Taxa,]
GenoMat<-GenoMat[match(filtered_pcs_nonduplicated$Taxa, GenoMat[,1]),]
GenoMat$Taxa<-filtered_pcs_nonduplicated$Taxa
GenoMat$Isolate<-filtered_pcs_nonduplicated$Isolate
GenoMat$TaxaShort<-filtered_pcs_nonduplicated$TaxaShort

GenoMat[1:5,1:5]
rownames(GenoMat)<-GenoMat$Isolate


dim(GenoMat)

GenoMat[1:5,1:5]

GenoMat<-GenoMat[,-1]

# read kinship matrix
kinship<-read.table("6_tassel_analysis/kinship_matrix_filtered.txt", header = FALSE, sep = "\t", as.is = TRUE, skip=3)
kinship[1:5,1:5]
rownames(kinship)<-kinship[,1]
kinship<-kinship[,-1]
colnames(kinship)<-rownames(kinship)
# get the kinship matrix for the isolates in the filtered_pcs_nonduplicated dataframe
kinship<-kinship[match(filtered_pcs_nonduplicated$Taxa, rownames(kinship)),]
kinship<-kinship[,match(filtered_pcs_nonduplicated$Taxa, colnames(kinship))]
rownames(kinship)<-filtered_pcs_nonduplicated$Isolate
colnames(kinship)<-filtered_pcs_nonduplicated$Isolate

kinship[1:5,1:5]
kinship<-as.matrix(kinship)

image(kinship)


# make a new column in the phenotypes dataframe that has the Isolate column of the filtered_pcs_nonduplicated dataframe
# this is done by finding the pictures that contain the Isolate names, but the picture names are longer than the Isolate names
# so use partial matching, grepl. cannot use match because the picture names are longer than the Isolate names
# maybe use strinfind

# go over each row of the phenotypes dataframe and find the Isolate name that is in the picture name
# put that name in the phenotype file in a new column
for (i in 1:nrow(phenotypes)){
  for (j in 1:nrow(filtered_pcs_nonduplicated)){
    if (grepl(filtered_pcs_nonduplicated$Isolate[j], phenotypes$Zt_Strain[i])){
      phenotypes$Isolate[i]<-filtered_pcs_nonduplicated$Isolate[j]
    }
  }
}
table(phenotypes$Isolate)

apply(phenotypes, 2, function(x) sum(is.na(x)))
colnames(phenotypes)
phenotypes$Date<-NULL
# boxplot of phenotypes$Area_Endophyte wrt phenotypes$Isolate
ggplot(phenotypes, aes(x=Isolate, y=Area_Endophyte)) + geom_boxplot() + theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))



write.csv(phenotypes, "phenotypes.csv", row.names = FALSE)
write.csv(filtered_pcs_nonduplicated, "filtered_pcs_nonduplicated.csv", row.names = FALSE)
write.csv(GenoMat, "GenoMat.csv", row.names = FALSE)
write.csv(kinship, "kinship.csv", row.names = FALSE)


# Now the data is ready for GWAS
# Gwas for PLACL in phenotypes, leave_id, REP  are fixed effects, also add PCA1, PCA2, PCA3 as fixed effects
colnames(phenotypes)
PhenoGWAS<-phenotypes[,c("Isolate", "BReplicate", "TReplicate","Area_Endophyte", "Area_halo", "Area_inhibition")]


# add the first three PCs as fixed effects
PhenoGWAS<-merge(PhenoGWAS, filtered_pcs_nonduplicated[,c("Isolate", "PC1", "PC2", "PC3")], by="Isolate")
PhenoGWAS[1:5,]

table(PhenoGWAS$BReplicate)
table(PhenoGWAS$TReplicate)



# Now a map file 
SNPMAP<- data.frame(
  rs = SNPSummary$Site.Name,
  chrom = SNPSummary$Chromosome,
  pos = SNPSummary$Physical.Position
)

SNPMAP<-SNPMAP[SNPMAP$rs%in%colnames(GenoMat),]


write.csv(SNPMAP, "SNPMAP.csv", row.names = FALSE)


library(rrBLUP)

geno <- data.frame(marker=SNPMAP$rs,chrom=SNPMAP$chrom,pos=SNPMAP$pos,(t(GenoMat[, SNPMAP$rs])-.5)*2)

geno[1:5,1:5]
dim(geno)


PhenoGWAS[1:5,]
colnames(PhenoGWAS)[1]<-"line"
PhenoGWAS$BReplicate<-as.factor(PhenoGWAS$BReplicate)
PhenoGWAS$TReplicate<-as.factor(PhenoGWAS$TReplicate)
PhenoGWAS$line<-as.factor(PhenoGWAS$line)

PhenoGWAS[1:5,]

# use lmer to fit a model for the phenotype
library(lme4)

kinship[1:5,1:5]
geno[1:5,1:5]
# remove leading X from colnames of geno
colnames(geno)<-gsub("X", "", colnames(geno))
colnames(kinship)<-rownames(kinship)<-gsub("X", "", rownames(kinship))
setdiff(colnames(kinship), colnames(geno))
setdiff(colnames(geno), colnames(kinship))

#  take any extra spaces from the names in all datasets
colnames(kinship)<-rownames(kinship)<-gsub(" ", ".", rownames(kinship))
colnames(geno)<-gsub(" ", ".", colnames(geno))


# remove the columns that are not in both geno and kinship
intersectnames<-intersect(PhenoGWAS$line,intersect(colnames(kinship), colnames(geno)))
kinship<-kinship[,intersectnames]
kinship<-kinship[intersectnames,]
geno<-geno[,c("marker", "chrom", "pos", intersectnames)]


geno<-geno[,c("marker", "chrom", "pos", colnames(kinship))]

dim(kinship)
dim(geno)



library(rrBLUP)
library(lme4)
library(car)

# Assuming PhenoGWAS, GenoMat, filtered_pcs_nonduplicated, kinship, and SNPMAP are already defined

# Define the traits to run GWAS on
traits <- c("Area_inhibition")



# Updated function to generate Manhattan plots using qqman
plot_gwas <- function(GWASOut, GenoMat, trait, annotatePval = 0.005, annotateTop = FALSE) {
  # Prepare the data for qqman
  GWASData <- data.frame(
    CHR = as.numeric(as.factor(GenoMat$chrom)), # Ensure chromosome numbers are numeric
    BP = GenoMat$pos,
    SNP = GenoMat$rs,
    P = 10^(-GWASOut[[trait]]) # Convert to p-values
  )
  
  # Make the Manhattan plot
  qqman::manhattan(GWASData, chr="CHR", bp="BP", snp="SNP", p="P", main=paste("Manhattan Plot for", trait),
            col=c("blue4", "orange3"), annotatePval =annotatePval, annotateTop = annotateTop,suggestiveline = F) # Customize colors as needed
}

# Updated function to generate QQ plots using qqman
plot_qq <- function(GWASOut, trait) {
  # Extract p-values directly from GWASOut
  pvals <- 10^(-GWASOut[[trait]])
  
  # Generate QQ plot
  qqman::qq(pvals, main=paste("QQ Plot for", trait))
}

run_gwas<-function(trait, PhenoGWAS, GenoMat, SNPMAP, kinship, dataName, annotatePval = 0.005, annotateTop = FALSE) {
  # Update PhenoGWAS for the current trait
  PhenoGWAS[,trait] <- as.numeric(PhenoGWAS[,trait])  
  # LMER model
  model <- lmer(as.formula(paste(trait, "~ BReplicate+TReplicate+ (1|line)", sep="")), data=PhenoGWAS)
  ranef <- ranef(model)$line[,1]
  ranefDF <- data.frame(Isolate=rownames(ranef(model)$line), ranef=ranef) 
  
  # Add PCs
  ranefDF <- merge(ranefDF, filtered_pcs_nonduplicated[,c("Isolate", "PC1", "PC2", "PC3")], by="Isolate")
  # GWAS analysis
  GWASOut <- GWAS(pheno=ranefDF[,c("Isolate", "ranef")], geno=GenoMat, fixed = c(), K = NULL, n.PC = 3, 
                  min.MAF=0.00, n.core=1, P3D=TRUE, plot=FALSE)
  
  # change the name from ranef to trait
  colnames(GWASOut)[colnames(GWASOut) == "ranef"] <- trait
  # Save results
  save(GWASOut, file=paste("GWASOut_", trait,"_",dataName, ".RData", sep=""))
  
  # Generate and save plots
  png(paste("GWAS_", trait,"_",dataName, ".png", sep=""))
  plot_gwas(GWASOut, GenoMat, trait, annotatePval, annotateTop)
  dev.off()
  # QQ plot
  png(paste("QQ_", trait,"_",dataName, ".png", sep=""))
  plot_qq(GWASOut, trait)
  dev.off()

}

# impute by mean 
genoMat<-as.matrix(geno[,-c(1:3)])
imputeColbyMean<-function(x){
  x[is.na(x)]<-mean(x, na.rm=TRUE)
  return(x)
}
ImputeGenomatByMean<-function(genoMat){
  genoMat<-apply(genoMat, 2, imputeColbyMean)
  return(genoMat)
}

genoMat<-ImputeGenomatByMean(genoMat)

sum(is.na(genoMat))
genoMat[1:5,1:10]
colnames(genoMat)

genoMat<-genoMat[,apply(genoMat,2,function(x){sum(is.na(x))})==0]

geno<-data.frame(rs=geno$marker, chrom=geno$chrom, pos=geno$pos, genoMat)
str(geno)
dim(geno)
# remove leading X from colnames of geno
colnames(geno)<-gsub("X", "", colnames(geno))

# remove any snps that have single allele
geno<-geno[apply(geno[,4:ncol(geno)], 1, function(x) length(unique(x))>1),]
dim(geno)
geno<-na.omit(geno)
dim(geno)
PhenoGWAS<-PhenoGWAS[PhenoGWAS$line%in%intersectnames,]
dim(PhenoGWAS)
PhenoGWAS$line

remlines<-c("22_ConilVit_L1", "22_CorCris_L1", "22_EcijaSecRic_L1")

PhenoGWAS<-PhenoGWAS[!(PhenoGWAS$line%in%remlines),]


# Main loop to run GWAS for each trait
for(trait in traits) {
  run_gwas(trait, PhenoGWAS, geno, SNPMAP, kinship, "GenoMat_PC3_only_R3E5",annotatePval = 0.1/nrow(geno), annotateTop = FALSE)
}



# haplotype based gwas

library(dplyr)


genoMA <- function(geno, windowsize, stepsize) {
  # Check for necessary columns
  if (!("chrom" %in% names(geno)) || !("pos" %in% names(geno))) {
    stop("geno must contain 'chrom' and 'pos' columns")
  }
  
  # Identify genotype columns, assuming they start with "geno"
  genotype_cols <- names(geno)[-c(1:3)]
  
  geno %>% 
    group_by(chrom) %>%
    mutate(window_start = floor((pos - min(pos)) / stepsize) * stepsize + min(pos),
           window_end = window_start + windowsize) %>%
    group_by(chrom, window_start, window_end) %>%
    summarise(across(all_of(genotype_cols), mean, na.rm = TRUE), .groups = 'drop') %>%
    ungroup()
}


GenoMAout<-genoMA(geno, 5000, 2500)

GenoMAout[1:5,1:5]

dim(GenoMAout)

# replace the window_start window_end columns with one column that is pos (mean of the two)

GenoMAGWAS<-data.frame(
  rs=paste("SNP" , 1:nrow(GenoMAout), sep=""),
  chrom = as.numeric(as.factor(GenoMAout$chrom)), # Ensure chromosome numbers are numeric
  pos = (GenoMAout$window_start+GenoMAout$window_end)/2,
  GenoMAout[,-c(1:3)]
)
dim(GenoMAGWAS)
GenoMAGWAS[1:5,1:5]

# run GWAS for the GenoMAout

SNPMAPMA<- data.frame(
  rs = GenoMAGWAS$rs,
  chrom = GenoMAGWAS$chrom,
  pos = GenoMAGWAS$pos
)
# remove leading X from colnames of geno
colnames(GenoMAGWAS)<-gsub("X", "", colnames(GenoMAGWAS))

# Main loop to run GWAS for each trait
for(trait in traits) {
  run_gwas(trait, PhenoGWAS, GenoMAGWAS, SNPMAPMA, kinship, "GenoMA_PC3_only_R3E3",annotatePval = 0.1/nrow(GenoMAGWAS), annotateTop = FALSE)
}


# PCA based GWAS
genoPC <- function(geno, windowsize, stepsize) {
  # Ensure geno has necessary columns
  if (!("chrom" %in% names(geno)) || !("pos" %in% names(geno))) {
    stop("geno must contain 'chrom' and 'pos' columns")
  }
  
  # Identify genotype columns, assuming they do not include "rs", "chrom", or "pos"
  genotype_cols <- setdiff(names(geno), c("rs", "chrom", "pos"))
  
  # Prepare data with window start and end
  geno$window_start <- ave(geno$pos, geno$chrom, FUN = function(x) floor((x - min(x)) / stepsize) * stepsize + min(x))
  geno$window_end <- geno$window_start + windowsize
  
  # Split data by chromosome, window start, and window end
  geno$group <- with(geno, paste(chrom, window_start, window_end, sep = "_"))
  split_geno <- split(geno, geno$group)
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = length(split_geno), style = 3)
  
  # Function to perform PCA and return the first principal component
  pca_first_pc <- function(df) {
    if(ncol(df) > 0) {
      pr <- prcomp(scale(t(df[, genotype_cols, drop = FALSE]),scale=FALSE, center=TRUE), scale. = FALSE, center = TRUE)
      return(pr$x[, 1])
    } else {
      return(NA)
    }
  }
  
  results_list <- lapply(seq_along(split_geno), function(i) {
    setTxtProgressBar(pb, i)
    chunk <- split_geno[[i]]
    chrom <- unique(chunk$chrom)
    window_start <- unique(chunk$window_start)
    window_end <- unique(chunk$window_end)
    PC1 <- pca_first_pc(chunk)
    return(data.frame(chrom = chrom, window_start = window_start, window_end = window_end, PC1 = t(PC1)))
  })
  
  close(pb)
  
  # Combine results and return
  results_df <- do.call(rbind, results_list)
  return(results_df)
}



# Example usage:
genoPCs <- genoPC(geno, windowsize = 500, stepsize = 250)
genoPCs[1:5,1:5]
dim(genoPCs)

# drop PC1. from the column names
colnames(genoPCs)<-gsub("PC1.", "", colnames(genoPCs))

# combine the genoPCs with the SNPMAPMA
SNPMAPPC<- data.frame(
  rs = paste("SNP" , 1:nrow(genoPCs), sep=""),
  chrom = as.numeric(as.factor(genoPCs$chrom)), # Ensure chromosome numbers are numeric
  pos = (genoPCs$window_start+genoPCs$window_end)/2
)

dim(SNPMAPPC)


GenoPCGWAS<-data.frame(
  rs=SNPMAPPC$rs,
  chrom = as.numeric(as.factor(genoPCs$chrom)), # Ensure chromosome numbers are numeric
  pos = (genoPCs$window_start+genoPCs$window_end)/2,
  genoPCs[,-c(1:3)]
)


dim(SNPMAPPC)
dim(genoPCs)
GenoPCGWAS[1:5,1:5]
# drop leading X from colnames of genoPCs
colnames(GenoPCGWAS)<-gsub("X", "", colnames(GenoPCGWAS))

# Main loop to run GWAS for each trait
for(trait in traits) {
  run_gwas(trait, PhenoGWAS, GenoPCGWAS, SNPMAPPC, kinship, "GenoPC_PC3_only_R3E3",annotatePval = 0.1/nrow(genoPCs), annotateTop = FALSE)
}

filtered_pcs_nonduplicated$Isolate

# Assuming phenotypes is your dataframe and Picture is the column of interest
filtered_pcs_nonduplicated <- filtered_pcs_nonduplicated %>%
  mutate(
    Year = str_extract(Isolate, "^[^_]+"), # Extracts everything before the first underscore
  )

filtered_pcs_nonduplicated$Year<-as.factor(filtered_pcs_nonduplicated$Year)

# principal component analysis using year as color
p<-ggplot(filtered_pcs_nonduplicated, aes(x=PC1, y=PC2, col=as.factor(Year))) + geom_point() + theme_minimal()

ggsave("PCA_year.png", p, width = 10, height = 10, units = "cm")


# heatmap of the kinship matrix
#install.packages("pheatmap")
library(pheatmap)
png("kinship.png")
pheatmap(kinship, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 8, cellwidth = 10, cellheight = 10, main = "Kinship Matrix")
dev.off()

png("kinship.png")
heatmap(kinship, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 8, cellwidth = 10, cellheight = 10, main = "Kinship Matrix")
dev.off()



### thinned geno
sampleGeno<-sample(1:nrow(geno), 30000)
genoThinned<-geno[sampleGeno,]
SNPMAPThinned<-SNPMAP[sampleGeno,]
# Main loop to run GWAS for each trait
for(trait in traits) {
  run_gwas(trait, PhenoGWAS, genoThinned, SNPMAPThinned, kinship, "GenoMatThinned_PC3",annotatePval = 0.1/nrow(geno), annotateTop = FALSE)
}


# a function that will thin every chromosome by position, every 1000bp
thinGeno<-function(geno, SNPMAP, stepsize){
  # thin the geno by stepsize
  genoThinned<-geno[geno$pos%%stepsize==0,]
  SNPMAPThinned<-SNPMAP[geno$pos%%stepsize==0,]
  return(list(genoThinned, SNPMAPThinned))
}

# thin the geno by stepsize
genoThinned<-thinGeno(geno, SNPMAP, 50)

SNPMAPThinned<-genoThinned[[2]]
genoThinned<-genoThinned[[1]]

# Main loop to run GWAS for each trait
for(trait in traits) {
  run_gwas(trait, PhenoGWAS, genoThinned, SNPMAPThinned, kinship, "GenoMatThinned_PC3",annotatePval = 0.1/nrow(geno), annotateTop = FALSE)
}



#########

require(tidyverse)
#install.packages("EMMREML")
#install.packages("pacman")
library("pacman")
pacman::p_load("genetics")
pacman::p_load("scatterplot3d")
pacman::p_load(tidyverse)


# use BiocManager
pacsBIO<-c("snpStats","rtracklayer","GenomicRanges","GenomInfoDb","IRanges")
#BiocManager::install(pacsBIO)
## Install the latest development version from GitHub with
#devtools::install_github("SFUStatgen/LDheatmap")


source("http://zzlab.net/GAPIT/gapit_functions.txt")
pacman::p_load(AssocTests)


# format data for GAPIT
# phenotypes
# we will get the BLUPS for each of the phenotypes for the isolates (taxa) using lmer
PhenoGWAS$line<-as.factor(PhenoGWAS$line)
ModelLMEArea_Endophyte<-lmer(Area_Endophyte~BReplicate+TReplicate+(1|line), data=PhenoGWAS)
ModelLMEArea_halo<-lmer(Area_halo~BReplicate+TReplicate+(1|line), data=PhenoGWAS)
ModelLMEArea_inhibition<-lmer(Area_inhibition~BReplicate+TReplicate+(1|line), data=PhenoGWAS)

# get the BLUPS
BLUPArea_Endophyte<-ranef(ModelLMEArea_Endophyte)$line[,1]
BLUPArea_halo<-ranef(ModelLMEArea_halo)$line[,1]
BLUPArea_inhibition<-ranef(ModelLMEArea_inhibition)$line[,1]

# put these in a dataframe
BLUPS<-data.frame(Isolate=rownames(ranef(ModelLMEArea_inhibition)$line), BLUPArea_inhibition, BLUPArea_halo)
BLUPS[1:5,]
traits
# first column is the isolate name, change this to Taxa
colnames(BLUPS)[1]<-"Taxa"

GenoMat<-geno[,-c(1:3)]
GenoMat<-t(GenoMat)
GenoMat<-as.data.frame(GenoMat)

dim(GenoMat)
# Now GenoMat
GenoMat[1:5,1:5]
# rownames as the first column as TAxa
GenoMat$Taxa<-rownames(GenoMat)
rownames(GenoMat)<-NULL
# make taxa the first column
GenoMat<-GenoMat[,c(ncol(GenoMat), 1:(ncol(GenoMat)-1))]
kinship[1:5,1:5]
kinshipGapit<-kinship
# put the rownames as first column
kinshipGapit<-cbind(rownames(kinshipGapit), as.data.frame(kinshipGapit))
# make taxa the first column
rownames(kinshipGapit)<-colnames(kinshipGapit)<-NULL
kinshipGapit[1:5,1:5]
mapping<-data.frame(Name=geno$rs, Chromosome=geno$chrom, Position=geno$pos)
models = c("FarmCPU","Blink","SUPER")



# impute each column of GenoMat using mean of that column
GenoMat[1:5,1:5]

GenoMat[, 2:ncol(GenoMat)]<-GenoMat[, 2:ncol(GenoMat)]+1
GenoMat[1:5,1:5]
str(GenoMat)
sum(is.na(GenoMat))
modK <- GAPIT(Y=BLUPS,
                GD=GenoMat,
                GM=mapping,
                KI=kinshipGapit,
                CV=NULL,
                PCA.total=0,
                model=models)

modK


# only 13 chromosomes, use mapping
GenoMat[1:5,1:5]
GenoMat[,-1]<-round(GenoMat[,-1])

modNULL <- GAPIT(Y=BLUPS,
                GD=GenoMat,
                GM=mapping,
                KI=NULL,
                CV=NULL,
                PCA.total=3,
                model=models)

modNULL


