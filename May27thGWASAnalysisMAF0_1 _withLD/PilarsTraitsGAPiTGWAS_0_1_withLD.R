# Load necessary data
load("Pheno_geno_Map_PillarTraits.RData")

#######################################
# Create genoImp data frame
genoImp <- data.frame(marker = SNPMAP$rs, chrom = SNPMAP$chrom, pos = SNPMAP$pos, t(genoImp))

#######################################
# Split phenotypes based on 'Removal' criteria
phenotypesSR <- phenotypes[phenotypes$Removal == "Astringent", ]
PhenoGWASSR <- phenotypesSR[, c("Isolate", "BReplicate", "TReplicate", "Endophyte", "Growth_endophyte", "Area_inhibition", "Radius_inhibition", "Ratio_inhibition", "Percentage_inhibition")]

phenotypesLoose <- phenotypes[phenotypes$Removal != "Astringent", ]
PhenoGWASLoose <- phenotypesLoose[, c("Isolate", "BReplicate", "TReplicate", "Endophyte", "Growth_endophyte", "Area_inhibition", "Radius_inhibition", "Ratio_inhibition", "Percentage_inhibition")]

# Calculate new traits: sqrt and log of Area_inhibition
PhenoGWASSR$Sqrt_Area_inhibition <- sqrt(PhenoGWASSR$Area_inhibition)
PhenoGWASSR$Log_Area_inhibition <- log(PhenoGWASSR$Area_inhibition + 1) # Adding 1 to avoid log(0)

PhenoGWASLoose$Sqrt_Area_inhibition <- sqrt(PhenoGWASLoose$Area_inhibition)
PhenoGWASLoose$Log_Area_inhibition <- log(PhenoGWASLoose$Area_inhibition + 1) # Adding 1 to avoid log(0)

# Rename column 'Isolate' to 'line'
colnames(PhenoGWASSR)[1] <- "line"
colnames(PhenoGWASLoose)[1] <- "line"

# Convert specific columns to factors
factor_columns <- c("BReplicate", "TReplicate", "line", "Endophyte")
PhenoGWASSR[factor_columns] <- lapply(PhenoGWASSR[factor_columns], as.factor)
PhenoGWASLoose[factor_columns] <- lapply(PhenoGWASLoose[factor_columns], as.factor)

#######################################
# Filter data based on 'Endophyte'
filter_endophyte <- function(data, endophyte) {
  subset_data <- data[data$Endophyte == endophyte, ]
  na.omit(subset_data)
}

PhenotypeDataList <- list(
  PhenoGWASSRR3E3 = filter_endophyte(PhenoGWASSR, "R3E3"),
  PhenoGWASLooseR3E3 = filter_endophyte(PhenoGWASLoose, "R3E3")
)

#######################################
# Process genotype data
GenoMat <- genoImp[,-c(1:3)]
GenoMat <- t(GenoMat) + 1
GenoMat <- as.data.frame(GenoMat)

# Set rownames as the first column 'Taxa'
GenoMat$Taxa <- rownames(GenoMat)
rownames(GenoMat) <- NULL
GenoMat <- GenoMat[, c(ncol(GenoMat), 1:(ncol(GenoMat) - 1))]

# Clean up 'Taxa' and match with 'pcs2022'
GenoMat$Taxa <- gsub("\\.", "-", GenoMat$Taxa)
GenoMat <- GenoMat[GenoMat$Taxa %in% pcs2022$Taxa, ]
GenoMat <- GenoMat[match(pcs2022$Taxa, GenoMat$Taxa), ]

# Replace 'Taxa' with 'Isolate'
GenoMat$Taxa <- pcs2022$Isolate

# Update PhenotypeDataList to match 'Taxa' levels
PhenotypeDataList <- lapply(PhenotypeDataList, function(x) {
  x$line <- factor(x$line, levels = GenoMat$Taxa)
  x
})

#######################################
# Filter SNPMAP and GenoMat based on LD pruning
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("SNPRelate")
library(SNPRelate)

# Convert to GDS format
snpgdsCreateGeno("geno.gds", genmat = as.matrix(GenoMat[, -1]), sample.id = GenoMat$Taxa, snp.id = colnames(GenoMat)[-1],
                 snp.chromosome = SNPMAP$chrom, snp.position = SNPMAP$pos, snpfirstdim = FALSE)

genofile <- snpgdsOpen("geno.gds")

# LD pruning
set.seed(1000) # For reproducibility
snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.95)
snpset.id <- unlist(snpset)
length(snpset.id)

# Close the GDS file
snpgdsClose(genofile)
file.remove("geno.gds")

# Filter SNPMAP and GenoMat based on LD pruning
SNPMAP <- SNPMAP[SNPMAP$rs %in% snpset.id, ]
GenoMat <- GenoMat[, c("Taxa", SNPMAP$rs)]
dim(GenoMat)
#######################################
# Load external GAPIT functions
source("gapit_functions.txt")

# Ensure phenotype and genotype data are aligned
GenoMat$Taxa <- as.character(GenoMat$Taxa)
GenoMat[1:5, 1:5]


Kmat <- rrBLUP::A.mat(GenoMat[, -1] - 1)
rownames(Kmat) <- colnames(Kmat) <- GenoMat$Taxa
gc()

#######################################
# Run analysis for each phenotype data in PhenotypeDataList
library(lme4)

calculate_blups <- function(data, trait) {
  model <- lmer(formula(paste(trait, "~ BReplicate + TReplicate + (1|line)")), data = data)
  return(ranef(model)$line[, 1])
}

run_analysis <- function(pheno_data) {
  intersectNames <- intersect(pheno_data$line, rownames(Kmat))
  Kmat_sub <- Kmat[intersectNames, intersectNames]
  GenoMat_sub <- GenoMat[GenoMat$Taxa %in% intersectNames, ]
  pheno_data <- pheno_data[pheno_data$line %in% intersectNames, ]
  pheno_data$line <- factor(pheno_data$line, levels = rownames(Kmat_sub))

  traits <- c("Growth_endophyte", "Radius_inhibition", "Area_inhibition", "Ratio_inhibition", "Percentage_inhibition", "Sqrt_Area_inhibition", "Log_Area_inhibition")
  blups <- sapply(traits, calculate_blups, data = pheno_data)

  GenoMat_sub <- as.data.frame(GenoMat_sub)

  BLUPS <- data.frame(Isolate = GenoMat_sub$Taxa, blups)

  list(BLUPS = BLUPS, GenoMat = GenoMat_sub, Kmat = Kmat_sub)
}
results <- lapply(PhenotypeDataList, run_analysis)

#######################################
# Save data and run GAPIT for each result set
new_dir <- "GWAS_resultsForPilarsTraits"
if (!dir.exists(new_dir)) dir.create(new_dir, recursive = TRUE)

initial_dir <- getwd()
pc_values <- 0:5
models <- c("FarmCPU", "Blink")

for (i in seq_along(results)) {
  result <- results[[i]]
  BLUPS <- result$BLUPS
  GenoMat_sub <- result$GenoMat
  Kmat_sub <- result$Kmat
  Kmat_sub <- as.data.frame(Kmat_sub)
  Kmat_sub$Taxa <- rownames(Kmat_sub)
  Kmat_sub <- Kmat_sub[, c(ncol(Kmat_sub), 1:(ncol(Kmat_sub) - 1))]

  phenotype_name <- names(PhenotypeDataList)[i]
  pheno_dir <- file.path(new_dir, phenotype_name)
  if (!dir.exists(pheno_dir)) dir.create(pheno_dir, recursive = TRUE)

  save(BLUPS, GenoMat_sub, SNPMAP, Kmat_sub, file = file.path(pheno_dir, "GWASdataBefGapit.RData"))

  for (pc in pc_values) {
    pc_dir_name <- paste(phenotype_name, "_", pc, "PC", sep = "")
    dir_name <- file.path(pheno_dir, pc_dir_name)

    if (!dir.exists(dir_name)) {
      dir.create(dir_name, recursive = TRUE)
    }

    setwd(dir_name)

    tmp <- capture.output({
      GAPIT(
        Y = BLUPS,
        GD = GenoMat_sub,
        GM = SNPMAP,
        KI = Kmat_sub,
        CV = NULL,
        PCA.total = pc,
        model = models,
        file.output = TRUE,
        cutOff = 0.1
      )
    })

    setwd(initial_dir)
  }
}

cat("GAPIT analysis complete for all principal component settings.\n")
