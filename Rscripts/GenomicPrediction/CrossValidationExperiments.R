setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/GenomicPrediction")

library(ggplot2)
library(tidyr)
library(dplyr)

load("DataforGenomicPrediction.RData")
ls()


# Load BGLR
if (!require(BGLR)) install.packages("BGLR")

Pheno12cm<-as.data.frame(Pheno12cm)
Pheno17cm<-as.data.frame(Pheno17cm)
Pheno12cm$Plant<-factor(Pheno12cm$Plant, levels=rownames(GenoWheat))
Pheno17cm$Plant<-factor(Pheno17cm$Plant, levels=rownames(GenoWheat))
Pheno12cm$Strain<-factor(Pheno12cm$Strain, levels=rownames(KmatSepMix))
Pheno17cm$Strain<-factor(Pheno17cm$Strain, levels=rownames(KmatSepMix))

sum(is.na(Pheno12cm))
sum(is.na(Pheno17cm))

# leafArea, necrosisArea, PLACL, pycnidiaCount, meanPycnidiaArea
# pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion, pycnidiaGreyValue are all numeric. 

NumVars<-c("leafArea", "necrosisArea", "PLACL", "pycnidiaCount", "meanPycnidiaArea", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion", "pycnidiaGreyValue")
Pheno12cm[NumVars]<-lapply(Pheno12cm[NumVars], as.numeric)
Pheno17cm[NumVars]<-lapply(Pheno17cm[NumVars], as.numeric)

# use missRanger to impute missing values
if (!require(missRanger)) install.packages("missRanger")

Pheno12cmImp<-missRanger(Pheno12cm, num.trees=100)
Pheno17cmImp<-missRanger(Pheno17cm, num.trees=100)

head(Pheno12cmImp)

str(Pheno12cmImp)


# aggregate over reps, leaf and set
Pheno12cmImpAgg<-aggregate(Pheno12cmImp[,NumVars], by=list(Plant=Pheno12cmImp$Plant, Strain=Pheno12cmImp$Strain), FUN=mean)

Pheno17cmImpAgg<-aggregate(Pheno17cmImp[,NumVars], by=list(Plant=Pheno17cmImp$Plant, Strain=Pheno17cmImp$Strain), FUN=mean)
dim(Pheno12cmImpAgg)
dim(Pheno17cmImpAgg)

hist(Pheno17cmImpAgg$PLACL)
############################################

# we will define cross validation startegies.
# the experiment consists of 4 mixes and 2000 wheat lines
# we would like to do
# 1- leave 20% of the wheat lines out for validation
# 2- Leave one mix out for validation
# 3- Leave one mix out for validation and 20% of the wheat lines out for validation
# this will be repeated 10 times for each strategy

KmatSepAll<-KmatSepMix
KmatSepMix<-KmatSepAll[c("mix1","mix2","mix3","mix4"),c("mix1","mix2","mix3","mix4") ]

KmatSepMixDiag<-diag(nrow(KmatSepMix))
rownames(KmatSepMixDiag)<-colnames(KmatSepMixDiag)<-rownames(KmatSepMix)
KmatWheatDiag<-diag(nrow(KmatWheat))
rownames(KmatWheatDiag)<-colnames(KmatWheatDiag)<-rownames(KmatWheat)

###################################################################



runGenomicPredictionStrategy <- function(traitName, kinshipMatrixWheat, kinshipMatrixMix, phenotypeData12cm, phenotypeData17cm) {
    wheatLines <- rownames(kinshipMatrixWheat)
    trainWheatLines <- sample(wheatLines, 0.8 * length(wheatLines))
    testWheatLines <- setdiff(wheatLines, trainWheatLines)

    phenotypeTrain12 <- phenotypeData12cm
    phenotypeTrain17 <- phenotypeData17cm
    phenotypeTrain12[phenotypeTrain12$Plant %in% testWheatLines, traitName] <- NA
    phenotypeTrain17[phenotypeTrain17$Plant %in% testWheatLines, traitName] <- NA
    print(sum(is.na(phenotypeTrain12[phenotypeTrain12$Plant %in% testWheatLines, traitName])))
    ZtrainWheat12 <- model.matrix(~0 + Plant, data = phenotypeTrain12)
    ZtrainWheat17 <- model.matrix(~0 + Plant, data = phenotypeTrain17)
    XtrainWheat12 <- model.matrix(~1, data = phenotypeTrain12)
    XtrainWheat17 <- model.matrix(~1, data = phenotypeTrain17)
    ZTrainMix12 <- model.matrix(~0 + Strain, data = phenotypeTrain12)
    ZTrainMix17 <- model.matrix(~0 + Strain, data = phenotypeTrain17)

    fm12 <- BGLR(y = as.numeric(phenotypeTrain12[[traitName]]), ETA = list(
        list(X = XtrainWheat12, model = "FIXED"),
        list(K = ZtrainWheat12 %*% kinshipMatrixWheat %*% t(ZtrainWheat12), model = "RKHS"),
        list(K = ZTrainMix12 %*% kinshipMatrixMix %*% t(ZTrainMix12), model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    fm17 <- BGLR(y = as.numeric(phenotypeTrain17[[traitName]]), ETA = list(
        list(X = XtrainWheat17, model = "FIXED"),
        list(K = ZtrainWheat17 %*% kinshipMatrixWheat %*% t(ZtrainWheat17), model = "RKHS"),
        list(K = ZTrainMix17 %*% kinshipMatrixMix %*% t(ZTrainMix17), model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    K12_mix <- ZTrainMix12 %*% kinshipMatrixMix %*% t(ZTrainMix12)
    K12_wheat <- ZtrainWheat12 %*% kinshipMatrixWheat %*% t(ZtrainWheat12)
    K12_combined <- K12_mix * K12_wheat

    K17_mix <- ZTrainMix17 %*% kinshipMatrixMix %*% t(ZTrainMix17)
    K17_wheat <- ZtrainWheat17 %*% kinshipMatrixWheat %*% t(ZtrainWheat17)
    K17_combined <- K17_mix * K17_wheat
    fm12_interaction <- BGLR(y = as.numeric(phenotypeTrain12[[traitName]]), ETA = list(
        list(X = XtrainWheat12, model = "FIXED"),
        list(K = K12_mix, model = "RKHS"),
        list(K = K12_wheat, model = "RKHS"),
        list(K = K12_combined, model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    fm17_interaction <- BGLR(y = as.numeric(phenotypeTrain17[[traitName]]), ETA = list(
        list(X = XtrainWheat17, model = "FIXED"),
        list(K = K17_mix, model = "RKHS"),
        list(K = K17_wheat, model = "RKHS"),
        list(K = K17_combined, model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)
    return(list(fm12 = fm12, fm17 = fm17, fm12_interaction = fm12_interaction, fm17_interaction = fm17_interaction, TestWheat = testWheatLines))
}




evaluateModelPredictions <- function(models, phenotypeData12, phenotypeData17, traitName) {
    # Generate predictions from models
    predictions12 <- predict(models$fm12)
    predictions17 <- predict(models$fm17)
    predictions12_interaction <- predict(models$fm12_interaction)
    predictions17_interaction <- predict(models$fm17_interaction)
    
    # Identify test locations
    testLocations12 <- phenotypeData12$Plant %in% models$TestWheat
    testLocations17 <- phenotypeData17$Plant %in% models$TestWheat
    
    # Extract test data
    phenotypeTest12 <- phenotypeData12[testLocations12,]
    phenotypeTest17 <- phenotypeData17[testLocations17,]
    
    # Calculate correlations
    correlation12 <- cor(predictions12[testLocations12], phenotypeTest12[[traitName]], use = "complete.obs")
    correlation17 <- cor(predictions17[testLocations17], phenotypeTest17[[traitName]], use = "complete.obs")
    correlation12_interaction <- cor(predictions12_interaction[testLocations12], phenotypeTest12[[traitName]], use = "complete.obs")
    correlation17_interaction <- cor(predictions17_interaction[testLocations17], phenotypeTest17[[traitName]], use = "complete.obs")

    # Calculate within-strain correlations for both interaction and non-interaction models
    withinStrainCorrelations <- function(predictions, phenotypeData, testLocations) {
        data <- data.frame(predictions = predictions[testLocations],
                           traits = phenotypeData[[traitName]][testLocations],
                           Strain = phenotypeData$Strain[testLocations])
        
        listData <- split(data, data$Strain)
        lapply(listData, function(df) cor(df$predictions, df$traits, use = "complete.obs"))
    }
    
    Cor12withinStrain <- withinStrainCorrelations(predictions12, phenotypeData12, testLocations12)
    Cor17withinStrain <- withinStrainCorrelations(predictions17, phenotypeData17, testLocations17)
    Cor12withinStrain_interaction <- withinStrainCorrelations(predictions12_interaction, phenotypeData12, testLocations12)
    Cor17withinStrain_interaction <- withinStrainCorrelations(predictions17_interaction, phenotypeData17, testLocations17)

    # Return results
    correlationResults <- list(
        Correlation12 = correlation12,
        Correlation17 = correlation17,
        Correlation12Interaction = correlation12_interaction,
        Correlation17Interaction = correlation17_interaction,
        Cor12withinStrain = Cor12withinStrain,
        Cor17withinStrain = Cor17withinStrain,
        Cor12withinStrain_interaction = Cor12withinStrain_interaction,
        Cor17withinStrain_interaction = Cor17withinStrain_interaction
    )

    return(list(CorrelationResults = correlationResults))
}


ST1Models<-runGenomicPredictionStrategy("PLACL", KmatWheatDiag, KmatSepMixDiag, Pheno12cmImpAgg, Pheno17cmImpAgg)
outcor<-evaluateModelPredictions(ST1Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")

outcor$CorrelationResults$Correlation12
outcor$CorrelationResults$Correlation17
outcor$CorrelationResults$Cor12withinStrain
outcor$CorrelationResults$Cor17withinStrain
outcor$CorrelationResults$Cor12withinStrain_interaction
outcor$CorrelationResults$Cor17withinStrain_interaction

CorResultsDFwithinStrains<-rbind(outcor$CorrelationResults$Cor12withinStrain, outcor$CorrelationResults$Cor17withinStrain, outcor$CorrelationResults$Cor12withinStrain_interaction, outcor$CorrelationResults$Cor17withinStrain_interaction)
rownames(CorResultsDFwithinStrains)<-c("12", "17", "12_interaction", "17_interaction")
CorResultsDFwithinStrains



ST1Models<-runGenomicPredictionStrategy("PLACL", KmatWheat, KmatSepMix, Pheno12cmImpAgg, Pheno17cmImpAgg)

outcor<-evaluateModelPredictions(ST1Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")


CorResultsDFwithinStrains<-rbind(outcor$CorrelationResults$Cor12withinStrain, outcor$CorrelationResults$Cor17withinStrain, outcor$CorrelationResults$Cor12withinStrain_interaction, outcor$CorrelationResults$Cor17withinStrain_interaction)
rownames(CorResultsDFwithinStrains)<-c("12", "17", "12_interaction", "17_interaction")
CorResultsDFwithinStrains



# for loop to repeat the experiment 10 times

resultsDiag<-list()
resultsKmat<-list()
for(i in 1:10){
    print(i)
    ST1ModelsDiag<-runGenomicPredictionStrategy("PLACL", KmatWheatDiag, KmatSepMixDiag, Pheno12cmImpAgg, Pheno17cmImpAgg)
    outcorDiag<-evaluateModelPredictions(ST1Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
    resultsDiag[[i]]<-outcorDiag
    ST1ModelsKmat<-runGenomicPredictionStrategy("PLACL", KmatWheat, KmatSepMix, Pheno12cmImpAgg, Pheno17cmImpAgg)
    outcorKmat<-evaluateModelPredictions(ST1Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
    resultsKmat[[i]]<-outcorKmat
}

# save results
save(resultsDiag, resultsKmat, file="resultsDiagKmat.RData")

###########Scenario2: Leave an entire mix out for validation

runGenomicPredictionStrategy2 <- function(traitName, kinshipMatrixWheat, kinshipMatrixMix, phenotypeData12cm, phenotypeData17cm, testMix) {
    phenotypeTrain12 <- phenotypeData12cm
    phenotypeTrain17 <- phenotypeData17cm
    phenotypeTrain12[phenotypeTrain12$Strain %in% testMix, traitName] <- NA
    phenotypeTrain17[phenotypeTrain17$Strain %in% testMix, traitName] <- NA
    print(sum(is.na(phenotypeTrain12[phenotypeTrain12$Strain %in% testMix, traitName])))
    ZtrainWheat12 <- model.matrix(~0 + Plant, data = phenotypeTrain12)
    ZtrainWheat17 <- model.matrix(~0 + Plant, data = phenotypeTrain17)
    XtrainWheat12 <- model.matrix(~1, data = phenotypeTrain12)
    XtrainWheat17 <- model.matrix(~1, data = phenotypeTrain17)
    ZTrainMix12 <- model.matrix(~0 + Strain, data = phenotypeTrain12)
    ZTrainMix17 <- model.matrix(~0 + Strain, data = phenotypeTrain17)

    fm12 <- BGLR(y = as.numeric(phenotypeTrain12[[traitName]]), ETA = list(
        list(X = XtrainWheat12, model = "FIXED"),
        list(K = ZtrainWheat12 %*% kinshipMatrixWheat %*% t(ZtrainWheat12), model = "RKHS"),
        list(K = ZTrainMix12 %*% kinshipMatrixMix %*% t(ZTrainMix12), model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    fm17 <- BGLR(y = as.numeric(phenotypeTrain17[[traitName]]), ETA = list(
        list(X = XtrainWheat17, model = "FIXED"),
        list(K = ZtrainWheat17 %*% kinshipMatrixWheat %*% t(ZtrainWheat17), model = "RKHS"),
        list(K = ZTrainMix17 %*% kinshipMatrixMix %*% t(ZTrainMix17), model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    K12_mix <- ZTrainMix12 %*% kinshipMatrixMix %*% t(ZTrainMix12)
    K12_wheat <- ZtrainWheat12 %*% kinshipMatrixWheat %*% t(ZtrainWheat12)
    K12_combined <- K12_mix * K12_wheat

    K17_mix <- ZTrainMix17 %*% kinshipMatrixMix %*% t(ZTrainMix17)
    K17_wheat <- ZtrainWheat17 %*% kinshipMatrixWheat %*% t(ZtrainWheat17)
    K17_combined <- K17_mix * K17_wheat
    fm12_interaction <- BGLR(y = as.numeric(phenotypeTrain12[[traitName]]), ETA = list(
        list(X = XtrainWheat12, model = "FIXED"),
        list(K = K12_mix, model = "RKHS"),
        list(K = K12_wheat, model = "RKHS"),
        list(K = K12_combined, model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    fm17_interaction <- BGLR(y = as.numeric(phenotypeTrain17[[traitName]]), ETA = list(
        list(X = XtrainWheat17, model = "FIXED"),
        list(K = K17_mix, model = "RKHS"),
        list(K = K17_wheat, model = "RKHS"),
        list(K = K17_combined, model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)
    return(list(fm12 = fm12, fm17 = fm17, fm12_interaction = fm12_interaction, fm17_interaction = fm17_interaction, TestMix = testMix))
}


evaluateModelPredictions2 <- function(models, phenotypeData12, phenotypeData17, traitName) {
    predictions12 <- predict(models$fm12)
    predictions17 <- predict(models$fm17)
    predictions12_interaction <- predict(models$fm12_interaction)
    predictions17_interaction <- predict(models$fm17_interaction)
    testLocations12 <- phenotypeData12$Strain %in% models$TestMix
    testLocations17 <- phenotypeData17$Strain %in% models$TestMix

    phenotypeTest12 <- phenotypeData12[testLocations12,]
    phenotypeTest17 <- phenotypeData17[testLocations17,]
    phenotypeTest12Trait<-phenotypeTest12[[traitName]]
    phenotypeTest17Trait<-phenotypeTest17[[traitName]]
    correlation12 <- cor(predictions12[testLocations12], phenotypeTest12Trait)
    correlation17 <- cor(predictions17[testLocations17], phenotypeTest17Trait)
    correlation12_interaction <- cor(predictions12_interaction[testLocations12], phenotypeTest12Trait)
    correlation17_interaction <- cor(predictions17_interaction[testLocations17], phenotypeTest17Trait)
    
    # Plotting results
    plot(predictions12_interaction[testLocations12], phenotypeTest12Trait, col = phenotypeTest12$Plant)
    plot(predictions17_interaction[testLocations17], phenotypeTest17Trait, col = phenotypeTest17$Plant)

    correlationResults <- list(
        Correlation12 = correlation12, 
        Correlation17 = correlation17, 
        Correlation12Interaction = correlation12_interaction, 
        Correlation17Interaction = correlation17_interaction
    )

    return(list(CorrelationResults = correlationResults))
}


ST2Models<-runGenomicPredictionStrategy2("PLACL", KmatWheatDiag, KmatSepMixDiag, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix1")
outcor<-evaluateModelPredictions2(ST2Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
ST2Models<-runGenomicPredictionStrategy2("PLACL", KmatWheatDiag, KmatSepMixDiag, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix2")
outcor<-evaluateModelPredictions2(ST2Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
ST2Models<-runGenomicPredictionStrategy2("PLACL", KmatWheatDiag, KmatSepMixDiag, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix3")
outcor<-evaluateModelPredictions2(ST2Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")



######## scenario 3: leave one mix out and 20% of the wheat lines out for validation

runGenomicPredictionStrategy3 <- function(traitName, kinshipMatrixWheat, kinshipMatrixMix, phenotypeData12cm, phenotypeData17cm, testMix, testWheat) {
    phenotypeTrain12 <- phenotypeData12cm
    phenotypeTrain17 <- phenotypeData17cm
    phenotypeTrain12[phenotypeTrain12$Strain %in% testMix, traitName] <- NA
    phenotypeTrain17[phenotypeTrain17$Strain %in% testMix, traitName] <- NA
    phenotypeTrain12[phenotypeTrain12$Plant %in% testWheat, traitName] <- NA
    phenotypeTrain17[phenotypeTrain17$Plant %in% testWheat, traitName] <- NA
    print(sum(is.na(phenotypeTrain12[phenotypeTrain12$Strain %in% testMix, traitName])))
    ZtrainWheat12 <- model.matrix(~0 + Plant, data = phenotypeTrain12)
    ZtrainWheat17 <- model.matrix(~0 + Plant, data = phenotypeTrain17)
    XtrainWheat12 <- model.matrix(~1, data = phenotypeTrain12)
    XtrainWheat17 <- model.matrix(~1, data = phenotypeTrain17)
    ZTrainMix12 <- model.matrix(~0 + Strain, data = phenotypeTrain12)
    ZTrainMix17 <- model.matrix(~0 + Strain, data = phenotypeTrain17)

    fm12 <- BGLR(y = as.numeric(phenotypeTrain12[[traitName]]), ETA = list(
        list(X = XtrainWheat12, model = "FIXED"),
        list(K = ZtrainWheat12 %*% kinshipMatrixWheat %*% t(ZtrainWheat12), model = "RKHS"),
        list(K = ZTrainMix12 %*% kinshipMatrixMix %*% t(ZTrainMix12), model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    fm17 <- BGLR(y = as.numeric(phenotypeTrain17[[traitName]]), ETA = list(
        list(X = XtrainWheat17, model = "FIXED"),
        list(K = ZtrainWheat17 %*% kinshipMatrixWheat %*% t(ZtrainWheat17), model = "RKHS"),
        list(K = ZTrainMix17 %*% kinshipMatrixMix %*% t(ZTrainMix17), model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    K12_mix <- ZTrainMix12 %*% kinshipMatrixMix %*% t(ZTrainMix12)
    K12_wheat <- ZtrainWheat12 %*% kinshipMatrixWheat %*% t(ZtrainWheat12)
    K12_combined <- K12_mix * K12_wheat

    K17_mix <- ZTrainMix17 %*% kinshipMatrixMix %*% t(ZTrainMix17)
    K17_wheat <- ZtrainWheat17 %*% kinshipMatrixWheat %*% t(ZtrainWheat17)
    K17_combined <- K17_mix * K17_wheat
    fm12_interaction <- BGLR(y = as.numeric(phenotypeTrain12[[traitName]]), ETA = list(
        list(X = XtrainWheat12, model = "FIXED"),
        list(K = K12_mix, model = "RKHS"),
        list(K = K12_wheat, model = "RKHS"),
        list(K = K12_combined, model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    fm17_interaction <- BGLR(y = as.numeric(phenotypeTrain17[[traitName]]), ETA = list(
        list(X = XtrainWheat17, model = "FIXED"),
        list(K = K17_mix, model = "RKHS"),
        list(K = K17_wheat, model = "RKHS"),
        list(K = K17_combined, model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    return(list(fm12 = fm12, fm17 = fm17, fm12_interaction = fm12_interaction, fm17_interaction = fm17_interaction, TestMix = testMix, TestWheat = testWheat))
}

evaluateModelPredictions3<-function(models, phenotypeData12, phenotypeData17, traitName){
    predictions12 <- predict(models$fm12)
    predictions17 <- predict(models$fm17)
    predictions12_interaction <- predict(models$fm12_interaction)
    predictions17_interaction <- predict(models$fm17_interaction)
    testLocations12 <- phenotypeData12$Strain %in% models$TestMix
    testLocations17 <- phenotypeData17$Strain %in% models$TestMix
    testLocations12Wheat<-phenotypeData12$Plant %in% models$TestWheat
    testLocations17Wheat<-phenotypeData17$Plant %in% models$TestWheat


    testLocations12<-testLocations12 & testLocations12Wheat
    testLocations17<-testLocations17 & testLocations17Wheat

    phenotypeTest12 <- phenotypeData12[testLocations12,]
    phenotypeTest17 <- phenotypeData17[testLocations17,]
    phenotypeTest12Trait<-phenotypeTest12[[traitName]]
    phenotypeTest17Trait<-phenotypeTest17[[traitName]]
    correlation12 <- cor(predictions12[testLocations12], phenotypeTest12Trait)
    correlation17 <- cor(predictions17[testLocations17], phenotypeTest17Trait)
    correlation12_interaction <- cor(predictions12_interaction[testLocations12], phenotypeTest12Trait)
    correlation17_interaction <- cor(predictions17_interaction[testLocations17], phenotypeTest17Trait)
    
    # Plotting results
    plot(predictions12_interaction[testLocations12], phenotypeTest12Trait, col = phenotypeTest12$Plant)
    plot(predictions17_interaction[testLocations17], phenotypeTest17Trait, col = phenotypeTest17$Plant)

    correlationResults <- list(
        Correlation12 = correlation12, 
        Correlation17 = correlation17, 
        Correlation12Interaction = correlation12_interaction, 
        Correlation17Interaction = correlation17_interaction
    )

    return(list(CorrelationResults = correlationResults))
}
testWheat<-sample(rownames(KmatWheat), ceiling(0.2*nrow(KmatWheat)))
TS3Models<-runGenomicPredictionStrategy3("PLACL", KmatWheatDiag, KmatSepMixDiag, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix1", testWheat = testWheat)

outcor<-evaluateModelPredictions3(TS3Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
outcor


TS3Models<-runGenomicPredictionStrategy3("PLACL", KmatWheat, KmatSepMix, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix1", testWheat = testWheat)

outcor<-evaluateModelPredictions3(TS3Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
outcor

#############################################

# for loop to repeat the experiment 10 times

# Define the kinship matrices to be used
kinshipWheat <- list(Normal = KmatWheat, Diagonal = KmatWheatDiag)
kinshipMix <- list(Normal = KmatSepMix, Diagonal = KmatSepMixDiag)

# Store results in a structured way
allResults <- list()

# Loop over 10 iterations
for (iteration in 1:10) {
  # Sample 20% of wheat lines randomly
  testWheatLines <- sample(rownames(KmatWheat), ceiling(0.2 * nrow(KmatWheat)))
  
  # Loop over each combination of kinship matrices
  for (wheatKey in names(kinshipWheat)) {
    for (mixKey in names(kinshipMix)) {
      # Store results for each iteration and matrix combination
      scenarioResults <- list()

      # Scenario 1: Leave 20% of wheat lines out for validation
      ST1Models <- runGenomicPredictionStrategy("PLACL", kinshipWheat[[wheatKey]], kinshipMix[[mixKey]], Pheno12cmImpAgg, Pheno17cmImpAgg)
      scenarioResults$Scenario1 <- evaluateModelPredictions(ST1Models, Pheno12cmImpAgg, Pheno17cmImpAgg, "PLACL")
      
      # Scenario 2: Leave one mix out for validation
      for (mix in unique(Pheno12cmImpAgg$Strain)) {
        ST2Models <- runGenomicPredictionStrategy2("PLACL", kinshipWheat[[wheatKey]], kinshipMix[[mixKey]], Pheno12cmImpAgg, Pheno17cmImpAgg, mix)
        scenarioResults[[paste("Scenario2", mix)]] <- evaluateModelPredictions2(ST2Models, Pheno12cmImpAgg, Pheno17cmImpAgg, "PLACL")
      }
      
      # Scenario 3: Leave one mix out and 20% of the wheat lines out for validation
      
      for (mix in unique(Pheno12cmImpAgg$Strain)) {
        ST3Models <- runGenomicPredictionStrategy3("PLACL", kinshipWheat[[wheatKey]], kinshipMix[[mixKey]], Pheno12cmImpAgg, Pheno17cmImpAgg, mix, testWheatLines)
        scenarioResults[[paste("Scenario3", mix)]] <- evaluateModelPredictions3(ST3Models, Pheno12cmImpAgg, Pheno17cmImpAgg, "PLACL")
      }
      
      # Store all results for this combination of matrices
      allResults[[paste("Iteration", iteration, "Wheat", wheatKey, "Mix", mixKey)]] <- scenarioResults
    }
  }
}

# Save all results to a file
save(allResults, file = "AllGenomicPredictionResults.RData")


################## Results Gathering
# Load the results
setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/GenomicPrediction")

load("AllGenomicPredictionResults.RData")


# Extract the correlation results for Scenario 1
# only within mixes

correlationsScenario1 <- sapply(allResults, function(iterationResults) {
  scenario1Results <- iterationResults$Scenario1
  corwithinStrains<-scenario1Results$CorrelationResults
  return(corwithinStrains)
})
dim(correlationsScenario1)
resultsScenario1_12<-correlationsScenario1[5,]
resultsScenario1_17<-correlationsScenario1[6,]
resultsScenario1_12_int<-correlationsScenario1[7,]
resultsScenario1_17_int<-correlationsScenario1[8,]



resultsScenario1_12_1<-sapply(resultsScenario1_12, function(x) x)
resultsScenario1_17_1<-sapply(resultsScenario1_17, function(x) x)
resultsScenario1_12_1_int<-sapply(resultsScenario1_12_int, function(x) x)
resultsScenario1_17_1_int<-sapply(resultsScenario1_17_int, function(x) x)
resultsScenario1_12_1[1:4,1:6]

# Assuming resultsScenario1_12_1 is already loaded in your environment
# Let's reshape this into a 3D array.

# Number of mixes, methods, and iterations
num_mixes <- 4
num_methods <- 4
num_iterations <- 10

# We need to extract data for each iteration and method combination correctly.
# Reshape the data from 2D matrix to 3D array
array_data_12 <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_17 <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_12_int <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_17_int <- array(dim = c(num_mixes, num_methods, num_iterations))
# Filling the array
for (iteration in 1:num_iterations) {
  for (method in 1:num_methods) {
    # Calculate column indices for the current iteration and method
    start_col <- (iteration - 1) * num_methods + method
    array_data_12[, method, iteration] <- unlist(resultsScenario1_12_1[, start_col])
    array_data_17[, method, iteration] <- unlist(resultsScenario1_17_1[, start_col])
    array_data_12_int[, method, iteration] <- unlist(resultsScenario1_12_1_int[, start_col])
    array_data_17_int[, method, iteration] <- unlist(resultsScenario1_17_1_int[, start_col])
  }
}

# Adding dimension names for clarity
dimnames(array_data_12)<- dimnames(array_data_17) <- list(
  Mix = paste0("Mix ", 1:num_mixes),
  Method = c("Normal Mix Normal", "Normal Mix Diagonal", "Diagonal Mix Normal", "Diagonal Mix Diagonal"),
  Iteration = paste0("Iteration ", 1:num_iterations)
)

dimnames(array_data_12_int)<- dimnames(array_data_17_int) <- list(
  Mix = paste0("Mix ", 1:num_mixes),
  Method = c("Normal Mix Normal", "Normal Mix Diagonal", "Diagonal Mix Normal", "Diagonal Mix Diagonal"),
  Iteration = paste0("Iteration ", 1:num_iterations)
)

# Display the structure of the created 3D array
print(array_data_12)

dimnames(array_data_12)

# means for each method and mix
means12Scenario1<-apply(array_data_12, c(1,2), function(x){mean(x, na.rm = TRUE)})
means17Scenario1<-apply(array_data_17, c(1,2), function(x){mean(x, na.rm = TRUE)})
means12Scenario1_int<-apply(array_data_12_int, c(1,2), function(x){mean(x, na.rm = TRUE)})
means17Scenario1_int<-apply(array_data_17_int, c(1,2), function(x){mean(x, na.rm = TRUE)})


se12Scenario1<-apply(array_data_12, c(1,2), function(x){sd(x, na.rm = TRUE)/sqrt(length(x))})   
se17Scenario1<-apply(array_data_17, c(1,2), function(x){sd(x, na.rm = TRUE)/sqrt(length(x))})
se12Scenario1_int<-apply(array_data_12_int, c(1,2), function(x){sd(x, na.rm = TRUE)/sqrt(length(x))})
se17Scenario1_int<-apply(array_data_17_int, c(1,2), function(x){sd(x, na.rm = TRUE)/sqrt(length(x))})

# barplot of means with error bars
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)
means12Scenario1<-as.data.frame(means12Scenario1)
means17Scenario1<-as.data.frame(means17Scenario1)
means12Scenario1_int<-as.data.frame(means12Scenario1_int)
means17Scenario1_int<-as.data.frame(means17Scenario1_int)
means12Scenario1$Mix<-rownames(means12Scenario1)
means17Scenario1$Mix<-rownames(means17Scenario1)
means12Scenario1_int$Mix<-rownames(means12Scenario1_int)
means17Scenario1_int$Mix<-rownames(means17Scenario1_int)

# Reshape the means data frame from wide to long format
means_long_12 <- means12Scenario1 %>%
  pivot_longer(cols = -Mix, names_to = "Method", values_to = "Mean")
means_long_17 <- means17Scenario1 %>%
    pivot_longer(cols = -Mix, names_to = "Method", values_to = "Mean")
means_long_12_int <- means12Scenario1_int %>%
  pivot_longer(cols = -Mix, names_to = "Method", values_to = "Mean")  
means_long_17_int <- means17Scenario1_int %>% 
  pivot_longer(cols = -Mix, names_to = "Method", values_to = "Mean")  


se12Scenario1<-as.data.frame(se12Scenario1)
se12Scenario1$Mix<-rownames(se12Scenario1)
se17Scenario1<-as.data.frame(se17Scenario1)
se17Scenario1$Mix<-rownames(se17Scenario1)
se12Scenario1_int<-as.data.frame(se12Scenario1_int)
se12Scenario1_int$Mix<-rownames(se12Scenario1_int)
se17Scenario1_int<-as.data.frame(se17Scenario1_int)
se17Scenario1_int$Mix<-rownames(se17Scenario1_int)



sd_long_12 <- se12Scenario1 %>%
  pivot_longer(cols = -Mix, names_to = "Method", values_to = "SE")
sd_long_17 <- se17Scenario1 %>%
  pivot_longer(cols = -Mix, names_to = "Method", values_to = "SE")
sd_long_12_int <- se12Scenario1_int %>%
  pivot_longer(cols = -Mix, names_to = "Method", values_to = "SE")
sd_long_17_int <- se17Scenario1_int %>%
  pivot_longer(cols = -Mix, names_to = "Method", values_to = "SE")



# Combine the means and standard errors into a single data frame
means_long_12$SE <- sd_long_12$SE
means_long_17$SE <- sd_long_17$SE
means_long_12_int$SE <- sd_long_12_int$SE
means_long_17_int$SE <- sd_long_17_int$SE

# Plot the means with error bars
p1_12<-ggplot(means_long_12, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()

p1_17<-ggplot(means_long_17, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()

p1_int_12<-ggplot(means_long_12_int, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()

p1_int_17<-ggplot(means_long_17_int, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()


##################

correlationsScenario2 <- lapply(allResults, function(iterationResults) {
  scenario2Results <- iterationResults[grep("Scenario2", names(iterationResults))][1:4]
  return(scenario2Results)
})


# Create an array to store the correlation results for Scenario 2
num_mixes <- 4
num_methods_wheat <- 2
num_methods_mix <- 2
num_methods <- 4
num_iterations <- 10

array_data_12_scenario2 <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_17_scenario2 <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_12_int_scenario2 <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_17_int_scenario2 <- array(dim = c(num_mixes, num_methods, num_iterations))

# Fill the arrays with correlation results
for (iteration in 1:num_iterations) {
  method=0
  for (method_wheat in 1:num_methods_wheat) {
    for (method_mix in 1:num_methods_mix) {
      method=method+1
    for (mix in 1:num_mixes) {
      cor_results <- correlationsScenario2[[paste("Iteration", iteration, "Wheat", names(kinshipWheat)[method_wheat], "Mix", names(kinshipMix)[method_mix])]][[paste("Scenario2 mix", mix, sep = "")]]$CorrelationResults
      array_data_12_scenario2[mix, method, iteration] <- cor_results$Correlation12
      array_data_17_scenario2[mix, method,iteration] <- cor_results$Correlation17
      array_data_12_int_scenario2[mix, method, iteration] <- cor_results$Correlation12Interaction
      array_data_17_int_scenario2[mix, method, iteration] <- cor_results$Correlation17Interaction
    }
  }
  }
}

# Add dimension names for clarity
dimnames(array_data_12_scenario2) <- dimnames(array_data_17_scenario2) <- dimnames(array_data_12_int_scenario2) <- dimnames(array_data_17_int_scenario2) <- list(
  Mix = paste0("Mix ", 1:num_mixes),
  Method = c("Normal Mix Normal", "Normal Mix Diagonal", "Diagonal Mix Normal", "Diagonal Mix Diagonal"),
  Iteration = paste0("Iteration ", 1:num_iterations)
)

# Calculate means and standard errors for each method and mix
means12Scenario2 <- apply(array_data_12_scenario2, c(1, 2), mean, na.rm = TRUE)
means17Scenario2 <- apply(array_data_17_scenario2, c(1, 2), mean, na.rm = TRUE)
means12Scenario2_int <- apply(array_data_12_int_scenario2, c(1, 2), mean, na.rm = TRUE)
means17Scenario2_int <- apply(array_data_17_int_scenario2, c(1, 2), mean, na.rm = TRUE)

se12Scenario2 <- apply(array_data_12_scenario2, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
se17Scenario2 <- apply(array_data_17_scenario2, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
se12Scenario2_int <- apply(array_data_12_int_scenario2, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
se17Scenario2_int <- apply(array_data_17_int_scenario2, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

# Convert means and standard errors to data frames for plotting
means12Scenario2 <- as.data.frame(means12Scenario2)
means17Scenario2 <- as.data.frame(means17Scenario2)
means12Scenario2_int <- as.data.frame(means12Scenario2_int)
means17Scenario2_int <- as.data.frame(means17Scenario2_int)

se12Scenario2 <- as.data.frame(se12Scenario2)
se17Scenario2 <- as.data.frame(se17Scenario2)
se12Scenario2_int <- as.data.frame(se12Scenario2_int)
se17Scenario2_int <- as.data.frame(se17Scenario2_int)

# Add mix names to the data frames
means12Scenario2$Mix<-means17Scenario2$Mix<-means12Scenario2_int$Mix<-means17Scenario2_int$Mix<-rownames(means12Scenario2) <- rownames(means17Scenario2) <- rownames(means12Scenario2_int) <- rownames(means17Scenario2_int) <- paste0("Mix ", 1:num_mixes)
se12Scenario2$Mix<-se17Scenario2$Mix<-se12Scenario2_int$Mix<-se17Scenario2_int$Mix<-rownames(se12Scenario2) <- rownames(se17Scenario2) <- rownames(se12Scenario2_int) <- rownames(se17Scenario2_int) <- paste0("Mix ", 1:num_mixes)

# Reshape the data frames for plotting
means_long_12_scenario2 <- means12Scenario2 %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario2),"Mix"), names_to = "Method", values_to = "Mean")
means_long_17_scenario2 <- means17Scenario2 %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario2),"Mix"), names_to = "Method", values_to = "Mean")
means_long_12_int_scenario2 <- means12Scenario2_int %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario2),"Mix"), names_to = "Method", values_to = "Mean")
means_long_17_int_scenario2 <- means17Scenario2_int %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario2),"Mix"), names_to = "Method", values_to = "Mean")

sd_long_12_scenario2 <- se12Scenario2 %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario2),"Mix"), names_to = "Method", values_to = "SE")
sd_long_17_scenario2 <- se17Scenario2 %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario2),"Mix"), names_to = "Method", values_to = "SE")
sd_long_12_int_scenario2 <- se12Scenario2_int %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario2),"Mix"), names_to = "Method", values_to = "SE")
sd_long_17_int_scenario2 <- se17Scenario2_int %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario2),"Mix"), names_to = "Method", values_to = "SE")

# Combine means and standard errors for plotting
means_long_12_scenario2$SE <- sd_long_12_scenario2$SE
means_long_17_scenario2$SE <- sd_long_17_scenario2$SE
means_long_12_int_scenario2$SE <- sd_long_12_int_scenario2$SE
means_long_17_int_scenario2$SE <- sd_long_17_int_scenario2$SE

# Plot the means with error bars for Scenario 2
ggplot(means_long_12_scenario2, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 2 (12 cm)", x = "Method", y = "Mean Correlation") +
  theme_minimal() +
  facet_wrap(~ Mix)

ggplot(means_long_17_scenario2, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 2 (17 cm)", x = "Method", y = "Mean Correlation") +
  theme_minimal() +
  facet_wrap(~ Mix)

p2_12<-ggplot(means_long_12_scenario2, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()

p2_17<-ggplot(means_long_17_scenario2, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()


p2_12_int<-ggplot(means_long_12_int_scenario2, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()


p2_17_int<-ggplot(means_long_17_int_scenario2, aes(x = Mix, y = Mean, fill = Method))   +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()

###########################

correlationsScenario3 <- lapply(allResults, function(iterationResults) {
  scenario3Results <- iterationResults[grep("Scenario3", names(iterationResults))][1:4]
  return(scenario3Results)
})


# Create an array to store the correlation results for Scenario 2
num_mixes <- 4
num_methods_wheat <- 2
num_methods_mix <- 2
num_methods <- 4
num_iterations <- 10

array_data_12_scenario3 <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_17_scenario3 <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_12_int_scenario3 <- array(dim = c(num_mixes, num_methods, num_iterations))
array_data_17_int_scenario3 <- array(dim = c(num_mixes, num_methods, num_iterations))

# Fill the arrays with correlation results
for (iteration in 1:num_iterations) {
  method=0
  for (method_wheat in 1:num_methods_wheat) {
    for (method_mix in 1:num_methods_mix) {
      method=method+1
      for (mix in 1:num_mixes) {
        cor_results <- correlationsScenario3[[paste("Iteration", iteration, "Wheat", names(kinshipWheat)[method_wheat], "Mix", names(kinshipMix)[method_mix])]][[paste("Scenario3 mix", mix, sep = "")]]$CorrelationResults
        array_data_12_scenario3[mix, method, iteration] <- cor_results$Correlation12
        array_data_17_scenario3[mix, method,iteration] <- cor_results$Correlation17
        array_data_12_int_scenario3[mix, method, iteration] <- cor_results$Correlation12Interaction
        array_data_17_int_scenario3[mix, method, iteration] <- cor_results$Correlation17Interaction
      }
    }
  }
}

# Add dimension names for clarity
dimnames(array_data_12_scenario3) <- dimnames(array_data_17_scenario3) <- dimnames(array_data_12_int_scenario3) <- dimnames(array_data_17_int_scenario3) <- list(
  Mix = paste0("Mix ", 1:num_mixes),
  Method = c("Normal Mix Normal", "Normal Mix Diagonal", "Diagonal Mix Normal", "Diagonal Mix Diagonal"),
  Iteration = paste0("Iteration ", 1:num_iterations)
)

# Calculate means and standard errors for each method and mix
means12Scenario3 <- apply(array_data_12_scenario3, c(1, 2), mean, na.rm = TRUE)
means17Scenario3 <- apply(array_data_17_scenario3, c(1, 2), mean, na.rm = TRUE)
means12Scenario3_int <- apply(array_data_12_int_scenario3, c(1, 2), mean, na.rm = TRUE)
means17Scenario3_int <- apply(array_data_17_int_scenario3, c(1, 2), mean, na.rm = TRUE)

se12Scenario3 <- apply(array_data_12_scenario3, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
se17Scenario3 <- apply(array_data_17_scenario3, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
se12Scenario3_int <- apply(array_data_12_int_scenario3, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
se17Scenario3_int <- apply(array_data_17_int_scenario3, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

# Convert means and standard errors to data frames for plotting
means12Scenario3 <- as.data.frame(means12Scenario3)
means17Scenario3 <- as.data.frame(means17Scenario3)
means12Scenario3_int <- as.data.frame(means12Scenario3_int)
means17Scenario3_int <- as.data.frame(means17Scenario3_int)

se12Scenario3 <- as.data.frame(se12Scenario3)
se17Scenario3 <- as.data.frame(se17Scenario3)
se12Scenario3_int <- as.data.frame(se12Scenario3_int)
se17Scenario3_int <- as.data.frame(se17Scenario3_int)

# Add mix names to the data frames
means12Scenario3$Mix<-means17Scenario3$Mix<-means12Scenario3_int$Mix<-means17Scenario3_int$Mix<-rownames(means12Scenario3) <- rownames(means17Scenario3) <- rownames(means12Scenario3_int) <- rownames(means17Scenario3_int) <- paste0("Mix ", 1:num_mixes)
se12Scenario3$Mix<-se17Scenario3$Mix<-se12Scenario3_int$Mix<-se17Scenario3_int$Mix<-rownames(se12Scenario3) <- rownames(se17Scenario3) <- rownames(se12Scenario3_int) <- rownames(se17Scenario3_int) <- paste0("Mix ", 1:num_mixes)

# Reshape the data frames for plotting
means_long_12_scenario3 <- means12Scenario3 %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario3),"Mix"), names_to = "Method", values_to = "Mean")
means_long_17_scenario3 <- means17Scenario3 %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario3),"Mix"), names_to = "Method", values_to = "Mean")
means_long_12_int_scenario3 <- means12Scenario3_int %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario3),"Mix"), names_to = "Method", values_to = "Mean")
means_long_17_int_scenario3 <- means17Scenario3_int %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario3),"Mix"), names_to = "Method", values_to = "Mean")

sd_long_12_scenario3 <- se12Scenario3 %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario3),"Mix"), names_to = "Method", values_to = "SE")
sd_long_17_scenario3 <- se17Scenario3 %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario3),"Mix"), names_to = "Method", values_to = "SE")
sd_long_12_int_scenario3 <- se12Scenario3_int %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario3),"Mix"), names_to = "Method", values_to = "SE")
sd_long_17_int_scenario3 <- se17Scenario3_int %>%
  pivot_longer(cols = setdiff(colnames(means12Scenario3),"Mix"), names_to = "Method", values_to = "SE")

# Combine means and standard errors for plotting
means_long_12_scenario3$SE <- sd_long_12_scenario3$SE
means_long_17_scenario3$SE <- sd_long_17_scenario3$SE
means_long_12_int_scenario3$SE <- sd_long_12_int_scenario3$SE
means_long_17_int_scenario3$SE <- sd_long_17_int_scenario3$SE

# Plot the means with error bars for Scenario 2
ggplot(means_long_12_scenario3, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 2 (12 cm)", x = "Method", y = "Mean Correlation") +
  theme_minimal() +
  facet_wrap(~ Mix)

ggplot(means_long_17_scenario3, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 2 (17 cm)", x = "Method", y = "Mean Correlation") +
  theme_minimal() +
  facet_wrap(~ Mix)



p3_12<-ggplot(means_long_12_scenario3, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()

p3_17<-ggplot(means_long_17_scenario3, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()

p3_12_int<-ggplot(means_long_12_int_scenario3, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()

p3_17_int<-ggplot(means_long_17_int_scenario3, aes(x = Mix, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean Correlation for Scenario 1", x = "Mix", y = "Mean Correlation") +
  theme_minimal()




#############Another scenario where we leave one mix out we make predictions for the left out mix for the wheat, store the predictions.
# After going over all mixes, we calculate the correlation between the predictions and the actual phenotypes.

runGenomicPredictionStrategy4 <-function(traitName, kinshipMatrixWheat, kinshipMatrixMix, testMix) {
    phenotypeTrain12 <- phenotypeData12cm
    phenotypeTrain17 <- phenotypeData17cm
    phenotypeTrain12[phenotypeTrain12$Strain %in% testMix, traitName] <- NA
    phenotypeTrain17[phenotypeTrain17$Strain %in% testMix, traitName] <- NA
    print(sum(is.na(phenotypeTrain12[phenotypeTrain12$Strain %in% testMix, traitName])))
    ZtrainWheat12 <- model.matrix(~0 + Plant, data = phenotypeTrain12)
    ZtrainWheat17 <- model.matrix(~0 + Plant, data = phenotypeTrain17)
    XtrainWheat12 <- model.matrix(~1, data = phenotypeTrain12)
    XtrainWheat17 <- model.matrix(~1, data = phenotypeTrain17)
    ZTrainMix12 <- model.matrix(~0 + Strain, data = phenotypeTrain12)
    ZTrainMix17 <- model.matrix(~0 + Strain, data = phenotypeTrain17)

    fm12 <- BGLR(y = as.numeric(phenotypeTrain12[[traitName]]), ETA = list(
        list(X = XtrainWheat12, model = "FIXED"),
        list(K = ZtrainWheat12 %*% kinshipMatrixWheat %*% t(ZtrainWheat12), model = "RKHS"),
        list(K = ZTrainMix12 %*% kinshipMatrixMix %*% t(ZTrainMix12), model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    fm17 <- BGLR(y = as.numeric(phenotypeTrain17[[traitName]]), ETA = list(
        list(X = XtrainWheat17, model = "FIXED"),
        list(K = ZtrainWheat17 %*% kinshipMatrixWheat %*% t(ZtrainWheat17), model = "RKHS"),
        list(K = ZTrainMix17 %*% kinshipMatrixMix %*% t(ZTrainMix17), model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    K12_mix <- ZTrainMix12 %*% kinshipMatrixMix %*% t(ZTrainMix12)
    K12_wheat <- ZtrainWheat12 %*% kinshipMatrixWheat %*% t(ZtrainWheat12)
    K12_combined <- K12_mix * K12_wheat

    K17_mix <- ZTrainMix17 %*% kinshipMatrixMix %*% t(ZTrainMix17)
    K17_wheat <- ZtrainWheat17 %*% kinshipMatrixWheat %*% t(ZtrainWheat17)
    K17_combined <- K17_mix * K17_wheat
    fm12_interaction <- BGLR(y = as.numeric(phenotypeTrain12[[traitName]]), ETA = list(
        list(X = XtrainWheat12, model = "FIXED"),
        list(K = K12_mix, model = "RKHS"),
        list(K = K12_wheat, model = "RKHS"),
        list(K = K12_combined, model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)

    fm17_interaction <- BGLR(y = as.numeric(phenotypeTrain17[[traitName]]), ETA = list(
        list(X = XtrainWheat17, model = "FIXED"),
        list(K = K17_mix, model = "RKHS"),
        list(K = K17_wheat, model = "RKHS"),
        list(K = K17_combined, model = "RKHS")
    ), nIter = 10000, burnIn = 2000, verbose = FALSE)
    return(list(fm12 = fm12, fm17 = fm17, fm12_interaction = fm12_interaction, fm17_interaction = fm17_interaction, TestMix = testMix))
}


GetModelPredictions4 <- function(models, phenotypeData12, phenotypeData17, traitName) {
  predictions12 <- predict(models$fm12)
  predictions17 <- predict(models$fm17)
  predictions12_interaction <- predict(models$fm12_interaction)
  predictions17_interaction <- predict(models$fm17_interaction)
  testLocations12 <- phenotypeData12$Strain %in% models$TestMix
  testLocations17 <- phenotypeData17$Strain %in% models$TestMix
  
  predResults <- list(
    pred12 = predictions12[testLocations12], 
    pred17 = predictions17[testLocations17], 
    pred12_int = predictions12_interaction[testLocations12], 
    pred17_int = predictions12_interaction[testLocations12],
    mix=models$TestMix,
    testLocations12=testLocations12,
    testLocations17=testLocations17
  )
  return(predResults)
}


ST2ModelsMix1<-runGenomicPredictionStrategy4("PLACL", KmatWheatDiag, KmatSepMix, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix1")
predsMix1<-GetModelPredictions4(ST2ModelsMix1, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
ST2ModelsMix2<-runGenomicPredictionStrategy4("PLACL", KmatWheatDiag, KmatSepMix, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix2")
predsMix2<-GetModelPredictions4(ST2ModelsMix2, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
ST2ModelsMix3<-runGenomicPredictionStrategy4("PLACL", KmatWheatDiag, KmatSepMix, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix3")
predsMix3<-GetModelPredictions4(ST2ModelsMix3, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")
ST2ModelsMix4<-runGenomicPredictionStrategy4("PLACL", KmatWheatDiag, KmatSepMix, Pheno12cmImpAgg, Pheno17cmImpAgg, "mix4")
predsMix4<-GetModelPredictions4(ST2ModelsMix4, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")


cor(c(predsMix1$pred12,predsMix2$pred12,predsMix3$pred12,predsMix4$pred12),c(Pheno12cmImpAgg[predsMix1$testLocations12,"PLACL"],
                                                                         Pheno12cmImpAgg[predsMix2$testLocations12,"PLACL"],
                                                                         Pheno12cmImpAgg[predsMix3$testLocations12,"PLACL"],
                                                                         Pheno12cmImpAgg[predsMix4$testLocations12,"PLACL"]))

colvec<-rep(0,length(c(Pheno12cmImpAgg[predsMix1$testLocations12,"PLACL"],
                       Pheno12cmImpAgg[predsMix2$testLocations12,"PLACL"],
                       Pheno12cmImpAgg[predsMix3$testLocations12,"PLACL"],
                       Pheno12cmImpAgg[predsMix4$testLocations12,"PLACL"])))


colvec[predsMix1$testLocations12]<-1
sum(predsMix1$testLocations12)
colvec[predsMix2$testLocations12]<-2
sum(predsMix2$testLocations12)
colvec[predsMix3$testLocations12]<-3
sum(predsMix3$testLocations12)
colvec[predsMix4$testLocations12]<-4
sum(predsMix4$testLocations12)
plot(c(predsMix1$pred12,predsMix2$pred12,predsMix3$pred12,predsMix4$pred12),c(Pheno12cmImpAgg[predsMix1$testLocations12,"PLACL"],
                                                                              Pheno12cmImpAgg[predsMix2$testLocations12,"PLACL"],
                                                                              Pheno12cmImpAgg[predsMix3$testLocations12,"PLACL"],
                                                                              Pheno12cmImpAgg[predsMix4$testLocations12,"PLACL"]), col= colvec)


cor(c(mean(predsMix1$pred12),mean(predsMix2$pred12),mean(predsMix3$pred12),mean(predsMix4$pred12)),c(mean(Pheno12cmImpAgg[predsMix1$testLocations12,"PLACL"]),
                                                                                                      mean(Pheno12cmImpAgg[predsMix2$testLocations12,"PLACL"]),
                                                                                                           mean(Pheno12cmImpAgg[predsMix3$testLocations12,"PLACL"]),
                                                                                                                mean(Pheno12cmImpAgg[predsMix4$testLocations12,"PLACL"])))

alldata<-cbind(predsMix1$pred12,predsMix2$pred12,predsMix3$pred12,predsMix4$pred12,Pheno12cmImpAgg[predsMix1$testLocations12,"PLACL"],
                                                                               Pheno12cmImpAgg[predsMix2$testLocations12,"PLACL"],
                                                                               Pheno12cmImpAgg[predsMix3$testLocations12,"PLACL"],
                                                                               Pheno12cmImpAgg[predsMix4$testLocations12,"PLACL"])

corsoverwheat<-apply(alldata,1,function(x)cor(x[1:4],x[5:8]))

boxplot(corsoverwheat)



alldata<-cbind(predsMix1$pred17,predsMix2$pred17,predsMix3$pred17,predsMix4$pred17,Pheno17cmImpAgg[predsMix1$testLocations17,"PLACL"],
               Pheno17cmImpAgg[predsMix2$testLocations17,"PLACL"],
               Pheno17cmImpAgg[predsMix3$testLocations17,"PLACL"],
               Pheno17cmImpAgg[predsMix4$testLocations17,"PLACL"])

corsoverwheat<-apply(alldata,1,function(x)cor(x[1:4],x[5:8]))





##############################################GWAS for this phenotype. 
source("../gapit_functions.txt")  # GAPIT functions
models=c("Blink")


GenoWheat[1:5,1:5]

tGenoWheat<-t(GenoWheat)
tGenoWheat[1:5,1:5]
tGenoWheat<-tGenoWheat[!duplicated(tGenoWheat),]
dim(tGenoWheat)
GenoWheat<-t(tGenoWheat)
GenoWheat[1:5,1:5]
Pheno12cmImpAgg[1:5,]
KmatGapit<-rrBLUP::A.mat(GenoWheat-1)
KmatGapit[1:5,1:5]
KmatGapit<-data.frame(taxa=rownames(KmatGapit), as.data.frame(KmatGapit))
colnames(KmatGapit)<-NULL
rownames(KmatGapit)<-NULL
GenoWforGAPIT<-cbind(data.frame(taxa=rownames(GenoWheat)),as.data.frame(GenoWheat))
MapWheatGapit<-MapWheat
MapWheatGapit[1:5,]
colnames(MapWheatGapit)<-c("Name","Chromosome", "Position")
MapWheatGapit<-MapWheatGapit[MapWheatGapit$Name%in%colnames(GenoWheat),]
Pheno12cmGAPIT<-Pheno12cmImpAgg[,-2]
colnames(Pheno12cmGAPIT)[1]<-"taxa"
Pheno12cmImpAgg$Strain<-drop.levels(Pheno12cmImpAgg$Strain)
Pheno12cmGAPITFixed<-as.data.frame(model.matrix(~Strain, data=Pheno12cmImpAgg))[,-1]
Pheno12cmGAPITFixed<-cbind(data.frame(taxa=Pheno12cmGAPIT$taxa), Pheno12cmGAPITFixed)
Pheno12cmGAPITFixed[1:5,]
# Assuming you start in the root directory where you want to create subdirectories
initial_dir <- getwd()

# Array to hold the number of principal components to be used in each analysis
pc_values <- c(0,2, 5)

# Model types to be used
models <- c("MLMM", "Blink")

# Loop through each principal component setting
for (pc in pc_values) {
  # Loop through model types
  for (model in models) {
    # Loop through with and without Kinship matrix (K)
    for (use_K in c(TRUE, FALSE)) {
      # Create a directory name based on the principal component and Kinship usage
      K_label <- ifelse(use_K, "withK", "withoutK")
      dir_name <- paste("GWAS_results_WheatJune25_", pc, "PC_", model, "_", K_label, sep = "")
      
      # Check if the directory exists, if not create it
      if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE)
      }
      
      # Set the working directory to the new directory
      setwd(file.path(initial_dir, dir_name))
      
      # Determine the Kinship matrix input based on the use_K flag
      if (use_K){Kmat_input <-KmatGapit } else {Kmat_input<-NULL}
      
      # Run the GAPIT function
      tmp <- capture.output({
        GAPIT(
          Y = Pheno12cmGAPIT,
          GD = GenoWforGAPIT,
          GM = MapWheatGapit,
          KI = Kmat_input,
          CV = Pheno12cmGAPITFixed,
          PCA.total = pc,
          model = model,
          file.output = TRUE,
          cutOff = 0.1
        )
      })
      
      # Go back to the initial directory
      setwd(initial_dir)
    }
  }
}

# Optionally print a message when the loop is complete
cat("GAPIT analysis complete for all principal component settings and models.\n")

boxplot(corsoverwheat)