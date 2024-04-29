setwd("/Users/denizakdemir/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/Rscripts/GenomicPrediction")

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


KmatSepMixDiag<-diag(nrow(KmatSepMix))
rownames(KmatSepMixDiag)<-colnames(KmatSepMixDiag)<-rownames(KmatSepMix)
KmatWheatDiag<-diag(nrow(KmatWheat))
rownames(KmatWheatDiag)<-colnames(KmatWheatDiag)<-rownames(KmatWheat)


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
    predictions12 <- predict(models$fm12)
    predictions17 <- predict(models$fm17)
    predictions12_interaction <- predict(models$fm12_interaction)
    predictions17_interaction <- predict(models$fm17_interaction)
    testLocations12 <- phenotypeData12$Plant %in% models$TestWheat
    testLocations17 <- phenotypeData17$Plant %in% models$TestWheat

    phenotypeTest12 <- phenotypeData12[testLocations12,]
    phenotypeTest17 <- phenotypeData17[testLocations17,]
    phenotypeTest12Trait<-phenotypeTest12[[traitName]]
    phenotypeTest17Trait<-phenotypeTest17[[traitName]]
    correlation12 <- cor(predictions12[testLocations12], phenotypeTest12Trait)
    correlation17 <- cor(predictions17[testLocations17], phenotypeTest17Trait)
    correlation12_interaction <- cor(predictions12_interaction[testLocations12], phenotypeTest12Trait)
    correlation17_interaction <- cor(predictions17_interaction[testLocations17], phenotypeTest17Trait)
    
    # Plotting results
    plot(predictions12_interaction[testLocations12], phenotypeTest12Trait, col = phenotypeTest12$Strain)
    plot(predictions17_interaction[testLocations17], phenotypeTest17Trait, col = phenotypeTest17$Strain)

    correlationResults <- list(
        Correlation12 = correlation12, 
        Correlation17 = correlation17, 
        Correlation12Interaction = correlation12_interaction, 
        Correlation17Interaction = correlation17_interaction
    )

    # Prepare data for ggplot
    dataFrameForPlot12 <- data.frame(
        Predictions = predictions12_interaction[testLocations12], 
        Trait = phenotypeTest12Trait, 
        Strain = phenotypeTest12$Strain
    )
    ggplot12 <- ggplot(dataFrameForPlot12, aes(x = Predictions, y = Trait, color = Strain)) +
        geom_point() + geom_smooth(method = "lm")

    dataFrameForPlot17 <- data.frame(
        Predictions = predictions17_interaction[testLocations17], 
        Trait = phenotypeTest17Trait, 
        Strain =phenotypeTest17$Strain
    )
    ggplot17 <- ggplot(dataFrameForPlot17, aes(x = Predictions, y = Trait, color = Strain)) +
        geom_point() + geom_smooth(method = "lm")

    plotsResults <- list(ggplot12 = ggplot12, ggplot17 = ggplot17)
    # Calculate correlations within strains
# Combine predictions and actual trait values into a new data frame
data12 <- data.frame(predictions = predictions12_interaction[testLocations12],
                     traits = phenotypeTest12Trait,
                     Strain = phenotypeTest12$Strain[testLocations12])

data17 <- data.frame(predictions = predictions17_interaction[testLocations17],
                     traits = phenotypeTest17Trait,
                     Strain = phenotypeTest17$Strain[testLocations17])

# Split the data by 'Strain'
listData12 <- split(data12, data12$Strain)
listData17 <- split(data17, data17$Strain)

# Calculate correlation for each strain
Cor12withinStrain <- lapply(listData12, function(df) cor(df$predictions, df$traits, use = "complete.obs"))
Cor17withinStrain <- lapply(listData17, function(df) cor(df$predictions, df$traits, use = "complete.obs"))
    correlationResults$Cor12withinStrain<-Cor12withinStrain
    correlationResults$Cor17withinStrain<-Cor17withinStrain
    return(list(CorrelationResults = correlationResults, PlotsResults = plotsResults))
}


ST1Models<-runGenomicPredictionStrategy("PLACL", KmatWheatDiag, KmatSepMixDiag, Pheno12cmImpAgg, Pheno17cmImpAgg)
outcor<-evaluateModelPredictions(ST1Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")

outcor$CorrelationResults
outcor$PlotsResults$ggplot12
outcor$PlotsResults$ggplot17


ST1Models<-runGenomicPredictionStrategy("PLACL", KmatWheat, KmatSepMix, Pheno12cmImpAgg, Pheno17cmImpAgg)

outcor<-evaluateModelPredictions(ST1Models, Pheno12cmImpAgg, Pheno17cmImpAgg,"PLACL")

outcor$CorrelationResults
outcor$PlotsResults$ggplot12
outcor$PlotsResults$ggplot17

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

correlationsScenario1 <- lapply(allResults, function(iterationResults) {
  scenario1Results <- iterationResults$Scenario1
  cor12withinStrain<-scenario1Results$CorrelationResults$Cor12withinStrain
  cor17withinStrain<-scenario1Results$CorrelationResults$Cor17withinStrain
  return(list(Cor12withinStrain=cor12withinStrain, Cor17withinStrain=cor17withinStrain))
})

# convert into data frame
correlationsScenario1DF <- do.call(rbind, correlationsScenario1)
correlationsScenario1DF <- as.data.frame(correlationsScenario1DF)
correlationsScenario1DF$Iteration <- 1:nrow(correlationsScenario1DF)
correlationsScenario1DF$experiment<-rownames(correlationsScenario1DF)
# table of mean correlations and standard errors
colnames(correlationsScenario1DF)
meanCorrelations <- aggregate(. ~ Iteration+experiment, data = correlationsScenario1DF, FUN = function(x) c(Mean = mean(x), SE = sd(x) / sqrt(length(x))))
meanCorrelations
