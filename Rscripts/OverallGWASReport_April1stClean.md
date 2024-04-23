---
title: "Genome-Wide Association Study (GWAS) Methodology"
author: "Deniz Akdemir"
date: "April 1rst, 2024"
output:
  html_document:
    keep_md: true
    toc: true
  pdf_document: default
  word_document: default  
---





## Abstract

This document outlines the methodology employed in our Genome-Wide Association Study (GWAS) to dissect the genetic architecture underlying both agronomic traits and a specialized set of "endophyte traits" in wheat. By utilizing a robust two-step process, our GWAS framework is designed to provide comprehensive insights essential for breeding programs focused on crop improvement.

## Introduction

The genetic exploration of phenotypic traits in crops like wheat is critical for the advancement of targeted breeding strategies. Our study leverages a GWAS framework to systematically investigate the genetic basis of diverse traits, including agronomic and unique "endophyte traits," aiming to enhance our understanding of wheat's genetic potential.

## Methodology

## Phenotypic Data Description for Endophyte Traits

These experiments were designed to explore the interaction between wheat strains and two specific treatments of endophytes, identified as R3E3 and R3E5. The objective was to investigate the underlying genetic factors influencing wheat's response to these endophytes.

### Experimental Design

Our phenotypic dataset encapsulates various measurements captured from the interaction of Zymoseptoria tritici (Zt) strains with the endophytes across biological and technical replicates. The dataset columns are defined as follows:

-   **Endophyte:** Represents the two treatments under investigation, R3E3 and R3E5.
-   **Breplicate:** Indicates the biological replicate number.
-   **Treplicate:** Denotes the technical replicate number.
-   **Zt_Strain:** Specifies the Zymoseptoria tritici strain used in each experiment.
-   **Area_endophyte:** Measures the area occupied by the endophyte within the Petri dish.
-   **Area_inhibition:** Quantifies the area showing no growth of Zt, indicative of inhibition by the endophyte.
-   **Ratio_inhibition:** Represents the normalized area of inhibition relative to the endophyte-occupied area.
-   **Ratio_Scales1 and Ratio_Scales2:** Categorical scales for the inhibition ratio, with different granularity (every cm and every 0.5 cm, respectively).
-   **Radius_inhibition:** Calculated as the halo radius minus the endophyte radius, based on perimeter measurements. This metric is considered reliable only for the R3E3 treatment due to its approximate circular shape, whereas R3E5 exhibits more irregular patterns.
-   **Percentage_inhibition:** The percentage of the Petri dish area covered by the inhibition zone.
-   **Melanization:** A binary indicator (1 for melanized, 0 for not melanized) reflecting whether the Zt strain exhibits melanization when in contact with the endophyte, compared to a control without endophyte.

Additional notes include the identification of data points considered for exclusion due to less dense lawns (`Remove_astringent`) or difficulty in distinguishing the inhibition area (`Remove_loose`).

## Phenotypic Data Description for Greenhouse Traits

### Data Overview

The greenhouse phenotypic dataset incorporates several key attributes:

-   **Set:** Groups the data into distinct experimental sets for systematic analysis.
-   **Plant:** A unique identifier for each plant subject to measurement.
-   **Strain:** Indicates the specific strain of wheat under observation.
-   **Leaf:** Assigns an identifier to individual leaves, allowing for detailed trait analysis at the foliar level.
-   **Rep:** Represents the replication number, underlining the repeatability of measurements for robustness.
-   **leafArea:** Quantifies the total area of the leaf in square centimeters, providing insights into growth patterns.
-   **necrosisArea:** Measures the extent of necrosis, offering a direct indicator of disease impact or stress response.
-   **PLACL:** Calculates the percentage of the leaf area affected by lesions, crucial for assessing disease severity.
-   **pycnidiaCount:** Enumerates the total pycnidia present, reflecting fungal disease prevalence on the leaf surface.
-   **meanPycnidiaArea:** Averages the area occupied by each pycnidium, indicating the extent of fungal colonization.
-   **pycnidiaPerCm2Leaf:** Determines the density of pycnidia relative to the leaf area, useful for quantifying infection intensity.
-   **pycnidiaPerCm2Lesion:** Assesses the concentration of pycnidia within lesioned areas, offering insights into disease progression.
-   **pycnidiaGreyValue:** Assigns a grey value to the pycnidia, potentially correlating with spore viability or maturity.
-   **Picture:** Provides photographic evidence of the leaf, essential for visual confirmation of phenotypic traits.
-   **cm:** Denotes the height at which measurements were taken, ensuring consistency across data collection.

Total Observations: The dataset comprises 9,619 individual leaf measurements, providing a robust foundation for our GWAS analysis

Wheat Varieties: The dataset encapsulates the responses of 4 unique wheat varieties to pathogen exposure. This diversity is pivotal for understanding the genetic breadth of disease resistance and susceptibility across different wheat genotypes.

Data Coverage: Phenotypic measurements span a total of 2 unique leaf identifiers (leave_id) across 2 experimental replications (REP). This comprehensive data collection ensures robustness and reliability in our phenotypic assessments.

## Genomic Data Preprocessing and SNP Selection

### SNP Preprocessing

Following initial data acquisition, our genomic dataset underwent preprocessing using TASSEL. This step was essential for standardizing data formats and facilitating subsequent analyses. Through this preprocessing phase, we aimed to refine our dataset to include only the most informative genetic markers. To enhance the robustness of our GWAS, we applied stringent filtering criteria to our SNP dataset. Specifically, we filtered SNPs based on missingness and Minor Allele Frequency (MAF), setting thresholds at 0.3 and 0.1, respectively. Analysis was confined to genetic markers located on chromosomes 1 to 13. These criteria were designed to exclude SNPs with high levels of missing data and those with low allele frequencies, which might otherwise obscure true genetic associations due to insufficient statistical power or poor data quality.

After applying our filtering criteria, our final SNP dataset comprised approximately 400,000 SNPs. This refined set of genetic markers represents a comprehensive basis for our GWAS, capturing a wide array of genetic variations within the Z. tritici strains under study.

### SNP Thinning for Enhanced Analysis

In addition to analyzing the entire set of 400,000 SNPs, we explored a reduced SNP dataset to further refine our genetic association mapping. This reduction involved selecting the SNP that exhibited the highest correlation with its neighbors within 100kb windows across the genome. This thinning approach aimed to distill our SNP dataset to those markers most representative of genetic variation in each genomic region, potentially enhancing the resolution and interpretability of our GWAS findings and increasing the statistical power to detect genetic associations.

### Population Structure Analysis in Zymoseptoria tritici Strains

Understanding the population structure is crucial for our GWAS to minimize false associations due to stratification. In this analysis we use the entire set of 400,000 SNPs to infer the population structure of Z. tritici strains. We employed principal component analysis (PCA) to visualize the genetic relationships among strains and identify potential subpopulations within our dataset. By incorporating PCA results as covariates in our GWAS models, we aim to account for population structure and reduce spurious associations.

![](FiguresForReport/Combined_GenoMat.png)

### Supplementary GWAS Analysis Informed by Structure Analysis

To further refine our insights, an additional GWAS analysis was conducted, informed by a principal component (PC) analysis. This PC plot highlighted a distinct cluster of 10 genotypes, prompting their exclusion to assess the impact on the genetic associations identified. This targeted analysis enhances the accuracy of our conclusions, ensuring they reflect the comprehensive genetic diversity of wheat.

## GWAS methods

### Step 1: Estimation of Genomic Values Using Mixed-Effects Models

#### Agronomic Traits

For agronomic traits, the estimation of genomic values involves adjusting for known fixed effects and accounting for random effects associated with isolates. The linear mixed-effects models (LMMs) are specified as:

$$
y_{ijklm} = \mu + {\alpha_{rep}}_i + {\beta_{geno}}_j + {\gamma_{year}}_k + (1 | {isolate}_l) + \epsilon_{ijklm}
$$

Where:

-   $y_{ijklm}$ represents the observed phenotype for the $i^{th}$ replicate, $j^{th}$ genotype, $k^{th}$ year, and $l^{th}$ isolate.
-   $\mu$ is the overall mean phenotype.
-   ${\alpha_{rep}}_i$, ${\beta_{geno}}_j$, and ${\gamma_{year}}_k$ are fixed effects for replicate, genotype, and year, respectively.
-   $1 | {isolate}_l$ captures the random effects associated with the $l^{th}$ isolate.
-   $\epsilon_{ijklm}$ is the residual error term, incorporating an additional index $m$ for more precise error modeling.

#### Endophyte Traits

For "endophyte traits", the model adjusts to reflect the experimental design, specified as:

$$
y_{ijkl} = \mu + {\alpha_{BReplicate}}_i + {\beta_{TReplicate}}_j + (1 | {isolate}_l) + \epsilon_{ijlm}
$$

Where:

-   $y_{ijlm}$ denotes the observed value for each trait, such as Radius_inhibition, Area_Endophyte, etc.
-   ${\alpha_{BReplicate}}_i$ and ${\beta_{TReplicate}}_j$ are the fixed effects for biological and technical replicates, respectively.
-   $1 | {isolate}_l$ represents the random effects associated with the $l^{th}$ isolate.
-   $\epsilon_{ijlm}$ is the residual error term for each observation.

BLUPs for the random effects are then extracted for each phenotypic trait within both datasets.

### Step 2: Genetic Association Mapping Using GAPIT

This phase is identical for both sets of traits, employing GAPIT in R for association mapping to identify SNP markers associated with the adjusted phenotypic values. The analysis incorporates models like FarmCPU and BLINK, designed to manage the complexities of genomic data effectively. FarmCPU (Fixed and Random Model Circulating Probability Unification) significantly enhances our analysis by addressing computational challenges, optimizing marker selection, and refining adjustments for population structure and kinship. By replacing the REM with a more efficient FEM and utilizing Bayesian information criteria for model selection, BLINK addresses the limitations of previous GWAS methodologies.

## GWAS results

To present the GWAS results methodically and cohesively, the findings will be organized by phenotype experiment, followed by detailed results for each trait.

### Agronomic Traits Analysis

Manhattan and QQ Plots: For each agronomic trait, a Manhattan plot to visualize significant SNPs and a QQ plot to assess the overall significance level and potential inflation due to population structure. Significant Markers Table: A table listing all significant SNPs identified for each agronomic trait, including their chromosomal location, effect size, and p-value.

Results for whaole data and aLso results are shown for each trait using the whole data or using the second and third leaves only.

#### Analysis 1: Complete SNP Set (634156 SNPs), use both leaves 2 and 3







|   |       X|SNP         | Chr|     Pos| P.value|       MAF|traits                       |
|:--|-------:|:-----------|---:|-------:|-------:|---------:|:----------------------------|
|4  | 1192521|S2_835783   |   2|  835783|   0e+00| 0.2188209|BLINK.BLUPpycnidiaCount      |
|7  | 5732161|S11_1490008 |  11| 1490008|   0e+00| 0.1433318|BLINK.BLUPmeanPycnidiaArea   |
|6  |  596212|S12_966299  |  12|  966299|   0e+00| 0.1531515|BLINK.BLUPpycnidiaCount      |
|5  |  448409|S8_1998041  |   8| 1998041|   0e+00| 0.1546240|BLINK.BLUPpycnidiaCount      |
|2  |  573216|S11_1490008 |  11| 1490008|   0e+00| 0.1433318|FarmCPU.BLUPmeanPycnidiaArea |
|3  |  573217|S11_1490011 |  11| 1490011|   1e-07| 0.1866747|FarmCPU.BLUPmeanPycnidiaArea |
|1  |  119252|S2_835783   |   2|  835783|   1e-07| 0.2188209|FarmCPU.BLUPpycnidiaCount    |



```
##  Marker density plotting.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-3-1.png)<!-- -->







```
##  (warning: all phenotypes will use the same ylim.)
##  Rectangular Manhattan plotting pycnidiaCount.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
##  Rectangular Manhattan plotting meanPycnidiaArea.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-4-2.png)<!-- -->



```
##  Q-Q plotting pycnidiaCount.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```
##  Q-Q plotting meanPycnidiaArea.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-5-2.png)<!-- -->


#### Analysis 2: Thinned SNP Set (167938 SNPs), use both leaves 2 and 3



|   |      X|SNP        | Chr|     Pos| P.value|       MAF|traits                           |
|:--|------:|:----------|---:|-------:|-------:|---------:|:--------------------------------|
|6  | 102602|S2_830502  |   2|  830502|   0e+00| 0.4032468|FarmCPU.BLUPpycnidiaCount        |
|10 | 102605|S2_830502  |   2|  830502|   0e+00| 0.4032468|BLINK.BLUPpycnidiaCount          |
|9  | 102604|S2_830502  |   2|  830502|   0e+00| 0.4032468|BLINK.BLUPpycnidiaPerCm2Lesion   |
|8  | 102603|S2_830502  |   2|  830502|   0e+00| 0.4032468|BLINK.BLUPpycnidiaPerCm2Leaf     |
|7  |  18627|S3_2969080 |   3| 2969080|   0e+00| 0.2310873|FarmCPU.BLUPpycnidiaCount        |
|1  |  10260|S2_830502  |   2|  830502|   5e-07| 0.4032468|FarmCPU.BLUPpycnidiaPerCm2Leaf   |
|3  | 102601|S2_830502  |   2|  830502|   5e-07| 0.4032468|FarmCPU.BLUPpycnidiaPerCm2Lesion |
|4  |  10270|S2_837377  |   2|  837377|   6e-07| 0.2685285|FarmCPU.BLUPpycnidiaPerCm2Lesion |
|5  | 102711|S2_837741  |   2|  837741|   6e-07| 0.2657498|FarmCPU.BLUPpycnidiaPerCm2Lesion |
|2  |  10271|S2_837741  |   2|  837741|   1e-06| 0.2657498|FarmCPU.BLUPpycnidiaPerCm2Leaf   |




```
##  Marker density plotting.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-8-1.png)<!-- -->





```
##  (warning: all phenotypes will use the same ylim.)
##  Rectangular Manhattan plotting pycnidiaCount.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```
##  Rectangular Manhattan plotting pycnidiaPerCm2Leaf.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

```
##  Rectangular Manhattan plotting pycnidiaPerCm2Lesion.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-9-3.png)<!-- -->



```
##  Q-Q plotting pycnidiaCount.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```
##  Q-Q plotting pycnidiaPerCm2Leaf.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

```
##  Q-Q plotting pycnidiaPerCm2Lesion.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-10-3.png)<!-- -->



#### Analysis 3: Thinned SNP Set (167938 SNPs), use Leaf 2 only



|   |     X|SNP        | Chr|     Pos| P.value|       MAF|traits                         |
|:--|-----:|:----------|---:|-------:|-------:|---------:|:------------------------------|
|2  | 98711|S2_566639  |   2|  566639|   0e+00| 0.2774916|BLINK.BLUPpycnidiaPerCm2Leaf   |
|3  | 38280|S9_389950  |   9|  389950|   0e+00| 0.1114674|BLINK.BLUPpycnidiaPerCm2Leaf   |
|4  | 39280|S9_1072639 |   9| 1072639|   1e-07| 0.2819307|BLINK.BLUPpycnidiaPerCm2Leaf   |
|1  |  9871|S2_566639  |   2|  566639|   5e-07| 0.2774916|FarmCPU.BLUPpycnidiaPerCm2Leaf |







```
##  Rectangular Manhattan plotting pycnidiaPerCm2Leaf.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-13-1.png)<!-- -->



```
##  Q-Q plotting pycnidiaPerCm2Leaf.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-14-1.png)<!-- -->




#### Analysis 4: Thinned SNP Set (167938 SNPs), use Leaf 3 only



|   |      X|SNP       | Chr|   Pos| P.value|       MAF|traits            |
|:--|------:|:---------|---:|-----:|-------:|---------:|:-----------------|
|2  | 435401|S11_92810 |  11| 92810|   0e+00| 0.1014031|BLINK.BLUPPLACL   |
|1  |  43540|S11_92810 |  11| 92810|   9e-07| 0.1014031|FarmCPU.BLUPPLACL |







```
##  Rectangular Manhattan plotting PLACL.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```
##  Rectangular Manhattan plotting PLACL.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-17-2.png)<!-- -->



```
##  Q-Q plotting PLACL.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-18-1.png)<!-- -->



### Endophyte Traits Analysis


#### Analysis 1: Complete SNP Set (634156 SNPs)







|   |       X|SNP         | Chr|     Pos| P.value|       MAF|traits                                |
|:--|-------:|:-----------|---:|-------:|-------:|---------:|:-------------------------------------|
|8  | 1498483|S2_3164588  |   2| 3164588|   0e+00| 0.1133621|BLINK.BLUPPercentage_inhibitionR3E3   |
|9  |  597188|S12_1035520 |  12| 1035520|   0e+00| 0.4604587|BLINK.BLUPPercentage_inhibitionR3E3   |
|3  | 1498481|S2_3164588  |   2| 3164588|   1e-07| 0.1133621|FarmCPU.BLUPPercentage_inhibitionR3E3 |
|5  | 1498501|S2_3164618  |   2| 3164618|   1e-07| 0.1133621|FarmCPU.BLUPPercentage_inhibitionR3E3 |
|1  |  149848|S2_3164588  |   2| 3164588|   1e-07| 0.1133621|FarmCPU.BLUPArea_inhibitionR3E3       |
|2  |  149850|S2_3164618  |   2| 3164618|   1e-07| 0.1133621|FarmCPU.BLUPArea_inhibitionR3E3       |
|6  | 1498482|S2_3164588  |   2| 3164588|   1e-07| 0.1133621|BLINK.BLUPArea_inhibitionR3E3         |
|7  | 1498502|S2_3164618  |   2| 3164618|   1e-07| 0.1133621|BLINK.BLUPArea_inhibitionR3E3         |
|4  |  149849|S2_3164606  |   2| 3164606|   1e-07| 0.1236207|FarmCPU.BLUPPercentage_inhibitionR3E3 |



```
##  Marker density plotting.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-21-1.png)<!-- -->







```
##  (warning: all phenotypes will use the same ylim.)
##  Rectangular Manhattan plotting Area_inhibition.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```
##  Rectangular Manhattan plotting Percentage_inhibition.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-22-2.png)<!-- -->



```
##  Q-Q plotting Area_inhibition.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```
##  Q-Q plotting Percentage_inhibition.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-23-2.png)<!-- -->



#### Analysis 2: Filtered SNP Set (167938 SNPs)







|   |      X|SNP         | Chr|     Pos| P.value|       MAF|traits                                |
|:--|------:|:-----------|---:|-------:|-------:|---------:|:-------------------------------------|
|4  | 449163|S2_3164606  |   2| 3164606|   0e+00| 0.1236207|BLINK.BLUPPercentage_inhibitionR3E3   |
|5  | 159597|S12_1035503 |  12| 1035503|   0e+00| 0.4676106|BLINK.BLUPPercentage_inhibitionR3E3   |
|2  | 449161|S2_3164606  |   2| 3164606|   1e-07| 0.1236207|FarmCPU.BLUPPercentage_inhibitionR3E3 |
|1  |  44916|S2_3164606  |   2| 3164606|   2e-07| 0.1236207|FarmCPU.BLUPArea_inhibitionR3E3       |
|3  | 449162|S2_3164606  |   2| 3164606|   2e-07| 0.1236207|BLINK.BLUPArea_inhibitionR3E3         |



```
##  Marker density plotting.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-26-1.png)<!-- -->







```
##  (warning: all phenotypes will use the same ylim.)
##  Rectangular Manhattan plotting Area_inhibition.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

```
##  Rectangular Manhattan plotting Percentage_inhibition.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-27-2.png)<!-- -->



```
##  Q-Q plotting Area_inhibition.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

```
##  Q-Q plotting Percentage_inhibition.
```

![](OverallGWASReport_April1stClean_files/figure-html/unnamed-chunk-28-2.png)<!-- -->



![FiguresForReport/BLUPPercentage_inhibitionR3E3_vs_SNP_S12_1035503.png](FiguresForReport/BLUPPercentage_inhibitionR3E3_vs_SNP_S12_1035503.png)

![FiguresForReport/BLUPArea_inhibitionR3E3_vs_SNP_S2_3164606.png](FiguresForReport/BLUPArea_inhibitionR3E3_vs_SNP_S2_3164606.png)

![FiguresForReport/BLUPPercentage_inhibitionR3E3_vs_SNP_S2_3164606.png](FiguresForReport/BLUPPercentage_inhibitionR3E3_vs_SNP_S2_3164606.png)



|Taxa                | SNPS2_3164606| SNPS12_1035503|
|:-------------------|-------------:|--------------:|
|21_CasTuj_L4.3      |          1.00|           0.00|
|21_ConRic_L3.1      |          1.00|           1.00|
|21_EciAmi_L1.3      |          1.00|           1.00|
|21_EciCon_L4.3      |          1.00|           0.00|
|21_EciRic_L4.1      |          1.00|           0.00|
|21_EciTej_L1.1      |          1.00|           1.00|
|21_EscAmi_L1.1      |          0.00|           0.00|
|21_EscRic_L1.2      |          1.00|           1.00|
|21_JerAmi_L1.3      |          1.00|           1.00|
|21_JerRic_L4.4      |          1.00|           0.00|
|21_ManTuj_L3.1      |          1.00|           0.00|
|22_Conil3758_L1     |          1.00|           1.00|
|22_Conil3806_L1     |          1.00|           0.00|
|22_Conil3911_L1     |          1.00|           0.00|
|22_ConilAmi_L1      |          1.00|           0.00|
|22_ConilCris_L1     |          1.00|           1.00|
|22_ConilFel_L1      |          1.00|           0.00|
|22_ConilFer_L1      |          0.00|           1.00|
|22_ConilMax_L1      |          1.00|           1.00|
|22_ConilOrt_L1      |          0.00|           0.00|
|22_ConilScul_L1     |          1.00|           1.00|
|22_ConilVit_L1      |          0.88|           0.54|
|22_Cor3806_L1       |          0.00|           0.00|
|22_Cor3911_L1       |          1.00|           1.00|
|22_Cor3927_L1       |          1.00|           1.00|
|22_CorAmi_L1        |          1.00|           0.00|
|22_CorAth_L1        |          1.00|           1.00|
|22_CorCale_L1       |          1.00|           1.00|
|22_CorCar_L1        |          1.00|           0.54|
|22_CorCris_L1       |          0.88|           0.54|
|22_CorFel_L1        |          1.00|           0.00|
|22_CorFer_L1        |          0.00|           0.00|
|22_CorIca200_L1     |          1.00|           1.00|
|22_CorKiko_L1       |          1.00|           0.00|
|22_CorMax_L1        |          1.00|           0.00|
|22_CorOrt_L1        |          0.00|           0.00|
|22_CorRic_L1        |          1.00|           1.00|
|22_CorScul_L1       |          0.00|           1.00|
|22_CorSim_L1        |          1.00|           0.00|
|22_CorVal_L1        |          1.00|           0.00|
|22_EcijaRegAmi_L1   |          0.00|           1.00|
|22_EcijaRegBon_L1   |          1.00|           0.00|
|22_EcijaRegKiko_L1  |          1.00|           1.00|
|22_EcijaRegRum_L1   |          1.00|           1.00|
|22_EcijaRegTej_L1   |          1.00|           0.00|
|22_EcijaSec3758_L1  |          1.00|           0.00|
|22_EcijaSec3806_L1  |          1.00|           0.00|
|22_EcijaSec3911_L1  |          1.00|           1.00|
|22_EcijaSec3927_L1  |          1.00|           1.00|
|22_EcijaSec83Ica_L2 |          1.00|           1.00|
|22_EcijaSecAmi_L1   |          1.00|           1.00|
|22_EcijaSecAth_L2   |          1.00|           1.00|
|22_EcijaSecBon_L1   |          1.00|           1.00|
|22_EcijaSecCris_L1  |          1.00|           0.00|
|22_EcijaSecEuro_L1  |          1.00|           1.00|
|22_EcijaSecFel_L1   |          0.00|           0.54|
|22_EcijaSecFer_L1   |          1.00|           0.00|
|22_EcijaSecJos_L1   |          1.00|           0.00|
|22_EcijaSecMax_L1   |          1.00|           0.00|
|22_EcijaSecNor_L2   |          1.00|           0.00|
|22_EcijaSecOrt_L1   |          1.00|           1.00|
|22_EcijaSecRic_L1   |          0.88|           0.54|
|22_EcijaSecRum_L1   |          1.00|           0.00|
|22_EcijaSecSah_L1   |          1.00|           1.00|
|22_EcijaSecSeb_L1   |          1.00|           0.00|
|22_EcijaSecSim_L2   |          1.00|           1.00|
|22_EcijaSecTej_L1   |          1.00|           1.00|
|22_EcijaSecVal_L1   |          0.00|           1.00|
|22_EcijaSecVit_L2   |          1.00|           0.54|
|22_Esc3911_L2       |          1.00|           1.00|
|22_Esc3992_L1       |          1.00|           1.00|
|22_EscCale_L1       |          1.00|           1.00|
|22_EscCar_L1        |          1.00|           1.00|
|22_EscCris_L1       |          1.00|           1.00|
|22_EscEuro_L1       |          1.00|           0.00|
|22_EscFer_L1        |          1.00|           0.00|
|22_EscIca_L1        |          1.00|           1.00|
|22_EscJos_L1        |          1.00|           1.00|
|22_EscMax_L1        |          0.00|           1.00|
|22_EscOrt_L1        |          1.00|           0.00|
|22_EscRic_L1.1      |          1.00|           1.00|
|22_EscSim_L1        |          1.00|           1.00|
|22_EscVal_L1        |          1.00|           1.00|
|22_EscVit_L2        |          1.00|           1.00|
|22_Jerez3911_L1     |          1.00|           0.00|
|22_Jerez3927_L1     |          1.00|           0.00|
|22_JerezEuro_L1     |          1.00|           1.00|
|22_JerezFel_L1      |          1.00|           0.00|
|22_JerezFer_L1      |          0.00|           1.00|
|22_JerezFue_L1      |          1.00|           1.00|
|22_JerezJos_L1      |          1.00|           1.00|
|22_JerezMax_L1      |          1.00|           0.00|
|22_JerezNor_L1      |          1.00|           0.00|
|22_JerezRic_L1      |          1.00|           1.00|
|22_JerezRum_L1      |          1.00|           0.00|
|22_JerezSim_L1      |          1.00|           0.00|
|22_JerezVal_L1      |          1.00|           0.00|
|22_MontanaRic_L1    |          1.00|           0.00|
|22_MonteraRic_L1    |          1.00|           1.00|
|22_ViveroRic_L1     |          1.00|           0.00|

#### Analysis 3: Filtered SNP Set (167938 SNPs) and filtered genotypes







|   |       X|SNP         | Chr|     Pos| P.value|       MAF|traits                           |
|:--|-------:|:-----------|---:|-------:|-------:|---------:|:--------------------------------|
|3  | 1595971|S12_1035503 |  12| 1035503|   0e+00| 0.4467535|BLINK.BLUPRatio_inhibitionR3E3   |
|4  |   45411|S1_988927   |   1|  988927|   0e+00| 0.0596668|BLINK.BLUPRatio_Scales1R3E3      |
|2  |    4541|S1_988927   |   1|  988927|   2e-07| 0.0596668|FarmCPU.BLUPRatio_Scales1R3E3    |
|1  |  159597|S12_1035503 |  12| 1035503|   3e-07| 0.4467535|FarmCPU.BLUPRatio_inhibitionR3E3 |


#### Conclusion

A brief synthesis of how the identified genetic markers contribute to our understanding of wheat's phenotypic variation and implications for breeding strategies.

## Conclusion

Our GWAS methodology, addressing both agronomic and "endophyte traits," employs a rigorous analytical framework to uncover the genetic underpinnings of phenotypic variation in wheat. The findings from this study are poised to contribute significantly to breeding strategies, advancing our genetic understanding and cultivation of wheat.
