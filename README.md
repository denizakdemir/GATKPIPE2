# GATKPIPE2

A bioinformatics pipeline for analyzing *Zymoseptoria tritici* sequencing data, closely following the GATK variant calling best practices (Van der Auwera & O'Connor, 2020; Poplin et al., 2017). 

## Pipeline Overview

1. **Quality Control and Trimming:**
   - **Decompression and Initial QC:** Raw FASTQ files were decompressed and subjected to quality control using FastQC, with results compiled via MultiQC.
   - **Trimming:** Using Trimmomatic, we trimmed adapters and low-quality sequences, followed by additional quality control.

2. **Read Alignment:**
   - **Alignment:** Reads were aligned to the *Z. tritici* reference genome (IPO323) from EnsemblGenomes.org using BWA-MEM.
   - **SAM to BAM Conversion:** SAM files were converted to sorted BAM files with SAMtools.
   - **Duplicate Marking:** Duplicate reads were marked, and metrics were collected using Picard.

3. **Variant Calling:**
   - **Raw VCF Generation:** GATK's HaplotypeCaller was employed to generate raw VCF files, applying stringent filters to identify high-confidence SNPs.
   - **Base Quality Score Recalibration (BQSR):** Systematic errors were corrected, and recalibrated BAM files were used to produce GVCF files.
   - **Joint Genotyping:** GVCF files were merged into a GenomicsDB using GATK's GenomicsDBImport tool, and joint genotyping across all samples was performed to create a consolidated VCF.
   - **SNP Annotation:** SnpEff was used for SNP annotation.

4. **Downstream Analysis:**
   - **Conversion and Haplotype Analysis:** TASSEL (Bradbury et al., 2007) was utilized to convert annotated VCFs into hapmap format, perform haplotype finding and imputation, and calculate kinship matrices and PCA for population structure assessment.
   - **Consensus Sequence Generation:** Consensus sequences for specific genomic regions were generated, providing insights into the genetic variability of *Z. tritici* strains.

This pipeline integrates robust tools and methodologies, adhering to bioinformatics best practices for high-throughput sequencing data analysis (McKenna et al., 2010; DePristo et al., 2011).

## Directory Structure

Folders other than `slurmscripts` are empty and are present to establish the directory structure. All the code is located in the `slurmscripts` folder.

## References

- Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120. [DOI: 10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170)

- McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297-1303. [DOI: 10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110)

- DePristo, M. A., Banks, E., Poplin, R., Garimella, K. V., Maguire, J. R., Hartl, C., Philippakis, A. A., del Angel, G., Rivas, M. A., Hanna, M., McKenna, A., Fennell, T. J., Kernytsky, A. M., Sivachenko, A. Y., Cibulskis, K., Gabriel, S. B., Altshuler, D., & Daly, M. J. (2011). A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics, 43(5), 491-498. [DOI: 10.1038/ng.806](https://doi.org/10.1038/ng.806)

- Poplin, R., Ruano-Rubio, V., DePristo, M. A., Fennell, T. J., Carneiro, M. O., Van der Auwera, G. A., Kling, D. E., Gauthier, L. D., Levy-Moonshine, A., Roazen, D., Shakir, K., Thibault, J., Chandran, S., Whelan, C., Lek, M., Gabriel, S., Daly, M. J., Neale, B., MacArthur, D. G., & Banks, E. (2017). Scaling accurate genetic variant discovery to tens of thousands of samples. bioRxiv, 201178. [DOI: 10.1101/201178](https://doi.org/10.1101/201178)

- Van der Auwera, G. A., & O'Connor, B. D. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.

- Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. (2007). TASSEL: Software for association mapping of complex traits in diverse samples. Bioinformatics 23:2633-2635. [DOI: 10.1093/bioinformatics/btm308](https://doi.org/10.1093/bioinformatics/btm308)
