Tassell Manual:

1- open filtered_snps.vcf.gz using files>open
2- Save filtered_snps files>save
3- thin by position
data>thin by position, filtered_snps_thinned_1000
3- Snp summary
    data>Geno summary
    
4- Filter Taxa
    filter>genotype table taxa 0.6, 0.0, 1.0
5- Filter snps:
    filter>genotype table locations  0.05, .95 maf maf
6- Numeric Geno
    file>Numeric Genotype
7- Impute

