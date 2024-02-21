
setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1EOdeUkPja2XuqDeJsIRWboaDvxqGTk52/Julian_Garcia_AbadilloVelasco/PROJECTS/PAOLA_Grape")
## libraries
require(tidyverse)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
#require(GAPIT3)
require(AssocTests)

# load data2gwas
#load("data/geno_beagle/pre_GWAS_geno.RData")
#load("data/BLUEs_v2.RData")

# # numeric K
# taxas = paste("G",rownames(G_grape),sep="_")
# newK = data.frame(Taxa = taxas, G_grape)
# rownames(newK) = 1:nrow(newK)
# colnames(newK) = c("Taxa",newK$Taxa)
# G_grape = newK
# rm(newK)


# fix Q
#Q$Taxa = paste("G",Q$Taxa,sep="_")

#
# # rearrange map
# map = map[,c(3,1,2)]
# map$POS = as.numeric(map$POS)
# sum(map$ID == colnames(snp_matrix)[2:ncol(snp_matrix)])
# # fix chr
# fix_chr = function(x){
#   output = c()
#   for (i in x){
#     if (is.na(i)){
#       i = NA
#     }
#     else{
#     if (nchar(i)==1){
#       i = paste("0",i,sep="")
#       }
#     }
#     output = c(output,i)
#   }
#   return(output)
# }
# #
# mapping$CHROM = fix_chr(as.character(mapping$CHROM))
# mapping$ID = paste("M",mapping$CHROM,mapping$POS, sep="_")
#mapping$CHROM[is.na(mapping$CHROM)] = "20"
#mapping$CHROM = as.numeric(mapping$CHROM)

# marker matrix
# newM = data.frame(Taxa = taxas, M_grape)
# rownames(newM) = 1:nrow(newM)
# colnames(newM) = c("Taxa",mapping$ID)
# M_grape = newM
# rm(newM)
# colnames(phenoT[[1]])[1] = colnames(phenoT[[2]])[1] = colnames(phenoT[[3]])[1] ="Taxa"
#colnames(M_grape) = c("Taxa",mapping$ID)

load("data/GWAS_input_v2.RData")
#order check
#sum(G_grape$Taxa == colnames(G_grape)[2:ncol(G_grape)]) #G names are correct
#sum(G_grape$Taxa == M_grape$Taxa) # G and M names match
#sum(G_grape$Taxa == Q$Taxa) # G and Q names match
#sum(mapping$Name == colnames(M_grape)[2:ncol(M_grape)]) # SNP ID match in map and M
#match(mapping$Name, colnames(M_grape)[2:ncol(M_grape)])
#identical(mapping$Name, colnames(M_grape)[2:ncol(M_grape)])

#mapping = data.frame(mapping)
#na.indices = which(is.na(mapping$CHROM))
#mapping[na.indices,]$CHROM = 20
#mapping$ID = paste("M",mapping$CHROM, mapping$POS, sep = "_")


#
# load("Data/geno_beagle/curated_lines.RData")
# curated_lines = paste("G",unique(final$Taxa),sep="_")
#

# seedless_lines = unique(c(new.data.pre %>% filter(seedless) %>% .$varieties,
#                     new.data.berry %>% filter(seedless) %>% .$varieties,
#                     new.data.post %>% filter(seedless) %>% .$varieties))

# save(G_grape, M_grape, mapping, phenoT, Q,
#       curated_lines, seedless_lines,
#         file = "data/GWAS_input_v2.RData")

models = c("Blink") #
set.names = c("harvest","berry","post-harvest")
seed.traits = c("seed_number","seed_fresh_weight") #remove seedless GIDs for these

gapit_output <- list() # overwritten take a look!

for (k in 1:length(phenoT)){
  data = phenoT[[k]] %>%
    filter(Taxa %in% curated_lines)
  for(i in 1:length(models)){
    for(j in 1:length(colnames(data)[-1])){ #
       s <- getwd()
       f <- paste0(getwd(), paste0("/GWAS_PC6_output/"))
       y2 <- colnames(data)[1+j]
       print(y2)

       y <- strsplit(y2,":")[[1]][1]
       dir.create(paste0(f,'/GAPIT_',models[i],'_',y,"/"), recursive = TRUE)
       setwd(paste0(f,'/GAPIT_',models[i],'_',y,"/"))
       # Y
       Y <- data.frame(data[c(1,(j+1))])
       if (y2 %in% seed.traits){
         Y = filter(Y,!Taxa %in% seedless_lines)
       }
       # model
       mod <- GAPIT(Y=Y,
                GD=M_grape,
                GM=mapping,
                KI=G_grape,
                CV=NULL,
                PCA.total=6,
                model=models[i])

      gapit_output[[y2]] <- mod$GWAS
      setwd(s)
    }
  }
}

save(gapit_output, file = "./GWAS_PC6_output/GWAS_object.RData")

