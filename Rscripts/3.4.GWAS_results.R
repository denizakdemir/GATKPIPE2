
setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1EOdeUkPja2XuqDeJsIRWboaDvxqGTk52/Julian_Garcia_AbadilloVelasco/PROJECTS/PAOLA_Grape")
## libraries
library(tidyverse)
# load("GWAS_output/GWAS_object.RData")
#
# #this object is too big and redundant, will keep just the important (GWAS)
# sets = c("harvest","berry","post-harvest")
#
# clean.GWAS = list()
# for (set in sets){
#   data = gapit_output[[set]]$Blink
#   for (trait in names(data)){
#     clean.GWAS[[trait]] = data[[trait]]$GWAS
#   }
# }
# save(clean.GWAS, file = "GWAS_output/output.clean.GWAS.RData")

require(ggrepel)
source("./code/GWAS_humberto/fun/manhattanfunction.R")
source("./code/GWAS_humberto/fun/FDR_function.R")
source("./code/GWAS_humberto/fun/chromossomes_changes.R")

load("GWAS_output/output.clean.GWAS.RData")

SNP_table = list()
for (i in 1:length(clean.GWAS)){
  trait = names(clean.GWAS)[i]
  mod <- clean.GWAS[[i]][,c(1,2,3,4)]
  colnames(mod)<-c("SNP","Chromosome","Position","P.value")
  fdr_value<-FDR(mod$P.value, 0.05)
  mod$LOD <- -log10(mod$P.value)
  if(any(is.infinite(mod$LOD)) == TRUE){
    print("Infinity value was found for LOD, it will replace by 'LOD = max(LOD)+10', in order to plot!")
    pos_inf<-which(is.infinite(mod$LOD))
    max_val<-max(mod$LOD[is.finite(mod$LOD)])
    mod2<-mod
    mod2$LOD[pos_inf]<-max_val+10
    #
    plotGWAS <- manhattan_fun(df = mod2[,c(1,2,3,5)], ylims = c(0,(3+max(mod2[,5]))),
                              FDR=-log10(fdr_value), titleplot = paste0("A) Manhattan plot - ",trait))
    #
  } else {
    #
    plotGWAS <- manhattan_fun(df = mod[,c(1,2,3,5)], ylims = c(0,(3+max(mod[,5]))),
                              FDR=-log10(fdr_value), titleplot = paste0("A) Manhattan plot - ",trait))
  }
  #
  ggsave(filename = paste0("Manhattan_",trait,".pdf"),
         plot = plotGWAS, path = paste0("./GWAS_plots/"), scale = 1.8, dpi = 150, limitsize = FALSE)
  # QQ
  h<-length(mod$P.value)
  x<--log10(c(1:h)/h)
  x<-x[order(x)]
  y<--log10(mod$P.value[1:h])
  y<-y[order(y)]
  dfe <- cbind.data.frame(x,y)
  #
  qq <- ggplot2::ggplot(dfe, aes(x, y))+
    geom_point(alpha=0.8, size=2)+
    geom_smooth(method = lm, color = 'gray', fill = "#69b3a2", se = TRUE)+
    geom_hline(yintercept = -log10(0.05/dim(dfe)[1]), color = "blue", linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_value), color = "red", linetype = "dashed") +
    # stat_smooth(method = 'auto', color = 'red', fill = "#69b3a2", se = TRUE, n = dim(dfe)[1], size = 1)+
    labs(title=paste0("B) QQplot - ",trait), x = 'Theoretical [-log10(p)]', y = 'Observed [-log10(p)]') +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(size = 14),
      #plot.subtitle = element_text(size = 8, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
    )
  ggsave(filename=paste0("QQ_",trait,".pdf"),
         plot=qq,
         path="GWAS_plots",
         width = 12, height = 12, units = "cm", scale = 1.5, dpi = 150)
  ## List of significant SNPs
  snpdf <- mod[mod$LOD>=-log10(fdr_value), ]
  if(dim(snpdf)[1]==0){
    snpdf[1,] <- NA
    snpdf$FDR <- -log10(fdr_value)
    snpdf$Bonferroni <- -log10(0.05/dim(dfe)[1])
    snpdf$trait<-trait
  } else {
    snpdf$FDR <- -log10(fdr_value)
    snpdf$Bonferroni <- -log10(0.05/dim(dfe)[1])
    snpdf$trait<-trait
  }
  #
  SNP_table[[trait]]<-snpdf
}

save(SNP_table, file = paste0("./GWAS_output/associations.RData"))

assos = do.call("rbind",SNP_table)
assos = assos[order(assos$Chromosome,assos$Position),]
assos$bonfe_6.114818 =  assos$LOD > assos$Bonferroni
assos$Bonferroni = NULL
assos$P.value = NULL
assos = assos[,c(6,1,2,3,5,4,7)]
colnames(assos) = c("Trait","Marker","Chr","Pos","FDR","LOD","Bonfe>6.114818")
save(assos, file = "./GWAS_output/associations_clean.RData")

