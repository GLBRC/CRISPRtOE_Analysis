#Working on Box Plot for CRISPRtOE
#2023-11-15

library(ggpubr)

#set working directory. will have to change
setwd("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/counting_analysis/")

data <- read.table(file = "box_plot_topGenes.txt", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
ggbarplot(data, x = "gene", y="value", title = "Average Top 10 Genes", ylab = "TMP vs DMSO Log2 Fold Change",
          xlab = "Gene (>2 Spacers per Gene)", add=c("mean_sd", "jitter"),
          color = "black", fill = "grey", ylim=c(0,12))

#################################
#Promoter and No Promoter Values#
#################################

#set working directory. will have to change
setwd("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/counting_analysis/")

data2 <- read.table(file = "box_plot_topGenes_prom_noProm.txt", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
ggbarplot(data2, x = "gene", y="value", title = "Average Top 10 Genes", ylab = "TMP vs DMSO Log2 Fold Change",
          xlab = "Gene (>2 Spacers per Gene)", add=c("mean_sd", "jitter"),
          pallette = c("red", "blue"), position = position_dodge(),
          color = "group", fill = "lightgrey", ggtheme = theme_classic()) +
  scale_y_continuous("TMP vs DMSO Log2 Fold Change",limits=c(-4,12), breaks = c(-4,-2,0,2,4,6,8,10,12)) +
  geom_hline(yintercept=0, color = "black")

ggsave("Prom_vs_noProm_top_genes.eps")
ggsave("Prom_vs_noProm_top_genes.pdf") 
