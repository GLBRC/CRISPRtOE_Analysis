#Working on Box Plot for CRISPRtOE
#2023-11-15

library(ggpubr)
library(patchwork)

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

data2 <- read.table(file = "filtered_promoter_noPromoter_outlier_gene_spacer_data.txt", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
ggbarplot(data2, x = "gene", y="value", title = "Average Top 10 Genes", ylab = "TMP vs DMSO Log2 Fold Change",
          xlab = "Gene (>2 Spacers per Gene)", add=c("mean_sd", "jitter"),
          pallette = c("red", "blue"), position = position_dodge(),
          color = "group", fill = "lightgrey", ggtheme = theme_classic()) +
  scale_y_continuous("TMP vs DMSO Log2 Fold Change",limits=c(-4,12), breaks = c(-4,-2,0,2,4,6,8,10,12)) +
  geom_hline(yintercept=0, color = "black")

ggsave("Prom_vs_noProm_top_genes_updated.eps")
ggsave("Prom_vs_noProm_top_genes_updated.pdf") 

data_firstGene <- read.table(file = "filtered_promoter_noPromoter_outlier_gene_spacer_data_firstGene_Operon.tx", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
ggbarplot(data2, x = "gene", y="value", title = "Average Top 10 First Genes in Operon", ylab = "TMP vs DMSO Log2 Fold Change",
          xlab = "Gene (>2 Spacers per Gene)", add=c("mean_sd", "jitter"),
          pallette = c("red", "blue"), position = position_dodge(),
          color = "group", fill = "lightgrey", ggtheme = theme_classic()) +
  scale_y_continuous("TMP vs DMSO Log2 Fold Change",limits=c(-4,12), breaks = c(-4,-2,0,2,4,6,8,10,12)) +
  geom_hline(yintercept=0, color = "black")

ggsave("Prom_vs_noProm_top_genes_updated_firstGene.eps")
ggsave("Prom_vs_noProm_top_genes_updated_firstGene.pdf") 

########################################
#Promoter and No Promoter Values MEDIAN#
########################################

#set working directory. will have to change
setwd("/Users/kevinmyers/OneDrive - UW-Madison/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/counting_analysis/")

data2 <- read.table(file = "filtered_promoter_noPromoter_outlier_gene_spacer_data_medians.txt", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
ggbarplot(data2, x = "gene", y="value", title = "Median Top 10 Genes", ylab = "TMP vs DMSO Log2 Fold Change",
          xlab = "Gene (>2 Spacers per Gene)", add=c("median_mad", "jitter"),
          pallette = c("red", "blue"), position = position_dodge(),
          color = "group", fill = "lightgrey", ggtheme = theme_classic()) +
  scale_y_continuous("TMP vs DMSO Log2 Fold Change",limits=c(-4,12), breaks = c(-4,-2,0,2,4,6,8,10,12)) +
  geom_hline(yintercept=0, color = "black")

ggsave("Prom_vs_noProm_top_genes_MEDIAN.eps")
ggsave("Prom_vs_noProm_top_genes_MEDIAN.pdf") 

data_firstGene <- read.table(file = "filtered_FIRST_GENE_promoter_noPromoter_outlier_gene_spacer_data_medians.txt", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
ggbarplot(data_firstGene, x = "gene", y="value", title = "Median Top 10 First Genes in Operon", ylab = "TMP vs DMSO Log2 Fold Change",
          xlab = "Gene (>2 Spacers per Gene)", add=c("median_mad", "jitter"),
          pallette = c("red", "blue"), position = position_dodge(),
          color = "group", fill = "lightgrey", ggtheme = theme_classic()) +
  scale_y_continuous("TMP vs DMSO Log2 Fold Change",limits=c(-4,12), breaks = c(-4,-2,0,2,4,6,8,10,12)) +
  geom_hline(yintercept=0, color = "black")

ggsave("Prom_vs_noProm_top_genes_updated_MEDIAN_firstGene.eps")
ggsave("Prom_vs_noProm_top_genes_updated_MEDIAN_firstGene.pdf") 

panelA <- read.table(file = "PanelA_Supp_Figure_genes.txt", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
PanelA_plot <- ggbarplot(panelA, x = "gene", y="value", title = "Panel A", ylab = "TMP vs DMSO Log2 Fold Change",
          xlab = "Gene (>2 Spacers per Gene)", add=c("median_mad", "jitter"),
          pallette = c("red", "blue"), position = position_dodge(),
          color = "group", fill = "lightgrey", ggtheme = theme_classic()) +
  scale_y_continuous("TMP vs DMSO Log2 Fold Change",limits=c(-4,8), breaks = c(-4,-2,0,2,4,6,8)) +
  geom_hline(yintercept=0, color = "black")
PanelA_plot

ggsave("PanelA_sup_fig.eps")
ggsave("PanelA_sup_fig.pdf") 

panelB <- read.table(file = "PanelB_Supp_Figure_genes.txt", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
PanelB_plot <- ggbarplot(panelB, x = "gene", y="value", title = "Panel B", ylab = "TMP vs DMSO Log2 Fold Change",
                         xlab = "Gene (>2 Spacers per Gene)", add=c("median_mad", "jitter"),
                         pallette = c("red", "blue"), position = position_dodge(),
                         color = "group", fill = "lightgrey", ggtheme = theme_classic()) +
  scale_y_continuous("TMP vs DMSO Log2 Fold Change",limits=c(-6,4), breaks = c(-6,-4,-2,0,2,4)) +
  geom_hline(yintercept=0, color = "black")
PanelB_plot

ggsave("PanelB_sup_fig.eps")
ggsave("PanelB_sup_fig.pdf") 
