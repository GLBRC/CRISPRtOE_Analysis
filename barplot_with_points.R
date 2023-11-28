#Working on Box Plot for CRISPRtOE
#2023-11-15

library(ggpubr)

#set working directory. will have to chagne
setwd("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/counting_analysis/")

data <- read.table(file = "box_plot_topGenes.txt", header = T, sep = "\t")

#use ggbarplot to make bargraph with dots and error bars
ggbarplot(data, x = "gene", y="value", title = "Average Top 10 Genes", ylab = "TMP vs DMSO Log2 Fold Change",
          xlab = "Gene (>2 Spacers per Gene)", add=c("mean_sd", "dotplot"),
          color = "black", fill = "grey", ylim=c(0,12))
