#Constructing FDR vs Fold Change Plot
#Notes on how to produce the scatter plots for CRISPRtOE analysis
#CRISPRtOE Project
#2023-11-13

library(ggrepel)
library(tidyverse)

#adjust working directory
setwd("/Users/kevinmyers/OneDrive - UW-Madison/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/counting_analysis/")

#read in table for plotting
prom5 <- read.table(file = "promoter_minCount5_updated_for_scatter_plot.txt", header = TRUE, sep = "\t")

#use ggplot for plotting the scatter plot
p <- ggplot(prom5, aes(prom5$negLogFDR, prom5$FC_linear)) +
  geom_point(aes(color=prom5$group2), size=1) +
  scale_color_manual(values=c("red", "gray")) +
  labs(title = "TMP vs DMSO, Promoter Constructs, Minmum 5 counts per spacer") +
  geom_text(data=subset(prom5, Gene == "b0048"), aes(negLogFDR, FC_linear, label = Gene), color = "black", nudge_x = -0.5) +
  scale_x_continuous(
    name = "-log(FDR)",
    limits = c(0.00001,100),
    breaks = c(0.00001,0.0001,0.001,0.01,0.1,1,10,100),
    expand = c(0,0),
    trans='log10',
  ) +
  scale_y_continuous(
    name = "Spacer Count Fold Change (TMP/DMSO)",
    limits = c(0.01,10000),
    breaks = c(0.01,0.1,1,10,100,1000,10000),
    expand = c(0,0),
    trans='log10'
  ) +
  theme_classic() +
  theme(legend.title=element_blank())
p
