#Plotting histograms for Tn5 vs CRISPRtOE data
#R commands to construct the histogram comparing Tn5 and CRISPRtOE distance around gene start
#2023-11-15

library(RColorBrewer)

#adjust working directory as needed
setwd("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/distance_measurements/")

#read in data
data <- read.table(file = "crisprtoe_up_and_down_gene_start.txt", head = F)
data2 <- read.table(file = "wholegenome_tn5_up_and_down_gene_start.txt", head = F)

#use specific color pallete
brewer.pal(n=5,"Set1")

#adjust colors
adjust_red <- adjustcolor("#E41A1C", alpha.f = 0.5)
adjust_green <- adjustcolor("#4DAF4A", alpha.f = 0.5)

#make histogram
hist(data$V1, breaks = 100, main = "Distance from Gene Start", xlab = "Distance (bp)", col=adjust_red, freq = F)
hist(data2$V1, breaks = 100, main = "Distance from Gene Start", xlab = "Distance (bp)", col=adjust_green, add=T, freq = F)
