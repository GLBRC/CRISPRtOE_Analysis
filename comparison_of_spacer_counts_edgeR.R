#Comparing spacer counts with edgeR
#CRISPRtOE Project
#2023-11-06

library(edgeR)
library(tidyverse)
#Set working directory
setwd("~/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/counting_analysis/edgeR/")
#Read in HTseq (raw count) data and targets file
all <- read.table(file = "spacer_count_output_file_first.txt", header = TRUE, row.names = 1)
targets <- read.table(file = "targets.txt", header = TRUE, sep = "\t")
#Group by Strain and Time Point since unpaired and add to the targets file
group <- factor(targets$Strain)
cbind(targets, group = group)
#Construct the design matrix using the group and no intercept
design <- model.matrix(~0 + group)
#filter out spacer IDs where the DMSO samples (33 and 34 OR 35 and 36) are 0s
filtered_zeros_DMSO <- filter(all, CtOE.EC.33 > 0 & CtOE.EC.34 > 0
                   & CtOE.EC.35 > 0 & CtOE.EC.36 > 0)

filtered_minCount_5_noProm <- filter(all, CtOE.EC.33 > 5 & CtOE.EC.34 > 5
                                & CtOE.EC.45 > 5 & CtOE.EC.46 > 5
                                  )

filtered_minCount_5_Prom <- filter(all, CtOE.EC.35 > 5 & CtOE.EC.36 > 5
                                     & CtOE.EC.47 > 5 & CtOE.EC.48 > 5
                                      )

#Construct the DGE for zeros in DMSO
dge <- DGEList(counts = filtered_zeros_DMSO, group=group)
dge <- calcNormFactors(dge)
#Plot updated MDS plot after all that - Write to file
pdf(file = "./MDS_Corrected_Reps_CRISPRtOE_DMSOzerosRemoved.pdf")
plotMDS(dge, main = "MDS Plot of Replicates for DMSO and TMP CRISPRtOE\n(DMSO Zeros Removed)")
dev.off()
#Estimate the dispersions - the heart of edgeR
x<-estimateGLMCommonDisp(dge, design)
y<-estimateGLMTrendedDisp(x, design)
z<-estimateGLMTagwiseDisp(y, design)
fit<-glmQLFit(z, design)
#Make the contrasts, or comparisons to make
my.contrasts<-makeContrasts(NoPromoter = groupTMP_NoPromoter - groupDMSO_NoPromoter,
                            Promoter = groupTMP_Promoter - groupDMSO_Promoter,
                            levels = design)

#Set up the DE test using the desired comparison
NoPromoter <- glmQLFTest(fit, contrast = my.contrasts[, "NoPromoter"])
#Determine DE genes and use BH as the method to determine FDR
de <- decideTestsDGE(NoPromoter, adjust.method="BH", p.value = 0.05)
detags <- rownames(NoPromoter)[as.logical(de)]
#Plot PDF of DE genes and Fold Change boundaries - Write to file
pdf(file = "./NoPromoter_DMSOzerosRemoved.pdf")
plotSmear(NoPromoter, de.tags=detags, main = "NoPromoter TMP vs DMSO (DMSO Zeros Removed)", xlab = "Log Counts Per Million (CPM)", ylab = "Log Fold Change (FC)", ylim=c(-15,15))
abline(h=c(-1,1), col = "blue")
abline(h=c(-2,2), col = "green")
legend(x = "topright", c("FDR <= 0.05", "FDR >= 0.05", "LogFC = 1", "LogFC = 2"), cex = 0.8, col = c("red", "black", "blue", "green"), lty = c(NA, NA, 1, 1), pch = c(20, 20, NA, NA))
dev.off()
#Print the top DE by FDR to the screen
topTags(NoPromoter)
#To export the data with the FDR, calculate the FDR separately and add it to the table of results
FDR <- p.adjust(NoPromoter$table$PValue, method = "BH")
table <- cbind(NoPromoter$table, FDR)
sum(table[,"FDR"]<= 0.05)
#Write results table to file
write.table(table, file = "NoPromoter_edgeR_results_DMSOzerosRemoved.txt", quote = FALSE, sep = "\t", col.names = NA)

Promoter <- glmQLFTest(fit, contrast = my.contrasts[, "Promoter"])
de <- decideTestsDGE(Promoter, adjust.method="BH", p.value = 0.05)
detags <- rownames(Promoter)[as.logical(de)]
pdf(file = "./Promoter_DMSOzerosRemoved.pdf")
plotSmear(Promoter, de.tags=detags, main = "Promoter TMP vs DMSO (DMSO Zeros Removed)", xlab = "Log Counts Per Million (CPM)", ylab = "Log Fold Change (FC)", ylim = c(-15,15))
abline(h=c(-1,1), col = "blue")
abline(h=c(-2,2), col = "green")
legend(x = "bottomright", c("FDR <= 0.05", "FDR >= 0.05", "LogFC = 1", "LogFC = 2"), cex = 0.8, col = c("red", "black", "blue", "green"), lty = c(NA, NA, 1, 1), pch = c(20, 20, NA, NA))
dev.off()
topTags(Promoter)
FDR <- p.adjust(Promoter$table$PValue, method = "BH")
table <- cbind(Promoter$table, FDR)
sum(table[,"FDR"]<= 0.05)
write.table(table, file = "Promoter_edgeR_results_DMSOzerosRemoved.txt", quote = FALSE, sep = "\t", col.names = NA)

cpmData <- cpm(filtered_zeros_DMSO)  #normalize counts for counts per million (same as normalized DGE above)
write.table(cpmData, file = "normalized_spacer_count_output_file_DMSOzerosRemoved.txt", sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)  #write normalized counts results to file

#########################################################################
#####       Set up the DE test using the desired comparison         #####
#########################################################################

dge_noProm <- DGEList(counts = filtered_minCount_5_noProm, group=group)
dge_noProm <- calcNormFactors(dge_noProm)
#Plot updated MDS plot after all that - Write to file
pdf(file = "./MDS_Corrected_Reps_CRISPRtOE_minCount5.pdf")
plotMDS(dge, main = "MDS Plot of Replicates for DMSO and TMP CRISPRtOE\n(Minimum 5 Counts per Spacer)")
dev.off()
#Estimate the dispersions - the heart of edgeR
x<-estimateGLMCommonDisp(dge_noProm, design)
y<-estimateGLMTrendedDisp(x, design)
z<-estimateGLMTagwiseDisp(y, design)
fit<-glmQLFit(z, design)
#Make the contrasts, or comparisons to make
my.contrasts<-makeContrasts(NoPromoter = groupTMP_NoPromoter - groupDMSO_NoPromoter,
                            levels = design)
NoPromoter <- glmQLFTest(fit, contrast = my.contrasts[, "NoPromoter"])
#Determine DE genes and use BH as the method to determine FDR
de <- decideTestsDGE(NoPromoter, adjust.method="BH", p.value = 0.05)
detags <- rownames(NoPromoter)[as.logical(de)]
#Plot PDF of DE genes and Fold Change boundaries - Write to file
pdf(file = "./NoPromoter_minCount5.pdf")
plotSmear(NoPromoter, de.tags=detags, main = "NoPromoter TMP vs DMSO\n(Minimum 5 Counts per Spacer)", xlab = "Log Counts Per Million (CPM)", ylab = "Log Fold Change (FC)", ylim=c(-15,15))
abline(h=c(-1,1), col = "blue")
abline(h=c(-2,2), col = "green")
legend(x = "topright", c("FDR <= 0.05", "FDR >= 0.05", "LogFC = 1", "LogFC = 2"), cex = 0.8, col = c("red", "black", "blue", "green"), lty = c(NA, NA, 1, 1), pch = c(20, 20, NA, NA))
dev.off()
#Print the top DE by FDR to the screen
topTags(NoPromoter)
#To export the data with the FDR, calculate the FDR separately and add it to the table of results
FDR <- p.adjust(NoPromoter$table$PValue, method = "BH")
table <- cbind(NoPromoter$table, FDR)
sum(table[,"FDR"]<= 0.05)
#Write results table to file
write.table(table, file = "NoPromoter_edgeR_results_minCount5.txt", quote = FALSE, sep = "\t", col.names = NA)

dge_Prom <- DGEList(counts = filtered_minCount_5_Prom, group=group)
dge_Prom <- calcNormFactors(dge_Prom)
#Plot updated MDS plot after all that - Write to file
pdf(file = "./MDS_Corrected_Reps_CRISPRtOE_minCount5.pdf")
plotMDS(dge, main = "MDS Plot of Replicates for DMSO and TMP CRISPRtOE\n(Minimum 5 Counts per Spacer)")
dev.off()
#Estimate the dispersions - the heart of edgeR
x<-estimateGLMCommonDisp(dge_Prom, design)
y<-estimateGLMTrendedDisp(x, design)
z<-estimateGLMTagwiseDisp(y, design)
fit<-glmQLFit(z, design)
#Make the contrasts, or comparisons to make
my.contrasts<-makeContrasts(Promoter = groupTMP_Promoter - groupDMSO_Promoter,
                            levels = design)
Promoter <- glmQLFTest(fit, contrast = my.contrasts[, "Promoter"])
de <- decideTestsDGE(Promoter, adjust.method="BH", p.value = 0.05)
detags <- rownames(Promoter)[as.logical(de)]
pdf(file = "./Promoter_minCount5.pdf")
plotSmear(Promoter, de.tags=detags, main = "Promoter TMP vs DMSO\n(Minimum 5 Counts per Spacer)", xlab = "Log Counts Per Million (CPM)", ylab = "Log Fold Change (FC)", ylim = c(-15,15))
abline(h=c(-1,1), col = "blue")
abline(h=c(-2,2), col = "green")
legend(x = "bottomright", c("FDR <= 0.05", "FDR >= 0.05", "LogFC = 1", "LogFC = 2"), cex = 0.8, col = c("red", "black", "blue", "green"), lty = c(NA, NA, 1, 1), pch = c(20, 20, NA, NA))
dev.off()
topTags(Promoter)
FDR <- p.adjust(Promoter$table$PValue, method = "BH")
table <- cbind(Promoter$table, FDR)
sum(table[,"FDR"]<= 0.05)
write.table(table, file = "Promoter_edgeR_results_minCount5_onlyProm.txt", quote = FALSE, sep = "\t", col.names = NA)

cpmData_Prom <- cpm(filtered_minCount_5_Prom)  #normalize counts for counts per million (same as normalized DGE above)
write.table(cpmData_Prom, file = "normalized_spacer_count_output_file_minCount5_promoter.txt", sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)  #write normalized counts results to file

cpmData_noProm <- cpm(filtered_minCount_5_noProm)  #normalize counts for counts per million (same as normalized DGE above)
write.table(cpmData_noProm, file = "normalized_spacer_count_output_file_minCount5_noPromoter.txt", sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)  #write normalized counts results to file
