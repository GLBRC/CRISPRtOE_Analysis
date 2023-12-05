#Commands for scatter plot for comparing replicates for CRISPRtOE
#R commands for plotting comparisons between replicates for both protospacers and gene
#2023-11-16

library(ggplot2)
library(ggpubr)
library(patchwork)

#adjust working directory as needed
setwd("/Users/kevinmyers/OneDrive - UW-Madison/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/counting_analysis/Scatter_Plots/")

############################################################################
##########            SPACERS            ###################################
############################################################################

#read in data
data <- read.table(file = "FC_NoProm_vs_Prom.txt", header = T, sep = "\t")
plot(data)
summary(data)

#use ggplot to plot data as needed
p <- ggplot(data, aes(x=NoProm, y=Prom)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "FC TMP vs DMSO Comparing Spacers for CRISPRtOE With and Without Promoters", x = "Log2 Fold Change TMP vs DMSO - NO Promoter CRISPRtOE",
       y = "Log2 Fold Change TMP vs DMSO - WITH Promoter CRISPRtOE") +
  xlim(-6,8) +
  ylim(-6,15) +
  stat_regline_equation(label.y = 15, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 14, aes(label = ..rr.label..)) +
  theme_classic()
p
#use ggsave to save files
ggsave("FC_Spacers_Prom_vs_NoProm.eps")
ggsave("FC_Spacers_Prom_vs_NoProm.pdf")

rep_data <- read.table(file = "spacer_replicates.txt", header = T, sep = "\t")
colnames(rep_data)
summary(rep_data)
log_rep_data <- log2(rep_data)
summary(log_rep_data)

fc_prom <- ggplot(data=subset(log_rep_data, !is.na(FC_Prom_A)), aes(x=FC_Prom_A, y=FC_Prom_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "FC TMP vs DMSO Comparing Spacers WITH PROMOTER CRISPRtOE Replicates", x = "Log2 Fold Change TMP vs DMSO A",
       y = "Log2 Fold Change TMP vs DMSO B") +
  xlim(-8,12) +
  ylim(-8,12) +
  stat_regline_equation(label.y = 11, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 10, aes(label = ..rr.label..)) +
  theme_classic()
fc_prom
ggsave("FC_Spacers_PromA_vs_PromB.eps")
ggsave("FC_Spacers_PromA_vs_PromB.pdf")

fc_noprom <- ggplot(data=subset(log_rep_data, !is.na(FC_NoProm_A)), aes(x=FC_NoProm_A, y=FC_NoProm_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "FC TMP vs DMSO Comparing Spacers NO PROMOTER CRISPRtOE Replicates", x = "Log2 Fold Change TMP vs DMSO A",
       y = "Log2 Fold Change TMP vs DMSO B") +
  xlim(-8,12) +
  ylim(-8,12) +
  stat_regline_equation(label.y = 11, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 10, aes(label = ..rr.label..)) +
  theme_classic()
fc_noprom
ggsave("FC_Spacers_NoPromA_vs_NoPromB.eps")
ggsave("FC_Spacers_NoPromA_vs_NoPromB.pdf")

dmso_noprom <- ggplot(data=subset(log_rep_data, !is.na(DMSO_NoProm_A)), aes(x=DMSO_NoProm_A, y=DMSO_NoProm_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "DMSO Comparing Spacers NO PROMOTER CRISPRtOE Replicates", x = "Log2 Spacer Counts DMSO A",
       y = "Log2 Spacer Counts DMSO B") +
  xlim(-2,18) +
  ylim(-2,18) +
  stat_regline_equation(label.y = 13, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 12, aes(label = ..rr.label..)) +
  theme_classic()
dmso_noprom
ggsave("DMSO_Spacers_NoPromA_vs_NoPromB.eps")
ggsave("DMSO_Spacers_NoPromA_vs_NoPromB.pdf")

dmso_prom <- ggplot(data=subset(log_rep_data, !is.na(DMSO_Prom_A)), aes(x=DMSO_Prom_A, y=DMSO_Prom_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "DMSO Comparing Spacers WITH PROMOTER CRISPRtOE Replicates", x = "Log2 Spacer Counts DMSO A",
       y = "Log2 Spacer Counts DMSO B") +
  xlim(-2,18) +
  ylim(-2,18) +
  stat_regline_equation(label.y = 13, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 12, aes(label = ..rr.label..)) +
  theme_classic()
dmso_prom
ggsave("DMSO_Spacers_PromA_vs_PromB.eps")
ggsave("DMSO_Spacers_PromA_vs_PromB.pdf")

tmp_noprom <- ggplot(data=subset(log_rep_data, !is.na(TMP_NoProm_A)), aes(x=TMP_NoProm_A, y=TMP_NoProm_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "TMP Comparing Spacers NO PROMOTER CRISPRtOE Replicates", x = "Log2 Spacer Counts TMP A",
       y = "Log2 Spacer Counts TMP B") +
  xlim(-2,18) +
  ylim(-2,18) +
  stat_regline_equation(label.y = 13, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 12, aes(label = ..rr.label..)) +
  theme_classic()
tmp_noprom
ggsave("TMP_Spacers_noPromA_vs_noPromB.eps")
ggsave("TMP_Spacers_noPromA_vs_noPromB.pdf")

tmp_prom <- ggplot(data=subset(log_rep_data, !is.na(TMP_Prom_A)), aes(x=TMP_Prom_A, y=TMP_Prom_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "TMP Comparing Spacers WITH PROMOTER CRISPRtOE Replicates", x = "Log2 Spacer Counts TMP A",
       y = "Log2 Spacer Counts TMP B") +
  xlim(-2,18) +
  ylim(-2,18) +
  stat_regline_equation(label.y = 13, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 12, aes(label = ..rr.label..)) +
  theme_classic()
tmp_prom
ggsave("TMP_Spacers_PromA_vs_PromB.eps")
ggsave("TMP_Spacers_PromA_vs_PromB.pdf")

fc_prom2 <- ggplot(data=subset(log_rep_data, !is.na(FC_Prom_A)), aes(x=FC_Prom_A, y=FC_Prom_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "With Promoter", x = "Log2 FC TMP vs DMSO A",
       y = "Log2 FC TMP vs DMSO B") +
  xlim(-8,12) +
  ylim(-8,12) +
  stat_regline_equation(label.y = 11, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 10, aes(label = ..rr.label..)) +
  theme_classic()
fc_prom2

fc_noprom2 <- ggplot(data=subset(log_rep_data, !is.na(FC_NoProm_A)), aes(x=FC_NoProm_A, y=FC_NoProm_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "No Promoter", x = "Log2 FC TMP vs DMSO A",
       y = "Log2 FC TMP vs DMSO B") +
  xlim(-8,12) +
  ylim(-8,12) +
  stat_regline_equation(label.y = 11, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 10, aes(label = ..rr.label..)) +
  theme_classic()
fc_noprom2

fc_prom2 + fc_noprom2 + plot_annotation(title = "Comparing TMP vs DMSO FC With or Without Promoter - Spacers")
ggsave("Combined_FC_Prom_noProm_SPACERS.eps")
ggsave("Combined_FC_Prom_noProm_SPACERS.pdf")


dmso_noprom2 <- ggplot(data=subset(log_rep_data, !is.na(DMSO_NoProm_A)), aes(x=DMSO_NoProm_A, y=DMSO_NoProm_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "DMSO NO PROMOTER", x = "Log2 Spacer Counts DMSO A",
       y = "Log2 Spacer Counts DMSO B") +
  xlim(-2,18) +
  ylim(-2,18) +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..)) +
  theme_classic()

dmso_prom2 <- ggplot(data=subset(log_rep_data, !is.na(DMSO_Prom_A)), aes(x=DMSO_Prom_A, y=DMSO_Prom_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "DMSO WITH PROMOTER", x = "Log2 Spacer Counts DMSO A",
       y = "Log2 Spacer Counts DMSO B") +
  xlim(-2,18) +
  ylim(-2,18) +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..)) +
  theme_classic()

tmp_noprom2 <- ggplot(data=subset(log_rep_data, !is.na(TMP_NoProm_A)), aes(x=TMP_NoProm_A, y=TMP_NoProm_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "TMP NO PROMOTER", x = "Log2 Spacer Counts TMP A",
       y = "Log2 Spacer Counts TMP B") +
  xlim(-2,18) +
  ylim(-2,18) +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..)) +
  theme_classic()

tmp_prom2 <- ggplot(data=subset(log_rep_data, !is.na(TMP_Prom_A)), aes(x=TMP_Prom_A, y=TMP_Prom_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "TMP WITH PROMOTER", x = "Log2 Spacer Counts TMP A",
       y = "Log2 Spacer Counts TMP B") +
  xlim(-2,18) +
  ylim(-2,18) +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..)) +
  theme_classic()

dmso_noprom2 + dmso_prom2 + tmp_noprom2 + tmp_prom2 + plot_annotation(title = "Comparing Replicates With or Without Promoter - Spacers")
ggsave("Grid_Spacer_Replicates.eps")
ggsave("Grid_Spacer_Replicates.pdf")


############################################################################
##########            GENES              ###################################
############################################################################


genes <- read.table(file = "median_genes_FC_promoter_vs_noPromoter.txt", header = T, sep = "\t", row.names = 1)
#log2_genes <- log2(genes)
head(genes)
summary(genes)

FC_gene_med <- ggplot(genes, aes(x=Log2_Median_FC_TMP_vs_DMSO_No_Promoter, y=Log2_Median_FC_TMP_vs_DMSO_Promoter)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "FC TMP vs DMSO Comparing Gene Medians for CRISPRtOE With and Without Promoters") +
  scale_x_continuous(
    name = "Log2 FC Gene Median TMP vs DMSO\nNO Promoter CRISPRtOE",
    limits = c(-6,8),
    breaks = c(-6, -4, -2, 0, 2, 4, 6, 8),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Log2 FC Gene Median TMP vs DMSO\nWITH Promoter CRISPRtOE",
    limits = c(-6,14),
    breaks = c(-6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14),
    expand = c(0,0)
  ) +
  stat_regline_equation(label.y = 10, label.x = -5, aes(label = ..rr.label..)) +
  theme_classic()
FC_gene_med
ggsave("FC_Genes_Prom_vs_NoProm_MEDIAN.eps")
ggsave("FC_Genes_Prom_vs_NoProm_MEDIAN.pdf")

noProm_spacers <- read.table(file = "Spacer_NO_Promoter_for_Median.txt", header = T, sep = "\t")
noProm_median <- aggregate(list(noProm_spacers$DMSO_NoPromoter_A_EC33,
                                noProm_spacers$DMSO_NoPromoter_B_EC34,
                                noProm_spacers$TMP_NoPromoter_A_EC45,
                                noProm_spacers$TMP_NoPromoter_A_EC45),
                           by=list(noProm_spacers$Locus_Tag), FUN = median)
colnames(noProm_median) <- c("Locus_Tag", "DMSO_NoPromoter_A_EC33_Median",
                             "DMSO_NoPromoter_B_EC34_Median",
                             "TMP_NoPromoter_A_EC45_Median",
                             "TMP_NoPromoter_A_EC45_Median")

write.table(noProm_median, file = "NoPromoter_Gene_Median_Values.txt", col.names = NA, sep = "\t", quote = F)

Prom_spacers <- read.table(file = "Spacer_Promoter_for_Median.txt", header = T, sep = "\t")

Prom_median <- aggregate(list(Prom_spacers$DMSO_Promoter_A_EC35,
                                Prom_spacers$DMSO_Promoter_B_EC36,
                                Prom_spacers$TMP_Promoter_A_EC47,
                                Prom_spacers$TMP_Promoter_B_EC48),
                           by=list(Prom_spacers$Locus_Tag), FUN = median)
colnames(Prom_median) <- c("Locus_Tag", "DMSO_Promoter_A_EC35_Median",
                             "DMSO_Promoter_B_EC36_Median",
                             "TMP_Promoter_A_EC47_Median",
                             "TMP_Promoter_B_EC48_Median")

write.table(Prom_median, file = "Promoter_Gene_Median_Values.txt", col.names = NA, sep = "\t", quote = F)

prom_values <- read.table(file = "Promoter_Gene_Median_Values.txt", header = T, sep = "\t", row.names = "Locus_Tag")
log2_prom_values <- log2(prom_values)
noProm_values <- read.table(file = "NoPromoter_Gene_Median_Values.txt", header = T, sep = "\t", row.names="Locus_Tag")
log2_noProm_values <- log2(noProm_values)

fc_gene_prom <- ggplot(data=log2_prom_values, aes(x=FC_Prom_A, y=FC_Prom_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "FC TMP vs DMSO Comparing Gene Medians WITH PROMOTER CRISPRtOE Replicates", x = "Log2 FC Gene Median\nTMP vs DMSO A",
       y = "Log2 FC Gene Median\nTMP vs DMSO B") +
  xlim(-8,12) +
  ylim(-8,12) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 9, aes(label = ..rr.label..)) +
  theme_classic()
fc_gene_prom
ggsave("FC_Genes_PromA_vs_PromB_MEDIAN.eps")
ggsave("FC_Genes_PromA_vs_PromB_MEDIAN.pdf")

fc_gene_noprom <- ggplot(data=log2_noProm_values, aes(x=FC_noProm_A, y=FC_noProm_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "FC TMP vs DMSO Comparing Gene Medians NO PROMOTER CRISPRtOE Replicates", x = "Log2 FC Gene Median\nTMP vs DMSO A",
       y = "Log2 FC Gene Median\nTMP vs DMSO B") +
  xlim(-8,12) +
  ylim(-8,12) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 9, aes(label = ..rr.label..)) +
  theme_classic()
fc_gene_noprom
ggsave("FC_Genes_NoPromA_vs_NoPromB_MEDIAN.eps")
ggsave("FC_Genes_NoPromA_vs_NoPromB_MEDIAN.pdf")

dmso_gene_noprom <- ggplot(data=log2_noProm_values, aes(x=DMSO_NoPromoter_A_EC33_Median, y=DMSO_NoPromoter_B_EC34_Median)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "DMSO Comparing Gene Medians NO PROMOTER CRISPRtOE Replicates", x = "Log2 Median Gene Counts DMSO A",
       y = "Log2 Median Gene Counts DMSO B") +
  xlim(-2,8) +
  ylim(-2,8) +
  stat_regline_equation(label.y = 7, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 6, aes(label = ..rr.label..)) +
  theme_classic()
dmso_gene_noprom
ggsave("DMSO_Genes_NoPromA_vs_NoPromB_MEDIAN.eps")
ggsave("DMSO_Genes_NoPromA_vs_NoPromB_MEDIAN.pdf")

dmso_gene_prom <- ggplot(data=log2_prom_values, aes(x=DMSO_Promoter_A_EC35_Median, y=DMSO_Promoter_B_EC36_Median)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "DMSO Comparing Gene Medians WITH PROMOTER CRISPRtOE Replicates", x = "Log2 Median Gene Counts DMSO A",
       y = "Log2 Median Gene Counts DMSO B") +
  xlim(-2,8) +
  ylim(-2,8) +
  stat_regline_equation(label.y = 7, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 6, aes(label = ..rr.label..)) +
  theme_classic()
dmso_gene_prom
ggsave("DMSO_Genes_PromA_vs_PromB_MEDIAN.eps")
ggsave("DMSO_Genes_PromA_vs_PromB_MEDIAN.pdf")

tmp_gene_noprom <- ggplot(data=log2_noProm_values, aes(x=TMP_NoPromoter_A_EC45_Median, y=TMP_NoPromoter_B_EC45_Median)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "TMP Comparing Gene Medians NO PROMOTER CRISPRtOE Replicates", x = "Log2 Median Gene Counts TMP A",
       y = "Log2 Median Gene Counts TMP B") +
  xlim(-2,10) +
  ylim(-2,10) +
  stat_regline_equation(label.y = 9, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 8, aes(label = ..rr.label..)) +
  theme_classic()
tmp_gene_noprom
ggsave("TMP_Genes_noPromA_vs_noPromB_MEDIAN.eps")
ggsave("TMP_Genes_noPromA_vs_noPromB_MEDIAN.pdf")

tmp_gene_prom <- ggplot(data=log2_prom_values, aes(x=TMP_Promoter_A_EC47_Median, y=TMP_Promoter_B_EC48_Median)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "TMP Comparing Gene Medians WITH PROMOTER CRISPRtOE Replicates", x = "Log2 Median Gene Counts TMP A",
       y = "Log2 Median Gene Counts TMP B") +
  xlim(-2,16) +
  ylim(-2,16) +
  stat_regline_equation(label.y = 13, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 12, aes(label = ..rr.label..)) +
  theme_classic()
tmp_gene_prom
ggsave("TMP_Genes_PromA_vs_PromB_MEDIAN.eps")
ggsave("TMP_Genes_PromA_vs_PromB_MEDIAN.pdf")

fc_gene_prom2 <- ggplot(data=log2_prom_values, aes(x=FC_Prom_A, y=FC_Prom_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "With Promoter", x = "Log2 FC TMP vs DMSO A",
       y = "Log2 FC TMP vs DMSO B") +
  xlim(-8,12) +
  ylim(-8,12) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 9, aes(label = ..rr.label..)) +
  theme_classic()
fc_gene_prom2

fc_gene_noprom2 <- ggplot(data=log2_noProm_values, aes(x=FC_noProm_A, y=FC_noProm_B)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "No Promoter", x = "Log2 FC TMP vs DMSO A",
       y = "Log2 FC TMP vs DMSO B") +
  xlim(-8,12) +
  ylim(-8,12) +
  stat_regline_equation(label.y = 11, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 10, aes(label = ..rr.label..)) +
  theme_classic()
fc_gene_noprom2

fc_gene_prom2 + fc_gene_noprom2 + plot_annotation(title = "Comparing TMP vs DMSO FC With or Without Promoter - Gene Medians")
ggsave("Combined_FC_Prom_noProm_GENES_MEDIAN.eps")
ggsave("Combined_FC_Prom_noProm_GENES_MEDIAN.pdf")


dmso_gene_noprom2 <- ggplot(data=log2_noProm_values, aes(x=DMSO_NoPromoter_A_EC33_Median, y=DMSO_NoPromoter_B_EC34_Median)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "DMSO NO PROMOTER", x = "Log2 Gene Counts DMSO A",
       y = "Log2 Gene Counts DMSO B") +
  xlim(-2,16) +
  ylim(-2,16) +
  stat_regline_equation(label.y = 15, aes(label = ..rr.label..)) +
  theme_classic()

dmso_gene_prom2 <- ggplot(data=log2_prom_values, aes(x=DMSO_Promoter_A_EC35_Median, y=DMSO_Promoter_B_EC36_Median)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "DMSO WITH PROMOTER", x = "Log2 Gene Counts DMSO A",
       y = "Log2 Gene Counts DMSO B") +
  xlim(-2,16) +
  ylim(-2,16) +
  stat_regline_equation(label.y = 15, aes(label = ..rr.label..)) +
  theme_classic()

tmp_gene_noprom2 <- ggplot(data=log2_noProm_values, aes(x=TMP_NoPromoter_A_EC45_Median, y=TMP_NoPromoter_B_EC45_Median)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "TMP NO PROMOTER", x = "Log2 Gene Counts TMP A",
       y = "Log2 Gene Counts TMP B") +
  xlim(-2,16) +
  ylim(-2,16) +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..)) +
  theme_classic()

tmp_gene_prom2 <- ggplot(data=log2_prom_values, aes(x=TMP_Promoter_A_EC47_Median, y=TMP_Promoter_B_EC48_Median)) +
  geom_point(size = 2, shape = 20) +
  geom_smooth(method = "lm", se = F, level=0.95, color="#E41A1C") +
  labs(title = "TMP WITH PROMOTER", x = "Log2 Gene Counts TMP A",
       y = "Log2 Gene Counts TMP B") +
  xlim(-2,16) +
  ylim(-2,16) +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..)) +
  theme_classic()

dmso_gene_noprom2 + dmso_gene_prom2 + tmp_gene_noprom2 + tmp_gene_prom2 + plot_annotation(title = "Comparing Replicates With or Without Promoter - Gene Medians")
ggsave("Grid_Gene_Replicates_MEDIAN.eps")
ggsave("Grid_Gene_Replicates_MEDIAN.pdf")
