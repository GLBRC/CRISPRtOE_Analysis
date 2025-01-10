#Histogram of distance of CRISPRtOE
#Code for graphing histograms for CRISPRtOE and distance from insertion site to gene start
#2023-10-19

library(optparse)

#use optparse to import directory to work from
parser <- OptionParser()

option_list <- list(
  make_option(c("-d", "--directory"), help = "Directory containing SAM files", type = "character", action = 'store', dest = "dir")
)

opt <- parse_args(OptionParser(option_list=option_list))

setwd(opt$dir)

#find all SAM files in working directory
SAM_files <- Sys.glob(file.path(opt$dir, "*position_list.txt"))
print(SAM_files)

#process each SAM file and make the histogram plots
for(i in SAM_files){
  print(i)
  sample_name <- strsplit(i, split='_out')[[1]][1]
  data <- read.table(file = i, header = TRUE, sep = "\t")
  data_1000 <- subset(data, data$distance_between_reads <= 1000)
  data_100 <- subset(data, data$distance_between_reads <= 100)
  first_pdf <- paste(sample_name,"all_interval_values.pdf",sep="_")
  second_pdf <- paste(sample_name,"less_than_1000bp_interval_values.pdf",sep="_")
  third_pdf <- paste(sample_name,"less_than_1000bp_focused_interval_values.pdf",sep="_")
  data_out <- paste(sample_name,"less_than_100.txt",sep="_")
  title_ID <- paste(sample_name,"Combined Histogram of End of Spacer\nto Transposon Insertion Site (≤ 1000 bp)",sep=" ")
  pdf(file = first_pdf)
  hist(data$distance_between_reads, breaks = 50000, freq = T, main = "Histogram of End of Spacer to Transposon Insertion Site", xlab = "Distance (bp)")
  dev.off()
  pdf(file = second_pdf)
  hist(data_1000$distance_between_reads, breaks = 100, xlim = c(0,1000), freq = T, main = "", xlab = "Distance (bp)")
  title(main = title_ID)
  axis(side=1, at=seq(0,1000,50), labels = seq(0,1000,50))
  dev.off()  
  pdf(file = third_pdf)
  hist(data_1000$distance_between_reads, breaks = 1000, xlim = c(0,100), freq = T, xlab = "Distance (bp)", main = "")
  title(main = title_ID)
  axis(side=1, at=seq(0,100,5), labels = seq(0,100,5))
  dev.off()
  write.table(data_100, file = data_out, col.names = NA, quote = FALSE, sep = "\t")
}