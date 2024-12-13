# CRISPRtOE Scripts and Commands

## Description and Reference

Repository of scripts and commands used to process and analyze sequencing files for CRISPRtOE (CRISPR transposon over-expression) experiments.

This repository maybe updated to streamline the analysis.

This is associated with the manuscript MANUSCRIPT TITLE AND LINK

## Requirements

There are four [Anaconda](https://anaconda.org) environments required for different steps of the analysis. Each environment contains the specific versions of software used for each step as listed below. The environment files are:

- `cutadaptenv.yml`
- `samtools.yml`
- `r-3.6.yml`
- `bedops_env.yml`
- `deeptools_env.yml`

These are found in the `environments` directory and can be installed as follows:

`conda env create -f cutadaptenv.yml`

The specific R libraries required for plotting are:

- `edgeR`
- `ggplot2`
- `tidyverse`
- `ggrepel`

## Usage and Steps

### Combine and trim FASTQ files

1. Combine FASTQ files if run over multiple lanes, keeping the R1 and R2 files separate:<br>
`cat Sample1*R1*fastq > Sample1_R1_combined.fastq`<br>
`cat Sample1*R2*fastq > Sample1_R2_combined.fastq`<br>

2. Remove all non-genomic sequence from the R1 and R2 FASTQ files with [Cutadapt (v3.4)](https://journal.embnet.org/index.php/embnetjournal/article/view/200):<br>
`conda activate cutadaptenv`<br>
`cutadapt -a TGTTGGAACAACCATAAAATGATAATTACACCCATAAA -o Sample1_R1_combined_TnMatchFiltered.fastq Sample1_R1_combined.fastq`<br>
`cutadapt -g GGATCCGTTATCAGCTACCTACTCGGCAGTTCAC -o Sample1_R2_combined_filter1.fastq Sample1_R2_combined.fastq`<br>
`cutadapt -a GTTATCAGCTACCTACTCGGC -o Sample1_R2_combined_filter2.fastq Sample1_R2_combined_filter1.fastq`<br>

### Align filtered FASTQ files to genome with [Bowtie2 (v2.4.4)](http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html)

Use [Bowtie2 (v2.4.4)](http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html) and the [_E. coli_ K-12 MG1655 genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)<br>

`conda activate cutadaptenv`<br>
`bowtie2 -N 0 -x ecoli -1 Sample1_R1_combined_TnMatchFiltered.fastq -2 Sample1_R2_combined_filter2.fastq --un Sample1_unmapped.sam -S Sample1_out.sam`<br>

### Filter the alignment files for length and same strand alignment of both reads

Use [Samtools (v1.13)](https://academic.oup.com/bioinformatics/article/25/16/2078/204688) and [BBMap (v38.32)](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/) for the filtering steps

1. Filter for only mapped reads:<br>
`conda activate cutadaptenv`<br>
`samtools view -F 8 -h -O sam Sample1_out.sam -o Sample1_out_mapped.sam`<br>

2. Filter for reads that map to the same strand:<br>
`samtools view -F 0x2 -h -O sam Sample1_out_mapped.sam -o Sample1_out_mapped_filtered.sam`

3. Filter for reads with a minimum length of 2 nt and a maximum length of 50 nt using [BBMap (v38.32)]:<br>
`bbmap-38.32/reformat.sh in=Sample1_out_mapped_filtered.sam out=Sample1_out_mapped_filtered_lengthFiltered.sam minlength=2 maxlength=50`

4. Generate a historgram from the final filtered SAM file using `parsing_sam_for_histogram_single_file.py` with a text file listing the SAM files to process and then the Rscript `Histogram_CRISPRtOE.R`:<br>
`python parsing_sam_for_histogram_single_file.py -f sam_files.txt`<br>
`Rscript Histogram_CRISPRtOE.R -d <directory with parsed SAM files>`

The histogram will show the distance between the spacer alignment and the transposon insertion location.

### Filter to include reads where the distance between the spacer and the transposon insertion site is <100 nts

One of the files from the `Histogram_CRISPRtOE.R` Rscript is `SampleA_less_than_100.txt`. This file is a list of the reads with a distance <100 nts. We can use this to filter the SAM file.

1. Convert the SAM files to BAM files, then use AWK to refine the text file:<br>
`conda activate samtools`<br>
`samtoosl view -b Sample1_out_mapped_filtered_lengthFiltered.sam -o Sample1_out_mapped_filtered_lengthFiltered.bam`<br>
`awk -F'\t' '{print $2}' SampleA_less_than_100.txt > SampleA_less_than_100_readName.txt`<br>

2. Use [Picard tools (v2.20)](https://broadinstitute.github.io/picard/) to filter the BAM file and include ONLY reads where the spacer sequence is less than 100 nts from the transposon insertion site:<br>
`conda activate samtools`<br>
`picard -Xmx124g FilterSamReads I=Sample1_out_mapped_filtered_lengthFiltered.bam O=Sample1_out_mapped_filtered_lengthFiltered_100bpInt_filtered.bam READ_LIST_FILE=SampleA_less_than_100_readName.txt FILTER=includeReadList`<br>

### Construct BigWig and Bed Files using [DeepTools (v3.5.1)](https://academic.oup.com/nar/article/44/W1/W160/2499308)

1. Sort and index the BAM file:<br>
`conda activate samtools`<br>
`samtools sort -O bam -o Sample1_out_mapped_filtered_lengthFiltered_sort.bam Sample1_out_mapped_filtered_lengthFiltered_100bpInt_filtered.bam`<br>
`samtools index Sample1_out_mapped_filtered_lengthFiltered_sort.bam`<br>

2. Normalized and write BigWig and Bed files:<br>
`conda activate deeptools_env`<br>
`bamCoverage -bs 1 -b Sample1_out_mapped_filtered_lengthFiltered_sort.bam -o Sample1_out_mapped_filtered_lengthFiltered_sort_BinSize1_CPM_norm.bw --effectiveGenomeSize 4641652 --normalizeUsing CPM`<br>
`bamCoverage -bs 1 -b Sample1_out_mapped_filtered_lengthFiltered_sort.bam -o $Sample1_out_mapped_filtered_lengthFiltered_sort_BinSize1.bed --outFileFormat bedgraph`<br>

### Filter BAM file to determine exact location of transposon insertion site

1. Filter the BAM file to only include reads that are <25 nt long. This will remove all the spacer reads because they are 32 nt long. Then the resulting reads can be used to find the location of the transposon insertion across the genome:<br>
`conda activate samtools`<br>
`samtools view -e 'length(seq)<25' -O BAM -o Sample1_out_mapped_filtered_lengthFiltered_sort_lessThan25.bam Sample1_out_mapped_filtered_lengthFiltered_sort.bam`<br>

2. Combine the BAM and corresponding BED files into a file (tab delimited) called `bam_files.txt` and then run the following python script:<br>
`ls *bam > bam_files.txt`<br>
`python site_of_insertion.py -f bam_files.txt`<br>

This will result in a WIG file for visualization and a table.txt file with the location of each insertion site and the count of the reads at that location.

### Determine the distance from transposon insertion site to gene downstream

Use [BEDOPS (v2.4.41)](https://academic.oup.com/bioinformatics/article/28/14/1919/218826) and the `closest-features` function for this.

1. Convert the table.txt file from the previous step into a BED file:<br>
`conda activate bedops_env`<br>
`ls *site_table.txt > input_file.txt`<br>
`python table_to_bed.py -f input_file.txt`<br>

2. Sort BED file and run `closest-features` with _E. coli_ gene BED file:<br>
`sort_bed Sample1_out_mapped_filtered_lengthFiltered_sort_lessThan25.bed > Sample1_out_mapped_filtered_lengthFiltered_sort_lessThan25_sort.bed`<br>
`closest-features --dist --closest Sample1_out_mapped_filtered_lengthFiltered_sort_lessThan25_sort.bed NC_000913.3.bed > Sample1_out_mapped_filtered_lengthFiltered_sort_lessThan25_sort_dist.bed`<br>

### Using DIST.BED files to construct histogram of distances between transposon insertion site and gene start site

1. Combine sample DIST.BED files into one and make a table to construct a histogram:<br>
`cat *dist.bed > combined_crisprtoe.bed`<br>
`python make_table_for_histogram.py`<br>

Use Tn-seq data from [Goodall et al](https://journals.asm.org/doi/10.1128/mbio.02096-17) from _E. coli_ K-12 BW25113 as a comparison.

2. Use R to construct a histogram comparing CRISPRtOE distance upstream of genes to traditional Tn-seq data:<br>
`CRISPRtOE <- read.table(file = "combined_crisprtoe_distance.txt", head = F)`<br>
`tnseq <- read.table(file = "whole_genome_tnseq_sort_dist.txt", header = F)`<br>
`hist(CRISPRtOE$V1, breaks = 50, main = "Distance from Gene Start", xlab = "Distance (bp)", col=rgb(1,0,0,0.5), freq = F)`<br>
`hist(tnseq$V1, breaks = 50, main = "Distance from Gene Start", xlab = "Distance (bp)", col=rgb(0,0,1,0.5), add=T, freq = F)`<br>

### Count the number of times each spacer appears in the filtered alignment files

Use `SampleA_combined_out_mapped_lengthFiltered_100bpInt_filtered_sort.bam` files. Combine if more than one:<br>
`ls *bam > bam_files.txt`<br>
`python counting_spacers.py -f bam_files.txt -t spacers_to_search.txt`<br>

Use the `comparison_of_spacer_counts_edgeR.R` Rscript to compare control and treated samples.<br>
Use the `plotting_scatter_plot.R` Rscript to generate a plot comparing fold change and FDR values.<br>
Use the `barplot_with_points.R` Rscript to generate bar plots of gene fold change with individual data points superimposed.<br>
Use the `scatter_plot_commands.R` Rscript to generate correlation scatter plots for replicates and other comparisons.

### Repository Maintenance

Repo maintained by Kevin Myers (kmyers2@wisc.edu) at https://gitpub.wei.wisc.edu/kmyers/crisprtoe_analysis (Access restricted for primary repo)
