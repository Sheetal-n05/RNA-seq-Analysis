
# Project - RNA-seq data analysis and visualization

# Author - Sheetal Gulab Nighut
# comments - used anonymous cell types like A, B and E


# install essential package and libraries
install.packages("DESeq2")
install.packages("cowplot")
install.packages("devtools")
install.packages("pheatmap")
install.packages("dplyr")
install.packages("tidyr")
install.packages("tidyverse")

library(cowplot)
library(devtools)
library(dplyr)
library(ggplot2)
library(Matrix) 
library(pheatmap)
library(tidyr)
library(tidyverse)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)


View(counts(dds))

#The object class used by the DESeq2 package to store the read counts and the intermediate estimated quantities during 
#statistical analysis is the DESeqDataSet, which will usually be represented in the code here as an object dds .

colData(dds)
type(dds)
head(colData)

# A, B, C, E AND T ARE DIFFRENT CELL TYPES

#check dds colnames
colData(dds) 

#add dds column for sample number
colData(dds)$sample_number <- c(1:21)

# added sample_type 4rd column in colData(dds)
colData(dds)$sample_type <- c("A",	"A",	"A",	"B",	"B",	"B",	"T",	"T",	"T",	"D",	
                              "D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D", "E",	"F", "F") 

# added sample_subtype 5th column in colData(dds)
colData(dds)$sample_subtype <- c("X",	"X",	"X",	"Y",	"Y",	"Y",	"T",	"T",	"T",	"Z",
                                 "Z",	"Z",	"Z",	"Z",	"Z",	"Z",	"Z",	"Z",	"W",	"W",	"W") 
head(sample_type)


# rlog =This function transforms the count data to the log2 scale in a way which minimizes differences between samples 

# blind=TRUE should be used for comparing samples in an manner unbiased by prior information on samples, 

rld <- rlog(dds, blind=TRUE) #run rld (pca), default uses top 500 most variable genes
colData(rld)

# Plot Principle component analysis for sample similarity among all cell types

#PC1 vs PC2 colored by samples number
plotPCA(rld, intgroup="sample_number") + theme_classic()

#PC1 vs PC2 colored by samples type
plotPCA(rld, intgroup="sample_type") + theme_classic() 

#PC1 vs PC2 colored by samples subtype
plotPCA(rld, intgroup="sample_subtype") + theme_classic()

#PC1 vs PC2 colored by samples subtype
plotPCA(rld, intgroup="sample_subtype") + theme_classic() + xlim(70,80) + ylim(10,20)


#####
#Pearson correlation (r) heatmap
rld_mat <- assay(rld) #extract rlog matrix from object

#pearson pairwise correlation (default), could also choose spearman for rank-based measure of association
rld_cor <- cor(rld_mat, method = "pearson") 

#Pearson correlation: evaluates linear relationship between two continuous variables, most appropriate for measurements taken from an interval scale
#Spearman correlation: evaluates monotonic relationship, based on the ranked values for each variable rather than raw data,more appropriate for measurements taken from ordinal scales

head(rld_cor) #check corr table output

# keeping name pf colum annotation from colData(dds)$sample_type
annotation <- as.data.frame(colData(dds)$sample_type)

# putting rowname of annotation as colmnames in rld_cor i.e sample name i.e 21
colnames(annotation) <- "annotation"
rownames(annotation) <- colnames(rld_cor)


#plot heatmap of r values (i.e., pearson correlation) with hierarchical clustering by cell type
pheatmap(rld_cor, fontsize = 10, fontsize_row = 10, 
         height=20, border_color=NA, 
         annotation = annotation) 

#####
#Differential expression analysis

#build DESeqDataSet from count table plus metadata
#read in raw counts (rsem.merged.gene_counts.tsv form S3), not dds above, so can specify design by cell type

counts <- read.table(file = "/Users/snighut/Desktop/genesdata/rsem.merged.gene_counts.tsv", sep = '\t', header = TRUE)

# Showing roes 29,972 entries, 23 total columns ( geneid, transcriptid, sample 1-21)
head(counts)

rownames(counts) <- counts[,1] #genes as rownames
head(counts)

dim(counts)
counts <- cbind(counts[3:22]) #remove non-count columns
head(counts)

# import metadata file
metadata <- read.table(file = "/Users/snighut/Desktop/genesdata/metadata.csv", sep = '\t', header = TRUE)

#use dds colData above fore metadata
metadata <- colData(dds)

# check dim
dim(metadata)

#subset the metadata for all 18 samples
counts_subset <- counts[,1:18]
metadata_subset <- metadata[1:18,]

## differential expressed genes for cell type A vs T
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = round(counts_subset), #add round if integer error (even though raw counts are in intergers...)
  colData = metadata_subset,
  design =~ A) #construct data object from counts and metadata

ddsFullCountTable
colData(ddsFullCountTable)#check colData 
#run differential expression
dds_DEA_A <- DESeq(ddsFullCountTable)

#extract log2 fold changes and p-values for last variable in design formula (here only cell_type)
results_A <- results(dds_DEA_A, contrast = c("A", "A", "T")) 

#view table
results_A 


## differential expressed genes for cell type B vs E
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = round(counts_subset), #add round if integer error (even though raw counts are in intergers...)
  colData = metadata_subset,
  design =~ B) #construct data object from counts and metadata


ddsFullCountTable
colData(ddsFullCountTable) #check colData 
#run differential expression
dds_B <- DESeq(ddsFullCountTable) 
#extract log2 fold changes and p-values for last variable in design formula (here only cell_type)
results_B <- results(dds_DEA_B, contrast = c("B", "B", "E")) 
#view table
results_B 

# store results as dataframe
results_B <- as.data.frame(results_B)

results_B_DEA_table <- write.csv(results_B_DEA_table, file = "/Users/snighut/Desktop/genesdata/results_B_DEA_table.csv", row.names = TRUE)


###
##subset DEA tables to ID top DE genes (up and down regulated)
results_A_subset <- results_A[!is.na(results_A$padj) & results_A$padj<0.05,] #keep rows with padj not = NA or padj<0.05
results_B_subset <- results_B[!is.na(results_B$padj) & results_B$padj<0.05,]

#check dimensions
dim(results_A_subset)
dim(results_B_subset)


#keep rows with log2FoldChange not = NA or |log2FoldChange|>1
results_A_subset <- results_A_subset[!is.na(results_A_subset$log2FoldChange) & abs(results_A_subset$log2FoldChange)>1,]
results_B_subset <- results_B_subset[!is.na(results_B_subset$log2FoldChange) & abs(results_B_subset$log2FoldChange)>1,]

dim(results_A_subset) #check dimensions
dim(results_B_subset)

#add column of absolute value log2FoldChange = log2FoldChange_abs (all postive values)
results_A_subset$log2FoldChange_abs <- abs(results_A_subset$log2FoldChange)
results_B_subset$log2FoldChange_abs <- abs(results_B_subset$log2FoldChange)

#rank tables high to low log2FoldChange_abs
results_A_subset <- results_A_subset[order(-results_A_subset$log2FoldChange_abs),]
results_B_subset <- results_B_subset[order(-results_B_subset$log2FoldChange_abs),]


#rank tables high to low log2FoldChange
results_A_subset_up <- results_A_subset[order(-results_A_subset$log2FoldChange),] 
results_B_subset_up <- results_B_subset[order(-results_B_subset$log2FoldChange),]


#rank tables low to high log2FoldChange
results_A_subset_down <- results_A_subset[order(results_A_subset$log2FoldChange),]
results_B_subset_down <- results_B_subset[order(results_B_subset$log2FoldChange),]


#save results_xxx_subset tables as .csv files
write.csv(results_A_subset,"/Users/snighut/Desktop/genesdata/DEA_A.csv", row.names = TRUE)
write.csv(results_B_subset,"/Users/snighut/Desktop/genesdata/DEA_B.csv", row.names = TRUE)


#ID top 10 up-regulated diff exp genes per cell type
results_A_subset_up_10 <-  results_A_subset_up[1:10,]
results_B_subset_up_10 <- results_B_subset_up[1:10,]


#merge up tables and add column to ID cell type
up_combo <- as.data.frame(rbind(results_A_subset_up_10,results_B_subset_up_10))
up_combo$source <- c(replicate(10,"A"),replicate(10,"B"))

length(rownames(up_combo)) #check if rownames are unique --> YES
length(unique(rownames(up_combo))) 

#Normalize full dataset (counts) and subset to top up-regulated genes
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = metadata,
  design =~sample_type) #construct data object from counts and metadata
ddsFullCountTable <- estimateSizeFactors(ddsFullCountTable)
normalized_counts <- counts(ddsFullCountTable, normalized=TRUE)

#subset
normalized_counts_up_combo <- normalized_counts[rownames(normalized_counts)%in%rownames(up_combo),] 
dim(normalized_counts_up_combo) #check dims


#######

#heatmap up-regulated
annotation_up <- as.data.frame(up_combo$source)
rownames(annotation_up) <- rownames(normalized_counts_up_combo) #annotation of source

annotation_up2 <- as.data.frame(metadata$sample_type)
rownames(annotation_up2) <- colnames(normalized_counts_up_combo) #annotation of sample type


#plot heatmap, annotation_row=annotation_up
pheatmap(normalized_counts_up_combo, fontsize = 10, fontsize_row = 10, 
         height=20, border_color=NA, annotation=annotation_up2,cluster_rows=TRUE)


#ID top 10 down-regulated diff exp genes per cell type
results_A_subset_down_10 <-  results_A_subset_down[1:10,]
results_B_subset_down_10 <- results_BC_subset_down[1:10,]

#merge down tables and add column to ID cell type
down_combo <- as.data.frame(rbind(results_A_subset_down_10,results_B_subset_down_10))
down_combo$source <- c(replicate(10,"A"),replicate(10,"B"))

length(rownames(down_combo)) #check if rownames are unique --> YES
length(unique(rownames(down_combo))) 

#subset to top down-regulated genes
normalized_counts_down_combo <- normalized_counts[rownames(normalized_counts)%in%rownames(down_combo),] 

#check dims
dim(normalized_counts_down_combo)

#heatmap down-regulated
##annotation_down <- as.data.frame(down_combo$source)
##rownames(annotation_down) <- rownames(normalized_counts_down_combo) #annotation of source

annotation_down2 <- as.data.frame(metadata$sample_type)
rownames(annotation_down2) <- colnames(normalized_counts_down_combo) #annotation of sample type

##pheatmap(normalized_counts_down_combo, fontsize = 10, fontsize_row = 10, height=20, border_color=NA, annotation=annotation_down2,annotation_row=annotation_down,cluster_rows=FALSE) #plot heatmap
pheatmap(normalized_counts_down_combo, fontsize = 10, fontsize_row = 10, 
         height=20, border_color=NA, annotation=annotation_down2,cluster_rows=TRUE) 
#plot heatmap


#####
# Volcano Plot using DEA table

res <- as.data.frame (results_B)

P1 <- EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'padj',
                          title = 'DEG for celll type A vs T ',
                          subtitle = "cutoff padj<0.05, LogFC>1.5",
                          pCutoff = 0.05, #note padjvalue cutoff of 0.05
                          FCcutoff = 1.5, #log fold change cutoff of 1.5
                          pointSize = 2.0,
                          labSize = 3.0,
                          col=c('black', 'black', 'black', 'red3'), #show Differentially expressed genes as red dots on plot
                          colAlpha = 0.4,
                          vline =  c(-1.5, 1.5),
                          vlineCol = c('red3', 'red3'),
                          legendPosition = 'none', gridlines.major = FALSE, gridlines.minor = FALSE) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size=12))


#now plot if logfold change cutoff of 2.0
P2 <- EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'padj',
                      title = 'DEG for celll type A vs T',
                      subtitle = "cutoff padj<0.05, LogFC>2.0",
                      pCutoff = 0.05,
                      FCcutoff = 2.0,
                      pointSize = 2.0,
                      labSize = 3.0,
                      col=c('black','black','black','blue3'),
                      colAlpha = 0.4,
                      vline =  c(-2.0, 2.0),
                      vlineCol = c('blue3','blue3'),
                      legendPosition = 'none', gridlines.major = FALSE, gridlines.minor = FALSE) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size=12))



#now plot if logfold change cutoff of 3.0
P3 <- EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'padj',
                      title = 'DEG for celll type A vs T',
                      subtitle = "cutoff padj<0.05, LogFC>3.0",
                      pCutoff = 0.05,
                      FCcutoff = 3.0,
                      pointSize = 2.0,
                      labSize = 3.0,
                      col=c('black', 'black', 'black', 'green3'),
                      colAlpha = 0.4,
                      vline =  c(-3.0, 3.0),
                      vlineCol = c('green3', 'green3'),
                      legendPosition = 'none', gridlines.major = FALSE, gridlines.minor = FALSE) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size=12))

#combo plot
library(cowplot)

plot_grid(P1, P2, P3, ncol=3)

