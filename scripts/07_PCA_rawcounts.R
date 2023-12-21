#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly=TRUE))
  BiocManager::install("DESeq2")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("SummarizedExperiment")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("ggpubr")

library(ggplot2)

# load raw count files
print("Loading files...")
raw_counts <- read.csv("03_mapped_reads/gene_count_matrix.csv", row.names = 1)

# load metadata and order the "Age" column
metadata <- read.table("00_input/SRA_metadata.tsv", header = TRUE, sep = "\t")
metadata$Age <- factor(metadata$Age, levels = c("0_hpf", "4_hpf", "8_hpf", "12_hpf", "16_hpf", "20_hpf", "24_hpf",
                                                "28_hpf", "32_hpf", "36_hpf", "40_hpf", "44_hpf", "48_hpf",
                                                "52_hpf", "72_hpf"))

# generate DESeq data object
DESeq.ds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ Age)

# investigate library sizes
# colSums(DESeq2::counts(DESeq.ds))

# read count normalization (rlog transformation)
print("Normalizing read count with rlog transformation...")
rld <- DESeq2::rlog(DESeq.ds, blind = FALSE)

# read count normalization (vst transformation)
print("Normalizing read count with vst transformation...")
vst <- DESeq2::vst(DESeq.ds, blind = FALSE)

# PCA plot
print("Plotting results...")
PCA.rld <- DESeq2::plotPCA(rld, intgroup = "Age", returnData = TRUE)
PCA.rld.plot <- ggplot(PCA.rld, aes(PC1, PC2, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_bw() +
  ggtitle("Rlog transformed counts")

PCA.vst <- DESeq2::plotPCA(vst, intgroup = "Age", returnData = TRUE)
PCA.vst.plot <- ggplot(PCA.vst, aes(PC1, PC2, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_bw() +
  ggtitle("Vst transformed counts")


panel <- ggpubr::ggarrange(PCA.rld.plot, PCA.vst.plot,
                           ncol = 1, nrow = 2,
                           common.legend = TRUE, legend = "right")

print("Saving plot to file...")
ggsave("03_mapped_reads/PCA_rawcounts.pdf", plot = panel, device = "pdf",
       dpi = 300, height = 6, width = 5.5, units = ("in"))

print("BYE BYE!!")
