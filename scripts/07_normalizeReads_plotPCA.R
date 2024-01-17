#!/usr/bin/env Rscript

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")

# BiocManager::install("DESeq2")
# install.packages("SummarizedExperiment")
# install.packages("ggplot2")
# install.packages("ggpubr")

library(ggplot2)

# define function to normalize read counts and get a PCA analysis
normalizeReads_plotPCA <- function(raw_counts_file) {
  
  # load raw count files
  write("     Loading files...", stdout())
  raw_counts <- read.csv(raw_counts_file, row.names = 1)
  
  # load metadata and order the "Age" column
  metadata <- read.table("00_input/SRA_metadata.tsv", header = TRUE, sep = "\t")
  metadata$Age <- factor(metadata$Age, levels = c("0_hpf", "4_hpf", "8_hpf", "12_hpf", "16_hpf", "20_hpf", "24_hpf",
                                                  "28_hpf", "32_hpf", "36_hpf", "40_hpf", "44_hpf", "48_hpf",
                                                  "52_hpf", "72_hpf"))
  
  # generate DESeq data object
  DESeq.ds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ Age)
  
  # keep genes with more than 10 reads in at least 2 runs
  keep <- rowSums(DESeq2::counts(DESeq.ds) >= 20) >= 2
  DESeq.ds <- DESeq.ds[keep,]
  
  # normalize read count (rlog transformation) and export table
  write("     Normalizing read count with rlog transformation...", stdout())
  suppressMessages(rld <- DESeq2::rlog(DESeq.ds, blind = FALSE))
  rld_table <- SummarizedExperiment::assay(rld)
  
  # normalize read count (vst transformation) and export table
  write("     Normalizing read count with vst transformation...", stdout())
  vst <- DESeq2::vst(DESeq.ds, blind = FALSE)
  vst_table <- SummarizedExperiment::assay(vst)
  
  # perform PCA with rld transformed data
  write("     Plotting results...", stdout())
  PCA.rld <- DESeq2::plotPCA(rld, intgroup = "Age", returnData = TRUE)
  
  #extract PC variance
  variance <- attr(PCA.rld, "percentVar")
  
  # plot PCA
  PCA.rld.plot <- ggplot(PCA.rld, aes(PC1, PC2, color = group)) +
    geom_point(size = 2, alpha = 0.7) +
    xlab(paste("PC1 (", round(variance[1]*100, 2), "%)", sep = "")) +
    ylab(paste("PC2 (", round(variance[2]*100, 2), "%)", sep = "")) +
    theme_bw() +
    ggtitle("Rlog transformed counts")
  
  # perform PCA with vst transformed data
  PCA.vst <- DESeq2::plotPCA(vst, intgroup = "Age", returnData = TRUE)
  
  #extract PC variance
  variance <- attr(PCA.vst, "percentVar")
  
  # plot PCA
  PCA.vst.plot <- ggplot(PCA.vst, aes(PC1, PC2, color = group)) +
    geom_point(size = 2, alpha = 0.7) +
    xlab(paste("PC1 (", round(variance[1]*100, 2), "%)", sep = "")) +
    ylab(paste("PC2 (", round(variance[2]*100, 2), "%)", sep = "")) +
    theme_bw() +
    ggtitle("Vst transformed counts")
  
  
  panel <- ggpubr::ggarrange(PCA.rld.plot, PCA.vst.plot,
                             ncol = 1, nrow = 2,
                             common.legend = TRUE, legend = "right")
  
  return(list(rld_table = rld_table, vst_table = vst_table, panel = panel, pca = PCA.rld))
}


write("", stdout())
write("--- Analyzing gene raw counts ---", stdout())
gene_counts_results <- normalizeReads_plotPCA("03_mapped_reads/gene_count_matrix.csv")

write("     Saving normalized read counts to files...", stdout())
write.table(gene_counts_results$rld_table, file = "03_mapped_reads/gene_counts_rlogTransformed.tsv", sep = "\t", quote = FALSE)
write.table(gene_counts_results$vst_table, file = "03_mapped_reads/gene_counts_vstTransformed.tsv", sep = "\t", quote = FALSE)

write("     Saving plot to file...", stdout())
ggsave("03_mapped_reads/PCA_generawcounts.pdf", plot = gene_counts_results$panel, device = "pdf",
       dpi = 300, height = 6, width = 5.5, units = ("in"))

write("DONE", stdout())

write("", stdout())
write("--- Analyzing transcript raw counts ---", stdout())
transcript_counts_results <- normalizeReads_plotPCA("03_mapped_reads/transcript_count_matrix.csv")

write("     Saving normalized read counts to files...", stdout())
write.table(transcript_counts_results$rld_table, file = "03_mapped_reads/transcript_counts_rlogTransformed.tsv", sep = "\t", quote = FALSE)
write.table(transcript_counts_results$vst_table, file = "03_mapped_reads/transcript_counts_vstTransformed.tsv", sep = "\t", quote = FALSE)

write("     Saving plot to file...", stdout())
ggsave("03_mapped_reads/PCA_transcriptrawcounts.pdf", plot = transcript_counts_results$panel, device = "pdf",
       dpi = 300, height = 6, width = 5.5, units = ("in"))

write("DONE", stdout())
