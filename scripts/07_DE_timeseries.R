#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly=TRUE))
  BiocManager::install("maSigPro")

library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)


#####################
#     FUNCTIONS     #
#####################

# define function to normalize read counts and get a PCA analysis

normalizeReads_plotPCA <- function(raw_counts_file, runs_metadata_file) {

  # load raw count files
  write("     Loading files...", stdout())
  raw_counts <- read.csv(raw_counts_file, row.names = 1)

  # load metadata and order the "Age" column
  metadata <- read.table(runs_metadata_file, header = TRUE, sep = "\t")
  metadata$Age <- factor(metadata$Age, levels = c("0_hpf", "4_hpf", "8_hpf", "12_hpf", "16_hpf", "20_hpf", "24_hpf",
                                                  "28_hpf", "32_hpf", "36_hpf", "40_hpf", "44_hpf", "48_hpf",
                                                  "52_hpf", "72_hpf"))

  # generate DESeq data object
  DESeq.ds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ Age)

  # normalize read counts
  write("     Normalizing read count with DESeq2 median of ratios...", stdout())
  DESeq.ds <- DESeq2::estimateSizeFactors(DESeq.ds)

  return(DESeq.ds)
}

# function to run maSigPro on single timeseries data
timeseries_maSigPro <- function(norm_counts,edesign) {
  
  # load normalized read counts (output from 07_PCA_rawcounts.R)
  normalized.counts <- read.table(norm_counts, header = TRUE, row.names = 1, sep ="\t")
  
  # load edesign object
  edesign_object <- read.table(edesign, header = TRUE, sep = "\t", row.names = 1)
  
  # define regression model
  design <- maSigPro::make.design.matrix(edesign_object)
  
  # find significant genes
  fit <- maSigPro::p.vector(normalized.counts, design, counts = TRUE)
  fit <- maSigPro::p.vector(normalized.counts, design)

  # find significant differences
  tstep <- maSigPro::T.fit(fit)

  # obtain lists of significant genes
  get <- maSigPro::get.siggenes(tstep, vars = "all")

  # visualization
  seegenes <- maSigPro::see.genes(get$sig.genes, k.mclust = TRUE)
  
  return(list(counts = normalized.counts,
              edesign = edesign_object,
              sig_genes_list = get,
              visualization_genes = seegenes
              ))
}


# function to plot expression values of selected genes
plot_expression_values <- function(counts,gene_list,gene_names,edesign) {
  
  # get gene expression values
  expr <- as_tibble(counts[gene_list,], rownames = NA)
  
  # create tibble to plot results
  expr_to_plot <- rownames_to_column(expr, var = "cds") %>%                     # include rownames in tibble
    pivot_longer(cols = -c(cds), names_to = "run") %>%                          # make pivot longer
    left_join(rownames_to_column(edesign, var = "run")) %>%                     # add information on hours
    select(-c("Replicate","Group")) %>%                                         # drop unwanted columns
    left_join(gene_names, keep = FALSE) %>%                                     # add information on gene names
    drop_na()                                                                   # drop rows with NAs
  
  expr_to_plot$name <- factor(expr_to_plot$name, levels = c("Vasa", "Dmrt1L", "SoxH", "FoxL2", "Wnt8a", "FoxB2"))
  
  plot <- ggplot(data = expr_to_plot, aes(x = as.factor(Time), y = value, color = name, linetype = name)) +
    stat_summary(aes(group = gene), fun = mean, geom = "line", linewidth = 0.8) +
    geom_point(size = 2, alpha = 0.5) +
    xlab("Hours post fertilization (hpf)") +
    ylab("Normalized expression value")
  
  return(plot)
}


#------ACTUAL CODE------------------------------------------------------------

########################################
#     RUN READ COUNT NORMALIZATION     #
########################################

write("", stdout())

write("--- Analyzing gene raw counts ---", stdout())

normalized_reads <- normalizeReads_plotPCA("03b_mapped_reads_BOWTIE/ALL.rawmapping.stats.csv", "00_input/SRA_metadata.tsv")


########################
#     RUN MASIGPRO     #
########################

# run maSigPro on vst normalized gene counts from all timepoints
DE_alltimepoints <- timeseries_maSigPro(normalized_reads, "00_input/edesign_object.tsv")

# run maSigPro on vst normalized gene counts from timepoints up to 24hpf
DE_24hpf <- timeseries_maSigPro(normalized_reads, "00_input/edesign_object_24hpf.tsv")


####################################
#     CHECK DE OF TARGET GENES     #
####################################

# vst norm, all timepoints
vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B017427"] # vasa:    no
vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B093608"] # dmrt1L:  no
vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B014180"] # soxH:    no
vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B094018"] # foxL2:   YES, cluster 3
vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B085403"] # wnt8a:   YES, cluster 9

# vst norm, all timepoints
vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI03912.1"] # vasa:    YES, cluster 3
vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI03798.1"] # dmrt1L:  no
vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI30824.1"] # soxH:    no
vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI49864.1"] # foxL2:   YES, cluster 4
vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI54402.1"] # wnt8a:   no
vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI35942.1"] # foxB2:   no

# vst norm, 24hpf
vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B017427"] # vasa:    YES, cluster 6
vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B093608"] # dmrt1L:  no
vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B014180"] # soxH:    no
vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B094018"] # foxL2:   no
vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B085403"] # wnt8a:   YES, cluster 2

# rlog norm, all timepoints
rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B017427"] # vasa:   no
rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B093608"] # dmrt1L: no
rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B014180"] # soxH:   no
rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B094018"] # foxL2:  YES, cluster 6
rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B085403"] # wnt8a:  YES, cluster 1

# rlog norm, 24hpf
rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B017427"] # vasa:   YES, cluster 7
rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B093608"] # dmrt1L: no
rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B014180"] # soxH:   no
rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B094018"] # foxL2:  no
rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B085403"] # wnt8a:  YES, cluster 4


################################
#     PLOT GENE EXPRESSION     #
################################

# generate gene conversion dataframe
gene_names <- data.frame(gene = c("gene-MGAL_10B017427","gene-MGAL_10B093608","gene-MGAL_10B014180","gene-MGAL_10B094018","gene-MGAL_10B085403","gene-MGAL_10B093191"),
                         stringtie = c("MSTRG.9718", "MSTRG.9631", "MSTRG.26773", "MSTRG.38958", "MSTRG.41862", "MSTRG.30071"),
                         cds = c("Mgal_VDI03912.1", "Mgal_VDI03798.1", "Mgal_VDI30824.1", "Mgal_VDI49864.1", "Mgal_VDI54402.1", "Mgal_VDI35942.1"),
                         name = c("Vasa", "Dmrt1L", "SoxH", "FoxL2", "Wnt8a", "FoxB2"))


vst.expr.plot <- plot_expression_values(vstnorm_alltimepoints$counts,
                                        gene_names$cds,
                                        gene_names,
                                        vstnorm_alltimepoints$edesign) +
  scale_colour_manual(values = c("Red", "Cyan", "Magenta", "Yellow", "Grey", "Grey")) +
  scale_linetype_manual(values = c(rep("solid", 4), rep("dashed", 2))) +
  ggtitle("Vst transformed gene expression") +
  labs(color = "Gene", linetype = "Gene") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        axis.line = element_line(linewidth = 0.8),
        axis.ticks = element_line(color = "black", size = 0.8),
        axis.text = element_text(color = "black"),
        legend.background = element_rect(fill = 'transparent', color = NA),
        legend.text = element_text(face = "bold"))

vst.expr.plot

ggsave("03b_mapped_reads_BOWTIE/vstgene_expression.pdf", plot = vst.expr.plot, device = "pdf",
       dpi = 300, height = 4.7, width = 6, units = ("in"), bg = 'transparent')

ggsave("03b_mapped_reads_BOWTIE/vstgene_expression.png", plot = vst.expr.plot, device = "png",
       dpi = 300, height = 4.7, width = 6, units = ("in"), bg = 'transparent')



rlog.expr.plot <- plot_expression_values(rlognorm_alltimepoints$counts,
                                        gene_names$cds,
                                        gene_names,
                                        rlognorm_alltimepoints$edesign) +
  ggtitle("Rlog transformed gene expression")


panel <- ggpubr::ggarrange(rlog.expr.plot, vst.expr.plot,
                           ncol = 2, nrow = 1,
                           common.legend = TRUE, legend = "right")
panel
