#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  BiocManager::install("maSigPro")
}

library(tidyverse)
library(maSigPro)

library(GOSemSim)
library(smacof)

setwd("/data/evassvis/fn76/mgal_DE_development")


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
timeseries_maSigPro <- function(norm_counts, edesign) {
  
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
plot_expression_values <- function(counts, gene_list, gene_names, edesign) {
  
  # get gene expression values
  expr <- as_tibble(counts[gene_list,], rownames = NA)
  
  # create tibble to plot results
  expr_to_plot <- rownames_to_column(expr, var = "cds") %>%                     # include rownames in tibble
    pivot_longer(cols = -c(cds), names_to = "run") %>%                          # make pivot longer
    left_join(rownames_to_column(edesign, var = "run")) %>%                     # add information on hours
    select(-c("Replicate", "Group")) %>%                                        # drop unwanted columns
    left_join(gene_names, keep = FALSE)                                         # drop rows with NAs
  
  # expr_to_plot$name <- factor(expr_to_plot$name, levels = c("Vasa", "Piwia", "Piwib", "Nanos", "spPHI",
  #                                                           "Dmrt1L", "SoxH", "FoxL2",
  #                                                           "Wnt8a", "FoxB2"))
  
  plot <- ggplot(data = expr_to_plot, aes(x = as.factor(Time), y = value, color = name, linetype = name)) +
    stat_summary(aes(group = gene), fun = mean, geom = "line", linewidth = 0.8) +
    geom_point(size = 2, alpha = 0.5) +
    xlab("Hours post fertilization (hpf)") +
    ylab("Normalized expression value")
  
  return(plot)
}


############################
#     THEMES FOR PLOTS     #
############################

theme_for_plots <- theme(
  # aspect.ratio = 1,
  plot.background = element_rect(fill = "transparent", colour = NA), 
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey90", lineend = "round"),
  panel.border = element_rect(colour = "black", linewidth = .6),
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.key = element_rect(fill = "transparent", colour = NA),
  legend.key.width = unit(.4, "cm"),
  legend.key.height = unit(.4, "cm"),
  legend.position = "right",
  legend.title = element_text(face = "bold"),
  # plot.title = element_text(size = 13, hjust = 0.0, vjust = 1.75, face = "bold"),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = "black", linewidth = .4),
  axis.ticks.length = unit(0.10, "cm"),
  axis.text.x = element_text(color = "black",
                             margin = margin(t = 4, r = 0, b = 0, l = 0)),
  axis.text.y = element_text(color = "black",
                             margin = margin(t = 0, r = 4, b = 0, l = 0)),
  axis.title.y = element_text(angle = 90, size = 13,
                              margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(angle = 0, size = 13,
                              margin = margin(t = 10, r = 0, b = 0, l = 0)),
  strip.text = element_text(color = "black", face = "bold", hjust = 0),
  strip.placement = "outside"
  # strip.background = element_rect(color = "black", linewidth = .6, linetype = "solid"),
  # strip.text.x = element_text(color = "black")
)
  

########################################
#     RUN READ COUNT NORMALIZATION     #
########################################

write("", stdout())
write("--- Analyzing gene raw counts from Bowtie output---", stdout())

normalized_reads_bowtie <- normalizeReads_plotPCA("03b_mapped_reads_BOWTIE/ALL.rawmapping.stats.csv", "00_input/SRA_metadata.tsv")

write("", stdout())
write("--- Analyzing gene raw counts from STAR output---", stdout())

normalized_reads_star <- normalizeReads_plotPCA("03a_mapped_reads_STAR/gene_count_matrix.csv", "00_input/SRA_metadata.tsv")


########################
#     RUN MASIGPRO     #
########################

write("", stdout())
write("--- maSigPro with vstnorm all timepoints ---", stdout())

# run maSigPro on vst normalized gene counts from all timepoints
vstnorm_alltimepoints <- timeseries_maSigPro("04_PCA_readcounts/gene_counts_vstTransformed_norm.tsv", "00_input/edesign_object.tsv")

write("", stdout())
write("--- maSigPro with vstnorm 24h ---", stdout())

# run maSigPro on vst normalized gene counts from timepoints up to 24hpf
vstnorm_24hpf <- timeseries_maSigPro("04_PCA_readcounts/gene_counts_vstTransformed_norm.tsv", "00_input/edesign_object_24hpf.tsv")

write("", stdout())
write("--- maSigPro with rlognorm all timepoints ---", stdout())

write("", stdout())
write("--- Saving RData ---", stdout())

save.image(file = "05_masigpro_analysis/DE_timeseries.Rdata")

# run maSigPro on vst normalized gene counts from all timepoints
rlog_alltimepoints <- timeseries_maSigPro("04_PCA_readcounts/gene_counts_rlogTransformed_norm.tsv", "00_input/edesign_object.tsv")

write("", stdout())
write("--- maSigPro with rlognorm 24h ---", stdout())

# run maSigPro on vst normalized gene counts from timepoints up to 24hpf
rlog_24hpf <- timeseries_maSigPro("04_PCA_readcounts/gene_counts_rlogTransformed_norm.tsv", "00_input/edesign_object_24hpf.tsv")

write("", stdout())
write("--- Saving RData ---", stdout())

save.image(file = "05_masigpro_analysis/DE_timeseries.Rdata")


####################################
#     CHECK DE OF TARGET GENES     #
####################################

# # vst norm, all timepoints
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B017427"] # vasa:    no
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B086491"] # nanos:    no
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B073040"] # piwia:    no
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B064020"] # piwib:    no
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B067558"] # spPHI:    no
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B093608"] # dmrt1L:  no
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B014180"] # soxH:    no
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B094018"] # foxL2:   YES, cluster 3
# vstnorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B085403"] # wnt8a:   YES, cluster 9

# # vst norm, all timepoints
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI03912.1"] # vasa:    YES, cluster 3
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI80212.1"] # nanos:    YES, cluster 3
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI60324.1"] # piwia:    no
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI83778.1"] # piwib:    YES, cluster 2
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI56617.1"] # spPHI:    YES, cluster 2
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI03798.1"] # dmrt1L:  no
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI30824.1"] # soxH:    no
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI49864.1"] # foxL2:   YES, cluster 4
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI54402.1"] # wnt8a:   no
# vstnorm_alltimepoints$visualization_genes$cut["Mgal_VDI35942.1"] # foxB2:   no

# # vst norm, 24hpf
# vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B017427"] # vasa:    YES, cluster 6
# vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B093608"] # dmrt1L:  no
# vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B014180"] # soxH:    no
# vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B094018"] # foxL2:   no
# vstnorm_24hpf$visualization_genes$cut["gene-MGAL_10B085403"] # wnt8a:   YES, cluster 2

# # rlog norm, all timepoints
# rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B017427"] # vasa:   no
# rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B093608"] # dmrt1L: no
# rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B014180"] # soxH:   no
# rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B094018"] # foxL2:  YES, cluster 6
# rlognorm_alltimepoints$visualization_genes$cut["gene-MGAL_10B085403"] # wnt8a:  YES, cluster 1

# # rlog norm, 24hpf
# rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B017427"] # vasa:   YES, cluster 7
# rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B093608"] # dmrt1L: no
# rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B014180"] # soxH:   no
# rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B094018"] # foxL2:  no
# rlognorm_24hpf$visualization_genes$cut["gene-MGAL_10B085403"] # wnt8a:  YES, cluster 4


################################
#     PLOT GENE EXPRESSION     #
################################

load(file = "05_masigpro_analysis/DE_timeseries.Rdata")

# generate gene conversion dataframe
gene_names <- data.frame(gene = c("gene-MGAL_10B017427", "gene-MGAL_10B086491", "gene-MGAL_10B073040", "gene-MGAL_10B064020", "gene-MGAL_10B067558", "gene-MGAL_10B093608","gene-MGAL_10B014180","gene-MGAL_10B094018","gene-MGAL_10B085403","gene-MGAL_10B093191"),
                         stringtie = c("MSTRG.9718", NA, NA, NA, NA, "MSTRG.9631", "MSTRG.26773", "MSTRG.38958", "MSTRG.41862", "MSTRG.30071"),
                         cds = c("Mgal_VDI03912.1", "Mgal_VDI80212.1", "Mgal_VDI60324.1", "Mgal_VDI83778.1", "Mgal_VDI56617.1", "Mgal_VDI03798.1", "Mgal_VDI30824.1", "Mgal_VDI49864.1", "Mgal_VDI54402.1", "Mgal_VDI35942.1"),
                         name = c("Vasa", "Nanos", "Piwia", "Piwib", "spPHI", "Dmrt1L", "SoxH", "FoxL2", "Wnt8a", "FoxB2"))


vst.expr.plot <- plot_expression_values(vstnorm_alltimepoints$counts,
                                        c("Mgal_VDI03912.1", "Mgal_VDI80212.1", "Mgal_VDI56617.1"),
                                        # gene_names$cds,
                                        gene_names,
                                        vstnorm_alltimepoints$edesign) +
  # scale_colour_manual(values = c("Red", "Cyan", "Magenta", "Yellow", "Grey", "Grey")) +
  # scale_linetype_manual(values = c(rep("solid", 4), rep("dashed", 2))) +
  ggtitle("Vst transformed gene expression") +
  labs(color = "Gene", linetype = "Gene") +
  theme_bw(base_size = 12) +
  theme_for_plots

vst.expr.plot

# ggsave("03b_mapped_reads_BOWTIE/vstgene_expression.pdf", plot = vst.expr.plot, device = "pdf",
       # dpi = 300, height = 4.7, width = 6, units = ("in"), bg = 'transparent')

# ggsave("03b_mapped_reads_BOWTIE/vstgene_expression.png", plot = vst.expr.plot, device = "png",
       # dpi = 300, height = 4.7, width = 6, units = ("in"), bg = 'transparent')



# rlog.expr.plot <- plot_expression_values(rlognorm_alltimepoints$counts,
                                        # gene_names$cds,
                                        # gene_names,
                                        # rlognorm_alltimepoints$edesign) +
  # ggtitle("Rlog transformed gene expression")


# panel <- ggpubr::ggarrange(rlog.expr.plot, vst.expr.plot,
                           # ncol = 2, nrow = 1,
                           # common.legend = TRUE, legend = "right")
# panel


##################################
#     PLOT MASIGPRO CLUSTERS     #
##################################

# labels for plot strips
labeller_cluster <- vstnorm_alltimepoints$visualization_genes$cut %>%
  tibble(group = .) %>%
  group_by(group) %>%
  summarise(n_genes = n()) %>%
  ungroup() %>%
  mutate(label = paste0("Cluster ", group,
                        " (", format(n_genes, big.mark = ","), " genes)")) %>%
  select(-n_genes) %>%
  deframe() %>%
  as_labeller()

# create dataframe to be plotted
data_to_plot <- vstnorm_alltimepoints$visualization_genes$cut %>%
  tibble(gene = names(.),
         group = .) %>%
  left_join(vstnorm_alltimepoints$counts %>%
              rownames_to_column(var = "gene")) %>%
  pivot_longer(-c(gene, group), names_to = "sample", values_to = "profile") %>%
  left_join(vstnorm_alltimepoints$sig_genes_list$sig.genes$edesign %>%
              rownames_to_column(var = "sample") %>%
              select(sample, Time))

# calculate mean expression values for each gene (grey lines)
gene_mean <- data_to_plot %>%
  group_by(gene, group, Time) %>%
  summarise(profile_mean = mean(profile), .groups = "drop") %>%
  mutate(line_type = case_when(gene == "Mgal_VDI03912.1" ~ "Vasa",
                               gene == "Mgal_VDI80212.1" ~ "Nanos",
                               gene == "Mgal_VDI83778.1" ~ "Piwi-b",
                               gene == "Mgal_VDI56617.1" ~ "spPHI",
                               gene == "Mgal_VDI49864.1" ~ "Fox-L2",
                               TRUE ~ "Other genes")) %>%
  mutate(line_type = factor(line_type,
                            levels = c("Vasa", "Nanos", "Piwi-b", "spPHI", "Fox-L2", "Other genes", "Mean expression")))

# calculate mean expression values for each cluster (dashed lines)
group_mean <- data_to_plot %>%
  group_by(group, Time) %>%
  summarise(mean_profile = mean(profile),
            sd_profile   = sd(profile),
            .groups = "drop") %>%
  mutate(ymin = mean_profile - sd_profile,
         ymax = mean_profile + sd_profile,
         line_type = "Mean expression")

# create plot
masigpro_groups_panel <- ggplot() +
  
  # grey gene lines
  geom_line(data = gene_mean,
            aes(x = Time, y = profile_mean, group = gene, color = line_type),
            alpha = 0.3, linewidth = 0.3, lineend = "round") +
  
  # # ribbon around group mean
  # geom_ribbon(data = group_mean,
  #             aes(x = Time, ymin = ymin, ymax = ymax, group = group),
  #             fill = "azure4", alpha = 0.2) +
  
  # mean line
  geom_line(data = group_mean,
            aes(x = Time, y = mean_profile, group = group, color = line_type),
            linewidth = 0.3, linetype = "dashed", lineend = "round") +
  
  # vasa line
  geom_line(data = gene_mean %>%
              filter(gene == "Mgal_VDI03912.1"),
            aes(x = Time, y = profile_mean, group = gene, color = line_type),
            linewidth = 0.5, lineend = "round") +
  
  # nanos line
  geom_line(data = gene_mean %>%
              filter(gene == "Mgal_VDI80212.1"),
            aes(x = Time, y = profile_mean, group = gene, color = line_type),
            linewidth = 0.5, lineend = "round") +
  
  # piwib line
  geom_line(data = gene_mean %>%
              filter(gene == "Mgal_VDI83778.1"),
            aes(x = Time, y = profile_mean, group = gene, color = line_type),
            linewidth = 0.5, lineend = "round") +
  
  # spPHI line
  geom_line(data = gene_mean %>%
              filter(gene == "Mgal_VDI56617.1"),
            aes(x = Time, y = profile_mean, group = gene, color = line_type),
            linewidth = 0.5, lineend = "round") +

  # FoxL2 line
  geom_line(data = gene_mean %>%
              filter(gene == "Mgal_VDI49864.1"),
            aes(x = Time, y = profile_mean, group = gene, color = line_type),
            linewidth = 0.5, lineend = "round") +
  
  scale_color_manual(name = "Gene",
                     values = c("Other genes" = "azure3",
                                "Vasa" = "firebrick1",
                                "Nanos" = "chocolate1",
                                "Piwi-b" = "darkorchid1",
                                "spPHI" = "deeppink1",
                                "Fox-L2" = "gold1",
                                "Mean expression" = "azure4")) +
  
  scale_x_continuous(breaks = unique(gene_mean$Time), expand = c(0, 0)) +
  
  facet_wrap(~group, labeller = labeller_cluster) +
  
  xlab("Hours post fertilization (hpf)") +
  ylab("Normalized expression value") +
  
  theme_minimal(base_size = 8) +
  theme_for_plots +
  theme(panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(0.4, "lines"))

masigpro_groups_panel

# save panel
ggsave("05_masigpro_analysis/vstnorm_alltimepoints_maSigPro_clusters.pdf",
       plot = masigpro_groups_panel, device = "pdf",
       dpi = 300, height = 6, width = 9, units = ("in"), bg = "white")

ggsave("05_masigpro_analysis/vstnorm_alltimepoints_maSigPro_clusters.png",
       plot = masigpro_groups_panel, device = "png",
       dpi = 300, height = 6, width = 10, units = ("in"), bg = "white")
  
save.image(file = "05_masigpro_analysis/DE_timeseries.Rdata")
load("05_masigpro_analysis/DE_timeseries.Rdata")


#######################################
#     SAVE CLUSTER GENES TO FILES     #
#######################################

# get the list of genes per cluster
gene_lists <- data_to_plot %>%
  select(group, gene) %>%
  distinct(group, gene) %>%
  group_by(group) %>%
  summarise(genes = list(gene), .groups = "drop") %>%
  deframe()
  { setNames(.$group, .$genes) }

# write the list of genes per cluster to separate files
for (cluster in names(gene_lists)) {
  file_path <- paste0("05_masigpro_analysis/01_genes_per_cluster/cluster_", cluster, "_genes.ls")
  writeLines(gene_lists[[cluster]],
             file_path)
}

vstnorm_alltimepoints$counts %>%
  rownames_to_column(var = "genes") %>%
  pull(genes) %>%
  writeLines("05_masigpro_analysis/01_genes_per_cluster/gene_universe.ls")