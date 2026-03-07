#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  BiocManager::install("maSigPro")
}

library(tidyverse)
library(GOSemSim)
library(smacof)

setwd("/data/evassvis/fn76/mgal_DE_development")


#####################
#     FUNCTIONS     #
#####################

plot_module_expression_values <- function(data_to_plot, cluster) {
  
  # create plot
  masigpro_expression <- ggplot() +
    
    # grey gene lines
    geom_line(data = data_to_plot %>% filter(group %in% cluster),
              aes(x = Time, y = profile_mean, group = gene, color = line_type),
              alpha = 0.3, linewidth = 0.3, lineend = "round") +
    
    # mean line
    # geom_line(data = group_mean %>% filter(group %in% cluster),
    #           aes(x = Time, y = mean_profile, group = group, color = line_type),
    #           linewidth = 0.3, linetype = "dashed", lineend = "round") +
    
    # vasa line
    geom_line(data = data_to_plot %>% filter(group %in% cluster) %>%
                filter(gene == "Mgal_VDI03912.1"),
              aes(x = Time, y = profile_mean, group = gene, color = line_type),
              linewidth = 1, lineend = "round") +
    
    # nanos line
    geom_line(data = data_to_plot %>% filter(group %in% cluster) %>%
                filter(gene == "Mgal_VDI80212.1"),
              aes(x = Time, y = profile_mean, group = gene, color = line_type),
              linewidth = 0.5, linetype = "dashed", lineend = "round") +
    
    # piwib line
    geom_line(data = data_to_plot %>% filter(group %in% cluster) %>%
                filter(gene == "Mgal_VDI83778.1"),
              aes(x = Time, y = profile_mean, group = gene, color = line_type),
              linewidth = 0.5, linetype = "dashed", lineend = "round") +
    
    # spPHI line
    geom_line(data = data_to_plot %>% filter(group %in% cluster) %>%
                filter(gene == "Mgal_VDI56617.1"),
              aes(x = Time, y = profile_mean, group = gene, color = line_type),
              linewidth = 0.5, linetype = "dashed", lineend = "round") +
    
    # FoxL2 line
    geom_line(data = data_to_plot %>% filter(group %in% cluster) %>%
                filter(gene == "Mgal_VDI49864.1"),
              aes(x = Time, y = profile_mean, group = gene, color = line_type),
              linewidth = 0.5, lineend = "round") +
    
    scale_color_manual(name = "Genes",
                       values = c("Other genes" = "azure3",
                                  "vasa" = "firebrick1",
                                  "nanos" = "chocolate1",
                                  "piwi-b" = "darkorchid1",
                                  "spPHI" = "deeppink1",
                                  "FoxL2" = "gold1",
                                  "Mean expression" = "azure4")) +
    
    # scale_x_continuous(breaks = unique(data_to_plot$Time), expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0), labels = c("0" = "Oocyte")) +
    
    xlab("Hours post fertilization (hpf)") +
    ylab("Normalized expression value") +
    
    theme_minimal(base_size = 12) +
    theme_for_plots +
    theme(panel.spacing.x = unit(1, "lines"),
          panel.spacing.y = unit(0.4, "lines"))
  
  return(masigpro_expression)
}

get_GOterm_semantic_mds <- function(goterm_enrich_file) {
  goerms_bp_elim <- read.table(goterm_enrich_file, sep = "\t", header = TRUE, quote = "") %>%
    filter(classicFisher < 0.05)
  
  annoDb_sycon <- read.table("00_input/mgal_proteome_GOannotation.tsv", sep = "\t", skip = 4) %>%
    select(V2, V5, V9) %>%
    mutate(V9 = case_when(V9 == "C" ~ "CC",
                          V9 == "P" ~ "BP",
                          V9 == "F" ~ "MF",
                          TRUE ~ "")) %>%
    rename(c("V2" = "GENE", "V5" = "GO", "V9" = "ONTOLOGY"))
  annoDb_sycon
  
  semantic_data <- godata(annoDb = annoDb_sycon, ont = "BP", computeIC = FALSE)
  
  semantic_similarity <- mgoSim(goerms_bp_elim$GO.ID, goerms_bp_elim$GO.ID,
                                semData = semantic_data, measure = "Wang", combine = NULL)
  
  semantic_dissimilarity <- sim2diss(semantic_similarity, method = "confusion", to.dist = TRUE)
  
  mds <- mds(semantic_dissimilarity)
  
  tibble_to_plot <- mds$conf %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "GO") %>%
    left_join(goerms_bp_elim, by = join_by(GO == GO.ID)) %>%
    arrange(desc(classicFisher))
  
  return(tibble_to_plot)
}


############################
#     THEMES FOR PLOTS     #
############################

theme_for_MDS <- theme(
  # aspect.ratio = 1,
  plot.background = element_blank(), 
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 10, face = "bold"),
  plot.title = element_text(size = 13, hjust = 0.5, vjust = 1.75, face = "bold"),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank()
)

md_arrows <- list(
  
  annotation_custom(grob = grid::segmentsGrob(
    x0 = unit(0, "mm"), x1 = unit(12, "mm"),
    y0 = unit(0, "mm"), y1 = unit(0, "mm"),
    arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
    gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  
  annotation_custom(grob = grid::segmentsGrob(
    x0 = unit(0, "mm"), x1 = unit(0, "mm"),
    y0 = unit(0, "mm"), y1 = unit(12, "mm"),
    arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
    gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  
  annotation_custom(grob = grid::textGrob(label = "MD1",
                                          x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))),
  
  annotation_custom(grob = grid::textGrob(label = "MD2",
                                          x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
  
  coord_cartesian(clip = "off"))


theme_for_plots <- theme(
  # aspect.ratio = 1,
  plot.background = element_rect(fill = "transparent", colour = NA), 
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey90", lineend = "round"),
  panel.border = element_rect(colour = "black", linewidth = .6),
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.key = element_rect(fill = "transparent", colour = NA),
  # legend.key.width = unit(.4, "cm"),
  # legend.key.height = unit(.4, "cm"),
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


################################
#     PLOT GENE EXPRESSION     #
################################

load("05_masigpro_analysis/DE_timeseries.Rdata")

cluster2_3_expression <- plot_module_expression_values(gene_mean %>% mutate(Time = as_factor(Time)), c("2", "3"))
cluster2_3_expression


#########################
#     PLOT GO TERMS     #
#########################

cluster2_mds <- get_GOterm_semantic_mds("05_masigpro_analysis/01_genes_per_cluster/cluster_2_genes_GOterms_topGO_BP_elim.txt")

cluster2_mds_plot <- cluster2_mds %>%
  ggplot(aes(D1, D2)) +
  geom_point(aes(size = Significant, fill = -log(classicFisher)),
             shape = 21, color = "#08306B", alpha = 0.7) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0000398", "GO:0006397", "GO:0005944")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "mRNA\nsplicing", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0031324", "GO:0010558", "GO:0009890")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Negative regulation\nof biosynthetis", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0003047", "GO:0001301", "GO:0071103")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Mitosis and\nchromatin remodelling", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0006403", "GO:0000165", "GO:0032012")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Cell proliferation", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +
  
  # ggrepel::geom_text_repel(data = . %>%
  #                            filter(Significant >= 5),
  #                          aes(label = GO)) +
   
  labs(#title = "Enriched GOs",
    x = "MD1", y = "MD2",
    fill = "−log(p-val)", size = "Sig. genes") +
  
  scale_size_continuous(range = c(1,14), breaks = c(30, 60, 90)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue",
                      breaks = c(5, 15, 25)) +
  
  guides(size = guide_legend(override.aes = list(color = "#706db8")),
         color = guide_colorbar(barheight = 4, barwidth = 1)) +
  
  coord_cartesian(clip = "off") +
  
  md_arrows +
  
  theme_bw(base_size = 12) +
  theme_for_MDS +
  theme(plot.margin = margin(0, 0, 6, 6, "mm"))
cluster2_mds_plot

cluster3_mds <- get_GOterm_semantic_mds("05_masigpro_analysis/01_genes_per_cluster/cluster_3_genes_GOterms_topGO_BP_elim.txt")

cluster3_mds_plot <- cluster3_mds %>%
  ggplot(aes(D1, D2)) +
  geom_point(aes(size = Significant, fill = -log(classicFisher)),
             shape = 21, color = "#08306B", alpha = 0.7) +
  
  geom_label(data = . %>%
               filter(GO %in% c("GO:0006325", "GO:0006338")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Chromatin\nremodelling", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0002402", "GO:0000278")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Mitosis\nand cell cyle", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0031327", "GO:0009890")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Negative regulation\nof biosynthetis", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0006260", "GO:0006261")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "DNA replication", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0000398", "GO:0006397", "GO:0006567", "GO:0006511")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "mRNA splicing\nprotein modifications", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +
  
  # ggrepel::geom_text_repel(data = . %>%
  #                            filter(Significant >= 50),
  #                          aes(label = GO)) +
  
  labs(#title = "Enriched GOs",
    x = "MD1", y = "MD2",
    fill = "−log(p-val)", size = "Sig. genes") +
  
  scale_size_continuous(range = c(1,14), breaks = c(30, 60, 90)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue",
                      breaks = c(5, 15, 25)) +
  
  guides(size = guide_legend(override.aes = list(color = "#706db8")),
         color = guide_colorbar(barheight = 4, barwidth = 1)) +
  
  coord_cartesian(clip = "off") +
  
  md_arrows +
  
  theme_bw(base_size = 12) +
  theme_for_MDS +
  theme(plot.margin = margin(0, 0, 6, 6, "mm"))
cluster3_mds_plot


#################
#     PANEL     #
#################

panel <- ggpubr::ggarrange(cluster2_3_expression + theme(axis.title.x = element_blank()),
                           cluster3_mds_plot,
                           widths = c(1,0.6), labels = "AUTO")
panel


ggsave("06_figures_for_ms/fig1_panel.pdf",
       plot = panel, device = "pdf",
       dpi = 300, height = 6/1.5, width = 20/1.5, units = ("in"), bg = "white")
ggsave("06_figures_for_ms/fig1_panel.png",
       plot = panel, device = "png",
       dpi = 300, height = 6/1.5, width = 20/1.5, units = ("in"), bg = "white")
