#!/usr/bin/env Rscript


library(ggtree)
library(tidyverse)


####################
#     FUNCTION     #
####################

# function to load, plot, and save tree
plot.tree <- function(input_filename) {
  
  tree <- read.tree(file = input_filename) %>%
    ape::drop.tip("Mgal_VDI03911.1")
  
  tree$tip.label <- ifelse(grepl("Drer|Mmus|Hsap|Cele|Dmel", tree$tip.label),
                           paste0("*", tree$tip.label),
                           tree$tip.label)
  
  tree.plot <- ggtree(tree, size = 0.4) +
    
    geom_tiplab(align = TRUE, size = 1.5, colour = "grey80", linesize = 0.3) +
    geom_tiplab(align = TRUE, size = 1.5, linetype = NULL) +
    geom_point(aes(color = as.numeric(label)), size = 1) +
  
    coord_cartesian(clip = 'off') +
    
    scale_color_gradient(name = "Bs values",
                          low = "firebrick3", high = "forestgreen",
                          limits = c(0,100), breaks = seq(0,100,25), na.value = "transparent") +
    
    ggtitle(expression(paste("A) ML phylogenetic tree of ", italic("Vasa/Ddx4"), " and ", italic("Ddx3"), "genes"))) +
    
    theme(legend.position = c(0.1,0.9),
          legend.text = element_text(size = 6),
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.3, "cm"),
          legend.background = element_rect(fill = "white", color = "grey60", linewidth = 0.2),
          legend.title.position = "top",
          legend.title = element_text(size = 8, face = "bold", hjust = 0.5),
          legend.ticks.length = unit(0.05, "cm"),
          title = element_text(size = 8))
  
  return(tree.plot %>% flip(89,132))
}


#####################
#     PLOT TREES    #
#####################

vasa.tree <- plot.tree("06_vasa_ML/vasa_bivalves_withRefOut_reduced_aligned_trim04.faa_rooted.treefile")
  
vasa.tree


########################################
#     PLOT ALIGNMENT DEAD/DEAH BOX     #
########################################

alignment.domain <- Biostrings::readAAStringSet("06_vasa_ML/vasa_bivalves_withRefOut_reduced_aligned_trim04_DEADDEAHdomain_treeOrder.faa",
                                                format = "fasta")

alignment.domain <- alignment.domain[names(alignment.domain) != "Mgal_VDI03911.1"]

names(alignment.domain) <- ifelse(grepl("Drer|Mmus|Hsap|Cele|Dmel", names(alignment.domain)),
                           paste0("*", names(alignment.domain)),
                           names(alignment.domain))

msa.domain <- ggmsa::ggmsa(alignment.domain, seq_name = TRUE, border = FALSE,
                           # color = "Clustal",
                           font = NULL) +
  
  geom_segment(aes(x = 270, xend = 270, y = 87, yend = 45),
               color = "grey60", linewidth = 0.5, lineend = "round") +
  
  geom_text(aes(x = 276, y = mean(c(87,45))), label = "Vasa/Ddx4",
            color = "grey30", size = 3.5, fontface = "bold.italic", angle = -90) +
  
  geom_segment(aes(x = 270, xend = 270, y = 44, yend = 1),
  color = "grey60", linewidth = 0.5, lineend = "round") +
  
  geom_text(aes(x = 276, y = mean(c(44,1))), label = "Ddx3",
            color = "grey30", size = 3.5, fontface = "bold.italic", angle = -90) +
  
  annotate("rect", xmin = 0.5, xmax = lengths(alignment.domain[1])+0.5, ymin = 64.5, ymax = 65.5,
           linewidth = 0.5, color = "grey30", alpha = 0) +
  
  geom_text(aes(x = 266, y = 64.7),
            color = "grey30", size = 2.3, fontface = "bold", label = "**") +
  
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,lengths(alignment.domain[1]), 25)) +
  
  ggtitle("B) Amino acid alignment of the DEAD/DEAH domain") +
  
  coord_cartesian(clip = 'off') +
  
  theme_minimal() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.ticks.x = element_line(linewidth = 0.3, color = "grey30"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8, color = "grey30")
        )

msa.domain


###################################
#     PLOT ALIGNMENT PROTEINS     #
###################################

# alignment.protein <- phangorn::read.aa("vasa_MgalDrerHsapDmel_aligned.faa", format = "fasta")
# 
# names(alignment.protein) <- ifelse(grepl("Drer|Mmus|Hsap|Cele|Dmel", names(alignment.protein)),
#                                   paste0("*", names(alignment.protein)),
#                                   names(alignment.protein))
# 
# 
# plot.protein.alignment <- function(alignment, start, stop) {
#   
#   ggmsa::ggmsa(alignment, seq_name = TRUE, border = FALSE,
#                # color = "Clustal",
#                # font = NULL,
#                start = start, end = stop,
#                char_width = 0.5
#                # position_highlight = TRUE
#                # consensus_views = TRUE,
#                # ref = "*Drer_NP.571132.1_ddx4"
#                # use_dot = TRUE
#                ) +
#     
#     coord_cartesian(clip = 'off') +
#     
#     theme_minimal() +
#     theme(panel.grid = element_blank(),
#           # axis.line.x = element_line(linewidth = 0.1),
#           axis.ticks.x = element_line(linewidth = 0.1, color = "grey30"),
#           axis.text.y = element_text(size = 4.2),
#           axis.text.x = element_text(size = 5.3, color = "grey30")
#           # axis.line = element_line(color = "black"),
#           # axis.title.x = element_text(hjust = -10)
#     )
# }
# 
# msa.title <- ggplot() +
#   ggtitle(expression("C) Protein alignment of Vasa from"~italic("Mytilus galloprovincialis")~"and reference species")) +
#   coord_cartesian(clip = 'off') +
#   theme_minimal() +
#   theme(title = element_text(size = 8))
# 
# msa.protein1 <- plot.protein.alignment(alignment.protein, 1, 206)
# msa.protein2 <- plot.protein.alignment(alignment.protein, 207, 413) +
#   annotate("rect", xmin = 346.5, xmax = 413.5, ymin = 0.5, ymax = 5.5,
#            linewidth = 0.6, color = "grey30", alpha = 0)
# msa.protein3 <- plot.protein.alignment(alignment.protein, 414, 620) +
#   annotate("rect", xmin = 413.5, xmax = 601.5, ymin = 0.5, ymax = 5.5,
#            linewidth = 0.6, color = "grey30", alpha = 0) +
#   annotate("rect", xmin = 605.5, xmax = 620.5, ymin = 0.5, ymax = 5.5,
#            linewidth = 0.6, color = "grey30", alpha = 0)
# msa.protein4 <- plot.protein.alignment(alignment.protein, 621, 828) +
#   annotate("rect", xmin = 620.5, xmax = 735.5, ymin = 0.5, ymax = 5.5,
#            linewidth = 0.6, color = "grey30", alpha = 0)
# 
# alignment.panel <- ggpubr::ggarrange(msa.title, msa.protein1, msa.protein2, msa.protein3, msa.protein4, nrow = 5)


######################
#     PLOT PANEL     #
######################

panel <- ggpubr::ggarrange(vasa.tree, NULL, msa.domain, align = "hv",
                           ncol = 3,
                           widths = c(0.6, 0.045, 0.4))
                           

panel

ggsave("07_figures_for_ms/Supp_Fig_S1.pdf",
       plot = panel, device = "pdf",
       dpi = 500, height = 7, width = 12, units = ("in"), bg = 'white',
       limitsize = FALSE)
ggsave("07_figures_for_ms/Supp_Fig_S1.png",
       plot = panel, device = "png",
       dpi = 500, height = 7, width = 12, units = ("in"), bg = 'white',
       limitsize = FALSE)
