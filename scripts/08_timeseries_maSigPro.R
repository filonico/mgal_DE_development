#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly=TRUE))
  BiocManager::install("maSigPro")


timeseries_maSigPro <- function(file_norm_counts,edesign) {
  
  # load normalized read counts (output from 07_PCA_rawcounts.R)
  normalized.counts <- read.table(file_norm_counts, header = TRUE, row.names = 1, sep ="\t")
  
  # load edesign object
  edesign_object <- read.table(edesign, header = TRUE, sep = "\t", row.names = 1)
  
  # define regression model
  design <- maSigPro::make.design.matrix(edesign_object)
  
  # find significant genes
  fit <- maSigPro::p.vector(normalized.counts, design, counts = TRUE)
  
  # find significant differences
  tstep <- maSigPro::T.fit(fit)
  
  # obtain lists of significant genes
  get <- maSigPro::get.siggenes(tstep, vars = "all")
  
  # visualization
  seegenes <- maSigPro::see.genes(get$sig.genes, k.mclust = TRUE)
  
  return(list(sig_genes_list = get, visualization_genes = seegenes))
}

vstnorm_alltimepoints <- timeseries_maSigPro("03_mapped_reads/gene_count_vstTransformed.tsv", "00_input/edesign_object.tsv")
vstnorm_24hpf <- timeseries_maSigPro("03_mapped_reads/gene_count_vstTransformed.tsv", "00_input/edesign_object_24hpf.tsv")
