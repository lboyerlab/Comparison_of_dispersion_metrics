library(dplyr)
library(Seurat)
library(patchwork)
library(DescTools)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggvenn)
library(patchwork)
library(tidyr)
library(ggrepel)
library(plyr)
library(ggpubr)

setwd("~")

## load in Seurat object from preprocessing & calculate noise ##
cm_only <- readRDS("T21/cm_only.rds")
conditions <- c("WT_A1", "WT_A7", "Dp1Tyb")
by_condition <- SplitObject(cm_only, split.by = "orig.ident")
Idents(cm_only) <- cm_only@meta.data$orig.ident
cm_only <- subset(cm_only, ident = c("WT_A1", "WT_A7", "Dp1Tyb"))

# merge the same conditions together, which are the control and Dp1TybDyrk1a
wt_merge <- JoinLayers(merge(by_condition[[1]], by_condition[[2]]), assay = "RNA")
merged_conds <- list(wt_merge, by_condition[[3]])

# Calculate noise for each condition
vmr <- function(x) var(x)/mean(x)

#' Calculates dispersion in a Seurat object using Gini index, VMR, variance, and entropy
#' @param seur Seurat object
#' @returns Matrix of transcriptional variability values 
calc_noise <- function(seur) {
  norm_counts <- NormalizeData(seur@assays$RNA@layers$counts)
  noise <- apply(norm_counts, 1, function(x) c(Gini = Gini(x), VMR = vmr(x), Variance = var(x), Entropy = Entropy(table(x), base = exp(1)))) # want to apply across the gene (across rows)
  noise <- t(noise)
  rownames(noise) <- rownames(seur) 
  return(noise)
}

noise_by_cond <- lapply(merged_conds, calc_noise)

# save noise calculations
m_str_conds <- c("WT", "Dp1Tyb")
path_to_data <- "T21/" # TODO: change this for each dataset 
for (i in 1:length(noise_by_cond)) {
  saveRDS(noise_by_cond[[i]], file = paste(path_to_data, m_str_conds[[i]], "_noise.rds", sep = ""))
}
noise_dfs <- lapply(noise_by_cond, as.data.frame)

# load files if needed
dp1tyb_noise <- as.data.frame(readRDS(paste(path_to_data, "Dp1Tyb_noise.rds", sep = "")))
wt_noise <- as.data.frame(readRDS(paste(path_to_data, "WT_noise.rds", sep = "")))
noise_dfs <- list(wt_noise, dp1tyb_noise) # merged_conds <- list(wt_merge, by_condition[[3]], dp1tybdyrk1a)

## Filter out ribosomal, mitochondrial, and lowly expressed genes ##
only_ribo_mito_regex <- "^Mt|^mt|^Rps|^Rpl"

#' Filters our ribosomal and mitochondrial genes
#' @param df data.frame of transcriptional variability of genes
#' @returns filtered data.frame
filter_genes <- function(df) {
  filtered_df <- as.data.frame(df[!(grepl(only_ribo_mito_regex, rownames(df))), , drop = F])
  return(filtered_df)
}

noise_dfs_filtered <- noise_dfs
for (i in 1:length(noise_dfs_filtered)) {
  noise_df <- noise_dfs_filtered[[i]]
  noise_df$mean_exp <- rowMeans(merged_conds[[i]]@assays$SCT@data) 
  noise_df <- noise_df[noise_df$mean_exp > 0.05, ]
  noise_dfs_filtered[[i]] <- filter_genes(noise_df)
}

metrics <- head(colnames(noise_dfs_filtered[[1]]), -1) # Gini, VMR, Variance, Entropy

## Calculate change in noise between conditions ##
#' Calculates the absolute change in transcriptional variability of a gene between two conditions
#' @param metric String name of the dispersion metric of interest
#' @param df1 data.frame of transcriptional variability of genes in condition 1 
#' @param df2 data.frame of transcriptional variability of genes in condition 2 
#' @returns data.frame of genes ordered by absolute change in transcriptional variability 
changing_genes <- function(metric, df1, df2) {
  comb <- merge(df1[metric], df2[metric], by = 0)
  rownames(comb) <- comb$Row.names
  comb$abs_change <- abs(comb[, 2] - comb[, 3])
  comb <- comb[order(comb$abs_change, decreasing = T), ]
  return(comb)
}

#' Calculates the change in transcriptional variability at both time points
#' @param metric String name of the dispersion metric of interest
#' @param noise_func function to use to calculate change in transcriptional variability (ie, absolute change vs signed change) 
#' @returns data.frame of genes ordered by change in transcriptional variability 
calc_noise_change <- function(metric, noise_func) {
  dp1tyb <- noise_func(metric, noise_dfs_filtered[[1]], noise_dfs_filtered[[2]])
  return(dp1tyb)
}

#' Calculates the change in transcriptional variability of a gene between two conditions
#' @param metric String name of the dispersion metric of interest
#' @param df1 data.frame of transcriptional variability of genes in condition 1 
#' @param df2 data.frame of transcriptional variability of genes in condition 2 
#' @returns data.frame of genes ordered by change in transcriptional variability 
raw_noise_change <- function(metric, df1, df2) {
  comb <- merge(df1[metric], df2[metric], by = 0)
  rownames(comb) <- comb$Row.names
  comb$change <- comb[, 2] - comb[, 3]
  comb <- comb[order(comb$change, decreasing = T), ]
  return(comb)
}

noise_change_dfs <- lapply(metrics, calc_noise_change, noise_func = changing_genes) # list of lists, first layer is by metrics (Gini, VMR, Variance, Entropy)
raw_noise_change_dfs <- lapply(metrics, calc_noise_change, noise_func = raw_noise_change) # 1st index is metric, 2nd is condition
ranked_abs_noise_dfs <- lapply(raw_noise_change_dfs, function(x) x[order(abs(x$change), decreasing = T), ])

#### correlation between change in VMR and gene-level factors #### (S6 Fig)
library(biomaRt)
library(ggpubr)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

#' Calculates the average normalized expression of a Seurat object 
#' @param seu Seurat object
#' @returns data.frame of average normalized expression 
calc_norm_avg_exp <- function(seu) {
  norm <- NormalizeData(GetAssayData(seu, layer = "counts"))
  norm_df <- data.frame(avg_exp = rowMeans(norm))
  rownames(norm_df) <- rownames(norm)
  return(norm_df)
}

#' Calculates the average normalized expression across the control and matHG conditions at a given time point
#' @param control data.frame of average normalized expression in the control condition
#' @param matHG data.frame of average normalized expression in the matHG condition
#' @returns data.frame of average normalized expression across both condtions
make_avg_exp_df <- function(sample1, sample2) {
  avg_df <- data.frame(s1 = as.numeric(sample1[, 1]), s2 = as.numeric(sample2[, 1]))
  rownames(avg_df) <- rownames(sample1)
  avg_df$avg_exp <- rowMeans(avg_df)
  return(avg_df)
}

#' Plots the average normalized expression at a time point against the change in VMR between control and matHG
#' @param noise data.frame of change in VMR 
#' @param avg_exp data.frame of average normalized expression
#' @param x x-coordinate of the correlation coefficient
#' @param y y-coordinate of the correlation coefficient
#' @returns list of the data.frame with combined average expression and change in VMR and a ggplot object
plot_avg_exp_noise_change <- function(noise, avg_exp, x, y) {
  comb <- merge(noise, avg_exp, by = 0)
  comb <- subset(comb, select = -c(Row.names))
  print(colnames(comb))
  p <- ggplot(comb, aes(x = avg_exp, y = change)) + geom_point() + 
    labs(x = "Average Normalized Expression", y = "Absolute Change in VMR") + 
    geom_smooth(se = T) + 
    stat_cor(aes(label = after_stat(r.label)), method = "kendall", cor.coef.name = "tau", label.x = x, label.y = y, size = 8, color = "darkblue") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  return(list(comb, p))
}

# find the mean normalized expression across control and Dp1Tyb (T21) conditions
norm_avg_exp <- lapply(by_condition, calc_norm_avg_exp)
wt_avg_exp_norm <- make_avg_exp_df(norm_avg_exp[[1]], norm_avg_exp[[2]])
dp1tyb_wt_avg_exp_norm <- make_avg_exp_df(wt_avg_exp_norm, norm_avg_exp[[3]])

# test for correlation & plot
abs_change_VMR <- raw_noise_change_dfs[[2]] %>% mutate(change = abs(change))
dp1tyb_exp_abs_vmr_comb_all <- plot_avg_exp_noise_change(abs_change_VMR, dp1tyb_wt_avg_exp_norm, 2.75, 0.9) 
cor.test(abs(dp1tyb_exp_abs_vmr_comb_all[[1]]$change), dp1tyb_exp_abs_vmr_comb_all[[1]]$avg_exp, method = "kendall")

# Test correlation between absolute change in VMR and gene length, gene GC content
#' Annotates genes with gene length and gene GC content
#' @param noise data.frame of genes 
#' @returns data.frame of genes annotated with gene length and gene GC content
annotat_length_gc <- function(noise) {
  annotat <- getBM(
    attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "percentage_gene_gc_content"),
    filters = "external_gene_name",
    values = rownames(noise),
    mart = ensembl
  )
  comb <- merge(noise, annotat, by.x = 0, by.y = "external_gene_name")
  comb <- subset(comb, select = -c(Row.names)) %>% mutate(gene_length = end_position - start_position + 1)
  return(comb)
}

#' Plots absolute change in VMR against gene length, with the correlation coefficient
#' @param df data.frame of genes and their gene length and absolute change in VMR
#' @param x x-coordinate of the correlation coefficient
#' @param y y-coordinate of the correlation coefficient
#' @returns ggplot object
plot_gene_length <- function(df, x, y) {
  return(ggplot(df, aes(x = log(gene_length), y = abs(change))) + geom_point() + 
           labs(x = "log(Gene Length)", y = "Absolute Change in VMR") + 
           geom_smooth(se = T) + 
           stat_cor(aes(label = after_stat(r.label)), method = "kendall", cor.coef.name = "tau", label.x = x, label.y = y, size = 8, color = "darkblue") + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

#' Plots absolute change in VMR against gene GC content, with the correlation coefficient
#' @param df data.frame of genes and their gene GC content and absolute change in VMR
#' @param x x-coordinate of the correlation coefficient
#' @param y y-coordinate of the correlation coefficient
#' @returns ggplot object
plot_gc_content <- function(df, x, y) {
  return(ggplot(df, aes(x = percentage_gene_gc_content, y = abs(change))) + geom_point() + 
           labs(x = "Gene GC Content %", y = "Absolute Change in VMR") + 
           geom_smooth(se = T) + 
           stat_cor(aes(label = after_stat(r.label)), method = "kendall", cor.coef.name = "tau", label.x = x, label.y = y, size = 8, color = "darkblue") + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

dp1tyb_vmr_annot <- annotat_length_gc(raw_noise_change_dfs[[2]])

# test for correlation & generate plots
cor.test(abs(dp1tyb_vmr_annot$change), dp1tyb_vmr_annot$gene_length, method = "kendall")
cor.test(abs(dp1tyb_vmr_annot$change), dp1tyb_vmr_annot$percentage_gene_gc_content, method = "kendall")

dp1tyb_gene_length <- plot_gene_length(dp1tyb_vmr_annot, 10, 1)
dp1tyb_gc_content <- plot_gc_content(dp1tyb_vmr_annot, 52, 1)

# Test correlation between absolute change in VMR and gene phylostrata
library(stringr)

Mus_musculus.PhyloMap <- read.csv("GenEra_Mus_musculus.csv", header = T, row.names = 1) # from Barrera-Redondo et al., generated using GenEra
vmr_change_w_genEra_age <- merge(raw_noise_change_dfs[[2]], Mus_musculus.PhyloMap, by.x = "Row.names", by.y = "Gene")

#' Plots absolute change in VMR against gene phylostrata, with the correlation coefficient
#' @param df data.frame of genes and their gene phylostrata and absolute change in VMR
#' @param cc correlation coefficient
#' @param x x-coordinate of the correlation coefficient
#' @param y y-coordinate of the correlation coefficient
#' @returns ggplot object
plot_genEra_boxplot <- function(df, cc, x_cord, y_cord) {
  return(ggplot(df, aes(x = as.numeric(Phylostratum), y = abs(change), group = Phylostratum)) + geom_boxplot() + 
           labs(x = "Phylostratum", y = "Absolute Change in VMR") +
           # annotate("text", label = bquote(tau == .(cc) * ", p = " * .(p_val)), x = x_cord, y = y_cord, size = 8, color = "darkblue") + 
           annotate("text", label = bquote(tau == .(cc)), x = x_cord, y = y_cord, size = 8, color = "darkblue") + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

# test for correlation & plot
cor.test(abs(vmr_change_w_genEra_age$change), vmr_change_w_genEra_age$Phylostratum, method = "kendall")
dp1tyb_genEra <- plot_genEra_boxplot(vmr_change_w_genEra_age, 0.027, 14, 0.5)

# combine plots (S6 Fig)
(dp1tyb_exp_abs_vmr_comb_all[[2]] | dp1tyb_gene_length) /
  (dp1tyb_gc_content | dp1tyb_genEra)  & 
  theme(plot.margin = unit(c(0.8, 0.8, 0.8, 0.8), "cm"), axis.title = element_text(size = 14), axis.text = element_text(size = 12))

#### Compare against differential expression ####
#' Filters out DEGs that are not significant or that do not have a sufficiently large logFC 
#' @param metric String name of the dispersion metric of interest
#' @returns data.frame of genes ordered by change in transcriptional variability 
filter_degs <- function(df) {return(filter_genes(df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.25, ]))}

#' Find DEGs between 3 comparisons: Dp1Tyb vs WT, Dp1TybDyrk1a vs WT, and Dp1Tyb vs Dp1TybDyrk1a
#' @param str_test String name of the test to use for DE
#' @returns list of data.frames of DEGs after filtering
get_degs <- function(str_test) {
  Dp1Tyb <- FindMarkers(cm_only, ident.1 = c("WT_A1", "WT_A7"), ident.2 = "Dp1Tyb", verbose = T, recorrect_umi = FALSE)
  return(filter_degs(Dp1Tyb))
}


wilcox_degs <- get_degs("wilcox")

## Venn Diagram: DEGs & change in noise ##
#' Takes the top 100 DEGs as ranked by a given metric
#' @param df data.frame of DEGs
#' @param col_name name of the column / metric to order the DEGs by 
#' @param is_pos logical; T if the DEGs should be ordered in decreasing order, F otherwise
#' @param is_abs logical; T if we are interested in the absolute FC, F otherwise 
#' @returns data.frame of the top 100 DEGs as ordered by the specified metric
most_change <- function(df, col_name, is_pos, is_abs) {
  if (is_abs) {df$abs_FC <- abs(df$avg_log2FC); col_name <- "abs_FC"}
  return(head(df[order(df[[col_name]], decreasing = is_pos), ], 100))
}

# separate out top DEGs
wilcox_degs_abs_change <- most_change(wilcox_degs, "avg_log2FC", T, T)
wilcox_degs_p_val <- most_change(wilcox_degs, "p_val_adj", F, F)

# plot Venn Diagrams of top DEGs vs top genes with largest absolute change in VMR, Dp1Tyb vs WT 
vmr_v_degs <- ggvenn(list(HVGs = rownames(head(ranked_abs_noise_dfs[[2]], 100)), DEGs = rownames(wilcox_degs_abs_change)), set_name_size = 4, auto_scale = TRUE) + 
  ggtitle('T21 vs WT') + theme(plot.margin = margin(2, 2, 10, 2, "pt")) + coord_cartesian(clip = "off")
vmr_v_degs_pval <- ggvenn(list(HVGs = rownames(head(ranked_abs_noise_dfs[[2]], 100)), DEGs = rownames(wilcox_degs_p_val)), set_name_size = 4, auto_scale = TRUE) + 
  ggtitle('T21 vs WT') + theme(plot.margin = margin(2, 2, 10, 2, "pt")) + coord_cartesian(clip = "off") 

# explore what's at the intersection of DEGs & genes with largest absolute change in VMR
#' Takes a list of top DEGs and list of top genes by change in VMR and formats it so it can be printed nicely
#' @param degs list of data.frames of DEGs
#' @param changes list of data.frames of genes w/ largest absolute change in VMR 
#' @returns data.frame of the genes that are only in the list of DEGs, only in the list of genes with change in VMR, and the genes that are in both
inter_outersect <- function(degs, changes) {
  union <- intersect(rownames(degs), rownames(changes))
  chas <- setdiff(rownames(changes), rownames(degs))
  des <- setdiff(rownames(degs), rownames(changes))
  max_l <- max(sapply(list(union, chas, des), length))
  df <- data.frame(degs = c(des, rep("", max_l - length(des))), union = c(union, rep("", max_l - length(union))), change_vmr = c(chas, rep("", max_l - length(chas))))
  return(df)
}

p_val_degs_v_vmr_hvgs <- inter_outersect(wilcox_degs_p_val, head(ranked_abs_noise_dfs[[2]], 100))
abs_fc_degs_v_vmr_hvgs <- inter_outersect(wilcox_degs_abs_change, head(ranked_abs_noise_dfs[[2]], 100))

## GSEA ##
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)

#' Performs GSEA using KEGG on a given set of genes 
#' @param markers data.frames of genes with some metric they can be ranked by
#' @param metric the metric to rank the genes with for GSEA
#' @returns gseaResult object 
gsea_for_seurat <- function(markers, metric) {
  filtered <- filter_genes(markers)
  
  # convert to Entrez IDs & remove NAs
  filtered$entrez <- mapIds(org.Mm.eg.db, keys=rownames(filtered), keytype="SYMBOL", column="ENTREZID")
  filtered <- filtered[!is.na(filtered$entrez), ]
  
  # create ranked gene list 
  gene_list <- filtered[[metric]]
  names(gene_list) <- filtered$entrez
  gene_list <- sort(gene_list, decreasing=TRUE)
  print(head(gene_list))
  
  kegg <-  gseKEGG(geneList = gene_list, organism = "mmu")
  return(kegg)
}

#' Performs GSEA using KEGG once using DEGs and once using change in VMR
#' @param degs data.frame of DEGs
#' @param noise_df data.frame of genes and their change in transcriptional variability
#' @returns gseaResult object 
compare_FC_noise_GSEA <- function(degs, noise_df) {
  using_degs <- gsea_for_seurat(degs, "avg_log2FC")
  using_noise <- gsea_for_seurat(noise_df, "change")
  return(list(using_degs, using_noise)) 
}

# using 1000 genes each
dp1tyb_gsea_vmr_matched <- compare_FC_noise_GSEA(head(wilcox_degs, 1000), head(ranked_abs_noise_dfs[[2]], 1000))

dp1tyb_gsea_vmr_matched[[1]]@result[, c("Description", "setSize", "enrichmentScore", "pvalue")] 

# visualize GSEA results
dp1tyb_gsea_vmr_matched[[1]]@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", dp1tyb_gsea_vmr_matched[[1]]@result$Description, fixed = T)
dp1tyb_gsea_vmr_matched[[2]]@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", dp1tyb_gsea_vmr_matched[[2]]@result$Description, fixed = T)

dotplot(dp1tyb_gsea_vmr_matched[[1]], showCategory=5) + ggtitle("Using DEGs") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  dotplot(dp1tyb_gsea_vmr_matched[[2]], showCategory=5) + ggtitle("Using Change in VMR") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

## TF motif enrichment ##
library(RcisTarget)

motif_rankings <- importRankings("mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather") # from Aerts Lab (https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/)
data(motifAnnotations_mgi) # this creates object called motifAnnotations

#' Performs TF motif enrichment
#' @param df data.frames of the genes to analyze 
#' @returns list of a list of enriched TFs by time point, and a data.table with enriched motifs 
motif_enrichment <- function(df) {
  geneList <- list(dp1tyb = rownames(head(df, 100)))
  enrichment <- cisTarget(geneList, motif_rankings, motifAnnot=motifAnnotations)
  TFs <- lapply(split(enrichment$TF_highConf, enrichment$geneSet),
                function(x) {
                  genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                  genesSplit <- unique(unlist(strsplit(genes, "; ")))
                  return(genesSplit)})
  return(list(TFs, enrichment))
}

deg_p_val_enrichment <- motif_enrichment(wilcox_degs)
vmr_largest_change_enrichment <- motif_enrichment(ranked_abs_noise_dfs[[2]])

# filter out TF motif enrichment results that are not expressed in CMs
expressed_counts <- lapply(by_condition, function(x) rownames(GetAssayData(x, layer = "counts")[rowSums(GetAssayData(x, layer = "counts")) > 0, ]))
dp1tyb_expressed <- c(expressed_counts[[1]], expressed_counts[[2]], expressed_counts[[3]])

#' Filters out TFs that are not expressed at their respective time points
#' @param vmr_TFs list of TFs
#' @param deg_TFs list of TFs
#' @returns list of lists of 
TF_expression_filter <- function(vmr_TFs, deg_TFs) {
  dp <- intersect(setdiff(vmr_TFs$dp1tyb, deg_TFs$dp1tyb), dp1tyb_expressed)
  return(dp)
}

vmr_abs_change_filtTFs <- TF_expression_filter(vmr_largest_change_enrichment[[1]], deg_p_val_enrichment[[1]])

# Plot TF motif enrichment analysis
library(igraph)
library(ggraph)

#' Plots dendrogram of enriched TFs and their genes, colored by whether they are a DEG or not, and edges weighted by change in VMR
#' @param motif_df data.table of enriched motifs
#' @param expr_genes list of genes that are expressed at the given time point
#' @param degs data.frame of DEGs
#' @param noise data.frame of genes and their change in transcriptional variability
#' @returns list of all the enriched TFs that are expressed, all the enriched genes that are expressed, and a dendrogram plot
tf_gene_dendrogram <- function(motif_df, expr_genes, degs, noise) {
  edges <- motif_df %>%
    mutate(TF = strsplit(as.character(TF_highConf), ";")) %>%
    mutate(Gene = strsplit(as.character(enrichedGenes), ";")) %>%
    unnest(TF) %>%
    unnest(Gene) %>% 
    filter(TF %in% expr_genes & Gene %in% expr_genes) %>% 
    dplyr::select(motif, TF, Gene, NES)
  edges$TF <- gsub(" \\(inferredBy_Orthology\\).| \\(directAnnotation\\).", "", edges$TF)
  edges$gene_VMR <- noise[edges$Gene, ]$change
  
  all_tfs <- unique(edges$TF)
  all_genes <- unique(edges$Gene)
  nodes <- data.frame(
    name = c(all_tfs, all_genes),
    type = rep(c("TF", "Gene"), times = c(length(all_tfs), length(all_genes)))
  )
  nodes$is_DEG <- ifelse(nodes$name %in% rownames(degs), "DEG", "Not a DEG")
  
  g <- graph_from_data_frame(
    d = edges[, c("TF", "Gene")],  # Edge list (TF --> Gene)
    vertices = nodes,              # Node metadata
    directed = FALSE
  )
  
  edge_attr(g, "Change_in_VMR") <- edges$gene_VMR
  
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(alpha = Change_in_VMR), color = "gray70") +
    geom_node_point(aes(color = is_DEG, size = type), alpha = 0.8) +
    geom_node_text(
      aes(label = name,
          color = is_DEG),
      size = 3,
      repel = TRUE,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("DEG" = "gray40", "Not a DEG" = "blue")) +
    scale_size_manual(values = c("TF" = 5, "Gene" = 2)) +
    scale_edge_alpha(name = "Change in VMR") + 
    theme_void() +
    labs(color = "DEG status", size = "Node Type")
  
  return(list(all_tfs, all_genes, p))
}

vmr_tf_genes_t21 <- tf_gene_dendrogram(vmr_largest_change_enrichment[[2]], dp1tyb_expressed, wilcox_degs, ranked_abs_noise_dfs[[2]])
vmr_tf_genes_t21[[3]] # see the plot 

deg_tf_genes_t21 <- tf_gene_dendrogram(deg_p_val_enrichment[[2]], dp1tyb_expressed, wilcox_degs, ranked_abs_noise_dfs[[2]])
deg_tf_genes_t21[[3]]

## Investigate biological pathways##
mmu_kegg <- download_KEGG("mmu")
kegg_genes <- mmu_kegg$KEGGPATHID2EXTID

dp1tyb_all_genes_fc_ <- FoldChange(cm_only, ident.1 = c("WT_A1", "WT_A7"), ident.2 = "Dp1Tyb", recorrect_umi = FALSE)

# S2 Table
merge_noise_fc <- subset(merge(raw_noise_change_dfs[[2]], dp1tyb_all_genes_fc_, by = 0), select = -c(Row.names))

#' Generates scatter plot of average FC vs change in VMR, colored by whether a gene is a DEG 
#' @param df data.frame of the genes of interest and their change in transcriptional variability
#' @param timepoint_fc data.frame of FC of all genes at a given time point
#' @param timepoint_degs data.frame of DEGs at a given time point
#' @param neg_cut X-axis cutoff for labeling genes
#' @param pos_cut X-axis cutoff for labeling genes
#' @returns ggplot object
add_fc_and_scatterplot <- function(df, fc, deg, neg_cut, pos_cut) {
  comb <- merge(df, fc, by = 0, all.x = T, no.dups = T)
  comb <- subset(comb, select = -c(Row.names))
  degs <- intersect(rownames(df), rownames(deg))
  comb$is_DEG <- with(comb, ifelse(Row.names %in% degs, T, F))
  p <- ggplot(comb) + geom_point(aes(x = change, y = avg_log2FC, colour = is_DEG)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_text_repel(aes(x = change, y = avg_log2FC), size = 5, force = 2,
                    label = ifelse(comb$is_DEG == F & (comb$change < neg_cut | comb$change > pos_cut), as.character(comb$Row.names),""),
                    nudge_y = 0.05) +
    scale_color_manual(name = "DEG status", values = c("TRUE" = "#619CFF", "FALSE" = "#F8766D"), labels = c("Not a DEG", "DEG")) + 
    labs(x = "Change in VMR", y = "Average Log FC") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text=element_text(size=rel(1.2)), legend.title = element_text(size=rel(1.2)),
          axis.title = element_text(size = 16))
  return(p)
}

# investigate PI3K signaling pathway 
pi3k_kegg_genes <- kegg_genes[kegg_genes$from == "mmu04151", "to"]
pi3k_kegg_genes_sym <- mapIds(org.Mm.eg.db, keys=pi3k_kegg_genes, keytype="ENTREZID", column="SYMBOL") # 367
pi3k_genes_noise <- raw_noise_change_dfs[[2]][rownames(raw_noise_change_dfs[[2]]) %in% pi3k_kegg_genes_sym, ] # 149 rows
pi3k_genes_noise <- pi3k_genes_noise[order(abs(pi3k_genes_noise$change), decreasing = T), ] # ordered by abs change in VMR
pi3k_dot_p1 <- DotPlot(cm_only, features = rownames(pi3k_genes_noise)[1:30], idents =  c("WT_A7", "WT_A1", "Dp1Tyb"), split.by = "orig.ident", cols = "RdBu", scale = F) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
add_fc_and_scatterplot(pi3k_genes_noise, dp1tyb_all_genes_fc_, wilcox_degs, -0.25, 0.05)

# investigate ECM-receptor interactions
ecm_kegg_genes <- kegg_genes[kegg_genes$from == "mmu04512", "to"]
ecm_kegg_genes_sym <- mapIds(org.Mm.eg.db, keys=ecm_kegg_genes, keytype="ENTREZID", column="SYMBOL") # 89
ecm_genes_noise <- raw_noise_change_dfs[[2]][rownames(raw_noise_change_dfs[[2]]) %in% ecm_kegg_genes_sym, ] # 36 rows
ecm_genes_noise <- ecm_genes_noise[order(abs(ecm_genes_noise$change), decreasing = T), ] # ordered by abs change in VMR
ecm_dot <- DotPlot(cm_only, features = rownames(ecm_genes_noise), idents =  c("WT_A7", "WT_A1", "Dp1Tyb"), split.by = "orig.ident", cols = "RdBu", scale = F) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
add_fc_and_scatterplot(ecm_genes_noise, dp1tyb_all_genes_fc_, wilcox_degs, -0.25, 0.05)

# investigate focal adhesion
focal_kegg_genes <- kegg_genes[kegg_genes$from == "mmu04510", "to"]
focal_kegg_genes_sym <- mapIds(org.Mm.eg.db, keys=focal_kegg_genes, keytype="ENTREZID", column="SYMBOL") # 202
focal_genes_noise <- raw_noise_change_dfs[[2]][rownames(raw_noise_change_dfs[[2]]) %in% focal_kegg_genes_sym, ] # 111 rows
focal_genes_noise <- focal_genes_noise[order(abs(focal_genes_noise$change), decreasing = T), ] # ordered by abs change in VMR
focal_dot_1 <- DotPlot(cm_only, features = rownames(focal_genes_noise)[1:30], idents =  c("WT_A7", "WT_A1", "Dp1Tyb"), split.by = "orig.ident", cols = "RdBu", scale = F) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
add_fc_and_scatterplot(focal_genes_noise, dp1tyb_all_genes_fc_, wilcox_degs, -0.25, 0.05)
