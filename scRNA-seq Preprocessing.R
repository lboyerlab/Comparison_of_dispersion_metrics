library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(gridExtra)

## set up to load & save files ##
# setwd("~/MatHG")
# conditions <- c("con_E9_5", "con_E11_5", "matHG_E9_5", "matHG_E11_5")

setwd("~/T21")
conditions <- c("WT_A1", "WT_A7", "Dp1Tyb")

# QC cutoffs 
nFeature_min <- 500
nFeature_max <- 7000
percent_mt_cutoff <- 20 
nCount_min <- 800 # UMI cutoff 
min_cells <- 10

conditions_seurat <- list()

# create individual Seurat objects
for (x in 1:length(conditions)) {
  condition <- conditions[x]
  data <- Read10X(data.dir = condition)
  condition <- CreateSeuratObject(counts = data, project = condition, min.cells = min_cells, min.features = nFeature_min)
  condition[["percent.mt"]] <- PercentageFeatureSet(condition, pattern = "^mt-")
  condition <- subset(condition, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt < percent_mt_cutoff & nCount_RNA > nCount_min)
  conditions_seurat[[x]] <- condition
}

# combine into one Seurat object
combined <- merge(x = conditions_seurat[[1]], y = tail(conditions_seurat, -1), add.cell.ids = conditions)

# check the combined object 
head(combined@meta.data, 5) 
tail(combined@meta.data, 5) 

# normalize using sctransform (normalizes using log1p, replaces normalization, FindVariableFeatures, and ScaleData)
combined <- SCTransform(combined, vars.to.regress = "percent.mt", verbose = TRUE)

# check the 10 most highly variable genes
top10 <- head(VariableFeatures(combined), 10)
print(top10)

# PCA
combined <- RunPCA(combined)
print(combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(combined, dims = 1:2, reduction = "pca")
DimPlot(combined, reduction = "pca")

# TODO: determine the dimensionality of data set (how many PCs to keep), based on heatmaps & elbow plot
DimHeatmap(combined, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(combined, dims = 1:20, cells = 500, balanced = TRUE)
DimHeatmap(combined, dims = 15:25, cells = 500, balanced = TRUE)
ElbowPlot(combined)

# cluster the cells
n_dims <- 19 # TODO: update this w/ number of dimensions to use
combined <- FindNeighbors(combined, dims = 1:n_dims)
combined <- FindClusters(combined, resolution = 0.4) # NOTE: can change resolution if needed 

# check results of clustering
head(Idents(combined), 5)
table(combined@meta.data$seurat_clusters)

# UMAP
combined <- RunUMAP(combined, dims = 1:n_dims)
DimPlot(combined, reduction = "umap", label = FALSE)
DimPlot(combined, reduction = "umap", label = FALSE, group.by = "orig.ident")

# otherwise, find cluster markers (compared to all other cells), report only the positive ones 
combined <- PrepSCTFindMarkers(combined)
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# marker genes for CMs
cm_markers <- c("Myh7","Ryr2","Ttn","Mybpc3","Actn2","Tnnc1","Actc1","Myh6","Tnnt2","Acta2", "Acta1", "Tnnc3")
ep_markers <- c("Wt1", "Alcam", "Dcn", "Upk3b", "Tbx18")
nc_markers <- c("Nr2f1", "Dlx5", "Dlx2", "Serpinf1", "Prrx1")
ec_markers <- c("Plvap", "Egfl7", "Klf2", "Emcm", "Cdh5")
fm_markers <- c("Vcan", "Postn", "Papss2", "Cthrc1", "Sox9")

# different ways to examine the cluster markers:
# finds top 10 markers for every cluster
combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  print(n=200) -> top10_markers

# look at markers of a specific cluster
combined.markers %>%
  group_by(cluster) %>%
  filter(cluster == 6) %>% # TODO: change this to get the specific cluster
  slice_max(n = 30, order_by = avg_log2FC) %>%
  print(n=60)

# look at expression of specific genes across clusters 
combined.markers %>%
  group_by(cluster) %>%
  filter(gene == "Col11a1") %>% 
  print(n=100)

# look at expression of CM marker genes across clusters 
combined.markers %>%
  group_by(cluster) %>%
  filter(gene %in% cm_markers) %>% 
  print(n=100) -> cm_markers_in_clusters

# look at expression of other marker genes 
combined.markers %>%
  group_by(cluster) %>%
  filter(gene %in% nc_markers) %>% 
  print(n=100) 

# use violin plots to visualize expression of markers across clusters
VlnPlot(combined, features = cm_markers)
VlnPlot(combined, features = c("Ryr2"), layer = "counts", log = TRUE) # plots the raw counts 

# Heat map with the top genes from each cluster
DoHeatmap(combined, features = top10_markers$gene) + NoLegend() 

# save 
path_to_data <- ""
write.csv(top10_markers, paste(path_to_data, "top10_cluster_markers.csv"))
write.csv(cm_markers_in_clusters, paste(path_to_data, "cm_markers_in_clusters.csv"))
saveRDS(combined, file = paste(path_to_data, "combined.rds"))
saveRDS(combined.markers, file = paste(path_to_data, "combined_markers.rds"))

# TODO: rename using cell types 
new.cluster.ids <- c("CM", "CM", "CM", "EC", "?", "FB",
                     "?", "EP", "FB", "EC", "MP")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)

# subset of only CM
cm_only <- subset(combined, ident = "CM")
saveRDS(cm_only, file = paste(path_to_data, "cm_only.rds"))


