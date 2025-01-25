### PACKAGE INSTALLATION ###

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)

### DATASET INITIALIZATION ###

immune_cells = readRDS("/Users/eunicelee153/Desktop/WORK/CLINICAL/Angelo Lab/DCIS Project/IMMrev.rds")
dcis_immune <- subset(immune_cells, subset = orig.ident %in% c("DCIS1", "DCIS2","DCIS3","DCIS4","DCIS5","DCIS6","DCIS7"))

### QUALITY CONTROL, NORMALIZATION, SCALING ###

dcis_immune[["percent.mt"]] <- PercentageFeatureSet(dcis_immune, pattern = "^MT-")
dcis_immune_qc <- subset(dcis_immune, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 20) # Tokura Paper
dcis_immune_norm <- NormalizeData(dcis_immune_qc)
dcis_immune_highvar <- FindVariableFeatures(dcis_immune_norm, selection.method = "vst", nfeatures = 2000)
all_genes <- rownames(dcis_immune_highvar)
dcis_immune_scaled <- ScaleData(dcis_immune_highvar, features = all_genes)

### DATA CLUSTERING, DIMENSIONALITY REDUCTION ###
dcis_immune_11 <- RunPCA(dcis_immune_scaled, features = VariableFeatures(object = dcis_immune_scaled))
VizDimLoadings(dcis_immune_11, dims = 1:2, reduction = "pca")
DimHeatmap(dcis_immune_11, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(dcis_immune_11)
dcis_immune_11 <- FindNeighbors(dcis_immune_11, dims = 1:12)
dcis_immune_11 <- FindClusters(dcis_immune_11, resolution = 4.0)
dcis_immune_11 <- RunUMAP(dcis_immune_11, dims = 1:12)
DimPlot(dcis_immune_11, group.by = "ident", reduction = "umap", label = TRUE)

DimPlot(dcis_immune_11, group.by = "IMMCLASS", reduction = "umap")

### MARKERS IDENTIFICATION ###

# Resolution 4.0
testcluster19 <- FindMarkers(dcis_immune_11, ident.1 = 19, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster1 <- FindMarkers(dcis_immune_11, ident.1 = 1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster3 <- FindMarkers(dcis_immune_11, ident.1 = 3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster17 <- FindMarkers(dcis_immune_11, ident.1 = 17, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster16 <- FindMarkers(dcis_immune_11, ident.1 = 16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster23 <- FindMarkers(dcis_immune_11, ident.1 = 23, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster21 <- FindMarkers(dcis_immune_11, ident.1 = 21, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster27 <- FindMarkers(dcis_immune_11, ident.1 = 27, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster36 <- FindMarkers(dcis_immune_11, ident.1 = 36, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster19, n = 100) #Unclassified
head(testcluster1, n = 100) #Unclassified
head(testcluster3, n = 50) #CD14+ monocyte
head(testcluster17, n = 50) #CD14+ monocyte
head(testcluster16, n = 50) #Proliferating immune cell
head(testcluster23, n = 50) #Macrophage
head(testcluster21, n = 50) #Macrophage
head(testcluster27, n = 50) #Macrophage
head(testcluster36, n = 50) #Neutrophil
testcluster25 <- FindMarkers(dcis_immune_11, ident.1 = 25, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster29 <- FindMarkers(dcis_immune_11, ident.1 = 29, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster26 <- FindMarkers(dcis_immune_11, ident.1 = 26, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster12 <- FindMarkers(dcis_immune_11, ident.1 = 12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster28 <- FindMarkers(dcis_immune_11, ident.1 = 28, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster25, n = 50) #CD8+ T cell
head(testcluster29, n = 70) #Unclassified
head(testcluster26, n = 70) #Unclassified
head(testcluster12, n = 70) #Unclassified
head(testcluster28, n = 70) #B
testcluster2 <- FindMarkers(dcis_immune_11, ident.1 = 2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster8 <- FindMarkers(dcis_immune_11, ident.1 = 8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster10 <- FindMarkers(dcis_immune_11, ident.1 = 10, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster5 <- FindMarkers(dcis_immune_11, ident.1 = 5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster34 <- FindMarkers(dcis_immune_11, ident.1 = 34, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster2, n = 20) #CD8+
head(testcluster8, n = 20) #CD8+
head(testcluster10, n = 20) #CD8+
head(testcluster5, n = 50) #CD4+
head(testcluster34, n = 50) #Unclassified
testcluster22 <- FindMarkers(dcis_immune_11, ident.1 = 22, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster7 <- FindMarkers(dcis_immune_11, ident.1 = 7, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster22, n = 50) #CD4+
head(testcluster7, n = 50) #CD4+
testcluster6 <- FindMarkers(dcis_immune_11, ident.1 = 6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster6, n = 20) #Treg
testcluster35 <- FindMarkers(dcis_immune_11, ident.1 = 35, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster18 <- FindMarkers(dcis_immune_11, ident.1 = 18, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster35, n = 70) #Resting T
head(testcluster18, n = 70) #Resting T
testcluster20 <- FindMarkers(dcis_immune_11, ident.1 = 20, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster14 <- FindMarkers(dcis_immune_11, ident.1 = 14, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster20, n = 20) #Erythrocyte
head(testcluster14, n = 20) #Erythrocyte
testcluster4 <- FindMarkers(dcis_immune_11, ident.1 = 4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster4, n = 20) #Plasma
testcluster31 <- FindMarkers(dcis_immune_11, ident.1 = 31, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster13 <- FindMarkers(dcis_immune_11, ident.1 = 13, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster31, n = 20) #NK/NKT
head(testcluster13, n = 20) #NK/NKT
testcluster0 <- FindMarkers(dcis_immune_11, ident.1 = 0, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster30 <- FindMarkers(dcis_immune_11, ident.1 = 30, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster24 <- FindMarkers(dcis_immune_11, ident.1 = 24, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster15 <- FindMarkers(dcis_immune_11, ident.1 = 15, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster11 <- FindMarkers(dcis_immune_11, ident.1 = 11, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster33 <- FindMarkers(dcis_immune_11, ident.1 = 33, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster32 <- FindMarkers(dcis_immune_11, ident.1 = 32, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster9 <- FindMarkers(dcis_immune_11, ident.1 = 9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster0, n = 20) #B
head(testcluster30, n = 50) #B
head(testcluster24, n = 20) #B
head(testcluster15, n = 20) #B
head(testcluster11, n = 20) #B
head(testcluster33, n = 50) #B
head(testcluster32, n = 50) #B
head(testcluster9, n = 50) #B

# Resolution 1.0
cluster0 <- FindMarkers(dcis_immune_11, ident.1 = 0, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster1 <- FindMarkers(dcis_immune_11, ident.1 = 1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster2 <- FindMarkers(dcis_immune_11, ident.1 = 2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster3 <- FindMarkers(dcis_immune_11, ident.1 = 3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster4 <- FindMarkers(dcis_immune_11, ident.1 = 4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster5 <- FindMarkers(dcis_immune_11, ident.1 = 5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster6 <- FindMarkers(dcis_immune_11, ident.1 = 6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster7 <- FindMarkers(dcis_immune_11, ident.1 = 7, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster8 <- FindMarkers(dcis_immune_11, ident.1 = 8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster9 <- FindMarkers(dcis_immune_11, ident.1 = 9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster10 <- FindMarkers(dcis_immune_11, ident.1 = 10, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster11 <- FindMarkers(dcis_immune_11, ident.1 = 11, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster12 <- FindMarkers(dcis_immune_11, ident.1 = 12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster13 <- FindMarkers(dcis_immune_11, ident.1 = 13, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster14 <- FindMarkers(dcis_immune_11, ident.1 = 14, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster15 <- FindMarkers(dcis_immune_11, ident.1 = 15, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster16 <- FindMarkers(dcis_immune_11, ident.1 = 16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster17 <- FindMarkers(dcis_immune_11, ident.1 = 17, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster18 <- FindMarkers(dcis_immune_11, ident.1 = 18, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(cluster0, n = 100) # CD14+ monocytes
head(cluster1, n = 20) # Resting T cell
head(cluster2, n = 20) # CD8+ T cell
head(cluster3, n = 50) # Macrophage & Neutrophil
head(cluster4, n = 50) # B cell
head(cluster5, n = 20) # CD4+ T cell
head(cluster6, n = 50) # B cell
head(cluster7, n = 100) # Unclassified (both on mine and theirs)
head(cluster8, n = 20) # Erythrocyte
head(cluster9, n = 100) # B cell
head(cluster10, n = 50) # NK/NKT
head(cluster11, n = 20) # Plasma
head(cluster12, n = 20) # Treg
head(cluster13, n = 50) # Resting T cell
head(cluster14, n = 20) # Proliferating immune cell
head(cluster15, n = 50) # CD8+ T cell
head(cluster16, n = 20) # B cell
head(cluster17, n = 100) # Unclassified (both on mine and theirs)
head(cluster18, n = 100) # B cell

# Resolution 2.0
testcluster0 <- FindMarkers(dcis_immune_11, ident.1 = 0, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster8 <- FindMarkers(dcis_immune_11, ident.1 = 8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster17 <- FindMarkers(dcis_immune_11, ident.1 = 17, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster14 <- FindMarkers(dcis_immune_11, ident.1 = 14, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster20 <- FindMarkers(dcis_immune_11, ident.1 = 20, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster22 <- FindMarkers(dcis_immune_11, ident.1 = 22, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster13 <- FindMarkers(dcis_immune_11, ident.1 = 13, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster15 <- FindMarkers(dcis_immune_11, ident.1 = 15, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster0, n = 100) #Unclassified
head(testcluster8, n = 100) #CD14
head(testcluster17, n = 100) #CD14
head(testcluster14, n = 100) #Macrophage
head(testcluster20, n = 100) #Macrophage
head(testcluster22, n = 100) #Neutrophil/Macrophage
head(testcluster13, n = 100) #Unclassified
head(testcluster15, n = 100) #Unclassified

# Resolution 3.0 
testcluster18 <- FindMarkers(dcis_immune_11, ident.1 = 22, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster22 <- FindMarkers(dcis_immune_11, ident.1 = 22, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testcluster28 <- FindMarkers(dcis_immune_11, ident.1 = 22, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(testcluster18, n = 50) #Macrophage
head(testcluster22, n = 75) #Neutrophil/Macrophage?
head(testcluster28, n = 75) #Neutrophil/Macrophage?

# Identify highly variable genes
head(VariableFeatures(dcis_immune_11), 10)

FeaturePlot(dcis_immune_11, features = c("CD14", "CD68"))
VlnPlot(dcis_immune_11, features = c("CD14", "CD68"))
FeaturePlot(dcis_immune_11, features = c("FCGR3A", "LYZ", "CSF3R", "S100A8", "S100A9", "CXCR1", "CD33"))
VlnPlot(dcis_immune_11, features = c("FCGR3A", "LYZ", "CSF3R", "S100A8", "S100A9", "CXCR1", "CD33"))
FeaturePlot(dcis_immune_11, features = c("MKI67"))
VlnPlot(dcis_immune_11, features = c("MKI67"))
FeaturePlot(dcis_immune_11, features = cd14monocyte)
VlnPlot(dcis_immune_11, features = cd14monocyte)
FeaturePlot(dcis_immune_11, features = c("CD4", "CD8A", "CD8B"))
VlnPlot(dcis_immune_11, features = c("CD4", "CD8A", "CD8B"))
FeaturePlot(dcis_immune_11, features = c("ITGAX", "CD68", "CD163", "MRC1"))
VlnPlot(dcis_immune_11, features = c("ITGAX", "CD68", "CD163", "MRC1"))
FeaturePlot(dcis_immune_11, features = c("ANXA1", "ITGAX", "IL3RA", "CPA3", "MS4A2"))
VlnPlot(dcis_immune_11, features = c("ITGAX", "CD68", "CD163", "MRC1"))

### MERGING ###

# B Cell = 40
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(0, 30, 24, 15, 11, 28, 33, 9, 32)), value = 40)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# CD4+ T Cell = 41
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(22, 7, 5)), value = 41)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# CD8+ T Cell = 42
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(2, 8, 10, 25)), value = 42)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# Treg = 43
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = 6), value = 43)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# NK/NKT = 44
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(13, 31)), value = 44)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# Resting T Cell = 45
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(35, 18)), value = 45)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# CD14+ Monocyte = 46
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(3, 17)), value = 46)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# Macrophage = 47
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(23, 21, 27)), value = 47)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# Erythrocyte = 48
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(20, 14)), value = 48)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# Neutrophil = 49
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = 36), value = 49)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# Proliferating Immune Cell = 50
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = 16), value = 50)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# Plasma = 51
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = 4), value = 51)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

# Unclassified = 52, 53, 54, 55
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(34, 26)), value = 52)
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = c(19, 1)), value = 53)
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = 12), value = 54)
dcis_immune_11 <- SetIdent(dcis_immune_11, cells = WhichCells(dcis_immune_11, idents = 29), value = 55)
dcis_immune_11$subcluster <- as.character(Idents(dcis_immune_11))
DimPlot(dcis_immune_11, reduction = "umap", group.by = "subcluster")

### FINAL PLOTTING ###

# Labeling
dcis_immune_11 <- SetIdent(dcis_immune_11, value = dcis_immune_11$subcluster)
levels(Idents(dcis_immune_11))

Idents(dcis_immune_11) <- factor(Idents(dcis_immune_11), 
                                 levels = c("41", "42", "44", "43", "45", "48",  "46", "47", "49", "50", "40", "51", "52", "53", "54", "55"))
new.cluster.ids <- c("CD4+ T Cell", "CD8+ T Cell", "NK/NKT", "Treg", "Resting T cell", "Erythrocyte", "CD14+ Monocyte", "Macrophage", "Neutrophil", "Proliferating Immune Cell", "B Cell", "Plasma", "Unclassified1", "Unclassified2", "Unclassified3", "Unclassified4")
names(new.cluster.ids) <- levels(dcis_immune_11)
dcis_immune_11_label <- RenameIdents(dcis_immune_11, new.cluster.ids)
Idents(dcis_immune_11_label)[Idents(dcis_immune_11_label) == "Unclassified1"] <- NA
Idents(dcis_immune_11_label)[Idents(dcis_immune_11_label) == "Unclassified2"] <- NA
Idents(dcis_immune_11_label)[Idents(dcis_immune_11_label) == "Unclassified3"] <- NA
Idents(dcis_immune_11_label)[Idents(dcis_immune_11_label) == "Unclassified4"] <- NA

DimPlot(dcis_immune_11_label, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 4)

# Save Plot
finalplot <- DimPlot(dcis_immune_11_label, reduction = "umap", label = TRUE, label.size = 4) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "immune_umap_4_EL.jpg", height = 7, width = 12, plot = finalplot, quality = 50)

### SCRATCH PAD ###

subset_markers <- FindAllMarkers(dcis_immune_11, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DimPlot(dcis_immune_scaled, group.by = "IMMCLASS")
table(dcis_immune_scaled@meta.data$IMMCLASS)
length(unique(dcis_immune_scaled@meta.data$IMMCLASS))

immclass_cluster <- c("B Cells", "Plasma Cells", "CD14+ Monocytes", "CD8+ T Cells", "Erythrocytes", "Macrophages", "Neutrophiles", "NK/NKT", "Proliferating Immune Cells", "T Cells", "Treg", "Unclassified")
names(immclass_cluster) <- levels(dcis_immune_12)
dcis_immune_rename <- RenameIdents(dcis_immune_12, immclass_cluster)
DimPlot(dcis_immune_rename, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Sub-clustering for Treg: Method 1
dcis_immune_11[["main.cluster"]] <- Idents(dcis_immune_11);
dcis_immune_11_sub1 <- subset(dcis_immune_11, subset = IMMCLASS == "Treg");
Idents(dcis_immune_11) <- "main.cluster";
Idents(dcis_immune_11, cells=colnames(dcis_immune_11)[colnames(dcis_immune_11) %in% colnames(dcis_immune_11_sub1)]) <- paste0("sub1_", Idents(dcis_immune_11_sub1))
dcis_immune_11[["merged.cluster"]] <- Idents(dcis_immune_11);
DimPlot(dcis_immune_11, group.by = "ident")
