# packages - with sttkit
library(sttkit)
library(org.Hs.eg.db)
library(patchwork)
library(scater)
library(Seurat)
# 
orgdb <- org.Hs.eg.db

# Define where you want the cell-phone output to be saved
prefix <- "/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tokura et al DCIS SC dataset/cellphoneDB"

# Prepare data and export for CellPhoneDB
cellphone_for_seurat(DCIS, orgdb, prefix)


# Define the output path where the `cellphonedb` analysis results locate
cellphone_outpath <- "/Users/innaa/cellphonedb/results/method1"

# Import and rerank CellPhoneDB Output
res <- import_cellphone(DCIS, cellphone_outpath, orgdb)


# Plot cell-cell interactions
# 'id' is the line of the significant pair in 'significant_means_ranked_spatial.csv' file
plot <- plot_cellphone(DCIS, id = 1, cellphone_outpath)
print(plot)

# packages - with ktplots
library(ktplots)
DCIS_sc = as.SingleCellExperiment(DCIS)
p1 <- plotExpression(DCIS_sc, features = "MS4A1", x = "ARTCLASS") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                   hjust = 1))
p2 <- plotPCA(DCIS_sc, colour_by = "ARTCLASS")
p1 + p2

pvals <- read.delim("/Users/innaa/cellphonedb/results/method2/statistical_analysis_pvalues_06_18_2024_140357.txt", check.names = FALSE)
means <- read.delim("/Users/innaa/cellphonedb/results/method2/statistical_analysis_means_06_18_2024_140357.txt", check.names = FALSE)
decon = read.delim("/Users/innaa/cellphonedb/results/method2/statistical_analysis_deconvoluted_06_18_2024_140357.txt", sep="\t")

plot_cpdb_heatmap(pvals = pvals, cellheight = 10, cellwidth = 10, symmetrical = FALSE)

plot_cpdb(
  scdata=DCIS_sc,
  cell_type1="Stroma",
  cell_type2="Macrophage",
  celltype_key="ARTCLASS",
  means=means,
  pvals=pvals,
)
