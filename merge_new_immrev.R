# Install packages
devtools::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)

# Load objects
IMMrev_new = dcis_immune_11_label
sDCISrev <- readRDS("sDCISrev.rds")
head(IMMrev_new@meta.data)
head(sDCISrev@meta.data$subcluster)
unique(sDCISrev@meta.data$ARTCLASS)

# Add new subcluster column from IMMrev_new to sDCISrev
immune_metadata <- IMMrev_new@meta.data[ , "subcluster", drop = FALSE]
sDCISrev@meta.data <- cbind(
  sDCISrev@meta.data,
  subcluster = immune_metadata[rownames(sDCISrev@meta.data), "subcluster"]
)

# Create a vector that maps cluster numbers to cell types
cluster_mapping <- c(
  "40" = "B Cell",
  "41" = "CD4+ T Cell",
  "42" = "CD8+ T Cell",
  "43" = "Treg",
  "44" = "NK/NKT",
  "45" = "Resting T Cell",
  "46" = "CD14+ Monocyte",
  "47" = "Macrophage",
  "48" = "Erythrocyte",
  "49" = "Neutrophil",
  "50" = "Proliferating Immune Cell",
  "51" = "Plasma",
  "52" = "Unclassified1",
  "53" = "Unclassified2",
  "54" = "Unclassified3",
  "55" = "Unclassified4"
)

# Label change
sDCISrev@meta.data$subcluster <- cluster_mapping[as.character(sDCISrev@meta.data$subcluster)]
sDCISrev@meta.data$subcluster[is.na(sDCISrev@meta.data$subcluster)] <- "Non-Immune"
table(sDCISrev@meta.data$subcluster)

# Merge epithelial cell types from ARTCLASS into subcluster
immune_labels <- c(
  "B Cell", "CD4+ T Cell", "CD8+ T Cell", "Treg", "NK/NKT", "Resting T Cell",
  "CD14+ Monocyte", "Macrophage", "Erythrocyte", "Neutrophil",
  "Proliferating Immune Cell", "Plasma", "Unclassified1", "Unclassified2",
  "Unclassified3", "Unclassified4"
)
epithelial_labels <- c("LuminalEpithelia", "Basal", "Stroma", "ProliferatingBasal", "ProliferatingLuminal")

sDCISrev@meta.data$subcluster[!(sDCISrev@meta.data$subcluster %in% immune_labels) &
                                sDCISrev@meta.data$ARTCLASS %in% epithelial_labels] <- 
  sDCISrev@meta.data$ARTCLASS[!(sDCISrev@meta.data$subcluster %in% immune_labels) &
                                sDCISrev@meta.data$ARTCLASS %in% epithelial_labels]

# Replace "Non-Immune" labels with NA in the subcluster column
sDCISrev@meta.data$subcluster[sDCISrev@meta.data$subcluster == "Non-Immune"] <- NA

# Verify the updated subcluster column
table(sDCISrev@meta.data$subcluster)
unique(sDCISrev@meta.data$subcluster)

# Save new RDS file
saveRDS(sDCISrev, file = "sDCISrev_new.rds")

# Save as AnnData file
SaveH5Seurat(sDCISrev, filename = "sDCISrev_new.h5seurat", overwrite = TRUE)
Convert("sDCISrev_new.h5seurat", dest = "h5ad")

### MEMORY ERROR ###

## Convert counts and data matrices to sparse format
sDCISrev[["RNA"]]@counts <- as.sparse(sDCISrev[["RNA"]]@counts)
sDCISrev[["RNA"]]@data <- as.sparse(sDCISrev[["RNA"]]@data)
IMMrev_new[["RNA"]]@counts <- as.sparse(IMMrev_new[["RNA"]]@counts)
IMMrev_new[["RNA"]]@data <- as.sparse(IMMrev_new[["RNA"]]@data)

# Inspect metadata
head(IMMrev_new@meta.data)
head(sDCISrev@meta.data)

# Transfer updated metadata from IMMrev to sDCISrev
common_cells <- intersect(rownames(IMMrev_new@meta.data), rownames(sDCISrev@meta.data))
sDCISrev@meta.data[common_cells, ] <- IMMrev_new@meta.data[common_cells, ]

# Subset immune cells in sDCISrev
new_immune_cells <- colnames(IMMrev_new)
sDCISrev <- subset(sDCISrev, cells = setdiff(colnames(sDCISrev), new_immune_cells)) #ERROR

##Trial 1 Divide into batches
cell_batches <- split(colnames(IMMrev_new), ceiling(seq_along(colnames(IMMrev_new)) / 5000))  # Adjust batch size
for (batch in cell_batches) {
  batch_subset <- subset(IMMrev_new, cells = batch)  # Subset IMMrev_new by batch
  sDCISrev <- merge(sDCISrev, y = batch_subset)  # Merge into sDCISrev
  gc()  # Clear memory
}

##Trial 2 Identify cells to keep
cells_to_keep <- setdiff(colnames(sDCISrev), new_immune_cells)
chunk_size <- 5000  
cell_chunks <- split(cells_to_keep, ceiling(seq_along(cells_to_keep) / chunk_size))
filtered_sDCISrev <- NULL
for (chunk in cell_chunks) {
  chunk_subset <- sDCISrev[, chunk]
  if (is.null(filtered_sDCISrev)) {
    filtered_sDCISrev <- chunk_subset
  } else {
    filtered_sDCISrev <- merge(filtered_sDCISrev, y = chunk_subset)
  }
  gc()  
}
sDCISrev <- sDCISrev[, cells_to_keep]

# Merge immune cells back
sDCISrev_new <- merge(sDCISrev, y = IMMrev_new, project = "sDCISrev_new")
head(sDCISrev_new@meta.data) # Verify metadata
sDCISrev_new # Verify dimensions

# Save new RDS file
saveRDS(sDCISrev_new, file = "sDCISrev_new.rds")
