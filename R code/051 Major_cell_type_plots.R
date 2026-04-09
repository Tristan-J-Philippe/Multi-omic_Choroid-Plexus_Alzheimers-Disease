###
### 061 Spatial_proxmity.R
### This script was written by Dr. Denis R Avey
# Purpose: Calculate spatial proximity and distance-based features between cell types.
# Dependencies:
source("/050 ST_specific_plots_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder

# Example of checking if points are within polygons
ChP <- readRDS("ChP_combined.rds")
annot_meta <- ChP@meta.data
rownames(annot_meta) <- names(ChP@active.ident) # non-unique values detected in @meta.data, so use rownames of active.ident instead

### cell polygon coordinates
# Filter cell IDs to include only those ending with '1' (slide 1) or '2' (slide 2)
cell_ids_slide1 <- annot_meta$unique_cellID[grepl("1$", annot_meta$unique_cellID)]
cell_ids_slide2 <- annot_meta$unique_cellID[grepl("2$", annot_meta$unique_cellID)]

# Create a list to store the polygon coordinates for each cell
cell_boundaries1 <- list()

# Iterate over all cell IDs for Slide 1
for (cell_id in cell_ids_slide1) {
  # Construct the full path to the polygon for the given cell
  cell_polygon <- ChP@images[["Slide.1"]]@boundaries[["segmentation"]]@polygons[[cell_id]]@Polygons[[1]]@coords
  
  # Store the polygon coordinates in the list
  cell_boundaries1[[cell_id]] <- cell_polygon
}

# Create a list to store the polygon coordinates for each cell
cell_boundaries2 <- list()

# Iterate over all cell IDs for Slide 2
for (cell_id in cell_ids_slide2) {
  # Construct the full path to the polygon for the given cell
  cell_polygon <- ChP@images[["Slide.2"]]@boundaries[["segmentation"]]@polygons[[cell_id]]@Polygons[[1]]@coords
  
  # Store the polygon coordinates in the list
  cell_boundaries2[[cell_id]] <- cell_polygon
}

### molecule coordinates
# Get the list of all gene names
gene_names <- names(ChP@images[["Slide.1"]]@molecules[["molecules"]])

# Initialize an empty list to store molecule coordinates for all genes
all_molecule_coords <- list()

# Iterate through each gene and extract its coordinates
for (gene in gene_names) {
  # Extract coordinates for the current gene
  gene_coords <- ChP@images[["Slide.1"]]@molecules[["molecules"]][[gene]]@coords
  
  # Store the coordinates in the list with the gene name as the key
  all_molecule_coords[[gene]] <- gene_coords
}

# Combine all coordinates into a single data frame (with a gene name column)
all_molecules <- do.call(rbind, lapply(names(all_molecule_coords), function(gene) {
  data.frame(gene = gene, x = all_molecule_coords[[gene]][, 1], y = all_molecule_coords[[gene]][, 2])
}))

# Repeat for slide 2
gene_names2 <- names(ChP@images[["Slide.2"]]@molecules[["molecules"]])

# Initialize an empty list to store molecule coordinates for all genes
all_molecule_coords2 <- list()

# Iterate through each gene and extract its coordinates
for (gene in gene_names2) {
  # Extract coordinates for the current gene
  gene_coords <- ChP@images[["Slide.2"]]@molecules[["molecules"]][[gene]]@coords
  
  # Store the coordinates in the list with the gene name as the key
  all_molecule_coords2[[gene]] <- gene_coords
}

# Combine all coordinates into a single data frame (with a gene name column)
all_molecules2 <- do.call(rbind, lapply(names(all_molecule_coords2), function(gene) {
  data.frame(gene = gene, x = all_molecule_coords2[[gene]][, 1], y = all_molecule_coords2[[gene]][, 2])
}))

# save cell and molecule coordinate files for future use
saveRDS(cell_boundaries1, file = "ChP_cell_boundaries_slide1.rds")
saveRDS(cell_boundaries2, file = "ChP_cell_boundaries_slide2.rds")
saveRDS(all_molecules, file = "ChP_molecule_coords_slide1.rds")
saveRDS(all_molecules2, file = "ChP_molecule_coords_slide2.rds")

cell_boundaries1 <- readRDS("ChP_cell_boundaries_slide1.rds")
cell_boundaries2 <- readRDS("ChP_cell_boundaries_slide2.rds")
all_molecules <- readRDS("ChP_molecule_coords_slide1.rds")
all_molecules2 <- readRDS("ChP_molecule_coords_slide2.rds")

## Slide 1 ##
# Convert cell boundaries to an sf object
cell_sf <- lapply(names(cell_boundaries1), function(cell_id) {
  st_polygon(list(cell_boundaries1[[cell_id]])) %>% 
    st_sfc(crs = 32633) %>% 
    st_sf(cell_id = cell_id)
}) %>%
  bind_rows()

# Convert molecule coordinates to an sf object
molecule_sf <- st_as_sf(all_molecules, coords = c("x", "y"), crs = 32633)

# Perform spatial join to match molecules with cells
filtered_molecules_sf <- st_join(molecule_sf, cell_sf, left = FALSE) # from cell_sf to here takes ~ 100 GB of RAM; ~1.5 hrs

# Convert back to a data frame (optional)
filtered_molecules1 <- as.data.frame(filtered_molecules_sf)

# Now 'filtered_molecules1' contains only the molecules within cell boundaries for Slide 1
saveRDS(filtered_molecules1, file = "ChP_filtered_molecules_slide1.rds")

## Slide 2 ##
cell_sf <- lapply(names(cell_boundaries2), function(cell_id) {
  st_polygon(list(cell_boundaries2[[cell_id]])) %>% 
    st_sfc(crs = 32633) %>% 
    st_sf(cell_id = cell_id)
}) %>%
  bind_rows()

# Convert molecule coordinates to an sf object
molecule_sf <- st_as_sf(all_molecules2, coords = c("x", "y"), crs = 32633)

# Perform spatial join to match molecules with cells
filtered_molecules_sf <- st_join(molecule_sf, cell_sf, left = FALSE)

# Convert back to a data frame (optional)
filtered_molecules2 <- as.data.frame(filtered_molecules_sf)

# Now 'filtered_molecules2' contains only the molecules within cell boundaries for Slide 2
saveRDS(filtered_molecules2, file = "ChP_filtered_molecules_slide2.rds")

# Now, to filter the Seurat object:
# First, retrieve the molecules that correspond to the filtered molecule coordinates
molecules1 <- readRDS(file = "ChP_filtered_molecules_slide1.rds")
molecules2 <- readRDS(file = "ChP_filtered_molecules_slide2.rds")

# Subset the Seurat object by selecting only the filtered molecules (those that are inside cell boundaries)
filtered_seurat <- readRDS("ChP_combined.rds")

# single gene test:
test_gene <- "A1BG"
test_coords <- filter_molecules_for_gene(test_gene, molecules1, filtered_seurat@images[["Slide.1"]])
if (!is.null(test_coords)) {
  cat("Filtered coordinates for gene:", test_gene, "\n")
  print(head(test_coords))
} else {
  cat("No matching coordinates found for gene:", test_gene, "\n")
}

cat("Checking for duplicates in test_coords...\n")
cat("Number of unique coordinates in gene_coords:", nrow(unique(test_coords)), "\n")
cat("Total rows in gene_coords:", nrow(test_coords), "\n")



# Apply to both slides
filtered_seurat@images[["Slide.1"]] <- filter_molecules(molecules1, filtered_seurat@images[["Slide.1"]])
filtered_seurat@images[["Slide.2"]] <- filter_molecules(molecules2, filtered_seurat@images[["Slide.2"]])

# Save the filtered Seurat object
saveRDS(filtered_seurat, file = "ChP_combined_filtered_coordinates.rds")


### repeat for specific FOV(s) ###
ChP <- readRDS(file = "ChP_combined_filtered_coordinates.rds")
annot_meta <- ChP@meta.data
rownames(annot_meta) <- names(ChP@active.ident) # non-unique values detected in @meta.data, so use rownames of active.ident instead

### cell polygon coordinates
# Filter cell IDs to include only those ending with '1' (slide 1) or '2' (slide 2)
cell_ids_slide1 <- annot_meta$unique_cellID[
  grepl("^c_1_89", annot_meta$unique_cellID) & grepl("1$", annot_meta$unique_cellID)
]

# Create a list to store the polygon coordinates for each cell
cell_boundaries1 <- list()

# Iterate over all cell IDs for Slide 1
for (cell_id in cell_ids_slide1) {
  # Construct the full path to the polygon for the given cell
  cell_polygon <- ChP@images[["Slide.1.89"]]@boundaries[["segmentation"]]@polygons[[cell_id]]@Polygons[[1]]@coords
  
  # Store the polygon coordinates in the list
  cell_boundaries1[[cell_id]] <- cell_polygon
}

### molecule coordinates
# Get the list of all gene names
gene_names <- names(ChP@images[["Slide.1.89"]]@molecules[["molecules"]])

# Initialize an empty list to store molecule coordinates for all genes
all_molecule_coords <- list()

# Iterate through each gene and extract its coordinates
for (gene in gene_names) {
  # Extract coordinates for the current gene
  gene_coords <- ChP@images[["Slide.1.89"]]@molecules[["molecules"]][[gene]]@coords
  
  # Store the coordinates in the list with the gene name as the key
  all_molecule_coords[[gene]] <- gene_coords
}

# Combine all coordinates into a single data frame (with a gene name column)
all_molecules <- do.call(rbind, lapply(names(all_molecule_coords), function(gene) {
  data.frame(gene = gene, x = all_molecule_coords[[gene]][, 1], y = all_molecule_coords[[gene]][, 2])
}))

# save cell and molecule coordinate files for future use
saveRDS(cell_boundaries1, file = "ChP_cell_boundaries_slide1.89.rds")
saveRDS(all_molecules, file = "ChP_molecule_coords_slide1.89.rds")

cell_boundaries1 <- readRDS("ChP_cell_boundaries_slide1.89.rds")
all_molecules <- readRDS("ChP_molecule_coords_slide1.89.rds")

# filter for molecules in spots #
## Slide 1 ##
# Convert cell boundaries to an sf object
cell_sf <- lapply(names(cell_boundaries1), function(cell_id) {
  st_polygon(list(cell_boundaries1[[cell_id]])) %>% 
    st_sfc(crs = 32633) %>% 
    st_sf(cell_id = cell_id)
}) %>%
  bind_rows()

# Convert molecule coordinates to an sf object
molecule_sf <- st_as_sf(all_molecules, coords = c("x", "y"), crs = 32633)


# Perform spatial join to match molecules with cells
filtered_molecules_sf <- st_join(molecule_sf, cell_sf, left = FALSE)

# Convert back to a data frame (optional)
filtered_molecules1 <- as.data.frame(filtered_molecules_sf)

# Now 'filtered_molecules1' contains only the molecules within cell boundaries for Slide 1
saveRDS(filtered_molecules1, file = "ChP_filtered_molecules_slide1.89.rds")

# Now, to filter the Seurat object:
# First, retrieve the molecules that correspond to the filtered molecule coordinates

# Subset the Seurat object by selecting only the filtered molecules (those that are inside cell boundaries)
filtered_seurat <- readRDS("ChP_combined.rds")



test_gene <- "TTR"
test_coords <- filter_molecules_for_gene(test_gene, filtered_molecules1, ChP@images[["Slide.1.89"]])

if (!is.null(test_coords)) {
  cat("Filtered coordinates for gene:", test_gene, "\n")
  print(head(test_coords))
} else {
  cat("No coordinates matched for gene:", test_gene, "\n")
}


# Apply to both slides
ChP@images[["Slide.1.89"]] <- filter_molecules(filtered_molecules1, ChP@images[["Slide.1.89"]])

# Save the filtered Seurat object
saveRDS(ChP, file = "ChP_combined_filtered.rds")


# load cleaned Seurat object
ChP <- readRDS("data_objects/ChP_full.rds") # version with images

# define colors
major <- c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="white", "low-quality"="lightgrey")
sub <- c("Epi_1a"="#CD8500","Epi_1b"="#FFA500","Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382","Fib_2"="#27D7EF", "Fib_3"="#2a91ea","BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#6F845A", "Endo_1"="#c875c4", "Endo_2" = "#976697", "Endo_3" = "#c5aac5", "Mural_1"="#261132", "Mural_2"="#6d207b", "Mural_3"="#8f2aa2", "Parenchyma"="grey", "low-quality"="lightgrey")

DimPlot(ChP, reduction = "umap")

# remove transcripts from low-quality cells for visualization
fov <- "Slide.1.96" # repeat for specific FOVs, or all
gene_mols <- ChP@images[[fov]]@molecules[["molecules"]]

for (gene in names(gene_mols)) {
  mol_data <- gene_mols[[gene]]
  
  # Molecule coords rownames are numeric indices (as strings)
  meta_indices <- as.integer(rownames(mol_data@coords))
  
  # Subset metadata by numeric row index, NOT by rowname
  meta_sub <- ChP@meta.data[meta_indices, , drop = FALSE]
  
  # Keep only molecules from cells NOT labeled "low-quality"
  keep_idx <- meta_sub$Majorcelltype != "low-quality"
  
  # Subset coords & reset rownames accordingly (keep original numeric indices)
  mol_data@coords <- mol_data@coords[keep_idx, , drop = FALSE]
  rownames(mol_data@coords) <- rownames(meta_sub)[keep_idx]
  
  # Save back
  gene_mols[[gene]] <- mol_data
}

# Update FOV molecule data in the object
ChP@images[[fov]]@molecules[["molecules"]] <- gene_mols

# visualize FOV-level data with cells segmented (replace 'Slide.X.XX' with FOV of interest)

# ImageDimPlot of single FOV for Figure 1
DefaultBoundary(ChP@images$Slide.1.89) <- 'segmentation'
DefaultBoundary(ChP@images$Slide.1.96) <- 'segmentation'
ChP$Majorcelltype <- as.factor(ChP$Majorcelltype)
ChP@active.ident <- ChP$Majorcelltype

# major cell type on FOVs of interest (Figure 1)
ImageDimPlot(ChP, fov = "Slide.1.89", group.by = "Majorcelltype", cols=major, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2)
ImageDimPlot(ChP, fov = "Slide.1.96", group.by = "Majorcelltype", cols=major, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2)

# custom transcript colors
t_cols1 <- c("#E41A1CFF", "#377EB8FF", "#4DAF4AFF") # default colors for 3 transcripts
t_cols2 <- c("#E41A1CFF", "#FFFF33FF","#4DAF4AFF") # custom transcript colors, e.g. because fibroblasts are blue

### Transcript plots (Figure 1)
# Epithelial
ImageDimPlot(ChP, fov = "Slide.1.96", group.by = "Majorcelltype", cols=major, alpha = 0.8, molecules = c("TTR", "HTR2C", "GPX3"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)
ImageDimPlot(ChP, fov = "Slide.1.89", group.by = "Majorcelltype", cols=major, alpha = 0.8, molecules = c("TTR", "HTR2C", "GPX3"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)

# Fibroblasts
ImageDimPlot(ChP, fov = "Slide.1.96", cols=major, mols.cols = t_cols2, alpha = 0.8, group.by = "Majorcelltype", molecules = c("DCN", "FN1", "LEPR"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)
ImageDimPlot(ChP, fov = "Slide.1.89", cols=major, mols.cols = t_cols2, alpha = 0.8, group.by = "Majorcelltype", molecules = c("DCN", "FN1", "LEPR"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)

# Endothelial
ImageDimPlot(ChP, fov = "Slide.1.89", cols=major, alpha = 0.8, group.by = "Majorcelltype", molecules = c("VWF", "PECAM1", "INSR"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)
ImageDimPlot(ChP, fov = "Slide.1.96", cols=major, alpha = 0.8, group.by = "Majorcelltype", molecules = c("VWF", "PECAM1", "INSR"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)

# Mural
ImageDimPlot(ChP, fov = "Slide.1.89", cols=major, alpha = 0.8, group.by = "Majorcelltype", molecules = c("MYH11", "TAGLN", "ACTA2"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)
ImageDimPlot(ChP, fov = "Slide.1.96", cols=major, alpha = 0.8, group.by = "Majorcelltype", molecules = c("MYH11", "TAGLN", "ACTA2"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.31)

# Immune
ImageDimPlot(ChP, fov = "Slide.1.89", cols=major, alpha = 0.8, group.by = "Majorcelltype", molecules = c("STAB1", "CD163", "MSR1"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)
ImageDimPlot(ChP, fov = "Slide.1.96", cols=major, alpha = 0.8, group.by = "Majorcelltype", molecules = c("STAB1", "CD163", "MSR1"), nmols = 10000, dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2, mols.size=0.3)