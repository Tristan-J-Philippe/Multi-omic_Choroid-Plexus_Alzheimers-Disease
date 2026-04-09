###
### 050 ST_specific_plots_util.R
### Some of these functions were written by Dr. Denis R Avey.
# Purpose: Spatial molecular expression visualizations
# Dependencies:
library(Matrix)
library(data.table)
library(plyr)
library(dplyr)
library(GenomicRanges)
library(Seurat)
library(Signac)
library(future)
library(ggplot2)
library(cowplot)
library(reticulate)
library(GenomeInfoDb)
library(gridExtra)
library(patchwork)
library(tibble)
library(tidyr)
library(sf)
set.seed(123)
options(future.globals.maxSize = 10*1024^3)

#set cols
major_cols <- c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="brown", "low-quality"="lightgrey")
# Define a function to get projID based on FOV name
get_projID <- function(FOV) {
  if (startsWith(FOV, "Slide.1.")) {
    fov_num <- as.numeric(sub("Slide.1\\.", "", FOV))
    if (fov_num >= 1 && fov_num <= 46) {
      return("35395295")
    } else if (fov_num >= 47 && fov_num <= 104) {
      return("52052115")
    } else if (fov_num >= 105 && fov_num <= 146) {
      return("67364573")
    } else if (fov_num >= 147 && fov_num <= 212) {
      return("49297367")
    }
  } else if (startsWith(FOV, "Slide.2.")) {
    fov_num <- as.numeric(sub("Slide.2\\.", "", FOV))
    if (fov_num >= 1 && fov_num <= 47) {
      return("96074993")
    } else if (fov_num >= 48 && fov_num <= 108) {
      return("08536593")
    } else if (fov_num >= 109 && fov_num <= 170) {
      return("07761236")
    } else if (fov_num >= 171 && fov_num <= 225) {
      return("98625132")
    }
  }
  return(NA) # Return NA if FOV doesn't match any range
}


# Function to identify FOVs meeting specified criteria, e.g. for a given cell type, specifiy minimum cell # per FOV, as well as total cell #
get_valid_FOVs <- function(ChP, gene_list, major_cell_type = "Epithelial",  major_cell_type2 = NULL,  major_cell_type3 = NULL, min_cells = 100, max_cells = 200, min_total_cells = 400, percent_in_cells=10) {
  valid_FOVs <- list()
  
  # Loop through unique projIDs
  for (proj in unique(ChP$projID)) {
    FOVs <- names(ChP@images)
    
    for (FOV in FOVs) {
      # Determine projID for the current FOV
      fov_projID <- get_projID(FOV)
      if (is.na(fov_projID) || fov_projID != proj) next
      
      cat("Checking FOV:", FOV, "ProjID:", fov_projID, "\n")
      
      # Get the total number of cells in the FOV
      total_cells <- length(ChP@images[[FOV]]@boundaries[["centroids"]]@cells)
      if (total_cells < min_total_cells) {
        cat("Skipping FOV:", FOV, "- Total cells less than", min_total_cells, "\n")
        next
      }
      
      # Filter cells of the specified major cell type within the FOV
      fov_cells <- ChP@images[[FOV]]@boundaries[["centroids"]]@cells
      major_cells <- which(ChP$Majorcelltype == major_cell_type & ChP$projID == proj)
      major_cells_in_fov <- intersect(fov_cells, names(major_cells))
      
      if (length(major_cells_in_fov) < min_cells) {
        cat("Skipping FOV:", FOV, "- Major cell type less than", min_cells, "\n")
        next}
      
      if(!is.null(major_cell_type2)){
        major_cells2 <- which(ChP$Majorcelltype == major_cell_type2 & ChP$projID == proj)
        major_cells_in_fov2 <- intersect(fov_cells, names(major_cells2))
        if (length(major_cells_in_fov2) < min_cells) {
          cat("Skipping FOV:", FOV, "- Major cell type 2 less than", min_cells, "\n")
          next}
      }
      
      if(!is.null(major_cell_type3)){
        major_cells3 <- which(ChP$Majorcelltype == major_cell_type3 & ChP$projID == proj)
        major_cells_in_fov3 <- intersect(fov_cells, names(major_cells3))
        if (length(major_cells_in_fov3) < min_cells) {
          cat("Skipping FOV:", FOV, "- Major cell type 3 less than", min_cells, "\n")
          next}
      }
      
      if(length(major_cells_in_fov) > max_cells){
        cat("Skipping FOV:", FOV, "- Major cell type more than", max_cells, "\n")
        next}
      
      
      # Check if gene list is detected in >10% of major cells in the FOV; adjust '0.1' as needed for genes with lower/higher expression
      gene_counts <- rowSums(ChP@assays$SCT@counts[gene_list, major_cells_in_fov, drop = FALSE] > 0)
      gene_detected <- sum(gene_counts > (percent_in_cells/100 * length(major_cells_in_fov)))
      if (gene_detected < length(gene_list)) {
        cat("Skipping FOV:", FOV, "- Not all genes detected in >", percent_in_cells, "% of major cells\n")
        next
      }
      
      # Add FOV to the valid list
      valid_FOVs <- append(valid_FOVs, FOV)
    }
  }
  
  return(valid_FOVs)
}


### Apply further filtering and select candidate FOV pairs ###
# function to help select FOV pairs
get_valid_FOVs_for_category <- function(ChP, category_genes, valid_FOVs, decreasing=T) {
  # Create a list to store the average expressions for each FOV, categorized by AD and NCI
  AD_FOVs_avg_expr <- list()
  NCI_FOVs_avg_expr <- list()
  
  # Loop through valid FOVs in the given category
  for (FOV in valid_FOVs) {
    # Ensure FOV is a valid name (as string) within ChP@images
    if (!FOV %in% names(ChP@images)) {
      cat("Warning: FOV", FOV, "not found in ChP@images\n")
      next  # Skip invalid FOVs
    }
    
    # Get the cell IDs for the current FOV
    fov_cells <- ChP@images[[FOV]]@boundaries[["centroids"]]@cells
    
    # Debugging: Check if fov_cells are valid
    if (length(fov_cells) == 0) {
      cat("Warning: No cells found for FOV:", FOV, "\n")
      next
    }
    
    # Check if all category genes are present in the data
    missing_genes <- setdiff(category_genes, rownames(ChP@assays$SCT@data))
    if (length(missing_genes) > 0) {
      cat("Warning: These genes are missing from the dataset:", missing_genes, "\n")
    }
    
    # Get the average expression of the category genes for all cells in the FOV
    valid_genes <- intersect(category_genes, rownames(ChP@assays$SCT@data))
    
    if (length(valid_genes) == 0) {
      cat("Warning: No valid genes to calculate average expression for FOV:", FOV, "\n")
      next
    }
    
    # Check if fov_cells are valid columns in the expression data
    valid_cells <- intersect(fov_cells, colnames(ChP@assays$SCT@data))
    if (length(valid_cells) == 0) {
      cat("Warning: No valid cells found for FOV:", FOV, "\n")
      next
    }
    
    # Now, calculate the average expression for the valid cells and genes
    avg_expr <- colMeans(ChP@assays$SCT@data[valid_genes, valid_cells, drop = FALSE])
    
    # Determine the AD status for each cell in the FOV
    projID <- ChP$projID[valid_cells]
    
    # Calculate avg expression of specified genes in each FOV
    if (all(projID == "52052115" | projID == "67364573" | projID == "98625132" | projID == "96074993")){
      AD_avg_expr <- mean(avg_expr, na.rm = TRUE)
      AD_FOVs_avg_expr[[FOV]] <- AD_avg_expr
    } else if (all(projID == "08536593" | projID == "07761236" | projID == "35395295" | projID == "49297367")) {
      NCI_avg_expr <- mean(avg_expr, na.rm = TRUE)
      NCI_FOVs_avg_expr[[FOV]] <- NCI_avg_expr
    } else {
      cat("Warning: FOV", FOV, "contains mixed AD/NCI cells\n")
    }
  }
  
  # Now, compare each AD FOV with each NCI FOV
  valid_FOV_pairs <- list()
  
  for (AD_FOV in names(AD_FOVs_avg_expr)) {
    for (NCI_FOV in names(NCI_FOVs_avg_expr)) {
      # Get the average expression for AD and NCI FOVs
      AD_avg_expr <- AD_FOVs_avg_expr[[AD_FOV]]
      NCI_avg_expr <- NCI_FOVs_avg_expr[[NCI_FOV]]
      
      # Perform comparison only if both values are valid (not NA)
      if (!is.na(AD_avg_expr) && !is.na(NCI_avg_expr)) {
        # Calculate the absolute difference
        diff_expr <- AD_avg_expr - NCI_avg_expr
        
        # Add the pair with its difference to the list
        valid_FOV_pairs[[paste(AD_FOV, NCI_FOV, sep = "_vs_")]] <- list(
          AD_FOV = AD_avg_expr, 
          NCI_FOV = NCI_avg_expr, 
          Difference = diff_expr
        )
      }
    }
  }
  
  # Sort the pairs by the magnitude of the difference
  sorted_FOV_pairs <- valid_FOV_pairs[order(sapply(valid_FOV_pairs, function(x) x$Difference), decreasing = decreasing)]
  
  sorted_FOV_pairs <- data.frame(do.call(rbind, sorted_FOV_pairs))
  sorted_FOV_pairs$ad_fov_names <- sapply(strsplit(rownames(sorted_FOV_pairs), "_vs_"), function(x) x[1])  # First FOV (AD)
  sorted_FOV_pairs$nci_fov_names <- sapply(strsplit(rownames(sorted_FOV_pairs), "_vs_"), function(x) x[2])  # Second FOV (NCI)
  # Return the sorted pairs
  return(sorted_FOV_pairs)
}

### Different way of plotting for given gene/FOV
scaled_imagedimplot <- function(srt, gene, fov1, fov2, major_cell_type){
  #Makes Image dimplots scaled for two FOVs
  #srt: seurat object
  #gene of interest
  #chosen fov1
  #chosen fov2
  
  # Subset to cells in the second FOV (FOV 1)
  DefaultAssay(srt) <- "SCT"
  expr_values_fov1 <- FetchData(
    object = srt,
    vars = gene,
    layer = "data",  # Use "data" for SCT or "counts" for Nanostring
    cells = srt@images[[fov1]]@boundaries[["centroids"]]@cells
  )
  
  # Subset to cells in the second FOV (FOV 2)
  expr_values_fov2 <- FetchData(
    object = srt,
    vars = gene,
    layer = "data",  # Use "data" for SCT or "counts" for Nanostring
    cells = srt@images$Slide.1.183@boundaries$centroids@cells
  )
  
  combined_expr_values <- FetchData(object = srt, vars = gene, layer = "data")
  
  mdat <- data.frame(srt@meta.data)
  mdat <- mdat[mdat$Majorcelltype==major_cell_type,]
  
  # Drop NA values for both FOVs (not necessary if using 'Nanostring' as assay, but SCT scaled exp is more reasonable)
  combined_expr_values <- combined_expr_values[!is.na(combined_expr_values[[gene]]), , drop = FALSE]
  expr_values_fov1 <- expr_values_fov1[rownames(expr_values_fov1)%in%rownames(mdat), , drop = FALSE]
  expr_values_fov1 <- expr_values_fov1[!is.na(expr_values_fov1[[gene]]), , drop = FALSE]
  expr_values_fov2 <- expr_values_fov2[rownames(expr_values_fov2)%in%rownames(mdat), , drop = FALSE]
  expr_values_fov2 <- expr_values_fov2[!is.na(expr_values_fov2[[gene]]), , drop = FALSE]
  
  # Combine expression values from both FOVs to get global min and max
  #combined_expr_values <- rbind(expr_values_fov1, expr_values_fov2)
  
  # Get the global min and max
  global_min <- min(combined_expr_values[1])
  global_max <- max(combined_expr_values[1])
  
  # Scale expression values for FOV 1 based on global range
  expr_values_fov1_scaled <- (expr_values_fov1[[gene]] - global_min) / (global_max - global_min)
  expr_values_fov2_scaled <- (expr_values_fov2[[gene]] - global_min) / (global_max - global_min)
  combined_expr_values_scaled <- (combined_expr_values[[gene]] - global_min) / (global_max - global_min)
  
  # Sort scaled expression values and factorize them
  unique_scaled_values <- sort(unique(c(expr_values_fov1_scaled, expr_values_fov2_scaled)))
  
  # Factorize scaled values for FOV 1
  expr_fov1_factorized <- factor(expr_values_fov1_scaled, levels = unique_scaled_values)
  
  # Factorize scaled values for FOV 2
  expr_fov2_factorized <- factor(expr_values_fov2_scaled, levels = unique_scaled_values)
  
  # Generate a smooth color gradient
  library(viridis)
  color_gradient <- inferno(length(unique_scaled_values))  # Unique color for each factor level
  
  # Map factor levels to colors from the color pallete
  color_mapping_fov1 <- color_gradient[as.integer(expr_fov1_factorized)]
  color_mapping_fov2 <- color_gradient[as.integer(expr_fov2_factorized)]
  
  # Add factorized expression values as metadata for both FOVs
  metadata_fov1 <- data.frame(expr_fov1_factorized, row.names = rownames(expr_values_fov1))
  metadata_fov2 <- data.frame(expr_fov2_factorized, row.names = rownames(expr_values_fov2))
  
  # Add metadata to Seurat object
  srt <- AddMetaData(srt, metadata_fov1, col.name = paste0(fov1, "_", gene, "_scaled_SCT"))
  srt@meta.data[[paste0(fov1, "_", gene, "_scaled_SCT")]][which(srt$Majorcelltype != major_cell_type)]  <- NA
  srt <- AddMetaData(srt, metadata_fov2, col.name = paste0(fov2, "_", gene, "_scaled_SCT"))
  srt@meta.data[[paste0(fov2, "_", gene, "_scaled_SCT")]][which(srt$Majorcelltype != major_cell_type)]  <- NA
  
  ## Plot
  # Plot for FOV 1
  DefaultBoundary(srt@images[[fov1]]) <- 'segmentation'
  pl1 <- ImageDimPlot(
    object = srt,
    fov = fov1,
    group.by = paste0(fov1, "_", gene, "_scaled_SCT"),  # Use the factorized expression for color
    cols = color_gradient,  # Map to colors based on factorized values
    size = 0.5,
    alpha = 1,
    dark.background = FALSE,
    border.color = "black",
    coord.fixed = FALSE,
    border.size = 0.2
  )
  #Fix coordinates because of weird bug
  pl <- ImageDimPlot(
    object = srt,
    fov = fov1,
    group.by = "Majorcelltype",  # Use the factorized expression for color
    cols = color_gradient,  # Map to colors based on factorized values
    size = 0.5,
    alpha = 1,
    dark.background = FALSE,
    border.color = "black",
    coord.fixed = FALSE,
    border.size = 0.2
  )
  pl1$data$x <- pl$data$x
  pl1$data$y <- pl$data$y
  
  # Plot for FOV 2
  DefaultBoundary(srt@images$Slide.1.183) <- 'segmentation'
  pl2 <- ImageDimPlot(
    object = srt,
    fov = "Slide.1.183",
    group.by = paste0(fov2, "_", gene, "_scaled_SCT"),  # Use the factorized expression for color
    cols = color_gradient,  # Map to colors based on factorized values
    size = 0.5,
    alpha = 1,
    dark.background = FALSE,
    border.color = "black",
    coord.fixed = FALSE,
    border.size = 0.2
  )
  #Fix coordinates because of weird bug
  pl <- ImageDimPlot(
    object = srt,
    fov = "Slide.1.183",
    group.by = "Majorcelltype",  # Use the factorized expression for color
    cols = color_gradient,  # Map to colors based on factorized values
    size = 0.5,
    alpha = 1,
    dark.background = FALSE,
    border.color = "black",
    coord.fixed = FALSE,
    border.size = 0.2
  )
  pl2$data$x <- pl$data$x
  pl2$data$y <- pl$data$y
  
  pdf(paste0(curr_dir, "/ST_images/", gene, fov1, fov2, ".pdf"))
  #plot
  plot(pl1)
  plot(pl2)
  dev.off()
  plot(pl1)
  plot(pl2)
}


filter_molecules_for_gene <- function(gene, molecule_data, seurat_image) {
  cat("\n--- Filtering for gene:", gene, "---\n")
  
  # Step 1: Extract molecules for this gene from the filtered data
  gene_molecules <- molecule_data[molecule_data$gene == gene, ]
  if (nrow(gene_molecules) == 0) {
    cat("No filtered molecules found for gene:", gene, "\n")
    return(NULL)
  }
  cat("Filtered molecule rows for gene:", gene, "->", nrow(gene_molecules), "\n")
  
  # Extract x and y coordinates from 'geometry' using sf::st_coordinates
  coords_matrix <- sf::st_coordinates(gene_molecules$geometry)
  gene_molecules$x <- coords_matrix[, 1]
  gene_molecules$y <- coords_matrix[, 2]
  
  # Step 2: Extract existing coordinates for this gene in the Seurat object
  gene_coords <- seurat_image@molecules[["molecules"]][[gene]]@coords
  
  if (is.null(gene_coords)) {
    cat("No existing coordinates in Seurat for gene:", gene, "\n")
    return(NULL)
  }
  if (!is.data.frame(gene_coords)) {
    cat("Converting gene_coords to a data frame.\n")
    gene_coords <- as.data.frame(gene_coords)
    colnames(gene_coords) <- c("x", "y")  # Ensure column names are correct
  }
  
  cat("Existing Seurat coordinate rows for gene:", gene, "->", nrow(gene_coords), "\n")
  
  # Step 3: Match based on paired x-y coordinates
  molecule_xy <- paste(gene_molecules$x, gene_molecules$y, sep = "_")
  coord_xy <- paste(gene_coords$x, gene_coords$y, sep = "_")
  
  matched_indices <- which(coord_xy %in% molecule_xy)
  
  if (length(matched_indices) == 0) {
    cat("No matching coordinates found for gene:", gene, "\n")
    return(NULL)
  }
  
  filtered_coords <- gene_coords[matched_indices, , drop = FALSE]
  cat("Number of filtered coordinates for gene:", gene, "->", nrow(filtered_coords), "\n")
  
  return(filtered_coords)
}

filter_molecules <- function(molecule_data, seurat_image) {
  for (gene in names(seurat_image@molecules[["molecules"]])) {
    cat("\nProcessing gene:", gene, "\n")
    filtered_coords <- filter_molecules_for_gene(gene, molecule_data, seurat_image)
    
    # Convert to matrix if not NULL
    if (!is.null(filtered_coords)) {
      filtered_coords <- as.matrix(filtered_coords)
      seurat_image@molecules[["molecules"]][[gene]]@coords <- filtered_coords
    } else {
      # If no matching coordinates, remove the molecule entry
      seurat_image@molecules[["molecules"]][[gene]] <- NULL
      cat("Removed molecule entry for gene:", gene, "\n")
    }
  }
  return(seurat_image)
}
