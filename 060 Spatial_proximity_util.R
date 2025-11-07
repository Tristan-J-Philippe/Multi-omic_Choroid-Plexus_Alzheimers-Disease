###
### 060 Spatial_proximity_util.R
### These functions were written by Dr. Denis R Avey.
# Purpose: Calculate spatial proximity and distance-based features between cell types.
# Dependencies:
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(purrr)
library(readr)
library(limma)
library(lme4)
library(emmeans)

# define function
get_closest_endo <- function(seu, mural_subtypes, endo_subtypes) {
  meta <- seu@meta.data
  results <- list()
  
  for (mural in mural_subtypes) {
    mural_cells <- meta %>% filter(Subcelltype_v2 == mural)
    
    res <- mural_cells %>%
      group_by(projID) %>%
      group_map(~{
        id <- .y$projID
        source_cells <- .x
        target_cells <- meta %>% filter(projID == id, Subcelltype_v2 %in% endo_subtypes)
        if (nrow(target_cells) == 0) return(NULL)
        
        source_coords <- source_cells[, c("x_um", "y_um")]
        target_coords <- target_cells[, c("x_um", "y_um")]
        target_labels <- target_cells$Subcelltype_v2
        
        # Find closest endothelial subtype for each mural cell
        closest <- apply(source_coords, 1, function(s) {
          dists <- sqrt((target_coords$x_um - s[1])^2 + (target_coords$y_um - s[2])^2)
          target_labels[which.min(dists)]
        })
        
        data.frame(projID = id,
                   mural_subtype = mural,
                   closest_endo = closest)
      }) %>% bind_rows()
    
    results[[mural]] <- res
  }
  
  bind_rows(results)
}


get_nearest_distance_to_celltype <- function(seu, type, target_types) {
  all_meta <- seu@meta.data
  meta <- all_meta %>% filter(newMajorcelltype == type)
  
  meta %>%
    group_by(individual_ID) %>%
    group_map(~{
      id <- .y$individual_ID
      source_cells <- .x
      target_cells <- all_meta %>%
        filter(individual_ID == id, newMajorcelltype %in% target_types)
      
      if (nrow(target_cells) == 0) {
        return(data.frame(individual_ID = id, avg_distance_um = NA, newMajorcelltype = type))
      }
      
      # Extract coordinates
      source_coords <- source_cells[, c("x_um", "y_um")]
      target_coords <- target_cells[, c("x_um", "y_um")]
      
      # Compute distances while avoiding self-distance (exact coord match)
      nearest_dists <- mapply(function(ix, iy) {
        dists <- sqrt((target_coords$x_um - ix)^2 + (target_coords$y_um - iy)^2)
        # Exclude any zero-distance matches (i.e., self or same spatial location)
        dists <- dists[dists > 0]
        if (length(dists) == 0) {
          return(NA)  # fallback if only self existed
        }
        min(dists)
      },
      ix = source_coords$x_um,
      iy = source_coords$y_um)
      
      data.frame(individual_ID = id, avg_distance_um = mean(nearest_dists, na.rm = TRUE), newMajorcelltype = type)
    }) %>%
    bind_rows()
}