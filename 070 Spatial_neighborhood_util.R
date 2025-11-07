###
### 070 Spatial_neighborhood_util.R
### These functions were written by Dr. Denis R Avey.
# Purpose: Compute and cluster spatial neighborhoods based on neighbor-averaged expression.
# Dependencies:
library(devtools)
library(CellularNeighborhoods) # devtools::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space@Main", subdir = "_code/cellular-neighborhoods")
library(Seurat)
library(Matrix)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(tibble)
library(lme4)
library(emmeans)
library(edgeR)
library(limma)
library(fgsea)
library(msigdbr)
library(pheatmap)
library(harmony)
library(ggplot2)
library(rlang)
library(ggrepel)
library(forcats)
library(tidyverse)
set.seed(123)


# helper function to plot only specific subtype proportions/numbers
plot_subtype_proportions <- function(df, subtype_list, palette, title) {
  df %>%
    filter(Subcelltype_v2 %in% subtype_list) %>%
    ggplot(aes(x = niche_k3, fill = Subcelltype_v2)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = palette[names(palette) %in% subtype_list]) +
    theme_minimal() +
    labs(title = title,
         x = "Spatial Cluster",
         y = "Proportion",
         fill = "Sub Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
plot_subtype_numbers <- function(df, subtype_list, palette, title) {
  df %>%
    filter(Subcelltype_v2 %in% subtype_list) %>%
    ggplot(aes(x = niche_k3, fill = Subcelltype_v2)) +
    geom_bar(position = "stack") +
    scale_fill_manual(values = palette[names(palette) %in% subtype_list]) +
    theme_minimal() +
    labs(title = title,
         x = "Spatial Cluster",
         y = "Proportion",
         fill = "Sub Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# define function for generating pseudobulk, performing limma DEG, and running fGSEA; each cluster vs. all others
run_fgsea <- function(seurat_obj, cluster_col, assay = "Nanostring", gene_set_category) {
  print("Step 1: Retrieving gene sets from MSigDB")
  genesets <- msigdbr(species = "Homo sapiens", collection = gene_set_category) %>%
    split(x = .$gene_symbol, f = .$gs_name)
  
  print("Step 2: Filtering cells with non-NA cluster assignments")
  # Optionally, filter cells from Seurat object based on non-NA cluster values
  
  print("Step 3: Creating sample_info df from metadata")
  meta_filt <- seurat_obj@meta.data
  sample_info <- data.frame(
    cell_id = meta_filt$unique_cellID,
    sample_id = paste0(
      meta_filt[[cluster_col]], ".", 
      meta_filt$arb_ID
    ),
    cluster = as.character(meta_filt[[cluster_col]]),
    individual = meta_filt$arb_ID
  )
  
  print("Step 4: Creating pseudobulk matrix and filtering out genes detected in <5% of cells")
  pseudobulk_mat <- seurat_obj@assays$Nanostring@layers$counts
  threshold <- 0.05 * ncol(pseudobulk_mat)  # 5% of total cells
  pseudobulk_mat <- pseudobulk_mat[rowSums(pseudobulk_mat > 0) > threshold, ] # removes 134 (of 6175) genes with low expression
  
  print("Step 5: Aggregating counts by cluster and individual")
  # Step 5: Aggregating counts by cluster and individual
  pseudobulk_mat_cluster <- pseudobulk_mat %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "cell_id", values_to = "count") %>%
    left_join(sample_info, by = "cell_id") %>%
    group_by(gene, sample_id) %>%  # Now group by sample_id, which already includes cluster and individual
    summarize(count = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = sample_id, values_from = count)
  pseudobulk_mat_cluster <- pseudobulk_mat_cluster[, -1]  # Remove the 'gene' column
  rownames(pseudobulk_mat_cluster) <- rownames(pseudobulk_mat)
  
  print("Step 6: Create design matrix")
  # step 6
  # Create a data frame of experimental factors for each sample
  sample_info_design <- data.frame(
    Sample = colnames(pseudobulk_mat_cluster),
    # Extract Cluster (first number before the first dot)
    Cluster = sub("^([0-9]+)\\..*", "\\1", colnames(pseudobulk_mat_cluster)),
    # Extract Individual (the text after the second dot)
    Individual = sub("^[0-9]+\\.(.*)", "\\1", colnames(pseudobulk_mat_cluster))
  )
  
  # Assuming 'sample_info_design' is your original data frame
  # Convert categorical variables to factors (this step is correct)
  sample_info_design$Cluster <- as.factor(sample_info_design$Cluster)
  sample_info_design$Individual <- as.factor(sample_info_design$Individual)
  
  # Create the design matrix but remove the intercept to include all levels
  design <- model.matrix(~ 0 + Cluster + Individual, data = sample_info_design) # niche-enriched, accounting for individual
  
  # Retain the 'Sample' column as row names for the design matrix
  rownames(design) <- sample_info_design$Sample
  
  # stop if design and pseudobulk matrices are misaligned
  print("Checking pseudobulk/design alignment")
  stopifnot(all(colnames(pseudobulk_mat_cluster) == rownames(design)))
  
  print("Step 7: DGEList and voom for DE analysis")
  # DGEList and voom
  dge <- DGEList(pseudobulk_mat_cluster)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  
  print("Step 8: Running fGSEA")
  # Run fGSEA
  fgsea_results_all <- list()
  de_results <- list()
  
  clusters <- unique(sample_info$cluster)
  
  for (clust in clusters) {
    print(paste("Analyzing cluster:", clust))
    
    cluster_cols <- grep("^Cluster", colnames(design), value = TRUE)
    target <- paste0("Cluster", clust)
    others <- setdiff(cluster_cols, target)
    
    if (length(others) == 0) {
      stop("No other clusters to compare against for contrast.")
    }
    
    # Construct contrast string
    contrast_formula <- paste0(target, " - (", paste(others, collapse = " + "), ")/", length(others))
    contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # Differential expression table for cluster
    de_table <- topTable(fit2, number = Inf, sort.by = "none")
    de_table$gene <- rownames(de_table)
    de_table$cluster <- clust
    de_results[[clust]] <- de_table
    
    # Gene rankings for fGSEA
    ranks <- setNames(sign(de_table$logFC) * -log10(de_table$adj.P.Val + 1e-10), de_table$gene)
    ranks <- ranks[!is.na(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    
    # fGSEA
    set.seed(123)
    # Using fgseaMultilevel instead of fgseaSimple
    fgsea_out <- fgseaMultilevel(pathways = genesets, stats = ranks)
    fgsea_out$cluster <- clust
    fgsea_results_all[[clust]] <- fgsea_out
  }
  
  print("Step 9: Saving results")
  # Save results
  de_results_combined <- bind_rows(de_results)
  # Modify the file name to include cluster_col and gene_set_category
  de_file_name <- paste0("DE_results_", cluster_col, "_", gene_set_category, ".csv") 
  write.csv(de_results_combined, de_file_name, row.names = FALSE)
  
  fgsea_results_combined <- bind_rows(fgsea_results_all)
  # Flatten the 'leadingEdge' column to a single string per row
  fgsea_results_combined$leadingEdge <- sapply(fgsea_results_combined$leadingEdge, function(x) paste(x, collapse = ", "))
  # Now save the data as a CSV
  # Modify the file name to include cluster_col and gene_set_category
  fgsea_file_name <- paste0("fGSEA_results_", cluster_col, "_", gene_set_category, ".csv")
  write.csv(fgsea_results_combined, fgsea_file_name, row.names = FALSE)
  
  return(list(
    fgsea = fgsea_results_all,
    de_results = de_results
  ))
}

# version with pairwise cluster comparisons (used for Figure 6)
run_fgsea_single_cluster <- function(seurat_obj, cluster_col, assay = "Nanostring", gene_set_category) {
  print("Step 1: Retrieving gene sets from MSigDB")
  genesets <- msigdbr(species = "Homo sapiens", collection = gene_set_category) %>%
    split(x = .$gene_symbol, f = .$gs_name)
  
  print("Step 2: Filtering cells with non-NA cluster assignments")
  # Optionally, filter cells from Seurat object based on non-NA cluster values
  
  print("Step 3: Creating sample_info df from metadata")
  meta_filt <- seurat_obj@meta.data
  sample_info <- data.frame(
    cell_id = meta_filt$unique_cellID,
    sample_id = paste0(
      meta_filt[[cluster_col]], ".", 
      meta_filt$arb_ID
    ),
    cluster = as.character(meta_filt[[cluster_col]]),
    individual = meta_filt$arb_ID
  )
  
  print("Step 4: Creating pseudobulk matrix and filtering out genes detected in <5% of cells")
  # Assuming ChP_filt has a 'counts' layer in the 'Nanostring' assay
  pseudobulk_mat <- seurat_obj@assays$Nanostring@layers$counts
  threshold <- 0.05 * ncol(pseudobulk_mat)  # 5% of total cells
  pseudobulk_mat <- pseudobulk_mat[rowSums(pseudobulk_mat > 0) > threshold, ] # removes 134 (of 6175) genes with low expression
  
  print("Step 5: Aggregating counts by cluster and individual")
  # Step 5: Aggregating counts by cluster and individual
  pseudobulk_mat_cluster <- pseudobulk_mat %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "cell_id", values_to = "count") %>%
    left_join(sample_info, by = "cell_id") %>%
    group_by(gene, sample_id) %>%  # Now group by sample_id, which already includes cluster and individual
    summarize(count = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = sample_id, values_from = count)
  pseudobulk_mat_cluster <- pseudobulk_mat_cluster[, -1]  # Remove the 'gene' column
  rownames(pseudobulk_mat_cluster) <- rownames(pseudobulk_mat)
  
  print("Step 6: Create design matrix")
  # step 6
  # Create a data frame of experimental factors for each sample
  sample_info_design <- data.frame(
    Sample = colnames(pseudobulk_mat_cluster),
    # Extract Cluster (first number before the first dot)
    Cluster = sub("^([0-9]+)\\..*", "\\1", colnames(pseudobulk_mat_cluster)),
    # Extract Individual (the text after the second dot)
    Individual = sub("^[0-9]+\\.(.*)", "\\1", colnames(pseudobulk_mat_cluster))
  )
  
  # Assuming 'sample_info_design' is your original data frame
  # Convert categorical variables to factors (this step is correct)
  sample_info_design$Cluster <- as.factor(sample_info_design$Cluster)
  sample_info_design$Individual <- as.factor(sample_info_design$Individual)
  
  # Create the design matrix but remove the intercept to include all levels
  design <- model.matrix(~ 0 + Cluster + Individual, data = sample_info_design) # niche-enriched, accounting for individual
  
  # Retain the 'Sample' column as row names for the design matrix
  rownames(design) <- sample_info_design$Sample
  
  # stop if design and pseudobulk matrices are misaligned
  print("Checking pseudobulk/design alignment")
  stopifnot(all(colnames(pseudobulk_mat_cluster) == rownames(design)))
  
  print("Step 7: DGEList and voom for DE analysis")
  # DGEList and voom
  dge <- DGEList(pseudobulk_mat_cluster)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  
  print("Step 8: Running fGSEA")
  # Step 8: Run DE and fGSEA for each pairwise cluster comparison
  print("Step 8: Running pairwise fGSEA")
  cluster_cols <- grep("^Cluster", colnames(design), value = TRUE)
  cluster_levels <- gsub("^Cluster", "", cluster_cols)
  pairwise_comparisons <- combn(cluster_levels, 2, simplify = FALSE)
  
  fgsea_results_all <- list()
  de_results <- list()
  
  for (pair in pairwise_comparisons) {
    cluster1 <- pair[1]
    cluster2 <- pair[2]
    print(paste("Analyzing pair:", cluster1, "vs", cluster2))
    
    contrast_formula <- paste0("Cluster", cluster1, " - Cluster", cluster2)
    contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # Differential expression table for the pair
    de_table <- topTable(fit2, number = Inf, sort.by = "none")
    de_table$gene <- rownames(de_table)
    de_table$comparison <- paste0(cluster1, "_vs_", cluster2)
    de_results[[paste0(cluster1, "_vs_", cluster2)]] <- de_table
    
    # Gene rankings for fGSEA
    ranks <- setNames(sign(de_table$logFC) * -log10(de_table$adj.P.Val + 1e-10), de_table$gene)
    ranks <- ranks[!is.na(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    
    set.seed(123)
    fgsea_out <- fgseaMultilevel(pathways = genesets, stats = ranks)
    fgsea_out$comparison <- paste0(cluster1, "_vs_", cluster2)
    fgsea_results_all[[paste0(cluster1, "_vs_", cluster2)]] <- fgsea_out
  }
  
  print("Step 9: Saving results")
  # Save results
  de_results_combined <- bind_rows(de_results, .id = "comparison")
  # Modify the file name to include cluster_col and gene_set_category
  de_file_name <- paste0("DE_results_pairwise", cluster_col, "_", gene_set_category, ".csv") 
  write.csv(de_results_combined, de_file_name, row.names = FALSE)
  
  # save fGSEA
  fgsea_results_combined <- bind_rows(fgsea_results_all, .id = "comparison")
  # Flatten the 'leadingEdge' column to a single string per row
  fgsea_results_combined$leadingEdge <- sapply(fgsea_results_combined$leadingEdge, function(x) paste(x, collapse = ", "))
  # Now save the data as a CSV
  # Modify the file name to include cluster_col and gene_set_category
  fgsea_file_name <- paste0("fGSEA_results_pairwise", cluster_col, "_", gene_set_category, ".csv")
  write.csv(fgsea_results_combined, fgsea_file_name, row.names = FALSE)
  
  return(list(
    fgsea = fgsea_results_all,
    de_results = de_results
  ))
}

