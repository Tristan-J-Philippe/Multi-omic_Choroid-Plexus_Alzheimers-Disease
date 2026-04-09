###
### 072 Spatial_neighborhood_proportions.R
### This script was written by Dr. Denis R Avey
# Purpose: Calculate proportions of cells within spatial neighborhoods.
# Dependencies:
source("/070 Spatial_neighborhood_util", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder


# subcelltype color palette
sub <- c("Epi_1a"="#CD8500","Epi_1b"="#FFA500","Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382","Fib_2"="#27D7EF", "Fib_3"="#2a91ea","BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#6F845A", "Endo_1"="#c875c4", "Endo_2" = "#976697", "Endo_3" = "#c5aac5", "Mural_1"="#261132", "Mural_2"="#6d207b", "Mural_3"="#8f2aa2", "Parenchyma"="grey", "low-quality"="lightgrey")

# load CosMx Seurat object
ChP_filt <- readRDS("data_objects/ChP_filt_nano.rds")

# Extract metadata as a dataframe
metadata_df <- ChP_filt@meta.data

# Define cell types of interest
cell_types <- c("Endothelial", "Epithelial", "Fibroblast", "Mural", "Immune")
epi_sub <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b")
fib_sub <- c("Fib_1","Fib_2", "Fib_3")
immune_sub <- c("BAM_1", "BAM_2", "BAM_3", "T_Cell")
endo_sub <- c("Endo_1", "Endo_2", "Endo_3")
mural_sub <- c("Mural_1", "Mural_2", "Mural_3")
all_sub <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b","Fib_1","Fib_2", "Fib_3","BAM_1", "BAM_2", "BAM_3", "T_Cell","Endo_1", "Endo_2", "Endo_3", "Mural_1", "Mural_2", "Mural_3")

# Plot grouped by neighborhood cluster (Proportions)
ggplot(metadata_df, aes(x = niche_k3, fill = newMajorcelltype)) +
  geom_bar(position = "fill") +  # Normalize to proportions
  scale_fill_manual(values = major) +
  theme_minimal() +
  labs(title = "Proportion of Major Cell Types within Niches",
       x = "Niche",
       y = "Proportion",
       fill = "Major Cell Type (stringent)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot all subtypes grouped by neighborhood cluster (Proportions)
ggplot(metadata_df, aes(x = niche_k3, fill = Subcelltype_v2)) +
  geom_bar(position = "fill") +  # Normalize to proportions
  scale_fill_manual(values = sub) +
  theme_minimal() +
  labs(title = "Proportion of Sub-cell-types within Niches",
       x = "Spatial Cluster",
       y = "Proportion",
       fill = "Sub Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# helper function to plot only specific subtype proportions/numbers

# Example usage:
plot_subtype_proportions(metadata_df, epi_sub, sub, "Proportion of Epithelial Subtypes")
plot_subtype_proportions(metadata_df, fib_sub, sub, "Proportion of Fibroblast Subtypes")
plot_subtype_proportions(metadata_df, immune_sub, sub, "Proportion of Immune Subtypes")
plot_subtype_proportions(metadata_df, endo_sub, sub, "Proportion of Endothelial Subtypes")
plot_subtype_proportions(metadata_df, mural_sub, sub, "Proportion of Mural Subtypes")

plot_subtype_numbers(metadata_df, epi_sub, sub, "Number of Epithelial Subtypes")
plot_subtype_numbers(metadata_df, fib_sub, sub, "Number of Fibroblast Subtypes")
plot_subtype_numbers(metadata_df, immune_sub, sub, "Number of Immune Subtypes")
plot_subtype_numbers(metadata_df, endo_sub, sub, "Proportion of Endothelial Subtypes")
plot_subtype_numbers(metadata_df, mural_sub, sub, "Proportion of Mural Subtypes")


per_indiv_props <- metadata_df %>%
  group_by(individual_ID, AD_status, niche_k3, Subcelltype_v2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(individual_ID, AD_status, niche_k3) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Plot
ggplot(per_indiv_props, aes(x = niche_k3, y = prop, fill = Subcelltype_v2)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_jitter(aes(color = Subcelltype_v2),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = sub) +
  scale_color_manual(values = sub) +
  labs(y = "Proportion of cells per individual",
       x = "Niche (k=3)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## limma stats

# Load participant metadata
participant_meta <- read_csv("data_csv/CosMX-projID metadata.csv") # <-- Update path if needed

# major cell type level
results_list <- list()

# Loop over major cell types
for (celltype in cell_types) {
  
  # 1. Logical column for cell type presence
  metadata_df$target_type <- metadata_df$newMajorcelltype == celltype
  
  # 2. Proportion of target cell type per niche × individual
  prop_df <- metadata_df %>%
    group_by(projID, niche_k3) %>%
    summarise(
      prop_target = mean(target_type),
      .groups = "drop"
    )
  
  # 3. Global proportion per individual
  individual_baseline <- metadata_df %>%
    group_by(projID) %>%
    summarise(global_prop = mean(target_type), .groups = "drop")
  
  # 4. Merge and compute deviation
  prop_df <- left_join(prop_df, individual_baseline, by = "projID") %>%
    mutate(prop_deviation = prop_target - global_prop)
  
  # 5. Add metadata (AD_status, etc.)
  prop_df <- left_join(prop_df, participant_meta, by = c("projID"))
  
  # 6. Design matrix: niche_k3 and AD_status ('Cog_Path')
  # Keep only complete cases (to avoid NA issues in model matrix)
  design_df <- prop_df %>%
    select(prop_deviation, niche_k3, Cog_Path, projID, msex, age_death, pmi) %>%
    filter(complete.cases(.))
  
  design_df$niche_k3 <- factor(design_df$niche_k3)
  design_df$Cog_Path <- factor(design_df$Cog_Path, levels = c("NCI", "AD"))
  
  design_df <- design_df %>%
    mutate(sample = paste0("s_", projID, "_", niche_k3))
  
  # Optional: include covariates, but may be too much for n = 8
  # Assign sample as rownames before model matrix
  rownames(design_df) <- design_df$sample
  design <- model.matrix(~ 0 + niche_k3 + Cog_Path + msex, data = design_df) # AD status and sex; version used for figures
  
  # 7. Format proportions as a matrix (genes = niches, samples = individuals)
  # limma expects rows = features, columns = samples
  prop_matrix <- design_df %>%
    select(sample, prop_deviation) %>%
    column_to_rownames("sample") %>%
    as.matrix() %>%
    t()  # matrix with 1 row (cell type), columns = samples
  
  # 8. Apply limma
  fit <- lmFit(prop_matrix, design)
  
  # 9. Contrast: deviation from global for each niche
  contrast_matrix <- makeContrasts(
    niche1_vs_baseline = niche_k31,
    niche2_vs_baseline = niche_k32,
    niche3_vs_baseline = niche_k33,
    levels = design
  )
  
  # 10. Fit contrasts and get p-values
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # 11. Store one row per contrast with p-values
  results_list[[celltype]] <- data.frame(
    contrast = colnames(contrast_matrix),
    t = as.numeric(fit2$t),
    p.value = as.numeric(fit2$p.value),
    adj.p.value = p.adjust(as.numeric(fit2$p.value), method = "BH")
  )
  # 11B. Save topTable - ANOVA-style p-value (overall)
  #results_list[[celltype]] <- result
}

# sub cell type level
results_list <- list()

# Loop over sub cell types (e.g. BAM_1 / all subtypes)
for (celltype in all_sub) {
  
  # 1. Logical column for cell type presence
  metadata_df$target_type <- metadata_df$Subcelltype_v2 == celltype
  
  # 2. Proportion of target cell type per niche × individual
  prop_df <- metadata_df %>%
    group_by(projID, niche_k3) %>%
    summarise(
      prop_target = mean(target_type),
      .groups = "drop"
    )
  
  # 3. Global proportion per individual
  individual_baseline <- metadata_df %>%
    group_by(projID) %>%
    summarise(global_prop = mean(target_type), .groups = "drop")
  
  # 4. Merge and compute deviation
  prop_df <- left_join(prop_df, individual_baseline, by = "projID") %>%
    mutate(prop_deviation = prop_target - global_prop)
  
  # 5. Add metadata (AD_status, pmi, age_death, etc.)
  prop_df <- left_join(prop_df, participant_meta, by = c("projID"))
  
  # 6. Design matrix: niche_k3 and AD_status ('Cog_Path')
  # Keep only complete cases (to avoid NA issues in model matrix)
  design_df <- prop_df %>%
    select(prop_deviation, niche_k3, Cog_Path, projID, msex, age_death, pmi) %>%
    filter(complete.cases(.))
  
  design_df$niche_k3 <- factor(design_df$niche_k3)
  design_df$Cog_Path <- factor(design_df$Cog_Path, levels = c("NCI", "AD"))
  
  design_df <- design_df %>%
    mutate(sample = paste0("s_", projID, "_", niche_k3))
  
  # Optional: include covariates, but may be too much for n = 8
  # Assign sample as rownames before model matrix
  rownames(design_df) <- design_df$sample
  design <- model.matrix(~ 0 + niche_k3 + Cog_Path + msex, data = design_df) # AD status and sex; version used for figures
  
  # 7. Format proportions as a matrix (genes = niches, samples = individuals)
  # limma expects rows = features, columns = samples
  prop_matrix <- design_df %>%
    select(sample, prop_deviation) %>%
    column_to_rownames("sample") %>%
    as.matrix() %>%
    t()  # matrix with 1 row (cell type), columns = samples
  
  # 8. Apply limma
  fit <- lmFit(prop_matrix, design)
  
  # 9. Contrast: deviation from global for each niche
  contrast_matrix <- makeContrasts(
    niche1_vs_baseline = niche_k31,
    niche2_vs_baseline = niche_k32,
    niche3_vs_baseline = niche_k33,
    levels = design
  )
  
  # 10. Fit contrasts and get p-values
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # 11. Store one row per contrast with p-values
  results_list[[celltype]] <- data.frame(
    contrast = colnames(contrast_matrix),
    t = as.numeric(fit2$t),
    p.value = as.numeric(fit2$p.value),
    adj.p.value = p.adjust(as.numeric(fit2$p.value), method = "BH")
  )
  # 11B. Save topTable - ANOVA-style p-value (overall)
  #results_list[[celltype]] <- result
}

# version restricted to specific subtypes (e.g. BAM_1 / all immune subtypes)
results_list <- list()

# Loop over sub cell types
for (celltype in epi_sub) {
  
  # 1. Define the subset of subtypes to consider (e.g., epi_sub for epithelial)
  subset_df <- metadata_df %>% filter(Subcelltype_v2 %in% epi_sub)  # replace with the subset of interest
  
  # 2. Logical column for cell type presence within subset
  subset_df$target_type <- subset_df$Subcelltype_v2 == celltype
  
  # 3. Proportion of target cell type per niche × individual (subset denominator)
  prop_df <- subset_df %>%
    group_by(projID, niche_k3) %>%
    summarise(prop_target = mean(target_type), .groups = "drop")
  
  # 4. Baseline proportion per individual (within same subset)
  individual_baseline <- subset_df %>%
    group_by(projID) %>%
    summarise(global_prop = mean(target_type), .groups = "drop")
  
  # 5. Merge and compute deviation
  prop_df <- left_join(prop_df, individual_baseline, by = "projID") %>%
    mutate(prop_deviation = prop_target - global_prop)
  
  # 6. Add metadata
  prop_df <- left_join(prop_df, participant_meta, by = c("projID"))
  
  # 7. Design matrix
  design_df <- prop_df %>%
    select(prop_deviation, niche_k3, Cog_Path, projID, msex, age_death, pmi) %>%
    filter(complete.cases(.))
  
  design_df$niche_k3 <- factor(design_df$niche_k3)
  design_df$Cog_Path <- factor(design_df$Cog_Path, levels = c("NCI", "AD"))
  
  design_df <- design_df %>%
    mutate(sample = paste0("s_", projID, "_", niche_k3))
  
  rownames(design_df) <- design_df$sample
  design <- model.matrix(~ 0 + niche_k3 + Cog_Path + msex, data = design_df)
  
  # 8. Format matrix for limma
  prop_matrix <- design_df %>%
    select(sample, prop_deviation) %>%
    column_to_rownames("sample") %>%
    as.matrix() %>%
    t()
  
  # 9. Apply limma
  fit <- lmFit(prop_matrix, design)
  
  # 10. Contrasts
  contrast_matrix <- makeContrasts(
    niche1_vs_baseline = niche_k31,
    niche2_vs_baseline = niche_k32,
    niche3_vs_baseline = niche_k33,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # 11. Store results
  results_list[[celltype]] <- data.frame(
    contrast = colnames(contrast_matrix),
    t = as.numeric(fit2$t),
    p.value = as.numeric(fit2$p.value),
    adj.p.value = p.adjust(as.numeric(fit2$p.value), method = "BH")
  )
} # edit subtype group in step 1 of function code (e.g. 'epi_sub') as well
