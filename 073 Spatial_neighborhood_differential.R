###
### 073 Spatial_neighborhood_differential.R
### This script was written by Dr. Denis R Avey
# Purpose: Calculate spatial proximity and distance-based features between cell types.
# Dependencies:
source("/070 Spatial_neighborhood_util", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder


# load object
ChP_filt_nano <- readRDS("ChP_filt_nano.rds") # small object, without images or SCT assay

# to run AD-differential within niche, follow the next 4 steps; for niche-enrichment, skip these steps
# 1. subset seurat object by niche
ChP_filt_nano_1 <- subset(ChP_filt_nano, subset = niche_k3 == "1")
ChP_filt_nano_2 <- subset(ChP_filt_nano, subset = niche_k3 == "2")
ChP_filt_nano_3 <- subset(ChP_filt_nano, subset = niche_k3 == "3")

# add dimnames back to seurat objects; subsetting can cause them to become null
# 2. Get the gene and cell names
genes <- rownames(ChP_filt_nano)
cells1 <- colnames(ChP_filt_nano_1)
cells2 <- colnames(ChP_filt_nano_2)
cells3 <- colnames(ChP_filt_nano_3)

# 3. Assign gene and cell names to the counts layer
ChP_filt_nano_1@assays[["Nanostring"]]@layers[["counts"]]@Dimnames <- list(genes, cells1)
ChP_filt_nano_2@assays[["Nanostring"]]@layers[["counts"]]@Dimnames <- list(genes, cells2)
ChP_filt_nano_3@assays[["Nanostring"]]@layers[["counts"]]@Dimnames <- list(genes, cells3)


# Example usage - pairwise comparison of neighborhood clusters
result_C5_10_paired <- run_fgsea_single_cluster(
  seurat_obj = ChP_filt_nano,
  cluster_col = "niche_k3",
  gene_set_category = "C5"
)

# volcano plot

# Version with pairwise comparisons
df_pair <- read_csv("DE_results_pairwiseniche_k3_C5.csv") %>%
  mutate(comparison = as.character(comparison)) %>%
  filter(comparison == "1_vs_3")

# Identify top genes to label (you can change the criteria)
top_genes <- df_pair %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)  # top 10 by significance

gene <- c() #define manual list of genes to label

# Volcano plot with gene labels
ggplot(df_pair, aes(x = t, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  max.overlaps = 20,
                  size = 3,
                  box.padding = 0.3,
                  point.padding = 0.2) +
  theme_minimal() +
  labs(
    title = "Niche1 vs. Niche3 DEGs",
    x = "t-value (Niche1 vs. Niche3)",
    y = "-log10(adjusted P-value)",
    color = "adj.P.Val < 0.05"
  )

# visualization of fGSEA with barplot (Fig. 6h)

df <- read_csv("fGSEA_results_niche_k3_C5.csv")
filtered_terms <- c() # manually define gsea terms to include

# Step 1: Filter for positive NES and significant terms (optional: p.adj < 0.05)
results1 <- df[df$cluster=="1_vs_3",]
top_pathways1 <- results1 %>%
  filter(NES > 0, padj < 0.05) %>%
  filter(pathway %in% filtered_terms) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%  # optionally, filter by padj
  ungroup()

results3 <- df[df$cluster=="1_vs_3",]
top_pathways3 <- results3 %>%
  filter(NES < 0, padj < 0.05) %>%
  filter(pathway %in% filtered_terms) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%  # optionally, filter by padj
  ungroup()

# repeat for other comparisons of interest
# Step 2: Plot
ggplot(top_pathways1, aes(x = fct_reorder(pathway, -log10(padj)), y = -log10(padj), fill = -log10(padj))) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Gene Set",
    y = "-log10(adj. p-value)",
    fill = "-log10(adj. p-value)",
    title = "Top Enriched GSEA Terms (Niche 1)"
  ) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  theme_minimal(base_size = 12)

ggplot(top_pathways3, aes(x = fct_reorder(pathway, -log10(padj)), y = -log10(padj), fill = -log10(padj))) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Gene Set",
    y = "-log10(adj. p-value)",
    fill = "-log10(adj. p-value)",
    title = "Top Enriched GSEA Terms (Niche 3)"
  ) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  theme_minimal(base_size = 12)


### boxplot of gene expression, grouped by niche, averaged by individual (Extended Data Figure 6)

ChP_filt_nano <- readRDS("ChP_filt_nano.rds")
cluster_col <- 'niche_k3'
meta_filt <- ChP_filt_nano@meta.data
assay_counts <- ChP_filt_nano@assays$Nanostring@layers$counts

# Filter for TTR only
TTR_counts <- assay_counts["TTR", , drop = FALSE]

# Build metadata frame for cell-level info
cell_metadata <- meta_filt %>%
  select(unique_cellID, arb_ID, all_of(cluster_col), AD_status = AD_status) %>%
  rename(cell_id = unique_cellID,
         individual = arb_ID,
         cluster = !!sym(cluster_col))

# Join expression with metadata
TTR_long <- as.data.frame(TTR_counts) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "cell_id", values_to = "count") %>%
  left_join(cell_metadata, by = "cell_id") %>%
  filter(!is.na(cluster))

# Average expression per cell (i.e. mean count per cell per individual/cluster)
TTR_avg <- TTR_long %>%
  group_by(cluster, individual) %>%
  summarize(avg_expr = mean(count), .groups = "drop")

# niche color palette
niche_col <- c('1'    = "#33a02c",   
               '2'    = "#1f78b4",
               '3'    = "#ffb000")

ggplot(TTR_avg, aes(x = cluster, y = avg_expr, fill = cluster)) +
  geom_boxplot(outlier.color = "red") +
  geom_jitter(width = 0.1, alpha = 0.6) +
  scale_fill_manual(values = niche_col) +
  labs(x = "Niche Cluster", y = "Avg TTR Counts per Cell", title = "TTR Expression by Cluster") +
  theme_minimal()

# stats - ANOVA and Tukey's
# Perform ANOVA
anova_model <- aov(avg_expr ~ cluster, data = TTR_avg)

# Tukey's Honest Significant Difference test
tukey_results <- TukeyHSD(anova_model)

# View all pairwise niche comparisons
print(tukey_results)

# Optionally convert to data frame for easier plotting or filtering
tukey_df <- as.data.frame(tukey_results$cluster) %>%
  rownames_to_column("comparison") %>%
  arrange(p.adj)

print(tukey_df)
