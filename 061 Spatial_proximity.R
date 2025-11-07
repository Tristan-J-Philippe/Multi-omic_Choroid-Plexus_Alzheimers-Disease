###
### 061 Spatial_proximity.R
### This script was written by Dr. Denis R Avey
# Purpose: Calculate spatial proximity and distance-based features between cell types.
# Dependencies:
source("/060 Spatial_proximity_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder

# custom color palettes
sub <- c("Epi_1a"="#CD8500","Epi_1b"="#FFA500","Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382","Fib_2"="#27D7EF", "Fib_3"="#2a91ea","BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#6F845A", "Endo_1"="#c875c4", "Endo_2" = "#976697", "Endo_3" = "#c5aac5", "Mural_1"="#261132", "Mural_2"="#6d207b", "Mural_3"="#8f2aa2", "Parenchyma"="grey", "low-quality"="lightgrey")
major <- c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B", "Parenchyma"="grey", "low-quality" = "lightgrey")
niche_col <- c(
  "1"    = "#33a02c",   
  "2"    = "#1f78b4",
  "3"    = "#ffb000"
)

### nearest Endothelial subtype among Mural subtypes (Fig. 5e)
# load seurat object
ChP_filt <- readRDS("data_objects/ChP_filt_082025.rds")

# Step 1: Define subtypes
endo_subtypes <- c("Endo_1", "Endo_2", "Endo_3")
mural_subtypes <- c("Mural_1", "Mural_2", "Mural_3")

# Step 2: Run function
closest_endo_df <- get_closest_endo(ChP_filt, mural_subtypes, endo_subtypes)

# Step 3: Aggregate proportions per mural subtype
prop_df <- closest_endo_df %>%
  group_by(mural_subtype, closest_endo) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mural_subtype) %>%
  mutate(prop = n / sum(n))

# Step 4: Plot
ggplot(prop_df, aes(x = mural_subtype, y = prop, fill = closest_endo)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Endo_1"="#c875c4", "Endo_2"="#976697", "Endo_3"="#c5aac5")) +
  labs(x = "Mural subtype", y = "Proportion of cells with closest endothelial", fill = "Closest endothelial") +
  theme_minimal()

# stats
# Compute proportions per individual
indiv_prop_df <- closest_endo_df %>%
  group_by(projID = projID, mural_subtype, closest_endo) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(projID, mural_subtype) %>%
  mutate(prop_target = n_cells / sum(n_cells)) %>%  # proportion within mural subtype
  ungroup()

# Compute individual-level global null proportion for each endothelial subtype
indiv_baseline <- closest_endo_df %>%
  group_by(projID = projID, closest_endo) %>%
  summarise(global_prop = n() / sum(nrow(closest_endo_df[closest_endo_df$projID==projID,])),
            .groups = "drop")

# Merge with per-mural-type proportions
prop_df <- indiv_prop_df %>%
  left_join(indiv_baseline, by = c("projID", "closest_endo")) %>%
  mutate(prop_deviation = prop_target - global_prop)

# Load participant metadata
participant_meta <- read_csv("data_csv/CosMX-projID metadata.csv") # <-- Update path if needed

prop_df <- prop_df %>%
  left_join(participant_meta, by = "projID")

# Example: model with mural subtype, AD status, and sex
design_df <- prop_df %>%
  select(prop_deviation, mural_subtype, Cog_Path, projID, msex, closest_endo) %>%
  filter(complete.cases(.))

design_df <- design_df %>%
  mutate(sample = paste0("s_", projID, "_", mural_subtype, "_", closest_endo))
rownames(design_df) <- design_df$sample

design <- model.matrix(~ 0 + mural_subtype + Cog_Path + msex, data = design_df)

# Format prop_deviation as matrix (features = mural×endothelial pairs, samples = individuals)
prop_matrix <- design_df %>%
  select(sample, prop_deviation) %>%
  column_to_rownames("sample") %>%
  as.matrix() %>%
  t()

# Fit limma
fit <- lmFit(prop_matrix, design)

contrast_matrix <- makeContrasts(
  mural1_vs_baseline = mural_subtypeEndo_1,
  mural2_vs_baseline = mural_subtypeEndo_2,
  mural3_vs_baseline = mural_subtypeEndo_3,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number = Inf, sort.by = "none")

# binomial tests
# Step 1: compute global proportions of each endothelial subtype
global_props <- closest_endo_df %>%
  group_by(closest_endo) %>%
  summarise(global_prop = n() / nrow(closest_endo_df), .groups = "drop")

# Step 2: compute counts per mural subtype
mural_counts <- closest_endo_df %>%
  group_by(mural_subtype, closest_endo) %>%
  summarise(n_success = n(), .groups = "drop") %>%
  # compute total number of cells per mural subtype
  group_by(mural_subtype) %>%
  mutate(n_total = sum(n_success)) %>%
  ungroup() %>%
  left_join(global_props, by = "closest_endo")

# Step 3: run binomial test for each pair
mural_counts <- mural_counts %>%
  rowwise() %>%
  mutate(
    binom_p = binom.test(n_success, n_total, p = global_prop, alternative = "greater")$p.value
  ) %>%
  ungroup()

# Step 4: view results
mural_counts

### Function to compute and store nearest distances from one major cell type to another ###

### distance to nearest Epithelial cell (Fig. 6f)
# load seurat object
ChP_filt <- readRDS("data_objects/ChP_filt_082025.rds")


# Define cell types of interest (e.g. all)
types <- unique(ChP_filt$newMajorcelltype)

# Create a unique mapping of individual_ID to AD_status
id_status <- ChP_filt@meta.data %>%
  distinct(individual_ID, AD_status)

# Run distance calculations for all types
all_distance_df <- lapply(types, function(itype) {
  get_nearest_distance_to_celltype(ChP_filt, type = itype, target_types = c("Epithelial"))
}) %>%
  bind_rows()

# Join with metadata
all_distance_df <- left_join(all_distance_df, id_status, by = "individual_ID")

# Drop NAs
global_df <- all_distance_df %>% filter(!is.na(avg_distance_um))

# Make sure the factor is ordered
global_df$newMajorcelltype <- factor(global_df$newMajorcelltype, levels = names(major))

# Boxplot: Compare distances between subtypes
ggplot(global_df, aes(x = newMajorcelltype, y = avg_distance_um, fill = newMajorcelltype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = sub, drop = FALSE) +
  labs(title = "Distance to nearest Epithelial cell by cell type", y = "Avg. Distance (um)", x = "cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Kruskal-Wallis test (nonparametric equivalent of ANOVA)
kruskal.test(avg_distance_um ~ newMajorcelltype, data = global_df)

# export df
write.csv(global_df, file="dist_to_epi.csv")

# repeat with target_types as c("Endothelial","Mural") for distance to vasculature (Fig. 6e)
