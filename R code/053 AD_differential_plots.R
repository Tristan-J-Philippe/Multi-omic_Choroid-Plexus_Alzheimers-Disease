###
### 053 AD_differential_plots.R
###
# Purpose: Produce spatial molecular expression visualizations
# Dependencies:
source("/050 ST_specific_plots_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder

###Epithelial
####Cilia

major_cols <- c("Epithelial"="#C48F45", "Fibroblast"="white", "Endothelial"="white", "Mural"="white", "Immune"="white", "Parenchyma"="white", "low-quality"="white")
DefaultAssay(Cos_ChP) <- "Nanostring"

genelist <- c("TUBB4B", "FOXJ1")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Epithelial", genelist, min_cells=400, max_cells=600, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs, decreasing=T)

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs, "Slide.1.89")

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Epi_Cilia.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF", "green3", "lightgreen"), 
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Cilia; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device


scaled_imagedimplot(Cos_ChP, gene="TUBB4B", fov1="Slide.2.23", fov2="Slide.2.124", major_cell_type="Epithelial")
scaled_imagedimplot(Cos_ChP, "FOXJ1", "Slide.2.23", "Slide.2.124", "Epithelial")


####Transport

genelist <- c("CPE", "KCNN2")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Epithelial", genelist, min_cells=400, max_cells=700, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs, decreasing=T)
DefaultAssay(Cos_ChP) <- "Nanostring"

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:10])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Epi_Transport.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#F781BF", "#4DAF4A", "#984EA3", "#A65628"), 
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Transport; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, gene="CPE", fov1="Slide.1.113", fov2="Slide.2.124", major_cell_type = "Epithelial")
scaled_imagedimplot(Cos_ChP, gene="KCNN2", fov1="Slide.1.113", fov2="Slide.2.124", major_cell_type = "Epithelial")


####Metabolism

genelist <- c("NQO1", "PLTP", "UBB")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Epithelial", genelist, min_cells=400, max_cells=600, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)
DefaultAssay(Cos_ChP) <- "Nanostring"

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:5])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:15])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Epi_Metabolism.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#F781BF", "#4DAF4A", "#984EA3", "#A65628"), 
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Metabolism; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device


scaled_imagedimplot(Cos_ChP, "NQO1", "Slide.2.23", "Slide.1.150", "Epithelial")
scaled_imagedimplot(Cos_ChP, "PLTP", "Slide.2.23", "Slide.1.150", "Epithelial")
scaled_imagedimplot(Cos_ChP, "UBB", "Slide.2.23", "Slide.1.150", "Epithelial")


####Immune

genelist <- c("CHI3L1", "IL6R", "SERPINA3")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Epithelial", genelist, min_cells=400, max_cells=600, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)
DefaultAssay(Cos_ChP) <- "Nanostring"

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Epi_Immune.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#F781BF", "#4DAF4A", "#984EA3", "#A65628"), 
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Immune; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device


scaled_imagedimplot(Cos_ChP, gene="CHI3L1", fov1="Slide.1.84", fov2="Slide1.183", major_cell_type="Epithelial")
scaled_imagedimplot(Cos_ChP, gene="IL6R", fov1="Slide.1.84", fov2="Slide1.183", major_cell_type="Epithelial")
scaled_imagedimplot(Cos_ChP, gene="SERPINA3", fov1="Slide.1.84", fov2="Slide1.183", major_cell_type="Epithelial")



####Adhesion

DefaultAssay(Cos_ChP) <- "Nanostring"

genelist <- c("CLDN5", "COL9A3", "LAMB1")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Epithelial", genelist, min_cells=500, max_cells=700, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Epi_Adhesion.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#F781BF", "#4DAF4A", "#984EA3", "#A65628"),,
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Adhesion; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, gene="CLDN5", fov1 = "Slide.1.184", fov2 = "Slide.1.174", major_cell_type = "Epithelial")
scaled_imagedimplot(Cos_ChP, "COL9A3", "Slide.1.184", "Slide.1.174", "Epithelial")
scaled_imagedimplot(Cos_ChP, "LAMB1", "Slide.1.184", "Slide.1.174", "Epithelial")



###Fibroblast

####Adhesion

major_cols <- c("Epi_1a"="#C48F45", "Epi_1b"="#C48F45", "Epi_2a"="#C48F45", "Epi_2b"="#C48F45", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Fib_4"="#3249A6", "Endothelial"="#926392", "Endo_1"="#926392", "Endo_3"="#926392", "Endo_2"="#926392", "Mural_2"="#261132", "Mural_3"="#261132", "Mural_1"="#261132", "BAM_1"="#26532B", "BAM_2"="#26532B", "BAM_3"="#26532B", "T_Cell"="#26532B", "Parenchyma"="brown", "low-quality"="lightgrey")
DefaultAssay(Cos_ChP) <- "Nanostring"

genelist <- c("COL8A1", "COL15A1", "COL1A2", "ADAMTS9", "LAMA2")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Fibroblast", genelist, min_cells=10, max_cells=400, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)
all_fovs <- c("Slide.1.80", "Slide.1.86", "Slide.2.119", "Slide.1.194")

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Fib_Adhesion.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Subcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Adhesion; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

pdf(paste0(curr_dir, "/ST_images/Fib_Adhesion_alltypes.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="lightgrey", "low-quality"="lightgrey"),
    dark.background = FALSE, 
    border.color = NA, 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Subcelltype"))
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, gene="COL8A1", fov1 = "Slide.1.80", fov2 = "Slide.2.119", major_cell_type = "Fibroblast")
scaled_imagedimplot(Cos_ChP, "COL15A1", "Slide.1.80", "Slide.2.119", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "ADAMTS9", "Slide.1.80", "Slide.2.119", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "CD96", "Slide.1.80", "Slide.2.119", "Fibroblast")


####Collagen

major_cols <- c("Epi_1a"="lightgrey", "Epi_1b"="lightgrey", "Epi_2a"="lightgrey", "Epi_2b"="lightgrey", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Endo_1"="lightgrey", "Endo_3"="lightgrey", "Endo_2"="lightgrey", "Mural_1"="lightgrey", "Mural_2"="lightgrey", "Mural_3"="lightgrey","BAM_1"="lightgrey", "BAM_2"="lightgrey", "BAM_3"="lightgrey", "T_Cell"="lightgrey", "Parenchyma"="lightgrey", "low-quality"="lightgrey")
DefaultAssay(Cos_ChP) <- "Nanostring"

genelist <- c("COL15A1", "COL18A1", "COL27A1", "COL4A1", "COL4A2", "COL4A4", "COL5A3", "COL7A1", "COL8A1")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Fibroblast", genelist, min_cells=10, max_cells=400, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)
all_fovs <- c("Slide.1.80", "Slide.1.86", "Slide.2.119", "Slide.1.194")

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Fib_Collagen.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Subcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF",  "#984EA3", "lightgreen", "#4DAF4A", "darkgreen", "orange", "#FF7F00", "#A65628"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Adhesion; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

####Adhesion CC

genelist <- c("LAMA2", "COL4A2", "COL1A2", "JAM3")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Fibroblast", genelist, min_cells=10, max_cells=400, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)
DefaultAssay(Cos_ChP) <- "Nanostring"

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Fib_Adhesion_CC.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Adhesion CC; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

pdf(paste0(curr_dir, "/ST_images/Fib_Adhesion_CC_alltypes.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="lightgrey", "low-quality"="lightgrey"),
    dark.background = FALSE, 
    border.color = NA, 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Subcelltype"))
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, "LAMA2", "Slide.2.126", "Slide.1.183", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "COL4A2", "Slide.2.126", "Slide.1.183", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "COL1A2", "Slide.2.126", "Slide.1.183", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "JAM3", "Slide.2.126", "Slide.1.183", "Fibroblast")



####Transport

genelist <- c("KCNMA1", "SLC38A2", "TFRC")
valid_FOVs <- get_valid_FOVs(Cos_ChP, genelist, major_cell_type="Fibroblast", min_cells=6, max_cells=400, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)
DefaultAssay(Cos_ChP) <- "Nanostring"

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- c(unique(interest_FOVs$nci_fov_names[1:4]), unique(interest_FOVs$nci_fov_names[4:8]))

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Fib_Transport.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF", "#4DAF4A"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Transport; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

pdf(paste0(curr_dir, "/ST_images/Fib_Transport_alltypes.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="lightgrey", "low-quality"="lightgrey"),
    dark.background = FALSE, 
    border.color = NA, 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Subcelltype"))
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, "KCNMA1", "Slide.1.185", "Slide.2.12", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "SLC38A2", "Slide.1.185", "Slide.2.12", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "TFRC", "Slide.1.185", "Slide.2.12", "Fibroblast")


####Metabolism

genelist <- c("APOE", "ACSL4", "APOLD1")
valid_FOVs <- get_valid_FOVs(Cos_ChP, genelist, major_cell_type="Fibroblast", min_cells=6, max_cells = 400, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)
DefaultAssay(Cos_ChP) <- "Nanostring"

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Fib_Metabolism.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Metabolism; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

pdf(paste0(curr_dir, "/ST_images/Fib_Metabolism_alltypes.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="lightgrey", "low-quality"="lightgrey"),
    dark.background = FALSE, 
    border.color = NA, 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Subcelltype"))
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, "APOE", "Slide.1.118", "Slide.1.174", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "ALDH1L1", "Slide.1.118", "Slide.1.174", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "ACSL4", "Slide.1.118", "Slide.1.174", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "APOLD1", "Slide.1.118", "Slide.1.174", "Fibroblast")


####Immune

genelist <- c("CD9", "GPX4", "C7", "GPX3")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type="Fibroblast", genelist, min_cells=6, max_cells = 400, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:40])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

DefaultAssay(Cos_ChP) <- "Nanostring"
# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Fib_Immune.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Immune; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

pdf(paste0(curr_dir, "/ST_images/Fib_Imm_alltypes.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="lightgrey", "low-quality"="lightgrey"),
    dark.background = FALSE, 
    border.color = NA, 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Subcelltype"))
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, "CD9", "Slide.1.118", "Slide.2.133", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "GPX4", "Slide.1.118", "Slide.2.133", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "C7", "Slide.1.118", "Slide.2.133", "Fibroblast")
scaled_imagedimplot(Cos_ChP, "GPX3", "Slide.1.118", "Slide.2.133", "Fibroblast")


###Immune
####HSP

genelist <- c("HSPA1A", "HSP90B1", "HSPA6", "HSPA1B", "DNAJA1", "HSP90AA1")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Immune", genelist, min_cells=10, max_cells=400, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)
DefaultAssay(Cos_ChP) <- "Nanostring"

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:6])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:15])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Imm_HSP.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Subcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF", "#984EA3", "#FF7F00", "#A65628", "lightblue"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("HSP; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

pdf(paste0(curr_dir, "/ST_images/Imm_HSP_alltypes.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="lightgrey", "low-quality"="lightgrey"),
    dark.background = FALSE, 
    border.color = NA, 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Subcelltype"))
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, gene="HSPA1A", fov1 = "Slide.2.36", fov2 = "Slide.2.56", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="HSP90B1", fov1 = "Slide.2.36", fov2 = "Slide.2.56", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="HSPA6", fov1 = "Slide.2.36", fov2 = "Slide.2.56", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="HSPA1B", fov1 = "Slide.2.36", fov2 = "Slide.2.56", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="DNAJA1", fov1 = "Slide.2.36", fov2 = "Slide.2.56", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="HSP90AA1", fov1 = "Slide.2.36", fov2 = "Slide.2.56", major_cell_type = "Immune")


####Inflammation

genelist <- c("CD68", "HLA.DRA", "HLA.DRB", "HLA.DQA1", "HLA.E")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Immune", genelist, min_cells=10, max_cells=400, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs)
DefaultAssay(Cos_ChP) <- "Nanostring"

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:20])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:20])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs)
all_fovs <- all_fovs[!(is.na(all_fovs))]
all_fovs <- c("Slide.2.24", "Slide.2.181", "Slide.2.23", "Slide.1.200", "Slide.1.185") #These have T cells

genelist <- c("CD68", "HLA-DRA", "HLA-DRB", "HLA-DQA1", "HLA-E") #They are a dash in the @images
# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Imm_Inflammation.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Subcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    mols.cols = c("#E41A1C", "#F781BF", "#984EA3", "#FF7F00", "#A65628"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Inflammation; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

pdf(paste0(curr_dir, "/ST_images/Imm_Inflammation_alltypes.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B","Parenchyma"="brown", "low-quality"="lightgrey"),
    dark.background = FALSE, 
    border.color = NA, 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Subcelltype"))
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device

scaled_imagedimplot(Cos_ChP, gene="CD68", fov1 = "Slide.2.20", fov2 = "Slide.1.200", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="HLA.DRA", fov1 = "Slide.2.20", fov2 = "Slide.1.200", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="HLA.DRB", fov1 = "Slide.2.20", fov2 = "Slide.1.200", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="HLA.DQA1", fov1 = "Slide.2.20", fov2 = "Slide.1.200", major_cell_type = "Immune")
scaled_imagedimplot(Cos_ChP, gene="HLA.E", fov1 = "Slide.2.20", fov2 = "Slide.1.200", major_cell_type = "Immune")
