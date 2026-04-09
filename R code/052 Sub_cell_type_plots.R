###
### 052 Sub_cell_type_plots.R
###
# Purpose: Produce spatial molecular expression visualizations
# Dependencies:
source("/050 ST_specific_plots_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder

###All subs

DefaultAssay(Cos_ChP) <- "Nanostring"
Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
DefaultBoundary(Cos_ChP@images[["Slide.1"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Sub_Slide1.pdf"))
ImageDimPlot(
  Cos_ChP, 
  fov = "Slide.1", 
  cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Fib_4"="#002fff", "Endo_1"="#c875c4", "Endo_3"="#c5aac5", "Endo_2"="#976697", "Mural"="#261132", "Oligo"="lightsalmon", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97", "Parenchyma"="brown", "low-quality"="lightgray"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Subcelltype"))
dev.off()

DefaultBoundary(Cos_ChP@images[["Slide.2"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Sub_Slide2.pdf"))
ImageDimPlot(
  Cos_ChP, 
  fov = "Slide.2", 
  cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Fib_4"="#002fff", "Endo_1"="#c875c4", "Endo_3"="#c5aac5", "Endo_2"="#976697", "Mural"="#261132", "Oligo"="lightsalmon", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97", "Parenchyma"="brown", "low-quality"="lightgray"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Subcelltype"))
dev.off()



###Epithelial

#Epi <- subset(Cos_ChP, Majorcelltype=="Epithelial")
Idents(Epi) <- Epi@meta.data$Subcelltype

DefaultBoundary(Epi@images[["Slide.1"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide1_Epi.pdf"))
ImageDimPlot(
  Epi, 
  fov = "Slide.1", 
  cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Epi sub; FOV:"))
dev.off()

DefaultBoundary(Epi@images[["Slide.2"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide2_Epi.pdf"))
ImageDimPlot(
  Epi, 
  fov = "Slide.2", 
  cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Epi sub; FOV:"))
dev.off()


####Celltype

major_cols <- c("Epi_1a"="#FFA500", "Epi_1b"="#CD8500", "Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#3249A6", "Fib_2"="#3249A6", "Fib_3"="#3249A6", "Fib_4"="#3249A6", "Endothelial"="#926392", "Endo_1"="#926392", "Endo_3"="#926392", "Endo_2"="#926392", "Mural_2"="#261132", "Mural_3"="#261132", "Mural_1"="#261132", "BAM_1"="#26532B", "BAM_2"="#26532B", "BAM_3"="#26532B", "T_Cell"="#26532B", "Parenchyma"="brown", "low-quality"="lightgrey")
genelist <- c("KRT8", "KRT18", "FOXJ1", "DNAH6", "COX1", "COX2")
valid_FOVs <- get_valid_FOVs(Cos_ChP, genelist, min_cells=400, max_cells=600, min_total_cells=600)
all_fovs <- c(valid_FOVs[[1]], valid_FOVs[[2]], valid_FOVs[[3]], valid_FOVs[[4]], "Slide.1.89")
DefaultAssay(Cos_ChP) <- "Nanostring"

pdf(paste0(curr_dir, "/ST_images/Epi_sub.pdf"))  # Open a PDF device to save the plots
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
    nmols = 300, 
    mols.cols = c("blue", "lightblue", "#E41A1C", "#F781BF", "green3", "lightgreen"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Epi sub; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device



###Fibroblast
####Celltype

#Fib <- subset(Cos_ChP, Majorcelltype=="Fibroblast")
Idents(Fib) <- Fib@meta.data$Subcelltype

DefaultBoundary(Fib@images[["Slide.1"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide1_Fib.pdf"))
ImageDimPlot(
  Fib, 
  fov = "Slide.1", 
  cols=c("Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Fib sub; FOV:"))
dev.off()

DefaultBoundary(Fib@images[["Slide.2"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide2_Fib.pdf"))
ImageDimPlot(
  Fib, 
  fov = "Slide.2", 
  cols=c("Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Fib sub; FOV:"))
dev.off()


####Celltype markers

genelist <- c("KCNMA1", "PLXNA4", "ROBO1", "TRPM3", "LAMA4", "COL4A2", "SLC2A3")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Fibroblast", genelist, min_cells=10, max_cells=300, min_total_cells=200)
all_fovs <- c(valid_FOVs[[1]], valid_FOVs[[2]], valid_FOVs[[3]], valid_FOVs[[4]], valid_FOVs[[5]])
DefaultAssay(Cos_ChP) <- "Nanostring"


major_cols <- c("Epi_1"="#C48F45", "Epi_2"="#C48F45", "Epi_3"="#C48F45", "Epi_4"="#C48F45", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Endo_1"="#926392", "Endo_2"="#926392", "Endo_3"="#926392", "Mural_1"="#261132", "Mural_2"="#261132", "Mural_3"="#261132", "BAM_1"="#26532B", "BAM_2"="#26532B", "BAM_3"="#26532B", "T_Cell"="#26532B", "Parenchyma"="brown", "low-quality"="lightgrey")

pdf(paste0(curr_dir, "/ST_images/Fibroblast_sub_paren.pdf"))  # Open a PDF device to save the plots
for (fov in c("Slide.2.160", "Slide.2.149", "Slide.1.169", "Slide.2.97", "Slide.2.133", "Slide.2.150", "Slide.2.161")) {
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
    mols.cols = c("#E41A1C", "#F781BF", "green3", "lightgreen", "darkgreen", "#FF7F00", "#A65628"),
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Fibroblast sub; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device



###Immune
####Celltype

#Imm <- subset(Cos_ChP, Majorcelltype=="Immune")
Idents(Imm) <- Imm@meta.data$Subcelltype

DefaultBoundary(Imm@images[["Slide.1"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide1_Imm.pdf"))
ImageDimPlot(
  Imm, 
  fov = "Slide.1", 
  cols=c("BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Imm sub; FOV:"))
dev.off()

DefaultBoundary(Imm@images[["Slide.2"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide2_Imm.pdf"))
ImageDimPlot(
  Imm, 
  fov = "Slide.2", 
  cols=c("BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Imm sub; FOV:"))
dev.off()



genelist <- c("STAB1", "CD163", "HSPA1B", "DNAJA1", "HSP90AA1", "HSPA8", "MKI67", "CENPF", "THEMIS", "CD3E")

major_cols <- c("Epi_1"="#C48F45", "Epi_2"="#C48F45", "Epi_3"="#C48F45", "Epi_4"="#C48F45", "Fib_1"="#3249A6", "Fib_2"="#3249A6", "Fib_3"="#3249A6", "Fib_4"="#3249A6", "Endo_1"="#926392", "Endo_3"="#926392", "Endo_2"="#926392", "Mural_2"="#261132", "Mural_3"="#261132", "Mural_1"="#261132", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97", "Parenchyma"="brown", "low-quality"="lightgray")

valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Immune", genelist, min_cells=10, max_cells=200, min_total_cells=400)
all_fovs <- c("Slide.1.83", "Slide.1.71", "Slide.1.183", "Slide.1.201")
DefaultAssay(Cos_ChP) <- "Nanostring"

pdf(paste0(curr_dir, "/ST_images/Immune_sub.pdf"))  # Open a PDF device to save the plots
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
    mols.cols = c("#B41A1C", "#E41A1C", "#F781BF", "blue2", "blue4", "cyan2", "darkcyan", "cadetblue2", "chocolate1", "darkorange", "chocolate3", "darkorange4", "bisque4", "burlywood4", "darkolivegreen"), 
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Immune sub; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device



###Endo
####Celltype

#Endo <- subset(Cos_ChP, Majorcelltype=="Endothelial")

Idents(Endo) <- Endo@meta.data$Subcelltype

DefaultBoundary(Endo@images[["Slide.1"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide1_Endo.pdf"))
ImageDimPlot(
  Endo, 
  fov = "Slide.1", 
  cols=c("Endo_1"="#c875c4", "Endo_3"="#c5aac5", "Endo_2"="#976697"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Endo sub; FOV:"))
dev.off()

DefaultBoundary(Endo@images[["Slide.2"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide2_Endo.pdf"))
ImageDimPlot(
  Endo, 
  fov = "Slide.2", 
  cols=c("Endo_1"="#c875c4", "Endo_3"="#c5aac5", "Endo_2"="#976697"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Endo sub; FOV:"))
dev.off()


####Celltype markers

genelist <- c("CLDN5", "VEGFC", "COL8A1",
              "TSHZ2", "IL1R1", "IL4R",
              "PLVAP", "FLT1") %>% unique()

major_cols <- c("Epi_1"="#C48F45", "Epi_2"="#C48F45", "Epi_3"="#C48F45", "Epi_4"="#C48F45", "Fib_1"="#3249A6", "Fib_2"="#3249A6", "Fib_3"="#3249A6", "Fib_4"="#3249A6", "Endo_1"="#c875c4", "Endo_3"="#c5aac5", "Endo_2"="#976697", "Mural_2"="#8f2aa2", "Mural_3"="#6d207b", "Mural_1"="#261132", "BAM_1"="#26532B", "BAM_2"="#26532B", "BAM_3"="#26532B", "T_Cell"="#26532B", "Parenchyma"="brown", "low-quality"="lightgray")

valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Endo_2", genelist, min_cells=1, max_cells=200, min_total_cells=200, percent_in_cells=5)
all_fovs <- c("Slide.1.186", "Slide.1.64", "Slide.1.168")
DefaultAssay(Cos_ChP) <- "Nanostring"

pdf(paste0(curr_dir, "/ST_images/Endo_sub.pdf"))  # Open a PDF device to save the plots
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
    mols.cols = c("#B41A1C", "#E41A1C", "#F781BF", "cyan2", "darkcyan", "cadetblue2", "chocolate1", "darkorange"), 
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Endo sub; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device



###Mural
####Celltype

Mural <- subset(Cos_ChP, Majorcelltype=="Mural")

Idents(Mural) <- Mural@meta.data$Subcelltype

DefaultBoundary(Mural@images[["Slide.1"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide1_Mural.pdf"))
ImageDimPlot(
  Mural, 
  fov = "Slide.1", 
  cols=c("Mural_2"="#8f2aa2", "Mural_3"="#6d207b", "Mural_1"="#261132"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Mural sub; FOV:"))
dev.off()

DefaultBoundary(Mural@images[["Slide.2"]]) <- 'segmentation'
pdf(paste0(curr_dir, "/ST_images/Slide2_Mural.pdf"))
ImageDimPlot(
  Mural, 
  fov = "Slide.2", 
  cols=c("Mural_2"="#8f2aa2", "Mural_3"="#6d207b", "Mural_1"="#261132"),
  dark.background = FALSE, 
  border.color = NA, 
  coord.fixed = FALSE, 
  border.size = 0.2, 
  mols.size = 0.3
) + ggtitle(paste("Mural sub; FOV:"))
dev.off()


####Celltype markers

genelist <- c("ENPP2", "MID1", "GLIS3",
              "ADAMTS9", "HIPK2", "COL18A1",
              "MYH11", "ACTA2", "TPM1") %>% unique()

major_cols <- c("Epi_1"="#C48F45", "Epi_2"="#C48F45", "Epi_3"="#C48F45", "Epi_4"="#C48F45", "Fib_1"="#3249A6", "Fib_2"="#3249A6", "Fib_3"="#3249A6", "Fib_4"="#3249A6", "Endo_1"="#926392", "Endo_3"="#926392", "Endo_2"="#926392", "Mural_2"="#8f2aa2", "Mural_3"="#6d207b", "Mural_1"="#261132", "BAM_1"="#26532B", "BAM_2"="#26532B", "BAM_3"="#26532B", "T_Cell"="#26532B", "Parenchyma"="brown", "low-quality"="lightgray")

valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Mural", genelist, min_cells=10, max_cells=200, min_total_cells=200, percent_in_cells=5)
all_fovs <- c("Slide.1.201")
DefaultAssay(Cos_ChP) <- "Nanostring"

pdf(paste0(curr_dir, "/ST_images/Mural_sub.pdf"))  # Open a PDF device to save the plots
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
    mols.cols = c("hotpink", "#E41A1C", "pink", "cyan2", "darkcyan", "cadetblue2", "green2", "green4", "darkgreen"), 
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Mural sub; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device


