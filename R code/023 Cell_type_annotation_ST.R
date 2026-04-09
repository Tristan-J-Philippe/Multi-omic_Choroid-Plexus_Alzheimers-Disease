###
### 023 Cell_type_annotation_ST.R
###
# Purpose:  Cell type annotation and cleaning of spatial transcriptomics data
# Dependencies:
source("/ChP_ST/020 Cell_type_annotation_util.R", verbose = FALSE)
library(Azimuth)
library(SpatialExperiment)
library(SummarizedExperiment)
library(spacexr)
library(stringr)
library(RhpcBLASctl)
library(ComplexHeatmap)
Sys.setenv(OMP_NUM_THREADS = "20", MKL_NUM_THREADS = "20",
           OPENBLAS_NUM_THREADS = "20", NUMEXPR_NUM_THREADS = "20", R_MAX_NUM_THREADS = "20")
blas_set_num_threads(20)
omp_set_num_threads(20)

curr_dir <- "/ChP_ST" #output results to specific folder
set.seed(1234)
options(future.globals.maxSize = 5*1024^3)

#Load data and prep data
##cosmx_1

#CosMX_1 <- LoadNanostring(data.dir= "/Ex_ChP/Atmx_export/CosChP_ctrl_FLAT/", fov="Slide.1") #Test only

#Load data
mat <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_ctrl_FLAT/Control_ChP1_14_exprMat_file.csv.gz"))
mdat <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_ctrl_FLAT/Control_ChP1_14_metadata_file.csv.gz"))
molecs <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_ctrl_FLAT/Control_ChP1_14_tx_file.csv.gz"))
segs <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_ctrl_FLAT/Control_ChP1_14-polygons.csv.gz"))
fov <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_ctrl_FLAT/Control_ChP1_14_fov_positions_file.csv.gz"))

#choose metadata
mdata <- mdat[,c("cell", "Area", "fov", "AspectRatio", "Mean.PanCK", "Max.PanCK", "Mean.CD68_CK8_18", "Max.CD68_CK8_18", "Mean.CD298_B2M", "Max.CD298_B2M", "Mean.CD45", "Max.CD45", "Mean.DAPI", "Max.DAPI", "cell_id", "nCount_RNA", "nFeature_RNA", "Area.um2", "propNegative", "complexity", "errorCtEstimate", "percOfDataFromError", "qcFlagsCellCounts", "qcFlagsCellPropNeg", "qcFlagsCellComplex", "qcFlagsCellArea", "qcCellsFlagged", "nCell", "nCount", "nCountPerCell", "nFeaturePerCell", "propNegativeCellAvg", "complexityCellAvg", "errorCtPerCellEstimate", "percOfDataFromErrorPerCell")]

#make seurat object
rownames(mat) <- paste("c_1", mat$fov, mat$cell_ID, sep="_")
mat <- t(mat)
mat <- Matrix(mat, sparse=T)
mat <- mat[-grep("Negative*|SystemControl*|fov|cell_ID", mat@Dimnames[[1]]),]
mdata <- column_to_rownames(mdata, "cell")
CosMX_1 <- CreateSeuratObject(counts = mat, assay = "Nanostring", meta.data=mdata)

#add global spatial
seg <- CreateSegmentation(data.frame(cell=segs$cell, x=segs$x_global_px, y=segs$y_global_px, stringsAsFactors = F))
cent <- CreateCentroids(data.frame(x=mdat$CenterX_global_px, y=mdat$CenterY_global_px, stringsAsFactors = F))
molec <- data.frame(x=molecs$x_global_px, y=molecs$y_global_px, gene=molecs$target, stringsAsFactors = F)
cent@cells <- mdat$cell
segmentations.data <- list(centroids = cent, segmentation = seg)
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                          "centroids"), molecules = molec, assay = "Nanostring")
cells <- intersect(Cells(x = coords, boundary = "segmentation"), 
                   Cells(x = coords, boundary = "centroids"))
cells <- intersect(Cells(CosMX_1), cells)
coords <- subset(x = coords, cells = cells)
CosMX_1[["Slide.1"]] <- coords

#add local spatial
for(i in unique(mdat$fov))
{
  #select FOV
  seg <- segs[segs$fov==i,]
  mdati <- mdat[mdat$fov==i,]
  molec <- molecs[molecs$fov==i,]
  
  #local
  seg <- data.frame(cell=seg$cell, x=seg$x_local_px, y=seg$y_local_px, stringsAsFactors = F)
  cent <- data.frame(x=mdati$CenterX_local_px, y=mdati$CenterY_local_px, stringsAsFactors = F)
  molec <- data.frame(x=molec$x_local_px, y=molec$y_local_px, gene=molec$target, stringsAsFactors = F)
  
  #make coords
  seg <- CreateSegmentation(seg)
  cent <- CreateCentroids(cent)
  cent@cells <- mdati$cell
  segmentations.data <- list(centroids = cent, segmentation = seg)
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                            "centroids"), molecules = molec, assay = "Nanostring")
  cells <- intersect(Cells(x = coords, boundary = "segmentation"), 
                     Cells(x = coords, boundary = "centroids"))
  cells <- intersect(Cells(CosMX_1), cells)
  coords <- subset(x = coords, cells = cells)
  CosMX_1[[paste0("Slide.1.", i)]] <- coords #add fov
}

#Verify slide
DefaultBoundary(CosMX_1@images[["Slide.1"]]) <- 'segmentation'
ImageDimPlot(CosMX_1, fov = "Slide.1", axes = TRUE)
DefaultBoundary(CosMX_1@images[["Slide.1"]]) <- 'centroid'
ImageDimPlot(CosMX_1, fov = "Slide.1", axes = TRUE)
ImageDimPlot(CosMX_1, molecules=c("COX2"), fov="Slide.1", axes = TRUE, nmols=10000)
ImageFeaturePlot(CosMX_1, features=c("COX2", "TTR"), fov="Slide.1")
#Verify fovs
DefaultBoundary(CosMX_1@images[["Slide.1.1"]]) <- 'segmentation'
ImageDimPlot(CosMX_1, fov = "Slide.1.1", axes = TRUE)
ImageFeaturePlot(CosMX_1, features=c("COX2", "TTR"), fov="Slide.1.1")
DefaultBoundary(CosMX_1@images[["Slide.1.1"]]) <- 'centroid'
ImageDimPlot(CosMX_1, fov = "Slide.1.1", axes = TRUE)
ImageDimPlot(CosMX_1, molecules=c("COX2","TTR"), fov="Slide.1.1", axes = TRUE) 
ImageFeaturePlot(CosMX_1, features=c("COX2", "TTR"), fov="Slide.1.1")

##choroid 6k ctrl slide
#note that the slide is flipped vertically
#note that projIDs were replaced with fkID
# upper left 1-46 fkID #127 NCI F
# upper right 47-104 fkID #144 AD F
# lower left 105-146  fkID #157 AD M
# lower right 147-212 fkID #139 NCI M
CosMX_1@meta.data$projID <- CosMX_1@meta.data$fov
CosMX_1@meta.data$projID <- replace(CosMX_1@meta.data$projID, CosMX_1@meta.data$fov<=46 & CosMX_1@meta.data$fov>=1, "127")
CosMX_1@meta.data$projID <- replace(CosMX_1@meta.data$projID, CosMX_1@meta.data$fov<=104 & CosMX_1@meta.data$fov>=47, "144")
CosMX_1@meta.data$projID <- replace(CosMX_1@meta.data$projID, CosMX_1@meta.data$fov<=146 & CosMX_1@meta.data$fov>=105, "157")
CosMX_1@meta.data$projID <- replace(CosMX_1@meta.data$projID, CosMX_1@meta.data$fov<=212 & CosMX_1@meta.data$fov>=147, "139")
#sanitycheck
test <- data.frame(CosMX_1@meta.data$fov, CosMX_1@meta.data$projID)
test <- test[!duplicated(test$CosMX_1.meta.data.fov),]
CosMX_1@meta.data$Slide <- 1
saveRDS(CosMX_1, paste0(curr_dir, "/Sample_preprocess/CosMX_1.RDS"))


##cosmx_2

#Load data
mat <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_AD_flat/AD_ChP2_9_exprMat_file.csv.gz"))
mdat <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_AD_flat/AD_ChP2_9_metadata_file.csv.gz"))
molecs <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_AD_flat/AD_ChP2_9_tx_file.csv.gz"))
segs <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_AD_flat/AD_ChP2_9-polygons.csv.gz"))
fov <- read.csv(gzfile("/Ex_ChP/Atmx_export/CosChP_AD_flat/AD_ChP2_9_fov_positions_file.csv.gz"))

#choose metadata
mdata <- mdat[,c("cell", "Area", "fov", "AspectRatio", "Mean.PanCK", "Max.PanCK", "Mean.CD68_CK8_18", "Max.CD68_CK8_18", "Mean.CD298_B2M", "Max.CD298_B2M", "Mean.CD45", "Max.CD45", "Mean.DAPI", "Max.DAPI", "cell_id", "nCount_RNA", "nFeature_RNA", "Area.um2", "propNegative", "complexity", "errorCtEstimate", "percOfDataFromError", "qcFlagsCellCounts", "qcFlagsCellPropNeg", "qcFlagsCellComplex", "qcFlagsCellArea", "qcCellsFlagged", "nCell", "nCount", "nCountPerCell", "nFeaturePerCell", "propNegativeCellAvg", "complexityCellAvg", "errorCtPerCellEstimate", "percOfDataFromErrorPerCell")]

#make seurat object
rownames(mat) <- paste("c_2", mat$fov, mat$cell_ID, sep="_")
mat <- t(mat)
mat <- Matrix(mat, sparse=T)
mat <- mat[-grep("Negative*|SystemControl*|fov|cell_ID", mat@Dimnames[[1]]),]
mdata <- column_to_rownames(mdata, "cell")
CosMX_2 <- CreateSeuratObject(counts = mat, assay = "Nanostring", meta.data=mdata)

#add global spatial
seg <- CreateSegmentation(data.frame(cell=segs$cell, x=segs$x_global_px, y=segs$y_global_px, stringsAsFactors = F))
cent <- CreateCentroids(data.frame(x=mdat$CenterX_global_px, y=mdat$CenterY_global_px, stringsAsFactors = F))
molec <- data.frame(x=molecs$x_global_px, y=molecs$y_global_px, gene=molecs$target, stringsAsFactors = F)
cent@cells <- mdat$cell
segmentations.data <- list(centroids = cent, segmentation = seg)
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                          "centroids"), molecules = molec, assay = "Nanostring")
cells <- intersect(Cells(x = coords, boundary = "segmentation"), 
                   Cells(x = coords, boundary = "centroids"))
cells <- intersect(Cells(CosMX_2), cells)
coords <- subset(x = coords, cells = cells)
CosMX_2[["Slide.2"]] <- coords

#add local spatial
for(i in unique(mdat$fov))
{
  #select FOV
  seg <- segs[segs$fov==i,]
  mdati <- mdat[mdat$fov==i,]
  molec <- molecs[molecs$fov==i,]
  
  #local
  seg <- data.frame(cell=seg$cell, x=seg$x_local_px, y=seg$y_local_px, stringsAsFactors = F)
  cent <- data.frame(x=mdati$CenterX_local_px, y=mdati$CenterY_local_px, stringsAsFactors = F)
  molec <- data.frame(x=molec$x_local_px, y=molec$y_local_px, gene=molec$target, stringsAsFactors = F)
  
  #make coords
  seg <- CreateSegmentation(seg)
  cent <- CreateCentroids(cent)
  cent@cells <- mdati$cell
  segmentations.data <- list(centroids = cent, segmentation = seg)
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                            "centroids"), molecules = molec, assay = "Nanostring")
  cells <- intersect(Cells(x = coords, boundary = "segmentation"), 
                     Cells(x = coords, boundary = "centroids"))
  cells <- intersect(Cells(CosMX_2), cells)
  coords <- subset(x = coords, cells = cells)
  CosMX_2[[paste0("Slide.2.", i)]] <- coords #add fov
}

#Verify slide
DefaultBoundary(CosMX_2@images[["Slide.2"]]) <- 'segmentation'
ImageDimPlot(CosMX_2, fov = "Slide.2", axes = TRUE)
DefaultBoundary(CosMX_2@images[["Slide.2"]]) <- 'centroid'
ImageDimPlot(CosMX_2, fov = "Slide.2", axes = TRUE)
ImageDimPlot(CosMX_2, molecules=c("COX2"), fov="Slide.2", axes = TRUE, nmols=10000)
ImageFeaturePlot(CosMX_2, features=c("COX2", "TTR"), fov="Slide.2")
#Verify fovs
DefaultBoundary(CosMX_2@images[["Slide.2.1"]]) <- 'segmentation'
ImageDimPlot(CosMX_2, fov = "Slide.2.1", axes = TRUE)
DefaultBoundary(CosMX_2@images[["Slide.2.1"]]) <- 'centroid'
ImageDimPlot(CosMX_2, fov = "Slide.2.1", axes = TRUE)
ImageDimPlot(CosMX_2, molecules=c("COX2","TTR"), fov="Slide.2.1", axes = TRUE) 
ImageFeaturePlot(CosMX_2, features=c("COX2", "TTR"), fov="Slide.2.1")

## choroid 6k ad slide
# upper left 1-47 projID #181 AD F
# upper right 48-108 projID #106 NCI F
# lower left 109-170 projID #105 NCI M
# lower right 171-225 projID #183 AD M

CosMX_2@meta.data$projID <- CosMX_2@meta.data$fov
CosMX_2@meta.data$projID <- replace(CosMX_2@meta.data$projID, CosMX_2@meta.data$fov<=47 & CosMX_2@meta.data$fov>=1, "181")
CosMX_2@meta.data$projID <- replace(CosMX_2@meta.data$projID, CosMX_2@meta.data$fov<=108 & CosMX_2@meta.data$fov>=48, "106")
CosMX_2@meta.data$projID <- replace(CosMX_2@meta.data$projID, CosMX_2@meta.data$fov<=170 & CosMX_2@meta.data$fov>=109, "105")
CosMX_2@meta.data$projID <- replace(CosMX_2@meta.data$projID, CosMX_2@meta.data$fov<=225 & CosMX_2@meta.data$fov>=171, "183")

test <- data.frame(CosMX_2@meta.data$fov, CosMX_2@meta.data$projID)
test <- test[!duplicated(test$CosMX_2.meta.data.fov),]
CosMX_2@meta.data$Slide <- 2
saveRDS(CosMX_2, paste0(curr_dir, "/Sample_preprocess/CosMX_2.RDS"))



Cos_ChP <- merge(CosMX_1, CosMX_2)
DefaultAssay(Cos_ChP) <- "Nanostring"
Cos_ChP <- JoinLayers(Cos_ChP)

hist(Cos_ChP@meta.data[["nCount_Nanostring"]], breaks=40)
hist(Cos_ChP@meta.data[["nFeature_Nanostring"]], breaks=50)
Cos_ChP <- subset(Cos_ChP, subset = nCount_Nanostring>100 & nCount_Nanostring<3500 & nFeature_Nanostring>50 & nFeature_Nanostring<2000 & qcCellsFlagged==FALSE & Area.um2>10)

saveRDS(Cos_ChP, paste(curr_dir, "/CellType/Firstcelltype/Cos_ChP_image.RDS", sep=""))




#Celltyping
##Wholeset
###preprocess no correction

hist(Cos_ChP@meta.data[["nCount_Nanostring"]], breaks=40)
hist(Cos_ChP@meta.data[["nFeature_Nanostring"]], breaks=50)
Cos_ChP <- subset(Cos_ChP, subset = nCount_Nanostring>100 & nCount_Nanostring<3500 & nFeature_Nanostring>50 & nFeature_Nanostring<2000 & qcCellsFlagged==FALSE & Area.um2>10)

SCTPCAfun(Cos_ChP, assay="Nanostring")#, spikein=ChP_genes) #first pass

use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Cos_ChP <- RunHarmony(object = Cos_ChP, group.by.vars = c("Slide", "projID"), dims.use = use_dims, assay="SCT", max_iter = 50, plot_convergence = TRUE) 
CNDRfun(Cos_ChP, PCA_dims = use_dims, assay="SCT", reduction = "harmony")

#Pick res
res="0.3"
Idents(Cos_ChP) <- Cos_ChP@meta.data[[paste("SCT_snn_res.", res, sep="")]]
SCNDRdimfun(Cos_ChP, reduction = "harmony", res=res, split.by=c("Slide", "projID", "qcCellsFlagged"))

#Reference markers
Idents(Cos_ChP) <- Cos_ChP@meta.data[["predicted.Majorcelltype"]]
plts <- DimPlot(Cos_ChP, raster=T, reduction='umap')
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Cos_ChP/ctrl_Ex_ChP_markers.png"), plot=plts)

Idents(Cos_ChP) <- Cos_ChP@meta.data[[paste("SCT_snn_res.", res, sep="")]]
Cos_ChP <- JoinLayers(Cos_ChP)
return_markersfun(Cos_ChP, assay="SCT", features = Cos_ChP_Nanostring_var_feat)

Featurefun(Cos_ChP, c("HTR2C", "TTR", "CLIC6", "AQP1", "TUBB4B", "FOXJ1", "ARL13B", "KRT18", "COX2", "TRPM3", "VIM", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "FOXP2", "COL6A3", "SLC6A13", "LUM", "DCN", "FN1", "ACTA2", "KCNMA1", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "CD3E", "HSPE1", "HSPH1"), "Key Markers",  assay="SCT")


###clusters

curr_clust <- 0
Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_BAM, known_peric, known_Mural_1), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Cos_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Cos_ChP, curr_clust)

search_markersfun(Cos_ChP, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Cos_ChP, curr_clust, known_fibro_prog)

search_markersfun(Cos_ChP, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Cos_ChP, curr_clust, c(known_peric, known_Mural_1) %>% unique())

search_markersfun(Cos_ChP, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Cos_ChP, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Cos_ChP, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Cos_ChP, curr_clust, c(known_light) %>% unique())
search_markersfun(Cos_ChP, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Cos_ChP, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Cos_ChP, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Cos_ChP, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Cos_ChP, curr_clust, known_dividing %>% unique())

search_markersfun(Cos_ChP, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Cos_ChP, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Cos_ChP, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

search_markersfun(Cos_ChP, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Cos_ChP, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Cos_ChP, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Cos_ChP, "LAMA4", raster=F)


####Major types
#0 Epi TTR, HTR2C, KRT18, AQP4, RSPO2, FOXJ1, TUBB4B, KRT18, GPX3 lower key markers Epi_2
#1 Epi TTR, HTR2C, KRT18, AQP4, RSPO2,  COXs Epi_4
#2 Epi TTR, HTR2C, TRPM3, AQP1, CLDN5, FOXJ1, TUBB4B, KRT18, GPX3 Epi_1
#3 Epi TTR, HTR2C, TRPM3, CLDN5, FOXJ1, DNAHs, TUBB4B, SLC12A2, GPX3, COXs Epi_3
#4 Fibro-Endo-Mural mix LAMA2, DCN, LUM, ACTA2, KCNMA1, PECAM1, VWF, TAGLN, COL4A1, COL6A3, COLs, CALD1,
#5 Imm CD163, LYVE1, APOE, PTPRC, F13A1, P2RY12, MS4A7, MS4A6A, MS4A4A, SPP1, STABB1, (BAM)
#6 Paren MBP, GFAP, MOBP (Astro+Olig+Olig prog)

###Rename all identities

Idents(Cos_ChP) <- Cos_ChP@meta.data[[paste("SCT_snn_res.", res, sep="")]]
Cos_ChP <- RenameIdents(Cos_ChP, "0"="Epi", "1"="Epi", "2"="Epi", "3"="Epi", "4"="Fibro", "5"="Imm", "6"="Paren")
Cos_ChP@meta.data$Majorcelltype <- Idents(Cos_ChP)
Idents(Cos_ChP) <- Cos_ChP@meta.data[[paste("SCT_snn_res.", res, sep="")]]
Cos_ChP <- RenameIdents(Cos_ChP, "0"="Epi_2", "1"="Epi_4", "2"="Epi_1", "3"="Epi_3", "4"="Fibro", "5"="Imm", "6"="Paren")
Cos_ChP@meta.data$Subcelltype <- Idents(Cos_ChP)


###Subtypes

saveRDS(Cos_ChP, paste(curr_dir, "/CellType/Firstcelltype/Cos_ChP.RDS", sep=""))

for(i in unique(Cos_ChP@meta.data$Majorcelltype)){
  sub <- subset(Cos_ChP, Majorcelltype==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/", i, ".RDS", sep=""))
}


###Overlap between SN and ST variable genes

plts <- ggVennDiagram(list(Cos_ChP_Nanostring_var_feat, Ex_ChP_soupXcounts_var_feat), category.names = c("CosMX", "SN")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2png(plts, paste0(curr_dir, "/CellType/SNvCosMX/Cos_ChP_variable_features_overlap.png"))


##Epi
###preprocess

Epi <- readRDS("/Ex_ChP/CosMX2/CellType/Firstcelltype/Epi.RDS")
Epi_meta <- Epi@meta.data
Epi_mat <- Epi@assays[["Nanostring"]]@layers[["counts"]]
colnames(Epi_mat) <- colnames(Epi)
rownames(Epi_mat) <- rownames(Epi)
Epi <- CreateSeuratObject(counts=Epi_mat, assay="Nanostring", meta.data=Epi_meta)
SCTPCAfun(Epi, assay="Nanostring")

use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Epi <- RunHarmony(object = Epi, group.by.vars = c("Slide", "projID"), dims.use = use_dims, assay="SCT", max.iter.harmony = 50, plot_convergence = TRUE) 
CNDRfun(Epi, PCA_dims = use_dims, assay="SCT", reduction = "harmony")

#Pick res
res="0.3"
Idents(Epi) <- Epi@meta.data[[paste("SCT_snn_res.", res, sep="")]]
SCNDRdimfun(Epi, reduction = "harmony", res=res, split.by= c("Slide", "projID"))

#Reference markers
Idents(Epi) <- Epi@meta.data[["Predicted.Subcelltype"]]
plts <- DimPlot(Epi, raster=T, reduction='umap')
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Epi/ctrl_Ex_ChP_Submarkers.png"), plot=plts)

Idents(Epi) <- Epi@meta.data[[paste("SCT_snn_res.", res, sep="")]]
return_markersfun(Epi, assay="SCT", features = Epi_Nanostring_var_feat)

Featurefun(Epi, c("HTR2C", "TTR", "CLIC6", "AQP1", "TUBB4B", "FOXJ1", "ARL13B", "KRT18", "COX2", "TRPM3", "VIM", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "FOXP2", "COL6A3", "SLC6A13", "LUM", "DCN", "FN1", "ACTA2", "KCNMA1", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "GFAP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "CD3E", "HSPE1", "HSPH1", "GPX3", "KRT17"), "Key Markers",  assay="SCT")



###Overlap between SN and ST variable genes

plts <- ggVennDiagram(list(Epi_Nanostring_var_feat, Epi_soupXcounts_var_feat), category.names = c("CosMX", "SN")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2png(plts, paste0(curr_dir, "/CellType/SNvCosMX/Epithelial_variable_features_overlap.png"))


###clusters

curr_clust <- 3
Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_BAM, known_peric, known_Mural_1), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Epi, curr_clust)

search_markersfun(Epi, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Epi, curr_clust, known_fibro_prog)

search_markersfun(Epi, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Epi, curr_clust, c(known_BAM, known_peric, known_Mural_1) %>% unique())

search_markersfun(Epi, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Epi, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Epi, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Epi, curr_clust, c(known_light) %>% unique())
search_markersfun(Epi, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Epi, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Epi, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Epi, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Epi, curr_clust, known_dividing %>% unique())

search_markersfun(Epi, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Epi, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Epi, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

search_markersfun(Epi, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Epi, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Epi, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Epi, "LAMA4", raster=F)


###Subtypes
#0 TTR, HTR2C, KRT18 low classical Epi_2
#2 TTR, HTR2C, KRT18 low classical Epi_2
#1 TTR, GPX3, CLDN, APP, FOXJ1, TUBB4B, ARL13B, KRT18 Epi_1
#3 TTR, GPX3, CLDN, APP, FOXJ1, TUBB4B, ARL13B, KRT18 Epi_1
#5 TTR, HTR2C, COX high Epi_4
#6 TTR, HTR2C, COX high Epi_4
#4 TTR, GPX3, CLDN, APP, FOXJ1, TUBB4B, ARL13B, KRT18, COX high Epi_3

###Rename all identities

Epi <- RenameIdents(Epi, "0"="Epi_2", "1"="Epi_1", "2"="Epi_2", "3"="Epi_1", "4"="Epi_3", "5"="Epi_4", "6"="Epi_4")
Epi@meta.data$Subcelltype <- Idents(Epi)
DimPlot(Epi)
saveRDS(Epi, "/Ex_ChP/CosMX2/CellType/Finalcelltype/Epi.RDS")


###Remove doublets

Epi <- subset(Epi, SCT_snn_res.0.5==5, invert=T)
Epi <- subset(Epi, SCT_snn_res.0.5==6, invert=T)
Epi <- subset(Epi, SCT_snn_res.0.5==7, invert=T)
Epi <- subset(Epi, SCT_snn_res.0.5==8, invert=T)
#Remove cells that don't belong
#Epi <- subset(Epi, HTR2C < 3)
#Epi <- subset(Epi, TTR < 3)
#Epi <- subset(Epi, CLIC6 < 3)
#Epi <- subset(Epi, AQP1 < 3)
#Epi <- subset(Epi, FOXJ1 < 3)
#Epi <- subset(Epi, FAM227A < 3)
#Epi <- subset(Epi, ARL13B < 3)
#Epi <- subset(Epi, TRPM3 < 3)
Epi <- subset(Epi, SLC4A4 < 3)
Epi <- subset(Epi, LAMA2 < 3)
Epi <- subset(Epi, LAMA4 < 3)
Epi <- subset(Epi, PRICKLE1 < 3)
Epi <- subset(Epi, MECOM < 3)
Epi <- subset(Epi, VWF < 3)
Epi <- subset(Epi, TAGLN < 3)
Epi <- subset(Epi, MYH11 < 3)
Epi <- subset(Epi, NOTCH3 < 3)
Epi <- subset(Epi, MAP2 < 3)
#Epi <- subset(Epi, CLU < 3)
#Epi <- subset(Epi, AQP4 < 3)
Epi <- subset(Epi, MOBP < 3)
Epi <- subset(Epi, MBP < 3)
Epi <- subset(Epi, GFAP < 3)
Epi <- subset(Epi, CD68 < 3)
Epi <- subset(Epi, F13A1 < 3)
Epi <- subset(Epi, STAB1 < 3)
Epi <- subset(Epi, CD163 < 3)
Epi <- subset(Epi, LYVE1 < 3)
Epi <- subset(Epi, P2RY12 < 3)
Epi <- subset(Epi, CD3E < 3)

Idents(Epi) <- Epi@meta.data$Subcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi1")
Epi_1 <- subset(Epi, ident="Epi1")
Epi_1 <- subset(Epi_1, Subcelltype=="Epi_1")

Idents(Epi) <- Epi@meta.data$Subcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi2")
Epi_2 <- subset(Epi, ident="Epi2")
Epi_2 <- subset(Epi_2, Subcelltype=="Epi_2")

Idents(Epi) <- Epi@meta.data$Subcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi3")
Epi_3 <- subset(Epi, ident="Epi3")
Epi_3 <- subset(Epi_3, Subcelltype=="Epi_3")

Idents(Epi) <- Epi@meta.data$Subcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi4")
Epi_4 <- subset(Epi, ident="Epi4")
Epi_4 <- subset(Epi_4, Subcelltype=="Epi_4")

Epi_1_meta <- Epi_1@meta.data
Epi_2_meta <- Epi_2@meta.data
Epi_3_meta <- Epi_3@meta.data
Epi_4_meta <- Epi_4@meta.data
merged_meta <- rbind(Epi_1_meta, Epi_2_meta, Epi_3_meta, Epi_4_meta)

Epi@meta.data$Subcelltype <- NULL
Epi_meta <- Epi@meta.data

Epi_meta <- rownames_to_column(Epi_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Epi_meta <- left_join(Epi_meta, merged_meta, by="rowname")
Epi_meta <- column_to_rownames(Epi_meta, "rowname")
Epi_meta$Subcelltype <- as.character(Epi_meta$Subcelltype)
Epi_meta[is.na(Epi_meta)==TRUE] <- "Doublet"
Epi@meta.data <- Epi_meta
Epi <- subset(Epi, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Epi) <- Epi@meta.data[["Subcelltype"]] 
levels(Epi) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4")
Epi@meta.data[["Subcelltype"]] <- Idents(Epi)
DimPlot(Epi, raster=F)




##Fibro
###preprocess

Fibro <- readRDS("/Ex_ChP/CosMX2/CellType/Firstcelltype/Fibro.RDS")
Fibro_meta <- Fibro@meta.data
Fibro_mat <- Fibro@assays[["Nanostring"]]@layers[["counts"]]
colnames(Fibro_mat) <- colnames(Fibro)
rownames(Fibro_mat) <- rownames(Fibro)
Fibro <- CreateSeuratObject(counts=Fibro_mat, assay="Nanostring", meta.data=Fibro_meta)

#Fibro <- readRDS("/Ex_ChP/CosMX2/CellType/Secondcelltype/Fib.RDS")
SCTPCAfun(Fibro, assay="Nanostring")#, spikein=ChP_genes) #first pass to remove doublets
use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Fibro <- RunHarmony(object = Fibro, group.by.vars =  c("Slide", "projID"), dims.use = use_dims, assay="SCT", max.iter.harmony = 50, plot_convergence = TRUE) 
CNDRfun(Fibro, PCA_dims = use_dims, assay="SCT", reduction = "harmony")

#Pick res
#res="0.7" #second pass
res="0.3" #last pass
Idents(Fibro) <- Fibro@meta.data[[paste("SCT_snn_res.", res, sep="")]]
SCNDRdimfun(Fibro, reduction = "harmony", res=res, split.by= c("Slide", "projID"))

#Reference markers
Idents(Fibro) <- Fibro@meta.data[["predicted.Subcelltype"]]
plts <- DimPlot(Fibro, raster=T, reduction='umap')
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Fibro/ctrl_Ex_ChP_Submarkers.png"), plot=plts)

Idents(Fibro) <- Fibro@meta.data[[paste("SCT_snn_res.", res, sep="")]]
return_markersfun(Fibro, assay="SCT", features = Fibro_Nanostring_var_feat)

Featurefun(Fibro, c("HTR2C", "TTR", "CLIC6", "AQP1", "TUBB4B", "FOXJ1", "ARL13B", "KRT18", "COX2", "TRPM3", "VIM", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "FOXP2", "COL6A3", "SLC6A13", "LUM", "DCN", "FN1", "ACTA2", "KCNMA1", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "CD3E", "HSPE1", "HSPH1"), "Key Markers",  assay="SCT")


###Overlap between SN and ST

plts <- ggVennDiagram(list(Fibro_Nanostring_var_feat, Fibro_soupXcounts_var_feat), category.names = c("CosMX", "SN")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2png(plts, paste0(curr_dir, "/CellType/SNvCosMX/Fibro_Cos_ChP_variable_features_overlap.png"))


###clusters

curr_clust <- 3
Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_BAM, known_peric, known_Mural_1), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Fibro, curr_clust)

search_markersfun(Fibro, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Fibro, curr_clust, known_fibro_prog)

search_markersfun(Fibro, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_light) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Fibro, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Fibro, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Fibro, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Fibro, curr_clust, known_dividing %>% unique())

search_markersfun(Fibro, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Fibro, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Fibro, "LAMA4", raster=F)


###Rename all identities

#Second pass res=0.7
Idents(Fibro) <- Fibro@meta.data[[paste("SCT_snn_res.", res, sep="")]]
Fibro <- RenameIdents(Fibro, "0"="Fibroblast", "1"="Mural", "2"="Fibroblast", "3"="Endo", "4"="Mural", "5"="Mural", "6"="Fibroblast")
Fibro@meta.data$Majorcelltype <- Idents(Fibro)
Idents(Fibro) <- Fibro@meta.data[[paste("SCT_snn_res.", res, sep="")]]
Fibro <- RenameIdents(Fibro, "0"="Fibroblast", "1"="Mural", "2"="Fibroblast", "3"="Endo", "4"="Mural", "5"="Mural", "6"="Fibroblast")
Fibro@meta.data$Subcelltype <- Idents(Fibro)

for(i in unique(Fibro@meta.data$Majorcelltype)){
  sub <- subset(Fibro, Majorcelltype==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Secondcelltype/", i, ".RDS", sep=""))
}

saveRDS(Fibro, "/Ex_ChP/CosMX2/CellType/Secondcelltype/Fibro.RDS")



###Remove doublets

#Remove cells that don't belong
Fibro <- subset(Fibro, HTR2C < 3)
Fibro <- subset(Fibro, TTR < 3)
Fibro <- subset(Fibro, CLIC6 < 3)
Fibro <- subset(Fibro, AQP1 < 3)
Fibro <- subset(Fibro, FOXJ1 < 3)
Fibro <- subset(Fibro, FAM227A < 3)
Fibro <- subset(Fibro, ARL13B < 3)
#Fibro <- subset(Fibro, TRPM3 < 3)
#Fibro <- subset(Fibro, SLC4A4 < 3)
#Fibro <- subset(Fibro, LAMA2 < 3)
#Fibro <- subset(Fibro, LAMA4 < 3)
#Fibro <- subset(Fibro, PRICKLE1 < 3)
#Fibro <- subset(Fibro, MECOM < 3) #!!!! remove Endo and Mural cells first
#Fibro <- subset(Fibro, VWF < 3) #!!!! remove Endo and Mural cells first
#Fibro <- subset(Fibro, TAGLN < 3) #!!!! remove Endo and Mural cells first
#Fibro <- subset(Fibro, MYH11 < 3) #!!!! remove Endo and Mural cells first
#Fibro <- subset(Fibro, NOTCH3 < 3) #!!!! remove Endo and Mural cells first
Fibro <- subset(Fibro, MAP2 < 3)
Fibro <- subset(Fibro, CLU < 3)
Fibro <- subset(Fibro, AQP4 < 3)
Fibro <- subset(Fibro, MOBP < 3)
Fibro <- subset(Fibro, MBP < 3)
Fibro <- subset(Fibro, CD68 < 3)
Fibro <- subset(Fibro, F13A1 < 3)
Fibro <- subset(Fibro, STAB1 < 3)
Fibro <- subset(Fibro, CD163 < 3)
Fibro <- subset(Fibro, LYVE1 < 3)
Fibro <- subset(Fibro, P2RY12 < 3)
Fibro <- subset(Fibro, CD3E < 3)

###Supervised clustering using RCTD
####load genes
#rename genes
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))

#load gene sets
Epithelial <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Epithelial/Spreadsheet__Epithelial_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Epithelial <- str_replace(Epithelial, pairs$HG38[i], pairs$CosMX[i])}
Epithelial <- Epithelial[Epithelial %in% rownames(Fib)]

Epi_1_genes <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Epi_1/Spreadsheet__Epi_1_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Epi_1_genes <- str_replace(Epi_1_genes, pairs$HG38[i], pairs$CosMX[i])}
Epi_1_genes <- Epi_1_genes[Epi_1_genes %in% rownames(Fib)]

Epi_2_genes <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Epi_2/Spreadsheet__Epi_2_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Epi_2_genes <- str_replace(Epi_2_genes, pairs$HG38[i], pairs$CosMX[i])}
Epi_2_genes <- Epi_2_genes[Epi_2_genes %in% rownames(Fib)]

Epi_3_genes <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Epi_3/Spreadsheet__Epi_3_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Epi_3_genes <- str_replace(Epi_3_genes, pairs$HG38[i], pairs$CosMX[i])}
Epi_3_genes <- Epi_3_genes[Epi_3_genes %in% rownames(Fib)]

Epi_4_genes <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Epi_4/Spreadsheet__Epi_4_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Epi_4_genes <- str_replace(Epi_4_genes, pairs$HG38[i], pairs$CosMX[i])}
Epi_4_genes <- Epi_4_genes[Epi_4_genes %in% rownames(Fib)]

Fibroblast <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Fibroblast/Spreadsheet__Fibroblast_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Fibroblast <- str_replace(Fibroblast, pairs$HG38[i], pairs$CosMX[i])}
Fibroblast <- Fibroblast[Fibroblast %in% rownames(Fib)]

Fib_1_genes <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Fib_1/Spreadsheet__Fib_1_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Fib_1_genes <- str_replace(Fib_1_genes, pairs$HG38[i], pairs$CosMX[i])}
Fib_1_genes <- Fib_1_genes[Fib_1_genes %in% rownames(Fib)]

Fib_2_genes <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Fib_2/Spreadsheet__Fib_2_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Fib_2_genes <- str_replace(Fib_2_genes, pairs$HG38[i], pairs$CosMX[i])}
Fib_2_genes <- Fib_2_genes[Fib_2_genes %in% rownames(Fib)]

Fib_3_genes <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Fib_3/Spreadsheet__Fib_3_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Fib_3_genes <- str_replace(Fib_3_genes, pairs$HG38[i], pairs$CosMX[i])}
Fib_3_genes <- Fib_3_genes[Fib_3_genes %in% rownames(Fib)]

Fib_4_genes <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Fib_4/Spreadsheet__Fib_4_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Fib_4_genes <- str_replace(Fib_4_genes, pairs$HG38[i], pairs$CosMX[i])}
Fib_4_genes <- Fib_4_genes[Fib_4_genes %in% rownames(Fib)]

Endo_1 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Endo_1/Spreadsheet__Endo_1_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Endo_1 <- str_replace(Endo_1, pairs$HG38[i], pairs$CosMX[i])}
Endo_1 <- Endo_1[Endo_1 %in% rownames(Fib)]

Endo_3 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Endo_3/Spreadsheet__Endo_3_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Endo_3 <- str_replace(Endo_3, pairs$HG38[i], pairs$CosMX[i])}
Endo_3 <- Endo_3[Endo_3 %in% rownames(Fib)]

Endo_2 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Endo_2/Spreadsheet__Endo_2_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Endo_2 <- str_replace(Endo_2, pairs$HG38[i], pairs$CosMX[i])}
Endo_2 <- Endo_2[Endo_2 %in% rownames(Fib)]

Mural_2 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Mural_2/Spreadsheet__Mural_2_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Mural_2 <- str_replace(Mural_2, pairs$HG38[i], pairs$CosMX[i])}
Mural_2 <- Mural_2[Mural_2 %in% rownames(Fib)]

Mural_3 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Mural_3/Spreadsheet__Mural_3_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Mural_3 <- str_replace(Mural_3, pairs$HG38[i], pairs$CosMX[i])}
Mural_3 <- Mural_3[Mural_3 %in% rownames(Fib)]

Mural_1 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Mural_1/Spreadsheet__Mural_1_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Mural_1 <- str_replace(Mural_1, pairs$HG38[i], pairs$CosMX[i])}
Mural_1 <- Mural_1[Mural_1 %in% rownames(Fib)]

Immune <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Immune/Spreadsheet__Immune_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Immune <- str_replace(Immune, pairs$HG38[i], pairs$CosMX[i])}
Immune <- Immune[Immune %in% rownames(Fib)]

BAM_1 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/BAM_1/Spreadsheet__BAM_1_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {BAM_1 <- str_replace(BAM_1, pairs$HG38[i], pairs$CosMX[i])}
BAM_1 <- BAM_1[BAM_1 %in% rownames(Fib)]

BAM_2 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/BAM_2/Spreadsheet__BAM_2_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {BAM_2 <- str_replace(BAM_2, pairs$HG38[i], pairs$CosMX[i])}
BAM_2 <- BAM_2[BAM_2 %in% rownames(Fib)]

BAM_3 <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/BAM_3/Spreadsheet__BAM_3_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {BAM_3 <- str_replace(BAM_3, pairs$HG38[i], pairs$CosMX[i])}
BAM_3 <- BAM_3[BAM_3 %in% rownames(Fib)]

T_Cell <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/T_Cell/Spreadsheet__T_Cell_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {T_Cell <- str_replace(T_Cell, pairs$HG38[i], pairs$CosMX[i])}
T_Cell <- T_Cell[T_Cell %in% rownames(Fib)]

Parenchyma <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Parenchyma/Spreadsheet__Parenchyma_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Parenchyma <- str_replace(Parenchyma, pairs$HG38[i], pairs$CosMX[i])}
Parenchyma <- Parenchyma[Parenchyma %in% rownames(Fib)]

Oligo <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Oligo/Spreadsheet__Oligo_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Oligo <- str_replace(Oligo, pairs$HG38[i], pairs$CosMX[i])}
Oligo <- Oligo[Oligo %in% rownames(Fib)]

Astro <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Astro/Spreadsheet__Astro_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Astro <- str_replace(Astro, pairs$HG38[i], pairs$CosMX[i])}
Astro <- Astro[Astro %in% rownames(Fib)]

Neuro <- c(read_csv(paste0(curr_dir, "/CellType/sn-limma/Neuro/Spreadsheet__Neuro_corrected.csv"))$gene)
for (i in 1:length(pairs$HG38)) {Neuro <- str_replace(Neuro, pairs$HG38[i], pairs$CosMX[i])}
Neuro <- Neuro[Neuro %in% rownames(Fib)]


####Prep and run
Fib <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Fibro.RDS")
Fibro_SN <- readRDS("/Ex_ChP/ref_sets/Fibro.rds")

Ex_ChP_RNA <- Fibro_SN@assays$soupXcounts$counts
colnames(Ex_ChP_RNA) <- colnames(Fibro_SN@assays$soupXcounts)
Ex_ChP_RNA <- as.matrix(Ex_ChP_RNA)
Ex_ChP_meta <- Fibro_SN@meta.data

#replace rownames
Ex_rname <- rownames(Ex_ChP_RNA)
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))

for (i in 1:length(pairs$HG38)) {
  Ex_rname <- str_replace(Ex_rname, pairs$HG38[i], pairs$CosMX[i])
}
rownames(Ex_ChP_RNA) <- Ex_rname
Ex_ChP_RNA <- Ex_ChP_RNA[rownames(Ex_ChP_RNA) %in% rownames(Fib), , drop = FALSE]

cell_types <- Ex_ChP_meta$Subcelltype
names(cell_types) <- rownames(Ex_ChP_meta) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- Ex_ChP_meta$nCount_soupXcounts;
names(nUMI) <- rownames(Ex_ChP_meta) # create nUMI named list

reference <- Reference(Ex_ChP_RNA, cell_types, n_max_cells=50000, require_int=F)
rm(Fibro_SN)
gc()

#query object
cellcoords <- rbind(Fib@images[["Slide.1"]]@boundaries[["centroids"]]@coords+5000, Fib@images[["Slide.2"]]@boundaries[["centroids"]]@coords+100000)
cellcoords <- data.frame(cellcoords)
rownames(cellcoords) <- colnames(Fib)
colnames(cellcoords) <- c("xcoord", "ycoord")
counts <- Fib@assays$Nanostring$counts
puck <- SpatialRNA(coords=cellcoords, counts = counts, nUMI=colSums(counts))

RCTD <- create.RCTD(puck, reference, max_cores=20)
RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet') #Hangs at choose sigma choose your own

results <- RCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD@spatialRNA
resultsdir <- paste0(curr_dir, "/CellType/LT_RCTD")
dir.create(resultsdir)

plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
RCTD_call <- data.frame(results[["results_df"]])
tabyl(RCTD_call, "spot_class")

#for(i in names(Fib@images)){Fib@images[[i]]<- NULL}
Fib_meta <- Fib@meta.data
Fib_meta <- rownames_to_column(Fib_meta, "rowname")
merged_meta <- rownames_to_column(RCTD_call, "rowname")
Fib_meta <- left_join(Fib_meta, merged_meta, by="rowname")
Fib_meta <- column_to_rownames(Fib_meta, "rowname")
Fib_meta[is.na(Fib_meta)==TRUE] <- "Other"
Fib@meta.data <- Fib_meta
#Fib <- subset(Fib, spot_class=="singlet" | spot_class=="doublet_uncertain") #too stringent

merge(Fib1, merge(Fib2, Fib3))


####bubble plot

highlights_genes <- c("DCN", "LEPR", "FN1", "CDH11", "STXBP6", "FBLN1", "KCNMA1", "SLC4A4", "PLXNA4", "ALPL", "NRG1", "ROBO1", "TRPM3", "MYRIP", "GPM6A", "CLDN11", "TJP1", "TJP2", "ECM2", "SLC6A20", "SLC13A3", "LAMA4", "ABCA9", "LAMA2", "COL4A2", "COL3A1", "COL6A3", "SLC4A7", "SLC2A3")
highlights_genes <- highlights_genes[highlights_genes %in% rownames(Fib)]

Idents(Fib) <- Fib$first_type
unique(Fib$first_type)
levels(Fib) <- c("Fib_1", "Fib_2", "Fib_3")
Idents(Fib) -> Fib$first_type
plts <- DotPlot(Fib, features = highlights_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/LT_RCTD/bubble_Fib_AB_all.pdf", sep=""), width=3+length(highlights_genes)*.15, height=1.5+length(unique(Fib@meta.data[["first_type"]]))*.4)


md <- data.frame(Fib@active.ident, Fib@meta.data[["first_type"]])
colnames(md) <- c("active.ident", "group.by")
tabyl(md, active.ident)



####correlation heatmap
Ex_ChP_NS <- CreateSeuratObject(Ex_ChP_RNA, meta.data = Ex_ChP_meta)
NFSP_fun(Ex_ChP_NS, assay="RNA")

my_genes <- unique(c(Fib_1_genes, Fib_2_genes, Fib_3_genes))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP_NS, group.by = "Subcelltype", assay = "RNA", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Fib_1_SN", "Fib_2_SN", "Fib_3_SN")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Cos_cts <- data.frame(AggregateExpression(Fib, group.by = "first_type", assay = "Nanostring", slot = "data", return.seurat = FALSE))
colnames(Cos_cts) <- c("Fib_1_ST", "Fib_2_ST", "Fib_3_ST")
Cos_cts <- Cos_cts[,c("Fib_1_ST", "Fib_2_ST", "Fib_3_ST")]
Cos_cts <- data.frame(t(scale(t(Cos_cts))))
Cos_cts <- Cos_cts %>% rownames_to_column("gene")
Cos_cts <- Cos_cts[Cos_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, Cos_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Fib_1_ST", "Fib_1_SN", "Fib_2_ST", "Fib_2_SN", "Fib_3_ST", "Fib_3_SN")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("ST", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/LT_RCTD/Correlation_heatmap_Fib_all.pdf", sep = ""), width=7, height=7)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Fib.pdf", sep = ""), width=7, height=7)

Fib@meta.data$Subcelltype <- Fib@meta.data$first_type
Fib@meta.data$spot_class <- NULL
Fib@meta.data$first_type <- NULL
Fib@meta.data$second_type <- NULL
Fib@meta.data$first_class <- NULL
Fib@meta.data$second_class <- NULL
Fib@meta.data$min_score <- NULL
Fib@meta.data$singlet_score <- NULL
Fib@meta.data$conv_all <- NULL
Fib@meta.data$conv_doublet <- NULL
saveRDS(Fib, paste0(curr_dir, "/CellType/Finalcelltype/Fib_123_LT_RCTD.RDS"))



##Endo
###Remove doublets
#Remove cells that don't belong
Endo <- subset(Endo, HTR2C < 3)
#Endo <- subset(Endo, TTR < 3)
Endo <- subset(Endo, AQP1 < 3)
Endo <- subset(Endo, FOXJ1 < 3)
Endo <- subset(Endo, ARL13B < 3)
Endo <- subset(Endo, TRPM3 < 3)
Endo <- subset(Endo, SLC4A4 < 3)
Endo <- subset(Endo, LAMA2 < 3)
Endo <- subset(Endo, LAMA4 < 3)
Endo <- subset(Endo, PRICKLE1 < 3)
#Endo <- subset(Endo, MECOM < 3)
#Endo <- subset(Endo, VWF < 3)
Endo <- subset(Endo, TAGLN < 3)
Endo <- subset(Endo, MYH11 < 3)
Endo <- subset(Endo, NOTCH3 < 3)
Endo <- subset(Endo, MAP2 < 3)
#Endo <- subset(Endo, CLU < 3)
#Endo <- subset(Endo, AQP4 < 3)
Endo <- subset(Endo, MOBP < 3)
Endo <- subset(Endo, MBP < 3)
Endo <- subset(Endo, CD68 < 3)
Endo <- subset(Endo, F13A1 < 3)
Endo <- subset(Endo, STAB1 < 3)
Endo <- subset(Endo, CD163 < 3)
Endo <- subset(Endo, LYVE1 < 3)
Endo <- subset(Endo, P2RY12 < 3)
Endo <- subset(Endo, CD3E < 3)


###preprocess

Endo <- readRDS("/Ex_ChP/CosMX2/CellType/Secondcelltype/Endo.RDS")
SCTPCAfun(Endo, assay="Nanostring")

use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Endo <- RunHarmony(object = Endo, group.by.vars =  c("Slide", "projID"), dims.use = use_dims, assay="SCT", max.iter.harmony = 50, plot_convergence = TRUE) 
CNDRfun(Endo, PCA_dims = use_dims, assay="SCT", reduction = "harmony")

#Pick res
res="0.3"
Idents(Endo) <- Endo@meta.data[[paste("SCT_snn_res.", res, sep="")]]
SCNDRdimfun(Endo, reduction = "harmony", res=res, split.by= c("Slide", "projID"))

Idents(Endo) <- Endo@meta.data[[paste("SCT_snn_res.", res, sep="")]]
return_markersfun(Endo, assay="SCT", features = Endo_Nanostring_var_feat)

###Supervised
####Prep and run

Endo <- readRDS("/Ex_ChP/CosMX2/CellType/Secondcelltype/Endo.RDS")
Endo_NS <- readRDS("/Ex_ChP/ref_sets/Endo.RDS")

Ex_ChP_RNA <- Endo_NS@assays$soupXcounts$counts
colnames(Ex_ChP_RNA) <- colnames(Endo_NS@assays$soupXcounts)
Ex_ChP_RNA <- as.matrix(Ex_ChP_RNA)
Ex_ChP_meta <- Endo_NS@meta.data

#replace rownames
Ex_rname <- rownames(Endo_NS@assays$soupXcounts)
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))

for (i in 1:length(pairs$HG38)) {
  Ex_rname <- str_replace(Ex_rname, pairs$HG38[i], pairs$CosMX[i])
}
rownames(Ex_ChP_RNA) <- Ex_rname
Ex_ChP_RNA <- Ex_ChP_RNA[rownames(Ex_ChP_RNA) %in% rownames(Endo), , drop = FALSE]

cell_types <- Ex_ChP_meta$Subcelltype
names(cell_types) <- rownames(Ex_ChP_meta) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- Ex_ChP_meta$nCount_soupXcounts;
names(nUMI) <- rownames(Ex_ChP_meta) # create nUMI named list

reference <- Reference(Ex_ChP_RNA, cell_types, n_max_cells=50000, require_int=F)
rm(Endo_NS)
gc()

#query object
cellcoords <- rbind(Endo@images[["Slide.1"]]@boundaries[["centroids"]]@coords+5000, Endo@images[["Slide.2"]]@boundaries[["centroids"]]@coords+100000)
cellcoords <- data.frame(cellcoords)
rownames(cellcoords) <- colnames(Endo)
colnames(cellcoords) <- c("xcoord", "ycoord")
counts <- Endo@assays$Nanostring$counts
puck <- SpatialRNA(coords=cellcoords, counts = counts, nUMI=colSums(counts))

RCTD <- create.RCTD(puck, reference, max_cores=20)
RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet') #Hangs at choose sigma choose your own

results <- RCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD@spatialRNA
resultsdir <- paste0(curr_dir, "/CellType/LT_RCTD")
dir.create(resultsdir)

plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
RCTD_call <- data.frame(results[["results_df"]])
tabyl(RCTD_call, "spot_class")

for(i in names(Endo@images)){Endo@images[[i]]<- NULL}
Endo_meta <- Endo@meta.data
Endo_meta <- rownames_to_column(Endo_meta, "rowname")
merged_meta <- rownames_to_column(RCTD_call, "rowname")
Endo_meta <- left_join(Endo_meta, merged_meta, by="rowname")
Endo_meta <- column_to_rownames(Endo_meta, "rowname")
Endo_meta[is.na(Endo_meta)==TRUE] <- "Other"
Endo@meta.data <- Endo_meta
#Endo <- subset(Endo, spot_class=="singlet" | spot_class=="doublet_uncertain") #too stringent


####bubble plot
highlights_genes <- c("PECAM1", "VWF", "INSR", "CLDN5", "CDH13", "ARL15", "MECOM", "VEGFC", "COL8A1", "FBN1", "THSD7A", 
                      "TSHZ2", "CEMIP2", "IL1R1", "IL4R", "FLI1", "SORBS2", 
                      "PLVAP", "FLT1", "ADGRF5", "KCNQ3", "SLC2A3", "ITGA1")
highlights_genes <- highlights_genes[highlights_genes %in% rownames(Endo)]

Idents(Endo) <- Endo$first_type
unique(Endo$first_type)
levels(Endo) <- c("Endo_1", "Endo_2", "Endo_3")
Idents(Endo) -> Endo$first_type
plts <- DotPlot(Endo, features = highlights_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/LT_RCTD/bubble_Endo_all.pdf", sep=""), width=3+length(highlights_genes)*.15, height=1.5+length(unique(Endo@meta.data[["first_type"]]))*.4)


md <- data.frame(Endo@active.ident, Endo@meta.data[["first_type"]])
colnames(md) <- c("active.ident", "group.by")
tabyl(md, active.ident)



####correlation heatmap
Ex_ChP_NS <- CreateSeuratObject(Ex_ChP_RNA, meta.data = Ex_ChP_meta)
NFSP_fun(Ex_ChP_NS, assay="RNA")

my_genes <- unique(c(Endo_1, Endo_2, Endo_3))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP_NS, group.by = "Subcelltype", assay = "RNA", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Endo_1_SN", "Endo_3_SN", "Endo_2_SN")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Cos_cts <- data.frame(AggregateExpression(Endo, group.by = "first_type", assay = "Nanostring", slot = "data", return.seurat = FALSE))
colnames(Cos_cts) <- c("Endo_1_ST", "Endo_2_ST", "Endo_3_ST")
Cos_cts <- Cos_cts[,c("Endo_1_ST", "Endo_2_ST", "Endo_3_ST")]
Cos_cts <- data.frame(t(scale(t(Cos_cts))))
Cos_cts <- Cos_cts %>% rownames_to_column("gene")
Cos_cts <- Cos_cts[Cos_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, Cos_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Endo_1_ST", "Endo_1_SN", "Endo_2_ST", "Endo_2_SN", "Endo_3_ST", "Endo_3_SN")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("ST", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/LT_RCTD/Correlation_heatmap_Endo_all.pdf", sep = ""), width=7, height=7)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Endo.pdf", sep = ""), width=7, height=7)

Endo@meta.data$Subcelltype <- Endo@meta.data$first_type
Endo@meta.data$spot_class <- NULL
Endo@meta.data$first_type <- NULL
Endo@meta.data$second_type <- NULL
Endo@meta.data$first_class <- NULL
Endo@meta.data$second_class <- NULL
Endo@meta.data$min_score <- NULL
Endo@meta.data$singlet_score <- NULL
Endo@meta.data$conv_all <- NULL
Endo@meta.data$conv_doublet <- NULL
saveRDS(Endo, "/Ex_ChP/CosMX2/CellType/Finalcelltype/Endo_LT_RCTD.RDS")




##Mural
###preprocess
Mural <- readRDS("/Ex_ChP/CosMX2/CellType/Secondcelltype/Mural.RDS")

SCTPCAfun(Mural, assay="Nanostring")

use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Mural <- RunHarmony(object = Mural, group.by.vars = c("Slide", "projID"), dims.use = use_dims, assay="SCT", max.iter.harmony = 50, plot_convergence = TRUE) 
CNDRfun(Mural, PCA_dims = use_dims, assay="SCT", reduction = "harmony")

#Pick res
res="0.5"
Idents(Mural) <- Mural@meta.data[[paste("SCT_snn_res.", res, sep="")]]
SCNDRdimfun(Mural, reduction = "harmony", res=res, split.by= c("Slide", "projID"))
return_markersfun(Mural, assay="SCT", features = Mural_Nanostring_var_feat)

Featurefun(Mural, c("HTR2C", "TTR", "CLIC6", "AQP1", "TUBB4B", "FOXJ1", "ARL13B", "KRT18", "COX2", "TRPM3", "VIM", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "FOXP2", "COL6A3", "SLC6A13", "LUM", "DCN", "FN1", "ACTA2", "KCNMA1", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "CD3E", "HSPE1", "HSPH1"), "Key Markers",  assay="SCT")
Featurefun(Mural, c("ENPP2","PLXDC2","PALMD", "MID1", "CEMIP", "GLIS3", "PDZD2", 
                    "COL4A4", "CRISPLD2", "ADAMTS9", "TLE1", "PIPNC1", "ADCY3", "HIPK2", "COL18A1", "SCL22A3",
                    "SOX6", "MYO1D", "MYH11", "ACTA2", "PRUNE2", "RYR2", "TPM1", "LMOD1"), "Mural", assay="SCT")


###Rename all identities

Mural <- RenameIdents(Mural, "0"="Mural_2", "1"="Mural_1", "2"="Mural_2", "3"="Mural_3") #last pass
Mural@meta.data$Subcelltype <- Mural@active.ident
DimPlot(Mural)

saveRDS(Mural, paste0(curr_dir, "/CellType/Finalcelltype/Mural.RDS"))

Ncellfun(Mural, paste0(curr_dir, "/CellType/Mural/Final per projid.csv", sep=""), "projID")


###Supervised
####Prep and run

Mural <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Mural.RDS")
Mural_SN <- readRDS("/Ex_ChP/ref_sets/Mural.RDS")

Ex_ChP_RNA <- Mural_SN@assays$soupXcounts$counts
colnames(Ex_ChP_RNA) <- colnames(Mural_SN@assays$soupXcounts)
Ex_ChP_RNA <- as.matrix(Ex_ChP_RNA)
Ex_ChP_meta <- Mural_SN@meta.data

#replace rownames
Ex_rname <- rownames(Mural_SN@assays$soupXcounts)
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))

for (i in 1:length(pairs$HG38)) {
  Ex_rname <- str_replace(Ex_rname, pairs$HG38[i], pairs$CosMX[i])
}
rownames(Ex_ChP_RNA) <- Ex_rname
Ex_ChP_RNA <- Ex_ChP_RNA[rownames(Ex_ChP_RNA) %in% rownames(Mural), , drop = FALSE]

cell_types <- Ex_ChP_meta$Subcelltype
names(cell_types) <- rownames(Ex_ChP_meta) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- Ex_ChP_meta$nCount_soupXcounts;
names(nUMI) <- rownames(Ex_ChP_meta) # create nUMI named list

reference <- Reference(Ex_ChP_RNA, cell_types, n_max_cells=50000, require_int=F)
rm(Mural_SN)
gc()

#query object
cellcoords <- rbind(Mural@images[["Slide.1"]]@boundaries[["centroids"]]@coords+5000, Mural@images[["Slide.2"]]@boundaries[["centroids"]]@coords+100000)
cellcoords <- data.frame(cellcoords)
rownames(cellcoords) <- colnames(Mural)
colnames(cellcoords) <- c("xcoord", "ycoord")
counts <- Mural@assays$Nanostring$counts
puck <- SpatialRNA(coords=cellcoords, counts = counts, nUMI=colSums(counts))

RCTD <- create.RCTD(puck, reference, max_cores=20)
RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet') #Hangs at choose sigma choose your own

results <- RCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD@spatialRNA
resultsdir <- paste0(curr_dir, "/CellType/LT_RCTD")
dir.create(resultsdir)

plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
RCTD_call <- data.frame(results[["results_df"]])
tabyl(RCTD_call, "spot_class")

for(i in names(Mural@images)){Mural@images[[i]]<- NULL}
Mural_meta <- Mural@meta.data
Mural_meta <- rownames_to_column(Mural_meta, "rowname")
merged_meta <- rownames_to_column(RCTD_call, "rowname")
Mural_meta <- left_join(Mural_meta, merged_meta, by="rowname")
Mural_meta <- column_to_rownames(Mural_meta, "rowname")
Mural_meta[is.na(Mural_meta)==TRUE] <- "Other"
Mural@meta.data <- Mural_meta
#Mural <- subset(Mural, spot_class=="singlet" | spot_class=="doublet_uncertain") #too stringent


####bubble plot
highlights_genes <- c("ENPP2","PLXDC2","PALMD", "MID1", "CEMIP", "GLIS3", "PDZD2", 
                      "COL4A4", "CRISPLD2", "ADAMTS9", "TLE1", "PIPNC1", "ADCY3", "HIPK2", "COL18A1", "SCL22A3",
                      "SOX6", "MYO1D", "MYH11", "ACTA2", "PRUNE2", "RYR2", "TPM1", "LMOD1") %>% unique()

highlights_genes <- highlights_genes[highlights_genes %in% rownames(Mural)]

Idents(Mural) <- Mural$first_type
unique(Idents(Mural))
levels(Mural) <- c("Mural_1", "Mural_2", "Mural_3")
Idents(Mural) -> Mural$first_type
plts <- DotPlot(Mural, features = highlights_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/LT_RCTD/bubble_Mural_all.pdf", sep=""), width=3+length(highlights_genes)*.15, height=1.5+length(unique(Mural@meta.data[["first_type"]]))*.4)

md <- data.frame(Mural@active.ident, Mural@meta.data[["first_type"]])
colnames(md) <- c("active.ident", "group.by")
tabyl(md, active.ident)



####correlation heatmap
Ex_ChP_NS <- CreateSeuratObject(Ex_ChP_RNA, meta.data = Ex_ChP_meta)
NFSP_fun(Ex_ChP_NS, assay="RNA")

my_genes <- unique(c(Mural_1, Mural_2, Mural_3))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP_NS, group.by = "Subcelltype", assay = "RNA", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Mural_1_SN", "Mural_3_SN", "Mural_2_SN")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Cos_cts <- data.frame(AggregateExpression(Mural, group.by = "first_type", assay = "Nanostring", slot = "data", return.seurat = FALSE))
colnames(Cos_cts) <- c("Mural_1_ST", "Mural_2_ST", "Mural_3_ST")
Cos_cts <- Cos_cts[,c("Mural_1_ST", "Mural_2_ST", "Mural_3_ST")]
Cos_cts <- data.frame(t(scale(t(Cos_cts))))
Cos_cts <- Cos_cts %>% rownames_to_column("gene")
Cos_cts <- Cos_cts[Cos_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, Cos_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Mural_1_ST", "Mural_1_SN", "Mural_2_ST", "Mural_2_SN", "Mural_3_ST", "Mural_3_SN")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("ST", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/LT_RCTD/Correlation_heatmap_Mural_all.pdf", sep = ""), width=7, height=7)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Mural.pdf", sep = ""), width=7, height=7)

Mural@meta.data$Subcelltype <- Mural@meta.data$first_type
Mural@meta.data$spot_class <- NULL
Mural@meta.data$first_type <- NULL
Mural@meta.data$second_type <- NULL
Mural@meta.data$first_class <- NULL
Mural@meta.data$second_class <- NULL
Mural@meta.data$min_score <- NULL
Mural@meta.data$singlet_score <- NULL
Mural@meta.data$conv_all <- NULL
Mural@meta.data$conv_doublet <- NULL
saveRDS(Mural, "/Ex_ChP/CosMX2/CellType/Finalcelltype/Mural_LT_RCTD.RDS")




##Imm
###preprocess
Imm <- readRDS("/Ex_ChP/CosMX2/CellType/Firstcelltype/Imm.RDS")
Imm_meta <- Imm@meta.data
Imm_mat <- Imm@assays[["SCT"]]@layers[["counts"]]
colnames(Imm_mat) <- colnames(Imm)
rownames(Imm_mat) <- rownames(Imm)
Imm <- CreateSeuratObject(counts=Imm_mat, assay="SCT", meta.data=Imm_meta)
SCTPCAfun(Imm, assay="Nanostring") #, spikein=ChP_genes) # first pass

use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Imm <- RunHarmony(object = Imm, group.by.vars = c("Slide", "projID"), dims.use = use_dims, assay="Nanostring", max.iter.harmony = 50, plot_convergence = TRUE) 
CNDRfun(Imm, PCA_dims = use_dims, assay="SCT", reduction = "harmony")

#Pick res
res="0.3"
Idents(Imm) <- Imm@meta.data[[paste("SCT_snn_res.", res, sep="")]]
SCNDRdimfun(Imm, reduction = "harmony", res=res, split.by=c("Slide", "projID"))

Idents(Imm) <- Imm@meta.data[[paste("SCT_snn_res.", res, sep="")]]
return_markersfun(Imm, assay="SCT", features = Imm_Nanostring_var_feat)


###Overlap between SN and ST

plts <- ggVennDiagram(list(Imm_Nanostring_var_feat, Imm_soupXcounts_var_feat), category.names = c("CosMX", "SN")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2png(plts, paste0(curr_dir, "/CellType/SNvCosMX/Imm_Cos_ChP_variable_features_overlap.png"))


###Remove doublets
Imm <- RenameIdents(Imm, "0"="Imm", "1"="Imm", "2"="Imm", "3"="Imm") #res=0.3
Imm@meta.data$Subcelltype <- Idents(Imm)
Imm <- subset(Imm, Subcelltype=="Doublet", invert=T)
#Remove cells that don't belong
Imm <- subset(Imm, HTR2C < 2)
Imm <- subset(Imm, TTR < 2)
Imm <- subset(Imm, AQP1 < 3)
Imm <- subset(Imm, FOXJ1 < 3)
Imm <- subset(Imm, ARL13B < 3)
Imm <- subset(Imm, TRPM3 < 3)
Imm <- subset(Imm, SLC4A4 < 3)
Imm <- subset(Imm, LAMA2 < 3)
Imm <- subset(Imm, LAMA4 < 3)
Imm <- subset(Imm, PRICKLE1 < 3)
Imm <- subset(Imm, MECOM < 3)
Imm <- subset(Imm, VWF < 3)
Imm <- subset(Imm, TAGLN < 3)
Imm <- subset(Imm, MYH11 < 3)
Imm <- subset(Imm, NOTCH3 < 3)
Imm <- subset(Imm, MAP2 < 3)
#Imm <- subset(Imm, CLU < 3)
Imm <- subset(Imm, AQP4 < 3)
Imm <- subset(Imm, MOBP < 3)
Imm <- subset(Imm, MBP < 3)
#Imm <- subset(Imm, CD68 < 3)
#Imm <- subset(Imm, F13A1 < 3)
#Imm <- subset(Imm, STAB1 < 3)
#Imm <- subset(Imm, CD163 < 3)
#Imm <- subset(Imm, LYVE1 < 3)
#Imm <- subset(Imm, P2RY12 < 3)
#Imm <- subset(Imm, CD3E < 3)


###Rename all identities

Imm <- RenameIdents(Imm, "0"="Immune", "1"="Immune", "2"="Immune", "3"="Immune") #res=0.3
Imm@meta.data$Subcelltype <- Idents(Imm)
saveRDS(Imm, "/Ex_ChP/CosMX2/CellType/Finalcelltype/Imm.RDS")


###Supervised clustering
Imm <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Imm.RDS")
Imm_NS <- readRDS("/Ex_ChP/ref_sets/Imm_SCT.rds")

Ex_ChP_RNA <- Imm_NS@assays$soupXcounts$counts
colnames(Ex_ChP_RNA) <- colnames(Imm_NS@assays$soupXcounts)
Ex_ChP_RNA <- as.matrix(Ex_ChP_RNA)
Ex_ChP_meta <- Imm_NS@meta.data

#replace rownames
Ex_rname <- rownames(Imm_NS@assays$soupXcounts)
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))

for (i in 1:length(pairs$HG38)) {
  Ex_rname <- str_replace(Ex_rname, pairs$HG38[i], pairs$CosMX[i])
}
rownames(Ex_ChP_RNA) <- Ex_rname
Ex_ChP_RNA <- Ex_ChP_RNA[rownames(Ex_ChP_RNA) %in% rownames(Imm), , drop = FALSE]

cell_types <- Ex_ChP_meta$Subcelltype
names(cell_types) <- rownames(Ex_ChP_meta) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- Ex_ChP_meta$nCount_soupXcounts;
names(nUMI) <- rownames(Ex_ChP_meta) # create nUMI named list

reference <- Reference(Ex_ChP_RNA, cell_types, n_max_cells=50000, require_int=F)
rm(Imm_NS)
gc()

#query object
cellcoords <- rbind(Imm@images[["Slide_1"]]@boundaries[["centroids"]]@coords+5000, Imm@images[["Slide_2"]]@boundaries[["centroids"]]@coords+100000)
cellcoords <- data.frame(cellcoords)
rownames(cellcoords) <- colnames(Imm)
colnames(cellcoords) <- c("xcoord", "ycoord")
counts <- Imm@assays$Nanostring$counts
puck <- SpatialRNA(coords=cellcoords, counts = counts, nUMI=colSums(counts))

RCTD <- create.RCTD(puck, reference, max_cores=20)
RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet') #Hangs at choose sigma choose your own

results <- RCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD@spatialRNA
resultsdir <- paste0(curr_dir, "/CellType/LT_RCTD")
dir.create(resultsdir)

plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
RCTD_call <- data.frame(results[["results_df"]])
tabyl(RCTD_call, "spot_class")

for(i in names(Imm@images)){Imm@images[[i]]<- NULL}
Imm_meta <- Imm@meta.data
Imm_meta <- rownames_to_column(Imm_meta, "rowname")
merged_meta <- rownames_to_column(RCTD_call, "rowname")
Imm_meta <- left_join(Imm_meta, merged_meta, by="rowname")
Imm_meta <- column_to_rownames(Imm_meta, "rowname")
Imm_meta[is.na(Imm_meta)==TRUE] <- "Other"
Imm@meta.data <- Imm_meta
#Imm <- subset(Imm, spot_class=="singlet" | spot_class=="doublet_uncertain")


####bubble plot
highlights_genes <- c("MSR1", "F13A1", "LYVE1", "MRC1", "STAB1", "CD163", "MERTK", "HCK", "IL4R", "CD74", "HLA-DRA", "CSF1R", "JAK2", "CD86", "CD83", "DNAJA1", "HSPE1", "HSPD1", "HSPB1", "MKI67", "CENPF", "TOP2A", "ITGAM", "S100A8", "CCR2", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11", "THEMIS", "CD3E", "PDE7A", "CAMK4", "HLA-B", "CD8A")
highlights_genes <- highlights_genes[highlights_genes %in% rownames(Imm)]

Idents(Imm) <- Imm$first_type
unique(Imm$first_type)
levels(Imm) <- c("BAM_1", "BAM_2", "BAM_3", "T_Cell")
Idents(Imm) -> Imm$first_type
plts <- DotPlot(Imm, features = highlights_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/LT_RCTD/bubble_Imm_all.pdf", sep=""), width=3+length(highlights_genes)*.15, height=1.5+length(unique(Imm@meta.data[["first_type"]]))*.4)


md <- data.frame(Imm@active.ident, Imm@meta.data[["first_type"]])
colnames(md) <- c("active.ident", "group.by")
tabyl(md, active.ident)

####correlation heatmap
Ex_ChP_NS <- CreateSeuratObject(Ex_ChP_RNA, meta.data = Ex_ChP_meta)
NFSP_fun(Ex_ChP_NS, assay="RNA")

my_genes <- unique(c(Epi_1_genes, Epi_2_genes, Epi_3_genes, Epi_4_genes, Fib_1_genes, Fib_2_genes, Fib_3_genes, Fib_4_genes, Endo, Mural, BAM_1, BAM_2, BAM_3, T_Cell, Oligo, Astro, Neuro))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP_NS, group.by = "Subcelltype", assay = "RNA", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("BAM_1_SN", "BAM_2_SN", "BAM_3_SN", "T_Cell_SN")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Cos_cts <- data.frame(AggregateExpression(Imm, group.by = "first_type", assay = "Nanostring", slot = "data", return.seurat = FALSE))
colnames(Cos_cts) <- c("BAM_1_ST", "BAM_2_ST", "BAM_3_ST", "T_Cell_ST")
Cos_cts <- Cos_cts[,c("BAM_1_ST", "BAM_2_ST", "BAM_3_ST", "T_Cell_ST")]
Cos_cts <- data.frame(t(scale(t(Cos_cts))))
Cos_cts <- Cos_cts %>% rownames_to_column("gene")
Cos_cts <- Cos_cts[Cos_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, Cos_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c( "BAM_1_ST", "BAM_1_SN", "BAM_2_ST", "BAM_2_SN", "BAM_3_ST", "BAM_3_SN", "T_Cell_ST", "T_Cell_SN")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("ST", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/LT_RCTD/Correlation_heatmap_Imm_all.pdf", sep = ""), width=7, height=7)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Imm.pdf", sep = ""), width=7, height=7)

#clean up metadata so it's consistent prior to save
Imm@meta.data$Subcelltype <- Imm@meta.data$first_type
Imm@meta.data$spot_class <- NULL
Imm@meta.data$first_type <- NULL
Imm@meta.data$second_type <- NULL
Imm@meta.data$first_class <- NULL
Imm@meta.data$second_class <- NULL
Imm@meta.data$min_score <- NULL
Imm@meta.data$singlet_score <- NULL
Imm@meta.data$conv_all <- NULL
Imm@meta.data$conv_doublet <- NULL
saveRDS(Imm, "/Ex_ChP/CosMX2/CellType/Finalcelltype/Imm_LT_RCTD.RDS")



##Paren
###preprocess
Paren <- readRDS("/Ex_ChP/CosMX2/CellType/Firstcelltype/Paren.RDS")
Paren_meta <- Paren@meta.data
Paren_mat <- Paren@assays[["Nanostring"]]@layers[["counts"]]
colnames(Paren_mat) <- colnames(Paren)
rownames(Paren_mat) <- rownames(Paren)
Paren <- CreateSeuratObject(counts=Paren_mat, assay="Nanostring", meta.data=Paren_meta)
SCTPCAfun(Paren, assay="Nanostring") #, spikein=ChP_genes # firstpass

use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Paren <- RunHarmony(object = Paren, group.by.vars = c("Slide", "projID"), dims.use = use_dims, assay="SCT", max.iter.harmony = 50, plot_convergence = TRUE) 
CNDRfun(Paren, PCA_dims = use_dims, assay="SCT", reduction = "harmony")

#Pick res
res="0.3"
Idents(Paren) <- Paren@meta.data[[paste("SCT_snn_res.", res, sep="")]]
SCNDRdimfun(Paren, reduction = "harmony", res=res, split.by=c("Slide", "projID"))

#Reference markers
Idents(Paren) <- Paren@meta.data[["predicted.Subcelltype"]]
plts <- DimPlot(Paren, raster=T, reduction='umap')
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Paren/ctrl_Ex_ChP_Submarkers.png"), plot=plts)

Idents(Paren) <- Paren@meta.data[[paste("SCT_snn_res.", res, sep="")]]
return_markersfun(Paren, assay="SCT", features = Paren_Nanostring_var_feat)

Featurefun(Paren, c("HTR2C", "TTR", "CLIC6", "AQP1", "TUBB4B", "FOXJ1", "ARL13B", "KRT18", "COX2", "TRPM3", "VIM", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "FOXP2", "COL6A3", "SLC6A13", "LUM", "DCN", "FN1", "ACTA2", "KCNMA1", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "CD3E", "HSPE1", "HSPH1"), "Key Markers",  assay="SCT")


###clusters

curr_clust <- 1
Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_BAM, known_peric, known_Mural_1), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Paren_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Paren, curr_clust)

search_markersfun(Paren, curr_clust, c(known_Mesen, known_Paren, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Paren, curr_clust, known_Paren_prog)

search_markersfun(Paren, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Paren, curr_clust, c(known_BAM, known_peric, known_Mural_1) %>% unique())

search_markersfun(Paren, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Paren, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Paren, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Paren, curr_clust, c(known_light) %>% unique())
search_markersfun(Paren, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Paren, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Paren, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Paren, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Paren, curr_clust, known_dividing %>% unique())

search_markersfun(Paren, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Paren, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Paren, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

search_markersfun(Paren, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Paren, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Paren, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Paren, "LAMA4", raster=F)



###Remove doublets

Paren <- RenameIdents(Paren, "0"="Olig", "1"="Olig", "2"="Doublet", "3"="Doublet") #res=0.3
Paren@meta.data$Subcelltype <- Idents(Paren)
Paren <- subset(Paren, Subcelltype=="Doublet", invert=T)
#Remove cells that don't belong
Paren <- subset(Paren, HTR2C < 3)
Paren <- subset(Paren, TTR < 3)
Paren <- subset(Paren, AQP1 < 3)
Paren <- subset(Paren, FOXJ1 < 3)
Paren <- subset(Paren, ARL13B < 3)
Paren <- subset(Paren, TRPM3 < 3)
Paren <- subset(Paren, SLC4A4 < 3)
Paren <- subset(Paren, LAMA2 < 3)
Paren <- subset(Paren, LAMA4 < 3)
Paren <- subset(Paren, PRICKLE1 < 3)
Paren <- subset(Paren, SLC6A3 < 3)
Paren <- subset(Paren, MECOM < 3)
Paren <- subset(Paren, VWF < 3)
Paren <- subset(Paren, TAGLN < 3)
Paren <- subset(Paren, MYH11 < 3)
Paren <- subset(Paren, NOTCH3 < 3)
#Paren <- subset(Paren, CADM2 < 3)
#Paren <- subset(Paren, MAP2 < 3)
#Paren <- subset(Paren, CLU < 3)
#Paren <- subset(Paren, AQP4 < 3)
#Paren <- subset(Paren, MOBP < 3)
#Paren <- subset(Paren, MBP < 3)
Paren <- subset(Paren, CD68 < 3)
Paren <- subset(Paren, F13A1 < 3)
Paren <- subset(Paren, STAB1 < 3)
Paren <- subset(Paren, CD163 < 3)
Paren <- subset(Paren, LYVE1 < 3)
Paren <- subset(Paren, P2RY12 < 3)
Paren <- subset(Paren, CD3E < 3)


###Rename all identities

saveRDS(Paren, "/Ex_ChP/CosMX2/CellType/Secondcelltype/Paren.RDS")



##Wholeset
###merge

Epi <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Epi.RDS")
Fibro <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Fib_LT_RCTD.RDS")
Endo <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Endo.RDS")
Mural <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Mural.RDS")
Imm <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Imm_LT_RCTD.RDS")
Paren <- readRDS("/Ex_ChP/CosMX2/CellType/Finalcelltype/Paren.RDS")
Cos_ChP <- readRDS("/Ex_ChP/CosMX2/CellType/Firstcelltype/Cos_ChP.RDS")

Epi_meta <- Epi@meta.data
Fibro_meta <- Fibro@meta.data
Endo_meta <- Endo@meta.data
Mural_meta <- Mural@meta.data
Paren_meta <- Paren@meta.data
Imm_meta <- Imm@meta.data
#rm(Epi, Fibro, Endo, Mural, Paren, Imm)

merged_meta <- rbind(Epi_meta, Fibro_meta, Endo_meta, Mural_meta, Paren_meta, Imm_meta)
#rm(Epi_meta, Fibro_meta, Endo_meta, Mural_meta, Paren_meta, Imm_meta)

Cos_ChP@meta.data$Subcelltype <- NULL
Cos_ChP_meta <- Cos_ChP@meta.data

Cos_ChP_meta <- rownames_to_column(Cos_ChP_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Cos_ChP_meta <- left_join(Cos_ChP_meta, merged_meta, by="rowname")
Cos_ChP_meta <- column_to_rownames(Cos_ChP_meta, "rowname")
Cos_ChP_meta$Subcelltype <- as.character(Cos_ChP_meta$Subcelltype)
Cos_ChP_meta[is.na(Cos_ChP_meta)==TRUE] <- "Doublet"
Cos_ChP@meta.data <- Cos_ChP_meta
Cos_ChP <- subset(Cos_ChP, Subcelltype=="Doublet", invert=T)

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
DimPlot(Cos_ChP)


###rerun

SCTPCAfun(Cos_ChP, assay="Nanostring")
use_dims <- 1:15
set.seed(1234)
Cos_ChP <- RunHarmony(object = Cos_ChP, group.by.vars = c("Slide", "projID"), dims.use = use_dims, assay="SCT", max.iter.harmony = 50, plot_convergence = TRUE) 
CNDRfun_quick(Cos_ChP, PCA_dims = use_dims, assay="SCT", reduction = "harmony", resolution=0.1)
res="0.1"
Idents(Cos_ChP) <- Cos_ChP@meta.data[[paste("SCT_snn_res.", res, sep="")]]
SCNDRdimfun(Cos_ChP, reduction = "harmony", res=res, split.by=c("Slide", "projID"))

Idents(Cos_ChP) <- Cos_ChP@meta.data[["Subcelltype"]]
Cos_ChP  <- RenameIdents(Cos_ChP, "Epi_1"="Epithelial", "Epi_2"="Epithelial", "Epi_3"="Epithelial", "Epi_4"="Epithelial", "Fib_1"="Fibroblast", "Fib_2"="Fibroblast", "Fib_3"="Fibroblast", "Fib_4"="Fibroblast", "Endo_1"="Endothelial", "Endo_2"="Endothelial", "Endo_3"="Endothelial", "Mural_2"="Mural", "Mural_3"="Mural", "Mural_1"="Mural", "BAM_1"="Immune", "BAM_2"="Immune", "BAM_3"="Immune", "T_Cell"="Immune")
levels(Cos_ChP) <- c("Epithelial", "Fibroblast", "Endothelial", "Mural", "Immune", "Parenchyma")
Cos_ChP@meta.data[["Majorcelltype"]] <- Idents(Cos_ChP)
DimPlot(Cos_ChP, raster=T)


###Umap vs nearest neighbour cleanup

Cos_ChP_hold <- Cos_ChP
for(i in names(Cos_ChP@images)){Cos_ChP@images[[i]]<- NULL}

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
Cos_ChP <- CellSelector(DimPlot(Cos_ChP), Cos_ChP, ident = "Epi")
Epi <- subset(Cos_ChP, ident="Epi")
Epi <- subset(Epi, c(Subcelltype=="Epi_1" | Subcelltype=="Epi_2" | Subcelltype=="Epi_3" | Subcelltype=="Epi_4"))

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
Cos_ChP <- CellSelector(DimPlot(Cos_ChP), Cos_ChP, ident = "Fibro")
Fibro <- subset(Cos_ChP, ident="Fibro")
Fibro <- subset(Fibro, Subcelltype=="Fibroblast")

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
Cos_ChP <- CellSelector(DimPlot(Cos_ChP), Cos_ChP, ident = "End")
Endo <- subset(Cos_ChP, ident="End")
Endo <- subset(Cos_ChP, Subcelltype=="Endo")

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
Cos_ChP <- CellSelector(DimPlot(Cos_ChP), Cos_ChP, ident = "Mur")
Mural <- subset(Cos_ChP, ident="Mur")
Mural <- subset(Mural, Subcelltype=="Mural")

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
Cos_ChP <- CellSelector(DimPlot(Cos_ChP), Cos_ChP, ident = "BAMs")
Imm <- subset(Cos_ChP, ident="BAMs")
Imm <- subset(Imm, Subcelltype=="Immune")

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
Cos_ChP <- CellSelector(DimPlot(Cos_ChP), Cos_ChP, ident = "Oligo")
Oligo <- subset(Cos_ChP, ident="Oligo")
Oligo <- subset(Oligo, Subcelltype=="Olig")

Epi_meta <- Epi@meta.data
Fibro_meta <- Fibro@meta.data
Endo_meta <- Endo@meta.data
Mural_meta <- Mural@meta.data
Imm_meta <- Imm@meta.data
Oligo_meta <- Oligo@meta.data
merged_meta <- rbind(Epi_meta, Fibro_meta, Endo_meta, Mural_meta, Imm_meta, Oligo_meta)

Cos_ChP@meta.data$Subcelltype <- NULL
Cos_ChP_meta <- Cos_ChP@meta.data

Cos_ChP_meta <- rownames_to_column(Cos_ChP_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Cos_ChP_meta <- left_join(Cos_ChP_meta, merged_meta, by="rowname")
Cos_ChP_meta <- column_to_rownames(Cos_ChP_meta, "rowname")
Cos_ChP_meta$Subcelltype <- as.character(Cos_ChP_meta$Subcelltype)
Cos_ChP_meta[is.na(Cos_ChP_meta)==TRUE] <- "Doublet"
Cos_ChP@meta.data <- Cos_ChP_meta
Cos_ChP <- subset(Cos_ChP, Subcelltype=="Doublet", invert=T)

#Reorder
#Reorder
Idents(Cos_ChP) <- Cos_ChP@meta.data[["Subcelltype"]]
levels(Cos_ChP) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b", "Fib_1", "Fib_2", "Fib_3", "Endo_1", "Endo_2", "Endo_3", "Mural_1", "Mural_2", "Mural_3", "BAM_1", "BAM_2", "BAM_3", "T_Cell", "Parenchyma")
Cos_ChP@meta.data[["Subcelltype"]] <- Idents(Cos_ChP)
DimPlot(Cos_ChP)
Idents(Cos_ChP) <- Cos_ChP@meta.data[["Majorcelltype"]]
levels(Cos_ChP) <- c("Epithelial", "Fibroblast", "Endothelial", "Mural", "Immune", "Parenchyma")
Cos_ChP@meta.data[["Majorcelltype"]] <- Idents(Cos_ChP)
DimPlot(Cos_ChP)

Cos_ChP_meta <- Cos_ChP@meta.data

Cos_ChP_hold@meta.data$Subcelltype <- NULL
Cos_ChP_meta_hold <- Cos_ChP_hold@meta.data

Cos_ChP_meta_hold <- rownames_to_column(Cos_ChP_meta_hold, "rowname")
Cos_ChP_meta <- rownames_to_column(Cos_ChP_meta, "rowname")
Cos_ChP_meta <- data.frame(Cos_ChP_meta$rowname, Cos_ChP_meta$Subcelltype)
colnames(Cos_ChP_meta) <- c("rowname", "Subcelltype")
Cos_ChP_meta_hold <- left_join(Cos_ChP_meta_hold, Cos_ChP_meta, by="rowname")
Cos_ChP_meta_hold <- column_to_rownames(Cos_ChP_meta_hold, "rowname")
Cos_ChP_meta_hold$Subcelltype <- as.character(Cos_ChP_meta_hold$Subcelltype)
Cos_ChP_meta_hold[is.na(Cos_ChP_meta_hold)==TRUE] <- "Doublet"
Cos_ChP_hold@meta.data <- Cos_ChP_meta_hold
Cos_ChP_hold <- subset(Cos_ChP_hold, Subcelltype=="Doublet", invert=T)

DimPlot(Cos_ChP_hold)
Cos_ChP <- Cos_ChP_hold
rm(Cos_ChP_hold, Endo, Endo_, Epi, Fibro, BAM, Olig, Mural)
gc()

saveRDS(Cos_ChP, paste0(curr_dir, "/CellType/Finalcelltype_images/Cos_ChP_images.RDS"))
