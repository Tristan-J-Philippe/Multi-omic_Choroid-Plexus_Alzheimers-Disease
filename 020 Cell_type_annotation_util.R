###
### 020 Cell_type_annotation_util.R
###
# Purpose: Cell type annotation and cleaning functions
# Dependencies:
library(data.table)
library(plyr)
library(dplyr)
library(Seurat)
library(Signac)
library(future)
library(ggplot2)
library(Matrix)
library(harmony)
library(tidyverse)
library(ggVennDiagram)
library(janitor)
library(cowplot)
library(clustree)
library(igraph)
library(igraphdata)
library(export)

curr_dir <- getwd() #Get current directory; dependency for most functions to output in folders

#SCTPCA
SCTPCAfun <- function(srt, assay="Nanostring", spikein=NULL){
  #Performs SCTransform, PCA, and regular plots.
  #srt: seurat object
  #assay: Assay, default is RNA
  
  nameobj <- deparse(substitute(srt)) #make char
  nameass <- deparse(substitute(assay)) #make char
  path <- paste(curr_dir, "/Sample_preprocess/", nameobj, "/Norm_Clustering/", sep="") #specify dir
  if(dir.exists(path)==FALSE){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  DefaultAssay(srt) <- assay
  
  srt <- SCTransform(srt, assay=assay, return.only.var.genes = FALSE, variable.features.n = 2000, min_cells = 0.5) #!!!!setting return.only.var.genes = FALSE and min_cells = 1 makes it run on nearly all the genes, that is useful for heatmaps and data later but can take a while and requires a lot of memory change to TRUE if you have a large dataset. #!!!!variable features set to 2000 for COSMX, set to 3000 for other datasets
  var_feat <- VariableFeatures(srt)
  saveRDS(var_feat, file = paste(path, nameobj, "_SCT_", assay, "_var_feat.RDS", sep=""))
  
  #clean up var_feat
  var_feat <- VariableFeatures(srt)[0:2000] #only use top 2000
  if(is.null(spikein)==FALSE){ #spikein
    var <- VariableFeatures(srt)
    spikein <- var[var %in% spikein] #select spikein that are variable but note that the order of var is maintained.
    var_feat <- c(var_feat, spikein) %>% unique()}
  var_feat <-  var_feat[var_feat %in% rownames(srt)]
  saveRDS(var_feat, file = paste(path, nameobj, "_", assay, "_var_feat.RDS", sep=""))
  assign(paste(nameobj, "_", assay, "_var_feat", sep=""), var_feat, envir = .GlobalEnv)
  
  
  srt <- RunPCA(object = srt, assay="SCT", features=var_feat)
  
  plts <- VariableFeaturePlot(srt)
  plot(plts)
  ggsave2(paste(path, nameobj, "_", nameass, "_VariableFeatures.png", sep=""), plot = plts, width=4, height=2.5)
  
  head(VariableFeatures(srt), 20)
  
  plts <- ElbowPlot(srt) #you need to cut off the number of PCs before the line flattens out aka cut it off at the elbow
  plot(plts)
  ggsave2(paste(path, nameobj, "_", nameass, "_Elbowplot.png", sep=""), plot = plts, width=4, height=2.5)
  
  assign(paste(nameobj), srt, envir = .GlobalEnv)
  assign(paste(nameobj, "_", assay, "_var_feat", sep=""), var_feat, envir = .GlobalEnv)
  saveRDS(var_feat, file = paste(path, nameobj, "_", assay, "_var_feat.RDS", sep=""))
}



NFSP_fun <- function(srt, assay="RNA", spikein=NULL){
  #Performs ScaleData, PCA, and regular plots.
  #srt: seurat object
  
  nameobj <- deparse(substitute(srt)) #make char
  path <- paste(curr_dir, "/Sample_preprocess/", nameobj, "/Norm_Clustering/", sep="") #specify dir
  if(dir.exists(path)==FALSE){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  DefaultAssay(srt) <- assay
  srt <- NormalizeData(srt, assay=assay)
  srt <- FindVariableFeatures(srt, assay=assay, nfeatures=10000)
  
  #clean up var_feat
  var_feat <- VariableFeatures(srt)[0:3000] #only use top 3000
  # var_feat <-  var_feat[-grep("^RP.*|XIST", var_feat)] !!!! don't remove things unless you need to
  if(is.null(spikein)==FALSE){ #spikein
    var <- VariableFeatures(srt)
    spikein <- var[var %in% spikein] #select spikein that are variable but note that the order of var is maintained.
    var_feat <- c(var_feat, spikein) %>% unique()}
  var_feat <-  var_feat[var_feat %in% rownames(srt)]
  saveRDS(var_feat, file = paste(path, nameobj, "_", assay, "_var_feat.RDS", sep=""))
  assign(paste(nameobj, "_", assay, "_var_feat", sep=""), var_feat, envir = .GlobalEnv)
  
  srt <- ScaleData(srt, assay=assay, features=var_feat)
  srt <- RunPCA(srt, assay=assay, features=var_feat)
  plts <- VariableFeaturePlot(srt)
  plot(plts)
  ggsave2(paste(path, nameobj, "_VariableFeatures.png", sep=""), plot = plts, width=4, height=2.5)
  
  head(VariableFeatures(srt), 20)
  
  plts <- ElbowPlot(srt) #you need to cut off the number of PCs before the line flattens out aka cut it off at the elbow
  plot(plts)
  ggsave2(paste(path, nameobj, "_Elbowplot.png", sep=""), plot = plts, width=4, height=2.5)
  
  assign(paste(nameobj), srt, envir = .GlobalEnv)
}


#Cluster and nonlinear dimensionality reduction function

CNDRfun <- function(srt, PCA_dims, assay = "soupXcounts", reduction = "pca"){
  #Performs FindNeighbours, RunUMAP, FindClusters with wide ranges to the use clustree to make a dendogram of celltypes.
  #srt: seurat object
  #assay: Specifies the assay to use (default=soupXcounts)
  #PCA_dims: Dimensions (1:X) based on PCA elbow plot, start with 1:10 then check elbow plot.
  #reduction: Specifies batch correction method in FindNeighbours, RunUMAP (default=pca)
  #Returns modified seurat object and some plots
  
  #Check if there are any NAs in clusters, which causes errors in clustree
  md <- srt@meta.data
  res <- grep("res", names(md), value = TRUE)
  for(i in res){
    if(any(is.na(md[[i]]))==TRUE){
      stop(paste("Error found NAs in", i))}}
  
  nameobj <- deparse(substitute(srt)) #make char
  path <- paste(curr_dir, "/CellType/", nameobj, "/", sep="") #specify dir
  if(dir.exists(path)==FALSE){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  DefaultAssay(srt) <- assay
  srt <- FindNeighbors(object = srt, dims = PCA_dims, assay = assay, reduction = reduction)
  srt <- RunUMAP(object = srt, dims = PCA_dims, assay = assay, reduction =  reduction)
  
  res_range <- c(seq(from = 0.03, to = 0.09, by = 0.02), seq(from = .1, to = .7, by = 0.2)) #!!!!Modify if necessary, it's usually between there though
  srt <- FindClusters(srt, resolution = res_range, assay=assay)
  plts <- clustree(srt)
  plot(plts)
  ggsave2(paste(path, nameobj, "_cluster_tree.png", sep=""), plot = plts, width=14, height=6)
  
  assign(paste(nameobj), srt, envir = .GlobalEnv)
}

CNDRfun_quick <- function(srt, PCA_dims, assay = "soupXcounts", reduction = "pca", resolution=0.1){
  #Performs FindNeighbours, RunUMAP, FindClusters with specific resolution
  #srt: seurat object
  #assay: Specifies the assay to use (default=soupXcounts)
  #PCA_dims: Dimensions (1:X) based on PCA elbow plot, start with 1:10 then check elbow plot.
  #reduction: Specifies batch correction method in FindNeighbours, RunUMAP (default=pca)
  #Returns modified seurat object and some plots
  
  #Check if there are any NAs in clusters, which causes errors in clustree
  md <- srt@meta.data
  res <- grep("res", names(md), value = TRUE)
  for(i in res){
    if(any(is.na(md[[i]]))==TRUE){
      stop(paste("Error found NAs in", i))}}
  
  nameobj <- deparse(substitute(srt)) #make char
  path <- paste(curr_dir, "/CellType/", nameobj, "/", sep="") #specify dir
  if(dir.exists(path)==FALSE){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  DefaultAssay(srt) <- assay
  srt <- FindNeighbors(object = srt, dims = PCA_dims, assay = assay, reduction = reduction)
  srt <- RunUMAP(object = srt, dims = PCA_dims, assay = assay, reduction =  reduction)
  
  srt <- FindClusters(srt, resolution = resolution, assay = assay)
  
  assign(paste(nameobj), srt, envir = .GlobalEnv)
}


#Normalization of peaks

#Fix ElbowPlot so that dim 1 can be excluded
ElbowPlotfun <- function(object, ndims = 1:20, reduction = "pca") 
{
  data.use <- Stdev(object = object, reduction = reduction)
  if (length(x = data.use) == 0) {
    stop(paste("No standard deviation info stored for", 
               reduction))
  }
  stdev <- "Standard Deviation"
  plot <- ggplot(data = data.frame(dims = ndims, stdev = data.use[ndims])) + 
    geom_point(mapping = aes_string(x = "dims", y = "stdev")) + 
    labs(x = gsub(pattern = "_$", replacement = "", x = Key(object = object[[reduction]])), 
         y = stdev) + theme_cowplot()
  return(plot)
}

peaknormfun <- function(srt){
  #Performs standard peak data normalization, PCA and plots (note the first PC is not plotted in the elbowplot)
  #srt: seurat object
  
  nameobj <- deparse(substitute(srt)) #make char
  path <- paste(curr_dir, "/", nameobj, "/Norm_Clustering/", sep="") #specify dir
  if(dir.exists(path)==FALSE){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  srt <- RunTFIDF(srt, assay = "ATAC") #normalize
  srt <- FindTopFeatures(srt, assay = "ATAC", min.cutoff = "q5") #Top features
  srt <- RunSVD(srt, assay = "ATAC") #Dimension reduction with SVD
  var_feat <- VariableFeatures(srt, assay = "ATAC")
  save(var_feat, file = paste(curr_dir, "/", nameobj, "/", nameobj, "_ATAC_var_feat.RDS", sep=""))
  
  plts <- ElbowPlotfun(srt, ndims=2:20, reduction = "lsi") #you need to cut off the number of PCs before the line flattens out aka cut it off at the elbow
  plot(plts)
  ggsave2(paste(path, nameobj, "_Elbowplot.png", sep=""), plot = plts, width=4, height=2.5)
  
  plts <- DepthCor(srt, n=40) #Check if any PCs correlate to the depth of sequencing. Usually the first one will.
  plot(plts)
  ggsave2(paste(path, nameobj, "_Correlation_to_seqdepth.png", sep=""), plot = plts, width=4, height=2.5)
  
  assign(paste(nameobj), srt, envir = .GlobalEnv)
  assign(paste(nameobj, "_ATAC_var_feat", sep=""), var_feat, envir = .GlobalEnv)
}



CNDRdimfun <- function(srt, reduction="pca", res=NULL, assay = "soupXcounts", split.by = NULL, group.by = NULL, cols=NULL){
  #Plots standard Dim plots and a couple Feature and Violin plots to check cell quality
  #srt: seurat object
  #assay: Specifies the assay to use (default=SCT)
  #reduction: Specifies batch correction method to save some plots (default = pca)
  #res: Specifies the resolution paramater to save in different folders (default = NULL)
  #split.by and group.by: Can take lists of metadata variables to split or group the data by when making plots.
  
  nameobj <- deparse(substitute(srt)) #make char
  path <- paste0(curr_dir, "/CellType/", nameobj, "/", reduction, res, "/") #specify dir
  if(dir.exists(path)==FALSE){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  plts <- DimPlot(srt, raster=FALSE, cols=cols) #Make dimplots
  plot(plts)
  ggsave2(paste(path, nameobj, "_umap.png", sep=""), plot = plts, height = 4, width = 7)
  
  if(is.null(split.by)==FALSE){
    for(var in split.by){
      plts <- DimPlot(srt, split.by = var, raster=FALSE, cols=cols) #Make dimplots and split them based on group
      plot(plts)
      width_max <- 4*length(unique(srt@meta.data[[var]]))
      if(width_max > 166){width_max <- 166} #max width to save
      ggsave2(paste(path, nameobj, "_umap_split_", var, ".png", sep=""), plot = plts, limitsize = FALSE, height = 4, width = width_max)
    }}
  
  if(is.null(group.by)==FALSE){
    for(var in group.by){
      plts <- DimPlot(srt, group.by = var, raster=FALSE) #Make dimplots and group them based on group
      plot(plts)
      ggsave2(paste(path, nameobj, "_umap_", var, ".png", sep=""), plot = plts, width = 5+floor(length(unique(srt@meta.data[[var]]))/17)*2, height=4)
    }}
  
  #Check if low quality cells are clustering together
  if(assay=="peaks" | assay=="ATAC"){
    plts <- FeaturePlot(srt, "TSS.enrichment", raster=FALSE, max.cutoff=6)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_TSS.enrichment", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "nCount_ATAC", raster=FALSE, max.cutoff=1000)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_nCount_ATAC", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "blacklist_ratio", raster=FALSE)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_blacklist_ratio", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "nucleosome_signal", raster=FALSE, max.cutoff=1.5)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_nucleosome_signal", ".pdf", sep=""), plot = plts, width=4, height=4)
    plts <- VlnPlot(srt, "TSS.enrichment", raster=FALSE, cols=cols)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_TSS.enrichment", ".pdf", sep=""), plot = plts, width=4, height=4)
    plts <- VlnPlot(srt, "nCount_ATAC", raster=FALSE, cols=cols)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_nCount_ATAC", ".pdf", sep=""), plot = plts, width=4, height=4)
    plts <- VlnPlot(srt, "blacklist_ratio", raster=FALSE, cols=cols)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_blacklist_ratio", ".pdf", sep=""), plot = plts, width=4, height=4)
    plts <- VlnPlot(srt, "nucleosome_signal", raster=FALSE, cols=cols)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_nucleosome_signal", ".pdf", sep=""), plot = plts, width=4, height=4)
  }else{
    plts <- FeaturePlot(srt, "nFeature_RNA", raster=FALSE)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_nFeatureRNA", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "nFeature_RNA", raster=FALSE, max.cutoff=1000)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_nFeatureRNA_max1000", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "nFeature_RNA", raster=FALSE, min.cutoff=8000)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_nFeatureRNA_min9000", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "percent.rb", raster=F, max.cutoff = 5)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_percentrb_5", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "percent.mt", raster=FALSE)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_percentmt", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- VlnPlot(srt, "nFeature_RNA", pt.size=F, cols=cols)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_Vln_nFeatureRNA", ".pdf", sep=""), plot = plts, width=4, height=4)
    plts <- VlnPlot(srt, "percent.mt", pt.size=F, cols=cols)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_Vln_percentmt", ".pdf", sep=""), plot = plts, width=4, height=4)
    plts <- VlnPlot(srt, "nCount_RNA", pt.size=F, cols=cols, y.max=10000)
    plot(plts)
    ggsave2(paste(path, nameobj, "_Vln_nCount_RNA", ".pdf", sep=""), plot = plts, width=4, height=4)
    #plts <- FeaturePlot(srt, "scrublet_score", raster=FALSE, order=TRUE)
    #plot(plts)
    #ggsave2(paste(path, nameobj, "_cluster_FP_scrublet", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "DoubletFinder_score", raster=FALSE, order=TRUE)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_doubletFinder", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "scDblFinder_doublet_score", raster=FALSE, order=TRUE)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_scDbltFinder", ".png", sep=""), plot = plts, width=4, height=4)
    plts <- FeaturePlot(srt, "scds_hybrid_score", raster=FALSE, order=TRUE)
    plot(plts)
    ggsave2(paste(path, nameobj, "_cluster_FP_scds", ".png", sep=""), plot = plts, width=4, height=4)
  }
}


SSCNDRdimfun <- function(srt, reduction="pca", res=NULL, assay = "soupXcounts", split.by = NULL, group.by = NULL, cols=NULL){
  #Plots standard Dim plots and a couple Feature and Violin plots to check cell quality
  #srt: seurat object
  #assay: Specifies the assay to use (default=SCT)
  #reduction: Specifies batch correction method to save some plots (default = pca)
  #res: Specifies the resolution paramater to save in different folders (default = NULL)
  #split.by and group.by: Can take lists of metadata variables to split or group the data by when making plots.
  
  nameobj <- deparse(substitute(srt)) #make char
  path <- paste(curr_dir, "/CellType/", nameobj, "/", reduction, res, "/", sep="") #specify dir
  if(dir.exists(path)==F){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  plts <- DimPlot(srt, raster=T, cols=cols) #Make dimplots
  plot(plts)
  ggsave2(paste(path, nameobj, "_umap.pdf", sep=""), plot = plts, height = 4, width = 7)
  
  if(is.null(split.by)==F){
    for(var in split.by){
      plts <- DimPlot(srt, split.by = var, raster=T, cols=cols) #Make dimplots and split them based on group
      plot(plts)
      width_max <- 4*length(unique(srt@meta.data[[var]]))
      if(width_max > 166){width_max <- 166} #max width to save
      ggsave2(paste(path, nameobj, "_umap_split_", var, ".pdf", sep=""), plot = plts, limitsize = T, height = 4, width = width_max)   
    }}
  
  if(is.null(group.by)==F){
    for(var in group.by){
      plts <- DimPlot(srt, group.by = var, raster=T) #Make dimplots and group them based on group
      plot(plts)
      ggsave2(paste(path, nameobj, "_umap_", var, ".pdf", sep=""), plot = plts, width = 7, height=4)
    }}
  
  #Check if low quality cells are clustering together
  plts <- FeaturePlot(srt, "Area.um2", raster=T)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Area.um2", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- FeaturePlot(srt, "complexity", raster=T, max.cutoff=1000)
  plot(plts)
  ggsave2(paste(path, nameobj, "_complexity", ".pdf", sep=""), plot = plts, width=4, height=4)
  
  plts <- FeaturePlot(srt, "nFeature_Nanostring", raster=T)
  plot(plts)
  ggsave2(paste(path, nameobj, "_cluster_FP_nFeature_Nanostring", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- FeaturePlot(srt, "nFeature_Nanostring", raster=T, max.cutoff=1000)
  plot(plts)
  ggsave2(paste(path, nameobj, "_cluster_FP_nFeature_Nanostring_max1000", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- FeaturePlot(srt, "nFeature_Nanostring", raster=T, min.cutoff=1500)
  plot(plts)
  ggsave2(paste(path, nameobj, "_cluster_FP_nFeature_Nanostring_min1500", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- VlnPlot(srt, "nFeature_Nanostring", pt.size=F, cols=cols)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Vln_nFeature_Nanostring", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- VlnPlot(srt, "nCount_Nanostring", pt.size=F, cols=cols)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Vln_nCount_Nanostring", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- VlnPlot(srt, "Area.um2", pt.size=F, cols=cols)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Vln_Area.um2", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- VlnPlot(srt, "complexity", pt.size=F, cols=cols)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Vln_complexity", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- VlnPlot(srt, "Mean.PanCK", pt.size=F, cols=cols)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Vln_Mean.PanCK", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- VlnPlot(srt, "Mean.CD68_CK8_18", pt.size=F)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Vln_Mean.CD68_CK8_18", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- VlnPlot(srt, "Mean.CD45", pt.size=F, cols=cols)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Vln_Mean.CD45", ".pdf", sep=""), plot = plts, width=4, height=4)
  plts <- VlnPlot(srt, "Mean.DAPI", pt.size=F, cols=cols)
  plot(plts)
  ggsave2(paste(path, nameobj, "_Vln_Mean.DAPI", ".pdf", sep=""), plot = plts, width=4, height=4)
}



#Function to label missing cells from a dataset
label_doubletsfun <- function(srt_no_doublet, srt_all){
  #labels cells that were removed in the cleaner dataset as doublets
  nameobj <- deparse(substitute(srt_all)) #make char
  
  srt_no_doublet@meta.data$Doublet <- "Singlet"
  srt_no_doublet_meta <- srt_no_doublet@meta.data
  srt_all_meta <- srt_all@meta.data
  
  srt_all_meta <- rownames_to_column(srt_all_meta, "rowname")
  srt_no_doublet_meta <- rownames_to_column(srt_no_doublet_meta, "rowname")
  srt_no_doublet_meta <- data.frame(srt_no_doublet_meta$rowname, srt_no_doublet_meta$Doublet)
  colnames(srt_no_doublet_meta) <- c("rowname", "Doublet")
  srt_all_meta <- left_join(srt_all_meta, srt_no_doublet_meta, by="rowname")
  srt_all_meta <- column_to_rownames(srt_all_meta, "rowname")
  srt_all_meta$Doublet <- as.character(srt_all_meta$Doublet)
  srt_all_meta$Doublet[is.na(srt_all_meta$Doublet)==TRUE] <- "Doublet"
  srt_all@meta.data <- srt_all_meta
  
  assign(paste(nameobj), srt_all, envir = .GlobalEnv)
}



#Find markers function
get_markersfun <- function(seurat_object, cluster, assay="soupXcounts", features=NULL){
  #Function to extract markers per different conditions (basic FindMarkers and more advanced FindConservedMarkers)
  #seurat_object
  #cluster number to analyse
  #assay: default is soupXcounts
  #features: Specify features to test
  #Returns a list of markers per condition
  c <- FindMarkers(seurat_object, ident.1 = cluster, only.pos = TRUE, logfc.threshold = 0.25, assay = assay, features = features, recorrect_umi = FALSE) #Specify how to use FindMarkers or other
  c <- c %>% rownames_to_column("gene")
  c <- c[order(c$avg_log2FC, decreasing = TRUE),] #Orders markers in descending order by logFC
  
  #!!!!Modify all FindConservedMarkers code chunks
  #cgt <- FindConservedMarkers(seurat_object, ident.1 = cluster, grouping.var = "group_time", only.pos = TRUE, logfc.threshold = 0.25, assay = assay)
  #cgt <- cgt %>% mutate(Total_avg_fc = rowMeans(cgt[,str_detect(names(cgt), "log2FC")]) %>% rownames_to_column("gene") #Makes a single avg_logFC column for next step
  #cgt <- cgt[order(ctg$avg_fc, decreasing = TRUE),]
  
  #cg <- FindConservedMarkers(seurat_object, ident.1 = cluster, grouping.var = "group", only.pos = TRUE, logfc.threshold = 0.25, assay = assay)
  #cg <- cg %>% mutate(Total_avg_fc = rowMeans(cg[,str_detect(names(cg), "log2FC")])) %>% rownames_to_column("gene") #Makes a single avg_logFC column for next step
  #cg <- cg[order(cg$avg_fc, decreasing = TRUE),]
  
  #ct <- FindConservedMarkers(seurat_object, ident.1 = cluster, grouping.var = "time", only.pos = TRUE, logfc.threshold = 0.25, assay = assay)
  #ct <- ct %>% mutate(Total_avg_fc = rowMeans(ct[,str_detect(names(ct), "log2FC")])) %>% rownames_to_column("gene") #Makes a single avg_logFC column for next step
  #ct <- ct[order(ct$avg_fc, decreasing = TRUE),]
  
  list(c) %>% 
    list(cluster_id = cluster, .)
}


#Aggregate clusters in a list

return_markersfun <- function(seurat_object, assay="soupXcounts", features=NULL, num_clust=NULL){
  #Runs get_markersfun in a forloop for every cluster in the sample
  #Needs a get_markersfun tailored to your experimental design in your global environment
  #seurat_object
  #features: Specify features to test, NULL by default
  #num_clust: Specify number of clusters 0:x, NULL by default will determine the total range of clusters.
  #Returns and saves a list of genes more highly expressed in each cluster as lists.
  nameobj <- deparse(substitute(seurat_object)) #make char
  DefaultAssay(seurat_object) <- assay
  path <- paste(curr_dir, "/CellType/", nameobj, "/", sep="")
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  if(is.null(num_clust)){num_clust <- c(0:(nlevels(seurat_object@active.ident)-1))}
  markers_list <- empty.dump() #create empty list for the for loop to use
  for(i in num_clust){
    markers_i <- get_markersfun(seurat_object = seurat_object, cluster = i, assay=assay, features = features)
    if(!is.null(markers_i)){
      markers_list[[paste(i)]] <- markers_i
    }else{
      markers_list[[paste(i)]] <- "Too few cells in cluster to analyse"}
  }
  assign(c(paste(nameobj, "_markers_list", sep = "")), markers_list, envir = .GlobalEnv) #Save in Global environment
  saveRDS(markers_list, paste(path, nameobj, "_markers_list.RDS", sep = ""))
}


#Feature plots for top 20 genes

Top20_Featuresfun <- function(seurat_object, Cluster_number){
  #Loops through the top 20 markers in the first list (from return_markersfun) and outputs and saves their Featureplot and outputs the result for the cluster
  #seurat_object
  #Cluster_number: To analyse
  tryCatch({
    Cluster_numb <- paste(Cluster_number)
    nameobj <- deparse(substitute(seurat_object)) #make char
    markers_list <- get(paste(nameobj, "_markers_list", sep=""), envir = .GlobalEnv)
    for(i in 1:20){ #I am limiting it to the top 20, some do 10 #I also limited it to the first dataframe of markers found since there wasn't a huge difference between them.
      gene <- markers_list[[Cluster_numb]][[2]][[1]][["gene"]][i]
      plts <- FeaturePlot(object = seurat_object, features = c(gene), order = TRUE, min.cutoff = 'q10', label = TRUE, repel = TRUE, raster=FALSE) #make a feature plot for every gene
      ggsave2(paste(curr_dir, "/CellType/", nameobj, "/Cluster/", Cluster_number, "/", gene, ".png", sep=""), plot = plts, width=4, height=2.5)
      print(plts)
    }
    print(head(markers_list[[Cluster_numb]][[2]][[1]], n=20)) #showing the top 20
    #print(head(markers_list[[Cluster_numb]][[2]][[2]], n=20)) #showing the top 20 !!!!Add however many as you have markers lists within each cluster. The name is [[Cluster_number]][[leave_as_a_2]][[number_of_markers_dataframe_based_on_your_get_markersfun]]
  }, error=function(e){cat("Warning : You may have fewer than 20 Top Features.",conditionMessage(e), "\n")})
}



#Feature plot of known markers

search_markersfun <- function(seurat_object, Cluster_number, known_markers, MakeFigs=F){
  #Loops through the known_markers list and outputs and saves their Featureplot and outputs the result for the cluster (if they are in the markers list)
  #seurat_object
  #known_markers: List of known markers as characters (repeats are fine, extra commas are not)
  #Cluster_number: To analyse
  #MakeFigs: Make FeaturePlots, default is FALSE to get an idea of what's there. Turn to TRUE to make sure the distribution makes sense.
  
  nameobj <- deparse(substitute(seurat_object)) #make char
  known_markers <- unique(known_markers)
  markers_list <- get(paste(nameobj, "_markers_list", sep=""), envir = .GlobalEnv)
  
  #Make lists
  sig <- c(known_markers[known_markers %in% markers_list[[as.character(Cluster_number)]][[2]][[1]][["gene"]]])
  print("Siginificant")
  print(sig)
  print("Not Significant")
  print(known_markers[!(known_markers %in% markers_list[[as.character(Cluster_number)]][[2]][[1]][["gene"]])])
  
  #Make feature plots
  if(MakeFigs==TRUE){
    for(i in sig){
      tryCatch({
        plts <- FeaturePlot(object = seurat_object, features = i, order = TRUE, min.cutoff = 'q10', label = TRUE, repel = TRUE, raster=FALSE) #make a feature plot for every gene
        ggsave2(paste(curr_dir, "/CellType/", nameobj, "/Cluster/", Cluster_number, "/", i, ".png", sep=""), plot = plts, width=4, height=2.5)
        plot(plts)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
  
  assign(c(paste(nameobj, Cluster_number, "sig_list", sep = "_")), sig, envir = .GlobalEnv) #Save in Global environment
}


search_insig_markersfun <- function(srt, Cluster_number, known_markers, assay="soupXcounts"){
  #Loops through the known_markers list and outputs and saves their Featureplot and outputs the result for the cluster (if they are in the markers list)
  #srt
  #known_markers: List of known markers as characters (repeats are fine, extra commas are not)
  #Cluster_number: To analyse
  #assay: char, "soupXcounts" default
  nameobj <- deparse(substitute(srt)) #make char
  known_markers <- unique(known_markers)
  markers_list <- get(paste(nameobj, "_markers_list", sep=""), envir = .GlobalEnv)
  sub_seurat <- subset(srt, idents = as.character(Cluster_number)) #remove anything that is lowly expressed in cluster
  df <- data.frame(sub_seurat@assays[[assay]]@data)
  known_markers <- known_markers[known_markers %in% srt@assays[[assay]]@data@Dimnames[[1]]] #makes sure genes are in assay.
  
  for(i in known_markers){
    tryCatch({
      if(max(df[i,]) < 1){
        print(paste("Warning:", i, " is not different accross clusters"))
      }else{
        plts <- FeaturePlot(object = sub_seurat, features = i, order = TRUE, min.cutoff = 'q10', label = TRUE, repel = TRUE, raster=FALSE) #make a feature plot for every gene
        print(plts)
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

#Function to check if two genes belong in the same cluster

OneThesefun <- function(seurat_object, Cluster_number, gene1, gene2, assay="soupXcounts"){
  #One of these things (is not like the other), checks discordant genes and returns a venn diagram (set_1 = gene1; set_2 = gene2, set_3 = all cells in cluster). This allows you to see if there is overlap between the two genes or if there might be a few different cell types within this cluster.
  #seurat_object
  #Cluster_number to subset and check
  #gene1 typically found in one cell type (char)
  #gene2 typically found in another cell type (char)
  sub_seurat <- subset(seurat_object, idents = as.character(Cluster_number))
  total <- GetAssayData(sub_seurat, assay = assay, slot = "data")
  g1 <- total[gene1,]
  g1 <- names(which(g1>max(g1)/2))
  g2 <- GetAssayData(sub_seurat, assay = assay, slot = "data")[gene2,]
  g2 <- names(which(g2>max(g2)/2))
  total <- total@Dimnames[[2]]
  ggVennDiagram(list(g1,g2, total))
}


#Function to make individual FeaturePlots and Vlnplots of genes of interest, useful for a quick search

Featurefun <- function(srt, gene_set, gene_set_name, max.cutoff= c(NA), assay="soupXcounts"){
  #Function to make feature plots and default violin plots of genes in gene_set at different resolutions (max.cutoff)
  #seurat_object
  #gene_set: list of genes (char) to plot
  #gene_set_name: char, Name of gene set folder to save under
  #max.cutoff: c(char) list of cutoffs, default is NA
  #assay: char, default is "soupXcounts"
  #Returns and saves plots in curr_dir/name of srt/Featureplots/gene_set_name
  
  #Set up
  nameobj <- deparse(substitute(srt))
  path = paste(curr_dir, "/CellType/", nameobj, "/Featureplots/", gene_set_name, "/", sep="")
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  gene_set <- gene_set %>% unique()
  gene_set <- gene_set[gene_set %in% srt@assays[[assay]]@data@Dimnames[[1]]] #makes sure genes are in assay.
  
  ##FeaturePlot  
  for(maxcutoff in max.cutoff){
    for(gene in gene_set){
      plts <- FeaturePlot(srt, features = gene, pt.size=.02, max.cutoff=maxcutoff, label = FALSE, raster=FALSE)
      if(is.na(maxcutoff)==TRUE){maxcutoff<-"Default"}
      ggsave2(paste(path, maxcutoff, "/", gene, ".png", sep=""), plot = plts, width=7, height=8)
      print(plts)
    }}
  
  
  ##Violin
  for(gene in gene_set){
    data <- data.frame(gene = srt[[assay]]@data[gene,], cluster = srt@active.ident)
    #data <- data[order(data[,1]),]
    #data2 <- dplyr::filter(data, gene>0)#!!!!This removes cells with 0 counts and may be innapropriate
    #data3 <- data[1:5000,]
    #data <- rbind(data2, data3)
    plts <- ggplot(data, aes(x = cluster, y = gene)) + geom_violin(aes(fill = cluster), trim=TRUE, scale = "width") + theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                            panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust=0.5)) + ggtitle(gene) + ylab("") + xlab("")
    ggsave2(paste(path, "VLN/", gene, ".png", sep=""), plot = plts, width=4, height=6)
    print(plts)
  }
  
}

Featurefun_pdf <- function(srt, gene_set, max.cutoff= c(NA), assay="soupXcounts"){
  #Function to make feature plots and default violin plots of genes in gene_set at different resolutions (max.cutoff)
  #seurat_object
  #gene_set: list of genes (char) to plot
  #max.cutoff: c(char) list of cutoffs, default is NA
  #assay: char, default is "soupXcounts"
  #Returns and saves plots in curr_dir/name of srt/Cluster/Markers of interest/gene_set_name
  
  #Set up
  nameobj <- deparse(substitute(srt))
  path = paste0(curr_dir, "/CellType/Publish/")
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  gene_set <- gene_set %>% unique()
  gene_set <- gene_set[gene_set %in% srt@assays[[assay]]@data@Dimnames[[1]]] #makes sure genes are in assay.
  
  ##FeaturePlot  
  for(maxcutoff in max.cutoff){
    for(gene in gene_set){
      plts <- FeaturePlot(srt, features = gene, pt.size=.02, max.cutoff=maxcutoff, label = FALSE, raster=FALSE) + theme(text=element_text(size=7))
      graph2pdf(plts, paste0(path, gene, maxcutoff, ".pdf"), width=4, height=4)
      plot(plts)
    }}
  
  
  ##Violin
  for(gene in gene_set){
    data <- data.frame(gene = srt[[assay]]@data[gene,], cluster = srt@active.ident)
    #data <- data[order(data[,1]),]
    #data2 <- dplyr::filter(data, gene>0)#!!!!This removes cells with 0 counts and may be innapropriate
    #data3 <- data[1:5000,]
    #data <- rbind(data2, data3)
    plts <- ggplot(data, aes(x = cluster, y = gene)) + geom_violin(aes(fill = cluster), trim=TRUE, scale = "width") + theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust=0.5), text=element_text(size=7)) + ggtitle(gene) + ylab("") + xlab("")
    graph2pdf(plts, paste(path, "VLN/", gene, ".pdf", sep=""), width=4, height=6)
    plot(plts)
  }
  
}

Ncellfun <- function(seq_sample, file_name, group.by=NULL){
  #Compares number and percent of cells per active.ident; returns result in a data.table
  #seq_sample: seurat_object (seurat object)
  #file_name: Save the result as CSV (char)
  #group.by: subdivide by a group (char)
  md <- data.frame(seq_sample@active.ident, seq_sample@meta.data[[group.by]])
  colnames(md) <- c("active.ident", "group.by")
  Ncell <- tabyl(md, active.ident, group.by) #tabulate (but better) number of cells per X and Y
  Ncellrow <- Ncell %>% adorn_percentages("row") %>% adorn_pct_formatting(digits = 2) #row percentages
  
  Ncellcol <- Ncell %>% adorn_percentages("col") %>% adorn_pct_formatting(digits = 2) #col percentages
  
  suppressWarnings(Ncell <- rbind(Ncell, "% by rows", Ncellrow, "% by cols", Ncellcol))
  write.csv(Ncell, file_name) #saving
  print(Ncell)
}
#Ncellfun(Duras2, "clusters per region.csv", "region")


#Known markers lists
known_Mesen <- c("SLC4A4", "LAMA2", "DLK1", "LUM", "DCN", "SNAI", "CDH1", "AKAP1", "PEAR1", "CXCL12", "DCN", "ADAM12", "EMILIN1", "LUM", "COL3A1", "COL1A2", "COL1A1", "ROR1", "EFNA5", "IGF2", "IGF1R", "PDGFRA", "COL1A1", "COL1A2", "LUM", "DCN", "ACTA2", "RGS5", "LUM", "DCN", "DLK1", "COL1A1", "COL3A1", "COL1A1", "MGP", "COL1A2", "LGALS1", "IGF2", "VCAN", "POSTN", "DCN", "FN1", "COL6A2", "DLK1", "COL6A3", "COL5A2", "COL5A1", "CDH11", "CEMIP")
known_Fibro <- c("PDGRFB", "PDGRFA", "COL1A1", "ERTR7", "LAMA1", "CD13", "FN", "ECM1", "LGALS3", "COL1A2", "TGFB2", "LUM", "COL11A1", "NESTIN", "PDGRFB", "PDGRFA", "COL1A1", "ERTR7", "LAMA1", "CD13", "FN", "CRABP2", "ALDH1A2", "SLC6A13", "S100A6", "NGFR", "COL1A1", "COL1A2", "LUM", "DCN", "ACTA2", "RGS5", "PDGRFB", "PDGRFA", "COL1A1", "ERTR7", "LAMA1", "CD13", "FN", "FXYD5", "FOXP1", "SIX1", "PDGRFB", "PDGRFA", "COL1A1", "ERTR7", "LAMA1", "CD13", "FN",  "CEMIP")
known_barrier <- c("CLDN1", "CLDN3", "CLDN4", "CLDN5", "AQP1", "ATP1A1", "CLDN11", "CLDN5", "CLDN4", "CLDN5", "JAMA", "PGP", "LRP1", "RAGE", "GLUT1", "PECAM", "ABCG2", "CDH5", "CGN", "SLC38A5", "ABCC2", "VWF", "SLC7A5")
known_periv <- c("LAMA2", "ZFP36L1", "SVIL", "ADD3", "EEF2", "STAT6", "ENOSF1", "OGA", "BABAM2", "VPS41", "TTC41", "ADK", "CAB39L", "SLC2A3", "FTL", "FLRT2", "MYEF2", "EGFR", "LGR4", "TTC14", "RANBP2", "CAB39L", "TNXB", "KSR1", "MID1", "SLK", "SPAG9", "CEMIP", "BICC1", "FLRT2", "LAMA2", "FBLN1", "APOD", "NTRK3", "PRICKLE1", "FGF13", "CYP1B1", "DCN", "COL12A1", "FHL2", "ANTXR2", "LAMC1", "PID1", "KCNK2", "LAM1", "COL6A3", "LAMA4", "KCNIP1", "FOXP2", "SRPX2", "NR2F1", "SVIL", "GSN", "PTPRS", "HIVEP3", "CDYL", "SOBP")
known_menin <- c("SLC4A4", "KCNMA1", "SLC4A4", "TRPM3") #Have SLC4A4; NO CLDN11, TUJ1, LAMA2 
known_arach <- c("CLDN11", "TJP1", "PTDGS")
known_fibro_prog <- c("SLC6A13", "SLIT2", "FN1", "KCNB2", "CRYM", "NRXN1", "SLC1A3", "SRGAP1", "ROBO2", "PTGDS", "ROBO1", "ROBO2")

known_Mural <- c("TAGLN", "ACTA2", "COL4A1", "MYH11", "SGIP1", "NESTIN", "ACTA2", "CD248", "DLK1", "PDGFRB", "DES", "APOE", "ACTA2", "MYH11", "CNN1", "RGS5", "PDGFRB", "NOTCH3", "MCAM", "CSPG4","TAGLN", "ACTA2", "MYH11", "SLC20A2", "SLC30A10", "GRM8")
known_peric <- c("LAMA2", "GPC5", "PDE7B", "PDGFRB", "DLC1", "LGALS3", "COL1A2", "TGFB2", "LUM", "COL11A1", "PDGFRB", "PDGFRB", "NG2", "DES", "KCNJ8", "ABCC9", "CD13", "GRM8", "PDG-FRB")
known_VSM <- c("MYH11", "TAGLN", "ACTA2", "LPP", "MYL9", "PDGFRB", "NG2", "DES", "CD13", "CD146", "ACTA2", "GRM8", "FRMD3")

known_Endothelial <- c("TTR", "CLDN5", "CLDN4", "ANO2", "ST6GALNAC3", "ELOVL7", "FLT1", "ABCB1", "PLVAP", "PECAM1", "CDH5", "KDR", "SELE", "VWF", "ARL15", "MECOM", "VWF", "HSP90AA1", "ANO2", "LDB2", "FLT1", "SERPINE1", "DOCK9", "HSPH1", "HSPD1", "PECAM1", "MSN", "TCF4", "PELI1", "PTGES3", "PIK3R3", "P4HA1", "RBMS3", "PLS3", "CRIM1", "IL3RA", "HSPA5", "RALGAPA2", "CLIC4", "MOB4", "DNAJB6", "PTPRB", "PICALM", "PLPP1", "SERPINH1", "ATP8B1", "MYH9", "TFPI", "BMPR2", "ST8SIA6", "UBC", "TEK", "CMIP", "GPBP1", "ST6GALNAC3", "HSPB1", "GRB10", "MCTP1", "MLKL", "ARHGAP26", "HSP90AB1", "CDK17", "DKK2", "BMP6", "PLEKHG1", "CDH13", "UTRN", "ITGA5", "PODXL", "ELOVL7", "MRTFB", "FNIP2", "MBNL1", "DOCK4", "HSPG2", "SPAG9", "PTGS2", "ERG", "SETBP1", "GRK5", "EGFL7", "IRAK2", "NOTCH4", "FNDC3B", "SASH1", "ARHGAP29", "THSD7A", "AKT3", "VEGFC", "SEC14L1", "HSPE1", "ADAMTS9", "ADGRL4", "BAG3", "SLC2A3", "ELMO1", "PALM2-AKAP2", "IL4R", "TPM4", "CLIC2", "CTTNBP2NL", "COL8A1", "TANC1", "NEDD9", "ATP10A", "INSR", "ERO1A", "APBB2", "NFIB", "CCNY", "RAPGEF2", "SHANK3", "MAST4", "RUNX1T1", "PRKCH", "STOM", "A2M", "PDE4D", "ARAP3", "PIK3C2B", "LIMCH1", "ABCB1", "FLNB", "DNAJB1", "ECE1", "RAPGEF5", "UGCG", "ENG", "PRKY", "GFOD1", "NKAIN2", "KMT2E", "PKHD1L1", "ARID4B", "SWAP70", "RND1", "QKI", "ASAP1", "NUDT4", "IGFBP3", "CAPZA2", "ADAMTS1", "LRMDA", "STIP1", "KIAA0355", "FCHSD2", "HNRNPU", "FLI1", "ZEB1", "STC1", "RPS6KA3", "PARVB", "HSPA8", "DNAJA1", "SORBS1", "ETV6", "DOCK6", "ARID5A", "EFNB2", "MALAT1", "GRASP", "PCSK5", "IQSEC1", "HRH1", "FRYL", "CALCRL", "AC010737.1", "FGD5", "CALD1", "ADIPOR2", "TCF12", "CHORDC1", "PTPRM", "JUN", "PDLIM5", "PTK2", "ITSN2", "TMTC1", "SAT1", "EMCN", "DUSP5", "KAT6A", "PREX2", "ETS2", "SNAP23", "TRA2B", "WSB1", "CLK1", "TSHZ2", "CAMK2D", "PRRC2C", "ADD1", "COL4A1", "RAPGEF4", "ELOVL5", "ICAM1", "SVIL", "PTPN12", "LRIG1", "PRSS23", "KYNU", "ARHGAP31", "CPNE8", "ELF2", "AC008050.1", "JMJD1C", "VMP1", "MEF2A", "MACF1", "PLCB1", "ABR", "HSP90B1", "HIF3A", "MTCH1", "CBLB", "TSC22D2", "PDLIM1", "PTGIS", "HMCN1", "GNB4", "CDC42BPB", "GALNT18", "PLCG2", "ATP11A", "HIF1A", "ELL2", "FBLIM1", "B2M", "ZFAND3", "DIPK2B", "PIK3C2A", "FYN", "ATP13A3", "AHCTF1", "FBXO11", "KIF13B", "ATF3", "SRGN", "MED13L", "AAK1", "RAP1B", "ABI3BP", "GNAQ", "EXT1", "NECTIN2", "IRAK3", "ARHGEF7", "SMAD6", "GPCPD1", "CFH", "CTNND1", "GNB1", "XAF1", "TRIO", "PTMA", "DGKH", "DENND5A", "RBMS2", "AC007681.1", "MRTFA", "TCF7L1", "STK38L", "RNF144B", "CFLAR", "PLA2G4A", "ARHGEF2", "COL25A1", "BAZ2B", "TDRD10", "SLC20A1", "ZNF423", "APP", "SESTD1", "NUAK1", "DOCK1", "NDRG1", "WAC", "NUCB1", "LIMS2", "OPHN1", "BACE2", "HSPA1A", "TLL1", "ZFAND2A", "MAP4", "WWTR1", "NAMPT", "PCDH11X", "NOVA2", "LINC01578", "HERC2", "HLA-E", "HIP1", "TOP1", "HNRNPA2B1", "GOLPH3", "HLA-B", "ADGRF5", "UBE2D3", "NEDD4", "AC007952.4", "PPP1CB", "TGM2", "PLSCR4", "PCDH11Y", "MEF2C", "AC020916.1", "ESYT2", "CNOT6L", "SLCO4A1", "TSPAN14", "FAM155A", "HSPA1B", "RIMKLB", "ANKRD28", "SAMD12", "STAT3", "TNFRSF10B", "ACTN1", "RBPMS", "COL4A2", "TMEM140", "MYO1E", "RNF145", "RBFOX2", "TJP2", "ZNF83", "CHRM3", "SH3RF3", "CLEC1A", "LAMA3", "C22orf34", "ATP2C1", "UBR2", "KMT2C", "ITPR2", "PPP2R2A", "NOP58", "NCOA3", "FOSL2", "CYTH1", "ICA1", "CCNL1", "CCN2", "RICTOR", "DNAJA4", "EXOC6", "GABARAPL1", "MGP", "TANC2", "ARHGAP23", "ADGRL2", "LINC02147", "CNTNAP3B", "NPDC1", "ARHGAP17", "NR4A1", "NEAT1", "ARHGEF12", "SHROOM4", "EWSR1", "CTNNB1", "GIPC2", "GASK1B", "FXR1", "HSPA4", "ABL2", "NFKB1", "EML4", "USPL1", "RHOJ", "ITGB1", "ACTR3", "PPP3CA", "ATP1A1", "STK3", "GBP2", "ITGAV", "PSMA6", "HMBOX1", "AZIN1", "HECW2", "GMDS", "AC119674.1", "DNAJC1", "SOX5", "EMP1", "FERMT2", "STAG1", "EPAS1", "ZSWIM6", "TNFRSF10D", "IL1R1", "LIFR", "SQSTM1", "MSRB3", "PEAK1", "NLRP1", "KPNA1", "TBCD", "ARPC2", "FMNL3", "VCAM1", "TIMP3", "MAP4K4", "CASC3", "SLC45A4", "ARL5B", "FOS", "IFI16", "STX12", "AHR", "TBK1", "HERPUD1", "YES1", "CHIC2", "AFDN", "TOX", "PDE3A", "CADPS2", "CNOT2", "PRKD1", "FAT4", "HBP1", "COL5A2", "LDLRAD3", "ARGLU1", "ADAM17", "RPGR", "NRIP1", "SESN3", "ST13", "PPP1R12A", "B4GALT1", "CWF19L2", "RIT1", "ZFPM2", "SULF2", "ACTB", "PPP1R15A", "FKBP1A", "IL15", "MTHFD1L", "PIP5K1C", "IQCJ-SCHIP1", "SUSD6", "SMURF2", "PDE10A", "KCNT2", "IPO7", "ERBIN", "UBE2B", "PGM5", "SNRK", "SOS1", "NLK", "SRSF7", "PPFIBP1", "CX3CL1", "TAF1", "NCKAP1", "NFKBIA", "NPIPB5", "HSPA4L", "ITGA1", "IFNAR2", "CACYBP", "CCN1", "TLE1", "KLF6", "TMC7", "HIVEP2", "VPS13C", "NAV1", "CMTM8", "YTHDC1", "TBC1D4", "CALR", "ATL2", "TGFBR3", "CEMIP2", "CARMIL1", "CASP4", "MOSMO", "RDX", "BCAS2", "UBA6", "CHD1", "RASA2", "PITPNB", "DENND3", "FRMD4B", "LDLR", "ABLIM3", "PAPSS2", "ID2", "CCDC88A", "PPP1R3B", "GAREM1", "FAM241A", "DYSF", "WARS", "KLHL2", "RHEB", "KDM7A", "YPEL2", "PER1", "MCC", "ACSL3", "ELN", "FUT8", "SMAD1", "DISC1", "TM4SF1", "PRKD2", "LDHA", "PIP4K2A", "TPM3", "MIR222HG", "INSIG1", "RPH3AL", "SFT2D2", "AHSA1", "MGAT5", "BHLHE40", "TSHZ1", "TSPYL2", "PHF14", "ZC3HAV1", "GALNT15", "LIMA1", "AC004889.1", "RELB", "KLF7", "ARHGEF28", "ZNF654", "DDX3X", "MAN1A1", "ADAM9", "RGS3", "PTBP3", "ROCK2", "FN1", "GCNT1", "TMEM87B", "AKAP12", "PREX1", "PLPP3", "JUNB", "CPAMD8", "MET", "CDKN1A", "TMTC2", "B4GALT5", "BICD1", "CHN1", "TXNIP", "BACH2", "SLCO2A1", "DDX58", "EIF4A1", "ALPL", "CORO1C", "PBX3", "CCDC85A", "AUTS2", "PARP14", "ITGA9", "PLA2G4C", "PLCXD3", "SRGAP1", "RANGAP1", "TBC1D8", "SLC4A7", "CDK14", "ZFP36", "USP53", "EPS8", "ANXA1", "SMAD9", "PFKFB3", "EBF1", "ABCA1", "RNF220", "FLT1", "CLDN5", "SLC38A5")
known_veinous <- c("TSHZ2", "AFF3", "SNTG2", "MCTP1", "KCTD8", "ADAMTS6", "MKL2", "RALCAPA2", "PDE4D", "MTUS2", "PLA2G4A", "FLI1", "NAMPT", "CEMIP2", "ELOVL7", "PDLIM1", "CFH", "CLDN5", "THSD7A", "SLC26A2", "TBC1D4", "THSD4", "PECAM1", "SLC16A1", "VEGFR-2", "TSHZ2", "ADGRG6", "MYOCD", "CD74", "TSHZ2", "LRRC1", "BNC2", "ETV6")
known_capilary <- c("FLT1", "ATP10A", "CLDN5", "SLC7A5", "SLC7A1", "ANO2", "CGNL1", "VWF", "MECOM", "BSG", "PTPRB", "ELOVL7", "THSD4", "ERG", "PODXL", "SLC16A1", "PECAM1", "SLC2A1", "VEGFR-2", "MFSD2A", "SLCA5", "TFRC")
known_arterial <- c("ARL15", "VEGFC", "MECOM", "THSD7A", "BMP6", "CP", "THSD4", "MCTP1", "ALPL", "ELOVL7", "SMAD6", "GFOD1", "FLT1", "CRIM1", "COL8A1", "CLDN5", "PLCG2", "ERG", "PECAM1", "LAMA3", "PELI1", "ADAMTS6", "PTPRB", "VEGFR-2", "VEGFC", "ARL15", "BMX", "EFNB2", "HEY1", "GATA2")
known_lymphatic <- c("PDPN", "LYVE1", "PROX1", "AQP1", "CCBE1", "PDPN", "PECAM1", "CD31", "PROX1", "LYVE1", "FLT4", "CCL21", "CD45", "PTPRC", "GATA2", "EFNB2", "FOXC2", "ITGA9", "FN1", "GJA1", "GJC2", "CCBE1", "VEGFR-2", "CD96")

known_Epithelial <- c("FLT1", "NKCC1", "NBCE2", "TTR", "ZO1", "CLIC6", "ARL13B", "GLUT1", "CTNNB1", "FOXJ1", "KCNJ13", "ENPP2", "HTR2C", "PDGFA", "PDGFB", "IGF2", "IGF1R", "ROR1", "EFNA5", "LINC00276", "HTR2C", "TMEM72-AS1", "GRM8", "PCAT1", "CP", "SNTB1", "PIP5K1B", "HDAC9", "CHST11", "ABCC5", "PLD5", "MTUS2", "TRPM3", "OTX2-AS1", "DYNC1I1", "LINC00880", "CLIC6", "SLC24A4", "CDH12", "ATP2B2", "CRPPA", "ROR2", "AMMECR1", "PLXNA2", "GNA14", "FYB2", "DOCK8", "SHANK2", "CRYBG1", "RYR3", "DNAH11", "LARGE1", "POU2AF1", "TRMT9B", "DGKI", "KIAA1217", "CUX2", "SLC4A10", "AGBL4", "NWD1", "PRDM16", "GLIS3", "AL008633.1", "ANKRD30BL", "GALNT17", "ABCA4", "COL9A1", "TMEM117", "ESRRG", "SLC4A5", "C2orf92", "AC079793.1", "LINC00598", "AC019197.1", "AL161757.4", "FBXO32", "HTR2C")
known_myoepi <- c("TAGLN", "TPM1", "MYL9", "MYL12A", "KRT17", "GBP1", "CALD1", "ANXA1", "ACTA2", "CTGF", "TUBB4B", "FAM183A", "RSPH1", "TPPP3", "MUC12", "CDK2AP2", "CCNO", "AQP1")
known_ependymal <- c( "GFAP", "FOXJ1", "PIFO", "CCDC15", "DYNLRB", "CFAP44", "ZMYND1", "LRRC23", "ZNF207", "ZNF800", "RPS6KC1", "EPHA4", "PHF20L1", "BCAS3", "USP32", "BTBD7", "GNAI2", "EGLN1", "PDXK", "SIK3", "DICER1", "HBS1L", "TRPC1", "TRPC4AP", "RASA2", "WASF2", "BCL2L1", "DPP10", "CTNNA2", "NPAS3", "FMN2", "CFAP299", "ADCY2", "SORBS1", "CASC15", "LINC00609", "ERBB4", "NEBL", "DTHD1", "CRB1", "RFX4", "DNAH9", "TTC6", "PTPRZ1", "ADGRL3", "DTNA", "CFAP47", "AC012405.1", "AGBL1", "AC020718.1", "ANKFN1", "LINC01088", "NRCAM", "PAMR1", "TNIK", "KCNMB2", "PLEKHA5", "ADGB", "TRDN", "TMEM232", "CSMD3", "LMNTD1", "MIR4300HG", "C8orf34", "LINC02055", "PPP2R2B", "ST8SIA1", "AL357093.2", "DPYSL3", "NAV2", "RPH3A", "HYDIN", "TTC34", "TSHZ2", "CFAP54", "MIR99AHG", "AC019330.1", "DNAH3", "KLHL32", "KIAA2012", "TTC29", "KCNMA1", "FYN", "GFAP", "GPM6B", "ASTN2", "PALLD", "CFAP46", "ANK2", "AC091078.1", "NTRK2", "DNAH10", "AC110023.1", "CHST9", "CCDC146", "GRAMD2A", "FUT9", "HPSE2", "PHACTR1", "CFAP44", "DNAH12", "ANXA1", "MUSK", "LINC00907", "WLS", "CFAP52", "SLC1A3", "AC013470.2", "KIAA0825", "SOX6", "BBOX1", "NEK10", "ULK4", "PPM1E", "AL022068.1", "CDH2", "PACRG", "SLC1A2", "FAM184A", "ERICH3", "FRMD5", "PPFIA2", "ZBBX", "DST", "FAT3", "NEK5", "SERPINI2", "ID4", "PLCB1", "SPAG16", "C8orf37-AS1", "PCDH9", "SGIP1", "DNAAF1", "KIAA1549L", "JAKMIP2", "RIPOR2", "PTCHD1-AS", "LRRIQ1", "LAMA1", "MYLK3", "TIAM1", "WDR49", "DYNC2H1", "PREX2", "CCDC200", "PALM2-AKAP2", "DCDC1", "KIF6", "PZP", "DCLK2", "NR2F1-AS1", "KCNN3", "DNAH7", "PANTR1", "LINC00511", "TTYH1", "AL157886.1", "SLC7A11", "LMO2", "AC093689.1", "RASSF9", "DCLK1", "BASP1-AS1", "PSD3", "CDH20", "DOC2A", "DNAH2", "SOX5", "FRMPD2", "CLVS1", "COL8A1", "LAMA2", "AC114971.1", "QKI", "NKAIN3", "C1orf87", "AC004949.1", "GYG2", "CFAP43", "VWA3B", "TMEM178B", "MAP3K19", "UBE3D", "LINC00461", "LINC02649", "FAM135B", "RSPH1", "ASB3", "VWA3A", "SPAG17", "ARAP2", "RNF43", "GMPR", "TMEM163", "EFHB", "FBXO15", "DNAH6", "MMP16", "ZBTB7C", "FABP6", "AC073941.1", "CHL1", "LMO3", "STK32A", "CFAP70", "CFAP206", "ZNF709", "ANGPT1", "EMP1", "LIFR", "TOGARAM2", "IQGAP2", "KIAA0319", "CTNND2", "PLEKHG1", "NFASC", "BMPR1B", "AL133304.3", "ITGA2", "STON2", "FGF14", "CD109", "SCN1A", "MLIP", "MAP1B", "FHAD1", "PARD3", "RFTN2", "SLC14A1", "MEIS2", "LRRC9", "NELFA", "PHYHIPL", "NEK11", "GPC6", "CCDC30", "SGCD", "CFAP100", "CCDC88C", "CFAP299")

known_ChPMature <- c("EMX2", "OTX1", "SIX3", "HOXA", "HOXB", "EN2", "MEIS1", "KRT18", "TTR", "KRT16", "NME5", "CLDN1", "CLDN3", "CLDN5", "AQP1", "COL3A1", "COL1A1", "MGP", "COL1A2", "LGALS1", "IGF2", "VCAN", "LUM", "POSTN", "DCN", "FN1", "COL6A2", "DLK1", "COL6A3", "COL5A2", "COL5A1", "CDH11")
known_dark <- c("NOS3", "AMPK", "ENDOG", "TFAM", "SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "PPARGC1A", "GABPA", "NRF1", "NRF2", "TFAM", "PPARG", "ADH1B", "ALDH2", "ENO1", "EPM2A", "GAPDH", "LDHA", "MPC1", "PDHA2", "PFKFB4", "PHKG2", "MT-CO1", "MT-ND4", "MT-ND1", "RARRES2", "IGFBP7", "MT-CO3", "STRIP2", "CXCL14", "AQP1", "CLDN5", "SLC26A3", "PRLR", "GPX3", "MT-CO3", "MT-ND4", "MT-CO2", "MT-ND3", "MT-ND1", "MT-ATP6", "MT-CO1", "MT-ND2", "MT-CYB", "MT-ND5", "TTR", "MT-ND4L", "CTSD", "MT-ATP8", "CLDN5")
known_light <- c("MEIG1", "FOXJ1", "DYNLRB2", "PIFO", "SHISA8", "CCDC67", "CCNO", "MCIDAS", "DEUP1", "C20orf85", "TUBB4B", "FAM183A", "RSPH1", "TPPP3", "ROPN1L", "MUC12", "CDK2AP2", "CCNO", "CARD19", "FOXJ1", "CCDC67", "KRT18")
other_dark <- c("CLDN5", "RARRES2", "FOCAD", "HTR2C", "STRIP2", "ATP2B3", "PRNP", "SLC16A10", "AQP1", "CXCL14", "SLC5A3", "FBLN1", "CP", "WFIKKN2", "IGFBP7", "FOLR1", "CYP27A1", "KCNE2", "KCNK1", "SERPINF1", "SLC13A4", "STIM2", "MT-ND1", "LAMB1", "NTRK2", "SLC4A5", "SOSTDC1", "PLBD1", "GPM6A", "ANXA4", "GALNT11", "RAMP3", "CA12", "ACOT11", "SLC2A1", "CLIC6", "MT-ND4", "VAT1L", "SLC2A12", "SLC20A2", "ABHD2", "SLC16A6", "ENPP2", "CYP1B1", "VEGFA", "MT-CO1", "GDF15", "FXYD1", "SYNGR1", "MT-ND2", "NEAT1", "PIP5K1B", "SLC1A3", "ATP2B2", "TTLL4", "MT-CO3", "C2orf40", "ETS2", "ITM2A", "MT-ND6", "MT-ND4L", "SYNE2", "RAB3GAP2", "CGNL1", "MT-ND3", "MT-ND5", "DNER", "SLC15A4", "ECEL1", "WFDC2", "FZD8", "SYNE1", "MARC2", "TSPO", "TMEM117", "SOD3", "KCNJ13", "ART3", "SPTBN1", "KIAA1456", "EPAS1", "COBLL1", "ATP2C1", "LGALS3BP", "MT-CO2", "MT-ATP6", "HLA-C", "TMEM185A", "ZADH2", "NENF", "WIF1", "GPNMB", "SLC23A2", "LPGAT1", "CLSTN1", "ARL4A", "TNFAIP8", "TYR", "SLC22A17", "IGFBP5", "SLC16A12", "HPGD", "ATP1B1", "SLC39A12", "PLA2G16", "MTRNR2L1", "MT-CYB", "SLC12A2", "DPP6", "SLC4A10", "ROR2", "SPOCK3", "GALNT10", "PMCH", "SHISA5", "F5", "DDR2", "TRPV4", "PXDN", "NT5DC1", "RNASET2", "RBP1", "DNAH11", "TMEM14A", "CLDN1", "ATP1B2", "TSPAN11", "MTRNR2L8", "PPP1R1B", "DR1", "CST3", "PHACTR2", "PLAC9", "SLC6A15", "CA14", "PTPRG", "CTSD", "TANC2", "NDUFA4L2", "DKK3", "SLC4A2", "CHCHD10", "SLC7A8", "ISPD", "TIMP3", "GPR155", "SMC6", "AP1S2", "B2M", "LINGO1", "S100B", "PTGDS", "PLP2", "ZFYVE16", "IGF1", "TPD52L1", "LAMA5", "WNK1", "DPY19L3", "SLC41A1", "CTNNAL1", "SULF1", "FUOM", "SLC17A8", "LINC01088", "NWD1", "PFKL", "HOMER3", "QSER1", "SLC19A1", "PDE11A", "CCK", "ARRDC3", "RBBP8", "TMC5", "PRND", "IQGAP2", "ZDBF2", "RABGAP1L", "OGDHL", "CPXM2", "ASPH", "DPP7", "CLDN4", "LRRN2", "LHFP", "TMED5", "TMEM38A", "CD9", "MDM4", "KTN1", "TTLL7", "TMEM164", "ANKRD37", "COL18A1", "DCDC2", "PRDX4", "BACE2", "ELN", "LAMB2")

known_dividing <- c("CRABP2", "MDK", "LIM2", "MARCKSL1", "PLS3", "TPBG", "HMGN2", "APOE", "HIST1H4C", "MCM2", "MCM4", "MCM6", "MCM7", "CDC47", "CDK1", "CDK2", "DDK", "CDT1", "CDC25C")

known_ChPImmature <- c("AQP1", "PAX6", "OTX2", "RSPO3", "MSX1", "CLDN3", "DCT", "PMEL", "PTGDS", "MDK", "ID3", "TPBG", "GNG11", "APOE", "WNT2B", "APCDD1", "WNT4", "RSPO3")
known_epi_prog <- c("AQP1", "CCDC67", "HES1", "HES5", "HES3")
known_epi_neuro_prog <- c("PAX6", "RSPO2", "RSPO3", "WNT5A", "WNT3A", "WNT1", "FABP7", "MEIS1", "NEGR1", "OTX2", "WNT1", "LMX1A", "GLI1")
known_neuro_prog <- c("LAMA2", "TUJ1", "TUBB3", "SOX2", "SOX2OT", "NEUROD6", "TUBB3", "PAX6", "MKI67", "RSPO2", "RSPO3", "RTN1", "EOMES", "NEUROD1", "NEUROD2", "NGN2", "EMX2")

known_Neuron <- c("SNAP25", "SOX2", "TBR2", "SOX2", "MAP2", "CADM2", "SYP", "NETO1", "MINDY4B", "AMPH", "SV2B", "NAPB", "HCN1", "NR2E1", "ATP8A2", "PPFIA2", "LRRC4C", "CVOP", "MCF2L2", "RTN1", "SOX2", "SOX9", "RBFOX3", "DLX1", "DLX2", "HOXA2", "HOXB2", "SCG2", "VGL", "SCG3", "CHGB", "CHGA", "MEF2C", "SATB2", "DCX", "GRID2", "GRIK2")
known_neuro_Exc <- c("GRID2", "GRIK2", "MAPK1", "MAP2K1", "MAP2K2", "RAB21", "SNX12", "ALG6", "ALG3", "SLC17A7", "CUX2", "RSPO1", "SCUBE1", "FEZF2", "NDST4", "NXPH4", "HRST4", "TSHZ2", "CHAT", "PTGRU" , "SLC17A7", "CUX2", "COL5A2", "RORB", "GAL", "THEMIS", "NVIA", "PCP4", "NR4A2", "NTNG2", "CTGF", "FEZF2", "SYT6", "SEMA3D", "VSNL1", "SNAP25", "DNM1")
known_neuro_Exc_subS <- c("SATB2", "TYRO3", "ARPP21", "SLC17A7", "TBR1", "CAMK2A", "ITPKA")
known_neuro_Inh <- c("SEMA3A", "SEMA3D", "GABRP", "PDGFD", "LHX6", "PVALB", "LHX6", "ADARB2", "LAMP5", "CALB2", "GAD1", "GAD2", "CXCL14", "CALB2", "VIP", "SCUBE3", "SLC6A1")
known_neuro_inh_subS <- c("SST", "PVALB", "LAMP5", "VIP", "GAD2", "SYT6", "CPNE7", "LHX6", "ADRB2", "CALB2", "LHX6", "PVALB", "LHX6", "SST", "ADARB2", "LAMP5", "CALB2", "GAD1", "GAD2", "CXCL14", "CALB2", "VIP", "SCUBE3")
known_neuro_inh_CGE <- c("ADARB2", "PROX1", "SV2C")
known_neuro_inh_MGE <- c("SOX6", "RELN", "LHX6")

known_astro <- c("AQP4", "MOBP", "SLC1A2", "SLC1A3", "FABP7", "ALDH1L1", "ALDOC", "S100B", "LAMA2", "SLC4A4", "AQP4", "GLUT1", "RBFOX2", "H2AZ1", "LRRN3", "BTBD17", "ZFYVE21", "GPR37L1", "ACSL6", "MMD2", "RFX4", "P3H4", "SOX11", "DKK3", "CCDC80", "SCD", "CD9", "SOAT1", "VIM", "PLCD1", "SPARCL1", "C2ORF72", "MT1E", "TECPR2", "RAMP3", "IFITM3", "CPNE5", "PLEKHB1", "SERPINE2", "ANOS1", "PON2", "LRP4", "PRCP", "ABAT", "AQP4", "TFPI", "ARNTL", "GNG5", "PREX1", "GNG7", "DAG1", "S100A13", "S100A16", "APBB2", "HES1", "B2M", "HES4", "HES5", "DTNA", "PAQR8", "BCL11A", "HLA-C", "HLA-A", "HLA-B", "GPR137B", "AIF1L", "MLLT11", "FERMT2", "CAN", "KLHDC8A", "FAT1", "SFXN5", "TAGLN3", "ACSBG1", "ITM2C", "GPM6B", "LRRC10B", "IFI27L1", "SPON1", "PPP1R17", "BEX1", "TEX264", "TTYH2", "TTYH1", "TRIL", "MLC1", "PHYHIPL", "PTTG1IP", "SLC6A1", "ADCYAP1R1", "LIPA", "MT2A", "SOX4", "LGALS3", "GNG11", "SOX2", "GJA1", "NKAIN4", "LGALS1", "ALCAM", "RASSF4", "STMN1", "TSPAN7", "TSPAN6", "SOX9", "GLUL", "APOE", "TSPAN3", "SOX8", "LYN", "BBS2", "AASS", "JUN", "TSC22D4", "ANXA5", "TUBB", "GNG12", "RGMA", "RHOC", "GLUD1", "FXYD1", "CLDN10", "PTMA", "ALDH6A1", "GPRC5B", "P4HA1", "ID1", "TMSB15A", "VEPH1", "GNPTG", "ID2", "ID4", "FAM181B", "ALDOC", "CD24", "MGLL", "DDR1", "SDC3", "ABHD3", "ELN", "PEA15", "LITAF", "TOR1AIP1", "SNX5", "SNX3", "UCHL1", "MAGT1", "CSRP2", "NEUROD6", "CTSL", "TIMP3", "TIMP1", "PHGDH", "IL33", "ELOVL2", "CBR1", "ELOVL5", "SLC39A12", "F3", "PRSS35", "SERPINB6", "HADHB", "RIT2", "GATM", "FABP5", "LRATD2", "BDH2", "FABP7", "PLIN3", "SCRG1", "RAMP1", "PNPLA3", "ROBO3", "OAT", "CD63", "NCAN", "LRRC17", "STON2", "GDPD2", "PTN", "CST3", "UBL3", "HEY1", "DNER", "MRC2", "PLAGL1", "BCAP29", "METTL7A", "TBC1D10A", "FADS2", "CA12", "PDGFRB", "BCHE", "NPNT", "PPP1R16A", "DIO2", "C1ORF122", "SHISA4", "LIFR", "ATP1B2", "VCAM1", "AGT", "YIF1A", "OLFM2", "RAB34", "RAB31", "DOK5", "TP53I3", "PSAT1", "SYPL1", "LIX1", "CRYL1", "PMP22", "CSPG5", "DHCR7", "MASP1", "CD82", "SLC1A2", "SLC1A3", "ZFP36L1", "TNC", "SLC1A4", "CLU", "FBLN2", "ATP1A2", "CXCL14", "COMT", "NDRG2", "ZFP36L2", "NT5E", "PDPN", "EMC2", "APOE", "ALDOC", "GFAP")
known_oligo <- c("VCAN", "LAMA2", "MOBP", "MOG", "LIF", "OMG", "OLIG1", "IGF1", "OLIG2", "MYT1", "SOX10", "BMP4", "IL1B", "BMP2", "MAG", "PLP1", "CNTF", "CNP", "PDGFB", "CXCL2", "TNF", "CXCL1", "FGF2", "GLI2", "ASCL1", "SHH", "NKX2-2", "SOX9", "MBP", "NKX2-6", "SOX8", "SOX5", "SOX6", "PLP1", "PLP1")
known_OPC <- c("CLK1", "FTO", "BCAS3", "IWS1", "MAP3K3", " PRKAR2A", " TTC19", "CDK19", "FOXK2", "SIK23", "RAB7A", "CNTLN", "HIF1A", "PAG1", "PACS2", "PDGFRA", "SOX10", "NKX2.2", "FGF12", "RTN1", "NKX2-2", "CALM1", "RAB31")


known_Immune <- c("CD163", "CD247", "CD45", "CD11B", "CD206", "CD163", "PU1", "IBA1", "LYVE1", "APOE", "LY86", "P2RY6", "HMHA1", "PTPRC", "PTPRC", "CD3E", "SP11", "CD14", "CD45", "RNASE1", "C1QA", "C1QC", "FCAR", "GPNMB")
known_Macrophage <- c("PTPRC", "CD68", "LYVE1", "PECAM1", "F13A1", "MERTK", "IGF1", "CD163", "SOLR1", "C3", "GRID2", "FCGRB2", "CR1", "P2RY12", "MRC2", "SIGLEC1")
known_mglia_M0 <- c("LAMA2", "TMEM119", "APOE", "CD33", "CD45", "CD43", "CD68", "MS4A10", "MS4A4A", "MS4A4E", "MS4A6A", "MS4A5", "MS4A8", "MS4A12", "MS4A6E", "MS4A3", "MS4A1", "IBA1", "CD14", "CXCR1", "GPR34", "MERTK", "C1QA", "C1QB", "C1QC", "C1QBP", "SBBB or BCSFB1", "GLUT5", "GAS6", "CD68", "TMEM119", "P2RY12", "AIF1", "CX3CR1", "CST1R", "MERTK", "SPARC", "IRF8", "SBBB or BCSFB1", "SPARC", "CD11B", "CD45", "CRYBB1", "SALL1", "TMEM119", "FCRLS", "APBB1IP", "PTPRC", "CX3CR1", "SALL1", "FCRLS", "SLC2A5", "CSF1R")
known_mglia_prog <- c("CD45", "CD43", "CD34", "SPI1", "MKI67", "IFIT3")
known_mglia_M1 <- c("CX3CR1", "IL1B", "TNF", "STAT3", "NFKB", "GSK", "ROS", "RNS", "CCL3", "CCL4", "CXCL16", "MIF", "H2D1", "APOE", "LAGALS3BP", "TOP2A", "UHRF1", "RRM2", "RAD51", "CHAF1B", "POLE", "RAD1AP1", "CCNE1", "HELLS", "GMNN", "CHAF1A", "CCNF", "CDC45", "PKMYT1", "SPC25", "PLK1", "NUSAP1", "NDC80", "KNSTRN", "KIF23", "KIF22", "CENPE", "CDCA8", "CDCA2", "AURKA", "ANLN", "CXD10", "OASL1", "ATF3", "EGR1")
known_mglia_M2 <- c("SELPLG", "MED12L", "BIN1", "CX3CR1", "P2RY13", "TMEM119", "CD44")
known_mglia_Dis <- c("TREM2", "H2-AA", "H2-AB1", "CD74", "IFITM3", "IRF7", "OAS1A", "RSAD2", "ZBP1", "IFIT2", "MX1", "CD40", "OASL1", "CXCL10", "H2-D1", "H2-Q4", "H2-Q5", "CTSB", "CTSZ", "CCL12", "CCL3", "CCL4", "CXCL16", "TLR2", "LILRB4", "PKM", "PGK1", "GAPDH", "PGAM1", "CLEC7A", "LGALS3", "AXL", "CFB", "C4B", "C3", "CD83", "TNFAIP3", "C5AR1", "CD68", "GPNMB", "ABCA1", "CXCR4", "CLEC7A", "CST7", "LPL", "APOE", "CD274", "CLEC7A", "ITGAX", "CD274", "IL1B", "CCL5", "CXCL10", "MMP14", "CD68", "GPNMB", "RSP26", "FTH1", "SPP1", "CCL5", "CXCL14", "CCL6", "CD68")
known_BAM <- c("CD163", "PTPRC", "AIF1", "CD68", "ADGRE1", "P2RX7", "CD206", "FOLR2", "NRP1", "CD63", "MHCII", "CLEC12A", "CD45", "ADGRE1", "FCGR1", "AIF16", "APOE", "MS4A7", "MS4A6C", "LYZ2", "TGFBI", "MS4A7 ", "MS4A6C", "LYZ2", "CLEC4A1", "LFITM2", "TGFBI", "CYBB", "MINDAL", "IFI27I2A", "IFITM3", "CD36", "ALDH2", "PLA2G7", "PF4", "MSRF", "AOAH", "DAB2", "CLEC12A", "MS4A7", "MRC1", "LGALS3", "CD163", "PF4", "CD68", "CD163", "CD206", "LYVE1", "ADGRE1", "CD74", "PF4", "MRC1", "CD68", "APOE", "MS4A7", "CD38", "CX3CR1")
known_BAM_prog <- c("S100A6", "HSPA1A")
known_BAM_strom <- c("CCND2", "LILRA5", "TTR", "ENPP2", "COL14A1", "CST7", "GM1673", "CLEC7A", "APOE",  "MS4A7", "MS4A6C", "LYZ2", "TGFBI")
known_BAM_epi <- c("LYVE1", "CST7", "GM1673", "ST8SIA6", "ICAM1", "LPL", "CLEC7A", "PSAT1", "CD63", "CD33", "CD45", "CD43", "CD68", "MS4A10", "MS4A4A", "MS4A4E", "MS4A6A", "MS4A5", "MS4A8", "MS4A12", "MS4A6E", "MS4A3", "MS4A1")
known_BAM_MHCIIlo <- c("CD163", "MRC1", "LYVE1", "MSR1", "SIGLEC1", "IRF8", "IFITM2", "IFITM6", "H2-AA", "CD74", "KLRA2", "NRP1", "BLVRB", "SMAGP", "VCAM1", "MRC1", "DSE", "CD163", "CD38", "GPX3", "CP", "STARD8", "MAF", "STAB1", "CBR2", "IGF1", "IGBP4", "FCG", "SNX2", "F13A1", "GAS6", "NINJ1", "CD163", "MRC1", "LYVE1", "MSR1", "SIGLEC1", "IRF8", "IFITM2", "IFITM6", "CLEC4N", "CLEC10A", "FOLR2", "CD7", "CD209F", "FOLR2", "TSLP", "CCL24", "LRP6", "FCNA", "CCL2", "SLC40A1")
known_BAM_MHCIIhi <- c("MHCII", "CD38", "MRC1", "CD163", "GAS6", "CD74", "AXL", "GPR65", "CLEC12A", "CD52", "CXCL16", "H2-AA", "H2-EB1", "H2-DMB1", "FXYD5", "FGL2", "SCIMP", "CD72", "FCGR4", "PLBD1", "KLRA2", "MHCII", "CD38", "LYVE1", "P2RX7", "EGFL7", "IRF7", "CRIP1", "CCL6", "CCL9")
known_Macro_M0 <- c("MMP9", "MARCO", "CCL20", "CD16")
known_Macro_M1 <- c("CD36", "STAT1", "CSF2RB,", "IL15", "IL2RA" ,"IL6R", "CD38", "CD69", "CD97", "ICAM1", "ITGAL", "ITGA4", "ITGB7", "MUC1", "SIAT1", "JAK2", "IRF1", "PTPRC", "PTPRO", "NMI", "CISH", "SOCS1", "CXCL10", "CXCL11", "CSF2", "STAT5", "ERK", "AKT", "IRF5", "NFKB1", "CXCL8", "IL6", "CSF3", "CSF1", "TNF", "IL1B", "CD14", "CD64", "SOCS1", "NOS2", "ARG1", "IL10", "IL4", "TNF", "IL17A", "CD274", "TGFB", "STAT6", "MMP9", "CXCL3", "IL18", "FCGR1A", "CD64", "CCL3", "CXCR4", "CD274", "CXCL10", "IFT3", "IL1B", "MMP14", "CCL5")
known_Macro_M2 <- c("HLA-DPB1", "HLA-DPA1", "HLA-DRA", "CCL3", "CCL4", "CCL5", "C1QA", "C1QB", "C1QC", "FTL", "SPP1", "TGM2", "MRC1", "CH25H", "PTGS1", "KLF4", "CISH", "SOCS1", "CD64", "CD32", "CD16A", "PI3K", "CD16B", "C1QA", "DSIPI", "THBS1", "CD163", "CXCL13", "CXCL4", "FPR1", "TLR1", "TLR8", "FOXP3", "TGFB", "CD86", "TREM2", "FCGR1", "FCGR2B", "FCGR3", "ROCK1", "MERTK", "CD14", "LRP1")
known_Macro_Dis <- c("GBP1", "IFI6", "CD40", "TNF", "CCL3", "CXCL8", "SLAM7", "CD14", "RPL36A", "FCG3A", "IGHG1", "LAIR2", "ATF3", "RGS1", "CCR2", "LYZ", "CD86", "RFX5", "RFXANK", "RFXAP", "SPI1", "TUNX1", "NFKB1", "RELA", "FOXL2", "FOXP1", "MXI1", "RUNX2", "MEEIS1", "KLF5", "PBX3", "AHR", "ARNT", "SOX2", "TCF12", "MEF2A", "TEAD1", "GPNMB", "FABP5", "MMP14", "CCL5", "LGALS3", "H2AA", "CD74", "CD274", "MMP14")
known_monocyte <- c("CD4", "FCER1A", "CD68", "CD14", "CD16")

known_plasma_dendritic <- c("CD4", "FCER1A", "CD68", "CD14", "CD16", "SIGLECH", "CCR9", "PASCIN", "ITGAX", "THBD", "IL3RA") #NO ZBTB46
known_conv_dendritic <- c("CD4", "FCER1A", "CD68", "CD14", "CD16", "FLT3", "IRF8", "XCR1", "CD209", "ZBTB46", "CD14", "ITGAX", "THBD", "IL3RA")
known_mast <- c("GATA2", "KIT", "HPGD5")

known_memory_T <- c("CAMK4", "USP36", "H3F3B", "NR4A2", "HSP00AA1", "SGCD", "KLRG1", "LEF1", "PVT1", "SELL", "CCR7", "LEF1", "TCF7", "KLF2", "CD3")
known_effector_T <- c("THEMIS", "CD3D", "CD4", "CD8A", "CD8B", "FOXP3", "NKG7", "CD44", "CD4", "SELL", "CCR7", "IL7R", "CD8", "CD8A", "CD8B", "GZMK", "CXCR3")
known_resident_T <- c("CD3D", "CD4", "CD8A", "CD8B", "FOXP3", "NKG7", "CD44", "CD69", "NR4A2", "IL7R", "CD8", "PRF1", "NKG7", "ZNF683", "GZMB", "CD8A", "CD8B")
known_NTK <- c("KLRF1", "NCAM1", "GNLY", "NKG7", "CD44", "PRF1", "NKG1", "GZMB", "KLRD1", "KLRF1")
known_plasma_B  <- c("CD19", "JCHAIN", "CD44", "CD68", "CD4", "CD79A", "MS4A1", "MHC", "IGHG3", "IGHA1", "DERL1", "FKBP11")


Cell_Stress <- c("CXCL8", "CCL2", "EIF2AK4", "HERPUD1", "HSPA5", "NF2L2", "EIF2S1", "IFBP1", "ATF4", "PGK1", "GORASP2", "ARCN1", "HSPB7", "HSPG2", "HSPB11",  "HSPA6", "HSPD1", "HSPE1", "HSPBAP1", "HSPA4L",  "HSPB3", "HSPA4", "HSPA9", "HSPA1L",  "HSPA1A",  "HSPA1B",  "HSP90AB1", "HSPB1", "HSPA5", "HSPA14",  "HSPA14.1", "HSPA12A", "HSPB2", "HSPA8", "HSP90B1", "HSPB8", "HSPH1", "HSPA2", "HSP90AA1", "CARHSP1", "HSPB6", "HSPBP1", "HSPA12B", "HSPA13", "PGK1", "ARCN1", "GORASP2", "HIST1H1C", "HIST1H2BC", "UBC", "JUND", "RGS1", "HSPA1A", "HSP90AA1", "CCL4", "DUSP1", "HSPA1B", "CCL3", "RHOB", "JUN", "ZFP36", "KLF2", "JUNB", "FOS", "TXNIP", "EGR1", "ADAMTS1", "BTG2", "WFDC21", "IER5", "ATF3", "HIST1H4I", "GEM", "IER2", "IER3", "HIST1H2BR", "HIST1H1E", "IL1B", "SERPINE1", "NFKBID", "FOSB", "EGR2", "CITED2", "KLF6", "NFKBIZ", "HIST2H2AA1", "HIST1H4D", "GM26532", "HSPA1B", "HSPA1A", "HIST1H2AC", "RHOB", "UBC", "RGS1", "SLC2A3", "HSPH1", "DNAJB1", "DUSP1", "HMOX1", "GLUL", "HSP90AA1", "BCAS2", "SRGN", "GPR183", "RGS16", "GADD45B", "PPP1R15A", "HIST1H3H", "HIST1H1D", "HIST1H2BH", "HIST1H2BJ", "HIST1H2BG", "ZFP36", "HSPB1", "HIST1H2BF", "HIST2H2BE", "P4HA1", "HIST1H2BD", "ZFP36L1", "BTG2", "DDIT4", "FOS", "HIST2H2AA4", "HSPE1", "JUN", "HIST2H2BF", "GFAP", "SPARC", "ID3", "ATF3", "CRYAB", "LINC01094", "JUNB", "FOS", "NR4A2", "DCLK1", "HSP90AA1", "IGFBP5", "MIR4300HG", "SLC38A1", "NR4A1", "JUN", "NEAT1", "PPP1R3C", "ID2", "HSPA1B", "HSPB8", "SFRP2", "UBC", "HSPB1", "DUSP1", "DDIT4", "HSPA1A")
whole_cells <- c("KCNJ8", "KCNJ11", "ABCC8", "ABCC9")

Whyss_Correy_Epithelial <- c("LINC00276", "HTR2C", "TMEM72-AS1", "GRM8", "PCAT1", "CP", "SNTB1", "PIP5K1B", "HDAC9", "CHST11", "ABCC5", "PLD5", "MTUS2", "TRPM3", "OTX2-AS1", "DYNC1I1", "LINC00880", "CLIC6", "SLC24A4", "CDH12", "ATP2B2", "CRPPA", "ROR2", "AMMECR1", "PLXNA2", "GNA14", "FYB2", "DOCK8", "SHANK2", "CRYBG1", "RYR3", "DNAH11", "LARGE1", "POU2AF1", "TRMT9B", "DGKI", "KIAA1217", "CUX2", "SLC4A10", "AGBL4", "NWD1", "PRDM16", "GLIS3", "AL008633.1", "ANKRD30BL", "GALNT17", "ABCA4", "COL9A1", "TMEM117", "ESRRG", "SLC4A5", "C2orf92", "AC079793.1", "LINC00598", "AC019197.1", "AL161757.4", "FBXO32", "XIST")

Whyss_Correy_Epithelial_highMT <- c("MT-CO3", "MT-ND4", "MT-CO2", "MT-ND3", "MT-ND1", "MT-ATP6", "MT-CO1", "MT-ND2", "MT-CYB", "MT-ND5", "TTR", "GPX3", "MT-ND4L", "CTSD", "MT-ATP8", "CLDN5")

Whyss_Correy_Epithelial_highGPX3 <- c("TTR", "GPX3", "CTSD", "BSG", "SOD3", "CST3", "SELENOW", "RPL13", "SFRP1", "TMBIM6", "EEF1A1", "RPS28", "SERPINF1", "RPL28", "FBLN1", "CLDN5", "RPL18A", "SERF2", "SLC25A6", "PPP1R1B", "MT-CO2", "RPL13A", "PSAP", "JUND", "GSTP1", "ACTB", "NPC2", "CD59", "FOLR1", "METRN", "AQP1", "PCSK1N", "ATP6V0C", "SLC4A2", "RPS11", "GAPDH", "IER3", "CA2", "WFIKKN2", "RPL3", "HPD", "SNHG25", "PEBP1", "RPS19", "RPL36", "ECRG4", "RPS5", "RPS16", "METTL7A", "COX8A", "GPX4", "CLDN2", "SLC22A17", "RPS27", "PBXIP1", "PDK4", "PRNP", "B2M", "RPL8", "C4orf48", "RPL21", "RPS17", "GUK1", "RPS15", "FXYD1", "EEF2", "ATP5F1D", "NFE2L1", "CHCHD10", "CLU", "LYPLA2", "RNASEK", "RPLP1", "UBA52", "TUBB4B", "RPS14", "CALM1", "MSX1", "RPS27A", "UQCR11", "RPL15", "RARRES2", "MIF", "BEX3", "UQCRB", "RPL6", "RPL41", "ATP5F1E", "COX7C", "FABP5")

Whyss_Correy_Epithelial_highSLC26A3 <- c("LINGO1", "SLC26A3", "AC060765.2", "RAPGEF6", "RASGEF1B", "LINC00486", "MT-CO1", "MT-CO2", "GPX3", "FP236383.3", "MT-ND2", "CTSD", "GREM1", "BCOR", "MT-ND6", "HSBP1L1", "INO80D", "LINC00893", "DHFR", "AC007906.2", "LINC01358")

Whyss_Correy_Mesenchymal <- c("SLC4A4", "PTGDS", "CEMIP", "THSD4", "KCNMA1", "AC008014.1", "PHLDB2", "BNC2", "PDE7B", "SLC7A2", "GPC6", "CDH11", "LEPR", "SLC24A3", "EYA2", "CACNA2D3", "LAMA2", "GMDS", "TMOD1", "UNC5C", "APOE", "STEAP4", "FN1", "CHSY3", "BICC1", "SMOC1", "FAM20A", "ALPL", "SLC38A2", "SPTLC3", "ZNF804B", "PDE10A", "IGFBP6", "FOXP1", "EBF1", "AKAP12", "SLC26A2", "CLMP", "CYP1B1", "ADAM12", "IGFBP4", "GNAL", "LRP1B", "NAMPT", "TIMP1", "RARRES1", "FOXP2", "CAMK2D", "USP53", "PTPRE", "IGFBP5", "MGP", "PDZRN4", "KCNT2", "C3", "SORBS1", "ABCA8", "FOXC1", "C8orf34", "NTRK3", "TGFBR3", "MIR100HG", "NR2F1", "AC079298.3", "COL11A1", "SLC1A1", "TMTC1", "NEBL", "POU6F2", "RNF24", "BMPR1B", "PLXNA4", "SLC22A23", "OSMR", "LINC01482", "CFH", "FTL", "WWTR1", "MPP6", "NCOA7", "NDRG1", "DLG2", "SLC47A1", "MAST4", "EGFR", "PLA2R1", "CLU", "LRP1", "LDB2", "AUTS2", "LRFN5", "PIEZO2", "PLPP3", "ZIC1", "TNXB", "EPB41L2", "PDE1A", "SVEP1", "WDR86", "ZFPM2", "ITGA8", "KLF9", "NTN1", "JAM2", "COL1A2", "PRKCH", "LDLRAD3", "ALCAM", "SULT1E1", "CLDN11", "ID1", "MID1", "SH3RF3", "MAP3K5", "PKP2", "NNMT", "SESTD1", "GSN", "SLC2A3", "GLCE", "UACA", "NGF", "PLD1", "ENPP6", "RIPOR2", "RDH10", "BMPER", "IFITM3", "TSC22D3", "C1S", "BMP5", "NUDT4", "NCAM1", "SLC9B2", "NFKBIZ", "BCL6", "FILIP1", "SULF2", "WBP1L", "MCC", "ITGA9", "ARHGAP21", "CPED1", "ADAMTSL3", "LPAR1", "FAT4", "SH3BP5", "SH3RF1", "OGN", "SFRP2", "LINGO2", "PRRX1", "TJP1", "TBX15", "CEBPD", "NRP2", "SDC2", "SCUBE1", "EHBP1", "RHBDD2", "EVA1C", "IL4R", "ARHGEF28", "MEG3", "KLF5", "FP236383.3", "IL6ST", "PITPNC1", "LINC01197", "PLEKHA5", "ATP13A3", "AFAP1L1", "RBMS3", "STAT3", "ACADL", "NPNT", "SCARA5", "FMNL2", "SRGN", "PRKG1", "FBXL7", "PARD3", "RBMS1", "MALAT1", "TTC28", "PCDH9", "FTX", "ZBTB16", "PARD3B", "C7", "CLIC4", "TGFB2", "SVIL", "TBC1D8", "SERPING1", "DNASE1L3", "PTGDR", "JUNB", "PTN", "ATP1B3", "C1R", "LMO4", "TXNIP", "ROCK2", "PAM", "DUSP1", "FCHSD2", "MYLK", "SLC7A14-AS1", "FNDC3A", "CCL26", "KLF2", "STXBP6", "GAS1", "JCAD", "NR2F1-AS1", "ABCA9", "CD55", "ISLR", "ABCB1", "PRDM6", "ENPEP", "LAMA4", "SNX9", "SETBP1", "FYN", "FTH1", "NAV3", "ZFHX3", "SLPI", "NBL1", "TBX18", "PLCB1", "CPQ", "APCDD1", "MAP1B", "CALD1", "AHI1", "SPIDR", "SLIT3", "BTG2", "EYS", "AFDN", "CYP1B1-AS1", "NFE2L2", "COL6A2", "RASAL2", "ANO6", "PDE8B", "MEG8", "FXYD5", "HECW1", "A4GALT", "PRELP", "ARNT2", "KCNE4", "EYA1", "MEDAG", "APOD", "IGF2R", "ELL2", "LIMA1", "KANSL1L", "MRC2", "TRABD2B", "MCTP2", "NMNAT2", "SH3PXD2A", "PDLIM1", "FBN1", "NALCN", "PAQR5", "ME1", "GRK5", "GJB6", "EPHA7", "AEBP1", "MDGA2", "MN1", "TPT1", "AOX1", "MSC-AS1", "FSTL3", "PPFIBP1", "BX284613.2", "CFB", "PPM1H", "PRICKLE1", "ZFP36L2", "CRTC3", "GABARAPL1", "ACKR3", "AC018742.1", "TRIB2", "BACH2", "CYBRD1", "AC027288.3", "PID1", "TSPAN18", "NHS", "EGFLAM", "CCDC200", "SLC39A11", "GLI3", "PIAS1", "TLE2", "ELF2", "ITGBL1", "PLEKHA6", "ADAMTS1", "SIK2", "CDH23", "CEBPB", "ENG", "DPYD", "HIF3A", "FOSL2", "SMG6", "SSTR2", "RPS8", "SLC38A1", "COL8A1", "EPS8", "GASK1A", "AC092155.1", "PBX1", "OSMR-AS1", "ALDH1A2", "RBPMS", "B4GALT5", "ASS1", "B4GALT1", "TACC1", "SLC9A9", "AC092957.1", "URB1", "STK3", "MAML3", "ACTB", "LIFR", "SMCHD1", "NCOA3", "ZCCHC24", "MYO3A", "CCND2", "CNTLN", "PRKAG2", "ADCY5", "ITPR2", "PAWR", "MTMR7", "UTRN", "NID1", "RHOA", "UST", "SYBU", "SLC7A11", "ADAMTS3", "ANKS1B", "SEMA3C", "NUCKS1", "HIP1", "KLRD1", "CPD", "FGL2", "GNG4", "ABCD3", "HES1", "CFI", "BACH1", "SBNO2", "VIM", "PHACTR1", "NAV1", "ADIRF", "ZNF423", "SYNPO2", "ARHGAP10", "NR2F2", "C1RL", "GLUL", "GALNT1", "ACTN1", "ARHGAP29", "NEAT1", "SGCD", "ECE1", "SYNJ2", "TPM4", "RNF125", "C16orf95", "CADPS2", "YBX3", "ZNF516", "CDON", "AC002460.2", "ROBO1", "KANK2", "WSB1", "RANBP3L", "DSP", "KCNQ1OT1", "AF117829.1", "FGD6", "FKBP5", "ACTG1", "ARHGAP6", "GK5", "PAPSS2", "GJA1", "MAN1A1", "ADGRL2", "SMAD9", "TXNRD1", "TPST1", "NIPAL2", "NOVA1", "ZBTB7C", "CAMTA1", "SUCLG2-AS1", "ARHGAP28", "TLE5", "CYB5A", "DUSP5", "SLC13A3", "IFITM2", "ACSL3", "PPP3CA", "FLG-AS1", "AC087564.1", "TMSB10", "GNA12", "AC092691.1", "IL1R1", "ERBIN", "EDA", "F3", "PBX3", "ADAM9", "PTPRU", "ABCA1", "DOCK6", "GRIP1", "EPHX1", "RAB11FIP3", "SOCS3", "TSPYL2", "SLC35G1", "EPB41L3", "TSBP1-AS1", "FIGN", "TFRC", "KSR1", "ARHGAP31", "MEIS2", "MET", "SAT1", "SLC43A2", "RPS6KA3", "ETS2", "RHOB", "PC", "PGAP1", "ADAMTSL1", "EEF1A1", "CNTNAP4", "RHEB", "CCDC3", "RPL34", "EIF4A1", "MGST1", "TMSB4X", "PDGFC", "CCDC14", "BZW1", "CACNA2D1", "AL365295.1", "CDK14", "RAP1A", "RPLP1", "RPS6", "CPAMD8", "ITGA1", "RPLP2", "MIDN", "GLI2", "KLHL13", "NRG1", "PNRC1", "TBC1D16", "CHD3", "IKZF2", "RPS24", "ELOVL5", "RPL37", "RPL32", "DENND5A", "ID2", "HMGB1", "RPL12", "GPC5", "JAK2", "ADGRL3", "PLSCR4", "HIVEP2", "RPS14", "NRP1", "TTN-AS1", "H3F3B", "AAMDC", "CENPP", "CDC37L1", "STK38", "ADCY2", "NR4A2", "SYNDIG1", "SRGAP1", "TJP2", "IFNGR1", "DNAJA1", "RPS27", "PCBP3", "CMTM8", "RPL37A", "ADAMTS9", "TIPARP", "RPS18", "DAB1", "EEPD1", "CNTN4")

Whyss_Correy_Macrophage <- c("LRMDA", "DOCK4", "ARHGAP24", "SLC8A1", "RUNX1", "ELMO1", "FRMD4A", "TBXAS1", "CD163", "DOCK2", "CTSB", "MYO1F", "ACSL1", "HIF1A", "RNF149", "SSH2", "ETV6", "PLXDC2", "ARHGAP15", "AL163541.1", "SLC11A1", "STAB1", "MEF2A", "SLCO2B1", "HSP90AA1", "DISC1", "INPP5D", "EPB41L3", "FMN1", "APBB1IP", "BMP2K", "DPYD", "FYB1", "SPP1", "RAB31", "LNCAROD", "PTPRC", "RGS1", "SMAP2", "QKI", "ATP8B4", "SRGN", "ADAM28", "ARHGAP26", "MSR1", "ABR", "TLR2", "LYN", "MS4A6A", "PREX1", "CMIP", "SLC9A9", "FCGR2A", "KCNQ3", "TNFRSF1B", "MEF2C", "AOAH", "CYTH1", "LPAR6", "SFMBT2", "CELF2", "FKBP5", "SAT1", "ITPR2", "SORL1", "FAM49B", "PIK3R5", "RIN3", "ZSWIM6", "ZEB2", "LCP2", "MAML3", "F13A1", "DENND3", "SRGAP1", "P4HA1", "SRGAP2", "GNA13", "PFKFB3", "GRB2", "SLC1A3", "NCKAP5", "MBNL1", "HSPH1", "SAMSN1", "MGAT4A", "DOCK10", "IGSF21", "TAOK3", "RBPJ", "RHBDF2", "RNASET2", "PLAUR", "PICALM", "CIITA", "TCF12", "DLEU1", "VAV1", "IRAK3", "MYO9B", "SLA", "KYNU", "SYK", "MYO5A", "GNAQ", "FGD4", "LAPTM5", "SLC2A3", "CSF2RA", "AC093895.1", "ANKRD44", "PPARD", "HCLS1", "ARHGAP22", "FCGBP", "SKAP2", "GPRIN3", "NAIP", "IKZF1", "TBC1D14", "TRPM2", "CD53", "REL", "TBC1D22A", "FRMD4B", "MGAT5", "CD83", "ALOX5", "FLI1", "TBC1D5", "POU2F2", "AP2A2", "MERTK", "NCKAP1L", "PADI2", "CSF3R", "OLR1", "ADAP2", "PALD1", "SRGAP2B", "ZNF710", "UGCG", "CD86", "ANKRD11", "PIK3AP1", "ARHGAP18", "MRC1", "FCHSD2", "USP15", "EPB41L2", "UBE2E2", "GAB2", "RYR1", "WDFY4", "FMNL1", "CSGALNACT1", "USP36", "ELL2", "JARID2", "ADGRE2", "DOCK11", "HLA-DRB1", "CD163L1", "OSBPL3", "SH2B3", "SIPA1L1", "KMT2E", "RAPGEF1", "VMP1", "ITGAX", "HLA-DRA", "CYTH4", "ZC3HAV1", "PDE3B", "FAM49A", "SLC2A5", "ZFP36L1", "C22orf34", "SNX29", "KMT2C", "LIMS1", "IL10RA", "MAF", "TFRC", "RANBP2", "GPBP1", "MANBA", "IFI16", "RNF144B", "SNX6", "GRK3", "CEP170", "PRKAG2", "NAMPT", "THEMIS2", "BCL2", "UTRN", "HCK", "TNS3", "FNDC3A", "TGFBR1", "ZFAND3", "PEAK1", "AGFG1", "WDFY3", "RNF213", "DNM2", "CTSS", "NUMB", "SLC16A10", "ATM", "NIPBL", "SUSD6", "KLHL6", "EPS15", "MALAT1", "PTPN2", "PARP8", "PELI1", "CLEC7A", "IQGAP2", "RNF130", "EMB", "SLC7A5", "CHD1", "CYBB", "CSF1R", "TCF4", "GLUL", "MFSD1", "MED13L", "TANC2", "LAT2", "GNB1", "COP1", "FCHO2", "SDCCAG8", "TYMP", "IFNGR2", "ST8SIA4", "PAG1", "TFEC", "PHF20", "JAZF1", "MB21D2", "SLC4A7", "PARVG", "STAG1", "ARHGAP25", "RAB1A", "RAP1A", "ITGAM", "FMNL3", "RAB7A", "CCDC88A", "LILRB1", "DNAJB6", "TET2", "CAB39", "SP100", "MYCBP2", "SH3BGRL", "PRKCB", "ANKS1A", "KDM2A", "HSPD1", "HSP90AB1", "PTGES3", "STK10", "CTNND1", "GAB3", "MS4A4A", "LRRK2", "PIP4K2A", "BCAT1", "AAK1", "RUBCNL", "VPS13C", "RB1", "STK17B", "SEC14L1", "FNBP1", "ZFHX3", "ATF6", "INPP4A", "RASA2", "CPVL", "HERPUD1", "SH3TC1", "SPI1", "PSME4", "HERC4", "DLEU2", "BID", "HLA-DPB1", "TMEM131L", "ROCK1", "ARMH3", "STK4", "GNB4", "GPR183", "PIK3R1", "SMCHD1", "DNAJC5", "MAP3K2", "C20orf194", "SMURF2", "ENTPD1", "ARRB2", "ANKRD12", "UBE2D3", "CYFIP1", "LINC02798", "CD74", "ST6GAL1", "FAM102B", "PCED1B", "MAN2A1", "CEMIP2", "BACH1", "SETX", "TTN", "USP4", "ATG7", "KCNQ1", "MIS18BP1", "GPCPD1", "BAZ1A", "PIAS1", "C1QA", "LCOR", "NCK2", "DNAJB1", "BAZ2B", "MAN1A1", "LRRK1", "CARD11", "SH3RF3", "TSPAN14", "CTTNBP2NL", "ITSN2", "SLC16A3", "HAVCR2", "PRRC2C", "LITAF", "WSB1", "MSN", "PARP14", "LHFPL2", "CTNNB1", "HNRNPU", "SLC25A37", "MTSS1", "RUNX2", "MARCH1", "PTBP3", "SPRED1", "ARHGAP6", "PECAM1", "TRA2B", "CTSZ", "SRRM2", "RAC1", "GSK3B", "ATP1B3", "LINC02649", "ZNF331", "PRDM1", "NEDD4L", "MKNK1", "LPCAT2", "RCSD1", "PSAP", "CMTM7", "MORC3", "ARPC2", "TGFBI", "SUMF1", "ACSL4", "IL17RA", "ADAM9", "IQGAP1", "FRYL", "PAPOLG", "RREB1", "SPTLC2", "PLXND1", "PLXNC1", "PLCG2", "MAP2K1", "CLK1", "CBL", "GTDC1", "SLC43A2", "EHBP1L1", "ATRN", "C3", "CPM", "MTHFD1L", "HSPB1", "RESF1", "ADCY7", "SKI", "IL1RAP", "APAF1", "ELF1", "RPS6KA2", "TRAF3", "CNOT6L", "SWAP70", "ABHD12", "RILPL2", "PACSIN2", "WIPF1", "NLRP1", "SRGAP2C", "AP3B1", "RAD51B", "UBC", "PABPC1", "ELF2", "SNAP23", "CDC42SE2", "RAB20", "ATP11A", "RCOR1", "FAM168A", "JDP2", "SOAT1", "HSPA1A", "IFNAR2", "MXD1", "PABPC4", "RGS2", "GK", "FNIP2", "LGMN", "CAPZA1", "B3GNT5", "SERPINB9", "AL035446.2", "OSBPL8", "MIR181A1HG", "HEXA", "STX7", "HIST1H2AC", "AP1B1", "RIPK2", "DENND5A", "TTYH3", "SKIL", "BAG3", "EAF2", "ZNF124", "TGFB1", "TET3", "ATP6V0A1", "ELOVL5", "CRLF3", "TMEM51", "PLEKHA2", "NIBAN1", "KIF16B", "USP9X", "P2RX7", "TANK", "AGTPBP1", "PIK3CD", "ERP44", "AC068587.4", "USPL1", "TMEM163", "PHIP", "JAK2", "CARD8", "MGAT1", "ABI1", "MALT1", "SLC2A9", "ABHD3", "RAPH1", "NFKB1", "INSIG1", "IL13RA1", "PLCL2", "ACER3", "PPP1CB", "E2F3", "HIVEP3", "ATP6V1B2", "PRKCH", "CDK6", "PILRA", "DGKD", "SDCBP", "RASSF2", "PARP4", "CD58", "CERS6", "RAP1GAP2", "HSPA1B", "LINC01374", "NCOA3", "SNCA", "UBXN2B", "PDE4A", "PRIM2", "ARHGAP21", "B2M", "SEPTIN6", "ETS2", "MAT2A", "DOP1B", "CDK19", "TPK1", "RPS6KA3", "EIF4A1", "TXNIP", "UBASH3B", "TMCC3", "ANP32A", "IFNGR1", "GSAP", "OGFRL1", "HEATR3", "PLSCR1", "PDE4B", "CNTRL", "KDM7A", "ERBIN", "GNG7", "HEATR5B", "MAP4K4", "GRAMD4", "DAGLB", "ZNF804A", "SNX30", "CHORDC1", "TBC1D4", "NRP1", "IL4R", "STIP1", "ELL", "PBX3", "OSTF1", "SNX8", "OXR1", "PRPF4B", "MBP", "SGPP1", "B4GALT1", "MAN2B1", "FTH1", "CLEC2D", "PSMA1", "PARVB", "CHD7", "PPP3CA", "PSTPIP2", "ARRB1", "HMOX1", "COTL1", "LAMP1", "PRPF38B", "RASA3", "XKR6", "IVNS1ABP", "NPL", "AXL", "TMSB4X", "EPS8", "GPD2", "AMPD3", "MARS", "PTPRE", "ITGB1", "ARAP2", "PLTP", "CD44", "APBA1", "NSMAF", "AUTS2", "KLHL5", "DRAM1", "KLHL2", "CPED1", "ARHGEF6", "ARHGAP31", "SIPA1L2", "ZFP36L2", "ATP2B1", "GRN", "EDA", "STAT6", "USP53", "ABCA1", "PPM1L", "PLA2G4A", "DDX21", "ARHGEF2", "FOSL2", "ID2", "CDK14", "ZFP36", "MARCH3", "CADM1", "MCTP1", "HLA-B", "TIAM1", "FAM20A", "GRK5", "TPM4", "NR4A2", "JUNB")

Whyss_Correy_Ependymal <- c("DPP10", "CTNNA2", "NPAS3", "FMN2", "CFAP299", "ADCY2", "SORBS1", "CASC15", "LINC00609", "ERBB4", "NEBL", "DTHD1", "CRB1", "RFX4", "DNAH9", "TTC6", "PTPRZ1", "ADGRL3", "DTNA", "CFAP47", "AC012405.1", "AGBL1", "AC020718.1", "ANKFN1", "GRIN2A", "LINC01088", "NRCAM", "PAMR1", "TNIK", "KCNMB2", "PLEKHA5", "ADGB", "TRDN", "TMEM232", "CSMD3", "LMNTD1", "MIR4300HG", "C8orf34", "LINC02055", "PPP2R2B", "ST8SIA1", "AL357093.2", "DPYSL3", "NAV2", "RPH3A", "HYDIN", "TTC34", "MAP2", "TSHZ2", "CFAP54", "MIR99AHG", "AC019330.1", "DNAH3", "RSPO2", "KLHL32", "KIAA2012", "TTC29", "NCAM2", "GRIA1", "KCNMA1", "FYN", "GFAP", "GPM6B", "ASTN2", "PALLD", "CFAP46", "ANK2", "AC091078.1", "NTRK2", "DNAH10", "AC110023.1", "CHST9", "CCDC146", "GRAMD2A", "FUT9", "HPSE2", "PHACTR1", "CFAP44", "DNAH12", "ANXA1", "MUSK", "LINC00907", "WLS", "CFAP52", "GABRG1", "SLC1A3", "AC013470.2", "KIAA0825", "SOX6", "BBOX1", "NEK10", "ULK4", "PPM1E", "AL022068.1", "CDH2", "CADM1", "PACRG", "SLC1A2", "FAM184A", "ERICH3", "FRMD5", "PPFIA2", "ZBBX", "DST", "FAT3", "NEK5", "SERPINI2", "ID4", "PLCB1", "SPAG16", "C8orf37-AS1", "PCDH9", "SGIP1", "DNAAF1", "KIAA1549L", "JAKMIP2", "RIPOR2", "PTCHD1-AS", "LRRIQ1", "LAMA1", "MYLK3", "TIAM1", "NRG3", "WDR49", "DYNC2H1", "PREX2", "CCDC200", "PALM2-AKAP2", "DCDC1", "KIF6", "PZP", "DCLK2", "NR2F1-AS1", "KCNN3", "DNAH7", "PANTR1", "LINC00511", "TTYH1", "AL157886.1", "SLC7A11", "LMO2", "AC093689.1", "RASSF9", "DCLK1", "BASP1-AS1", "PSD3", "CDH20", "DOC2A", "DNAH2", "SOX5", "FRMPD2", "CLVS1", "COL8A1", "LAMA2", "AC114971.1", "QKI", "NKAIN3", "C1orf87", "AC004949.1", "GYG2", "CFAP43", "VWA3B", "TMEM178B", "MAP3K19", "UBE3D", "LINC00461", "LINC02649", "GABRB1", "FAM135B", "RSPH1", "ASB3", "VWA3A", "SPAG17", "ARAP2", "RNF43", "GMPR", "TMEM163", "EFHB", "FBXO15", "DNAH6", "MMP16", "ZBTB7C", "FABP6", "AC073941.1", "CHL1", "NEGR1", "LMO3", "STK32A", "CFAP70", "CFAP206", "ZNF709", "ANGPT1", "NRG4", "EMP1", "LIFR", "TOGARAM2", "IQGAP2", "KIAA0319", "CTNND2", "PLEKHG1", "NFASC", "BMPR1B", "AL133304.3", "ITGA2", "STON2", "FGF14", "CD109", "SCN1A", "MLIP", "MAP1B", "FHAD1", "PARD3", "RFTN2", "SLC14A1", "MEIS2", "LRRC9", "NELFA", "PHYHIPL", "NEK11", "GPC6", "CCDC30", "SGCD", "CFAP100", "CCDC88C", "BICD1", "SDK2", "CTXND1", "LMCD1-AS1", "AC097518.2", "SPATA17", "AL109809.5", "SPATA6", "TSGA10", "ARMC3", "CAPS", "ADGRL2", "LRP2BP", "GALNT8", "AC073050.1", "SYNPO2", "MKX", "GJA1", "GASK1B", "GABRG3", "ANTXR2", "ARMC4", "CMIP", "ODF2L", "STRBP", "AUTS2", "SLC9C2", "BMP2K", "CFAP221", "PPP1R42", "CNTN1", "KIF5C", "CFAP97D2", "RFX3", "AGBL4", "ADGRV1", "AK8", "MAPK10", "RMST", "CEP112", "PPP1R12B", "SMOC1", "NR2F1", "CASC2", "OPHN1", "AC002460.2", "OSBPL6", "VCAN", "TRPS1", "CFAP73", "KIAA1211", "FRY", "UNC79", "AC092957.1", "CCDC39", "EFCAB6", "GREB1L", "HKDC1", "LIMCH1", "PSD2", "AC097515.1", "CD36", "WWOX", "COL16A1", "AQP4", "HMGB1", "ARHGAP5", "AC005999.2", "P3H2", "CREB5", "GRIP1", "USP43", "UBL3", "SERPINE2", "CCDC170", "CD44", "GNG7", "HIF3A", "PRRT2", "ATP1A2", "C1orf198", "CFAP58", "DNAI2", "ARHGEF26-AS1", "IQCK", "CELSR1", "AL137139.2", "HSP90AA1", "KIF1B", "AC046134.2", "CFAP157", "TEKT1", "TAPT1-AS1", "SRGAP3", "AP002026.1", "IQCG", "GAREM1", "AC087564.1", "PRUNE2", "NRXN1", "AL390957.1", "IGSF11", "COL4A6", "WDR11-AS1", "AQP4-AS1", "WWTR1", "FIGN", "PCED1B", "LRGUK", "IFT88", "MGAT5", "SDK1", "PBX1", "BCL2", "ANKMY1", "IQCH", "MGST1", "PRKD1", "BAALC", "ZNRF3", "RAD51B", "SETBP1", "ZFHX4", "KNDC1", "CDH4", "CDC14A", "TNFRSF19", "ROPN1L", "EFCAB2", "MYO15B", "AL589740.1", "PCSK5", "DACH1", "LRRC74B", "ZMAT1", "DOCK7", "TMEM67", "TCTEX1D1", "TRDMT1", "ZFHX2", "CASC1", "PFKFB3", "AP001831.1", "CCN2", "DRC3", "MOK", "DGKB", "MYO16", "PDE4DIP", "MIR100HG", "C6orf118", "AC063979.2", "AC090578.1", "PTPRK", "PCGF5", "TRIM9", "NOVA1", "CDK14", "DYNLRB2", "ARMH1", "APC", "OLMALINC", "ARHGAP26", "MATN2", "ECT2L", "AGBL2", "CUBN", "MYCBP2", "RFX3-AS1", "MACF1", "CSPP1", "CSGALNACT1", "WDFY3", "INTU", "SDC2", "HULC", "NAV3", "DNAI1", "SRGAP1", "CABCOCO1", "RGS22", "DPYSL2", "WDR63", "FOXP2", "ARHGAP32", "RTN4", "NME9", "AKT3", "MAOB", "ARNT2", "MARCH3", "AZIN1-AS1", "MEGF9", "FREM1", "MAP4", "NFIB", "AC138627.1", "RPS6KA5", "BBOF1", "LDLRAD3", "WDR60", "GON4L", "ZFP36L1", "AK7", "SOX9", "GRAMD2B", "ERBIN", "CFAP61", "LTBP1", "OSBPL3", "SLC35F1", "LINC00271", "ZNF423", "LINC01138", "TTC23L", "FAM227A", "MPDZ", "CEP83", "CCDC81", "CROCC2", "BCL6", "DLEC1", "EML6", "SFXN5", "CMSS1", "C9orf135", "KCTD1", "USP2-AS1", "CCDC173", "CHD7", "CLU", "PLCH1", "TENM4", "CCDC18-AS1", "RNF19A", "CFAP69", "KREMEN1", "LRIG1", "SPAG6", "FAM27C", "CDKL3", "LINC00894", "SCD5", "TSC22D2", "LEKR1", "KCNN2", "AC108519.1", "ITPR2", "LMNA", "MTSS2", "GOLIM4", "RUFY3", "DOK5", "MYO9A", "CELF2", "ANKDD1B", "FAM169A", "AL513166.1", "LINGO2", "FAM81A", "TCTEX1D2", "PRTFDC1", "MAMLD1", "CACNA1A", "HSPB8", "DZANK1", "NALCN", "ENAH", "CEP290", "SOX2-OT", "COL28A1", "TTC21A", "CNTRL", "RPGR", "IQCA1", "AC007389.1", "MAP4K4", "DAW1", "FAM189A2", "RGS6", "AC092691.1", "OSBPL11", "COBL", "KATNAL1", "FGD6", "ABCA1", "BACH2", "USP54", "LRRC8D", "HEATR5A", "GPC5", "PPP3CA", "CAMTA1", "ATP2B1", "MARCH10", "NEK1", "DNAJB1", "GAN", "AL645568.1", "MYLK", "C20orf194", "AC020911.2", "SGO1-AS1", "EFR3B", "TRPC6", "MFHAS1", "NCALD", "GRIK4", "ALMS1", "RPS6KA2", "CCDC171", "AC239809.3", "ITGB4", "ABR", "FOS", "EXPH5", "CCDC88A", "LINC01876", "PDE3A", "NTM", "AC068587.4", "CPAMD8")

Whyss_Correy_Neural <- c("KCNIP4", "KCND2", "STXBP5L", "ADGRB3", "NRXN1", "RBFOX1", "CDH18", "GRIK2", "FGF14", "FSTL5", "RIMS2", "RALYL", "ZNF385D", "CADM2", "SYT1", "RIMS1", "RBFOX3", "UNC13C", "CA10", "ANKS1B", "GRIA4", "GRID2", "LINC00632", "CACNA1A", "CNTN1", "CAMK4", "NEGR1", "ARPP21", "CALN1", "NKAIN2", "EPHA6", "RYR2", "TIAM1", "RGS7", "SNAP25", "MEG3", "CNTNAP2", "ATP8A2", "ANK3", "LRRC4C", "FAM155A", "ZFPM2", "DLGAP1", "GALNT13", "NTM", "FGF12", "OPCML", "CNKSR2", "SGCZ", "SH3GL2", "CADM1", "PCLO", "CADPS2", "SNTG1", "MYT1L", "NRG1", "GABRB2", "SCN2A", "RELN", "RIT2", "ADAM22", "TENM1", "TRIM9", "ZNF385B", "SRRM3", "GRM4", "NRXN3", "SNHG14", "SPHKAP", "FRY", "GRIA2", "CACNA1B", "SRRM4", "CADPS", "SYN2", "SV2B", "UNC80", "GRM1", "MSRA", "DLG2", "CNTN4", "CELF4", "PDE1A", "MAP2", "SYNPR", "GABBR2", "ELAVL3", "KCNH1", "OLFM3", "GPR158", "MAPT", "SLIT3", "PTPRD", "SLC35F3", "PDE3B", "GRIN2A", "UNC79", "CDH22", "EPHB1", "MCTP1", "RUNX1T1", "FUT9", "AC073225.1", "ZNF804A", "KSR2", "KCNH7", "KIF5C", "TLL1", "VWC2", "BRINP3", "PPFIA2", "PDE10A", "JAKMIP2", "PTPRR", "FAT2", "SEZ6L", "NRCAM", "LINC00599", "CCSER1", "SLC35F4", "ADGRL3", "CACNG2", "PPP3CA", "SPOCK1", "MTCL1", "GABRB3", "MAP1B", "SCN8A", "CLVS2", "DAB1", "PLCXD3", "MACROD2", "AMPH", "PHACTR3", "CSMD3", "NMNAT2", "TRHDE", "CSMD1", "CDH10", "LINC02798", "LRRTM3", "ADD2", "NCAM1", "STMN2", "KCND3", "SLC8A2", "MMP16", "ERC2", "SPTBN4", "RASGRF1", "AC007614.1", "KIRREL3", "ASTN1", "RAB37", "CNTN6", "OXR1", "SCN3A", "ANK1", "VSNL1", "CTNNA2", "GRIN1", "NDST3", "PLCB1", "OTUD7A", "RUNDC3B", "HCN1", "TMEM178A", "IL1RAPL1", "IGSF21", "NAV3", "JPH4", "TMEM266", "NRXN2", "ETV1", "NLGN4X", "PAK5", "XKR6", "NDRG4", "SYT14", "SYN3", "CNTNAP4", "ERBB4", "BRINP1", "CAMK2B", "ADCY1", "GPHN", "ANK2", "HECW2", "KCNJ6", "TMEM63C", "ST8SIA5", "SCN1A", "ASTN2", "PCSK2", "RASGEF1C", "TMEM163", "STXBP5", "DGCR9", "LINC01619", "LIMA1", "ABLIM3", "KCNJ3", "RAP1GAP2", "TMEFF2", "MEG8", "HECW1", "UBASH3B", "TMEM178B", "NOL4", "KCNC1", "MDGA1", "NTRK3", "LHFPL4", "MEGF11", "CHD7", "PRR16", "TMEM130", "KCNQ2", "TNR", "NEURL1", "MCF2L2", "PPP2R2C", "AC093765.3", "XKR4", "PCDH9", "LRRC7", "DNM1", "ADAMTS16", "ASIC2", "NFASC", "PTCHD1-AS", "CACNA1D", "GABRG2", "GPR137C", "FRMD5", "CLVS1", "ICA1", "PHACTR1", "GRM5", "PRUNE2", "AC007563.2", "GAP43", "GUCY1A2", "CDK14", "GABRA1", "BICD1", "SGIP1", "EPB41", "CADM3", "HTR1E", "DPF3", "LINC00461", "NOVA1", "GRIN2B", "CPLX2", "TCF4", "WSCD2", "AUTS2", "SYP", "LRRTM4", "RAPGEF4", "AC073050.1", "TMEM132B", "PTCHD1", "FMN2", "JMJD1C", "BSN", "SLITRK5", "ST8SIA1", "SYNDIG1", "ERC1", "PRAG1", "LINC01798", "COL19A1", "LRFN5", "GNAO1", "CDH7", "ADAMTS18", "FAM135B", "PKIA", "FSTL4", "RAPGEF5", "SNAP91", "PCBP3", "LRRC3B", "KIF26B", "CELF5", "SLC4A4", "CNTNAP5", "GFOD1", "ST18", "CCDC144A", "CLASP2", "LRP1B", "GRM7", "GLRA2", "ANKRD18B", "GPRIN3", "UNC13A", "PDE4DIP", "LUZP2", "LHX1-DT", "ANKRD36", "PTPRN2", "APBA2", "PRKCB", "EBF1", "PPFIA4", "DPP10", "SLC7A14", "JPH3", "SYT2", "SETBP1", "CAMK2N1", "AKAP12", "ATP2B1", "CAMTA1", "PPM1E", "DSCAM", "ANKRD36C", "HIVEP3", "CSMD2", "BMPER", "CACNA1E", "EPB41L3", "DIRAS2", "ABLIM1", "DGKB", "MDGA2", "DNM3", "CNTN5", "APBA1", "EPHA7", "PTPRO", "AC092691.1", "AC005906.2", "ARFGEF3", "NBEA", "MEIS2", "IQSEC3", "IQCJ-SCHIP1", "RCAN2", "NSG2", "LRFN2", "GALNT9", "PPM1L", "C2CD5", "DCC", "AFF3", "RNF144A", "NLGN1", "SORL1", "SLC8A1", "SRGAP3", "AC016821.1", "DNAJC6", "AC013287.1", "DOCK3", "NALCN", "GRIA3", "DOK6", "PAUPAR", "LINC01151", "SEL1L3", "LRCH1", "ENOX1", "AC006148.1", "SIDT1", "MYT1", "MFSD6", "KCNQ5", "GABRB1", "GALNT7", "SEMA6A", "CAMK2D", "NEBL", "GRIA1", "MMP24", "MTMR7", "CAMKK2", "SORCS1", "KCNB1", "FUT8", "SUSD4", "ZNF704", "NEDD4L", "CACNA2D3", "FAM106A", "TSPAN5", "NRG2", "KIF1B", "CERS6", "KCNMB4", "CCNJL", "ROBO2", "SOGA1", "LDLRAD3", "TSPAN18", "PKIB", "FAM153CP", "KLHL3", "GLCE", "KIF3C", "KIAA1211", "SLC24A2", "FAT3", "SVEP1", "UNC5C", "DOP1B", "PTPN4", "SYBU", "NR2F1-AS1", "ST7", "PRKACB", "CACNA2D1", "PPP2R2B", "LSAMP", "AC105411.1", "TTC7B", "RTN1", "GDAP1", "CHL1", "AC090578.1", "GAREM1", "MAML3", "SRCIN1", "KCTD16", "AC068587.4", "CDYL2", "AKT3", "CDON", "AP001347.1", "DLGAP2", "SUPT3H", "PPM1H", "ELP4", "FYN", "PGBD5", "PHYHIP", "GRIP1", "GLCCI1", "ADAM23", "TRANK1", "LIN7A", "TMEM132D", "AKAP6", "WASF1", "MEIS1", "CDC42BPA", "DIPK2A", "BEND6", "SRGAP2", "COL27A1", "SHANK3", "RBFOX2", "IL16", "CCDC88A", "ARHGAP21", "GRID1", "EML5", "SLC44A5", "PAXIP1", "FRRS1L", "MAP7", "AC011246.1", "DCLK1", "ROBO1", "AC025159.1", "NLGN4Y", "LIMCH1", "RNF130", "LINC00342", "PLXNA4", "GARNL3", "AL513166.1", "ROCK2", "FRMD4A", "PLCL2", "NCAM2", "KCTD8", "HS6ST3", "CDH2", "MYCBP2", "PLA2R1", "CLIP1", "ATP9A", "RNF150", "ZBTB18", "CLSTN3", "TAFA1", "KCNMA1", "PSMA1", "DGKH", "CCDC200", "PER3", "NIN", "NIPAL3", "SLC35F1", "TMTC2", "POU2F1", "INPP4A", "LINGO2", "MAP3K5", "ADGRL1", "CD47", "KCNQ3", "AGTPBP1", "EVL", "PHF14", "DGKG", "C8orf34", "ADGRL2", "LINC01266", "DGKD", "GNG7", "HIVEP2", "AC245060.5", "SHISA9", "PDE4B", "ZNF804B", "AC008014.1", "ADAM12", "PDE7B", "POU6F2", "PHLDB2")

Whyss_Correy_Glial <- c("PCDH9", "IL1RAPL1", "CTNNA3", "ST18", "PLP1", "SLC24A2", "RNF220", "PIP4K2A", "MAN2A1", "MOBP", "PDE4B", "QKI", "SLC44A1", "MBP", "C10orf90", "EDIL3", "NCAM2", "PEX5L", "PDE1C", "SPOCK3", "FRMD5", "PPP2R2B", "TF", "TMTC2", "PRUNE2", "ELMO1", "TMEM144", "DOCK10", "PTPRD", "SLC5A11", "SLAIN1", "GPM6B", "KIRREL3", "AC008571.2", "MARCH1", "ZEB2", "DLG2", "TTLL7", "CDH20", "NCAM1", "ATP8A1", "NTM", "ARAP2", "ANLN", "OTUD7A", "LIMCH1", "PHLPP1", "SH3GL3", "NCKAP5", "SPOCK1", "CNTN2", "UGT8", "BCAS1", "NKAIN2", "FRYL", "AGAP1", "GRID1", "PTPRK", "DNAJC6", "PLEKHH1", "CCDC88A", "SLC7A14-AS1", "SHTN1", "ANO4", "PALM2-AKAP2", "MAP7", "TMEM165", "FRMD4B", "ERBIN", "FOLH1", "RAPGEF5", "PXK", "ANK3", "LRP2", "LINC01608", "LINC00609", "SHROOM4", "PLCL1", "UNC5C", "DAAM2", "QDPR", "AMER2", "MOG", "DNM3", "NFASC", "GRM3", "AK5", "CD22", "KCNH8", "SEPTIN7", "MAP4K4", "CNDP1", "POLR2F", "SH3TC2", "MYRF", "SRCIN1", "DPYD", "DST", "DPYSL5", "TMEFF2", "RTN4", "KLHL32", "LINC00320", "CADM2", "SOX2-OT", "ENOX1", "SCD5", "PCSK6", "KCNMB4", "AC009063.2", "HDAC2-AS2", "SLCO1A2", "PHACTR3", "CDK19", "AC090015.1", "AC108721.2", "TRIM2", "HHIP", "SILC1", "CREB5", "CLASP2", "AC012494.1", "ANK2", "GNAO1", "MEGF10", "RFTN2", "NPAS3", "GLDN", "DIP2B", "AKAP6", "CLDN11", "GPR37", "ADGRB3", "ARHGAP21", "DSCAML1", "MAPT", "LINC00639", "CLMN", "DLC1", "FUT8", "AATK", "GPRC5B", "CNTNAP4", "COLGALT2", "MAGI2", "ANKS1B", "SCD", "TJP1", "MAP4K5", "GAB1", "ELAVL3", "SIK3", "MAG", "AGTPBP1", "AC016597.1", "TCF12", "CERCAM", "CLDND1", "MAL", "IGSF11", "LAMP2", "SORT1", "C5orf64", "AL033523.1", "DLG1", "ABCA2", "LPAR1", "HECW2", "EPCAM-DT", "LINC01630", "PRR5L", "ADARB2", "RAB40B", "ERBB3", "SEMA6A", "ALCAM", "LINC01170", "NRXN3", "DOCK5", "KLHL4", "DPYSL2", "FA2H", "PIEZO2", "CPB2-AS1", "CADM1", "PREX1", "CARNS1", "CDH19", "P2RX7", "RNF130", "HAPLN2", "KIF5C", "HHATL", "ZDHHC20", "AUTS2", "SLC39A11", "CNP", "MIR181A1HG", "ASPA", "LDB3", "HIP1", "ADIPOR2", "ATG4C", "APBB2", "DIP2C", "KIF6", "PPP1R16B", "AL359091.1", "NALCN", "PTBP2", "SPP1", "AC079352.1", "HIPK2", "FGF1", "AL392086.3", "AL133304.3", "SYNJ2", "USP54", "LRRC63", "ANKIB1", "LURAP1L-AS1", "LANCL1", "PHYHIPL", "GNAI1", "GPHN", "PDE8A", "ANKRD18A", "DOCK1", "AC006148.1", "DOCK9", "ARHGAP23", "NEO1", "SH3TC2-DT", "AC004690.2", "PLD1", "TPPP", "ZNF536", "HS3ST5", "SEMA4D", "LGR5", "VRK2", "ATP10B", "KIF1B", "ZSWIM6", "DOCK4", "TTYH2", "CNTNAP2", "SORCS2", "ARFGEF3", "PLPPR1", "COL4A5", "CCP110", "ZNF638", "FCHO1", "CHRM5", "LRRC7", "YPEL2", "ITGA2", "AC096564.1", "CDC42BPA", "ABCA8", "ERBB4", "TRIM59", "C12orf76", "SGIP1", "CACNA2D1", "SECISBP2L", "PLAAT3", "KIF13B", "ZNF704", "LRRC8D", "JAKMIP2", "B3GAT1", "ACER3", "TMEFF1", "PDE1A", "APLP1", "OSBP2", "JAKMIP3", "PRIMA1", "FEZ1", "COBL", "AL354809.1", "FAM222A", "PWRN1", "BIN1", "CDK18", "TMCC2", "HPN-AS1", "CPM", "OLMALINC", "FAM13C", "ROBO1", "PLLP", "PTMA", "RFFL", "SLC4A8", "HEPACAM", "GRIA2", "TP53TG5", "AC087564.1", "GAREM1", "TMTC4", "SYNDIG1", "NUTM2A-AS1", "EFHD1", "PLA2G4C", "SLC22A15", "LINC01505", "TCF4", "SLC6A1", "RYBP", "PTGDS", "ARHGAP22", "MBNL2", "SLCO3A1", "AGPAT4", "DICER1", "LIPA", "BBS2", "SSH2", "PICALM", "NINJ2", "RASGRF1", "ZNF365", "RNF13", "GALNT13", "STXBP6", "PRRG1", "PPP2R2A", "LIPE", "DYSF", "SEC14L5", "NXPE3", "GFAP", "TFEB", "IQCJ-SCHIP1", "ZKSCAN1", "FMNL2", "DUBR", "CPQ", "BAZ2B", "KIAA1324L", "UBE2E2", "LRRC1", "HSP90AA1", "DEPTOR", "ZCCHC24", "PPP1R21", "EEPD1", "IQGAP1", "MYO1D", "GRIA4", "RASGRP3", "ADAMTS18", "SPATA6", "NDE1", "EPS15", "ZNF652", "SEPTIN8", "PRKACB", "KIAA0930", "TLE4", "SEPTIN4", "VPS35L", "TSPAN5", "PROX1", "SLC22A23", "CHD7", "GNG7", "FAM171A1", "RAP1GDS1", "MEGF9", "RNF144A", "FNBP1", "PADI2", "MIR100HG", "RICTOR", "CAMK2N1", "GTDC1", "BACH2", "GLTP", "LSS", "HTRA1", "FNIP1", "AC106791.1", "NAV2", "MEIS2", "CEP170", "FUT9", "FAM53B", "NCOA7", "CALM2", "PPM1H", "SNX30", "PACS2", "L3MBTL4", "SPTLC2", "SERPINI1", "MAPRE2", "MOB3B", "ATL1", "CDK6", "TMC7", "GREB1L", "NLK", "PHLDB1", "TECPR2", "SLC13A3", "ZFHX4", "TMEM63A", "GATM", "RDX", "RAB11FIP4", "SYT9", "AMD1", "RASSF2", "VPS13C", "GARNL3", "CDH8", "PPFIA2", "CAPN3", "MAP1B", "WIPF1", "KLHL3", "FP236383.3", "ASTN2", "ZNF708", "AC007566.1", "DLEU2", "NAV1", "SNX29", "AL589740.1", "ANKRD36C", "KIAA1755", "RCAN2", "SCARB2", "TUBA1A", "AC016821.1", "HBS1L", "PANTR1", "STRN", "TMOD2", "EPHB2", "SH3PXD2A", "AC006059.1", "CORO2B", "PIK3C2B", "MOSPD2", "TRAPPC10", "SGCD", "CAMTA1", "RAB30", "SYT11", "EPB41L2", "ADCY5", "AF165147.1", "SELENOP", "MTURN", "FNIP2", "TBC1D12", "FIGN", "MFSD6", "SLF2", "RYR2", "FANCL", "C21orf62-AS1", "DENND5A", "GPSM2", "LRRC8B", "DCC", "AC011346.1", "CDK17", "MYO9B", "LIFR", "SLC9A9", "DTNB", "PHIP", "C1orf198", "CALN1", "CDKL1", "FAM168A", "LINC00877", "KNOP1", "CRYBG3", "ATP11A", "HMGB1", "WDFY2", "GAB2", "FOXN2", "SEMA3C", "NAPEPLD", "SUN2", "ABCA6", "TNIK", "DNAJB2", "CCDC144A", "GALNT15", "GSN", "NPC1", "PRTFDC1", "BOK", "PHF14", "FKBP5", "ABHD17B", "KLHL2", "SEMA3B", "MTMR7", "ZBTB37", "RHOBTB1", "ENTPD3-AS1", "TAPT1-AS1", "ELL2", "PCBP1-AS1", "FAM102A", "DYNC1LI2", "CSMD2", "TARSL2", "ZNF146", "DAPK2", "NEGR1", "MGAT5", "ENPP6", "TIAM1", "HULC", "ANKRD36", "AC090517.5", "PREX2", "NDRG1", "FYN", "NDRG2", "NDUFAF6", "DPF3", "LRRC4C", "RHOBTB3", "CRYAB", "TESK2", "RASSF4", "CD47", "NASP", "DIPK2A", "RTTN", "SESTD1", "SORL1", "DDHD1", "PTN", "RNASEH2B", "CENPP", "LDLRAD3", "PMP22", "INF2", "MYO18A", "AC007364.1", "FDFT1", "TRIM9", "AFMID", "SNCA", "MT-ND2", "CALD1", "IKZF2", "CCDC200", "GRID2", "SLC1A3", "CLIC4", "GLUL", "ADGRL3", "MDGA2", "TXNIP", "LAMA2")

Whyss_Correy_Endothelial <- c("ARL15", "MECOM", "VWF", "HSP90AA1", "ANO2", "LDB2", "FLT1", "SERPINE1", "DOCK9", "HSPH1", "HSPD1", "PECAM1", "MSN", "TCF4", "PELI1", "PTGES3", "PIK3R3", "P4HA1", "RBMS3", "PLS3", "CRIM1", "IL3RA", "HSPA5", "RALGAPA2", "CLIC4", "MOB4", "DNAJB6", "PTPRB", "PICALM", "PLPP1", "SERPINH1", "ATP8B1", "MYH9", "TFPI", "BMPR2", "ST8SIA6", "UBC", "TEK", "CMIP", "GPBP1", "ST6GALNAC3", "HSPB1", "GRB10", "MCTP1", "MLKL", "ARHGAP26", "HSP90AB1", "CDK17", "DKK2", "BMP6", "PLEKHG1", "CDH13", "UTRN", "ITGA5", "PODXL", "ELOVL7", "MRTFB", "FNIP2", "MBNL1", "DOCK4", "HSPG2", "SPAG9", "PTGS2", "ERG", "SETBP1", "GRK5", "EGFL7", "IRAK2", "NOTCH4", "FNDC3B", "SASH1", "ARHGAP29", "THSD7A", "AKT3", "VEGFC", "SEC14L1", "HSPE1", "ADAMTS9", "ADGRL4", "BAG3", "SLC2A3", "ELMO1", "PALM2-AKAP2", "IL4R", "TPM4", "CLIC2", "CTTNBP2NL", "COL8A1", "TANC1", "NEDD9", "ATP10A", "INSR", "ERO1A", "APBB2", "NFIB", "CCNY", "RAPGEF2", "SHANK3", "MAST4", "RUNX1T1", "PRKCH", "STOM", "A2M", "PDE4D", "ARAP3", "PIK3C2B", "LIMCH1", "ABCB1", "FLNB", "DNAJB1", "ECE1", "RAPGEF5", "UGCG", "ENG", "PRKY", "GFOD1", "NKAIN2", "KMT2E", "PKHD1L1", "ARID4B", "SWAP70", "RND1", "QKI", "ASAP1", "NUDT4", "IGFBP3", "CAPZA2", "ADAMTS1", "LRMDA", "STIP1", "KIAA0355", "FCHSD2", "HNRNPU", "FLI1", "ZEB1", "STC1", "RPS6KA3", "PARVB", "HSPA8", "DNAJA1", "SORBS1", "ETV6", "DOCK6", "ARID5A", "EFNB2", "MALAT1", "GRASP", "PCSK5", "IQSEC1", "HRH1", "FRYL", "CALCRL", "AC010737.1", "FGD5", "CALD1", "ADIPOR2", "TCF12", "CHORDC1", "PTPRM", "JUN", "PDLIM5", "PTK2", "ITSN2", "TMTC1", "SAT1", "EMCN", "DUSP5", "KAT6A", "PREX2", "ETS2", "SNAP23", "TRA2B", "WSB1", "CLK1", "TSHZ2", "CAMK2D", "PRRC2C", "ADD1", "COL4A1", "RAPGEF4", "ELOVL5", "ICAM1", "SVIL", "PTPN12", "LRIG1", "PRSS23", "KYNU", "ARHGAP31", "CPNE8", "ELF2", "AC008050.1", "JMJD1C", "VMP1", "MEF2A", "MACF1", "PLCB1", "ABR", "HSP90B1", "HIF3A", "MTCH1", "CBLB", "TSC22D2", "PDLIM1", "PTGIS", "HMCN1", "GNB4", "CDC42BPB", "GALNT18", "PLCG2", "ATP11A", "HIF1A", "ELL2", "FBLIM1", "B2M", "ZFAND3", "DIPK2B", "PIK3C2A", "FYN", "ATP13A3", "AHCTF1", "FBXO11", "KIF13B", "ATF3", "SRGN", "MED13L", "AAK1", "RAP1B", "ABI3BP", "GNAQ", "EXT1", "NECTIN2", "IRAK3", "ARHGEF7", "SMAD6", "GPCPD1", "CFH", "CTNND1", "GNB1", "XAF1", "TRIO", "PTMA", "DGKH", "DENND5A", "RBMS2", "AC007681.1", "MRTFA", "TCF7L1", "STK38L", "RNF144B", "CFLAR", "PLA2G4A", "ARHGEF2", "COL25A1", "BAZ2B", "TDRD10", "SLC20A1", "ZNF423", "APP", "SESTD1", "NUAK1", "DOCK1", "NDRG1", "WAC", "NUCB1", "LIMS2", "OPHN1", "BACE2", "HSPA1A", "TLL1", "ZFAND2A", "MAP4", "WWTR1", "NAMPT", "PCDH11X", "NOVA2", "LINC01578", "HERC2", "HLA-E", "HIP1", "TOP1", "HNRNPA2B1", "GOLPH3", "HLA-B", "ADGRF5", "UBE2D3", "NEDD4", "AC007952.4", "PPP1CB", "TGM2", "PLSCR4", "PCDH11Y", "MEF2C", "AC020916.1", "ESYT2", "CNOT6L", "SLCO4A1", "TSPAN14", "FAM155A", "HSPA1B", "RIMKLB", "ANKRD28", "SAMD12", "STAT3", "TNFRSF10B", "ACTN1", "RBPMS", "COL4A2", "TMEM140", "MYO1E", "RNF145", "RBFOX2", "TJP2", "ZNF83", "CHRM3", "SH3RF3", "CLEC1A", "LAMA3", "C22orf34", "ATP2C1", "UBR2", "KMT2C", "ITPR2", "PPP2R2A", "NOP58", "NCOA3", "FOSL2", "CYTH1", "ICA1", "CCNL1", "CCN2", "RICTOR", "DNAJA4", "EXOC6", "GABARAPL1", "MGP", "TANC2", "ARHGAP23", "ADGRL2", "LINC02147", "CNTNAP3B", "NPDC1", "ARHGAP17", "NR4A1", "NEAT1", "ARHGEF12", "SHROOM4", "EWSR1", "CTNNB1", "GIPC2", "GASK1B", "FXR1", "HSPA4", "ABL2", "NFKB1", "EML4", "USPL1", "RHOJ", "ITGB1", "ACTR3", "PPP3CA", "ATP1A1", "STK3", "GBP2", "ITGAV", "PSMA6", "HMBOX1", "AZIN1", "HECW2", "GMDS", "AC119674.1", "DNAJC1", "SOX5", "EMP1", "FERMT2", "STAG1", "EPAS1", "ZSWIM6", "TNFRSF10D", "IL1R1", "LIFR", "SQSTM1", "MSRB3", "PEAK1", "NLRP1", "KPNA1", "TBCD", "ARPC2", "FMNL3", "VCAM1", "TIMP3", "MAP4K4", "CASC3", "SLC45A4", "ARL5B", "FOS", "IFI16", "STX12", "AHR", "TBK1", "HERPUD1", "YES1", "CHIC2", "AFDN", "TOX", "PDE3A", "CADPS2", "CNOT2", "PRKD1", "FAT4", "HBP1", "COL5A2", "LDLRAD3", "ARGLU1", "ADAM17", "RPGR", "NRIP1", "SESN3", "ST13", "PPP1R12A", "B4GALT1", "CWF19L2", "RIT1", "ZFPM2", "SULF2", "ACTB", "PPP1R15A", "FKBP1A", "IL15", "MTHFD1L", "PIP5K1C", "IQCJ-SCHIP1", "SUSD6", "SMURF2", "PDE10A", "KCNT2", "IPO7", "ERBIN", "UBE2B", "PGM5", "SNRK", "SOS1", "NLK", "SRSF7", "PPFIBP1", "CX3CL1", "TAF1", "NCKAP1", "NFKBIA", "NPIPB5", "HSPA4L", "ITGA1", "IFNAR2", "CACYBP", "CCN1", "TLE1", "KLF6", "TMC7", "HIVEP2", "VPS13C", "NAV1", "CMTM8", "YTHDC1", "TBC1D4", "CALR", "ATL2", "TGFBR3", "CEMIP2", "CARMIL1", "CASP4", "MOSMO", "RDX", "BCAS2", "UBA6", "CHD1", "RASA2", "PITPNB", "DENND3", "FRMD4B", "LDLR", "ABLIM3", "PAPSS2", "ID2", "CCDC88A", "PPP1R3B", "GAREM1", "FAM241A", "DYSF", "WARS", "KLHL2", "RHEB", "KDM7A", "YPEL2", "PER1", "MCC", "ACSL3", "ELN", "FUT8", "SMAD1", "DISC1", "TM4SF1", "PRKD2", "LDHA", "PIP4K2A", "TPM3", "MIR222HG", "INSIG1", "RPH3AL", "SFT2D2", "AHSA1", "MGAT5", "BHLHE40", "TSHZ1", "TSPYL2", "PHF14", "ZC3HAV1", "GALNT15", "LIMA1", "AC004889.1", "RELB", "KLF7", "ARHGEF28", "ZNF654", "DDX3X", "MAN1A1", "ADAM9", "RGS3", "PTBP3", "ROCK2", "FN1", "GCNT1", "TMEM87B", "AKAP12", "PREX1", "PLPP3", "JUNB", "CPAMD8", "MET", "CDKN1A", "TMTC2", "B4GALT5", "BICD1", "CHN1", "TXNIP", "BACH2", "SLCO2A1", "DDX58", "EIF4A1", "ALPL", "CORO1C", "PBX3", "CCDC85A", "AUTS2", "PARP14", "ITGA9", "PLA2G4C", "PLCXD3", "SRGAP1", "RANGAP1", "TBC1D8", "SLC4A7", "CDK14", "ZFP36", "USP53", "EPS8", "ANXA1", "SMAD9", "PFKFB3", "EBF1", "ABCA1", "RNF220")


known_dark <- c("NOS3", "AMPK", "ENDOG", "TFAM", "SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "PPARGC1A", "GABPA", "NRF1", "NRF2", "TFAM", "PPARG", "ADH1B", "ALDH2", "ENO1", "EPM2A", "GAPDH", "LDHA", "MPC1", "PDHA2", "PFKFB4", "PHKG2", "COX1", "MT-ND4", "MT-ND1", "RARRES2", "IGFBP7", "COX3", "STRIP2", "CXCL14", "AQP1", "CLDN5", "SLC26A3", "PRLR", "GPX3", "COX3", "MT-ND4", "COX2", "MT-ND3", "MT-ND1", "MT-ATP6", "COX1", "MT-ND2", "MT-CYB", "MT-ND5", "TTR", "MT-ND4L", "CTSD", "MT-ATP8", "CLDN5")

ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)