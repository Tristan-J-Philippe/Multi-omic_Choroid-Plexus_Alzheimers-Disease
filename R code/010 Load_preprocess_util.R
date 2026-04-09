###
### 010 Load_preprocess_util.R
###
# Purpose: Load, preprocess, and QC hashed multi omic single nuclei RNA sequencing functions


# Dependencies:
#function to make python environment modified from sctk to take care of numpy-numba dependency issue.
library(reticulate)
sctkPythonInstallVirtualEnv <- function (envname = "sctk-reticulate", packages = c("scrublet"), selectEnvironment = TRUE, python = NULL) #scipy", "numpy", "astroid", "six", "scrublet", "scanpy", "scanorama", "bbknn", "anndata", "Cython", "pytz", "llvmlite", "importlib-resources", "Pillow", "sortedcontainers", "fbpca", "annoy", "zipp", "threadpoolctl", "six", "pyparsing", "pillow", "packaging", "kiwisolver", "joblib", "intervaltree", "fonttools", "cycler", "python-dateutil", "contourpy", "scikit-learn", "matplotlib", "geosketch", "pandas", "umap_learn", "scikit-image", "python_dateutil", "tzdata", "networkx", "imageio", "seaborn"
{
  path <- reticulate::virtualenv_create(envname = envname)
  for (i in packages) {
    reticulate::virtualenv_install(envname = envname, packages = i, ignore_installed=F, force=F)
  }
reticulate::virtualenv_remove(envname = envname, packages="numpy==1.26.0", confirm = FALSE) #remove new numpy version
reticulate::virtualenv_install(envname = envname, packages = "numpy==1.23.0", ignore_installed = FALSE) #use slightly older numpy version for numba to work, necessary for scrublet
  if (isTRUE(selectEnvironment))
    sctkPythonInstallVirtualEnv(envname = envname)
  invisible(path)
}

sctkPythonInstallVirtualEnv() #make python environment !!!!comment out if already made!!!!

use_python(python="/home/tristan-j-philippe/.virtualenvs/sctk-reticulate/bin/python3.9", required=TRUE) #!!!!change path to python

reticulate::py_config() #Confirm correct virtual environment and dependencies.
reticulate::py_available()
reticulate::py_module_available('cellbender') #sanity check
reticulate::py_module_available('scrublet')

#Load required pacages
library(Matrix)
library(data.table)
library(plyr)
library(dplyr)
library(GenomicRanges)
library(singleCellTK)
library(Seurat)
library(Signac)
library(future)
library(ggplot2)
library(cowplot)
library(SoupX)
library(DoubletFinder)
library(qlcMatrix)
library(FastCAR)
library(gridExtra)
library(dsb)
library(org.Hs.eg.db)
library(GenomeInfoDb)
library(cluster)
library(fitdistrplus)


#For ATAC
library(EnsDb.Hsapiens.v86)
#Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC" #!!!!

#Set current working direct.
curr_dir <- getwd() #Get current directory; dependency for most functions to output in folders


##function to load and separate RNA, HT, and ATAC data
loadfun <- function(nameobj, GEX_data_dir, ATAC_data_dir=NULL, meta=NULL, exp_cmo=NULL, tfidfMin=0.9, forceAccept = FALSE, ContaminationFraction=NULL){
  #Description: Loads data matrices, filters, and qcs raw Cell ranger input matrices
  #Inputs
  ##nameobjs: Sample name, if multiple use c(), must be char
  ##GEX_data_dir: Directory to sample; up to and not including [nameobj]/outs)
  ##ATAC_data_dir: Directory to ATAC sample; up to and not including [nameobj]/outs)
  ##meta: Two column metadata table of metadata to add for this sample (e.g. batches)
  ##exp_cmo: Must specify CMO range since not all experiments have all 16, otherwise you get errors (limfit #100, cells with zero counts exist in clusters, etc.) in HTODemux. Also used to determine if CMOs are expected. Default is NULL (no CMOs expected)
  ##tfidfMin: SoupX variable default is 0.8, can modify if <100 or >500 genes per cluster. Usually more genes per cluster doesn't hurt, but <100 makes rho unstable
  ##forceAccept: SoupX variable default is FALSE but if sample is very dirty set to TRUE.
  ##ContaminationFraction: Manually set soupX contamination fraction
  #Returns: Seurat object with RNA, ATAC and CMO as assays with all QC information in metadata.
  capture.output({ #capture entire output
    
    path <- paste(curr_dir, "/Sample_preprocess/", sep="") #specify dir
    if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
    
    mtx.RNA <- Read10X_h5(paste(GEX_data_dir, nameobj, "_raw_feature_bc_matrix.h5", sep = ""))
    
    #create sce for RNA mtx !!!!You need to repeat sce/sctk steps with each assay separately if you want to use on CMO, antibodies, ATAC etc.
    if(is.null(exp_cmo)==FALSE){
      sce <- SingleCellExperiment(mtx.RNA$`Gene Expression`)
      mtx.RNA_GEX <- mtx.RNA$`Gene Expression`
    }else{
      sce <- SingleCellExperiment(mtx.RNA)
      mtx.RNA_GEX <- mtx.RNA
    }
    
    ##knee plot
    srt <- CreateSeuratObject(mtx.RNA_GEX)
    srt <- CalculateBarcodeInflections(srt, barcode.column="nFeature_RNA")
    plts <- BarcodeInflectionsPlot(srt)
    plot(plts)
    ggsave2(paste(path, nameobj, "/nfeature_kneeplot.png", sep=""), plot = plts, width=4, height=2.5)
    
    plts <- ggplot(srt@meta.data, aes(nFeature_RNA)) + geom_histogram()
    plot(plts)
    ggsave2(paste(path, nameobj, "/pre_cut_nFeat_RNA.png", sep=""), plot = plts, width=4, height=2.5)
    
    ##Remove empty droplets:
    #RunemptydropQCs
    names(assays(sce)) <- 'counts'
    sce <- runDropletQC(sce)
    plotEmptyDropsResults(inSCE = sce,  axisLabelSize = 20, sample = NULL, fdrCutoff = 0.01, dotSize = 0.5, defaultTheme = TRUE)
    plotBarcodeRankScatter(inSCE = sce, title = "Rank Plot", legendSize = 14)
    
    #check dims pre
    plotSCEViolinColData(sce, coldata = "sum", summary = "mean", title = "Pre-filter", ylab = "Total counts")
    
    #filter
    sce_cell <- subsetSCECols(sce, colData = c("!is.na(dropletUtils_emptyDrops_fdr)", "detected < 10000", "detected > 500")) #remove: NA, Total Features per cell <10K, Total Features per cell >500 !!!!

    #check dims post
    print("Post drop filter:")
    print(dim(sce_cell))
    plotSCEViolinColData(sce_cell, coldata = "sum", summary = "mean", title = "Post-filter", ylab = "Total counts")
    
    
    ##Get clusters to feed soupX (and DoubletFinder)
    srt <- NFSPUC_fun(sce_cell)
    
    
    ##soupX
    sc <- SoupChannel(sce@assays@data@listData[["counts"]], sce_cell@assays@data@listData[["counts"]]) #make sX object
    sc <- setClusters(sc, srt@meta.data[["seurat_clusters"]])
    sc <- autoEstCont(sc, tfidfMin=tfidfMin, forceAccept = forceAccept) #!!!! You want ~100 genes to pass tf-idf cutoffs so change the min appropriately. More won't hurt but less really hurts rho (aka contamination) estimation
    print("Rho, aka contamination:")
    print(sc[["metaData"]][["rho"]][1])
    if(is.null(ContaminationFraction)==FALSE){sc <- setContaminationFraction(sc, ContaminationFraction)} #Manually set contaminationFraction
    dirty_feat <- sc[["fit"]][["markersUsed"]][["gene"]]
    sc <- adjustCounts(sc, roundToInt=TRUE)
    
    
    ##doublet finder
    #pK
    min.pc <- min_pc_fun(srt)
    sweep <- paramSweep_v3(srt, PCs = 1:min.pc, sct = FALSE)
    sweep <- summarizeSweep(sweep, GT = FALSE)
    bcmvn <- find.pK(sweep)
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
    pK <- bcmvn.max$pK
    pK <- as.numeric(levels(pK))[pK]
    print("Doublet Finder pK:")
    print(pK)
    
    nExp_poi <- round(0.25 * nrow(srt@meta.data)) #Expected doublets !!!! Assuming 25% doublet formation rate - tailor for your dataset
    annotations <- as.character(srt@meta.data$seurat_clusters) #annotations
    
    srt <- doubletFinder_v3(srt, PCs = 1:min.pc, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE, annotations = annotations)
    
    #Just generalizing the name of the score, call, and annotation
    srt@meta.data[["DoubletFinder_call"]] <- srt@meta.data[,grepl("DF.classifications", names(srt@meta.data))]
    srt@meta.data[,grepl("DF.classifications", names(srt@meta.data))] <- NULL
    srt@meta.data[["DoubletFinder_score"]] <- srt@meta.data[,grepl("pANN", names(srt@meta.data))]
    srt@meta.data[,grepl("pANN", names(srt@meta.data))] <- NULL
    srt@meta.data[["DoubletFinder_DFclusters"]] <- srt@meta.data[,grepl("DF.doublet.contributors", names(srt@meta.data))]
    srt@meta.data[,grepl("DF.doublet.contributors", names(srt@meta.data))] <- NULL
    
    
    ##FastCAR
    ambProfile <- describe.ambient.RNA.sequence(mtx.RNA_GEX, start = 50, stop = 2000, by = 50, contaminationChanceCutoff = 0.05)
    plot.ambient.profile(ambProfile)
    ambientProfile <- determine.background.to.remove(mtx.RNA_GEX, emptyDropletCutoff=500, contaminationChanceCutoff=0.2) #Contamination varies a lot and this package doesn't have an objective way of determining contamination.
    sc_FC <- remove.background(mtx.RNA_GEX, ambientProfile)
    
    
    ##General QCs (doublet, RNA contamination, MT): runCellQC or addPerCellQC to assign to matrix? Unsure what the difference is.
    sce <- sce_cell
    sce <- runCellQC(sce, algorithms = c("QCMetrics", "scDblFinder", "cxds", "bcds", "cxds_bcds_hybrid"), sample = colData(sce)$sample, collectionName = "mito", mitoRef = "human", mitoIDType = "symbol", mitoGeneLocation = "rownames") #"scrublet", "decontX", 
    
    
    ##Merge metadata with seurat
    meta_srt <- data.frame(srt@meta.data)
    meta_sce <- data.frame(sce@colData)
    meta_comb <- cbind(meta_srt, meta_sce)
    srt@meta.data <- meta_comb
    srt@meta.data[["nCount_originalexp"]] <- NULL
    srt@meta.data[["nFeature_originalexp"]] <- NULL
    srt[["soupXcounts"]] <- CreateAssayObject(sc) #add soupXcounts]
    srt <- subset(srt, nFeatures_soupXcounts>200)

    ##Add antibody capture
    if(is.null(exp_cmo)==FALSE){
      mtx.AC <- mtx.RNA$'Antibody Capture'
      mtx.AC <- mtx.AC[exp_cmo,]
      
      mtx_cells <- mtx.AC[,colMaxs(mtx.AC,)>10]
      mtx_cells <-  mtx_cells[,mtx_cells@Dimnames[[2]] %in% colnames(srt)] #select cells based on srt
      srt <- srt[,which(colnames(srt) %in% mtx_cells@Dimnames[[2]])] #select cells that have hashtags.
      mtx_cells <-  mtx_cells[,mtx_cells@Dimnames[[2]] %in% colnames(srt)] #select cells based on srt
      mtx_cells@Dimnames[[2]] <-  mtx_cells@Dimnames[[2]][order(colnames(srt)), drop=FALSE] #Ensure cells are ordered the same way  
      srt[["CMO"]] <- CreateAssayObject(counts=mtx_cells)
      srt <- NucleosomeSignal(object = srt, assay="ATAC")
      srt <- TSSEnrichment(object = srt, assay="ATAC")
      
    }else{print('No HTO/CMO/AbC')}
    
    
    
    
    ##Add ATAC
    if(is.null(ATAC_data_dir)==TRUE){print('No ATAC added')}else{
      mtx.ATAC <- Read10X_h5(paste(ATAC_data_dir, nameobj, "_raw_feature_bc_matrix.h5", sep = ""))
      fragpath <- paste(ATAC_data_dir, nameobj, "_atac_fragments.tsv.gz", sep = "")
      mtx.ATAC <- mtx.ATAC$'Peaks'
      mtx.ATAC <-  mtx.ATAC[,mtx.ATAC@Dimnames[[2]] %in% colnames(srt)] #select cells based on srt
      srt <- srt[,which(colnames(srt) %in% mtx.ATAC@Dimnames[[2]])] #select cells that have ATAC 
      mtx.ATAC <-  mtx.ATAC[,mtx.ATAC@Dimnames[[2]] %in% colnames(srt)] #select cells based on srt
      mtx.ATAC@Dimnames[[2]] <-  mtx.ATAC@Dimnames[[2]][order(colnames(srt)), drop=FALSE] #Ensure cells are ordered the same way 
      srt[["ATAC"]] <- CreateChromatinAssay(counts = mtx.ATAC, sep = c(":", "-"), fragments = fragpath, annotation=annotation) 
    }
    
    
    ##MT and RB funs
    DefaultAssay(srt) <- "soupXcounts"
    srt <- mt_fun(srt, nameobj) #mitochondria function
    srt <- rb_fun(srt, nameobj) #ribosome function
    
    plts <- ggplot(srt@meta.data, aes(nFeature_RNA)) + geom_histogram()
    plot(plts)
    ggsave2(paste(path, nameobj, "/post_cut_nFeat_RNA.png", sep=""), plot = plts, width=4, height=2.5)
    
    ##Add metadata from meta
    srt <- AddMetaData(object=srt, list(rep(x = nameobj, times = ncol(srt))), col.name = "Lane") #add nameobj as sample name !!!!
    
    if(is.null(meta)==TRUE){print('No metadata added')}else{
      for(r in 1:nrow(meta)){
        name <- as.character(meta[r,1])
        info <- as.character(meta[r,2])
        srt <- AddMetaData(object = srt, metadata = list(rep(x = info, times = ncol(srt))), col.name = name)
      }
    }
    
    
    #print some stats
    print("min. ncounts per nuclei")
    print(min(srt@meta.data$nCount_RNA))
    print("avg ncounts per nuclei")
    print(mean(srt@meta.data$nCount_RNA))
    print("avg nfeat per nuclei")
    print(mean(srt@meta.data$nFeature_RNA))
    
    
    assign(paste('dirty_feat_', nameobj, sep=""), dirty_feat, envir = .GlobalEnv) #to global environment
    saveRDS(dirty_feat, file = paste(path, nameobj, "/dirty_feat_", nameobj, "_", ".RDS", sep=""))
    assign(paste(nameobj, sep=""), srt, envir = .GlobalEnv) #to global environment
  }, file = paste(curr_dir, "/Sample_preprocess/", nameobj, "_log_load.csv", sep=""))
  
}


#Helper functions
##MT fun
mt_fun <- function(srt, nameobj){
  #Description: Calculates percent of mitochondrial features.
  #Inputs:
  ##srt: seurat object
  #Returns: Seurat object with percent.mt and standard plots
  path <- paste(curr_dir, "/Sample_preprocess/", nameobj, "/Mitochondria_QC/", sep="") #specify dir
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "^MT-") #Check genes that start with MT (mitochondrial genes), usually the higher this is the lower quality the data is
  
  plts <- VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(plts)
  ggsave2(paste(path, "with_mitochondial_genes_Vln.png", sep=""), plot = plts, width=4, height=2.5)
  
  plts <- ggplot(srt@meta.data, aes(percent.mt)) + geom_histogram(bins=5)
  plot(plts)
  ggsave2(paste(path, "with_mitochondial_genes_Hist_percentMT.png", sep=""), plot = plts, width=4, height=2.5)
  
  plts <- FeatureScatter(srt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot(plts)
  ggsave2(paste(path, "with_mitochondial_genes_Sctr.png", sep=""), plot = plts, width=4, height=2.5)
  
  plts <- FeatureScatter(srt, feature1 = "percent.mt", feature2 = "nFeature_RNA")
  plot(plts)
  ggsave2(paste(path, "with_mitochondrial_genes_SctrMT.png", sep=""), plot = plts, width=4, height=2.5)
  
  #Remove mitochondrial genes
  srt <- subset(x = srt, subset = percent.mt < 30) #normally 5%, 30% only for Epi_1!!!!
  plts <- VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(plts)
  ggsave2(paste(path, "without_mitochondial_genes_Vln.png", sep=""), plot = plts, width=4, height=2.5)
  
  plts <- FeatureScatter(srt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot(plts)
  ggsave2(paste(path, "without_mitochondial_genes_sctr.png", sep=""), plot = plts, width=4, height=2.5)
  return(srt)  
  
  plts <- FeatureScatter(srt, feature1 = "percent.mt", feature2 = "nFeature_RNA")
  plot(plts)
  ggsave2(paste(path, "without_mitochondrial_genes_SctrMT.png", sep=""), plot = plts, width=4, height=2.5)
}

##RB_fun
rb_fun <- function(srt, nameobj){  
  #Description: Calculates percent of ribosomal features.
  #Inputs:
  ##srt: seurat object
  #Returns: Seurat object with percent.rb and standard plots
  path <- paste(curr_dir, "/Sample_preprocess/", nameobj, "/Ribosomal_QC/", sep="") #specify dir
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  srt[["percent.rb"]] <- PercentageFeatureSet(object = srt, pattern = "^RP.*") #Check genes that start with MT (mitochondrial genes), usually the higher this is the lower quality the data is
  
  plts <- VlnPlot(srt, features = c("percent.rb"), ncol = 3)
  plot(plts)
  ggsave2(paste(path, "with_ribosomal_genes_VlnRB.png", sep=""), plot = plts, width=4, height=2.5)
  
  plts <- ggplot(srt@meta.data, aes(percent.rb)) + geom_histogram(bins=5)
  plot(plts)
  ggsave2(paste(path, "with_ribosomal_genes_Hist_percentRB.png", sep=""), plot = plts, width=4, height=2.5)
  
  plts <- FeatureScatter(srt, feature1 = "percent.rb", feature2 = "nFeature_RNA")
  plot(plts)
  ggsave2(paste(path, "with_ribosomal_genes_SctrRB.png", sep=""), plot = plts, width=4, height=2.5)
  return(srt)
  
}

##min_pc_fun
min_pc_fun <- function(srt){
  #Description: Calculates the minimum number of PCs
  #Inputs:
  ##srt: seurat object
  #Returns: Returns minimum number of PCs required for Soupx
  stdv <- srt[["pca"]]@stdev
  sum.stdv <- sum(srt[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  return(min.pc)
  print(min.pc)
}

#NSPUC_fun
NFSPUC_fun<- function(sce){
  #Description: Creates a seurat object, normalizes, scales, identifies PCs, neighbours, clusters, and reduces dimensions
  #Inputs:
  ##srt: sce object
  #Returns: Seurat object
  srt <- CreateSeuratObject(sce@assays@data@listData[["counts"]])
  srt <- NormalizeData(srt)
  srt <- FindVariableFeatures(srt, nfeatures = 3000)
  srt <- ScaleData(srt)
  srt <- RunPCA(srt)
  min.pc <- min_pc_fun(srt)
  srt <- FindNeighbors(srt, dims=1:min.pc)
  srt <- RunUMAP(srt, dims=1:min.pc)
  srt <- FindClusters(srt, resolution = 0.1)
  return(srt)
}


#metafun
metafun <- function(srt, projid, Metadata=Metadata, subset_projid=TRUE){
  #Description: Function to add metadata from a CSV and subset individual project IDs.
  #Input
  ##srt: seurat object
  ##projid: projid of sample
  #Returns: Subset projid from seurat object, default=TRUE
  nameobj_out <- paste("ChP_", projid, sep="")
  if(subset_projid==TRUE){
    srt <- subset(srt, sample == projid)
  }
  #addmetadata
  meta <- srt@meta.data
  metaadd <- Metadata[Metadata$projid==projid,]
  metaadd <- as.data.frame(lapply(metaadd, rep, nrow(meta)))
  meta <- cbind(meta, metaadd)
  srt@meta.data <- meta
  #save
  saveRDS(srt, file=paste(curr_dir, "/Sample_preprocess/", nameobj_out, ".RDS", sep=""))
  assign(paste(nameobj_out), srt, envir = .GlobalEnv)  #to global environment
}

#metafun
metafunB <- function(srt, projid, subset_projid=TRUE){
  #Description: Function to add metadata from a CSV and subset individual project IDs. !!!!For samples to which we are adding cells.
  #Input
  ##srt: seurat object
  ##projid: projid of sample
  #Returns: Subset projid from seurat object, default=TRUE
  nameobj_out <- paste("ChP_", projid, sep="")
  if(subset_projid==TRUE){
    srt <- subset(srt, sample == projid)
  }
  #addmetadata
  meta <- srt@meta.data
  metaadd <- Metadata[Metadata$projid==projid,]
  metaadd <- as.data.frame(lapply(metaadd, rep, nrow(meta)))
  meta <- cbind(meta, metaadd)
  srt@meta.data <- meta
  #save
  saveRDS(srt, file=paste(curr_dirF, "/Sample_preprocess/", nameobj_out, "B.RDS", sep=""))
  assign(paste(nameobj_out), srt, envir = .GlobalEnv)  #to global environment
}

##demultiplexfun
demultiplexfun <- function(srt, assay="CMO", p.quant_CMO1=0.99, p.quant_CMO2=0.99){
  #Description: Demultiplexes data with specific parameters and makes standard plots.
  #Inputs:
  ##seq_sample: seurat object with CMOs
  ##assay: Specify assay to demultiplex, default is CMO
  ##p.quant_CMO1: Affects the minimum threshold for the top CMO. AKA Higher number decreases singlet/doublet, increases negative. Passed to HTODemux_low_highfun; default is 0.99 as per Seurat's original HTODemux. 
  ##p.quant_CMO2: Affects the minimum threshold for the second top CMO. AKA Higher number decreases doublets, increases singlets, no effect on negatives. Passed to HTODemux_low_highfun; default is 0.99 as per Seurat's original HTODemux.
  #Returns: Modified seurat object and plots
  
  nameobj <- deparse(substitute(srt)) #make char
  path <- paste(curr_dir, "/Sample_preprocess/", nameobj, "/Demultiplex/", sep="") #specify dir
  DefaultAssay(srt) <- assay
  
  if(dir.exists(path)==FALSE){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  capture.output({ #capture entire output
    #Filter low and high expressing cells
    plts <- VlnPlot(srt, features = "nFeature_RNA")
    plot(plts)
    ggsave2(paste(path, "features.png", sep=""), plot = plts, width=4, height=2.5)
    
    plts <- ggplot(srt@meta.data, aes(nFeature_RNA)) + geom_histogram(bins=100)
    plot(plts)
    ggsave2(paste(path, "nFeatures_Hist.png", sep=""), plot = plts, width=4, height=2.5)
    
    #Demultiplex
    if(assay!="CMO_dsb"){srt <- NormalizeData(srt, assay = assay, normalization.method = "CLR")} #Normalize CMO #Don't renormalize if you used dsb!!!!
    
    #srt <- HTODemux(srt, assay = assay, positive.quantile=p.quant) #!!!!
    srt <- HTODemux_low_highfun(srt, assay = assay, p.quant_CMO1=p.quant_CMO1, p.quant_CMO2=p.quant_CMO2) #use modified HTODemux function (below) !!!!
    
    #Check number of singlets, doublets, negative
    nSDN <- table(srt[[paste(assay, "_classification.global", sep="")]])
    print(nSDN)
    write.csv(nSDN, paste(path, "Number_of_SDNs.csv", sep=""))
    
    #violin of above
    Idents(srt) <- paste(assay, "_classification.global", sep="")
    plts <- VlnPlot(srt, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
    plot(plts)
    ggsave2(paste(path, assay, "_vln_SND.png", sep=""), plot = plts, width=4, height=2.5)
    
    #Check CMO expression
    plts <- VlnPlot(srt, features = paste("nCount_", assay, sep=""), split.by = paste(assay, "_maxID", sep=""))
    plot(plts)
    ggsave2(paste(path, assay, "_all_ncounts.png", sep=""), plot = plts, width=4, height=2.5)
    
    plts <- VlnPlot(srt, features = paste("nCount_", assay, sep=""), split.by = paste(assay, "_classification.global", sep=""))
    plot(plts)
    ggsave2(paste(path, assay, "_SND_ncounts.png", sep=""), plot = plts, width=4, height=2.5)
    
    plts <- ggplot(srt@meta.data, aes(x=srt@meta.data[[paste("nCount_", assay, sep="")]])) + geom_histogram(bins=100)
    plot(plts)
    ggsave2(paste(path, assay, "_nCount_Hist.png", sep=""), plot = plts, width=4, height=2.5)
    
    #Per CMO
    #ridge plots per CMO
    Idents(srt) <- paste(assay, "_maxID", sep="")
    for(i in rownames(srt[[assay]]))
    {
      plts <- RidgePlot(srt, assay = assay, features = i, ncol = 1)
      plot(plts)
      ggsave2(paste(path, assay, "_ridge_", i, ".png", sep=""), plot = plts, width=4, height=2.5)
    }
    
    #Scatter plot CMOvCMO
    for(i in 1:length(rownames(srt[[assay]]))){
      for(j in 1:length(rownames(srt[[assay]]))){
        if(i<=j)
          next
        feat1 = rownames(srt[[assay]])[i]
        feat2 = rownames(srt[[assay]])[j]
        plts <- FeatureScatter(srt, feature1 = feat1, feature2 = feat2, pt.size=0.01)
        plot(plts)
        ggsave2(paste(path, assay, "_feature_", feat1, "_vs_", feat2, ".png", sep=""), plot = plts, width=4, height=2.5)
      }
    }
    
    plts <- HTOHeatmap(srt, assay = assay)
    plot(plts)
    ggsave2(paste(path, assay, "_heatmap.png", sep=""), plot = plts, width=4, height=2.5)
    
    if(ncol(srt)<100){print("Object contains too few cells to RunTSNE ... skipping tsne")
    }else{if(ncol(srt)<40000){
      hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = srt, assay = assay))))
      srt <- RunTSNE(srt, distance.matrix = hto.dist.mtx)
      Idents(srt) <- paste(assay, "_classification.global", sep="")
      plot(DimPlot(srt))
      ggsave2(paste(path, assay, "_tsne_NSD.png", sep=""), plot = DimPlot(srt), width=4, height=2.5)
      
      Idents(srt) <- paste(assay, "_maxID", sep="")#CMO plot
      plot(DimPlot(srt))
      ggsave2(paste(path, assay, "_tsne_CMO.png", sep=""), plot = DimPlot(srt), width=4, height=2.5)
      
      meta <- srt@meta.data #extract metadata
      #Correlate CMO to Doublet Finder
      plts <- ggplot(meta, aes(x=meta[[paste(assay, "_classification.global", sep="")]], y=DoubletFinder_score)) + 
        geom_violin(trim=FALSE) +
        geom_text(x = 2, y = 0.8, label = lm_eqn(meta[[paste(assay, "_classification.global", sep="")]], meta$DoubletFinder_score), parse = TRUE)
      plot(plts)
      ggsave2(paste(path, assay, "vsDoubletFinder.png", sep=""), plot = plts, width=4, height=2.5)
      
    }else{print("Object contains too many cells to RunTSNE ... skipping tsne")}}
  
    ##remove doublets and negatives
    Idents(srt) <- paste(assay, "_classification.global", sep="")
    #srt <- subset(srt, idents="Singlet")
    
    DefaultAssay(srt) <- "RNA"
    #to global environment
    assign(paste(nameobj), srt, envir = .GlobalEnv)
  }, file = paste(curr_dir, "/Sample_preprocess/", nameobj, "_log_demultiplex.csv", sep=""))
}

##Linear model equation to add to ggplot
lm_eqn <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

##HTOdemux dependencies
MaxN <- function(x, N = 2){
  len <- length(x)
  if (N > len) {
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
}


##Modified HTO demux
HTODemux_low_highfun <- function (object, assay = "HTO", p.quant_CMO1 = 0.99, p.quant_CMO2 = 0.99, init = NULL, 
                                  nstarts = 100, kfunc = "clara", nsamples = 100, seed = 42, 
                                  verbose = TRUE) 
{
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(object = object, assay = assay, slot = "counts")[, 
                                                                          colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(EXPR = kfunc, kmeans = {
    init.clusters <- kmeans(x = t(x = GetAssayData(object = object, 
                                                   assay = assay)), centers = ncenters, nstart = nstarts)
    Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
  }, clara = {
    init.clusters <- clara(x = t(x = GetAssayData(object = object, 
                                                  assay = assay)), k = ncenters, samples = nsamples)
    Idents(object = object, cells = names(x = init.clusters$clustering), 
           drop = TRUE) <- init.clusters$clustering
  }, stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'"))
  average.expression <- AverageExpression(object = object, 
                                          assays = assay, verbose = FALSE)[[assay]]
  discrete_weak <- GetAssayData(object = object, assay = assay)
  discrete_weak[discrete_weak > 0] <- 0
  discrete_strong <- discrete_weak
  
  ###start modified section
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, 
    ])]])]
    fit <- suppressWarnings(expr = fitdist(data = values.use, 
                                           distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = p.quant_CMO1)$quantiles[1])
    discrete_weak[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, 
                     " reads"))
    }
  }
  
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, 
    ])]])]
    fit <- suppressWarnings(expr = fitdist(data = values.use, 
                                           distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = p.quant_CMO2)$quantiles[1])
    discrete_strong[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, 
                     " reads"))
    }
  }
  
  npositive_weak <- colSums(x = discrete_weak)
  npositive_weak[npositive_weak >= 1] <- "Singlet"        
  npositive_weak[npositive_weak == 0] <- "Negative"
  
  npositive_strong <- colSums(x = discrete_strong)
  npositive_strong[npositive_strong >= 2] <- "Doublet"
  npositive_strong[npositive_strong == 1 | npositive_strong == 0] <- "Singlet"
  
  
  classification.global <- data.frame(npositive_weak, npositive_strong)
  classification.global$npositive_strong[classification.global$npositive_weak == "Negative"] <- "Negative"
  classification.global <- classification.global$npositive_strong
  ###End modififed section
  
  donor.id = rownames(x = data)
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
  hash.maxID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data), 
                                                 FUN = function(x) {
                                                   return(which(x = data[, x] == hash.max[x])[1])
                                                 })])
  hash.secondID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data), 
                                                    FUN = function(x) {
                                                      return(which(x = data[, x] == hash.second[x])[1])
                                                    })])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == 
                                                                           "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == 
                                                                           "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID, 
                                        hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, c("maxID", 
                                                          "secondID", "margin", "classification", "classification.global"), 
                                                 sep = "_")
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, "_classification")
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, 
                                                            "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  object$hash.ID <- Idents(object = object)
  return(object)
}
