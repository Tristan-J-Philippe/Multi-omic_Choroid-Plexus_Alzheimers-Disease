###
### 021 Cell_type_annotation_RNA.R
###
# Purpose: Cell type annotation and cleaning of single nuclei RNA dataset.
# Dependencies:
source("/ChP_RNA/020 Cell_type_annotation_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(1234)


##Load and merge samples (note all projIDs replaced with fkID)
ChP_101 <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_101.RDS", sep=""))
...
ChP_154 <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_154.RDS", sep=""))

#merge first half due to limits on total number of cells and computational resources
Ex_ChP <- merge(x=ChP_101, y=c(ChP_104, ChP_106, ChP_107, ChP_108, ChP_109, ChP_110, ChP_111, ChP_112, ChP_113, ChP_114, ChP_116, ChP_117, ChP_118, ChP_119, ChP_120, ChP_121, ChP_122, ChP_122B, ChP_124, ChP_125, ChP_126, ChP_128, ChP_129, ChP_131, ChP_131B, ChP_132, ChP_133, ChP_134, ChP_135, ChP_136, ChP_137, ChP_138, ChP_140, ChP_141B, ChP_144, ChP_145))
rm(list = grep("ChP_", ls(), value=TRUE))
gc()

#second half
Ex_ChPB <- merge(x=ChP_145B, y=c(ChP_147, ChP_148, ChP_149, ChP_150, ChP_151, ChP_152, ChP_153, ChP_154, ChP_155, ChP_156, ChP_157, ChP_158, ChP_159, ChP_160, ChP_161, ChP_162, ChP_163, ChP_163B, ChP_164, ChP_165, ChP_166, ChP_167, ChP_168B, ChP_168, ChP_170, ChP_172, ChP_174, ChP_175, ChP_176, ChP_177, ChP_178, ChP_178B, ChP_179, ChP_180, ChP_181, ChP_182, ChP_183, ChP_184, ChP_185))
#rm(list = grep("ChP_", ls(), value=TRUE))
gc()

Ex_ChP <- merge(Ex_ChP, Ex_ChP2)
rm(Ex_ChP2)
gc()

Ex_ChP@meta.data$Lane <- Idents(Ex_ChP)
unique(Ex_ChP@meta.data$Lane) #sanity check
Ex_ChP@meta.data$Lane_projid <- paste(Ex_ChP@meta.data$projid, Ex_ChP@meta.data$Lane, sep="_")
Ex_ChP@meta.data$sample <- Ex_ChP@meta.data$projid
gc()
Ex_ChP@meta.data[is.na(Ex_ChP@meta.data)==TRUE] <- "unk"

L_p <- unique(Ex_ChP@meta.data$Lane_projid)
L_p[order(L_p)] #Check this when making metadata for DEG later





#Celltyping
##Wholeset for each part
###NFSP, Harmonize, CNDR
ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)
NFSP_fun(Ex_ChP, assay="soupXcounts", spikein = ChP_genes)
use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)

Ex_ChP <- RunHarmony(object = Ex_ChP, group.by.vars = c("Lane"), dims.use = use_dims, assay="soupXcounts", max.iter.harmony = 50, plot_convergence = TRUE) #Correct for batches !!!!
#CNDRfun_quick(Ex_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony", resolution=0.1)
#Ex_ChP <- FindClusters(Ex_ChP, resolution = c(0.01), assay=assay)
#plts <- clustree(Ex_ChP)
#  plot(plts)
#  graph2png(plts, paste(curr_dir, "/CellType/", "Ex_ChP/Ex_ChP_cluster_tree.png", sep=""), width=14, height=6)
CNDRfun(Ex_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")


###Choose resolution
res="0.1" #Pick res
Idents(Ex_ChP) <- Ex_ChP@meta.data[[paste("soupXcounts_snn_res.", res, sep="")]]
CNDRdimfun(Ex_ChP, reduction = "harmony", res=res, split.by = c("X10X_batch", "projid", "Lane", "CMO_classification.global", "DoubletFinder_call"))
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per Lane.csv", sep=""), "Lane")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per 10X_batch.csv", sep=""), "X10X_batch")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per projid.csv", sep=""), "projid")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
DimPlot(Ex_ChP, raster=T)

return_markersfun(Ex_ChP, assay="soupXcounts", features = Ex_ChP_soupXcounts_var_feat)###Find markers

###Dividing cells spend more time in S phase
Ex_ChP <- CellCycleScoring(Ex_ChP, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident=FALSE)
DimPlot(Ex_ChP, group.by = "Phase", raster=F)
DimPlot(Ex_ChP, split.by = "Phase", raster=F)
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/CMO per phase.csv", "Phase"))



curr_clust <- 0 #iterate over top 30 markers, number of markers per set and which markers
Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

Top20_Featuresfun(Ex_ChP, curr_clust)

search_markersfun(Ex_ChP, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Ex_ChP, curr_clust, known_fibro_prog)

search_markersfun(Ex_ChP, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_light) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Ex_ChP, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Ex_ChP, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Ex_ChP, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Ex_ChP, curr_clust, known_dividing %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Ex_ChP, "LAMA4", raster=F)


#Ex_ChPA
#0: Epi_1
#1: Epi_2
#2: Fibro
#3: Oligo
#4: Neuro?
#5: Imm
#6: Endo
#7: Astro

#Ex_ChPB
#0: Epi_2
#1: Epi_1
#2: Fibro
#3: Epi_1
#4: Imm
#5: Endo
#6: Neuro
#7: Oligo
#8: Astro
#9: Neuro
#10: Oligo

###Rename all identities
Ex_ChP <- RenameIdents(Ex_ChP, "0"="Epi_2", "1"="Epi_1", "2"="Fibro", "3"="Epi_1", "4"="Imm", "5"="Endo", "6"="Neuro", "7"="Oligo", "8"="Astro", "9"="Neuro", "10"="Oligo")
Ex_ChP@meta.data$Firstcelltype <- Idents(Ex_ChP)
saveRDS(Ex_ChP, paste(curr_dir, "/CellType/Ex_ChP.RDS", sep=""))


###Subtypes
for(i in unique(Ex_ChP@meta.data$Firstcelltype)){
  sub <- subset(Ex_ChP, Firstcelltype==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/", i, "A.RDS", sep=""))
}



##Epi_1
###PCA
Epi_1A <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_1A.RDS"))
Epi_1B <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_1B.RDS"))
Epi_1 <- merge(Epi_1A, Epi_1B) #merge both batches
rm(Epi_1A, Epi_1B)
gc()

ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)
NFSP_fun(Epi_1, assay="soupXcounts", spikein=ChP_genes)

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Epi_1 <- RunHarmony(object = Epi_1, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Epi_1, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

Idents(Epi_1) <- Epi_1@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Epi_1, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global", "Doublet"))
return_markersfun(Epi_1, assay="soupXcounts", features = Epi_1_soupXcounts_var_feat)

for(i in c(0:3)){
  sub <- subset(Epi_1, soupXcounts_snn_res.0.1==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/Epi_1.", i, ".RDS", sep=""))
}

saveRDS(Epi_1, paste0(curr_dir, "/CellType/Firstcelltype/Epi_1.RDS"))



curr_clust <- 0
Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

Top20_Featuresfun(Epi_1, curr_clust)

search_markersfun(Epi_1, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Epi_1, curr_clust, known_Fibro_prog)

search_markersfun(Epi_1, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_light) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Epi_1, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Epi_1, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Epi_1, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Epi_1, curr_clust, known_dividing %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Epi_1, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Epi_1, "TTR", raster=F)


###Label
Epi_1 <- RenameIdents(Epi_1, "0"="Fib_3", "2"="Fib_3", "1"="Fib_2", "3"="Fib_2")
Epi_1@meta.data$Subcelltype <- Epi_1@active.ident
DoubletE2 <- subset(Epi_2, Subcelltype=="Doublet")

#remove doublets
Epi_2 <- subset(Epi_2, CMO_classification.global=="Singlet")
Epi_2 <- subset(Epi_2, percent.mt<5)

Epi_2 <- subset(Epi_2, Subcelltype=="Epi_2", spikein=known_dark)

saveRDS(Epi_1, paste0(curr_dir, "/CellType/Secondcelltype/Epi_1.RDS"))


###Venn of agreement
plts <- ggVennDiagram(list(colnames(Epi_1)[Epi_1@meta.data$DoubletFinder_call=="Doublet"], colnames(Epi_1)[Epi_1@meta.data$CMO_classification.global=="Doublet"], colnames(Epi_1)[Epi_1@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Epi_1)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Epi_1/Epi_1 DoubletFinder-CMO-Manual venn.png", sep=""))


###Remove doublets

Epi_1 <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_1.RDS"))
Epi_1 <- subset(Ex_ChPB, subset=c(Firstcelltype=="Epi_1"))
Epi_1_all <- Epi_1

Featurefun(Epi_1, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")

#Remove cells that don't belong
#Epi_1 <- subset(Epi_1, HTR2C < 3)
#Epi_1 <- subset(Epi_1, TTR < 3)
#Epi_1 <- subset(Epi_1, CLIC6 < 3)
#Epi_1 <- subset(Epi_1, AQP1 < 3)
#Epi_1 <- subset(Epi_1, CFAP54 < 3)
#Epi_1 <- subset(Epi_1, FAM227A < 3)
#Epi_1 <- subset(Epi_1, ARL13B < 3)
#Epi_1 <- subset(Epi_1, TRPM3 < 3)
Epi_1 <- subset(Epi_1, SLC4A4 < 3)
Epi_1 <- subset(Epi_1, LAMA2 < 3)
Epi_1 <- subset(Epi_1, LAMA4 < 3)
Epi_1 <- subset(Epi_1, CLDN1 < 3)
Epi_1 <- subset(Epi_1, CLDN11 < 3)
Epi_1 <- subset(Epi_1, PRICKLE1 < 3)
Epi_1 <- subset(Epi_1, CEMIP < 3)
Epi_1 <- subset(Epi_1, FOXP2 < 3)
Epi_1 <- subset(Epi_1, COL6A3 < 3)
Epi_1 <- subset(Epi_1, SLC6A13 < 3)
Epi_1 <- subset(Epi_1, MECOM < 3)
Epi_1 <- subset(Epi_1, VWF < 3)
Epi_1 <- subset(Epi_1, TAGLN < 3)
Epi_1 <- subset(Epi_1, MYH11 < 3)
Epi_1 <- subset(Epi_1, NOTCH3 < 3)
Epi_1 <- subset(Epi_1, CADM2 < 3)
Epi_1 <- subset(Epi_1, MAP2 < 3)
Epi_1 <- subset(Epi_1, GRIK2 < 3)
#Epi_1 <- subset(Epi_1, CLU < 3)
#Epi_1 <- subset(Epi_1, AQP4 < 3)
Epi_1 <- subset(Epi_1, MOBP < 3)
Epi_1 <- subset(Epi_1, MBP < 3)
Epi_1 <- subset(Epi_1, CDK19 < 3)
Epi_1 <- subset(Epi_1, CD68 < 3)
Epi_1 <- subset(Epi_1, F13A1 < 3)
Epi_1 <- subset(Epi_1, STAB1 < 3)
Epi_1 <- subset(Epi_1, CD163 < 3)
Epi_1 <- subset(Epi_1, LYVE1 < 3)
Epi_1 <- subset(Epi_1, P2RY12 < 3)
Epi_1 <- subset(Epi_1, THEMIS < 3)

Epi_1 <- CellSelector(DimPlot(Epi_1), Epi_1, ident = "Epi_1")
Epi_1 <- CellSelector(DimPlot(Epi_1), Epi_1, ident = "Doublet")
#Epi_1@meta.data$Immtemp <- Idents(Epi_1)
Epi_1 <- subset(Epi_1, ident="Epi_1")
#Epi_1@meta.data$Immtemp <- NULL
label_doubletsfun(Epi_1, Epi_1_all)
Epi_1 <- Epi_1_all


##Epi_2
###PCA

Epi_2A <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_2A.RDS"))
Epi_2B <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_2B.RDS"))
Epi_2 <- merge(Epi_2A, Epi_2B)
rm(Epi_2A, Epi_2B)
gc()

ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)
NFSP_fun(Epi_2, assay="soupXcounts", spikein = ChP_genes) #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Epi_2 <- RunHarmony(object = Epi_2, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Epi_2, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Epi_2) <- Epi_2@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Epi_2, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Epi_2, assay="soupXcounts", features = Epi_2_soupXcounts_var_feat)

for(i in c(0:3)){
  sub <- subset(Epi_2, soupXcounts_snn_res.0.1==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/Epi_2.", i, ".RDS", sep=""))
}

saveRDS(Epi_2, paste0(curr_dir, "/CellType/Firstcelltype/Epi_2.RDS"))



curr_clust <- 3
Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Epi_2_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Epi_2, curr_clust)

search_markersfun(Epi_2, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Epi_2, curr_clust, known_Fibro_prog)

search_markersfun(Epi_2, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Epi_2, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Epi_2, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Epi_2, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Epi_2, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Epi_2, curr_clust, c(known_light) %>% unique())
search_markersfun(Epi_2, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Epi_2, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Epi_2, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Epi_2, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Epi_2, curr_clust, known_dividing %>% unique())

search_markersfun(Epi_2, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Epi_2, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Epi_2, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Epi_2, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Epi_2, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Epi_2, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Epi_2, "TTR", raster=F)

###Label

Epi_2 <- RenameIdents(Epi_2, "0"="Fib_3", "2"="Fib_3", "1"="Fib_2", "3"="Fib_2")
Epi_2@meta.data$Subcelltype <- Epi_2@active.ident
DoubletE2 <- subset(Epi_2, Subcelltype=="Doublet")

#remove doublets
Epi_2 <- subset(Epi_2, CMO_classification.global=="Singlet")
Epi_2 <- subset(Epi_2, percent.mt<5)

Epi_2 <- subset(Epi_2, Subcelltype=="Epi_2")

saveRDS(Epi_2, paste0(curr_dir, "/CellType/Secondcelltype/Epi_2.RDS"))


###Venn of agreement

plts <- ggVennDiagram(list(colnames(Epi_2)[Epi_2@meta.data$DoubletFinder_call=="Doublet"], colnames(Epi_2)[Epi_2@meta.data$CMO_classification.global=="Doublet"], colnames(Epi_2)[Epi_2@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Epi_2)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Epi_2/Epi_2 DoubletFinder-CMO-Manual venn.png", sep=""))


###remove doublets

Epi_2 <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_2.RDS"))
Epi_2_all <- Epi_2

Featurefun(Epi_2, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong
Epi_2 <- subset(Epi_2, HTR2C < 4)
Epi_2 <- subset(Epi_2, TTR < 4)
Epi_2 <- subset(Epi_2, CLIC6 < 4)
#Epi_2 <- subset(Epi_2, AQP1 < 3)
Epi_2 <- subset(Epi_2, CFAP54 < 3)
Epi_2 <- subset(Epi_2, FAM227A < 3)
Epi_2 <- subset(Epi_2, ARL13B < 3)
#Epi_2 <- subset(Epi_2, TRPM3 < 3)
#Epi_2 <- subset(Epi_2, SLC4A4 < 3)
#Epi_2 <- subset(Epi_2, LAMA2 < 3)
#Epi_2 <- subset(Epi_2, LAMA4 < 3)
#Epi_2 <- subset(Epi_2, CLDN1 < 3)
#Epi_2 <- subset(Epi_2, CLDN11 < 3)
#Epi_2 <- subset(Epi_2, PRICKLE1 < 3)
#Epi_2 <- subset(Epi_2, CEMIP < 3)
#Epi_2 <- subset(Epi_2, FOXP2 < 3)
#Epi_2 <- subset(Epi_2, COL6A3 < 3)
#Epi_2 <- subset(Epi_2, SLC6A13 < 3)
Epi_2 <- subset(Epi_2, MECOM < 3)
Epi_2 <- subset(Epi_2, VWF < 3)
Epi_2 <- subset(Epi_2, TAGLN < 3)
Epi_2 <- subset(Epi_2, MYH11 < 3)
Epi_2 <- subset(Epi_2, NOTCH3 < 3)
Epi_2 <- subset(Epi_2, CADM2 < 3)
Epi_2 <- subset(Epi_2, MAP2 < 3)
Epi_2 <- subset(Epi_2, GRIK2 < 3)
#Epi_2 <- subset(Epi_2, CLU < 3)
#Epi_2 <- subset(Epi_2, AQP4 < 3)
Epi_2 <- subset(Epi_2, MOBP < 3)
Epi_2 <- subset(Epi_2, MBP < 3)
Epi_2 <- subset(Epi_2, CDK19 < 3)
Epi_2 <- subset(Epi_2, CD68 < 3)
Epi_2 <- subset(Epi_2, F13A1 < 3)
Epi_2 <- subset(Epi_2, STAB1 < 3)
Epi_2 <- subset(Epi_2, CD163 < 3)
Epi_2 <- subset(Epi_2, LYVE1 < 3)
Epi_2 <- subset(Epi_2, P2RY12 < 3)
Epi_2 <- subset(Epi_2, THEMIS < 3)

Epi_2 <- CellSelector(DimPlot(Epi_2), Epi_2, ident = "Epi_2")
Epi_2 <- CellSelector(DimPlot(Epi_2), Epi_2, ident = "Doublet")
#Epi_2@meta.data$Immtemp <- Idents(Epi_2)
Epi_2 <- subset(Epi_2, ident="Epi_2")
#Epi_2@meta.data$Immtemp <- NULL
label_doubletsfun(Epi_2, Epi_2_all)
Epi_2 <- Epi_2_all

##Epi
###PCA

Epi_1 <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_1.RDS"))
Epi_2 <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_2.RDS"))
Epi <- merge(Epi_1, Epi_2)
rm(Epi_1, Epi_2)
gc()

Epi_soupXcounts_var_feat <- VariableFeatures(Epi)
Epi <- RunPCA(Epi, assay="soupXcounts", features=Epi_soupXcounts_var_feat)

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Epi <- RunHarmony(object = Epi, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Epi, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Epi, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
Ncellfun(Epi, paste0(curr_dir, "/CellType/Epi/Final per AD.csv", sep=""), "Cog_Path")
#return markers
return_markersfun(Epi, assay="soupXcounts", features = Epi_soupXcounts_var_feat)
Epi <- RenameIdents(Epi, "0"="Epi_1", "1"="Epi_2", "2"="Epi_3", "3"="Epi_4")
Epi@meta.data$Subcelltype <- Idents(Epi)

Epi <- RenameIdents(Epi, "Epi_0"="0", "Epi_1"="1", "Epi_2"="2", "Epi_3"="3")
return_markersfun(Epi, assay="soupXcounts", features = Epi_soupXcounts_var_feat)

Featurefun(Epi, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "GPX3", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "XIST", "XACT"), "KeyMarkers",  assay="soupXcounts")

saveRDS(Epi, paste0(curr_dir, "/CellType/Finalcelltype/Epi.RDS"))

#1 Epi_2 few clear defining features
#3 Epi_3 MTs
#0 Epi_1 CFAPs
#2,4,5,6 Epi_4 CFAPs, MTs, GPX3, TUBB4B, SPAGs, DNAs


###Remove interdoublets

Idents(Epi) <- Epi@meta.data$Subcelltype
DimPlot(Epi, raster=F)
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi1")
Epi_1 <- subset(Epi, ident="Epi1")
Epi_1 <- subset(Epi_1, Subcelltype=="Epi_1")


Idents(Epi) <- Epi@meta.data$Subcelltype
DimPlot(Epi)
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi2")
Epi_2 <- subset(Epi, Subcelltype=="Epi_2")
Epi_2 <- subset(Epi_2, ident="Epi2")


Idents(Epi) <- Epi@meta.data$Subcelltype
DimPlot(Epi)
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi3")
Epi_3 <- subset(Epi, Subcelltype=="Epi_3")
Epi_3 <- subset(Epi_3, ident="Epi3")


Idents(Epi) <- Epi@meta.data$Subcelltype
DimPlot(Epi)
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi4")
Epi_4 <- subset(Epi, Subcelltype=="Epi_4")
Epi_4 <- subset(Epi_4, ident="Epi4")

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

Epi <- RunUMAP(object = Epi, dims = 1:20, assay = "soupXcounts", reduction =  "harmony")


###Subtypes

Epi_1 <- subset(Epi, Finalcelltype=="Epi_2")



NFSP_fun(Epi_1, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Epi_1 <- RunHarmony(object = Epi_1, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Epi_1, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Epi_1) <- Epi_1@meta.data[["soupXcounts_snn_res.0.7"]]
CNDRdimfun(Epi_1, reduction = "harmony", res="_0.7")
#return markers
return_markersfun(Epi_1, assay="soupXcounts", features = Epi_1_soupXcounts_var_feat, num_clust=0:3)
saveRDS(Epi_1, paste0(curr_dir, "/CellType/Finalcelltype/Epi_1__Ependymal.RDS"))



curr_clust <- 3
Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:50]
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Epi_1_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Epi_1, curr_clust)

search_markersfun(Epi_1, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Epi_1, curr_clust, known_Fibro_prog)

search_markersfun(Epi_1, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_light) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Epi_1, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Epi_1, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Epi_1, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Epi_1, curr_clust, known_dividing %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Epi_1, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Epi_1, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Epi_1, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Epi_1, "MT-CO2", raster=F)


##Fibro
###PCA

Fibro <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Fibro.RDS"))
Fibro <- subset(Ex_ChP, subset=Firstcelltype=="Fibro")
Fibro_all <- Fibro

FibroA <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/FibroA.RDS"))
FibroB <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/FibroB.RDS"))
Fibro <- merge(FibroA, FibroB)
rm(FibroA, FibroB)
gc()

ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)
NFSP_fun(Fibro, assay="soupXcounts", spikein= ChP_genes)

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Fibro <- RunHarmony(object = Fibro, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Fibro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
#Idents(Fibro) <- Fibro@meta.data$Subcelltype
Idents(Fibro) <- Fibro@meta.data[["soupXcounts_snn_res.0.3"]]
CNDRdimfun(Fibro, reduction = "harmony", res="0.3", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Fibro, assay="soupXcounts", features = Fibro_soupXcounts_var_feat)
saveRDS(Fibro, paste0(curr_dir, "/CellType/Firstcelltype/Fibro.RDS"))



curr_clust <- 0
Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Fibro, curr_clust)

search_markersfun(Fibro, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, "CLDN11", "TJP1", "TJP2") %>% unique())
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
FeaturePlot(Fibro, "NRXN3", raster=F)

#Fib_2: LAMA2 TRPM3: 1 keep 7
#Fib_1: LAMA2: 3, 11
#Fib_3: SLC4A4: 0, 4 keep 2, 5 10
#Fib_4: SLC4A4 TRPM3: 6, 9
#Mural: 8

#remove doublets and Mural res 0.3
#Mural: 8
#Doublets: 7, 2, 5, 10

#2nd pass
#Fib_2: 0
#Fib_4: 7
#Fib_3: 1,2,3,6
#Fib_1: 4,5
#8 Doublet

###remove doublets

Featurefun(Fibro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "XIST", "XACT"), "Test2",  assay="soupXcounts")

#Remove cells that don't belong
Fibro <- subset(Fibro, HTR2C < 3)
Fibro <- subset(Fibro, TTR < 3)
Fibro <- subset(Fibro, CLIC6 < 3)
#Fibro <- subset(Fibro, AQP1 < 3)
Fibro <- subset(Fibro, CFAP54 < 3)
Fibro <- subset(Fibro, FAM227A < 3)
Fibro <- subset(Fibro, ARL13B < 3)
#Fibro <- subset(Fibro, TRPM3 < 3)
#Fibro <- subset(Fibro, SLC4A4 < 3)
#Fibro <- subset(Fibro, LAMA2 < 3)
#Fibro <- subset(Fibro, LAMA4 < 3)
#Fibro <- subset(Fibro, CLDN1 < 3)
#Fibro <- subset(Fibro, CLDN11 < 3)
#Fibro <- subset(Fibro, PRICKLE1 < 3)
#Fibro <- subset(Fibro, CEMIP < 3)
#Fibro <- subset(Fibro, FOXP2 < 3)
#Fibro <- subset(Fibro, COL6A3 < 3)
#Fibro <- subset(Fibro, SLC6A13 < 3)
Fibro <- subset(Fibro, MECOM < 3)
Fibro <- subset(Fibro, VWF < 3)
Fibro <- subset(Fibro, TAGLN < 3)
Fibro <- subset(Fibro, MYH11 < 3)
Fibro <- subset(Fibro, NOTCH3 < 3)
Fibro <- subset(Fibro, CADM2 < 3)
Fibro <- subset(Fibro, MAP2 < 3)
Fibro <- subset(Fibro, GRIK2 < 3)
#Fibro <- subset(Fibro, CLU < 3)
#Fibro <- subset(Fibro, AQP4 < 3)
Fibro <- subset(Fibro, MOBP < 3)
Fibro <- subset(Fibro, MBP < 3)
Fibro <- subset(Fibro, CDK19 < 3)
Fibro <- subset(Fibro, CD68 < 3)
Fibro <- subset(Fibro, F13A1 < 3)
Fibro <- subset(Fibro, STAB1 < 3)
Fibro <- subset(Fibro, CD163 < 3)
Fibro <- subset(Fibro, LYVE1 < 3)
Fibro <- subset(Fibro, P2RY12 < 3)
Fibro <- subset(Fibro, THEMIS < 3)


###Label

#Firstpass
Fibro <- RenameIdents(Fibro, "0"="Fib", "1"="Fib", "2"="Fib", "3"="Fib", "4"="Fib", "5"="Doublet", "6"="Doublet", "7"="Doublet", "8"="Fib", "9"="Mural", "10"="Fib") #res 0.3
Fibro@meta.data$Subcelltype <- Fibro@active.ident
MuralB <- subset(Fibro, Subcelltype=="Mural")
Fibro <- subset(Fibro, Subcelltype=="Fib")

Fibro <- subset(Fibro, CMO_classification.global=="Singlet")
Fibro <- subset(Fibro, percent.mt<5)
Fibro <- subset(Fibro, nFeature_soupXcounts>300)

saveRDS(Fibro, paste0(curr_dir, "/CellType/Secondcelltype/Fibro.RDS"))


###Venn of agreement

plts <- ggVennDiagram(list(colnames(Fibro)[Fibro@meta.data$DoubletFinder_call=="Doublet"], colnames(Fibro)[Fibro@meta.data$CMO_classification.global=="Doublet"], colnames(Fibro)[Fibro@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Fibro)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Fibro/Fibro DoubletFinder-CMO-Manual venn.png", sep=""))



###re-PCA

Fibro <- readRDS(paste0(curr_dir, "/CellType/Secondcelltype/Fibro.RDS"))

NFSP_fun(Fibro, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Fibro <- RunHarmony(object = Fibro, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Fibro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Fibro) <- Fibro@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Fibro, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
Ncellfun(Fibro, paste0(curr_dir, "/CellType/Fibro/Final per AD.csv", sep=""), "Cog_Path")
#return markers
return_markersfun(Fibro, assay="soupXcounts", features = Fibro_soupXcounts_var_feat)

Featurefun(Fibro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "KeyMarkers",  assay="soupXcounts")

Fibro <- RenameIdents(Fibro, "0"="Fib_1", "1"="Fib_2", "2"="Fib_3", "3"="Fib_1")
Idents(Fibro) -> Fibro@meta.data$Subcelltype

saveRDS(Fibro, paste0(curr_dir, "/CellType/Finalcelltype/Fibro.RDS"))

Fibro <- RunUMAP(object = Fibro, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist = 1.4)
DimPlot(Fibro, raster=F)

#res.01
#0 Fib_1 SLC4A4, KCNMA1, TJP2,
#1 Fib_2 LAMA2, COLs, SLCs, TRPM3, CLDN, ROBO2
#2 Fib_3 LAMA2 COLs
#3 Fib_1 SLC4A4, KCNMA1, TJP2,

###Remove interdoublets

Idents(Fibro) <- Fibro@meta.data$Subcelltype
Fibro <- CellSelector(DimPlot(Fibro), Fibro, ident = "Fib1")
Fib_1 <- subset(Fibro, Subcelltype=="Fib_1")
Fib_1 <- subset(Fib_1, ident="Fib1")

Idents(Fibro) <- Fibro@meta.data$Subcelltype
Fibro <- CellSelector(DimPlot(Fibro), Fibro, ident = "Fib2")
Fib_2 <- subset(Fibro, Subcelltype=="Fib_2")
Fib_2 <- subset(Fib_2, ident="Fib2")

Idents(Fibro) <- Fibro@meta.data$Subcelltype
Fibro <- CellSelector(DimPlot(Fibro), Fibro, ident = "Fib3")
Fib_3 <- subset(Fibro, Subcelltype=="Fib_3")
Fib_3 <- subset(Fib_3, ident="Fib3")

Fib_1_meta <- Fib_1@meta.data
Fib_2_meta <- Fib_2@meta.data
Fib_3_meta <- Fib_3@meta.data
merged_meta <- rbind(Fib_1_meta, Fib_2_meta, Fib_3_meta)
rm(Fib_1, Fib_2, Fib_3)
gc()

Fibro@meta.data$Subcelltype <- NULL
Fibro_meta <- Fibro@meta.data

Fibro_meta <- rownames_to_column(Fibro_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Fibro_meta <- left_join(Fibro_meta, merged_meta, by="rowname")
Fibro_meta <- column_to_rownames(Fibro_meta, "rowname")
Fibro_meta$Subcelltype <- as.character(Fibro_meta$Subcelltype)
Fibro_meta[is.na(Fibro_meta)==TRUE] <- "Doublet"
Fibro@meta.data <- Fibro_meta
Fibro <- subset(Fibro, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Fibro) <- Fibro@meta.data[["Subcelltype"]] 
levels(Fibro) <- c("Fib_1", "Fib_2", "Fib_3")
Fibro@meta.data[["Subcelltype"]] <- Idents(Fibro)
DimPlot(Fibro, raster=F)

Fibro <- RunUMAP(object = Fibro, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist = 1.4)


##Imm
###PCA

ImmA <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/ImmA.RDS"))
ImmB <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/ImmB.RDS"))
Imm <- merge(ImmA, ImmB)
rm(ImmA, ImmB)
gc()

NFSP_fun(Imm, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Imm <- RunHarmony(object = Imm, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Imm, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony
Imm <- RunUMAP(Imm, dims=1:20, assay="soupXcounts", reduction="harmony", min.dist=1, n.neighbors=50L)
#Pick res
Idents(Imm) <- Imm@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Imm, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Imm, assay="soupXcounts", features = Imm_soupXcounts_var_feat)
saveRDS(Imm, paste0(curr_dir, "/CellType/Firstcelltype/Imm.RDS"))



curr_clust <- 5
Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Imm, curr_clust)

search_markersfun(Imm, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Imm, curr_clust, known_Fibro_prog)

search_markersfun(Imm, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Imm, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Imm, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Imm, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Imm, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Imm, curr_clust, c(known_light) %>% unique())
search_markersfun(Imm, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Imm, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Imm, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Imm, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Imm, curr_clust, known_dividing %>% unique())

search_markersfun(Imm, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Imm, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Imm, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Imm, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Imm, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Imm, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Imm, "TTR", raster=F)

###Label
Imm <- RenameIdents(Imm, "0"="BAM", "2"="BAM", "5"="BAM", "3"="T_cell", "1"="Doublet", "6"="Doublet", "4"="Doublet")
Imm@meta.data$Subcelltype <- Imm@active.ident
DoubletI <- subset(Imm, Subcelltype=="Doublet")
#remove doublets
#Imm <- subset(Imm, CMO_classification.global=="Singlet")
#Imm <- subset(Imm, percent.mt<5)

BAM <- subset(Imm, Subcelltype=="BAM")
T_cell <- subset(Imm, Subcelltype=="T_cell")

#saveRDS(Imm, paste0(curr_dir, "/CellType/Secondcelltype/Imm.RDS"))


###Venn of agreement

plts <- ggVennDiagram(list(colnames(Imm)[Imm@meta.data$DoubletFinder_call=="Doublet"], colnames(Imm)[Imm@meta.data$CMO_classification.global=="Doublet"], colnames(Imm)[Imm@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Imm)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Imm/Imm DoubletFinder-CMO-Manual venn.png", sep=""))


###remove doublets

Imm <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Imm.RDS"))
Imm <- subset(Ex_ChP, subset=Firstcelltype=="Imm")
Imm_all <- Imm

Featurefun(Imm, c("HSPH1"), "Imm_markers",  assay="soupXcounts")
#Remove cells that don't belong
Imm <- subset(Imm, HTR2C < 4)
Imm <- subset(Imm, TTR < 4)
Imm <- subset(Imm, CLIC6 < 4)
Imm <- subset(Imm, AQP1 < 3)
Imm <- subset(Imm, CFAP54 < 3)
Imm <- subset(Imm, FAM227A < 3)
Imm <- subset(Imm, ARL13B < 3)
Imm <- subset(Imm, TRPM3 < 3)
Imm <- subset(Imm, SLC4A4 < 3)
Imm <- subset(Imm, LAMA2 < 3)
Imm <- subset(Imm, LAMA4 < 3)
Imm <- subset(Imm, CLDN1 < 3)
Imm <- subset(Imm, CLDN11 < 3)
Imm <- subset(Imm, PRICKLE1 < 3)
Imm <- subset(Imm, CEMIP < 3)
Imm <- subset(Imm, FOXP2 < 3)
Imm <- subset(Imm, COL6A3 < 3)
Imm <- subset(Imm, SLC6A13 < 3)
Imm <- subset(Imm, MECOM < 3)
Imm <- subset(Imm, VWF < 3)
Imm <- subset(Imm, TAGLN < 3)
Imm <- subset(Imm, MYH11 < 3)
Imm <- subset(Imm, NOTCH3 < 3)
Imm <- subset(Imm, CADM2 < 3)
Imm <- subset(Imm, MAP2 < 3)
Imm <- subset(Imm, GRIK2 < 3)
Imm <- subset(Imm, CLU < 3)
Imm <- subset(Imm, AQP4 < 3)
Imm <- subset(Imm, MOBP < 3)
Imm <- subset(Imm, MBP < 3)
Imm <- subset(Imm, CDK19 < 3)
#Imm <- subset(Imm, CD68 < 3)
#Imm <- subset(Imm, F13A1 < 3)
#Imm <- subset(Imm, STAB1 < 3)
#Imm <- subset(Imm, CD163 < 3)
#Imm <- subset(Imm, LYVE1 < 3)
#Imm <- subset(Imm, P2RY12 < 3)
#Imm <- subset(Imm, THEMIS < 3)

DimPlot(Imm)

###Remove inter doublets

Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM1")
BAM1 <- subset(Imm, ident="BAM1")
BAM1 <- subset(BAM1, Subcelltype=="BAM_1")

Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM2")
BAM2 <- subset(Imm, ident="BAM2")
BAM2 <- subset(BAM2, Subcelltype=="BAM_2")

Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM3")
BAM3 <- subset(Imm, ident="BAM3")
BAM3 <- subset(BAM3, Subcelltype=="BAM_3")


Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "T_c")
T_cell <- subset(Imm, ident="T_c")
T_cell <- subset(T_cell, Subcelltype=="T_Cell")

BAM1_meta <- BAM1@meta.data
BAM2_meta <- BAM2@meta.data
BAM3_meta <- BAM3@meta.data
T_cell_meta <- T_cell@meta.data
merged_meta <- rbind(BAM1_meta, BAM2_meta, BAM3_meta, T_cell_meta)

Imm@meta.data$Subcelltype <- NULL
Imm_meta <- Imm@meta.data

Imm_meta <- rownames_to_column(Imm_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Imm_meta <- left_join(Imm_meta, merged_meta, by="rowname")
Imm_meta <- column_to_rownames(Imm_meta, "rowname")
Imm_meta$Subcelltype <- as.character(Imm_meta$Subcelltype)
Imm_meta[is.na(Imm_meta)==TRUE] <- "Doublet"
Imm@meta.data <- Imm_meta
Imm <- subset(Imm, Subcelltype=="Doublet", invert=T)



##T_cell
###PCA

NFSP_fun(T_cell, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
T_cell <- RunHarmony(object = T_cell, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(T_cell, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(T_cell) <- T_cell@meta.data[["soupXcounts_snn_res.0.09"]]
CNDRdimfun(T_cell, reduction = "harmony", res="_0.09", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(T_cell, assay="soupXcounts", features = T_cell_soupXcounts_var_feat)



curr_clust <- 0
T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature , known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(T_cell_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mgliaM1, known_mgliaM2, known_BAM,  known_BAM_MCHIIlo, known_BAM_MCHIIlo_dura, known_BAM_MCHIIhi, known_BAM_MCHIIhi_dura, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(T_cell, curr_clust)

search_markersfun(T_cell, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(T_cell, curr_clust, known_Fibro_prog)

search_markersfun(T_cell, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(T_cell, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(T_cell, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(T_cell, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(T_cell, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(T_cell, curr_clust, c(known_light) %>% unique())
search_markersfun(T_cell, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(T_cell, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(T_cell, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(T_cell, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(T_cell, curr_clust, known_dividing %>% unique())

search_markersfun(T_cell, curr_clust, c(known_ChPT_cellature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(T_cell, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(T_cell, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(T_cell, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mgliaM1, known_mgliaM2, known_BAM, known_BAM_MCHIIlo, known_BAM_MCHIIlo_dura, known_BAM_MCHIIhi, known_BAM_MCHIIhi_dura, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(T_cell, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(T_cell, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(T_cell, "TTR", raster=F)

###Label
T_cell <- RenameIdents(T_cell, "1"="T_Cell", "2"="T_Cell", "5"="T_Cell", "0"="Doublet", "3"="Doublet", "4"="Doublet")
T_cell@meta.data$Subcelltype <- T_cell@active.ident

T_cell <- subset(T_cell, Subcelltype=="T_Cell")

T_cell <- subset(T_cell, CMO_classification.global=="Singlet")
T_cell <- subset(T_cell, percent.mt<5)

saveRDS(T_cell, paste0(curr_dir, "/CellType/Secondcelltype/T_cell.RDS"))


###Venn of agreement

plts <- ggVennDiagram(list(colnames(T_cell)[T_cell@meta.data$DoubletFinder_call=="Doublet"], colnames(T_cell)[T_cell@meta.data$CMO_classification.global=="Doublet"], colnames(T_cell)[T_cell@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(T_cell)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/T_cell/T_cell DoubletFinder-CMO-Manual venn.png", sep=""))


###Remove doublets

T_cell_all <- T_cell

Featurefun(T_cell, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong
T_cell <- subset(T_cell, HTR2C < 4)
T_cell <- subset(T_cell, TTR < 4)
T_cell <- subset(T_cell, CLIC6 < 4)
#T_cell <- subset(T_cell, AQP1 < 3)
T_cell <- subset(T_cell, CFAP54 < 3)
T_cell <- subset(T_cell, FAM227A < 3)
T_cell <- subset(T_cell, ARL13B < 3)
#T_cell <- subset(T_cell, TRPM3 < 3)
#T_cell <- subset(T_cell, SLC4A4 < 3)
#T_cell <- subset(T_cell, LAMA2 < 3)
#T_cell <- subset(T_cell, LAMA4 < 3)
#T_cell <- subset(T_cell, CLDN1 < 3)
#T_cell <- subset(T_cell, CLDN11 < 3)
#T_cell <- subset(T_cell, PRICKLE1 < 3)
#T_cell <- subset(T_cell, CEMIP < 3)
#T_cell <- subset(T_cell, FOXP2 < 3)
#T_cell <- subset(T_cell, COL6A3 < 3)
#T_cell <- subset(T_cell, SLC6A13 < 3)
T_cell <- subset(T_cell, MECOM < 3)
T_cell <- subset(T_cell, VWF < 3)
T_cell <- subset(T_cell, TAGLN < 3)
T_cell <- subset(T_cell, MYH11 < 3)
T_cell <- subset(T_cell, NOTCH3 < 3)
T_cell <- subset(T_cell, CADM2 < 3)
T_cell <- subset(T_cell, MAP2 < 3)
T_cell <- subset(T_cell, GRIK2 < 3)
#T_cell <- subset(T_cell, CLU < 3)
T_cell <- subset(T_cell, AQP4 < 3)
T_cell <- subset(T_cell, MOBP < 3)
T_cell <- subset(T_cell, MBP < 3)
T_cell <- subset(T_cell, CDK19 < 3)
T_cell <- subset(T_cell, CD68 < 3)
T_cell <- subset(T_cell, F13A1 < 3)
T_cell <- subset(T_cell, STAB1 < 3)
T_cell <- subset(T_cell, CD163 < 3)
T_cell <- subset(T_cell, LYVE1 < 3)
T_cell <- subset(T_cell, P2RY12 < 3)
T_cell <- subset(T_cell, THEMIS < 3)

T_cell <- CellSelector(DimPlot(T_cell), T_cell, ident = "T_cell")
T_cell <- CellSelector(DimPlot(T_cell), T_cell, ident = "Doublet")
#T_cell@meta.data$T_celltemp <- Idents(T_cell)
T_cell <- subset(T_cell, ident="T_cell")
#T_cell@meta.data$T_celltemp <- NULL
label_doubletsfun(T_cell, T_cell_all)
T_cell <- T_cell_all


##BAM
###PCA

#BAM <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/BAM.RDS"))
NFSP_fun(BAM, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
BAM <- RunHarmony(object = BAM, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(BAM, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(BAM) <- BAM@meta.data$soupXcounts_snn_res.0.1
CNDRdimfun(BAM, reduction = "harmony", res="_Final", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global", "apoe_genotype", "Cog_Path"))

#return markers
return_markersfun(BAM, assay="soupXcounts", features = BAM_soupXcounts_var_feat)

Ncellfun(BAM, paste0(curr_dir, "/CellType/BAM/Final per APOE.csv", sep=""), "apoe_genotype")
Ncellfun(BAM, paste0(curr_dir, "/CellType/BAM/Final per AD.csv", sep=""), "Cog_Path")
Ncellfun(BAM, paste0(curr_dir, "/CellType/BAM/Final per projid.csv", sep=""), "projid")




curr_clust <- 4
BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature , known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(BAM_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(BAM, curr_clust)

search_markersfun(BAM, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(BAM, curr_clust, known_Fibro_prog)

search_markersfun(BAM, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(BAM, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(BAM, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(BAM, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(BAM, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(BAM, curr_clust, c(known_light) %>% unique())
search_markersfun(BAM, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(BAM, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(BAM, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(BAM, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(BAM, curr_clust, known_dividing %>% unique())

search_markersfun(BAM, curr_clust, c(known_ChPBAMature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(BAM, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(BAM, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(BAM, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(BAM, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(BAM, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(BAM, c("ABCA1", "AOAH", "CYBB"), raster=F)

#0: BAM activated
#1: BAM activated
#2: BAM activated
#3: BAM hHSP
#4: BAM hHSP
#5: BAM activated
#6: BAM hHSP
#7: BAM hHSP
#8: BAM activated
#9: BAM hHSP-dividing
#10: BAM activated

###Label

#First pass
BAM <- RenameIdents(BAM, "0"="Doublet", "1"="BAM", "2"="Doublet?", "3"="Doublet", "4"="BAM", "5"="BAM", "6"="BAM", "7"="BAM", "8"="BAM", "9"="BAM", "10"="Doublet") #resw=0.1
BAM@meta.data$Subcelltype <- BAM@active.ident

BAM@meta.data$Subcelltype <- BAM@active.ident
DoubletM <- subset(BAM, Subcelltype=="Doublet")
BAM <- subset(BAM, Subcelltype=="BAM")
#remove doublets
BAM <- subset(BAM, CMO_classification.global=="Singlet")
BAM <- subset(BAM, percent.mt<5)

#Second pass
BAM <- RenameIdents(BAM, "0"="BAM_1", "1"="BAM_2", "2"="BAM_3") #res=0.03
BAM@meta.data$Subcelltype <- BAM@active.ident
BAM <- subset(BAM, c(Subcelltype=="BAM_1" | Subcelltype=="BAM_1" | Subcelltype=="BAM_1" | Subcelltype=="BAM_2" | Subcelltype=="BAM_3"))
levels(BAM) <- c("BAM_1", "BAM_1", "BAM_1", "BAM_2", "BAM_3")
BAM@meta.data[["Subcelltype"]] <- Idents(BAM)

saveRDS(BAM, paste0(curr_dir, "/CellType/Secondcelltype/BAM.RDS"))


###Venn of agreement

plts <- ggVennDiagram(list(colnames(BAM)[BAM@meta.data$DoubletFinder_call=="Doublet"], colnames(BAM)[BAM@meta.data$CMO_classification.global=="Doublet"], colnames(BAM)[BAM@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(BAM)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/BAM/BAM DoubletFinder-CMO-Manual venn.png", sep=""))


###remove doublets

BAM_all <- BAM

Featurefun(BAM, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong
BAM <- subset(BAM, HTR2C < 4)
BAM <- subset(BAM, TTR < 4)
BAM <- subset(BAM, CLIC6 < 4)
#BAM <- subset(BAM, AQP1 < 3)
BAM <- subset(BAM, CFAP54 < 3)
BAM <- subset(BAM, FAM227A < 3)
BAM <- subset(BAM, ARL13B < 3)
#BAM <- subset(BAM, TRPM3 < 3)
#BAM <- subset(BAM, SLC4A4 < 3)
#BAM <- subset(BAM, LAMA2 < 3)
#BAM <- subset(BAM, LAMA4 < 3)
#BAM <- subset(BAM, CLDN1 < 3)
#BAM <- subset(BAM, CLDN11 < 3)
#BAM <- subset(BAM, PRICKLE1 < 3)
#BAM <- subset(BAM, CEMIP < 3)
#BAM <- subset(BAM, FOXP2 < 3)
#BAM <- subset(BAM, COL6A3 < 3)
#BAM <- subset(BAM, SLC6A13 < 3)
BAM <- subset(BAM, MECOM < 3)
BAM <- subset(BAM, VWF < 3)
BAM <- subset(BAM, TAGLN < 3)
BAM <- subset(BAM, MYH11 < 3)
BAM <- subset(BAM, NOTCH3 < 3)
BAM <- subset(BAM, CADM2 < 3)
BAM <- subset(BAM, MAP2 < 3)
BAM <- subset(BAM, GRIK2 < 3)
#BAM <- subset(BAM, CLU < 3)
#BAM <- subset(BAM, AQP4 < 3)
BAM <- subset(BAM, MOBP < 3)
BAM <- subset(BAM, MBP < 3)
BAM <- subset(BAM, CDK19 < 3)
BAM <- subset(BAM, CD68 < 3)
BAM <- subset(BAM, F13A1 < 3)
BAM <- subset(BAM, STAB1 < 3)
BAM <- subset(BAM, CD163 < 3)
BAM <- subset(BAM, LYVE1 < 3)
BAM <- subset(BAM, P2RY12 < 3)
BAM <- subset(BAM, THEMIS < 3)

Idents(BAM) <- BAM@meta.data$Subcelltype
BAM1 <- subset(BAM, Subcelltype=="BAM_1")

Idents(BAM) <- BAM@meta.data$Subcelltype
BAM <- CellSelector(DimPlot(BAM), BAM, ident = "BAM2a")
BAM2a <- subset(BAM, ident="BAM2a")
BAM2a <- subset(BAM2a, Subcelltype=="BAM_1")

Idents(BAM) <- BAM@meta.data$Subcelltype
BAM <- CellSelector(DimPlot(BAM), BAM, ident = "BAM2b")
BAM2b <- subset(BAM, ident="BAM2b")
BAM2b <- subset(BAM2b, Subcelltype=="BAM_1")

Idents(BAM) <- BAM@meta.data$Subcelltype
BAM <- CellSelector(DimPlot(BAM), BAM, ident = "BAM3")
BAM3 <- subset(BAM, ident="BAM3")
BAM3 <- subset(BAM3, Subcelltype=="BAM_2")

Idents(BAM) <- BAM@meta.data$Subcelltype
BAM <- CellSelector(DimPlot(BAM), BAM, ident = "BAM4")
BAM4 <- subset(BAM, ident="BAM4")
BAM4 <- subset(BAM4, Subcelltype=="BAM_3")

BAM0_meta <- BAM0@meta.data
BAM1_meta <- BAM1@meta.data
BAM2a_meta <- BAM2a@meta.data
BAM2b_meta <- BAM2b@meta.data
BAM3_meta <- BAM3@meta.data
BAM4_meta <- BAM4@meta.data
merged_meta <- rbind(BAM0_meta, BAM1_meta, BAM2a_meta, BAM2b_meta, BAM3_meta, BAM4_meta)

BAM@meta.data$Subcelltype <- NULL
BAM_meta <- BAM@meta.data

BAM_meta <- rownames_to_column(BAM_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
BAM_meta <- left_join(BAM_meta, merged_meta, by="rowname")
BAM_meta <- column_to_rownames(BAM_meta, "rowname")
BAM_meta$Subcelltype <- as.character(BAM_meta$Subcelltype)
BAM_meta[is.na(BAM_meta)==TRUE] <- "Doublet"
BAM@meta.data <- BAM_meta
BAM <- subset(BAM, Subcelltype=="Doublet", invert=T)

DimPlot(BAM, raster=F)

BAM <- CellSelector(DimPlot(BAM), BAM, ident = "Doublet")
#BAM@meta.data$BAMtemp <- Idents(BAM)
BAM <- subset(BAM, ident="BAM")
#BAM@meta.data$BAMtemp <- NULL
label_doubletsfun(BAM, BAM_all)
BAM <- BAM_all


###Imm

Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- RenameIdents(Imm, "BAM_1"="0", "BAM_2"="1", "BAM_3"="4", "T_Cell"="5")
return_markersfun(Imm, assay="soupXcounts", features = Imm_soupXcounts_var_feat)

unk_BAM_Phag <- c(Imm_markers_list[["0"]][[2]][[1]][["gene"]][0:100])
unk_BAM_2 <- c(Imm_markers_list[["1"]][[2]][[1]][["gene"]][0:100])
unk_BAM_1 <- c(Imm_markers_list[["3"]][[2]][[1]][["gene"]][0:100])
unk_BAM_Div <- c(Imm_markers_list[["4"]][[2]][[1]][["gene"]][0:100])
unk_T_cell <- c(Imm_markers_list[["5"]][[2]][[1]][["gene"]][0:100])

Idents(Imm) <- Imm@meta.data$Subcelltype
Featurefun(Imm, unk_BAM_Phag, "BAM_Phag", assay="soupXcounts")
Featurefun(Imm, unk_BAM_2, "BAM_2", assay="soupXcounts")
Featurefun(Imm, unk_BAM_Div, "BAM_Dividing", assay="soupXcounts")
Featurefun(Imm, unk_T_cell, "BAM_T_Cell", assay="soupXcounts")


##Neural

Neural <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Neural.RDS"))
Neural <- subset(Ex_ChP, subset=Firstcelltype=="Neural")
Neural_all <- Neural

Featurefun(Neural, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong

#Remove cells that don't belong
Neuro <- subset(Neuro, HTR2C < 3)
Neuro <- subset(Neuro, TTR < 3)
Neuro <- subset(Neuro, CLIC6 < 3)
#Neuro <- subset(Neuro, AQP1 < 3)
Neuro <- subset(Neuro, CFAP54 < 3)
Neuro <- subset(Neuro, FAM227A < 3)
Neuro <- subset(Neuro, ARL13B < 3)
Neuro <- subset(Neuro, TRPM3 < 3)
#Neuro <- subset(Neuro, SLC4A4 < 3)
Neuro <- subset(Neuro, LAMA2 < 3)
Neuro <- subset(Neuro, LAMA4 < 3)
Neuro <- subset(Neuro, CLDN1 < 3)
Neuro <- subset(Neuro, CLDN11 < 3)
Neuro <- subset(Neuro, PRICKLE1 < 3)
Neuro <- subset(Neuro, CEMIP < 3)
Neuro <- subset(Neuro, FOXP2 < 3)
Neuro <- subset(Neuro, COL6A3 < 3)
#Neuro <- subset(Neuro, SLC6A13 < 3)
Neuro <- subset(Neuro, MECOM < 3)
Neuro <- subset(Neuro, VWF < 3)
Neuro <- subset(Neuro, TAGLN < 3)
Neuro <- subset(Neuro, MYH11 < 3)
Neuro <- subset(Neuro, NOTCH3 < 3)
#Neuro <- subset(Neuro, CADM2 < 3)
#Neuro <- subset(Neuro, MAP2 < 3)
#Neuro <- subset(Neuro, GRIK2 < 3)
#Neuro <- subset(Neuro, CLU < 3)
#Neuro <- subset(Neuro, AQP4 < 3)
#Neuro <- subset(Neuro, MOBP < 3)
#Neuro <- subset(Neuro, MBP < 3)
#Neuro <- subset(Neuro, CDK19 < 3)
Neuro <- subset(Neuro, CD68 < 3)
Neuro <- subset(Neuro, F13A1 < 3)
Neuro <- subset(Neuro, STAB1 < 3)
Neuro <- subset(Neuro, CD163 < 3)
Neuro <- subset(Neuro, LYVE1 < 3)
Neuro <- subset(Neuro, P2RY12 < 3)
Neuro <- subset(Neuro, THEMIS < 3)

Neural <- CellSelector(DimPlot(Neural), Neural, ident = "Oligo")
Neural <- CellSelector(DimPlot(Neural), Neural, ident = "Neuron")
Neural <- CellSelector(DimPlot(Neural), Neural, ident = "Astro")
Neural <- CellSelector(DimPlot(Neural), Neural, ident = "Doublet")
Neural@meta.data$Doublet <- Neural@active.ident
Neural <- subset(Neural, ident="Neural")
label_doubletsfun(Neural, Neural_all)
Neural <- Neural_all


###PCA

NFSP_fun(Neural, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Neural <- RunHarmony(object = Neural, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Neural, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Neural) <- Neural@meta.data[["soupXcounts_snn_res.0.05"]]
CNDRdimfun(Neural, reduction = "harmony", res="_0.05", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global", "Doublet"))
#return markers
return_markersfun(Neural, assay="soupXcounts", features = Neural_soupXcounts_var_feat)
#saveRDS(Neural, paste0(curr_dir, "/CellType/Firstcelltype/Neural.RDS"))



curr_clust <- 5
Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

Top20_Featuresfun(Neural, curr_clust)

search_markersfun(Neural, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Neural, curr_clust, known_Fibro_prog)

search_markersfun(Neural, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Neural, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Neural, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Neural, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Neural, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Neural, curr_clust, c(known_light) %>% unique())
search_markersfun(Neural, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Neural, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Neural, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Neural, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Neural, curr_clust, known_dividing %>% unique())

search_markersfun(Neural, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Neural, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Neural, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Neural, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Neural, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Neural, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Neural, "TTR", raster=F)

###Label
Neural <- RenameIdents(Neural, "0"="Oligo", "3"="Neuron", "4"="Astro", "1"="Doublet", "2"="Doublet")
Neural@meta.data$Subcelltype <- Neural@active.ident
Neuro <- subset(Neural, Subcelltype=="Neuron")
Neuro <- subset(Neuro, Doublet=="Neuron")
Oligo <- subset(Neural, Subcelltype=="Oligo")
Oligo <- subset(Oligo, Doublet=="Oligo")
Astro <- subset(Neural, Subcelltype=="Astro")
Astro <- subset(Astro, Doublet=="Astro")

DoubletN <- subset(Neural, Subcelltype=="Doublet")


###Venn of agreement

plts <- ggVennDiagram(list(colnames(Neural)[Neural@meta.data$DoubletFinder_call=="Doublet"], colnames(Neural)[Neural@meta.data$CMO_classification.global=="Doublet"], colnames(Neural)[Neural@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Neural)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Neural/Neural DoubletFinder-CMO-Manual venn.png", sep=""))



###Neuro
Neuro <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Neuro.RDS"))
Neuro_all <- Neuro

Featurefun(Neuro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong
Neuro <- subset(Neuro, HTR2C < 4)
Neuro <- subset(Neuro, TTR < 4)
Neuro <- subset(Neuro, CLIC6 < 4)
#Neuro <- subset(Neuro, AQP1 < 3)
Neuro <- subset(Neuro, CFAP54 < 3)
Neuro <- subset(Neuro, FAM227A < 3)
Neuro <- subset(Neuro, ARL13B < 3)
#Neuro <- subset(Neuro, TRPM3 < 3)
#Neuro <- subset(Neuro, SLC4A4 < 3)
#Neuro <- subset(Neuro, LAMA2 < 3)
#Neuro <- subset(Neuro, LAMA4 < 3)
#Neuro <- subset(Neuro, CLDN1 < 3)
#Neuro <- subset(Neuro, CLDN11 < 3)
#Neuro <- subset(Neuro, PRICKLE1 < 3)
#Neuro <- subset(Neuro, CEMIP < 3)
#Neuro <- subset(Neuro, FOXP2 < 3)
#Neuro <- subset(Neuro, COL6A3 < 3)
#Neuro <- subset(Neuro, SLC6A13 < 3)
Neuro <- subset(Neuro, MECOM < 3)
Neuro <- subset(Neuro, VWF < 3)
Neuro <- subset(Neuro, TAGLN < 3)
Neuro <- subset(Neuro, MYH11 < 3)
Neuro <- subset(Neuro, NOTCH3 < 3)
Neuro <- subset(Neuro, CADM2 < 3)
Neuro <- subset(Neuro, MAP2 < 3)
Neuro <- subset(Neuro, GRIK2 < 3)
#Neuro <- subset(Neuro, CLU < 3)
Neuro <- subset(Neuro, AQP4 < 3)
Neuro <- subset(Neuro, MOBP < 3)
Neuro <- subset(Neuro, MBP < 3)
Neuro <- subset(Neuro, CDK19 < 3)
Neuro <- subset(Neuro, CD68 < 3)
Neuro <- subset(Neuro, F13A1 < 3)
Neuro <- subset(Neuro, STAB1 < 3)
Neuro <- subset(Neuro, CD163 < 3)
Neuro <- subset(Neuro, LYVE1 < 3)
Neuro <- subset(Neuro, P2RY12 < 3)
Neuro <- subset(Neuro, THEMIS < 3)

Neuro <- CellSelector(DimPlot(Neuro), Neuro, ident = "Neuro")
Neuro@meta.data$Doublet <- Neuro@active.ident
Neuro <- subset(Neuro, ident="Neuro")
#Neuro@meta.data$Immtemp <- NULL
label_doubletsfun(Neuro, Neuro_all)
Neuro <- Neuro_all


####PCA
NeuroA <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/NeuroA.RDS"))
NeuroB <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/NeuroB.RDS"))
Neuro <- merge(NeuroA, NeuroB)
rm(NeuroA, NeuroB)
gc()

NFSP_fun(Neuro, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Neuro <- RunHarmony(object = Neuro, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Neuro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Neuro) <- Neuro@meta.data[["soupXcounts_snn_res.0.03"]]
CNDRdimfun(Neuro, reduction = "harmony", res="_0.03", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Neuro, assay="soupXcounts", features = Neuro_soupXcounts_var_feat)
saveRDS(Neuro, paste0(curr_dir, "/CellType/Firstcelltype/Neuro.RDS"))



curr_clust <- 0
Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Neuro, curr_clust)

search_markersfun(Neuro, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Neuro, curr_clust, known_Fibro_prog)

search_markersfun(Neuro, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Neuro, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_light) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Neuro, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Neuro, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Neuro, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Neuro, curr_clust, known_dividing %>% unique())

search_markersfun(Neuro, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Neuro, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Neuro, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Neuro, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Neuro, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Neuro, "TTR", raster=F)


####Label

Neuro <- RenameIdents(Neuro, "1"="Neuro", "2"="Neuro", "3"="Neuro", "4"="Neuro", "5"="Neuro", "7"="Neuro", "8"="Neuro", "0"="Doublet", "6"="Doublet")
Neuro@meta.data$Subcelltype <- Neuro@active.ident
DoubletN <- subset(Neuro, Subcelltype=="Doublet")
Neuro <- subset(Neuro, Subcelltype=="Neuro")
#remove HTO doublets and %MTO
Neuro <- subset(Neuro, CMO_classification.global=="Singlet")
Neuro <- subset(Neuro, percent.mt<5)

saveRDS(Neuro, paste0(curr_dir, "/CellType/Secondcelltype/Neuro.RDS"))


####Venn of agreement

plts <- ggVennDiagram(list(colnames(Neuro)[Neuro@meta.data$DoubletFinder_call=="Doublet"], colnames(Neuro)[Neuro@meta.data$CMO_classification.global=="Doublet"], colnames(Neuro)[Neuro@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Neuro)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Neuro/Neuro DoubletFinder-CMO-Manual venn.png", sep=""))



###Astro
Astro <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Astro.RDS"))
Astro_all <- Astro

Featurefun(Astro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong

#Remove cells that don't belong
Astro <- subset(Astro, HTR2C < 3)
Astro <- subset(Astro, TTR < 3)
Astro <- subset(Astro, CLIC6 < 3)
#Astro <- subset(Astro, AQP1 < 3)
Astro <- subset(Astro, CFAP54 < 3)
Astro <- subset(Astro, FAM227A < 3)
Astro <- subset(Astro, ARL13B < 3)
Astro <- subset(Astro, TRPM3 < 3)
#Astro <- subset(Astro, SLC4A4 < 3)
Astro <- subset(Astro, LAMA2 < 3)
Astro <- subset(Astro, LAMA4 < 3)
Astro <- subset(Astro, CLDN1 < 3)
Astro <- subset(Astro, CLDN11 < 3)
Astro <- subset(Astro, PRICKLE1 < 3)
Astro <- subset(Astro, CEMIP < 3)
Astro <- subset(Astro, FOXP2 < 3)
Astro <- subset(Astro, COL6A3 < 3)
#Astro <- subset(Astro, SLC6A13 < 3)
Astro <- subset(Astro, MECOM < 3)
Astro <- subset(Astro, VWF < 3)
Astro <- subset(Astro, TAGLN < 3)
Astro <- subset(Astro, MYH11 < 3)
Astro <- subset(Astro, NOTCH3 < 3)
#Astro <- subset(Astro, CADM2 < 3)
#Astro <- subset(Astro, MAP2 < 3)
#Astro <- subset(Astro, GRIK2 < 3)
#Astro <- subset(Astro, CLU < 3)
#Astro <- subset(Astro, AQP4 < 3)
#Astro <- subset(Astro, MOBP < 3)
#Astro <- subset(Astro, MBP < 3)
#Astro <- subset(Astro, CDK19 < 3)
Astro <- subset(Astro, CD68 < 3)
Astro <- subset(Astro, F13A1 < 3)
Astro <- subset(Astro, STAB1 < 3)
Astro <- subset(Astro, CD163 < 3)
Astro <- subset(Astro, LYVE1 < 3)
Astro <- subset(Astro, P2RY12 < 3)
Astro <- subset(Astro, THEMIS < 3)

Astro <- CellSelector(DimPlot(Astro), Astro, ident = "Astro")
Astro@meta.data$Doublet <- Astro@active.ident
#Astro@meta.data$Immtemp <- Idents(Astro)
Astro <- subset(Astro, ident="Astro")
#Astro@meta.data$Immtemp <- NULL
label_doubletsfun(Astro, Astro_all)
Astro <- Astro_all


####PCA
AstroA <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/AstroA.RDS"))
AstroB <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/AstroB.RDS"))
Astro <- merge(AstroA, AstroB)
rm(AstroA, AstroB)
gc()

NFSP_fun(Astro, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Astro <- RunHarmony(object = Astro, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Astro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Astro) <- Astro@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Astro, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Astro, assay="soupXcounts", features = Astro_soupXcounts_var_feat)
saveRDS(Astro, paste0(curr_dir, "/CellType/Firstcelltype/Astro.RDS"))



curr_clust <- 2
Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Astro, curr_clust)

search_markersfun(Astro, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Astro, curr_clust, known_Fibro_prog)

search_markersfun(Astro, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Astro, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Astro, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Astro, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Astro, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Astro, curr_clust, c(known_light) %>% unique())
search_markersfun(Astro, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Astro, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Astro, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Astro, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Astro, curr_clust, known_dividing %>% unique())

search_markersfun(Astro, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Astro, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Astro, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Astro, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Astro, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Astro, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Astro, "VCAN", raster=F)

0: Astro

####Label

Astro <- RenameIdents(Astro, "2"="Astro", "0"="Doublet", "1"="Doublet", "3"="Doublet", "4"="Doublet", "5"="Doublet", "6"="Doublet")
Astro@meta.data$Subcelltype <- Astro@active.ident
DoubletA <- subset(Astro, Subcelltype=="Doublet")
Astro <- subset(Astro, Subcelltype=="Astro")
#remove doublets
Astro <- subset(Astro, CMO_classification.global=="Singlet")
Astro <- subset(Astro, percent.mt<5)

saveRDS(Astro, paste0(curr_dir, "/CellType/Secondcelltype/Astro.RDS"))


####Venn of agreement

plts <- ggVennDiagram(list(colnames(Astro)[Astro@meta.data$DoubletFinder_call=="Doublet"], colnames(Astro)[Astro@meta.data$CMO_classification.global=="Doublet"], colnames(Astro)[Astro@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Astro)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Astro/Astro DoubletFinder-CMO-Manual venn.png", sep=""))



###Oligo

Oligo <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Oligo.RDS"))
Oligo_all <- Oligo

Featurefun(Oligo, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong

#Remove cells that don't belong
Oligo <- subset(Oligo, HTR2C < 3)
Oligo <- subset(Oligo, TTR < 3)
Oligo <- subset(Oligo, CLIC6 < 3)
#Oligo <- subset(Oligo, AQP1 < 3)
Oligo <- subset(Oligo, CFAP54 < 3)
Oligo <- subset(Oligo, FAM227A < 3)
Oligo <- subset(Oligo, ARL13B < 3)
Oligo <- subset(Oligo, TRPM3 < 3)
#Oligo <- subset(Oligo, SLC4A4 < 3)
Oligo <- subset(Oligo, LAMA2 < 3)
Oligo <- subset(Oligo, LAMA4 < 3)
Oligo <- subset(Oligo, CLDN1 < 3)
Oligo <- subset(Oligo, CLDN11 < 3)
Oligo <- subset(Oligo, PRICKLE1 < 3)
Oligo <- subset(Oligo, CEMIP < 3)
Oligo <- subset(Oligo, FOXP2 < 3)
Oligo <- subset(Oligo, COL6A3 < 3)
#Oligo <- subset(Oligo, SLC6A13 < 3)
Oligo <- subset(Oligo, MECOM < 3)
Oligo <- subset(Oligo, VWF < 3)
Oligo <- subset(Oligo, TAGLN < 3)
Oligo <- subset(Oligo, MYH11 < 3)
Oligo <- subset(Oligo, NOTCH3 < 3)
#Oligo <- subset(Oligo, CADM2 < 3)
#Oligo <- subset(Oligo, MAP2 < 3)
#Oligo <- subset(Oligo, GRIK2 < 3)
#Oligo <- subset(Oligo, CLU < 3)
#Oligo <- subset(Oligo, AQP4 < 3)
#Oligo <- subset(Oligo, MOBP < 3)
#Oligo <- subset(Oligo, MBP < 3)
#Oligo <- subset(Oligo, CDK19 < 3)
Oligo <- subset(Oligo, CD68 < 3)
Oligo <- subset(Oligo, F13A1 < 3)
Oligo <- subset(Oligo, STAB1 < 3)
Oligo <- subset(Oligo, CD163 < 3)
Oligo <- subset(Oligo, LYVE1 < 3)
Oligo <- subset(Oligo, P2RY12 < 3)
Oligo <- subset(Oligo, THEMIS < 3)

Oligo <- CellSelector(DimPlot(Oligo), Oligo, ident = "Olig")
Oligo@meta.data$Doublet <- Oligo@active.ident
#Oligo@meta.data$Immtemp <- Idents(Oligo)
Oligo <- subset(Oligo, ident="Olig")
#Oligo@meta.data$Immtemp <- NULL
label_doubletsfun(Oligo, Oligo_all)
Oligo <- Oligo_all


####PCA

OligoA <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/OligoA.RDS"))
OligoB <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/OligoB.RDS"))
Oligo <- merge(OligoA, OligoB)
rm(OligoA, OligoB)
gc()

NFSP_fun(Oligo, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Oligo <- RunHarmony(object = Oligo, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Oligo, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Oligo) <- Oligo@meta.data[["soupXcounts_snn_res.0.05"]]
CNDRdimfun(Oligo, reduction = "harmony", res="_0.05", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global", "Doublet"))
#return markers
return_markersfun(Oligo, assay="soupXcounts", features = Oligo_soupXcounts_var_feat)
saveRDS(Oligo, paste0(curr_dir, "/CellType/Firstcelltype/Oligo.RDS"))



curr_clust <- 0
Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Oligo, curr_clust)

search_markersfun(Oligo, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Oligo, curr_clust, known_Fibro_prog)

search_markersfun(Oligo, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Oligo, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_light) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Oligo, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Oligo, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Oligo, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Oligo, curr_clust, known_dividing %>% unique())

search_markersfun(Oligo, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Oligo, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Oligo, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Oligo, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Oligo, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Oligo, "TTR", raster=F)

####Label
Oligo <- RenameIdents(Oligo, "0"="Doublet", "1"="Oligo", "2"="Doublet")
Oligo@meta.data$Subcelltype <- Oligo@active.ident
DoubletO <- subset(Oligo, Subcelltype=="Doublet")
Oligo <- subset(Oligo, Subcelltype=="Oligo")
#remove doublets
Oligo <- subset(Oligo, CMO_classification.global=="Singlet")
Oligo <- subset(Oligo, percent.mt<5)

saveRDS(Oligo, paste0(curr_dir, "/CellType/Secondcelltype/Oligo.RDS"))


####Venn of agreement

plts <- ggVennDiagram(list(colnames(Oligo)[Oligo@meta.data$DoubletFinder_call=="Doublet"], colnames(Oligo)[Oligo@meta.data$CMO_classification.global=="Doublet"], colnames(Oligo)[Oligo@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Oligo)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Oligo/Oligo DoubletFinder-CMO-Manual venn.png", sep=""))




##Endo
###PCA
EndoA <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/EndoA.RDS"))
EndoB <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/EndoB.RDS"))
Endo <- merge(EndoA, EndoB)
rm(EndoA, EndoB)
gc()

NFSP_fun(Endo, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Endo <- RunHarmony(object = Endo, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Endo, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Endo) <- Endo@meta.data[["soupXcounts_snn_res.0.07"]]
CNDRdimfun(Endo, reduction = "harmony", res="_0.07", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Endo, assay="soupXcounts", features = Endo_soupXcounts_var_feat)
saveRDS(Endo, paste0(curr_dir, "/CellType/Firstcelltype/Endo.RDS"))

Featurefun(Endo, c("HTR2C", "TTR", "CLIC6", "AQP1", "TUBB4B", "FOXJ1", "ARL13B", "KRT18", "COX2", "TRPM3", "VIM", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "FOXP2", "COL6A3", "SLC6A13", "LUM", "DCN", "FN1", "ACTA2", "KCNMA1", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "CD3E", "HSPE1", "HSPH1"), "Key Markers",  assay="soupXcounts")


curr_clust <- 2
Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Endo, curr_clust)

search_markersfun(Endo, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
#search_markersfun(Endo, curr_clust, known_Fibro_prog)

search_markersfun(Endo, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Endo, curr_clust, c(known_veinous) %>% unique())
search_markersfun(Endo, curr_clust, c(known_capilary) %>% unique())
search_markersfun(Endo, curr_clust, c(known_Endo_1) %>% unique())
search_markersfun(Endo, curr_clust, c(known_lymphatic) %>% unique())

search_markersfun(Endo, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())


search_markersfun(Endo, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Endo, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Endo, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Endo, curr_clust, c(known_light) %>% unique())
search_markersfun(Endo, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Endo, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Endo, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Endo, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Endo, curr_clust, known_dividing %>% unique())

search_markersfun(Endo, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Endo, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Endo, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Endo, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Endo, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Endo, curr_clust, c(Cell_Stress, whole_cells) %>% unique())

search_markersfun(Endo, 5, c("CLDN5", "IL1R1", "VEGFC", "RGCC", "PLVAP", "SGK1"))
FeaturePlot(Endo, c("CLDN5", "IL1R1"), raster=F)
FeaturePlot(Endo, c("VEGFC", "RGCC"), raster=F)
FeaturePlot(Endo, c("PLVAP", "SGK1"), raster=F)

VlnPlot(Endo, c("CLDN5", "IL1R1", "VEGFC", "RGCC", "PLVAP", "SGK1"))

#last pass
#0 Endo_3
#1 Endo_1
#2 Endo_2
#3 Endo_2
#4 Endo_1

###Label

#Endo <- RenameIdents(Endo, "0"="Doublet", "1"="Doublet", "2"="Endo", "4"="Endo","5"="Endo", "6"="Endo", "9"="Endo", "10"="Endo", "3"="Doublet", "7"="Doublet", "8"="Doublet") #First pass
Endo <- RenameIdents(Endo, "0"="Endo", "1"="Doublet", "2"="Mural") #second pass
Endo@meta.data$Subcelltype <- Endo@active.ident
DoubletE <- subset(Endo, Subcelltype=="Doublet")
Mural <- subset(Endo, Subcelltype=="Mural")
Endo <- subset(Endo, Subcelltype=="Endo")
#remove doublets
Endo <- subset(Endo, CMO_classification.global=="Singlet")
Endo <- subset(Endo, percent.mt<5)

saveRDS(Endo, paste0(curr_dir, "/CellType/Secondcelltype/Endo.RDS"))



###Venn of agreement

plts <- ggVennDiagram(list(colnames(Endo)[Endo@meta.data$DoubletFinder_call=="Doublet"], colnames(Endo)[Endo@meta.data$CMO_classification.global=="Doublet"], colnames(Endo)[Endo@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Endo)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Endo/Endo DoubletFinder-CMO-Manual venn.png", sep=""))


###Remove interdoublets

Endo <- RenameIdents(Endo, "0"="Endo_3", "1"="Endo_1", "2"="Endo_2") #Last pass
Endo@meta.data$Subcelltype <- Endo@active.ident

Idents(Endo) <- Endo@meta.data$Subcelltype
Endo <- CellSelector(DimPlot(Endo), Endo, ident = "Cap")
Cap <- subset(Endo, ident="Cap")
Cap <- subset(Cap, Subcelltype=="Endo_3")


Idents(Endo) <- Endo@meta.data$Subcelltype
Endo <- CellSelector(DimPlot(Endo), Endo, ident = "Art")
Art <- subset(Endo, Subcelltype=="Endo_1")
Art <- subset(Art, ident="Art")


Idents(Endo) <- Endo@meta.data$Subcelltype
Endo <- CellSelector(DimPlot(Endo), Endo, ident = "Ven")
Ven <- subset(Endo, Subcelltype=="Endo_2")
Ven <- subset(Ven, ident="Ven")


Cap_meta <- Cap@meta.data
Art_meta <- Art@meta.data
Ven_meta <- Ven@meta.data
merged_meta <- rbind(Cap_meta, Art_meta, Ven_meta)

Endo@meta.data$Subcelltype <- NULL
Endo_meta <- Endo@meta.data

Endo_meta <- rownames_to_column(Endo_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Endo_meta <- left_join(Endo_meta, merged_meta, by="rowname")
Endo_meta <- column_to_rownames(Endo_meta, "rowname")
Endo_meta$Subcelltype <- as.character(Endo_meta$Subcelltype)
Endo_meta[is.na(Endo_meta)==TRUE] <- "Doublet"
Endo@meta.data <- Endo_meta
Endo <- subset(Endo, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Endo) <- Endo@meta.data[["Subcelltype"]] 
levels(Endo) <- c("Endo_1", "Endo_3", "Endo_2")
Endo@meta.data[["Subcelltype"]] <- Idents(Endo)
DimPlot(Endo, raster=F)

Endo <- RunUMAP(object = Endo, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist = 1.2)
saveRDS(Endo, paste0(curr_dir, "/CellType/Finalcelltype/Endo.RDS"))

Ncellfun(Endo, paste0(curr_dir, "/CellType/Endo/Final per APOE.csv", sep=""), "apoe_genotype")
Ncellfun(Endo, paste0(curr_dir, "/CellType/Endo/Final per AD.csv", sep=""), "Cog_Path")
Ncellfun(Endo, paste0(curr_dir, "/CellType/Endo/Final per projid.csv", sep=""), "projid")





##Mural
###PCA

NFSP_fun(Mural, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Mural <- RunHarmony(object = Mural, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Mural, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Mural) <- Mural@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Mural, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global", "Cog_Path"))
#return markers
return_markersfun(Mural, assay="soupXcounts", features = Mural_soupXcounts_var_feat)



curr_clust <- 2
Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Mural, curr_clust)

search_markersfun(Mural, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Mural, curr_clust, known_Fibro_prog)

search_markersfun(Mural, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Mural, curr_clust, c(known_Mural) %>% unique())
search_markersfun(Mural, curr_clust, c(known_peric) %>% unique())
search_markersfun(Mural, curr_clust, c(known_Mural_1) %>% unique())


search_markersfun(Mural, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Mural, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Mural, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Mural, curr_clust, c(known_light) %>% unique())
search_markersfun(Mural, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Mural, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Mural, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Mural, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Mural, curr_clust, known_dividing %>% unique())

search_markersfun(Mural, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Mural, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Mural, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Mural, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Mural, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Mural, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
Featurefun(Mural, c("ENPP2","PLXDC2","PALMD", "MID1", "CEMIP", "GLIS3", "PDZD2", 
                    "COL4A4", "CRISPLD2", "ADAMTS9", "TLE1", "PIPNC1", "ADCY3", "HIPK2", "COL18A1", "SCL22A3",
                    "SOX6", "MYO1D", "MYH11", "ACTA2", "PRUNE2", "RYR2", "TPM1", "LMOD1"), "Mural")

#0 : Mural_2
#1 : Mural_1
#2 : Mural_3
#3 : Mural_2
#4 : Mural_2


###Label

Mural <- RenameIdents(Mural, "0"="Mural", "1"="Mural", "3"="Mural", "4"="Mural", "2"="Doublet") #first pass
Mural <- subset(Mural, Subcelltype=="Mural")

#remove doublets
Mural <- subset(Mural, CMO_classification.global=="Singlet")
Mural <- subset(Mural, percent.mt<5)
saveRDS(Mural, paste0(curr_dir, "/CellType/Secondcelltype/Mural.RDS"))


###Venn of agreement

plts <- ggVennDiagram(list(colnames(Mural)[Mural@meta.data$DoubletFinder_call=="Doublet"], colnames(Mural)[Mural@meta.data$CMO_classification.global=="Doublet"], colnames(Mural)[Mural@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Mural)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Mural/Mural DoubletFinder-CMO-Manual venn.png", sep=""))


###Remove interdoublets

Mural <- RenameIdents(Mural, "0"="Mural_2", "1"="Mural_1", "2"="Mural_3", "3"="Mural_2", "4"="Mural_2") #last pass
Mural@meta.data$Subcelltype <- Mural@active.ident
DimPlot(Mural)

Idents(Mural) <- Mural@meta.data$Subcelltype
Mural <- CellSelector(DimPlot(Mural), Mural, ident = "Peri1")
Peri1 <- subset(Mural, ident="Peri1")
Peri1 <- subset(Peri1, Subcelltype=="Mural_2")

Idents(Mural) <- Mural@meta.data$Subcelltype
Mural <- CellSelector(DimPlot(Mural), Mural, ident = "Peri2")
Peri2 <- subset(Mural, ident="Peri2")
Peri2 <- subset(Peri2, Subcelltype=="Mural_3")

Idents(Mural) <- Mural@meta.data$Subcelltype
Mural <- CellSelector(DimPlot(Mural), Mural, ident = "Mural_1s")
Mural_1 <- subset(Mural, ident="Mural_1s")
Mural_1 <- subset(Mural_1, Subcelltype=="Mural_1")


Peri1_meta <- Peri1@meta.data
Peri2_meta <- Peri2@meta.data
Mural_1_meta <- Mural_1@meta.data
merged_meta <- rbind(Peri1_meta, Peri2_meta, Mural_1_meta)

Mural@meta.data$Subcelltype <- NULL
Mural_meta <- Mural@meta.data

Mural_meta <- rownames_to_column(Mural_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Mural_meta <- left_join(Mural_meta, merged_meta, by="rowname")
Mural_meta <- column_to_rownames(Mural_meta, "rowname")
Mural_meta$Subcelltype <- as.character(Mural_meta$Subcelltype)
Mural_meta[is.na(Mural_meta)==TRUE] <- "Doublet"
Mural@meta.data <- Mural_meta
Mural <- subset(Mural, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Mural) <- Mural@meta.data[["Subcelltype"]] 
levels(Mural) <- c("Mural_2", "Mural_3", "Mural_1")
Mural@meta.data[["Subcelltype"]] <- Idents(Mural)
DimPlot(Mural, raster=F)

Mural <- RunUMAP(object = Mural, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist = 1.5)
saveRDS(Mural, paste0(curr_dir, "/CellType/Finalcelltype/Mural.RDS"))

Ncellfun(Mural, paste0(curr_dir, "/CellType/Mural/Final per APOE.csv", sep=""), "apoe_genotype")
Ncellfun(Mural, paste0(curr_dir, "/CellType/Mural/Final per AD.csv", sep=""), "Cog_Path")
Ncellfun(Mural, paste0(curr_dir, "/CellType/Mural/Final per projid.csv", sep=""), "projid")




##Whole dataset reprocess
###load and merge

Epi <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Epi.RDS", sep=""))
Fibro <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Fibro.RDS", sep=""))
Mural <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Mural.RDS", sep=""))
Endo <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Endo.RDS", sep=""))
Oligo <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Oligo.RDS", sep=""))
Astro <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Astro.RDS", sep=""))
Neuro <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Neuro.RDS", sep=""))
Imm <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Imm.RDS", sep=""))
#Ex_ChP <- readRDS(paste(curr_dir, "/CellType/Ex_ChP.RDS", sep=""))

Ex_ChP <- merge(Epi, c(Fibro, Mural, Endo, Oligo, Astro, Neuro, Imm))
rm (Epi, Fibro, Endo, Mural, Oligo, Astro, Neuro, Imm)


###merge

Fibro$Majorcelltype1 <- Fibro$Majorcelltype
Fibro$Majorcelltype <- NULL
Fibro$Majorcelltype <- Fibro$Majorcelltype1
Fibro$Majorcelltype1 <- NULL
Fibro$Majorcelltype <-"Fibroblast"
Fibro$dropletUtils_emptyDrops_total <- NULL
Fibro$dropletUtils_emptyDrops_logprob <- NULL
Fibro$dropletUtils_emptyDrops_pvalue <- NULL
Fibro$dropletUtils_emptyDrops_limited <- NULL
Fibro$dropletUtils_emptyDrops_fdr <- NULL
Fibro$dropletUtils_BarcodeRank_Inflection <- NULL
Fibro$dropletUtils_BarcodeRank_Knee <- NULL

Ex_ChP$Finalcelltype <- NULL
Ex_ChP$Phase <- NULL
Ex_ChP$S.Score <- NULL
Ex_ChP$G2M.Score <- NULL




Fibro_meta <- Fibro@meta.data
Ex_ChP_meta <- Ex_ChP@meta.data
unique(colnames(Fibro_meta) == colnames(Ex_ChP_meta))

Fibro_meta <- Fibro_meta[(rownames(Fibro_meta) %in% rownames(Ex_ChP_meta)),]
Ex_ChP_meta_less <- Ex_ChP_meta[!(rownames(Ex_ChP_meta) %in% rownames(Fibro_meta)),]
Ex_ChP_meta_Imm_Fib <- rbind(Ex_ChP_meta_less, Fibro_meta)
Ex_ChP_meta_Imm_Fib <- Ex_ChP_meta_Imm_Fib[order(match(rownames(Ex_ChP_meta_Imm_Fib), rownames(Ex_ChP_meta))),]
Ex_ChP_meta_Imm_Fib -> Ex_ChP@meta.data
Idents(Ex_ChP) <- Ex_ChP$Subcelltype
DimPlot(Ex_ChP, raster=F)



gc()
Ex_ChP <- subset(Ex_ChP, Majorcelltype!="Fibroblast")
gc()
Ex_ChP <- merge(Ex_ChP, Fibro)

rm(Fibro)
gc()


###process

NFSP_fun(Ex_ChP, assay="soupXcounts")
use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)

Ex_ChP <- RunHarmony(object = Ex_ChP, group.by.vars = c("Lane"), dims.use = use_dims, assay="soupXcounts", max.iter.harmony = 50, plot_convergence = TRUE) #Correct for batches !!!!
CNDRfun_quick(Ex_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony", resolution=0.1)
#Ex_ChP <- FindClusters(Ex_ChP, resolution = c(0.01), assay=assay)
#plts <- clustree(Ex_ChP)
#  plot(plts)
#  graph2png(plts, paste(curr_dir, "/CellType/", "Ex_ChP/Ex_ChP_cluster_tree.png", sep=""), width=14, height=6)
#CNDRfun(Ex_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
res="0.1"
Idents(Ex_ChP) <- Ex_ChP@meta.data[[paste("soupXcounts_snn_res.", res, sep="")]]
CNDRdimfun(Ex_ChP, reduction = "harmony", res=res, split.by = c("X10X_batch", "projid", "Lane"))
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per Lane.csv", sep=""), "Lane")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per 10X_batch.csv", sep=""), "X10X_batch")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per projid.csv", sep=""), "projid")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
DimPlot(Ex_ChP, raster=F)
ggsave2(paste0(curr_dir, "/CellType/Ex_ChP/Finalmarkers.png", plot=plts))
#return_markersfun(Ex_ChP, assay="soupXcounts", features = Ex_ChP_soupXcounts_var_feat)

Ex_ChP <- CellCycleScoring(Ex_ChP, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident=FALSE)
#DimPlot(Ex_ChP, group.by = "Phase", raster=F)
#DimPlot(Ex_ChP, split.by = "Phase", raster=F)
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/CMO per phase.csv", "Phase"))
saveRDS(Ex_ChP, paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP_EMF.RDS"))


###Remove stragglers

Ex_ChP <- RunUMAP(Ex_ChP, dims=1:20, assay="soupXcounts", reduction="harmony")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Epi")
Epi <- subset(Ex_ChP, ident="Epi")
Epi <- subset(Epi, Subcelltype=="Epi_1" | Subcelltype=="Epi_2" | Subcelltype=="Epi_3")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Fib")
Fibro <- subset(Ex_ChP, ident="Fib")
Fibro <- subset(Fibro, Subcelltype=="Fib_1" | Subcelltype=="Fib_2" | Subcelltype=="Fib_3" | Subcelltype=="Fib_4")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "End")
Endo <- subset(Ex_ChP, ident="End")
Endo <- subset(Endo, Subcelltype=="Endo")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Mur")
Mural <- subset(Ex_ChP, ident="Mur")
Mural <- subset(Mural, Subcelltype=="Mural")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Neur")
Neuro <- subset(Ex_ChP, ident="Neur")
Neuro <- subset(Neuro, Subcelltype=="Neuro")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Ast")
Astro <- subset(Ex_ChP, ident="Ast")
Astro <- subset(Astro, Subcelltype=="Astro")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Oli")
Oligo <- subset(Ex_ChP, ident="Oli")
Oligo <- subset(Oligo, Subcelltype=="Oligo")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Imm")
Imm <- subset(Ex_ChP, ident="Imm")
Imm <- subset(Imm, Subcelltype=="BAM_1" | Subcelltype=="BAM_2" | Subcelltype=="BAM_3" | Subcelltype=="T_Cell")

Epi_meta <- Epi@meta.data
Fibro_meta <- Fibro@meta.data
Endo_meta <- Endo@meta.data
Mural_meta <- Mural@meta.data
Neuro_meta <- Neuro@meta.data
Astro_meta <- Astro@meta.data
Oligo_meta <- Oligo@meta.data
Imm_meta <- Imm@meta.data
merged_meta <- rbind(Epi_meta, Fibro_meta, Endo_meta, Mural_meta, Neuro_meta, Astro_meta, Oligo_meta, Imm_meta)

Ex_ChP@meta.data$Subcelltype <- NULL
Ex_ChP_meta <- Ex_ChP@meta.data

Ex_ChP_meta <- rownames_to_column(Ex_ChP_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Ex_ChP_meta <- left_join(Ex_ChP_meta, merged_meta, by="rowname")
Ex_ChP_meta <- column_to_rownames(Ex_ChP_meta, "rowname")
Ex_ChP_meta$Subcelltype <- as.character(Ex_ChP_meta$Subcelltype)
Ex_ChP_meta[is.na(Ex_ChP_meta)==TRUE] <- "Doublet"
Ex_ChP@meta.data <- Ex_ChP_meta

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
DimPlot(Ex_ChP, label=TRUE, raster=F) #sanity check 

Ex_ChP <- subset(Ex_ChP, subset=Subcelltype=="Doublet", invert=TRUE)


###save and remove srt

saveRDS(Epi, paste(curr_dir, "/CellType/Finalcelltype/Epi.RDS", sep=""))
saveRDS(Fibro, paste(curr_dir, "/CellType/Finalcelltype/Fibro.RDS", sep=""))
saveRDS(Mural, paste(curr_dir, "/CellType/Finalcelltype/Mural.RDS", sep=""))
saveRDS(Endo, paste(curr_dir, "/CellType/Finalcelltype/Endo.RDS", sep=""))
saveRDS(Oligo, paste(curr_dir, "/CellType/Finalcelltype/Oligo.RDS", sep=""))
saveRDS(Astro, paste(curr_dir, "/CellType/Finalcelltype/Astro.RDS", sep=""))
saveRDS(Neuro, paste(curr_dir, "/CellType/Finalcelltype/Neuro.RDS", sep=""))
saveRDS(Imm, paste(curr_dir, "/CellType/Finalcelltype/Imm.RDS", sep=""))
rm(Epi, Fibro, Mural, Endo, Neuro, Astro, Oligo, BAM, T_cell, Imm)
gc()

saveRDS(Ex_ChP, paste(curr_dir, "/CellType/Finalcelltype/Ex_ChP.RDS", sep=""))



###Ncells

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
#Cell numbers
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_projid_Lane.csv", sep=""), group.by="Lane_projid")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_projid.csv", sep=""), group.by="projid")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_cogdx.csv", sep=""), group.by="cogdx")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Cog_Path.csv", sep=""), group.by="Cog_Path")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_apoe.csv", sep=""), group.by="apoe_genotype")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_CMO per phase.csv", "Phase"))

Idents(Ex_ChP) <- Ex_ChP@meta.data$Majorcelltype
#Cell numbers
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Majorcelltype_projid.csv", sep=""), group.by="projid")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Majorcelltype_cogdx.csv", sep=""), group.by="cogdx")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Majorcelltype_Cog_Path.csv", sep=""), group.by="Cog_Path")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Majorcelltype_apoe.csv", sep=""), group.by="apoe_genotype")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_CMO per phase.csv", "Phase"))

