###
### 022 Cell_type_annotation_ATAC.R
###
# Purpose: Cell type annotation and cleaning of single nuclei ATAC dataset.
# Dependencies:
source("/ChP_ATAC/020 Cell_type_annotation_util.R", verbose = FALSE)
curr_dir <- "/ChP_ATAC" #output results to specific folder
set.seed(1234)

##Load and merge samples

ChP_37085610B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_131B.RDS", sep=""))
ChP_50405592B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_141B.RDS", sep=""))
ChP_52052115B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_144B.RDS", sep=""))
ChP_52064033B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_145B.RDS", sep=""))
ChP_78511269B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_163B.RDS", sep=""))
ChP_80487613B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_167B.RDS", sep=""))
Ex_ChP_ATAC <- merge(x=ChP_131B, y=c(ChP_141B, ChP_144B, ChP_145B, ChP_163B, ChP_167B), add.cell.ids=c("26_1", "1_3", "35_1", "1_2", "19_2", "24_2"))
#rm(list = grep("ChP_", ls(), value=TRUE))
gc()


Idents(Ex_ChP_ATAC) <- "Lane"
Ex_ChP_ATAC <- RenameIdents(Ex_ChP_ATAC, 'C1_GEX'='ChP11_1', 'C2_GEX'='ChP11_2', 'C3_GEX'='ChP11_3', 'C4_GEX'='ChP11_4')

Ex_ChP_ATAC@meta.data$Lane <- Idents(Ex_ChP_ATAC)
unique(Ex_ChP_ATAC@meta.data$Lane) #sanity check
Ex_ChP_ATAC@meta.data$Lane_projid <- paste(Ex_ChP_ATAC@meta.data$projid, Ex_ChP_ATAC@meta.data$Lane, sep="_")
Ex_ChP_ATAC@meta.data$sample <- Ex_ChP_ATAC@meta.data$projid
gc()
Ex_ChP_ATAC@meta.data[is.na(Ex_ChP_ATAC@meta.data)==TRUE] <- "unk"

L_p <- unique(Ex_ChP_ATAC@meta.data$Lane_projid)
L_p[order(L_p)] #Check this when making metadata for DEG later
saveRDS(Ex_ChP_ATAC, paste(curr_dir, "/CellType/Firstcelltype/Ex_ChP_ATAC.RDS", sep=""))


##QC

Ex_ChP_ATAC <- subset(Ex_ChP_ATAC, CMO_classification.global=="Singlet")

# add blacklist ratio
Ex_ChP_ATAC$blacklist_ratio <- FractionCountsInRegion(object = Ex_ChP_ATAC, assay = 'ATAC', regions = blacklist_hg38_unified)
Ex_ChP_ATAC$nucleosome_group <- ifelse(Ex_ChP_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
Ex_ChP_ATAC@meta.data$nucleosome_signal <- as.numeric(Ex_ChP_ATAC@meta.data$nucleosome_signal)

VlnPlot(object = Ex_ChP_ATAC, assay="ATAC",
        features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
        pt.size = 0.1, ncol = 4)



Ex_ChP_ATAC <- subset(x = Ex_ChP_ATAC, subset = nCount_ATAC > 100 & nCount_ATAC < 5000 & blacklist_ratio < 0.1 & nucleosome_signal < 4 & TSS.enrichment > 3)


#Celltyping
##Wholeset
###SVD UMAP

DefaultAssay(Ex_ChP_ATAC) <- "ATAC"
peaknormfun(Ex_ChP_ATAC)
dims_use <- c(2,4:30) #avoid dims that correlate to read depth.

Ex_ChP_ATAC@meta.data$Lane <- droplevels(Ex_ChP_ATAC@meta.data$Lane)
Ex_ChP_ATAC@meta.data$nCount_ATAC_round <- round(Ex_ChP_ATAC@meta.data$nCount_ATAC, digits=-2)
Ex_ChP_ATAC@meta.data$nCount_ATAC_round <- pmin(Ex_ChP_ATAC@meta.data$nCount_ATAC_round, 2500)
Ex_ChP_ATAC@meta.data$nCount_ATAC_round <- pmax(Ex_ChP_ATAC@meta.data$nCount_ATAC_round, 500)
Ex_ChP_ATAC@meta.data$nCount_ATAC_round <- as.factor(Ex_ChP_ATAC@meta.data$nCount_ATAC_round)
table(Ex_ChP_ATAC$nCount_ATAC_round)
Ex_ChP_ATAC <- RunHarmony(Ex_ChP_ATAC, group.by.vars=c("Lane", "nCount_ATAC_round"), reduction="lsi", dims.use=dims_use, assay.use="ATAC", max.iter.harmony=50, plot_convergence=TRUE, project.dim=F)
Ex_ChP_ATAC <- RunUMAP(Ex_ChP_ATAC, reduction = "harmony", dims = dims_use)
Ex_ChP_ATAC <- FindNeighbors(object = Ex_ChP_ATAC, reduction = 'harmony', dims = dims_use)

res_range <- c(seq(from = 0.03, to = 0.09, by = 0.02), seq(from = .1, to = .7, by = 0.2)) #!!!!Modify if necessary, it's usually between there though
Ex_ChP_ATAC <- FindClusters(object = Ex_ChP_ATAC, resolution=res_range, verbose=F, algorithm =3)
plts <- clustree(Ex_ChP_ATAC)
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Ex_ChP_ATAC/Ex_ChP_ATAC_cluster_tree.png"), plot = plts, width=14, height=6)

res="0.1"
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data[[paste("ATAC_snn_res.", res, sep="")]]
CNDRdimfun(Ex_ChP_ATAC, reduction="harmony", assay="ATAC", res=res, split.by = c("X10X_batch", "projid", "Lane"))
md <- data.frame(Ex_ChP_ATAC@meta.data$Lane, Ex_ChP_ATAC@meta.data$Lane_projid, Ex_ChP_ATAC@meta.data$projid)
colnames(md) <- c("Lane", "Lane_projid", "projid")
tabyl(md, Lane_projid)
tabyl(md, Lane)

saveRDS(Ex_ChP_ATAC, paste(curr_dir, "/CellType/Firstcelltype/Ex_ChP_ATAC_preprocess.RDS", sep=""))


###Label concurrent cells

Ex_ChP_ATAC <- RenameCells(Ex_ChP_ATAC, new.names = paste0(substr(colnames(Ex_ChP_ATAC), 6, nchar(colnames(Ex_ChP_ATAC))), "_", substr(colnames(Ex_ChP_ATAC), 1, nchar(colnames(Ex_ChP_ATAC))-21)))

Ex_ChP <- readRDS("/Ex_ChP/Ex_ChP4-5/CellType/Finalcelltype/Ex_ChP.RDS")
Ex_ChP_meta <- Ex_ChP@meta.data

Ex_ChP_ATAC@meta.data$Subcelltype <- NULL
Ex_ChP_ATAC_meta <- Ex_ChP_ATAC@meta.data

Ex_ChP_ATAC_meta <- rownames_to_column(Ex_ChP_ATAC_meta, "rowname")
Ex_ChP_meta <- rownames_to_column(Ex_ChP_meta, "rowname")
Ex_ChP_meta$rowname <- paste0(substr(Ex_ChP_meta$rowname, 1, nchar(Ex_ChP_meta$rowname)-2)) #Removing Mjr type number

Ex_ChP_meta <- data.frame(Ex_ChP_meta$rowname, Ex_ChP_meta$Subcelltype)
colnames(Ex_ChP_meta) <- c("rowname", "Subcelltype")
Ex_ChP_ATAC_meta <- left_join(Ex_ChP_ATAC_meta, Ex_ChP_meta, by="rowname")
Ex_ChP_ATAC_meta <- column_to_rownames(Ex_ChP_ATAC_meta, "rowname")
Ex_ChP_ATAC_meta$Subcelltype <- as.character(Ex_ChP_ATAC_meta$Subcelltype)
Ex_ChP_ATAC_meta[is.na(Ex_ChP_ATAC_meta)==TRUE] <- "Unknown"
Ex_ChP_ATAC@meta.data <- Ex_ChP_ATAC_meta

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
plts <- DimPlot(Ex_ChP_ATAC, label=F, raster=F, order=c("Epi_1", "Epi_2", "Epi_3", "Epi_4", "Fib_1", "Fib_2", "Fib_3", "Fib_4", "Endo", "Mural", "BAM_1", "BAM_2", "BAM_3", "T_Cell", "Oligo", "Neuro", "Unknown")) #sanity check 
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Ex_ChP_ATAC/Concurrent cells to SN.png"), plot = plts, width=14, height=6)


###Integrate

Ex_ChP <- FindTopFeatures(Ex_ChP, min.cutoff = 10)
Ex_ChP <- RunTFIDF(Ex_ChP)
Ex_ChP <- RunSVD(Ex_ChP)
Ex_ChP <- RunUMAP(Ex_ChP, reduction = "lsi", dims = 2:30, return.model=TRUE)
Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
DimPlot(Ex_ChP, label=T, raster=F)
saveRDS(Ex_ChP, paste(curr_dir, "/CellType/Firstcelltype/Ex_ChP_lsi.RDS", sep=""))
#!!!!lsi integration doesn't match labelling of concurrent cells!!!!

Ex_ChP <- RunUMAP(Ex_ChP, reduction = "pca", dims = 2:30, return.model=TRUE)
DefaultAssay(Ex_ChP_ATAC) <- "ACTIVITY"
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = Ex_ChP,
  query = Ex_ChP_ATAC,
  reference.reduction = "pca",
  reduction = "pcaproject",
  dims = 2:30
)

# map query onto the reference dataset
Ex_ChP_ATAC <- MapQuery(
  anchorset = transfer.anchors,
  reference = Ex_ChP,
  query = Ex_ChP_ATAC,
  refdata = Ex_ChP$Subcelltype,
  reference.reduction = "pca",
  new.reduction.name = "ref.pca",
  reduction.model = 'umap'
)

DimPlot(Ex_ChP, label=T, raster=F)
DimPlot(Ex_ChP_ATAC, label=T, raster=F, reduction="ref.umap")

saveRDS(Ex_ChP_ATAC, paste(curr_dir, "/CellType/Firstcelltype/Ex_ChP_ATAC_integrated.RDS", sep=""))


###Label remaining cells with activity assay

# normalize gene activities
DefaultAssay(Ex_ChP_ATAC) <- "ACTIVITY"
Ex_ChP_ATAC <- NormalizeData(Ex_ChP_ATAC)
Ex_ChP_ATAC <- ScaleData(Ex_ChP_ATAC, features = rownames(Ex_ChP_ATAC))

res="0.1"
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data[[paste("ATAC_snn_res.", res, sep="")]]
return_markersfun(Ex_ChP_ATAC, assay="ACTIVITY")


###Resolution 0.1
###Dividing cells spend more time in S phase

Ex_ChP_ATAC <- CellCycleScoring(Ex_ChP_ATAC, assay="ACTIVITY", g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident=FALSE)
DimPlot(Ex_ChP_ATAC, group.by = "Phase", raster=F)
DimPlot(Ex_ChP_ATAC, split.by = "Phase", raster=F)
Ncellfun(Ex_ChP_ATAC, paste0(curr_dir, "/Ex_ChP_ATAC/CMO per phase.csv", "Phase"))



curr_clust <- 5
Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Ex_ChP_ATAC_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Ex_ChP_ATAC, curr_clust)

search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, known_fibro_prog)

search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_light) %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Ex_ChP_ATAC, curr_clust, known_dividing %>% unique())

search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Ex_ChP_ATAC, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Ex_ChP_ATAC, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
Featurefun(Ex_ChP_ATAC, c("TTR", "HTR2C", "GPX3", "DCN", "LEPR", "FN1", "VWF", "PECAM1", "INSR", "MYH11", "TAGLN", "ACTA2", "STAB1", "CD163", "MSR1", "CDH12", "OTX2", "JUND", "CDH1", "KRT8", "KRT18", "LHX2", "GFAP", "KCNMA1", "FOXJ1", "DNAH6", "CFAP43", "CFAP54", "SPAG17", "TUBB4B", "CLDN1", "CLDN3", "CLDN5", "OCLN", "TJP1", "TJP2", "TRPM3", "GPX3", "MT-CO1", "MT-CO2", "MT-CO3", "AQP1", "ATP1A1", "SLC2A12", "SLC25A6", "SLC25A23", "SLC22A17", "SLC7A2", "SLC25A3", "SLC4A2", "SLC4A4", "FOLR1", "PIFO", "IGFBP2", "LUM", "VIM", "S100B", "CD2", "DCN", "FN1", "LEPR", "CDH11", "STXBP6", "FBLN1", "KCNMA1", "SLC4A4", "PLXNA4", "ALPL", "DPEP1", "DLK1", "NRG1", "ROBO1", "TRPM3", "MYRIP", "GPM6A", "CLDN11", "ECM2", "SLC6A20", "SLC13A3", "ABCA9", "LAMA2", "LAMA4", "COL18A1", "COL4A2", "COL1A2", "COL4A4", "COL5A2",  "COL3A1", "COL6A3", "PDGFRA", "SLC4A7", "SLC2A3", "LUM"), "highlights_genes",  assay="ACTIVITY")
plts <-VlnPlot(Ex_ChP_ATAC, features=c("TTR", "HTR2C", "GPX3", "DCN", "LEPR", "FN1", "VWF", "PECAM1", "INSR", "MYH11", "TAGLN", "ACTA2", "STAB1", "CD163", "MSR1", "OTX2", "JUND", "CDH1", "KRT8", "KRT18", "LHX2", "GFAP", "FOXJ1", "DNAH6", "CFAP43", "CFAP54", "SPAG17", "TUBB4B", "CLDN1", "CLDN3", "CLDN5", "OCLN", "TJP1", "TJP2", "GPX3", "AQP1", "ATP1A1", "SLC2A12", "SLC25A6", "SLC25A23", "SLC22A17", "SLC7A2", "SLC25A3", "SLC4A2", "SLC4A4", "FOLR1", "PIFO", "IGFBP2", "LUM", "VIM", "S100B", "CD2", "DCN", "FN1", "LEPR", "CDH11", "STXBP6", "FBLN1", "SLC4A4", "ALPL", "DPEP1", "DLK1", "NRG1", "MYRIP", "GPM6A", "CLDN11", "ECM2", "SLC6A20", "SLC13A3", "ABCA9", "LAMA4", "COL18A1", "COL4A2", "COL1A2", "COL4A4", "COL5A2",  "COL3A1", "COL6A3", "PDGFRA", "SLC4A7", "SLC2A3", "LUM"), stack=T, raster=T, same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
graph2png(plts, paste(curr_dir, "/CellType/vln.png", sep=""), width=3+82*.4, height=1.5+20*1)
remove(plts)


#Ex_ChP_ATAC
#0: Epi
#1: Epi
#2: Epi
#3: Fib
#4: Epi
#5: Epi
#6: Imm
#7: Paren
#8: Epi
#9: Fib
#10: Mural

  
  
  ###Rename all identities
  
Ex_ChP_ATAC <- RenameIdents(Ex_ChP_ATAC, "0"="Epi", "1"="Epi", "2"="Epi", "3"="Fib", "4"="Epi", "5"="Epi", "6"="Imm", "7"="Paren", "8"="Epi", "9"="Fib", "10"="Mural")
Ex_ChP_ATAC@meta.data$Firstcelltype <- Idents(Ex_ChP_ATAC)
saveRDS(Ex_ChP_ATAC, paste(curr_dir, "/CellType/Firstcelltype/Ex_ChP_ATAC_Firsttype.RDS", sep=""))


##Subtypes

for(i in unique(Ex_ChP_ATAC@meta.data$Firstcelltype)){
  sub <- subset(Ex_ChP_ATAC, Firstcelltype==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/", i, ".RDS", sep=""))
}



##Epi
###PCA

Epi <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi.RDS"))
gc()

Featurefun(Epi, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="ACTIVITY")
#Remove cells that don't belong
#Epi <- subset(Epi, HTR2C < 3)
#Epi <- subset(Epi, TTR < 3)
#Epi <- subset(Epi, CLIC6 < 3)
#Epi <- subset(Epi, AQP1 < 3)
#Epi <- subset(Epi, CFAP54 < 3)
#Epi <- subset(Epi, FAM237A < 3)
#Epi <- subset(Epi, ARL13B < 3)
#Epi <- subset(Epi, TRPM3 < 3)
Epi <- subset(Epi, SLC4A4 < 3)
Epi <- subset(Epi, LAMA4 < 3)
Epi <- subset(Epi, CLDN1 < 3)
Epi <- subset(Epi, CLDN11 < 3)
Epi <- subset(Epi, PRICKLE1 < 3)
Epi <- subset(Epi, CEMIP < 3)
Epi <- subset(Epi, COL6A3 < 3)
Epi <- subset(Epi, SLC6A13 < 3)
Epi <- subset(Epi, VWF < 3)
Epi <- subset(Epi, TAGLN < 3)
Epi <- subset(Epi, MYH11 < 3)
Epi <- subset(Epi, NOTCH3 < 3)
Epi <- subset(Epi, MAP2 < 3)
#Epi <- subset(Epi, CLU < 3)
#Epi <- subset(Epi, AQP4 < 3)
Epi <- subset(Epi, MOBP < 3)
Epi <- subset(Epi, MBP < 3)
Epi <- subset(Epi, CDK19 < 3)
Epi <- subset(Epi, CD68 < 3)
Epi <- subset(Epi, F13A1 < 3)
Epi <- subset(Epi, STAB1 < 3)
Epi <- subset(Epi, CD163 < 3)
Epi <- subset(Epi, LYVE1 < 3)
Epi <- subset(Epi, P2RY12 < 3)
Epi <- subset(Epi, THEMIS < 3)

DefaultAssay(Epi) <- "ATAC"
peaknormfun(Epi)
dims_use <- c(3:30) #avoid dims that correlate to read depth.

table(Epi$nCount_ATAC_round)
Epi <- RunHarmony(Epi, group.by.vars=c("Lane", "nCount_ATAC_round"), reduction="lsi", dims.use=dims_use, assay.use="ATAC", max.iter.harmony=50, plot_convergence=TRUE, project.dim=F)
Epi <- RunUMAP(Epi, reduction = "harmony", dims = dims_use, min.dist=2)
Epi <- FindNeighbors(object = Epi, reduction = 'harmony', dims = dims_use)

res_range <- c(seq(from = 0.03, to = 0.09, by = 0.02), seq(from = .1, to = .7, by = 0.2)) #!!!!Modify if necessary, it's usually between there though
Epi <- FindClusters(object = Epi, resolution=res_range, verbose=F, algorithm=3)
plts <- clustree(Epi)
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Epi/Epi_cluster_tree.png"), plot = plts, width=14, height=6)

res="0.03"
Idents(Epi) <- Epi@meta.data[["ATAC_snn_res.0.03"]]
CNDRdimfun(Epi, reduction="harmony", assay="ATAC", res=res, split.by = c("X10X_batch", "projid", "Lane", "CMO_classification.global"))
md <- data.frame(Epi@meta.data$Lane, Epi@meta.data$Lane_projid, Epi@meta.data$projid, Epi@meta.data$Cog_Path)
colnames(md) <- c("Lane", "Lane_projid", "projid", "Cog_Path")
tabyl(md, Lane_projid)
tabyl(md, Lane)
tabyl(md, Cog_Path)
Ncellfun(Epi, "Cog_Path", group.by="Cog_Path")


# normalize gene activities
DefaultAssay(Epi) <- "ACTIVITY"
Epi <- NormalizeData(Epi)
Epi <- ScaleData(Epi, features = rownames(Epi))

res="0.03"
Idents(Epi) <- Epi@meta.data[[paste("ATAC_snn_res.", res, sep="")]]
return_markersfun(Epi, assay="ACTIVITY")



curr_clust <- 4
Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

Top20_Featuresfun(Epi, curr_clust)

search_markersfun(Epi, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Epi, curr_clust, known_fibro_prog)

search_markersfun(Epi, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Epi, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

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
Featurefun(Epi, c("TTR", "HTR2C", "CROCC", "OTX2", "JUND", "CDH1", "KRT8", "KRT18", "LHX2", "GFAP", "KCNMA1", "FOXJ1", "DNAH6", "CFAP43", "CFAP54", "SPAG17", "TUBB4B", "CLDN1", "CLDN3", "CLDN5", "OCLN", "TJP1", "TJP2", "TRPM3", "GPX3", "MT-CO1", "MT-CO2", "MT-CO3", "AQP1", "ATP1A1", "SLC2A12"), "highlights_genes",  assay="ACTIVITY")
DimPlot(Epi)
plts <- VlnPlot(Epi, features=c("TTR", "HTR2C", "CROCC", "OTX2", "JUND", "CDH1", "KRT8", "KRT18", "LHX2", "GFAP", "KCNMA1", "FOXJ1", "DNAH6", "CFAP43", "CFAP54", "SPAG17", "TUBB4B", "CLDN1", "CLDN3", "CLDN5", "OCLN", "TJP1", "TJP2", "TRPM3", "GPX3", "MT-CO1", "MT-CO2", "MT-CO3", "AQP1", "ATP1A1", "SLC2A12"), stack=T, raster=T, same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/vln.png", sep=""), width=3+82*.4, height=1.5+10*1)
remove(plts)

#0: Cilia+Secretion
#1: Low markers
#2: Cilia++
#3: Secretion

###Label
Epi <- RenameIdents(Epi, "0"="Epi_3", "1"="Epi_1", "2"="Epi_2", "3"="Epi_4", "4"="Epi_2", "5"="Epi_2", "6"="Epi_2", "7"="Epi_2", "8"="Epi_2", "9"="Epi_2", "10"="Epi_2", "11"="Epi_2")
Epi@meta.data$Subcelltype <- Epi@active.ident

saveRDS(Epi, paste0(curr_dir, "/CellType/Finalcelltype/Epi.RDS"))


##Fib
###PCA

Fib <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Fib.RDS"))
gc()

Featurefun(Fib, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="ACTIVITY")

#Remove cells that don't belong
Fib <- subset(Fib, HTR2C < 3)
Fib <- subset(Fib, TTR < 3)
Fib <- subset(Fib, CLIC6 < 3)
#Fib <- subset(Fib, AQP1 < 3)
Fib <- subset(Fib, CFAP54 < 3)
Fib <- subset(Fib, FAM227A < 3)
Fib <- subset(Fib, ARL13B < 3)
#Fib <- subset(Fib, TRPM3 < 3)
#Fib <- subset(Fib, SLC4A4 < 3)
#Fib <- subset(Fib, LAMA2 < 3)
#Fib <- subset(Fib, LAMA4 < 3)
#Fib <- subset(Fib, CLDN1 < 3)
#Fib <- subset(Fib, CLDN11 < 3)
#Fib <- subset(Fib, PRICKLE1 < 3)
#Fib <- subset(Fib, CEMIP < 3)
#Fib <- subset(Fib, FOXP2 < 3)
#Fib <- subset(Fib, COL6A3 < 3)
#Fib <- subset(Fib, SLC6A13 < 3)
Fib <- subset(Fib, PECAM1 < 3)
Fib <- subset(Fib, VWF < 3)
Fib <- subset(Fib, TAGLN < 3)
Fib <- subset(Fib, MYH11 < 3)
Fib <- subset(Fib, NOTCH3 < 3)
Fib <- subset(Fib, MAP2 < 3)
Fib <- subset(Fib, CLU < 3)
Fib <- subset(Fib, AQP4 < 3)
Fib <- subset(Fib, MOBP < 3)
Fib <- subset(Fib, MBP < 3)
Fib <- subset(Fib, CDK19 < 3)
Fib <- subset(Fib, CD68 < 3)
Fib <- subset(Fib, F13A1 < 3)
Fib <- subset(Fib, STAB1 < 3)
Fib <- subset(Fib, CD163 < 3)
Fib <- subset(Fib, LYVE1 < 3)
Fib <- subset(Fib, P2RY12 < 3)
Fib <- subset(Fib, THEMIS < 3)

DefaultAssay(Fib) <- "ATAC"
peaknormfun(Fib)
dims_use <- c(2:30) #avoid dims that correlate to read depth.

Fib@meta.data$Lane <- droplevels(Fib@meta.data$Lane)
Fib@meta.data$nCount_ATAC_round <- pmax(as.numeric(Fib@meta.data$nCount_ATAC_round), 5)
Fib@meta.data$nCount_ATAC_round <- as.factor(Fib@meta.data$nCount_ATAC_round)
table(Fib$nCount_ATAC_round)
Fib <- RunHarmony(Fib, group.by.vars=c("Lane", "nCount_ATAC_round"), reduction="lsi", dims.use=dims_use, assay.use="ATAC", max.iter.harmony=50, plot_convergence=TRUE, project.dim=F)
Fib <- RunUMAP(Fib, reduction = "harmony", dims = dims_use, min.dist=2)
Fib <- FindNeighbors(object = Fib, reduction = 'harmony', dims = dims_use)

res_range <- c(seq(from = 0.03, to = 0.09, by = 0.02), seq(from = .1, to = .7, by = 0.2)) #!!!!Modify if necessary, it's usually between there though
Fib <- FindClusters(object = Fib, resolution=res_range, verbose=F, algorithm =3)
plts <- clustree(Fib)
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Fib/Fib_Vasc_cluster_tree.png"), plot = plts, width=14, height=6)

res="0.03"
Idents(Fib) <- Fib@meta.data[[paste("ATAC_snn_res.", res, sep="")]]
CNDRdimfun(Fib, reduction="harmony", assay="ATAC", res=res, split.by = c("X10X_batch", "projid", "Lane", "CMO_classification.global"))
md <- data.frame(Fib@meta.data$Lane, Fib@meta.data$Lane_projid, Fib@meta.data$projid)
colnames(md) <- c("Lane", "Lane_projid", "projid")
tabyl(md, Lane_projid)
tabyl(md, "count")



# normalize gene activities
DefaultAssay(Fib) <- "ACTIVITY"
Fib <- NormalizeData(Fib)
Fib <- ScaleData(Fib, features = rownames(Fib))

res="0.03"
Idents(Fib) <- Fib@meta.data[[paste("ATAC_snn_res.", res, sep="")]]
return_markersfun(Fib, assay="ACTIVITY")



curr_clust <- 1
Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

Top20_Featuresfun(Fib, curr_clust)

search_markersfun(Fib, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Fib, curr_clust, known_fibro_prog)

search_markersfun(Fib, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Fib, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Fib, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Fib, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Fib, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Fib, curr_clust, c(known_light) %>% unique())
search_markersfun(Fib, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Fib, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Fib, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Fib, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Fib, curr_clust, known_dividing %>% unique())

search_markersfun(Fib, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Fib, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Fib, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

search_markersfun(Fib, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Fib, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Fib, curr_clust, c(Cell_Stress, whole_cells) %>% unique())


Featurefun(Fib, unique(c("DCN", "LEPR", "FN1", "VWF", "PECAM1", "INSR", "MYH11", "TAGLN", "ACTA2", "ATP1A1", "SLC2A12", "SLC25A6", "SLC25A23", "SLC22A17", "SLC7A2", "SLC25A3", "SLC4A2", "SLC4A4", "FOLR1", "PIFO", "IGFBP2", "LUM", "VIM", "S100B", "CD2", "DCN", "FN1", "LEPR", "CDH11", "STXBP6", "FBLN1", "KCNMA1", "SLC4A4", "PLXNA4", "ALPL", "DPEP1", "DLK1", "NRG1", "ROBO1", "TRPM3", "MYRIP", "GPM6A", "CLDN11", "ECM2", "SLC6A20", "SLC13A3", "ABCA9", "LAMA2", "LAMA4", "COL18A1", "COL4A2", "COL1A2", "COL4A4", "COL5A2",  "COL3A1", "COL6A3", "PDGFRA", "SLC4A7", "SLC2A3", "LUM")), "highlights_genes",  assay="ACTIVITY")
plts <- VlnPlot(Fib, features=unique(c("DCN", "LEPR", "FN1", "VWF", "PECAM1", "INSR", "MYH11", "TAGLN", "ACTA2", "ATP1A1", "SLC2A12", "SLC25A6", "SLC25A23", "SLC22A17", "SLC7A2", "SLC25A3", "SLC4A2", "SLC4A4", "FOLR1", "PIFO", "IGFBP2", "LUM", "VIM", "S100B", "CD2", "DCN", "FN1", "LEPR", "CDH11", "STXBP6", "FBLN1", "KCNMA1", "SLC4A4", "PLXNA4", "ALPL", "DPEP1", "DLK1", "NRG1", "ROBO1", "TRPM3", "MYRIP", "GPM6A", "CLDN11", "ECM2", "SLC6A20", "SLC13A3", "ABCA9", "LAMA2", "LAMA4", "COL18A1", "COL4A2", "COL1A2", "COL4A4", "COL5A2",  "COL3A1", "COL6A3", "PDGFRA", "SLC4A7", "SLC2A3", "LUM")), stack=T, raster=T, same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/vln.png", sep=""), width=3+82*.4, height=1.5+20*1)
remove(plts)

###Label
Fib@meta.data$Subcelltype <- "Fibroblast"
saveRDS(Fib, paste0(curr_dir, "/CellType/Finalcelltype/Fibro.RDS"))


##Imm
###PCA
Imm <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Imm.RDS"))
gc()

Featurefun(Imm, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="ACTIVITY")

#Remove cells that don't belong
Imm <- subset(Imm, HTR2C < 4)
Imm <- subset(Imm, TTR < 4)
Imm <- subset(Imm, CLIC6 < 4)
#Imm <- subset(Imm, AQP1 < 3)
Imm <- subset(Imm, CFAP54 < 3)
Imm <- subset(Imm, FAM227A < 3)
Imm <- subset(Imm, ARL13B < 3)
Imm <- subset(Imm, LAMA4 < 3)
Imm <- subset(Imm, CLDN1 < 3)
Imm <- subset(Imm, CLDN11 < 3)
Imm <- subset(Imm, PRICKLE1 < 3)
Imm <- subset(Imm, CEMIP < 3)
Imm <- subset(Imm, COL6A3 < 3)
Imm <- subset(Imm, SLC6A13 < 3)
Imm <- subset(Imm, VWF < 3)
Imm <- subset(Imm, TAGLN < 3)
Imm <- subset(Imm, MYH11 < 3)
Imm <- subset(Imm, NOTCH3 < 3)
Imm <- subset(Imm, MAP2 < 3)
#Imm <- subset(Imm, CLU < 3)
Imm <- subset(Imm, AQP4 < 3)
Imm <- subset(Imm, MOBP < 3)
Imm <- subset(Imm, MBP < 3)
#Imm <- subset(Imm, CDK19 < 3)
#Imm <- subset(Imm, CD68 < 3)
#Imm <- subset(Imm, F13A1 < 3)
#Imm <- subset(Imm, STAB1 < 3)
#Imm <- subset(Imm, CD163 < 3)
#Imm <- subset(Imm, LYVE1 < 3)
#Imm <- subset(Imm, P2RY12 < 3)
#Imm <- subset(Imm, THEMIS < 3)

DefaultAssay(Imm) <- "ATAC"
peaknormfun(Imm)
dims_use <- c(3:30) #avoid dims that correlate to read depth.

Imm@meta.data$Lane <- droplevels(Imm@meta.data$Lane)
table(Imm$nCount_ATAC_round)
Imm <- RunHarmony(Imm, group.by.vars=c("Lane", "nCount_ATAC_round"), reduction="lsi", dims.use=dims_use, assay.use="ATAC", max.iter.harmony=50, plot_convergence=TRUE, project.dim=F)
Imm <- RunUMAP(Imm, reduction = "harmony", dims = dims_use, min.dist=2)
Imm <- FindNeighbors(object = Imm, reduction = 'harmony', dims = dims_use)

res_range <- c(seq(from = 0.03, to = 0.09, by = 0.02), seq(from = .1, to = .7, by = 0.2)) #!!!!Modify if necessary, it's usually between there though
Imm <- FindClusters(object = Imm, resolution=res_range, verbose=F, algorithm =3)
plts <- clustree(Imm)
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Imm/Imm_cluster_tree.png"), plot = plts, width=14, height=6)

res="0.7"
Idents(Imm) <- Imm@meta.data[[paste("ATAC_snn_res.", res, sep="")]]
CNDRdimfun(Imm, reduction="harmony", assay="ATAC", res=res, split.by = c("X10X_batch", "projid", "Lane", "CMO_classification.global"))
md <- data.frame(Imm@meta.data$Lane, Imm@meta.data$Lane_projid, Imm@meta.data$projid, Imm@meta.data$Phase)
colnames(md) <- c("Lane", "Lane_projid", "projid", "Phase")
tabyl(md, Lane)
tabyl(md, Phase)
Ncellfun(Imm, "Imm.csv", group.by="Phase")


# normalize gene activities
DefaultAssay(Imm) <- "ACTIVITY"
Imm <- NormalizeData(Imm)
Imm <- ScaleData(Imm, features = rownames(Imm))
res="0.7"
Idents(Imm) <- Imm@meta.data[[paste("ATAC_snn_res.", res, sep="")]]
return_markersfun(Imm, assay="ACTIVITY")



curr_clust <- 0
Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

Top20_Featuresfun(Imm, curr_clust)

search_markersfun(Imm, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Imm, curr_clust, known_fibro_prog)

search_markersfun(Imm, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Imm, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

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
Featurefun(Imm, c("F13A1", "LYVE1", "MRC1", "STAB1", "CD163", "CD68", "MERTK", "HCK", "IL4R", "CD74", "HLA-DRA", "S100A9", "CSF1R", "JAK2", "CD86", "VCAM1", "CD83", "LAMP1", "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1", "MKI67", "CENPF", "TOP2A", "ITGAM", "S100A8", "CCR2", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11", "THEMIS", "CD3E", "PDE7A", "CAMK4", "HLA-B", "CD8A"), "highlights_genes",  assay="ACTIVITY")
plts <- VlnPlot(Imm, features=c("F13A1", "LYVE1", "MRC1", "STAB1", "CD163", "CD68", "MERTK", "HCK", "IL4R", "CD74", "HLA-DRA", "S100A9", "CSF1R", "JAK2", "CD86", "VCAM1", "CD83", "LAMP1", "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1", "MKI67", "CENPF", "TOP2A", "ITGAM", "S100A8", "CCR2", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11", "THEMIS", "CD3E", "PDE7A", "CAMK4", "HLA-B", "CD8A"), stack=T, raster=T, same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
graph2png(plts, paste(curr_dir, "/CellType/vln.png", sep=""), width=3+82*.4, height=1.5+20*1)
remove(plts)


###Label

Imm <- RenameIdents(Imm, "0"="BAM_2", "1"="BAM_1", "2"="BAM_2")
Imm@meta.data$Subcelltype <- Imm@active.ident

saveRDS(Imm, paste0(curr_dir, "/CellType/Finalcelltype/Imm.RDS"))


##Whole dataset reprocess
###load and merge

Epi <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Epi.RDS", sep=""))
Fibro <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Fib.RDS", sep=""))
Imm <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Imm.RDS", sep=""))

Epi_meta <- Epi@meta.data
Fibro_meta <- Fibro@meta.data
Imm_meta <- Imm@meta.data
merged_meta <- rbind(Epi_meta, Fibro_meta, Imm_meta)

Ex_ChP_ATAC@meta.data$Subcelltype.concurrent <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC@meta.data$Subcelltype <- NULL
Ex_ChP_ATAC_meta <- Ex_ChP_ATAC@meta.data

Ex_ChP_ATAC_meta <- rownames_to_column(Ex_ChP_ATAC_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Ex_ChP_ATAC_meta <- left_join(Ex_ChP_ATAC_meta, merged_meta, by="rowname")
Ex_ChP_ATAC_meta <- column_to_rownames(Ex_ChP_ATAC_meta, "rowname")
Ex_ChP_ATAC_meta$Subcelltype <- as.character(Ex_ChP_ATAC_meta$Subcelltype)
Ex_ChP_ATAC_meta[is.na(Ex_ChP_ATAC_meta)==TRUE] <- "Doublet"
Ex_ChP_ATAC@meta.data <- Ex_ChP_ATAC_meta
Ex_ChP_ATAC <- subset(Ex_ChP_ATAC, Subcelltype!="Doublet")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
DimPlot(Ex_ChP_ATAC, label=TRUE, raster=F) #sanity check 


###SVD UMAP

DefaultAssay(Ex_ChP_ATAC) <- "ATAC"
peaknormfun(Ex_ChP_ATAC)
dims_use <- c(2,4:30) #avoid dims that correlate to read depth.

table(Ex_ChP_ATAC$nCount_ATAC_round)
Ex_ChP_ATAC <- RunHarmony(Ex_ChP_ATAC, group.by.vars=c("Lane", "nCount_ATAC_round"), reduction="lsi", dims.use=dims_use, assay.use="ATAC", max.iter.harmony=50, plot_convergence=TRUE, project.dim=F)
Ex_ChP_ATAC <- RunUMAP(Ex_ChP_ATAC, reduction = "harmony", dims = dims_use)
Ex_ChP_ATAC <- FindNeighbors(object = Ex_ChP_ATAC, reduction = 'harmony', dims = dims_use)

res_range <- c(seq(from = 0.03, to = 0.09, by = 0.02), seq(from = .1, to = .7, by = 0.2)) #!!!!Modify if necessary, it's usually between there though
Ex_ChP_ATAC <- FindClusters(object = Ex_ChP_ATAC, resolution=res_range, verbose=F, algorithm =3)
plts <- clustree(Ex_ChP_ATAC)
plot(plts)
ggsave2(paste0(curr_dir, "/CellType/Ex_ChP_ATAC/Ex_ChP_ATAC_cluster_tree.png"), plot = plts, width=14, height=6)

res="0.1"
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data[[paste("ATAC_snn_res.", res, sep="")]]
CNDRdimfun(Ex_ChP_ATAC, reduction="harmony", assay="ATAC", res=res, split.by = c("X10X_batch", "projid", "Lane"))

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
CNDRdimfun(Ex_ChP_ATAC, reduction="Sub", assay="ATAC", res=res, split.by = c("X10X_batch", "projid", "Lane"), cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700", "Fibroblast"="#3249A6", "BAM_1"="#70C25B", "BAM_2"="#0E9554"))

Ex_ChP_ATAC <- RenameIdents(Ex_ChP_ATAC, "Epi_1"="Epithelial", "Epi_2"="Epithelial", "Epi_3"="Epithelial", "Epi_4"="Epithelial", "BAM_1"="Immune", "BAM_2"="Immune")
Idents(Ex_ChP_ATAC) -> Ex_ChP_ATAC@meta.data$Majorcelltype
CNDRdimfun(Ex_ChP_ATAC, reduction="Final", assay="ATAC", res=res, split.by = c("X10X_batch", "projid", "Lane"), cols=c("#C48F45", "#3249A6", "#26532B"))

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype.concurrent
plts <- DimPlot(Ex_ChP_ATAC, label=TRUE, raster=F, order=c("Epi_1", "Epi_2", "Epi_3", "Epi_4", "Fib_1", "Fib_2", "Fib_3", "Fib_4", "Endo", "Mural", "BAM_1", "BAM_2", "BAM_3", "Oligo", "Neuro", "Unknown")) #sanity check 
plot(plts)

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype.concurrent
plts <- DimPlot(Ex_ChP_ATAC, label=TRUE, raster=F, order=c("Epi_1", "Epi_2", "Epi_3", "Epi_4", "Fib_1", "Fib_2", "Fib_3", "Fib_4", "Endo", "Mural", "BAM_1", "BAM_2", "BAM_3", "Oligo", "Neuro", "Unknown")) #sanity check 
plot(plts)

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ncellfun(Ex_ChP_ATAC, "CellType/Finalcelltype/Cog_Path.csv", "Cog_Path")
Ncellfun(Ex_ChP_ATAC, "CellType/Finalcelltype/projid.csv", "projid")
Ncellfun(Ex_ChP_ATAC, "CellType/Finalcelltype/Lane.csv", "Lane")

saveRDS(Ex_ChP_ATAC, paste(curr_dir, "/CellType/Finalcelltype/Ex_ChP_ATAC.RDS", sep=""))


###Remove stragglers
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC <- CellSelector(DimPlot(Ex_ChP_ATAC), Ex_ChP_ATAC, ident = "Epi")
Epi_1 <- subset(Ex_ChP_ATAC, ident="Epi")
Epi_1 <- subset(Epi_1, Subcelltype=="Epi_1")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC <- CellSelector(DimPlot(Ex_ChP_ATAC), Ex_ChP_ATAC, ident = "Epi")
Epi_2 <- subset(Ex_ChP_ATAC, ident="Epi")
Epi_2 <- subset(Epi_2, Subcelltype=="Epi_2")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC <- CellSelector(DimPlot(Ex_ChP_ATAC), Ex_ChP_ATAC, ident = "Epi")
Epi_3 <- subset(Ex_ChP_ATAC, ident="Epi")
Epi_3 <- subset(Epi_3, Subcelltype=="Epi_3")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC <- CellSelector(DimPlot(Ex_ChP_ATAC), Ex_ChP_ATAC, ident = "Epi")
Epi_4 <- subset(Ex_ChP_ATAC, ident="Epi")
Epi_4 <- subset(Epi_4, Subcelltype=="Epi_4")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC <- CellSelector(DimPlot(Ex_ChP_ATAC), Ex_ChP_ATAC, ident = "Fib")
Fibro <- subset(Ex_ChP_ATAC, ident="Fib")
Fibro <- subset(Fibro, Subcelltype=="Fib_1" | Subcelltype=="Fib_2" | Subcelltype=="Fib_3" | Subcelltype=="Fib_4")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC <- CellSelector(DimPlot(Ex_ChP_ATAC), Ex_ChP_ATAC, ident = "Mur")
Mural <- subset(Ex_ChP_ATAC, ident="Mur")
Mural <- subset(Mural, Subcelltype=="Mural")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC <- CellSelector(DimPlot(Ex_ChP_ATAC), Ex_ChP_ATAC, ident = "Imm")
Imm <- subset(Ex_ChP_ATAC, ident="Imm")
Imm <- subset(Imm, Subcelltype=="BAM_1" | Subcelltype=="BAM_2")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC <- CellSelector(DimPlot(Ex_ChP_ATAC), Ex_ChP_ATAC, ident = "Paren")
Paren <- subset(Ex_ChP_ATAC, ident="Paren")
Paren <- subset(Paren, Subcelltype=="Parenchyma")

Epi1_meta <- Epi_1@meta.data
Epi2_meta <- Epi_2@meta.data
Epi3_meta <- Epi_3@meta.data
Epi4_meta <- Epi_4@meta.data
Fibro_meta <- Fibro@meta.data
Mural_meta <- Mural@meta.data
Imm_meta <- Imm@meta.data
Paren_meta <- Paren@meta.data
merged_meta <- rbind(Epi1_meta, Epi2_meta, Epi3_meta, Epi4_meta, Fibro_meta, Mural_meta, Imm_meta, Paren_meta)

Ex_ChP_ATAC@meta.data$Subcelltype <- NULL
Ex_ChP_ATAC_meta <- Ex_ChP_ATAC@meta.data

Ex_ChP_ATAC_meta <- rownames_to_column(Ex_ChP_ATAC_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Ex_ChP_ATAC_meta <- left_join(Ex_ChP_ATAC_meta, merged_meta, by="rowname")
Ex_ChP_ATAC_meta <- column_to_rownames(Ex_ChP_ATAC_meta, "rowname")
Ex_ChP_ATAC_meta$Subcelltype <- as.character(Ex_ChP_ATAC_meta$Subcelltype)
Ex_ChP_ATAC_meta[is.na(Ex_ChP_ATAC_meta)==TRUE] <- "Doublet"
Ex_ChP_ATAC@meta.data <- Ex_ChP_ATAC_meta
Ex_ChP_ATAC <- subset(Ex_ChP_ATAC, Subcelltype!="Doublet")

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
DimPlot(Ex_ChP_ATAC, label=TRUE, raster=F) #sanity check 


###Rename cells and add fragments
Ex_ChP_ATAC <- RenameCells(Ex_ChP_ATAC, new.names = paste0(substr(colnames(Ex_ChP_ATAC), 1, nchar(colnames(Ex_ChP_ATAC))-6), Ex_ChP_ATAC@meta.data$Lane))
colnames(Ex_ChP_ATAC)
Fragments(Ex_ChP_ATAC@assays$ATAC) <- CreateFragmentObject(path = "/Ex_ChP/ATAC_Frags/ChP11_atac_fragments.tsv.gz", cells = colnames(Ex_ChP_ATAC), validate.fragments = TRUE)
