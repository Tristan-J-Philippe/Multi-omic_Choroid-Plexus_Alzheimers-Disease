###
### 035 Cell_type_annotation_differential_ATAC.R 
###
# Purpose: Cell type annotation differential analysis and gsea of single nuclei ATAC dataset.
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_ATAC" #output results to specific folder
set.seed(42)
pcutoff <- 0.05
LFcutoff <- 0.25

###Subcelltype
####Epi

###Heatmap by pseudobulk
cts_fun <- function(celltype){
  sub <- subset(Ex_ChP_ATAC, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projid", assay = "ACTIVITY", slot = "counts", return.seurat = FALSE)
  cts <- cts[["ACTIVITY"]]
  assign(c(paste("cts_ATAC", celltype, sep = "_")), cts, envir = .GlobalEnv)
}

cts_fun("Epi_1")
gc()
cts_fun("Epi_2")
gc()
cts_fun("Epi_3")
gc()
cts_fun("Epi_4")
gc()

dim(cts_Epi_1)
dim(cts_Epi_2)
dim(cts_Epi_3)
dim(cts_Epi_4)

cts_sub <- cbind(cts_ATAC_Epi_1, cts_ATAC_Epi_2, cts_ATAC_Epi_3, cts_ATAC_Epi_4)
colnames(cts_sub) <- c(rep("Epi", 6), rep("Epi_2", 6), rep("Epi_3", 6), rep("Epi_4", 6))

celltype <- factor(paste(colnames(cts_sub)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub <- cts_sub[,colnames(cts_sub)%in%rownames(Maindm)]
cts_sub <- cts_sub[,colnames(cts_sub)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Epi_1=celltypeEpi-(celltypeEpi_2+celltypeEpi_3+celltypeEpi_4)/3, 
                        Epi_2=celltypeEpi_2-(celltypeEpi+celltypeEpi_3+celltypeEpi_4)/3, 
                        Epi_3=celltypeEpi_3-(celltypeEpi_4+celltypeEpi_2+celltypeEpi_4)/3, 
                        Epi_4=celltypeEpi_4-(celltypeEpi+celltypeEpi_2+celltypeEpi_3)/3)

cts_sub <- cts_sub[rowSums(cts_sub)>2*ncol(cts_sub),]
hvps <- FindVariableFeatures(object = cts_sub, selection.method = "vst")
plot(sort(hvps$vst.variance.standardized, decreasing = TRUE), type = "l", ylim=c(1,2))
hvps <- hvps %>% rownames_to_column("gene")
hvps <- hvps[hvps$vst.variance.standardized>1,]
cts_sub <- cts_sub[rownames(cts_sub) %in% hvps$gene,]

vdat <- voom(cts_sub, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


####BAM

###Heatmap by pseudobulk
cts_fun <- function(celltype){
  sub <- subset(Ex_ChP_ATAC, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projid", assay = "ACTIVITY", slot = "counts", return.seurat = FALSE)
  cts <- cts[["ACTIVITY"]]
  assign(c(paste("cts_ATAC", celltype, sep = "_")), cts, envir = .GlobalEnv)
}

cts_fun("BAM_1")
cts_fun("BAM_2")

dim(cts_BAM_1)
dim(cts_BAM_2)

cts_sub <- cbind(cts_ATAC_BAM_1, cts_ATAC_BAM_2)
colnames(cts_sub) <- c(rep("BAM_1", 6), rep("BAM_2", 6))

celltype <- factor(paste(colnames(cts_sub)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub <- cts_sub[,colnames(cts_sub)%in%rownames(Maindm)]
cts_sub <- cts_sub[,colnames(cts_sub)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, BAM_1=celltypeBAM_1-celltypeBAM_2, BAM_2=celltypeBAM_2-celltypeBAM_1)

cts_sub <- cts_sub[rowSums(cts_sub)>2*ncol(cts_sub),]
hvps <- FindVariableFeatures(object = cts_sub, selection.method = "vst")
plot(sort(hvps$vst.variance.standardized, decreasing = TRUE), type = "l", ylim=c(1,2))
hvps <- hvps %>% rownames_to_column("gene")
hvps <- hvps[hvps$vst.variance.standardized>1,]
cts_sub <- cts_sub[rownames(cts_sub) %in% hvps$gene,]

vdat <- voom(cts_sub, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Majorcelltype

###Heatmap by pseudobulk
cts_fun <- function(celltype){
  sub <- subset(Ex_ChP_ATAC, Majorcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projid", assay = "ACTIVITY", slot = "counts", return.seurat = FALSE)
  cts <- cts[["ACTIVITY"]]
  assign(c(paste("cts_ATAC", celltype, sep = "_")), cts, envir = .GlobalEnv)
}

cts_fun("Epithelial")
cts_fun("Fibroblast")
cts_fun("Immune")

dim(cts_Epithelial)
dim(cts_Fibroblast)
dim(cts_Immune)

cts_major <- cbind(cts_ATAC_Epithelial, cts_ATAC_Fibroblast, cts_ATAC_Immune)
colnames(cts_major) <- c(rep("Epithelial", 6), rep("Fibroblast", 6), rep("Immune", 6))

celltype <- factor(paste(colnames(cts_major)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_major)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_major <- cts_major[,colnames(cts_major)%in%rownames(Maindm)]
cts_major <- cts_major[,colnames(cts_major)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, Epithelial=celltypeEpithelial-(celltypeFibroblast+celltypeImmune)/2,
                        Fibroblast=celltypeFibroblast-(celltypeEpithelial+celltypeImmune)/2,
                        Immune=celltypeImmune-(celltypeEpithelial+celltypeFibroblast)/2
)

cts_major <- cts_major[rowSums(cts_major)>2*ncol(cts_major),]
hvps <- FindVariableFeatures(object = cts_major, selection.method = "vst")
plot(sort(hvps$vst.variance.standardized, decreasing = TRUE), type = "l", ylim=c(1,2))
hvps <- hvps %>% rownames_to_column("gene")
hvps <- hvps[hvps$vst.variance.standardized>1,]
cts_major <- cts_major[rownames(cts_major) %in% hvps$gene,]

vdat <- voom(cts_major, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###FGSEA celltype
Path_fun(celltypeDEGs_Epi_1)
Path_fun(celltypeDEGs_Epi_2)
Path_fun(celltypeDEGs_Epi_3)
Path_fun(celltypeDEGs_Epi_4)
Path_fun(celltypeDEGs_BAM_1)
Path_fun(celltypeDEGs_BAM_2)

Path_fun(celltypeDEGs_Epithelial)
Path_fun(celltypeDEGs_Fibroblast)
Path_fun(celltypeDEGs_Immune)



###In RNA
####Load RNA

RNA_Epi_1_corrected <- read_csv("/Ex_RNA/Epi_1/Spreadsheet__Epi_1_corrected.csv")
RNA_Epi_2_corrected <- read_csv("/Ex_RNA/Epi_2/Spreadsheet__Epi_2_corrected.csv")
RNA_Epi_3_corrected <- read_csv("/Ex_RNA/Epi_3/Spreadsheet__Epi_3_corrected.csv")
RNA_Epi_4_corrected <- read_csv("/Ex_RNA/Epi_4/Spreadsheet__Epi_4_corrected.csv")

RNA_Fib_1_corrected <- read_csv("/Ex_RNA/Fib_1/Spreadsheet__Fib_1_corrected.csv")
RNA_Fib_2_corrected <- read_csv("/Ex_RNA/Fib_2/Spreadsheet__Fib_2_corrected.csv")
RNA_Fib_3_corrected <- read_csv("/Ex_RNA/Fib_3/Spreadsheet__Fib_3_corrected.csv")
RNA_Fib_4_corrected <- read_csv("/Ex_RNA/Fib_4/Spreadsheet__Fib_4_corrected.csv")

RNA_BAM_1_corrected <- read_csv("/Ex_RNA/BAM_1/Spreadsheet__BAM_1_corrected.csv")
RNA_BAM_2_corrected <- read_csv("/Ex_RNA/BAM_2/Spreadsheet__BAM_2_corrected.csv")

RNA_Epithelial_corrected <- read_csv("/Ex_RNA/Epithelial/Spreadsheet__Epithelial_corrected.csv")
RNA_Fibroblast_corrected <- read_csv("/Ex_RNA/Fibroblast/Spreadsheet__Fibroblast_corrected.csv")
RNA_Mural_corrected <- read_csv("/Ex_RNA/Mural/Spreadsheet__Mural_corrected.csv")
RNA_Immune_corrected <- read_csv("/Ex_RNA/Immune/Spreadsheet__Immune_corrected.csv")
RNA_Parenchyma_corrected <- read_csv("/Ex_RNA/Parenchyma/Spreadsheet__Parenchyma_corrected.csv")


####Overlap venns

RNA_atac_overlap_fun(RNA_Epi_1_corrected, celltypeDEGs_Epi_1_sig)
RNA_atac_overlap_fun(RNA_Epi_2_corrected, celltypeDEGs_Epi_2_sig)
RNA_atac_overlap_fun(RNA_Epi_3_corrected, celltypeDEGs_Epi_3_sig)
RNA_atac_overlap_fun(RNA_Epi_4_corrected, celltypeDEGs_Epi_4_sig)

RNA_atac_overlap_fun(RNA_Fib_1_corrected, celltypeDEGs_Fib_1_sig)
RNA_atac_overlap_fun(RNA_Fib_2_corrected, celltypeDEGs_Fib_2_sig)
RNA_atac_overlap_fun(RNA_Fib_3_corrected, celltypeDEGs_Fib_3_sig)
RNA_atac_overlap_fun(RNA_Fib_4_corrected, celltypeDEGs_Fib_4_sig)

RNA_atac_overlap_fun(RNA_BAM_1_corrected, celltypeDEGs_BAM_1_sig)
RNA_atac_overlap_fun(RNA_BAM_2_corrected, celltypeDEGs_BAM_2_sig)


RNA_atac_overlap_fun(RNA_Epithelial_corrected, celltypeDEGs_Epithelial_sig)
RNA_atac_overlap_fun(RNA_Fibroblast_corrected, celltypeDEGs_Fibroblast_sig)
RNA_atac_overlap_fun(RNA_Mural_corrected, celltypeDEGs_Mural_sig)
RNA_atac_overlap_fun(RNA_Immune_corrected, celltypeDEGs_Immune_sig)


##Publish figures
###UMAP

Ex_ChP_ATAC <- RunUMAP(Ex_ChP_ATAC, dims=c(2:20), assay="ATAC", reduction="harmony", min.dist=2.5)
#Major types
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Majorcelltype
levels(Ex_ChP_ATAC) <- c("Epithelial", "Fibroblast", "Immune")
Ex_ChP_ATAC@meta.data[["Majorcelltype"]] <- Idents(Ex_ChP_ATAC)
plts <- DimPlot(Ex_ChP_ATAC, raster=F, raster.dpi=c(600,600), cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Immune"="#26532B")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_ATAC_UMAP_major.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_ATAC_UMAP_major.png", sep=""), width=10, height=5)
remove(plts)


#Subcelltypes
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
levels(Ex_ChP_ATAC) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b", "Fibroblast", "BAM_1", "BAM_2")
Ex_ChP_ATAC@meta.data[["Subcelltype"]] <- Idents(Ex_ChP_ATAC)
plts <- DimPlot(Ex_ChP_ATAC, raster=F, raster.dpi=c(600,600), cols=c("Epi_1"="#CD8500", "Epi_2"="#FFA500", "Epi_3"="#FF4500", "Epi_4"="#CD3700", "Fibroblast"="#3249A6", "BAM_1"="#70C25B", "BAM_2"="#0E9554"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP.png", sep=""), width=10, height=5)
remove(plts)


#Sub type individual figures
Idents(Epi) <- Epi@meta.data[["Subcelltype"]] 
levels(Epi) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4")
Epi@meta.data[["Subcelltype"]] <- Idents(Epi)
Epi <- RunUMAP(Epi, dims=c(2:20), assay="ATAC", reduction="harmony", min.dist=2.5)
plts <- DimPlot(Epi, raster=F, cols=c("Epi_1"="#CD8500", "Epi_2"="#FFA500", "Epi_3"="#FF4500", "Epi_4"="#CD3700"), pt.size=3) + theme(text=element_text(size=7)) #make publishable
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Epi_UMAP.png", sep=""), width=10, height=5)
remove(plts)

Idents(Imm) <- Imm@meta.data[["Subcelltype"]]
levels(Imm) <- c("BAM_1", "BAM_2")
Imm@meta.data[["Subcelltype"]] <- Idents(Imm)
plts <- DimPlot(Imm, raster=F, cols=c("BAM_1"="#70C25B", "BAM_2"="#0E9554"), pt.size=4) + theme(text=element_text(size=7)) # make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Imm_UMAP.png", sep=""), width=10, height=5)
remove(plts)


###Heatmap

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
gc()

#select list
Epi_genes <- c("AQP1", "PRLR", "ABCA4", "ABCA5", "GALNT11", "ENPP2", "TPD52L1", "TRPV4", "GPM6A", "SFRP1", "SGK1", "NWD1", "CDH12", "FBLN1", "GLIS3", "HTR2C", "OTX2-AS1", "TRMT9B", "VAT1L", "WDR49", "ZBBX", "CLIC6", "ATP2B2", "STEAP1B", "SLC5A3", "DNAH12", "DNAH9", "DNAH6", "DNAH7", "DNAH11", "RMST", "ROR2", "RFX3", "ARL13B", "SPAG17", "TMEM72-AS1", "TMEM117", "RYR3", "SLC24A4", "SNTB1", "SLC12A2", "SLC39A12", "EYA4", "F5", "AGBL4", "IGFBP4", "RIMBP2", "FYB2", "SLC4A10", "CFAP43", "CRPPA", "OCA2", "GNA14", "NEK11", "CFAP54", "MTUS2", "CFAP44", "DGKI", "GPR155", "DOCK8", "POU2AF1", "PACRG", "HYDIN", "PRDM16", "SHANK2", "CFAP100", "TTR", "ESRRG", "EPHA4", "LRRIQ1", "CHST11", "COL9A1", "CTNND2", "CUX2", "DCDC1", "LRRC9", "MRPS6", "SDK1", "MYO1E", "SLC15A4", "SOX6", "SPAG16")
Epi34_genes <- c("MT-CO1", "MT-ND4", "GPX3", "MT-CO1", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB", "MT-ND1", "MT-ND2", "SLC26A3", "MT-ND3", "MT-ND4L", "MT-ND5", "MTRNR2L8")
Fib_genes <- c("PPM1H", "PDE7B", "PBX1", "NEBL", "ITGA9", "KLF9", "KCNT2", "ITGA8", "LRP1", "LRP1B", "PLCB1", "SESTD1", "SETBP1", "MID1", "GPC6", "MIR100HG", "GMDS", "PHLDB2", "CEMIP", "ABCA9", "APOD", "FBN1", "FN1", "EDA", "EBF1", "EGFR", "EYA2", "FAM20A", "FOXP2", "ADAMTSL3", "ARHGAP10", "EGFR", "ATP13A3", "ARHGAP21", "CDH11", "ADK", "SLC22A23", "BCL6", "THSD4", "TGFBR3", "BICC1", "BNC2", "CHSY3", "SLC7A2", "CLDN11", "ATP1B3", "ABCA8", "SPAG9", "BMPR1B", "AKAP12", "ADAM12", "ADD3", "ALCAM", "ALPL", "ARHGAP6", "ATP13A3", "WWTR1", "TJP1", "TJP2", "TMOD1", "UNC5C", "TMTC1", "USP53", "CADPS2", "CDH23", "CDK14", "CLIC4", "COL4A2", "COL11A1", "CPED1", "EPS8", "GLCE", "GNAL", "GRK5", "HIPK2", "HIVEP2", "LAMA4", "LDLRAD3", "LPAR1", "MAP3K5", "MAST4", "NID1", "NR2F1", "NR2F1-AS1", "PLD1", "POU6F2", "PRKCH", "SLC1A1", "SLC22A3", "SLC38A2", "SMOC1", "ACADL", "JAK2", "ALDH1A2", "COL1A2", "COL3A1", "LEPR", "SLC24A3", "STXBP6", "TBC1D8")
Fib1_genes <- c("PLXNA4", "KCNMA1", "CACNA2D3", "SPTLC3", "SLC4A4", "ARL15")
Fib2_genes <- c("COCH", "NTM", "COLEC12","LAMA2")
Fib3_genes <- c("VPS41", "TRPM3", "MDGA2", "PKP2", "PRKAG2", "PTGDS", "PTN", "SLC7A14-AS1", "FAM110B", "HS3ST5", "SVIL", "GREB1", "ATP1A2", "CACNA1A", "ROBO1", "MAPK4", "KCNQ3", "SLC6A1", "SLC6A13", "SLC7A11", "SLC13A3", "KCNK2", "KCNT1", "MYRIP", "NCKAP5", "RANBP3L", "ROBO2", "SLC1A3", "ARNT2", "EPHA7", "HIF3A", "SNCAIP", "PAK3", "NKD1", "CAB39L", "COL5A1", "COL5A2", "ADCY5")
Fib4_genes <- c("ANKS1B", "ABCA10", "CBLB", "ABCA6", "P4HA1", "PDGFRA", "HIF1A", "PRKD1", "PTGES3", "RIPOR3", "UBC")
Endo_genes <- c("CCNY", "CLIC4", "MECOM", "ARHGAP26", "ARHGAP29", "ATP13A3", "PODXL", "PTPRB", "SHANK3", "THSD7A", "SLC2A3", "RPGR", "CMIP", "EPAS1", "IQSEC1", "MSN", "IVNS1ABP", "GALNT18", "PLEKHG1", "PLPP1", "PECAM1", "CRIM1", "FLT1", "BMPR2", "GRB10", "DOCK9", "UBC", "INSR", "LDB2", "NFIB", "ARL15", "NOTCH4", "LYVE1", "PIK3R3", "STGALNAC3", "UGCG", "ADGRL4", "ELOVL7", "PRUNE2")
Mural_genes <- c("RBPMS", "ACTA2", "RCAN2", "RGS6", "MYH11", "CARMN", "CALD1", "SYNPO2", "DLC1", "SGIP1", "MCAM", "LTBP1", "NOTCH3", "ADAMTS1", "FRY", "JAG1", "LMOD1", "COL4A1", "COL4A2", "EBF1", "EPS8", "LPP", "SLIT3")
BAM_genes <- c("STAB1", "MERTK", "CD163", "AIF1", "CSF1R", "APBB1IP", "C1QA", "C1QB", "C1QC", "CLEC7A", "ENPP2", "MERTK", "MS4A4E", "MS4A6A", "MS4A7", "PTPRC")
BAM1_genes <- c("ABCA1", "AOAH", "STARD13", "ADAM28", "MEF2A", "MERTK", "MGAT4A", "MRC1", "PDE4B", "PLCG2", "RAB31", "SFMBT2", "SLC1A3", "SPRED1", "SRGN")
BAM2_genes <- c("DNAJB6", "HSP90AA1", "HSPH1", "DNAJB1", "HSPD1", "HSPE1", "DNAJA4", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA4L", "HSPA5", "HSPA6", "HSPA8", "HSPB1", "ATP2C1", "CD83", "CHORDC1", "P4HA1", "PTGES3", "SPECC1", "SPP1", "TRA2B", "UBC", "UGCG", "USP36", "USPL1")
BAM3_genes <- c("MKI67", "ANLN", "ANXA2", "APOLD1", "ARHGAP11", "ARHGAP1B", "ASPM", "ATAD2", "ATAD5", "BARD1", "BRCA1", "BRCA2", "BRIP1", "BUB1", "BUB1B", "CCDC18", "CIP2A", "CIT", "CKAP2", "CKAP2L", "CLSPN", "DEPDC1B", "DIAPH3", "DTL", "ECT2", "EZH2", "FAM111B", "FANCI", "GINS1", "GPSM2", "GTSE1", "HJURP", "HMMR", "KIF4A", "KIF11", "KIF14", "KIF15", "KIF18A", "KIF18B", "KIF20B", "KIF23", "KIFC1", "KNL1", "KNTC1", "LMNB1", "LNCAROD", "MELK", "NCAPG", "NCAPG2", "MCAPH", "NUF2", "NUSAP1", "PARPBP", "POLQ", "PRR1", "PVT1", "RFC3", "RRM2", "SCLT1", "SGO2", "TOP2A", "SMC4", "TPX2", "VRK1")
Tcell_genes <- c("THEMIS", "IL7R", "CANK4", "CBLB", "BCL11B", "ANK3", "ETS1", "FYN", "GLCCI1", "INPP4B", "IQGAP2", "ITGA4", "ITK", "KLF12", "PARP8", "PRKCH", "PPP1R16B", "PPP2R2B", "PRKACB", "PTPN22", "RASA3", "RESF1", "SKAP1", "SLC38A1", "SLFN12L", "STAT4", "SYNE2", "SYTL3", "TC2N", "TOX", "TRERF1", "TSPAN5", "ZEB1", "ZNF831")
paren_genes <- c("MOBP", "NAV2", "GFAP", "CLU", "SLC1A3", "FGF14", "FRMD5", "LSAMP", "GLUT1", "ADGRB3", "ADGRL3", "ANK2", "NRXN1", "PCDH9", "SLC9A9", "ANK3", "KCNIP4", "NRXN3", "RBFOX1", "NRG3", "CADM2")

###Heatmap by pseudobulk
cts_fun <- function(celltype){
  sub <- subset(Ex_ChP_ATAC, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projid", assay = "ACTIVITY", slot = "counts", return.seurat = FALSE)
  cts <- cts[["ACTIVITY"]]
  assign(c(paste("cts", celltype, sep = "_")), cts, envir = .GlobalEnv)
}

cts_fun("Epi_1")
gc()
cts_fun("Epi_2")
gc()
cts_fun("Epi_3")
gc()
cts_fun("Epi_4")
gc()

cts_fun("BAM_1")
cts_fun("BAM_2")

cts_fun <- function(celltype){
  sub <- subset(Ex_ChP_ATAC, Majorcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projid", assay = "ACTIVITY", slot = "counts", return.seurat = FALSE)
  cts <- cts[["ACTIVITY"]]
  assign(c(paste("cts", celltype, sep = "_")), cts, envir = .GlobalEnv)
}

cts_fun("Epithelial")
cts_fun("Fibroblast")
cts_fun("Immune")

#major no paren
cts <- cbind(cts_Epithelial, cts_Fibroblast, cts_Immune)
#curated findmarkers
my_genes <- c(Epi_genes, Fib_genes, BAM_genes) %>% unique()

cts <- data.frame(scale(cts))
cts <- cts[rownames(cts) %in% my_genes,] #!!!!no parenchyma in my genes
cts <- cts[order(match(rownames(cts), my_genes)),]
cts <- t(scale(t(cts)))
htg <- Heatmap(cts, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_major_no_paren.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_major_no_paren.pdf", sep = ""), width=7, height=7)


#sub with paren
cts <- cbind(cts_Epi_1, cts_Epi_2, cts_Epi_3, cts_Epi_4, cts_Fibroblast, cts_BAM_1, cts_BAM_2)
cts <- data.frame(scale(cts))

my_genes <- c(Epi_genes, Epi34_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, Fib4_genes, BAM_genes, BAM1_genes, BAM2_genes, BAM3_genes, Tcell_genes) %>% unique()

cts <- data.frame(scale(cts))
cts <- cts[rownames(cts) %in% my_genes,] #!!!!no parenchyma in my genes
cts <- cts[order(match(rownames(cts), my_genes)),]
cts <- t(scale(t(cts)))
htg <- Heatmap(cts, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_noparen.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_noparen.pdf", sep = ""), width=7, height=7)



###Bubble VLN Feature QC plots
####Main

DefaultAssay(Ex_ChP_ATAC) <- "ACTIVITY"
highlights_genes <- c("TTR", "HTR2C", "GPX3", "DCN", "LEPR", "FN1", "STAB1", "CD163", "MSR1")
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Majorcelltype
plts <- DotPlot(Ex_ChP_ATAC, features = highlights_genes, split.by="Majorcelltype", assay="ACTIVITY", cols=c("#C48F45", "#3249A6", "#26532B")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_ATAC_bubble.pdf", sep=""), width=5+length(highlights_genes)*.15, height=1.5+length(unique(Ex_ChP_ATAC@meta.data[["Majorcelltype"]]))*.4)

Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Majorcelltype
plts <- DotPlot(Ex_ChP_ATAC, features = highlights_genes, assay="ACTIVITY") + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_ATAC_bubble_scale.pdf", sep=""), width=5+length(highlights_genes)*.15, height=1.5+length(unique(Ex_ChP_ATAC@meta.data[["Majorcelltype"]]))*.4)

plts <- VlnPlot(Ex_ChP_ATAC, features = highlights_genes, assay="ACTIVITY", split.by='Majorcelltype', stack=T, raster=T, cols=c("#C48F45", "#3249A6", "#26532B"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_ATAC_vln.pdf", sep=""), width=3+length(highlights_genes)*.4, height=1.5+length(unique(Ex_ChP_ATAC@meta.data[["Majorcelltype"]]))*1)

Featurefun(Ex_ChP_ATAC, highlights_genes, "highlights_genes",  assay="ATAC", max.cutoff=c(3,4))


####Epi sub

Epi <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Epithelial.RDS"))

DefaultAssay(Epi) <- "ACTIVITY"
Epi_genes <- c("TTR", "HTR2C", "CDH12", "OTX2", "JUND", "CDH1", "KRT8", "KRT18", "LHX2", "GFAP", "KCNMA1", "FOXJ1", "DNAH6", "CFAP43", "CFAP54", "SPAG17", "TUBB4B", "CLDN1", "CLDN3", "CLDN5", "OCLN", "TJP1", "TJP2", "TRPM3", "GPX3", "MT-CO1", "MT-CO2", "MT-CO3", "AQP1", "ATP1A1", "SLC2A12", "SLC25A6", "SLC25A23", "SLC22A17", "SLC7A2", "SLC25A3", "SLC4A2", "SLC4A4", "FOLR1", "PIFO", "IGFBP2", "LUM", "VIM", "S100B", "CD2")

plts <- DotPlot(Epi, features = Epi_genes, split.by="Subcelltype", cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble.pdf", sep=""), width=3+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Epi, features = Epi_genes, cols=c("grey", "#C48F45")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble_scale.pdf", sep=""), width=5+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Epi, features = Epi_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble_default.pdf", sep=""), width=5+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Epi, features = Epi_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700"), same.y.lims=F) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_vln.pdf", sep=""), width=3+length(Epi_genes)*.4, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*1)
remove(plts)

plts <- VlnPlot(Epi, features = Epi_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))  + geom_boxplot(outliers=F)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_vln_box.pdf", sep=""), width=3+length(Epi_genes)*.4, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Epi, Epi_genes, "highlights_genes",  assay="ACTIVITY", max.cutoff=c(3,4))


####Imm sub

Imm <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Immune.RDS"))

DefaultAssay(Imm) <- "ACTIVITY"
Imm_genes <- c("F13A1", "LYVE1", "MRC1", "STAB1", "CD163", "CD68", "MERTK", "HCK", "IL4R", "CD74", "HLA-DRA", "S100A9", "CSF1R", "JAK2", "CD86", "VCAM1", "CD83", "LAMP1", "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1", "MKI67", "CENPF", "TOP2A", "ITGAM", "S100A8", "CCR2", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11", "THEMIS", "CD3E", "PDE7A", "CAMK4", "HLA-B", "CD8A") %>% unique()

Idents(Imm) <- Imm@meta.data$Subcelltype
plts <- DotPlot(Imm, features = Imm_genes, split.by="Subcelltype", cols=c("#70C25B", "#0E9554")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Imm, features = Imm_genes, cols=c("grey", "#26532B")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble_scale.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Imm, features = Imm_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble_default.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Imm, features = Imm_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#70C25B", "#0E9554"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_vln.pdf", sep=""), width=3+length(Imm_genes)*.4, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Imm, Imm_genes, "highlights_genes",  assay="ACTIVITY", max.cutoff=c(3,4))


###Correlation between sn RNAseq and ATAC
####major

my_genes <- c(sn_Epithelial_corrected$gene, sn_Fibroblast_corrected$gene, sn_Immune_corrected$gene) %>% unique()

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Epithelial_RNA", "Fibroblast_RNA", "Endothelial_RNA", "Mural_RNA", "Immune_RNA", "Parenchyma_RNA")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]
Ex_cts$Endothelial_RNA <- NULL
Ex_cts$Mural_RNA <- NULL
Ex_cts$Parenchyma_RNA <- NULL

atac_cts <- data.frame(AggregateExpression(Ex_ChP_ATAC, group.by = "Majorcelltype", assay = "ACTIVITY", slot = "data", return.seurat = FALSE))
colnames(atac_cts) <- c("Epithelial_ATAC", "Fibroblast_ATAC", "Immune_ATAC")
atac_cts <- data.frame(t(scale(t(atac_cts))))
atac_cts <- atac_cts %>% rownames_to_column("gene")
atac_cts <- atac_cts[atac_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, atac_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epithelial_ATAC", "Epithelial_RNA", "Fibroblast_ATAC", "Fibroblast_RNA", "Immune_ATAC", "Immune_RNA")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("ATAC", rownames(cts_cor)), grepl("RNA", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Major.pdf", sep = ""), width=7, height=7)


####within Epi

my_genes <- c(sn_Epithelial_corrected$gene, sn_Epi_1_corrected$gene, sn_Epi_2_corrected$gene, sn_Epi_3_corrected$gene, sn_Epi_4_corrected$gene) %>% unique()

#within overlapping genes
Ex_Epi_cts <- data.frame(AggregateExpression(Epi_sn, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_Epi_cts) <- c("Epi_1_SN", "Epi_2_SN", "Epi_3_SN", "Epi_4_SN")
Ex_Epi_cts <- data.frame(t(scale(t(Ex_Epi_cts))))
Ex_Epi_cts <- Ex_Epi_cts %>% rownames_to_column("gene")
Ex_Epi_cts[is.na(Ex_Epi_cts)] <- 0
Ex_Epi_cts <- Ex_Epi_cts[Ex_Epi_cts$gene %in% my_genes,]

atac_Epi_cts <- data.frame(AggregateExpression(Epi, group.by = "Subcelltype", assay = "ACTIVITY", slot = "data", return.seurat = FALSE))
colnames(atac_Epi_cts) <- c("Epi_1_ATAC", "Epi_2_ATAC", "Epi_3_ATAC", "Epi_4_ATAC")
atac_Epi_cts <- data.frame(t(scale(t(atac_Epi_cts))))
atac_Epi_cts <- atac_Epi_cts %>% rownames_to_column("gene")
atac_Epi_cts <- atac_Epi_cts[atac_Epi_cts$gene %in% my_genes,]

cts <- left_join(Ex_Epi_cts, atac_Epi_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_1_ATAC", "Epi_1_SN", "Epi_2_ATAC", "Epi_2_SN", "Epi_3_ATAC", "Epi_3_SN", "Epi_4_ATAC", "Epi_4_SN")]
cts <- cts[complete.cases(cts), ]

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("ATAC", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
Heatmap(cts_cor)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Epi.pdf", sep = ""), width=7, height=7)


####within BAM

my_genes <- c(sn_Immune_corrected$gene, sn_BAM_1_corrected$gene, sn_BAM_2_corrected$gene) %>% unique()

#within overlapping genes
Ex_Imm_cts <- data.frame(AggregateExpression(Imm_sn, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_Imm_cts) <- c("BAM_1_SN", "BAM_2_SN")
Ex_Imm_cts <- data.frame(t(scale(t(Ex_Imm_cts))))
Ex_Imm_cts <- Ex_Imm_cts %>% rownames_to_column("gene")
Ex_Imm_cts[is.na(Ex_Imm_cts)] <- 0
Ex_Imm_cts <- Ex_Imm_cts[Ex_Imm_cts$gene %in% my_genes,]

atac_Imm_cts <- data.frame(AggregateExpression(Imm, group.by = "Subcelltype", assay = "ACTIVITY", slot = "data", return.seurat = FALSE))
colnames(atac_Imm_cts) <- c("BAM_1_ATAC", "BAM_2_ATAC")
atac_Imm_cts <- data.frame(t(scale(t(atac_Imm_cts))))
atac_Imm_cts <- atac_Imm_cts %>% rownames_to_column("gene")
atac_Imm_cts <- atac_Imm_cts[atac_Imm_cts$gene %in% my_genes,]

cts <- left_join(Ex_Imm_cts, atac_Imm_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("BAM_1_ATAC", "BAM_1_SN", "BAM_2_ATAC", "BAM_2_SN")]
cts <- cts[complete.cases(cts), ]

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("ATAC", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
Heatmap(cts_cor)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_BAM.pdf", sep = ""), width=7, height=7)



###CoveragePlot of select markers
####Epithelial

DefaultAssay(Ex_ChP_ATAC) <- "ATAC"
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC$Subcelltype

Epi_genes <- c("TTR", "HTR2C", "CDH12", "OTX2", "JUND", "CDH1", "KRT8", "KRT18", "LHX2", "GFAP", "KCNMA1", "FOXJ1", "DNAH6", "CFAP43", "CFAP54", "CFAP45", "CFAP47", "CFAP58", "CFAP69", "CFAP70", "SPAG17", "TUBB4B", "CLDN1", "CLDN3", "CLDN5", "OCLN", "TJP1", "TJP2", "TRPM3", "GPX3", "AQP1", "ATP1A1", "SLC2A12", "SLC25A6", "SLC25A23", "SLC22A17", "SLC7A2", "SLC25A3", "SLC4A2", "SLC4A4", "FOLR1", "PIFO", "IGFBP2", "LUM", "VIM", "S100B", "CD2")

test <- celltypeDEGs_Epithelial_sig$gene
closest_totest <- ClosestFeature(Ex_ChP_ATAC, regions=test)
closest_totest <- left_join(closest_totest, Gene.search, join_by("gene_name"=="gene"))
view(closest_totest)
closest_totest[closest_totest$gene_name %in% Epi_genes,]

for(gene in Epi_genes){
  highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
  plts <- CoveragePlot(Ex_ChP_ATAC, region=gene, region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#3249A6", "#70C25B", "#0E9554"))
  plot(plts)
  graph2pdf(plts, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/", gene, ".pdf"))
}

gene <- "CFAP100"
highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
plts <- CoveragePlot(Ex_ChP_ATAC, region="chr3-126394939-126394978", region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#3249A6", "#70C25B", "#0E9554"))
plts
graph2pdf(plts, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/", gene, "_sub.pdf"))

gene <- "CFAP69"
highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
plts <- CoveragePlot(Ex_ChP_ATAC, region="chr7-90245545-90255422", region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#3249A6", "#70C25B", "#0E9554"))
plts
graph2pdf(plts, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/", gene, "_sub.pdf"))


####Fibroblast

Fibro_genes <- c("DCN", "FN1", "LEPR", "CDH11", "STXBP6", "FBLN1", "KCNMA1", "SLC4A4", "PLXNA4", "ALPL", "DPEP1", "DLK1", "NRG1", "ROBO1", "TRPM3", "MYRIP", "GPM6A", "CLDN11", "ECM2", "SLC6A20", "SLC13A3", "ABCA9", "LAMA2", "LAMA4", "COL18A1", "COL4A2", "COL1A2", "COL4A4", "COL5A2",  "COL3A1", "COL6A3", "PDGFRA", "SLC4A7", "SLC2A3", "LUM")

test <- celltypeDEGs_Fibroblast_sig$gene
closest_totest <- ClosestFeature(Ex_ChP_ATAC, regions=test)
closest_totest <- left_join(closest_totest, Gene.search, join_by("gene_name"=="gene"))
view(closest_totest)
closest_totest[closest_totest$gene_name %in% Fibro_genes,]

for(gene in Fibro_genes){
  highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
  plts <- CoveragePlot(Ex_ChP_ATAC, region=gene, region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#3249A6", "#70C25B", "#0E9554"))
  plot(plts)
  graph2pdf(plts, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/", gene, ".pdf"))
}

gene <- "SLC4A4"
highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
plts <- CoveragePlot(Ex_ChP_ATAC, region="chr4-71186757-71187401", region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#3249A6", "#70C25B", "#0E9554"))
plts
graph2pdf(plts, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/", gene, "_sub.pdf"))


####Immune

Imm_genes <- c("F13A1", "LYVE1", "MRC1", "STAB1", "CD163", "CD68", "MERTK", "HCK", "IL4R", "CD74", "HLA-DRA", "S100A9", "CSF1R", "JAK2", "CD86", "VCAM1", "CD83", "LAMP1", "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1", "MKI67", "CENPF", "TOP2A", "ITGAM", "S100A8", "CCR2", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11", "THEMIS", "CD3E", "PDE7A", "CAMK4", "HLA-B", "CD8A") %>% unique()

test <- celltypeDEGs_Immune_sig$gene
closest_totest <- ClosestFeature(Ex_ChP_ATAC, regions=test)
closest_totest <- left_join(closest_totest, Gene.search, join_by("gene_name"=="gene"))
view(closest_totest)
closest_totest[closest_totest$gene_name %in% Imm_genes,]

for(gene in Imm_genes){
  highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
  plts <- CoveragePlot(Ex_ChP_ATAC, region=gene, region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#3249A6", "#70C25B", "#0E9554"))
  plot(plts)
  graph2pdf(plts, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/", gene, ".pdf"))
}



###Transcription factors and motifs

#!!!!Ranges previously included non-standard contigs that need to removed for AddMotifs to work
peaks.keep <- seqnames(granges(Ex_ChP_ATAC)) %in% standardChromosomes(granges(Ex_ChP_ATAC))
Ex_ChP_ATAC_filt <- Ex_ChP_ATAC
Ex_ChP_ATAC_OG <- Ex_ChP_ATAC
Ex_ChP_ATAC_filt@assays$ATAC <- Ex_ChP_ATAC_filt@assays$ATAC[as.vector(peaks.keep), ]
DefaultAssay(Ex_ChP_ATAC) <- "ACTIVITY"
Ex_ChP_ATAC[["ATAC"]] <- NULL
Ex_ChP_ATAC[["ATAC"]] <- CreateChromatinAssay(counts = Ex_ChP_ATAC_filt@assays[["ATAC"]], sep = c(":", "-"), fragments = Ex_ChP_ATAC_OG@assays[["ATAC"]]@fragments[[1]], annotation=annotation)
DefaultAssay(Ex_ChP_ATAC) <- "ATAC"
rm(Ex_ChP_ATAC_filt, Ex_ChP_ATAC_OG)

library(motifmatchr)
pwm <- getMatrixSet(JASPAR2020, opts=list(species=9606, all_versions=FALSE))
Ex_ChP_ATAC <- AddMotifs(Ex_ChP_ATAC, genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pwm)
library(BSgenome.Hsapiens.UCSC.hg38)
saveRDS(Ex_ChP_ATAC, paste(curr_dir, "/CellType/Finalcelltype/Ex_ChP_ATAC_noextrachromosomes.RDS", sep=""))


####Differential motifs?

DefaultAssay(Ex_ChP_ATAC) <- "ATAC"
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC@meta.data$Subcelltype
Ex_ChP_ATAC_enriched.motifs <- FindMotifs(Ex_ChP_ATAC, features=unique(celltypeDEGs_Epithelial$gene), only.pos=F, verbose=T)


#TF footprinting
Ex_ChP_ATAC <- Footprint(Ex_ChP_ATAC, motif.name ="MA0599.1", genome = BSgenome.Hsapiens.UCSC.hg38, compute.expected = FALSE)
PlotFootprint(Ex_ChP_ATAC, features = c("MA0599.1"))
Ex_ChP_ATAC <- Footprint(Ex_ChP_ATAC, motif.name ="MA0481.3", genome = BSgenome.Hsapiens.UCSC.hg38, compute.expected = FALSE)
PlotFootprint(Ex_ChP_ATAC, features = c("MA0481.3"))
Ex_ChP_ATAC <- Footprint(Ex_ChP_ATAC, motif.name ="MA1648.1", genome = BSgenome.Hsapiens.UCSC.hg38, compute.expected = FALSE)
PlotFootprint(Ex_ChP_ATAC, features = c("MA1648.1"))

#Not different
for(i in Ex_ChP_ATAC_enriched.motifs$motif.name[1:5]){
  print(i)
  Ex_ChP_ATAC <- Footprint(Ex_ChP_ATAC, motif.name =i, genome = BSgenome.Hsapiens.UCSC.hg38, compute.expected = FALSE)
  plot(PlotFootprint(Ex_ChP_ATAC, features = c(i)))
}
