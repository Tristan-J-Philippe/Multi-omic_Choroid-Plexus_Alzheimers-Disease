###
### 037 Cell_type_annotation_differential_ST.R 
###
# Purpose: Cell type annotation differential analysis and gsea of spatial transcriptomics dataset.
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder
set.seed(42)
FDRcutoff <- 0.05
LFcutoff <- 0.25
pct_cutoff <- 0.1

###Heatmap by pseudobulk
cts_fun <- function(celltype){
  sub <- subset(Cos_ChP, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projID", assay = "Nanostring", slot = "counts", return.seurat = FALSE)
  cts <- cts[["Nanostring"]]
  assign(c(paste("cts", celltype, sep = "_")), cts, envir = .GlobalEnv)
  pct <- pctincell(sub)
  assign(c(paste("pct", celltype, sep = "_")), pct, envir = .GlobalEnv)
}


###Subcelltype

for(i in names(Fibro@images)){Fibro@images[[i]]<- NULL} #remove images before subbing
cts_fun("Epi_1")
gc()
cts_fun("Epi_2")
gc()
cts_fun("Epi_3")
gc()
cts_fun("Epi_4")
gc()
cts_fun("Fib_1")
cts_fun("Fib_2")
cts_fun("Fib_3")
cts_fun("Endo_1")
cts_fun("Endo_3")
cts_fun("Endo_2")
cts_fun("Mural_2")
cts_fun("Mural_3")
cts_fun("Mural_1")
cts_fun("BAM_1")
cts_fun("BAM_2")
cts_fun("Olig")


cts_sub_paren <- cbind(cts_Epi_1, cts_Epi_2, cts_Epi_3, cts_Epi_4, cts_Fibroblast, cts_Endo_1, cts_Endo_3, cts_Endo_2, cts_Mural, cts_Immune, cts_Olig)
colnames(cts_sub_paren) <- c(rep("Epi_1", 8), rep("Epi_2", 8), rep("Epi_3", 8), rep("Epi_4", 8), rep("Fibroblast", 8), rep("Endo_1", 8), rep("Capilary", 6), rep("Endo_2", 8), rep("Mural", 8), rep("Immune", 8), rep("Olig", 5))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, Epi_1=celltypeEpi_1-(celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFibroblast+celltypeEndo+celltypeMural+celltypeImmune+celltypeOlig)/8, Epi_2=celltypeEpi_2-c(celltypeEpi_1+celltypeEpi_3+celltypeEpi_4+celltypeFibroblast+celltypeEndo+celltypeMural+celltypeImmune+celltypeOlig)/8, Epi_3=celltypeEpi_3-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_4+celltypeFibroblast+celltypeEndo+celltypeMural+celltypeImmune+celltypeOlig)/8, Epi_4=celltypeEpi_4-c(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeFibroblast+celltypeEndo+celltypeMural+celltypeImmune+celltypeOlig)/8, Fibroblast=celltypeFibroblast-(celltypeEpi_1 +celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeEndo+celltypeMural+celltypeImmune+celltypeOlig)/8,  Endo=celltypeEndo-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFibroblast+celltypeMural+celltypeImmune+celltypeOlig)/8, Mural=celltypeMural-c(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFibroblast+celltypeEndo+celltypeImmune+celltypeOlig)/8, Immune=celltypeImmune-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFibroblast+celltypeEndo+celltypeMural+celltypeOlig)/8, Olig=celltypeOlig-c(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFibroblast+celltypeEndo+celltypeMural+celltypeImmune)/8)

vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Majorcelltype

for(i in names(Cos_ChP@images)){Cos_ChP@images[[i]]<- NULL} #remove images before subbing
cts_fun("Epithelial")
cts_fun("Fibroblast")
cts_fun("Endothelial")
cts_fun("Mural")
cts_fun("Immune")
cts_fun("Parenchyma")

cts_major_paren <- cbind(cts_Epithelial, cts_Fibroblast, cts_Endothelial, cts_Mural, cts_Immune, cts_Parenchyma)
colnames(cts_major_paren) <- c(rep("Epithelial", 8), rep("Fibroblast", 8), rep("Endothelial", 8), rep("Mural", 8), rep("Immune", 8), rep("Parenchyma", 5))

celltype <- factor(paste(colnames(cts_major_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_major_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_major_paren <- cts_major_paren[,colnames(cts_major_paren)%in%rownames(Maindm)]
cts_major_paren <- cts_major_paren[,colnames(cts_major_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, Epithelial=celltypeEpithelial-(celltypeFibroblast+celltypeEndothelial+celltypeMural+celltypeImmune+celltypeParenchyma)/5, Fibroblast=celltypeFibroblast-(celltypeEpithelial+celltypeEndothelial+celltypeMural+celltypeImmune+celltypeParenchyma)/5, Endothelial=celltypeEndothelial-(celltypeEpithelial+celltypeFibroblast+celltypeMural+celltypeImmune+celltypeParenchyma)/5, Mural=celltypeMural-(celltypeEpithelial+celltypeFibroblast+celltypeEndothelial+celltypeImmune+celltypeParenchyma)/5, Immune=celltypeImmune-(celltypeEpithelial+celltypeFibroblast+celltypeEndothelial+celltypeMural+celltypeParenchyma)/5, Parenchyma=celltypeParenchyma-(celltypeEpithelial+celltypeFibroblast+celltypeEndothelial+celltypeMural+celltypeImmune)/5)

cts_major_paren <- cts_major_paren[rowSums(cts_major_paren)>2*ncol(cts_major_paren),]
vdat <- voom(cts_major_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}



###Within Epi

cts_sub_paren <- cbind(cts_Epi_1, cts_Epi_2, cts_Epi_3, cts_Epi_4)
colnames(cts_sub_paren) <- c(rep("Epi_1", 8), rep("Epi_2", 8), rep("Epi_3", 8), rep("Epi_4", 8))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Epi_1=celltypeEpi_1-(celltypeEpi_2+celltypeEpi_3+celltypeEpi_4)/3, 
                        Epi_2=celltypeEpi_2-(celltypeEpi_1+celltypeEpi_3+celltypeEpi_4)/3, 
                        Epi_3=celltypeEpi_3-(celltypeEpi_4+celltypeEpi_2+celltypeEpi_4)/3, 
                        Epi_4=celltypeEpi_4-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3)/3)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Within Mural

cts_sub_paren <- cbind(cts_Mural_1, cts_Mural_2, cts_Mural_3)
colnames(cts_sub_paren) <- c(rep("Mural_1", 8), rep("Mural_2", 8),  rep("Mural_3", 8))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Mural_1=celltypeMural_1-(celltypeMural_2+celltypeMural_3)/2,
                        Mural_2=celltypeMural_2-(celltypeMural_1+celltypeMural_3)/2,
                        Mural_3=celltypeMural_3-(celltypeMural_1+celltypeMural_2)/2)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}




###In sn
####Load sn

sn_Epi_1_corrected <- read_csv("/Ex_RNA/Epi_1/Spreadsheet__Epi_1_corrected.csv")
sn_Epi_2_corrected <- read_csv("/Ex_RNA/Epi_2/Spreadsheet__Epi_2_corrected.csv")
sn_Epi_3_corrected <- read_csv("/Ex_RNA/Epi_3/Spreadsheet__Epi_3_corrected.csv")
sn_Epi_4_corrected <- read_csv("/Ex_RNA/Epi_4/Spreadsheet__Epi_4_corrected.csv")
sn_Mural_2_corrected <- read_csv("/Ex_RNA/Mural_2/Spreadsheet__Mural_2_corrected.csv")
sn_Mural_3_corrected <- read_csv("/Ex_RNA/Mural_3/Spreadsheet__Mural_3_corrected.csv")
sn_Mural_1_corrected <- read_csv("/Ex_RNA/Mural_1/Spreadsheet__Mural_1_corrected.csv")
sn_Epithelial_corrected <- read_csv("/Ex_RNA/Epithelial/Spreadsheet__Epithelial_corrected.csv")
sn_Fibroblast_corrected <- read_csv("/Ex_RNA/Fibroblast/Spreadsheet__Fibroblast_corrected.csv")
sn_Endothelial_corrected <- read_csv("/Ex_RNA/Endothelial/Spreadsheet__Endothelial_corrected.csv")
sn_Mural_corrected <- read_csv("/Ex_RNA/Mural/Spreadsheet__Mural_corrected.csv")
sn_Immune_corrected <- read_csv("/Ex_RNA/Immune/Spreadsheet__Immune_corrected.csv")
sn_Parenchyma_corrected <- read_csv("/Ex_RNA/Parenchyma/Spreadsheet__Parenchyma_corrected.csv")

sn_Epi_1 <- read_csv("/Ex_RNA/Epi_1/Spreadsheet__Epi_1_uncorrected.csv")
sn_Epi_2 <- read_csv("/Ex_RNA/Epi_2/Spreadsheet__Epi_2_uncorrected.csv")
sn_Epi_3 <- read_csv("/Ex_RNA/Epi_3/Spreadsheet__Epi_3_uncorrected.csv")
sn_Epi_4 <- read_csv("/Ex_RNA/Epi_4/Spreadsheet__Epi_4_uncorrected.csv")
sn_Mural_2 <- read_csv("/Ex_RNA/Mural_2/Spreadsheet__Mural_2_uncorrected.csv")
sn_Mural_3 <- read_csv("/Ex_RNA/Mural_3/Spreadsheet__Mural_3_uncorrected.csv")
sn_Mural_1 <- read_csv("/Ex_RNA/Mural_1/Spreadsheet__Mural_1_uncorrected.csv")
sn_Epithelial <- read_csv("/Ex_RNA/Epithelial/Spreadsheet__Epithelial_uncorrected.csv")
sn_Fibroblast <- read_csv("/Ex_RNA/Fibroblast/Spreadsheet__Fibroblast_uncorrected.csv")
sn_Endothelial <- read_csv("/Ex_RNA/Endothelial/Spreadsheet__Endothelial_uncorrected.csv")
sn_Mural <- read_csv("/Ex_RNA/Mural/Spreadsheet__Mural_uncorrected.csv")
sn_Immune <- read_csv("/Ex_RNA/Immune/Spreadsheet__Immune_uncorrected.csv")
sn_Parenchyma <- read_csv("/Ex_RNA/Parenchyma/Spreadsheet__Parenchyma_uncorrected.csv")


####Overlap venns

pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX2 = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))
sn_cosmx_overlap_fun(sn_Epi_1_corrected, celltypeDEGs_Epi_1_sig)
sn_cosmx_overlap_fun(sn_Epi_2_corrected, celltypeDEGs_Epi_2_sig)
sn_cosmx_overlap_fun(sn_Epi_3_corrected, celltypeDEGs_Epi_3_sig)
sn_cosmx_overlap_fun(sn_Epi_4_corrected, celltypeDEGs_Epi_4_sig)

sn_cosmx_overlap_fun(sn_Epithelial_corrected, celltypeDEGs_Epithelial_sig)
sn_cosmx_overlap_fun(sn_Fibroblast_corrected, celltypeDEGs_Fibroblast_sig)
sn_cosmx_overlap_fun(sn_Endothelial_corrected, celltypeDEGs_Endothelial_sig)
sn_cosmx_overlap_fun(sn_Immune_corrected, celltypeDEGs_Immune_sig)
sn_cosmx_overlap_fun(sn_Parenchyma_corrected, celltypeDEGs_Parenchyma_sig)

Featurefun(Cos_ChP, c("IBA1", "F13A1", "HLA", "TREM2", "LYVE1", "CD45", "CX3CR1", "P2RY12", "CD68", "TMEM119", "SPI1", "PF4", "MRC1", "MRC2"), "test")


##Figures Celltype
####UMAP

DimPlot(Cos_ChP, raster=T)

#Major types
Idents(Cos_ChP) <- Cos_ChP@meta.data[["Majorcelltype"]]
plts <- DimPlot(Cos_ChP, raster=F,raster.dpi=c(600,600), cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Parenchyma"="sienna", "Immune"="#26532B")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP_major.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP_major.png", sep=""), width=10, height=5)
remove(plts)

#Subcelltypes
Idents(Cos_ChP) <- Cos_ChP@meta.data[["Subcelltype"]] 
plts <- DimPlot(Cos_ChP, raster=F, raster.dpi=c(600,600), cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Endo_1"="#c875c4", "Endo_2"="#976697", "Endo_3"="#c5aac5", "Mural_1"="#261132", "Mural_2"="#6d207b", "Mural_3"="#8f2aa2", "Parenchyma"="brown", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP.png", sep=""), width=10, height=5)
remove(plts)

Cos_ChP <- subset(Cos_ChP, Majorcelltype!="Parenchyma")
Idents(Cos_ChP) <- Cos_ChP@meta.data[["Subcelltype"]] 
plts <- DimPlot(Cos_ChP, raster=T,raster.dpi=c(600,600), cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Endo_1"="#c875c4", "Endo_2"="#976697", "Endo_3"="#c5aac5", "Mural_1"="#261132", "Mural_2"="#6d207b", "Mural_3"="#8f2aa2", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP_noParenchyma.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP_noParenchyma.png", sep=""), width=10, height=5)
remove(plts)

#Major types
Idents(Cos_ChP) <- Cos_ChP@meta.data$Majorcelltype
plts <- DimPlot(Cos_ChP, raster=T,raster.dpi=c(600,600), cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP_major_noParenchyma.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP_major_noParenchyma.png", sep=""), width=10, height=5)
remove(plts)

#Spatial
plts <- ImageDimPlot(Cos_ChP, fov = "Slide.1", group.by = "Subcelltype", dark.background=F, size=0.1, cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Endo_1"="#c875c4", "Endo_2"="#976697", "Endo_3"="#c5aac5", "Mural_1"="#261132", "Mural_2"="#6d207b", "Mural_3"="#8f2aa2", "Parenchyma"="brown", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Spatial_S1.pdf", sep=""), width=5, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Spatial_S1.png", sep=""), width=5, height=5)
plts <- ImageDimPlot(Cos_ChP, fov = "Slide.2", group.by = 'Subcelltype', dark.background=F, size=0.1, cols=c("Epi_1"="#FFA500", "Epi_2"="#CD8500", "Epi_3"="#FF4500", "Epi_4"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Endo_1"="#c875c4", "Endo_2"="#976697", "Endo_3"="#c5aac5", "Mural_1"="#261132", "Mural_2"="#6d207b", "Mural_3"="#8f2aa2", "Parenchyma"="brown", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Spatial_S2.pdf", sep=""), width=5, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Spatial_S2.png", sep=""), width=5, height=5)


###Heatmap

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
gc()

#select list
Epi_genes <- c("AQP1", "PRLR", "ABCA4", "ABCA5", "GALNT11", "ENPP2", "TPD52L1", "TRPV4", "GPM6A", "SFRP1", "SGK1", "NWD1", "CDH12", "FBLN1", "GLIS3", "HTR2C", "OTX2-AS1", "TRMT9B", "VAT1L", "WDR49", "ZBBX", "CLIC6", "ATP2B2", "STEAP1B", "SLC5A3", "DNAH12", "DNAH9", "DNAH6", "DNAH7", "DNAH11", "RMST", "ROR2", "RFX3", "ARL13B", "SPAG17", "TMEM72-AS1", "TMEM117", "RYR3", "SLC24A4", "SNTB1", "SLC12A2", "SLC39A12", "EYA4", "F5", "AGBL4", "IGFBP4", "RIMBP2", "FYB2", "SLC4A10", "CFAP43", "CRPPA", "OCA2", "GNA14", "NEK11", "CFAP54", "MTUS2", "CFAP44", "DGKI", "GPR155", "DOCK8", "POU2AF1", "PACRG", "HYDIN", "PRDM16", "SHANK2", "CFAP100", "TTR", "ESRRG", "EPHA4", "LRRIQ1", "CHST11", "COL9A1", "CTNND2", "CUX2", "DCDC1", "LRRC9", "MRPS6", "SDK1", "MYO1E", "SLC15A4", "SOX6", "SPAG16")
Epi34_genes <- c("COX2", "MT-ND4", "GPX3", "COX1", "COX3", "MT-ATP6", "MT-ATP8", "MT-CYB", "MT-ND1", "MT-ND2", "SLC26A3", "MT-ND3", "MT-ND4L", "MT-ND5", "MTRNR2L8")
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
  sub <- subset(Cos_ChP, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projID", assay = "Nanostring", slot = "counts", return.seurat = FALSE)
  cts <- cts[["Nanostring"]]
  assign(c(paste("cts", celltype, sep = "_")), cts, envir = .GlobalEnv)
}

for(i in names(Cos_ChP@images)){Cos_ChP@images[[i]]<- NULL} #remove images before subbing
cts_fun("Epi_1")
gc()
cts_fun("Epi_2")
gc()
cts_fun("Epi_3")
gc()
cts_fun("Epi_4")
gc()
cts_fun("Fib_1")
cts_fun("Fib_2")
cts_fun("Fib_3")
cts_fun("Fib_4")
cts_fun("Endo")
cts_fun("Mural")
cts_fun("BAM_1")
cts_fun("BAM_2")
cts_fun("Olig")


cts_fun <- function(celltype){
  sub <- subset(Cos_ChP, Majorcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projID", assay = "Nanostring", slot = "counts", return.seurat = FALSE)
  cts <- cts[["Nanostring"]]
  assign(c(paste("cts", celltype, sep = "_")), cts, envir = .GlobalEnv)
}

cts_fun("Epithelial")
cts_fun("Fibroblast")
cts_fun("Immune")
cts_fun("Parenchyma")

#sub
cts_sub <- cbind(cts_Epi_1, cts_Epi_2, cts_Epi_3, cts_Epi_4, cts_Fib_1, cts_Fib_2, cts_Fib_3, cts_Fib_4, cts_Endo, cts_Mural, cts_BAM_1, cts_BAM_2)
my_genes <- c(Epi_genes, Epi34_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, Fib4_genes, Endo_genes, Mural_genes, BAM_genes, BAM1_genes, BAM2_genes) %>% unique()
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX2 = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))
for (i in 1:length(pairs$HG38)) {
  my_genes <- str_replace(my_genes, pairs$HG38[i], pairs$CosMX2[i])
}
my_genes <- my_genes[my_genes %in% rownames(Cos_ChP)]

cts_sub <- data.frame(scale(cts_sub))
cts_sub <- cts_sub[rownames(cts_sub) %in% my_genes,] #!!!!no parenchyma in my genes
cts_sub <- cts_sub[order(match(rownames(cts_sub), my_genes)),]
cts_sub <- t(scale(t(cts_sub)))
htg <- Heatmap(cts_sub, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_noparen.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_noparen.pdf", sep = ""), width=7, height=7)

#major no paren
cts_mjr <- cbind(cts_Epithelial, cts_Fibroblast, cts_Endo, cts_Mural, cts_Immune)
#curated from findmarkers
my_genes <- c(Epi_genes, Fib_genes, Endo_genes, Mural_genes, BAM_genes) %>% unique()
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX2 = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))
for (i in 1:length(pairs$HG38)) {
  my_genes <- str_replace(my_genes, pairs$HG38[i], pairs$CosMX2[i])
}
my_genes <- my_genes[my_genes %in% rownames(Cos_ChP)]
#from limmma
my_genes <- c(celltypeDEGs_Epithelial_sig$gene, celltypeDEGs_Fibroblast_sig$gene, celltypeDEGs_Endothelial_sig$gene, celltypeDEGs_Mural_sig$gene, celltypeDEGs_Immune_sig$gene) %>% unique()
#overlap
my_genes <- c(overlap_Epithelial$gene, overlap_Fibroblast$gene, overlap_Endothelial$gene, overlap_Mural$gene, overlap_Immune$gene) %>% unique()

cts_mjr <- data.frame(scale(cts_mjr))
cts_mjr <- cts_mjr[rownames(cts_mjr) %in% my_genes,] #!!!!no parenchyma in my genes
cts_mjr <- cts_mjr[order(match(rownames(cts_mjr), my_genes)),]
cts_mjr <- t(scale(t(cts_mjr)))
htg <- Heatmap(cts_mjr, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_major_noparen.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_major_noparen.pdf", sep = ""), width=7, height=7)

#major with paren
cts_mjr <- cbind(cts_Epithelial, cts_Fibroblast, cts_Endo, cts_Mural, cts_Immune, cts_Parenchyma)
#curated from findmarkers
my_genes <- c(Epi_genes, Fib_genes, Endo_genes, Mural_genes, BAM_genes) %>% unique()
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX2 = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))
for (i in 1:length(pairs$HG38)) {
  my_genes <- str_replace(my_genes, pairs$HG38[i], pairs$CosMX2[i])
}
my_genes <- my_genes[my_genes %in% rownames(Cos_ChP)]
#from limmma
my_genes <- c(celltypeDEGs_Epithelial_sig$gene, celltypeDEGs_Fibroblast_sig$gene, celltypeDEGs_Endothelial_sig$gene, celltypeDEGs_Mural_sig$gene, celltypeDEGs_Immune_sig$gene, celltypeDEGs_Parenchyma_sig$gene) %>% unique()
my_genes <- c(overlap_Epithelial$gene, overlap_Fibroblast$gene, overlap_Endothelial$gene, overlap_Mural$gene, overlap_Immune$gene, overlap_Parenchyma$gene) %>% unique()

cts_mjr <- data.frame(scale(cts_mjr))
cts_mjr <- cts_mjr[rownames(cts_mjr) %in% my_genes,] #!!!!no parenchyma in my genes
cts_mjr <- cts_mjr[order(match(rownames(cts_mjr), my_genes)),]
cts_mjr <- t(scale(t(cts_mjr)))
htg <- Heatmap(cts_mjr, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_major.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid_major.pdf", sep = ""), width=7, height=7)


#currated findmarkers
my_genes <- c(Epi_genes, Epi34_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, Endo_genes, Mural_genes, BAM_genes, BAM1_genes, BAM2_genes, paren_genes) %>% unique()
pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX2 = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))
for (i in 1:length(pairs$HG38)) {
  my_genes <- str_replace(my_genes, pairs$HG38[i], pairs$CosMX2[i])
}
my_genes <- my_genes[my_genes %in% rownames(Cos_ChP)]
#Limma
my_genes <- c(celltypeDEGs_Epi_1_sig$gene, celltypeDEGs_Epi_2_sig$gene, celltypeDEGs_Epi_3_sig$gene, celltypeDEGs_Epi_4_sig$gene, celltypeDEGs_Fib_1_sig$gene, celltypeDEGs_Fib_2_sig$gene, celltypeDEGs_Fib_3_sig$gene, celltypeDEGs_Fib_4_sig$gene, celltypeDEGs_Endo_sig$gene, celltypeDEGs_Mural_sig$gene, celltypeDEGs_BAM_1_sig$gene, celltypeDEGs_BAM_2_sig$gene, celltypeDEGs_BAM_3_sig$gene, celltypeDEGs_T_Cell_sig$gene, celltypeDEGs_Oligo_sig$gene, celltypeDEGs_Astro_sig$gene, celltypeDEGs_Neuro_sig$gene) %>% unique()
#overlap
my_genes <- c(overlap_Epithelial$gene, overlap_Epi_1$gene, overlap_Epi_2$gene, overlap_Epi_3$gene, overlap_Epi_4$gene, overlap_Fibroblast$gene, overlap_Fib_1$gene, overlap_Fib_2$gene, overlap_Fib_3$gene, overlap_Fib_4$gene, overlap_Endo$gene, overlap_Mural$gene, overlap_Immune$gene, overlap_BAM_1$gene, overlap_BAM_2$gene, overlap_Oligo$gene) %>% unique()

cts_sub_paren <- cbind(cts_Epi_1, cts_Epi_2, cts_Epi_3, cts_Epi_4, cts_Fib_1, cts_Fib_2, cts_Fib_3, cts_Fib_4, cts_Endo, cts_Mural, cts_BAM_1, cts_BAM_2, cts_Olig)
colnames(cts_sub_paren) <- c(rep("Epi_1_", 8), rep("Epi_2", 8), rep("Epi_3", 8), rep("Epi_4", 8), rep("Fib_1", 8), rep("Fib_2", 8), rep("Fib_3", 8), rep("Fib_4", 8), rep("Endo", 8), rep("Mural", 8), rep("BAM_1", 8), rep("BAM_2", 8), rep("Oligo", 5))
cts_sub_paren <- data.frame(scale(cts_sub_paren))
cts_sub_paren <- cts_sub_paren[rownames(cts_sub_paren) %in% my_genes,]
cts_sub_paren <- cts_sub_paren[order(match(rownames(cts_sub_paren), my_genes)),]
cts_sub_paren <- t(scale(t(cts_sub_paren)))
htg <- Heatmap(cts_sub_paren, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_projid.pdf", sep = ""), width=7, height=7)



###Bubble VLN Feature QC plots
####Major

highlights_genes <- c("TTR", "HTR2C", "GPX3", "DCN", "LEPR", "FN1", "VWF", "PECAM1", "INSR", "MYH11", "TAGLN", "ACTA2", "STAB1", "CD163", "MSR1")
Idents(Cos_ChP) <- Cos_ChP@meta.data$Majorcelltype
levels(Cos_ChP) <- c("Epithelial", "Fibroblast", "Endothelial", "Mural", "Immune", "Parenchyma")

plts <- DotPlot(Cos_ChP, features = highlights_genes, split.by="AUC_Subcelltype", cols=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "brown")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_bubble.pdf", sep=""), width=5+length(highlights_genes)*.15, height=1.5+length(unique(Cos_ChP@meta.data[["Majorcelltype"]]))*.4)

plts <- VlnPlot(Cos_ChP, features = highlights_genes, split.by='Majorcelltype', stack=T, raster=T, cols=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "brown"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_vln.pdf", sep=""), width=3+length(highlights_genes)*.4, height=1.5+length(unique(Cos_ChP@meta.data[["Majorcelltype"]]))*1)

plts <- VlnPlot(Cos_ChP, features = highlights_genes, split.by='Majorcelltype', stack=T, raster=T, cols=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "brown"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5)) + geom_boxplot(outliers=F)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_vln_box.pdf", sep=""), width=3+length(highlights_genes)*.4, height=1.5+length(unique(Cos_ChP@meta.data[["Majorcelltype"]]))*1)

Featurefun(Cos_ChP, highlights_genes, "highlights_genes",  assay="Nanostring", max.cutoff=c(1))

Idents(Cos_ChP) <- Cos_ChP@meta.data$Majorcelltype
SCNDRdimfun(Cos_ChP, reduction="harmony", res="Majorcelltype", split.by =  c("Slide", "projID"), cols=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "brown"))
Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
SCNDRdimfun(Cos_ChP, reduction="harmony", res="Subcelltype", split.by = c("Slide", "projID"), cols=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#002fff", "#c875c4", "#c5aac5", "#976697", "#8f2aa2", "#6d207b", "#261132", "#70C25B", "#0E9554", "#266417", "#3ded97", "brown"))


####Epi sub

Epi <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Epi.RDS"))
Idents(Epi) <- Epi@meta.data$Subcelltype
levels(Epi) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4")
Idents(Epi) -> Epi@meta.data$Subcelltype

Epi_genes <- c("TTR", "HTR2C", "CDH1", "KRT8", "KRT18", "FOXJ1", "DNAH6", "TUBB4B", "CFAP100", "CFAP69", "CLDN5", "TJP1", "GPX3", "COX1", "COX2", "AQP1", "ATP1A1", "SLC2A12", "SLC25A6", "SLC25A23", "FOLR1", "IGFBP2", "VIM")

plts <- DotPlot(Epi, features = Epi_genes, split.by="Subcelltype", cols=c("#FFA500", "#CD8500", "#FF4500", "#CD3700")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble.pdf", sep=""), width=3+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Epi, features = Epi_genes, cols=c("grey", "#C48F45")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble_scale.pdf", sep=""), width=5+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Epi, features = Epi_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble_default.pdf", sep=""), width=5+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Epi, features = Epi_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#FFA500", "#CD8500", "#FF4500", "#CD3700"), same.y.lims=F) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_vln.pdf", sep=""), width=3+length(Epi_genes)*.4, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*1)
remove(plts)

plts <- VlnPlot(Epi, features = Epi_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#FFA500", "#CD8500", "#FF4500", "#CD3700"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))  + geom_boxplot(outliers=F)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_vln_box.pdf", sep=""), width=3+length(Epi_genes)*.4, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Epi, Epi_genes, "highlights_genes",  assay="soupXcounts", max.cutoff=c(3,4))

Idents(Epi) <- Epi@meta.data$Subcelltype
SCNDRdimfun(Epi, reduction="harmony", res="Subcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#FFA500", "#CD8500", "#FF4500", "#CD3700"))


####Fibro sub

#Fibro <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Fibro.RDS"))

Idents(Fibro) <- Fibro$Subcelltype
Fibro_genes <- c("DCN", "FN1", "LEPR", "CDH11", "STXBP6", "FBLN1", "KCNMA1", "SLC4A4", "PLXNA4", "ALPL", "DPEP1", "DLK1", "NRG1", "ROBO1", "TRPM3", "MYRIP", "GPM6A", "CLDN11", "ECM2", "SLC6A20", "SLC13A3", "ABCA9", "LAMA2", "LAMA4", "COL18A1", "COL4A2", "COL1A2", "COL4A4", "COL5A2",  "COL3A1", "COL6A3", "PDGFRA", "SLC4A7", "SLC2A3", "LUM")
plts <- DotPlot(Fibro, features = Fibro_genes, split.by='Subcelltype', cols=c("#050382", "#27D7EF", "#2a91ea")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_bubble.pdf", sep=""), width=3+length(Fibro_genes)*.15, height=1.5+length(unique(Fibro@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Fibro, features = Fibro_genes, cols=c("grey", "#3249A6")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_bubble_scale.pdf", sep=""), width=5+length(Fibro_genes)*.15, height=1.5+length(unique(Fibro@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Fibro, features = Fibro_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_bubble_default.pdf", sep=""), width=5+length(Fibro_genes)*.15, height=1.5+length(unique(Fibro@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Fibro, features = Fibro_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#050382", "#27D7EF", "#2a91ea"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_vln.pdf", sep=""), width=3+length(Fibro_genes)*.4, height=1.5+length(unique(Fibro@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Fibro, Fibro_genes, "highlights_genes",  assay="SCT", max.cutoff=c(3,4))

Idents(Fibro) <- Fibro@meta.data$Subcelltype
SCNDRdimfun(Fibro, reduction="harmony", res="Subcelltype", split.by = c("projid", "Slide"), cols=c("#050382", "#27D7EF", "#2a91ea"))


####Imm sub

#Imm <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Imm.RDS"))

Imm_genes <- c("F13A1", "LYVE1", "MRC1", "STAB1", "CD163", "CD68", "MERTK", "HCK", "IL4R", "CD74", "HLA-DRA", "S100A9", "CSF1R", "JAK2", "CD86", "VCAM1", "CD83", "LAMP1", "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1", "MKI67", "CENPF", "TOP2A", "ITGAM", "S100A8", "CCR2", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11", "THEMIS", "CD3E", "PDE7A", "CAMK4", "HLA-B", "CD8A") %>% unique()

Idents(Imm) <- Imm@meta.data$Subcelltype
plts <- DotPlot(Imm, features = Imm_genes, split.by="Subcelltype", cols=c("#70C25B", "#0E9554", "#266417", "#3ded97")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Imm, features = Imm_genes, cols=c("grey", "#26532B")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble_scale.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Imm, features = Imm_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble_default.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Imm, features = Imm_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#70C25B", "#0E9554", "#266417", "#3ded97"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_vln.pdf", sep=""), width=3+length(Imm_genes)*.4, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Imm, Imm_genes, "highlights_genes",  assay="soupXcounts", max.cutoff=c(3,4))

Idents(Imm) <- Imm@meta.data$Subcelltype
SCNDRdimfun(Imm, reduction="harmony", res="Subcelltype", cols=c("#70C25B", "#0E9554", "#266417", "#3ded97"))


####Endo sub

Endo <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Endo.RDS"))

Endo_genes <- c("PECAM1", "VWF", "INSR", "CLDN5", "CDH13", "MECOM", "VEGFC", "COL8A1", "FBN1",
                "TSHZ2", "CEMIP2", "IL1R1", "IL4R",
                "PLVAP", "FLT1", "ADGRF5", "KCNQ3", "SLC2A3", "ITGA1") %>% unique()

Idents(Endo) <- Endo@meta.data$Subcelltype
plts <- DotPlot(Endo, features = Endo_genes, split.by="Subcelltype", cols=c("#c875c4", "#c5aac5", "#976697")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_bubble.pdf", sep=""), width=5+length(Endo_genes)*.15, height=1.5+length(unique(Endo@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Endo, features = Endo_genes, cols=c("grey", "#926392")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_bubble_scale.pdf", sep=""), width=5+length(Endo_genes)*.15, height=1.5+length(unique(Endo@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Endo, features = Endo_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_bubble_default.pdf", sep=""), width=5+length(Endo_genes)*.15, height=1.5+length(unique(Endo@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Endo, features = Endo_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#c875c4", "#c5aac5", "#976697"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_vln.pdf", sep=""), width=3+length(Endo_genes)*.4, height=1.5+length(unique(Endo@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Endo, Endo_genes, "highlights_genes",  assay="SCT", max.cutoff=c(3,4))

Idents(Endo) <- Endo@meta.data$Subcelltype
SCNDRdimfun(Endo, reduction="harmony", res="Subcelltype", split.by = c("projID", "Slide"), cols=c("#c875c4", "#c5aac5", "#976697"))
plts <- DimPlot(Endo, raster=F, cols=c("#c875c4", "#c5aac5", "#976697")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_UMAP.pdf", sep=""), width=10, height=5)


####Mural sub

Mural <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Mural.RDS"))

Mural_genes <- c("ENPP2","PLXDC2","PALMD", "MID1", "CEMIP", "GLIS3", "PDZD2", 
                 "COL4A4", "CRISPLD2", "ADAMTS9", "TLE1", "PIPNC1", "ADCY3", "HIPK2", "COL18A1", "SCL22A3",
                 "SOX6", "MYO1D", "MYH11", "ACTA2", "PRUNE2", "RYR2", "TPM1", "LMOD1") %>% unique()

Idents(Mural) <- Mural@meta.data$Subcelltype
plts <- DotPlot(Mural, features = Mural_genes, split.by="Subcelltype", cols=c("#8f2aa2", "#6d207b", "#261132")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_bubble.pdf", sep=""), width=5+length(Mural_genes)*.15, height=1.5+length(unique(Mural@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Mural, features = Mural_genes, cols=c("grey", "#261132")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_bubble_scale.pdf", sep=""), width=5+length(Mural_genes)*.15, height=1.5+length(unique(Mural@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Mural, features = Mural_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_bubble_default.pdf", sep=""), width=5+length(Mural_genes)*.15, height=1.5+length(unique(Mural@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Mural, features = Mural_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#8f2aa2", "#6d207b", "#261132"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_vln.pdf", sep=""), width=3+length(Mural_genes)*.4, height=1.5+length(unique(Mural@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Mural, Mural_genes, "highlights_genes",  assay="SCT", max.cutoff=c(3,4))

Idents(Mural) <- Mural@meta.data$Subcelltype
SCNDRdimfun(Mural, reduction="harmony", res="Subcelltype", split.by = c("projID", "Slide"), cols=c("#8f2aa2", "#6d207b", "#261132"))
plts <- DimPlot(Mural, raster=F, cols=c("#8f2aa2", "#6d207b", "#261132")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_UMAP.pdf", sep=""), width=10, height=5)

#Major types
Idents(Cos_ChP) <- Cos_ChP@meta.data$Majorcelltype
plts <- DimPlot(Cos_ChP, raster=T,raster.dpi=c(600,600), cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP_major_noParenchyma.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Cos_ChP_UMAP_major_noParenchyma.png", sep=""), width=10, height=5)
remove(plts)




####QC

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
cols <- c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#6d207b", "#8f2aa2", "#70C25B", "#0E9554", "#266417", "#3ded97", "brown")
plts <- VlnPlot(Cos_ChP, "nFeature_Nanostring", pt.size=F, cols=cols)
plot(plts)
ggsave2(paste(curr_dir, "/CellType/Publish/QC_Vln_nFeature_Nanostring", ".pdf", sep=""), plot = plts, width=4, height=4)
plts <- VlnPlot(Cos_ChP, "nCount_Nanostring", pt.size=F, cols=cols)
plot(plts)
ggsave2(paste(curr_dir, "/CellType/Publish/QC_Vln_nCount_Nanostring", ".pdf", sep=""), plot = plts, width=4, height=4)
plts <- VlnPlot(Cos_ChP, "Area.um2", pt.size=F, cols=cols)
plot(plts)
ggsave2(paste(curr_dir, "/CellType/Publish/QC_Vln_Area.um2", ".pdf", sep=""), plot = plts, width=4, height=4)
plts <- VlnPlot(Cos_ChP, "complexity", pt.size=F, cols=cols)
plot(plts)
ggsave2(paste(curr_dir, "/CellType/Publish/QC_Vln_complexity", ".pdf", sep=""), plot = plts, width=4, height=4)


###Correlation between sn RNAseq and CosMX2
####major

my_genes <- c(overlap_Epithelial$gene, overlap_Fibroblast$gene, overlap_Endothelial$gene, overlap_Mural$gene, overlap_Immune$gene, overlap_Parenchyma$gene) %>% unique()

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Epithelial_SN", "Fibroblast_SN", "Endothelial_SN", "Mural_SN", "Immune_SN", "Parenchyma_SN")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Cos_cts <- data.frame(AggregateExpression(Cos_ChP, group.by = "Majorcelltype", assay = "SCT", slot = "data", return.seurat = FALSE))
colnames(Cos_cts) <- c("Epithelial_ST", "Fibroblast_ST", "Endothelial_ST", "Mural_ST", "Immune_ST", "Parenchyma_ST")
Cos_cts <- data.frame(t(scale(t(Cos_cts))))
Cos_cts <- Cos_cts %>% rownames_to_column("gene")
Cos_cts <- Cos_cts[Cos_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, Cos_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epithelial_ST", "Epithelial_SN", "Fibroblast_ST", "Fibroblast_SN", "Endothelial_ST", "Endothelial_SN", "Mural_ST", "Mural_SN", "Immune_ST", "Immune_SN", "Parenchyma_ST", "Parenchyma_SN")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("ST", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
Heatmap(cts_cor)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Major.pdf", sep = ""), width=7, height=7)


####within Epi

my_genes <- c(overlap_Epithelial$gene, overlap_Epi_1$gene, overlap_Epi_2$gene, overlap_Epi_3$gene, overlap_Epi_4$gene) %>% unique()

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Epi_sn, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Epi_1_SN", "Epi_2_SN", "Epi_3_SN", "Epi_4_SN")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Cos_cts <- data.frame(AggregateExpression(Epi, group.by = "Subcelltype", assay = "SCT", slot = "data", return.seurat = FALSE))
colnames(Cos_cts) <- c("Epi_1_Cos", "Epi_2_Cos", "Epi_3_Cos", "Epi_4_Cos")
Cos_cts <- data.frame(t(scale(t(Cos_cts))))
Cos_cts <- Cos_cts %>% rownames_to_column("gene")
Cos_cts <- Cos_cts[Cos_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, Cos_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_1_Cos", "Epi_1_SN", "Epi_2_Cos", "Epi_2_SN", "Epi_3_Cos", "Epi_3_SN", "Epi_4_Cos", "Epi_4_SN")]

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Cos", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
Heatmap(cts_cor)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Epi.pdf", sep = ""), width=7, height=7)



####within Mural

my_genes <- c(sn_Mural_corrected$gene, sn_Mural_2_corrected$gene, sn_Mural_3_corrected$gene, sn_Mural_1_corrected$gene) %>% unique()

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Mural_sn, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Mural_2_SN", "Mural_3_SN", "Mural_1_3_SN")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Cos_cts <- data.frame(AggregateExpression(Epi, group.by = "Subcelltype", assay = "SCT", slot = "data", return.seurat = FALSE))
colnames(Cos_cts) <- c("Mural_2_Cos", "Mural_3_Cos", "Mural_1_Cos")
Cos_cts <- data.frame(t(scale(t(Cos_cts))))
Cos_cts <- Cos_cts %>% rownames_to_column("gene")
Cos_cts <- Cos_cts[Cos_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, Cos_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Mural_2_ST", "Mural_2_SN", "Mural_3_ST", "Mural_3_SN", "Mural_1_ST", "Mural_1_SN")]

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Cos", rownames(cts_cor)), grepl("SN", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F)
Heatmap(cts_cor)
plts
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Mural.pdf", sep = ""), width=7, height=7)


###FGSEA celltype

Path_fun(celltypeDEGs_Epi_1)
Path_fun(celltypeDEGs_Epi_2)
Path_fun(celltypeDEGs_Epi_3)
Path_fun(celltypeDEGs_Epi_4)

Path_fun(celltypeDEGs_Epithelial)
Path_fun(celltypeDEGs_Fibroblast)
Path_fun(celltypeDEGs_Mural)
Path_fun(celltypeDEGs_Endothelial)
Path_fun(celltypeDEGs_Immune)
Path_fun(celltypeDEGs_Parenchyma)



###Ncells

Idents(Cos_ChP) <- Cos_ChP@meta.data$Subcelltype
#Cell numbers
Ncellfun(Cos_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Cos_ChP_projID.csv", sep=""), group.by="projID")
Ncellfun(Cos_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Cos_ChP_cogdx.csv", sep=""), group.by="cogdx")
Ncellfun(Cos_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Cos_ChP_Cog_Path.csv", sep=""), group.by="Cog_Path")
Ncellfun(Cos_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Cos_ChP_apoe.csv", sep=""), group.by="apoe_genotype")

Idents(Cos_ChP) <- Cos_ChP@meta.data$Majorcelltype
#Cell numbers
Ncellfun(Cos_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Cos_ChP_Majorcelltype_projID.csv", sep=""), group.by="projID")
Ncellfun(Cos_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Cos_ChP_Majorcelltype_cogdx.csv", sep=""), group.by="cogdx")
Ncellfun(Cos_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Cos_ChP_Majorcelltype_Cog_Path.csv", sep=""), group.by="Cog_Path")
Ncellfun(Cos_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Cos_ChP_Majorcelltype_apoe.csv", sep=""), group.by="apoe_genotype")
