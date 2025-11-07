###
### 035 cogdx_AD_differential_ATAC.R
###
# Purpose: AD differential in single nuclei ATAC dataset
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_ATAC" #output results to specific folder
set.seed(42)

pcutoff <- 0.05
LFcutoff <- 0.25
pct_cutoff <- 1 #Use standardized variable peaks cutoff not pct in cells (it's much lower than RNA). Anything >1 is more variable than chance.


##Subset
Ex_ChP_ATAC <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP_ATAC.RDS"))

#!!!!Change ATAC/ACTIVITY and Major/Sub
datprepfun_reg(Ex_ChP_ATAC, "Epithelial")
datprepfun_reg(Ex_ChP_ATAC, "Fibroblast")
datprepfun_reg(Ex_ChP_ATAC, "Mural")
datprepfun_reg(Ex_ChP_ATAC, "Immune")

datprepfun_reg(Ex_ChP_ATAC, "BAM_1")
datprepfun_reg(Ex_ChP_ATAC, "BAM_2")
datprepfun_reg(Ex_ChP_ATAC, "Epi_1")
datprepfun_reg(Ex_ChP_ATAC, "Epi_2")
datprepfun_reg(Ex_ChP_ATAC, "Epi_3")
datprepfun_reg(Ex_ChP_ATAC, "Epi_4")


##limma runs
###factors set up
#discrete var
msex = factor(paste(mdat$msex))
X10X_batch = factor(paste(mdat$X10X_batch))
library_batch = factor(paste(mdat$library_batch))
apoe = factor(paste(mdat$apoe_genotype))
apoe_sex = factor(paste(mdat$apoe_genotype, mdat$msex, sep="."))
cog_path = factor(paste(mdat$Cog_Path))
ad_adnc <- factor(paste(mdat$ad_adnc))

#continuous
cogdx_c = as.numeric(mdat$cogdx)
projid = as.numeric(mdat$projid)
age = as.numeric(mdat$age_death)
pmi = as.numeric(mdat$pmi)
educ = as.numeric(mdat$educ)
gpath = as.numeric(mdat$gpath)
Braaksc = as.numeric(mdat$braaksc)
tdp = as.numeric(mdat$tdp_st4)
tangles = as.numeric(mdat$tangsqrt_est_8reg)
amyl = as.numeric(mdat$amylsqrt_est_8reg)

arteriol = as.numeric(mdat$arteriol_scler)
caa = as.numeric(mdat$caa_4gp)
cvda = as.numeric(mdat$cvda_4gp2)
cogn_global = as.numeric(mdat$cogn_global_last)


###cogdx
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cog_path + age)

rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdx_MCIvNCI=cog_pathMCI-cog_pathNCI, cogdx_ADvNCI=cog_pathAD-cog_pathNCI, cogdx_ADvMCI=cog_pathAD-cog_pathMCI)

limmafun(cts_Epithelial)
limmafun(cts_Epi_1)
limmafun(cts_Epi_2)
limmafun(cts_Epi_3)
limmafun(cts_Epi_4)
limmafun(cts_Fibroblast)
limmafun(cts_Immune)

limmafun(cts_ATAC_Epithelial)
limmafun(cts_ATAC_Epi_1)
limmafun(cts_ATAC_Epi_2)
limmafun(cts_ATAC_Epi_3)
#limmafun(cts_ATAC_Fibroblast) #Sanity checks don't look right
limmafun(cts_ATAC_Immune)


###AD_adnc (Epis, Fibs, and hHSP BAMs) 
options(na.action='na.pass')
Maindm <- model.matrix(~0 + ad_adnc + age)

rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, adnc_ADvNCI=ad_adnc1-ad_adnc0)

limmafun(cts_Epithelial)
limmafun(cts_Epi_1)
limmafun(cts_Epi_2)
limmafun(cts_Epi_3)
limmafun(cts_Epi_4)
#limmafun(cts_ATAC_Fibroblast) #Sanity checks don't look right
limmafun(cts_Immune)

limmafun(cts_ATAC_Epithelial)
limmafun(cts_ATAC_Epi_1)
limmafun(cts_ATAC_Epi_2)
limmafun(cts_ATAC_Epi_3)
limmafun(cts_ATAC_Fibroblast)
limmafun(cts_ATAC_Immune)


###continuous gpath
options(na.action='na.pass')
Maindm <- model.matrix(~0 + gpath + age)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, gpath)

limmafun(cts_Epithelial)


###continuous Braaksc (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + Braaksc + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, Braaksc)

limmafun(cts_Epithelial)


###continuous tdp (Epis, Fibs, and hHSP BAMs) 
options(na.action='na.pass')
Maindm <- model.matrix(~0 + tdp + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, tdp)

limmafun(cts_Epithelial)


###continuous tangles (Epis, Fibs, and hHSP BAMs) 
options(na.action='na.pass')
Maindm <- model.matrix(~0 + tangles + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, tangles)

limmafun(cts_Epithelial)


###continuous amyloid burden (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + amyl + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, amyl)

limmafun(cts_Epithelial)


###continuous arteriol (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + arteriol + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, arteriol)

limmafun(cts_Epithelial)


###continuous caa (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + caa + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, caa)

limmafun(cts_Epithelial)


###continuous cvda (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cvda + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cvda)

limmafun(cts_Epithelial)

###continuous last cog score

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogn_global + age)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogn_global)

limmafun(cts_Epithelial)
limmafun(cts_Epi_1)
limmafun(cts_Epi_2)
limmafun(cts_Epi_3)
limmafun(cts_Epi_4)
#limmafun(cts_ATAC_Fibroblast) #Sanity checks don't look right
limmafun(cts_Immune)

limmafun(cts_ATAC_Epithelial)
limmafun(cts_ATAC_Epi_1)
limmafun(cts_ATAC_Epi_2)
limmafun(cts_ATAC_Epi_3)
limmafun(cts_ATAC_Fibroblast)
limmafun(cts_ATAC_Immune)




#FGSEA
Path_fun(limma_cogn_global_Epithelial, invert=T)
Path_fun(limma_cogn_global_Fibroblast, invert=T)

Path_fun(limma_cogdx_ADvNCI_Epithelial)
Path_fun(limma_cogdx_ADvNCI_Epi_1)
Path_fun(limma_cogdx_ADvNCI_Epi_2)
Path_fun(limma_cogdx_ADvNCI_Epi_3)
Path_fun(limma_cogdx_ADvNCI_Epi_4)
Path_fun(limma_cogdx_ADvNCI_Fibroblast)
Path_fun(limma_cogdx_ADvNCI_Immune)

Path_fun(limma_cogdx_MCIvNCI_Epithelial)
Path_fun(limma_cogdx_MCIvNCI_Epi_1)
Path_fun(limma_cogdx_MCIvNCI_Epi_2)
Path_fun(limma_cogdx_MCIvNCI_Epi_3)
Path_fun(limma_cogdx_MCIvNCI_Epi_4)
Path_fun(limma_cogdx_MCIvNCI_Fibroblast)
Path_fun(limma_cogdx_MCIvNCI_Immune)

Path_fun(limma_cogdx_ADvMCI_Epithelial)
Path_fun(limma_cogdx_ADvMCI_Epi_1)
Path_fun(limma_cogdx_ADvMCI_Epi_2)
Path_fun(limma_cogdx_ADvMCI_Epi_3)
Path_fun(limma_cogdx_ADvMCI_Epi_4)
Path_fun(limma_cogdx_ADvMCI_Fibroblast)
Path_fun(limma_cogdx_ADvMCI_Immune)

Path_fun(Pattern_Epithelial_TM_late)

Path_fun(limma_adnc_ADvNCI_Epithelial)
Path_fun(limma_adnc_ADvNCI_Fibroblast)
Path_fun(limma_adnc_ADvNCI_Immune)


Path_fun(limma_gpath_Epithelial)

Path_fun(limma_tdp_Epithelial)

Path_fun(limma_tangles_Epithelial)

Path_fun(limma_amyl_Epithelial)



#Pretty figs for pub
##GO heatmap of overlap
col_fun <- colorRamp2(c(-2, 0, 2), c("blue","white", "red"))


###Epi
selection <- c("GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOBP_CILIUM_MOVEMENT", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT", "GOBP_EXPORT_FROM_CELL", "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOBP_INFLAMMATORY_RESPONSE", "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "GOCC_ANCHORING_JUNCTION", "GOCC_CELL_CELL_JUNCTION", "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS", "GOBP_REGULATION_OF_RESPONSE_TO_STRESS")

##per contrast
GSEA_Epithelial_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Epithelial, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Epithelial, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogdx_ADvMCI_Epithelial, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogn_global_Epithelial, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Epithelial, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Epithelial, by="pathway", suffix=c("", "_gpath"))

GSEA_Epithelial_join <- GSEA_Epithelial_join %>% column_to_rownames("pathway")
GSEA_Epithelial_join <- GSEA_Epithelial_join[selection,]
GSEA_Epithelial_join <- GSEA_Epithelial_join[match(selection, rownames(GSEA_Epithelial_join)),]
mat <- GSEA_Epithelial_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Epithelial_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
colnames(mat2) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat2[is.na(mat2)] <- 1
mat2 <- t(mat2)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(7, "mm"), height = nrow(mat)*unit(7, "mm"), col = colorRamp2(c(-2, 0, 2), c("blue","white", "red")), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y)
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Epithelial.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_MCIvNCI), col=col_fun(mat$cogdx_MCIvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvMCI), col=col_fun(mat$cogdx_ADvMCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()



###CoveragePlot of select markers
####Epithelial
DefaultAssay(Ex_ChP_ATAC) <- "ATAC"
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC$Cog_Path
levels(Ex_ChP_ATAC) <- c("NCI", "MCI", "AD")
SN_Epithelial_cogdx_ADvNCI_corrected <- read_csv("limma/cogdx_ADvNCI/Spreadsheet_Epithelial_cogdx_ADvNCI_corrected.csv")

Epi_genes <- c("CACNA2D2", "SLC9A3R1", "TUBA4A", "DYNC1H1", "CROCC", "KLF4", "TCF12", "FOXJ1")

test <- limma_cogdx_ADvNCI_ATAC_Epithelial$gene
closest_totest <- ClosestFeature(Ex_ChP_ATAC, regions=test)
closest_totest <- left_join(closest_totest, Gene.search, join_by("gene_name"=="gene"))
view(closest_totest)
closest_totest[closest_totest$gene_name %in% SN_Epithelial_cogdx_ADvNCI_corrected$gene,]
closest_totest[closest_totest$gene_name %in% Epi_genes,]

for(gene in Epi_genes){
  highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
  plts <- CoveragePlot(Ex_ChP_ATAC, region=gene, region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("skyblue", "#F86247", "red3"))
  plot(plts)
  graph2pdf(plts, file=paste0(curr_dir, "/limma/CoveragePlots/", gene, ".pdf"))
}

gene <- "DYNC1H1"
highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
plts <- CoveragePlot(Ex_ChP_ATAC, region="chr14-101960000-102000000", region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("skyblue", "#F86247", "red3"))
plts
graph2pdf(plts, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/", gene, "_sub.pdf"))


####Fibroblast
DefaultAssay(Ex_ChP_ATAC) <- "ATAC"
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC$Cog_Path
SN_Fibroblast_cogdx_ADvNCI_corrected <- read_csv("limma/cogdx_ADvNCI/Spreadsheet_Fibroblast_cogdx_ADvNCI_corrected.csv")

Fib_genes <- c("COL8A1", "COL15A1", "ADAMTS9", "CD96", "APOE", "ACSL4", "APOLD1", "KCNMA1", "SLC38A2", "TFRC", "CD9", "GPX4", "C7", "GPX3")

test <- limma_cogdx_ADvNCI_ATAC_Fibroblast$gene
closest_totest <- ClosestFeature(Ex_ChP_ATAC, regions=test)
closest_totest <- left_join(closest_totest, Gene.search, join_by("gene_name"=="gene"))
view(closest_totest)
closest_totest[closest_totest$gene_name %in% SN_Fibroblast_cogdx_ADvNCI_corrected$gene,]
closest_totest[closest_totest$gene_name %in% Fib_genes,]

for(gene in Fib_genes){
  highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
  plts <- CoveragePlot(Ex_ChP_ATAC, region=gene, region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("skyblue", "#F86247", "red3"))
  plot(plts)
  graph2pdf(plts, file=paste0(curr_dir, "/limma/CoveragePlots/", gene, ".pdf"))
}


####Immune
DefaultAssay(Ex_ChP_ATAC) <- "ATAC"
Idents(Ex_ChP_ATAC) <- Ex_ChP_ATAC$Cog_Path
SN_Immune_cogdx_ADvNCI_corrected <- read_csv("limma/cogdx_ADvNCI/Spreadsheet_Immune_cogdx_ADvNCI_corrected.csv")

Imm_genes <- c("HSP90B1", "HSPA6", "HSPA1B", "HSP90AA1", "HLA-DRA", "HLA-DMA", "HLA-DPA1", "HLA-E")

test <- limma_cogdx_ADvNCI_ATAC_Immune_sig$gene
closest_totest <- ClosestFeature(Ex_ChP_ATAC, regions=test)
closest_totest <- left_join(closest_totest, Gene.search, join_by("gene_name"=="gene"))
view(closest_totest)
closest_totest[closest_totest$gene_name %in% SN_Immune_cogdx_ADvNCI_corrected$gene,]
closest_totest[closest_totest$gene_name %in% Imm_genes,]

for(gene in Imm_genes){
  highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
  plts <- CoveragePlot(Ex_ChP_ATAC, region=gene, region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("skyblue", "#F86247", "red3"))
  plot(plts)
  graph2pdf(plts, file=paste0(curr_dir, "/limma/CoveragePlots/", gene, ".pdf"))
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
Epi_enriched.motifs <- FindMotifs(Epi, features=unique(c(limma_cogdx_ADvMCI_ATAC_Epithelial_sig$gene, limma_cogdx_ADvNCI_ATAC_Epithelial_sig$gene, limma_cogdx_MCIvNCI_ATAC_Epithelial_sig$gene)), only.pos=F, verbose=T)

Imm_enriched.motifs <- FindMotifs(Imm, features=unique(c(limma_cogdx_ADvMCI_ATAC_Immune_sig$gene, limma_cogdx_ADvNCI_ATAC_Immune_sig$gene, limma_cogdx_MCIvNCI_ATAC_Immune_sig$gene)), only.pos=F, verbose=T)


#TF footprinting
Epi <- Footprint(Epi, motif.name ="MA1513.1", genome = BSgenome.Hsapiens.UCSC.hg38, compute.expected = FALSE)
PlotFootprint(Epi, features = c("MA1513.1"))
#Not different
for(i in Epi_enriched.motifs$motif.name[1:5]){
  print(i)
  Epi <- Footprint(Epi, motif.name =i, genome = BSgenome.Hsapiens.UCSC.hg38, compute.expected = FALSE)
  plot(PlotFootprint(Epi, features = c(i)))
}


####Epithelial
Epi <- subset(Ex_ChP_ATAC, Majorcelltype=="Epithelial")
match_names <- data.frame(matches@colData@listData[["name"]])
colnames(match_names) <- "Gene"
match_names <- rownames_to_column(match_names, "group_name")

for(gene in c("TUBA4A", "DYNC1H1", "CROCC", "CFAP100", "CFAP69", "FOXJ1", "KLF4", "TCF12")){
  test <- limma_cogdx_ADvNCI_ATAC_Epithelial$gene
  highlights <- subsetByOverlaps(StringToGRanges(test), LookupGeneCoords(Ex_ChP_ATAC, gene))
  
  matches <- matchMotifs(pwm, subject=highlights, genome=BSgenome.Hsapiens.UCSC.hg38, out="scores")
  head(motifMatches(matches))
  head(motifScores(matches))
  
  Matched <-data.frame(matches@assays@data@listData[["motifScores"]])
  colnames(Matched) <- matches@colData$name
  Matched_sig <- Matched[, colSums(Matched)>0]
  view(Matched_sig)
  
  matches_pos <- matchMotifs(pwm, subject=highlights, genome=BSgenome.Hsapiens.UCSC.hg38, out="positions")
  motif.df <- as.data.frame(matches_pos)
  motif.df <- merge(motif.df, match_names, by="group_name")
  
  region <- FindRegion(object = Ex_ChP_ATAC, region = gene, sep = c("-", "-"), 
                       assay = "ATAC", extend.upstream = 1000, 
                       extend.downstream = 1000)
  start.pos <- start(x=region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  plts <- ggplot(data = motif.df) +
    geom_segment(aes(x = start, y = 0, xend = end, yend = 0), size = 2, data = motif.df) + theme_classic() + ylab(label = "Motif") + 
    geom_label(aes(label=str_wrap(Gene,12), x=(start + end)/2, y=0), size=3, position="jitter") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))
  Cov_plts <- CoveragePlot(Ex_ChP_ATAC, region=gene, region.highlight=highlights, extend.downstream=1000, extend.upstream=1000, features=gene, annotation=T, peaks=T, links=T) & scale_fill_manual(values=c("skyblue", "#F86247", "red3"))
  heights <- c(26,40)
  widths <- c(10, length(x = features))
  comb_plts <- CombineTracks(plotlist = list(Cov_plts, plts), heights = heights, widths = widths)
  plot(comb_plts)
  graph2pdf(comb_plts, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/Motifs_", gene, ".pdf"))
  write_csv(motif.df, file=paste0(curr_dir, "/CellType/Publish/CoveragePlots/Motifs_", gene, ".csv"))
}

test <- merge(Motifs_FOXJ1, merge(Motifs_CROCC, merge(Motifs_DYNC1H1, Motifs_TUBA4A, by="Gene"), by="Gene"), by="Gene") #Check for TFs in common

for(i in colnames(Matched_sig)[1:5]){
  print(i)
  Epi <- Footprint(Epi, motif.name =i, genome = BSgenome.Hsapiens.UCSC.hg38, compute.expected = FALSE)
  plot(PlotFootprint(Epi, features = c(i)))
}

for(i in c("KLF4", "TCF12(var.2)", "FOXP2")){
  print(i)
  Epi <- Footprint(Epi, motif.name =i, genome = BSgenome.Hsapiens.UCSC.hg38, compute.expected = FALSE)
  plts <- PlotFootprint(Epi, features = c(i))
  plot(plts)
  graph2pdf(comb_plts, file=paste0(curr_dir, "/limma/CoveragePlots/Footprint_", i, ".pdf"))
}

my_region <- GRanges(seqnames = "chr2", ranges= IRanges(start=219266735, end=219266745))
getSeq(BSgenome.Hsapiens.UCSC.hg38, my_region)
