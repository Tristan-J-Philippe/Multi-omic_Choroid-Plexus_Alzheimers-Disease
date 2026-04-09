###
### 038 cogdx_AD_differential_ST.R
###
# Purpose: AD differential in single nuclei ST dataset
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder
set.seed(42)

pvalcutoff <- 0.05
LFcutoff <- 0.5
pct_cutoff <- 0.05

mdat <- read_csv("/Ex_ChP/CosMX2/CosMX-projID metadata.csv") %>% column_to_rownames("projID")


##Subset
###Final cell type
Cos_ChP <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Cos_ChP.RDS", sep=""))
for(i in names(Cos_ChP@images)){Cos_ChP@images[[i]]<- NULL} #remove images before subbing

#Subcelltype
datprepfun(Cos_ChP, "Epi_1")
datprepfun(Cos_ChP, "Epi_2")
datprepfun(Cos_ChP, "Epi_3")
datprepfun(Cos_ChP, "Epi_4")

###Major cell type
#!!!!Change to Majorcelltype in datprepfun
datprepfun(Cos_ChP, "Epithelial")
datprepfun(Cos_ChP, "Endothelial")
datprepfun(Cos_ChP, "Mural")
datprepfun(Cos_ChP, "Fibroblast")
datprepfun(Cos_ChP, "Immune")

pcts <- cbind(pct_Epi_1, pct_Epi_2, pct_Epi_3, pct_Epi_4, pct_Epithelial, pct_Endo, pct_Mural, pct_Fibroblast, pct_BAM)
colnames(pcts) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4","Epithelial", "Endo", "Mural","Fibroblast", "BAM")
write.csv(pcts, paste0(curr_dir, "/limma/Percent_genes_in_cells.csv"))

##limma runs
###factors set up

#discrete var
msex = factor(paste(mdat$msex))
cogdx = factor(paste(mdat$cogdx))
ad_adnc <- factor(paste(mdat$ad_adnc))
projID <- factor(paste(rownames(mdat)))

#continuous
cogdx_c = as.numeric(mdat$cogdx)
age = as.numeric(mdat$age_death)
pmi = as.numeric(mdat$pmi)
educ = as.numeric(mdat$educ)
gpath = as.numeric(mdat$gpath)
adnc = as.numeric(mdat$adnc)
tdp = as.numeric(mdat$tdp_st4)
arteriol = as.numeric(mdat$arteriol_scler)
caa = as.numeric(mdat$caa_4gp)
cvda = as.numeric(mdat$cvda_4gp2)


###cogdx
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi)

rownames(Maindm) <- rownames(mdat) #colnames(cts_Epithelial) new versions of R/Seurat add a g infront of the participant ID
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1)

limmafun(cts_Epi_1, pct_Epi_1)
limmafun(cts_Epi_2, pct_Epi_2)
limmafun(cts_Epi_3, pct_Epi_3)
limmafun(cts_Epi_4, pct_Epi_4)

limmafun(cts_Epithelial, pct_Epithelial)
limmafun(cts_Endo, pct_Endothelial)
limmafun(cts_Mural, pct_Mural)
limmafun(cts_Fibroblast, pct_Fibroblast)
limmafun(cts_Immune, pct_Immune)



#FGSEA
#cogdx
Path_fun(limma_cogdx_ADvNCI_Epi_1)
Path_fun(limma_cogdx_ADvNCI_Epi_2)
Path_fun(limma_cogdx_ADvNCI_Epi_3)
Path_fun(limma_cogdx_ADvNCI_Epi_4)
Path_fun(limma_cogdx_ADvNCI_Endo)
Path_fun(limma_cogdx_ADvNCI_Mural)
Path_fun(limma_cogdx_ADvNCI_Epithelial)
Path_fun(limma_cogdx_ADvNCI_Fibroblast)
Path_fun(limma_cogdx_ADvNCI_BAM)

#adnc
Path_fun(limma_adnc_ADvNCI_Epi_1)
Path_fun(limma_adnc_ADvNCI_Epi_2)
Path_fun(limma_adnc_ADvNCI_Epi_3)
Path_fun(limma_adnc_ADvNCI_Epi_4)
Path_fun(limma_adnc_ADvNCI_Endo)
Path_fun(limma_adnc_ADvNCI_Mural)
Path_fun(limma_adnc_ADvNCI_Epithelial)
Path_fun(limma_adnc_ADvNCI_Fibroblast)
Path_fun(limma_adnc_ADvNCI_BAM)

#gpath
Path_fun(limma_gpath_Epi_1)
Path_fun(limma_gpath_Epi_2)
Path_fun(limma_gpath_Epi_3)
Path_fun(limma_gpath_Epi_4)
Path_fun(limma_gpath_Endo)
Path_fun(limma_gpath_Mural)
Path_fun(limma_gpath_Epithelial)
Path_fun(limma_gpath_Fibroblast)
Path_fun(limma_gpath_BAM)


#Figures DEG-GSEA
FDRcutoff <- 0.01
Pvalcutoff <- 0.05
LFcutoff <- 0.2
Epi_genes <- c("TUBB4B", "CFAP43", "DNAH12", "FOXJ1", "SERPINA3", "IL6R", "CHI3L1", "UBB", "PLTP", "NQO1", "KCNN2", "CPE", "LAMB1", "CLDN5", "COL9A3")
Fib_genes <- c("COL8A1", "COL15A1", "ADAMTS9", "CD96", "APOE", "ACSL4", "APOLD1", "KCNMA1", "SLC38A2", "TFRC", "CD9", "GPX4", "C7", "GPX3")


###Load sn

##custom Volcano
Volcano_fun(limma_cogdx_ADvNCI_Epithelial, Epi_genes)
Volcano_fun(limma_cogdx_ADvNCI_Epi_1, Epi_genes)
Volcano_fun(limma_cogdx_ADvNCI_Epi_2, Epi_genes)
Volcano_fun(limma_cogdx_ADvNCI_Epi_3, Epi_genes)
Volcano_fun(limma_cogdx_ADvNCI_Epi_4, Epi_genes)

Volcano_fun(limma_cogdx_ADvNCI_Fibroblast, Fib_genes)


##Boxplots
#Set up to rename data by disease
#Epithelial
Epi_genes <- Epi_genes[Epi_genes %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in Epi_genes){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_All.pdf", sep=""), plt, width=.7, height=.5*length(Epi_genes))

gene_list <- c("COL9A1", "EFNA5", "ITGAV")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_cellchat.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

gene_list <- c("ARL2BP", "CEP164", "CFAP58", "CEP43", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA", "DNAJB1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_Proteomic.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

#Fibroblast
Fib_genes <- Fib_genes[Fib_genes %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in Fib_genes){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_All.pdf", sep=""), plt, width=.7, height=.5*length(Fib_genes))

gene_list <- c("COL15A1", "COL18A1", "COL27A1", "COL4A1", "COL4A2", "COL4A4", "COL5A3", "COL7A1", "COL8A1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Collagen.pdf", sep=""), plt, width=.7, height=.5*length(Fib_genes))

gene_list <- c("LAMA2", "COL4A2", "COL1A2", "COL4A4", "SEMA3C", "BMP5", "JAM3", "CLDN11", "NAMPT")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Cellchat.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

gene_list <- c("ARL2BP", "CEP164", "CFAP58", "CEP43", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA", "DNAJB1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Proteomic.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

#Immune
gene_list <- c("HSPA1A", "HSP90B1", "HSPA6", "HSPA1B", "DNAJA1", "HSP90AA1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_BAM$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_BAM, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune_HSP.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

gene_list <- c("CD68", "HLA.DRA", "HLA.DRB", "HLA.DQA1", "HLA.E")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_BAM$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_BAM, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune_inflammation.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

gene_list <- c(gene_list <- c("ARL2BP", "CEP164", "CFAP58", "CEP43", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA", "DNAJB1"))
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Immune$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Immune, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune_proteomic.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

#Endothelial
gene_list <- c("LAMA5", "EFNB2", "PECAM1", "ITGA1", "ITGAV")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Endo$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Endo, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Endothelial_cellchat.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

#Mural
gene_list <- c("FN1", "COL4A2", "ITGA1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Mural$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Mural, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Mural_cellchat.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))


##GO heatmap of overlap
col_fun <- colorRamp2(c(-2, 0, 2), c("blue","white", "red"))


###Epi

selection <- c("GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOBP_CILIUM_MOVEMENT", "GOCC_CYTOPLASMIC_MICROTUBULE", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOMF_TUBULIN_BINDING", "GOCC_AXONEMAL_DYNEIN_COMPLEX", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT", "GOCC_INTRACILIARY_TRANSPORT_PARTICLE_B",
               "GOBP_EXPORT_FROM_CELL", "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOBP_EXOCYTOSIS", "GOBP_POSITIVE_REGULATION_OF_SECRETION", "GOBP_SECRETION", "GOBP_REGULATION_OF_PROTEIN_SECRETION", "GOBP_NEGATIVE_REGULATION_OF_HORMONE_SECRETION", "GOBP_REGULATION_OF_TRANSPORT", "GOCC_CATION_CHANNEL_COMPLEX", "GOCC_ION_CHANNEL_COMPLEX", "GOMF_MOLECULAR_TRANSDUCER_ACTIVITY", 
               "GOBP_REGULATION_OF_INSULIN_SECRETION", "GOBP_CELLULAR_GLUCOSE_HOMEOSTASIS", "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS", "GOBP_POSITIVE_REGULATION_OF_PROTEIN_PHOSPHORYLATION", "GOBP_GLYCOPROTEIN_BIOSYNTHETIC_PROCESS", "GOBP_REGULATION_OF_STEROL_TRANSPORT", "GOBP_LIPID_LOCALIZATION",
               "GOBP_IMMUNE_SYSTEM_PROCESS", "GOBP_INFLAMMATORY_RESPONSE", "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "GOBP_ADAPTIVE_IMMUNE_RESPONSE", "GOBP_RESPONSE_TO_CYTOKINE", "GOBP_REGULATION_OF_RESPONSE_TO_STRESS", "GOBP_REGULATION_OF_RESPONSE_TO_STRESS", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_REGULATION_OF_T_CELL_ACTIVATION", "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE",
               "GOBP_CELL_CELL_ADHESION", "GOCC_CELL_JUNCTION", "GOBP_LEUKOCYTE_CELL_CELL_ADHESION", "GOBP_CELL_MATRIX_ADHESION", "GOBP_CELL_SUBSTRATE_ADHESION", "GOCC_TIGHT_JUNCTION", "GOBP_TIGHT_JUNCTION_ORGANIZATION",
               "GOMF_SIGNALING_RECEPTOR_BINDING", "GOBP_RESPONSE_TO_GROWTH_FACTOR", "GOBP_ACTIVATION_OF_GTPASE_ACTIVITY")

##per contrast
GSEA_Epithelial_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Epithelial, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_adnc_ADvNCI_Epithelial, by="pathway", suffix=c("", "_adnc"))

GSEA_Epithelial_join <- GSEA_Epithelial_join %>% column_to_rownames("pathway")
GSEA_Epithelial_join <- GSEA_Epithelial_join[selection,]
GSEA_Epithelial_join <- GSEA_Epithelial_join[match(selection, rownames(GSEA_Epithelial_join)),]
mat <- GSEA_Epithelial_join[,c("NES", "NES_adnc")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Epithelial_join[,c("padj", "padj_adnc")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
mat2[is.na(mat2)] <- 1
mat2 <- t(mat2)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(7, "mm"), height = nrow(mat)*unit(7, "mm"), col = colorRamp2(c(-2, 0, 2), c("blue","white", "red")), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.01) {
    grid.text("*", x, y)
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Epithelial.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

#per celltype
GSEA_EpiType_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Epi_1, by="pathway", suffix=c("", "_Epi_1")) %>% full_join(GSEA_cogdx_ADvNCI_Epi_2, by="pathway", suffix=c("", "_Epi_2")) %>% full_join(GSEA_cogdx_ADvNCI_Epi_3, by="pathway", suffix=c("", "_Epi_3")) %>% full_join(GSEA_cogdx_ADvNCI_Epi_4, by="pathway", suffix=c("", "_Epi_4"))

GSEA_EpiType_join <- GSEA_EpiType_join %>% column_to_rownames("pathway")
GSEA_EpiType_join <- GSEA_EpiType_join[selection,]
GSEA_EpiType_join <- GSEA_EpiType_join[match(selection, rownames(GSEA_EpiType_join)),]
mat <- GSEA_EpiType_join[,c("NES", "NES_Epi_2", "NES_Epi_3", "NES_Epi_4")]
colnames(mat) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_EpiType_join[,c("padj", "padj_Epi_2", "padj_Epi_3", "padj_Epi_4")]
colnames(mat2) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_EpiTypes.pdf", sep = ""), font="Arial", width=w, height=h)


###Fib

selection <- c("GOBP_TIGHT_JUNCTION_ORGANIZATION", "GOCC_TIGHT_JUNCTION", "GOBP_CELL_MATRIX_ADHESION", "GOBP_CELL_ADHESION_MEDIATED_BY_INTEGRIN", "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT", "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", "GOCC_POSTSYNAPTIC_DENSITY_MEMBRANE",
               "GOMF_PASSIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOMF_TRANSMEMBRANE_TRANSPORTER_BINDING", "GOBP_ION_TRANSMEMBRANE_TRANSPORT", "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY","GOMF_INORGANIC_MOLECULAR_ENTITY_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOBP_REGULATION_OF_ION_TRANSMEMBRANE_TRANSPORT", "GOCC_ION_CHANNEL_COMPLEX", "GOCC_POTASSIUM_CHANNEL_COMPLEX", "GOCC_CATION_CHANNEL_COMPLEX", "GOBP_REGULATION_OF_ION_TRANSPORT", "GOBP_INORGANIC_ION_TRANSMEMBRANE_TRANSPORT", "GOBP_LIGAND_GATED_ION_CHANNEL_SIGNALING_PATHWAY",
               "GOBP_AMIDE_BIOSYNTHETIC_PROCESS", "GOBP_CELLULAR_AMIDE_METABOLIC_PROCESS", "GOBP_LIPID_LOCALIZATION", "GOBP_REGULATION_OF_RNA_SPLICING", "GOMF_RRNA_BINDING")

#Fibroblast
GSEA_Fibroblast_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Fibroblast, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_adnc_ADvNCI_Fibroblast, by="pathway", suffix=c("", "_adnc"))

GSEA_Fibroblast_join <- GSEA_Fibroblast_join %>% column_to_rownames("pathway")
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[selection,]
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[match(selection, rownames(GSEA_Fibroblast_join)),]
mat <- GSEA_Fibroblast_join[,c("NES", "NES_adnc")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fibroblast_join[,c("padj", "padj_adnc")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fibroblast.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()


###Fib_Neg

#Fibroblast
GSEA_Fibroblast_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Fibroblast_Neg, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_adnc_ADvNCI_Fibroblast_Neg, by="pathway", suffix=c("", "_adnc"))

GSEA_Fibroblast_join <- GSEA_Fibroblast_join %>% column_to_rownames("pathway")
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[selection,]
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[match(selection, rownames(GSEA_Fibroblast_join)),]
mat <- GSEA_Fibroblast_join[,c("NES", "NES_adnc")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fibroblast_join[,c("padj", "padj_adnc")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fibroblast_Neg.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_Neg_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()



###Fib_Pos

#Fibroblast
GSEA_Fibroblast_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Fibroblast_Pos, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_adnc_ADvNCI_Fibroblast_Pos, by="pathway", suffix=c("", "_adnc"))

GSEA_Fibroblast_join <- GSEA_Fibroblast_join %>% column_to_rownames("pathway")
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[selection,]
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[match(selection, rownames(GSEA_Fibroblast_join)),]
mat <- GSEA_Fibroblast_join[,c("NES", "NES_adnc")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fibroblast_join[,c("padj", "padj_adnc")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fibroblast_Pos.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_Pos_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

##upset plots

library(dplyr)
library(ComplexHeatmap)
library(ComplexUpset)

presence = ComplexUpset:::get_mode_presence('exclusive_intersection')
summarise_values = function(df) {
  aggregate(
    as.formula(paste0(presence, '~ intersection')),
    df,
    FUN=sum
  )
}


###Epi_2

limma_cogn_global_Epi_2_sig <- limma_cogn_global_Epi_2[limma_cogn_global_Epi_2$adj.P.Val<0.01 & abs(limma_cogn_global_Epi_2$logFC>0.15),]

cm <- make_comb_mat(list_to_matrix(list('AD(p)vNCI' = c(limma_adnc_NCIvAD_Epi_2_sig$gene), 'AD(d)vNCI' = c(limma_cogdx_ADvNCI_Epi_2_sig$gene), 'AD(d)vMCI' = c(limma_cogdx_ADvMCI_Epi_2_sig$gene), 'MCIvNCI' = c(limma_cogdx_MCIvNCI_Epi_2_sig$gene))))
plts <- UpSet(cm, set_order = rownames(cm), pt_size = unit(2, "mm"), top_annotation = upset_top_annotation(cm, height = unit(2.5, "inch"), annotation_name_rot = 90))
plts
graph2pdf(plts, file=paste(curr_dir, "/limma/upset_Epi_2.pdf", sep=""), width=7, height=5)


cm <- data.frame(list_to_matrix(list('AD(p)vNCI' = c(limma_adnc_NCIvAD_Epi_2_sig$gene), 'AD(d)vNCI' = c(limma_cogdx_ADvNCI_Epi_2_sig$gene), 'AD(d)vMCI' = c(limma_cogdx_ADvMCI_Epi_2_sig$gene), 'MCIvNCI' = c(limma_cogdx_MCIvNCI_Epi_2_sig$gene))))
plts <- upset(cm, c(colnames(cm)), name="Contrast",
              base_annotations=list('Intersection size'=(ggplot() + geom_bar() + ylab('Intersection size') + scale_y_log10())),
              width_ratio=0.1, sort_sets=F, themes=upset_default_themes(text=element_text(size=7)))
plts
graph2pdf(plts, paste(curr_dir, "/limma/upset_Epi_2_logscale.pdf", sep=""), width=7, height=5)


###Fib_3

limma_cogn_global_Fib_3_sig <- limma_cogn_global_Fib_3[limma_cogn_global_Fib_3$adj.P.Val<0.01 & abs(limma_cogn_global_Fib_3$logFC>0.1),]

cm <- make_comb_mat(list_to_matrix(list('AD(p)vNCI' = c(limma_adnc_NCIvAD_Fib_3_sig$gene), 'AD(d)vNCI' = c(limma_cogdx_ADvNCI_Fib_3_sig$gene), 'AD(d)vMCI' = c(limma_cogdx_ADvMCI_Fib_3_sig$gene), MCIvNCI = c(limma_cogdx_MCIvNCI_Fib_3_sig$gene))))
plts <- UpSet(cm, set_order = rownames(cm), pt_size = unit(2, "mm"), top_annotation = upset_top_annotation(cm, height = unit(2.5, "inch"), annotation_name_rot = 90))
plts
graph2pdf(plts, file=paste(curr_dir, "/limma/upset_Fib_3.pdf", sep=""), width=7, height=5)


cm <- data.frame(list_to_matrix(list('AD(p)vNCI' = c(limma_adnc_NCIvAD_Fib_3_sig$gene), 'AD(d)vNCI' = c(limma_cogdx_ADvNCI_Fib_3_sig$gene), 'AD(d)vMCI' = c(limma_cogdx_ADvMCI_Fib_3_sig$gene), MCIvNCI = c(limma_cogdx_MCIvNCI_Fib_3_sig$gene))))
plts <- upset(cm, c(colnames(cm)), name="Contrast",
              base_annotations=list('Intersection size'=(ggplot() + geom_bar() + ylab('Intersection size') + scale_y_log10()))
              , width_ratio=0.1, sort_sets=F, themes=upset_default_themes(text=element_text(size=7)))
plts
graph2pdf(plts, paste(curr_dir, "/limma/upset_Fib_3_logscale.pdf", sep=""), width=7, height=5)



###BAM

selection <- c("GOBP_INFLAMMATORY_RESPONSE", "GOBP_CELL_CHEMOTAXIS", "GOBP_B_CELL_MEDIATED_IMMUNITY", "GOBP_LEUKOCYTE_PROLIFERATION", "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION", "GOCC_MHC_PROTEIN_COMPLEX_BINDING",  "GOMF_MHC_CLASS_II_PROTEIN_COMPLEX_BINDING", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",  "GOMF_PEPTIDE_ANTIGEN_BINDING", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN", "GOMF_UNFOLDED_PROTEIN_BINDING")

##per contrast
GSEA_BAM_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_BAM, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_adnc_ADvNCI_BAM, by="pathway", suffix=c("", "_adnc"))

GSEA_BAM_join <- GSEA_BAM_join %>% column_to_rownames("pathway")
GSEA_BAM_join <- GSEA_BAM_join[selection,]
GSEA_BAM_join <- GSEA_BAM_join[match(selection, rownames(GSEA_BAM_join)),]
mat <- GSEA_BAM_join[,c("NES", "NES_adnc")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_BAM_join[,c("padj", "padj_adnc")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_BAM.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))

pdf(paste(curr_dir, "/FGSEA/Bar_BAM_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()


###Endo

selection <- c("GOMF_STRUCTURAL_MOLECULE_ACTIVITY", "GOMF_CELL_ADHESION_MOLECULE_BINDING", "GOCC_CELL_SUBSTRATE_JUNCTION", "GOCC_ANCHORING_JUNCTION", "GOCC_ACTIN_CYTOSKELETON", "GOBP_ACTIN_FILAMENT_BASED_PROCESS", "GOCC_ACTIN_FILAMENT_BUNDLE", "GOBP_ACTIN_FILAMENT_ORGANIZATION", "GOMF_ACTIN_BINDING", "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", "GOBP_VASCULAR_TRANSPORT", "GOCC_ION_CHANNEL_COMPLEX", "GOCC_TRANSPORT_VESICLE_MEMBRANE", "GOCC_CATION_CHANNEL_COMPLEX", "GOBP_INNATE_IMMUNE_RESPONSE", "GOBP_DEFENSE_RESPONSE", "GOMF_ANTIGEN_BINDING", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN", "GOBP_IMMUNE_RESPONSE", "GOBP_RESPONSE_TO_CYTOKINE", "GOBP_T_CELL_CYTOKINE_PRODUCTION", "GOBP_REGULATION_OF_IMMUNE_RESPONSE", "GOBP_IMMUNE_EFFECTOR_PROCESS", "GOBP_T_CELL_MEDIATED_IMMUNITY")

##per contrast
GSEA_Endo_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Endo, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_adnc_ADvNCI_Endo, by="pathway", suffix=c("", "_adnc"))

GSEA_Endo_join <- GSEA_Endo_join %>% column_to_rownames("pathway")
GSEA_Endo_join <- GSEA_Endo_join[selection,]
GSEA_Endo_join <- GSEA_Endo_join[match(selection, rownames(GSEA_Endo_join)),]
mat <- GSEA_Endo_join[,c("NES", "NES_adnc")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Endo_join[,c("padj", "padj_adnc")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Endo.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Endo_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()


###Mural

selection <- c("GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOBP_ION_TRANSPORT", "GOMF_ACTIVE_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOBP_TRANSMEMBRANE_TRANSPORT", "GOBP_METAL_ION_TRANSPORT", "GOMF_TRANSPORTER_ACTIVITY", "GOBP_CATION_TRANSPORT", "GOBP_EXTRACELLULAR_TRANSPORT", "GOBP_CATION_TRANSMEMBRANE_TRANSPORT", "GOBP_INORGANIC_ION_TRANSMEMBRANE_TRANSPORT", "GOMF_CATION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY")

##per contrast
GSEA_Mural_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Mural, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_adnc_ADvNCI_Mural, by="pathway", suffix=c("", "_adnc"))

GSEA_Mural_join <- GSEA_Mural_join %>% column_to_rownames("pathway")
GSEA_Mural_join <- GSEA_Mural_join[selection,]
GSEA_Mural_join <- GSEA_Mural_join[match(selection, rownames(GSEA_Mural_join)),]
mat <- GSEA_Mural_join[,c("NES", "NES_adnc")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Mural_join[,c("padj", "padj_adnc")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Mural.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Mural_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()


###DEP
####AD

ChP_DEP_corrected <- read_csv("/Ex_ChP/Ex_ChP5/ChP_DEP_corrected.csv")
CSF_ChP_DEP_corrected <- read_csv("/Ex_ChP/Ex_ChP5/CSF_ChP_DEP_corrected.csv")

cm <- make_comb_mat(list_to_matrix(list(Epi_1 = c(limma_cogdx_ADvNCI_Epi_1_sig$gene), Epi_2 = c(limma_cogdx_ADvNCI_Epi_2_sig$gene), Fib_1 = c(limma_cogdx_ADvNCI_Fib_1_sig$gene), Fib_2 = c(limma_cogdx_ADvNCI_Fib_2_sig$gene), Fib_3 = c(limma_cogdx_ADvNCI_Fib_3_sig$gene), DEP_CSF = c(CSF_ChP_DEP_corrected$gene), DEP_ChP = c(ChP_DEP_corrected$gene))))
plts <- UpSet(cm, set_order = rownames(cm), pt_size = unit(2, "mm"), top_annotation = upset_top_annotation(cm, height = unit(10.5, "inch"), annotation_name_rot = 90))
plts
graph2pdf(plts, file=paste(curr_dir, "/limma/upset_DEP_AD.pdf", sep=""), width=7, height=5)


cm <- data.frame(list_to_matrix(list(Epi_1 = c(limma_cogdx_ADvNCI_Epi_1_sig$gene), Epi_2 = c(limma_cogdx_ADvNCI_Epi_2_sig$gene), Fib_1 = c(limma_cogdx_ADvNCI_Fib_1_sig$gene), Fib_2 = c(limma_cogdx_ADvNCI_Fib_2_sig$gene), Fib_3 = c(limma_cogdx_ADvNCI_Fib_3_sig$gene), DEP_CSF = c(CSF_ChP_DEP_corrected$gene), DEP_ChP = c(ChP_DEP_corrected$gene))))
plts <- upset(cm, c(colnames(cm)), name="Contrast",
              base_annotations=list('Intersection size'=(ggplot() + geom_bar() + ylab('Intersection size') + scale_y_log10()))
              , width_ratio=0.1, sort_sets=F, themes=upset_default_themes(text=element_text(size=7)))
plts
graph2pdf(plts, paste(curr_dir, "/limma/upset_DEP_AD_logscale.pdf", sep=""), width=10, height=5)

