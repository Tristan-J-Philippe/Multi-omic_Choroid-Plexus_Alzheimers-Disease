###
### 033 AD_differential_RNA.R
###
# Purpose: AD differential in single nuclei RNA dataset
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(42)

FDRcutoff <- 0.01
LFcutoff <- 0.25
pct_cutoff <- 0.05

#metadata
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical).csv"))
mdat <- column_to_rownames(mdat, "projid")

mdatF <- read.csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) Female.csv"), row.names=1)
mdatM <- read.csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) Male.csv"), row.names=1)

##Subset
###Final cell type
Ex_ChP <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP.RDS"))
#Subcelltype (note change to SUbcelltype in datprepfun_lane)
datprepfun_lane(Ex_ChP, "Epi_1")
datprepfun_lane(Ex_ChP, "Epi_2")
datprepfun_lane(Ex_ChP, "Epi_3")
datprepfun_lane(Ex_ChP, "Epi_4")
datprepfun_lane(Ex_ChP, "Fib_1")
datprepfun_lane(Ex_ChP, "Fib_2")
datprepfun_lane(Ex_ChP, "Fib_3")
datprepfun_lane(Ex_ChP, "Fib_4")
datprepfun_lane(Ex_ChP, "Endo_1")
datprepfun_lane(Ex_ChP, "Endo_3")
datprepfun_lane(Ex_ChP, "Endo_2")
datprepfun_lane(Ex_ChP, "Mural_1")
datprepfun_lane(Ex_ChP, "Mural_2")
datprepfun_lane(Ex_ChP, "Mural_3")
datprepfun_lane(Ex_ChP, "BAM_1")
datprepfun_lane(Ex_ChP, "BAM_2")
datprepfun_lane(Ex_ChP, "BAM_3")
datprepfun_lane(Ex_ChP, "T_Cell")


###Major cell type

#note change to Majorcelltype in datprepfun_lane
datprepfun_lane(Ex_ChP, "Epithelial")
datprepfun_lane(Ex_ChP, "Fibroblast")
datprepfun_lane(Ex_ChP, "Endothelial")
datprepfun_lane(Ex_ChP, "Mural")
datprepfun_lane(Ex_ChP, "Immune")

pcts <- cbind(pct_Epi_1, pct_Epi_2, pct_Epi_3, pct_Epi_4, pct_Fib_1, pct_Fib_2, pct_Fib_3, pct_Fib_4, pct_Endo_1, pct_Endo_3, pct_Endo_2, pct_Mural, pct_BAM_1, pct_BAM_2, pct_BAM_3, pct_T_Cell, pct_Epithelial, pct_Fibroblast, pct_Immune)
colnames(pcts) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4", "Fib_1", "Fib_2", "Fib_3", "Fib_4", "Endo_1", "Endo_3", "Endo_2", "Mural", "BAM_1", "BAM_2", "BAM_3", "T_Cell", "Epithelial", "Fibroblast", "Immune")
write.csv(pcts, paste0(curr_dir, "/Percent_genes_in_cells.csv"))


##limma runs
###factors set up

#discrete var
msex = factor(paste(mdat$msex))
X10X_batch = factor(paste(mdat$X10X_batch))
library_batch = factor(paste(mdat$library_batch))
apoe = factor(paste(mdat$apoe_genotype))
apoe_sex = factor(paste(mdat$apoe_genotype, mdat$msex, sep="."))
cogdx = factor(paste(mdat$cogdx))
cogdx_sex = factor(paste(mdat$cogdx, mdat$msex, sep="."))
apoe_adnc <- factor(paste(mdat$apoe_genotype, mdat$ad_adnc, sep="."))
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

#sex
X10X_batch_F = factor(paste(mdatF$X10X_batch))
library_batch_F = factor(paste(mdatF$library_batch))
apoe_F = factor(paste(mdatF$apoe_genotype))
Cog_Path_F = factor(paste(mdatF$Cog_Path))
age_F = as.numeric(mdatF$age_death)

X10X_batch_M = factor(paste(mdatM$X10X_batch))
library_batch_M = factor(paste(mdatM$library_batch))
apoe_M = factor(paste(mdatM$apoe_genotype))
Cog_Path_M = factor(paste(mdatM$Cog_Path))
age_M = as.numeric(mdatM$age_death)

#cogdxAD
msex_cogdxAD = factor(paste(mdatAD$msex))
X10X_batch_cogdxAD = factor(paste(mdatAD$X10X_batch))
library_batch_cogdxAD = factor(paste(mdatAD$library_batch))
apoe_cogdxAD = factor(paste(mdatAD$apoe_genotype))
age_cogdxAD = as.numeric(mdatAD$age_death)
pmi_cogdxAD = as.numeric(mdatAD$pmi)

msex_NMCI = factor(paste(mdatNMCI$msex))
X10X_batch_NMCI = factor(paste(mdatNMCI$X10X_batch))
library_batch_NMCI = factor(paste(mdatNMCI$library_batch))
apoe_NMCI = factor(paste(mdatNMCI$apoe_genotype))
age_NMCI = as.numeric(mdatNMCI$age_death)
pmi_NMCI = as.numeric(mdatNMCI$pmi)

msex_MCI = factor(paste(mdatMCI$msex))
X10X_batch_MCI = factor(paste(mdatMCI$X10X_batch))
library_batch_MCI = factor(paste(mdatMCI$library_batch))
apoe_MCI = factor(paste(mdatMCI$apoe_genotype))
age_MCI = as.numeric(mdatMCI$age_death)
pmi_MCI = as.numeric(mdatMCI$pmi)

msex_NCI = factor(paste(mdatNCI$msex))
X10X_batch_NCI = factor(paste(mdatNCI$X10X_batch))
library_batch_NCI = factor(paste(mdatNCI$library_batch))
apoe_NCI = factor(paste(mdatNCI$apoe_genotype))
age_NCI = as.numeric(mdatNCI$age_death)
pmi_NCI = as.numeric(mdatNCI$pmi)

#ADNC
msex_ADNC0 = factor(paste(mdatADNC0$msex))
X10X_batch_ADNC0 = factor(paste(mdatADNC0$X10X_batch))
library_batch_ADNC0 = factor(paste(mdatADNC0$library_batch))
apoe_ADNC0 = factor(paste(mdatADNC0$apoe_genotype))
age_ADNC0 = as.numeric(mdatADNC0$age_death)
pmi_ADNC0 = as.numeric(mdatADNC0$pmi)

msex_ADNC1 = factor(paste(mdatADNC1$msex))
X10X_batch_ADNC1 = factor(paste(mdatADNC1$X10X_batch))
library_batch_ADNC1 = factor(paste(mdatADNC1$library_batch))
apoe_ADNC1 = factor(paste(mdatADNC1$apoe_genotype))
age_ADNC1 = as.numeric(mdatADNC1$age_death)
pmi_ADNC1 = as.numeric(mdatADNC1$pmi)

###cogdx

####AD (Epis, Fibs, and hHSP BAMs)
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + X10X_batch)

rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx4+cogdx5)/2-(cogdx2+cogdx3)/2)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Immune, pct_Immune)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)
limmafun_lane(cts_Endo, pct_Endo)
limmafun_lane(cts_Mural, pct_Mural)

#!!!!No batch correction or split
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + X10X_batch)

rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1, cogdx_ADvNCI=cogdx1-(cogdx4+cogdx5)/2, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)
limmafun(cts_T_Cell, pct_T_Cell)
limmafun_lane(cts_Immune, pct_Immune)

####within Female
options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path_F + age_F + X10X_batch_F)

rownames(Maindm) <- rownames(mdatF)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, NCIvMCI_F=Cog_Path_FNCI-Cog_Path_FMCI, ADvNCI_F=Cog_Path_FNCI-Cog_Path_FAD, MCIvAD_F=Cog_Path_FMCI-Cog_Path_FAD, N_MCIvAD_F=(Cog_Path_FNCI+Cog_Path_FMCI)/2-Cog_Path_FAD, NvMCI_AD_F=(Cog_Path_FNCI)-(Cog_Path_FMCI+Cog_Path_FAD)/2)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)
limmafun(cts_Immune, pct_Immune)

####within Male (Endo unique to males)
options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path_M + age_M + X10X_batch_M)

rownames(Maindm) <- rownames(mdatM)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, NCIvMCI_M=Cog_Path_MNCI-Cog_Path_MMCI, ADvNCI_M=Cog_Path_MNCI-Cog_Path_MAD, MCIvAD_M=Cog_Path_MMCI-Cog_Path_MAD, N_MCIvAD_M=(Cog_Path_MNCI+Cog_Path_MMCI)/2-Cog_Path_MAD, NvMCI_AD_M=(Cog_Path_MNCI)-(Cog_Path_MMCI+Cog_Path_MAD)/2)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)



###AD_adnc (Epis, Fibs, and hHSP BAMs) 
####All

options(na.action='na.pass')
Maindm <- model.matrix(~0 + ad_adnc + msex + age + X10X_batch)

rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, adnc_ADvNCI=ad_adnc1-ad_adnc0)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)


###continuous gpath (Epis, Fibs, BAM_2) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + gpath + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, gpath)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)


###continuous Braaksc (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + Braaksc + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, Braaksc)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)


###continuous tdp (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + tdp + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, tdp)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)


###continuous tangles (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + tangles + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, tangles)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)


###continuous amyloid burden (Epis, Fibs, and hHSP BAMs) 

options(na.action='na.pass')
Maindm <- model.matrix(~0 + amyl + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, amyl)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)


###continuous arteriol (Epis, Fibs, and hHSP BAMs) 
options(na.action='na.pass')
Maindm <- model.matrix(~0 + arteriol + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, arteriol)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)


###continuous caa (Epis, Fibs, and hHSP BAMs) 
options(na.action='na.pass')
Maindm <- model.matrix(~0 + caa + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, caa)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)



###continuous cvda (Epis, Fibs, and hHSP BAMs) 
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cvda + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cvda)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)


###continuous last cog score
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogn_global + msex + age + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogn_global)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Mural, pct_Mural)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)

#!!!!No batch correction or split by lane
limmafun(cts_Immune, pct_Immune)
limmafun(cts_BAM, pct_BAM)
limmafun(cts_BAM_1, pct_BAM_1)
limmafun(cts_BAM_2, pct_BAM_2)




#Pattern/Template matching

mdat <- read.csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) - X.csv"), row.names=1)
#Make factors
Cog_Path = factor(paste(mdat$Cog_Path))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path)
rownames(Maindm) <- rownames(mdat)
Maindm <- Maindm[complete.cases(Maindm),]



##Epithelial sort
#Get data
limma_cogdx_MCIvNCI_Epithelial <- limma_cogdx_MCIvNCI_Epithelial[!duplicated(limma_cogdx_MCIvNCI_Epithelial$gene),]
limma_cogdx_ADvMCI_Epithelial <- limma_cogdx_ADvMCI_Epithelial[!duplicated(limma_cogdx_ADvMCI_Epithelial$gene),]
limma_cogdx_ADvNCI_Epithelial <- limma_cogdx_ADvNCI_Epithelial[!duplicated(limma_cogdx_ADvNCI_Epithelial$gene),]

SN <- limma_cogdx_ADvNCI_Epithelial[,c(1,9:(length(limma_cogdx_ADvNCI_Epithelial)-2))]
rownames(SN) <- NULL
SN <- column_to_rownames(SN, "gene")

SN <- data.frame(t(SN))
temp <- list(substr(rownames(SN),1,9))
SN <- aggregate(SN, temp, mean)
SN <- column_to_rownames(SN, "Group.1")
SN <- data.frame(t(SN))
rownames(SN) <- limma_cogdx_ADvNCI_Epithelial$gene #put correct rownames, aggregate changes '-' to '.', thanks but no thanks?

vdatweights <- SN

#Templates
dm <- Maindm[rownames(Maindm) %in% colnames(vdatweights),]
dm <- data.frame(dm)
Temp_early <- dm$Cog_PathMCI
Temp_late <- dm$Cog_PathAD
Temp_early_late <- dm$Cog_PathNCI

#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.
SN_MCIvNCI <- data.frame(limma_cogdx_MCIvNCI_Epithelial$gene, limma_cogdx_MCIvNCI_Epithelial$P.Value, limma_cogdx_MCIvNCI_Epithelial$logFC)
colnames(SN_MCIvNCI) <- c("gene", "P.Value_MCIvNCI", "LFC_MCIvNCI")

SN_ADvMCI <- data.frame(limma_cogdx_ADvMCI_Epithelial$gene, limma_cogdx_ADvMCI_Epithelial$P.Value, limma_cogdx_ADvMCI_Epithelial$logFC)
colnames(SN_ADvMCI) <- c("gene", "P.Value_ADvMCI", "LFC_ADvMCI")

SN_ADvNCI <- data.frame(limma_cogdx_ADvNCI_Epithelial$gene, limma_cogdx_ADvNCI_Epithelial$P.Value, limma_cogdx_ADvNCI_Epithelial$logFC)
colnames(SN_ADvNCI) <- c("gene", "P.Value_ADvNCI", "LFC_ADvNCI")

stvals <- merge(merge(SN_MCIvNCI, SN_ADvMCI), SN_ADvNCI)
stvals <- column_to_rownames(stvals, "gene")
stvals <- t(stvals)

Epithelial_TM_early <- sort(stvals[1, stvals[1,]<Pvalcutoff & stvals[3,]<Pvalcutoff & stvals[5,]>Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])<LFcutoff])
Cor.fun(TM_select = Epithelial_TM_early, Template = Temp_early, tval = data.frame(limma_cogdx_MCIvNCI_Epithelial$gene, limma_cogdx_MCIvNCI_Epithelial$t))

Epithelial_TM_late <- sort(stvals[5, stvals[1,]>Pvalcutoff & stvals[3,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[2,])<LFcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])>LFcutoff])
Cor.fun(TM_select = Epithelial_TM_late, Template = Temp_late, tval = data.frame(limma_cogdx_ADvNCI_Epithelial$gene, limma_cogdx_ADvNCI_Epithelial$t))

Epithelial_TM_early_late <- sort(stvals[5, stvals[1,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[6,])>LFcutoff])
limma_cogdx_MCIvNCI_Epithelial <- limma_cogdx_MCIvNCI_Epithelial[order(match(limma_cogdx_MCIvNCI_Epithelial$gene,limma_cogdx_ADvNCI_Epithelial$gene)),]
Cor.fun(TM_select = Epithelial_TM_early_late, Template = Temp_early_late, tval = data.frame(limma_cogdx_ADvNCI_Epithelial$gene, (abs(limma_cogdx_ADvNCI_Epithelial$t)+abs(limma_cogdx_MCIvNCI_Epithelial$t))/2))

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Epithelial_TM_late_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansNCI <- data.frame(rowMeans(Sc[,c(grep("NCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansNCI) <- c("NCI")
rowMeansMCI <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansMCI) <- c("MCI")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansNCI, rowMeansMCI, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$AD),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Epithelial Late.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Epithelial_TM_late_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
as.data.frame(GO_results)

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Epithelial_TM_early_late_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansNCI <- data.frame(rowMeans(Sc[,c(grep("NCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansNCI) <- c("NCI")
rowMeansMCI <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansMCI) <- c("MCI")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansNCI, rowMeansMCI, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$AD),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Epithelial Prog.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Epithelial_TM_early_late_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
as.data.frame(GO_results)


##Fibroblast sort
#Get data
limma_cogdx_MCIvNCI_Fibroblast <- limma_cogdx_MCIvNCI_Fibroblast[!duplicated(limma_cogdx_MCIvNCI_Fibroblast$gene),]
limma_cogdx_ADvMCI_Fibroblast <- limma_cogdx_ADvMCI_Fibroblast[!duplicated(limma_cogdx_ADvMCI_Fibroblast$gene),]
limma_cogdx_ADvNCI_Fibroblast <- limma_cogdx_ADvNCI_Fibroblast[!duplicated(limma_cogdx_ADvNCI_Fibroblast$gene),]

SN <- limma_cogdx_ADvNCI_Fibroblast[,c(1,9:(length(limma_cogdx_ADvNCI_Fibroblast)-2))]
rownames(SN) <- NULL
SN <- column_to_rownames(SN, "gene")

SN <- data.frame(t(SN))
temp <- list(substr(rownames(SN),1,9))
SN <- aggregate(SN, temp, mean)
SN <- column_to_rownames(SN, "Group.1")
SN <- data.frame(t(SN))
rownames(SN) <- limma_cogdx_ADvNCI_Fibroblast$gene #put correct rownames, aggregate changes '-' to '.', thanks but no thanks?
vdatweights <- SN

#Templates
dm <- Maindm[rownames(Maindm) %in% colnames(SN),]
dm <- data.frame(dm)
Temp_early <- dm$Cog_PathMCI
Temp_late <- dm$Cog_PathAD
Temp_early_late <- dm$Cog_PathNCI

#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.
SN_MCIvNCI <- data.frame(limma_cogdx_MCIvNCI_Fibroblast$gene, limma_cogdx_MCIvNCI_Fibroblast$P.Value, limma_cogdx_MCIvNCI_Fibroblast$logFC)
colnames(SN_MCIvNCI) <- c("gene", "P.Value_MCIvNCI", "LFC_MCIvNCI")

SN_ADvMCI <- data.frame(limma_cogdx_ADvMCI_Fibroblast$gene, limma_cogdx_ADvMCI_Fibroblast$P.Value, limma_cogdx_ADvMCI_Fibroblast$logFC)
colnames(SN_ADvMCI) <- c("gene", "P.Value_ADvMCI", "LFC_ADvMCI")

SN_ADvNCI <- data.frame(limma_cogdx_ADvNCI_Fibroblast$gene, limma_cogdx_ADvNCI_Fibroblast$P.Value, limma_cogdx_ADvNCI_Fibroblast$logFC)
colnames(SN_ADvNCI) <- c("gene", "P.Value_ADvNCI", "LFC_ADvNCI")

stvals <- merge(merge(SN_MCIvNCI, SN_ADvMCI), SN_ADvNCI)
stvals <- column_to_rownames(stvals, "gene")
stvals <- t(stvals)

Fibroblast_TM_early <- sort(stvals[1, stvals[1,]<Pvalcutoff & stvals[3,]<Pvalcutoff & stvals[5,]>Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])<LFcutoff])
Cor.fun(TM_select = Fibroblast_TM_early, Template = Temp_early, tval = data.frame(limma_cogdx_MCIvNCI_Fibroblast$gene, limma_cogdx_MCIvNCI_Fibroblast$t))

Fibroblast_TM_late <- sort(stvals[5, stvals[1,]>Pvalcutoff & stvals[3,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[2,])<LFcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])>LFcutoff])
Cor.fun(TM_select = Fibroblast_TM_late, Template = Temp_late, tval = data.frame(limma_cogdx_ADvNCI_Fibroblast$gene, limma_cogdx_ADvNCI_Fibroblast$t))

Fibroblast_TM_early_late <- sort(stvals[5, stvals[1,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[6,])>LFcutoff])
limma_cogdx_MCIvNCI_Fibroblast <- limma_cogdx_MCIvNCI_Fibroblast[order(match(limma_cogdx_MCIvNCI_Fibroblast$gene,limma_cogdx_ADvNCI_Fibroblast$gene)),]
Cor.fun(TM_select = Fibroblast_TM_early_late, Template = Temp_early_late, tval = data.frame(limma_cogdx_ADvNCI_Fibroblast$gene, limma_cogdx_ADvNCI_Fibroblast$t*abs(limma_cogdx_MCIvNCI_Fibroblast$t)/2))

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Fibroblast_TM_early_late_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansNCI <- data.frame(rowMeans(Sc[,c(grep("NCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansNCI) <- c("NCI")
rowMeansMCI <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansMCI) <- c("MCI")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansNCI, rowMeansMCI, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$AD),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Fibroblast Prog.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Fibroblast_TM_early_late_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
as.data.frame(GO_results)


##Immune sort
#Get data
limma_cogdx_MCIvNCI_Immune <- limma_cogdx_MCIvNCI_Immune[!duplicated(limma_cogdx_MCIvNCI_Immune$gene),]
limma_cogdx_ADvMCI_Immune <- limma_cogdx_ADvMCI_Immune[!duplicated(limma_cogdx_ADvMCI_Immune$gene),]
limma_cogdx_ADvNCI_Immune <- limma_cogdx_ADvNCI_Immune[!duplicated(limma_cogdx_ADvNCI_Immune$gene),]

SN <- limma_cogdx_MCIvNCI_Immune[,c(1,9:(length(limma_cogdx_MCIvNCI_Immune)-2))]
rownames(SN) <- NULL
SN <- column_to_rownames(SN, "gene")

SN <- data.frame(t(SN))
temp <- list(substr(rownames(SN),1,9))
SN <- aggregate(SN, temp, mean)
SN <- column_to_rownames(SN, "Group.1")
SN <- data.frame(t(SN))
rownames(SN) <- limma_cogdx_MCIvNCI_Immune$gene #put correct rownames, aggregate changes '-' to '.', thanks but no thanks?
vdatweights <- SN

#Templates
dm <- Maindm[rownames(Maindm) %in% colnames(SN),]
dm <- data.frame(dm)
Temp_early <- dm$Cog_PathMCI
Temp_late <- dm$Cog_PathAD
Temp_early_late <- dm$Cog_PathNCI

#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.
SN_MCIvNCI <- data.frame(limma_cogdx_MCIvNCI_Immune$gene, limma_cogdx_MCIvNCI_Immune$P.Value, limma_cogdx_MCIvNCI_Immune$logFC)
colnames(SN_MCIvNCI) <- c("gene", "P.Value_MCIvNCI", "LFC_MCIvNCI")
SN_MCIvNCI <- SN_MCIvNCI[order(SN_MCIvNCI$gene),]

SN_ADvMCI <- data.frame(limma_cogdx_ADvMCI_Immune$gene, limma_cogdx_ADvMCI_Immune$P.Value, limma_cogdx_ADvMCI_Immune$logFC)
colnames(SN_ADvMCI) <- c("gene", "P.Value_ADvMCI", "LFC_ADvMCI")
SN_ADvMCI <- SN_ADvMCI[order(SN_ADvMCI$gene),]

SN_ADvNCI <- data.frame(limma_cogdx_ADvNCI_Immune$gene, limma_cogdx_ADvNCI_Immune$P.Value, limma_cogdx_ADvNCI_Immune$logFC)
colnames(SN_ADvNCI) <- c("gene", "P.Value_ADvNCI", "LFC_ADvNCI")
SN_ADvNCI <- SN_ADvNCI[order(SN_ADvNCI$gene),]

stvals <- merge(merge(SN_MCIvNCI, SN_ADvMCI), SN_ADvNCI)
stvals <- column_to_rownames(stvals, "gene")
stvals <- t(stvals)

Immune_TM_early <- sort(stvals[1, stvals[1,]<Pvalcutoff & stvals[3,]<Pvalcutoff & stvals[5,]>Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])<LFcutoff])
Cor.fun(TM_select = Immune_TM_early, Template = Temp_early, tval = data.frame(limma_cogdx_MCIvNCI_Immune$gene, limma_cogdx_MCIvNCI_Immune$t))

Immune_TM_late <- sort(stvals[5, stvals[1,]>Pvalcutoff & stvals[3,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[2,])<LFcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])>LFcutoff])
Cor.fun(TM_select = Immune_TM_late, Template = Temp_late, tval = data.frame(limma_cogdx_ADvNCI_Immune$gene, limma_cogdx_ADvNCI_Immune$t))

Immune_TM_early_late <- sort(stvals[1, stvals[1,]<Pvalcutoff  & stvals[5,]<Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[6,])>LFcutoff])
limma_cogdx_MCIvNCI_Immune <- limma_cogdx_MCIvNCI_Immune[order(match(limma_cogdx_MCIvNCI_Immune$gene,limma_cogdx_ADvNCI_Immune$gene)),]
Cor.fun(TM_select = Immune_TM_early_late, Template = Temp_early_late, tval = data.frame(limma_cogdx_ADvNCI_Immune$gene, (abs(limma_cogdx_ADvNCI_Immune$t)+limma_cogdx_MCIvNCI_Immune$t)/2))


###early

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]
rownames(SN_mdat) <- rownames(SN)
#select
genes <- unique(c(Pattern_Immune_TM_early_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansNCI <- data.frame(rowMeans(Sc[,c(grep("NCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansNCI) <- c("NCI")
rowMeansMCI <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansMCI) <- c("MCI")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansNCI, rowMeansMCI, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$MCI),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Immune Early.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Immune_TM_early_late_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
as.data.frame(GO_results)

###earlylate

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]
rownames(SN_mdat) <- rownames(SN)

#select
genes <- unique(c(Pattern_Immune_TM_early_late_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansNCI <- data.frame(rowMeans(Sc[,c(grep("NCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansNCI) <- c("NCI")
rowMeansMCI <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansMCI) <- c("MCI")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansNCI, rowMeansMCI, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$AD),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Immune Early_Late.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Immune_TM_early_late_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
as.data.frame(GO_results)

###late
SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]
rownames(SN_mdat) <- rownames(SN)

#select
genes <- unique(c(Pattern_Immune_TM_late_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansNCI <- data.frame(rowMeans(Sc[,c(grep("NCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansNCI) <- c("NCI")
rowMeansMCI <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansMCI) <- c("MCI")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansNCI, rowMeansMCI, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$AD),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Immune Late.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Immune_TM_early_late_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
as.data.frame(GO_results)



##Endothelial sort
#Get data
limma_cogdx_MCIvNCI_Endothelial <- limma_cogdx_MCIvNCI_Endothelial[!duplicated(limma_cogdx_MCIvNCI_Endothelial$gene),]
limma_cogdx_ADvMCI_Endothelial <- limma_cogdx_ADvMCI_Endothelial[!duplicated(limma_cogdx_ADvMCI_Endothelial$gene),]
limma_cogdx_ADvNCI_Endothelial <- limma_cogdx_ADvNCI_Endothelial[!duplicated(limma_cogdx_ADvNCI_Endothelial$gene),]

SN <- limma_cogdx_MCIvNCI_Endothelial[,c(1,9:(length(limma_cogdx_MCIvNCI_Endothelial)-2))]
rownames(SN) <- NULL
SN <- column_to_rownames(SN, "gene")

SN <- data.frame(t(SN))
temp <- list(substr(rownames(SN),1,9))
SN <- aggregate(SN, temp, mean)
SN <- column_to_rownames(SN, "Group.1")
SN <- data.frame(t(SN))
rownames(SN) <- limma_cogdx_MCIvNCI_Endothelial$gene #put correct rownames, aggregate changes '-' to '.', thanks but no thanks?
vdatweights <- SN

#Templates
dm <- Maindm[rownames(Maindm) %in% colnames(SN),]
dm <- data.frame(dm)
Temp_early <- dm$Cog_PathMCI
Temp_late <- dm$Cog_PathAD
Temp_early_late <- dm$Cog_PathNCI



#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.
SN_MCIvNCI <- data.frame(limma_cogdx_MCIvNCI_Endothelial$gene, limma_cogdx_MCIvNCI_Endothelial$P.Value, limma_cogdx_MCIvNCI_Endothelial$logFC)
colnames(SN_MCIvNCI) <- c("gene", "P.Value_MCIvNCI", "LFC_MCIvNCI")
SN_MCIvNCI <- SN_MCIvNCI[order(SN_MCIvNCI$gene),]

SN_ADvMCI <- data.frame(limma_cogdx_ADvMCI_Endothelial$gene, limma_cogdx_ADvMCI_Endothelial$P.Value, limma_cogdx_ADvMCI_Endothelial$logFC)
colnames(SN_ADvMCI) <- c("gene", "P.Value_ADvMCI", "LFC_ADvMCI")
SN_ADvMCI <- SN_ADvMCI[order(SN_ADvMCI$gene),]

SN_ADvNCI <- data.frame(limma_cogdx_ADvNCI_Endothelial$gene, limma_cogdx_ADvNCI_Endothelial$P.Value, limma_cogdx_ADvNCI_Endothelial$logFC)
colnames(SN_ADvNCI) <- c("gene", "P.Value_ADvNCI", "LFC_ADvNCI")
SN_ADvNCI <- SN_ADvNCI[order(SN_ADvNCI$gene),]

stvals <- merge(merge(SN_MCIvNCI, SN_ADvMCI), SN_ADvNCI)
stvals <- column_to_rownames(stvals, "gene")
stvals <- t(stvals)

Endothelial_TM_early <- sort(stvals[1, stvals[1,]<Pvalcutoff & stvals[3,]<Pvalcutoff & stvals[5,]>Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])<LFcutoff])
Cor.fun(TM_select = Endothelial_TM_early, Template = Temp_early, tval = data.frame(limma_cogdx_MCIvNCI_Endothelial$gene, limma_cogdx_MCIvNCI_Endothelial$t))

Endothelial_TM_late <- sort(stvals[5, stvals[1,]>Pvalcutoff & stvals[3,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[2,])<LFcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])>LFcutoff])
Cor.fun(TM_select = Endothelial_TM_late, Template = Temp_late, tval = data.frame(limma_cogdx_ADvNCI_Endothelial$gene, limma_cogdx_ADvNCI_Endothelial$t))

Endothelial_TM_early_late <- sort(stvals[1, stvals[1,]<Pvalcutoff  & stvals[5,]<Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[6,])>LFcutoff])
limma_cogdx_MCIvNCI_Endothelial <- limma_cogdx_MCIvNCI_Endothelial[order(match(limma_cogdx_MCIvNCI_Endothelial$gene,limma_cogdx_ADvNCI_Endothelial$gene)),]
Cor.fun(TM_select = Endothelial_TM_early_late, Template = Temp_early_late, tval = data.frame(limma_cogdx_ADvNCI_Endothelial$gene, (abs(limma_cogdx_ADvNCI_Endothelial$t)+abs(limma_cogdx_MCIvNCI_Endothelial$t))/2))

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(names(Endothelial_TM_early_late)))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansNCI <- data.frame(rowMeans(Sc[,c(grep("NCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansNCI) <- c("NCI")
rowMeansMCI <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansMCI) <- c("MCI")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansNCI, rowMeansMCI, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$MCI),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Endothelial Early_Late.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Endothelial_TM_early_late_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
as.data.frame(GO_results)


SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Endothelial_TM_late_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansNCI <- data.frame(rowMeans(Sc[,c(grep("NCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansNCI) <- c("NCI")
rowMeansMCI <- data.frame(rowMeans(Sc[,c(grep("MCI", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansMCI) <- c("MCI")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansNCI, rowMeansMCI, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$MCI),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Endothelial Late.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Endothelial_TM_early_late_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
as.data.frame(GO_results)




#FGSEA
###cogn_global
Path_fun(limma_cogn_global_Epi_1, invert=T)
Path_fun(limma_cogn_global_Epi_2, invert=T)
Path_fun(limma_cogn_global_Epi_3, invert=T)
Path_fun(limma_cogn_global_Epi_4, invert=T)
Path_fun(limma_cogn_global_Fib_1, invert=T)
Path_fun(limma_cogn_global_Fib_2, invert=T)
Path_fun(limma_cogn_global_Fib_3, invert=T)
Path_fun(limma_cogn_global_BAM_1, invert=T)
Path_fun(limma_cogn_global_BAM_2, invert=T)
Path_fun(limma_cogn_global_Mural, invert=T)
Path_fun(limma_cogn_global_Epithelial, invert=T)
Path_fun(limma_cogn_global_Fibroblast, invert=T)
Path_fun(limma_cogn_global_Endothelial, invert=T)
Path_fun(limma_cogn_global_Immune, invert=T)
Path_fun(limma_cogn_global_BAM, invert=T)

###cogdx
Path_fun(limma_cogdx_ADvNCI_Epi_1)
Path_fun(limma_cogdx_ADvNCI_Epi_2)
Path_fun(limma_cogdx_ADvNCI_Epi_3)
Path_fun(limma_cogdx_ADvNCI_Epi_4)
Path_fun(limma_cogdx_ADvNCI_Fib_1)
Path_fun(limma_cogdx_ADvNCI_Fib_2)
Path_fun(limma_cogdx_ADvNCI_Fib_3)
Path_fun(limma_cogdx_ADvNCI_BAM_1)
Path_fun(limma_cogdx_ADvNCI_BAM_2)
Path_fun(limma_cogdx_ADvNCI_Mural)
Path_fun(limma_cogdx_ADvNCI_Epithelial)
Path_fun(limma_cogdx_ADvNCI_Fibroblast)
Path_fun(limma_cogdx_ADvNCI_Endothelial)
Path_fun(limma_cogdx_ADvNCI_Immune)

Path_fun(limma_cogdx_MCIvNCI_Epi_1)
Path_fun(limma_cogdx_MCIvNCI_Epi_2)
Path_fun(limma_cogdx_MCIvNCI_Epi_3)
Path_fun(limma_cogdx_MCIvNCI_Epi_4)
Path_fun(limma_cogdx_MCIvNCI_Fib_1)
Path_fun(limma_cogdx_MCIvNCI_Fib_2)
Path_fun(limma_cogdx_MCIvNCI_Fib_3)
Path_fun(limma_cogdx_MCIvNCI_BAM_1)
Path_fun(limma_cogdx_MCIvNCI_BAM_2)
Path_fun(limma_cogdx_MCIvNCI_BAM_2)
Path_fun(limma_cogdx_MCIvNCI_Mural)
Path_fun(limma_cogdx_MCIvNCI_Epithelial)
Path_fun(limma_cogdx_MCIvNCI_Fibroblast)
Path_fun(limma_cogdx_MCIvNCI_Endothelial)
Path_fun(limma_cogdx_MCIvNCI_Immune)

Path_fun(limma_cogdx_ADvMCI_Epi_1)
Path_fun(limma_cogdx_ADvMCI_Epi_2)
Path_fun(limma_cogdx_ADvMCI_Epi_3)
Path_fun(limma_cogdx_ADvMCI_Epi_4)
Path_fun(limma_cogdx_ADvMCI_Fib_1)
Path_fun(limma_cogdx_ADvMCI_Fib_2)
Path_fun(limma_cogdx_ADvMCI_Fib_3)
Path_fun(limma_cogdx_ADvMCI_BAM_1)
Path_fun(limma_cogdx_ADvMCI_BAM_2)
Path_fun(limma_cogdx_ADvMCI_Mural)
Path_fun(limma_cogdx_ADvMCI_Epithelial)
Path_fun(limma_cogdx_ADvMCI_Fibroblast)
Path_fun(limma_cogdx_ADvMCI_Endothelial)
Path_fun(limma_cogdx_ADvMCI_Immune)


###Pattern
Path_fun(Pattern_Epithelial_TM_early)
Path_fun(Pattern_Epithelial_TM_early_late)
Path_fun(Pattern_Epithelial_TM_late)
Path_fun(Pattern_Fibroblast_TM_early)
Path_fun(Pattern_Fibroblast_TM_early_late)
Path_fun(Pattern_Fibroblast_TM_late)
Path_fun(Pattern_BAM_TM_early)
Path_fun(Pattern_BAM_TM_early_late)
Path_fun(Pattern_BAM_TM_early_prog)
Path_fun(Pattern_BAM_TM_late)

###adnc
Path_fun(limma_adnc_ADvNCI_Epi_1)
Path_fun(limma_adnc_ADvNCI_Epi_2)
Path_fun(limma_adnc_ADvNCI_Epi_3)
Path_fun(limma_adnc_ADvNCI_Epi_4)
Path_fun(limma_adnc_ADvNCI_Fib_1)
Path_fun(limma_adnc_ADvNCI_Fib_2)
Path_fun(limma_adnc_ADvNCI_Fib_3)
Path_fun(limma_adnc_ADvNCI_BAM_1)
Path_fun(limma_adnc_ADvNCI_BAM_2)
Path_fun(limma_adnc_ADvNCI_Endo_1)
Path_fun(limma_adnc_ADvNCI_Endo_3)
Path_fun(limma_adnc_ADvNCI_Endo_2)
Path_fun(limma_adnc_ADvNCI_Mural)
Path_fun(limma_adnc_ADvNCI_Epithelial)
Path_fun(limma_adnc_ADvNCI_Fibroblast)
Path_fun(limma_adnc_ADvNCI_Endothelial)
Path_fun(limma_adnc_ADvNCI_Immune)
Path_fun(limma_adnc_ADvNCI_BAM)

#gpath
Path_fun(limma_gpath_Epi_1)
Path_fun(limma_gpath_Epi_2)
Path_fun(limma_gpath_Epi_3)
Path_fun(limma_gpath_Epi_4)
Path_fun(limma_gpath_Fib_1)
Path_fun(limma_gpath_Fib_2)
Path_fun(limma_gpath_Fib_3)
Path_fun(limma_gpath_BAM_1)
Path_fun(limma_gpath_BAM_2)
Path_fun(limma_gpath_Endo_1)
Path_fun(limma_gpath_Endo_3)
Path_fun(limma_gpath_Endo_2)
Path_fun(limma_gpath_Mural)
Path_fun(limma_gpath_Epithelial)
Path_fun(limma_gpath_Fibroblast)
Path_fun(limma_gpath_Endothelial)
Path_fun(limma_gpath_Immune)
Path_fun(limma_gpath_BAM)

#tdp
Path_fun(limma_tdp_Epi_1)
Path_fun(limma_tdp_Epi_2)
Path_fun(limma_tdp_Epi_3)
Path_fun(limma_tdp_Epi_4)
Path_fun(limma_tdp_Fib_1)
Path_fun(limma_tdp_Fib_2)
Path_fun(limma_tdp_Fib_3)
Path_fun(limma_tdp_Epithelial)
Path_fun(limma_tdp_Fibroblast)
Path_fun(limma_tdp_Immune)

#tangles
Path_fun(limma_tangles_Epi_1)
Path_fun(limma_tangles_Epi_2)
Path_fun(limma_tangles_Epi_3)
Path_fun(limma_tangles_Epi_4)
Path_fun(limma_tangles_Fib_1)
Path_fun(limma_tangles_Fib_2)
Path_fun(limma_tangles_Fib_3)
Path_fun(limma_tangles_Epithelial)
Path_fun(limma_tangles_Fibroblast)
Path_fun(limma_tangles_Immune)

#amyl
Path_fun(limma_amyl_Epi_1)
Path_fun(limma_amyl_Epi_2)
Path_fun(limma_amyl_Epi_3)
Path_fun(limma_amyl_Epi_4)
Path_fun(limma_amyl_Fib_1)
Path_fun(limma_amyl_Fib_2)
Path_fun(limma_amyl_Fib_3)
Path_fun(limma_amyl_Epithelial)
Path_fun(limma_amyl_Fibroblast)
Path_fun(limma_amyl_Immune)


###Proteomics
#Prep CSF proteomics
limma_cogdx_ADvNCI_CSF <- dep_postmortemCSF[,c(1, 15:17)]
colnames(limma_cogdx_ADvNCI_CSF) <- c("gene", "logFC", "t", "P.Value")
limma_cogdx_MCIvNCI_CSF <- dep_postmortemCSF[,c(1, 19:21)]
colnames(limma_cogdx_MCIvNCI_CSF) <- c("gene", "logFC", "t", "P.Value")
limma_cogdx_ADvMCI_CSF <- dep_postmortemCSF[,c(1, 23:25)]
colnames(limma_cogdx_ADvMCI_CSF) <- c("gene", "logFC", "t", "P.Value")

#run
Path_fun(limma_ad_adnc_ADvNCI_CSF)
Path_fun(limma_cogdx_ADvNCI_CSF)
Path_fun(limma_cogdx_MCIvNCI_CSF)
Path_fun(limma_cogdx_ADvMCI_CSF)

#Prep ChP proteomics
limma_cogdx_ADvNCI_ChP <- dep_postmortemChP[,c(1, 15:17)]
colnames(limma_cogdx_ADvNCI_ChP) <- c("gene", "logFC", "t", "P.Value")
limma_cogdx_MCIvNCI_ChP <- dep_postmortemChP[,c(1, 19:21)]
colnames(limma_cogdx_MCIvNCI_ChP) <- c("gene", "logFC", "t", "P.Value")
limma_cogdx_ADvMCI_ChP <- dep_postmortemChP[,c(1, 23:25)]
colnames(limma_cogdx_ADvMCI_ChP) <- c("gene", "logFC", "t", "P.Value")

#run
Path_fun(limma_ad_adnc_ADvNCI_ChP)
Path_fun(limma_cogdx_ADvNCI_ChP)
Path_fun(limma_cogdx_MCIvNCI_ChP)
Path_fun(limma_cogdx_ADvMCI_ChP)


#Figures DEGs
FDRcutoff <- 0.01
LFcutoff <- 0.25
Epi_genes <- c("TUBB4B", "CFAP43", "DNAH12", "FOXJ1", "SERPINA3", "IL6R", "CHI3L1", "UBB", "PLTP", "NQO1", "KCNN2", "CPE", "LAMB1", "CLDN5", "COL9A3")
Fib_genes <- c("COL8A1", "COL15A1", "LAMA2", "TJP1", "COL8A1", "COL1A2", "ADAMTS9", "CD96", "APOE", "ACSL4", "APOLD1", "KCNMA1", "SLC38A2", "TFRC")
BAM_genes <- c("HSP90B1", "HSPA6", "HSPA1B", "HSP90AA1")

ChP_proteins <- c("CLN6", "CFAP58", "CEP43", "ARL2BP", "CEP164", "CFAP58",
                  "APOA4", "HSD11B1", "CHAC2", "MAN2B2", "FABP5", "ALDH9A1", "IP6K1", "KDSR", "SLC37A4", "PTGR1", "A2M", "DDO",
                  "SLC7A8", "UAP1L1", "STX2", "IPO4", "RGL1", "ITPA", "TMEM63A", "CHMP1A",
                  "IL1RAP", "HSD11B1", "PP1R1A", "GABRA3", "PARP14", "GRIA2", "HLA-B", "VRK2", "IFIH1", "SERPINF1", "ORAI3", "DNAJB1",
                  "APC", "MOB4", "TMEM98", "GATAD2A", "B3GALNT2", "EGFL8", "ICAM3")
ChP_proteins[ChP_proteins %in% limma_cogdx_ADvNCI_Fibroblast$gene]
ChP_proteins[ChP_proteins %in% limma_cogdx_ADvNCI_Epithelial_sig$gen]

CSF_proteins <- c("BCAN", "CSH1", "SCG3", "NGF", "TGM2", "RPE", "RELN", "CDH3", "CPE", "PTK7", "VASH1", "EPGN", "ARHGAP5", "AMOT",
                  "HLA-E", "MIF", "IL20RA", "B2M", "PTGES2", "IL13RA2")
CSF_proteins[CSF_proteins %in% limma_cogdx_ADvNCI_Fibroblast$gene]
CSF_proteins[CSF_proteins %in% limma_cogdx_ADvNCI_Epithelial_sig$gene]


##Venn diagram within celltype
#Epis
plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_Epi_1_sig$gene, limma_cogdx_ADvNCI_Epi_2_sig$gene, limma_cogdx_ADvNCI_Epi_3_sig$gene, limma_cogdx_ADvNCI_Epi_4_sig$gene), category.names = c("Epi_1", "Epi_2", "Epi_3", "Epi_4")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Epi_s_cogx_ADvNCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvMCI_Epi_1_sig$gene, limma_cogdx_ADvMCI_Epi_2_sig$gene, limma_cogdx_ADvMCI_Epi_3_sig$gene, limma_cogdx_ADvMCI_Epi_4_sig$gene), category.names = c("Epi_1", "Epi_2", "Epi_3", "Epi_4")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Epi_s_cogx_ADvMCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_MCIvNCI_Epi_1_sig$gene, limma_cogdx_MCIvNCI_Epi_2_sig$gene, limma_cogdx_MCIvNCI_Epi_3_sig$gene, limma_cogdx_MCIvNCI_Epi_4_sig$gene), category.names = c("Epi_1", "Epi_2", "Epi_3", "Epi_4")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Epi_s_cogx_MCIvNCI.pdf", sep=""), width=4, height=2.5)

#Fibs
plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_Fib_1_sig$gene, limma_cogdx_ADvNCI_Fib_2_sig$gene, limma_cogdx_ADvNCI_Fib_3_sig$gene), category.names = c("Fib_1", "Fib_2", "Fib_3")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Fib_s_cogx_ADvNCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvMCI_Fib_1_sig$gene, limma_cogdx_ADvMCI_Fib_2_sig$gene, limma_cogdx_ADvMCI_Fib_3_sig$gene), category.names = c("Fib_1", "Fib_2", "Fib_3")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Fib_s_cogx_ADvMCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_MCIvNCI_Fib_1_sig$gene, limma_cogdx_MCIvNCI_Fib_2_sig$gene, limma_cogdx_MCIvNCI_Fib_3_sig$gene), category.names = c("Fib_1", "Fib_2", "Fib_3")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Fib_s_cogx_MCIvNCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_Fib_1_sig$gene, limma_cogdx_ADvMCI_Fib_1_sig$gene, limma_cogdx_MCIvNCI_Fib_1_sig$gene), category.names = c("Overall", "Late", "Early")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Fib_cogdx_Fib_1.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_Fib_2_sig$gene, limma_cogdx_ADvMCI_Fib_2_sig$gene, limma_cogdx_MCIvNCI_Fib_2_sig$gene), category.names = c("Overall", "Late", "Early")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Fib_cogdx_Fib_2.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_Fib_3_sig$gene, limma_cogdx_ADvMCI_Fib_3_sig$gene, limma_cogdx_MCIvNCI_Fib_3_sig$gene), category.names = c("Overall", "Late", "Early")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Fib_cogdx_Fib_3.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_Fibroblast_sig$gene, limma_cogdx_ADvMCI_Fibroblast_sig$gene, limma_cogdx_MCIvNCI_Fibroblast_sig$gene), category.names = c("Overall", "Late", "Early")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_Fib_cogdx_Fibroblast.pdf", sep=""), width=4, height=2.5)

#DEPs
limma_cogdx_ADvNCI_ChP_sig <- limma_cogdx_ADvNCI_ChP[limma_cogdx_ADvNCI_ChP$P.Value<0.05,]
limma_cogdx_ADvMCI_ChP_sig <- limma_cogdx_ADvMCI_ChP[limma_cogdx_ADvMCI_ChP$P.Value<0.05,]
limma_cogdx_MCIvNCI_ChP_sig <- limma_cogdx_MCIvNCI_ChP[limma_cogdx_MCIvNCI_ChP$P.Value<0.05,]

limma_cogdx_ADvNCI_CSF_sig <- limma_cogdx_ADvNCI_CSF[limma_cogdx_ADvNCI_CSF$P.Value<0.05,]
limma_cogdx_ADvMCI_CSF_sig <- limma_cogdx_ADvMCI_CSF[limma_cogdx_ADvMCI_CSF$P.Value<0.05,]
limma_cogdx_MCIvNCI_CSF_sig <- limma_cogdx_MCIvNCI_CSF[limma_cogdx_MCIvNCI_CSF$P.Value<0.05,]

#SNepi-SNfib-CSF-ChP
plts <- ggVennDiagram(list(rownames(Ex_ChP), limma_cogdx_ADvNCI_ChP$gene, limma_cogdx_ADvNCI_CSF$gene), category.names = c("SN", "ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP-CSF_cogx_genes-proteins_available.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_Epithelial$gene, limma_cogdx_ADvNCI_Fibroblast$gene, limma_cogdx_ADvNCI_ChP$gene, limma_cogdx_ADvNCI_CSF$gene), category.names = c("Epithelial", "Fibroblast", "ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP-CSF_cogx_genes-proteins_expressed.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_Epithelial_sig$gene, limma_cogdx_ADvNCI_Fibroblast_sig$gene, limma_cogdx_ADvNCI_ChP_sig$gene, limma_cogdx_ADvNCI_CSF_sig$gene), category.names = c("Epithelial", "Fibroblast", "ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP_cogx_ADvNCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvMCI_Epithelial_sig$gene, limma_cogdx_ADvMCI_Fibroblast_sig$gene, limma_cogdx_ADvMCI_ChP_sig$gene, limma_cogdx_ADvMCI_CSF_sig$gene), category.names = c("Epithelial", "Fibroblast", "ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP_cogx_ADvMCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_MCIvNCI_Epithelial_sig$gene, limma_cogdx_MCIvNCI_Fibroblast_sig$gene, limma_cogdx_MCIvNCI_ChP_sig$gene, limma_cogdx_MCIvNCI_CSF_sig$gene), category.names = c("Epithelial", "Fibroblast", "ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP_cogx_MCIvNCI.pdf", sep=""), width=4, height=2.5)

#CSF-ChP only
plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_ChP$gene, limma_cogdx_ADvNCI_CSF$gene), category.names = c("ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP-CSF_cogx_proteins_available.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_ChP$gene, limma_cogdx_ADvNCI_CSF$gene), category.names = c("ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP-CSF_cogx_proteins_expressed.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvNCI_ChP_sig$gene, limma_cogdx_ADvNCI_CSF_sig$gene), category.names = c("ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP-CSF_cogx_ADvNCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_ADvMCI_ChP_sig$gene, limma_cogdx_ADvMCI_CSF_sig$gene), category.names = c("ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP-CSF_cogx_ADvMCI.pdf", sep=""), width=4, height=2.5)

plts <- ggVennDiagram(list(limma_cogdx_MCIvNCI_ChP_sig$gene, limma_cogdx_MCIvNCI_CSF_sig$gene), category.names = c("ChP", "CSF")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/limma/Venn_DEP-CSF_cogx_MCIvNCI.pdf", sep=""), width=4, height=2.5)

all_intersections <- attr(venn(list(limma_cogdx_MCIvNCI_Epithelial_sig$gene, limma_cogdx_MCIvNCI_Fibroblast_sig$gene, limma_cogdx_MCIvNCI_ChP_sig$gene, limma_cogdx_MCIvNCI_CSF_sig$gene), names =  c("Epithelial", "Fibroblast", "ChP", "CSF")), "intersections")


##custom Volcano
Volcano_fun(limma_cogdx_ADvNCI_Epithelial, Epi_genes)
Volcano_fun(limma_cogdx_ADvMCI_Epithelial, Epi_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Epithelial, Epi_genes)

Volcano_fun(limma_cogdx_ADvNCI_Epi_1, Epi_genes)
Volcano_fun(limma_cogdx_ADvMCI_Epi_1, Epi_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Epi_1, Epi_genes)

Volcano_fun(limma_cogdx_ADvNCI_Epi_2, Epi_genes)
Volcano_fun(limma_cogdx_ADvMCI_Epi_2, Epi_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Epi_2, Epi_genes)

Volcano_fun(limma_cogdx_ADvNCI_Epi_3, Epi_genes)
Volcano_fun(limma_cogdx_ADvMCI_Epi_3, Epi_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Epi_3, Epi_genes)

Volcano_fun(limma_cogdx_ADvNCI_Epi_4, Epi_genes)
Volcano_fun(limma_cogdx_ADvMCI_Epi_4, Epi_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Epi_4, Epi_genes)

Volcano_fun(limma_cogdx_ADvNCI_Fib_1, Fib_genes)
Volcano_fun(limma_cogdx_ADvMCI_Fib_1, Fib_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Fib_1, Fib_genes)

Volcano_fun(limma_cogdx_ADvNCI_Fib_2, Fib_genes)
Volcano_fun(limma_cogdx_ADvMCI_Fib_2, Fib_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Fib_2, Fib_genes)

Volcano_fun(limma_cogdx_ADvNCI_Fib_3, Fib_genes)
Volcano_fun(limma_cogdx_ADvMCI_Fib_3, Fib_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Fib_3, Fib_genes)

Volcano_fun(limma_cogdx_ADvNCI_Fib_4, Fib_genes)
Volcano_fun(limma_cogdx_ADvMCI_Fib_4, Fib_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Fib_4, Fib_genes)

Volcano_fun(limma_cogdx_ADvNCI_Fibroblast, Fib_genes)
Volcano_fun(limma_cogdx_ADvMCI_Fibroblast, Fib_genes)
Volcano_fun(limma_cogdx_MCIvNCI_Fibroblast, Fib_genes)

Volcano_fun(limma_cogdx_ADvNCI_BAM, BAM_genes)
Volcano_fun(limma_cogdx_ADvMCI_BAM, BAM_genes)
Volcano_fun(limma_cogdx_MCIvNCI_BAM, BAM_genes)


Volcano_fun(limma_cogdx_ADvNCI_ChP, ChP_proteins)
Volcano_fun(limma_cogdx_ADvMCI_ChP, ChP_proteins)
Volcano_fun(limma_cogdx_MCIvNCI_ChP, ChP_proteins)

Volcano_fun(limma_cogdx_ADvNCI_CSF, CSF_proteins)
Volcano_fun(limma_cogdx_ADvMCI_CSF, CSF_proteins)
Volcano_fun(limma_cogdx_MCIvNCI_CSF, CSF_proteins)


##Custom scatterplot
tcutoff=3.2

#Epithelial
tScat_fun(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_ADvMCI_Epithelial, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_MCIvNCI_Epithelial, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epithelial, limma_cogdx_ADvMCI_Epithelial, Epi_genes)

#Epi sub
tScat_fun(limma_cogdx_ADvNCI_Epi_1, limma_cogdx_ADvMCI_Epi_1, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epi_1, limma_cogdx_MCIvNCI_Epi_1, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epi_1, limma_cogdx_ADvMCI_Epi_1, Epi_genes)

tScat_fun(limma_cogdx_ADvNCI_Epi_2, limma_cogdx_ADvMCI_Epi_2, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epi_2, limma_cogdx_MCIvNCI_Epi_2, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epi_2, limma_cogdx_ADvMCI_Epi_2, Epi_genes)

tScat_fun(limma_cogdx_ADvNCI_Epi_3, limma_cogdx_ADvMCI_Epi_3, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epi_3, limma_cogdx_MCIvNCI_Epi_3, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epi_3, limma_cogdx_ADvMCI_Epi_3, Epi_genes)

tScat_fun(limma_cogdx_ADvNCI_Epi_4, limma_cogdx_ADvMCI_Epi_4, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epi_4, limma_cogdx_MCIvNCI_Epi_4, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epi_4, limma_cogdx_ADvMCI_Epi_4, Epi_genes)

#Fibroblast
tScat_fun(limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_ADvMCI_Fibroblast, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_MCIvNCI_Fibroblast, unique(Fib_genes))
tScat_fun(limma_cogdx_MCIvNCI_Fibroblast, limma_cogdx_ADvMCI_Fibroblast, unique(Fib_genes))

#Fibro sub
tScat_fun(limma_cogdx_MCIvNCI_Fib_1, limma_cogdx_ADvMCI_Fib_1, unique(Fib_genes))
tScat_fun(limma_cogdx_MCIvNCI_Fib_2, limma_cogdx_ADvMCI_Fib_2, unique(Fib_genes))
tScat_fun(limma_cogdx_MCIvNCI_Fib_3, limma_cogdx_ADvMCI_Fib_3, unique(Fib_genes))

tScat_fun(limma_cogdx_ADvNCI_Fib_1, limma_cogdx_MCIvNCI_Fib_1, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_2, limma_cogdx_MCIvNCI_Fib_2, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_3, limma_cogdx_MCIvNCI_Fib_3, unique(Fib_genes))

tScat_fun(limma_cogdx_ADvNCI_Fib_1, limma_cogdx_ADvMCI_Fib_1, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_2, limma_cogdx_ADvMCI_Fib_2, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_3, limma_cogdx_ADvMCI_Fib_3, unique(Fib_genes))

tScat_fun(limma_cogdx_ADvNCI_Fib_1, limma_cogdx_ADvNCI_Fib_2, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_3, limma_cogdx_ADvNCI_Fib_2, unique(Fib_genes))

#Proteomics CHP
tcutoff=2.064
tScat_fun(limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_ADvNCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_ADvMCI_Fibroblast, limma_cogdx_ADvMCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_MCIvNCI_Fibroblast, limma_cogdx_MCIvNCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_ADvNCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_ADvMCI_Epithelial, limma_cogdx_ADvMCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_MCIvNCI_Epithelial, limma_cogdx_MCIvNCI_ChP, unique(ChP_proteins))

#Proteomics CSF
tScat_fun(limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_ADvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_Fibroblast, limma_cogdx_ADvMCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_MCIvNCI_Fibroblast, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_ADvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_Epithelial, limma_cogdx_ADvMCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_MCIvNCI_Epithelial, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))

#Proteomics CHPvCSF
tScat_fun(limma_cogdx_ADvNCI_ChP, limma_cogdx_ADvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_ChP, limma_cogdx_ADvMCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_MCIvNCI_ChP, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))

#Proteomics cogdx
tScat_fun(limma_cogdx_ADvNCI_ChP, limma_cogdx_MCIvNCI_ChP, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvNCI_ChP, limma_cogdx_ADvMCI_ChP, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_ChP, limma_cogdx_MCIvNCI_ChP, unique(CSF_proteins))

tScat_fun(limma_cogdx_ADvNCI_CSF, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvNCI_CSF, limma_cogdx_ADvMCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_CSF, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))


##Boxplots
#Set up to rename data by disease
#Epithelial
Epi_genes <- Epi_genes[Epi_genes %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in Epi_genes){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_All.pdf", sep=""), plt, width=1, height=.5*length(Epi_genes))

gene_list <- c("COL9A1", "EFNA5", "ITGAV")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_cellchat.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

gene_list <- c("CEP164", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA-B", "DNAJB1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_proteomic.pdf", sep=""), plt, width=1, height=.5*length(gene_list))


#Fibroblast
gene_list <- Fib_genes
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_All.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

gene_list <- c("LAMA2", "COL4A2", "COL1A2", "COL4A4", "SEMA3C", "BMP5", "JAM3", "CLDN11", "NAMPT")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Cellchat.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

gene_list <- c("COL15A1", "COL18A1", "COL27A1", "COL4A1", "COL4A2", "COL4A4", "COL5A3", "COL7A1", "COL8A1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue","#F86247",  "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Collagen.pdf", sep=""), plt, width=.7, height=.5*length(Fib_genes))

gene_list <- c("CEP164", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA-B", "DNAJB1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Proteomic.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

#Immune
gene_list <- c("HSPA1A", "HSP90B1", "HSPA6", "HSPA1B", "DNAJA1", "HSP90AA1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_BAM$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_BAM, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

gene_list <- c("CD68", "HLA-DRA", "HLA-DRB", "HLA-DQA1", "HLA-E")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_BAM$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_BAM, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune_inflammation.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

gene_list <- c("CEP164", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA-B", "DNAJB1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_BAM$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_BAM, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune_Proteomic.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

#Endothelial
gene_list <- c("LAMA5", "EFNB2", "PECAM1", "ITGA1", "ITGAV")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Endo$gene]

plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Endo, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Endothelial_cellchat.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

#Mural
gene_list <- c("FN1", "COL4A2", "ITGA1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Mural$gene]

plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Mural, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Mural_cellchat.pdf", sep=""), plt, width=1, height=.5*length(gene_list))


#Epi2
selectPath_fun(limma_cogdx_ADvNCI_Epi_2, c("GOCC_CILIUM", "GOCC_CILIARY_TIP", "GOCC_DYNEIN_COMPLEX", "GOBP_MICROTUBULE_BUNDLE_FORMATION", "GOBP_REGULATION_OF_CILIUM_MOVEMENT", "GOBP_REGULATION_OF_MICROTUBULE_BASED_MOVEMENT", "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOBP_INTRACILIARY_TRANSPORT", "GOBP_EXPORT_FROM_CELL", "GOBP_REGULATION_OF_SECRETION", "GOCC_CATION_CHANNEL_COMPLEX", "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOMF_GATED_CHANNEL_ACTIVITY", "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS", "GOBP_PHOSPHORYLATION", "GOBP_INFLAMMATORY_RESPONSE", "GOBP_CELL_CELL_SIGNALING", "GOBP_CELL_ADHESION"))

#Fib
selectPath_fun(limma_cogdx_ADvNCI_Fib_3, c("GOCC_CELL_SUBSTRATE_JUNCTION", "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", "GOBP_CELL_ADHESION", "GOBP_REGULATION_OF_BINDING", "GOBP_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS", "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS", "GOBP_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS", "GOBP_IMMUNE_RESPONSE", "GOBP_MACROAUTOPHAGY", "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS", "GOMF_CATION_CHANNEL_ACTIVITY", "GOMF_GATED_CHANNEL_ACTIVITY", "GOCC_TRANSPORTER_COMPLEX", "GOBP_REGULATION_OF_INTRACELLULAR_TRANSPORT", "GOCC_ION_CHANNEL_COMPLEX"))

#BAM
selectPath_fun(limma_cogdx_ADvNCI_BAM_1, c("GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE", "GOBP_REGULATION_OF_STEROL_TRANSPORT", "GOBP_CELL_MOTILITY", "GOBP_CELLULAR_MACROMOLECULE_CATABOLIC_PROCESS"))


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
GSEA_Fibroblast_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fibroblast, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fibroblast, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fibroblast, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fibroblast, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fibroblast, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fibroblast, by="pathway", suffix=c("", "_gpath"))

GSEA_Fibroblast_join <- GSEA_Fibroblast_join %>% column_to_rownames("pathway")
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[selection,]
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[match(selection, rownames(GSEA_Fibroblast_join)),]
mat <- GSEA_Fibroblast_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fibroblast_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fibroblast.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_MCIvNCI), col=col_fun(mat$cogdx_MCIvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvMCI), col=col_fun(mat$cogdx_ADvMCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

#per celltype
GSEA_FibType_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Fib_1, by="pathway", suffix=c("", "_Fib_1")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_2, by="pathway", suffix=c("", "_Fib_2")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_3, by="pathway", suffix=c("", "_Fib_3"))

GSEA_FibType_join <- GSEA_FibType_join %>% column_to_rownames("pathway")
GSEA_FibType_join <- GSEA_FibType_join[selection,]
GSEA_FibType_join <- GSEA_FibType_join[match(selection, rownames(GSEA_FibType_join)),]
mat <- GSEA_FibType_join[,c("NES", "NES_Fib_2", "NES_Fib_3")]
colnames(mat) <- c("Fib_1", "Fib_2", "Fib_3")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_FibType_join[,c("padj", "padj_Fib_2", "padj_Fib_3")]
colnames(mat2) <- c("Fib_1", "Fib_2", "Fib_3")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_FibTypes.pdf", sep = ""), font="Arial", width=w, height=h)


#Fib_1
GSEA_Fib_1_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_1, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_1, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_1, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_1, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_1, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_1, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_1_join <- GSEA_Fib_1_join %>% column_to_rownames("pathway")
GSEA_Fib_1_join <- GSEA_Fib_1_join[selection,]
GSEA_Fib_1_join <- GSEA_Fib_1_join[match(selection, rownames(GSEA_Fib_1_join)),]
mat <- GSEA_Fib_1_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_1_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_1.pdf", sep = ""), font="Arial", width=w, height=h)

#Fib_2
GSEA_Fib_2_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_2, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_2, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_2, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_2, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_2, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_2, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_2_join <- GSEA_Fib_2_join %>% column_to_rownames("pathway")
GSEA_Fib_2_join <- GSEA_Fib_2_join[selection,]
GSEA_Fib_2_join <- GSEA_Fib_2_join[match(selection, rownames(GSEA_Fib_2_join)),]
mat <- GSEA_Fib_2_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_2_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_2.pdf", sep = ""), font="Arial", width=w, height=h)

#Fib_3
GSEA_Fib_3_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_3, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_3, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_3, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_3, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_3, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_3, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_3_join <- GSEA_Fib_3_join %>% column_to_rownames("pathway")
GSEA_Fib_3_join <- GSEA_Fib_3_join[selection,]
GSEA_Fib_3_join <- GSEA_Fib_3_join[match(selection, rownames(GSEA_Fib_3_join)),]
mat <- GSEA_Fib_3_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_3_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_3.pdf", sep = ""), font="Arial", width=w, height=h)

#Fib_4
GSEA_Fib_4_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_4, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_4, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_4, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_4, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_4, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_4, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_4_join <- GSEA_Fib_4_join %>% column_to_rownames("pathway")
GSEA_Fib_4_join <- GSEA_Fib_4_join[selection,]
GSEA_Fib_4_join <- GSEA_Fib_4_join[match(selection, rownames(GSEA_Fib_4_join)),]
mat <- GSEA_Fib_4_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_4_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_4.pdf", sep = ""), font="Arial", width=w, height=h)


###Fib_Neg

#Fibroblast
GSEA_Fibroblast_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fibroblast_Neg, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fibroblast_Neg, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fibroblast_Neg, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fibroblast_Neg, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fibroblast_Neg, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fibroblast_Neg, by="pathway", suffix=c("", "_gpath"))

GSEA_Fibroblast_join <- GSEA_Fibroblast_join %>% column_to_rownames("pathway")
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[selection,]
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[match(selection, rownames(GSEA_Fibroblast_join)),]
mat <- GSEA_Fibroblast_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fibroblast_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fibroblast_Neg.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_Neg_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_MCIvNCI), col=col_fun(mat$cogdx_MCIvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_Neg_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvMCI), col=col_fun(mat$cogdx_ADvMCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_Neg_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

#per celltype
GSEA_FibType_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Fib_1_Neg, by="pathway", suffix=c("", "_Fib_1")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_2_Neg, by="pathway", suffix=c("", "_Fib_2")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_3_Neg, by="pathway", suffix=c("", "_Fib_3"))

GSEA_FibType_join <- GSEA_FibType_join %>% column_to_rownames("pathway")
GSEA_FibType_join <- GSEA_FibType_join[selection,]
GSEA_FibType_join <- GSEA_FibType_join[match(selection, rownames(GSEA_FibType_join)),]
mat <- GSEA_FibType_join[,c("NES", "NES_Fib_2", "NES_Fib_3")]
colnames(mat) <- c("Fib_1", "Fib_2", "Fib_3")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_FibType_join[,c("padj", "padj_Fib_2", "padj_Fib_3")]
colnames(mat2) <- c("Fib_1", "Fib_2", "Fib_3")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_FibTypes_Neg.pdf", sep = ""), font="Arial", width=w, height=h)


#Fib_1
GSEA_Fib_1_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_1_Neg, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_1_Neg, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_1_Neg, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_1_Neg, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_1_Neg, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_1_Neg, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_1_join <- GSEA_Fib_1_join %>% column_to_rownames("pathway")
GSEA_Fib_1_join <- GSEA_Fib_1_join[selection,]
GSEA_Fib_1_join <- GSEA_Fib_1_join[match(selection, rownames(GSEA_Fib_1_join)),]
mat <- GSEA_Fib_1_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_1_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_1_Neg.pdf", sep = ""), font="Arial", width=w, height=h)

#Fib_2_Neg
GSEA_Fib_2_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_2_Neg, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_2_Neg, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_2_Neg, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_2_Neg, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_2_Neg, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_2_Neg, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_2_join <- GSEA_Fib_2_join %>% column_to_rownames("pathway")
GSEA_Fib_2_join <- GSEA_Fib_2_join[selection,]
GSEA_Fib_2_join <- GSEA_Fib_2_join[match(selection, rownames(GSEA_Fib_2_join)),]
mat <- GSEA_Fib_2_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_2_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_2_Neg.pdf", sep = ""), font="Arial", width=w, height=h)

#Fib_3
GSEA_Fib_3_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_3_Neg, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_3_Neg, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_3_Neg, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_3_Neg, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_3_Neg, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_3_Neg, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_3_join <- GSEA_Fib_3_join %>% column_to_rownames("pathway")
GSEA_Fib_3_join <- GSEA_Fib_3_join[selection,]
GSEA_Fib_3_join <- GSEA_Fib_3_join[match(selection, rownames(GSEA_Fib_3_join)),]
mat <- GSEA_Fib_3_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_3_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_3_Neg.pdf", sep = ""), font="Arial", width=w, height=h)



###Fib_Pos

#Fibroblast
GSEA_Fibroblast_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fibroblast_Pos, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fibroblast_Pos, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fibroblast_Pos, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fibroblast_Pos, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fibroblast_Pos, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fibroblast_Pos, by="pathway", suffix=c("", "_gpath"))

GSEA_Fibroblast_join <- GSEA_Fibroblast_join %>% column_to_rownames("pathway")
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[selection,]
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[match(selection, rownames(GSEA_Fibroblast_join)),]
mat <- GSEA_Fibroblast_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fibroblast_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fibroblast_Pos.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_Pos_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_MCIvNCI), col=col_fun(mat$cogdx_MCIvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_Pos_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvMCI), col=col_fun(mat$cogdx_ADvMCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_Pos_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()


#per celltype
GSEA_FibType_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Fib_1_Pos, by="pathway", suffix=c("", "_Fib_1")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_2_Pos, by="pathway", suffix=c("", "_Fib_2")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_3_Pos, by="pathway", suffix=c("", "_Fib_3"))

GSEA_FibType_join <- GSEA_FibType_join %>% column_to_rownames("pathway")
GSEA_FibType_join <- GSEA_FibType_join[selection,]
GSEA_FibType_join <- GSEA_FibType_join[match(selection, rownames(GSEA_FibType_join)),]
mat <- GSEA_FibType_join[,c("NES", "NES_Fib_2", "NES_Fib_3")]
colnames(mat) <- c("Fib_1", "Fib_2", "Fib_3")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_FibType_join[,c("padj", "padj_Fib_2", "padj_Fib_3")]
colnames(mat2) <- c("Fib_1", "Fib_2", "Fib_3")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_FibTypes_Pos.pdf", sep = ""), font="Arial", width=w, height=h)


#Fib_1
GSEA_Fib_1_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_1_Pos, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_1_Pos, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_1_Pos, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_1_Pos, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_1_Pos, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_1_Pos, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_1_join <- GSEA_Fib_1_join %>% column_to_rownames("pathway")
GSEA_Fib_1_join <- GSEA_Fib_1_join[selection,]
GSEA_Fib_1_join <- GSEA_Fib_1_join[match(selection, rownames(GSEA_Fib_1_join)),]
mat <- GSEA_Fib_1_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_1_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_1_Pos.pdf", sep = ""), font="Arial", width=w, height=h)

#Fib_2_Pos
GSEA_Fib_2_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_2_Pos, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_2_Pos, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_2_Pos, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_2_Pos, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_2_Pos, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_2_Pos, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_2_join <- GSEA_Fib_2_join %>% column_to_rownames("pathway")
GSEA_Fib_2_join <- GSEA_Fib_2_join[selection,]
GSEA_Fib_2_join <- GSEA_Fib_2_join[match(selection, rownames(GSEA_Fib_2_join)),]
mat <- GSEA_Fib_2_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_2_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_2_Pos.pdf", sep = ""), font="Arial", width=w, height=h)

#Fib_3
GSEA_Fib_3_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Fib_3_Pos, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvNCI_Fib_3_Pos, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvMCI_Fib_3_Pos, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Fib_3_Pos, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Fib_3_Pos, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Fib_3_Pos, by="pathway", suffix=c("", "_gpath"))

GSEA_Fib_3_join <- GSEA_Fib_3_join %>% column_to_rownames("pathway")
GSEA_Fib_3_join <- GSEA_Fib_3_join[selection,]
GSEA_Fib_3_join <- GSEA_Fib_3_join[match(selection, rownames(GSEA_Fib_3_join)),]
mat <- GSEA_Fib_3_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Fib_3_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fib_3_Pos.pdf", sep = ""), font="Arial", width=w, height=h)


###Immune
selection <- c("GOBP_INFLAMMATORY_RESPONSE", "GOBP_CELL_CHEMOTAXIS", "GOBP_B_CELL_MEDIATED_IMMUNITY", "GOBP_LEUKOCYTE_PROLIFERATION", "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION", "GOCC_MHC_PROTEIN_COMPLEX",  "GOMF_MHC_CLASS_II_PROTEIN_COMPLEX_BINDING", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",  "GOMF_PEPTIDE_ANTIGEN_BINDING", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN", "GOMF_UNFOLDED_PROTEIN_BINDING")

##per contrast
GSEA_BAM_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Immune, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvMCI_Immune, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvNCI_Immune, by="pathway", suffix=c("", "_ADvNCI"))

GSEA_BAM_join <- GSEA_BAM_join %>% column_to_rownames("pathway")
GSEA_BAM_join <- GSEA_BAM_join[selection,]
GSEA_BAM_join <- GSEA_BAM_join[match(selection, rownames(GSEA_BAM_join)),]
mat <- GSEA_BAM_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_BAM_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI")]
colnames(mat2) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Immune.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Immune_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_MCIvNCI), col=col_fun(mat$cogdx_MCIvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Immune_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvMCI), col=col_fun(mat$cogdx_ADvMCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Immune_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()


###Endo
selection <- c("GOMF_STRUCTURAL_MOLECULE_ACTIVITY", "GOMF_CELL_ADHESION_MOLECULE_BINDING", "GOCC_CELL_SUBSTRATE_JUNCTION", "GOCC_ANCHORING_JUNCTION", "GOCC_ACTIN_CYTOSKELETON", "GOBP_ACTIN_FILAMENT_BASED_PROCESS", "GOCC_ACTIN_FILAMENT_BUNDLE", "GOBP_ACTIN_FILAMENT_ORGANIZATION", "GOMF_ACTIN_BINDING", "GOCC_CATION_CHANNEL_COMPLEX")

##per contrast
GSEA_Endo_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Endo, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvMCI_Endo, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvNCI_Endo, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Endo, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Endo, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Endo, by="pathway", suffix=c("", "_gpath"))

GSEA_Endo_join <- GSEA_Endo_join %>% column_to_rownames("pathway")
GSEA_Endo_join <- GSEA_Endo_join[selection,]
GSEA_Endo_join <- GSEA_Endo_join[match(selection, rownames(GSEA_Endo_join)),]
mat <- GSEA_Endo_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Endo_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Endo.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Endo_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_MCIvNCI), col=col_fun(mat$cogdx_MCIvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Endo_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvMCI), col=col_fun(mat$cogdx_ADvMCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Endo_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()


#per celltype
GSEA_EndoType_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Endo_1, by="pathway", suffix=c("", "_Endo_1")) %>% full_join(GSEA_cogdx_ADvNCI_Endo_2, by="pathway", suffix=c("", "_Endo_2")) %>% full_join(GSEA_cogdx_ADvNCI_Endo_3, by="pathway", suffix=c("", "_Endo_3"))

GSEA_EndoType_join <- GSEA_EndoType_join %>% column_to_rownames("pathway")
GSEA_EndoType_join <- GSEA_EndoType_join[selection,]
GSEA_EndoType_join <- GSEA_EndoType_join[match(selection, rownames(GSEA_EndoType_join)),]
mat <- GSEA_EndoType_join[,c("NES", "NES_Endo_2", "NES_Endo_3")]
colnames(mat) <- c("Endo_1", "Endo_2", "Endo_3")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_EndoType_join[,c("padj", "padj_Endo_2", "padj_Endo_3")]
colnames(mat2) <- c("Endo_1", "Endo_2", "Endo_3")
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_EndoTypes.pdf", sep = ""), font="Arial", width=w, height=h)


###Mural
selection <- c("GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOBP_ION_TRANSPORT", "GOMF_ACTIVE_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOBP_TRANSMEMBRANE_TRANSPORT", "GOBP_METAL_ION_TRANSPORT", "GOMF_TRANSPORTER_ACTIVITY", "GOBP_CATION_TRANSPORT", "GOBP_EXTRACELLULAR_TRANSPORT", "GOBP_CATION_TRANSMEMBRANE_TRANSPORT", "GOBP_INORGANIC_ION_TRANSMEMBRANE_TRANSPORT", "GOMF_CATION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY")

##per contrast
GSEA_Mural_join <- full_join(GO_labelled, GSEA_cogdx_MCIvNCI_Mural, by="pathway", suffix=c("", "_MCIvNCI")) %>% full_join(GSEA_cogdx_ADvMCI_Mural, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_cogdx_ADvNCI_Mural, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogn_global_Mural, by="pathway", suffix=c("", "_CogDecl")) %>% full_join(GSEA_adnc_ADvNCI_Mural, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_gpath_Mural, by="pathway", suffix=c("", "_gpath"))

GSEA_Mural_join <- GSEA_Mural_join %>% column_to_rownames("pathway")
GSEA_Mural_join <- GSEA_Mural_join[selection,]
GSEA_Mural_join <- GSEA_Mural_join[match(selection, rownames(GSEA_Mural_join)),]
mat <- GSEA_Mural_join[,c("NES", "NES_ADvMCI", "NES_ADvNCI", "NES_CogDecl", "NES_adnc", "NES_gpath")]
colnames(mat) <- c("cogdx_MCIvNCI", "cogdx_ADvMCI", "cogdx_ADvNCI", "Cog_Dcl", "ADNC_ADvNCI", "gpath")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_Mural_join[,c("padj", "padj_ADvMCI", "padj_ADvNCI", "padj_CogDecl", "padj_adnc", "padj_gpath")]
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
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Mural.pdf", sep = ""), font="Arial", width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_Mural_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_MCIvNCI), col=col_fun(mat$cogdx_MCIvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Mural_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvMCI), col=col_fun(mat$cogdx_ADvMCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Mural_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$cogdx_ADvNCI), col=col_fun(mat$cogdx_ADvNCI), horiz=T, xlim=c(0,3), border=NA)
dev.off()


###CSF-DEP
selection <- c("GOCC_CILIUM", "GOMF_ACTIN_BINDING", "GOBP_AMIDE_BIOSYNTHETIC_PROCESS", "GOBP_PEPTIDE_METABOLIC_PROCESS", "GOBP_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS", "GOBP_REGULATION_OF_CATABOLIC_PROCESS", "GOMF_CHEMOKINE_ACTIVITY", "GOBP_MONOCYTE_CHEMOTAXIS")

GSEA_celltype_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_CSF, by="pathway", suffix=c("", "_ADvNCI_CSF")) %>% full_join(GSEA_cogdx_ADvMCI_CSF, by="pathway", suffix=c("", "_ADvMCI_CSF")) %>% full_join(GSEA_cogdx_MCIvNCI_CSF, by="pathway", suffix=c("", "_MCIvNCI_CSF"))

GSEA_celltype_join <- GSEA_celltype_join[!is.na(GSEA_celltype_join$pathway),]
GSEA_celltype_join <- GSEA_celltype_join %>% rownames_to_column("XX")
GSEA_celltype_join <- GSEA_celltype_join %>% column_to_rownames("pathway")
GSEA_celltype_join <- GSEA_celltype_join[selection,]
GSEA_celltype_join <- GSEA_celltype_join[match(selection, rownames(GSEA_celltype_join)),]

mat <- GSEA_celltype_join[,c(grep(pattern = "NES", x = colnames(GSEA_celltype_join), value = T))]
colnames(mat) <- c("ADvNCI_CSF", "ADvMCI_CSF", "MCIvNCI_CSF")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_celltype_join[,c(grep(pattern = "padj", x = colnames(GSEA_celltype_join), value = T))]
colnames(mat2) <- c("ADvNCI_CSF", "ADvMCI_CSF", "MCIvNCI_CSF")
mat2[is.na(mat2)] <- 1
mat2 <- t(mat2)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(7, "mm"), height = nrow(mat)*unit(7, "mm"), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y)
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/CSF heatmap.pdf", font="Arial", sep = ""), width=w, height=h)

#barplots
mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_CSF_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$MCIvNCI_CSF), col=col_fun(mat$MCIvNCI_CSF), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_CSF_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADvMCI_CSF), col=col_fun(mat$ADvMCI_CSF), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_CSF_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADvNCI_CSF), col=col_fun(mat$ADvNCI_CSF), horiz=T, xlim=c(0,3), border=NA)
dev.off()


###CHP-DEP
selection <- c(
  "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", 
  "GOBP_B_CELL_MEDIATED_IMMUNITY", "GOBP_LEUKOCYTE_PROLIFERATION", "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION", 
  "GOCC_MHC_PROTEIN_COMPLEX",  "GOMF_MHC_CLASS_II_PROTEIN_COMPLEX_BINDING", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",  "GOMF_PEPTIDE_ANTIGEN_BINDING", "GOBP_REGULATION_OF_T_CELL_ACTIVATION", "GOBP_T_CELL_CYTOKINE_PRODUCTION",
  "GOBP_RESPONSE_TO_CYTOKINE", "GOBP_CELL_CHEMOTAXIS", 
  "GOMF_UNFOLDED_PROTEIN_BINDING",
  
  "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
  "GOBP_CELL_MATRIX_ADHESION", "GOMF_CELL_ADHESION_MOLECULE_BINDING", "GOMF_STRUCTURAL_MOLECULE_ACTIVITY", 
  "GOBP_CELL_ADHESION_MEDIATED_BY_INTEGRIN",
  "GOCC_ACTIN_CYTOSKELETON", "GOBP_ACTIN_FILAMENT_ORGANIZATION",
  "GOBP_TIGHT_JUNCTION_ORGANIZATION", "GOCC_TIGHT_JUNCTION", "GOCC_CELL_SUBSTRATE_JUNCTION", "GOCC_ANCHORING_JUNCTION", 
  "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",
  
  "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOBP_REGULATION_OF_PROTEIN_SECRETION",
  "GOMF_PASSIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
  "GOCC_ION_CHANNEL_COMPLEX", "GOBP_REGULATION_OF_ION_TRANSPORT",
  "GOBP_VASCULAR_TRANSPORT",
  
  "GOBP_REGULATION_OF_INSULIN_SECRETION", "GOBP_CELLULAR_GLUCOSE_HOMEOSTASIS",
  "GOBP_REGULATION_OF_STEROL_TRANSPORT", "GOBP_LIPID_LOCALIZATION",
  "GOBP_AMIDE_BIOSYNTHETIC_PROCESS",  
  "GOBP_RESPONSE_TO_GROWTH_FACTOR", "GOBP_REGULATION_OF_RNA_SPLICING",
  
  "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOBP_CILIUM_MOVEMENT", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT")

GSEA_celltype_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_ChP, by="pathway", suffix=c("", "_ADvNCI_ChP")) %>% full_join(GSEA_cogdx_ADvMCI_ChP, by="pathway", suffix=c("", "_ADvMCI_ChP")) %>% full_join(GSEA_cogdx_MCIvNCI_ChP, by="pathway", suffix=c("", "_MCIvNCI_ChP")) %>% full_join(GSEA_cogdx_ADvNCI_Immune, by="pathway", suffix=c("", "Immune")) %>% full_join(GSEA_cogdx_ADvNCI_Fibroblast, by="pathway", suffix=c("", "Fibroblast")) %>% full_join(GSEA_cogdx_ADvNCI_Endo, by="pathway", suffix=c("", "Endothelial")) %>% full_join(GSEA_cogdx_ADvNCI_Mural, by="pathway", suffix=c("", "Mural")) %>% full_join(GSEA_cogdx_ADvNCI_Epithelial, by="pathway", suffix=c("", "Epithelial"))
write_csv(GSEA_celltype_join, paste(curr_dir, "/FGSEA/Celltype_fulljoin.csv", sep=""))

GSEA_celltype_join <- GSEA_celltype_join[!is.na(GSEA_celltype_join$pathway),]
GSEA_celltype_join <- GSEA_celltype_join %>% rownames_to_column("XX")
GSEA_celltype_join <- GSEA_celltype_join %>% column_to_rownames("pathway")
GSEA_celltype_join <- GSEA_celltype_join[selection,]
GSEA_celltype_join <- GSEA_celltype_join[match(selection, rownames(GSEA_celltype_join)),]

mat <- GSEA_celltype_join[,c(grep(pattern = "NES", x = colnames(GSEA_celltype_join), value = T))]
colnames(mat) <- c("ADvNCI_ChP", "ADvMCI_ChP", "MCIvNCI_ChP", "Immune", "Fibroblast", "Endothelial", "Mural", "Epithelial")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_celltype_join[,c(grep(pattern = "padj", x = colnames(GSEA_celltype_join), value = T))]
colnames(mat2) <- c("ADvNCI_ChP", "ADvMCI_ChP", "MCIvNCI_ChP", "Immune", "Fibroblast", "Endothelial", "Mural", "Epithelial")
mat2[is.na(mat2)] <- 1
mat2 <- t(mat2)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(7, "mm"), height = nrow(mat)*unit(7, "mm"), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y)
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_DEP-celltype.pdf", font="Arial", sep = ""), width=w, height=h)

#barplots
mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_ChP_MCIvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$MCIvNCI_ChP), col=col_fun(mat$MCIvNCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_ChP_ADvMCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADvMCI_ChP), col=col_fun(mat$ADvMCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_ChP_ADvNCI.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADvNCI_ChP), col=col_fun(mat$ADvNCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()


####Positive
GSEA_celltype_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_ChP_Pos, by="pathway", suffix=c("", "_ADvNCI_ChP")) %>% full_join(GSEA_cogdx_ADvMCI_ChP_Pos, by="pathway", suffix=c("", "_ADvMCI_ChP")) %>% full_join(GSEA_cogdx_MCIvNCI_ChP_Pos, by="pathway", suffix=c("", "_MCIvNCI_ChP")) %>% full_join(GSEA_cogdx_ADvNCI_Immune_Pos, by="pathway", suffix=c("", "Immune")) %>% full_join(GSEA_cogdx_ADvNCI_Fibroblast_Pos, by="pathway", suffix=c("", "Fibroblast")) %>% full_join(GSEA_cogdx_ADvNCI_Endo_Pos, by="pathway", suffix=c("", "Endothelial")) %>% full_join(GSEA_cogdx_ADvNCI_Mural_Pos, by="pathway", suffix=c("", "Mural")) %>% full_join(GSEA_cogdx_ADvNCI_Epithelial_Pos, by="pathway", suffix=c("", "Epithelial"))
write_csv(GSEA_celltype_join, paste(curr_dir, "/FGSEA/Celltype_fulljoin.csv", sep=""))

GSEA_celltype_join <- GSEA_celltype_join[!is.na(GSEA_celltype_join$pathway),]
GSEA_celltype_join <- GSEA_celltype_join %>% rownames_to_column("XX")
GSEA_celltype_join <- GSEA_celltype_join %>% column_to_rownames("pathway")
GSEA_celltype_join <- GSEA_celltype_join[selection,]
GSEA_celltype_join <- GSEA_celltype_join[match(selection, rownames(GSEA_celltype_join)),]

mat <- GSEA_celltype_join[,c(grep(pattern = "NES", x = colnames(GSEA_celltype_join), value = T))]
colnames(mat) <- c("ADvNCI_ChP", "ADvMCI_ChP", "MCIvNCI_ChP", "Immune", "Fibroblast", "Endothelial", "Mural", "Epithelial")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_celltype_join[,c(grep(pattern = "padj", x = colnames(GSEA_celltype_join), value = T))]
colnames(mat2) <- c("ADvNCI_ChP", "ADvMCI_ChP", "MCIvNCI_ChP", "Immune", "Fibroblast", "Endothelial", "Mural", "Epithelial")
mat2[is.na(mat2)] <- 1
mat2 <- t(mat2)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(7, "mm"), height = nrow(mat)*unit(7, "mm"), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y)
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_DEP-celltype_Pos.pdf", font="Arial", sep = ""), width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_ChP_MCIvNCI_Pos.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$MCIvNCI_ChP), col=col_fun(mat$MCIvNCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_ChP_ADvMCI_Pos.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADvMCI_ChP), col=col_fun(mat$ADvMCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_ChP_ADvNCI_Pos.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADvNCI_ChP), col=col_fun(mat$ADvNCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()


####Negative
GSEA_celltype_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_ChP_Neg, by="pathway", suffix=c("", "_ADvNCI_ChP")) %>% full_join(GSEA_cogdx_ADvMCI_ChP_Neg, by="pathway", suffix=c("", "_ADvMCI_ChP")) %>% full_join(GSEA_cogdx_MCIvNCI_ChP_Neg, by="pathway", suffix=c("", "_MCIvNCI_ChP")) %>% full_join(GSEA_cogdx_ADvNCI_Immune_Neg, by="pathway", suffix=c("", "Immune")) %>% full_join(GSEA_cogdx_ADvNCI_Fibroblast_Neg, by="pathway", suffix=c("", "Fibroblast")) %>% full_join(GSEA_cogdx_ADvNCI_Endo_Neg, by="pathway", suffix=c("", "Endothelial")) %>% full_join(GSEA_cogdx_ADvNCI_Mural_Neg, by="pathway", suffix=c("", "Mural")) %>% full_join(GSEA_cogdx_ADvNCI_Epithelial_Neg, by="pathway", suffix=c("", "Epithelial"))
write_csv(GSEA_celltype_join, paste(curr_dir, "/FGSEA/Celltype_fulljoin.csv", sep=""))

GSEA_celltype_join <- GSEA_celltype_join[!is.na(GSEA_celltype_join$pathway),]
GSEA_celltype_join <- GSEA_celltype_join %>% rownames_to_column("XX")
GSEA_celltype_join <- GSEA_celltype_join %>% column_to_rownames("pathway")
GSEA_celltype_join <- GSEA_celltype_join[selection,]
GSEA_celltype_join <- GSEA_celltype_join[match(selection, rownames(GSEA_celltype_join)),]

mat <- GSEA_celltype_join[,c(grep(pattern = "NES", x = colnames(GSEA_celltype_join), value = T))]
colnames(mat) <- c("ADvNCI_ChP", "ADvMCI_ChP", "MCIvNCI_ChP", "Immune", "Fibroblast", "Endothelial", "Mural", "Epithelial")
mat[is.na(mat)] <- 0
mat <- t(mat)

mat2 <- GSEA_celltype_join[,c(grep(pattern = "padj", x = colnames(GSEA_celltype_join), value = T))]
colnames(mat2) <- c("ADvNCI_ChP", "ADvMCI_ChP", "MCIvNCI_ChP", "Immune", "Fibroblast", "Endothelial", "Mural", "Epithelial")
mat2[is.na(mat2)] <- 1
mat2 <- t(mat2)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(7, "mm"), height = nrow(mat)*unit(7, "mm"), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y)
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_DEP-celltype_Neg.pdf", font="Arial", sep = ""), width=w, height=h)

mat <- data.frame(t(mat))
pdf(paste(curr_dir, "/FGSEA/Bar_ChP_MCIvNCI_Neg.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$MCIvNCI_ChP), col=col_fun(mat$MCIvNCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_ChP_ADvMCI_Neg.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADvMCI_ChP), col=col_fun(mat$ADvMCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_ChP_ADvNCI_Neg.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADvNCI_ChP), col=col_fun(mat$ADvNCI_ChP), horiz=T, xlim=c(0,3), border=NA)
dev.off()

#Parenchyma modules
###Load modules
Ast <- read_csv(paste0(curr_dir, "/GSVA_modules/Ast.csv"))
Ext <- read_csv(paste0(curr_dir, "/GSVA_modules/Ext.csv"))
Inh <- read_csv(paste0(curr_dir, "/GSVA_modules/Inh.csv"))
Mic <- read_csv(paste0(curr_dir, "/GSVA_modules/Mic.csv"))
Oli <- read_csv(paste0(curr_dir, "/GSVA_modules/Oli.csv"))
ST_mods <- read_csv(paste0(curr_dir, "/GSVA_modules/ST modules.csv"))

Paren_mods <- list("ast_m19"=c(Ast[Ast$module_assignment=="ast_m19",]$gene_name),
                   "ext_m3"=c(Ext[Ext$module_assignment=="ext_m3",]$gene_name),
                   "ext_m23"=c(Ext[Ext$module_assignment=="ext_m23",]$gene_name),
                   "ext_m8"=c(Ext[Ext$module_assignment=="ext_m8",]$gene_name),
                   "ext_m29"=c(Ext[Ext$module_assignment=="ext_m29",]$gene_name),
                   "ext_m28"=c(Ext[Ext$module_assignment=="ext_m28",]$gene_name),
                   "inh_m28"=c(Inh[Inh$module_assignment=="inh_m28",]$gene_name),
                   "inh_m31"=c(Inh[Inh$module_assignment=="inh_m31",]$gene_name),
                   "mic_m46"=c(Mic[Mic$module_assignment=="mic_m46",]$gene_name),
                   "oli_m26"=c(Oli[Oli$module_assignment=="oli_m26",]$gene_name),
                   "oli_m14"=c(Oli[Oli$module_assignment=="oli_m14",]$gene_name),
                   "DAM"=c(ST_mods[ST_mods$module=="DAM",]$symbol))

GSVA_fun(limma_cogdx_ADvNCI_Epithelial, pathways=c(Epi_modules, Paren_mods))
GSVA_fun(limma_cogdx_ADvNCI_Immune, pathways=c(Imm_modules, Paren_mods))

###Factors set up
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical).csv"))
mdat <- column_to_rownames(mdat, "projid")
msex = factor(paste(mdat$msex))
X10X_batch = factor(paste(mdat$X10X_batch))
library_batch = factor(paste(mdat$library_batch))
apoe = factor(paste(mdat$apoe_genotype))
apoe_sex = factor(paste(mdat$apoe_genotype, mdat$msex, sep="."))
cogdx = factor(paste(mdat$cogdx))
cogdx_sex = factor(paste(mdat$cogdx, mdat$msex, sep="."))
apoe_adnc <- factor(paste(mdat$apoe_genotype, mdat$ad_adnc, sep="."))
ad_adnc <- factor(paste(mdat$ad_adnc))

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


####AD cogdx
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + X10X_batch)

rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx4+cogdx5)/2-(cogdx2+cogdx3)/2)

GSVA_limmafun_lane(GSVA_Epithelial)
GSVA_limmafun_lane(GSVA_Immune)