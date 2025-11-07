###
### 036 Cell_type_proportion_ST.R
###
# Purpose: Cell type proportion analysis in spatial transcriptomics dataset
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder

##Cell proportion
###Per subtype
#load up to date meta
mdat <- read_csv(paste0(curr_dir, "/ST-projID metadata.csv", col_type=cols(.default="c")))

#cell counts
md <- data.frame(Cos_ChP@meta.data$Subcelltype, Cos_ChP@meta.data$projID)
colnames(md) <- c("active.ident", "projID")
prop <- tabyl(md, active.ident, projID)
prop <- prop %>% column_to_rownames("active.ident")
parenchyma <- prop[c("Parenchyma"),]
parenchyma <- colSums(parenchyma)
prop <- prop[!rownames(prop)%in%c("Parenchyma"),]

#discrete metadata
Cos_ChP_md <- left_join(md, mdat, by="projID")
Cos_ChP_md$active.ident <- NULL
Cos_ChP_md <- distinct(Cos_ChP_md, .keep_all = T)

#Factors
msex = factor(paste(Cos_ChP_md$msex))
X10X_batch = factor(paste(Cos_ChP_md$X10X_batch))
library_batch = factor(paste(Cos_ChP_md$library_batch))
apoe = factor(paste(Cos_ChP_md$apoe_genotype))
cogdx = factor(paste(Cos_ChP_md$cogdx))
Cog_Path = factor(paste(Cos_ChP_md$Cog_Path))
ad_adnc <- factor(paste(Cos_ChP_md$ad_adnc))

#continuous
age = as.numeric(Cos_ChP_md$age_death)
pmi = as.numeric(Cos_ChP_md$pmi)
educ = as.numeric(Cos_ChP_md$educ)
gpath = as.numeric(Cos_ChP_md$gpath)
cogn_global = as.numeric(Cos_ChP_md$cogn_global_last)
parenchyma = as.numeric(parenchyma)
tdp_st4 = as.numeric(Cos_ChP_md$tdp_st4)

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)
rownames(Maindm) <- colnames(prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
prop <- prop[,colnames(prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1)

vdat <- voom(prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Proportions/cogdx_ADvNCI.csv", sep=""))


###per subtype within major
####Epi

md <- data.frame(Cos_ChP@meta.data$Subcelltype, Cos_ChP@meta.data$projID)
colnames(md) <- c("active.ident", "projID")
prop <- tabyl(md, active.ident, projID)
prop <- prop %>% column_to_rownames("active.ident")
parenchyma <- prop[c("Parenchyma"),]
parenchyma <- colSums(parenchyma)
prop <- prop[!rownames(prop)%in%c("Parenchyma"),]

Epi_prop <- prop[c(0:4),]

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)
rownames(Maindm) <- colnames(Epi_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Epi_prop <- Epi_prop[,colnames(Epi_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1)

vdat <- voom(Epi_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_Episub.csv", sep=""))


####Fib

md <- data.frame(Cos_ChP@meta.data$Subcelltype, Cos_ChP@meta.data$projID)
colnames(md) <- c("active.ident", "projID")
prop <- tabyl(md, active.ident, projID)
prop <- prop %>% column_to_rownames("active.ident")
parenchyma <- prop[c("Parenchyma"),]
parenchyma <- colSums(parenchyma)
prop <- prop[!rownames(prop)%in%c("Parenchyma"),]


Fib_prop <- prop[c(5:7),]

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)
rownames(Maindm) <- colnames(Fib_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Fib_prop <- Fib_prop[,colnames(Fib_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1)

vdat <- voom(Fib_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_Fibsub.csv", sep=""))


####Endo

Endo_prop <- prop[c(8:10),]
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(Endo_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Endo_prop <- Endo_prop[,colnames(Endo_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1)

vdat <- voom(Endo_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_sub_Endo.csv", sep=""))


####Mural

Mural_prop <- prop[c(11:13),]
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi)

rownames(Maindm) <- colnames(Mural_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Mural_prop <- Mural_prop[,colnames(Mural_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1)

vdat <- voom(Mural_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_sub_Mural.csv", sep=""))


####Imm

Imm_prop <- prop[c(14:17),]

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)
rownames(Maindm) <- colnames(Imm_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Imm_prop <- Imm_prop[,colnames(Imm_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1)

vdat <- voom(Imm_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_sub.csv", sep=""))


####Per major celltype
#load up to date meta
mdat <- read_csv(paste0(curr_dir, "/CosMX-projID metadata.csv", col_type=cols(.default="c")))
#load up to date meta

#cell counts
md <- data.frame(Cos_ChP@meta.data$Majorcelltype, Cos_ChP@meta.data$projID)
colnames(md) <- c("active.ident", "projID")
prop <- tabyl(md, active.ident, projID)
prop <- prop %>% column_to_rownames("active.ident")
parenchyma <- prop[c("Parenchyma"),]

#discrete metadata
Cos_ChP_md <- left_join(md, mdat, by="projID")
Cos_ChP_md$active.ident <- NULL
Cos_ChP_md <- distinct(Cos_ChP_md, .keep_all = T)

#factors
msex = factor(paste(Cos_ChP_md$msex))
X10X_batch = factor(paste(Cos_ChP_md$X10X_batch))
library_batch = factor(paste(Cos_ChP_md$library_batch))
apoe = factor(paste(Cos_ChP_md$apoe_genotype))
cogdx = factor(paste(Cos_ChP_md$cogdx))
ad_adnc <- factor(paste(Cos_ChP_md$ad_adnc))

#continuous
age = as.numeric(Cos_ChP_md$age_death)
pmi = as.numeric(Cos_ChP_md$pmi)
educ = as.numeric(Cos_ChP_md$educ)
cogn_global = as.numeric(Cos_ChP_md$cogn_global_last)
parenchyma = as.numeric(parenchyma)

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
prop <- prop[,colnames(prop)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1)

vdat <- voom(prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)
#dir.create(paste(curr_dir, "/CellType/Proportions/", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_mjr.csv", sep=""))
