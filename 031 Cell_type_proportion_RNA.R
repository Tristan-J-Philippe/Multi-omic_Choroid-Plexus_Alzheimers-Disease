###
### 031 Cell_type_proportion_RNA.R
###
# Purpose: Cell type proportion analysis in single nuclei RNA dataset
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(1234)
##Cell proportion
####Per subtype

#load up to date meta
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical).csv", col_type=cols(.default="c")))

#cell counts
md <- data.frame(Ex_ChP@meta.data$Subcelltype, Ex_ChP@meta.data$Lane_projid, Ex_ChP@meta.data$projid)
colnames(md) <- c("active.ident", "Lane_projid", "projid")
prop <- tabyl(md, active.ident, Lane_projid)
prop <- prop %>% column_to_rownames("active.ident")
parenchyma <- prop[c("Neuro", "Astro"),]
parenchyma <- colSums(parenchyma)
prop <- prop[!rownames(prop)%in%c("Neuro", "Astro", "Oligo"),]

#discrete metadata
Ex_ChP_md <- left_join(md, mdat, by="projid")
Ex_ChP_md$active.ident <- NULL
Ex_ChP_md <- distinct(Ex_ChP_md, .keep_all = T)

#Factors
msex = factor(paste(Ex_ChP_md$msex))
X10X_batch = factor(paste(Ex_ChP_md$X10X_batch))
library_batch = factor(paste(Ex_ChP_md$library_batch))
apoe = factor(paste(Ex_ChP_md$apoe_genotype))
cogdx = factor(paste(Ex_ChP_md$cogdx))
Cog_Path = factor(paste(Ex_ChP_md$Cog_Path))
ad_adnc <- factor(paste(Ex_ChP_md$ad_adnc))

#continuous
age = as.numeric(Ex_ChP_md$age_death)
pmi = as.numeric(Ex_ChP_md$pmi)
educ = as.numeric(Ex_ChP_md$educ)
cogn_global = as.numeric(Ex_ChP_md$cogn_global_last)
parenchyma = as.numeric(parenchyma)
tdp_st4 = as.numeric(Ex_ChP_md$tdp_st4)

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
prop <- prop[,colnames(prop)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1)
#Maincm <- makeContrasts(levels=Maindm, tdp_st4)
#Maincm <- makeContrasts(levels=Maindm, ad_adnc_ADvNCI=(ad_adnc1-ad_adnc0))
#Maincm <- makeContrasts(levels=Maindm, cogn_global)

vdat <- voom(prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Lane <- factor(paste(sapply(str_split(colnames(prop), "_",),'[', 2), sapply(str_split(colnames(prop),"_",),'[', 3), sep="_")) #get Lane names
corfit <- duplicateCorrelation(vdat, block=Lane)
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm, block=Lane, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_MCIvNCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvMCI.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvMCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_MCIvNCI.csv", sep=""))

prop <- data.frame(t(prop))
prop <- prop/rowSums(prop)*100
#prop <- data.frame(t(vdata))
#rownames(prop) <- gsub("X", "", rownames(prop))
prop$projid <- sapply(str_split(rownames(prop), "_",),'[', 1)
prop <- left_join(prop, mdat, by="projid")

plts <- ggplot(prop, aes(Epi_1, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epi_1_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Epi_2, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epi_2_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Epi_3, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epi_3_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Epi_4, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epi_4_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Fib_1, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fib_1_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Fib_2, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fib_2_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Fib_3, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fib_3_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Fib_4, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fib_4_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Mural, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Mural_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Endo, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Endo_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(BAM_1, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/BAM_1_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(BAM_2, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/BAM_2_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(BAM_3, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/BAM_3_cogdx.pdf", sep=""), width=.8, height=1.2)

prop_heat <- as.matrix(prop[,0:13])
rownames(prop_heat) <- prop$cogn_global_last
prop_heat <- prop_heat[order(rownames(prop_heat), decreasing=T),] 
prop_heat <- (scale(prop_heat))
Heatmap(t(prop_heat), heatmap_legend_param = list(title = "Proportion"), cluster_rows = T, cluster_columns = F, show_column_names = T, show_row_names = TRUE, col = magma(256))


###Imm

Imm_prop <- data.frame(t(prop[c(15:18),]))
Imm_prop <- Imm_prop/colSums(Imm_prop)*100

plts <- ggplot(Imm_prop, aes(BAM_1, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/BAM_1_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(Imm_prop, aes(BAM_2, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/BAM_2_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(Imm_prop, aes(BAM_3, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/BAM_3_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

Imm_prop <- prop[c(15:18),]
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(Imm_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Imm_prop <- Imm_prop[,colnames(Imm_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1)

vdat <- voom(Imm_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Lane <- factor(paste(sapply(str_split(colnames(Imm_prop), "_",),'[', 2), sapply(str_split(colnames(Imm_prop),"_",),'[', 3), sep="_")) #get Lane names
corfit <- duplicateCorrelation(vdat, block=Lane)
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm, block=Lane, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

#dir.create(paste(curr_dir, "/CellType/Publish/Proportions/", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_sub_Imm.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_MCIvNCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvMCI_sub_Imm.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvMCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_MCIvNCI_sub_Imm.csv", sep=""))

###Epi

Epi_prop <- data.frame(t(prop[c(5:8),]))
Epi_prop <- Epi_prop/colSums(Epi_prop)*100

plts <- ggplot(Epi_prop, aes(Epi_1, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epi_1_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(Epi_prop, aes(Epi_2, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epi_2_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(Epi_prop, aes(Epi_3, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epi_3_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(Epi_prop, aes(Epi_4, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epi_4_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

Epi_prop <- prop[c(5:8),]
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(Epi_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Epi_prop <- Epi_prop[,colnames(Epi_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, tdp_st4)
Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1)
Maincm <- makeContrasts(levels=Maindm, ad_adnc_ADvNCI=(ad_adnc1-ad_adnc0))
Maincm <- makeContrasts(levels=Maindm, cogn_global)

vdat <- voom(Epi_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Lane <- factor(paste(sapply(str_split(colnames(Epi_prop), "_",),'[', 2), sapply(str_split(colnames(Imm_prop),"_",),'[', 3), sep="_")) #get Lane names
corfit <- duplicateCorrelation(vdat, block=Lane)
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm, block=Lane, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

#dir.create(paste(curr_dir, "/CellType/Publish/Proportions/", sep=""))
df <- topTable(Mainf, coef = "tdp_st4", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/tdp_st4_sub_Epi.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_sub_Epi.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_MCIvNCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvMCI_sub_Epi.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvMCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_MCIvNCI_sub_Epi.csv", sep=""))
df <- topTable(Mainf, coef = "ad_adnc_ADvNCI", p.value=1, number=Inf) #BAM_2
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/ad_adnc_ADvNCI_sub_Epi.csv", sep=""))
df <- topTable(Mainf, coef = "cogn_global", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/_cogn_global_sub_Epi.csv", sep=""))


###Fib

Fib_prop <- data.frame(t(prop[c(5:7),]))
Fib_prop <- Fib_prop/colSums(Fib_prop)*100
plts <- ggplot(Fib_prop, aes(Fib_1, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fib_1_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(Fib_prop, aes(Fib_2, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fib_2_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(Fib_prop, aes(Fib_3, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fib_3_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(Fib_prop, aes(Fib_4, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fib_4_cogdx_sub.pdf", sep=""), width=.8, height=1.2)

Fib_prop <- prop[c(5:7),]
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(Fib_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Fib_prop <- Fib_prop[,colnames(Fib_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1)

vdat <- voom(Fib_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Lane <- factor(paste(sapply(str_split(colnames(Fib_prop), "_",),'[', 2), sapply(str_split(colnames(Imm_prop),"_",),'[', 3), sep="_")) #get Lane names
corfit <- duplicateCorrelation(vdat, block=Lane)
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm, block=Lane, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

#dir.create(paste(curr_dir, "/CellType/Publish/Proportions/", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_sub_Fib.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_MCIvNCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvMCI_sub_Fib.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvMCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_MCIvNCI_sub_Fib.csv", sep=""))

###Endo

Endo_prop <- prop[c(8:10),]
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(Endo_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Endo_prop <- Endo_prop[,colnames(Endo_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1)

vdat <- voom(Endo_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Lane <- factor(paste(sapply(str_split(colnames(Endo_prop), "_",),'[', 2), sapply(str_split(colnames(Endo_prop),"_",),'[', 3), sep="_")) #get Lane names
corfit <- duplicateCorrelation(vdat, block=Lane)
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm, block=Lane, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

#dir.create(paste(curr_dir, "/CellType/Publish/Proportions/", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_sub_Endo.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_MCIvNCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvMCI_sub_Endo.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvMCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_MCIvNCI_sub_Endo.csv", sep=""))


###Mural

Mural_prop <- prop[c(12:14),]
options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(Mural_prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Mural_prop <- Mural_prop[,colnames(Mural_prop)%in%rownames(Maindm)]

Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1)

vdat <- voom(Mural_prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Lane <- factor(paste(sapply(str_split(colnames(Mural_prop), "_",),'[', 2), sapply(str_split(colnames(Mural_prop),"_",),'[', 3), sep="_")) #get Lane names
corfit <- duplicateCorrelation(vdat, block=Lane)
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm, block=Lane, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

#dir.create(paste(curr_dir, "/CellType/Publish/Proportions/", sep=""))
df <- topTable(Mainf, coef = "tdp_st4", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/tdp_st4_sub_Mural.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_sub_Mural.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_MCIvNCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvMCI_sub_Mural.csv", sep=""))


####Per major celltype
#load up to date meta
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical).csv", col_type=cols(.default="c")))
#load up to date meta

#cell counts
md <- data.frame(Ex_ChP@meta.data$Majorcelltype, Ex_ChP@meta.data$Lane_projid, Ex_ChP@meta.data$projid)
colnames(md) <- c("active.ident", "Lane_projid", "projid")
prop <- tabyl(md, active.ident, Lane_projid)
prop <- prop %>% column_to_rownames("active.ident")
parenchyma <- prop[c("Parenchyma"),]

#discrete metadata
Ex_ChP_md <- left_join(md, mdat, by="projid")
Ex_ChP_md$active.ident <- NULL
Ex_ChP_md <- distinct(Ex_ChP_md, .keep_all = T)

#factors
msex = factor(paste(Ex_ChP_md$msex))
X10X_batch = factor(paste(Ex_ChP_md$X10X_batch))
library_batch = factor(paste(Ex_ChP_md$library_batch))
apoe = factor(paste(Ex_ChP_md$apoe_genotype))
cogdx = factor(paste(Ex_ChP_md$cogdx))
ad_adnc <- factor(paste(Ex_ChP_md$ad_adnc))

#continuous
age = as.numeric(Ex_ChP_md$age_death)
pmi = as.numeric(Ex_ChP_md$pmi)
educ = as.numeric(Ex_ChP_md$educ)
cogn_global = as.numeric(Ex_ChP_md$cogn_global_last)
parenchyma = as.numeric(parenchyma)

options(na.action='na.pass')
Maindm <- model.matrix(~0 + cogdx + msex + age + pmi + parenchyma)

rownames(Maindm) <- colnames(prop)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
prop <- prop[,colnames(prop)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1)

vdat <- voom(prop, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Lane <- factor(paste(sapply(str_split(colnames(prop), "_",),'[', 2), sapply(str_split(colnames(prop),"_",),'[', 3), sep="_")) #get Lane names
corfit <- duplicateCorrelation(vdat, block=Lane)
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm, block=Lane, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)
#dir.create(paste(curr_dir, "/CellType/Publish/Proportions/", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvNCI_mjr.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_MCIvNCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_ADvMCI_mjr.csv", sep=""))
df <- topTable(Mainf, coef = "cogdx_ADvMCI", p.value=1, number=Inf)
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions/cogdx_MCIvNCI_mjr.csv", sep=""))

#prop <- data.frame(t(rbind(props$Proportions))*100)
prop <- data.frame(t(vdata))
rownames(prop) <- gsub("X", "", rownames(prop))
prop$projid <- sapply(str_split(rownames(prop), "_",),'[', 1)
prop <- left_join(prop, mdat, by="projid")

plts <- ggplot(prop, aes(Fibroblast, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Fibroblast_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Epithelial, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Epithelial_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Immune, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Immune_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Endothelial, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Endothelial_cogdx.pdf", sep=""), width=.8, height=1.2)

plts <- ggplot(prop, aes(Mural, Cog_Path)) + geom_boxplot() + geom_point() + coord_flip() + scale_y_discrete(limit=c("NCI","MCI","AD")) + labs(x="Normalized Cell Count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7))
plts
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Proportions/Mural_cogdx.pdf", sep=""), width=.8, height=1.2)



#Get metadata
Ex_ChP_md <- data.frame(Imm@meta.data$Subcelltype, Imm@meta.data$Lane_projid, Imm@meta.data$projid)
colnames(Imm_md) <- c("Subcelltype", "Lane_projid", "projid")
Imm_md <- left_join(Imm_md, mdat, by="projid") #add up to date metadata

#propeller
propeller(cluster=Imm_md$Subcelltype, sample=Imm_md$projid, group=Imm_md$Cog_Path)
propeller(cluster=Imm_md$Subcelltype, sample=Imm_md$Lane_projid, group=Imm_md$Cog_Path) #Lane as replicate

#better to get it out and do a proper limmma model
props <- getTransformedProps(clusters = Imm_md$Subcelltype, sample = Imm_md$Lane_projid)

Imm_md$Subcelltype <- NULL
Imm_md <- distinct(Imm_md, Lane_projid, .keep_all = T)
#discrete var
msex = factor(paste(Imm_md$msex))
X10X_batch = factor(paste(Imm_md$X10X_batch))
library_batch = factor(paste(Imm_md$library_batch))
apoe = factor(paste(Imm_md$apoe_genotype))
cogdx = factor(paste(Imm_md$cogdx))
ad_adnc <- factor(paste(Imm_md$ad_adnc))

#continuous
age = as.numeric(Imm_md$age_death)
pmi = as.numeric(Imm_md$pmi)
educ = as.numeric(Imm_md$educ)
cogn_global = as.numeric(Imm_md$cogn_global_last)

Maindm <- model.matrix(~0 + cogdx + msex + age + X10X_batch)
prop <- props$TransformedProps

rownames(Maindm) <- rownames(Imm_md)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdx_ADvNCI=(cogdx4+cogdx5)/2-cogdx1, cogdx_ADvMCI=(cogdx2+cogdx3)/2-(cogdx4+cogdx5)/2, cogdx_MCIvNCI=(cogdx2+cogdx3)/2-cogdx1)
#Maincm <- makeContrasts(levels=Maindm, ad_adnc_ADvNCI=(ad_adnc1-ad_adnc0))
#Maincm <- makeContrasts(levels=Maindm, cogn_global)
#prop <- prop[prop$sample %in% cogn_global,]
#Maincm <- makeContrasts(levels=Maindm, apoe34v33=(apoe34-apoe33))


#limma
Lane <- factor(paste(sapply(str_split(Imm_md$Lane_projid,"_",),'[', 2), sapply(str_split(Imm_md$Lane_projid,"_",),'[', 3), sep="_")) #get Lane names
corfit <- duplicateCorrelation(prop, block=Lane)
Mainf <- lmFit(prop, Maindm, block=Lane, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)
topTable(Mainf, coef = "cogdx_ADvNCI", p.value=1)
topTable(Mainf, coef = "cogdx_MCIvNCI", p.value=1) #BAM_1 trend
topTable(Mainf, coef = "cogdx_ADvMCI", p.value=1)