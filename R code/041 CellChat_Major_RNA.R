###
### 041 CellChat_Major_RNA.R
###
# Purpose: AD differential in single nuclei ST dataset
# Dependencies:
source("/040 CellChat_util_RNA.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder

Ex_ChP <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Ex_ChP_sans_paren_EMF.RDS", sep=""))
#Ex_ChP <- subset(Ex_ChP, Majorcelltype!="Parenchyma")
Idents(Ex_ChP) <- Ex_ChP@meta.data$Majorcelltype
levels(Ex_ChP) <- c("Epithelial", "Fibroblast", "Mural", "Endothelial", "Immune")
Idents(Ex_ChP) -> Ex_ChP@meta.data$Majorcelltype

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
levels(Ex_ChP) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4", "Fib_1", "Fib_2", "Fib_3", "Endo_1", "Endo_2", "Endo_3", "Mural_1", "Mural_2", "Mural_3", "BAM_1", "BAM_2", "BAM_3", "T_Cell")
Idents(Ex_ChP) -> Ex_ChP@meta.data$Subcelltype

#Cellchat
##RunCC
CCfun(Ex_ChP, "Majorcelltype")
gc()

###per group
Ex_ChP@meta.data$cogdx[which(Ex_ChP@meta.data$cogdx == "5")] <- "4"
Ex_ChP@meta.data$cogdx[which(Ex_ChP@meta.data$cogdx == "3")] <- "2"
Ex_ChPNCI <- subset(Ex_ChP, subset=(cogdx=="1"))
Ex_ChPMCI <- subset(Ex_ChP, subset=(cogdx=="2" | cogdx=="3"))
Ex_ChPAD <- subset(Ex_ChP, subset=(cogdx=="4" | cogdx=="5"))
rm(Ex_ChP)
gc()

###RunCC
CCfun(Ex_ChPNCI, "Majorcelltype")
rm(Ex_ChPNCI)
gc()
CCfun(Ex_ChPMCI, "Majorcelltype")
rm(Ex_ChPMCI)
gc()
CCfun(Ex_ChPAD, "Majorcelltype")
rm(Ex_ChPAD)
gc()


###Venns
####Interaction

interact_NCI <- unique(paste(Ex_ChPNCI_df.net$interaction_name))
interact_MCI <- unique(paste(Ex_ChPMCI_df.net$interaction_name))
interact_AD <- unique(paste(Ex_ChPAD_df.net$interaction_name))
plts <- ggVennDiagram(list(interact_NCI, interact_MCI, interact_AD), category.names = c("NCI", "MCI", "AD")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_Interaction_Comparison.png", sep=""), width=4, height=2.5)

NCI_MCI_AD <- Reduce(intersect, list(interact_NCI, interact_MCI, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_MCI_AD)]
A_MCI <- interact_MCI[!(interact_MCI %in% NCI_MCI_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_MCI_AD)]
Intersects <- list(NCI_MCI_AD, A_NCI, A_MCI, A_AD)
names(Intersects) <- c("NCI_MCI_AD", "NCI", "MCI", "AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_Interaction_NCI_AD.txt", sep=""))


####Pathway

interact_NCI <- unique(paste(Ex_ChPNCI_df.net$pathway_name))
interact_MCI <- unique(paste(Ex_ChPMCI_df.net$pathway_name))
interact_AD <- unique(paste(Ex_ChPAD_df.net$pathway_name))
plts <- ggVennDiagram(list(interact_NCI, interact_MCI, interact_AD), category.names = c("NCI", "MCI", "AD")) + theme(legend.position = "none", text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_pathway_Comparison.png", sep=""), width=4, height=2.5)

NCI_MCI_AD <- Reduce(intersect, list(interact_NCI, interact_MCI, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_MCI_AD)]
A_MCI <- interact_MCI[!(interact_MCI %in% NCI_MCI_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_MCI_AD)]
Intersects <- list(NCI_MCI_AD, A_NCI, A_MCI, A_AD)
names(Intersects) <- c("NCI_MCI_AD", "NCI", "MCI", "AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_pathway_NCI_AD.txt", sep=""))


####ligand-receptor

interact_NCI <- unique(c(paste(Ex_ChPNCI_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
interact_MCI <- unique(c(paste(Ex_ChPMCI_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
interact_AD <- unique(c(paste(Ex_ChPAD_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
plts <- ggVennDiagram(list(interact_NCI, interact_MCI, interact_AD), category.names = c("NCI", "MCI", "AD")) + theme(legend.position = "none", text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_ligand-receptor_Comparison.png", sep=""), width=4, height=2.5)

NCI_MCI_AD <- Reduce(intersect, list(interact_NCI, interact_MCI, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_MCI_AD)]
A_MCI <- interact_MCI[!(interact_MCI %in% NCI_MCI_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_MCI_AD)]
Intersects <- list(NCI_MCI_AD, A_NCI, A_MCI, A_AD)
names(Intersects) <- c("NCI_MCI_AD", "A_NCI", "A_MCI", "A_AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_ligand-receptor_NCI_AD.txt", sep=""))



##Merge, compare and DE
###Merge and compare per condition

#merge results
Ex_ChP_merge_CC <- mergeCellChat(list(NCI=Ex_ChPNCI_CC, MCI=Ex_ChPMCI_CC, AD=Ex_ChPAD_CC), add.names=names(list(NCI=Ex_ChPNCI_CC, MCI=Ex_ChPMCI_CC, AD=Ex_ChPAD_CC)))

# perform differential expression analysis
Ex_ChP_merge_CC <- identifyOverExpressedGenes(Ex_ChP_merge_CC, group.dataset = "datasets", pos.dataset = "AD", features.name = "AD", only.pos = FALSE) #!!!!pos.dataset set AD as default dataset for positive fold change #Keep all results
net <- netMappingDEG(Ex_ChP_merge_CC, features.name = "AD") #select DE

#select up and down regulated pathways with a logFC of .1
net.up <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "AD", ligand.logFC = 0.1, ligand.pvalues=0.01, receptor.logFC = 0.1, receptor.pvalues=0.01) #!!!!up in "ad"
net.down <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "AD", ligand.logFC = -0.1, ligand.pvalues=0.01, receptor.logFC = -0.1, receptor.pvalues=0.01) #!!!!up in "NCI" = down in AD
net_sig <- rbind(net.up, net.down)

# perform differential expression analysis
Ex_ChP_merge_CC <- identifyOverExpressedGenes(Ex_ChP_merge_CC, group.dataset = "datasets", pos.dataset = "MCI", features.name = "MCI", only.pos = FALSE) #!!!!pos.dataset set AD as default dataset for positive fold change 
net <- netMappingDEG(Ex_ChP_merge_CC, features.name = "MCI") #select DE

#select up and down regulated pathways with a logFC of .1
net.up <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "MCI", ligand.logFC = 0.1, ligand.pvalues=0.01, receptor.logFC = 0.1, receptor.pvalues=0.01) #!!!!up in "ad"
gene.up <- extractGeneSubsetFromPair(net.up, Ex_ChP_merge_CC) #get genes
net.down <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "MCI", ligand.logFC = -0.1, ligand.pvalues=0.01, receptor.logFC = -0.1, receptor.pvalues=0.01) #!!!!up in "NCI" = down in AD
gene.down <- extractGeneSubsetFromPair(net.down, Ex_ChP_merge_CC)
net_sig <- rbind(net_sig, net.up, net.down)

saveRDS(Ex_ChP_merge_CC, paste(curr_dir, "/CellChat/Ex_ChP_merge_CC_ADNCI.RDS", sep=""))



##DEG with limma instead
#get net
net <- subsetCommunication(Ex_ChP_merge_CC)


net <- rbind(net$NCI, net$MCI, net$AD)
net <- net[,c("source", "target", "ligand", "receptor", "prob", "pval", "interaction_name", "pathway_name", "annotation", "evidence")]

#get genes
net <- bind_cols(net, data.frame(t(sapply(str_split(net$receptor,"_",),'[',1:2))))
colnames(net) <- c("source", "target", "ligand", "receptor", "prob", "pval", "interaction_name", "pathway_name", "annotation", "evidence", "receptorA", "receptorB")


genes <- c(net$ligand, net$receptorA, net$receptorB)
genes <- unique(genes)
genes <- genes[!is.na(genes)]

cutoff_pval <- 0.05
cutoff_logFC <- 0.2


###NCIvAD

#make limma_res
limma_cogdx_ADvNCI_Epithelial <- cbind(Celltype="Epithelial", limma_cogdx_ADvNCI_Epithelial)
limma_cogdx_ADvNCI_Fibroblast <- cbind(Celltype="Fibroblast", limma_cogdx_ADvNCI_Fibroblast)
limma_cogdx_ADvNCI_Endothelial <- cbind(Celltype="Endothelial", limma_cogdx_ADvNCI_Endothelial)
limma_cogdx_ADvNCI_Mural <- cbind(Celltype="Mural", limma_cogdx_ADvNCI_Mural)
limma_cogdx_ADvNCI_Immune <- cbind(Celltype="Immune", limma_cogdx_ADvNCI_Immune)

limma_res <- bind_rows(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_ADvNCI_Endothelial, limma_cogdx_ADvNCI_Mural, limma_cogdx_ADvNCI_Immune)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdx_ADvNCI <- DEG_source_fun(limma_res, net)
net_cogdx_ADvNCI_sig <- net_cogdx_ADvNCI[net_cogdx_ADvNCI$P.Value_Ligand<cutoff_pval & abs(net_cogdx_ADvNCI$logFC_Ligand)>cutoff_logFC,]  #& c(net_cogdx_ADvMCI$P.Value_ReceptorA<cutoff_pval | net_cogdx_ADvMCI$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdx_ADvMCI$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdx_ADvMCI$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdx_ADvNCI_sig <- net_cogdx_ADvNCI_sig[complete.cases(net_cogdx_ADvNCI_sig$source),]
view(net_cogdx_ADvNCI_sig)



###NCIvMCI

#make limma_res
limma_cogdx_MCIvNCI_Epithelial <- cbind(Celltype="Epithelial", limma_cogdx_MCIvNCI_Epithelial)
limma_cogdx_MCIvNCI_Fibroblast <- cbind(Celltype="Fibroblast", limma_cogdx_MCIvNCI_Fibroblast)
limma_cogdx_MCIvNCI_Endothelial <- cbind(Celltype="Endothelial", limma_cogdx_MCIvNCI_Endothelial)
limma_cogdx_MCIvNCI_Mural <- cbind(Celltype="Mural", limma_cogdx_MCIvNCI_Mural)
limma_cogdx_MCIvNCI_Immune <- cbind(Celltype="Immune", limma_cogdx_MCIvNCI_Immune)

limma_res <- bind_rows(limma_cogdx_MCIvNCI_Epithelial, limma_cogdx_MCIvNCI_Fibroblast, limma_cogdx_MCIvNCI_Endothelial, limma_cogdx_MCIvNCI_Mural, limma_cogdx_MCIvNCI_Immune)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdx_MCIvNCI <- DEG_source_fun(limma_res, net)
net_cogdx_MCIvNCI_sig <- net_cogdx_MCIvNCI[net_cogdx_MCIvNCI$P.Value_Ligand<cutoff_pval & abs(net_cogdx_MCIvNCI$logFC_Ligand)>cutoff_logFC,] #& c(net_cogdx_ADvMCI$P.Value_ReceptorA<cutoff_pval | net_cogdx_ADvMCI$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdx_ADvMCI$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdx_ADvMCI$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdx_MCIvNCI_sig <- net_cogdx_MCIvNCI_sig[complete.cases(net_cogdx_MCIvNCI_sig$source),]
view(net_cogdx_MCIvNCI_sig)



###MCIvAD

#make limma_res
limma_cogdx_ADvMCI_Epithelial <- cbind(Celltype="Epithelial", limma_cogdx_ADvMCI_Epithelial)
limma_cogdx_ADvMCI_Fibroblast <- cbind(Celltype="Fibroblast", limma_cogdx_ADvMCI_Fibroblast)
limma_cogdx_ADvMCI_Endothelial <- cbind(Celltype="Endothelial", limma_cogdx_ADvMCI_Endothelial)
limma_cogdx_ADvMCI_Mural <- cbind(Celltype="Mural", limma_cogdx_ADvMCI_Mural)
limma_cogdx_ADvMCI_Immune <- cbind(Celltype="Immune", limma_cogdx_ADvMCI_Immune)

limma_res <- bind_rows(limma_cogdx_ADvMCI_Epithelial, limma_cogdx_ADvMCI_Fibroblast, limma_cogdx_ADvMCI_Endothelial, limma_cogdx_ADvMCI_Mural, limma_cogdx_ADvMCI_Immune)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdx_ADvMCI <- DEG_source_fun(limma_res, net)
net_cogdx_ADvMCI_sig <- net_cogdx_ADvMCI[net_cogdx_ADvMCI$P.Value_Ligand<cutoff_pval & abs(net_cogdx_ADvMCI$logFC_Ligand)>cutoff_logFC ,] #& c(net_cogdx_ADvMCI$P.Value_ReceptorA<cutoff_pval | net_cogdx_ADvMCI$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdx_ADvMCI$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdx_ADvMCI$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdx_ADvMCI_sig <- net_cogdx_ADvMCI_sig[complete.cases(net_cogdx_ADvMCI_sig$source),]
view(net_cogdx_ADvMCI_sig)



###Aggregate net_sigs

net <- rbind(net_cogdx_ADvNCI, net_cogdx_ADvMCI, net_cogdx_MCIvNCI)
write.csv(net, paste0(curr_dir, "/CellChat/Comparison/net_pathways.csv"))
net_sig <- rbind(net_cogdx_ADvNCI_sig, net_cogdx_ADvMCI_sig, net_cogdx_MCIvNCI_sig)
write.csv(net_sig, paste0(curr_dir, "/CellChat/Comparison/net_sig_pathways.csv"))


###Aggregate by pathway

ligandlogFC <- aggregate(net_cogdx_ADvNCI$logFC_Ligand, by=list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdx_ADvNCI$logFC_ReceptorA, list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdx_ADvNCI$logFC_ReceptorB, list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_ReceptorB), mean, all=T)
combinedlogFC_ADvNCI <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_ADvNCI <- cbind(combinedlogFC_ADvNCI[,1:2], rowMeans(combinedlogFC_ADvNCI[3:5], na.rm = TRUE))
view(combinedlogFC_ADvNCI)

ligandlogFC <- aggregate(net_cogdx_ADvMCI$logFC_Ligand, by=list(net_cogdx_ADvMCI$pathway_name, net_cogdx_ADvMCI$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdx_ADvMCI$logFC_ReceptorA, list(net_cogdx_ADvMCI$pathway_name, net_cogdx_ADvMCI$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdx_ADvMCI$logFC_ReceptorB, list(net_cogdx_ADvMCI$pathway_name, net_cogdx_ADvMCI$Celltype_ReceptorB), mean, all=T)
combinedlogFC_ADvMCI <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_ADvMCI <- cbind(combinedlogFC_ADvMCI[,1:2], rowMeans(combinedlogFC_ADvMCI[3:5], na.rm = TRUE))
view(combinedlogFC_ADvMCI)

ligandlogFC <- aggregate(net_cogdx_MCIvNCI$logFC_Ligand, by=list(net_cogdx_MCIvNCI$pathway_name, net_cogdx_MCIvNCI$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdx_MCIvNCI$logFC_ReceptorA, list(net_cogdx_MCIvNCI$pathway_name, net_cogdx_MCIvNCI$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdx_MCIvNCI$logFC_ReceptorB, list(net_cogdx_MCIvNCI$pathway_name, net_cogdx_MCIvNCI$Celltype_ReceptorB), mean, all=T)
combinedlogFC_MCIvNCI <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_MCIvNCI <- cbind(combinedlogFC_MCIvNCI[,1:2], rowMeans(combinedlogFC_MCIvNCI[3:5], na.rm = TRUE))
view(combinedlogFC_MCIvNCI)


##Plot differences
###heatmap comps

#barplots
plts <- compareInteractions(Ex_ChP_merge_CC, show.legend = F, group=c(1,2,3), color.use=c("skyblue", "#F86247", "red3")) + theme(text=element_text(size=7)) #Barplot of counts
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_bar.pdf", sep=""), width=2, height=2.5)

plts <- compareInteractions(Ex_ChP_merge_CC, show.legend = F, group=c(1,2,3), measure="weight", color.use=c("skyblue", "#F86247", "red3")) + theme(text=element_text(size=7)) #Barplot of weights
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_bar.pdf", sep=""), width=2, height=2.5)

#Network
plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,2)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_AD-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,2), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_AD-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(2,3)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_AD-MCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(2,3), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_AD-MCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,3)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_MCI-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,3), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_MCI-NCI.pdf", sep=""))

signaling_NCI <- dimnames(Ex_ChP_merge_CC@net$NCI$prob)[[3]]
signaling_MCI <- dimnames(Ex_ChP_merge_CC@net$MCI$prob)[[3]]
signaling_AD <- dimnames(Ex_ChP_merge_CC@net$AD$prob)[[3]]
signaling_shared <- intersect(signaling_AD, intersect(signaling_MCI, signaling_NCI))
signaling_sig <- unique(net_cogdx_ADvNCI_sig$interaction_name) # 284 sig pathways generated using CC_limma.R script
signaling_final <- intersect(signaling_sig, signaling_shared)

Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("NCI","AD"), signaling=signaling_final, n_sample=21, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"), group.color=c("skyblue", "red3"), col.max=0.1)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("NCI","AD"), signaling=signaling_final, n_sample=21, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"), group.color=c("skyblue", "red3"), col.max=8)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("NCI","MCI"), signaling=signaling_final, n_sample=21, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"), group.color=c("skyblue", "#F86247"), col.max=0.1)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("NCI","MCI"), signaling=signaling_final, n_sample=21, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"), group.color=c("skyblue", "#F86247"), col.max=8)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("MCI","AD"), signaling=signaling_final, n_sample=21, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"), group.color=c("#F86247", "red3"), col.max=0.1)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("MCI","AD"), signaling=signaling_final, n_sample=21, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"), group.color=c("#F86247", "red3"), col.max=8)

Compare_interaction_multibar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("NCI", "MCI","AD"), cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"), group.color=c("skyblue","#F86247", "red3"))
Compare_interaction_multibar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("NCI", "MCI","AD"), cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"), group.color=c("skyblue","#F86247", "red3"))


###stackedbar sig paths

net_sig <- rbind(net_cogdx_ADvNCI_sig, net_cogdx_MCIvNCI_sig, net_cogdx_ADvMCI_sig)
net_sig <- net_sig[order(net_sig$source, decreasing=T),]
view(net_sig)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_count.pdf", sep=""), width=3, height=6)
remove(plts)
Count_ADNC_sig <- rankNet(Ex_ChP_merge_CC, mode = "comparison", measure = "count", do.stat = TRUE, return.data=T)$signaling.contribution

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_weight.pdf", sep=""), width=3, height=6)
remove(plts)
Weight_ADNC_sig <- rankNet(Ex_ChP_merge_CC, mode = "comparison", measure = "weight", do.stat = TRUE, return.data=T)$signaling.contribution

#ADvNCI
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvNCI_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvNCI_weight.pdf", sep=""), width=3, height=6)
remove(plts)

#ADvMCI
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(2,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvMCI_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(2,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvMCI_weight.pdf", sep=""), width=3, height=6)
remove(plts)

#MCIvNCI
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_MCIvNCI_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_MCIvNCI_weight.pdf", sep=""), width=3, height=6)
remove(plts)


###Heatmap of sig paths

#preprocess
mat1 <- prep_mat(object=Ex_ChPNCI_CC, slot.name="netP", thresh=0.05, pattern="all")
mat2 <- prep_mat(object=Ex_ChPMCI_CC, slot.name="netP", thresh=0.05, pattern="all")
mat3 <- prep_mat(object=Ex_ChPAD_CC, slot.name="netP", thresh=0.05, pattern="all")

mat <- rbind.fill(mat1, mat2, mat3)

netAnalysis_signalingRole_heatmap_combofun(mat, pattern = "all", filename = "NCI-MCI-AD", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "#C48F45", "#3249A6", "#926392", "#261132", "#26532B","#C48F45", "#3249A6", "#926392", "#261132", "#26532B"))


###Communication bubble plot of sig paths

pairLR.use <- rbind(net_sig[, "interaction_name", drop = F])
plts <- netVisual_bubble(Ex_ChP_merge_CC, pairLR.use = pairLR.use, sources.use = net_sig$source, comparison = c(1, 2, 3),  angle.x = 90, remove.isolate = T, color.heatmap = c("white", "black"), title.name = paste("Differential signaling in ", "AD")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Dot_LRpair.pdf", sep=""), width=3+length(unique(net_sig$source))*30*.2, height=3+length(unique(pairLR.use$interaction_name))*.2)


###Chord plot

net_path <- unique(c(net$pathway_name))
chords1 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_NCI_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_aggregate(Ex_ChPNCI_CC, signaling = i, layout = "chord", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"))
  dev.off()
  netVisual_aggregate(Ex_ChPNCI_CC, signaling = i, layout = "chord", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"))
}
chords2 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_MCI_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_aggregate(Ex_ChPMCI_CC, signaling = i, layout = "chord", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"))
  dev.off()
  netVisual_aggregate(Ex_ChPMCI_CC, signaling = i, layout = "chord", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"))
}
chords3 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_AD_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_aggregate(Ex_ChPAD_CC, signaling = i, layout = "chord", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"))
  dev.off()
  netVisual_aggregate(Ex_ChPAD_CC, signaling = i, layout = "chord", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"))
}


for (i in net_path) {
  skip<-FALSE
  tryCatch(chords1(i), error=function(e){skip<-TRUE})
  if(skip){next}
  tryCatch(chords2(i), error=function(e){skip<-TRUE})
  if(skip){next}
  tryCatch(chords3(i), error=function(e){skip<-TRUE})
  if(skip){next}
}
dev.off()


###Circle plots

#NCI
Ex_ChPNCI_CC <- aggregateNet(Ex_ChPNCI_CC,  signaling=unique(net_sig$pathway_name))
pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Weight_NCI.pdf"))
netVisual_circle(Ex_ChPNCI_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
netVisual_circle(Ex_ChPNCI_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Count_NCI.pdf"))
netVisual_circle(Ex_ChPNCI_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
netVisual_circle(Ex_ChPNCI_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#MCI
Ex_ChPMCI_CC <- aggregateNet(Ex_ChPNCI_CC,  signaling=unique(net_sig$pathway_name))
pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Weight_MCI.pdf"))
netVisual_circle(Ex_ChPMCI_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
netVisual_circle(Ex_ChPMCI_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Count_MCI.pdf"))
netVisual_circle(Ex_ChPMCI_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
netVisual_circle(Ex_ChPMCI_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#AD
Ex_ChPNCI_CC <- aggregateNet(Ex_ChPAD_CC,  signaling=unique(net_sig$pathway_name))
pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Weight_AD.pdf"))
netVisual_circle(Ex_ChPAD_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
netVisual_circle(Ex_ChPAD_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Count_AD.pdf"))
netVisual_circle(Ex_ChPAD_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
netVisual_circle(Ex_ChPAD_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



#Cellchat on sub
##RunCC
Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- RenameIdents(Ex_ChP, "Epi_1"="Epithelial", "Epi_2"="Epithelial", "Epi_3"="Epithelial", "Epi_4"="Epithelial", "Endo"="Endothelial", "Mural"="Mural", "BAM_1"="Immune", "BAM_2"="Immune", "BAM_3"="Immune", "T_Cell"="Immune")
levels(Ex_ChP) <- c("Epithelial", "Fib_1", "Fib_2", "Fib_3", "Endothelial", "Mural", "Immune")
Idents(Ex_ChP) -> Ex_ChP@meta.data$Subcelltype
gc()
CCfun(Ex_ChP, "Subcelltype")
gc()


##RunCC per
###Preprocess
Ex_ChP@meta.data$cogdx[which(Ex_ChP@meta.data$cogdx == "5")] <- "4"
Ex_ChP@meta.data$cogdx[which(Ex_ChP@meta.data$cogdx == "3")] <- "2"
Ex_ChPNCI <- subset(Ex_ChP, subset=(cogdx=="1"))
Ex_ChPMCI <- subset(Ex_ChP, subset=(cogdx=="2" | cogdx=="3"))
Ex_ChPAD <- subset(Ex_ChP, subset=(cogdx=="4" | cogdx=="5"))
rm(Ex_ChP)
gc()

###RunCC
CCfun(Ex_ChPNCI, "Subcelltype")
rm(Ex_ChPNCI)
gc()
CCfun(Ex_ChPMCI, "Subcelltype")
rm(Ex_ChPMCI)
gc()
CCfun(Ex_ChPAD, "Subcelltype")
rm(Ex_ChPAD)
gc()


###Venns
####Interaction

interact_NCI <- unique(paste(Ex_ChPNCI_df.net$interaction_name))
interact_MCI <- unique(paste(Ex_ChPMCI_df.net$interaction_name))
interact_AD <- unique(paste(Ex_ChPAD_df.net$interaction_name))
plts <- ggVennDiagram(list(interact_NCI, interact_MCI, interact_AD), category.names = c("NCI", "MCI", "AD")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_Interaction_Comparison.png", sep=""), width=4, height=2.5)

NCI_MCI_AD <- Reduce(intersect, list(interact_NCI, interact_MCI, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_MCI_AD)]
A_MCI <- interact_MCI[!(interact_MCI %in% NCI_MCI_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_MCI_AD)]
Intersects <- list(NCI_MCI_AD, A_NCI, A_MCI, A_AD)
names(Intersects) <- c("NCI_MCI_AD", "NCI", "MCI", "AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_Interaction_NCI_AD.txt", sep=""))


####Pathway

interact_NCI <- unique(paste(Ex_ChPNCI_df.net$pathway_name))
interact_MCI <- unique(paste(Ex_ChPMCI_df.net$pathway_name))
interact_AD <- unique(paste(Ex_ChPAD_df.net$pathway_name))
plts <- ggVennDiagram(list(interact_NCI, interact_MCI, interact_AD), category.names = c("NCI", "MCI", "AD")) + theme(legend.position = "none", text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_pathway_Comparison.png", sep=""), width=4, height=2.5)

NCI_MCI_AD <- Reduce(intersect, list(interact_NCI, interact_MCI, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_MCI_AD)]
A_MCI <- interact_MCI[!(interact_MCI %in% NCI_MCI_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_MCI_AD)]
Intersects <- list(NCI_MCI_AD, A_NCI, A_MCI, A_AD)
names(Intersects) <- c("NCI_MCI_AD", "NCI", "MCI", "AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_pathway_NCI_AD.txt", sep=""))


####ligand-receptor

interact_NCI <- unique(c(paste(Ex_ChPNCI_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
interact_MCI <- unique(c(paste(Ex_ChPMCI_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
interact_AD <- unique(c(paste(Ex_ChPAD_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
plts <- ggVennDiagram(list(interact_NCI, interact_MCI, interact_AD), category.names = c("NCI", "MCI", "AD")) + theme(legend.position = "none", text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_ligand-receptor_Comparison.png", sep=""), width=4, height=2.5)

NCI_MCI_AD <- Reduce(intersect, list(interact_NCI, interact_MCI, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_MCI_AD)]
A_MCI <- interact_MCI[!(interact_MCI %in% NCI_MCI_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_MCI_AD)]
Intersects <- list(NCI_MCI_AD, A_NCI, A_MCI, A_AD)
names(Intersects) <- c("NCI_MCI_AD", "A_NCI", "A_MCI", "A_AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_ligand-receptor_NCI_AD.txt", sep=""))



##Merge, compare and DE
###Merge and compare per condition

#merge results
Ex_ChP_merge_CC <- mergeCellChat(list(NCI=Ex_ChPNCI_CC, MCI=Ex_ChPMCI_CC, AD=Ex_ChPAD_CC), add.names=names(list(NCI=Ex_ChPNCI_CC, MCI=Ex_ChPMCI_CC, AD=Ex_ChPAD_CC)))

#barplots
plts <- compareInteractions(Ex_ChP_merge_CC, show.legend = F, group=c(1,2,3), color.use=c("skyblue", "#F86247", "red3")) + theme(text=element_text(size=7)) #Barplot of counts
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_bar.pdf", sep=""), width=2, height=2.5)

plts <- compareInteractions(Ex_ChP_merge_CC, show.legend = F, group=c(1,2,3), measure="weight", color.use=c("skyblue", "#F86247", "red3")) + theme(text=element_text(size=7)) #Barplot of weights
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_bar.pdf", sep=""), width=2, height=2.5)

#Network
plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,2)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_AD-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,2), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_AD-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(2,3)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_AD-MCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(2,3), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_AD-MCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,3)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_MCI-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,3), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_MCI-NCI.pdf", sep=""))



##Differential gene expression analysis

# perform differential expression analysis
Ex_ChP_merge_CC <- identifyOverExpressedGenes(Ex_ChP_merge_CC, group.dataset = "datasets", pos.dataset = "AD", features.name = "AD", only.pos = FALSE) #!!!!pos.dataset set AD as default dataset for positive fold change 
net <- netMappingDEG(Ex_ChP_merge_CC, features.name = "AD") #select DE

#select up and down regulated pathways with a logFC of .1
net.up <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "AD", ligand.logFC = 0.1, ligand.pvalues=0.01, receptor.pvalues=0.01) #!!!!up in "ad"
gene.up <- extractGeneSubsetFromPair(net.up, Ex_ChP_merge_CC) #get genes
net.down <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "NCI", ligand.logFC = -0.1, ligand.pvalues=0.05, receptor.pvalues=0.05) #!!!!up in "NCI" = down in AD
gene.down <- extractGeneSubsetFromPair(net.down, Ex_ChP_merge_CC)
net.mid <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "MCI", ligand.logFC = 0.1, ligand.pvalues=0.01, receptor.pvalues=0.01) #!!!!up in "MCI"
gene.mid <- extractGeneSubsetFromPair(net.mid, Ex_ChP_merge_CC) #get genes

net_sig <- rbind(net.up, net.mid, net.down)
write.csv(net_sig, paste0(curr_dir, "/CellChat/Ex_ChP_net_sig.csv"))
saveRDS(Ex_ChP_merge_CC, paste(curr_dir, "/CellChat/Ex_ChP_merge_CC_ADNCI.RDS", sep=""))




##DEG with limma instead
###function and prep CC

#get net
net <- subsetCommunication(Ex_ChP_merge_CC)


net <- rbind(net$NCI, net$MCI, net$AD)
net <- net[,c("source", "target", "ligand", "receptor", "prob", "pval", "interaction_name", "pathway_name", "annotation", "evidence")]

#get genes
net <- bind_cols(net, data.frame(t(sapply(str_split(net$receptor,"_",),'[',1:2))))
colnames(net) <- c("source", "target", "ligand", "receptor", "prob", "pval", "interaction_name", "pathway_name", "annotation", "evidence", "receptorA", "receptorB")


genes <- c(net$ligand, net$receptorA, net$receptorB)
genes <- unique(genes)
genes <- genes[!is.na(genes)]

cutoff_pval <- 0.05
cutoff_logFC <- 0.2


###NCIvAD

#make limma_res
limma_cogdx_ADvNCI_Epithelial <- cbind(Celltype="Epithelial", limma_cogdx_ADvNCI_Epithelial)
limma_cogdx_ADvNCI_Fib_1 <- cbind(Celltype="Fib_1", limma_cogdx_ADvNCI_Fib_1)
limma_cogdx_ADvNCI_Fib_2 <- cbind(Celltype="Fib_2", limma_cogdx_ADvNCI_Fib_2)
limma_cogdx_ADvNCI_Fib_3 <- cbind(Celltype="Fib_3", limma_cogdx_ADvNCI_Fib_3)
limma_cogdx_ADvNCI_Endothelial <- cbind(Celltype="Endothelial", limma_cogdx_ADvNCI_Endo)
limma_cogdx_ADvNCI_Mural <- cbind(Celltype="Mural", limma_cogdx_ADvNCI_Mural)
limma_cogdx_ADvNCI_Immune <- cbind(Celltype="Immune", limma_cogdx_ADvNCI_Immune)

limma_res <- bind_rows(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_ADvNCI_Fib_1, limma_cogdx_ADvNCI_Fib_2, limma_cogdx_ADvNCI_Fib_3, limma_cogdx_ADvNCI_Endothelial, limma_cogdx_ADvNCI_Mural, limma_cogdx_ADvNCI_Immune)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdx_ADvNCI <- DEG_source_fun(limma_res, net)
net_cogdx_ADvNCI_sig <- net_cogdx_ADvNCI[net_cogdx_ADvNCI$P.Value_Ligand<cutoff_pval & abs(net_cogdx_ADvNCI$logFC_Ligand)>cutoff_logFC,]  #& c(net_cogdx_ADvMCI$P.Value_ReceptorA<cutoff_pval | net_cogdx_ADvMCI$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdx_ADvMCI$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdx_ADvMCI$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdx_ADvNCI_sig <- net_cogdx_ADvNCI_sig[complete.cases(net_cogdx_ADvNCI_sig$source),]
view(net_cogdx_ADvNCI_sig)



###NCIvMCI

#make limma_res
limma_cogdx_MCIvNCI_Epithelial <- cbind(Celltype="Epithelial", limma_cogdx_MCIvNCI_Epithelial)
limma_cogdx_MCIvNCI_Fib_1 <- cbind(Celltype="Fib_1", limma_cogdx_MCIvNCI_Fib_1)
limma_cogdx_MCIvNCI_Fib_2 <- cbind(Celltype="Fib_2", limma_cogdx_MCIvNCI_Fib_2)
limma_cogdx_MCIvNCI_Fib_3 <- cbind(Celltype="Fib_3", limma_cogdx_MCIvNCI_Fib_3)
limma_cogdx_MCIvNCI_Endothelial <- cbind(Celltype="Endothelial", limma_cogdx_MCIvNCI_Endo)
limma_cogdx_MCIvNCI_Mural <- cbind(Celltype="Mural", limma_cogdx_MCIvNCI_Mural)
limma_cogdx_MCIvNCI_Immune <- cbind(Celltype="Immune", limma_cogdx_MCIvNCI_Immune)

limma_res <- bind_rows(limma_cogdx_MCIvNCI_Epithelial, limma_cogdx_MCIvNCI_Fib_1, limma_cogdx_MCIvNCI_Fib_2, limma_cogdx_MCIvNCI_Fib_3, limma_cogdx_MCIvNCI_Endothelial, limma_cogdx_MCIvNCI_Mural, limma_cogdx_MCIvNCI_Immune)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdx_MCIvNCI <- DEG_source_fun(limma_res, net)
net_cogdx_MCIvNCI_sig <- net_cogdx_MCIvNCI[net_cogdx_MCIvNCI$P.Value_Ligand<cutoff_pval & abs(net_cogdx_MCIvNCI$logFC_Ligand)>cutoff_logFC,] #& c(net_cogdx_ADvMCI$P.Value_ReceptorA<cutoff_pval | net_cogdx_ADvMCI$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdx_ADvMCI$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdx_ADvMCI$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdx_MCIvNCI_sig <- net_cogdx_MCIvNCI_sig[complete.cases(net_cogdx_MCIvNCI_sig$source),]
view(net_cogdx_MCIvNCI_sig)



###MCIvAD

#make limma_res
limma_cogdx_ADvMCI_Epithelial <- cbind(Celltype="Epithelial", limma_cogdx_ADvMCI_Epithelial)
limma_cogdx_ADvMCI_Fib_1 <- cbind(Celltype="Fib_1", limma_cogdx_ADvMCI_Fib_1)
limma_cogdx_ADvMCI_Fib_2 <- cbind(Celltype="Fib_2", limma_cogdx_ADvMCI_Fib_2)
limma_cogdx_ADvMCI_Fib_3 <- cbind(Celltype="Fib_3", limma_cogdx_ADvMCI_Fib_3)
limma_cogdx_ADvMCI_Endothelial <- cbind(Celltype="Endothelial", limma_cogdx_ADvMCI_Endo)
limma_cogdx_ADvMCI_Mural <- cbind(Celltype="Mural", limma_cogdx_ADvMCI_Mural)
limma_cogdx_ADvMCI_Immune <- cbind(Celltype="Immune", limma_cogdx_ADvMCI_Immune)

limma_res <- bind_rows(limma_cogdx_ADvMCI_Epithelial, limma_cogdx_ADvMCI_Fib_1, limma_cogdx_ADvMCI_Fib_2, limma_cogdx_ADvMCI_Fib_3, limma_cogdx_ADvMCI_Endothelial, limma_cogdx_ADvMCI_Mural, limma_cogdx_ADvMCI_Immune)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdx_ADvMCI <- DEG_source_fun(limma_res, net)
net_cogdx_ADvMCI_sig <- net_cogdx_ADvMCI[net_cogdx_ADvMCI$P.Value_Ligand<cutoff_pval & abs(net_cogdx_ADvMCI$logFC_Ligand)>cutoff_logFC ,] #& c(net_cogdx_ADvMCI$P.Value_ReceptorA<cutoff_pval | net_cogdx_ADvMCI$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdx_ADvMCI$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdx_ADvMCI$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdx_ADvMCI_sig <- net_cogdx_ADvMCI_sig[complete.cases(net_cogdx_ADvMCI_sig$source),]
view(net_cogdx_ADvMCI_sig)



###Aggregate net_sigs

net <- rbind(net_cogdx_ADvNCI, net_cogdx_ADvMCI, net_cogdx_MCIvNCI)
write.csv(net, paste0(curr_dir, "/CellChat/Comparison/net_pathways.csv"))
net_sig <- rbind(net_cogdx_ADvNCI_sig, net_cogdx_ADvMCI_sig, net_cogdx_MCIvNCI_sig)
write.csv(net_sig, paste0(curr_dir, "/CellChat/Comparison/net_sig_pathways.csv"))


###Aggregate by pathway

ligandlogFC <- aggregate(net_cogdx_ADvNCI$logFC_Ligand, by=list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdx_ADvNCI$logFC_ReceptorA, list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdx_ADvNCI$logFC_ReceptorB, list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_ReceptorB), mean, all=T)
combinedlogFC_ADvNCI <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_ADvNCI <- cbind(combinedlogFC_ADvNCI[,1:2], rowMeans(combinedlogFC_ADvNCI[3:5], na.rm = TRUE))
view(combinedlogFC_ADvNCI)

ligandlogFC <- aggregate(net_cogdx_ADvMCI$logFC_Ligand, by=list(net_cogdx_ADvMCI$pathway_name, net_cogdx_ADvMCI$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdx_ADvMCI$logFC_ReceptorA, list(net_cogdx_ADvMCI$pathway_name, net_cogdx_ADvMCI$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdx_ADvMCI$logFC_ReceptorB, list(net_cogdx_ADvMCI$pathway_name, net_cogdx_ADvMCI$Celltype_ReceptorB), mean, all=T)
combinedlogFC_ADvMCI <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_ADvMCI <- cbind(combinedlogFC_ADvMCI[,1:2], rowMeans(combinedlogFC_ADvMCI[3:5], na.rm = TRUE))
view(combinedlogFC_ADvMCI)

ligandlogFC <- aggregate(net_cogdx_MCIvNCI$logFC_Ligand, by=list(net_cogdx_MCIvNCI$pathway_name, net_cogdx_MCIvNCI$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdx_MCIvNCI$logFC_ReceptorA, list(net_cogdx_MCIvNCI$pathway_name, net_cogdx_MCIvNCI$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdx_MCIvNCI$logFC_ReceptorB, list(net_cogdx_MCIvNCI$pathway_name, net_cogdx_MCIvNCI$Celltype_ReceptorB), mean, all=T)
combinedlogFC_MCIvNCI <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_MCIvNCI <- cbind(combinedlogFC_MCIvNCI[,1:2], rowMeans(combinedlogFC_MCIvNCI[3:5], na.rm = TRUE))
view(combinedlogFC_MCIvNCI)



signaling_NCI <- dimnames(Ex_ChP_merge_CC@net$NCI$prob)[[3]]
signaling_MCI <- dimnames(Ex_ChP_merge_CC@net$MCI$prob)[[3]]
signaling_AD <- dimnames(Ex_ChP_merge_CC@net$AD$prob)[[3]]
signaling_shared <- intersect(signaling_AD, intersect(signaling_MCI, signaling_NCI))
signaling_sig <- unique(net_sig$interaction_name) # 284 sig pathways generated using CC_limma.R script
signaling_final <- intersect(signaling_sig, signaling_shared)


Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("NCI","AD"), signaling=signaling_final, n_sample=21, cell.color=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"), group.color=c("skyblue", "red3"), col.max=0.3)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("NCI","AD"), signaling=signaling_final, n_sample=21, cell.color=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"), group.color=c("skyblue", "red3"), col.max=10)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("NCI","MCI"), signaling=signaling_final, n_sample=21, cell.color=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"), group.color=c("skyblue", "#F86247"), col.max=0.3)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("NCI","MCI"), signaling=signaling_final, n_sample=21, cell.color=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"), group.color=c("skyblue", "#F86247"), col.max=10)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("MCI","AD"), signaling=signaling_final, n_sample=21, cell.color=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"), group.color=c("#F86247", "red3"), col.max=0.3)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("MCI","AD"), signaling=signaling_final, n_sample=21, cell.color=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"), group.color=c("#F86247", "red3"), col.max=10)

Compare_interaction_multibar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("NCI", "MCI","AD"), signaling=signaling_final, cell.color=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"), group.color=c("skyblue","#F86247", "red3"))
Compare_interaction_multibar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("NCI", "MCI","AD"), signaling=signaling_final, cell.color=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"), group.color=c("skyblue","#F86247", "red3"))


###stackedbar sig paths

net_sig <- net_sig[order(net_sig$source, decreasing=T),]
view(net_sig)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_weight.pdf", sep=""), width=3, height=6)
remove(plts)

#ADvNCI
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvNCI_count.pdf", sep=""), width=3, height=6)
remove(plts)
Weight_ADNC_sig <- rankNet(Ex_ChP_merge_CC, mode = "comparison", measure = "weight", do.stat = TRUE, return.data=T)$signaling.contribution

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvNCI_weight.pdf", sep=""), width=3, height=6)
remove(plts)

#ADvMCI
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(2,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvMCI_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(2,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvMCI_weight.pdf", sep=""), width=3, height=6)
remove(plts)

#MCIvNCI
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_MCIvNCI_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_MCIvNCI_weight.pdf", sep=""), width=3, height=6)
remove(plts)

###Heatmap of sig paths

#preprocess
mat1 <- prep_mat(object=Ex_ChPNCI_CC, slot.name="netP", thresh=0.05, pattern="all")
mat2 <- prep_mat(object=Ex_ChPMCI_CC, slot.name="netP", thresh=0.05, pattern="all")
mat3 <- prep_mat(object=Ex_ChPAD_CC, slot.name="netP", thresh=0.05, pattern="all")

mat <- rbind.fill(mat1, mat2, mat3)

netAnalysis_signalingRole_heatmap_combofun(mat, pattern = "all", signaling = unique(net$pathway_name), filename = "NCI-MCI-AD", color.use=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97", "#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97", "#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"))


###Communication bubble plot of sig paths

pairLR.use <- rbind(net_sig[, "interaction_name", drop = F])
plts <- netVisual_bubble(Ex_ChP_merge_CC, pairLR.use = pairLR.use, sources.use = net_sig$source, comparison = c(1, 2, 3),  angle.x = 90, remove.isolate = T, color.heatmap = c("white", "black"), title.name = paste("Differential signaling in ", "AD")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Dot_LRpair.pdf", sep=""), width=3+length(unique(net_sig$source))*30*.2, height=3+length(unique(pairLR.use$interaction_name))*.2)


###Chord plot

net_path <- unique(c(net_sig$pathway_name))
chords1 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_NCI_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_aggregate(Ex_ChPNCI_CC, signaling = i, layout = "chord", color.use=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"))
  dev.off()
  netVisual_aggregate(Ex_ChPNCI_CC, signaling = i, layout = "chord", color.use=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"))
}
chords2 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_MCI_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_aggregate(Ex_ChPMCI_CC, signaling = i, layout = "chord", color.use=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"))
  dev.off()
  netVisual_aggregate(Ex_ChPMCI_CC, signaling = i, layout = "chord", color.use=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"))
}
chords3 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_AD_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_aggregate(Ex_ChPAD_CC, signaling = i, layout = "chord", color.use=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"))
  dev.off()
  netVisual_aggregate(Ex_ChPAD_CC, signaling = i, layout = "chord", color.use=c("#FFA500", "#CD8500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#8f2aa2", "#6d207b", "#70C25B", "#0E9554", "#266417", "#3ded97"))
}


for (i in net_path) {
  skip<-FALSE
  tryCatch(chords1(i), error=function(e){skip<-TRUE})
  if(skip){next}
  tryCatch(chords2(i), error=function(e){skip<-TRUE})
  if(skip){next}
  tryCatch(chords3(i), error=function(e){skip<-TRUE})
  if(skip){next}
}
dev.off()



Idents(BAM) <- BAM$Subcelltype
DimPlot(BAM)



for(i in c("IGKC", "IGHG1.2", "JCHAIN", "IGLL5", "PPP2R5E", "ARSA", "LRR1", "JAG2", "SOX4", "PIM2", "KCTD7", "IGHA1", "FGF18", "MBOAT7", "ERN1", "HMGA2", "KCNJ6", "IGKV4.1", "SUPT20H", "LDB1", "IGHV1.69", "CD27", "TBK1", "PLOD2", "BATF", "SERINC3", "REPS1", "SAE1", "VPS4B", "ACSM3", "PDGFA", "TMTC4", "XPO1", "MARK3", "ID1", "LRRC3B", "PET117", "AHI1", "BTLA", "SIGIRR", "DLGAP2", "MAP2K3", "BIRC6", "PAG1", "VPS37A", "IGLL1", "MAPK9", "SETD5", "RGS6")){
  if(i %in% rownames(Ex_ChP)){
    plot(FeaturePlot(Imm, i))
    plot(FeaturePlot(BAM, i))
    plot(FeaturePlot(Ex_ChP, i, raster=F))
  }
}
