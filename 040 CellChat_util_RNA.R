###
### 040 CellChat_util_RNA.R
###
# Purpose: Cellchat analysis functions
# Dependencies:
library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(circlize)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(ggVennDiagram)
library(export)
library(tidyr)
options(stringsAsFactors = FALSE)
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
unique(CellChatDB[["interaction"]][["annotation"]])

future::plan("multisession", workers = 2) #!!!!CORES
options(future.globals.maxSize = 200* 1024^3) #!!!!RAM
future.seed=TRUE



CCfun <- function(srt, group.by){
  #Description: Preprocessing, network interactions and some visualizations
  #input:
  ##srt: Seurat object
  ##group.by: to createCellChat object default is group.by=idents, but can specify celltypes
  #output: Cellchat object
  nameobj <- deparse(substitute(srt)) #make char
  path <- paste(curr_dir, "/CellChat/", nameobj, "/", sep="") #specify dir
  if(dir.exists(path)==FALSE){dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  #Make CC object from Seurat object
  CC_obj <- createCellChat(srt, group.by=group.by) 
  CC_obj@DB <- CellChatDB #Add database to use to object
  CC_obj <- subsetData(CC_obj) #Add data - database intersections
  
  #Project gene expression data onto database
  CC_obj <- identifyOverExpressedGenes(CC_obj)
  CC_obj <- identifyOverExpressedInteractions(CC_obj)
  
  #Infer network
  CC_obj <- computeCommunProb(CC_obj)
  CC_obj <- computeCommunProbPathway(CC_obj, thresh=0.05)
  CC_obj <- filterCommunication(CC_obj, min.cells = 10) #Filter out cell-cell interactions with few cells
  df.net <- subsetCommunication(CC_obj) #Check which cell-cell interaction pathways are significant
  head(df.net)
  
  #Aggregate cell-cell networks
  CC_obj <- aggregateNet(CC_obj)
  groupSize <- as.numeric(table(CC_obj@idents))
  
  
  #Visualize networks
  pdf(paste(path, "Count_network.pdf", sep=""))
  netVisual_circle(CC_obj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  
  pdf(paste(path, "Weight_network.pdf", sep=""))
  netVisual_circle(CC_obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(CC_obj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(CC_obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  
  #Visualize individual cell types
  mat <- CC_obj@net$weight
  for (i in rownames(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    pdf(paste(path, "circle_", i, ".pdf", sep=""), width=480, height=480)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  
  #Chord plot
  path <- unique(df.net$pathway_name)
  chords <- function(i){
    png(paste(curr_dir, "/CellChat/", nameobj, "/", path, "chord_", i, ".png", sep=""), width=480, height=480, type="cairo")
    netVisual_aggregate(CC_obj, signaling = i, layout = "chord")
    dev.off()
    netVisual_aggregate(CC_obj, signaling = i, layout = "chord")
  }
  
  for (i in path) {
    skip<-FALSE
    tryCatch(chords(i), error=function(e){skip<-TRUE})
    if(skip){next}
  }
  
  assign(paste(nameobj, "_CC", sep=""), CC_obj, envir = .GlobalEnv)
  assign(paste(nameobj, "_df.net", sep=""), df.net, envir = .GlobalEnv)
  saveRDS(CC_obj, file = paste(curr_dir, "/CellChat/", nameobj, "_CC.rds", sep=""))
}


#Construct stats table that matches changes in gene expression (more stringent)

DEG_source_fun <- function(limma_res, net){
  #Description: Merges CC net results to limma results table to be more stringent
  #input:
  ##limma_res: limma result
  ##net: output from subsetCommunications
  #output: Data table summarizing net from cellchat and matching limma results
  
  
  #Add source limma results
  net_sig <- data.frame(t(data.frame(sapply(
    unique(net$source), function(j){
      net_S <- net[net$source==j,]
      limma_S <- limma_res[limma_res$Celltype==j,]
      colnames(limma_S) <- paste0(colnames(limma_S), "_Ligand")
      sapply(1:nrow(net_S), function(i) bind_cols(net_S[i,], limma_S[limma_S$gene==net_S[i,]$ligand,]))
    }
  ))))
  rownames(net_sig) <- NULL
  net_sig <- unnest(net_sig, cols = colnames(net_sig))
  
  #Add receptorA limma results
  net_sig <- data.frame(t(data.frame(sapply(
    unique(net_sig$target), function(j){
      net_S <- net_sig[net_sig$target==j,]
      limma_S <- limma_res[limma_res$Celltype==j,]
      colnames(limma_S) <- paste0(colnames(limma_S), "_ReceptorA")
      sapply(1:nrow(net_S), function(i) bind_cols(net_S[i,], limma_S[limma_S$gene==net_S[i,]$receptorA,]))
    }
  ))))
  rownames(net_sig) <- NULL
  net_sig <- unnest(net_sig, cols = colnames(net_sig))
  
  #Add receptorB limma results
  net_sigB <- data.frame(t(data.frame(sapply(
    unique(net_sig$target), function(j){
      net_S <- net_sig[net_sig$target==j,]
      limma_S <- limma_res[limma_res$Celltype==j,]
      colnames(limma_S) <- paste0(colnames(limma_S), "_ReceptorB")
      sapply(1:nrow(net_S), function(i) bind_cols(net_S[i,], limma_S[limma_S$gene==net_S[i,]$receptorB,]))
    }
  ))))
  rownames(net_sigB) <- NULL
  net_sigB <- unnest(net_sigB, cols = colnames(net_sigB))
  net_sigB  <- net_sigB[complete.cases(net_sigB), ]
  
  net_sig  <- net_sig[!complete.cases(net_sig), ]
  net_sig <- bind_rows(net_sig, net_sigB)
  return(net_sig)
}



#Make custom heatmap and bar plots of interaction strengths with sig stars.

Compare_interaction_heatbar_fun <- function(object, measure=c("count", "weight"), comparison = c(1, 2), slot.name = "net", signaling = NULL, n_sample=2, cell.color = NULL, group.color=NULL, color.heatmap = c("blue", "white", "red"), col.max=NULL, norm_factor=NULL){
  measure <- match.arg(measure)
  object.names <- names(methods::slot(object, slot.name))
  
  #extract probe data
  prob.list <- list()
  for (i in 1:length(comparison)) {
    object.list <- methods::slot(object, slot.name)[[comparison[[i]]]]
    prob <- object.list$prob
    if (!is.null(norm_factor)){
      norm_fac <- norm_factor[[i]]
      prob <- as.data.frame(as.table(prob)) %>%
        rename(Var1="Sender", Var2="Receiver", Freq="Comm_Prob") %>%
        left_join(norm_fac, by = c("Sender", "Receiver"))
      prob <- prob %>% mutate(Norm_Comm_Prob = Comm_Prob / Geo_Mean)
      prob <- prob[,c("Sender", "Receiver", "Var3", "Norm_Comm_Prob")]
      prob <- data.frame(prob)
      prob <- with(prob, array(data=prob$Norm_Comm_Prob, dim=c(length(unique(prob$Sender)),
                                                               length(unique(prob$Receiver)),
                                                               length(unique(prob$Var3))),
                               dimnames=list(unique(prob$Sender), unique(prob$Receiver), unique(prob$Var3))))
    }
    prob[object.list$pval > 0.05] <- 0
    if(!is.null(signaling)){prob <- prob[,,signaling]}
    if (measure == "count") {
      prob <- 1 * (prob > 0)
    }
    prob.list[[comparison[[i]]]] <- prob
  }
  
  cell.name <- c(rownames(prob.list[[1]]), rownames(prob.list[[2]]))
  cell.name.all <- as.character(unique(unlist(cell.name)))
  
  #extract and stat test interaction data per celltype
  #Sending
  df.send <- data.frame(c(cell.name.all)) #prime df to take in pval
  colnames(df.send) <- c("Celltype")
  
  for (i in 1:length(cell.name.all)) {
    prob1 <- as.vector(prob.list[[1]][cell.name.all[i], ,]) #note dim 1 (Sender) is used here
    prob2 <- as.vector(prob.list[[2]][cell.name.all[i], ,]) #note dim 1 (Sender) is used here
    
    #assign to df
    df.send[[comparison[[1]]]][df.send$Celltype == cell.name.all[i]] <- sum(prob1)
    df.send[[comparison[[2]]]][df.send$Celltype == cell.name.all[i]] <- sum(prob2)
    
    #make both same length
    maxlength <- max(length(prob1), length(prob2))
    prob1 <- c(prob1, rep(NA, maxlength - length(prob1)))
    prob2 <- c(prob2, rep(NA, maxlength - length(prob2)))
    prob.values <- cbind(prob1, prob2)
    prob.values[is.na(prob.values)] <- 0
    
    prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE] #remove pathways negative in both
  }
  
  #ttest
  for (i in 1:nrow(df.send)){
    q <- (df.send[i,2]-df.send[i,3])/sd(pivot_longer(df.send[1:3], cols=c(comparison[[1]], comparison[[2]]))$value)
    p_value <- pt(q, n_sample-1, lower.tail = F)
    df.send$pvalues[i] <- p_value
  }
  
  write.csv(df.send, paste(curr_dir, "/CellChat/Comparison/Interaction_bar_heatmap_", measure, "_Send_STATs", comparison[[1]], "-", comparison[[2]], ".csv", sep=""))
  print("Note: You need to manually add pvalues to the bar graph, the STATs file was saved.")
  df.send$pvalues <- NULL
  df.send <- pivot_longer(df.send, cols=2:3, names_to="Group", values_to="Values")
  df.send$Celltype <- factor(df.send$Celltype, levels = unique(df.send$Celltype))
  df.send$Group <- factor(df.send$Group, levels = unique(df.send$Group))
  plts <- ggplot(df.send, aes(x=Celltype, y=Values, fill=Celltype, color=Group)) +
    geom_bar(stat="identity", position=position_dodge(width=.9, preserve="single"), width=.8, linewidth=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7)) + scale_fill_manual(values=cell.color) + scale_color_manual(values=group.color) +
    labs(y="Interaction Strength") 
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Interaction_bar_heatmap_", measure, "_Send", comparison[[1]], "-", comparison[[2]], ".pdf", sep=""), width=4, height=2)
  
  #Receiving
  df.receive <- data.frame(c(cell.name.all)) #prime df to take in pval
  colnames(df.receive) <- c("Celltype")
  
  for (i in 1:length(cell.name.all)) {
    prob1 <- as.vector(prob.list[[1]][, cell.name.all[i],]) #note dim 2 (Receiver) is used here
    prob2 <- as.vector(prob.list[[2]][, cell.name.all[i],]) #note dim 2 (Receiver) is used here
    
    #assign to df
    df.receive[[comparison[[1]]]][df.receive$Celltype == cell.name.all[i]] <- sum(prob1)
    df.receive[[comparison[[2]]]][df.receive$Celltype == cell.name.all[i]] <- sum(prob2)
    
    #make both same length
    maxlength <- max(length(prob1), length(prob2))
    prob1 <- c(prob1, rep(NA, maxlength - length(prob1)))
    prob2 <- c(prob2, rep(NA, maxlength - length(prob2)))
    prob.values <- cbind(prob1, prob2)
    prob.values[is.na(prob.values)] <- 0
    
    prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE] #remove pathways negative in both
  }
  #ttest
  for (i in 1:nrow(df.receive)){
    q <- (df.receive[i,2]-df.receive[i,3])/sd(pivot_longer(df.receive[1:3], cols=c(comparison[[1]], comparison[[2]]))$value)
    p_value <- pt(q, n_sample-1, lower.tail = F)
    df.receive$pvalues[i] <- p_value
  }
  
  write.csv(df.receive, paste(curr_dir, "/CellChat/Comparison/Interaction_bar_heatmap_", measure, "_Receive_STATs", comparison[[1]], "-", comparison[[2]], ".csv", sep=""))
  print("Note: You need to manually add pvalues to the bar graph, the STATs file was saved.")
  df.receive$pvalues <- NULL
  df.receive <- pivot_longer(df.receive, cols=2:3, names_to="Group", values_to="Values")
  df.receive$Celltype <- factor(df.receive$Celltype, levels = unique(df.receive$Celltype))
  df.receive$Group <- factor(df.receive$Group, levels = unique(df.receive$Group))
  plts <- ggplot(df.receive, aes(x=Celltype, y=Values, fill=Celltype, color=Group)) +
    geom_bar(stat="identity", position=position_dodge(width=.9, preserve="single"), width=.8, linewidth=2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7)) + scale_fill_manual(values=cell.color) +
    scale_color_manual(values=group.color) +
    labs(y="Interaction Strength") 
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Interaction_bar_heatmap_", measure, "_Receive", comparison[[1]], "-", comparison[[2]], ".pdf", sep=""), width=4, height=2)
  
  #Both Receive and send
  df.both <- data.frame(matrix(nrow=length(cell.name.all), ncol=length(cell.name.all))) #prime df to take in data and pval
  colnames(df.both) <- cell.name.all
  rownames(df.both) <- cell.name.all
  df.pval <- df.both
  lprob1 <- list()
  lprob2 <- list()
  lpval <- list()
  
  for (i in 1:length(cell.name.all)) {
    for (j in 1:length(cell.name.all)) {
      prob1 <- as.vector(prob.list[[1]][cell.name.all[i], cell.name.all[j],]) #dim1 (sender) permutation with dim2 (receiver)
      prob2 <- as.vector(prob.list[[2]][cell.name.all[i], cell.name.all[j],])
      
      #make both same length
      maxlength <- max(length(prob1), length(prob2))
      prob1 <- c(prob1, rep(NA, maxlength - length(prob1)))
      prob2 <- c(prob2, rep(NA, maxlength - length(prob2)))
      prob.values <- cbind(prob1, prob2)
      prob.values[is.na(prob.values)] <- 0
      
      prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE] #remove pathways negative in both
      
      #assign
      lprob1[[j]] <- as.vector(prob.list[[1]][cell.name.all[i], cell.name.all[j],]) #dim1 (sender) permutation with dim2 (receiver)
      lprob2[[j]] <- as.vector(prob.list[[2]][cell.name.all[i], cell.name.all[j],])
    }
    #assign to df
    df.both[rownames(df.both) == paste(cell.name.all[i], sep="_"),] <- sapply(lprob2, sum)-sapply(lprob1, sum)
    
    df.i <- data.frame(sapply(lprob2, sum),sapply(lprob1, sum))
    colnames(df.i) <- c(comparison)
    for (j in 1:nrow(df.i)){
      q <- (df.i[j,1]-df.i[j,2])/sd(c(df.i[[1]],df.i[[2]]))
      p_value <- pt(q, n_sample-1, lower.tail = F)
      df.pval[i,j] <- p_value
    }
  }
  
  df.pval[is.na(df.pval)] <- 1
  
  
  if (measure == "count") {
    title.name = "Differential number of interactions"
  } else if (measure == "weight") {
    title.name = "Differential interaction strength"
  }
  legend.name = "Relative values"
  
  mat <- as.matrix(df.both)
  #colors
  if (is.null(cell.color)) {
    cell.color <- scPalette(ncol(mat))
  }
  names(cell.color) <- colnames(mat)
  
  #set max color
  if(is.null(col.max)){col.max <- max(abs(mat))}
  
  if (length(color.heatmap) == 4) {
    color.heatmap.use = colorRamp3(c(-col.max, -col.max/2, col.max, col.max/2), color.heatmap)
  }else if (length(color.heatmap) == 3) {
    color.heatmap.use = colorRamp3(c(-col.max, 0, col.max), color.heatmap)
  }else if (length(color.heatmap) == 2) {
    color.heatmap.use = colorRamp3(c(-col.max, col.max), color.heatmap)
  }else if (length(color.heatmap) == 1) {
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                               name = color.heatmap))))(100)
  }
  
  #annotations
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.color), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.color), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  
  ht1 = Heatmap(mat,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(df.pval[i, j] < 0.05) {
                    grid.text("*", x, y, gp = gpar(fontsize = 30))
                  }
                },
                col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation = row_annotation,
                cluster_rows = F, cluster_columns = F, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
                column_title = gt_render(paste0(title.name, "<br> Receiver"), r = unit(2, "pt"), padding = unit(c(2, 2), "pt")),
                column_title_gp = gpar(fontsize = 7), column_names_rot = 90, row_title = "Sender", row_title_gp = gpar(fontsize = 7),
                row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),
                                                                title_position = "leftcenter-rot", border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, "mm")))
  draw(ht1)
  graph2pdf(ht1, paste(curr_dir, "/CellChat/Comparison/Interaction_bar_heatmap_", measure, comparison[[1]], "-", comparison[[2]], ".pdf", sep=""), width=4, height=4)
}

Compare_interaction_multibar_fun <- function(object, measure=c("count", "weight"), comparison = c(1, 2, 3), slot.name = "net", signaling = NULL, cell.color = NULL, group.color=NULL, color.heatmap = c("blue", "white", "red")){
  measure <- match.arg(measure)
  object.names <- names(methods::slot(object, slot.name))
  
  #extract probe data
  prob.list <- list()
  for (i in 1:length(comparison)) {
    object.list <- methods::slot(object, slot.name)[[comparison[[i]]]]
    prob <- object.list$prob
    if (!is.null(norm_factor)){
      norm_fac <- norm_factor[[i]]
      prob <- as.data.frame(as.table(prob)) %>%
        rename(Var1="Sender", Var2="Receiver", Freq="Comm_Prob") %>%
        left_join(norm_fac, by = c("Sender", "Receiver"))
      prob <- prob %>% mutate(Norm_Comm_Prob = Comm_Prob / Geo_Mean)
      prob <- prob[,c("Sender", "Receiver", "Var3", "Norm_Comm_Prob")]
      prob <- data.frame(prob)
      prob <- with(prob, array(data=prob$Norm_Comm_Prob, dim=c(length(unique(prob$Sender)),
                                                               length(unique(prob$Receiver)),
                                                               length(unique(prob$Var3))),
                               dimnames=list(unique(prob$Sender), unique(prob$Receiver), unique(prob$Var3))))
    }
    prob[object.list$pval > 0.05] <- 0
    prob <- prob[,,signaling]
    if (measure == "count") {
      prob <- 1 * (prob > 0)
    }
    prob.list[[comparison[[i]]]] <- prob
  }
  
  cell.name <- c(rownames(prob.list[[1]]), rownames(prob.list[[2]], rownames(prob.list[[3]])))
  cell.name.all <- as.character(unique(unlist(cell.name)))
  
  #extract and stat test interaction data per celltype
  #Sending
  df.send <- data.frame(c(cell.name.all)) #prime df to take in pval
  colnames(df.send) <- c("Celltype")
  for (j in 1:length(comparison)) {
    for (i in 1:length(cell.name.all)) {
      prob <- as.vector(prob.list[[j]][cell.name.all[i], ,]) #note dim 1 (Sender) is used here
      #assign to df
      df.send[[object.names[[j]]]][df.send$Celltype == cell.name.all[i]] <- sum(prob)
    }
  }
  
  df.send <- pivot_longer(df.send, cols=2:ncol(df.send), names_to="Group", values_to="Values")
  df.send$Celltype <- factor(df.send$Celltype, levels = unique(df.send$Celltype))
  df.send$Group <- factor(df.send$Group, levels = unique(df.send$Group))
  plts <- ggplot(df.send, aes(x=Celltype, y=Values, fill=Celltype, color=Group)) +
    geom_bar(stat="identity", position=position_dodge(width=.9, preserve="single"), width=.8, linewidth=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7)) + scale_fill_manual(values=cell.color) + scale_color_manual(values=group.color) +
    labs(y="Interaction Strength") 
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Interaction_bar_", measure, "_Send", comparison[[1]], "-", comparison[[2]], "-", comparison[[3]], ".pdf", sep=""), width=4, height=2)
  
  #Receiving
  df.receive <- data.frame(c(cell.name.all)) #prime df to take in pval
  colnames(df.receive) <- c("Celltype")
  
  for (j in 1:length(comparison)) {
    for (i in 1:length(cell.name.all)) {
      prob <- as.vector(prob.list[[j]][, cell.name.all[i],]) #note dim 2 (Receiver) is used here
      #assign to df
      df.receive[[object.names[[j]]]][df.receive$Celltype == cell.name.all[i]] <- sum(prob)
    }
  }
  
  df.receive$pvalues <- NULL
  df.receive <- pivot_longer(df.receive, cols=2:ncol(df.receive), names_to="Group", values_to="Values")
  df.receive$Celltype <- factor(df.receive$Celltype, levels = unique(df.receive$Celltype))
  df.receive$Group <- factor(df.receive$Group, levels = unique(df.receive$Group))
  plts <- ggplot(df.receive, aes(x=Celltype, y=Values, fill=Celltype, color=Group)) +
    geom_bar(stat="identity", position=position_dodge(width=.9, preserve="single"), width=.8, linewidth=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=7)) + scale_fill_manual(values=cell.color) + scale_color_manual(values=group.color) +
    labs(y="Interaction Strength") 
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Interaction_bar_", measure, "_Receive", comparison[[1]], "-", comparison[[2]], "-", comparison[[3]], ".pdf", sep=""), width=4, height=2)
}


#broken code the cellchat authors will hopefully fix in the next release
netVisual_diffInteraction <- function (object, comparison = c(1, 2), measure = c("count", "weight", "count.merged", "weight.merged"), color.use = NULL, 
                                       color.edge = c("#b2182b", "#2166ac"), title.name = NULL, 
                                       sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
                                       top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
                                       vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
                                       edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
                                       label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
                                       edge.curved = 0.2, shape = "circle", layout = in_circle(), 
                                       margin = 0.2, arrow.width = 1, arrow.size = 0.2) 
{
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top, 
                                 na.rm = T)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, 
                       -atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g), 
                                                                        1]), pi - atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g), 
                                                                                                                                  1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle <- NA
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}



#Allow me to pick all heatmap colors
netVisual_heatmap <- function (object, comparison = c(1, 2), measure = c("count", 
                                                                         "weight"), signaling = NULL, slot.name = c("netP", "net"), 
                               color.use = NULL, color.heatmap = c("blue", "white", "red"), 
                               title.name = NULL, width = NULL, height = NULL, font.size = 8, 
                               font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE, 
                               sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
                               row.show = NULL, col.show = NULL) 
{
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    }
    else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  #Simplified
  
  #making it symmetrical
  if(abs(min(mat))>max(mat)){col.max <- abs(min(mat))}
  else{col.max <- max(mat)}
  if (length(color.heatmap) == 4) {
    color.heatmap.use = colorRamp3(c(-col.max, -col.max/2, col.max, col.max/2), color.heatmap)
  }
  if (length(color.heatmap) == 3) {
    color.heatmap.use = colorRamp3(c(-col.max, 0, col.max), color.heatmap)
  }
  else if (length(color.heatmap) == 2) {
    color.heatmap.use = colorRamp3(c(-col.max, col.max), color.heatmap)
  }
  else if (length(color.heatmap) == 1) {
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                               name = color.heatmap))))(100)
  }
  
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = legend.name, bottom_annotation = col_annotation, 
                left_annotation = row_annotation, top_annotation = ha2, 
                right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), column_title = title.name, 
                column_title_gp = gpar(fontsize = font.size.title), 
                column_names_rot = 90, row_title = "Sources (Sender)", 
                row_title_gp = gpar(fontsize = font.size.title), row_title_rot = 90, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                            fontface = "plain"), title_position = "leftcenter-rot", 
                                            border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), 
                                            grid_width = unit(2, "mm")))
  return(ht1)
}

#transparancy default set to 0
netVisual_chord_cell_internal <- function (net, color.use = NULL, group = NULL, cell.order = NULL, 
                                           sources.use = NULL, targets.use = NULL, lab.cex = 0.8, small.gap = 1, 
                                           big.gap = 10, annotationTrackHeight = c(0.03), remove.isolate = FALSE, 
                                           link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, 
                                           reduce = -1, link.border = NA, title.name = NULL, 
                                           show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,
                                           ...) 
{
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source", "target")
  }
  else if (is.data.frame(net)) {
    if (all(c("source", "target", "prob") %in% colnames(net)) == 
        FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source, net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)
  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  net <- subset(net, prob > 0)
  if (dim(net)[1] <= 0) {
    message("No interaction between those cells")
  }
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source, 
                                                             net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(cells.removed, cells.removed, 
                             1e-10 * sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if (nrow(net) > nrow(net.fake)) {
        link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      }
      scale = TRUE
    }
  }
  df <- net
  cells.use <- union(df$source, df$target)
  order.sector <- cell.levels[cell.levels %in% cells.use]
  if (is.null(color.use)) {
    color.use = scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  }
  else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }
  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector
  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }
  edge.color <- color.use[as.character(df$source)]
  if (directional == 0 | directional == 2) {
    link.arr.type = "triangle"
  }
  else {
    link.arr.type = "big.arrow"
  }
  circos.clear()
  chordDiagram(df, order = order.sector, col = edge.color, 
               grid.col = grid.col, transparency = 0, link.border = link.border, 
               directional = directional, direction.type = c("diffHeight", 
                                                             "arrows"), link.arr.type = link.arr.type, annotationTrack = "grid", 
               annotationTrackHeight = annotationTrackHeight, preAllocateTracks = list(track.height = max(strwidth(order.sector))), 
               small.gap = small.gap, big.gap = big.gap, link.visible = link.visible, 
               scale = scale, group = group, link.target.prop = link.target.prop, 
               reduce = reduce, ...)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                niceFacing = TRUE, adj = c(0, 0.5), cex = lab.cex)
  }, bg.border = NA)
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(grid.col), 
                                  type = "grid", legend_gp = grid::gpar(fill = grid.col), 
                                  title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc") - unit(legend.pos.x, 
                                                        "mm"), y = unit(legend.pos.y, "mm"), just = c("right", 
                                                                                                      "bottom"))
  }
  if (!is.null(title.name)) {
    text(-0, 1.02, title.name, cex = 1)
  }
  circos.clear()
  gg <- recordPlot()
  return(gg)
}

netVisual_aggregate <- function (object, signaling, signaling.name = NULL, color.use = NULL, 
                                 thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, 
                                 targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE, 
                                 vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL, 
                                 weight.scale = TRUE, edge.weight.max = NULL, edge.width.max = 8, 
                                 layout = c("circle", "hierarchy", "chord", "spatial"), pt.title = 12, 
                                 title.space = 6, vertex.label.cex = 0.8, alpha.image = 0.15, 
                                 point.size = 1.5, group = NULL, cell.order = NULL, small.gap = 1, 
                                 big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, 
                                 legend.pos.x = 20, legend.pos.y = 20, ...) 
{
  layout <- match.arg(layout)
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    }
    else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, 
                       key = "pathway_name", matching.exact = T, pair.only = T)
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 
                                         3, sum) != 0]
  }
  else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 
                                     0]
  }
  if (length(pairLR.name.use) == 0) {
    stop(paste0("There is no significant communication of ", 
                signaling.name))
  }
  else {
    pairLR <- pairLR[pairLR.name.use, ]
  }
  nRow <- length(pairLR.name.use)
  prob <- prob[, , pairLR.name.use]
  pval <- pval[, , pairLR.name.use]
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
    pval <- replicate(1, pval, simplify = "array")
  }
  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1, 2), sum)
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow = c(1, 2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, 
                         sources.use = sources.use, targets.use = targets.use, 
                         remove.isolate = remove.isolate, top = top, color.use = color.use, 
                         vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, 
                         vertex.size.max = vertex.size.max, weight.scale = weight.scale, 
                         edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, 
                         title.name = NULL, vertex.label.cex = vertex.label.cex, 
                         ...)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum), 
                                                             vertex.receiver), sources.use = sources.use, targets.use = targets.use, 
                         remove.isolate = remove.isolate, top = top, color.use = color.use, 
                         vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, 
                         vertex.size.max = vertex.size.max, weight.scale = weight.scale, 
                         edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, 
                         title.name = NULL, vertex.label.cex = vertex.label.cex, 
                         ...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), 
                    side = 3, outer = TRUE, cex = 1, line = -title.space)
    gg <- recordPlot()
  }
  else if (layout == "circle") {
    prob.sum <- apply(prob, c(1, 2), sum)
    gg <- netVisual_circle(prob.sum, sources.use = sources.use, 
                           targets.use = targets.use, idents.use = idents.use, 
                           remove.isolate = remove.isolate, top = top, color.use = color.use, 
                           vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, 
                           vertex.size.max = vertex.size.max, weight.scale = weight.scale, 
                           edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, 
                           title.name = paste0(signaling.name, " signaling pathway network"), 
                           vertex.label.cex = vertex.label.cex, ...)
  }
  else if (layout == "spatial") {
    prob.sum <- apply(prob, c(1, 2), sum)
    if (vertex.weight == "incoming") {
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$indeg
    }
    if (vertex.weight == "outgoing") {
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$outdeg
    }
    coordinates <- object@images$coordinates
    labels <- object@idents
    gg <- netVisual_spatial(prob.sum, coordinates = coordinates, 
                            labels = labels, alpha.image = alpha.image, point.size = point.size, 
                            sources.use = sources.use, targets.use = targets.use, 
                            idents.use = idents.use, remove.isolate = remove.isolate, 
                            top = top, color.use = color.use, vertex.weight = vertex.weight, 
                            vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                            weight.scale = weight.scale, edge.weight.max = edge.weight.max, 
                            edge.width.max = edge.width.max, title.name = paste0(signaling.name, 
                                                                                 " signaling pathway network"), vertex.label.cex = vertex.label.cex, 
                            ...)
  }
  else if (layout == "chord") {
    prob.sum <- apply(prob, c(1, 2), sum)
    gg <- netVisual_chord_cell_internal(prob.sum, color.use = color.use, 
                                        sources.use = sources.use, targets.use = targets.use, 
                                        remove.isolate = remove.isolate, group = group, 
                                        cell.order = cell.order, lab.cex = vertex.label.cex, 
                                        small.gap = small.gap, big.gap = big.gap, scale = scale, 
                                        reduce = reduce, title.name = paste0(signaling.name, 
                                                                             " signaling pathway network"), show.legend = show.legend, 
                                        legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
  }
  return(gg)
}


#color my way
netVisual_bubble <- function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
                              pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, 
                              sort.by.source.priority = TRUE, color.heatmap = c("darkolivegreen1", "darkolivegreen"), n.colors = 10, direction = 1, thresh = 0.05, 
                              comparison = NULL, group = NULL, remove.isolate = FALSE, 
                              max.dataset = NULL, min.dataset = NULL, min.quantile = 0, 
                              max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, 
                              color.text = NULL, title.name = NULL, font.size = 10, font.size.title = 10, 
                              show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
                              angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE) 
{
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  }
  else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle = c(0, 45, 90)
    hjust = c(0, 1, 1)
    vjust = c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }
  else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      pairLR.use$pathway_name <- as.character(pairLR.use$pathway_name)
    }
    else if ("interaction_name" %in% colnames(pairLR.use)) {
      pairLR.use$interaction_name <- as.character(pairLR.use$interaction_name)
    }
  }
  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net", 
                                  sources.use = sources.use, targets.use = targets.use, 
                                  signaling = signaling, pairLR.use = pairLR.use, 
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, 
                                  sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), 
                           targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                             ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, 
                              " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T) * 
                             1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                         levels(df.net$target), sep = " -> ")
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2), 
    ])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, 
                                        levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, 
                                   levels = cells.order)
    df <- df.net
  }
  else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net", 
                                      sources.use = sources.use, targets.use = targets.use, 
                                      signaling = signaling, pairLR.use = pairLR.use, 
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, 
                                    sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), 
                             targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, 
                                       unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                               ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, 
                                " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), 
                               each = length(levels(df.net$target))), levels(df.net$target), 
                           sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
                            ")")
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 
                      0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      }
      else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                                       ncol = 5))
        colnames(df.net) <- c("interaction_name_2", 
                              "source.target", "prob", "pval", "prob.original")
        df.net$source.target <- group.names0
      }
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, 
                                     " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T) * 
                             1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2), ])
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = unique(df$interaction_name_2))
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], 
                                             " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, 
                                dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
  max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        dataset.na <- c(df.i.j$dataset[is.na(values)], 
                        setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          }
          else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            }
            else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
                                                    dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  if (!is.null(pairLR.use)) {
    interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name, 
    ]$interaction_name_2, unique(df$interaction_name_2))
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = interaction_name_2.order)
  }
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
                                                                    unique(df$source.target)))
  if (sort.by.target & !sort.by.source) {
    if (!is.null(targets.use)) {
      df$target <- factor(df$target, levels = intersect(targets.use, 
                                                        df$target))
      df <- with(df, df[order(target, source), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & !sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      df <- with(df, df[order(source, target), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      if (!is.null(targets.use)) {
        df$target <- factor(df$target, levels = intersect(targets.use, 
                                                          df$target))
      }
      if (sort.by.source.priority) {
        df <- with(df, df[order(source, target), ])
      }
      else {
        df <- with(df, df[order(target, source), ])
      }
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, 
                      color = prob, size = pval)) + geom_point(pch = 16) + 
    theme_linedraw() + theme(panel.grid.major = element_blank()) + 
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, 
                                     vjust = vjust.x), axis.title.x = element_blank(), 
          axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
  values <- c(1, 2, 3)
  names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), 
                        breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
                                                                                 sort(unique(df$pval))], name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                    na.value = "white", limits = c(quantile(df$prob, 
                                                                            0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
                                    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 
                                                                                         1, na.rm = T)), labels = c("min", "max")) + 
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }
  else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                    na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5, 
                                                                                         title = "Commun. Prob."))
  }
  g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5 + length(dataset.name[comparison]), 
                       length(group.names0) * length(dataset.name[comparison]), 
                       by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype = "dashed", 
                          color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      }
      else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, 
                                               "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 
                                             2, stringr::str_length(dataset.name.order) - 
                                               1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  }
  else {
    return(g)
  }
}

subsetCommunication_internal <- function (net, LR, cells.level, slot.name = "net", sources.use = NULL, 
                                          targets.use = NULL, signaling = NULL, pairLR.use = NULL, 
                                          thresh = 0.05, datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL, 
                                          ligand.pct.1 = NULL, ligand.pct.2 = NULL, receptor.pvalues = NULL, 
                                          receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL) 
{
  if (!is.data.frame(net)) {
    prob <- net$prob
    pval <- net$pval
    prob[pval >= thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source", "target", "interaction_name")
    net.pval <- reshape2::melt(pval, value.name = "pval")
    net$pval <- net.pval$pval
    net <- subset(net, prob > 0)
  }
  if (!("ligand" %in% colnames(net))) {
    pairLR <- dplyr::select(LR, c("interaction_name_2", 
                                  "pathway_name", "ligand", "receptor", "annotation", 
                                  "evidence"))
    idx <- match(net$interaction_name, rownames(pairLR))
    net <- cbind(net, pairLR[idx, ])
  }
  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], 
                                 pairLR.use = LR, key = "pathway_name", matching.exact = T, 
                                 pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }
  if (!is.null(pairLR.use)) {
    net <- tryCatch({
      subset(net, interaction_name %in% pairLR.use$interaction_name)
    }, error = function(e) {
      subset(net, pathway_name %in% pairLR.use$pathway_name)
    })
  }
  if (!is.null(datasets)) {
    if (!("datasets" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before selecting 'datasets'")
    }
    net <- net[net$datasets %in% datasets, , drop = FALSE]
  }
  if (!is.null(ligand.pvalues)) {
    if (!("ligand.pvalues" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pvalues'")
    }
    net <- net[net$ligand.pvalues <= ligand.pvalues, , drop = FALSE]
  }
  if (!is.null(ligand.logFC)) {
    if (!("ligand.logFC" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.logFC'")
    }
    if (ligand.logFC >= 0) {
      net <- net[net$ligand.logFC >= ligand.logFC, , drop = FALSE]
    }
    else {
      net <- net[net$ligand.logFC <= ligand.logFC, , drop = FALSE]
    }
  }
  if (!is.null(ligand.pct.1)) {
    if (!("ligand.pct.1" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pct.1'")
    }
    net <- net[net$ligand.pct.1 >= ligand.pct.1, , drop = FALSE]
  }
  if (!is.null(ligand.pct.2)) {
    if (!("ligand.pct.2" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pct.2'")
    }
    net <- net[net$ligand.pct.2 >= ligand.pct.2, , drop = FALSE]
  }
  if (!is.null(receptor.pvalues)) {
    if (!("receptor.pvalues" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pvalues'")
    }
    net <- net[net$receptor.pvalues <= receptor.pvalues, 
               , drop = FALSE]
  }
  if (!is.null(receptor.logFC)) {
    if (!("receptor.logFC" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.logFC'")
    }
    if (receptor.logFC >= 0) {
      net <- net[net$receptor.logFC >= receptor.logFC, 
                 , drop = FALSE]
    }
    else {
      net <- net[net$receptor.logFC <= receptor.logFC, 
                 , drop = FALSE]
    }
  }
  if (!is.null(receptor.pct.1)) {
    if (!("receptor.pct.1" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pct.1'")
    }
    net <- net[net$receptor.pct.1 >= receptor.pct.1, , drop = FALSE]
  }
  if (!is.null(receptor.pct.2)) {
    if (!("receptor.pct.2" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pct.2'")
    }
    net <- net[net$receptor.pct.2 >= receptor.pct.2, , drop = FALSE]
  }
  net <- net[rowSums(is.na(net)) != ncol(net), , drop = FALSE]
  if (nrow(net) == 0) {
    warning("No significant signaling interactions are inferred based on the input!")
  }
  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source", "target", "pathway_name", 
                                "prob", "pval", "annotation"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    net.pval <- net %>% group_by(source_target, pathway_name) %>% 
      summarize(pval = mean(pval), .groups = "drop")
    net <- net %>% group_by(source_target, pathway_name) %>% 
      summarize(prob = sum(prob), .groups = "drop")
    a <- stringr::str_split(net$source_target, "sourceTotarget", 
                            simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net <- dplyr::select(net, -source_target)
    net$pval <- net.pval$pval
  }
  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  net <- BiocGenerics::as.data.frame(net, stringsAsFactors = FALSE)
  if (nrow(net) == 0) {
    warning("No significant signaling interactions are inferred!")
  }
  else {
    rownames(net) <- 1:nrow(net)
  }
  if (slot.name == "net") {
    if (("ligand.logFC" %in% colnames(net)) & ("datasets" %in% 
                                               colnames(net))) {
      net <- net[, c("source", "target", "ligand", "receptor", 
                     "prob", "pval", "interaction_name", "interaction_name_2", 
                     "pathway_name", "annotation", "evidence", "datasets", 
                     "ligand.logFC", "ligand.pct.1", "ligand.pct.2", 
                     "ligand.pvalues", "receptor.logFC", "receptor.pct.1", 
                     "receptor.pct.2", "receptor.pvalues")]
    }
    else if ("ligand.logFC" %in% colnames(net)) {
      net <- net[, c("source", "target", "ligand", "receptor", 
                     "prob", "pval", "interaction_name", "interaction_name_2", 
                     "pathway_name", "annotation", "evidence", "ligand.logFC", 
                     "ligand.pct.1", "ligand.pct.2", "ligand.pvalues", 
                     "receptor.logFC", "receptor.pct.1", "receptor.pct.2", 
                     "receptor.pvalues")]
    }
    else {
      net <- net[, c("source", "target", "ligand", "receptor", 
                     "prob", "pval", "interaction_name", "interaction_name_2", 
                     "pathway_name", "annotation", "evidence")]
    }
  }
  else if (slot.name == "netP") {
    net <- net[, c("source", "target", "pathway_name", "prob", 
                   "pval")]
  }
  return(net)
}

subsetCommunication <- function (object = NULL, net = NULL, slot.name = "net", sources.use = NULL, 
                                 targets.use = NULL, signaling = NULL, pairLR.use = NULL, 
                                 thresh = 0.05, datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL, 
                                 ligand.pct.1 = NULL, ligand.pct.2 = NULL, receptor.pvalues = NULL, 
                                 receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL) 
{
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }
  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }
  if (object@options$mode == "single") {
    if (is.null(net)) {
      net <- slot(object, "net")
    }
    LR <- object@LR$LRsig
    cells.level <- levels(object@idents)
    df.net <- subsetCommunication_internal(net, LR, cells.level, 
                                           slot.name = slot.name, sources.use = sources.use, 
                                           targets.use = targets.use, signaling = signaling, 
                                           pairLR.use = pairLR.use, thresh = thresh, datasets = datasets, 
                                           ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC, 
                                           ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2, 
                                           receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC, 
                                           receptor.pct.1 = receptor.pct.1, receptor.pct.2 = receptor.pct.2)
  }
  else if (object@options$mode == "merged") {
    if (is.null(net)) {
      net0 <- slot(object, "net")
      df.net <- vector("list", length(net0))
      names(df.net) <- names(net0)
      for (i in 1:length(net0)) {
        net <- net0[[i]]
        LR <- object@LR[[i]]$LRsig
        cells.level <- levels(object@idents[[i]])
        df.net[[i]] <- subsetCommunication_internal(net, 
                                                    LR, cells.level, slot.name = slot.name, sources.use = sources.use, 
                                                    targets.use = targets.use, signaling = signaling, 
                                                    pairLR.use = pairLR.use, thresh = thresh, 
                                                    datasets = datasets, ligand.pvalues = ligand.pvalues, 
                                                    ligand.logFC = ligand.logFC, ligand.pct.1 = ligand.pct.1, 
                                                    ligand.pct.2 = ligand.pct.2, receptor.pvalues = receptor.pvalues, 
                                                    receptor.logFC = receptor.logFC, receptor.pct.1 = receptor.pct.1, 
                                                    receptor.pct.2 = receptor.pct.2)
      }
    }
    else {
      LR <- data.frame()
      for (i in 1:length(object@LR)) {
        LR <- rbind(LR, object@LR[[i]]$LRsig)
      }
      LR <- unique(LR)
      cells.level <- levels(object@idents$joint)
      df.net <- subsetCommunication_internal(net, LR, 
                                             cells.level, slot.name = slot.name, sources.use = sources.use, 
                                             targets.use = targets.use, signaling = signaling, 
                                             pairLR.use = pairLR.use, thresh = thresh, datasets = datasets, 
                                             ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC, 
                                             ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2, 
                                             receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC, 
                                             receptor.pct.1 = receptor.pct.1, receptor.pct.2 = receptor.pct.2)
    }
  }
  return(df.net)
}

#Change the color of the dots and scale per row
netVisual_bubble <- function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
                              pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, 
                              sort.by.source.priority = TRUE, color.heatmap = c("lightblue1", "midnightblue"), n.colors = 10, direction = 1, thresh = 0.05, 
                              comparison = NULL, group = NULL, remove.isolate = FALSE, 
                              max.dataset = NULL, min.dataset = NULL, min.quantile = 0, 
                              max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, 
                              color.text = NULL, title.name = NULL, font.size = 10, font.size.title = 10, 
                              show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
                              angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE) 
{
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  }
  else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle = c(0, 45, 90)
    hjust = c(0, 1, 1)
    vjust = c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }
  else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      pairLR.use$pathway_name <- as.character(pairLR.use$pathway_name)
    }
    else if ("interaction_name" %in% colnames(pairLR.use)) {
      pairLR.use$interaction_name <- as.character(pairLR.use$interaction_name)
    }
  }
  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net", 
                                  sources.use = sources.use, targets.use = targets.use, 
                                  signaling = signaling, pairLR.use = pairLR.use, 
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, 
                                  sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), 
                           targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                             ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, 
                              " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T) * 
                             1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                         levels(df.net$target), sep = " -> ")
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2), 
    ])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, 
                                        levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, 
                                   levels = cells.order)
    df <- df.net
  }
  else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net", 
                                      sources.use = sources.use, targets.use = targets.use, 
                                      signaling = signaling, pairLR.use = pairLR.use, 
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, 
                                    sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), 
                             targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, 
                                       unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                               ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, 
                                " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), 
                               each = length(levels(df.net$target))), levels(df.net$target), 
                           sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
                            ")")
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 
                      0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      }
      else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                                       ncol = 5))
        colnames(df.net) <- c("interaction_name_2", 
                              "source.target", "prob", "pval", "prob.original")
        df.net$source.target <- group.names0
      }
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, 
                                     " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T) * 
                             1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2), ])
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = unique(df$interaction_name_2))
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], 
                                             " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, 
                                dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  
  #scale per interaction_name  
  for(i in unique(df$interaction_name)) {
    df.prob <- df[df$interaction_name==i, "prob"]
    if(length(df.prob)>1){
      df.prob <- scale(df.prob)
      df.prob -> df[df$interaction_name==i, "prob"]
    }
  }
  
  min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
  max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        dataset.na <- c(df.i.j$dataset[is.na(values)], 
                        setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          }
          else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            }
            else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
                                                    dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  if (!is.null(pairLR.use)) {
    interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name, 
    ]$interaction_name_2, unique(df$interaction_name_2))
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = interaction_name_2.order)
  }
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
                                                                    unique(df$source.target)))
  if (sort.by.target & !sort.by.source) {
    if (!is.null(targets.use)) {
      df$target <- factor(df$target, levels = intersect(targets.use, 
                                                        df$target))
      df <- with(df, df[order(target, source), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & !sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      df <- with(df, df[order(source, target), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      if (!is.null(targets.use)) {
        df$target <- factor(df$target, levels = intersect(targets.use, 
                                                          df$target))
      }
      if (sort.by.source.priority) {
        df <- with(df, df[order(source, target), ])
      }
      else {
        df <- with(df, df[order(target, source), ])
      }
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, 
                      color = prob, size = pval)) + geom_point(pch = 16) + 
    theme_linedraw() + theme(panel.grid.major = element_blank()) + 
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, 
                                     vjust = vjust.x), axis.title.x = element_blank(), 
          axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
  values <- c(1, 2, 3)
  names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), 
                        breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
                                                                                 sort(unique(df$pval))], name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                    na.value = "white", limits = c(quantile(df$prob, 
                                                                            0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
                                    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 
                                                                                         1, na.rm = T)), labels = c("min", "max")) + 
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }
  else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                    na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5, 
                                                                                         title = "Commun. Prob."))
  }
  g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5 + length(dataset.name[comparison]), 
                       length(group.names0) * length(dataset.name[comparison]), 
                       by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype = "dashed", 
                          color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      }
      else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, 
                                               "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 
                                             2, stringr::str_length(dataset.name.order) - 
                                               1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  }
  else {
    return(g)
  }
}

#order based on signaling, if provided
rankNet <- function (object, slot.name = "netP", measure = c("weight", 
                                                             "count"), mode = c("comparison", "single"), comparison = c(1, 
                                                                                                                        2), color.use = NULL, stacked = FALSE, sources.use = NULL, 
                     targets.use = NULL, signaling = NULL, pairLR = NULL, signaling.type = NULL, 
                     do.stat = FALSE, cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, 
                     show.raw = FALSE, return.data = FALSE, x.rotation = 90, 
                     title = NULL, bar.w = 0.75, font.size = 8, do.flip = TRUE, 
                     x.angle = NULL, y.angle = 0, x.hjust = 1, y.hjust = 1, axis.gap = FALSE, 
                     ylim = NULL, segments = NULL, tick_width = NULL, rel_heights = c(0.9, 
                                                                                      0, 0.1), norm_factor=NULL) 
{
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  options(warn = -1)
  object.names <- names(methods::slot(object, slot.name))
  if (measure == "weight") {
    ylabel = "Information flow"
  }
  else if (measure == "count") {
    ylabel = "Number of interactions"
  }
  if (mode == "single") {
    object1 <- methods::slot(object, slot.name)
    prob = object1$prob
    prob[object1$pval > thresh] <- 0
    if (measure == "count") {
      prob <- 1 * (prob > 0)
    }
    if (!is.null(sources.use)) {
      if (is.character(sources.use)) {
        if (all(sources.use %in% dimnames(prob)[[1]])) {
          sources.use <- match(sources.use, dimnames(prob)[[1]])
        }
        else {
          stop("The input `sources.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), sources.use)
      prob[idx.t, , ] <- 0
    }
    if (!is.null(targets.use)) {
      if (is.character(targets.use)) {
        if (all(targets.use %in% dimnames(prob)[[1]])) {
          targets.use <- match(targets.use, dimnames(prob)[[2]])
        }
        else {
          stop("The input `targets.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), targets.use)
      prob[, idx.t, ] <- 0
    }
    if (sum(prob) == 0) {
      stop("No inferred communications for the input!")
    }
    pSum <- apply(prob, 3, sum)
    pSum.original <- pSum
    if (measure == "weight") {
      pSum <- -1/log(pSum)
      pSum[is.na(pSum)] <- 0
      idx1 <- which(is.infinite(pSum) | pSum < 0)
      values.assign <- seq(max(pSum) * 1.1, max(pSum) * 
                             1.5, length.out = length(idx1))
      position <- sort(pSum.original[idx1], index.return = TRUE)$ix
      pSum[idx1] <- values.assign[match(1:length(idx1), 
                                        position)]
    }
    else if (measure == "count") {
      pSum <- pSum.original
    }
    pair.name <- names(pSum)
    df <- data.frame(name = pair.name, contribution = pSum.original, 
                     contribution.scaled = pSum, group = object.names[comparison[1]])
    idx <- with(df, order(df$contribution))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    for (i in 1:length(pair.name)) {
      df.t <- df[df$name == pair.name[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name[i]), ]
      }
    }
    if (!is.null(signaling.type)) {
      LR <- subset(object@DB$interaction, annotation %in% 
                     signaling.type)
      if (slot.name == "netP") {
        signaling <- unique(LR$pathway_name)
      }
      else if (slot.name == "net") {
        pairLR <- LR$interaction_name
      }
    }
    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    }
    else if ((slot.name == "netP") && (!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    }
    else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    gg <- ggplot(df, aes(x = name, y = contribution.scaled)) + 
      geom_bar(stat = "identity", width = bar.w) + theme_classic() + 
      theme(axis.text = element_text(size = 10), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), axis.title.y = element_text(size = 10)) + 
      xlab("") + ylab(ylabel) + coord_flip()
    if (!is.null(title)) {
      gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  else if (mode == "comparison") {
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    for (i in 1:length(comparison)) {
      object.list <- methods::slot(object, slot.name)[[comparison[i]]]
      prob <- object.list$prob
      if (!is.null(norm_factor)){
        norm_fac <- norm_factor[[i]]
        prob <- as.data.frame(as.table(prob)) %>%
          rename(Var1="Sender", Var2="Receiver", Freq="Comm_Prob") %>%
          left_join(norm_fac, by = c("Sender", "Receiver"))
        prob <- prob %>% mutate(Norm_Comm_Prob = Comm_Prob / Geo_Mean)
        prob <- prob[,c("Sender", "Receiver", "Var3", "Norm_Comm_Prob")]
        prob <- data.frame(prob)
        prob <- with(prob, array(data=prob$Norm_Comm_Prob, dim=c(length(unique(prob$Sender)),
                                                                 length(unique(prob$Receiver)),
                                                                 length(unique(prob$Var3))),
                                 dimnames=list(unique(prob$Sender), unique(prob$Receiver), unique(prob$Var3))))
      }
      
      prob[object.list$pval > thresh] <- 0
      if (measure == "count") {
        prob <- 1 * (prob > 0)
      }
      prob.list[[i]] <- prob
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% dimnames(prob)[[1]])) {
            sources.use <- match(sources.use, dimnames(prob)[[1]])
          }
          else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), sources.use)
        prob[idx.t, , ] <- 0
      }
      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% dimnames(prob)[[1]])) {
            targets.use <- match(targets.use, dimnames(prob)[[2]])
          }
          else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), targets.use)
        prob[, idx.t, ] <- 0
      }
      if (sum(prob) == 0) {
        stop("No inferred communications for the input!")
      }
      pSum.original[[i]] <- apply(prob, 3, sum)
      if (measure == "weight") {
        pSum[[i]] <- -1/log(pSum.original[[i]])
        pSum[[i]][is.na(pSum[[i]])] <- 0
        idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 
                            0)
        pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
      }
      else if (measure == "count") {
        pSum[[i]] <- pSum.original[[i]]
      }
      pair.name[[i]] <- names(pSum.original[[i]])
      object.names.comparison <- c(object.names.comparison, 
                                   object.names[comparison[i]])
    }
    if (measure == "weight") {
      values.assign <- seq(max(unlist(pSum)) * 1.1, max(unlist(pSum)) * 
                             1.5, length.out = length(unlist(idx)))
      position <- sort(pSum.original.all, index.return = TRUE)$ix
      for (i in 1:length(comparison)) {
        if (i == 1) {
          pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), 
                                                     position)]
        }
        else {
          pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i - 
                                                                         1])) + 1:length(unlist(idx[1:i])), position)]
        }
      }
    }
    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
      df[[i]] <- data.frame(name = pair.name.all, contribution = 0, 
                            contribution.scaled = 0, group = object.names[comparison[i]], 
                            row.names = pair.name.all)
      df[[i]][pair.name[[i]], 3] <- pSum[[i]]
      df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]
    }
    contribution.relative <- list()
    for (i in 1:(length(comparison) - 1)) {
      contribution.relative[[i]] <- as.numeric(format(df[[length(comparison) - 
                                                            i + 1]]$contribution/df[[1]]$contribution, digits = 1))
      contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    names(contribution.relative) <- paste0("contribution.relative.", 
                                           1:length(contribution.relative))
    for (i in 1:length(comparison)) {
      for (j in 1:length(contribution.relative)) {
        df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
      }
    }
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    if (length(comparison) == 2) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 contribution, -contribution.data2))
    }
    else if (length(comparison) == 3) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, contribution, -contribution.data2))
    }
    else if (length(comparison) == 4) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, -contribution.relative.3, 
                                 contribution, -contribution.data2))
    }
    else {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, -contribution.relative.3, 
                                 -contribution.relative.4, contribution, -contribution.data2))
    }
    for (i in 1:length(comparison)) {
      df[[i]] <- df[[i]][idx, ]
      df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL
    df <- do.call(rbind, df)
    df$group <- factor(df$group, levels = object.names.comparison)
    if (is.null(color.use)) {
      color.use = ggPalette(length(comparison))
    }
    df$group <- factor(df$group, levels = rev(levels(df$group)))
    color.use <- rev(color.use)
    if (do.stat & length(comparison) == 2) {
      for (i in 1:length(pair.name.all)) {
        if (nrow(prob.list[[j]]) != nrow(prob.list[[1]])) {
          stop("Statistical test is not applicable to datasets with different cellular compositions! Please set `do.stat = FALSE`")
        }
        prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * 
                                nrow(prob.list[[1]]), ncol = length(comparison))
        for (j in 1:length(comparison)) {
          if (pair.name.all[i] %in% comparison[[j]]) {
            prob.values[, j] <- as.vector(prob.list[[j]][, 
                                                         , pair.name.all[i]])
          }
          else {
            prob.values[, j] <- NA
          }
        }
        prob.values <- prob.values[rowSums(prob.values, 
                                           na.rm = TRUE) != 0, , drop = FALSE]
        if (nrow(prob.values) > 3 & sum(is.na(prob.values)) == 
            0) {
          pvalues <- wilcox.test(prob.values[, 1], prob.values[, 2], exact=F)$p.value
        }
        else {
          pvalues <- 0
        }
        pvalues[is.na(pvalues)] <- 0
        df$pvalues[df$name == pair.name.all[i]] <- pvalues
      }
    }
    if (length(comparison) == 2) {
      if (do.stat) {
        colors.text <- ifelse((df$contribution.relative < 
                                 1 - tol) & (df$pvalues < cutoff.pvalue), color.use[2], 
                              ifelse((df$contribution.relative > 1 + tol) & 
                                       df$pvalues < cutoff.pvalue, color.use[1], 
                                     "black"))
      }
      else {
        colors.text <- ifelse(df$contribution.relative < 
                                1 - tol, color.use[2], ifelse(df$contribution.relative > 
                                                                1 + tol, color.use[1], "black"))
      }
    }
    else {
      message("The text on the y-axis will not be colored for the number of compared datasets larger than 3!")
      colors.text = NULL
    }
    for (i in 1:length(pair.name.all)) {
      df.t <- df[df$name == pair.name.all[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name.all[i]), 
        ]
      }
    }
    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    }
    else if ((slot.name == "netP") && (!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    }
    else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    if(!is.null(signaling)){
      df$name <- factor(df$name, levels = signaling)
    }
    if (stacked) {
      gg <- ggplot(df, aes(x = name, y = contribution, 
                           fill = group)) + geom_bar(stat = "identity", 
                                                     width = bar.w, position = "fill")
      if (measure == "weight") {
        gg <- gg + xlab("") + ylab("Relative information flow")
      }
      else if (measure == "count") {
        gg <- gg + xlab("") + ylab("Relative number of interactions")
      }
      gg <- gg + geom_hline(yintercept = 0.5, linetype = "dashed", 
                            color = "grey50", size = 0.5)
    }
    else {
      if (show.raw) {
        gg <- ggplot(df, aes(x = name, y = contribution, 
                             fill = group)) + geom_bar(stat = "identity", 
                                                       width = bar.w, position = position_dodge(0.8)) + 
          xlab("") + ylab(ylabel)
      }
      else {
        gg <- ggplot(df, aes(x = name, y = contribution.scaled, 
                             fill = group)) + geom_bar(stat = "identity", 
                                                       width = bar.w, position = position_dodge(0.8)) + 
          xlab("") + ylab(ylabel)
      }
      if (axis.gap) {
        gg <- gg + theme_bw() + theme(panel.grid = element_blank())
        gg.gap::gg.gap(gg, ylim = ylim, segments = segments, 
                       tick_width = tick_width, rel_heights = rel_heights)
      }
    }
    gg <- gg + CellChat_theme_opts() + theme_classic()
    if (do.flip) {
      gg <- gg + coord_flip() + theme(axis.text.y = element_text(colour = colors.text))
      if (is.null(x.angle)) {
        x.angle = 0
      }
    }
    else {
      if (is.null(x.angle)) {
        x.angle = 45
      }
      gg <- gg + scale_x_discrete(limits = rev) + theme(axis.text.x = element_text(colour = rev(colors.text)))
    }
    gg <- gg + theme(axis.text = element_text(size = font.size), 
                     axis.title.y = element_text(size = font.size))
    gg <- gg + scale_fill_manual(name = "", values = color.use)
    gg <- gg + guides(fill = guide_legend(reverse = TRUE))
    gg <- gg + theme(axis.text.x = element_text(angle = x.angle, 
                                                hjust = x.hjust), axis.text.y = element_text(angle = y.angle, 
                                                                                             hjust = y.hjust))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  }
  else {
    return(gg)
  }
}


#Make true comparison heatmap (scale by row)

netAnalysis_signalingRole_heatmap_combofun <- function(mat, pattern = c("outgoing", 
                                                                        "incoming", "all"), signaling=NULL, color.use = NULL, color.heatmap.use = c("white", "black"), filename = NULL,
                                                       font.size = 7, font.size.title = 7, cluster.rows = FALSE, cluster.cols = FALSE) 
{
  pattern <- match.arg(pattern)
  if (pattern == "outgoing") {
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    legend.name <- "Overall"
  }
  
  title <- paste0(legend.name, " signaling patterns")
  
  mat$group <- make.unique(mat$group)
  mat <- column_to_rownames(mat, "group")
  if(!is.null(signaling)){
    signaling <- signaling[signaling %in% colnames(mat)]
    mat <- mat[colnames(mat) %in% signaling]
    mat <- mat[,match(c(rev(signaling)), colnames(mat))]
  }
  mat[mat == 0] <- NA
  mat <- t(scale(mat)) #scale
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  
  color.heatmap.use <- colorRamp2(c(ceiling(min(mat, na.rm=T)), floor(max(mat, na.rm=T))), color.heatmap.use)
  
  df <- data.frame(group = colnames(mat))
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = "Relative strength", bottom_annotation = col_annotation, cluster_rows = cluster.rows, 
                cluster_columns = cluster.cols, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), 
                column_names_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                                   fontface = "plain"), title_position = "leftcenter-rot", 
                                                                   border = NA, at = legend.break, legend_height = unit(20, 
                                                                                                                        "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                                                                                                                                                                                 "mm")))
  draw(ht1)
  graph2pdf(ht1, paste0(curr_dir, "/CellChat/Comparison/Interaction_DESP_heatmap_", pattern, "_", filename, ".pdf"), width=3, height=6) #!!!!change save location
}

prep_mat <- function(object, slot.name, thresh, pattern){
  object <- netAnalysis_computeCentrality(object, slot.name = slot.name, thresh=thresh)
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- outgoing
  }
  else if (pattern == "incoming") {
    mat <- incoming
  }
  else if (pattern == "all") {
    mat <- outgoing + incoming
  }
  mat <- data.frame(mat)
  mat <- rownames_to_column(mat, "group")
  return(mat)
}



rbindf <- function(...) {
  
  l <- list(...)
  if(length(l) == 1) l <- l[[1]]
  nn <- length(l)
  
  x <- l[[1]]
  if(length(l)>1){
    for(i in 2:nn) {
      y <- l[[i]]
      if(nrow(x) > 0 & nrow(y) > 0) {
        if(!all(yinx <- names(y) %in% names(x))) {
          x[, names(y)[!yinx]] <- NA
        } 
        if(!all(xiny <- names(x) %in% names(y))) {
          y[, names(x)[!xiny]] <- NA
        } 
      }
      x <- rbind(x, y)
    }
  }
  return(x)
}
