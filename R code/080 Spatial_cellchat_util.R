###
### 080 Spatial_cellchat_util.R
### These functions were written by Drs. Denis R Avey and Tristan J Philippe
# Purpose: Compare CellChat-derived signaling across three niches.
# Dependencies:
library(CellChat)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggplot2)
library(rlang)
library(circlize)
library(ComplexHeatmap)
library(export)
library(RColorBrewer)
library(stringr)
library(grid)
library(ggrepel)
library(gridExtra)
library(pheatmap)
library(future)
library(future.apply)
library(edgeR)
library(limma)
library(fgsea)
options(future.globals.maxSize = 6*1024^3)
options(future.seed = TRUE)
options(stringsAsFactors = FALSE)
set.seed(123)  # For reproducibility

# define function to calculate potential interacting cell partners in a given micron radius
ComputePotentialSpatialFractions <- function(cellchat,
                                             contact_thresh = 20,
                                             short_range_thresh = 100,
                                             long_range_thresh = 250, # default CellChat interaction.range
                                             max_cells = NULL,
                                             n_workers,
                                             log_file,
                                             output_csv) {
  if (future::supportsMulticore()) {
    plan(multicore, workers = n_workers)
  } else {
    plan(multisession, workers = n_workers)
  }
  
  coords <- cellchat@images[["coordinates"]]
  idents <- cellchat@idents
  names(idents) <- rownames(coords)
  coords <- coords[names(idents), , drop = FALSE]
  if (!all(names(idents) %in% rownames(coords))) stop("Mismatch between idents and coordinates.")
  
  # Convert coords to data.table
  coord_dt <- as.data.table(coords)
  coord_dt[, cell_id := rownames(coords)]
  coord_dt[, ident := idents[cell_id]]
  setnames(coord_dt, names(coord_dt)[1:2], c("x", "y"))
  
  # Precompute lookup by identity
  id_lookup <- split(coord_dt, coord_dt$ident)
  
  sender_types <- names(id_lookup)
  receiver_types <- names(id_lookup)
  
  # All possible sender-receiver combinations
  sender_receiver_grid <- expand.grid(sender = sender_types, receiver = receiver_types, stringsAsFactors = FALSE)
  total_pairs <- nrow(sender_receiver_grid)
  writeLines(paste0("Starting computation for ", total_pairs, " sender-receiver combinations..."), con = log_file)
  
  results <- future_lapply(seq_len(total_pairs), function(idx) {
    pair <- sender_receiver_grid[idx, ]
    send_type <- pair$sender
    recv_type <- pair$receiver
    
    send_cells <- id_lookup[[send_type]]
    recv_cells <- id_lookup[[recv_type]]
    if (nrow(send_cells) == 0 || nrow(recv_cells) == 0) return(NULL)
    
    if (!is.null(max_cells) && nrow(send_cells) > max_cells) {
      set.seed(123)
      send_cells <- send_cells[sample(.N, max_cells)]
    }
    
    pair_dists <- c()
    for (s in seq_len(nrow(send_cells))) {
      sx <- send_cells$x[s]; sy <- send_cells$y[s]; sid <- send_cells$cell_id[s]
      
      # Axis-aligned bounding box filtering
      rcandidates <- recv_cells[
        x >= (sx - long_range_thresh) & x <= (sx + long_range_thresh) &
          y >= (sy - long_range_thresh) & y <= (sy + long_range_thresh)
      ]
      
      if (nrow(rcandidates) == 0) next
      
      rcandidates[, dist := sqrt((x - sx)^2 + (y - sy)^2)]
      rcandidates <- rcandidates[dist <= long_range_thresh]
      
      if (send_type == recv_type) rcandidates <- rcandidates[cell_id != sid]
      
      pair_dists <- c(pair_dists, rcandidates$dist)
    }
    
    total <- length(pair_dists)
    if (total == 0) return(NULL)
    
    contact <- sum(pair_dists <= contact_thresh)
    short_range <- sum(pair_dists > contact_thresh & pair_dists < short_range_thresh)
    long_range <- sum(pair_dists >= short_range_thresh & pair_dists <= long_range_thresh)
    
    if (idx %% 100 == 0) {
      msg <- paste(Sys.time(), "Processed", idx, "of", total_pairs, "sender-receiver pairs")
      write(msg, file = log_file, append = TRUE)
    }
    
    data.frame(
      sender = send_type,
      receiver = recv_type,
      n_pairs = total,
      contact_fraction = contact / total,
      short_range_fraction = short_range / total,
      long_range_fraction = long_range / total
    )
  })
  
  write(paste(Sys.time(), "Finished computation."), file = log_file, append = TRUE)
  
  final_df <- do.call(rbind, results)
  write.csv(final_df, output_csv, row.names = FALSE)
  future::plan(future::sequential)
  return(final_df)
}


# normalize CC objects based on potential interacting pairs (within 250 um)
# define normalization function
normalize_cellchat_by_log_npairs <- function(cellchat_obj, potential_df) {
  # Step 1: Pivot to wide matrix (sender x receiver) of n_pairs
  n_pairs_wide <- potential_df %>%
    select(sender, receiver, n_pairs) %>%
    pivot_wider(names_from = receiver, values_from = n_pairs)
  
  # Convert to matrix
  row_names <- n_pairs_wide$sender
  n_pairs_matrix <- as.matrix(n_pairs_wide[, -1])
  rownames(n_pairs_matrix) <- row_names
  
  # Step 2: Take log2(n_pairs + 1)
  log_npairs_matrix <- log2(n_pairs_matrix + 1)
  mean_log_npairs <- mean(log_npairs_matrix, na.rm = TRUE)
  
  # Step 3: Normalize net$prob
  net_prob_dims <- dim(cellchat_obj@net$prob)
  norm_net_prob <- cellchat_obj@net$prob / array(rep(log_npairs_matrix, net_prob_dims[3]),
                                                 dim = net_prob_dims)
  norm_net_prob <- norm_net_prob * mean_log_npairs
  
  # Step 4: Normalize netP$prob
  netP_prob_dims <- dim(cellchat_obj@netP$prob)
  norm_netP_prob <- cellchat_obj@netP$prob / array(rep(log_npairs_matrix, netP_prob_dims[3]),
                                                   dim = netP_prob_dims)
  norm_netP_prob <- norm_netP_prob * mean_log_npairs
  
  # Step 5: Normalize count and weight
  norm_count <- cellchat_obj@net$count / log_npairs_matrix * mean_log_npairs
  norm_weight <- cellchat_obj@net$weight / log_npairs_matrix * mean_log_npairs
  
  # Step 6: Assign normalized values
  cellchat_obj@net$prob <- norm_net_prob
  cellchat_obj@netP$prob <- norm_netP_prob
  cellchat_obj@net$count <- norm_count
  cellchat_obj@net$weight <- norm_weight
  
  return(cellchat_obj)
}


# define function for computing net centrality from cellchat object
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

# define custom rankNet function
rankNet <- function (object, slot.name = "netP", measure = c("weight", 
                                                             "count"), mode = c("comparison", "single"), comparison = c(1, 
                                                                                                                        2), color.use = NULL, stacked = FALSE, sources.use = NULL, 
                     targets.use = NULL, signaling = NULL, pairLR = NULL, signaling.type = NULL, 
                     do.stat = FALSE, cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, 
                     show.raw = FALSE, return.data = FALSE, x.rotation = 90, 
                     title = NULL, bar.w = 0.75, font.size = 8, do.flip = TRUE, 
                     x.angle = NULL, y.angle = 0, x.hjust = 1, y.hjust = 1, axis.gap = FALSE, 
                     ylim = NULL, segments = NULL, tick_width = NULL, rel_heights = c(0.9, 
                                                                                      0, 0.1)) 
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
          if (pair.name.all[i] %in% pair.name[[j]]) {
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
          pvalues <- wilcox.test(prob.values[, 1], prob.values[, 
                                                               2], paired = TRUE)$p.value
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

# define custom signalingRole_heatmap function for plotting differential net centrality
netAnalysis_signalingRole_heatmap_diff_or_ratio_comp <- function(mat,
                                                                 pattern = c("outgoing", "incoming", "all"),
                                                                 signaling = NULL,
                                                                 color.use = NULL,
                                                                 color.heatmap.use = c("blue", "white", "red"),
                                                                 filename = NULL,
                                                                 font.size = 7,
                                                                 font.size.title = 10,
                                                                 cluster.rows = FALSE,
                                                                 cluster.cols = FALSE,
                                                                 comparison.type = c("diff", "ratio", "diff.rowscaled"),
                                                                 significance.matrix = NULL) {
  
  pattern <- match.arg(pattern)
  comparison.type <- match.arg(comparison.type)
  
  legend.name <- switch(pattern,
                        outgoing = "Outgoing",
                        incoming = "Incoming",
                        all = "Overall")
  title <- switch(comparison.type,
                  diff = paste0(legend.name, " signaling difference"),
                  ratio = paste0(legend.name, " signaling log2 ratio"),
                  diff.rowscaled = paste0(legend.name, " signaling difference (Row-scaled signed min-max)"))
  
  mat$group <- make.unique(mat$group)
  mat <- column_to_rownames(mat, "group")
  
  if (nrow(mat) %% 2 != 0) stop("Matrix must have an even number of rows for paired comparison.")
  n_half <- nrow(mat) / 2
  mat1 <- mat[1:n_half, , drop = FALSE]
  mat2 <- mat[(n_half + 1):(2 * n_half), , drop = FALSE]
  rownames(mat1) <- rownames(mat2)
  
  mat_diff <- if (comparison.type == "diff") {
    mat2 - mat1
  } else {
    log2(mat2 / mat1)
  }
  mat_diff[!is.finite(as.matrix(mat_diff))] <- NA
  
  # Standardize rownames (remove .1 if present)
  rownames(mat_diff) <- gsub("\\.1$", "", rownames(mat_diff))
  
  if (!is.null(signaling)) {
    missing <- setdiff(signaling, colnames(mat_diff))
    if (length(missing) > 0) warning("Missing signaling pathways: ", paste(missing, collapse = ", "))
    signaling_use <- intersect(rev(signaling), colnames(mat_diff))
    mat_diff <- mat_diff[, signaling_use, drop = FALSE]
  }
  
  if (comparison.type == "diff.rowscaled") {
    # Column-wise signed min-max scaling preserving zero = white
    mat_scaled <- apply(mat_diff, 2, function(x) {
      scaled <- numeric(length(x))
      
      pos_idx <- which(!is.na(x) & x > 0)
      neg_idx <- which(!is.na(x) & x < 0)
      
      max_pos <- if(length(pos_idx) > 0) max(x[pos_idx], na.rm = TRUE) else NA
      max_neg <- if(length(neg_idx) > 0) min(x[neg_idx], na.rm = TRUE) else NA
      
      if (!is.na(max_pos) && max_pos != 0) {
        scaled[pos_idx] <- x[pos_idx] / max_pos
      } else {
        scaled[pos_idx] <- 0
      }
      
      if (!is.na(max_neg) && max_neg != 0) {
        scaled[neg_idx] <- x[neg_idx] / abs(max_neg)
      } else {
        scaled[neg_idx] <- 0
      }
      
      zero_idx <- which(!is.na(x) & x == 0)
      scaled[zero_idx] <- 0
      
      # For NA in original vector, set scaled to 0 or NA (your choice)
      na_idx <- which(is.na(x))
      scaled[na_idx] <- 0
      
      scaled
    })
    mat_scaled <- as.matrix(mat_scaled)
    rownames(mat_scaled) <- rownames(mat_diff)
    colnames(mat_scaled) <- colnames(mat_diff)
    color.heatmap.use <- colorRamp2(c(-1, 0, 1), color.heatmap.use)
  } else {
    mat_scaled <- mat_diff
    max_val <- max(abs(mat_scaled), na.rm = TRUE)
    #color.heatmap.use <- colorRamp2(c(-max_val/2, 0, max_val/2), color.heatmap.use) # define color scale
    color.heatmap.use <- colorRamp2(c(-0.002, 0, 0.002), color.heatmap.use) # or set values manually to enable direct comparison across datasets
  }
  
  legend.label <- switch(comparison.type,
                         diff = "Normalized Centrality Score (AD - NCI)",
                         ratio = "Normalized Centrality Score log2(AD/NCI)",
                         diff.rowscaled = "Normalized Centrality Score (row-scaled)")
  
  # Row color annotations
  row_anno <- NULL
  if (!is.null(color.use)) {
    ct <- rownames(mat_scaled)
    ct_colors <- color.use[ct]
    if (any(is.na(ct_colors))) {
      stop("Some cell types in rownames(mat_scaled) are missing from color.use: ",
           paste(unique(ct[is.na(ct_colors)]), collapse = ", "))
    }
    row_anno <- rowAnnotation(CellType = rownames(mat_scaled),
                              col = list(CellType = color.use),
                              show_annotation_name = FALSE,
                              simple_anno_size = unit(0.2, "cm"))
  }
  
  # Significance annotation (asterisks)
  text_anno <- NULL
  # Align significance.matrix if provided
  if (!is.null(significance.matrix)) {
    sig_rows_use <- intersect(rownames(mat_scaled), rownames(significance.matrix))
    sig_cols_use <- intersect(colnames(mat_scaled), colnames(significance.matrix))
    significance.matrix <- significance.matrix[sig_rows_use, sig_cols_use, drop = FALSE]
    
    # Reorder to match mat_scaled
    significance.matrix <- significance.matrix[rownames(mat_scaled), colnames(mat_scaled), drop = FALSE]
    
    stopifnot(all(dim(significance.matrix) == dim(mat_scaled)))
    stopifnot(all(rownames(significance.matrix) == rownames(mat_scaled)))
    stopifnot(all(colnames(significance.matrix) == colnames(mat_scaled)))
    
    label_mat <- significance.matrix
    
    text_anno <- function(j, i, x, y, width, height, fill) {
      if (label_mat[i, j] != "") {
        grid.text(label = label_mat[i, j], x = x, y = y,
                  gp = gpar(fontsize = 10, fontface = "bold"))
      }
    }
  }
  
  # Ensure correct row ordering
  row_order <- rev(seq_len(nrow(mat_scaled)))
  
  ht1 <- Heatmap(mat_scaled,
                 col = color.heatmap.use,
                 name = legend.label,
                 right_annotation = row_anno,
                 top_annotation = NULL,
                 cluster_rows = cluster.rows,
                 cluster_columns = cluster.cols,
                 row_names_side = "left",
                 row_names_rot = 0,
                 row_names_gp = gpar(fontsize = font.size),
                 column_names_gp = gpar(fontsize = font.size),
                 column_title = title,
                 column_title_gp = gpar(fontsize = font.size.title),
                 column_names_rot = 90,
                 na_col = "white",
                 cell_fun = text_anno,
                 row_order = row_order,
                 width = unit(60, "mm"),   
                 height = unit(40, "mm"),
                 heatmap_legend_param = list(
                   title_gp = gpar(fontsize = 8, fontface = "plain"),
                   title_position = "leftcenter-rot",
                   border = NA,
                   legend_height = unit(20, "mm"),
                   labels_gp = gpar(fontsize = 8),
                   grid_width = unit(2, "mm")
                 )
  )
  
  draw(ht1, heatmap_legend_side = "right")
  
  if (!is.null(filename)) {
    graph2pdf(ht1, paste0(curr_dir, "/CellChat/Comparison/Interaction_DESP_heatmap_", comparison.type, "_", pattern, "_", filename, ".pdf"),
              width = 4, height = 3)
  }
}


# define function for generating pseudobulk by cell type
pseudobulk_by_celltype <- function(seurat_obj, celltype_col = "newMajorcelltype", assay = "Nanostring") {
  meta <- seurat_obj@meta.data
  counts <- seurat_obj@assays[[assay]]@layers$counts
  celltypes <- unique(meta[[celltype_col]])
  genenames <- rownames(seurat_obj)
  results <- list()
  
  for (ct in celltypes) {
    message("Processing: ", ct)
    ct_cells <- rownames(meta)[meta[[celltype_col]] == ct]
    ct_meta <- meta[ct_cells, ]
    ct_counts <- counts[, ct_cells]
    
    # Sample ID is individual
    ct_meta$sample_id <- ct_meta$arb_ID
    df <- data.frame(cell_id = colnames(ct_counts), sample_id = ct_meta$sample_id)
    
    # Sum across cells from the same individual
    pseudobulk_mat <- as.data.frame(as.matrix(ct_counts)) %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "cell_id", values_to = "count") %>%
      left_join(df, by = "cell_id") %>%
      group_by(gene, sample_id) %>%
      summarize(count = sum(count), .groups = "drop") %>%
      pivot_wider(names_from = sample_id, values_from = count) %>%
      column_to_rownames("gene")
    
    results[[ct]] <- pseudobulk_mat
  }
  
  return(results)
}


# define DE function (AD vs. NCI)
run_de_by_celltype <- function(pseudobulk_list, seurat_obj, celltype_col = "newMajorcelltype", assay = "Nanostring") {
  meta <- seurat_obj@meta.data
  raw_counts <- seurat_obj@assays[[assay]]@layers$counts
  
  all_de_results <- list()
  
  for (ct in names(pseudobulk_list)) {
    message("Running DE for: ", ct)
    pb <- pseudobulk_list[[ct]]  # genes x samples (individuals)
    
    # Metadata setup
    sample_ids <- colnames(pb)
    ct_meta <- meta %>% filter(arb_ID %in% sample_ids)
    ct_meta <- ct_meta[!duplicated(ct_meta$arb_ID), ]
    ct_meta <- ct_meta %>% filter(arb_ID %in% sample_ids)
    ct_meta <- ct_meta[match(sample_ids, ct_meta$arb_ID), ]
    
    # DE setup
    design <- model.matrix(~0 + AD_status, data = ct_meta)
    colnames(design) <- c("AD", "NCI")
    
    dge <- DGEList(pb)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design)
    fit <- lmFit(v, design)
    contrast <- makeContrasts(AD - NCI, levels = design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    
    de_table <- topTable(fit2, number = Inf, sort.by = "none")
    de_table$gene <- rownames(de_table)
    de_table$celltype <- ct
    
    # Add: total pseudobulk count across samples
    de_table$total_count <- rowSums(pb[rownames(de_table), , drop = FALSE], na.rm = TRUE)
    
    # Add: % of expressing cells (non-zero counts) in original data
    ct_cells <- rownames(meta)[meta[[celltype_col]] == ct]
    expr_mat <- raw_counts[, ct_cells, drop = FALSE]
    gene_expr_pct <- Matrix::rowMeans(expr_mat > 0) * 100
    de_table$pct_expressing_cells <- gene_expr_pct[rownames(de_table)]
    
    all_de_results[[ct]] <- de_table
  }
  
  all_de_results_df <- bind_rows(all_de_results)
  return(all_de_results_df)
}


# create net from de, including parsing of cellchat gene names
create_net_from_de <- function(net, de_results) {
  # Manual mapping for ligand/receptor complexes
  complex_gene_map <- list(
    "Activin AB" = c("INHBA", "INHBB"),
    "Inhibin A" = c("INHA", "INHBA"),
    "Inhibin B" = c("INHA", "INHBB"),
    "IL12AB" = c("IL12A", "IL12B"),
    "IL23 complex" = c("IL23A", "IL12B"),
    "IL27 complex" = c("EBI3", "IL27"),
    "LTa1b2" = c("LTA", "LTB"),
    "HT" = c("HDC"),       
    "Endorphin" = c("POMC"),
    "Ach-CHAT" = "CHAT",                   # Acetylcholine → Choline acetyltransferase
    "LTC4-LTC4S" = "LTC4S",                # Leukotriene C4 synthase
    "DOPA-DDC" = "DDC",                    # Dopamine → DOPA decarboxylase
    "LXA4-ALOX12" = "ALOX12",              # Lipoxin A4 synth via ALOX12
    "LXA4-ALOX5" = "ALOX5",                # Alternative route via ALOX5
    "GABA-GAD1" = "GAD1",                  # GABA via glutamate decarboxylase 1
    "GABA-GAD2" = "GAD2",                  # GABA via glutamate decarboxylase 2
    "Glu-SLC17A7" = "SLC17A7",             # VGLUT1 transporter
    "Glu-SLC1A1" = "SLC1A1",               # Glutamate reuptake transporter
    "Glu-SLC1A2" = "SLC1A2",               # EAAT2 (glial)
    "Glu-SLC1A3" = "SLC1A3",               # EAAT1 (glial)
    "Hist-HDC" = "HDC",                    # Histamine → histidine decarboxylase
    "5-HT-TPH2" = "TPH2",                  # Tryptophan hydroxylase 2 (neural serotonin)
    "b-Endorphin-POMC" = "POMC",           # β-Endorphin precursor
    "PGD2-PTGDS" = "PTGDS",                # Prostaglandin D2 synthase
    "PGE2-PTGES" = "PTGES",                # PGE2 synthase
    "PGE2-PTGES2" = "PTGES2",              # Alternate isoform
    "PGE2-PTGES3" = "PTGES3",              # Alternate isoform
    "PGI2-PTGIS" = "PTGIS",                # Prostacyclin synthase
    "TXA2-TBXAS1" = "TBXAS1",              # Thromboxane A2 synthase
    "RA-ALDH1A1" = "ALDH1A1",              # Retinoic acid from retinal
    "RA-ALDH1A3" = "ALDH1A3",              # Isoform
    "Cholesterol-DHCR24" = "DHCR24",       # Cholesterol biosynthesis enzyme
    "Cholesterol-LIPA" = "LIPA",           # Lysosomal acid lipase
    "DHT-SRD5A1" = "SRD5A1",               # Dihydrotestosterone from testosterone
    "DHT-SRD5A2" = "SRD5A2",
    "DHT-SRD5A3" = "SRD5A3",
    "E2-CYP19A1" = "CYP19A1",              # Estrogen synth (aromatase)
    "5-HT-DOPA-DDC" = "DDC",               # DDC used again here
    "Testosterone-HSD17B12" = "HSD17B12",  # Testosterone synth from androstenedione
    "TGFbR1" = "TGFBR1",
    "TGFbR2" = "TGFBR2",
    "TGFbR" = c("TGFBR1", "TGFBR2", "TGFBR3"),
    "R2" = "TGFBR2",                    # Likely alias in some CC networks
    "LXA4" = "FPR2",                    # Lipoxin A4 receptor
    "CD45" = "PTPRC",
    "NKG2D" = "KLRK1",
    "CD8 receptor" = c("CD8A", "CD8B")
  )
  
  # Enhanced parsing function
  extract_genes <- function(name) {
    name <- trimws(name)
    if (name %in% names(complex_gene_map)) {
      return(complex_gene_map[[name]])
    }
    
    # Fallback: split on '_' and filter out empty strings
    if (str_detect(name, "_")) {
      parts <- str_split(name, "_", simplify = TRUE)
      return(parts[parts != ""])
    }
    
    # Final fallback: assume single gene name
    return(name)
  }
  
  # Apply to ligand and receptor columns
  net <- net %>%
    mutate(
      ligand_genes = lapply(ligand, extract_genes),
      receptor_genes = lapply(receptor, extract_genes)
    )
  
  # DE lookup table
  de_lookup <- de_results %>%
    select(gene, celltype, logFC, adj.P.Val, total_count, pct_expressing_cells) %>%
    rename(
      de_gene = gene,
      de_celltype = celltype,
      de_logFC = logFC,
      de_adjP = adj.P.Val,
      de_count = total_count,
      de_pct = pct_expressing_cells
    )
  
  # Lookup helpers
  get_avg_logFC <- function(gene_list, celltype) {
    vals <- de_lookup %>%
      filter(de_gene %in% gene_list, de_celltype == celltype) %>%
      pull(de_logFC)
    if (length(vals) == 0) return(NA_real_)
    mean(vals, na.rm = TRUE)
  }
  
  get_min_adjP <- function(gene_list, celltype) {
    vals <- de_lookup %>%
      filter(de_gene %in% gene_list, de_celltype == celltype) %>%
      pull(de_adjP)
    if (length(vals) == 0) return(NA_real_)
    min(vals, na.rm = TRUE)
  }
  
  get_total_count <- function(gene_list, celltype) {
    vals <- de_lookup %>%
      filter(de_gene %in% gene_list, de_celltype == celltype) %>%
      pull(de_count)
    if (length(vals) == 0) return(NA_real_)
    sum(vals, na.rm = TRUE)
  }
  
  get_pct_expressing <- function(gene_list, celltype) {
    vals <- de_lookup %>%
      filter(de_gene %in% gene_list, de_celltype == celltype) %>%
      pull(de_pct)
    if (length(vals) == 0) return(NA_real_)
    mean(vals, na.rm = TRUE)
  }
  
  # Apply gene list lookups row-wise
  net_with_fc <- net %>%
    rowwise() %>%
    mutate(
      ligand.logFC = get_avg_logFC(ligand_genes, source),
      receptor.logFC = get_avg_logFC(receptor_genes, target),
      ligand.adjP = get_min_adjP(ligand_genes, source),
      receptor.adjP = get_min_adjP(receptor_genes, target),
      ligand.total_count = get_total_count(ligand_genes, source),
      receptor.total_count = get_total_count(receptor_genes, target),
      ligand.pct_expressing = get_pct_expressing(ligand_genes, source),
      receptor.pct_expressing = get_pct_expressing(receptor_genes, target),
      combined_logFC = mean(c(ligand.logFC, receptor.logFC), na.rm = TRUE),
      combined_adjP = pmin(ligand.adjP, receptor.adjP, na.rm = TRUE)
    ) %>%
    ungroup()
  
  return(net_with_fc)
}

# function to run pseudobulk, de, and create net
process_niche <- function(seurat_obj, net_obj, niche_name) {
  # 1. Generate pseudobulk matrix
  pb_list <- pseudobulk_by_celltype(seurat_obj)
  
  # 2. Differential expression analysis
  de_results <- run_de_by_celltype(pb_list, seurat_obj)
  
  # 3. Annotate interaction network
  updated_net <- create_net_from_de(net = net_obj, de_results = de_results)
  
  # 4. Plotting
  major <- c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B") 
  
  p_combined <- ggplot(updated_net, aes(x = combined_logFC, y = -log10(combined_adjP), color = source)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = major) +
    theme_minimal() +
    labs(x = "Combined logFC", y = "-log10(combined adj P)", title = paste0("Niche ", niche_name, ": Combined DE Summary"), color = "Sender")
  
  return(list(
    pseudobulk = pb_list,
    de_results = de_results,
    updated_net = updated_net,
    plot_combined = p_combined
  ))
}

# visualize distribution of ligand/receptor expression
# Assuming 'net' is your annotated network dataframe with the new columns
plot_pct_expressing_distribution <- function(net) {
  net_ligand <- net %>%
    select(interaction_name, ligand.pct_expressing) %>%
    filter(!is.na(ligand.pct_expressing)) %>%
    mutate(GeneType = "Ligand", PctExpressing = ligand.pct_expressing)
  
  net_receptor <- net %>%
    select(interaction_name, receptor.pct_expressing) %>%
    filter(!is.na(receptor.pct_expressing)) %>%
    mutate(GeneType = "Receptor", PctExpressing = receptor.pct_expressing)
  
  net_long <- bind_rows(net_ligand, net_receptor)
  
  ggplot(net_long, aes(x = PctExpressing, fill = GeneType)) +
    geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 5, color = "red", linetype = "dashed") +
    labs(
      title = "Distribution of % Expressing Cells (Pseudobulk)",
      subtitle = "Ligands and Receptors across interactions in 'net'",
      x = "% of Cells Expressing",
      y = "Count of Interactions",
      fill = "Gene Type"
    ) +
    theme_minimal(base_size = 12)
}

# revised plotting function
plot_interactions_delta <- function(df, pathway, major, axis_limits, label_LR = TRUE, dataset_name = "Full", delta_thresh = 0.0001) {
  df_pathway <- df %>%
    filter(pathway_name == pathway) %>%
    mutate(
      highlight = ifelse(abs(delta_prob) > delta_thresh & abs(combined_logFC) > 0.15, "High Activity", "Other"),
      source_color = ifelse(highlight == "High Activity", major[source], "gray70"),
      target_color = ifelse(highlight == "High Activity", major[target], "gray90")
    )
  
  p <- ggplot(df_pathway, aes(x = combined_logFC, y = delta_prob)) +
    geom_point(
      aes(fill = source_color, color = target_color),
      shape = 21,
      size = ifelse(df_pathway$highlight == "High Activity", 1.5, 1),
      stroke = ifelse(df_pathway$highlight == "High Activity", 0.6, 0.2),
      alpha = ifelse(df_pathway$highlight == "High Activity", 0.9, 0.5)
    ) +
    geom_hline(yintercept = c(-delta_thresh, delta_thresh), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-0.15, 0.15), linetype = "dotted", color = "gray50") +
    scale_fill_identity() +
    scale_color_identity() +
    labs(
      title = paste0(dataset_name, " - ", pathway),
      subtitle = paste0("AD-diff Prob (AD - NCI), threshold = ±", delta_thresh),
      x = "Combined logFC (AD - NCI)",
      y = "AD-diff Interaction Probability"
    ) +
    coord_cartesian(
      xlim = axis_limits$x,
      ylim = axis_limits$y
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
  
  if (label_LR) {
    p <- p + geom_text_repel(
      data = df_pathway %>% filter(highlight == "High Activity"),
      aes(label = interaction_name),
      size = 2.5,
      max.overlaps = 25,
      box.padding = 0.4
    )
  }
  
  return(p)
}

plot_pathways_pdf <- function(pathways, major, label_LR = FALSE, output_file = "interaction_grid.pdf") {
  nets <- list(
    Full = net_diff,
    Niche1 = net1_diff,
    Niche2 = net2_diff,
    Niche3 = net3_diff
  )
  
  all_plots <- list()
  
  for (p in pathways) {
    # Calculate delta prob threshold from FULL dataset only
    full_data <- net_diff %>% filter(pathway_name == p)
    delta_thresh <- quantile(full_data$delta_prob, 0.95, na.rm = TRUE)
    
    # Axis limits
    all_data <- bind_rows(lapply(nets, \(df) df %>% filter(pathway_name == p)))
    axis_limits <- list(
      x = range(all_data$combined_logFC, na.rm = TRUE) %>% scales::squish_infinite(),
      y = range(all_data$delta_prob, na.rm = TRUE) %>% scales::squish_infinite()
    )
    
    for (n in names(nets)) {
      p_plot <- plot_interactions_delta(
        df = nets[[n]], 
        pathway = p, 
        major = major, 
        delta_thresh = delta_thresh,
        axis_limits = axis_limits,
        label_LR = label_LR,
        dataset_name = n
      )
      all_plots[[paste(n, p)]] <- p_plot
    }
  }
  
  pdf(output_file, width = 12, height = 12)
  grid.arrange(grobs = all_plots, ncol = 4)
  dev.off()
}



# relative prob version
plot_relative_interactions <- function(df, pathway, major, axis_limits, label_LR = TRUE, dataset_name = "Full", rel_thresh = 0.25) {
  df_pathway <- df %>%
    filter(pathway_name == pathway) %>%
    mutate(
      relative_delta = delta_prob / max(abs(delta_prob), na.rm = TRUE),
      highlight = ifelse(abs(relative_delta) > rel_thresh & abs(combined_logFC) > 0.15, "High Activity", "Other"),
      source_color = ifelse(highlight == "High Activity", major[source], "gray70"),
      target_color = ifelse(highlight == "High Activity", major[target], "gray90"),
      point_size = ifelse(highlight == "High Activity", 1.5, 0.5),
      point_stroke = ifelse(highlight == "High Activity", 0.6, 0.2),
      point_alpha = ifelse(highlight == "High Activity", 0.9, 0.5)
    )
  
  p <- ggplot(df_pathway, aes(x = combined_logFC, y = relative_delta)) +
    geom_point(
      aes(fill = source_color, color = target_color, size = point_size, stroke = point_stroke, alpha = point_alpha),
      shape = 21
    ) +
    geom_hline(yintercept = c(-rel_thresh, rel_thresh), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-0.15, 0.15), linetype = "dotted", color = "gray50") +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_identity() +     
    scale_alpha_identity() + 
    labs(
      title = paste0(dataset_name, " - ", pathway),
      subtitle = paste0("Relative AD-diff Prob (AD - NCI), threshold = ±", rel_thresh),
      x = "Combined logFC (AD - NCI)",
      y = "Relative AD-diff Interaction Probability"
    ) +
    coord_cartesian(xlim = axis_limits$x, ylim = c(-1, 1)) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
  
  if (label_LR) {
    p <- p + geom_text_repel(
      data = df_pathway %>% filter(highlight == "High Activity"),
      aes(label = interaction_name),
      size = 2.5,
      max.overlaps = 25,
      box.padding = 0.4
    )
  }
  
  return(p)
}

plot_rel_pathways_pdf <- function(pathways, major, label_LR = FALSE, output_file = "interaction_grid_relative.pdf") {
  nets <- list(
    Full = net_diff,
    Niche1 = net1_diff,
    Niche2 = net2_diff,
    Niche3 = net3_diff
  )
  
  all_plots <- list()
  
  for (p in pathways) {
    all_data <- bind_rows(lapply(nets, \(df) df %>% filter(pathway_name == p)))
    axis_limits <- list(
      x = range(all_data$combined_logFC, na.rm = TRUE) %>% scales::squish_infinite()
    )
    
    for (n in names(nets)) {
      p_plot <- plot_relative_interactions(
        df = nets[[n]],
        pathway = p,
        major = major,
        axis_limits = axis_limits,
        label_LR = label_LR,
        dataset_name = n,
        rel_thresh = 0.25
      )
      all_plots[[paste(n, p)]] <- p_plot
    }
  }
  
  pdf(output_file, width = 12, height = 12)
  gridExtra::grid.arrange(grobs = all_plots, ncol = 4)
  dev.off()
}


