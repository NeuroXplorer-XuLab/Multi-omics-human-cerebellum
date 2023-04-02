library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(metap)
library(UpSetR)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(pheatmap)
library(plotly)
library(enrichplot)
library(DOSE)
library(gg.gap)
library(viridis)
library(CellChat)
library(ggalluvial)
library(proxyC)
library(monocle3)
library(topGO)
library(loomR)
library(SCENIC)
library(SeuratDisk)
library(AUCell)
library(GENIE3)
library(RcisTarget)
library(SCopeLoomR)
library(scales)
library(StemSC)
library(CytoTRACE)
library(slingshot)
library(grDevices)
library(destiny)
library(tradeSeq)
library(ggformula)
rm(list = ls())
rm(list = ls())
setwd("E:/cerebellum/result")

mockdim <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  panel.border = element_blank()
) & NoLegend()

trydot <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"),
                    col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                    idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE,
                    scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA,
                    idencolor) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = RColorBrewer::brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", 
              call. = FALSE, immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", 
                                         color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", 
                     no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  
  clusinfo <- unlist(lapply(strsplit(as.character(data.plot$id), split = '_'), function(x) x[1]))
  clusinfo <- as.numeric(clusinfo)
  sampleinfo <- unlist(lapply(strsplit(as.character(data.plot$id), split = '_'), function(x) x[2]))
  Tcolor <- idencolor[(clusinfo+1)]
  lineinfo <- paste(clusinfo, data.plot$features.plot, sep = '_')
  data.plot <- cbind(data.plot, clusinfo, sampleinfo, Tcolor, lineinfo)
  
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "sampleinfo", 
                                                        y = "avg.exp")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = 'Tcolor', shape = 'features.plot')) + 
    geom_line(mapping = aes_string(group = 'lineinfo', color = 'Tcolor')) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Sequential Samples", y = ifelse(test = is.null(x = split.by), 
                                              yes = "Identity", no = "Average Expression")) + 
    cowplot::theme_cowplot()
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  return(plot)
}

ccosg <- function(object, groups = "all", assay = "RNA", slot = "data", mu = 1, n_genes_user = 100){
  genexcell <- Seurat::GetAssayData(object = object[[assay]], 
                                    slot = slot)
  if (groups == "all") {
    group_info <- Seurat::Idents(object = object)
  }
  else {
    group_info <- factor(object[[groups]][, 1])
    names(group_info) <- colnames(object)
  }
  groups_order = sort(unique(group_info))
  n_cluster = length(groups_order)
  if (n_cluster == 1) {
    stop("Cannot perform marker gene identification on a single cluster.")
  }
  n_cell = ncol(genexcell)
  n_gene = nrow(genexcell)
  gene_name = rownames(genexcell)
  if (n_genes_user > n_gene) {
    n_genes_user = n_gene
  }
  cluster_mat = matrix(0, nrow = n_cluster, ncol = n_cell)
  order_i = 1
  for (group_i in groups_order) {
    idx_i = group_info == group_i
    cluster_mat[order_i, idx_i] = 1
    order_i = order_i + 1
  }
  cluster_mat_sparse = as(cluster_mat, "dgCMatrix")
  cosine_sim = proxyC::simil(genexcell, cluster_mat_sparse, 
                             method = "cosine", drop0 = TRUE)
  pos_nonzero = cosine_sim != 0
  pos_nonzero = which(as.matrix(pos_nonzero), arr.ind = TRUE)
  genexlambda = cosine_sim * cosine_sim
  e_power2_sum = Matrix::rowSums(genexlambda)
  if (mu == 1) {
    genexlambda[pos_nonzero] = genexlambda[pos_nonzero]/(replicate(ncol(genexlambda), 
                                                                   e_power2_sum)[as.matrix(pos_nonzero)])
  }
  else {
    genexlambda[pos_nonzero] = genexlambda[pos_nonzero]/(((1 - 
                                                             mu) * genexlambda[pos_nonzero] + mu * (replicate(ncol(genexlambda), 
                                                                                                              e_power2_sum)))[as.matrix(pos_nonzero)])
  }
  genexlambda = genexlambda * cosine_sim
  rank_stats_names = data.frame(matrix(matrix(), n_genes_user, 
                                       length(groups_order), dimnames = list(seq(1, n_genes_user), 
                                                                             groups_order)), stringsAsFactors = F)
  rank_stats_scores = data.frame(matrix(matrix(), n_genes_user, 
                                        length(groups_order), dimnames = list(seq(1, n_genes_user), 
                                                                              groups_order)), stringsAsFactors = F)
  order_i = 1
  for (group_i in groups_order) {
    idx_i = group_info == group_i
    scores = genexlambda[, order_i]
    global_indices = COSG:::select_top_n(scores, n_genes_user)
    rank_stats_names[, order_i] = gene_name[global_indices]
    rank_stats_scores[, order_i] = scores[global_indices]
    order_i = order_i + 1
  }
  colnames(rank_stats_names) <- groups_order
  colnames(rank_stats_scores) <- groups_order
  ranks_stats = list(names = rank_stats_names, scores = rank_stats_scores)
  return(ranks_stats)
}

modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

sDimplot <- function(object, color, group, bgcolor = '#eeeeee'){
  spx <- object[['x']]
  spy <- object[['y']]
  cluster <- object[[group]]
  mforsplot <- cbind(as.data.frame(cbind(spx, spy)), cluster)
  if('' %in% color){
    color = hue_pal()(length(table(spgw13$seurat_clusters)))
  }
  fmap <- ggplot(mforsplot, aes(x = mforsplot$x, y = mforsplot$y, color = mforsplot[,3])) + 
    geom_point(shape=19) + scale_color_manual(values = color) + 
    xlim(min(spx), max(spx)) + ylim(min(spy), max(spy)) + 
    labs(x = "", y = "", title = "") +
    theme_bw() + theme(panel.background = element_rect(fill = bgcolor)) +
    theme(panel.grid =element_blank()) +
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank())+
    theme(panel.border = element_blank())
  fmap
  return(fmap)
} 
sDimplot(spgw13, color = '', group = 'seurat_clusters')

sFeatureplot <- function(object, feature, slot = 'data', color = c('grey', 'navy'), bgcolor = '#eeeeee'){
  spx <- object[['x']]
  spy <- object[['y']]
  mfeature <- FetchData(object, feature, slot = slot)
  mforsplot <- cbind(as.data.frame(cbind(spx, spy)), mfeature)
  fmap <- ggplot(mforsplot, aes(x = mforsplot$x, y = mforsplot$y, color = mforsplot[,3])) + 
    geom_point(shape=19) + scale_color_gradientn(colours = color) + 
    xlim(min(spx), max(spx)) + ylim(min(spy), max(spy)) + 
    labs(x = "", y = "", title = feature) +
    theme_bw() + theme(panel.background = element_rect(fill = bgcolor)) +
    theme(panel.grid =element_blank()) +
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank())+
    theme(panel.border = element_blank())
  fmap
  return(fmap)
} 


