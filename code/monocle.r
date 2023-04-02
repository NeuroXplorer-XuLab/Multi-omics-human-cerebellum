load(file = "E:/cerebellum/result/human UBC two lineage.Rdata")
load(file = "E:/cerebellum/result/cebclean.Rdata")
load(file = 'E:/cerebellum/result/human pkc+pro+in.Rdata')
load(file = 'E:/cerebellum/result/human only GC.Rdata')
DefaultAssay(ubwithsimi) <- 'RNA'
DefaultAssay(cebclean) <- 'RNA'
DefaultAssay(ppicell) <- 'RNA'

VlnPlot(subset(ubwithsimi, seurat_clusters %in% c(12,0,4,1,2,5,8,9,6,7,14)), c('GRM1', 'CALB2'), group.by = 'celltype')

ubcfortraj <- subset(ubwithsimi, seurat_clusters %in% setdiff(0:18, c(11,10,18,13)))
data <- GetAssayData(ubcfortraj, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(ubcfortraj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
# 
cds <- cluster_cells(cds, k = 15)
p3 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p4 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p3 + p4
# 
cds <- learn_graph(cds, close_loop = F)
p6 <- plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE)
p6
cds <- order_cells(cds)
p7 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, 
                 label_leaves = F, label_branch_points = F, label_roots = F)+ 
  scale_color_gradientn(colors = c('#63669d', '#a9d6af'))
p7
# 
save(cds, file = 'human UBC monocle.Rdata')
# 
pdf(file = 'hsa-ceb UBC traj.pdf', height = 3, width = 3)
print(p7)
print(p7 & mockdim & labs(title = element_blank()))
dev.off()
# 
#
load(file = 'human UBC monocle.Rdata')
rrtoii <- graph_test(cds, neighbor_graph = "principal_graph")
changgg <- rrtoii %>% filter(q_value < 0.01)
changgg <- changgg[order(changgg$morans_I, decreasing = T), ]
# tempoubc <- c('PIK3C2G', 'MEGF10', 'PTPRT', 'ATP1A2', 'SLC1A3', 'GDF10', 'MFGE8', 
#                'ELN', 'WLS', 'CALCB', 'CRYAB', 'SULF1', 'CITED4', 'CALCA', 'MSX1', 
#                'TRPM3', 'TMTC2', 'EYS', 'GRIK1', 'LRP1B', 'PTPRM', 'GALNTL6', 
#                'MASP1', 'PLP1', 'ZFP36L1', 'NCKAP5', 'SPARC', 'SMOC1', 'RFX4', 
#                'STMN2', 'CCSER1', 'CACNA2D1')
# pseudot <- cds@principal_graph_aux$UMAP$pseudotime
# pseudot <- pseudot[order(pseudot, decreasing = F)]
# dtmatrix <- ubcfortraj@assays$RNA@data[tempoubc, names(pseudot)]
# dtmatrix <- as.data.frame(dtmatrix)
# # 
# pheatmap(dtmatrix, cluster_rows = F, cluster_cols = F, border = F, 
#          color = colorRampPalette(c("#f0eeec", '#b8d4e9', '#63669d'))(100),
#          show_rownames = T, show_colnames = F, fontsize = 250/length(tempoubc))
tempoubc <- c('PIK3C2G', 'CRYAB', 'TRPM3', 'EOMES', 'TMTC2', 'LRP1B')
pdf(file = 'hsa-ceb UBC traj change gene.pdf', height = 4, width = 9)
fprint <- plot_genes_in_pseudotime(cds[tempoubc, ], cell_size = 1, ncol = 3, 
                                   panel_order = tempoubc) + 
  scale_y_continuous(trans='log1p') + 
  scale_color_gradientn(colors = c('#63669d', '#a9d6af')) + theme_bw() & 
  theme(axis.title = element_blank()) + NoLegend()
print(fprint)
dev.off()


regumatseu <- regumatrix[colnames(ubcfortraj), ]
# ubcfortraj$regu_BCL11B <- regumatseu$BCL11B
ubcfortraj$regu_ATOH1 <- regumatseu$ATOH1
ubcfortraj$regu_TCF7L2 <- regumatseu$TCF7L2
ubcfortraj$regu_MYBL2 <- regumatseu$MYBL2
# 
ubrrrrr <- c('regu_BCL11B', 'regu_ATOH1', 'regu_TCF7L2', 'regu_MYBL2')
for(i in ubrrrrr){
  png(filename = paste0('UBC ', i, '.png'), width = 400, height = 400)
  f1 <- FeaturePlot(ubcfortraj, i, cols = c('#dddddd', '#4b5861'), max.cutoff = 0.3) & mockdim & 
    labs(title = element_blank())
  print(f1)
  dev.off()
}