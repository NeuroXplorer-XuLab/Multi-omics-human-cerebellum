load(file = 'human pkc+pro+in.Rdata')
VlnPlot(ppicell, 'nFeature_RNA', pt.size = 0)
pim17 <- FindMarkers(ppicell, ident.1 = 17, min.diff.pct = 0.3, logfc.threshold = 0.7)
pim16 <- FindMarkers(ppicell, ident.1 = 16, min.diff.pct = 0.3, logfc.threshold = 0.7)
pim11 <- FindMarkers(ppicell, ident.1 = 11, min.diff.pct = 0.3, logfc.threshold = 0.7)
pim1 <- FindMarkers(ppicell, ident.1 = 1, min.diff.pct = 0.3, logfc.threshold = 0.7)
pim56 <- FindMarkers(ppicell, ident.1 = c(5,6), 
                     min.diff.pct = 0.3, logfc.threshold = 0.5, only.pos = T)
pim5anti6 <- FindMarkers(ppicell, ident.1 = 5, ident.2 = 6)
pimpur <- FindMarkers(ppicell, ident.1 = c(2,3,4,7), min.diff.pct = 0.3, logfc.threshold = 0.7)
# 
ppicell <- subset(ppicell, seurat_clusters %in% setdiff(0:17, 16))
fppi <- DimPlot(ppicell, label = T, label.size = 5)
fppi
ppicell$celltype <- ppicell$seurat_clusters
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(1))] <- 'BG'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(8,10))] <- 'G2M'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(11))] <- 'GCP'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(5,6))] <- 'VZP'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(12))] <- 'EPP'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(15,2,3,4,7))] <- 'PKC'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(0,13,14,9))] <- 'IN'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(17))] <- 'OPC'
table(ppicell$celltype)


fintp <- DimPlot(subset(ppicell, celltype != 'GCP'), group.by = 'celltype', shuffle = T, 
                cols = c('#63669d', '#b7c2c8', '#916e8c', '#83949e', '#bab8a9', 
                         '#fee2d4', '#eac2e4')) & 
  mockdim & labs(title = element_blank())
fintp
pdf('hsa-ceb ppi umap.pdf', width = 3, height = 3)
print(fintp)
dev.off()
#
top2appi <- subset(ppicell, seurat_clusters %in% c(8,10))
DefaultAssay(top2appi) <- 'integrated'
top2appi <- ScaleData(top2appi)
top2appi <- RunPCA(top2appi, npcs = 50)
ElbowPlot(object = top2appi, ndims = 40, reduction = "pca") 
DefaultAssay(top2appi) <- 'integrated'
top2appi <- FindNeighbors(top2appi, reduction = "pca", dims = 1:11)
top2appi <- FindClusters(top2appi, resolution = 1)
top2appi <- RunUMAP(top2appi, reduction = "pca", dims = 1:11)
ftop <- DimPlot(top2appi, label = T, label.size = 5)
ftop
DefaultAssay(top2appi) <- 'RNA'
VlnPlot(top2appi, c('RGS16', 'GADD45G', 'PTF1A', 'KIRREL2', 'NRN1', 
                    'THSD7B', 'LMX1A', 'TRPM3', 'PIK3C2G'))
VlnPlot(top2appi, c('TOP2A', 'SOX2'))
tp2 <- FindMarkers(top2appi, ident.1 = 2, logfc.threshold = 0.7)
tp4 <- FindMarkers(top2appi, ident.1 = 4, logfc.threshold = 0.5)

ppiforvln <- subset(ppicell, seurat_clusters %in% c(1,5,6,12,8,10,17))
table(ppiforvln$celltype)
# 
progdifgene <- c('NPY', 'TNC', 'LINC01727', 'FST', 'MT2A', 'PTF1A', 'DHRS3', 'KIRREL2', 
                 'TOP2A', 'MKI67', 'AQP4', 'PIFO', 'RSPH1', 
                 'OLIG2', 'PDGFRA')
pdf(file = 'hsa-ceb prog violin.pdf', height = 15, width = 6)
fprint <- StackedVlnPlot(ppiforvln, features = progdifgene, 
                         group.by = 'celltype')
print(fprint)
dev.off()
# 
# 5/6 anti GO
DefaultAssay(ppicell) <- 'RNA'
mvzrgc111 <- FindMarkers(ppicell, ident.1 = 5, ident.2 = 6, 
                         min.diff.pct = 0.1, logfc.threshold = 0.3)
glist <- bitr(mvzrgc111 %>% filter(avg_log2FC<0) %>% rownames(), 
              fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gogogo <- enrichGO(glist$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', readable = T)
gogogo<- pairwise_termsim(gogogo)
clusterProfiler::dotplot(gogogo, showCategory = 10)+ 
        scale_color_gradientn(colors = c('#80325c', '#edeaf1'))


ppifortraj <- subset(ppicell, seurat_clusters %in% c(5,6,9,13,14,0))
data <- GetAssayData(ppifortraj, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(ppifortraj, reduction = "umap")
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
  scale_color_gradientn(colors = c('#68926a', '#b8d4e9'))
p7
# 
pdf(file = 'hsa-ceb prog traj to interneuron.pdf', height = 3, width = 3)
print(p7)
print(p7 & mockdim & labs(title = element_blank()))
dev.off()
#
save(cds, file = 'human pro monocle.Rdata')
# 
rrtoii <- graph_test(cds, neighbor_graph = "principal_graph")
changgg <- rrtoii %>% filter(q_value < 0.01)
changgg <- changgg[order(changgg$morans_I, decreasing = T), ]
# 
load(file = 'E:/cerebellum/result/human pro monocle.Rdata')
temporalg <- c('TTYH1', 'HES1', 'SMOC1', 'GAD2', 'NRXN3', 'XKR4', 'SLN', 'LRRTM4', 'C1QTNF3', 
               'HES6', 'TFPI', 'RGS16')
pdf(file = 'hsa-ceb RGC-IPC-inte traj change gene.pdf', height = 8, width = 6)
fprint <- plot_genes_in_pseudotime(cds[temporalg, ], cell_size = 1, ncol = 3, 
                                   panel_order = c('C1QTNF3', 'LRRTM4', 'SLN', 'TTYH1', 'HES1', 'SMOC1', 
                                                   'TFPI', 'HES6', 'RGS16', 'GAD2', 'NRXN3', 'XKR4')) + 
  scale_y_continuous(trans='log1p') + 
  scale_color_gradientn(colors = c('#68926a', '#b8d4e9')) + theme_bw() & 
  theme(axis.title = element_blank()) + NoLegend()
print(fprint)
dev.off()
#

test <- subset(ppicell, seurat_clusters %in% c(5,1,6,9))
DefaultAssay(test) <- 'RNA'
for(i in names(table(test$orig.ident))){
  assign(i, subset(test, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
rm(test); gc()
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30, k.filter = 140)
gc()
test <- IntegrateData(anchorset = anchor, dims = 1:30, k.weight = 140)
gc()
test <- ScaleData(test)
test <- RunPCA(test, npcs = 50)
ElbowPlot(object = test, ndims = 40, reduction = "pca") 
DefaultAssay(test) <- 'integrated'
test <- FindNeighbors(test, reduction = "pca", dims = 1:19)
test <- FindClusters(test, resolution = 1)
test <- RunUMAP(test, reduction = "pca", dims = 1:19)
ftest <- DimPlot(test, label = T, label.size = 5)
ftest
DefaultAssay(test) <- 'RNA'
# 
mtest <- FindMarkers(test, ident.1 = 9, logfc.threshold = 0.7, min.diff.pct = 0.2)
# 
# 
# 
DefaultAssay(test) <- 'RNA'
gcforgsea <- FindMarkers(test, ident.1 = c(2), ident.2 = c(0,9),
                         logfc.threshold = 0, min.diff.pct = 0)
gcforgsea <- cbind(gcforgsea, rownames(gcforgsea))
colnames(gcforgsea)[6] <- 'SYMBOL'
gcid <- bitr(rownames(gcforgsea), OrgDb = 'org.Hs.eg.db', 
             fromType = 'SYMBOL', toType = c('ENTREZID', 'SYMBOL'))
gcgsea <- merge(gcforgsea, gcid, by = 'SYMBOL')
gcgsea <- gcgsea[order(gcgsea$avg_log2FC, decreasing = T), ]
ggcclist <- gcgsea$avg_log2FC
names(ggcclist) <- gcgsea$ENTREZID
# 
gcgsego <- gseGO(geneList = ggcclist, OrgDb = org.Hs.eg.db, verbose = T)
options(clusterProfiler.download.method = "wininet")
gcgsekegg <- gseKEGG(geneList = ggcclist, organism = "hsa")
# 

test2 <- subset(test, seurat_clusters %in% c(5,1,6,9))
test <- ScaleData(test)
test <- RunPCA(test, npcs = 50)
ElbowPlot(object = test, ndims = 40, reduction = "pca") 
DefaultAssay(test) <- 'integrated'
test <- FindNeighbors(test, reduction = "pca", dims = 1:19)
test <- FindClusters(test, resolution = 1)
test <- RunUMAP(test, reduction = "pca", dims = 1:19)
ftest <- DimPlot(test, label = T, label.size = 5)
ftest



ppifortraj
DefaultAssay(ppifortraj) <- 'integrated'
ppifortraj <- FindNeighbors(ppifortraj, reduction = "pca", dims = 1:8)
ppifortraj <- FindClusters(ppifortraj, resolution = 0.5)
DimPlot(ppifortraj, label = T)
VlnPlot(ppifortraj, 'nFeature_RNA')
DefaultAssay(ppifortraj) <- 'RNA'
inhi6 <- FindMarkers(ppifortraj, ident.1 = 6, min.diff.pct = 0.2, logfc.threshold = 0.5)
subtraj <- subset(ppifortraj, seurat_clusters %in% c(0,4,7))
DefaultAssay(subtraj) <- 'RNA'
DimPlot(subtraj, label = T, label.size = 5)
mrgc2 <- FindMarkers(subtraj, ident.1 = 4, ident.2 = 0, logfc.threshold = 0.3, min.diff.pct = 0.2)
# 
subtraj$celltype <- subtraj$seurat_clusters
levels(subtraj$celltype)[which(levels(subtraj$celltype) %in% c(0))] <- 'RGC-1'
levels(subtraj$celltype)[which(levels(subtraj$celltype) %in% c(4))] <- 'RGC-2'
levels(subtraj$celltype)[which(levels(subtraj$celltype) %in% c(7))] <- 'IPC'
table(subtraj$celltype)
# 
proggene <- c('PTF1A', 'SLN', 'CDH4', 'C1QTNF3', 'ESRRG', 'FOXP2', 'GAS2L3', 
              'LUZP2', 'SPARCL1', 'CADM2', 'KCND2', 'TFPI', 'TSC22D4', 
              'HES6', 'RGS16')
# , 'AL137139.2', 'DACH1','GDF10', 'AL133346.1', 'BAALC', 'MSX2',, 'ZIC2' 
pdf(file = 'hsa-ceb RGC-IPC-inte progenitor diff gene.pdf', height = 4, width = 4)
print(DotPlot(subtraj, features = proggene, group.by = 'celltype', 
              cluster.idents = F, cols = c('#dddddd', '#006934')) + RotatedAxis() + coord_flip())
dev.off()
# 
# 
rgipip <- aggregate(subtraj@meta.data, by = list(subtraj@meta.data$timepoint, subtraj@meta.data$seurat_clusters), FUN=table)
write.csv(file = 'hsa-ceb prog RGC IPC ratio.csv', rgipip[, 1:3], quote = F)
# 
pdf(file = 'hsa-ceb prog umap.pdf', height = 3, width = 3)
fprint <- DimPlot(subtraj, group.by = 'seurat_clusters', shuffle = T, 
                 cols = c('#68926a', '#84b186', '#a9d6af')) & 
  mockdim & labs(title = element_blank())
print(fprint)
dev.off()
# 
mrgc2 <- FindMarkers(subtraj, ident.1 = 4, ident.2 = 0, logfc.threshold = 0.3, min.diff.pct = 0.1, only.pos = T)
mrgc1 <- FindMarkers(subtraj, ident.1 = 0, ident.2 = 4, logfc.threshold = 0.3, min.diff.pct = 0.2, only.pos = T)
#
rgcgo1 <- bitr(rownames(mrgc1), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
rgcgo1 <- enrichGO(rgcgo1$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', readable = T)
clusterProfiler::dotplot(rgcgo1, showCategory = 10)
rgcgo2 <- bitr(rownames(mrgc2), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
rgcgo2 <- enrichGO(rgcgo2$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', readable = T)
clusterProfiler::dotplot(rgcgo2, showCategory = 10)








# 
# 
#
# 
interneuron <- subset(ppicell, seurat_clusters %in% c(9,13,0,14))
#
DefaultAssay(interneuron) <- 'RNA'
for(i in names(table(interneuron$orig.ident))){
  assign(i, subset(interneuron, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20, k.filter = 49)
gc()
interneuron <- IntegrateData(anchorset = anchor, dims = 1:20, k.weight = 49)
#
DefaultAssay(interneuron) <- 'integrated'
interneuron <- ScaleData(interneuron)
interneuron <- RunPCA(interneuron, npcs = 50)
ElbowPlot(object = interneuron, ndims = 40, reduction = "pca") 
DefaultAssay(interneuron) <- 'integrated'
interneuron <- FindNeighbors(interneuron, reduction = "pca", dims = 1:9)
interneuron <- FindClusters(interneuron, resolution = 0.8)
interneuron <- RunUMAP(interneuron, reduction = "pca", dims = 1:9)
finhi <- DimPlot(interneuron, label = T, label.size = 5)
finhi
DefaultAssay(interneuron) <- 'RNA'
VlnPlot(interneuron, 'nFeature_RNA')
interneuron <- subset(interneuron, seurat_clusters != 10)
DefaultAssay(interneuron) <- 'integrated'
interneuron <- ScaleData(interneuron)
interneuron <- RunPCA(interneuron, npcs = 50)
ElbowPlot(object = interneuron, ndims = 40, reduction = "pca") 
DefaultAssay(interneuron) <- 'integrated'
interneuron <- FindNeighbors(interneuron, reduction = "pca", dims = 1:8)
interneuron <- FindClusters(interneuron, resolution = 1)
interneuron <- RunUMAP(interneuron, reduction = "pca", dims = 1:8)
finhi <- DimPlot(interneuron, label = T, label.size = 5)
finhi
DefaultAssay(interneuron) <- 'RNA'
# 
# 
VlnPlot(interneuron, c('MEG3', 'ZIC2', 'IGFBP5', 'PANTR1'))
interneuron$celltype <- interneuron$seurat_clusters
levels(interneuron$celltype)[which(levels(interneuron$celltype) %in% c(1,2,9,6))] <- 'LU'
levels(interneuron$celltype)[which(levels(interneuron$celltype) %in% c(0,3,4,5,7,8,10))] <- 'GO'
table(interneuron$celltype)
luggolgi <- aggregate(interneuron@meta.data, by = list(interneuron@meta.data$timepoint, interneuron@meta.data$celltype), FUN=table)
write.csv(file = 'hsa-ceb inter lugaro-golgi ratio.csv', luggolgi[, 1:3], quote = F)
# 
# 
pdf(file = 'hsa-ceb inter RGS16.pdf', height = 3, width = 3)
f1 <- FeaturePlot(interneuron, 'RGS16', cols = c('#dddddd', 'navy')) & mockdim & 
  labs(title = element_blank())
print(f1)
dev.off()
# 
VlnPlot(interneuron, c('MEG3', 'ZIC2'))
f1 <- FeaturePlot(interneuron, c('MEG3', 'ZIC2'), blend = T, blend.threshold = 0, 
            cols = c('#eeeeee', '#68926a', '#80325c')) & mockdim & labs(title = element_blank())
f2 <- FeaturePlot(interneuron, c('PANTR1', 'IGFBP5'), blend = T, blend.threshold = 0,
            cols = c('#eeeeee', '#68926a', '#80325c')) & mockdim & labs(title = element_blank())
pdf(file = 'hsa-ceb inter counter pair.pdf', height = 3, width = 12)
print(f1)
print(f2)
dev.off()
# fprint <- DimPlot(interneuron, label = F, group.by = 'timepoint', shuffle = T,
#                   cols = c('#68926a', '#a9d6af', '#b8d4e9', '#8eb0c9', '#916e8c', '#80325c')) &
#   mockdim & labs(title = element_blank())
# fprint
# png('hsa-ceb main diff time.png', height = 600, width = 600)
# 
# print(fprint)
# dev.off()
# 
# 
purkinje <- subset(ppicell, seurat_clusters %in% c(2,3,4,7,15))
DefaultAssay(purkinje) <- 'integrated'
purkinje <- ScaleData(purkinje)
purkinje <- RunPCA(purkinje, npcs = 50)
ElbowPlot(object = purkinje, ndims = 40, reduction = "pca") 
DefaultAssay(purkinje) <- 'integrated'
purkinje <- FindNeighbors(purkinje, reduction = "pca", dims = 1:14)
purkinje <- FindClusters(purkinje, resolution = 1)
purkinje <- RunUMAP(purkinje, reduction = "pca", dims = 1:14)
fpur <- DimPlot(purkinje, label = T, label.size = 5)
fpur
DefaultAssay(purkinje) <- 'RNA'
VlnPlot(purkinje, 'nFeature_RNA')
FeaturePlot(purkinje, 'LSAMP')
purkinje <- subset(purkinje, seurat_clusters != 9)
# 
purkinje$celltype <- purkinje$seurat_clusters
levels(purkinje$celltype)[which(levels(purkinje$celltype) %in% c(8))] <- 'p1'
levels(purkinje$celltype)[which(levels(purkinje$celltype) %in% c(1,7))] <- 'p2'
levels(purkinje$celltype)[which(levels(purkinje$celltype) %in% c(0,2,3,4,5,6))] <- 'p3'
table(purkinje$celltype)
# 
save(purkinje, file = 'human pkc.Rdata')
# 
# 
load(file = 'human pkc.Rdata')
# 
pkctype <- aggregate(purkinje@meta.data, by = list(purkinje@meta.data$timepoint, purkinje@meta.data$celltype), FUN=table)
write.csv(file = 'hsa-ceb pkc type ratio.csv', pkctype[, 1:3], quote = F)
# 
pp1 <- FindMarkers(purkinje, ident.1 = c(8), min.diff.pct = 0.2, logfc.threshold = 0.5)
pp2 <- FindMarkers(purkinje, ident.1 = c(1,7), min.diff.pct = 0.2, logfc.threshold = 0.5)
pp3 <- FindMarkers(purkinje, ident.1 = c(0,2,3,4,5,6))
pur2anti3 <- FindMarkers(purkinje, ident.1 = 'p2', ident.2 = 'p3', 
                         group.by = 'celltype', logfc.threshold = 0.5)
pdf(file = 'hsa-ceb pkc test p2p3.pdf', height = 4, width = 4)
for(i in rownames(pur2anti3)){
  print(FeaturePlot(purkinje, i)& mockdim)
}
dev.off()
# # 
pkcgenedot <- c('PCDH7', 'UNC5D', 'EBF2', 'IGF1', 'LINC02336', 'SCG2', 
                'LSAMP', 'NLGN1', 'DLGAP2', 'CACNA1C', 'AC092691.1', 
                'PBX3', 'LRRTM4', 'SLA', 'LMO4', 'KCTD16', 'FGF3')
pdf(file = 'hsa-ceb pkc marker dotplot.pdf', height = 5, width = 4.5)
DotPlot(purkinje, features = pkcgenedot, group.by = 'celltype', 
        cluster.idents = F, cols = c('#dddddd', '#49354e')) + RotatedAxis() + coord_flip()
dev.off()
# pkc的umap
pdf(file = 'hsa-ceb pkc umap.pdf', height = 2, width = 2)
fprint <- DimPlot(purkinje, group.by = 'celltype', shuffle = T, 
                  cols = c('#916e8c', '#80325c', '#c2b9ee')) & 
  mockdim & labs(title = element_blank())
print(fprint)
dev.off()
#
# pkc spatial
pdf(file = 'hsa-ceb pkc feature PLCB4 CADM1.pdf', height = 2, width = 4)
fprint <- FeaturePlot(purkinje, c('PLCB4', 'CADM1'), cols = c('#dedede', 'navy')) & 
  mockdim & labs(title = element_blank())
print(fprint)
dev.off()


#
# 
load(file = 'human pkc.Rdata')
# 
pkccoso <- ccosg(purkinje, groups = 'celltype', 
                assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 20)
pkccoso <- pkccoso[[1]]
pkccos <- c(pkccoso$p1, pkccoso$p2, pkccoso$p3)
pkccos <- pkccos[!grepl(x = pkccos, pattern = '^RPS|^RPL')]
grninter <- read.csv(file = "E:/cerebellum/purkinje forced.csv")

write.csv(pkccoso, file = 'purkinje DEG.csv')
pkpkpk <- grninter %>% filter(Downstream_Gene %in% pkccos)
pkpkpk <- pkpkpk[, c("Upstream_TF", "Downstream_Gene")]
pkpkpk <- pkpkpk[!duplicated(pkpkpk), ]
upupname <- pkpkpk$Upstream_TF
selename <- c()
for(i in 1:length(upupname)){
  if(max(AverageExpression(purkinje, assays = 'RNA', features = upupname[i], group.by = 'celltype')[[1]]) > 0.15){
    selename <- c(selename, upupname[i])
    print(i)
  }
}
pkpkpk <- pkpkpk[which(pkpkpk$Upstream_TF %in% selename), ]


# 
pkpkpk <- grn005 %>% filter(Downstream_Gene %in% pkccos)
pkpkpk <- pkpkpk[, c("Upstream_TF", "Downstream_Gene")]
pkpkpk <- pkpkpk[!duplicated(pkpkpk), ]
write.csv(pkpkpk, file = 'grn part of pkc net.csv', quote = F)
# 
subp1 <- pkpkpk %>% filter(Downstream_Gene %in% pkccoso$p1)
subp2 <- pkpkpk %>% filter(Downstream_Gene %in% pkccoso$p2)
subp3 <- pkpkpk %>% filter(Downstream_Gene %in% pkccoso$p3)
nodefile <- unique(c(pkpkpk$Upstream_TF, pkpkpk$Downstream_Gene))
nodefile <- as.data.frame(cbind(nodefile, rep(0, length(nodefile)), seq(1:length(nodefile))))
colnames(nodefile)[2:3] <- c('celltype', 'ifTF')
for(i in 1:nrow(nodefile)){
  if(nodefile[i, 1] %in% grninter$Upstream_TF){
    nodefile[i, 3] <- 1
  } else {
    nodefile[i, 3] <- 0
  }
  if(nodefile[i, 1] %in% c(subp1$Upstream_TF, subp1$Downstream_Gene)) {nodefile[i, 2] <- as.numeric(nodefile[i, 2])+1}
  if(nodefile[i, 1] %in% c(subp2$Upstream_TF, subp2$Downstream_Gene)) {nodefile[i, 2] <- as.numeric(nodefile[i, 2])+2}
  if(nodefile[i, 1] %in% c(subp3$Upstream_TF, subp3$Downstream_Gene)) {nodefile[i, 2] <- as.numeric(nodefile[i, 2])+4}
}
write.csv(nodefile, file = 'grn part of pkc node.csv', quote = F)
# 
# 
(grninter %>% filter(Downstream_Gene == 'CLSTN2'))[,1:6]
# 







# 
pk12cotf1 <- c('LBX1', 'IRF2', 'PBX3', 'MEIS2', 'LHX5', 'NFIX', 'NR4A2', 'PAX2', 'ESRRG', 
               'ZBTB18', 'NHLH2', 'BBX', 'ZIC1', 'SOX9', 'SOX5', 'MXI1', 'JUND', 'SP5', 
               'NR1H2', 'DACH1', 'PAX6')
pkcotf1 <- bitr(pk12cotf1, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
pkcotf1go <- enrichGO(pkcotf1$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', readable = T)
clusterProfiler::dotplot(pkcotf1go, showCategory = 10)
# 
pk12cotf2 <- c('NFKB1', 'MEFB', 'KLF7', 'MEF2A', 'TBP','LHX1', 'NFYC', 'FOXO6', 
               'NR2F1', 'KLF2', 'SP4', 'ZBTB14','MEF2C', 
               'ARID1A', 'UBP1','BCL11B', 'ASCL2', 'DBP', 'EBF3', 
               'RORA', 'SOX4', 'KLF5', 'HES4', 'ARNTL', 'RFX1', 'NFE2L3', 
               'SKOR2', 'TFEB', 'ZBTB38', 'HIVEP2', 'KDM5D', 
               'ARID4A', 'HMGXB4', 'FOXP1', 'STAT4', 'FOXP4', 
               'BCL11A', 'ARNT2')
#'ZNF566','ZNF785', 'ZNF71','ZNF12','ZNF281','ZNF529','ZNF82''ZNF337', 'ZNF793','ZNF519''ZNF33B', 
#', 'ZNF765', 'ZNF467''ZNF611', 'ZNF236','ZNF624', 'ZNF587','ZNF35','ZNF573', 'ZNF513','ZNF667',, , 
pkcotf2 <- bitr(pk12cotf2, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
pkcotf2go <- enrichGO(pkcotf2$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', readable = T)
clusterProfiler::dotplot(pkcotf2go, showCategory = 10)




# gc
load(file = 'cebclean.Rdata')
cebclean$timepoint <- cebclean$orig.ident
cebclean$timepoint[which(cebclean$timepoint == 'E083')] <- 'E102'
cebclean$timepoint[which(cebclean$timepoint == 'E075')] <- 'E083'
cebclean$timepoint[which(cebclean$timepoint == 'E078')] <- 'E108'
cebclean$timepoint[which(cebclean$timepoint == 'E085')] <- 'E117'
cebclean$timepoint[which(cebclean$timepoint == 'E101')] <- 'E101'
cebclean$timepoint[which(cebclean$timepoint == 'E104')] <- 'E093'
# 
DefaultAssay(cebclean) <- 'RNA'
cebclean$celltype <- cebclean$seurat_clusters
levels(cebclean$celltype)[which(levels(cebclean$celltype) %in% c(2,1,10,6))] <- 'DCNs'
levels(cebclean$celltype)[which(levels(cebclean$celltype) %in% c(8,5,3,4,12,7,19,24,15,14))] <- 'GCs'
levels(cebclean$celltype)[which(levels(cebclean$celltype) %in% c(17,9,0))] <- 'UBCs'
levels(cebclean$celltype)[which(levels(cebclean$celltype) %in% c(11,18))] <- 'PROs'
levels(cebclean$celltype)[which(levels(cebclean$celltype) %in% c(16))] <- 'INs'
levels(cebclean$celltype)[which(levels(cebclean$celltype) %in% c(13))] <- 'PKCs'
table(cebclean$celltype)
# 
DefaultAssay(cebclean) <- 'integrated'
cebclean <-  FindSubCluster(cebclean, graph.name = 'integrated_snn', cluster = 18, resolution = 0.5)
DimPlot(cebclean, group.by = 'sub.cluster', label = T)
DefaultAssay(cebclean) <- 'RNA'
prog2m <- subset(cebclean, sub.cluster %in% c('11','18_0', '18_1', '5', '3', '8',
                                              '4', '12', '7', '14', '15', '19', '2', 
                                              '17', '9', '0', '1', '10', '6'))
rm(cebclean)
DimPlot(prog2m)
gc()

# anchorgene
E083 <- subset(prog2m, orig.ident == 'E083')
E083 <- NormalizeData(E083, normalization.method = 'LogNormalize')
E083 <- FindVariableFeatures(E083, selection.method = "vst", nfeatures = 3000)
E083 <- ScaleData(E083)
E083 <- RunPCA(E083, npcs = 50)
ElbowPlot(object = E083, ndims = 40, reduction = "pca")
E083 <- FindNeighbors(E083, reduction = "pca", dims = 1:15)
E083 <- FindClusters(E083, resolution = 1)
E083 <- RunUMAP(E083, reduction = "pca", dims = 1:15)
DimPlot(E083, label = T, label.size = 6)
FeaturePlot(E083, c('SOX2', 'CRABP2', 'TOP2A', 'LMX1A'), label = T) & mockdim
cd14 <- FindMarkers(E083, ident.1 = 14, min.diff.pct = 0.3, logfc.threshold = 0.7)

# 
E083 <- FindSubCluster(E083, cluster = c(11), resolution = 0.5, graph.name = 'RNA_snn')
DimPlot(E083, group.by = 'sub.cluster', label = T)
VlnPlot(E083, 'ROR2', group.by = 'sub.cluster')
E083$sub11 <- E083$sub.cluster
# 11-0/2 ROR2-的
E083 <- FindSubCluster(E083, cluster = c(4), resolution = 0.2, graph.name = 'RNA_snn')
DimPlot(E083, group.by = 'sub.cluster', label = T)
VlnPlot(E083, 'ROR2', group.by = 'sub.cluster')
E083$sub4 <- E083$sub.cluster
# 4-2  ROR2-
E083 <- FindSubCluster(E083, cluster = c(5), resolution = 0.3, graph.name = 'RNA_snn')
DimPlot(E083, group.by = 'sub.cluster', label = T)
VlnPlot(E083, 'ROR2', group.by = 'sub.cluster')
E083$sub5 <- E083$sub.cluster
# 5-2 ROR2-
E083 <- FindSubCluster(E083, cluster = c(7), resolution = 0.3, graph.name = 'RNA_snn')
DimPlot(E083, group.by = 'sub.cluster', label = T)
VlnPlot(E083, 'ROR2', group.by = 'sub.cluster')
E083$sub7 <- E083$sub.cluster
# 7-1/2 ROR2-

twogc <- as.data.frame(cbind(E083$sub4, E083$sub5, E083$sub7, E083$sub11))
twogc <- cbind(twogc, seq(1:nrow(twogc)))
colnames(twogc) <- c('sub4', 'sub5', 'sub7', 'sub11', 'final')
for(i in 1:nrow(twogc)){
  if(TRUE %in% grepl('_', twogc[i, ])){
    twogc[i, 5] <- twogc[i, grepl('_', twogc[i, ])]
  } else {
    twogc[i, 5] <- twogc[i, 1]
  }
}
E083$finaltype <- twogc$final
DimPlot(E083, group.by = 'finaltype', label = T)
E083$fffftype <- E083$finaltype
E083$fffftype[which(E083$fffftype %in% c('16','12'))] <- 'Pros'
E083$fffftype[which(E083$fffftype %in% c('15'))] <- 'preUBC'
E083$fffftype[which(E083$fffftype %in% c('13', '11_2', '9', '1', '10'))] <- 'UBC'
E083$fffftype[which(E083$fffftype %in% c('14', '11_0', '4_2', '6', '0', '8', '5_2', '7_1', '7_2'))] <- 'GC'
E083$fffftype[which(E083$fffftype %in% c('11_1', '11_3', '4_0', '4_1', '7_0', '2', '3', '5_0', '5_1', '5_3'))] <- 'DCN'
table(E083$fffftype)
E083marker <- ccosg(E083, groups = 'fffftype', 
                    assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 50)
gcanchorset <- c(E083marker[[1]][, 1], E083marker[[1]][, 2], E083marker[[1]][, 3], 
                 E083marker[[1]][, 4], E083marker[[1]][, 5])
# 
# 
rm(prog2m)
DefaultAssay(prog2m) <- 'RNA'
for(i in names(table(prog2m$orig.ident))){
  assign(i, subset(prog2m, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 3000))
}
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = gcanchorset, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
prog2m <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
prog2m <- ScaleData(prog2m)
prog2m <- RunPCA(prog2m, npcs = 50)
ElbowPlot(object = prog2m, ndims = 40, reduction = "pca") 
DefaultAssay(prog2m) <- 'integrated'
prog2m <- FindNeighbors(prog2m, reduction = "pca", dims = 1:14)
prog2m <- FindClusters(prog2m, resolution = 1)
prog2m <- RunUMAP(prog2m, reduction = "pca", dims = 1:14)
fgm <- DimPlot(prog2m, label = T, label.size = 5)
fgm
DefaultAssay(prog2m) <- 'RNA'
save(prog2m, file = 'hsa exci lineage.Rdata')
FeaturePlot(prog2m, 'CRABP2')
VlnPlot(prog2m, 'nFeature_RNA', pt.size = 0)
gc18 <- FindMarkers(prog2m, ident.1 = 18, logfc.threshold = 0.7, only.pos = T)
gc16 <- FindMarkers(prog2m, ident.1 = 16, logfc.threshold = 0.7, only.pos = T)

prog2mnd <- subset(prog2m, seurat_clusters != 18)
rm(prog2m)
gc()
DefaultAssay(prog2mnd) <- 'integrated'
prog2mnd <- ScaleData(prog2mnd)
prog2mnd <- RunPCA(prog2mnd, npcs = 50)
ElbowPlot(object = prog2mnd, ndims = 40, reduction = "pca") 
DefaultAssay(prog2mnd) <- 'integrated'
prog2mnd <- FindNeighbors(prog2mnd, reduction = "pca", dims = 1:12)
prog2mnd <- FindClusters(prog2mnd, resolution = 1)
prog2mnd <- RunUMAP(prog2mnd, reduction = "pca", dims = 1:12)
fgmnd <- DimPlot(prog2mnd, label = T, label.size = 5)
fgmnd
DefaultAssay(prog2mnd) <- 'RNA'
save(prog2mnd, file = 'human exci lineage.Rdata')
FeaturePlot(prog2mnd, 'CRABP2')
VlnPlot(prog2mnd, 'nFeature_RNA', pt.size = 0)


load(file = 'human exci lineage.Rdata')
# 
ububub <- subset(prog2mnd, seurat_clusters %in% c(14,16,18,10,8,6))
DefaultAssay(ububub) <- 'RNA'
for(i in names(table(ububub$orig.ident))){
  assign(i, subset(ububub, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 1000))
}
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 1000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
gc()
ububub <- IntegrateData(anchorset = anchor, dims = 1:20)
gc()
ububub <- ScaleData(ububub)
ububub <- RunPCA(ububub, npcs = 50)
ElbowPlot(object = ububub, ndims = 40, reduction = "pca") 
DefaultAssay(ububub) <- 'integrated'
ububub <- FindNeighbors(ububub, reduction = "pca", dims = 1:12)
ububub <- FindClusters(ububub, resolution = 1)
ububub <- RunUMAP(ububub, reduction = "pca", dims = 1:12)
fubc <- DimPlot(ububub, label = T, label.size = 5)
fubc
DefaultAssay(ububub) <- 'RNA'
cub16 <- FindMarkers(ububub, ident.1 = 16, min.diff.pct = 0.3, logfc.threshold = 0.7)
cub11 <- FindMarkers(ububub, ident.1 = 11, ident.2 = 16, only.pos = T)
cub15 <- FindMarkers(ububub, ident.1 = 15, min.diff.pct = 0.3, logfc.threshold = 0.7)
cub15anti514 <- FindMarkers(ububub, ident.1 = 15, ident.2 = c(14,5), logfc.threshold = 0.5)
FeaturePlot(ububub, 'SOX2')
#
table(cebclean$seurat_clusters[names(ububub$seurat_clusters[ububub$seurat_clusters == 16])])
# 
onlyub <- subset(ububub, seurat_clusters %in% setdiff(0:17, c(8,15,11,16)))
DefaultAssay(onlyub) <- 'RNA'
for(i in names(table(onlyub$orig.ident))){
  assign(i, subset(onlyub, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 1000))
}
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 1000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
gc()
onlyub <- IntegrateData(anchorset = anchor, dims = 1:20)
gc()
onlyub <- ScaleData(onlyub)
onlyub <- RunPCA(onlyub, npcs = 50)
ElbowPlot(object = onlyub, ndims = 40, reduction = "pca") 
DefaultAssay(onlyub) <- 'integrated'
onlyub <- FindNeighbors(onlyub, reduction = "pca", dims = 1:10)
onlyub <- FindClusters(onlyub, resolution = 1)
onlyub <- RunUMAP(onlyub, reduction = "pca", dims = 1:10)
fubc <- DimPlot(onlyub, label = T, label.size = 5)
fubc
DefaultAssay(onlyub) <- 'RNA'
# 
onlyub <- FindSubCluster(onlyub, cluster = c(5), resolution = 0.3, graph.name = 'integrated_snn')
DimPlot(onlyub, group.by = 'sub.cluster', label = T)
VlnPlot(onlyub, 'PTF1A', group.by = 'sub.cluster')
onlyub <- subset(onlyub, `sub.cluster` %in% c('0', '1', '3', '4', '5_0', '6', '7', '8', 
                                              '9', '10', '11', '12', '13', '14', '15', '16'))
DefaultAssay(onlyub) <- 'integrated'
onlyub <- ScaleData(onlyub)
onlyub <- RunPCA(onlyub, npcs = 50)
ElbowPlot(object = onlyub, ndims = 40, reduction = "pca") 
DefaultAssay(onlyub) <- 'integrated'
onlyub <- FindNeighbors(onlyub, reduction = "pca", dims = 1:10)
onlyub <- FindClusters(onlyub, resolution = 1)
onlyub <- RunUMAP(onlyub, reduction = "pca", dims = 1:10)
fubc <- DimPlot(onlyub, label = T, label.size = 5)
fubc
DefaultAssay(onlyub) <- 'RNA'
save(onlyub, file = 'human UBC lineage.Rdata')
# 
onlyub$celltype <- onlyub$seurat_clusters
levels(onlyub$celltype)[which(levels(onlyub$celltype) %in% c(13))] <- 'UBCpro1'
levels(onlyub$celltype)[which(levels(onlyub$celltype) %in% c(15))] <- 'UBCpro2'
levels(onlyub$celltype)[which(levels(onlyub$celltype) %in% c(10,11))] <- 'UBCpro3'
levels(onlyub$celltype)[which(levels(onlyub$celltype) %in% c(7,6,12))] <- 'UBCpro3-cc'
levels(onlyub$celltype)[which(levels(onlyub$celltype) %in% c(4,5,1,3,9,0))] <- 'imUBC'
levels(onlyub$celltype)[which(levels(onlyub$celltype) %in% c(2,14,8))] <- 'mUBC'
table(onlyub$celltype)
pdf(file = 'hsa-ceb UBC umap.pdf', height = 3, width = 3)
fprint <- DimPlot(onlyub, group.by = 'celltype', shuffle = T, 
                  cols = c('#84b186', '#68926a', '#9db4c2', '#779ac2', '#63669d', '#83949e')) & 
  mockdim & labs(title = element_blank())
print(fprint)
dev.off()
# 


# 
# 
# 
gcpro <- subset(ububub, seurat_clusters %in% c(16))
prog2mnd <- subset(prog2mnd, seurat_clusters %in% setdiff(0:18, c(14,16,18,10,8,6)))
DefaultAssay(prog2mnd) <- 'RNA'
onlygc <- merge(gcpro, prog2mnd)
rm(prog2mnd, ububub)
gc()
# 
DefaultAssay(onlygc) <- 'RNA'
for(i in names(table(onlygc$orig.ident))){
  assign(i, subset(onlygc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
onlygc <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
onlygc <- ScaleData(onlygc)
onlygc <- RunPCA(onlygc, npcs = 50)
ElbowPlot(object = onlygc, ndims = 40, reduction = "pca") 
DefaultAssay(onlygc) <- 'integrated'
onlygc <- FindNeighbors(onlygc, reduction = "pca", dims = 1:7)
onlygc <- FindClusters(onlygc, resolution = 1)
onlygc <- RunUMAP(onlygc, reduction = "pca", dims = 1:7, min.dist = 0.1)
fgcc <- DimPlot(onlygc, label = T, label.size = 5)
fgcc
DefaultAssay(onlygc) <- 'RNA'
save(onlygc, file = 'human GC lineage.Rdata')

load(file = 'human GC lineage.Rdata')
#
table(ububub$seurat_clusters[names(onlygc$seurat_clusters[onlygc$seurat_clusters == 23])])
table(ppicell$seurat_clusters[names(onlygc$seurat_clusters[onlygc$seurat_clusters == 23])])
# 
gc22 <- FindMarkers(onlygc, ident.1 = 22, logfc.threshold = 0.7, min.diff.pct = 0.3)
gc422 <- FindMarkers(onlygc, ident.1 = c(4,22), ident.2 = c(1,2,3,5,6,7,8,12,13,16,17,18,19,21), 
                     logfc.threshold = 0., min.diff.pct = 0.3)
gc0 <- FindMarkers(onlygc, ident.1 = 0, logfc.threshold = 0.7, min.diff.pct = 0.3)
# 
FeaturePlot(onlygc, c('PCDH7', 'TLX3'), blend = T, blend.threshold = 0, 
            cols = c('#eeeeee', '#68926a', '#80325c')) & mockdim & labs(title = element_blank())
# 
DefaultAssay(onlygc)
DefaultAssay(onlyub)
# 
onlygc$celltype <- onlygc$seurat_clusters
levels(onlygc$celltype)[which(levels(onlygc$celltype) %in% c(4))] <- 'likeUBC'
ubwithsimi <- merge(onlyub, subset(onlygc, celltype == 'likeUBC'))
likeubm <- FindMarkers(ubwithsimi, ident.1 = 'likeUBC', ident.2 = c('imUBC', 'mUBC'), 
                       min.diff.pct = 0.3, logfc.threshold = 0.7, group.by = 'celltype')
VlnPlot(onlygc, 'nFeature_RNA', pt.size = 0)

# 
DefaultAssay(ubwithsimi) <- 'RNA'
for(i in names(table(ubwithsimi$orig.ident))){
  assign(i, subset(ubwithsimi, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 1000))
}
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 1000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
gc()
ubwithsimi <- IntegrateData(anchorset = anchor, dims = 1:20)
gc()
ubwithsimi <- ScaleData(ubwithsimi)
ubwithsimi <- RunPCA(ubwithsimi, npcs = 50)
ElbowPlot(object = ubwithsimi, ndims = 40, reduction = "pca") 
DefaultAssay(ubwithsimi) <- 'integrated'
ubwithsimi <- FindNeighbors(ubwithsimi, reduction = "pca", dims = 1:8)
ubwithsimi <- FindClusters(ubwithsimi, resolution = 1)
ubwithsimi <- RunUMAP(ubwithsimi, reduction = "pca", dims = 1:8)
fubc <- DimPlot(ubwithsimi, label = T, label.size = 5)
fubc
DefaultAssay(ubwithsimi) <- 'RNA'
ubwithsimi$celltype <- ubwithsimi$seurat_clusters
levels(ubwithsimi$celltype)[which(levels(ubwithsimi$celltype) %in% c(16))] <- 'UBCpro1'
levels(ubwithsimi$celltype)[which(levels(ubwithsimi$celltype) %in% c(17))] <- 'UBCpro2'
levels(ubwithsimi$celltype)[which(levels(ubwithsimi$celltype) %in% c(15,3))] <- 'UBCpro3'
levels(ubwithsimi$celltype)[which(levels(ubwithsimi$celltype) %in% c(11,10,13,18))] <- 'UBCpro3-cc'
levels(ubwithsimi$celltype)[which(levels(ubwithsimi$celltype) %in% c(0,4,12,1,2,8,5,9))] <- 'UBC-L1'
levels(ubwithsimi$celltype)[which(levels(ubwithsimi$celltype) %in% c(6,7,14))] <- 'UBC-L2'
table(ubwithsimi$celltype)
save(ubwithsimi, file = 'human UBC two lineage.Rdata')
# 



# 
load(file = 'human UBC two lineage.Rdata')
# 
DefaultAssay(ubwithsimi) <- 'RNA'
fub <- DimPlot(ubwithsimi, label = T)
fub
ubwithsimi$celltype <- factor(ubwithsimi$celltype, levels = c('UBCpro1', 'UBCpro2', 'UBCpro3', 'UBCpro3-cc', 'UBC-L1', 'UBC-L2'))
# 
ubm2 <- FindMarkers(ubwithsimi, ident.1 = c(6,7,14), min.diff.pct = 0.3, logfc.threshold = 0.7)
pdf(file = 'hsa-ceb UBC two lineage umap.pdf', height = 3, width = 3)
fprint <- DimPlot(ubwithsimi, group.by = 'celltype', shuffle = T, 
                  cols = c('#68926a', '#84b186', '#779ac2', '#9db4c2', '#63669d', '#83949e')) & 
  mockdim & labs(title = element_blank())
print(fprint)
dev.off()
# 
tempoubc <- c('CYP26B1','GDF10', 'COL2A1', 'PIK3C2G', 'MEGF10', 'PTPRT', 'ATP1A2', 'SLC1A3', 
              'SOX2', 'WLS', 'CRYAB', 'CALCB', 'SULF1', 'CITED4', 'MSX1',
              'TRPM3', 'NEUROG1', 'EYS', 'GRIK1', 'ATOH1', 'MKI67', 'TOP2A', 
               'LMX1A', 'DGKB', 'DCDC1', 'GALNTL6', 'CNTN1', 'CDH18', 
              'MASP1', 'PLP1', 'ZFP36L1', 'NCKAP5', 
              'STMN2', 'CACNA2D1', 'CCSER1', 'EOMES')
# , 'CALCA', 'SPARC', 'TMTC2''LRP1B',, 'SMOC1', 'RFX4'
pdf(file = 'hsa-ceb UBC two lineage dotplot.pdf', height = 3, width = 11)
print(DotPlot(ubwithsimi, features = tempoubc, group.by = 'celltype', 
              cluster.idents = F, cols = c('#dddddd', '#63669d')) + RotatedAxis())
dev.off()
# 
#
ubctwo <- subset(ubwithsimi, celltype %in% c('UBC-L1', 'UBC-L2'))
ubc2cosm <- ccosg(ubctwo, groups = 'celltype', 
                 assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 20)
ubc2cosm <- cbind(ubc2cosm$names[1], ubc2cosm$scores[1], ubc2cosm$names[2], ubc2cosm$scores[2])
ubgrninput <- c(ubc2cosm[,1][!grepl(ubc2cosm[, 1], pattern = '^MT-|^RPL|^RPS')], ubc2cosm[,3])
# 
excigrn <- read.csv(file = "E:/cerebellum/GCS+UBC grn.csv")

excigrn <- cbind(excigrn, updownadj, creupadj)
egrn005 <- excigrn %>% filter(updownadj <0.05 & creupadj <0.05)
# 
ubcubc <- egrn005 %>% filter(Downstream_Gene %in% ubgrninput)
ubcubc <- ubcubc[, c("Upstream_TF", "Downstream_Gene")]
ubcubc <- ubcubc[!duplicated(ubcubc), ]
write.csv(ubcubc, file = 'grn part of ubc net.csv', quote = F)
# 
ubc2cosm <- ccosg(ubctwo, groups = 'celltype', 
                  assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 100)
ubc2cosm <- cbind(ubc2cosm$names[1], ubc2cosm$scores[1], ubc2cosm$names[2], ubc2cosm$scores[2])
ubgrninput <- c(ubc2cosm[,1][!grepl(ubc2cosm[, 1], pattern = '^MT-|^RPL|^RPS')], ubc2cosm[,3])
# 
excigrn <- read.csv(file = "E:/cerebellum/GCS+UBC grn.csv")
# 

excigrn <- cbind(excigrn, updownadj, creupadj)
egrn005 <- excigrn %>% filter(updownadj <0.05 & creupadj <0.05)
# 
ubcubc <- egrn005 %>% filter(Downstream_Gene %in% ubgrninput & Upstream_TF %in% ubgrninput)
ubcubc <- ubcubc[, c("Upstream_TF", "Downstream_Gene")]
ubcubc <- ubcubc[!duplicated(ubcubc), ]
write.csv(ubcubc, file = 'grn part of ubc net.csv', quote = F)
# 
# correlation
ubcforcor <- unique(c(ubcubc$Upstream_TF, ubcubc$Downstream_Gene))
allubccor <- c()
for(i in ubcforcor){
  cor1 <- apply(ubwithsimi@assays$RNA@data[ubcforcor, ], 1, 
                function(x) cor.test(x, ubwithsimi@assays$RNA@data[i, ]))
  cor1 <- unlist(lapply(cor1, function(x) x$estimate))
  allubccor <- cbind(allubccor, cor1)
}
rownames(allubccor) <- unlist(lapply(strsplit(rownames(allubccor), split = '.', fixed = T), function(x) x[[1]]))
colnames(allubccor) <- rownames(allubccor)
for(i in 1:nrow(allubccor)){
  allubccor[i, i] <- 0
}
pdf(file = 'hsa-ceb UBC network correlation.pdf', height = 8, width = 8)
pheatmap(allubccor, cluster_rows = T, cluster_cols = T, border = F, scale = 'none', 
         color = colorRampPalette(c('#63669d', "#f0eeec", '#80325c'))(10),
         show_rownames = T, show_colnames = F, fontsize = 250/length(tempoubc), breaks = seq(-0.5, 0.5, 0.1))
dev.off()
#
#
pdf(file = 'hsa-ceb anti-pair CNTNAP5 and LMX1A.pdf', height = 3, width = 12)
fprint <- FeaturePlot(ubwithsimi, features = c('CNTNAP5', 'LMX1A'), blend = T, blend.threshold = 0.1, 
                      cols = c('#eeeeee', '#036eb8', '#D40041'), pt.size = 0.1)& mockdim& labs(title = element_blank())
print(fprint)
dev.off()
# 
# 
ubcfortraj <- subset(ubwithsimi, seurat_clusters %in% setdiff(0:18, c(11,10,18,13)))
data <- GetAssayData(ubcfortraj, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
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
# 
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
# 
tempoubc <- c('PIK3C2G', 'CRYAB', 'TRPM3', 'EOMES', 'TMTC2', 'LRP1B')
pdf(file = 'hsa-ceb UBC traj change gene.pdf', height = 4, width = 9)
fprint <- plot_genes_in_pseudotime(cds[tempoubc, ], cell_size = 1, ncol = 3, 
                                   panel_order = tempoubc) + 
  scale_y_continuous(trans='log1p') + 
  scale_color_gradientn(colors = c('#63669d', '#a9d6af')) + theme_bw() & 
  theme(axis.title = element_blank()) + NoLegend()
print(fprint)
dev.off()
# 
# 
ubcregu <- loomR::connect(filename = "E:/cerebellum/result/scenic_integrated-output.loom",
                          mode = 'r', skip.validate = TRUE)
metaloom <- ubcregu$get.attribute.df(MARGIN = 2)
DimPlot(ubwithsimi, cells.highlight = rownames(metaloom), label = T)
length(ubcregu[["row_attrs/Regulons"]][1]) # 317 regulon
regumatrix <- metaloom[, 13:329]
colreguname <- gsub(pattern = '_(+)', replacement = '', fixed = T, 
                    names(ubcregu[["row_attrs/Regulons"]][1]))
colnames(regumatrix) <- colreguname
pseudot <- cds@principal_graph_aux$UMAP$pseudotime
pseudot <- pseudot[rownames(regumatrix)]
FALSE %in% c(names(pseudot) == rownames(regumatrix))
scluster <- metaloom[,3]
regumatrix <- cbind(pseudot, regumatrix)
regumatrix <- cbind(scluster, regumatrix)
regumatrix <- regumatrix[order(regumatrix$pseudot, decreasing = F), ]
tryaggregu <- aggregate(regumatrix, by = list(regumatrix$scluster), FUN = mean)
tryaggregu <- tryaggregu[c(3,4,2,1), ]
# 
seleubcregu <- c('PRDM15', 'EN1', 'FOXO1', 'ELF1', 'ETV3', 'E4F1', 'GABPA', 'GMEB2', 
                 'YY2', 'PPARGC1A', 'SOX15', 'ESRRB', 'ETV3L', 'ZNF41', 'ZNF607', 
                 'HDAC1', 'NR2F1', 'ZNF596', 'ATOH7', 'IRF3', 'JUN', 'ZEB1', 'ATF4', 
                 'SAP30', 'ASCL1', 'ZNF66', 'ATOH1')
reguheat1 <- pheatmap(tryaggregu[, seleubcregu], cluster_rows = F, cluster_cols = T, scale = 'column', border = F, 
                      color = colorRampPalette(c('#f2f0ef', '#84b186', '#4b5861'))(10),
                      show_rownames = F, show_colnames = T, breaks = seq(-1,1, 0.2))
reguheat1
pdf(file = 'hsa-ceb two UBC regulon with pseudo heat.pdf')
print(reguheat1)
dev.off()
# reguheat1 <- pheatmap(regumatrix[, 3:290], cluster_rows = F, cluster_cols = T, scale = 'none', border = F, 
#                       color = colorRampPalette(c('#f2f0ef', '#84b186', '#4b5861'))(10),
#                       show_rownames = F, show_colnames = T, breaks = seq(0, 0.5, 0.05))
# pdf(file = 'test.pdf', height = 6, width = 30)
# print(reguheat1)
# dev.off()
# 

# 
load(file = 'human UBC two lineage.Rdata')
DefaultAssay(ubwithsimi) <- 'RNA'
gcforgsea <- FindMarkers(ubwithsimi, ident.1 = 'UBC-L1', ident.2 = 'UBC-L2',
                         group.by = 'celltype', logfc.threshold = 0, min.diff.pct = 0)
gcforgsea <- cbind(gcforgsea, rownames(gcforgsea))
colnames(gcforgsea)[6] <- 'SYMBOL'
gcid <- bitr(rownames(gcforgsea), OrgDb = 'org.Hs.eg.db', 
             fromType = 'SYMBOL', toType = c('ENTREZID', 'SYMBOL'))
gcgsea <- merge(gcforgsea, gcid, by = 'SYMBOL')
gcgsea <- gcgsea[order(gcgsea$avg_log2FC, decreasing = T), ]
ggcclist <- gcgsea$avg_log2FC
names(ggcclist) <- gcgsea$ENTREZID
# 
gcgsego <- gseGO(geneList = ggcclist, OrgDb = org.Hs.eg.db, verbose = T)
options(clusterProfiler.download.method = "wininet")
gcgsekegg <- gseKEGG(geneList = ggcclist, organism = "hsa")
# 
#
# lineage feature plot
pdf(file = 'hsa-ceb UBC feature of two sublinege.pdf', height = 4, width = 4)
f1 <- FeaturePlot(ubwithsimi, features = c('DGKB', 'PLCB4', 'KCNQ5', 'PLCB1'), 
                  cols = c('#dedede', 'navy')) & mockdim & labs(title = element_blank())
print(f1)
dev.off()
# 

load(file = 'human GC lineage.Rdata')
onlygc <- subset(onlygc, seurat_clusters != 4)
gc()
DefaultAssay(onlygc) <- 'RNA'
for(i in names(table(onlygc$orig.ident))){
  assign(i, subset(onlygc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
rm(onlygc)
gc()
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
# save(anchor, file = 'anchor.Rdata')
# #
load(file = 'anchor.Rdata')
onlygc <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
onlygc <- ScaleData(onlygc)
onlygc <- RunPCA(onlygc, npcs = 50)
ElbowPlot(object = onlygc, ndims = 40, reduction = "pca") 
DefaultAssay(onlygc) <- 'integrated'
onlygc <- FindNeighbors(onlygc, reduction = "pca", dims = 1:13)
onlygc <- FindClusters(onlygc, resolution = 1)
onlygc <- RunUMAP(onlygc, reduction = "pca", dims = 1:6, min.dist = 0.1)
fgcc <- DimPlot(onlygc, label = T, label.size = 5)
fgcc
DefaultAssay(onlygc) <- 'RNA'
# 23 
onlygc <- subset(onlygc, seurat_clusters != 23)
gc()
DefaultAssay(onlygc) <- 'RNA'
for(i in names(table(onlygc$orig.ident))){
  assign(i, subset(onlygc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
rm(onlygc)
gc()
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
onlygc <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
onlygc <- ScaleData(onlygc)
onlygc <- RunPCA(onlygc, npcs = 50)
ElbowPlot(object = onlygc, ndims = 40, reduction = "pca") 
DefaultAssay(onlygc) <- 'integrated'
onlygc <- FindNeighbors(onlygc, reduction = "pca", dims = 1:20)
onlygc <- FindClusters(onlygc, resolution = 1)
onlygc <- RunUMAP(onlygc, reduction = "pca", dims = 1:20)
fgcc <- DimPlot(onlygc, label = T, label.size = 5)
fgcc
DefaultAssay(onlygc) <- 'RNA'
mmm23 <- FindMarkers(onlygc, ident.1 = 23, min.diff.pct = 0.3, logfc.threshold = 0.7)
# 
onlygc <- subset(onlygc, seurat_clusters != 23)
rm(anchor)
gc()
DefaultAssay(onlygc) <- 'RNA'
for(i in names(table(onlygc$orig.ident))){
  assign(i, subset(onlygc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
rm(onlygc)
gc()
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
# save(anchor, file = 'anchor.Rdata')
# #
# load(file = 'anchor.Rdata')
onlygc <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
onlygc <- ScaleData(onlygc)
onlygc <- RunPCA(onlygc, npcs = 50)
ElbowPlot(object = onlygc, ndims = 40, reduction = "pca") 
DefaultAssay(onlygc) <- 'integrated'
onlygc <- FindNeighbors(onlygc, reduction = "pca", dims = 1:35)
onlygc <- FindClusters(onlygc, resolution = 1)
onlygc <- RunUMAP(onlygc, reduction = "pca", dims = 1:35)
fgcc <- DimPlot(onlygc, label = T, label.size = 5)
fgcc
DefaultAssay(onlygc) <- 'RNA'
# mmm20 <- FindMarkers(onlygc, ident.1 = 20, min.diff.pct = 0.3, logfc.threshold = 0.7)
onlygc <- CellCycleScoring(onlygc, s.features = cc.genes.updated.2019$s.genes, 
                           g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
pdf(file = 'hsa-ceb GC cellcycle.pdf', width = 4, height = 4)
f1 <- DimPlot(onlygc, label = T, group.by = 'Phase', cols = c('#916e8c', '#84b186', '#779ac2')) & 
  mockdim & labs(title = element_blank())
print(f1)
dev.off()
# 
# 
onlygc@active.ident <- onlygc$seurat_clusters
DefaultAssay(onlygc) <- 'integrated'
onlygc <-  FindSubCluster(onlygc, graph.name = 'integrated_snn', cluster = 7, resolution = 0.3)
DimPlot(onlygc, group.by = 'sub.cluster', label = T)
# 
onlygc$protype <- onlygc$sub.cluster
onlygc$protype[which(onlygc$protype %in% c('11','6','3','0','23','7_0','7_2','4','2','19','14','22'))] <- 'pro'
onlygc$protype[which(onlygc$protype %in% c('10','17','5','8','15','20'))] <- 'sg2m'
onlygc$protype[which(onlygc$protype %in% c('9','16','18','13','1','12','21','7_1','7_3'))] <- 'afcc'
table(onlygc$protype)
# 
pdf(file = 'hsa-ceb GC pre umap.pdf', width = 4, height = 4)
f1 <- DimPlot(onlygc, label = F, group.by = 'protype', 
              cols = c( '#779ac2', '#63669d','#83949e', '#80325c')) & mockdim & labs(title = element_blank())
print(f1)
dev.off()
#
pdf(file = 'hsa-ceb GC a-p anti-pair.pdf', width = 12, height = 3)
f1 <- FeaturePlot(onlygc, c('TLX3', 'BARHL1'), blend = T, blend.threshold = 0, 
                  cols = c('#eeeeee', '#68926a', '#63669d'), max.cutoff = 1.5) & mockdim & labs(title = element_blank())
print(f1)
dev.off()
# feature plot
DefaultAssay(onlygc) <- 'RNA'
gcmainlist <- c('HEY1', 'TRPM3', 'SOX2', 'ZFP36L1', 'PCDH7', 'ASTN2', 
                'ATOH1', 'NHLH1', 'TOP2A', 'PCNA', 'PRPH', 'ZFP36L1', 'SFRP1', 'NEUROD1')
for(i in gcmainlist){
  png(filename = paste0('GC ', i, '.png'), width = 600, height = 600)
  f1 <- FeaturePlot(onlygc, i, cols = c('#dddddd', 'navy')) & mockdim & 
    labs(title = element_blank())
  print(f1)
  dev.off()
}
DefaultAssay(onlygc) <- 'RNA'
gcmainlist <- c('PCDH7', 'PTPRK', 'NRG3', 'LINC02609',
                'KCNIP4', 'ADARB2', 'NTF3', 'PARM1')
for(i in gcmainlist){
  png(filename = paste0('GC ', i, '.png'), width = 600, height = 600)
  f1 <- FeaturePlot(onlygc, i, cols = c('#dddddd', 'navy')) & mockdim & 
    labs(title = element_blank())
  print(f1)
  dev.off()
}
## 
DefaultAssay(onlygc) <- 'RNA'
mmmprotiter <- FindMarkers(onlygc, ident.1 = 'pro', ident.2 = 'afcc', group.by = 'protype', 
                           min.diff.pct = 0.2, logfc.threshold = 0.5)
glist<- mmmprotiter %>% filter(avg_log2FC > 0)
symglist <- bitr(rownames(glist), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
upgo <- enrichGO(symglist$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', readable = T)
clusterProfiler::dotplot(upgo, showCategory = 10)
# terminated
glist<- mmmprotiter %>% filter(avg_log2FC < 0)
symglist <- bitr(rownames(glist), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
downgo <- enrichGO(symglist$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', readable = T)
clusterProfiler::dotplot(downgo, showCategory = 10)
# 
pdf(file = 'hsa-ceb GC pro and terminated.pdf', height = 6, width = 8)
print(clusterProfiler::dotplot(upgo, showCategory = 10)+ 
        scale_color_gradientn(colors = c('#63669d', '#edeaf1')))
print(clusterProfiler::dotplot(downgo, showCategory = 10)+ 
        scale_color_gradientn(colors = c('#779ac2', '#edeaf1')))
dev.off()
# 
#
save(onlygc, file = 'human only GC.Rdata')



# 
load(file = 'human only GC.Rdata')
progc <- subset(onlygc, protype %in% c('afcc', 'pro') & Phase == 'G1')
gc(); DefaultAssay(progc) <- 'RNA'
for(i in names(table(progc$orig.ident))){
  assign(i, subset(progc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
rm(progc); gc()
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
progc <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
progc <- ScaleData(progc)
progc <- RunPCA(progc, npcs = 50)
ElbowPlot(object = progc, ndims = 40, reduction = "pca") 
DefaultAssay(progc) <- 'integrated'
progc <- FindNeighbors(progc, reduction = "pca", dims = 1:21)
progc <- FindClusters(progc, resolution = 1)
progc <- RunUMAP(progc, reduction = "pca", dims = 1:21, min.dist = 0.1)
fgcp <- DimPlot(progc, label = T, label.size = 5)
fgcp
DefaultAssay(progc) <- 'RNA'
save(progc, file = 'human G1 GCs.Rdata')
# umap
colforumap <- c('#dbeddf', '#b8d4e9', '#84b186', '#d9c6da', '#80325c', '#c2b9ee', '#f8e4db', 
                '#63669d', '#bab8a9', '#4b5861', '#8f9a8a', '#d1caad', '#916e8c', '#83949e', 
                '#8f9a8a', '#00709f', '#5c9400', '#b7c2c8', '#edeaf1')
pdf(file = 'hsa-ceb GC only G1 umap.pdf', width = 4, height = 3)
f1 <- DimPlot(progc, label = F, 
              cols = colforumap) & mockdim & labs(title = element_blank())
print(f1)
dev.off()
# 
png('hsa-ceb GC only G1 selected from umap.png', height = 600, width = 600)
print(DimPlot(onlygc, cells.highlight = colnames(progc), cols.highlight = '#779ac2')
      & mockdim)
dev.off()
# featureplot
DefaultAssay(progc) <- 'RNA'
gcmainlist <- c('NR2F2', 'DOK6', 'EBF2', 'BARHL1',
                'GALNTL6', 'PRKCB', 'NRXN3', 'TLX3')
for(i in gcmainlist){
  png(filename = paste0('G1 GC ', i, '.png'), width = 500, height = 500)
  f1 <- FeaturePlot(progc, i, cols = c('#dddddd', 'navy')) & mockdim & 
    labs(title = element_blank())
  print(f1)
  dev.off()
}



ggid <- bitr(rownames(progc), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
progc <- progc[as.numeric(rownames(ggid)), ]
progc <- progc@assays$RNA@data
rownames(progc) <- progc$ENTREZID
# 
stemness <- StemSC(as.matrix(progc))
# totalstem <- c()
# for(i in 1:ncol(testsc)){
#   stemness <- StemSC(as.matrix(testsc[, i]))
#   print(i)
#   totalstem <- c(totalstem, stemness)
# }
# 1:3673, 5k一个bin 

# FeaturePlot(progc, 'stemindex', cols = c('#dddddd', 'navy', 'green', 'red', 'yellow'), min.cutoff = 0.4)
# test <- aggregate(progc$stemindex, by = list(progc$seurat_clusters), FUN=mean)
results <- CytoTRACE(as.matrix(test@assays$RNA@counts))
plotCytoTRACE(results, emb = test@reductions$umap@cell.embeddings, gene = "HEY1", 
              colors = c('#fcfdbf', '#BA3878', '#010005'))
# 
load(file = "E:/cerebellum/result/human G1 GCs stems score.Rdata")
load(file = 'human G1 GCs.Rdata')
plotCytoTRACE(results, emb = progc@reductions$umap@cell.embeddings, gene = "HEY1", 
              colors = c('#fcfdbf', '#BA3878', '#010005'))
# 
progc$stemindex <- scale(stemness)
VlnPlot(progc, 'stemindex', pt.size = 0)
test <- aggregate(progc$stemindex, by = list(progc$seurat_clusters), FUN=mean)
pdf(file = 'hsa-ceb GC stemSC result.pdf', height = 3, width = 3)
print(FeaturePlot(progc, 'stemindex', cols = c( '#010005', '#BA3878','#fcfdbf')) & 
        mockdim & labs(title = element_blank()))
dev.off()
#
DefaultAssay(progc) <- 'integrated'
data <- GetAssayData(progc, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data, gene_metadata = gene_annotation); rm(data, gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(progc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
# 
cds <- cluster_cells(cds, k = 13)
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
  scale_color_gradientn(colors = c('#a9d6af', '#63669d'))
p7
## monocle
pdf(file = 'hsa-ceb GC monocle traj.pdf', height = 3, width = 4)
print(p7)
print(p7 & mockdim & labs(title = element_blank()))
dev.off()
#  
save(cds, file = 'human GC monocle.Rdata')



#
load(file = 'human GC monocle.Rdata')
rrtoii <- graph_test(cds, neighbor_graph = "principal_graph")
write.csv(rrtoii, file = 'GC monocle pseudo temporal.csv')
# 
rrtoii <- read.csv(file = 'GC monocle pseudo temporal.csv', row.names = 1)
changgg <- rrtoii %>% filter(q_value < 0.01)
changgg <- changgg[order(changgg$morans_I, decreasing = T), ]
temporalg <- c('MGP','ZFP36L1', 'ATOH1', 'NRG1', 'CNTN1', 'NEUROD1', 'PRPH', 'PCDH9', 'NAV3')
plot_genes_in_pseudotime(cds[temporalg[7:8], ], cell_size = 1, ncol = 3) + 
  scale_y_continuous(trans='log1p') + 
  scale_color_gradientn(colors = c('#a9d6af', '#63669d')) + theme_bw() & 
  theme(axis.title = element_blank()) + NoLegend()
# 
pdf(file = 'hsa-ceb GC only G1 traj change gene.pdf', height = 4, width = 9)
for(i in c(1,4,7)){
  fprint <- plot_genes_in_pseudotime(cds[temporalg[i:(i+2)], ], cell_size = 1, ncol = 3) + 
    scale_y_continuous(trans='log1p') + 
    scale_color_gradientn(colors = c('#a9d6af', '#63669d')) + theme_bw() & 
    theme(axis.title = element_blank()) + NoLegend() 
  print(fprint)
}
dev.off()

# 
load(file = 'human G1 GCs.Rdata')
endgc <- subset(progc, seurat_clusters %in% c(4,12,7,14,13,16,2,10,15,9))
gc()
DefaultAssay(endgc) <- 'RNA'
for(i in names(table(endgc$orig.ident))){
  assign(i, subset(endgc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
rm(endgc); gc()
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
endgc <- IntegrateData(anchorset = anchor, dims = 1:30); gc()
endgc <- ScaleData(endgc)
endgc <- RunPCA(endgc, npcs = 50)
ElbowPlot(object = endgc, ndims = 40, reduction = "pca") 
DefaultAssay(endgc) <- 'integrated'
endgc <- FindNeighbors(endgc, reduction = "pca", dims = 1:13)
endgc <- FindClusters(endgc, resolution = 1)
endgc <- RunUMAP(endgc, reduction = "pca", dims = 1:13)
fend <- DimPlot(endgc, label = T, label.size = 5)
fend
DefaultAssay(endgc) <- 'RNA'
# dim 13
end8 <- FindMarkers(endgc, ident.1 = 8, min.diff.pct = 0.3, logfc.threshold = 0.7)
end010 <- FindMarkers(endgc, ident.1 = c(0,10), ident.2 = c(3,9,5), only.pos = T)
end15 <- FindMarkers(endgc, ident.1 = 15)
save(endgc, file = 'human GC after cc.Rdata')

# 
# 
# progenitor
DefaultAssay(progc) <- 'integrated'
atatgc <- merge(subset(progc, seurat_clusters %in% c(0,1,3,5,6,8,11,17,18)), 
                subset(endgc, seurat_clusters %in% c(8)))
DefaultAssay(atatgc) <- 'RNA'
for(i in names(table(atatgc$orig.ident))){
  assign(i, subset(atatgc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
rm(atatgc); gc()
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
rm(E075, E078, E083, E085, E101, E104); gc()
atatgc <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
atatgc <- ScaleData(atatgc)
atatgc <- RunPCA(atatgc, npcs = 50)
ElbowPlot(object = atatgc, ndims = 40, reduction = "pca") 
DefaultAssay(atatgc) <- 'integrated'
atatgc <- FindNeighbors(atatgc, reduction = "pca", dims = 1:15)
atatgc <- FindClusters(atatgc, resolution = 1)
atatgc <- RunUMAP(atatgc, reduction = "pca", dims = 1:15)
fggc <- DimPlot(atatgc, label = T, label.size = 5)
fggc
DefaultAssay(atatgc) <- 'RNA'
# 
atatgc$egltype <- atatgc$seurat_clusters
levels(atatgc$egltype)[which(levels(atatgc$egltype) %in% c(16,9,5,4,12,15,2,10))] <- 'ND'
levels(atatgc$egltype)[which(levels(atatgc$egltype) %in% c(14,7,11,6,1,8,13,3,0))] <- 'AT'
table(atatgc$egltype)
# 
save(atatgc, file = 'human GC before cc.Rdata')
# 
load(file = 'human GC before cc.Rdata')
# 
pdf(file = 'hsa-ceb GC before cc.pdf', width = 3, height = 3)
f1 <- DimPlot(atatgc, label = F, group.by = 'egltype', 
              cols = c('#80325c', '#779ac2')) & mockdim & labs(title = element_blank())
print(f1)
dev.off()
## 
matantind <- FindMarkers(atatgc, ident.1 = 'AT', group.by = 'egltype', 
                         min.diff.pct = 0.2, logfc.threshold = 0.5)
diffegl <- c('SOX2', 'ATOH1', 'BTBD11', 'ELAVL2','NEUROD1', 'PPP1R17')
pdf(file = 'hsa-ceb GC diff EGL.pdf', height = 10, width = 4)
fprint <- StackedVlnPlot(atatgc, features = diffegl, group.by = 'egltype')
print(fprint)
dev.off()
# 
#
DefaultAssay(atatgc) <- 'RNA'
gcforgsea <- FindMarkers(atatgc, ident.1 = 'AT', ident.2 = 'ND',
                         group.by = 'egltype', logfc.threshold = 0, min.diff.pct = 0)
write.csv(gcforgsea, file = 'GC before cc logFC for gsea.csv')
# 
gcforgsea <- read.csv(file = 'GC before cc logFC for gsea.csv', row.names = 1)
gcforgsea <- cbind(gcforgsea, rownames(gcforgsea))
colnames(gcforgsea)[6] <- 'SYMBOL'
gcid <- bitr(rownames(gcforgsea), OrgDb = 'org.Hs.eg.db', 
             fromType = 'SYMBOL', toType = c('ENTREZID', 'SYMBOL'))
gcgsea <- merge(gcforgsea, gcid, by = 'SYMBOL')
gcgsea <- gcgsea[order(gcgsea$avg_log2FC, decreasing = T), ]
# 
ggcclist <- gcgsea$avg_log2FC
names(ggcclist) <- gcgsea$ENTREZID
# 
gcgsego <- gseGO(geneList = ggcclist, OrgDb = org.Hs.eg.db, verbose = T)
options(clusterProfiler.download.method = "wininet")
gcgsekegg <- gseKEGG(geneList = ggcclist, organism = "hsa")
# 
pdf(file = 'hsa-ceb GC AT-ND gsea.pdf', width = 8, height = 5)
fprint1 <- gseaplot2(gcgsego, geneSetID = c('GO:0002181', 'GO:0042254'), 
                     color = c('#00AAAA', '#00AA55'))
fprint2 <- gseaplot2(gcgsego, geneSetID = c('GO:0007409', 'GO:0048667'), 
                     color = c('#5555FF', '#B088FF'))
fprint4 <- gseaplot2(gcgsekegg, geneSetID = c('hsa04015', 'hsa04062', 'hsa04014'), 
                     color = c('#00AAAA', '#00AA55', '#009FCC'))
fprint5 <- gseaplot2(gcgsekegg, geneSetID = c('hsa03010'), 
                     color = c('#00AAAA'))
print(fprint1)
print(fprint2)
print(fprint4)
print(fprint5)
dev.off()
# 
DefaultAssay(progc) <- 'RNA'
mmgcgc12 <- FindMarkers(progc, ident.1 = 12, ident.2 = c(14,16), min.diff.pct = 0.3, 
                        logfc.threshold = 0.7)
mmgcgc4 <- FindMarkers(progc, ident.1 = 4, ident.2 = c(7), min.diff.pct = 0.4, 
                        logfc.threshold = 0.3)
mmgcgc14 <- FindMarkers(progc, ident.1 = 14, ident.2 = c(4), min.diff.pct = 0.3, 
                       logfc.threshold = 0.7)
mmgcgc10 <- FindMarkers(progc, ident.1 = 10, ident.2 = c(9), min.diff.pct = 0.3, 
                        logfc.threshold = 0.7)
# lineage
DefaultAssay(progc) <- 'integrated'
progc <- FindSubCluster(progc, cluster = c(17), resolution = 0.4, graph.name = 'integrated_snn')
DimPlot(progc, group.by = 'sub.cluster', label = T)
VlnPlot(progc, 'TLX3', group.by = 'sub.cluster')
progc$sub17 <- progc$sub.cluster
# 17-1 TLX3-
DefaultAssay(progc) <- 'integrated'
progc <- FindSubCluster(progc, cluster = c(1), resolution = 1.5, graph.name = 'integrated_snn')
DimPlot(progc, group.by = 'sub.cluster', label = T)
DefaultAssay(progc) <- 'RNA'
VlnPlot(progc, 'TLX3', group.by = 'sub.cluster')
progc$sub1 <- progc$sub.cluster
# 0,1,2,4,6,10,11 TLX3-
DefaultAssay(progc) <- 'integrated'
progc <- FindSubCluster(progc, cluster = c(5), resolution = 0.2, graph.name = 'integrated_snn')
DimPlot(progc, group.by = 'sub.cluster', label = T)
VlnPlot(progc, 'TLX3', group.by = 'sub.cluster')
progc$sub5 <- progc$sub.cluster
# 5-1,5-2 TLX3-
# 
# subcluster
twogc <- as.data.frame(cbind(progc$sub17, progc$sub5, progc$sub1))
twogc <- cbind(twogc, seq(1:nrow(twogc)))
colnames(twogc) <- c('sub17', 'sub1', 'sub5', 'final')
for(i in 1:nrow(twogc)){
  if(TRUE %in% grepl('_', twogc[i, ])){
    twogc[i, 'final'] <- twogc[i, grepl('_', twogc[i, ])]
  } else {
    twogc[i, 'final'] <- twogc[i, 1]
  }
}
DefaultAssay(progc) <- 'RNA'
# 
progc$finaltype <- twogc$final
DimPlot(progc, group.by = 'finaltype', label = T)
progc$twotype <- progc$finaltype
progc$twotype[which(progc$twotype %in% c('3','12','4','7','13','14','17_0','17_2','5_0', '1_3', 
                                         '1_5', '1_7', '1_8', '1_9'))] <- 'post'
progc$twotype[which(progc$twotype %in% c('0','2','6','8','9','10','11','15','18','16', '17_1', 
                                         '5_1', '5_2', '1_0','1_1','1_2','1_4','1_6', '1_10', 
                                         '1_11'))] <- 'ant'
table(progc$twotype)
DimPlot(progc, group.by = 'twotype', label = T)
# 
# 
gc2cosm <- ccosg(progc, groups = 'twotype', 
                 assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 100)
gc2cosm <- cbind(gc2cosm$names[1], gc2cosm$scores[1], gc2cosm$names[2], gc2cosm$scores[2])
gcmpcdh <- FindMarkers(progc, ident.1 = 'ant', ident.2 = 'post',
                       group.by = 'twotype', features = c(gc2cosm[,1], gc2cosm[,3]))
gcmpcdh$change = ifelse(abs(gcmpcdh$avg_log2FC) >= 0.5 & abs(gcmpcdh$pct.1 - gcmpcdh$pct.2) >= 0.2, 
                        ifelse(gcmpcdh$avg_log2FC > 0.5 ,'Up','Down'),
                        'Stable')
p <- ggplot(
  gcmpcdh, aes(y = avg_log2FC, x = abs(`pct.1` - `pct.2`), colour=change)) +
  geom_point(alpha=1, size=3.5) +
  scale_color_manual(values=c("#68926a", "#d2dae2","#63669d"))+
  # 
  geom_hline(yintercept=c(-0.5,0.5),lty=1,col="black",lwd=0.8) +
  geom_vline(xintercept = 0.2,lty=1,col="black",lwd=0.8) +
  # 
  labs(y="log2(ant fold change of post)",
       x="abs(ant-post)")+ theme_bw() + 
  # 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()) + 
  geom_text(aes(y = avg_log2FC, x = abs(`pct.1` - `pct.2`), label = rownames(gcmpcdh)))
p
p2 <- gg.gap(p, segments = c(-0.2, 0.2), ylim = c(-2,2))
p2
pdf(file = 'hsa-ceb GC a-p volcano.pdf', height = 5, width = 8)
print(p2)
dev.off()
# 
#
# 的correlation
sppcor <- c('BARHL1', 'PCDH7', 'PTPRK', 'NCAM2', 'LSAMP', 'NRG3', 'TMTC2', 
            'TLX3', 'KCNIP4', 'ADARB2', 'PARM1', 'ASTN2','IGSF21', 'NTF3', 'ITM2C', 'MFAP2')
allspcor <- c()
for(i in sppcor){
  cor1 <- apply(onlygc@assays$RNA@data[sppcor, ], 1,
                function(x) cor.test(x, onlygc@assays$RNA@data[i, ]))
  cor1 <- unlist(lapply(cor1, function(x) x$estimate))
  allspcor <- cbind(allspcor, cor1)
}
write.csv(allspcor, file = 'corr of cerebellum spatial in URL GCP.csv', quote = F)
# ITM2C-NTF3, ASTN2-DCDC1, , 'ZNF804A', , 'DCDC1'
newspcor <- c('Barhl1', 'Lsamp', 'Tmsb4x', 'Sox4', 'Ckb', 'Runx1t1', 'Celf4', 
              'Tlx3', 'Ntf3', 'Mfap2', 'Igfbp2', 'Spsb4', 'Plagl1', 'Hpcal1', 'Ier5', 'Ephb2', 'Mdk', 'Zmiz1')
allspcor <- c()
for(i in toupper(newspcor)){
  cor1 <- apply(onlygc@assays$RNA@data[toupper(newspcor), ], 1,
                function(x) cor.test(x, onlygc@assays$RNA@data[i, ]))
  cor1 <- unlist(lapply(cor1, function(x) x$estimate))
  allspcor <- cbind(allspcor, cor1)
}
write.csv(allspcor, file = 'corr of cerebellum spatial in find mouse URL GCP human.csv', quote = F)
# 


# 
# PRR35 EBF2
pdf(file = 'hsa-ceb HEY1 EBF2.pdf', height = 2, width = 8)
f1 <- FeaturePlot(atatgc, features = c('HEY1', 'EBF2'), blend = T, blend.threshold = 0.1, 
                  cols = c('#eeeeee', '#009f35', '#D40041')) & mockdim & labs(title = element_blank())
print(f1)
dev.off()
mantiprr <- FindMarkers(atatgc, ident.1 = c(14,7,11), logfc.threshold = 0.7, min.diff.pct = 0.3)
pdf(file = 'hsa-ceb PRR35 EBF2 de PRR35.pdf', height = 2, width = 2)
f1 <- FeaturePlot(atatgc, features = c('PRR35'), max.cutoff = 1.5, 
                  cols = c('#eeeeee', '#036eb8')) & mockdim & labs(title = element_blank())
print(f1)
dev.off()
# 
# 

# 用GSEA看看ant-post区别？
# DefaultAssay(progc) <- 'RNA'
# gcforgsea <- FindMarkers(progc, ident.1 = 'ant', ident.2 = 'post',
#                          group.by = 'twotype', logfc.threshold = 0, min.diff.pct = 0)
# write.csv(gcforgsea, file = 'G1 GC a-p logFC for gsea.csv')
# # 
# gcforgsea <- read.csv(file = 'G1 GC a-p logFC for gsea.csv', row.names = 1)
# gcforgsea <- cbind(gcforgsea, rownames(gcforgsea))
# colnames(gcforgsea)[6] <- 'SYMBOL'
# gcid <- bitr(rownames(gcforgsea), OrgDb = 'org.Hs.eg.db', 
#              fromType = 'SYMBOL', toType = c('ENTREZID', 'SYMBOL'))
# gcgsea <- merge(gcforgsea, gcid, by = 'SYMBOL')
# gcgsea <- gcgsea[order(gcgsea$avg_log2FC, decreasing = T), ]
# # 
# ggcclist <- gcgsea$avg_log2FC
# names(ggcclist) <- gcgsea$ENTREZID
# # 
# gcgsego <- gseGO(geneList = ggcclist, OrgDb = org.Hs.eg.db, verbose = T)
# options(clusterProfiler.download.method = "wininet")
# gcgsekegg <- gseKEGG(geneList = ggcclist, organism = "hsa")
# # 
# pdf(file = 'hsa-ceb GC AT-ND gsea.pdf', width = 8, height = 5)
# fprint1 <- gseaplot2(gcgsego, geneSetID = c('GO:0002181', 'GO:0042254'), 
#                      color = c('#00AAAA', '#00AA55'))
# fprint2 <- gseaplot2(gcgsego, geneSetID = c('GO:0007409', 'GO:0048667'), 
#                      color = c('#5555FF', '#B088FF'))
# fprint4 <- gseaplot2(gcgsekegg, geneSetID = c('hsa04015', 'hsa04062', 'hsa04014'), 
#                      color = c('#00AAAA', '#00AA55', '#009FCC'))
# fprint5 <- gseaplot2(gcgsekegg, geneSetID = c('hsa03010'), 
#                      color = c('#00AAAA'))
# print(fprint1)
# print(fprint2)
# print(fprint4)
# print(fprint5)
# dev.off()
# # 

#
# osterior
DefaultAssay(progc) <- 'RNA'
# destiny
post <- subset(progc, twotype == 'post')
DefaultAssay(post) <- 'integrated'
ct <-GetAssayData(object = post)
ct <- ct[VariableFeatures(post),]
ct <- as.ExpressionSet(as.data.frame(t(ct)))
ct$clusters <- post@meta.data$seurat_clusters
ct$protype <- post@meta.data$protype
find_dm_k(nrow(ct)-1) # 1942
dm <- DiffusionMap(ct, k = 1942, n_pcs = 50, verbose = T)
# 
palette(c('#dbeddf', '#b8d4e9', '#84b186', '#d9c6da', '#80325c', '#c2b9ee', '#f8e4db', 
                   '#63669d', '#bab8a9', '#4b5861', '#8f9a8a', '#d1caad', '#916e8c', '#83949e', 
                   '#8f9a8a', '#00709f', '#5c9400', '#b7c2c8', '#edeaf1'))
plot(dm, 1:2, pch = 20, col_by = 'clusters', legend_main = 'Cell stage')
# slingshot
DefaultAssay(post) <- 'RNA'
sce <- as.SingleCellExperiment(post)
DefaultAssay(post) <- 'integrated'
sce2 <- as.SingleCellExperiment(post)
reducedDims(sce) <- reducedDims(sce2); rm(sce2)
FALSE %in% (rownames(sce@int_colData@listData$reducedDims$UMAP) == rownames(dm@eigenvectors))
sce@int_colData@listData$reducedDims$UMAP <- dm@eigenvectors
# 
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP', 
                 start.clus = 1, end.clus = c(7,4,13), approx_points = 45)
# summary(sce$slingPseudotime_1)
lin1 <- getLineages(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP', 
                    start.clus= 1, end.clus = c(7,4,13))
cuv1 <- getCurves(lin1)
plot(reducedDims(sce)[['UMAP']], col=post$seurat_clusters, asp = 1, pch = 16)
lines(SlingshotDataSet(cuv1), lwd = 3, col = 'black')
# 
pdf(file = 'hsa-ceb GCsublineage post desitiny.pdf', width = 7, height = 4)
plot(dm, 1:2, pch = 20, col_by = 'clusters', legend_main = 'Cell stage')
plot(reducedDims(sce)[['UMAP']], col=post$seurat_clusters, asp = 1, pch = 16)
lines(SlingshotDataSet(cuv1), lwd = 3, col = 'black')
dev.off()
# 
trageneplot <- function(sce, pseudo, feature, colors, xrange){
  totalexp <- c()
  for(i in 1:length(pseudo)){
    ttt <- sce[[pseudo[i]]]
    eee <- expddd <- as.matrix(sce@assays@data@listData$counts[feature, ])
    assign(paste0('m',i), as.data.frame(cbind(ttt, eee)) %>% filter(!is.na(eee)))
    traident <- rep(i, nrow(get(paste0('m',i))))
    assign(paste0('m',i), cbind(get(paste0('m',i)), traident))
    assign(rownames(get(paste0('m',i))), seq(1:nrow(get(paste0('m',i)))))
    totalexp <- rbind(totalexp, get(paste0('m',i)))
  }
  totalexp$traident <- as.factor(totalexp$traident)
  # 
  p <- ggplot(totalexp, aes(x = `ttt`, y = `V2`, col = `traident`)) + 
    scale_y_continuous(trans='log1p') + geom_smooth(method = 'loess') +
    xlim(min(xrange), max(xrange)) + labs(title = feature) +
    scale_color_manual(values = colors) + 
    theme_bw()
  return(p)
}
# 
postgene <- c('TLX3', 'NTF3', 'PRPH', 'GALNTL6', 'LINC02516', 'GALNT17', 'PRKCB', 'NRXN3', 'C1QTNF3')
postpseudo <- c('slingPseudotime_1', 'slingPseudotime_2', 'slingPseudotime_3')
postcol <- c('#83949e', '#80325c', '#63669d')
pdf(file = 'hsa-ceb GCsublineage post slingshot.pdf', width = 4, height = 3)
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[1], colors=postcol, xrange = c(0, 0.32))
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[2], colors=postcol, xrange = c(0, 0.32))
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[3], colors=postcol, xrange = c(0, 0.32))
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[4], colors=postcol, xrange = c(0, 0.32))
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[5], colors=postcol, xrange = c(0, 0.32))
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[6], colors=postcol, xrange = c(0, 0.32))
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[7], colors=postcol, xrange = c(0, 0.32))
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[8], colors=postcol, xrange = c(0, 0.32))
trageneplot(sce=sce, pseudo=postpseudo, feature=postgene[9], colors=postcol, xrange = c(0, 0.32))
dev.off()



# expttt <- sce$slingPseudotime_1
# expddd <- as.matrix(sce@assays@data@listData$counts['C1QTNF3', ])
# mmmexp <- as.data.frame(cbind(expttt, expddd)) %>% filter(!is.na(expttt))
# traident <- rep(1, nrow(mmmexp))
# mmmexp <- cbind(mmmexp, traident)
# rownames(mmmexp) <- seq(1:nrow(mmmexp))
# # 
# expttt <- sce$slingPseudotime_3
# expddd <- as.matrix(sce@assays@data@listData$counts['C1QTNF3', ])
# mmmexp2 <- as.data.frame(cbind(expttt, expddd)) %>% filter(!is.na(expttt))
# traident <- rep(2, nrow(mmmexp2))
# mmmexp2 <- cbind(mmmexp2, traident)
# rownames(mmmexp2) <- seq(1:nrow(mmmexp2))
# # 
# totalmmexp <- rbind(mmmexp, mmmexp2)
# totalmmexp$traident <- as.factor(totalmmexp$traident)
# # 
# ggplot(totalmmexp, aes(x = `expttt`, y = `V2`, col = `traident`)) + 
#   scale_y_continuous(trans='log1p') + 
#   geom_smooth(method = 'loess') + scale_color_manual(values = c('black', 'grey'))


## 再看anterior的
DefaultAssay(progc) <- 'RNA'
ant <- subset(progc, twotype == 'ant' & seurat_clusters %in% setdiff(0:18, c(8,18,17)))
DefaultAssay(ant) <- 'integrated'
ct <-GetAssayData(object = ant)
ct <- ct[VariableFeatures(ant),]
ct <- as.ExpressionSet(as.data.frame(t(ct)))
ct$clusters <- ant@meta.data$seurat_clusters
ct$protype <- ant@meta.data$protype
find_dm_k(nrow(ct)-1) # 1942
dm <- DiffusionMap(ct, k = 1942, n_pcs = 50, verbose = T)
# 
palette(c('#dbeddf', '#b8d4e9', '#84b186', '#d9c6da', '#80325c', '#c2b9ee', '#f8e4db', 
                   '#63669d', '#bab8a9', '#4b5861', '#8f9a8a', '#d1caad', '#916e8c', '#83949e', 
                   '#8f9a8a', '#00709f', '#5c9400', '#b7c2c8', '#edeaf1'))
plot(dm, c(1:2), pch = 20, col_by = 'clusters', 
     legend_main = 'Cell stage')
# lingshot
DefaultAssay(ant) <- 'RNA'
sce <- as.SingleCellExperiment(ant)
DefaultAssay(ant) <- 'integrated'
sce2 <- as.SingleCellExperiment(ant)
reducedDims(sce) <- reducedDims(sce2); rm(sce2)
FALSE %in% (rownames(sce@int_colData@listData$reducedDims$UMAP) == rownames(dm@eigenvectors))
sce@int_colData@listData$reducedDims$UMAP <- dm@eigenvectors
# 
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP', 
                 start.clus = 1, end.clus = c(9,10,15), approx_points = 45)
# summary(sce$slingPseudotime_1)
lin1 <- getLineages(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP', 
                    start.clus= 1, end.clus = c(9,10,15))
cuv1 <- getCurves(lin1)
plot(reducedDims(sce)[['UMAP']], col=ant$seurat_clusters, asp = 1, pch = 16)
lines(SlingshotDataSet(cuv1), lwd = 3, col = 'black')
# 
pdf(file = 'hsa-ceb GCsublineage anti desitiny.pdf', width = 7, height = 4)
plot(dm, 1:2, pch = 20, col_by = 'clusters', legend_main = 'Cell stage')
plot(reducedDims(sce)[['UMAP']], col=ant$seurat_clusters, asp = 1, pch = 16)
lines(SlingshotDataSet(cuv1), lwd = 3, col = 'black')
dev.off()
# 
DefaultAssay(ant) <- 'RNA'
mmmant9 <- FindMarkers(ant, ident.1 = 9, ident.2 = 10, 
                       logfc.threshold = 0.7, min.diff.pct = 0.3)
mmmant15 <- FindMarkers(ant, ident.1 = 15, ident.2 = c(9,10), 
                       logfc.threshold = 0.7, min.diff.pct = 0.3)
mmmant10 <- FindMarkers(ant, ident.1 = 10, ident.2 = c(9), 
                        logfc.threshold = 0.7, min.diff.pct = 0.3)
antgene <- c('PCDH7', 'PTPRK', 'PRPH', 'EBF2', 'MEIS2', 'DOK6', 'SLC8A1', 'NR2F2', 'CPNE4')
antpseudo <- c('slingPseudotime_1', 'slingPseudotime_2', 'slingPseudotime_3')
antcol <- c('#8f9a8a', '#4b5861', '#00709f')
# 
pdf(file = 'hsa-ceb GCsublineage antt slingshot.pdf', width = 4, height = 3)
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[1], colors=antcol, xrange = c(0, 0.25))
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[2], colors=antcol, xrange = c(0, 0.25))
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[3], colors=antcol, xrange = c(0, 0.25))
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[4], colors=antcol, xrange = c(0, 0.25))
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[5], colors=antcol, xrange = c(0, 0.25))
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[6], colors=antcol, xrange = c(0, 0.25))
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[7], colors=antcol, xrange = c(0, 0.25))
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[8], colors=antcol, xrange = c(0, 0.25))
trageneplot(sce=sce, pseudo=antpseudo, feature=antgene[9], colors=antcol, xrange = c(0, 0.25))
dev.off()


## 
load(file = 'human GC after cc.Rdata')
# 12,15,10,11
# 
DefaultAssay(endgc) <- 'integrated'
endgc <- FindSubCluster(endgc, cluster = c(4), resolution = 0.4, graph.name = 'integrated_snn')
DimPlot(endgc, group.by = 'sub.cluster', label = T)
VlnPlot(endgc, 'PRPH', group.by = 'sub.cluster')
endgc$sub4 <- endgc$sub.cluster
# 4-2
DefaultAssay(endgc) <- 'integrated'
endgc <- FindSubCluster(endgc, cluster = c(9), resolution = 1, graph.name = 'integrated_snn')
DimPlot(endgc, group.by = 'sub.cluster', label = T)
DefaultAssay(endgc) <- 'RNA'
VlnPlot(endgc, 'PRPH', group.by = 'sub.cluster')
endgc$sub9 <- endgc$sub.cluster
# 9-1/3/5
DefaultAssay(endgc) <- 'integrated'
endgc <- FindSubCluster(endgc, cluster = c(3), resolution = 0.5, graph.name = 'integrated_snn')
DimPlot(endgc, group.by = 'sub.cluster', label = T)
VlnPlot(endgc, 'PRPH', group.by = 'sub.cluster')
endgc$sub3 <- endgc$sub.cluster
# 3-4,3-3
# 
# 
twogc <- as.data.frame(cbind(endgc$sub4, endgc$sub9, endgc$sub3))
twogc <- cbind(twogc, seq(1:nrow(twogc)))
colnames(twogc) <- c('sub4', 'sub9', 'sub3', 'final')
for(i in 1:nrow(twogc)){
  if(TRUE %in% grepl('_', twogc[i, ])){
    twogc[i, 'final'] <- twogc[i, grepl('_', twogc[i, ])]
  } else {
    twogc[i, 'final'] <- twogc[i, 1]
  }
}
DefaultAssay(endgc) <- 'RNA'
# 
endgc$ffftype <- twogc$final
DimPlot(endgc, group.by = 'ffftype', label = T)
endgc$twotype <- endgc$finaltype
# endgc$twotype[which(endgc$twotype %in% c('3','12','4','7','13','14','17_0','17_2','5_0', '1_3', 
#                                          '1_5', '1_7', '1_8', '1_9'))] <- 'post'
# progc$twotype[which(progc$twotype %in% c('0','2','6','8','9','10','11','15','18','16', '17_1', 
#                                          '5_1', '5_2', '1_0','1_1','1_2','1_4','1_6', '1_10', 
#                                          '1_11'))] <- 'ant'
# table(progc$twotype)
# DimPlot(progc, group.by = 'twotype', label = T)
# 
eeegc <- subset(endgc, ffftype %in% c('12','15','10','11','3_3','3_4','4_2','9_1','9_3','9_5'))
for(i in names(table(eeegc$orig.ident))){
  assign(i, subset(eeegc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
rm(eeegc); gc()
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
eeegc <- IntegrateData(anchorset = anchor, dims = 1:30); gc()
eeegc <- ScaleData(eeegc)
eeegc <- RunPCA(eeegc, npcs = 50)
ElbowPlot(object = eeegc, ndims = 40, reduction = "pca") 
DefaultAssay(eeegc) <- 'integrated'
eeegc <- FindNeighbors(eeegc, reduction = "pca", dims = 1:17)
eeegc <- FindClusters(eeegc, resolution = 1)
eeegc <- RunTSNE(eeegc, reduction = "pca", dims = 1:17, verbose = T, )
feee <- DimPlot(eeegc, label = T, label.size = 5)
feee
DefaultAssay(eeegc) <- 'RNA'
# 
mmmeee4 <- FindMarkers(eeegc, ident.1 = 7, logfc.threshold = 0.7, min.diff.pct = 0.2)
# 
neeegc <- subset(eeegc, seurat_clusters %in% setdiff(0:11, c(4,7)))
DefaultAssay(neeegc) <- 'integrated'
neeegc <- ScaleData(neeegc)
neeegc <- RunPCA(neeegc, npcs = 50)
ElbowPlot(object = neeegc, ndims = 40, reduction = "pca") 
DefaultAssay(neeegc) <- 'integrated'
neeegc <- FindNeighbors(neeegc, reduction = "pca", dims = 1:16)
neeegc <- FindClusters(neeegc, resolution = 1)
neeegc <- RunTSNE(neeegc, reduction = "pca", dims = 1:16, verbose = T, )
fnee <- DimPlot(neeegc, label = T, label.size = 5)
fnee
DefaultAssay(neeegc) <- 'RNA'
# 
neeegc <- subset(neeegc, seurat_clusters != 9)
neeegc <- subset(neeegc, seurat_clusters != 6)
neeegc <- subset(neeegc, seurat_clusters != 10)
fnee <- DimPlot(neeegc, label = T, label.size = 5)
fnee
mmmeee6 <- FindMarkers(neeegc, ident.1 = 5, logfc.threshold = 0.5)
mmmeee6 <- FindMarkers(neeegc, ident.1 = 10, ident.2 = 2, logfc.threshold = 0.5)
# 
neeegc$detailtp <- neeegc$seurat_clusters
levels(neeegc$detailtp)[which(levels(neeegc$detailtp) %in% c(2))] <- 'An1'
levels(neeegc$detailtp)[which(levels(neeegc$detailtp) %in% c(1))] <- 'An2'
levels(neeegc$detailtp)[which(levels(neeegc$detailtp) %in% c(3))] <- 'An3'
levels(neeegc$detailtp)[which(levels(neeegc$detailtp) %in% c(8))] <- 'Po3'
levels(neeegc$detailtp)[which(levels(neeegc$detailtp) %in% c(7,0))] <- 'Po1'
levels(neeegc$detailtp)[which(levels(neeegc$detailtp) %in% c(4,5))] <- 'Po2'
table(neeegc$detailtp)
# 
# 
cahedpp <- rownames(endgc)[grepl(x = rownames(endgc), '^PCDH|^CDH|^GALNT|^EFN|^EPH')]
cahedpp <- setdiff(cahedpp, c('PCDH9-AS1', 'PCDH9-AS2', 'PCDH9-AS3', 'PCDH9-AS4'))
pdf(file = 'test.pdf', height = 4, width = 30)
DotPlot(neeegc, features = c('PRPH', cahedpp), group.by = 'detailtp', 
        cluster.idents = F, cols = c('#dddddd', '#63669d')) + RotatedAxis()
dev.off()
# 
cahedpp <- rownames(endgc)[grepl(x = rownames(endgc), '^CELSR|^CLSTN|^DCHS|^DSC|^DSG|^FAT')]
cahedpp <- setdiff(cahedpp, c('CLSTN2-AS1', 'DSG1-AS1', 'DSG2-AS1', 'DSCAM-AS1'))
pdf(file = 'test.pdf', height = 4, width = 10)
DotPlot(neeegc, features = c('PRPH', cahedpp), group.by = 'detailtp', 
        cluster.idents = F, cols = c('#dddddd', '#63669d')) + RotatedAxis()
dev.off()
# 
cahedpp <- rownames(endgc)[grepl(x = rownames(endgc), '^ALG|^B3GAL|^B4GAL|^CHSD|^FUT|^MAN|^MGAT')]
cahedpp <- setdiff(cahedpp, c('B4GALT4-AS1', 'MANEA-DT', 'B4GALT1-AS1', 'MAN1B1-DT', 'ALG9-IT1', 'B3GALT5-AS1', 'ALG13-AS1'))
pdf(file = 'test.pdf', height = 4, width = 20)
DotPlot(neeegc, features = c('PRPH', cahedpp), group.by = 'detailtp', 
        cluster.idents = F, cols = c('#dddddd', '#63669d')) + RotatedAxis()
dev.off()


# 
cadherin <- c('CDH9', 'CDH12', 'CDH18', 'CDH20', 'PCDH7', 'PCDH15', 'PCDH17', 'CLSTN2', 'DSCAM', 
              'EFNA3', 'EFNB2', 'EPHA3', 'EPHA5', 'EPHA8', 'EPHB1', 'EPHB6', 
              'GALNT9', 'GALNTL6', 'GALNT12', 'GALNT14', 'GALNT15', 'GALNT17', 'GALNT18', 
              'FUT8', 'FUT9', 'MAN1A1', 'MANBAL')
# AL157895.1, AC073225.1, 
neeegc$detailtp<-factor(neeegc$detailtp,levels = c('An1','An2','An3','Po1','Po2','Po3'))
pdf(file = 'hsa-ceb GC end cadherin+eph+NOglyco.pdf', height = 3, width = 10)
DotPlot(neeegc, features = c(cadherin), group.by = 'detailtp', 
        cluster.idents = F, cols = c('#dddddd', 'navy')) + RotatedAxis()
dev.off()
# 
#
pdf(file = 'hsa-ceb GC end umap.pdf', height = 2, width = 2)
print(DimPlot(neeegc, group.by = 'detailtp', cols = c('#8f9a8a', '#4b5861', '#00709f', '#83949e', '#80325c', '#63669d')) & 
        mockdim & labs(title = element_blank()))
dev.off()














#
excigrn <- read.csv(file = "E:/cerebellum/GCS+UBC grn.csv")
# 
updownadj <- p.adjust(excigrn$Up_Down_Pvalue, method = 'BH')
creupadj <- p.adjust(excigrn$CRE_Up_Pvalue, method = 'BH')
excigrn <- cbind(excigrn, updownadj, creupadj)
egrn005 <- excigrn %>% filter(updownadj <0.05 & creupadj <0.05)
# 
specigc <- rownames(gcmpcdh %>% filter(change != 'Stable'))
gcforgrn <- c(specigc, 'BARHL1', 'EPHB2', 'NRG3')
# gcforgrn <- gcmpcdh[which(rownames(gcmpcdh) %in% c(egrn005$Upstream_TF, egrn005$Downstream_Gene)), ]
gcforgrn <- gcforgrn[which(gcforgrn %in% c(egrn005$Upstream_TF, egrn005$Downstream_Gene))]
gcgcgc <- egrn005 %>% filter(Downstream_Gene %in% gcforgrn)
gcgcgc <- gcgcgc[, c("Upstream_TF", "Downstream_Gene")]
gcgcgc <- gcgcgc[!duplicated(gcgcgc), ]
write.csv(gcgcgc, file = 'grn part of two gc net.csv', quote = F)
# 
gcforgrn <- gcforgrn[which(gcforgrn %in% c(excigrn$Upstream_TF, excigrn$Downstream_Gene))]
gcgcgc <- excigrn %>% filter(Downstream_Gene %in% gcforgrn)
gcgcgc <- gcgcgc[, c("Upstream_TF", "Downstream_Gene")]
gcgcgc <- gcgcgc[!duplicated(gcgcgc), ]
write.csv(gcgcgc, file = 'grn part of two gc net 2.csv', quote = F)
# 
gcforgrn <- gcmpcdh[which(rownames(gcmpcdh) %in% c(excigrn$Upstream_TF, excigrn$Downstream_Gene)), ]
gcgcgc <- egrn005 %>% filter(Downstream_Gene %in% rownames(gcforgrn))
gcgcgc <- gcgcgc[, c("Upstream_TF", "Downstream_Gene")]
gcgcgc <- gcgcgc[!duplicated(gcgcgc), ]
write.csv(gcgcgc, file = 'grn part of two gc net 3.csv', quote = F)

endgc <- subset(endgc, seurat_clusters != 15)
fend <- DimPlot(endgc, label = T, label.size = 5)
fend
eeemarker <- ccosg(endgc, groups = 'seurat_clusters', 
                   assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 5)
eeemarker <- as.data.frame(eeemarker[[1]])
# 

endm1 <- FindMarkers(endgc, ident.1 = c(8,10,0,1,3,7,16,14,12), 
                     logfc.threshold = 0.5, min.diff.pct = 0.2)
endtlx1 <- FindMarkers(endgc, ident.1 = c(11,13), ident.2 = c(2,4,6,8,9,5),
                       logfc.threshold = 0.5, min.diff.pct = 0.2)
endtlx2 <- FindMarkers(endgc, ident.1 = c(6,9), ident.2 = c(2,4,5,11,13,8),
                       logfc.threshold = 0.5, min.diff.pct = 0.2)
endbar1 <- FindMarkers(endgc, ident.1 = c(14), ident.2 = c(8,10,0,1,3,7,16,12),
                       logfc.threshold = 0.5, min.diff.pct = 0.2)
endbar2 <- FindMarkers(endgc, ident.1 = c(3), ident.2 = c(8,10,0,1,14,7,16,12),
                       logfc.threshold = 0.5, min.diff.pct = 0.2)
endbar3 <- FindMarkers(endgc, ident.1 = c(1), ident.2 = c(8,10,0,3,14,7,16,12),
                       logfc.threshold = 0.5, min.diff.pct = 0.2)
# 
endgc$celltype <- endgc$seurat_clusters
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(2,4,5))] <- 'pGC-pre'
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(0,8,10,7,16))] <- 'aGC-pre'
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(13))] <- 'pGC-1'
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(11))] <- 'pGC-2'
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(6,9))] <- 'pGC-3'
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(3))] <- 'aGC-1'
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(14))] <- 'aGC-2'
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(12))] <- 'aGC-3'
levels(endgc$celltype)[which(levels(endgc$celltype) %in% c(1))] <- 'aGC-4'
table(endgc$celltype)

# 
load(file = 'human GC after cc.Rdata')
# 
DefaultAssay(endgc) <- 'RNA'
pdf(file = 'hsa-ceb GC end umap.pdf', width = 3, height = 3)
f1 <- DimPlot(endgc, label = F, group.by = 'celltype', 
              cols = c('#63669d', '#779ac2', '#68926a', '#916e8c', '#4b5861', '#8f9a8a', '#c2b9ee', '#d1caad', '#9ba5c0')) & 
  mockdim & labs(title = element_blank())
print(f1)
dev.off()
# 
ppro <- FindMarkers(endgc, ident.1 = c(2,4,5), ident.2 = c(11,13,6,6), only.pos = T, 
                    min.diff.pct = 0.3, logfc.threshold = 0.5)
apro <- FindMarkers(endgc, ident.1 = c(8,10,0,7,16), ident.2 = c(3,14,12,1), only.pos = T, 
                    min.diff.pct = 0.2, logfc.threshold = 0.3)
# engc7 <- FindMarkers(endgc, ident.1 = c(7), ident.2 = c(3,14,12,1), only.pos = T, 
#                      min.diff.pct = 0.2, logfc.threshold = 0.3)
diffgc <- c('ADAMTS6', 'TLL1', 'TMEM132B', 'EPHA3', 'NR2F2', 'AC015574.1', 
            'SEMA3C', 'ALCAM', 'C11orf96', 'MEIS2', 'EBF2', 
            'EYA2', 'RAB26', 'ITGBL1', 
            'GALNTL6','STC1', 'BRINP2', 'NTM', 'RPRM', 'CSMD3', 'DOK6', 
            'C1QTNF3', 'NRXN3', 'PPP1R1B', 'BHLHE41', 'GMDS')
endgc$celltype<-factor(endgc$celltype,levels = c('aGC-pre','aGC-1','aGC-2','aGC-3','aGC-4','pGC-pre','pGC-1','pGC-2','pGC-3'))
pdf(file = 'hsa-ceb GC end dotplot.pdf', height = 4, width = 10)
DotPlot(endgc, features = diffgc, group.by = 'celltype', 
        cluster.idents = F, cols = c('#dddddd', '#63669d')) + RotatedAxis()
dev.off()
# 
# 
cor1 <- apply(endgc@assays$RNA@data[rownames(endgc), ], 1, 
              function(x) cor.test(x, endgc@assays$RNA@data['PRPH', ]))
cor1 <- unlist(lapply(cor1, function(x) x$estimate))
cor1 <- as.data.frame(cbind(names(cor1), cor1))


load(file = 'human GC before cc.Rdata')
gcprotp <- progc$celltype %>% as.data.frame()
rm(progc)
# 
load(file = 'human GC after cc.Rdata')
gcendtp <- endgc$celltype %>% as.data.frame()
rm(endgc)
# 
onlygc <- subset(onlygc, seurat_clusters == 23)
onlygc$celltype <- onlygc$seurat_clusters
levels(onlygc$celltype)[which(levels(onlygc$celltype) %in% c(23))] <- 'GCstem'
gcstemtp <- onlygc$celltype %>% as.data.frame()
rm(onlygc)
# 
load(file = 'human UBC two lineage.Rdata')
ubctp <- ubwithsimi$celltype %>% as.data.frame()
rm(ubwithsimi)
# 
load(file = 'human pkc+pro+in.Rdata')
ppicell <- subset(ppicell, seurat_clusters %in% setdiff(0:17, 16))
ppicell$celltype <- ppicell$seurat_clusters
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(1))] <- 'BG'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(8,10))] <- 'G2M'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(11))] <- 'GCP'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(5,6))] <- 'VZP'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(12))] <- 'EPP'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(15,2,3,4,7))] <- 'PKC'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(0,13,14,9))] <- 'IN'
levels(ppicell$celltype)[which(levels(ppicell$celltype) %in% c(17))] <- 'OPC'
table(ppicell$celltype)
prointtp <- ppicell$celltype %>% as.data.frame()
rm(ppicell)
# 
load(file = 'human pkc.Rdata')
table(purkinje$celltype)
purktp <- purkinje$celltype %>% as.data.frame()
rm(purkinje)
# 
# 
gc()
c('gcprotp', 'gcendtp', 'gcstemtp', 'ubctp', 'prointtp', 'purktp')
gcprotp$. <- as.character(gcprotp$.) ; gcprotp$cellid <- rownames(gcprotp)
colnames(gcprotp)[1] <- 'gcpro'
gcendtp$. <- as.character(gcendtp$.); gcendtp$cellid <- rownames(gcendtp)
colnames(gcendtp)[1] <- 'gcend'
gcstemtp$. <- as.character(gcstemtp$.); gcstemtp$cellid <- rownames(gcstemtp)
colnames(gcstemtp)[1] <- 'gcstem'
ubctp$. <- as.character(ubctp$.); ubctp$cellid <- rownames(ubctp)
colnames(ubctp)[1] <- 'ubc'
prointtp$. <- as.character(prointtp$.); prointtp$cellid <- rownames(prointtp)
colnames(prointtp)[1] <- 'proint'
purktp$. <- as.character(purktp$.); purktp$cellid <- rownames(purktp)
colnames(purktp)[1] <- 'purkin'
# 
intersect(gcprotp$cellid, gcendtp$cellid)
intersect(ubctp$cellid, gcprotp$cellid)
intersect(ubctp$cellid, prointtp$cellid) # +
intersect(ubctp$cellid, purktp$cellid)
intersect(prointtp$cellid, purktp$cellid) # +
intersect(gcprotp$cellid, prointtp$cellid)
intersect(gcstemtp$cellid, prointtp$cellid) #+

# 
totalctp <- merge(gcprotp, gcendtp, by = 'cellid', all = T)
totalctp <- merge(totalctp, gcstemtp, by = 'cellid', all = T)
totalctp <- merge(totalctp, ubctp, by = 'cellid', all = T)
totalctp <- merge(totalctp, prointtp, by = 'cellid', all = T)
totalctp <- merge(totalctp, purktp, by = 'cellid', all = T)
#
totalctp[is.na(totalctp)] <- 0
# fififi <- c()
# for(i in 1:nrow(totalctp)){
#   if(length(totalctp[grepl(pattern = 0, totalctp[i, ]), ]) >= 2)
# }
totalctp$final <- seq(1,nrow(totalctp))
totalctp$final <- totalctp$proint
for(i in 1:nrow(totalctp)){
  if(totalctp[i, "purkin"] != 0){
    totalctp[i, "final"] <- totalctp[i, "purkin"]
  }
  if(totalctp[i, "ubc"] != 0){
    totalctp[i, "final"] <- totalctp[i, "ubc"]
  }
  if(totalctp[i, "gcend"] != 0){
    totalctp[i, "final"] <- totalctp[i, "gcend"]
  }
  if(totalctp[i, "gcpro"] != 0){
    totalctp[i, "final"] <- totalctp[i, "gcpro"]
  }
  if(totalctp[i, "gcstem"] != 0){
    totalctp[i, "final"] <- totalctp[i, "gcstem"]
  }
}
cellforchat <- totalctp %>% filter(final != 'GCP' & final != 'G2M')
write.csv(cellforchat, file = 'cell for cellchat with type.csv')

#
load(file = 'cebclean.Rdata')
cellforchat <- read.csv(file = 'cell for cellchat with type.csv')
rownames(cellforchat) <- cellforchat$cellid
cebclean <- cebclean[, cellforchat$cellid]
gc()
cellforchat <- cellforchat[names(cebclean$nCount_RNA), ]
cebclean$finaltype <- cellforchat$final
# 
# cellchat
datainput <- GetAssayData(cebclean, assay = 'RNA', slot = 'data')
metadt <- cebclean@meta.data
rm(cebclean, cellforchat)
cellcobj <- createCellChat(object = datainput, meta = metadt, group.by = 'finaltype')
levels(cellcobj@idents)
# 
cellcobj@DB <- CellChatDB.human
# 
cellcobj <- subsetData(cellcobj)
cellcobj <- identifyOverExpressedGenes(cellcobj, thresh.pc = 0, thresh.fc = 0.2, 
                                       thresh.p = 0.01)
gc()
cellcobj <- identifyOverExpressedInteractions(cellcobj)
# 
cellcobj <- computeCommunProb(object = cellcobj, type = 'triMean', trim = NULL, 
                              nboot = 100) # 15min 左右
cellcobj <- filterCommunication(cellcobj, min.cells = 10)
#
dfnet <- subsetCommunication(cellcobj)
dfnetpath <- subsetCommunication(cellcobj, slot.name = 'netP')
# 
cellcobj <- computeCommunProbPathway(object = cellcobj)
cellcobj <- aggregateNet(cellcobj)
pdf(file = 'all signaling in total.pdf', width = 8, height = 6)
for(i in unique(dfnet$pathway_name)){
  fpri1 <- netVisual_heatmap(cellcobj, signaling = i, color.heatmap = 'Reds')
  print(fpri1)
  fpri2 <- netAnalysis_contribution(cellcobj, signaling = i)
  print(fpri2)
}
dev.off()
# 
netVisual_aggregate(cellcobj, signaling = 'SEMA6', layout = 'hierarchy',
                    vertex.receiver = seq(1,10))
# 
netVisual_bubble(cellcobj, sources.use = 1:6, targets.use = 16:20, remove.isolate = FALSE)
#





ifdifspe <- c('AC010857.1','LINC02381','TULP3','AP002784.1',
              'AC104389.5','CALY','TRH','FGF8','SERPINH1','SLC6A13','CXCL12','ITIH5', 
              'DCN','ATP1A2','IGF2','CYP1B1','PMP22','PTGDS','ETV4','SACS',
              'LGR4','PDZD7',
              'RPS17P8','PPPIR17','AL450988.1',
              'AL162388.2','NEFL','PLSCR2','KIAA1755','RORB',
              'PRR35','NTRK3','COL1A1')
pdf(file = 'test.pdf')
for(i in intersect(rownames(spgw13), ifdifspe)){
  print(sFeatureplot(spgw13, i))
}
dev.off()


FeaturePlot(onlygc, c('SFRP1', 'ATOH1', 'PRPH', 'TOP2A'))



















