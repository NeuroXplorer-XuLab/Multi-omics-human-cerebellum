load(file = 'cebclean.Rdata')

timepoint_conversion <- list(
  'E083' = 'E102',
  'E075' = 'E083',
  'E078' = 'E108',
  'E085' = 'E117',
  'E101' = 'E101',
  'E104' = 'E093'
)

cebclean$timepoint <- as.character(cebclean$orig.ident)
cebclean$timepoint <- unlist(lapply(cebclean$timepoint, function(x) ifelse(x %in% names(timepoint_conversion), timepoint_conversion[[x]], x)))

DefaultAssay(cebclean) <- 'RNA'
cebclean$celltype <- as.factor(cebclean$seurat_clusters)

celltype_conversion <- list(
  GCs-1 = c(2, 1, 10, 6),
  GCs-2 = c(8, 5, 3, 4, 12, 7, 19, 24, 15, 14),
  UBCs = c(17, 9, 0),
  PROs = c(11, 18),
  INs = c(16),
  PKCs = c(13)
)

for (new_label in names(celltype_conversion)) {
  idx <- which(levels(cebclean$celltype) %in% celltype_conversion[[new_label]])
  levels(cebclean$celltype)[idx] <- new_label
}

fint <- DimPlot(cebclean, group.by = 'celltype', shuffle = T, 
                cols = c('#84b186', '#b8d4e9', '#8eb0c9', '#9ba5c0', '#916e8c', '#63669d')) & 
  mockdim & labs(title = element_blank())
fint
pdf('hsa-ceb main umap.pdf', width = 5, height = 5)
print(fint)
dev.off()
f1 <- DimPlot(cebclean, label = F, group.by = 'timepoint', shuffle = T,  
              cols = c('#68926a', '#a9d6af', '#b8d4e9', '#8eb0c9', '#916e8c', '#80325c')) & 
  mockdim & labs(title = element_blank())
f1
pdf('hsa-ceb main diff time.pdf', height = 7, width = 7)
print(f1)
dev.off()

mainratio <- aggregate(cebclean@meta.data, by = list(cebclean@meta.data$timepoint, cebclean@meta.data$celltype), FUN=table)
write.csv(file = 'hsa-ceb main ratio.csv', mainratio[, 1:3], quote = F)

bigheatmarker <- ccosg(cebclean, groups = 'celltype', 
                       assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 15)
heatgene <- c(bigheatmarker[[1]][, 1], bigheatmarker[[1]][, 2], bigheatmarker[[1]][, 3], 
              bigheatmarker[[1]][, 4], bigheatmarker[[1]][, 5], bigheatmarker[[1]][, 6])

for(i in 1:length(summary(cebclean$celltype))){
  assign(paste0('heat', i), sample(names(cebclean$celltype[which(cebclean$celltype == names(summary(cebclean$celltype)[i]))]), 
                                   summary(cebclean$celltype)[i]))
}
summary(cebclean$celltype)
DefaultAssay(cebclean) <- 'RNA'
cebclean <- ScaleData(cebclean, features = heatgene)
cebclean$celltype <- factor(cebclean$celltype, 
                            levels = c('UBCs', 'DCNs', 'GCs', 'PROs', 'PKCs', 'INs')) 
f14 <- DoHeatmap(object = cebclean, features = heatgene, 
                 cells = c(heat1, heat2, heat3, heat4, heat5, heat6), 
                 draw.lines = F, group.by = 'celltype', disp.min = 0, disp.max =2, slot = 'data') + 
  theme(axis.text.y = element_text(size = 250/length(heatgene))) + 
  scale_color_gradientn(colours = c("#f0eeec", '#b8d4e9', '#63669d'), aesthetics = "fill")
pdf(file = 'hsa-ceb main heatmap.pdf', width = 10, height = 6)
print(f14)
dev.off()
# umap featureplot
mainlist <- c('EOMES', 'LMX1A', 'PCDH7', 'PCSK2', 'TLX3', 'RASD1', 
              'SOX2', 'CYP26B1', 'CA8', 'WNT7B', 'PAX2', 'GAD2', 'HES1', 'ATOH1', 'PRPH', 
              'NHLH1', 'NRN1', 'TOP2A', 'MKI67')
for(i in mainlist){
  png(filename = paste0('main ', i, '.png'), width = 800, height = 800)
  f1 <- FeaturePlot(cebclean, i, cols = c('#dddddd', 'navy')) & mockdim & 
    labs(title = element_blank())
  print(f1)
  dev.off()
}

DefaultAssay(cebclean) <- 'integrated'
cebclean <-  FindSubCluster(cebclean, graph.name = 'integrated_snn', cluster = 18, resolution = 0.5)
DimPlot(cebclean, group.by = 'sub.cluster', label = T)
DefaultAssay(cebclean) <- 'RNA'

p180 <- FindMarkers(cebclean, ident.1 = '18_0', ident.2 = c('18_1', '18_2', '18_3'), 
                    group.by = 'sub.cluster', min.diff.pct = 0.3, logfc.threshold = 0.7)
p181 <- FindMarkers(cebclean, ident.1 = '18_1', ident.2 = c('18_0', '18_2', '18_3'), 
                    group.by = 'sub.cluster', min.diff.pct = 0.3, logfc.threshold = 0.7)
p182 <- FindMarkers(cebclean, ident.1 = '18_2', ident.2 = c('18_0', '18_1', '18_3'), 
                    group.by = 'sub.cluster', min.diff.pct = 0.3, logfc.threshold = 0.7)
p183 <- FindMarkers(cebclean, ident.1 = '18_3', ident.2 = c('18_0', '18_1', '18_2'), 
                    group.by = 'sub.cluster', min.diff.pct = 0.3, logfc.threshold = 0.7)

ppicell <- subset(cebclean, sub.cluster %in% c('11','18_2', '18_3','16', '13'))

DefaultAssay(ppicell) <- 'RNA'
for(i in names(table(ppicell$orig.ident))){
  assign(i, subset(ppicell, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
gc()
ppicell <- IntegrateData(anchorset = anchor, dims = 1:20)

DefaultAssay(ppicell) <- 'RNA'
ppicell <- CellCycleScoring(ppicell, s.features = cc.genes.updated.2019$s.genes, 
                            g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
ppicell$CC.Difference <- ppicell$S.Score - ppicell$G2M.Score
DefaultAssay(ppicell) <- 'integrated'
ppicell <- ScaleData(ppicell, vars.to.regress = c('S.Score', 'G2M.Score'))

gc()

ppicell <- RunPCA(ppicell, npcs = 50)
ElbowPlot(object = ppicell, ndims = 40, reduction = "pca") 
DefaultAssay(ppicell) <- 'integrated'
ppicell <- FindNeighbors(ppicell, reduction = "pca", dims = 1:9)
ppicell <- FindClusters(ppicell, resolution = 1)
ppicell <- RunUMAP(ppicell, reduction = "pca", dims = 1:9)
fppi <- DimPlot(ppicell, label = T, label.size = 5)
fppi
DefaultAssay(ppicell) <- 'RNA'
 
VlnPlot(ppicell, 'nFeature_RNA')
ppi6 <- FindMarkers(ppicell, ident.1 = 6, logfc.threshold = 0.5)
ppi11 <- FindMarkers(ppicell, ident.1 = 11, logfc.threshold = 0.5)

ppicell <- subset(ppicell, seurat_clusters %in% setdiff(0:15, c(6,11)))
 
DefaultAssay(ppicell) <- 'RNA'
for(i in names(table(ppicell$orig.ident))){
  assign(i, subset(ppicell, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
gc()
ppicell <- IntegrateData(anchorset = anchor, dims = 1:20)
ppicell <- ScaleData(ppicell)
ppicell <- RunPCA(ppicell, npcs = 50)
ElbowPlot(object = ppicell, ndims = 40, reduction = "pca") 
DefaultAssay(ppicell) <- 'integrated'
ppicell <- FindNeighbors(ppicell, reduction = "pca", dims = 1:8)
ppicell <- FindClusters(ppicell, resolution = 1)
ppicell <- RunUMAP(ppicell, reduction = "pca", dims = 1:8)
fppi <- DimPlot(ppicell, label = T, label.size = 5)
fppi
DefaultAssay(ppicell) <- 'RNA'
save(ppicell, file = 'human pkc+pro+in.Rdata')