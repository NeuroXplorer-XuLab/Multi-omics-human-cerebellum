mmue16 <- readRDS(file = "E:/cerebellum/GSE118068_RAW_rdata/GSM3318002_E16.Rdata")
mmue18 <- readRDS(file = "E:/cerebellum/GSE118068_RAW_rdata/GSM3318003_E18.Rdata")
mmup0 <- readRDS(file = "E:/cerebellum/GSE118068_RAW_rdata/GSM3318004_P0.Rdata")
mmup5 <- readRDS(file = "E:/cerebellum/GSE118068_RAW_rdata/GSM3318005_P5.Rdata")
mmup7 <- readRDS(file = "E:/cerebellum/GSE118068_RAW_rdata/GSM3318006_P7.Rdata")
# 
mgenelist <- c('Car8', 'Pcp4', 'Sox2', 'Gfap', 'Atoh1', 'Eomes', 'Pax2', 'Cx3cr1', 'Cldn5')
VlnPlot(mmue18, 'nFeature_RNA')
FeatureScatter(mmue18, 'nCount_RNA', 'nFeature_RNA')
mmue18 <- subset(mmue18, nFeature_RNA > 500)
mmue18 <- NormalizeData(mmue18, normalization.method = 'LogNormalize')
mmue18 <- FindVariableFeatures(mmue18, selection.method = "vst", nfeatures = 2000)
mmue18 <- ScaleData(mmue18)
mmue18 <- RunPCA(mmue18, npcs = 50)
ElbowPlot(mmue18, ndims = 40, reduction = "pca")
mmue18 <- FindNeighbors(mmue18, reduction = "pca", dims = 1:15)
mmue18 <- FindClusters(mmue18, resolution = 1)
mmue18 <- RunUMAP(mmue18, reduction = "pca", dims = 1:15)
fe18 <- DimPlot(mmue18, label = T, label.size = 6)
fe18
# 16-endo, 17,
# 
VlnPlot(mmup0, 'nFeature_RNA')
FeatureScatter(mmup0, 'nCount_RNA', 'nFeature_RNA')
mmup0 <- subset(mmup0, nFeature_RNA > 500)
mmup0 <- NormalizeData(mmup0, normalization.method = 'LogNormalize')
mmup0 <- FindVariableFeatures(mmup0, selection.method = "vst", nfeatures = 2000)
mmup0 <- ScaleData(mmup0)
mmup0 <- RunPCA(mmup0, npcs = 50)
ElbowPlot(mmup0, ndims = 40, reduction = "pca")
mmup0 <- FindNeighbors(mmup0, reduction = "pca", dims = 1:19)
mmup0 <- FindClusters(mmup0, resolution = 1)
mmup0 <- RunUMAP(mmup0, reduction = "pca", dims = 1:19)
fpp0 <- DimPlot(mmup0, label = T, label.size = 6)
fpp0
# 7/18 pkc, 3/16/19 in, 13/5/6 pro, 
m142 <- FindMarkers(mmup0, ident.1 = c(2,14), only.pos = T, 
                    logfc.threshold = 0.7, min.pct = 0.3)
# 
VlnPlot(mmup5, 'nFeature_RNA')
FeatureScatter(mmup5, 'nCount_RNA', 'nFeature_RNA')
mmup5 <- subset(mmup5, nFeature_RNA > 500)
mmup5 <- NormalizeData(mmup5, normalization.method = 'LogNormalize')
mmup5 <- FindVariableFeatures(mmup5, selection.method = "vst", nfeatures = 2000)
mmup5 <- ScaleData(mmup5)
mmup5 <- RunPCA(mmup5, npcs = 50)
ElbowPlot(mmup5, ndims = 40, reduction = "pca")
mmup5 <- FindNeighbors(mmup5, reduction = "pca", dims = 1:24)
mmup5 <- FindClusters(mmup5, resolution = 1)
mmup5 <- RunUMAP(mmup5, reduction = "pca", dims = 1:24)
fpp5 <- DimPlot(mmup5, label = T, label.size = 6)
fpp5
# 17 ubc, 9/11 in, 10 pro, 12/13 ast, 
m14 <- FindMarkers(mmup5, ident.1 = c(14), only.pos = T, 
                    logfc.threshold = 0.7, min.pct = 0.3)
# 
mean(mmue16$nFeature_RNA)
mmue16 <- subset(mmue16, nFeature_RNA > 500)
mmue16 <- NormalizeData(mmue16, normalization.method = 'LogNormalize')
mmue16 <- FindVariableFeatures(mmue16, selection.method = "vst", nfeatures = 2000)
mmue16 <- ScaleData(mmue16)
mmue16 <- RunPCA(mmue16, npcs = 50)
ElbowPlot(mmue16, ndims = 40, reduction = "pca")
mmue16 <- FindNeighbors(mmue16, reduction = "pca", dims = 1:17)
mmue16 <- FindClusters(mmue16, resolution = 1)
mmue16 <- RunUMAP(mmue16, reduction = "pca", dims = 1:17)
fe16 <- DimPlot(mmue16, label = T, label.size = 6)
fe16
# 2 pkc, 
# 
mean(mmup7$nFeature_RNA)
mmup7 <- subset(mmup7, nFeature_RNA > 500)
mmup7 <- NormalizeData(mmup7, normalization.method = 'LogNormalize')
mmup7 <- FindVariableFeatures(mmup7, selection.method = "vst", nfeatures = 2000)
mmup7 <- ScaleData(mmup7)
mmup7 <- RunPCA(mmup7, npcs = 50)
ElbowPlot(mmup7, ndims = 40, reduction = "pca")
mmup7 <- FindNeighbors(mmup7, reduction = "pca", dims = 1:24)
mmup7 <- FindClusters(mmup7, resolution = 1)
mmup7 <- RunUMAP(mmup7, reduction = "pca", dims = 1:24)
fpp7 <- DimPlot(mmup7, label = T, label.size = 6)
fpp7
# 8 pkc, 13/4 in
# 
me16 <- subset(mmue16, seurat_clusters %in% c(1,11,4,6,7,0,16,3))
me16 <- RenameCells(me16, add.cell.id = 'e16')
me18 <- subset(mmue18, seurat_clusters %in% c(2,6,3,12,5,0,13,1,4))
me18 <- RenameCells(me18, add.cell.id = 'e18')
mep0 <- subset(mmup0, seurat_clusters %in% c(1,4,9,11,0,8,10,12))
mep0 <- RenameCells(mep0, add.cell.id = 'p0')
mep5 <- subset(mmup5, seurat_clusters %in% c(17,3,0,1,7,2,6,4,5,8))
mep5 <- RenameCells(mep5, add.cell.id = 'p5')
mep7 <- subset(mmup7, seurat_clusters %in% c(0,2,3,10,1,5,6,7,9))
mep7 <- RenameCells(mep7, add.cell.id = 'p7')
# 
mouselist <- c('me16', 'me18', 'mep0', 'mep5', 'mep7')
for(i in 1:5){
  assign(mouselist[i], NormalizeData(get(mouselist[i])))
  assign(mouselist[i], FindVariableFeatures(get(mouselist[i]), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(me16, me18, mep0, mep5, mep7), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
gc(); rm(me16, me18, mep0, mep5, mep7)
mouseceb <- IntegrateData(anchorset = anchor, dims = 1:20)
gc()
mouseceb <- ScaleData(mouseceb)
DefaultAssay(mouseceb) <- 'integrated'
mouseceb <- RunPCA(mouseceb, npcs = 50)
ElbowPlot(object = mouseceb, ndims = 40, reduction = "pca") 
DefaultAssay(mouseceb) <- 'integrated'
mouseceb <- FindNeighbors(mouseceb, reduction = "pca", dims = 1:11)
mouseceb <- FindClusters(mouseceb, resolution = 1)
mouseceb <- RunUMAP(mouseceb, reduction = "pca", dims = 1:11)
fmceb <- DimPlot(mouseceb, label = T, label.size = 6)
fmceb
DefaultAssay(mouseceb) <- 'RNA'
nnmouse <- mouseceb


# 
# cross species PRPH
proantipost <- FindMarkers(progc, ident.1 = 'pro', ident.2 = 'afcc', group.by = 'protype', 
                           min.diff.pct = 0.3, logfc.threshold = 0.7)
pdf(file = 'hsa-mmu BOC PRPH diff.pdf', height = 3, width = 16)
print(FeaturePlot(nnmouse, features = c('Boc', 'Prph'), blend = T, blend.threshold = 0.1, 
                  cols = c('#eeeeee', '#036eb8', '#D40041')) & mockdim)
print(FeaturePlot(progc, features = c('BOC', 'PRPH'), blend = T, blend.threshold = 0.1, 
                  cols = c('#eeeeee', '#036eb8', '#D40041')) & mockdim)
dev.off()
# corr
sppcor <- c('BARHL1', 'PCDH7', 'PTPRK', 'NCAM2', 'LSAMP', 'NRG3', 'TMTC2', 
            'TLX3', 'KCNIP4', 'ADARB2', 'PARM1', 'ASTN2','IGSF21', 'NTF3', 'ITM2C', 'MFAP2')
msppcor <- stringr::str_to_title(tolower(sppcor))
allspcor <- c()
for(i in msppcor){
  cor1 <- apply(nnmouse@assays$RNA@data[msppcor, ], 1,
                function(x) cor.test(x, nnmouse@assays$RNA@data[i, ]))
  cor1 <- unlist(lapply(cor1, function(x) x$estimate))
  allspcor <- cbind(allspcor, cor1)
}
write.csv(allspcor, file = 'corr of cerebellum spatial in mouse URL GCP.csv', quote = F)
# 
# 
mousecebgc <- subset(mouseceb, seurat_clusters %in% c(10,4,22,13,17,1,3,23,2,16,20))
allspcor <- c()
for(i in msppcor){
  cor1 <- apply(mousecebgc@assays$RNA@data[msppcor, ], 1,
                function(x) cor.test(x, mousecebgc@assays$RNA@data[i, ]))
  cor1 <- unlist(lapply(cor1, function(x) x$estimate))
  allspcor <- cbind(allspcor, cor1)
}
write.csv(allspcor, file = 'corr of cerebellum spatial in mouse low URL GCP.csv', quote = F)
# 
# 
mouseexp <- AverageExpression(nnmouse, assays = 'RNA', group.by = 'orig.ident')[[1]]
mouseexp <- as.data.frame(mouseexp)
mouseexp <- mouseexp %>% filter(GSM3318003_E18>0.7 | GSM3318004_P0>0.7)
cor1 <- apply(nnmouse@assays$RNA@data[rownames(mouseexp), ], 1, 
              function(x) cor.test(x, nnmouse@assays$RNA@data['Barhl1', ]))
cor1 <- unlist(lapply(cor1, function(x) x$estimate))
cor1 <- as.data.frame(cbind(names(cor1), cor1))
# 
newspcor <- c('Barhl1', 'Lsamp', 'Tmsb4x', 'Sox4', 'Ckb', 'Runx1t1', 'Celf4', 
              'Tlx3', 'Ntf3', 'Mfap2', 'Igfbp2', 'Spsb4', 'Plagl1', 'Hpcal1', 'Ier5', 'Ephb2', 'Mdk', 'Zmiz1')
allspcor <- c()
for(i in newspcor){
  cor1 <- apply(nnmouse@assays$RNA@data[newspcor, ], 1,
                function(x) cor.test(x, nnmouse@assays$RNA@data[i, ]))
  cor1 <- unlist(lapply(cor1, function(x) x$estimate))
  allspcor <- cbind(allspcor, cor1)
}
write.csv(allspcor, file = 'corr of cerebellum spatial in find mouse URL GCP.csv', quote = F)
# 

# 
nnmouse$celltype <- nnmouse$seurat_clusters
levels(nnmouse$celltype)[which(levels(nnmouse$celltype) %in% c(9,3,2,8,6,14,19))] <- 'AT'
levels(nnmouse$celltype)[which(levels(nnmouse$celltype) %in% c(7,12,13,18,15,11))] <- 'ND'
levels(nnmouse$celltype)[which(levels(nnmouse$celltype) %in% c(16,1,5,10,0,4,17))] <- 'post'
mousepro <- as.data.frame(AverageExpression(nnmouse, assays = 'RNA', group.by = 'celltype')[[1]])
humanpro <- as.data.frame(AverageExpression(progc, assays = 'RNA', group.by = 'protype')[[1]])
humansele <- humanpro %>% filter(pro>5 & afcc<2)
# MGP,NPY
pdf(file = 'hsa-mmu cross diff.pdf', height = 6, width = 6)
print(FeaturePlot(nnmouse, features = c('Sfrp1', 'Mgp', 'Npy'), 
                  cols = c('#eeeeee', '#036eb8')) & mockdim & labs(title = element_blank()))
print(FeaturePlot(progc, features = c('SFRP1', 'MGP', 'NPY'),
                  cols = c('#eeeeee', '#D40041')) & mockdim & labs(title = element_blank()))
dev.off()



# 
# 
mmupost <- subset(nnmouse, seurat_clusters %in% c(4,0,5,1,16,11))
gc()
DefaultAssay(mmupost) <- 'integrated'
mmupost <- ScaleData(mmupost)
mmupost <- RunPCA(mmupost, npcs = 50)
ElbowPlot(object = mmupost, ndims = 40, reduction = "pca") 
DefaultAssay(mmupost) <- 'integrated'
mmupost <- FindNeighbors(mmupost, reduction = "pca", dims = 1:16)
mmupost <- FindClusters(mmupost, resolution = 1)
mmupost <- RunUMAP(mmupost, reduction = "pca", dims = 1:16)
fmgc <- DimPlot(mmupost, label = T, label.size = 6)
fmgc
DefaultAssay(mmupost) <- 'RNA'
mgc13 <- FindMarkers(mmupost, ident.1 = 13, logfc.threshold = 0.7, min.diff.pct = 0.3)
mgc11 <- FindMarkers(mmupost, ident.1 = 11, logfc.threshold = 0.7, min.diff.pct = 0.3)
# 
mmupost <- subset(mmupost, seurat_clusters %in% c(0:10,12))
gc()
DefaultAssay(mmupost) <- 'integrated'
mmupost <- ScaleData(mmupost)
mmupost <- RunPCA(mmupost, npcs = 50)
ElbowPlot(object = mmupost, ndims = 40, reduction = "pca") 
DefaultAssay(mmupost) <- 'integrated'
mmupost <- FindNeighbors(mmupost, reduction = "pca", dims = 1:14)
mmupost <- FindClusters(mmupost, resolution = 1)
mmupost <- RunUMAP(mmupost, reduction = "pca", dims = 1:14)
fmgc <- DimPlot(mmupost, label = T, label.size = 6)
fmgc
DefaultAssay(mmupost) <- 'RNA'
manti <- FindMarkers(mmupost, ident.1 = 4, ident.2 = 9, 
                     logfc.threshold = 0.5)





load(file = "E:/cerebellum/bior.Rdata")
pcw11 <- subset(seurat, Stage == '11 wpc')
DefaultAssay(pcw11) <- 'RNA'
samlist <- c('SN022', 'SN023', 'SN121')
for(i in samlist){
  assign(i, subset(pcw11, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(SN022, SN023, SN121), 
                                 anchor.features = 3000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
pcw11 <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
pcw11 <- ScaleData(pcw11)
pcw11 <- RunPCA(pcw11, npcs = 50)
ElbowPlot(object = pcw11, ndims = 40, reduction = "pca") 
DefaultAssay(pcw11) <- 'integrated'
pcw11 <- FindNeighbors(pcw11, reduction = "pca", dims = 1:35)
pcw11 <- FindClusters(pcw11, resolution = 1)
pcw11 <- RunUMAP(pcw11, reduction = "pca", dims = 1:35)
fp11 <- DimPlot(pcw11, label = T, label.size = 5)
fp11
DefaultAssay(pcw11) <- 'RNA'
# 
translist <- bitr(rownames(pcw11), fromType = 'ENSEMBL', toType =c('ENSEMBL', 'SYMBOL'), 
                  OrgDb = 'org.Hs.eg.db')
# 
# 
pcw17 <- subset(seurat, Stage == '17 wpc')
DefaultAssay(pcw17) <- 'RNA'
table(pcw17$orig.ident)
samlist <- c('SN097', 'SN100', 'SN133')
for(i in samlist){
  assign(i, subset(pcw17, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(SN097, SN100, SN133), 
                                 anchor.features = 3000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
pcw17 <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
pcw17 <- ScaleData(pcw17)
pcw17 <- RunPCA(pcw17, npcs = 50)
ElbowPlot(object = pcw17, ndims = 40, reduction = "pca") 
DefaultAssay(pcw17) <- 'integrated'
pcw17 <- FindNeighbors(pcw17, reduction = "pca", dims = 1:29)
pcw17 <- FindClusters(pcw17, resolution = 1)
pcw17 <- RunUMAP(pcw17, reduction = "pca", dims = 1:29)
fp17 <- DimPlot(pcw17, label = T, label.size = 5)
fp17
DefaultAssay(pcw17) <- 'RNA'
# 
translist <- bitr(rownames(pcw11), fromType = 'ENSEMBL', toType =c('ENSEMBL', 'SYMBOL'), 
                  OrgDb = 'org.Hs.eg.db')