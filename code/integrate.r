library(Seurat)
# read in data and rename timepoint
comceb <- readRDS(file = "E:/cerebellum/ceb.combined.rds")
timepoint_map <- c('cerebellum_12.5' = 'E075', 'cerebellum_13.1' = 'E078', 
                   'cerebellum_13.6' = 'E083', 'cerebellum_14.1' = 'E085', 
                   'cerebellum_16.3' = 'E101', 'cerebellum_16.6' = 'E104')
comceb$timepoint <- timepoint_map[comceb$orig.ident]

# loop through each timepoint and perform analysis
for (tp in unique(comceb$timepoint)) {
  sub <- subset(comceb, timepoint == tp)
  obj <- CreateSeuratObject(counts = GetAssayData(sub, slot = 'count'), 
                            project = tp)
  obj <- NormalizeData(obj, normalization.method = 'LogNormalize')
  features <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = features)
  obj <- RunPCA(obj, npcs = 50)
  obj <- FindNeighbors(obj, dims = 1:15)
  obj <- FindClusters(obj, resolution = 1)
  obj <- RunUMAP(obj, dims = 1:15)
  assign(paste0('sub', tp), obj)
}

# find marker genes using subE085
subE085 <- ScaleData(subE085)
subE085 <- RunPCA(subE085, npcs = 50)
ElbowPlot(object = subE085, ndims = 40, reduction = "pca") 
subE085 <- FindNeighbors(subE085, reduction = "pca", dims = 1:15)
subE085 <- FindClusters(subE085, resolution = 1)
subE085 <- RunUMAP(subE085, reduction = "pca", dims = 1:15, min.dist = 0.2)
fe <- DimPlot(subE085, label = T, label.size = 6)
mprogenitor <- FindMarkers(subE085, ident.1 = c(6, 4, 17), logfc.threshold = 0.8, 
                           min.diff.pct = 0.3, only.pos = T)
mubc <- FindMarkers(subE085, ident.1 = c(1, 11), ident.2 = c(0, 2, 13, 6, 14, 10), 
                    logfc.threshold = 0.8, min.diff.pct = 0.3, only.pos = T)
mgc1 <-FindMarkers(subE085, ident.1 = c(13, 0, 2), ident.2 = c(6, 14, 10, 1, 11), 
                   logfc.threshold = 0.8, min.diff.pct = 0.3, only.pos = T)
mgc2 <-FindMarkers(subE085, ident.1 = c(6, 14, 10), ident.2 = c(1, 11, 2, 0, 13), 
                   logfc.threshold = 0.8, min.diff.pct = 0.3, only.pos = T)

rm(comceb, origdata, obj)
gc()
# 用subE85来挑markergene
subE075 <- ScaleData(subE075)
subE075 <- RunPCA(subE075, npcs = 50)
ElbowPlot(object = subE075, ndims = 40, reduction = "pca") 
subE075 <- FindNeighbors(subE075, reduction = "pca", dims = 1:15)
subE075 <- FindClusters(subE075, resolution = 1)
subE075 <- RunUMAP(subE075, reduction = "pca", dims = 1:15)
fe75 <- DimPlot(subE075, label = T, label.size = 6)
fe75
mprogenitor <- FindMarkers(subE075, ident.1 = c(6,4,17), 
                           logfc.threshold = 0.8, min.diff.pct = 0.3, only.pos = T)
mubc <- FindMarkers(subE075, ident.1 = c(1,11), ident.2 = c(0,2,13,6,14,10), 
                    logfc.threshold = 0.8, min.diff.pct = 0.3, only.pos = T)
mgc1 <-FindMarkers(subE075, ident.1 = c(13,0,2), ident.2 = c(6,14,10,1,11), 
                   logfc.threshold = 0.8, min.diff.pct = 0.3, only.pos = T)
mgc2 <-FindMarkers(subE075, ident.1 = c(6,14,10), ident.2 = c(1,11,2,0,13), 
                   logfc.threshold = 0.8, min.diff.pct = 0.3, only.pos = T)
minter <- FindMarkers(subE075, ident.1 = 12, 
                      logfc.threshold = 0.8, min.diff.pct = 0.3, only.pos = T)
mpurkin <- FindMarkers(subE075, ident.1 = c(5,7,9), 
                       logfc.threshold = 1.3, min.diff.pct = 0.5, only.pos = T)
# m18 <- FindMarkers(subE075, ident.1 = 18, logfc.threshold = 0.8)
anchorgene <- unique(c(rownames(mprogenitor), rownames(mubc), rownames(mgc1), rownames(mgc2),
                       rownames(minter), rownames(mpurkin)))
# save(anchorgene, file = 'selected anchorgene.Rdata')
load(file = 'selected anchorgene.Rdata')
anchor <- FindIntegrationAnchors(object.list = list(subE075, subE078, subE083, subE085, subE101, subE104), 
                                 anchor.features = anchorgene, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
cebclean <- IntegrateData(anchorset = anchor, dims = 1:20)
gc()
cebclean <- ScaleData(cebclean)
DefaultAssay(cebclean) <- 'integrated'
cebclean <- RunPCA(cebclean, features = anchorgene, npcs = 50)
ElbowPlot(object = cebclean, ndims = 40, reduction = "pca") 
DefaultAssay(cebclean) <- 'integrated'
cebclean <- FindNeighbors(cebclean, reduction = "pca", dims = 1:7)
cebclean <- FindClusters(cebclean, resolution = 1)
cebclean <- RunUMAP(cebclean, reduction = "pca", dims = 1:7, min.dist = 0.2)
fceb <- DimPlot(cebclean, label = T, label.size = 6)
fceb
#
save(cebclean, file = 'cebclean.Rdata')
rm(subE075, subE078, subE083, subE085, subE101, subE104)
