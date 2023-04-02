load(file = "E:/cerebellum/cere_ns.Rdata")
for(i in names(table(cere_ns$age))){
  assign(i, subset(cere_ns, age == i))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 1000))
}
anchor <- FindIntegrationAnchors(object.list = list(`10 PCW`, `11 PCW`, `12 PCW`), 
                                 anchor.features = 1000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30)
gc()
cere_ns <- IntegrateData(anchorset = anchor, dims = 1:30)
gc()
cere_ns <- ScaleData(cere_ns)
cere_ns <- RunPCA(cere_ns, npcs = 50)
ElbowPlot(cere_ns, ndims = 40, reduction = "pca")
DefaultAssay(cere_ns) <- 'integrated'
cere_ns <- FindNeighbors(cere_ns, reduction = "pca", dims = 1:23)
cere_ns <- FindClusters(cere_ns, resolution = 1)
cere_ns <- RunUMAP(cere_ns, reduction = "pca", dims = 1:23)
fnns <- DimPlot(cere_ns, label = T, label.size = 6)
fnns
DefaultAssay(cere_ns) <- 'RNA'
mnns25 <- FindMarkers(cere_ns, ident.1 = 25, min.diff.pct = 0.3, logfc.threshold = 0.5)
VlnPlot(cere_ns, 'CA8', pt.size = 0)
# 0,1,2,4,6,7,8,9,10,13,14,16,18,21,22,24,25
nnspkc <- subset(cere_ns, seurat_clusters %in% c(0,1,2,4,6,7,8,9,10,13,14,16,18,21,22,24,25))
DefaultAssay(nnspkc) <- 'integrated'
nnspkc <- ScaleData(nnspkc)
nnspkc <- RunPCA(nnspkc, npcs = 50)
ElbowPlot(nnspkc, ndims = 40, reduction = "pca")
DefaultAssay(nnspkc) <- 'integrated'
nnspkc <- FindNeighbors(nnspkc, reduction = "pca", dims = 1:24)
nnspkc <- FindClusters(nnspkc, resolution = 1)
nnspkc <- RunUMAP(nnspkc, reduction = "pca", dims = 1:24)
fnpk <- DimPlot(nnspkc, label = T, label.size = 6)
fnpk
DefaultAssay(nnspkc) <- 'RNA'
# 
# 
nnsexi <- subset(cere_ns, seurat_clusters %in% c(3,5,11,12,15,17,19,20,23))
DefaultAssay(nnsexi) <- 'integrated'
nnsexi <- ScaleData(nnsexi)
nnsexi <- RunPCA(nnsexi, npcs = 50)
ElbowPlot(nnsexi, ndims = 40, reduction = "pca")
DefaultAssay(nnsexi) <- 'integrated'
nnsexi <- FindNeighbors(nnsexi, reduction = "pca", dims = 1:18)
nnsexi <- FindClusters(nnsexi, resolution = 1)
nnsexi <- RunUMAP(nnsexi, reduction = "pca", dims = 1:18)
fnex <- DimPlot(nnsexi, label = T, label.size = 6)
fnex
DefaultAssay(nnsexi) <- 'RNA'
# 
mnexi <- FindMarkers(nnsexi, ident.1 = c(8,14,17,0,9,4,2), 
                     min.diff.pct = 0.3, logfc.threshold = 0.5)
mnexi2 <- FindMarkers(nnsexi, ident.1 = c(13,5,7,1,6), 
                     min.diff.pct = 0.3, logfc.threshold = 0.5)
mnexi3 <- FindMarkers(nnsexi, ident.1 = c(6), 
                      min.diff.pct = 0.3, logfc.threshold = 0.5)
VlnPlot(cere_ns, 'NEUROD1', pt.size = 0)