load(file = 'human pkc.Rdata')
DimPlot(purkinje, label = T)
table(interneuron$celltype)
interneuron$celltype2 <- interneuron$celltype
#
purkinje$celltype <- purkinje$seurat_clusters
purkinje$celltype2 <- purkinje$seurat_clusters
levels(purkinje$celltype)[which(levels(purkinje$celltype) %in% c(8))] <- 'EBF2+'
levels(purkinje$celltype)[which(levels(purkinje$celltype) %in% c(9,0,2,3,4,5,6))] <- 'PBX3+'
levels(purkinje$celltype)[which(levels(purkinje$celltype) %in% c(1,7))] <- 'LSAMP+'
table(purkinje$celltype)
levels(purkinje$celltype2)[which(levels(purkinje$celltype2) %in% c(8))] <- 'EBF2+'
levels(purkinje$celltype2)[which(levels(purkinje$celltype2) %in% c(9,0,2,3,4,5,6,1,7))] <- 'EBF2-'
table(purkinje$celltype2)
# 
inpkc <- merge(purkinje, interneuron)
table(inpkc$celltype)
table(inpkc$celltype2)
save(inpkc, file = 'human IN+PKC for ATAC.Rdata')


# 
mouE13 <- readRDS(file = "E:/cerebellum/cerebellum_mouE13B.Rdata")
mouE14 <- readRDS(file = "E:/cerebellum/cerebellum_mouE14B.Rdata")
mouE15 <- readRDS(file = "E:/cerebellum/cerebellum_mouE15B.Rdata")
mouE17 <- readRDS(file = "E:/cerebellum/cerebellum_mouE17A.Rdata")
mouE18 <- readRDS(file = "E:/cerebellum/cerebellum_mouE18.Rdata")
mouselist <- c('mouE13', 'mouE14', 'mouE15', 'mouE17', 'mouE18')
for(i in 1:5){
  assign(mouselist[i], NormalizeData(get(mouselist[i])))
  assign(mouselist[i], FindVariableFeatures(get(mouselist[i]), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(mouE13, mouE14, mouE15, mouE17, mouE18), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
mouseceb <- IntegrateData(anchorset = anchor, dims = 1:20)
gc()
mouseceb <- ScaleData(mouseceb)
DefaultAssay(mouseceb) <- 'integrated'
mouseceb <- RunPCA(mouseceb, npcs = 50)
ElbowPlot(object = mouseceb, ndims = 40, reduction = "pca") 
DefaultAssay(mouseceb) <- 'integrated'
mouseceb <- FindNeighbors(mouseceb, reduction = "pca", dims = 1:14)
mouseceb <- FindClusters(mouseceb, resolution = 1)
mouseceb <- RunUMAP(mouseceb, reduction = "pca", dims = 1:14)
fmceb <- DimPlot(mouseceb, label = T, label.size = 6)
fmceb
# 
mouseceb <- subset(mouseceb, seurat_clusters %in% setdiff(0:26, c(22,23,26)))
rm(mouE13, mouE14, mouE15, mouE17, mouE18)
gc()
DefaultAssay(mouseceb) <- 'RNA'
for(i in 1:5){
  obj <- subset(mouseceb, orig.ident == rownames(table(mouseceb$orig.ident))[i])
  obj <- NormalizeData(obj, normalization.method = 'LogNormalize')
  assign(mouselist[i], FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))
}
rm(obj)
anchor <- FindIntegrationAnchors(object.list = list(mouE13, mouE14, mouE15, mouE17, mouE18), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
mouseceb <- IntegrateData(anchorset = anchor, dims = 1:20)
gc()
mouseceb <- ScaleData(mouseceb)
DefaultAssay(mouseceb) <- 'integrated'
mouseceb <- RunPCA(mouseceb, npcs = 50)
ElbowPlot(object = mouseceb, ndims = 40, reduction = "pca") 
DefaultAssay(mouseceb) <- 'integrated'
mouseceb <- FindNeighbors(mouseceb, reduction = "pca", dims = 1:15)
mouseceb <- FindClusters(mouseceb, resolution = 1)
mouseceb <- RunUMAP(mouseceb, reduction = "pca", dims = 1:15)
fmceb <- DimPlot(mouseceb, label = T, label.size = 6)
fmceb
DefaultAssay(mouseceb) <- 'RNA'
m25 <- FindMarkers(mouseceb, ident.1 = 25, min.diff.pct = 0.3, logfc.threshold = 0.8)
mouseceb <- subset(mouseceb, seurat_clusters %in% 0:24)
fmceb <- DimPlot(mouseceb, label = T, label.size = 6)
fmceb
save(mouseceb, file = 'mousecebE13-18.Rdata')
rm(mouE13, mouE14, mouE15, mouE17, mouE18)

load(file = 'mousecebE13-18.Rdata')


#
load(file = 'mousecebE13-18.Rdata')
# 
otholist <- read.csv(file = "E:/cerebellum/hsa-mmu-otholist v108.txt", 
                       header = T, encoding = 'utf-8', sep = '\t')
otholist <- otholist %>% filter(Gene.name != '' & Mouse.gene.name != '')
# 
dpkhsa <- otholist$Gene.name[duplicated(otholist$Gene.name)]
dpkmmu <- otholist$Mouse.gene.name[duplicated(otholist$Mouse.gene.name)]
dpkgene <- otholist %>% filter(Gene.name %in% dpkhsa | Mouse.gene.name %in% dpkmmu)
one2one <- anti_join(otholist, dpkgene, "Gene.name")
one2one <- rbind(one2one, dpkgene[which(dpkgene$Mouse.gene.name == 'Gapdh'), ])
one2one <- one2one %>% filter(Gene.name %in% rownames(cebclean) | 
                                Mouse.gene.name %in% rownames(mouseceb))
# 

mouseceb <- mouseceb[one2one$Mouse.gene.name, ];gc()
cebclean <- cebclean[one2one$Gene.name, ];gc()

mmuforint <- subset(mouseceb, seurat_clusters %in% c(0,11,14,24,19,9,21,12,5,6,7,18))
hsaforint <- ppicell
gc()
DefaultAssay(mmuforint) <- 'RNA'
DefaultAssay(hsaforint) <- 'RNA'
mmuforint <- mmuforint[one2one$Mouse.gene.name, ];gc()
o2oformmu <- one2one %>% filter(Mouse.gene.name %in% rownames(mmuforint))
rownames(o2oformmu) <- o2oformmu$Mouse.gene.name
o2oformmu <- o2oformmu[rownames(mmuforint), ]
# 
mmuforint$orig.ident[mmuforint$orig.ident == 'cerebellum_mouE13B'] <- 'E13'
mmuforint$orig.ident[mmuforint$orig.ident == 'cerebellum_mouE14B'] <- 'E14'
mmuforint$orig.ident[mmuforint$orig.ident == 'cerebellum_mouE15B'] <- 'E15'
mmuforint$orig.ident[mmuforint$orig.ident == 'cerebellum_mouE17A'] <- 'E17'
mmuforint$orig.ident[mmuforint$orig.ident == 'cerebellum_mouE18'] <- 'E18'

for(i in names(table(mmuforint$orig.ident))){
  obj <- subset(mmuforint, orig.ident == i)
  obj <- GetAssayData(obj, slot = 'count')
  rownames(obj) <- o2oformmu$Gene.name
  obj <- CreateSeuratObject(obj, project = i, min.cells = 0, min.features = 0)
  obj <- RenameCells(obj, add.cell.id = i)
  obj <- NormalizeData(obj, normalization.method = 'LogNormalize')
  assign(i, FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))
}
# 

hsaforint <- hsaforint[one2one$Gene.name, ];gc()
o2oforhsa <- one2one %>% filter(Gene.name %in% rownames(hsaforint))
rownames(o2oforhsa) <- o2oforhsa$Gene.name
o2oforhsa <- o2oforhsa[rownames(hsaforint), ]
# 
for(i in names(table(hsaforint$orig.ident))){
  obj <- subset(hsaforint, orig.ident == i)
  obj <- GetAssayData(obj, slot = 'count')
  rownames(obj) <- o2oforhsa$Gene.name
  obj <- CreateSeuratObject(obj, project = i, min.cells = 0, min.features = 0)
  obj <- RenameCells(obj, add.cell.id = i)
  obj <- NormalizeData(obj, normalization.method = 'LogNormalize')
  assign(i, FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))
}
# 
rm(obj)
#
anchor <- FindIntegrationAnchors(object.list = list(E13, E14, E15, E17, E18, E075, E078, E083, E085, E101, E104), 
                                 anchor.features = 2000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:20)
inteppi <- IntegrateData(anchorset = anchor, dims = 1:20)
gc()
inteppi <- ScaleData(inteppi)
DefaultAssay(inteppi) <- 'integrated'
inteppi <- RunPCA(inteppi, npcs = 50)
ElbowPlot(object = inteppi, ndims = 40, reduction = "pca") 
DefaultAssay(inteppi) <- 'integrated'
inteppi <- FindNeighbors(inteppi, reduction = "pca", dims = 1:9)
inteppi <- FindClusters(inteppi, resolution = 1)
inteppi <- RunUMAP(inteppi, reduction = "pca", dims = 1:9)
fhmp <- DimPlot(inteppi, label = T, label.size = 6)
fhmp
DefaultAssay(inteppi) <- 'RNA'
save(inteppi, file = 'int.hsa-mmu pro+pkc+in.Rdata')
# 
# 
load(file = 'int.hsa-mmu pro+pkc+in.Rdata')
DimPlot(inteppi, label = T)
FeaturePlot(inteppi, 'RGS16')
#
inte16 <- FindMarkers(inteppi, ident.1 = 16, logfc.threshold = 0.7, min.diff.pct = 0.3)
inteppi$species <- inteppi$orig.ident
inteppi$species[which(inteppi$species %in% c('E075','E078','E083','E085','E101','E104'))] <- 'hsa'
inteppi$species[which(inteppi$species %in% c('E13','E14','E15','E17','E18'))] <- 'mmu'
FeaturePlot(inteppi, 'WNT1', split.by = 'species')
# 
inte13 <- FindMarkers(inteppi, ident.1 = 13, ident.2 = 0, 
                      logfc.threshold = 0.5, min.diff.pct = 0.2)
# 
FeaturePlot(subset(inteppi, seurat_clusters %in% c(0,13)), 'EBF2', split.by = 'species')

# 
# 
DefaultAssay(inteppi) <- 'RNA'
intepkc <- subset(inteppi, seurat_clusters %in% c(2,4,0,13))
intepkc$oldcluster <- intepkc$seurat_clusters
for(i in names(table(intepkc$orig.ident))){
  assign(i, subset(intepkc, orig.ident == i))
  assign(i, NormalizeData(get(i), normalization.method = 'LogNormalize'))
  assign(i, FindVariableFeatures(get(i), selection.method = "vst", nfeatures = 2000))
}
anchor <- FindIntegrationAnchors(object.list = list(E13, E14, E15, E17, E18, E075, E078, E083, E101, E104), 
                                 anchor.features = 3000, 
                                 scale = T, 
                                 reduction = 'cca', 
                                 dims = 1:30, k.filter = 100)
gc()
intepkc <- IntegrateData(anchorset = anchor, dims = 1:30, k.weight = 100)
gc()
intepkc <- ScaleData(intepkc)
intepkc <- RunPCA(intepkc, npcs = 50)
ElbowPlot(object = intepkc, ndims = 40, reduction = "pca") 
DefaultAssay(intepkc) <- 'integrated'
intepkc <- FindNeighbors(intepkc, reduction = "pca", dims = 1:14)
intepkc <- FindClusters(intepkc, resolution = 1)
intepkc <- RunUMAP(intepkc, reduction = "pca", dims = 1:14)
fipk <- DimPlot(intepkc, label = T, label.size = 5)
fipk
DefaultAssay(intepkc) <- 'RNA'
ipk12 <- FindMarkers(intepkc, ident.1 = 12, ident.2 = c(1,8,14,15),
                     min.diff.pct = 0.2, logfc.threshold = 0.5)
ipk69 <- FindMarkers(intepkc, ident.1 = c(6,9),
                     min.diff.pct = 0.2, logfc.threshold = 0.5)
ipk24 <- FindMarkers(intepkc, ident.1 = c(2,4),
                     min.diff.pct = 0.2, logfc.threshold = 0.5)
DimPlot(intepkc, group.by = 'species')