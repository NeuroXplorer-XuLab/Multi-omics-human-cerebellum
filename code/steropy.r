setwd("E:/cerebellum/spatial")
load(file = 'E:/cerebellum/spatial/GW13_raw_out.Rdata')
spgw13 <- seurat_spatialObj;rm(seurat_spatialObj)
# 
VlnPlot(spgw13, 'nFeature_Spatial')
FeatureScatter(spgw13, 'nCount_Spatial', 'nFeature_Spatial')
spgw13 <- subset(spgw13, nFeature_Spatial > 500)
spgw13 <- NormalizeData(spgw13, normalization.method = 'LogNormalize')
spgw13 <- FindVariableFeatures(spgw13, selection.method = "vst", nfeatures = 2000)
spgw13 <- ScaleData(spgw13)
spgw13 <- RunPCA(spgw13, npcs = 50)
ElbowPlot(spgw13, ndims = 40, reduction = "pca")
spgw13 <- FindNeighbors(spgw13, reduction = "pca", dims = 1:8)
spgw13 <- FindClusters(spgw13, resolution = 1)
spgw13 <- RunUMAP(spgw13, reduction = "pca", dims = 1:8)
fsp13 <- DimPlot(spgw13, label = T, label.size = 6)
fsp13
# 
feamean <- mean(spgw13$nFeature_Spatial) # 2481.4

sFeatureplot(spgw13, feature = 'LSAMP')
# 
# 
sDimplot(spgw13, color = c('grey','grey','brown','green','grey','orange','red','navy','grey','grey','grey'), group = 'seurat_clusters')
#
sDimplot(spgw13, color = c('grey','grey','grey','grey','grey','grey','red','grey','grey','grey','grey'), group = 'seurat_clusters')

derem3 <- FindMarkers(spgw13, ident.1 = 3, min.diff.pct = 0.3, 
                      logfc.threshold = 0.5, only.pos = T)
derem2 <- FindMarkers(spgw13, ident.1 = 2, 
                      logfc.threshold = 0.5, only.pos = T)
derem7 <- FindMarkers(spgw13, ident.1 = 7, min.diff.pct = 0.3, 
                      logfc.threshold = 0.5, only.pos = T)
derem3anti7 <- FindMarkers(spgw13, ident.1 = 3, ident.2 = 7, 
                           min.diff.pct = 0.1, logfc.threshold = 0.5)
derem2anti3 <- FindMarkers(spgw13, ident.1 = 2, ident.2 = 3, 
                           min.diff.pct = 0.3, logfc.threshold = 0.5)
derem5 <- FindMarkers(spgw13, ident.1 = 5, min.diff.pct = 0.2, logfc.threshold = 0.5)
der1anti0 <- FindMarkers(spgw13, ident.1 = 1, ident.2 = 0)
derem2anti37 <- FindMarkers(spgw13, ident.1 = 2, ident.2 = c(3,7), 
                            min.diff.pct = 0.2, logfc.threshold = 0.5)

# 
# 2/5 
derem52 <- FindMarkers(spgw13, ident.1 = c(5,2), min.diff.pct = 0.3, 
                       logfc.threshold = 0.5)
derem2 <- FindMarkers(spgw13, ident.1 = 2, min.diff.pct = 0.3, 
                      logfc.threshold = 0.5, only.pos = T)
derem5 <- FindMarkers(spgw13, ident.1 = 5, only.pos = T)
derem2anti5 <- FindMarkers(spgw13, ident.1 = 2, ident.2 = 5, 
                           min.diff.pct = 0.3, logfc.threshold = 0.5)
# 
pdf('sp GW13 dimplot.pdf', width = 6.9, height = 8.9)
print(sDimplot(spgw13, color = '', group = 'seurat_clusters') & NoLegend() & 
        labs(title = element_blank()))
dev.off()
# 
gcmmm <- c('PRPH', 'ATOH1', 'PAX2', 'CA8', 'EOMES', 'HES1')
gcmmm <- c('WLS', 'NHLH1', 'SOX2')
for(i in gcmmm){
  png(paste0('sp GW13 feature ', i,  '.png'), width = 345, height = 445)
  print(sFeatureplot(spgw13, feature = i, bgcolor = '#ffffff', color = c('#dedede', 'navy'))
        &labs(title = element_blank())&NoLegend())
  dev.off()
}
#
# 
FeaturePlot(spgw13, features = c('BARHL1', 'TLX3'), blend = T, blend.threshold = 0.1, 
            cols = c('#eeeeee', '#036eb8', '#D40041'))
mcolor <- rbind(c('#e7e4e8', '#becfe1', '#85b3d8', '#4d97d0', '#036eb8'), 
                c('#e5bbcb', '#c6afce', '#999fd0', '#6d90d5', '#4683db'), 
                c('#e481a2', '#d081b1', '#af7dc0', '#907ad1', '#757ae3'), 
                c('#e4487b', '#dd5596', '#c85db1', '#b566ce', '#a572ec'), 
                c('#d40041', '#ea2d7e', '#e241a6', '#da56cf', '#d76ef9'))
mcolor <- as.data.frame(mcolor)
colnames(mcolor) <- c(1,2,3,4,5)
rownames(mcolor) <- c(1,2,3,4,5)
antipair <- FetchData(spgw13, c('PLCB4', 'CADM1'), slot = 'data')
# antipair[,2] <- 0
# 'BARHL1', 'TLX3'
cutbarhl <- cut(antipair[,1], 
                breaks = c(-Inf, max(antipair[,1])*0.2, max(antipair[,1])*0.4, 
                           max(antipair[,1])*0.6, max(antipair[,1])*0.8, max(antipair[,1])), 
                labels = c(1,2,3,4,5))
cuttlx3 <- cut(antipair[,2], 
                breaks = c(-Inf, max(antipair[,2])*0.2, max(antipair[,2])*0.4, 
                           max(antipair[,2])*0.6, max(antipair[,2])*0.8, max(antipair[,2])), 
                labels = c(1,2,3,4,5))
# cuttlx3 <- rep(1, length(cuttlx3))
colorcol <- c()
for(i in 1:nrow(antipair)){
  colorcol <- c(colorcol, mcolor[cutbarhl[i], cuttlx3[i]])
}
antipair <- cbind(antipair, colorcol, spgw13[['x']], spgw13[['y']])
ff <- ggplot(antipair, aes(x = antipair$x, y = antipair$y, color = antipair[,3])) + 
  geom_point(shape=19) + scale_color_manual(values = unique(antipair$colorcol)[order(unique(antipair$colorcol), decreasing = F)]) + 
  xlim(min(antipair$x), max(antipair$x)) + ylim(min(antipair$y), max(antipair$y)) + 
  labs(x = "", y = "", title = element_blank()) +
  theme_bw() + theme(panel.background = element_rect(fill = '#ffffff')) +
  theme(panel.grid =element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank()) + NoLegend()
ff
pdf('sp GW13 antipair PLCB4+CADM1.pdf', width = 6.9, height = 8.9)
print(ff)
dev.off()
# 
pdf('sp GW13 antipair ATOH1+0.pdf', width = 6.9, height = 8.9)
print(ff)
dev.off()

nnnspgw13 <- subset(spgw13, seurat_clusters %in% setdiff(0:8, c(4)))
for(i in 1:2){
  count <- seq((i-1)*100+1, i*100)
  pdf(file = paste0('hsa-ceb spatial ', i, ' feature', '.pdf'), width = 6.3, height = 6.8)
  for(j in count){
    print(sFeatureplot(nnnspgw13, feature = rownames(nnnspgw13)[j]))
  }
  dev.off()
}
# 
max(spgw13$y); min(spgw13$y) ; max(spgw13$x); min(spgw13$x)

locmatrix <- as.data.frame(cbind(spgw13$x, spgw13$y))
colnames(locmatrix) <- c('x', 'y')
sdiffventral <- locmatrix %>% filter((x<=13500 & y>=13000) | (x>=13500 & y>=12000)) %>% rownames()
sdiffdor1 <- locmatrix %>% filter((x<=13500 & y<13000)) %>% rownames()
sdiffdor2 <- locmatrix %>% filter((x>13500 & y<12000)) %>% rownames()
sdiffdor1 <- setdiff(sdiffdor1, sdiffventral)
# 
spgw13$threesp <- colnames(spgw13)
spgw13$threesp[which(spgw13$threesp %in% sdiffventral)] <- 'ven'
spgw13$threesp[which(spgw13$threesp %in% sdiffdor1)] <- 'dor1'
spgw13$threesp[which(spgw13$threesp %in% sdiffdor2)] <- 'dor2'
# 
sDimplot(spgw13, color = '', group = 'threesp')
# 
eglsp <- subset(spgw13, seurat_clusters %in% c(3))
mspvvv <- FindMarkers(eglsp, ident.1 = 'ven', group.by = 'threesp')
mspd1d1 <- FindMarkers(eglsp, ident.1 = 'dor1', group.by = 'threesp')
mspd2d2 <- FindMarkers(eglsp, ident.1 = 'dor2', group.by = 'threesp')


# 
# 
eglsp <- subset(spgw13, seurat_clusters %in% c(3,7))
pdf('sp GW13 PRR-EBF2-HEY1.pdf', width = 6.2, height = 6.4)
print(sFeatureplot(eglsp, feature = 'EBF2', color = c('grey', '#D40041', '#D40041'))&mockdim)
print(sFeatureplot(eglsp, feature = 'PRR35', color = c('grey', '#036eb8', '#036eb8'))&mockdim)
print(sFeatureplot(eglsp, feature = 'HEY1', color = c('grey', '#009f35', '#009f35'))&mockdim)
dev.off()



# 
coresp <- subset(spgw13, seurat_clusters %in% c(0,5,1))
sDimplot(coresp, color = c('orange','skyblue','firebrick'), group = 'seurat_clusters')
# max(coresp$y)-min(coresp$y)
pdf('sp GW13 PLCB4 CADM1.pdf', width = 4.8, height = 5.7)
print(sFeatureplot(coresp, feature = 'PLCB4', color = c('grey', 'grey', '#D40041', '#D40041'))&mockdim)
print(sFeatureplot(coresp, feature = 'CADM1', color = c('grey', 'grey', '#036eb8', '#036eb8'))&mockdim)
dev.off()


# 
sforprph <- subset(spgw13, x<=13500 & y>=13000)
sDimplot(sforprph, color = '', group = 'seurat_clusters')
sforprph <- FindVariableFeatures(sforprph, selection.method = "vst", nfeatures = 1500)
sforprph <- ScaleData(sforprph)
sforprph <- RunPCA(sforprph, npcs = 50)
ElbowPlot(sforprph, ndims = 40, reduction = "pca")
sforprph <- FindNeighbors(sforprph, reduction = "pca", dims = 1:20)
sforprph <- FindClusters(sforprph, resolution = 1)
sforprph <- RunUMAP(sforprph, reduction = "pca", dims = 1:20)
fsprp <- DimPlot(sforprph, label = T, label.size = 6)
fsprp
sDimplot(sforprph, color = '', group = 'seurat_clusters')
# 
pdf(file = 'test.pdf', width = 3, height = 3)
for(i in c('PCP4','PCDH17','SKOR2','WNT7B','SLC1A6','FGF3','SLC32A1','KITLG',
           'WFDC2','PNCK','CA8','PRPH','SFRP1','ZFP36L1')){
  print(sFeatureplot(sforprph, feature = i))
}
dev.off()
# 
pickedx <- colnames(sforprph)[sforprph$x >= 12100 & sforprph$x <= 12400]
pickedy <- colnames(sforprph)[sforprph$y == 14400 | sforprph$y == 14500]
picked <- intersect(pickedx, pickedy)
# 
sforprph$pickident <- colnames(sforprph)
sforprph$pickident[which(!(sforprph$pickident %in% picked))] <- 1
sforprph$pickident[which(sforprph$pickident %in% picked)] <- 2
mtest <- FindMarkers(sforprph, ident.1 = 2, group.by = 'pickident', 
                     logfc.threshold = 0.5, only.pos = T, min.pct = 0.3)
# 
sFeatureplot(sforprph, feature = 'VANGL1')


s23pg13 <- subset(spgw13, seurat_clusters %in% c(2,3))
s23pg13 <- FindVariableFeatures(s23pg13, selection.method = "vst", nfeatures = 1500)
s23pg13 <- ScaleData(s23pg13)
s23pg13 <- RunPCA(s23pg13, npcs = 50)
ElbowPlot(s23pg13, ndims = 40, reduction = "pca")
s23pg13 <- FindNeighbors(s23pg13, reduction = "pca", dims = 1:10)
s23pg13 <- FindClusters(s23pg13, resolution = 1)
s23pg13 <- RunUMAP(s23pg13, reduction = "pca", dims = 1:10)
fs222 <- DimPlot(s23pg13, label = T, label.size = 6)
fs222
msp22de1 <- FindMarkers(s23pg13, ident.1 = 1, ident.2 = 2)


# 
pdf('check sp GC endpoint.pdf')
for(i in c('ADAMTS6', 'TLL1', 'TMEM132B', 'EPHA3', 'NR2F2', 'AC015574.1', 
           'SEMA3C', 'ALCAM', 'C11orf96', 'MEIS2', 'EBF2', 
           'EYA2', 'RAB26', 'ITGBL1', 
           'GALNTL6','STC1', 'BRINP2', 'NTM', 'RPRM', 'CSMD3', 'DOK6', 
           'C1QTNF3', 'NRXN3', 'PPP1R1B', 'BHLHE41', 'GMDS')){
  print(sFeatureplot(spgw13, feature = i))
}
dev.off()
# 
# 
pdf(file = 'hsa-ceb UBC spatial four gene.pdf', width = 6.9, height = 8.9)
for(i in c('COL2A1', 'CRYAB', 'CALCB', 'EOMES')){
  fprint <- sFeatureplot(spgw13, feature = i, bgcolor = '#ffffff', color = c('#dddddd', 'navy'))& 
    NoLegend()
  print(fprint)
}
dev.off()