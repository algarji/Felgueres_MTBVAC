library(dplyr)
library(tidyverse)
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggthemes)
library(EnhancedVolcano)
library(sceasy)


# MERGE iBCG AND MTBVAC DATASETS

run1 <- readRDS("D:/blood/10x/ibcg.rds")	## Demultiplexed preprocessed run containing PBMCs from 3 healthy donors treated with iBCG
run2 <- readRDS("D:/blood/10x/mtbvac.rds")	## Demultiplexed preprocessed run containing PBMCs from 3 healthy donors treated with MTBVAC

pbmc.list <- list("run1" = run1, "run2" = run2)
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 2000)
immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)
pbmc.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(pbmc.combined) <- "integrated"
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, verbose = FALSE)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:30)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:30)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)

# FIGURE 1A

color <- c('grey', '#E26464', 'grey', '#967DB8', '#71A65D', 'grey', '#ADD8E4', 'grey', 'grey', 'grey', 'grey')
DimPlot(pbmc.combined, cols=color)


# IN DEPTH ANALYSIS OF MTBVAC DATASET 

mtbvac <- run2
mtbvac <- NormalizeData(mtbvac)
mtbvac <- FindVariableFeatures(mtbvac, selection.method = "mean.var.plot")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
mtbvac <- CellCycleScoring(mtbvac, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mtbvac <- ScaleData(mtbvac, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mtbvac), block.size = 100) 
mtbvac <- RunPCA(mtbvac, npcs = 20, verbose = FALSE)
mtbvac <- FindNeighbors(mtbvac, dims = 1:20)
mtbvac <- FindClusters(mtbvac, resolution = 0.5)
mtbvac <- RunUMAP(mtbvac, dims = 1:20)

# FIGURE 2A

DimPlot(mtbvac, cols=c(#43B072, #93881C, #254093, #37B9C6, #58B038, #EEA9B8, #F07E46, #DC173E, #D2A1CA, #9A9999, #B3BEB4, #EB5B98), label=TRUE)

# FIGURE 2B

DotPlot(mtbvac, features = c('TRDV2', 'TRGV9', 'CCR7', 'TCF7', 'IL7R', 'KLF2', 'CD4', 'TNFRSF4', 'CD8B', 'NELL2', 'CDC20', 'MKI67', 'TYROBP', 'NCAM1', 'ASPM', "HIST1H1B", "CD79B", "MS4A1", "FGFBP2", "FCRL6", "TRAV1-2", "KLRB1", "FCER1A", "LYZ")) + xlab('Gene') +  ylab('Cluster') + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + scale_colour_viridis(option="magma") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

# FIGURE 3A

cd4 <- subset(mtbvac, idents=c('2', '3', '4'))
VlnPlot(cd4, features=c('HOPX', 'CXCR3', 'CXCR6', 'IFNG', 'IL2RA', 'HLA-DRA', 'GNLY', 'SYTL2', 'CST7', 'CCL5', 'SELL'), cols=c(#93881C, #254093, #37B9C6))

# FIGURE 3B

cd4_subset <- subset(mtbvac, idents=c('4'))
cd4_subset <- DietSeurat(cd4_subset)
cd4_subset <- ScaleData(cd4_subset, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = rownames(cd4_subset))
cd4_subset <- RunPCA(cd4_subset, verbose = FALSE)
cd4_subset <- FindNeighbors(cd4_subset, dims = 1:15)
cd4_subset <- FindClusters(cd4_subset, resolution = 0.5)
cd4_subset <- RunUMAP(cd4_subset, dims = 1:15)
DimPlot(cd4_subset, label=TRUE)

# FIGURE 3C

cd4.markers <- FindAllMarkers(cd4_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cd4.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10_cd4
DoHeatmap(cd4_subset, features = top10_cd4$gene) + NoLegend() 

# FIGURE 3D

FeaturePlot(cd4_subset, features=c('IFNG', 'CSF2', 'IL26', 'TNF', 'IL21', 'IL17A'), order=TRUE)

# FIGURE 3E

cd4_comp <- subset(pbmc.combined, idents=c('CD4 activated'))
Idents(cd4_comp) <- 'condition'
DotPlot(cd4_comp, features = c('HLA-DMA', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DQA1', 'HLA-DQB1', 'RDH10', 'RHOC', 'PPIA', 'GAPDH', 'PLAAT4', 'SELL', 'TCF7')) + xlab('Gene') +  ylab('Cluster') + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + scale_colour_viridis(option="mako") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

# FIGURE 4A

gd <- subset(mtbvac, idents = c("1", "6", "8"))
gd <- DietSeurat(gd)
gd <- ScaleData(gd, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = rownames(gd), block.size = 100)
gd <- RunPCA(gd, verbose = FALSE)
gd <- FindNeighbors(gd, dims = 1:15)
gd <- FindClusters(gd, resolution = 0.5)
gd <- RunUMAP(gd, dims = 1:15, spread = 5, min.dist = 0.05)
DimPlot(gd, cols=c('#507AA8', '#F28E2D', '#E25859', '#77B7B2', '#58A14F', '#EDC947', '#B07AA1'), label=TRUE)

# FIGURE 4B

VlnPlot(gd, features=c('CD27), cols=c('#507AA8', '#F28E2D', '#E25859', '#77B7B2', '#58A14F', '#EDC947', '#B07AA1'))

# FIGURE 4C

FeaturePlot(gd, features=c('IFNG', 'IL2RA'), blend=TRUE, order=TRUE, cols=c('red', 'blue'))

# FIGURE 4D

gd_comp <- subset(pbmc.combined, idents=c('gd T cell'))
Idents(gd_comp) <- 'condition'
DotPlot(gd_comp, features = c('GZMM', 'CCL5', 'ISG20', 'IFI44L', 'IFITM1', 'RSAD2', 'HAVCR2', 'GZMB', 'GNLY', 'IFNG', 'TNFRSF18', 'TNFRSF9', 'IL2RA', 'IL21R', 'FURIN', 'CD82', 'BATF', 'DUSP4', 'TFRC')) + xlab('Gene') +  ylab('Cluster') + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + scale_colour_viridis(option="mako") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

# FIGURE 4E

VlnPlot(gd_comp, features=c('FCGR3A'), cols=c('#941F81', '#FCEA1C'))

# FIGURE 4F

poscells <- WhichCells(object = gd, expression = FCGR3A > 0)
gd$gd_logical <- ifelse(colnames(gd) %in% poscells, "Pos", "Neg")
DimPlot(gd, group.by='gd_logical', cols=c('grey', 'purple'))
Idents(gd) <- 'gd_logical'
gd_markers <- FindMarkers(gd, ident.1='Pos', ident.2='Neg', min.pct=0.25)
volcano_gd <- data.frame(row.names=row.names(gd_markers), log2FoldChange=gd_markers$avg_log2FC, pvalue=gd_markers$p_val_adj)
EnhancedVolcano(volcano_gd, lab = rownames(volcano_gd), x = 'log2FoldChange', y = 'pvalue', title = 'Neg vs Pos', pCutoff = 0.05, FCcutoff = 0.58, pointSize = 3.0, labSize = 2.0, col=c('black', 'black', 'black', 'red3'), colAlpha = 1)

# FIGURE 5A

nk <- subset(mtbvac, idents=c('7'))
nk <- DietSeurat(nk)
nk <- ScaleData(nk, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = rownames(nk))
nk <- RunPCA(nk, verbose = FALSE)
nk <- FindNeighbors(nk, dims = 1:15)
nk <- FindClusters(nk, resolution = 0.5)
nk <- RunUMAP(nk, dims = 1:15, spread = 1, min.dist = 0.05)
DimPlot(nk, cols=c('pink1', 'lightblue4'), label=TRUE)

# FIGURE 5B

nk_markers <- FindMarkers(nk, ident.1='1', ident.2='2', min.pct=0.25)
volcano_nk <- data.frame(row.names=row.names(nk_markers), log2FoldChange=nk_markers$avg_log2FC, pvalue=nk_markers$p_val_adj)
EnhancedVolcano(volcano_nk, lab = rownames(volcano_nk), x = 'log2FoldChange', y = 'pvalue', title = 'SC1 vs SC2', pCutoff = 0.05, FCcutoff = 0.58, pointSize = 3.0, labSize = 2.0, col=c('black', 'black', 'black', 'red3'), colAlpha = 1)

# FIGURE 5C

sceasy::convertFormat(nk, from="seurat", to="anndata", outFile='scanpy_nk.h5ad')
# Python in Spyder in Anaconda
import scanpy as sc
adata =sc.read_h5ad("scanpy_nk.h5ad")
markers=['FCER1G', 'NCAM1', 'NCR1', 'KLRB1', 'KLRC1', 'KLRC2', 'KLRD1', 'KLRK1', 'FCGR3A', 'CD226', 'KIR2DL1', 'KIR2DL4', 'IFNG', 'CSF1', 'IL2RA', 'IL12RB1', 'IL12RB2', 'IL21R', 'FASLG', 'TNFSF10', 'GZMA', 'GZMB', 'GZMH', 'GZMK', 'GZMM']
sc.pl.matrixplot(adata, markers, 'cluster', cmap='Reds', colorbar_title='scaled\nexpression', swap_axes=True, var_group_rotation=True, save='.pdf')

# FIGURE 5D

nk_comp <- subset(pbmc.combined, idents=c('NK cell'))
Idents(nk_comp) <- 'condition'
DotPlot(nk_comp, features = c('TNFRSF18', 'IGFBP7', 'GRINA', 'NFKBIA', 'BATF', 'XCL1', 'XCL2', 'GZMB', 'BZW1', 'GSN')) + xlab('Gene') +  ylab('Cluster') + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + scale_colour_viridis(option="mako") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

# FIGURE 6B

FeaturePlot(mtbvac, features=c('IL6', 'IFNG', 'TNF', 'IL17A', 'IL17F', 'CXCL10', 'IL2', 'IL22', 'IL4'), order=TRUE)
