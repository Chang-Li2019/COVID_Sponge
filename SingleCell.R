library(Seurat)
library(patchwork)
library(scCATCH)

mild1.data <- Read10X_h5('Data/GSM4339769_C141_filtered_feature_bc_matrix.h5')
mild2.data <- Read10X_h5('Data/GSM4339770_C142_filtered_feature_bc_matrix.h5')
severe1.data <- Read10X_h5('Data/GSM4339771_C143_filtered_feature_bc_matrix.h5')
mild3.data <- Read10X_h5('Data/GSM4339772_C144_filtered_feature_bc_matrix.h5')
severe2.data <- Read10X_h5('Data/GSM4339773_C145_filtered_feature_bc_matrix.h5')
severe3.data <- Read10X_h5('Data/GSM4339774_C146_filtered_feature_bc_matrix.h5')
healthy1.data <- Read10X_h5('Data/GSM4475048_C51_filtered_feature_bc_matrix.h5')
healthy2.data <- Read10X_h5('Data/GSM4475049_C52_filtered_feature_bc_matrix.h5')
healthy3.data <- Read10X_h5('Data/GSM4475050_C100_filtered_feature_bc_matrix.h5')
severe4.data <- Read10X_h5('Data/GSM4475051_C148_filtered_feature_bc_matrix.h5')
severe5.data <- Read10X_h5('Data/GSM4475052_C149_filtered_feature_bc_matrix.h5')
severe6.data <- Read10X_h5('Data/GSM4475053_C152_filtered_feature_bc_matrix.h5')

healthy4.data <- ReadMtx(mtx = 'Data/GSM3660650_SC249NORbal_fresh_matrix.mtx',
                         cells = 'Data/GSM3660650_SC249NORbal_fresh_barcodes.tsv',
                         features = 'Data/GSM3660650_SC249NORbal_fresh_genes.tsv')


mild1 <- CreateSeuratObject(counts=mild1.data)
mild2 <- CreateSeuratObject(counts=mild2.data)
severe1 <- CreateSeuratObject(counts=severe1.data)
mild3 <- CreateSeuratObject(counts=mild3.data)
severe2 <- CreateSeuratObject(counts=severe2.data)
severe3 <- CreateSeuratObject(counts=severe3.data)
healthy1 <- CreateSeuratObject(counts=healthy1.data)
healthy2 <- CreateSeuratObject(counts=healthy2.data)
healthy3 <- CreateSeuratObject(counts=healthy3.data)
severe4 <- CreateSeuratObject(counts=severe4.data)
severe5 <- CreateSeuratObject(counts=severe5.data)
severe6 <- CreateSeuratObject(counts=severe6.data)
healthy4 <- CreateSeuratObject(counts=healthy4.data,min.features = 100)

mild1$grp = "mild"
mild2$grp = "mild"
severe1$grp = "severe"
mild3$grp = "mild"
severe2$grp = "severe"
severe3$grp = "severe"
healthy1$grp = "healthy"
healthy2$grp = "healthy"
healthy3$grp = "healthy"
severe4$grp = "severe"
severe5$grp = "severe"
severe6$grp = "severe"
healthy4$grp = "healthy"

mild1 = subset(mild1,subset=nFeature_RNA>500)
mild2 = subset(mild2,subset=nFeature_RNA>500)
severe1 = subset(severe1,subset=nFeature_RNA>500)
mild3 = subset(mild3,subset=nFeature_RNA>500)
severe2 = subset(severe2,subset=nFeature_RNA>500)
severe3 = subset(severe3,subset=nFeature_RNA>500)
healthy1 = subset(healthy1,subset=nFeature_RNA>500)
healthy2 = subset(healthy2,subset=nFeature_RNA>500)
healthy3 = subset(healthy3,subset=nFeature_RNA>500)
severe4 = subset(severe4,subset=nFeature_RNA>500)
severe5 = subset(severe5,subset=nFeature_RNA>500)
severe6 = subset(severe6,subset=nFeature_RNA>500)
healthy4 = subset(healthy4,subset=nFeature_RNA>500)

mild1 = NormalizeData(mild1, verbose=F)
mild2 = NormalizeData(mild2, verbose=F)
severe1 = NormalizeData(severe1, verbose=F)
mild3 = NormalizeData(mild3, verbose=F)
severe2 = NormalizeData(severe2, verbose=F)
severe3 = NormalizeData(severe3, verbose=F)
healthy1 = NormalizeData(healthy1, verbose=F)
healthy2 = NormalizeData(healthy2, verbose=F)
healthy3 = NormalizeData(healthy3, verbose=F)
severe4 = NormalizeData(severe4, verbose=F)
severe5 = NormalizeData(severe5, verbose=F)
severe6 = NormalizeData(severe6, verbose=F)
healthy4 = NormalizeData(healthy4, verbose=F)

mild1 = FindVariableFeatures(mild1, selection.method="vst",nfeatures=2000)
mild2 = FindVariableFeatures(mild2, selection.method="vst",nfeatures=2000)
severe1 = FindVariableFeatures(severe1, selection.method="vst",nfeatures=2000)
mild3 = FindVariableFeatures(mild3, selection.method="vst",nfeatures=2000)
severe2 = FindVariableFeatures(severe2, selection.method="vst",nfeatures=2000)
severe3 = FindVariableFeatures(severe3, selection.method="vst",nfeatures=2000)
healthy1 = FindVariableFeatures(healthy1, selection.method="vst",nfeatures=2000)
healthy2 = FindVariableFeatures(healthy2, selection.method="vst",nfeatures=2000)
healthy3 = FindVariableFeatures(healthy3, selection.method="vst",nfeatures=2000)
severe4 = FindVariableFeatures(severe4, selection.method="vst",nfeatures=2000)
severe5 = FindVariableFeatures(severe5, selection.method="vst",nfeatures=2000)
severe6 = FindVariableFeatures(severe6, selection.method="vst",nfeatures=2000)
healthy4 = FindVariableFeatures(healthy4, selection.method="vst",nfeatures=2000)

features = SelectIntegrationFeatures(object.list = list(mild1, mild2, mild3, severe1, severe2, severe3, severe4, severe5, severe6, 
                                                        healthy1, healthy2, healthy3, healthy4))
combined.list = lapply(X=list(mild1, mild2, mild3, severe1, severe2, severe3, severe4, severe5, severe6, 
                              healthy1, healthy2, healthy3, healthy4), FUN=function(x){
                                x = ScaleData(x, features=features, verbose=FALSE)
                                x = RunPCA(x, features=features, verbose=FALSE)
                              })

# immune.anchors <- FindIntegrationAnchors(object.list = list(mild1, mild2, mild3, severe1, severe2, severe3, severe4, severe5, severe6, 
#                                                             healthy1, healthy2, healthy3, healthy4), dims = 1:20)

immune.anchors <- FindIntegrationAnchors(object.list = combined.list, anchor.features = features, reduction = "rpca")

immune.anchors_new <- FindIntegrationAnchors(object.list = combined.list, dims=1:50)
immune.combined_new <- IntegrateData(anchorset = immune.anchors_new)

immune.combined <- IntegrateData(anchorset = immune.anchors)

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

theme_set(theme_cowplot())
t.cells <- subset(immune.combined, idents = "0")
Idents(t.cells) <- "grp"
avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
avg.t.cells$gene <- rownames(avg.t.cells)

markers <- FindMarkers(object = immune.combined, ident.1 = "healthy", group.by = 'grp', ident.2 = "severe")

ggplot(as.data.frame(avg.t.cells), aes(severe, healthy)) + geom_point()

clu_markers <- findmarkergenes(immune.combined,
                               species = 'Human',
                               cluster = c(0:9),
                               match_CellMatch = T,
                               cancer = NULL,
                               tissue = "Lung",
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)


clu_ann <- scCATCH(clu_markers,
                   species = "Human",
                   cancer = NULL,
                   tissue = "Lung")

tcell = subset(x = immune.combined, idents = "4")
macrophage = subset(x = immune.combined, idents = "0")
tcell_markers <- FindMarkers(object = tcell, ident.1 = "severe", group.by = 'grp', ident.2 = "healthy")
macrophage_markers <-  FindMarkers(object = macrophage, ident.1 = "severe", group.by = 'grp', ident.2 = "healthy")
macrophage_markers_r <-  FindMarkers(object = macrophage, ident.1 = "healthy", group.by = 'grp', ident.2 = "severe")
fname = tcell_markers
fname[(fname$avg_log2FC>0) & (fname$avg_log2FC<1),]

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)