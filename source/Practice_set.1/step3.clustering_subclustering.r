# filteredSCEset
library(Seurat)
Seuratset <- as.seurat(filteredSCEset)
Seuratset@var.genes <- rownames(hvg)
Seuratset <- ScaleData(Seuratset)

##### Dimension reduction
PCA = 50
Seuratset <- RunPCA(object = Seuratset, pc.genes = Seuratset@var.genes, do.print = FALSE, pcs.compute = PCA, weight.by.var = FALSE)
plot(Seuratset@dr$pca@sdev)
Seuratset <- RunTSNE(Seuratset, dims.use = 1:PCA, do.fast = T, seed.use = 42, perplexity=100)
Seuratset <- RunUMAP(object = Seuratset, reduction.use = "pca", dims.use = 1:PCA, min_dist = 0.75)
Seuratset <- FindClusters(Seuratset, reduction.type="pca", dims.use = 1:PCA, save.SNN = TRUE, force.recalc = TRUE)

PCAPlot(Seuratset)
TSNEPlot(Seuratset, do.label = TRUE, pt.size = 0.05)
DimPlot(Seuratset, reduction.use = "umap", do.label = TRUE)

##### Reculstering
library(scater)
subset_Seuratset <- SubsetData(Seuratset, ident.use = c(0:12,14))
subset_sce <- Convert(from = subset_Seuratset, to = "sce")
rownames(subset_sce) <- uniquifyFeatureNames(rowData(subset_sce)$gene, rowData(subset_sce)$Symbol)

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rowData(subset_sce)$gene

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
subset_sce <- subset_sce[G_list$ensembl_gene_id, ]
rowData(subset_sce)$symbol <- G_list$hgnc_symbol
rownames(subset_sce) <- uniquifyFeatureNames(rowData(subset_sce)$gene, rowData(subset_sce)$symbol)
keep_feature <- rownames(subset_sce)[!grepl("ENSG", rownames(subset_sce))]
subset_sce <- subset_sce[keep_feature, ]

library(scran)
clusters <- quickCluster(subset_sce, method="igraph")
subset_sce <- computeSumFactors(subset_sce, cluster=clusters)
subset_sce <- normalize(subset_sce)

### Select Highly variable genes (feature selection)
var.fit <- trendVar(subset_sce, parametric=TRUE, use.spikes=FALSE)#, design = batchDesign)
var.out <- decomposeVar(subset_sce, var.fit)
hvg <- var.out[which(var.out$FDR < 0.05 & var.out$bio > .01),]
dim(hvg)

# subset_sce
library(Seurat)
Seuratset <- as.seurat(subset_sce)
Seuratset@var.genes <- rownames(hvg)
Seuratset@scale.data <- ScaleData(Seuratset)

##### Dimension reduction
PCA = 50
Seuratset <- RunPCA(Seuratset, pcs.compute = PCA, weight.by.var = FALSE)
plot(Seuratset@dr$pca@sdev)
Seuratset <- RunTSNE(Seuratset, dims.use = 1:PCA, do.fast = T, seed.use = 42, perplexity=100)
Seuratset <- RunUMAP(object = Seuratset, reduction.use = "pca", dims.use = 1:PCA, min_dist = 0.75)
Seuratset <- FindClusters(Seuratset, reduction.type="pca", dims.use = 1:PCA, save.SNN = TRUE, force.recalc = TRUE)

PCAPlot(Seuratset)
TSNEPlot(Seuratset, do.label = TRUE, pt.size = 0.05)
DimPlot(Seuratset, reduction.use = "umap", do.label = TRUE)
