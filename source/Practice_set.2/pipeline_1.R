# BIML data processing

# 1. load data (output of 10X genomics)
# This dataset passed all QC filters
# Just use the dataset for basic and downstream analysis

library(DropletUtils)

PATH_raw = "/BiO/data/scRNAseq_dataset/memory_Tcell/RawData/"

sce <- read10xCounts(PATH_raw)
colnames(sce) <- colData(sce)$Barcode

# Data Normalization
library(scran)
clusters <- quickCluster(sce, method="igraph") # for convergence.
table(clusters)

sce <- computeSumFactors(sce, cluster=clusters)
summary(sizeFactors(sce))

sce <- normalize(sce)

rownames(sce@assays$data$counts) <- rownames(sce)
colnames(sce@assays$data$counts) <- colnames(sce)

rownames(sce@assays$data$logcounts) <- rownames(sce)
colnames(sce@assays$data$logcounts) <- colnames(sce)

### filtering genes by rowmean==0
library(Matrix)
# keep_feature <- rowMeans(filteredSCEset@assays$data$logcounts)!=0
keep_feature <- rowSums(sce@assays$data$logcounts != 0) > 3

sce <- sce[keep_feature, ]

### Select Highly variable genes (feature selection)
var.fit <- trendVar(sce, parametric=TRUE, use.spikes=FALSE)#, design = batchDesign)
var.out <- decomposeVar(sce, var.fit)
hvg <- var.out[which(var.out$FDR < 0.05),] # var.out$bio > .01
dim(hvg)

plot(y=var.out$total, x=var.out$mean, pch=16, cex=0.3,
     ylab="variance of log expression", xlab="mean log expression")
o <- order(var.out$mean)
lines(y=var.out$tech[o], x=var.out$mean[o], col="dodgerblue", lwd=2)
points(y=var.out$total[var.out$FDR < 0.05],
       x=var.out$mean[var.out$FDR < 0.05],
       pch=16, cex=0.3, col="red")

library(Seurat)
Seuratset <- as.seurat(sce)
Seuratset@var.genes <- rownames(hvg)
Seuratset <- ScaleData(object = Seuratset)

##### Dimension reduction
PCA = 50
Seuratset <- RunPCA(object = Seuratset, pc.genes = Seuratset@var.genes, do.print = FALSE, pcs.compute = PCA, weight.by.var = FALSE)
plot(Seuratset@dr$pca@sdev)

Seuratset <- RunTSNE(Seuratset, dims.use = 1:15,
                     do.fast = T, seed.use = 42, perplexity=30)


Seuratset <- RunUMAP(object = Seuratset, reduction.use = "pca", dims.use = 1:15, min_dist = 0.75)

# Clustering
Seuratset <- FindClusters(Seuratset, reduction.type="pca",
                          dims.use = 1:15, save.SNN = TRUE,
                          force.recalc = TRUE)


# Visualization
PCAPlot(Seuratset)
TSNEPlot(Seuratset, do.label = TRUE, pt.size = 0.05)
DimPlot(Seuratset, reduction.use = "pca", do.label = TRUE)
DimPlot(Seuratset, reduction.use = "tsne", do.label = TRUE)
DimPlot(Seuratset, reduction.use = "umap", do.label = TRUE)

markerGenes <- FindAllMarkers(Seuratset)
markerGenes_reordered <- markerGenes[order(markerGenes[,"cluster"], markerGenes[,"p_val_adj"], -markerGenes[,"avg_logFC"]), ]

heatmap_features <- markerGenes_reordered %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)

DoHeatmap(Seuratset, genes.use = heatmap_features$gene)

# Violin and Ridge plots
VlnPlot(object = Seuratset, features.plot = c("ENSG00000148773"), point.size.use = 0)
VlnPlot(object = Seuratset, features.plot = c("ENSG00000111796"), point.size.use = 0)
VlnPlot(object = Seuratset, features.plot = c("ENSG00000126353"), point.size.use = 0)

VlnPlot(object = Seuratset, features.plot = c("ENSG00000148773","ENSG00000111796","ENSG00000126353"), point.size.use = 0)
RidgePlot(object = Seuratset, features = c("ENSG00000148773", "ENSG00000111796", "ENSG00000126353"))
Seuratset$groups <- Seuratset@active.ident
FeaturePlot(object = Seuratset, features = c("ENSG00000112486", "ENSG00000111796"), blend = TRUE)
DotPlot(object = Seuratset, genes.plot = heatmap_features$gene, x.lab.rot = TRUE, plot.legend = TRUE)

# marker gene expression pattern
library(ggplot2)
library(dplyr)
TGgene = "ENSG00000126353"
TGgeneExpression = Seuratset@data[TGgene, ]
tibble(x = Seuratset@dr$tsne@cell.embeddings[,1],
       y = Seuratset@dr$tsne@cell.embeddings[,2],
       TGgeneExpression = Seuratset@data[TGgene, ]) %>%
  ggplot(aes(x=x, y=y, colour=TGgeneExpression)) +
  geom_point(size=0.01) +
  ggtitle("CCR7") +
  # scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  # scale_colour_gradient(low = "grey", high = "red") +
  scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
                         breaks=c(0,max(TGgeneExpression)),
                         labels=c(0,round(as.numeric(max(TGgeneExpression)), digits = 2))) +
  ylab("Component 2") +
  xlab("Component 1") +
  theme_bw() +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        # legend.text=element_text(size=10),
        legend.title=element_blank(),
        # legend.key=element_blank(),
        axis.text.x = element_text(size=10)
  )

# Find markergenes by each clusters (find cluster specific marker genes)
markerGenes <- FindAllMarkers(Seuratset)
head(markerGenes[order(markerGenes[,"cluster"], markerGenes[,"p_val_adj"], -markerGenes[,"avg_logFC"]), ])
markerGenes_reordered <- markerGenes[order(markerGenes[,"cluster"], markerGenes[,"p_val_adj"], -markerGenes[,"avg_logFC"]), ]

# select gene for draw heatmap (top3 gene per each clusters)
heatmap_features <- markerGenes_reordered %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)

# various visualization method
DoHeatmap(Seuratset, features = heatmap_features$gene)
DimPlot(Seuratset, reduction = "pca", label = TRUE)
DimPlot(Seuratset, reduction = "tsne", do.label = TRUE)
DimPlot(Seuratset, reduction = "umap", do.label = TRUE)

# Violin and Ridge plots for express gene expression patterns
VlnPlot(object = Seuratset, features = c("ENSG00000112486", "ENSG00000111796"), pt.size = 0)
RidgePlot(object = Seuratset, features = c("ENSG00000112486", "ENSG00000111796"))
Seuratset$groups <- Seuratset@active.ident
FeaturePlot(object = Seuratset, features = c("ENSG00000112486", "ENSG00000111796"), blend = TRUE)
DotPlot(object = Seuratset, features = heatmap_features$gene)
