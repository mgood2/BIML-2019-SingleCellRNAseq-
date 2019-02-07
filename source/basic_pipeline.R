# BIML data processing

# 1. load data (output of 10X genomics)
library(DropletUtils)

PATH_raw = "A:/Conference/BIML_2019/data/RawData"

sce <- read10xCounts(PATH_raw)
my.counts <- sce@assays$data$counts

# 2. Filtering empty droplet
set.seed(100)
e.out <- emptyDrops(my.counts)
is.cell <- e.out$FDR < 0.01
sum(is.cell, na.rm=TRUE)


## Quality control on the cells ##
# MT gene QC
library(scater)
library(EnsDb.Hsapiens.v86)

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(sce)$CHR <- location
summary(location=="MT")
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")))

# par(mfrow=c(1,3))

hist(sce$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(sce$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(sce$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

# Dimensionality reduction plots
sce <- runPCA(sce, use_coldata=TRUE)
# reducedDimNames(sce)
plotReducedDim(sce, use_dimred = "PCA_coldata")

library(RColorBrewer)
ggplot(as.data.frame(sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = sce$pct_counts_Mito)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
  geom_point(size = 0.1) + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = sce$log10_total_counts)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
  geom_point(size = 0.1) + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

filtering <- numeric()
filtering[which(sce$log10_total_counts > 3 & sce$pct_counts_Mito <= 10)] <- 1
filtering[which(sce$log10_total_counts <= 3 | sce$pct_counts_Mito > 10)] <- 0
sce$filtering <- as.integer(filtering)
print (table(sce$filtering))
filteredSCEset <- sce[, which(sce$filtering==1)]

ggplot(as.data.frame(sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(sce$filtering))) +
  # scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
  geom_point(size = 1) + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

library(scran)
clusters <- quickCluster(filteredSCEset, method="igraph", min.mean=0.1,
                         irlba.args=list(maxit=1000)) # for convergence.
table(clusters)

filteredSCEset <- computeSumFactors(filteredSCEset, min.mean=0.1, cluster=clusters)
summary(sizeFactors(filteredSCEset))

filteredSCEset <- normalize(filteredSCEset)

rownames(filteredSCEset@assays$data$counts) <- rownames(filteredSCEset)
colnames(filteredSCEset@assays$data$counts) <- colnames(filteredSCEset)

rownames(filteredSCEset@assays$data$logcounts) <- rownames(filteredSCEset)
colnames(filteredSCEset@assays$data$logcounts) <- colnames(filteredSCEset)

### filtering genes by rowmean==0
library(Matrix)
# keep_feature <- rowMeans(filteredSCEset@assays$data$logcounts)!=0
keep_feature <- rowSums(filteredSCEset@assays$data$logcounts != 0) > 3

filteredSCEset <- filteredSCEset[keep_feature, ]

### Select Highly variable genes (feature selection)
var.fit <- trendVar(filteredSCEset, parametric=TRUE, use.spikes=FALSE)#, design = batchDesign)
var.out <- decomposeVar(filteredSCEset, var.fit)
hvg <- var.out[which(var.out$FDR < 0.05 & var.out$bio > .01),]
dim(hvg)

colnames(filteredSCEset) <- filteredSCEset$Barcode
colnames(filteredSCEset@assays$data$counts) <- filteredSCEset$Barcode
colnames(filteredSCEset@assays$data$logcounts) <- filteredSCEset$Barcode

saveRDS(filteredSCEset, file = "A:/Conference/BIML_2019/data/filteredSCEset/filteredSCEset.rds")

library(Seurat)
Seuratset <- CreateSeuratObject(filteredSCEset@assays$data$counts)
Seuratset@assays$RNA@var.features <- rownames(hvg)
Seuratset@assays$RNA@data <- filteredSCEset@assays$data$logcounts
Seuratset <- ScaleData(object = Seuratset, 
                       features = rownames(Seuratset@assays$RNA@data))

##### Dimension reduction


PCA = 50
Seuratset <- RunPCA(Seuratset, pcs.compute = PCA, weight.by.var = FALSE)
plot(Seuratset@reductions$pca@stdev,
     xlab = "PC",
     ylab = "Eigenvalue")



Seuratset <- FindNeighbors(object = Seuratset)
Seuratset <- FindClusters(Seuratset, random.seed = 123, 
                          verbose = TRUE, 
                          n.start = 10000, n.iter = 10000)

Seuratset <- RunTSNE(Seuratset, dims.use = 1:PCA, do.fast = T, 
                     features = rownames(hvg), seed.use = 12345, 
                     perplexity=100)


Seuratset <- RunUMAP(Seuratset, reduction.use = "pca", dims.use = 1:5)

markerGenes <- FindAllMarkers(Seuratset)
head(markerGenes[order(markerGenes[,"cluster"], markerGenes[,"p_val_adj"], -markerGenes[,"avg_logFC"]), ])
markerGenes_reordered <- markerGenes[order(markerGenes[,"cluster"], markerGenes[,"p_val_adj"], -markerGenes[,"avg_logFC"]), ]

heatmap_features <- markerGenes_reordered %>%
                    group_by(cluster) %>%
                    top_n(n = 5, wt = avg_logFC)

heatmap_features$gene

saveRDS(Seuratset, file = "A:/Conference/BIML_2019/data/Seurat/Seuratset.rds")

DoHeatmap(Seuratset, features = heatmap_features$gene)
DimPlot(Seuratset, reduction = "pca", label = TRUE)
DimPlot(Seuratset, reduction = "tsne", label = TRUE)
DimPlot(Seuratset, reduction = "umap", label = TRUE)

# Violin and Ridge plots
VlnPlot(object = Seuratset, features = c("ENSG00000112486", "ENSG00000111796"), pt.size = 0)
RidgePlot(object = Seuratset, features = c("ENSG00000112486", "ENSG00000111796"))
Seuratset$groups <- Seuratset@active.ident
FeaturePlot(object = Seuratset, features = c("ENSG00000112486", "ENSG00000111796"), blend = TRUE)
DotPlot(object = Seuratset, features = heatmap_features$gene)
