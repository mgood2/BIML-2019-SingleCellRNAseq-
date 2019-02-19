# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SingleCellExperiment", version = "3.8")
# BiocManager::install("DropletUtils", version = "3.8")
# BiocManager::install("scater", version = "3.8")
# BiocManager::install("EnsDb.Hsapiens.v86", version = "3.8")
# BiocManager::install("scran", version = "3.8")
#
# install.packages("irlba")
# install.packages("Seurat")

PATH_PBMC_dataset = "/BiO/data/scRNAseq_dataset/PBMC_dataset/PBMC8k"

library(SingleCellExperiment)
pbmc_rawcount <- readRDS(file = paste(PATH_PBMC_dataset, "pbmc8k_rawcount.rds", sep = "/"))
sce <- SingleCellExperiment(assays = list(counts = pbmc_rawcount))
rowData(sce)$ID <- rownames(sce@assays$data$counts)

library(DropletUtils)
set.seed(0)
my.counts <- pbmc_rawcount
br.out <- barcodeRanks(my.counts)

# Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))


set.seed(100)
e.out <- emptyDrops(my.counts)
is.cell <- e.out$FDR < 0.01
sum(is.cell, na.rm=TRUE)
sce <- sce[,which(e.out$FDR < 0.01)]

bcrank <- barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=bcrank$inflection, col="darkgreen", lty=2)
abline(h=bcrank$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"),
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)


## Quality control on the cells ##
# MT gene QC
library(scater)
library(EnsDb.Hsapiens.v86)

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID,
                   column="SEQNAME", keytype="GENEID")

rowData(sce)$CHR <- location
summary(location=="MT")
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
sce <- runPCA(sce, use_coldata=TRUE)
# reducedDimNames(sce)
plotReducedDim(sce, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))
dev.off()

hist(sce$log10_total_counts,
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(sce$log10_total_features_by_counts,
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(sce$pct_counts_Mito,
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(sce$log10_total_counts > 2 & sce$pct_counts_Mito <= 10 & sce$log10_total_features_by_counts > 2)] <- 1
filtering[which(sce$log10_total_counts <= 2 | sce$pct_counts_Mito > 10 | sce$log10_total_features_by_counts <= 2)] <- 0
sce$filtering <- as.integer(filtering)
table(sce$filtering)

library(dplyr)
library(ggplot2)
library(RColorBrewer)

sp<-ggplot(as.data.frame(sce@reducedDims$PCA_coldata),
           aes(x=PC1, y=PC2, color = sce$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

sp

sp<-ggplot(as.data.frame(sce@reducedDims$PCA_coldata),
           aes(x=PC1, y=PC2, color = sce$pct_counts_Mito)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

sp


sp<-ggplot(as.data.frame(sce@reducedDims$PCA_coldata),
           aes(x=PC1, y=PC2, color = sce$log10_total_features_by_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

sp

ggplot(as.data.frame(sce@reducedDims$PCA_coldata),
       aes(x=PC1, y=PC2, color = as.factor(sce$filtering))) +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())
