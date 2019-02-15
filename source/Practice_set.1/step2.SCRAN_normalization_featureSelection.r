# filteredSCEset <- readRDS(file = paste(PATH_PBMC_dataset, "filteredSCEset.rds", sep = "/"))

# filteredSCEset: QC completed dataset
filteredSCEset <- sce[, which(sce$filtering==1)]

# Normalization step
library(scran)
clusters <- quickCluster(filteredSCEset, method="igraph") # for large dataset (usually number of cells > 1000)
filteredSCEset <- computeSumFactors(filteredSCEset, cluster=clusters)
filteredSCEset <- normalize(filteredSCEset)

rownames(filteredSCEset@assays$data$logcounts) <- rownames(filteredSCEset)
colnames(filteredSCEset@assays$data$logcounts) <- colnames(filteredSCEset)

# plot(sizeFactors(filteredSCEset), (filteredSCEset$total_counts)/1000, log="xy",
#      ylab="Library size (kilo)", xlab = "Size factor")

### filtering genes
# remove lowly expressed genes
library(Matrix)
# keep_feature <- rowMeans(filteredSCEset@assays$data$logcounts)!=0
keep_feature <- rowSums((filteredSCEset@assays$data$logcounts) != 0) > 3
filteredSCEset <- filteredSCEset[keep_feature, ]

### Select Highly variable genes (feature selection)
var.fit <- trendVar(filteredSCEset, parametric=TRUE, use.spikes=FALSE)#, design = batchDesign)
var.out <- decomposeVar(filteredSCEset, var.fit)
hvg <- var.out[which(var.out$FDR < 0.05 & var.out$bio > .01),] # var.out$bio means biological variance (from decomposition of total variance)
dim(hvg)

saveRDS(filteredSCEset, file = "filteredSCEset.rds")
saveRDS(hvg, file = "hvg.rds")

plot(y= var.out$total, x=var.out$mean, pch=16, cex=0.3,
     ylab="Variance of log-expression", xlab="Mean log-expression")
o <- order(var.out$mean)
lines(y=var.out$tech[o], x=var.out$mean[o], col="dodgerblue", lwd=2)
points(y=var.out$total[var.out$FDR <=0.05 & var.out$bio > 0.1],
       x=var.out$mean[var.out$FDR <=0.05 & var.out$bio > 0.1],
       pch=16, cex=0.3, col="red")
