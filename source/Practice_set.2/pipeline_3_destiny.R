if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("destiny", version = "3.8")

library(destiny)
library(Matrix)

trajectory_seurat <- SubsetData(Seuratset, ident.use = c(0:3))
trajectory_seurat@raw.data <- trajectory_seurat@raw.data[, rownames(trajectory_seurat@meta.data)]

trajectory_sceset <- Convert(from = trajectory_seurat, to = "sce")

### Select Highly variable genes (feature selection)
var.fit <- trendVar(trajectory_sceset, parametric=TRUE, use.spikes=FALSE)#, design = batchDesign)
var.out <- decomposeVar(trajectory_sceset, var.fit)
trajectory_hvg <- var.out[which(var.out$FDR < 0.05 & var.out$bio > .01),] # var.out$bio > .01
dim(trajectory_hvg)

ordering_genes <- rownames(trajectory_hvg)
destinySCEset <- trajectory_sceset[ordering_genes, ]
dm <- DiffusionMap(t(as.matrix(logcounts(destinySCEset))))

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  clusters = trajectory_seurat@ident)

ggplot(tmp, aes(x = DC1, y = DC2, colour = clusters)) +
  geom_point(size=0.1) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic()

destinySCEset$pseudotime_diffusionmap <- rank(-eigenvectors(dm)[,1])
destinySCEset$clusters <- trajectory_seurat@ident

library(ggbeeswarm)
ggplot(as.data.frame(colData(destinySCEset)),
       aes(x = pseudotime_diffusionmap,
           y = clusters, colour = clusters)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        # legend.title=element_blank(),
        # legend.key=element_blank(),
  ) +
  # scale_color_tableau() +
  # theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("clusters") +
  ggtitle("Cells ordered by diffusion map pseudotime")
