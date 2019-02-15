# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("monocle", version = "3.8")

trajectory_seurat <- SubsetData(Seuratset, ident.use = c(0:3))
trajectory_seurat@raw.data <- trajectory_seurat@raw.data[, rownames(trajectory_seurat@meta.data)]

gene_annotation <- as.data.frame(rownames(trajectory_seurat@raw.data))
colnames(gene_annotation) <- "ENSG_ID"
rownames(gene_annotation) <- gene_annotation$ENSG_ID

library(monocle)
set.seed(12345)

pd <- new("AnnotatedDataFrame", data = trajectory_seurat@meta.data)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(trajectory_seurat@raw.data, phenoData = pd, featureData = fd,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 5))

set.seed(123)
cds <- reduceDimension(cds, max_components = 2, num_dim = 15, norm_method = 'log',
                       reduction_method = 'tSNE', verbose = T)

cds$clusters <- trajectory_seurat@ident

disp_table <- dispersionTable(cds)

trajectory_sceset <- Convert(from = trajectory_seurat, to = "sce")

### Select Highly variable genes (feature selection)
var.fit <- trendVar(trajectory_sceset, parametric=TRUE, use.spikes=FALSE)#, design = batchDesign)
var.out <- decomposeVar(trajectory_sceset, var.fit)
trajectory_hvg <- var.out[which(var.out$FDR < 0.05 & var.out$bio > .01),] # var.out$bio > .01
dim(trajectory_hvg)

ordering_genes <- rownames(trajectory_hvg)

cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "clusters")

df_monocle2 <- as.data.frame(cbind(cds$Pseudotime, cds$clusters))
colnames(df_monocle2) <- c("pseudotime", "clusters")

library(ggbeeswarm)
ggplot(df_monocle2,
       aes(x = pseudotime,
           y = cds$clusters, colour = cds$clusters)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.title=element_blank(),
        # legend.key=element_blank(),
  ) +
  # scale_color_tableau() +
  # theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("clusters") +
  ggtitle("Cells ordered by diffusion map pseudotime")


identical(colnames(cds@reducedDimS), rownames(Seuratset@meta.data))

DimPlot(Seuratset, reduction.use = "umap", do.label = TRUE)

library(ggplot2)
library(dplyr)
TGgene = "ENSG00000126353"
TGgeneExpression = Seuratset@data[TGgene, ]
tibble(x = Seuratset@dr$umap@cell.embeddings[,1],
       y = Seuratset@dr$umap@cell.embeddings[,2],
       TGgeneExpression = Seuratset@data[TGgene, ]) %>%
  ggplot(aes(x=x, y=y, colour=TGgeneExpression)) +
  geom_point(size=0.01) +
  ggtitle("KLRB1") +
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
        # legend.title=element_blank(),
        # legend.key=element_blank(),
        axis.text.x = element_text(size=10)
  )


library(ggplot2)
library(dplyr)

TGgene = "ENSG00000112486"
# TGgeneExpression = trajectory_seurat@data[TGgene, ]
tibble(x = cds@reducedDimS[1,],
       y = cds@reducedDimS[2,],
       TGgeneExpression = trajectory_seurat@data[TGgene, ]) %>%
  ggplot(aes(x=x, y=y, colour=TGgeneExpression)) +
  geom_point(size=0.01) +
  ggtitle("CCR6") +
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
        legend.text=element_text(size=10),
        legend.title=element_blank(),
        # legend.key=element_blank(),
        axis.text.x = element_text(size=10)
  )


pseudotime_de <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 8)

sig_gene_names <- row.names(subset(pseudotime_de, qval < 0.1))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 4,
                        cores = 1,
                        show_rownames = F)
