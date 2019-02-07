gene_annotation <- as.data.frame(rownames(Seuratset@assays$RNA@data))
colnames(gene_annotation) <- "ENSG_ID"
rownames(gene_annotation) <- gene_annotation$ENSG_ID

library(monocle)
set.seed(12345)

pd <- new("AnnotatedDataFrame", data = Seuratset@meta.data)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(Seuratset@assays$RNA@counts, phenoData = pd, featureData = fd,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 5))

# Log-transform each value in the expression matrix.
L <- log(exprs(cds[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
library(reshape)
library(reshape2)
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(UMI_count)") +
  ylab("Density")


set.seed(123)
cds <- reduceDimension(cds, max_components = 2, num_dim = 50, norm_method = 'log',
                       reduction_method = 'tSNE', verbose = T)

cds$clusters <- Seuratset@active.ident


markerGenes <- FindAllMarkers(Seuratset)

disp_table <- dispersionTable(cds)
ordering_genes <- markerGenes[markerGenes$p_val_adj < 0.05, ]$gene
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
cds$sample_name <- NULL

plot_cell_trajectory(cds, color_by = "clusters")

identical(colnames(cds@reducedDimS), rownames(Seuratset@meta.data))

library(ggplot2)

TGgene = "ENSG00000112486"
TGgeneExpression = Seuratset@assays$RNA@data[TGgene, ]
tibble(x = cds@reducedDimS[1,],
       y = cds@reducedDimS[2,], 
       TGgeneExpression = Seuratset@assays$RNA@data[TGgene, ]) %>% 
  ggplot(aes(x=x, y=y, colour=TGgeneExpression)) + 
  geom_point(size=0.01) + 
  # ggtitle(TGgene) +
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
        # legend.title=element_blank(), 
        # legend.key=element_blank(), 
        axis.text.x = element_text(size=10) 
  ) 

pseudotime_de <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 8)

sig_gene_names <- row.names(subset(pseudotime_de, qval < 0.1))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = F)
