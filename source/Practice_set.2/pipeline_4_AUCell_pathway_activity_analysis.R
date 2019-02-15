if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AUCell", version = "3.8")

library(AUCell)
library(GSEABase)

Immunologic_genesets <- getGmt("C:/Users/JUSUNG LEE/Project/mSigDB/Immunologic_signature.gmt")
kegg_genesets <- getGmt("C:/Users/JUSUNG LEE/Project/Human_BreastCancer_single_cell_analysis/data/Total_analysis/AUC_ref_gmt/c2.cp.kegg.v6.2.symbols.gmt")

sceset <- Convert(from = Seuratset, to = "sce")

BiocManager::install("biomaRt", version = "3.8")
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rowData(sceset)$gene

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=genes, mart= mart)

sceset <- sceset[G_list$ensembl_gene_id, ]
rowData(sceset)$symbol <- G_list$hgnc_symbol
library(scater)
rownames(sceset) <- uniquifyFeatureNames(rowData(sceset)$gene, rowData(sceset)$symbol)
keep_feature <- rownames(sceset)[!grepl("ENSG", rownames(sceset))]
sceset <- sceset[keep_feature, ]

# identical(rownames(Seuratset@meta.data), colnames(sceset))

rownames(sceset@assays$data$counts) <- rownames(sceset)
rownames(sceset@assays$data$logcounts) <- rownames(sceset)

cells_rankings <- AUCell_buildRankings(sceset@assays$data$logcounts)
cells_AUC_immunologic_SigDB <- AUCell_calcAUC(Immunologic_genesets, cells_rankings)

expAUC.immunologic <- getAUC(cells_AUC_immunologic_SigDB)
z_scaled_expAUC.immunologic <- t(scale(t(expAUC.immunologic)))

head(Seuratset@meta.data)
Seuratset@meta.data$clusters <- Seuratset@ident
metadata_ordering <- Seuratset@meta.data[order(Seuratset@meta.data$clusters), ]

library(pheatmap)
pheatmap(z_scaled_expAUC.immunologic[, match(rownames(metadata_ordering), colnames(z_scaled_expAUC.immunologic))],
         labels_col = "", annotation_col = metadata_ordering[, 6:7], cluster_cols = FALSE)
