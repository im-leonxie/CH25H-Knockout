library(Seurat)
library(SeuratWrappers)
remotes::install_github('satijalab/seurat-wrappers')

astrocytes <- readRDS("C:/Users/leonx/OneDrive/Documents/R Projects/LG253_AST_reclusted_res0.3.rds")
astrocytes <- UpdateSeuratObject(astrocytes)
DefaultAssay(astrocytes) <- "RNA"
microglia <- readRDS("C:/Users/leonx/Downloads/LG253_MG_reclusted_res0.3.rds")
microglia <- UpdateSeuratObject(microglia)
DefaultAssay(microglia) <- "RNA"

AST_MG <- merge(astrocytes, microglia, merge.data = T)

Idents(microglia) <- microglia@meta.data[["Condition"]]
Idents(astrocytes) <- astrocytes@meta.data[["Condition"]]
markers_tau <- FindMarkers(microglia, ident.1 = "CH25H: +/+; P301S tau: +", ident.2 = "CH25H: +/+; P301S tau: -", min.pct = 0.1)
markers_ch25h <- FindMarkers(microglia, ident.1 = "CH25H: -/-; P301S tau: +", ident.2 = "CH25H: +/+; P301S tau: +", min.pct = 0.1)

markers_tau <- markers_tau[markers_tau$avg_log2FC > 0,]
markers_ch25h <- markers_ch25h[markers_ch25h$avg_log2FC < 0,]

write.csv(intersect(rownames(markers_tau), rownames(markers_ch25h)), "fuhuhluhtoogan_mg.csv")

#-----------------------------------------------------------------------------------------------------------------------------------------

library(monocle3)
library(grr)
library(SeuratWrappers)
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))

remotes::install_github("bnprks/BPCells/r")
BiocManager::install('grr')
install.packages('BiocManager')

BiocManager::install("monocle3")

devtools::install_github('cole-trapnell-lab/monocle3')

cds_mg <- as.cell_data_set(microglia)

cds_mg <- preprocess_cds(cds_mg, num_dim = 100)
cds_mg <- reduce_dimension(cds_mg)
cds_mg <- cluster_cells(cds_mg)


marker_test_res <- top_markers(cds_mg, group_cells_by="Condition", 
                               reference_cells=1000, cores=8)
View(marker_test_res)
cds_mg <- learn_graph(cds_mg)
cds_mg <- order_cells(cds_mg)








cds_ast <- as.cell_data_set(astrocytes)

cds_ast <- preprocess_cds(cds_ast, num_dim = 100)
cds_ast <- reduce_dimension(cds_ast)
cds_ast <- cluster_cells(cds_ast)

plot_cells(cds_ast, color_cells_by = "Condition")

marker_test_res <- top_markers(cds_ast, group_cells_by="Condition", 
                               reference_cells=1000, cores=8)
View(marker_test_res)
cds_ast <- learn_graph(cds_ast)
cds_ast <- order_cells(cds_ast)

plot_cells(cds_ast, color_cells_by = "pseudotime")
