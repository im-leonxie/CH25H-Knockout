microglia <- readRDS("C:/Users/leonx/Downloads/working summer 2025/LG253_MG_reclusted_res0.3.rds")
library(Seurat)
library(SeuratWrappers)
library(dplyr)
microglia <- UpdateSeuratObject(microglia)
DefaultAssay(microglia) <- "RNA"

Idents(microglia) <-- microglia@meta.data[["Condition"]]

tau_wt <- subset(microglia, microglia@meta.data$Condition == "CH25H: +/+; P301S tau: +")
tau_ko <- subset(microglia, microglia@meta.data$Condition == "CH25H: -/-; P301S tau: +")
ctrl_wt <- subset(microglia, microglia@meta.data$Condition == "CH25H: +/+; P301S tau: -")
ctrl_ko <- subset(microglia, microglia@meta.data$Condition == "CH25H: -/-; P301S tau: -")

tau_markers <- FindMarkers(MG, ident.1 = "CH25H: +/+; P301S tau: +", ident.2 = "CH25H: +/+; P301S tau: -", min.pct = 0.1)
ch25h_markers <- FindMarkers(MG, ident.1 = "CH25H: -/-; P301S tau: +", ident.2 = "CH25H: +/+; P301S tau: +", min.pct = 0.1)

tau_markers <- tau_markers[tau_markers$avg_log2FC >= 0,]
tau_markers <- tau_markers[tau_markers$p_val <= 0.05,]
tau_markers <- tau_markers %>% filter(tau_markers, tau_markers$p_val_adj <= 0.05)
ch25h_markers <- ch25h_markers[ch25h_markers$avg_log2FC <= 0,]
ch25h_markers <- ch25h_markers[ch25h_markers$p_val <= 0.05,]


write.csv(intersect(rownames(tau_markers), rownames(ch25h_markers)), file = "corrected_mg_degs.csv")

write.csv(tau_markers, file = "tau_markers_microglia.csv")
write.csv(ch25h_markers, file = "ch25h_markers_microglia.csv")
View(ch25h_markers)
unloadNamespace("promises")
detachAll()

unlink("C:/Users/leonx/AppData/Local/R/win-library/4.5/00LOCK", recursive = TRUE)

tau_wt <- subset()