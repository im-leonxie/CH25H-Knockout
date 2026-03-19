install.packages("devtools")
library(devtools)
devtools::install_github("satijalab/seurat")
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(dplyr)
library(readr)

MG <- readRDS("C:/Users/leonx/Downloads/LG253_MG_reclusted_res0.3.rds") 
MG <- UpdateSeuratObject(MG)

EN <- readRDS("C:/Users/leonx/Downloads/2025_07_14/2025_07_14/LG253_integrated_MG_exitatoryNEURON_2025.rds")
EN <- UpdateSeuratObject(EN)

View(MG)

DefaultAssay(MG) <- "RNA"

DefaultAssay(EN) <- "RNA"

Idents(MG) <- MG@meta.data[["Condition"]]
Idents(EN) <- EN@meta.data[["Condition"]]

markers <- FindAllMarkers(EN, Condition == "CH25H: -/-; P301S tau: ")

tau_deg <- FindMarkers(MG, ident.1 = "CH25H: +/+; P301S tau: +", ident.2 = "CH25H: +/+; P301S tau: -", min.pct = 0.1)
tau_deg_EN <- FindMarkers(EN, ident.1 = "CH25H: +/+; P301S tau: +", ident.2 = "CH25H: +/+; P301S tau: -", min.pct = 0.1)
write.csv(tau_deg_EN, file = "C:/Users/leonx/Documents/nhsee/tau_diff_gene_EN_all.csv")
tau_deg_EN_up <- tau_deg_EN %>% filter(tau_deg_EN$avg_log2FC > 0)
tau_deg_EN_down <- tau_deg_EN %>% filter(tau_deg_EN$avg_log2FC < 0)

View(tau_deg_EN)
# tau_deg_EN <- tau_deg_EN %>% filter(tau_deg_EN$pct.1 >= 0.1)
# # tau_deg_EN <- tau_deg_EN %>% filter(tau_deg_EN$pct.1 >= tau_deg_EN$pct.2)
tau_deg_EN <- tau_deg_EN %>% filter(tau_deg_EN$p_val <= 0.05)
# 

rescue_deg_EN <- FindMarkers(EN, ident.1 = "CH25H: -/-; P301S tau: +", ident.2 = "CH25H: +/+; P301S tau: +", min.pct = 0.1)
write.csv(rescue_deg_EN, file = "C:/Users/leonx/Documents/nhsee/rescue_diff_gene_EN_all.csv")
rescue_deg_EN <- rescue_deg_EN %>% filter(rescue_deg_EN$avg_log2FC < 0)

rescue_deg_EN_up <- rescue_deg_EN %>% filter(rescue_deg_EN$avg_log2FC > 0)
rescue_deg_EN_down <- rescue_deg_EN %>% filter(rescue_deg_EN$avg_log2FC < 0)

deg_intersect_down_up <- as.data.frame(intersect(rownames(tau_deg_EN_down), rownames(rescue_deg_EN_up)))
write.csv(deg_intersect_down_up, file = "C:/Users/leonx/Documents/nhsee/deg_intersect_EN_down_up.csv")

View(rescue_deg_EN)
# rescue_deg_EN <- rescue_deg_EN %>% filter(rescue_deg_EN$pct.1 >= 0.1)
# # rescue_deg_EN <- rescue_deg_EN %>% filter(rescue_deg_EN$pct.1 <= rescue_deg_EN$pct.2)
rescue_deg_EN <- rescue_deg_EN %>% filter(rescue_deg_EN$p_val <= 0.05)
# View(rescue_deg_EN)



write.csv(tau_deg_EN, file = "C:/Users/leonx/Documents/nhsee/tau_diff_gene_EN.csv")
write.csv(rescue_deg_EN, file = "C:/Users/leonx/Documents/nhsee/rescue_diff_gene_EN.csv")

deg_intersect <- as.data.frame(intersect(rownames(tau_deg_EN), rownames(rescue_deg_EN)))
write.csv(deg_intersect, file = "C:/Users/leonx/Documents/nhsee/deg_intersect_EN.csv")
