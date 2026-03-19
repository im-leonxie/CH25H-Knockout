library(ggplot2)
# Install if not already installed
install.packages("ggrepel")

# Load the package
library(ggrepel)

# Install if not already installed
install.packages("ggpubr")

# Load the package
library(ggpubr)
#correlationmap
library(readr)
#comparison with DAM vs MG cluster1 (Figure 4F) ======
tau_deg <- read.csv("C:/Users/leonx/Documents/tau_markers_microglia.csv",header=T) # DAM genes published 3044obs
ch25h_deg <- read.csv("C:/Users/leonx/Documents/ch25h_markers_microglia.csv",header=T) # 1190obs
head(tau_deg)
head(ch25h_deg)
#damstage<-read.csv("Damstages.csv",header=T)
#mg <- LID_MG2_vs_MG1_DEGs

#dam<- subset(dam,dam$PVALUE>-log10(0.05)) # 1693 genes

#dam<- subset(dam,dam$PVALUE<0.05) # 136 observations
#mg<- subset(mg, mg$PVALUE<0.05 & abs(mg$LOG2FC) > 0.1375) # 141 genes 10% change
#damstage<-subset(damstage,X.log10.p.value..Mann.Whitney>-log10(0.05))
#mg <- mg[!duplicated(mg$X), ]

#rownames(mg) <- mg$X
tau_deg$gene <- tau_deg$X
ch25h_deg$gene <- ch25h_deg$X

#what happens if you just graph both?
# 
# tau_deg$ch25h_logFC <- ch25h_deg$avg_log2FC[match(tau_deg$gene, ch25h_deg$gene)]

ch25h_deg$tau_logFC <- tau_deg$avg_log2FC[match(ch25h_deg$gene, tau_deg$gene)]

# -> mg$damlogFC <- dam$avg_log2FC[match(mg$gene,dam$gene)] # match the DEGs in dam group to mg group
#mg$stagelogFC <- damstage$FC[match(mg$X,damstage$UNIQUD)]

tau_up <- tau_deg[tau_deg$avg_log2FC > 0.1, ]
ch25h_up <- ch25h_deg[ch25h_deg$avg_log2FC > 0.1, ]

tau_down <- tau_deg[tau_deg$avg_log2FC < -0.1,]
ch25h_down <- ch25h_deg[ch25h_deg$avg_log2FC < -0.1,]

# -> mg.up<-mg[mg$avg_log2FC>0.1,]
# -> dam.up<-dam[dam$avg_log2FC>0.1,]

#overlap.up.down<-mg.down[mg.down$gene %in% dam.up$gene,]

#ch25h is the first reading
overlap_down_up <- ch25h_down[ch25h_down$gene %in% tau_up$gene,]
overlap_up_down <- ch25h_up[ch25h_up$gene %in% tau_down$gene,]

# -> mg.down<-mg[mg$avg_log2FC<(-0.1),]
# -> dam.down<-dam[dam$avg_log2FC<(-0.1),]

#*does* investigate the up/down vs down/up

ch25h_deg$color <- "NC"

#red labeled: ch25h brings them down while tau up, blue labeled: ch25h brings them up while tau down; red are the ones we studied

ggplot(mg_inflam, aes(x = X, y = , fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "ggplot2 Heatmap")


ch25h_deg$color[ch25h_deg$gene %in% overlap_down_up$X] <- "red"
ch25h_deg$color[ch25h_deg$gene %in% overlap_up_down$X] <- "blue"

IR <- read.csv("C:/Users/leonx/Documents/inflammatory response.csv", header = T)
write.csv(ch25h_deg, "EN scatterplot raw data.csv")

inflam <- ch25h_deg[ch25h_deg$X %in% colnames(IR),]
inflam <- inflam[inflam$X %in% rownames(tau_deg),]

NEU_AP <- read.csv("C:/Users/leonx/Documents/neuron apoptosis.csv")

IFNALPHA <- read.csv("C:/Users/leonx/Downloads/HALLMARK_INTERFERON_ALPHA_RESPONSE.csv")
IFNGAMMA <- read.csv("C:/Users/leonx/Downloads/REACTOME_INTERFERON_GAMMA_SIGNALING.csv")
IFN <- read.csv("C:/Users/leonx/Downloads/REACTOME_INTERFERON_SIGNALING.csv")

mg_inflam <- ch25h_deg[ch25h_deg$X %in% colnames(IFNALPHA),]
mg_inflam <- mg_inflam[mg_inflam$X %in% rownames(tau_deg),]


#===========================================================
pdf(file = "C:/Users/leonx/Documents/nhsee/neuronal apoptosis reversed.pdf", width = 8, height = 8)
ggplot(data = mg_inflam, aes(x = avg_log2FC, y = tau_logFC)) +
  geom_point(aes(color = color, fill = color), shape = 16, size = 2) + #add pts
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("NC" = "black", "red" = "salmon", "blue" = "dodgerblue")) +
  #scale_fill_manual(values = c("NC" = "black", "red" = "salmon", "blue" = "dodgerblue"))
  #geom_text_repel(subset(ch25h_deg, color == "red"), mapping = aes(label = gene), max.overlaps = 20) +
  
  # geom_text_repel(data = three_genes, mapping = aes(label = gene), max.overlaps = 20, color = "indianred4") +
  # geom_point(data = three_genes, color = "indianred4", shape = 16, size = 2.5, alpha = 1) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +  
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+  
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")+
  labs(title = "Differentially Expressed Genes Plotted By Expression Change") +
  ylab("WT (no tau) vs. WT (tau)") + xlab("WT (tau) vs. KO (tau)") + xlim(-1, 1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "*`,`~")), method="pearson")   #R
dev.off()


#===========================================================

three_genes <- subset(ch25h_deg, ch25h_deg$gene %in% c("Rock1", "Smad7", "Jun", "Bax", "Fam162a", "Madd", "Dap3"))

pdf(file = "C:/Users/leonx/Documents/nhsee/reversed apoptosis fold change.pdf", width = 8, height = 8)
ggplot(data = mg_inflam, aes(x = avg_log2FC, y = tau_logFC)) +
  geom_point(aes(color = color), shape = 1, size = 2) + #add pts
  scale_color_manual(values = c("NC" = "black", "red" = "salmon", "blue" = "dodgerblue")) +
  #scale_fill_manual(values = c("NC" = "black", "red" = "salmon", "blue" = "dodgerblue"))
  #geom_text_repel(subset(ch25h_deg, color == "red"), mapping = aes(label = gene), max.overlaps = 20) +
  # geom_text_repel(data = three_genes, mapping = aes(label = gene), max.overlaps = 20, color = "indianred4") +
  # geom_point(data = three_genes, color = "indianred4", shape = 16, size = 2.5, alpha = 1) +
  theme_classic() +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +  
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+  
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")+
  ylab("Tau+_vs_Tau-") + xlab("CH25HKO_Tau+_vs_Tau+") + xlim(-1, 1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "*`,`~")), method="pearson")   #R
dev.off()

overlap.up.down<-mg.down[mg.down$gene %in% dam.up$gene,]
overlap.down.up<-mg.up[mg.up$gene %in% dam.down$gene,]
#scatterplot for correlation analysis
mg$color <- "NC"
mg$color[mg$gene %in% overlap.up.down$gene]<-"red"
mg$color[mg$gene %in% overlap.down.up$gene ]<-"blue"
write.csv(mg,"LG253_EN_scatterplot_raw_data.csv")
min(abs(mg$avg_log2FC))
mg$damlogFC
apoptosis <- read.csv("C:/Users/leonx/Documents/reversed_apoptosis_lg253.csv",header = T)
label_data <- subset(mg, mg$gene %in% apoptosis$mouse)
three_genes <- subset(mg, mg$gene %in% c("Rock1", "Smad7", "Jun"))

pdf(file="/Users/haochen/Desktop/LG253_Taudown_reverse_by_CH25H_deletion_correlation.pdf",width = 8, height = 6)
ggplot(data=mg_inflam, aes(x=avg_log2FC, y=damlogFC)) + 
  geom_point(aes(color= color),shape=1, size=3)+
  scale_color_manual(values = c("NC" = "black", "red" = "salmon", "blue" = "dodgerblue"))+
  geom_text_repel(subset(mg_inflam, color=='blue'), mapping = aes(label = gene), max.overlaps = 20)+
  # geom_text_repel(data = label_data[rownames(label_data) %notin% rownames(three_genes),], mapping = aes(label = gene), max.overlaps = 20)+
  # geom_text_repel(data = three_genes, mapping = aes(label = gene), max.overlaps = 20, color = "red")+
  # geom_point(data = label_data, color = "black", shape = 16, size = 3.5,alpha=0.8)+
  # geom_point(data = three_genes, color = "red", shape = 16, size = 3.5, alpha = 0.8) +
  #geom_smooth(method = "lm", colour='black') +
  theme_classic() +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +  
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+  
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")+
  ylab("Tau+_vs_Tau-") + xlab("CH25HKO_Tau+_vs_Tau+") + xlim(-1, 1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "*`,`~")), method="pearson")   #R
dev.off()

