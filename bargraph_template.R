library(stringr)
library(ggplot2)
install.packages("gghighlight")
library(gghighlight)

GSEA <- read_tsv("C:/Users/leonx/Documents/nhsee/rev_by_CH25H_deletion.tsv") #10obs
#GSEA <- GSEA[,-1]
colnames(GSEA)[1] <- c("Name")
colnames(GSEA)[7] <- c("FDR")
colnames(GSEA)[4] <- c("k")
GSEA$FDR <- as.numeric(GSEA$FDR)
#GSEA$FDR <- as.numeric(GSEA$FDR)
#GSEA <- GSEA[GSEA$X.log.p.value.>0,] #309 pathways
GSEA$FDR <- -log10(GSEA$FDR)
GSEA$enrich <- paste(GSEA$k, "/", GSEA$K, sep=" ")
#GSEA <- GSEA[order(-GSEA$logP),]
GSEA$Name <-  gsub("HALLMARK_", "", GSEA$Name)
GSEA$Name <- gsub("*_", " ", GSEA$Name)
GSEA$Name <- str_to_title(GSEA$Name)

GSEA <- GSEA[-c(1,3,7,10,12,13,15,18,19,22),]

#GSEA$z.score <- (GSEA$k.K - mean(GSEA$k.K))/sd(GSEA$k.K)
GSEA$NES <- GSEA$z.score
GSEA$threshold <- c("NC")
GSEA$threshold [GSEA$NES > 0 ] <- c("UP")
GSEA$threshold [GSEA$NES < 0 ] <- c("DOWN")
#marker$colours[marker$log2FoldChange >= 0.2 & marker$padj <= 0.1] <- c("UP")
GSEA$Name <-  factor(GSEA$Name, levels=rev(GSEA$Name))
# neg <- GSEA[order(GSEA$NES)[1:10],]
neg <- GSEA[order(GSEA$NES),]
neg <- GSEA[GSEA$NES<0,]
#choosing pathway
neg <- neg[neg$Name %in% c('Pathogen Induced Cytokine Storm Signaling Pathway','TREM1 Signaling','Neuroinflammation Signaling Pathway', 'Interferon Signaling', 'Toll-like Receptor Signaling'),]

neg <- neg[order(neg$X.log.p.value., decreasing = TRUE),]
pos <- GSEA [GSEA$NES>0,]#242 pathways
#choosing pathway
pos <- pos[pos$Name %in% c('Pathogen Induced Cytokine Storm Signaling Pathway','TREM1 Signaling','VDR/RXR Activation', 'Apoptosis Signaling', 'Activation of IRF by Cytosolic Pattern Recognition Receptors', 'Complement System',  'Interferon Signaling'),]
pos <- pos[order(pos$NES, decreasing = TRUE)[1:20],]

# GSEA bargraph
GSEA <- rbind(neg,pos)
GSEA <- neg
GSEA <- GSEA %>% arrange(threshold,logP) %>% mutate(Name=factor(Name, levels=Name))
GSEA<- head(GSEA,10)
#check numeric or character
class(neg$NES) 
neg$NES <- as.numeric (neg$NES)
neg$NES <- abs(neg$NES)

GSEA_select <- GSEA[GSEA$Name %in% c('Complement','Tgf Beta Signaling','Inflammatory Response','Interferon Gamma Response'),]
pathways_of_interest <- subset(GSEA, GSEA$Name %in% c("P53 Pathway", "Apoptosis", "Inflammatory Signaling", "Hypoxia"))
GSEA$colors <- "purple1"
GSEA$colors[GSEA$Name %in% c("P53 Pathway", "Apoptosis", "Inflammatory Response", "Hypoxia")] <- "cadetblue"

pdf("LG253_Taudown_reverse_by_CH25H_deletion!!.pdf", width=8, height=8)
ggplot(data = GSEA, aes(x = reorder(Name,FDR), y = FDR), fill = direction) + ylim(0 , 10)+
  theme_classic() +
  ylab("-log10(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  aes(width=0.8, fill=colors, alpha=0.5, color = "black" )) +
  scale_fill_manual(values = c("NC" = "black", "purple1" = "purple1", "cadetblue" = "cadetblue")) +
  scale_color_manual(values = c("NC" = "black", "purple1" = "purple1", "cadetblue" = "cadetblue")) + #legend
  # gghighlight(Name == "P53 Pathway" || Name == "Apoptosis" || Name == "Inflammatory Signaling") +
  #geom_col(aes(fill = enrichmentScore), color = 'black', lwd = 0.3) + 
  #scale_fill_gradient2(low = "blue",
  #             high = "red",
  #            mid = "white",
  #           guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.8)) +
  # coord_flip() +
  # geom_text(aes(label=enrich), vjust=0.4, 
  #           hjust=0, size=5, color="black", 
  #           stat="Identity", y= GSEA$FDR + 0.1) +
  theme(aspect.ratio = 0.8, axis.text.x = element_text(angle = -45, hjust = 0)) +  
  ggtitle("Pathways enriched by Differentially Expressed Genes") + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(colour = "black", size=16), axis.text.y = element_text(colour = "black", size=12)) + #changing Y and x Axis 
  geom_hline(yintercept=c(1),lty=4,col="black",lwd=0.5) 
dev.off()


pathways <- read_tsv("C:/Users/leonx/Documents/nhsee/rev_by_CH25H_deletion.tsv")
View(pathways)

pathways$Name <- gsub("HALLMARK", "", x = pathways$Name)

pathways$frac <- paste(pathways$k, "/", pathways$K, sep = " ")

pdf("LG253_rescued_by_CH25H.pdf", width = 8, height = 7) 
ggplot(data = pathways, aes(x = reorder(Pathway, FDR), y = FDR), fill = direction + ylim(0, 10) +
         theme_classic() + 
         ylab("-log10(FDR)") + xlab(NULL) +
         geom_bar(stat = "frac"), width = 0.8, fill = "purple1", alpha = 0.5, color = "black" +
         coord_flip() + geom_text(aes(label=frac), vjust = 0.4, hjust = 0, size = 5, color = "black",
                                  stat = "Pathway", y = pathways$FDR + 0.1) +
         theme(aspect.ratio = 0.8) + 
         ggtitle("CH25HKO_Tau_vs_Tau") + 
         theme(plot.title = element_Text(hjust = 0.5), axis.text.x = element_text(colour = "black", size = 12), axis.text.y = element_text(colour = "black", size = 12)) +
         geom_hline(yintercept = c(1), lty = 4, col = "black", lwd = 0.5))
ggsave("C:/Users/leonx/Documents/nhsee/LG253_rescued_by_CH25H.pdf", plot = bar, device = "pdf", width = 8, height = 6, units = "in")

bar
