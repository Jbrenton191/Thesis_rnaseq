 library(tximport)
library(DESeq2)
library(tidyverse)
library(magrittr)
library(DEGreport)

# # changed 11L to 11LL and 2R and 2RR in directory titles

txi.salmon <- readRDS("~/jb/post_viva_work/txi_salmon.rds")
sampleMetadata<-readRDS("~/jb/post_viva_work/sampleMetadata_deseq.rds")
mus_tx2id<-readRDS("~/jb/post_viva_work/mus_tx2id_deseq.rds")

dds <- DESeqDataSetFromTximport(txi.salmon,
                                colData=sampleMetadata,
                                # design=~Breeder+Age+Group)
                                # design=~Group+Age+Group:Age)
                                design=~Breeder+RIN+Group+Age+Group:Age)

dds$Group<-relevel(dds$Group, ref = "Low-Fat Diet")
dds$Age<-relevel(dds$Age, ref = "6 Months")
dds$Age_Group<-factor(x = dds$Age_Group, levels = c("6 Months-Low-Fat Diet", "6 Months-High-Fat Diet", "18 Months-Low-Fat Diet", "18 Months-High-Fat Diet"))
dds$Age_Group<-relevel(x = dds$Age_Group, ref = "6 Months-Low-Fat Diet")

filter_count <- DEGreport::degFilter(counts = counts(dds),
                                     metadata = as.data.frame(colData(dds)),
                                     group = "Age_Group",
                                     min = 1, # All samples in group must have more than expr > 0
                                     minreads = 0)
dds <- dds[rownames(filter_count),]
dds <- DESeq(dds)

res<-results(dds)

res2 <- DESeq2::results(dds, name = "Age_18.Months_vs_6.Months")
res3 <- DESeq2::results(dds, name = "Group_High.Fat.Diet_vs_Low.Fat.Diet")


res$gene <- row.names(res)
resOrdered<-res[order(res$padj),]
resOrdered$gene <- row.names(resOrdered)
resOrdered <- as.data.frame(resOrdered)

sigGenes <- subset(resOrdered, padj < 0.05)
sigGenes <- as.data.frame(sigGenes)
sigGenes$regulation<-"Upregulated"
sigGenes$regulation[sigGenes$log2FoldChange<0]<-"Downregulated"
sigGenes2<-mus_tx2id[which(mus_tx2id$GENE_ID %in% sigGenes$gene),]
sigGenes2 <- as.data.frame(sigGenes2)
sigGenes2<-dplyr::distinct(.data = sigGenes2, GENE_ID, .keep_all = T)
names(sigGenes2)[2]<-"gene"
sigGenes2

png(filename = "~/jb/Redo_graphs/MA_plot.png",
    type = "cairo", width = 4500, height = 3400, res = 500)

x<-plotMA(res, main="MA-plot for Diet*Age Interaction Results", returnData = T)

ggplot(x, aes(mean, lfc))+
  geom_point(size=2, alpha=0.1)+
  # scale_colour_manual(values = "")+
  xlab(paste0("Mean of Normalised Counts")) +
  ylab(paste0("Log2-Fold Change")) + 
  # coord_fixed()+
  xlim(0, 1e4)+
  # ylim(-0.5, 0.5)+
  ylim(-2, 2)+
  theme_bw()+
  # ggtitle(expression(underline("MA-plot: Genes Level Counts (Interaction Term Comparison)")))+
  theme(plot.title = element_text(size = 30,face = "bold.italic", hjust = 0.5, vjust = 0), 
        plot.subtitle = element_text(size = 18, face = "italic"),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(face = "bold", size = 18, margin = margin(3)),
        axis.title.y = element_text(face = "bold", size = 18, vjust=1.5, margin = margin(3)),
        axis.text.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_text(size=16),
        legend.spacing = unit(x = 0, units = "pt"), 
        legend.box.spacing = margin(0.0001),
        # axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        strip.text = element_text(size = 16))

dev.off()

vsd<-vst(dds, blind=TRUE)

png(filename = "~/jb/Redo_graphs/RNA_Seq_PCA.png",
    type = "cairo", width = 4500, height = 3400, res = 500)

colours <- c("#56A8CBFF", "#DA291CFF", "#E69F00", "mediumpurple3")

pcaData <- DESeq2::plotPCA(vsd, intgroup="Age_Group", ntop=500000, returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData %<>% separate(Age_Group, into=c("Age", "Group"), sep = "-", extra='merge')
order<-c("Low-Fat Diet", "High-Fat Diet")
pcaData$Group<-factor(pcaData$Group, levels = order)

ggplot(pcaData, aes(PC1, PC2, color=Group, shape=Age))+
  geom_point(size=4)+
  scale_colour_manual(values = colours)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  # ggtitle(expression(underline("PCA Plot")))+
  theme(plot.title = element_text(size = 30,face = "bold.italic", hjust = 0.5, vjust = 0), 
        plot.subtitle = element_text(size = 18, face = "italic"),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(face = "bold", size = 18, margin = margin(3)),
        axis.title.y = element_text(face = "bold", size = 18, vjust=1.5, margin = margin(3)),
        axis.text.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_text(size=16),
        legend.spacing = unit(x = 0, units = "pt"), 
        legend.box.spacing = margin(0.0001),
        # axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        strip.text = element_text(size = 16))

dev.off()

library(ggrepel)

names(resOrdered)[which(names(resOrdered)=="gene")]<-"GENE_ID"

resOrdered<-left_join(resOrdered, mus_tx2id, by="GENE_ID")

resOrdered$sig_reg<-"Unchanged"
resOrdered$sig_reg[which(resOrdered$padj<0.05 & resOrdered$log2FoldChange<0)]<-"Significantly Downregulated"
resOrdered$sig_reg[which(resOrdered$padj<0.05 & resOrdered$log2FoldChange>0)]<-"Significantly Upregulated"
resOrdered$sig_reg<-factor(resOrdered$sig_reg, levels=c("Unchanged", "Significantly Downregulated", "Significantly Upregulated"))

res_1gene<-resOrdered %>% distinct(GENENAME, .keep_all = T)

# colours<-c("grey", "purple", "red")

png(filename = "~/jb/Redo_graphs/volcano_plot_deseq_interaction.png",
    type = "cairo", width = 4500, height = 3400, res = 500)

ggplot(res_1gene, 
       aes(log2FoldChange, -log10(padj), color=sig_reg))+
  geom_point()+
  scale_color_manual(values=c("#666666","#3366FF", "#FF3333"))+
  labs(color = "Significantly Changed Genes\n(Adjusted P-value < 0.05)")+
  # scale_color_manual(values=colours)+
  # # scale_fill_brewer(palette="Spectral")+
  xlim(-20,10)+
   ylim(0,5)+
  # ggtitle(label = expression(underline("Volcano Plot: Diet and Age Interaction")))+
  xlab(label=expression('log'[2]*'Fold Change'))+
  ylab(label=expression(-'log'[10]*'(Adjusted P-Value)'))+
  theme(plot.title = element_text(size = 24, hjust = 0.5, vjust=1),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))+
  geom_text_repel(
    data = res_1gene %>% 
      filter(padj<0.05) %>% distinct(GENENAME, .keep_all = T),
    aes(label = GENENAME), nudge_x = 0.01, nudge_y = 0.01,
    max.overlaps = 20,
    size = 3,segment.colour = "grey",
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.2, "lines")
  )+
  theme_bw()+
  theme(plot.title = element_text(size = 26,face = "bold.italic", hjust = 0.5, vjust = 0), 
        plot.subtitle = element_text(size = 16, face = "italic"),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(face = "bold", size = 16, margin = margin(3)),
        axis.title.y = element_text(face = "bold", size = 16, vjust=1.5, margin = margin(3)),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(size=12),
        legend.spacing = unit(x = 0, units = "pt"), 
        legend.box.spacing = margin(0.0001),
        # axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        strip.text = element_text(size = 16))
dev.off()

sig_genes<-res_1gene %>% filter(sig_reg=="Significant")
sig_genes %<>% select(GENENAME, log2FoldChange,padj)
sig_genes$regulation<-"Upregulated"
sig_genes$regulation[sig_genes$log2FoldChange<0]<-"Downregulated"
knitr::kable(sig_genes, 
             format = "latex", digits = 2)

# Kegg ----
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichR)

sigGenes2
sig_vector<-as.character(sigGenes2$GENENAME)

anno <- AnnotationDbi::select(org.Mm.eg.db, 
                              keys=sig_vector, 
                              columns=c("GENENAME", "ENSEMBL", "ENTREZID"),
                              keytype="SYMBOL")
anSig <- as.data.frame(subset(anno, SYMBOL %in% sig_vector))
gene_list<-as.character(na.omit(anSig$ENTREZID))
kk <- enrichKEGG(gene = gene_list, organism = 'mmu', 
                 pvalueCutoff = 0.05)


