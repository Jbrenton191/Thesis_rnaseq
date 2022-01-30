# library(EnsDb.Hsapiens.v86)
# library(Biostrings)
# library(stringr)
library(tximport)
library(DESeq2)
# # library(plyr)
# # library(AnnotationDbi)
library(ggplot2)
# # library(forcats)
 library(tidyverse)
library(magrittr)
# library(dplyr)
library(DEGreport)
# # changed 11L to 11LL and 2R and 2RR in directory titles
# 
cat("Start")
# # setwd("~/jb/6mo_RNASeq_analysis/salmon_folder")
# setwd("/Volumes/SEAGATE/Post_viva_rnaseq/salmon")
# 
# dir<-getwd()
# ##files<-list.files(dir) lists all folders/files within subdirectory
# ###files<-list.dirs(dir) #lists all subfolders within directory
# files<-list.files(dir, recursive = TRUE)
# files<-files[grep("quant.sf", files)]
# filenames<-sub("(^.*)/.*", "\\1", files)
# files<-file.path(dir, files)
# names(files)<-filenames
# all(file.exists(files))
# print("############### got all folders, now importing")
# 
# # mus_tx2id<-read.delim(file = "/Volumes/SEAGATE/rnaseq/salmon_folder/mus_gencode_txid_to_genes.txt", sep = " ", header = FALSE)
# # mus_tx2id<-read.delim(file = "/Volumes/SEAGATE/LIDO_MAC_02_2021/jb/Thesis_RNASeq_Analysis/mus_gencode_txid_to_genes.txt", sep = " ", header = T)
# 
# names(mus_tx2id)<-c("TX_ID", "GENE_ID", "GENENAME")
# 
# #removing Chow samples:
# files<-files[-c(1:3)]
# # 708244, 708246, 71770 are chow, 5RR, 7NM, 15RR, 20LL, 5RR are LFD at 6mo
# ### JBE3RR,JBE5LL,JBE9L,JBE1NM are LFD, and JBE11L, JBE22NM, JBE17LL, JBE27L, at 18mo
# 
# meta_info<-read_csv(file = "../../LIDO_MAC_02_2021/jb/weight_and_blood_data-analysis/All_timepoints_full_mouse_experiments_list.csv", col_names = T) %>% dplyr::select(Animal, Group, Age)
# #
# names(files) %in% meta_info$Animal
# 
# meta_info<-meta_info[meta_info$Animal %in% names(files),]
# 
# condition<-c(1:length(files))
# for (i in 1:length(files)){
#   if(grepl("JBE3RR|5LL|9L|1NM", names(files[i]))){
#     condition[i]<-"LFD"
#   }else if(grepl("JBE11L|22NM|17LL|27L|5RR|7NM|15RR|20LL", names(files[i]))){
#     condition[i]<-"HFD"
#   } else{
#     condition[i]<-"HFD"
#   }
# }


# If SEAGATE HD not plugged in uncomment below and comment out the txi and metadata objects
txi.salmon <- readRDS("~/jb/post_viva_work/txi_salmon.rds")
sampleMetadata<-readRDS("~/jb/post_viva_work/sampleMetadata_deseq.rds")
mus_tx2id<-readRDS("~/jb/post_viva_work/mus_tx2id_deseq.rds")


# 
txi.salmon<-tximport(files, type = "salmon",
                         tx2gene = mus_tx2id,
                         ignoreTxVersion = T)
# 
# sampleMetadata <- data.frame(
#   samples = names(files),
#   condition=factor(condition))
# 
# reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9 10 Error in tximport(files, type = "salmon",
# tx2gene = mus_tx2id, ignoreTxVersion = T,  : 
#                                          all(txId == raw[[txIdCol]]) is not TRUE


# Age<-c(rep(6,9), rep(18,8))
# 
# sampleMetadata <- meta_info
# sampleMetadata$Breeder<-"BSU"
# sampleMetadata$Breeder[grep(x = sampleMetadata$Animal, pattern = "JBE1[5-9]|JBE20|JBE3[8-9]|JBE4[0-3]|JBE3[0|6-7]|JBE2[8-9]")]<-"ENVIGO"
# 
# 
# sampleMetadata$Age_Group<- factor(paste0(sampleMetadata$Age, "-", sampleMetadata$Group))
# 

# cat("End")


####################
#########################3

dds <- DESeqDataSetFromTximport(txi.salmon,
                                colData=sampleMetadata,
                                # design=~Breeder+Age+Group)
                                design=~Breeder+Group+Age+Group:Age)
                                # design=~Breeder+Age_Group)



                                
dds$Group<-relevel(dds$Group, ref = "Low-Fat Diet")
dds$Age<-relevel(dds$Age, ref = "6 Months")
# dds$Age_Group<-factor(x = dds$Age_Group, levels = c("6 Months-Low-Fat Diet", "6 Months-High-Fat Diet", "18 Months-Low-Fat Diet", "18 Months-High-Fat Diet"))
# dds$Age_Group<-relevel(x = dds$Age_Group, ref = "6 Months-Low-Fat Diet")

filter_count <- DEGreport::degFilter(counts = counts(dds),
                                     metadata = as.data.frame(colData(dds)),
                                     group = "Age_Group",
                                     min = 1, # All samples in group must have more than expr > 0
                                     minreads = 0)

# filter_count <- DEGreport::degFilter(counts = counts(dds), 
#                                      metadata = as.data.frame(colData(dds)),
#                                      group = "Group",
#                                      min = 1, # All samples in group must have more than expr > 0
#                                      minreads = 0)

dds <- dds[rownames(filter_count),]

     


                                # design=~Age+Group+Age:Group)
# keep <- rowSums(counts(dds))>=10
# dds <- dds[keep,]
# 
dds <- DESeq(dds)

res<-results(dds)

png(filename = "~/jb/Redo_graphs/MA_plot.png",
    type = "cairo", width = 4500, height = 3000, res = 300)

x<-plotMA(res, main="MA-plot for Diet*Age Interaction Results", returnData = T)

ggplot(x, aes(mean, lfc))+
  geom_point(size=2, alpha=0.2)+
  # scale_colour_manual(values = "")+
  xlab(paste0("Mean of Normalised Counts")) +
  ylab(paste0("Log2-Fold Change")) + 
  # coord_fixed()+
  # xlim(1, 10000)+
  # ylim(-1, 1)+
  theme_bw()+
  ggtitle(expression(underline("MA-plot: Diet*Age Interaction Results")))+
  theme(plot.title = element_text(size = 28, hjust = 0.5), 
        plot.subtitle = element_text(face = "italic", size = 12),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18, vjust=1.5),
        legend.text =element_text(size = 13),
        legend.title =element_text(size=14))


dev.off()

# dds <- DESeq(dds, test = "LRT", reduced = ~Age+Group)

# res <- DESeq2::results(dds, name = "Age_Group_18.Months.Low.Fat.Diet_vs_6.Months.Low.Fat.Diet")

# res <- DESeq2::results(dds, name ="Group_High.Fat.Diet_vs_Low.Fat.Diet")
# res <- DESeq2::results(dds, name ="Age_18.Months_vs_6.Months")
# res <- DESeq2::results(dds, name= "Group1")
# res <- DESeq2::results(dds, name = "Age1")


# res<-results(dds, contrast = c("Group", "Low-Fat Diet", "High-Fat Diet"))
# res<-results(dds, contrast = c("Age", "18 Months", "6 Months"))
# res<-results(dds, list ( c("Group_Low.Fat.Diet_vs_High.Fat.Diet", "GroupLow.Fat.Diet.Age6.Months")))
# res<-results(dds, name="GroupLow.Fat.Diet.Age6.Months")
# res<-results(dds, name="Age_Group1")

#to find contrasts:
# resultsNames(dds)

# res <- DESeq2::results(dds, name = "Age_18.Months_vs_6.Months")
res <- DESeq2::results(dds, name = "Group_High.Fat.Diet_vs_Low.Fat.Diet")
# res <- DESeq2::results(dds, name = "GroupHigh.Fat.Diet.Age18.Months")  


res$gene <- row.names(res)
resOrdered<-res[order(res$padj),]
resOrdered$gene <- row.names(resOrdered)
resOrdered <- as.data.frame(resOrdered)

#is this a good pvalue setting? No, it's unadjusted
# sigGenes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))
sigGenes <- subset(res, padj < 0.05)
sigGenes <- as.data.frame(sigGenes)
sigGenes2<-mus_tx2id[which(mus_tx2id$GENE_ID %in% sigGenes$gene),]
sigGenes2 <- as.data.frame(sigGenes2)
sigGenes2<-dplyr::distinct(.data = sigGenes2, GENE_ID, .keep_all = T)
names(sigGenes2)[2]<-"gene"


sigGenes_group<-left_join(sigGenes, sigGenes2, by="gene")

# sigGenes_group$Effect<-"Diet"
sigGenes_group



x<-plotCounts(dds, gene = "ENSMUSG00000026696", 
              intgroup = c("Group","Age"), returnData = TRUE)


x<-plotCounts(dds, gene = "ENSMUSG00000020538", 
              intgroup = c("Group","Age"), returnData = TRUE)

ggplot(x,aes(x = Age, y = count, color = Group, group = Group)) + 
  geom_point() + 
  stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()

# do again for age ----


# res<-results(dds, contrast = c("Group", "Low-Fat Diet", "High-Fat Diet"))
# res<-results(dds, contrast = c("Age", "18 Months", "6 Months"))
# res<-results(dds, list ( c("Group_Low.Fat.Diet_vs_High.Fat.Diet", "GroupLow.Fat.Diet.Age6.Months")))
# res<-results(dds, name="GroupLow.Fat.Diet.Age6.Months")
# res<-results(dds, name="Age_Group1")


res$gene <- row.names(res)
resOrdered<-res[order(res$padj),]
resOrdered$gene <- row.names(resOrdered)
resOrdered <- as.data.frame(resOrdered)

#is this a good pvalue setting? No, it's unadjusted
# sigGenes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))
sigGenes <- subset(res, padj < 0.05)
sigGenes <- as.data.frame(sigGenes)
sigGenes2<-mus_tx2id[which(mus_tx2id$GENE_ID %in% sigGenes$gene),]
sigGenes2 <- as.data.frame(sigGenes2)
sigGenes2<-dplyr::distinct(.data = sigGenes2, GENE_ID, .keep_all = T)
names(sigGenes2)[2]<-"gene"

sigGenes_age<-left_join(sigGenes, sigGenes2, by="gene")
sigGenes_age$Effect<-"Age"

# Combining the two ----

sigGenes_table<-rbind(sigGenes_age, sigGenes_group)
sigGenes_table %>% select(GENENAME, Effect, padj)

knitr::kable(sigGenes_table %>% select(GENENAME, Effect, padj), format = "latex")
###########

vsd<-vst(dds, blind=TRUE)

png(filename = "~/jb/Redo_graphs/RNA_Seq_PCA.png",
    type = "cairo", width = 4500, height = 3000, res = 300)

colours <- c("#56A8CBFF", "#DA291CFF", "#E69F00", "mediumpurple3")

pcaData <- DESeq2::plotPCA(vsd, intgroup="Age_Group", ntop=500000, returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData %<>% separate(Age_Group, into=c("Age", "Group"), sep = "-", extra='merge')
order<-c("Low-Fat Diet", "High-Fat Diet")
pcaData$Group<-factor(pcaData$Group, levels = order)

ggplot(pcaData, aes(PC1, PC2, color=Group, shape=Age))+
  geom_point(size=3)+
  scale_colour_manual(values = colours)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  ggtitle(expression(underline("PCA Plot")))+
  theme(plot.title = element_text(size = 28, hjust = 0.5), 
        plot.subtitle = element_text(face = "italic", size = 12),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18, vjust=1.5),
        legend.text =element_text(size = 13),
        legend.title =element_text(size=14))

dev.off()


library(ggrepel)

names(resOrdered)[which(names(resOrdered)=="gene")]<-"GENE_ID"

resOrdered<-left_join(resOrdered, mus_tx2id, by="GENE_ID")

resOrdered$sig_reg<-"Unchanged"
resOrdered$sig_reg[which(resOrdered$padj<0.05)]<-"Significant"
resOrdered$sig_reg<-factor(resOrdered$sig_reg, levels=c("Unchanged", "Significant"))

# res_1gene<-resOrdered_18mo %>% distinct(GENENAME, .keep_all = T)

# colours<-c("grey", "purple", "red")

png(filename = "~/jb/Redo_graphs/volcano_plot_deseq.png",
    type = "cairo", width = 4000, height = 2800, res = 350)

ggplot(resOrdered, 
       aes(log2FoldChange, -log10(padj), color=sig_reg))+
  geom_point()+
  scale_color_manual(values=c("#666666","#3366FF", "#FF3333"))+
  labs(color = "Significantly Changed Genes\n(Adjusted P-value < 0.05)")+
  # scale_color_manual(values=colours)+
  # # scale_fill_brewer(palette="Spectral")+
  xlim(-10,10)+
  ylim(0,5)+
  # ggtitle(label = expression(underline("Volcano Plot: 6 Month LFD and 18 Month HFD")))+
  xlab(label=expression('log'[2]*'Fold Change'))+
  ylab(label=expression(-'log'[10]*'(Adjusted P-Value)'))+
  theme(plot.title = element_text(size = 24, hjust = 0.5, vjust=1),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))+
  geom_text_repel(
    data = resOrdered %>% filter(padj<0.05) %>% distinct(GENENAME, .keep_all = T),
    aes(label = GENENAME), nudge_x = 0.4, nudge_y = 0.4,
    max.overlaps = 20,
    size = 4,segment.colour = "grey",
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.2, "lines")
    )+
  theme_bw()+
  theme(plot.title = element_text(size = 28,face = "bold.italic", hjust = 0.5, vjust = 0), 
        plot.subtitle = element_text(size = 14, face = "italic"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 16, margin = margin(3)),
        axis.title.y = element_text(face = "bold", size = 16, vjust=1.5, margin = margin(3)),
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_text(size=12),
        legend.spacing = unit(x = 0, units = "pt"), 
        legend.box.spacing = margin(0.0001))

dev.off()



# saveRDS(txi.salmon, "~/jb/post_viva_work/txi_salmon.rds")
# saveRDS(sampleMetadata , "~/jb/post_viva_work/sampleMetadata_deseq.rds")
# saveRDS(mus_tx2id, "~/jb/post_viva_work/mus_tx2id_deseq.rds")
# 



















# 
# 
# 
# ggplot(resOrdered, aes(log2FoldChange, -log10(padj)))+
#   geom_point()+
#   scale_color_manual(values=c("#666666","#3366FF", "#FF3333"))+
#   labs(color = "Significant up or\n down regulation?")+
#   # scale_color_manual(values=colours)+
#   # # scale_fill_brewer(palette="Spectral")+
#    xlim(-5,5)+
#   ylim(0,2)+
#   ggtitle(label = expression(underline("Volcano Plot: DESeq2 Analysis")))+
#   xlab(label=expression('log'[2]*'Fold Change'))+
#   ylab(label=expression(-'log'[10]*'(Adjusted P-Value)'))+
#   theme(plot.title = element_text(size = 24, hjust = 0.5, vjust=1),
#         axis.title.x = element_text(size=14),
#         axis.title.y = element_text(size=14))
# 








# End of original analysis












#####

filter_count <- DEGreport::degFilter(counts = counts(dds),
                     metadata = as.data.frame(colData(dds)),
                     group = "Group",
                     min = 1, # All samples in group must have more than expr > 0
                     minreads = 0)

filter_dds <- dds[rownames(filter_count),]

res<-results(filter_dds, name = "GroupLow.Fat.Diet.Age6.Months")
# x<-results(filter_dds, contrast=c("Group","Low-Fat Diet", "High-Fat Diet"))
# res<-results(filter_dds, contrast=c("Group","Low-Fat Diet", "High-Fat Diet"))


res$gene <- row.names(res)
resOrdered<-res[order(res$padj),]
resOrdered$gene <- row.names(resOrdered)
resOrdered <- as.data.frame(resOrdered)

sigGenes <- rownames(subset(resOrdered, padj < 0.05))

mus_tx2id[mus_tx2id$GENE_ID %in% sigGenes,]

# ###############
# dds_LRT <- DESeq(filter_dds, test = "LRT", reduced = ~Age)
# 
# res<-results(dds_LRT)
# 
# res$gene <- row.names(res)
# resOrdered<-res[order(res$padj),]
# resOrdered$gene <- row.names(resOrdered)
# resOrdered <- as.data.frame(resOrdered)
# 
# sigGenes <- rownames(subset(resOrdered, padj < 0.05))
# 

############------------------------------------------------------------------------
#  Reanalysis to look at lFD 6 months vs HFD 18 Months

x<-meta_info %>% filter(Group=="High-Fat Diet" & Age=='18 Months')
y<-meta_info %>% filter(Group=="Low-Fat Diet" & Age=='6 Months')

meta_info<-rbind(x,y)

files<-files[names(files) %in% meta_info$Animal]

txi.salmon<-tximport(files, type = "salmon",
                     tx2gene = mus_tx2id,
                     ignoreTxVersion = T,
                     countsFromAbundance = "lengthScaledTPM")

sampleMetadata <- meta_info

dds <- DESeqDataSetFromTximport(txi.salmon,
                                colData=sampleMetadata,
                                design=~Group)
dds <- DESeq(dds)



filter_count <- DEGreport::degFilter(counts = counts(dds),
                                     metadata = as.data.frame(colData(dds)),
                                     group = "Group",
                                     min = 1, # All samples in group must have more than expr > 0
                                     minreads = 0)

filter_dds <- dds[rownames(filter_count),]

res<-results(filter_dds)
# x<-results(filter_dds, contrast=c("Group","Low-Fat Diet", "High-Fat Diet"))

res$gene <- row.names(res)
resOrdered<-res[order(res$padj),]
resOrdered$gene <- row.names(resOrdered)
resOrdered <- as.data.frame(resOrdered)

sigGenes <- rownames(subset(resOrdered, padj < 0.05))

mus_tx2id[mus_tx2id$GENE_ID %in% sigGenes,]

---
  vsd<-vst(dds, blind=TRUE)
DESeq2::plotPCA(vsd, intgroup=c("Age", "Group"), returnData=T)

pca_vsd <- prcomp(t(assay(vsd)))


x<-tibble(PC = c(1:17),
sdev = pca_vsd$sdev) %>%
dplyr::mutate(d = sdev^2,
pev = d/sum(d),
cev = cumsum(d)/sum(d))
#missing step

y<-pca_vsd$x[,c(1,2)]
zz<-rownames(y)
y<-as.tibble(y)
y$Animal2<-zz
col_d<-as.data.frame(vsd@colData)
col_d$Animal2<-rownames(col_d)
bb<-left_join(x = y, y = col_d, by="Animal2")

ggplot(bb,mapping = aes(PC1, PC2, color=Group, shape=Age))+
  geom_point(size=3)


library(ggpubr)
library(pcaExplorer)
ggarrange(pcascree(pca_vsd, type = "pev"),
pcascree(pca_vsd, type = "cev"),
nrow = 2,
labels = c("a", "b"))

# first one better
# PC_corr_PCs <- correlatePCs(pcaobj = pca_vsd, coldata = colData(vsd)[,2:4], pcs = 1:17)
PC_corr_PCs <- correlatePCs(pcaobj = pca_vsd, coldata = colData(vsd), pcs = 1:17)

PC_corr_FDR_PCs <- PC_corr_PCs$pvalue %>% 
  as_tibble() %>% 
  dplyr::mutate(PC = str_c("PC_", c(1:17))) %>% 
  dplyr::select(-sample_id) %>% 
  tidyr::gather(key = covariate, value = p, -PC) %>% 
  dplyr::mutate(FDR = p.adjust(p, method = "fdr"))

# first PCA correlates with sizefactor
# sampleMetadata<-as.tibble(colData(dds)[,2:4])
sampleMetadata<-as.tibble(colData(dds)[,2:5])

dds <- DESeqDataSetFromTximport(txi.salmon,
                                colData=sampleMetadata,
                                # design=~Group+Age)
                                # design=~Age+Group+Age:Group)
                                # design=~Group+Age+sizeFactor)
                                design=~Group+Age+Breeder+Group:Age)

# dds <- DESeq(dds,  test="LRT", reduced = ~ Group+Age)
dds <- DESeq(dds,  test="LRT", reduced = ~ Age + Group + Breeder)

res <- DESeq2::results(dds)

res$gene <- row.names(res)
resOrdered<-res[order(res$padj),]
resOrdered$gene <- row.names(resOrdered)
resOrdered <- as.data.frame(resOrdered)

#is this a good pvalue setting? No, it's unadjusted
sigGenes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange > 1)))
sigGenes <- rownames(subset(res, padj < 0.05))


sigGenes<-mus_tx2id[which(mus_tx2id$GENE_ID %in% sigGenes),]

sigGenes<-unique(sigGenes$GENENAME)
sigGenes







