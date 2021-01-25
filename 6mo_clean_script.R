### done with Deseq
###### 6 months sig
library(EnsDb.Hsapiens.v86)
library(Biostrings)
library(stringr)
library(dplyr)
library(tximport)
library(plyr)
library(DESeq2)
library(plyr)
library(AnnotationDbi)
library(magrittr)
library(ggplot2)
library(forcats)

setwd("~/jb/6mo_RNASeq_analysis/salmon_folder")
dir<-getwd()
##files<-list.files(dir) lists all folders/files within subdirectory
###files<-list.dirs(dir) #lists all subfolders within directory
files<-list.files(dir, recursive = TRUE)
files<-files[grep("quant.sf", files)]
filenames<-sub("(^.*)/.*", "\\1", files)
files<-file.path(dir, files)
names(files)<-filenames
all(file.exists(files))
print("############### got all folders, now importing")

# mus_tx2id<-read.delim(file = "/Volumes/SEAGATE/rnaseq/salmon_folder/mus_gencode_txid_to_genes.txt", sep = " ", header = FALSE)
mus_tx2id<-read.delim(file = "~/jb/Thesis_RNASeq_Analysis/mus_gencode_txid_to_genes.txt", sep = " ", header = T)

names(mus_tx2id)<-c("TX_ID", "GENE_ID", "GENENAME")

txi.salmon_mus_6mo<-tximport(files, type = "salmon", tx2gene = mus_tx2id, ignoreTxVersion = T,countsFromAbundance = "lengthScaledTPM")


###708244, 708246, 71770 are chow, 5RR, 7NM, 15RR, 20LL, 5RR are LFD 
condition<-c(1:length(files))
for (i in 1:length(files)){
  if(grepl("^7.*", names(files[i]))){
    condition[i]<-"Control"
  } else if (grepl("5RR|7NM|15RR|20LL", names(files[i]))){
    condition[i]<-"Control"
  }
  else{
    condition[i]<-"HFD"
  }
}

sampleMetadata <- data.frame(
  samples = names(files),
  condition = condition
)
dds_6mo <- DESeqDataSetFromTximport(txi.salmon_mus_6mo,
                                    colData = sampleMetadata,
                                    design = ~ condition)
keep <- rowSums(counts(dds_6mo)) >= 10
dds_6mo <- dds_6mo[keep,]

dds_6mo <- DESeq(dds_6mo)
dds_6mo$condition <- relevel(dds_6mo$condition, ref = "Control")

res_6mo <- results(dds_6mo, contrast=c("condition","HFD","Control"))
res_6mo$gene <- row.names(res_6mo)
resOrdered_6mo <- res_6mo[order(res_6mo$padj),]
resOrdered_6mo$gene <- row.names(resOrdered_6mo)
resOrdered_6mo <- as.data.frame(resOrdered_6mo)

#is this a good pvalue setting? No, it's unadjusted
sigGenes <- rownames(subset(res_6mo, pvalue < 0.05 & abs(log2FoldChange > 1)))

sigGenes<-mus_tx2id[which(mus_tx2id$GENE_ID %in% sigGenes),]

unique(sigGenes$GENENAME)->deseq_sig_genes
###compare these with Edge_R
names(resOrdered_6mo)[which(names(resOrdered_6mo)=="gene")]<-"GENE_ID"

resOrdered_6mo<-left_join(resOrdered_6mo, mus_tx2id, by="GENE_ID")


#sigGenesOrdered <- subset(resOrdered_6mo, padj < 0.05 & abs(log2FoldChange > 1))

sigGenesOrdered_6mo <- subset(resOrdered_6mo, padj < 0.05 & abs(log2FoldChange) > 1)

sigGenesOrdered2_6mo <- subset(resOrdered_6mo, padj < 0.1)

# sigGenesOrdered_6mo <- subset(resOrdered_6mo, padj < 0.05)

# unique(sigGenesOrdered$GENENAME)

deseq_sig_genes<-unique(sigGenesOrdered_6mo$GENENAME)

sig_6mo<-unique(sigGenesOrdered2_6mo$GENENAME)

##############
########################

library(org.Mm.eg.db)
library(clusterProfiler)
library(DOSE)
library(genefilter)
library(geneplotter)
library(topGO)
library(limma)
library(enrichR)

# sig_vector<-as.character(unique(sigGenesOrdered2_6mo$GENENAME))
sig_vector<-as.character(sig_6mo)


anno <- AnnotationDbi::select(org.Mm.eg.db, 
                              keys=sig_vector, 
                              columns=c("GENENAME", "ENSEMBL", "ENTREZID"),
                              keytype="SYMBOL")

anSig <- as.data.frame(subset(anno, SYMBOL %in% sig_vector))

anSig<-na.omit(anSig)
anSig_6mo<-anSig
head(anSig)

#######################
######################
# 
# library(ggrepel)
# 
# resOrdered_6mo$sig_reg<-"none"
# resOrdered_6mo$sig_reg[which(resOrdered_6mo$log2FoldChange>=1 & resOrdered_6mo$padj<0.05)]<-"up"
# resOrdered_6mo$sig_reg[which(resOrdered_6mo$log2FoldChange<=-1 & resOrdered_6mo$padj<0.05)]<-"down"
# resOrdered_6mo$sig_reg<-factor(resOrdered_6mo$sig_reg, levels=c("none", "down", "up"))
# 
# res_1gene<-resOrdered_6mo %>% distinct(GENENAME, .keep_all = T)
# 
# # colours<-c("grey", "purple", "red")
# 
# png(filename = "/Users/Jonathan/jb/Thesis_RNASeq_Analysis/6mo_volcano_plot_deseq.png",
#     type = "cairo", width = 2500, height = 1500, res = 300)
# 
# ggplot(res_1gene, aes(log2FoldChange, -log10(padj), color=sig_reg))+
#   geom_point()+
#   scale_color_manual(values=c("#666666","#3366FF", "#FF3333"))+
#   labs(color = "Significant up or\n down regulation?")+
#   # scale_color_manual(values=colours)+
#   # # scale_fill_brewer(palette="Spectral")+
#   xlim(-10,10)+
#   ylim(0,5)+
#   ggtitle(label = expression(underline("6 Month Volcano Plot: DESeq2 Analysis")))+
#   xlab(label=expression('log'[2]*'Fold Change'))+
#   ylab(label=expression(-'log'[10]*'(Adjusted P-Value)'))+
#   theme(plot.title = element_text(size = 24, hjust = 0.5, vjust=1), 
#         axis.title.x = element_text(size=14),
#         axis.title.y = element_text(size=14))+
#   geom_text_repel(
#     data = filter(res_1gene, padj<0.05 & abs(log2FoldChange)>=1),
#     aes(label = GENENAME),
#     size = 4,segment.colour = "grey",
#     box.padding = unit(0.3, "lines"),
#     point.padding = unit(0.2, "lines")
#   )
# 
# dev.off()

#######################
###################### MA plot
# 
# DESeq2::plotMA(res_6mo)

###################
#################
# 
# kk <- enrichKEGG(gene = unique(as.vector(anSig$ENTREZID)), organism = 'mmu')
# 
# head(kk, n=10)
# 
# ##################
# #################
# 
# k <- kegga(unique(anSig$ENTREZID), species.KEGG = "mmu")
# 
# topKEGG(results = k)
# 

# x<-listEnrichrDbs()
# y<-enrichR::enrichr(unique(anSig$SYMBOL), databases = x$libraryName)

###########################
#######################
library(clusterProfiler)
library(rWikiPathways)
kegg_enrich <- enrichKEGG(gene = as.character(unique(na.omit(anSig$ENTREZID))), pAdjustMethod = "none", pvalueCutoff = 1, organism = 'mouse')
head(kegg_enrich)

y <- setReadable(kegg_enrich, 'org.Mm.eg.db', keyType="ENTREZID")
kegg_table<-head(summary(y),10)[,c(2,8)]
kegg_table<-tibble::as_tibble(kegg_table, rownames = NULL)
names(kegg_table)<-c("KEGG_Pathway", "Significantly Upregulated_Genes")

z<-enrichR::enrichr(unique(anSig$SYMBOL), databases = "KEGG_2019_Mouse")
head(z, 10)

y<-tolower(z$KEGG_2019_Mouse$Genes)
CapStr <- function(y) {
  c <- strsplit(y, ";")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=";")
}
x<-sapply(X = y, FUN = CapStr)
x<-unname(x)
x<-gsub(x = x, pattern = "H2-aa", replacement = "H2-Aa")
z$KEGG_2019_Mouse$Genes<-x

kegg_table<-tibble::as_tibble(z[[1]][,c(1,2,9)], rownames = NULL)
kegg_table<-kegg_table[1:10,]
names(kegg_table)<-c("KEGG_Pathway", "# of Genes Found In Pathway", "Gene Names")

knitr::kable(kegg_table[,c(1,2,3)], format = "simple", booktabs = TRUE, caption="6 Month KEGG Pathway Results")

knitr::kable(kegg_table[,c(1,3)], format = "latex", booktabs = TRUE, caption="6 Month KEGG Pathway Results")

# 

# 
# enrichplot::dotplot(kegg_enrich)
#   
#   
#   barplot(kegg_enrich, 
#         drop = TRUE, color = "pvalue",
#         showCategory = 15, 
#         title = "KEGG Enrichment Pathways: 6 Month DESeq2 Analysis",
#         font.size = 8)


z<-enrichR::enrichr(unique(anSig$SYMBOL), databases = "KEGG_2019_Mouse")

kegg<-head(as_tibble(z, validate = T),3)
kegg<-kegg$KEGG_2019_Mouse[,c(1,2,3,9)]

png(filename = "/Users/Jonathan/jb/Thesis_RNASeq_Analysis/6mo_kegg_graph.png",
    type = "cairo", width = 5000, height = 3000, res = 350)

# kegg %>% mutate(name = fct_reorder(Term, desc(P.value))) %>%
#   ggplot(aes(x=name, y=-log10(P.value), size=factor(Overlap))) +
#   geom_point(aes(color=P.value)) +
#   scale_color_continuous(type = "viridis", trans = 'reverse')+
#   coord_flip()+
#   ylab(label=expression(-'log'[10]*'(Uncorrected P-Value)'))+
#   xlab(label="KEGG Pathway")+
#   ggtitle(label = expression(underline("KEGG Enrichment Pathways: 6 Months")))+
#   labs(fill = "P-value (uncorrected)", size="Overlap (Genes Found  /\n Total Genes in Pathway)")+
#   theme(legend.title = element_text(size=12), axis.text.y = element_text(size=12), plot.title = element_text(size = 26),
#         axis.title.x = element_text(size=16),
#               axis.title.y = element_text(size=16),
#         legend.position = "top",
#         legend.direction = "horizontal")

kegg %>% mutate(name = fct_reorder(Term, desc(P.value))) %>%
  ggplot(aes(x=name, y=-log10(P.value))) +
  geom_point(aes(color=P.value, shape=factor(Overlap)),size=6) +
  scale_color_continuous(type = "viridis", trans = 'reverse', breaks=c(max(kegg$P.value),min(kegg$P.value)))+
  coord_flip()+
  ylab(label=expression(-'log'[10]*'(Uncorrected P-Value)'))+
  xlab(label="KEGG Pathway")+
  ggtitle(label = expression(underline("KEGG Enrichment Pathways: 18 Months")))+
  labs(fill = "P-value (uncorrected)", shape="Overlap (Genes Found  /\n Total Genes in Pathway)")+
  theme(legend.title = element_text(size=14), axis.text.y = element_text(size=12),
        plot.title = element_text(size = 26),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.position = "top",
        legend.direction = "horizontal")

dev.off()

############
############

BP<-head(as.data.frame(enrichR::enrichr(unique(anSig$SYMBOL), databases = "GO_Biological_Process_2018")), 5)
CC<-head(as.data.frame(enrichR::enrichr(unique(anSig$SYMBOL), databases = "GO_Cellular_Component_2018")), 5)
MF<-head(as.data.frame(enrichR::enrichr(unique(anSig$SYMBOL), databases = "GO_Molecular_Function_2018")), 5)

names(BP)<-gsub(names(BP), pattern = "^.*\\.(.*)$", replacement="\\1")
BP2<-BP[,c(1,3,9)]
BP2$Term<-gsub(BP2$Term, pattern = "^(.*) \\(.*\\)$", replacement="\\1")
BP2$ONT<-"BP"

names(CC)<-gsub(names(CC), pattern = "^.*\\.(.*)$", replacement="\\1")
CC2<-CC[,c(1,3,9)]
CC2$Term<-gsub(CC2$Term, pattern = "^(.*) \\(.*\\)$", replacement="\\1")
CC2$ONT<-"CC"

names(MF)<-gsub(names(MF), pattern = "^.*\\.(.*)$", replacement="\\1")
MF2<-MF[,c(1,3,9)]
MF2$Term<-gsub(MF2$Term, pattern = "^(.*) \\(.*\\)$", replacement="\\1")
MF2$ONT<-"MF"

xy<-rbind(BP2, CC2, MF2)
xy$Description<-factor(xy$Description)

xy<-head(xy %>% mutate(name = fct_reorder(Term, desc(value))), 3)

png(filename = "/Users/Jonathan/jb/Thesis_RNASeq_Analysis/6mo_GO_graph.png",
    type = "cairo", width = 5000, height = 3000, res = 350)

xy %>% mutate(name = fct_reorder(Term, desc(value))) %>%
  ggplot(aes(x=name, y=-log10(value), fill=ONT)) +
  geom_bar(stat="identity", alpha=.6, width=.4) +
  coord_flip()+
  ylab(label=expression(-'log'[10]*' P-Value (uncorrected)'))+
  xlab(label="GO Pathway")+
  scale_fill_discrete(name = "Gene Ontology", labels = c("Biological Process",
                                                         "Cellular Component",
                                                         "Molecular Function"))+
  ggtitle(label = expression(underline("GO Pathway Analysis: 6 Months")))+
  theme(legend.title = element_text(size=12), axis.text.y = element_text(size=12), 
        plot.title = element_text(size = 26),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.position = "top",
        legend.direction = "horizontal")

dev.off()
########################
################################
########################################

# gsea_list<-sigGenesOrdered %>% distinct(log2FoldChange, .keep_all = T) %>% select(GENENAME, log2FoldChange)
# names(gsea_list)[1]<-"SYMBOL"
# gsea_list<-left_join(gsea_list, anno, by='SYMBOL')
# gsea_list %<>% distinct(ENTREZID, .keep_all = T) %>% select(ENTREZID, log2FoldChange)
# gsea_list<-na.omit(gsea_list)
# e_id<-gsea_list$ENTREZID
# gsea_list<-as.vector(gsea_list$log2FoldChange)
# names(gsea_list)<-e_id
# gsea_list<-sort(gsea_list, decreasing = T)

# ewp2 <- GSEA(gsea_list, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=T, pvalueCutoff = 0.9)
# 
# go_enrich <- enrichGO(gene =anSig$ENTREZID,
#                       OrgDb = 'org.Mm.eg.db',
#                       readable = T,
#                       ont = "ALL",
#                       pvalueCutoff = 0.05)
# 
# head(go_enrich)
# 
# 
# y<-go_enrich@result[order(go_enrich@result$qvalue),]
# 
# BP<-y %>% filter(ONTOLOGY=="BP")
# CC<-y %>% filter(ONTOLOGY=="CC")
# MF<-y %>% filter(ONTOLOGY=="MF")
# 
# xy<-rbind(head(BP,4), head(CC,4), head(MF,4))
# xy$Description<-factor(xy$Description)
# 
# xy %>% mutate(name = fct_reorder(Description, desc(pvalue))) %>%
#   ggplot(aes(x=name, y=-log10(pvalue), fill=ONTOLOGY)) +
#   geom_bar(stat="identity", alpha=.6, width=.4) +
#   coord_flip()
# 

library(pathview)
library(KEGGREST)
org <- keggList("mmu", database = "pathway")
head(org)
y<-grep(x = org, pattern = "coagulation")
org[y]
setwd(dir = "~/jb/Thesis_RNASeq_Analysis/")
pathview(gene.data = anSig$ENTREZID, pathway.id = "04610", species = "mouse", , plot.col.key=F)

# Serpind1 is the gene for heparin cofactor 2- ie HCII

# To add 6mo and 18mo anSigs together to get F5 added to complement cascade picture

comb_sig<-rbind(anSig_6mo, anSig_18mo)

pathview(gene.data = comb_sig$ENTREZID, pathway.id = "04610", species = "mouse", , plot.col.key=F, res=400)

