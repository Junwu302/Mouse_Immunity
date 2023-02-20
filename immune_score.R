# load FPKM
FPKM2TPM = function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
df_tpm = read.delim2("matrix_gene.fpkm.annot.txt",header = T, sep='\t',stringsAsFactors = F)
df_tpm = df_tpm[,c(1:14,19,20)]
for(i in 2:14){
  df_tpm[,i] = as.numeric(df_tpm[,i])
}
ind = apply(df_tpm, 1, anyNA)
df_tpm = df_tpm[!ind,]
for(i in 2:14){
  df_tpm[,i] = FPKM2TPM(df_tpm[,i])
}
df_tpm = df_tpm[df_tpm$gene_biotype == "protein_coding",]
library(ggplot2)
data.pca <- prcomp(t(df_tpm[,2:14]))
id = c("C1","C2","C3","M1","M2","M3","N2","N3","N5")
df_pca = prcomp(t(df_tpm[,id]))
df_pca = data.frame(df_pca$x, Treatment = substr(colnames(df_tpm[,id]),1,1))  
ggplot(df_pca,aes(x=PC1,y=PC2,color=Treatment))+ geom_point(cex=3) +
  geom_vline(xintercept = 0, linetype=2) + geom_hline(yintercept = 0, linetype=2) +
  theme_bw()

# load count
df_count = read.delim2("matrix_gene.count.annot.txt",header = T, sep='\t',stringsAsFactors = F)
df_count = df_count[,c(1:14,19,20)]
df_count = df_count[df_count$gene_biotype == "protein_coding",]
id = c("C1","C2","C3","M1","M2","M3","N2","N3","N5")
df_count = df_count[,c("gene_id","gene_name",id)]
# DE analysis
library(DESeq2)
coldata = data.frame(Treat = substr(id,1,1),stringsAsFactors = F)
coldata$Treat = factor(coldata$Treat, levels = c("N","C","M"))
rownames(coldata) = id
dds = DESeq(DESeqDataSetFromMatrix(countData = df_count[,id], colData = coldata, design= ~Treat))
de_nc = as.data.frame(results(dds, contrast = c('Treat', 'N', 'C')))
de_nm = as.data.frame(results(dds, contrast = c('Treat', 'N', 'M')))
de_cm = as.data.frame(results(dds, contrast = c('Treat', 'C', 'M')))
de_cm = cbind(df_count[,1:2], de_cm)
de_cm = de_cm[!is.na(de_cm$log2FoldChange),]
#ind = !is.na(de_cm$padj) & de_cm$padj < 0.05 & abs(de_cm$log2FoldChange) >= log(1.5)
#de_cm = de_cm[ind,]

# 2.下载28种免疫细胞的参考基因集 (IRGs) Cell Rep. 2017 Jan 3;18(1):248-262
# HOM_MouseHumanSequence = read.table("HOM_MouseHumanSequence.rpt",header = T, sep='\t',stringsAsFactors = F)
# HOM_MouseHumanSequence = HOM_MouseHumanSequence[,c(1,2,4)]
# HOM_MouseHumanSequence = split(HOM_MouseHumanSequence, HOM_MouseHumanSequence$DB.Class.Key)
# IDtrans = data.frame()
# for(i in 1:length(HOM_MouseHumanSequence)){
#   df = HOM_MouseHumanSequence[[i]]
#   Human = unique(df$Symbol[df$Common.Organism.Name == "human"])
#   Mouse = unique(df$Symbol[df$Common.Organism.Name != "human"])
#   if(length(Human) == 0){Human = NA}
#   if(length(Mouse) == 0){Mouse = NA}
#   IDtrans = rbind(IDtrans, data.frame(Human = Human, Mouse= Mouse, stringsAsFactors = F))
# }
# save(IDtrans, file = "Human_Mouse_IDtrans.RData")
library("AnnotationDbi")
library("org.Hs.eg.db")
load("Human_Mouse_IDtrans.RData")
IRGs = read.csv("mmc3.csv",stringsAsFactors = F)
IDtrans = IDtrans[!is.na(IDtrans$Human) & !is.na(IDtrans$Mouse),]
IRGs = merge(IRGs, IDtrans, by.x = "Metagene", by.y = "Human")
IRGs$Metagene = mapIds(org.Hs.eg.db, IRGs$Metagene, 'ENTREZID', 'SYMBOL')

# 3. run ssGSEA
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
de_cm = merge(IDtrans, de_cm, by.x = "Mouse", by.y = "")
geneList = de_cm[!duplicated(de_cm$gene_name),]
geneList = geneList$log2FoldChange
names(geneList) = de_cm$gene_name[!duplicated(de_cm$gene_name)]
geneList = geneList[order(geneList, decreasing = T)]
#names(geneList) = mapIds(org.Hs.eg.db, names(geneList), 'ENTREZID', 'SYMBOL')
immune_gmt = IRGs[,c(2,4)]
colnames(immune_gmt) = c("term", "gene")

GSEA_Res = data.frame()
for(term in unique(immune_gmt$term)){
  set.seed(15)
  tmp = GSEA(geneList, TERM2GENE = immune_gmt[immune_gmt$term == term,],
              pvalueCutoff = 1)@result
  GSEA_Res = rbind(GSEA_Res, tmp)
}
write.csv(GSEA_Res, file = "GSEA_Res.csv", row.names = F, quote = F)
