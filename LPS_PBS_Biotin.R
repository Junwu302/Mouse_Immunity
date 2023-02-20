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
immune_gmt = IRGs[,c(2,4)]
colnames(immune_gmt) = c("term", "gene")

de_lps_pbs = read.csv("LPS_vs_PBS.csv", header = T, stringsAsFactors = F)
de_lps_pbs = merge(IDtrans, de_lps_pbs, by.x = "Mouse", by.y = "gene_name")
geneList = de_lps_pbs[!duplicated(de_lps_pbs$Mouse),]
geneList = geneList$log2FC
names(geneList) = de_lps_pbs$Mouse[!duplicated(de_lps_pbs$Mouse)]
geneList = geneList[order(geneList, decreasing = T)]

GSEA_Res = data.frame()
for(term in unique(immune_gmt$term)){
  set.seed(15)
  tmp = GSEA(geneList, TERM2GENE = immune_gmt[immune_gmt$term == term,],
             pvalueCutoff = 1)@result
  GSEA_Res = rbind(GSEA_Res, tmp)
}

write.csv(GSEA_Res, file = "LPS_PBS_GSEA_Res.csv", row.names = F, quote = F)

# Biotin_vs_LPS
de_biotin_lps = read.csv("Biotin_vs_LPS.csv", header = T, stringsAsFactors = F)
de_biotin_lps = merge(IDtrans, de_biotin_lps, by.x = "Mouse", by.y = "gene_name")
geneList = de_biotin_lps[!duplicated(de_biotin_lps$Mouse),]
geneList = geneList$log2FC
names(geneList) = de_biotin_lps$Mouse[!duplicated(de_biotin_lps$Mouse)]
geneList = geneList[order(geneList, decreasing = T)]

GSEA_Res = data.frame()
for(term in unique(immune_gmt$term)){
  set.seed(15)
  tmp = GSEA(geneList, TERM2GENE = immune_gmt[immune_gmt$term == term,],
             pvalueCutoff = 1)@result
  GSEA_Res = rbind(GSEA_Res, tmp)
}

write.csv(GSEA_Res, file = "Biotin_LPS_GSEA_Res.csv", row.names = F, quote = F)

