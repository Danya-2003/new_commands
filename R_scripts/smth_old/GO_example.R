
#library(clusterProfiler)
#library("org.Hs.eg.db")
#library(biomaRt)

#Retrieve the ensembl info
#ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
#dsets = listDatasets(ensembl)
#bm <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

#data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.csv"
#data2 <- "C:/ChIP-seq/fantom6_signatures/CTD-2587H24.5_ASO_G0267577_03.csv"
#out_mf_plus <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.mf.plus.csv"
#out_bp_plus <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.bp.plus.csv"
#out_cc_plus <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.cc.plus.csv"
#out_kegg_plus <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.kegg.plus.csv"
#out_mf_minus <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.mf.minus.csv"
#out_bp_minus <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.bp.minus.csv"
#out_cc_minus <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.cc.minus.csv"
#out_kegg_minus <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno_de/CTD-2587H24.5_03_sicer.both.de_only.kegg.minus.csv"

file <- read.table(data, header =TRUE, sep = "\t")
fantom <- read.table(data2, header =TRUE, sep = "\t")
gene <- file$feature
rownames(fantom) <- fantom$geneID
logfc <- fantom[gene, 8]
file <- cbind(file, logfc)
file_minus <- subset(file, file$logfc < 0)
file_plus <- subset(file, file$logfc > 0)

genes <- (file_minus$feature) # ENSEMBL ids

ego_MF <- enrichGO(gene         = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

write.table(data.frame(ego_MF),out_mf_minus,
            sep = '\t',  quote = F,  row.names = F)

ego_BP <- enrichGO(gene         = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

write.table(data.frame(ego_BP), out_bp_minus,
            sep = '\t',  quote = F,  row.names = F)

ego_CC <- enrichGO(gene         = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

write.table(data.frame(ego_CC),out_cc_minus,
            sep = '\t',  quote = F,  row.names = F)

# Filter the ensembl ID
ids <- getBM(filters = "ensembl_gene_id", 
             attributes = c("ensembl_gene_id", 'entrezgene_id', 'entrezgene_accession'), 
             values = genes, mart = bm)

kk <- enrichKEGG(gene         = ids$entrezgene_id,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 qvalueCutoff  = 0.05)

write.table(data.frame(kk), out_kegg_minus,
            sep = '\t',  quote = F,  row.names = F)


genes <- (file_plus$feature) # ? - ? ENSEMBL ids

ego_MF <- enrichGO(gene         = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

write.table(data.frame(ego_MF),out_mf_plus,
            sep = '\t',  quote = F,  row.names = F)

ego_BP <- enrichGO(gene         = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

write.table(data.frame(ego_BP), out_bp_plus,
            sep = '\t',  quote = F,  row.names = F)

ego_CC <- enrichGO(gene         = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

write.table(data.frame(ego_CC),out_cc_plus,
            sep = '\t',  quote = F,  row.names = F)

# Filter the ensembl ID
ids <- getBM(filters = "ensembl_gene_id", 
             attributes = c("ensembl_gene_id", 'entrezgene_id', 'entrezgene_accession'), 
             values = genes, mart = bm)

kk <- enrichKEGG(gene         = ids$entrezgene_id,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 
                 qvalueCutoff  = 0.05)

write.table(data.frame(kk), out_kegg_plus,
            sep = '\t',  quote = F,  row.names = F)

#3browseKEGG(kk, "hsa04110")