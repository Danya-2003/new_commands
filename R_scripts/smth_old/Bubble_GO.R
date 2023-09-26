#library(clusterProfiler)
#library("org.Hs.eg.db")
#library(biomaRt)
#library(GOplot)

#Retrieve the ensembl info
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
dsets = listDatasets(ensembl)
bm <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

#GO_analysis

peak_plus_logfc_minus <- read.table(in_peak_plus_logfc_minus, header =TRUE, sep = ",")


genes <- (peak_plus_logfc_minus$feature) # ENSEMBL ids

#ego <- enrichGO(gene         = genes,
#             OrgDb         = org.Hs.eg.db,
#            keyType       = 'ENSEMBL',
#           ont           = "ALL",
#          pAdjustMethod = "BH",
#         qvalueCutoff  = 0.05)

#write.table(data.frame(ego), out_peak_minus_logfc_minus,
#          sep = '\t',  quote = F,  row.names = F)

ego_MF <- enrichGO(gene         = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)


ego_BP <- enrichGO(gene         = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05)

#ego_CC <- enrichGO(gene         = genes,
#            OrgDb         = org.Hs.eg.db,
#            keyType       = 'ENSEMBL',
#            ont           = "CC",
#            pAdjustMethod = "BH",
#            qvalueCutoff  = 0.05)

genelist <- subset(peak_plus_logfc_minus, select = c("feature", "logfc"))
colnames(genelist)[colnames(genelist) == 'feature'] <- 'ID'
colnames(genelist)[colnames(genelist) == 'logfc'] <- 'logFC'

#data <- read.table(out_peak_plus_logfc_minus, header =TRUE, sep = "\t")
#geneid <- cSplit(subset(data, select = "geneID"), "geneID", "/")
#geneid <- geneid %>% unite(geneid, all_of(colnames(geneid)), sep = ", ", na.rm = TRUE)
#subdata <- subset(data, select = c("ONTOLOGY", "ID", "Description"))
#colnames(subdata)[colnames(subdata) == 'ONTOLOGY'] <- 'Category'
#colnames(subdata)[colnames(subdata) == 'Description'] <- 'Term'
#subdata <- cbind(subdata, geneid)
#colnames(subdata)[colnames(subdata) == 'geneid'] <- 'Genes'
#subdata <- cbind(subdata, data$p.adjust)
#colnames(subdata)[colnames(subdata) == 'data$p.adjust'] <- 'adj_pval'

data <- as.data.frame(cbind (ego_MF$ID, ego_MF$Description, ego_MF$p.adjust, ego_MF$geneID))
geneid <- cSplit(subset(data, select = "V4"), "V4", "/")
geneid <- geneid %>% unite(geneid, all_of(colnames(geneid)), sep = ", ", na.rm = TRUE)
subdata <- subset(data, select = c("V1", "V2"))
colnames(subdata)[colnames(subdata) == 'V1'] <- 'ID'
colnames(subdata)[colnames(subdata) == 'V2'] <- 'Term'
subdata <- cbind(subdata, geneid)
colnames(subdata)[colnames(subdata) == 'geneid'] <- 'Genes'
subdata <- cbind(subdata, data$V3)
colnames(subdata)[colnames(subdata) == 'data$V3'] <- 'adj_pval'
subdata <- cbind(Category = "MF", subdata)


bubble_table_MF <- circle_dat(subdata, genelist)
bubble_table_MF$adj_pval <- as.numeric(bubble_table_MF$adj_pval)

data <- as.data.frame(cbind (ego_BP$ID, ego_BP$Description, ego_BP$p.adjust, ego_BP$geneID))
geneid <- cSplit(subset(data, select = "V4"), "V4", "/")
geneid <- geneid %>% unite(geneid, all_of(colnames(geneid)), sep = ", ", na.rm = TRUE)
subdata <- subset(data, select = c("V1", "V2"))
colnames(subdata)[colnames(subdata) == 'V1'] <- 'ID'
colnames(subdata)[colnames(subdata) == 'V2'] <- 'Term'
subdata <- cbind(subdata, geneid)
colnames(subdata)[colnames(subdata) == 'geneid'] <- 'Genes'
subdata <- cbind(subdata, data$V3)
colnames(subdata)[colnames(subdata) == 'data$V3'] <- 'adj_pval'
subdata <- cbind(Category = "BP", subdata)

bubble_table_BP <- circle_dat(subdata, genelist)
bubble_table_BP$adj_pval <- as.numeric(bubble_table_BP$adj_pval)

