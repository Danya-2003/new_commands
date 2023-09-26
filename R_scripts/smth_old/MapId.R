library("org.Hs.eg.db")
library(biomaRt)
library(tibble)

#data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno_de/CTD-2587H24.5_01_sicer.both.de_only.csv"
#data2 <- "C:/ChIP-seq/fantom6_signatures/CTD-2587H24.5_ASO_G0267577_01.csv"
#out_minus <- "C:/ChIP-seq/BINGO/CTD-2587H24.5_01_sicer.both.de_only.minus.csv"
#out_plus <- "C:/ChIP-seq/BINGO/CTD-2587H24.5_01_sicer.both.de_only.plus.csv"


file <- read.table(data, header =TRUE, sep = "\t")
fantom <- read.table(data2, header =TRUE, sep = "\t")
gene <- file$feature
rownames(fantom) <- fantom$geneID
logfc <- fantom[gene, 8]
file <- cbind(file, logfc)
file_minus <- subset(file, file$logfc < 0)
file_plus <- subset(file, file$logfc > 0)

file_minus$SYMBOL <- mapIds (org.Hs.eg.db, keys = file_minus$feature, keytype = "ENSEMBL", column = "SYMBOL")
file_minus$ENTREZID <- mapIds (org.Hs.eg.db, keys = file_minus$feature, keytype = "ENSEMBL", column = "ENTREZID")
bingo_minus <- as.data.frame(cbind(file_minus$SYMBOL, HGNC = "HGNC:", ENSEMBL = "Ensembl:"))
bingo_minus$HGNC <- paste(bingo_minus$HGNC, file_minus$ENTREZID, sep = "")
bingo_minus$ENSEMBL <- paste(bingo_minus$ENSEMBL, file_minus$feature, sep = "")

write.table(data.frame(bingo_minus), out_minus, sep = '\t',  quote = F,  row.names = F)

file_plus$SYMBOL <- mapIds (org.Hs.eg.db, keys = file_plus$feature, keytype = "ENSEMBL", column = "SYMBOL")
file_plus$ENTREZID <- mapIds (org.Hs.eg.db, keys = file_plus$feature, keytype = "ENSEMBL", column = "ENTREZID")
bingo_plus <- as.data.frame(cbind(file_plus$SYMBOL, HGNC = "HGNC:", ENSEMBL = "Ensembl:"))
bingo_plus$HGNC <- paste(bingo_plus$HGNC, file_plus$ENTREZID, sep = "")
bingo_plus$ENSEMBL <- paste(bingo_plus$ENSEMBL, file_plus$feature, sep = "")

write.table(data.frame(bingo_plus), out_plus, sep = '\t',  quote = F,  row.names = F)

#keytypes(org.Hs.eg.db)
#column = "GENENAME", column = "ENTREZID", colum = "SYMBOL"

# map from one annotation to another using biomart

#m <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#listAttributes(m)
#listFilterOptions(m, "ensembl_gene_id")
#listFilters(m)
#file_minus$hgnc_symbol <- getBM(mart = m, attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = file_minus$feature, uniqueRows = FALSE)