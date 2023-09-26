#BiocManager::install()
#BiocManager::install("DESeq2")
#BiocManager::install("ChIPpeakAnno")
#BiocManager::install("EnsDb.Hsapiens.v86")

library(DESeq2)
library(ChIPpeakAnno)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86) ##GRCh38

cov <- read.table('~/Downloads/ctd/diffpeak/coverage.txt')
cov$id <- with(cov, paste0(V1,'_',V2,'_',V3))
group <- factor(c('Kd','Kd','NC','NC'))
cnts <- cov[,4:7]
dds <- DESeqDataSetFromMatrix(countData = cnts,colData = DataFrame(group),design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c('group','Kd','NC'))
results <- cbind(as.character(cov$V1),cov$V2,cov$V3,res$padj,res$log2FoldChange)
options(scipen = 999)
write.table(file='~/Downloads/ctd/diffpeak/deseq/deseq2_ctd.csv', results, quote = FALSE, row.names=F, col.names=F, sep ='\t')

cov <- read.table('~/Downloads/ctd/diffpeak/coverage_01.txt')
cov$id <- with(cov, paste0(V1,'_',V2,'_',V3))
group <- factor(c('Kd_01','NC','NC'))
cnts <- cov[,4:6]
dds <- DESeqDataSetFromMatrix(countData = cnts, colData = DataFrame(group), design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c('group','Kd_01','NC'))
results_01 <- cbind(as.character(cov$V1),cov$V2,cov$V3,res$padj,res$log2FoldChange)
options(scipen = 999)
write.table(file='~/Downloads/ctd/diffpeak/deseq/deseq2_ctd_01.csv', results, quote = FALSE, row.names=F, col.names=F, sep ='\t')

cov <- read.table('~/Downloads/ctd/diffpeak/coverage_03.txt')
cov$id <- with(cov, paste0(V1,'_',V2,'_',V3))
group <- factor(c('Kd_02','NC','NC'))
cnts <- cov[,4:6]
dds <- DESeqDataSetFromMatrix(countData = cnts,colData = DataFrame(group),design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c('group','Kd_02','NC'))
results_02 <- cbind(as.character(cov$V1),cov$V2,cov$V3,res$padj,res$log2FoldChange)
options(scipen = 999)
write.table(file='~/Downloads/ctd/diffpeak/deseq/deseq2_ctd_03.csv', results, quote = FALSE, row.names=F, col.names=F, sep ='\t')




results_df <- as.data.frame(results)
good <- subset(results_df, results_df$V4 < 0.05)
write.table(file='~/Downloads/ctd/diffpeak/deseq/DEpeaks_deseq_ctd.csv', good, quote = FALSE, row.names=F, col.names=F, sep ='\t')
good <- na.omit(good)
good_bed <- good [,c(1,2,3)]
write.table(file='~/Downloads/ctd/diffpeak/deseq/DEpeaks_deseq_ctd.bed', good_bed, quote = FALSE, row.names=F, col.names=F, sep ='\t')
#colnames(good) <- c('seqnames', 'start', 'end', 'fdr', 'logfold')
#peaks <- toGRanges(good, package="ChIPpeakAnno")

results_01_df <- as.data.frame(results)
good_01 <- subset(results_01_df, results_01_df$V4 < 0.05)
write.table(file='~/Downloads/ctd/diffpeak/deseq/DEpeaks_deseq_ctd_01.csv', good_01, quote = FALSE, row.names=F, col.names=F, sep ='\t')
good_01 <- na.omit(good_01)
good_bed_01 <- good_01 [,c(1,2,3)]
write.table(file='~/Downloads/ctd/diffpeak/deseq/DEpeaks_deseq_ctd_01.bed', good_bed_01, quote = FALSE, row.names=F, col.names=F, sep ='\t')
#colnames(good_01) <- c('seqnames', 'start', 'end', 'fdr', 'logfold')
#peaks_01 <- toGRanges(good_01, package="ChIPpeakAnno")

results_02_df <- as.data.frame(results)
good_02 <- subset(results_02_df, results_02_df$V4 < 0.05)
write.table(file='~/Downloads/ctd/diffpeak/deseq/DEpeaks_deseq_ctd_03.csv', good_02, quote = FALSE, row.names=F, col.names=F, sep ='\t')
good_02 <- na.omit(good_02)
good_bed_02 <- good_02 [,c(1,2,3)]
write.table(file='~/Downloads/ctd/diffpeak/deseq/DEpeaks_deseq_ctd_03.bed', good_bed_02, quote = FALSE, row.names=F, col.names=F, sep ='\t')
#colnames(good_02) <- c('seqnames', 'start', 'end', 'fdr', 'logfold')
#peaks_02 <- toGRanges(good_02, package="ChIPpeakAnno")

require(data.table)
require(annotate)
require(org.Hs.eg.db)
require(ChIPpeakAnno)
#data(TSS.human.GRCh38)
annoDataEnsDb <- toGRanges(EnsDb.Hsapiens.v86)
annotatedPeak = annotatePeakInBatch(peaks, AnnotationData=annoDataEnsDb)
annotatedPeak_01 = annotatePeakInBatch(peaks_01, AnnotationData=annoDataEnsDb)
annotatedPeak_02 = annotatePeakInBatch(peaks_02, AnnotationData=annoDataEnsDb)
write.table(file='~/Downloads/ctd/anno_deseq_ctd.csv', annotatedPeak, quote = FALSE, row.names=F, col.names=F, sep ='\t')
write.table(file='~/Downloads/ctd/anno_deseq_ctd_01.csv', annotatedPeak_01, quote = FALSE, row.names=F, col.names=F, sep ='\t')
write.table(file='~/Downloads/ctd/anno_deseq_ctd_02.csv', annotatedPeak_02, quote = FALSE, row.names=F, col.names=F, sep ='\t')
