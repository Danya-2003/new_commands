#BiocManager::install("MEDIPS")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#install.packages("~/Downloads/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz", repos = NULL, type = "source")

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-9_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-10.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
KD = c(kd_01,kd_02)

nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-19.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-20.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
NC = c(nc_01,nc_02)

kd_01_Input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-9_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_Input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-10.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_input = c(kd_01_Input, kd_02_Input)

nc_01_Input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-19.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_Input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-20.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_input = c(nc_01_Input,nc_02_Input)

CS = MEDIPS.couplingVector(pattern = 'CG', refObj = NC[[1]])

resultTable = MEDIPS.meth(MSet1 = NC, MSet2 = KD, CSet = CS, ISet1 = nc_input, ISet2 = kd_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)
options(scipen = 999)

resultTable_kd_01 = MEDIPS.meth(MSet1 = NC, MSet2 = kd_01, CSet = CS, ISet1 = nc_input, ISet2 = kd_01_Input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
sig_01 = MEDIPS.selectSig(results=resultTable_kd_01, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)
options(scipen = 999)

resultTable_kd_02 = MEDIPS.meth(MSet1 = NC, MSet2 = kd_02, CSet = CS, ISet1 = nc_input, ISet2 = kd_02_Input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
sig_02 = MEDIPS.selectSig(results=resultTable_kd_02, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)
options(scipen = 999)

write.table(file='~/Downloads/ctd/diffpeak/medips/medips_ctd.txt', sig, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/ctd/diffpeak/medips/medips_ctd_01.txt', sig_01, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/ctd/diffpeak/medips/medips_ctd_03.txt', sig_02, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


results_df <- sig[c(1,2,3,31,28)]
results_df <- as.data.frame(results_df)
good <- subset(results_df, results_df$edgeR.adj.p.value < 0.05)
write.table(file='~/Downloads/ctd/diffpeak/medips/DEpeaks_medips_ctd.csv', good, quote = FALSE, row.names=F, col.names=F, sep ='\t')
good <- na.omit(good)
good <- good[,c(1,2,3)]
write.table(file='~/Downloads/ctd/diffpeak/medips/DEpeaks_medips_ctd.bed', good, quote = FALSE, row.names=F, col.names=F, sep ='\t')
#colnames(good) <- c('seqnames', 'start', 'end', 'fdr', 'logfold')
#peaks <- toGRanges(good, package="ChIPpeakAnno")

results_01_df <- sig_01[c(1,2,3,25,22)]
results_01_df <- as.data.frame(results_01_df)
good_01 <- subset(results_01_df, results_01_df$edgeR.adj.p.value < 0.05)
write.table(file='~/Downloads/ctd/diffpeak/medips/DEpeaks_medips_ctd_01.csv', good_01, quote = FALSE, row.names=F, col.names=F, sep ='\t')
good_01 <- na.omit(good_01)
good_01 <- good_01[,c(1,2,3)]
write.table(file='~/Downloads/ctd/diffpeak/medips/DEpeaks_medips_ctd_01.bed', good_01, quote = FALSE, row.names=F, col.names=F, sep ='\t')
#colnames(good_01) <- c('seqnames', 'start', 'end', 'fdr', 'logfold')
#peaks_01 <- toGRanges(good_01, package="ChIPpeakAnno")

results_02_df <- sig_02[c(1,2,3,25,22)]
results_02_df <- as.data.frame(results_02_df)
good_02 <- subset(results_02_df, results_02_df$edgeR.adj.p.value < 0.05)
write.table(file='~/Downloads/ctd/diffpeak/medips/DEpeaks_medips_ctd_03.csv', good_02, quote = FALSE, row.names=F, col.names=F, sep ='\t')
good_02 <- na.omit(good_02)
good_02 <- good_02[,c(1,2,3)]
write.table(file='~/Downloads/ctd/diffpeak/medips/DEpeaks_medips_ctd_03.bed', good_02, quote = FALSE, row.names=F, col.names=F, sep ='\t')
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
write.table(file='~/Downloads/ctd/anno_medips_ctd.csv', annotatedPeak, quote = FALSE, row.names=F, col.names=F, sep ='\t')
write.table(file='~/Downloads/ctd/anno_medips_ctd_01.csv', annotatedPeak_01, quote = FALSE, row.names=F, col.names=F, sep ='\t')
write.table(file='~/Downloads/ctd/anno_medips_ctd_02.csv', annotatedPeak_02, quote = FALSE, row.names=F, col.names=F, sep ='\t')


