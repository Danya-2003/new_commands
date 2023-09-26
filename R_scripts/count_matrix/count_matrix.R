library("DESeq2")
library("edgeR")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

### JPX

BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-11.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-12.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-17.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-18.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-11.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-12.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-17.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-18.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)

resultTable_01 = MEDIPS.meth(MSet1 = nc_01, MSet2 = kd_01, ISet1 = nc_01_input, ISet2 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
resultTable_02 = MEDIPS.meth(MSet1 = nc_02, MSet2 = kd_02, ISet1 = nc_02_input, ISet2 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')

count_matrix <- as.data.frame(paste(resultTable_01$chr, resultTable_01$start, resultTable_01$stop, sep = "_"))
rownames(count_matrix) <- count_matrix[,1]
count_matrix$nc_01 <- resultTable_01$C.17.bam.counts
count_matrix$nc_02 <- resultTable_02$C.18.bam.counts
count_matrix$kd_01 <- resultTable_01$C.11.bam.counts
count_matrix$kd_02 <- resultTable_02$C.12.bam.counts
count_matrix$nc_01_input <- resultTable_01$I.17.bam.counts
count_matrix$nc_02_input <- resultTable_02$I.18.bam.counts
count_matrix$kd_01_input <- resultTable_01$I.11.bam.counts
count_matrix$kd_02_input <- resultTable_02$I.12.bam.counts
count_matrix <- count_matrix[,-1]

write.table(count_matrix, '~/count_matrix/jpx', quote = FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')

rm(list = ls())

### CTD

BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-9_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-10.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-19.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-20.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-9_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-10.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-19.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-20.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)

resultTable_01 = MEDIPS.meth(MSet1 = nc_01, MSet2 = kd_01, ISet1 = nc_01_input, ISet2 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
resultTable_02 = MEDIPS.meth(MSet1 = nc_02, MSet2 = kd_02, ISet1 = nc_02_input, ISet2 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')

count_matrix <- as.data.frame(paste(resultTable_01$chr, resultTable_01$start, resultTable_01$stop, sep = "_"))
rownames(count_matrix) <- count_matrix[,1]
count_matrix$nc_01 <- resultTable_01$C.9_.bam.counts
count_matrix$nc_02 <- resultTable_02$C.10.bam.counts
count_matrix$kd_01 <- resultTable_01$C.19.bam.counts
count_matrix$kd_02 <- resultTable_02$C.20.bam.counts
count_matrix$nc_01_input <- resultTable_01$I.9_.bam.counts
count_matrix$nc_02_input <- resultTable_02$I.10.bam.counts
count_matrix$kd_01_input <- resultTable_01$I.19.bam.counts
count_matrix$kd_02_input <- resultTable_02$I.20.bam.counts
count_matrix <- count_matrix[,-1]

write.table(count_matrix, '~/count_matrix/ctd', quote = FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')

rm(list = ls())

### FGD5-AS1
BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-1_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-2_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-21.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-22.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-1_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-2_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-21.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-22.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)

resultTable_01 = MEDIPS.meth(MSet1 = nc_01, MSet2 = kd_01, ISet1 = nc_01_input, ISet2 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
resultTable_02 = MEDIPS.meth(MSet1 = nc_02, MSet2 = kd_02, ISet1 = nc_02_input, ISet2 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')

count_matrix <- as.data.frame(paste(resultTable_01$chr, resultTable_01$start, resultTable_01$stop, sep = "_"))
rownames(count_matrix) <- count_matrix[,1]
count_matrix$nc_01 <- resultTable_01$C.1_.bam.counts
count_matrix$nc_02 <- resultTable_02$C.2_.bam.counts
count_matrix$kd_01 <- resultTable_01$C.21.bam.counts
count_matrix$kd_02 <- resultTable_02$C.22.bam.counts
count_matrix$nc_01_input <- resultTable_01$I.1_.bam.counts
count_matrix$nc_02_input <- resultTable_02$I.2_.bam.counts
count_matrix$kd_01_input <- resultTable_01$I.21.bam.counts
count_matrix$kd_02_input <- resultTable_02$I.22.bam.counts
count_matrix <- count_matrix[,-1]

write.table(count_matrix, '~/count_matrix/fgd5', quote = FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')

rm(list = ls())

### EMX2OS
BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-3_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-4_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-21.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-22.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-3_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-4_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-21.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-22.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)

resultTable_01 = MEDIPS.meth(MSet1 = nc_01, MSet2 = kd_01, ISet1 = nc_01_input, ISet2 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
resultTable_02 = MEDIPS.meth(MSet1 = nc_02, MSet2 = kd_02, ISet1 = nc_02_input, ISet2 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')

count_matrix <- as.data.frame(paste(resultTable_01$chr, resultTable_01$start, resultTable_01$stop, sep = "_"))
rownames(count_matrix) <- count_matrix[,1]
count_matrix$nc_01 <- resultTable_01$C.3_.bam.counts
count_matrix$nc_02 <- resultTable_02$C.4_.bam.counts
count_matrix$kd_01 <- resultTable_01$C.21.bam.counts
count_matrix$kd_02 <- resultTable_02$C.22.bam.counts
count_matrix$nc_01_input <- resultTable_01$I.3_.bam.counts
count_matrix$nc_02_input <- resultTable_02$I.4_.bam.counts
count_matrix$kd_01_input <- resultTable_01$I.21.bam.counts
count_matrix$kd_02_input <- resultTable_02$I.22.bam.counts
count_matrix <- count_matrix[,-1]

write.table(count_matrix, '~/count_matrix/emx2os', quote = FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')

rm(list = ls())


### AC005592.2
BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-5_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-6_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-17.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-18.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-5_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-6_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-17.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-18.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)

resultTable_01 = MEDIPS.meth(MSet1 = nc_01, MSet2 = kd_01, ISet1 = nc_01_input, ISet2 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
resultTable_02 = MEDIPS.meth(MSet1 = nc_02, MSet2 = kd_02, ISet1 = nc_02_input, ISet2 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')

count_matrix <- as.data.frame(paste(resultTable_01$chr, resultTable_01$start, resultTable_01$stop, sep = "_"))
rownames(count_matrix) <- count_matrix[,1]
count_matrix$nc_01 <- resultTable_01$C.5_.bam.counts
count_matrix$nc_02 <- resultTable_02$C.6_.bam.counts
count_matrix$kd_01 <- resultTable_01$C.17.bam.counts
count_matrix$kd_02 <- resultTable_02$C.18.bam.counts
count_matrix$nc_01_input <- resultTable_01$I.5_.bam.counts
count_matrix$nc_02_input <- resultTable_02$I.6_.bam.counts
count_matrix$kd_01_input <- resultTable_01$I.17.bam.counts
count_matrix$kd_02_input <- resultTable_02$I.18.bam.counts
count_matrix <- count_matrix[,-1]

write.table(count_matrix, '~/count_matrix/ac', quote = FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')

rm(list = ls())

### RP11-398K22.12
BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-7_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-8_.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-23.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-24.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-7_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-8_.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-23.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-24.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)

resultTable_01 = MEDIPS.meth(MSet1 = nc_01, MSet2 = kd_01, ISet1 = nc_01_input, ISet2 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
resultTable_02 = MEDIPS.meth(MSet1 = nc_02, MSet2 = kd_02, ISet1 = nc_02_input, ISet2 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')

count_matrix <- as.data.frame(paste(resultTable_01$chr, resultTable_01$start, resultTable_01$stop, sep = "_"))
rownames(count_matrix) <- count_matrix[,1]
count_matrix$nc_01 <- resultTable_01$C.7_.bam.counts
count_matrix$nc_02 <- resultTable_02$C.8_.bam.counts
count_matrix$kd_01 <- resultTable_01$C.23.bam.counts
count_matrix$kd_02 <- resultTable_02$C.24.bam.counts
count_matrix$nc_01_input <- resultTable_01$I.7_.bam.counts
count_matrix$nc_02_input <- resultTable_02$I.8_.bam.counts
count_matrix$kd_01_input <- resultTable_01$I.23.bam.counts
count_matrix$kd_02_input <- resultTable_02$I.24.bam.counts
count_matrix <- count_matrix[,-1]

write.table(count_matrix, '~/count_matrix/rp11', quote = FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')

rm(list = ls())

### MAPKAPK5-AS1_H3K27me3
BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-13.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-15.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-25.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-27.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-13.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-15.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-25.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-27.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)

resultTable_01 = MEDIPS.meth(MSet1 = nc_01, MSet2 = kd_01, ISet1 = nc_01_input, ISet2 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
resultTable_02 = MEDIPS.meth(MSet1 = nc_02, MSet2 = kd_02, ISet1 = nc_02_input, ISet2 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')

count_matrix <- as.data.frame(paste(resultTable_01$chr, resultTable_01$start, resultTable_01$stop, sep = "_"))
rownames(count_matrix) <- count_matrix[,1]
count_matrix$nc_01 <- resultTable_01$C.13.bam.counts
count_matrix$nc_02 <- resultTable_02$C.15.bam.counts
count_matrix$kd_01 <- resultTable_01$C.25.bam.counts
count_matrix$kd_02 <- resultTable_02$C.27.bam.counts
count_matrix$nc_01_input <- resultTable_01$I.13.bam.counts
count_matrix$nc_02_input <- resultTable_02$I.15.bam.counts
count_matrix$kd_01_input <- resultTable_01$I.25.bam.counts
count_matrix$kd_02_input <- resultTable_02$I.27.bam.counts
count_matrix <- count_matrix[,-1]

write.table(count_matrix, '~/count_matrix/mapkapk_h3k27me3', quote = FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')

rm(list = ls())

### MAPKAPK5-AS1_H3K4me3
BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-14.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-16.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-26.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/sambamba/C-28.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-14.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-16.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_01_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-26.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/sambamba/I-28.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)

resultTable_01 = MEDIPS.meth(MSet1 = nc_01, MSet2 = kd_01, ISet1 = nc_01_input, ISet2 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
resultTable_02 = MEDIPS.meth(MSet1 = nc_02, MSet2 = kd_02, ISet1 = nc_02_input, ISet2 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')

count_matrix <- as.data.frame(paste(resultTable_01$chr, resultTable_01$start, resultTable_01$stop, sep = "_"))
rownames(count_matrix) <- count_matrix[,1]
count_matrix$nc_01 <- resultTable_01$C.14.bam.counts
count_matrix$nc_02 <- resultTable_02$C.16.bam.counts
count_matrix$kd_01 <- resultTable_01$C.26.bam.counts
count_matrix$kd_02 <- resultTable_02$C.28.bam.counts
count_matrix$nc_01_input <- resultTable_01$I.14.bam.counts
count_matrix$nc_02_input <- resultTable_02$I.16.bam.counts
count_matrix$kd_01_input <- resultTable_01$I.26.bam.counts
count_matrix$kd_02_input <- resultTable_02$I.28.bam.counts
count_matrix <- count_matrix[,-1]

write.table(count_matrix, '~/count_matrix/mapkapk_h3k4me3', quote = FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')
rm(list = ls())


### JPX

count_matrix <- read.table("~/count_matrix/jpx", header = TRUE, sep = "\t")

png(filename="~/scatter/jpx_nc_01.png")
plot(count_matrix$nc_01)
dev.off()

png(filename="~/scatter/jpx_nc_02.png")
plot(count_matrix$nc_02)
dev.off()

png(filename="~/scatter/jpx_kd_01.png")
plot(count_matrix$kd_01)
dev.off()

png(filename="~/scatter/jpx_kd_02.png")
plot(count_matrix$kd_02)
dev.off()

png(filename="~/scatter/jpx_nc_01_input.png")
plot(count_matrix$nc_01_input)
dev.off()

png(filename="~/scatter/jpx_nc_02_input.png")
plot(count_matrix$nc_02_input)
dev.off()

png(filename="~/scatter/jpx_kd_01_input.png")
plot(count_matrix$kd_01_input)
dev.off()

png(filename="~/scatter/jpx_kd_02_input.png")
plot(count_matrix$kd_02_input)
dev.off()

hist(count_matrix$nc_01)
hist(count_matrix$nc_02)
hist(count_matrix$kd_01)
hist(count_matrix$kd_02)
hist(count_matrix$nc_01_input)
hist(count_matrix$nc_02_input)
hist(count_matrix$kd_01_input)
hist(count_matrix$kd_02_input)

count_matrix_1 <- count_matrix[(count_matrix$nc_01!="0") & 
                     (count_matrix$nc_02!="0") & 
                     (count_matrix$kd_01!="0") & 
                     (count_matrix$kd_02!="0") & 
                     (count_matrix$nc_01_input!="0") & 
                     (count_matrix$nc_02_input!="0") & 
                     (count_matrix$kd_01_input!="0") & 
                     (count_matrix$kd_02_input!="0"), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01>"200") & 
                                 (count_matrix$nc_02>"200") & 
                                 (count_matrix$kd_01>"200") & 
                                 (count_matrix$kd_02>"200") & 
                                 (count_matrix$nc_01_input>"200") & 
                                 (count_matrix$nc_02_input>"200") & 
                                 (count_matrix$kd_01_input>"200") & 
                                 (count_matrix$kd_02_input>"200"), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())

### CTD


count_matrix <- read.table("~/count_matrix/ctd", header = TRUE, sep = "\t")

png(filename="~/scatter/ctd_nc_01.png")
plot(count_matrix$nc_01)
dev.off()

png(filename="~/scatter/ctd_nc_02.png")
plot(count_matrix$nc_02)
dev.off()

png(filename="~/scatter/ctd_kd_01.png")
plot(count_matrix$kd_01)
dev.off()

png(filename="~/scatter/ctd_kd_02.png")
plot(count_matrix$kd_02)
dev.off()

png(filename="~/scatter/ctd_nc_01_input.png")
plot(count_matrix$nc_01_input)
dev.off()

png(filename="~/scatter/ctd_nc_02_input.png")
plot(count_matrix$nc_02_input)
dev.off()

png(filename="~/scatter/ctd_kd_01_input.png")
plot(count_matrix$kd_01_input)
dev.off()

png(filename="~/scatter/ctd_kd_02_input.png")
plot(count_matrix$kd_02_input)
dev.off()

hist(count_matrix$nc_01)
hist(count_matrix$nc_02)
hist(count_matrix$kd_01)
hist(count_matrix$kd_02)
hist(count_matrix$nc_01_input)
hist(count_matrix$nc_02_input)
hist(count_matrix$kd_01_input)
hist(count_matrix$kd_02_input)

count_matrix_1 <- count_matrix[(count_matrix$nc_01!="0") & 
                                 (count_matrix$nc_02!="0") & 
                                 (count_matrix$kd_01!="0") & 
                                 (count_matrix$kd_02!="0") & 
                                 (count_matrix$nc_01_input!="0") & 
                                 (count_matrix$nc_02_input!="0") & 
                                 (count_matrix$kd_01_input!="0") & 
                                 (count_matrix$kd_02_input!="0"), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01>"200") & 
                                 (count_matrix$nc_02>"200") & 
                                 (count_matrix$kd_01>"200") & 
                                 (count_matrix$kd_02>"200") & 
                                 (count_matrix$nc_01_input>"200") & 
                                 (count_matrix$nc_02_input>"200") & 
                                 (count_matrix$kd_01_input>"200") & 
                                 (count_matrix$kd_02_input>"200"), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())

### FGD5-AS1

count_matrix <- read.table("~/count_matrix/fgd5", header = TRUE, sep = "\t")

png(filename="~/scatter/fgd5_nc_01.png")
plot(count_matrix$nc_01)
dev.off()

png(filename="~/scatter/fgd5_nc_02.png")
plot(count_matrix$nc_02)
dev.off()

png(filename="~/scatter/fgd5_kd_01.png")
plot(count_matrix$kd_01)
dev.off()

png(filename="~/scatter/fgd5_kd_02.png")
plot(count_matrix$kd_02)
dev.off()

png(filename="~/scatter/fgd5_nc_01_input.png")
plot(count_matrix$nc_01_input)
dev.off()

png(filename="~/scatter/fgd5_nc_02_input.png")
plot(count_matrix$nc_02_input)
dev.off()

png(filename="~/scatter/fgd5_kd_01_input.png")
plot(count_matrix$kd_01_input)
dev.off()

png(filename="~/scatter/fgd5_kd_02_input.png")
plot(count_matrix$kd_02_input)
dev.off()

hist(count_matrix$nc_01)
hist(count_matrix$nc_02)
hist(count_matrix$kd_01)
hist(count_matrix$kd_02)
hist(count_matrix$nc_01_input)
hist(count_matrix$nc_02_input)
hist(count_matrix$kd_01_input)
hist(count_matrix$kd_02_input)

count_matrix_1 <- count_matrix[(count_matrix$nc_01!="0") & 
                                 (count_matrix$nc_02!="0") & 
                                 (count_matrix$kd_01!="0") & 
                                 (count_matrix$kd_02!="0") & 
                                 (count_matrix$nc_01_input!="0") & 
                                 (count_matrix$nc_02_input!="0") & 
                                 (count_matrix$kd_01_input!="0") & 
                                 (count_matrix$kd_02_input!="0"), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01>"200") & 
                                 (count_matrix$nc_02>"200") & 
                                 (count_matrix$kd_01>"200") & 
                                 (count_matrix$kd_02>"200") & 
                                 (count_matrix$nc_01_input>"200") & 
                                 (count_matrix$nc_02_input>"200") & 
                                 (count_matrix$kd_01_input>"200") & 
                                 (count_matrix$kd_02_input>"200"), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())


### EMX2OS

count_matrix <- read.table("~/count_matrix/emx2os", header = TRUE, sep = "\t")

png(filename="~/scatter/emx2os_nc_01.png")
plot(count_matrix$nc_01)
dev.off()

png(filename="~/scatter/emx2os_nc_02.png")
plot(count_matrix$nc_02)
dev.off()

png(filename="~/scatter/emx2os_kd_01.png")
plot(count_matrix$kd_01)
dev.off()

png(filename="~/scatter/emx2os_kd_02.png")
plot(count_matrix$kd_02)
dev.off()

png(filename="~/scatter/emx2os_nc_01_input.png")
plot(count_matrix$nc_01_input)
dev.off()

png(filename="~/scatter/emx2os_nc_02_input.png")
plot(count_matrix$nc_02_input)
dev.off()

png(filename="~/scatter/emx2os_kd_01_input.png")
plot(count_matrix$kd_01_input)
dev.off()

png(filename="~/scatter/emx2os_kd_02_input.png")
plot(count_matrix$kd_02_input)
dev.off()

hist(count_matrix$nc_01)
hist(count_matrix$nc_02)
hist(count_matrix$kd_01)
hist(count_matrix$kd_02)
hist(count_matrix$nc_01_input)
hist(count_matrix$nc_02_input)
hist(count_matrix$kd_01_input)
hist(count_matrix$kd_02_input)

count_matrix_1 <- count_matrix[(count_matrix$nc_01!="0") & 
                                 (count_matrix$nc_02!="0") & 
                                 (count_matrix$kd_01!="0") & 
                                 (count_matrix$kd_02!="0") & 
                                 (count_matrix$nc_01_input!="0") & 
                                 (count_matrix$nc_02_input!="0") & 
                                 (count_matrix$kd_01_input!="0") & 
                                 (count_matrix$kd_02_input!="0"), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01>"200") & 
                                 (count_matrix$nc_02>"200") & 
                                 (count_matrix$kd_01>"200") & 
                                 (count_matrix$kd_02>"200") & 
                                 (count_matrix$nc_01_input>"200") & 
                                 (count_matrix$nc_02_input>"200") & 
                                 (count_matrix$kd_01_input>"200") & 
                                 (count_matrix$kd_02_input>"200"), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())


### AC005592.2

count_matrix <- read.table("~/count_matrix/ac", header = TRUE, sep = "\t")

png(filename="~/scatter/ac_nc_01.png")
plot(count_matrix$nc_01)
dev.off()

png(filename="~/scatter/ac_nc_02.png")
plot(count_matrix$nc_02)
dev.off()

png(filename="~/scatter/ac_kd_01.png")
plot(count_matrix$kd_01)
dev.off()

png(filename="~/scatter/ac_kd_02.png")
plot(count_matrix$kd_02)
dev.off()

png(filename="~/scatter/ac_nc_01_input.png")
plot(count_matrix$nc_01_input)
dev.off()

png(filename="~/scatter/ac_nc_02_input.png")
plot(count_matrix$nc_02_input)
dev.off()

png(filename="~/scatter/ac_kd_01_input.png")
plot(count_matrix$kd_01_input)
dev.off()

png(filename="~/scatter/ac_kd_02_input.png")
plot(count_matrix$kd_02_input)
dev.off()


hist(count_matrix$nc_01)
hist(count_matrix$nc_02)
hist(count_matrix$kd_01)
hist(count_matrix$kd_02)
hist(count_matrix$nc_01_input)
hist(count_matrix$nc_02_input)
hist(count_matrix$kd_01_input)
hist(count_matrix$kd_02_input)

count_matrix_1 <- count_matrix[(count_matrix$nc_01!="0") & 
                                 (count_matrix$nc_02!="0") & 
                                 (count_matrix$kd_01!="0") & 
                                 (count_matrix$kd_02!="0") & 
                                 (count_matrix$nc_01_input!="0") & 
                                 (count_matrix$nc_02_input!="0") & 
                                 (count_matrix$kd_01_input!="0") & 
                                 (count_matrix$kd_02_input!="0"), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01>"200") & 
                                 (count_matrix$nc_02>"200") & 
                                 (count_matrix$kd_01>"200") & 
                                 (count_matrix$kd_02>"200") & 
                                 (count_matrix$nc_01_input>"200") & 
                                 (count_matrix$nc_02_input>"200") & 
                                 (count_matrix$kd_01_input>"200") & 
                                 (count_matrix$kd_02_input>"200"), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())


### RP11-398K22.12

count_matrix <- read.table("~/count_matrix/rp11", header = TRUE, sep = "\t")

png(filename="~/scatter/rp11_nc_01.png")
plot(count_matrix$nc_01)
dev.off()

png(filename="~/scatter/rp11_nc_02.png")
plot(count_matrix$nc_02)
dev.off()

png(filename="~/scatter/rp11_kd_01.png")
plot(count_matrix$kd_01)
dev.off()

png(filename="~/scatter/rp11_kd_02.png")
plot(count_matrix$kd_02)
dev.off()

png(filename="~/scatter/rp11_nc_01_input.png")
plot(count_matrix$nc_01_input)
dev.off()

png(filename="~/scatter/rp11_nc_02_input.png")
plot(count_matrix$nc_02_input)
dev.off()

png(filename="~/scatter/rp11_kd_01_input.png")
plot(count_matrix$kd_01_input)
dev.off()

png(filename="~/scatter/rp11_kd_02_input.png")
plot(count_matrix$kd_02_input)
dev.off()

hist(count_matrix$nc_01)
hist(count_matrix$nc_02)
hist(count_matrix$kd_01)
hist(count_matrix$kd_02)
hist(count_matrix$nc_01_input)
hist(count_matrix$nc_02_input)
hist(count_matrix$kd_01_input)
hist(count_matrix$kd_02_input)

count_matrix_1 <- count_matrix[(count_matrix$nc_01!="0") & 
                                 (count_matrix$nc_02!="0") & 
                                 (count_matrix$kd_01!="0") & 
                                 (count_matrix$kd_02!="0") & 
                                 (count_matrix$nc_01_input!="0") & 
                                 (count_matrix$nc_02_input!="0") & 
                                 (count_matrix$kd_01_input!="0") & 
                                 (count_matrix$kd_02_input!="0"), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01>"200") & 
                                 (count_matrix$nc_02>"200") & 
                                 (count_matrix$kd_01>"200") & 
                                 (count_matrix$kd_02>"200") & 
                                 (count_matrix$nc_01_input>"200") & 
                                 (count_matrix$nc_02_input>"200") & 
                                 (count_matrix$kd_01_input>"200") & 
                                 (count_matrix$kd_02_input>"200"), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())

### MAPKAPK5-AS1_H3K27me3

count_matrix <- read.table("~/count_matrix/mapkapk_h3k27me3", header = TRUE, sep = "\t")

png(filename="~/scatter/mapkapk_h3k27me3_nc_01.png")
plot(count_matrix$nc_01)
dev.off()

png(filename="~/scatter/mapkapk_h3k27me3_nc_02.png")
plot(count_matrix$nc_02)
dev.off()

png(filename="~/scatter/mapkapk_h3k27me3_kd_01.png")
plot(count_matrix$kd_01)
dev.off()

png(filename="~/scatter/mapkapk_h3k27me3_kd_02.png")
plot(count_matrix$kd_02)
dev.off()

png(filename="~/scatter/mapkapk_h3k27me3_nc_01_input.png")
plot(count_matrix$nc_01_input)
dev.off()

png(filename="~/scatter/mapkapk_h3k27me3_nc_02_input.png")
plot(count_matrix$nc_02_input)
dev.off()

png(filename="~/scatter/mapkapk_h3k27me3_kd_01_input.png")
plot(count_matrix$kd_01_input)
dev.off()

png(filename="~/scatter/mapkapk_h3k27me3_kd_02_input.png")
plot(count_matrix$kd_02_input)
dev.off()

hist(count_matrix$nc_01)
hist(count_matrix$nc_02)
hist(count_matrix$kd_01)
hist(count_matrix$kd_02)
hist(count_matrix$nc_01_input)
hist(count_matrix$nc_02_input)
hist(count_matrix$kd_01_input)
hist(count_matrix$kd_02_input)

count_matrix_1 <- count_matrix[(count_matrix$nc_01!="0") & 
                                 (count_matrix$nc_02!="0") & 
                                 (count_matrix$kd_01!="0") & 
                                 (count_matrix$kd_02!="0") & 
                                 (count_matrix$nc_01_input!="0") & 
                                 (count_matrix$nc_02_input!="0") & 
                                 (count_matrix$kd_01_input!="0") & 
                                 (count_matrix$kd_02_input!="0"), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01>"200") & 
                                 (count_matrix$nc_02>"200") & 
                                 (count_matrix$kd_01>"200") & 
                                 (count_matrix$kd_02>"200") & 
                                 (count_matrix$nc_01_input>"200") & 
                                 (count_matrix$nc_02_input>"200") & 
                                 (count_matrix$kd_01_input>"200") & 
                                 (count_matrix$kd_02_input>"200"), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())

### MAPKAPK5-AS1_H3K4me3

count_matrix <- read.table("~/count_matrix/mapkapk_h3k4me3", header = TRUE, sep = "\t")

png(filename="~/scatter/mapkapk_h3k4me3_nc_01.png")
plot(count_matrix$nc_01)
dev.off()

png(filename="~/scatter/mapkapk_h3k4me3_nc_02.png")
plot(count_matrix$nc_02)
dev.off()

png(filename="~/scatter/mapkapk_h3k4me3_kd_01.png")
plot(count_matrix$kd_01)
dev.off()

png(filename="~/scatter/mapkapk_h3k4me3_kd_02.png")
plot(count_matrix$kd_02)
dev.off()

png(filename="~/scatter/mapkapk_h3k4me3_nc_01_input.png")
plot(count_matrix$nc_01_input)
dev.off()

png(filename="~/scatter/mapkapk_h3k4me3_nc_02_input.png")
plot(count_matrix$nc_02_input)
dev.off()

png(filename="~/scatter/mapkapk_h3k4me3_kd_01_input.png")
plot(count_matrix$kd_01_input)
dev.off()

png(filename="~/scatter/mapkapk_h3k4me3_kd_02_input.png")
plot(count_matrix$kd_02_input)
dev.off()

hist(count_matrix$nc_01)
hist(count_matrix$nc_02)
hist(count_matrix$kd_01)
hist(count_matrix$kd_02)
hist(count_matrix$nc_01_input)
hist(count_matrix$nc_02_input)
hist(count_matrix$kd_01_input)
hist(count_matrix$kd_02_input)

count_matrix_1 <- count_matrix[(count_matrix$nc_01!="0") & 
                                 (count_matrix$nc_02!="0") & 
                                 (count_matrix$kd_01!="0") & 
                                 (count_matrix$kd_02!="0") & 
                                 (count_matrix$nc_01_input!="0") & 
                                 (count_matrix$nc_02_input!="0") & 
                                 (count_matrix$kd_01_input!="0") & 
                                 (count_matrix$kd_02_input!="0"), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01>"200") & 
                                 (count_matrix$nc_02>"200") & 
                                 (count_matrix$kd_01>"200") & 
                                 (count_matrix$kd_02>"200") & 
                                 (count_matrix$nc_01_input>"200") & 
                                 (count_matrix$nc_02_input>"200") & 
                                 (count_matrix$kd_01_input>"200") & 
                                 (count_matrix$kd_02_input>"200"), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())
