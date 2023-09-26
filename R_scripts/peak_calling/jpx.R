library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

BSgenome='BSgenome.Hsapiens.UCSC.hg38'

uniq=0
extend=300
shift=0
ws=1000

kd_01 = MEDIPS.createSet(file = '~/Downloads/jpx/bam/C-11.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02 = MEDIPS.createSet(file = '~/Downloads/jpx/bam/C-12.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
KD = c(kd_01,kd_02)

nc_01 = MEDIPS.createSet(file = '~/Downloads/jpx/bam/C-17.bam', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02 = MEDIPS.createSet(file = '~/Downloads/jpx/bam/C-18.bam' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
NC = c(nc_01,nc_02)

kd_01_input = MEDIPS.createSet(file = '~/Downloads/jpx/bam/I-11.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_02_input = MEDIPS.createSet(file = '~/Downloads/jpx/bam/I-12.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
kd_input = c(kd_01_Input, kd_02_Input)

nc_01_input = MEDIPS.createSet(file = '~/Downloads/jpx/bam/I-17.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_02_input = MEDIPS.createSet(file = '~/Downloads/jpx/bam/I-18.bam', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
nc_input = c(nc_01_Input,nc_02_Input)

CS = MEDIPS.couplingVector(pattern = 'CG', refObj = NC[[1]])
CS_01 = MEDIPS.couplingVector(pattern = 'CG', refObj = nc_01[[1]])
CS_02 = MEDIPS.couplingVector(pattern = 'CG', refObj = nc_02[[1]])

### FOR SCALE FACTOR

resultTable_nc_01 = MEDIPS.meth(MSet1 = nc_01, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_nc_02 = MEDIPS.meth(MSet1 = nc_02, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_kd_01 = MEDIPS.meth(MSet1 = kd_01, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_kd_02 = MEDIPS.meth(MSet1 = kd_02, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_nc_01_input = MEDIPS.meth(MSet1 = nc_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_nc_02_input = MEDIPS.meth(MSet1 = nc_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_kd_01_input = MEDIPS.meth(MSet1 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_kd_02_input = MEDIPS.meth(MSet1 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

### FOR SCALE FACTOR AGAINST INPUT

resultTable_nc_input_01 = MEDIPS.meth(MSet1 = nc_01, ISet1 = nc_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_nc_input_02 = MEDIPS.meth(MSet1 = nc_02, ISet1 = nc_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_kd_input_01 = MEDIPS.meth(MSet1 = kd_01, ISet1 = kd_01_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

resultTable_kd_input_02 = MEDIPS.meth(MSet1 = kd_02, ISet1 = kd_02_input, p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')
options(scipen = 999)

### WRITE TABLE

write.table(file='~/Downloads/jpx/peaks/medips_jpx_kd_01.txt', resultTable_kd_01, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_kd_02.txt', resultTable_kd_02, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_nc_01.txt', resultTable_nc_01, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_nc_02.txt', resultTable_nc_02, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_kd_01_input.txt', resultTable_kd_01_input, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_kd_02_input.txt', resultTable_kd_02_input, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_nc_01_input.txt', resultTable_nc_01_input, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_nc_02_input.txt', resultTable_nc_02_input, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(file='~/Downloads/jpx/peaks/medips_jpx_kd_input_01.txt', resultTable_kd_input_01, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_kd_input_02.txt', resultTable_kd_input_02, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_nc_input_01.txt', resultTable_nc_input_01, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(file='~/Downloads/jpx/peaks/medips_jpx_nc_input_02.txt', resultTable_nc_input_02, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
