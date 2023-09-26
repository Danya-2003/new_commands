library("DESeq2")
library("edgeR")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)



### FGD5-AS1

count_matrix <- read.table("~/count_matrix/fgd5", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01<162) & 
                                   (count_matrix$nc_02<162) &
                                   (count_matrix$kd_01<162) &
                                   (count_matrix$kd_02<162) &
                                   (count_matrix$nc_01_input<162) &
                                   (count_matrix$nc_02_input<162) &
                                   (count_matrix$kd_01_input<162) &
                                   (count_matrix$kd_02_input<162), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, regex('[:digit:]_[:digit:]', type = line_break), 3)
write.table(count_matrix, '~/filter_bam/fgd5_162.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

count_matrix <- read.table("~/count_matrix/fgd5", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01!=0) & 
                               (count_matrix$nc_02!=0) & 
                               (count_matrix$kd_01!=0) & 
                               (count_matrix$kd_02!=0) & 
                               (count_matrix$nc_01_input!=0) & 
                               (count_matrix$nc_02_input!=0) & 
                               (count_matrix$kd_01_input!=0) & 
                               (count_matrix$kd_02_input!=0), ]

count_matrix <- count_matrix[(count_matrix$nc_01<162) & 
                               (count_matrix$nc_02<162) &
                               (count_matrix$kd_01<162) &
                               (count_matrix$kd_02<162) &
                               (count_matrix$nc_01_input<162) &
                               (count_matrix$nc_02_input<162) &
                               (count_matrix$kd_01_input<162) &
                               (count_matrix$kd_02_input<162), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/fgd5_0_162.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')


rm(list = ls())


### EMX2OS

count_matrix <- read.table("~/count_matrix/emx2os", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01!=0) & 
                                 (count_matrix$nc_02!=0) & 
                                 (count_matrix$kd_01!=0) & 
                                 (count_matrix$kd_02!=0) & 
                                 (count_matrix$nc_01_input!=0) & 
                                 (count_matrix$nc_02_input!=0) & 
                                 (count_matrix$kd_01_input!=0) & 
                                 (count_matrix$kd_02_input!=0), ]

count_matrix <- count_matrix[(count_matrix$nc_01<166) & 
                                   (count_matrix$nc_02<166) & 
                                   (count_matrix$kd_01<166) & 
                                   (count_matrix$kd_02<166) & 
                                   (count_matrix$nc_01_input<166) & 
                                   (count_matrix$nc_02_input<166) & 
                                   (count_matrix$kd_01_input<166) & 
                                   (count_matrix$kd_02_input<166), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/emx2os_0_166.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

count_matrix <- read.table("~/count_matrix/emx2os", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01<166) & 
                               (count_matrix$nc_02<166) & 
                               (count_matrix$kd_01<166) & 
                               (count_matrix$kd_02<166) & 
                               (count_matrix$nc_01_input<166) & 
                               (count_matrix$nc_02_input<166) & 
                               (count_matrix$kd_01_input<166) & 
                               (count_matrix$kd_02_input<166), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/emx2os_166.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

rm(list = ls())


### AC005592.2

count_matrix <- read.table("~/count_matrix/ac", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01!=0) & 
                                 (count_matrix$nc_02!=0) & 
                                 (count_matrix$kd_01!=0) & 
                                 (count_matrix$kd_02!=0) & 
                                 (count_matrix$nc_01_input!=0) & 
                                 (count_matrix$nc_02_input!=0) & 
                                 (count_matrix$kd_01_input!=0) & 
                                 (count_matrix$kd_02_input!=0), ]

count_matrix <- count_matrix[(count_matrix$nc_01<127) & 
                                   (count_matrix$nc_02<127) & 
                                   (count_matrix$kd_01<127) & 
                                   (count_matrix$kd_02<127) & 
                                   (count_matrix$nc_01_input<127) & 
                                   (count_matrix$nc_02_input<127) & 
                                   (count_matrix$kd_01_input<127) & 
                                   (count_matrix$kd_02_input<127), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/ac_0_127.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

count_matrix <- read.table("~/count_matrix/ac", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01<127) & 
                               (count_matrix$nc_02<127) & 
                               (count_matrix$kd_01<127) & 
                               (count_matrix$kd_02<127) & 
                               (count_matrix$nc_01_input<127) & 
                               (count_matrix$nc_02_input<127) & 
                               (count_matrix$kd_01_input<127) & 
                               (count_matrix$kd_02_input<127), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/ac_127.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

rm(list = ls())


### RP11-398K22.12

count_matrix <- read.table("~/count_matrix/rp11", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01!=0) & 
                                 (count_matrix$nc_02!=0) & 
                                 (count_matrix$kd_01!=0) & 
                                 (count_matrix$kd_02!=0) & 
                                 (count_matrix$nc_01_input!=0) & 
                                 (count_matrix$nc_02_input!=0) & 
                                 (count_matrix$kd_01_input!=0) & 
                                 (count_matrix$kd_02_input!=0), ]


count_matrix <- count_matrix[(count_matrix$nc_01<67) & 
                                   (count_matrix$nc_02<67) & 
                                   (count_matrix$kd_01<67) & 
                                   (count_matrix$kd_02<67) & 
                                   (count_matrix$nc_01_input<67) & 
                                   (count_matrix$nc_02_input<67) & 
                                   (count_matrix$kd_01_input<67) & 
                                   (count_matrix$kd_02_input<67), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/rp11_0_67.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

count_matrix <- read.table("~/count_matrix/rp11", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01<67) & 
                               (count_matrix$nc_02<67) & 
                               (count_matrix$kd_01<67) & 
                               (count_matrix$kd_02<67) & 
                               (count_matrix$nc_01_input<67) & 
                               (count_matrix$nc_02_input<67) & 
                               (count_matrix$kd_01_input<67) & 
                               (count_matrix$kd_02_input<67), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/rp11_67.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

rm(list = ls())


### CTD


count_matrix <- read.table("~/count_matrix/ctd", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01!=0) & 
                                 (count_matrix$nc_02!=0) & 
                                 (count_matrix$kd_01!=0) & 
                                 (count_matrix$kd_02!=0) & 
                                 (count_matrix$nc_01_input!=0) & 
                                 (count_matrix$nc_02_input!=0) & 
                                 (count_matrix$kd_01_input!=0) & 
                                 (count_matrix$kd_02_input!=0), ]


count_matrix <- count_matrix[(count_matrix$nc_01<46) & 
                                   (count_matrix$nc_02<46) & 
                                   (count_matrix$kd_01<46) & 
                                   (count_matrix$kd_02<46) & 
                                   (count_matrix$nc_01_input<46) & 
                                   (count_matrix$nc_02_input<46) & 
                                   (count_matrix$kd_01_input<46) & 
                                   (count_matrix$kd_02_input<46), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/ctd_0_46.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

count_matrix <- read.table("~/count_matrix/ctd", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01<46) & 
                               (count_matrix$nc_02<46) & 
                               (count_matrix$kd_01<46) & 
                               (count_matrix$kd_02<46) & 
                               (count_matrix$nc_01_input<46) & 
                               (count_matrix$nc_02_input<46) & 
                               (count_matrix$kd_01_input<46) & 
                               (count_matrix$kd_02_input<46), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/ctd_46.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

rm(list = ls())

### JPX

count_matrix <- read.table("~/count_matrix/jpx", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01!=0) & 
                                 (count_matrix$nc_02!=0) & 
                                 (count_matrix$kd_01!=0) & 
                                 (count_matrix$kd_02!=0) & 
                                 (count_matrix$nc_01_input!=0) & 
                                 (count_matrix$nc_02_input!=0) & 
                                 (count_matrix$kd_01_input!=0) & 
                                 (count_matrix$kd_02_input!=0), ]

count_matrix <- count_matrix[(count_matrix$nc_01<53) & 
                                   (count_matrix$nc_02<53) & 
                                   (count_matrix$kd_01<53) & 
                                   (count_matrix$kd_02<53) & 
                                   (count_matrix$nc_01_input<53) & 
                                   (count_matrix$nc_02_input<53) & 
                                   (count_matrix$kd_01_input<53) & 
                                   (count_matrix$kd_02_input<53), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/jpx_0_53.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

count_matrix <- read.table("~/count_matrix/jpx", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01<53) & 
                               (count_matrix$nc_02<53) & 
                               (count_matrix$kd_01<53) & 
                               (count_matrix$kd_02<53) & 
                               (count_matrix$nc_01_input<53) & 
                               (count_matrix$nc_02_input<53) & 
                               (count_matrix$kd_01_input<53) & 
                               (count_matrix$kd_02_input<53), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/jpx_53.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

rm(list = ls())


### MAPKAPK5-AS1_H3K27me3

count_matrix <- read.table("~/count_matrix/mapkapk_h3k27me3", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01!=0) & 
                                 (count_matrix$nc_02!=0) & 
                                 (count_matrix$kd_01!=0) & 
                                 (count_matrix$kd_02!=0) & 
                                 (count_matrix$nc_01_input!=0) & 
                                 (count_matrix$nc_02_input!=0) & 
                                 (count_matrix$kd_01_input!=0) & 
                                 (count_matrix$kd_02_input!=0), ]

count_matrix <- count_matrix[(count_matrix$nc_01<227) & 
                                   (count_matrix$nc_02<227) & 
                                   (count_matrix$kd_01<227) & 
                                   (count_matrix$kd_02<227) & 
                                   (count_matrix$nc_01_input<227) & 
                                   (count_matrix$nc_02_input<227) & 
                                   (count_matrix$kd_01_input<227) & 
                                   (count_matrix$kd_02_input<227), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/mapkapk5_h3k27_0_227.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

count_matrix <- read.table("~/count_matrix/mapkapk_h3k27me3", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01<227) & 
                               (count_matrix$nc_02<227) & 
                               (count_matrix$kd_01<227) & 
                               (count_matrix$kd_02<227) & 
                               (count_matrix$nc_01_input<227) & 
                               (count_matrix$nc_02_input<227) & 
                               (count_matrix$kd_01_input<227) & 
                               (count_matrix$kd_02_input<227), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/mapkapk5_h3k27_227.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

rm(list = ls())

### MAPKAPK5-AS1_H3K4me3

count_matrix <- read.table("~/count_matrix/mapkapk_h3k4me3", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01!=0) & 
                               (count_matrix$nc_02!=0) & 
                               (count_matrix$kd_01!=0) & 
                               (count_matrix$kd_02!=0) & 
                               (count_matrix$nc_01_input!=0) & 
                               (count_matrix$nc_02_input!=0) & 
                               (count_matrix$kd_01_input!=0) & 
                               (count_matrix$kd_02_input!=0), ]

count_matrix <- count_matrix[(count_matrix$nc_01<53) & 
                               (count_matrix$nc_02<53) & 
                               (count_matrix$kd_01<53) & 
                               (count_matrix$kd_02<53) & 
                               (count_matrix$nc_01_input<53) & 
                               (count_matrix$nc_02_input<53) & 
                               (count_matrix$kd_01_input<53) & 
                               (count_matrix$kd_02_input<53), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/mapkapk5_h3k4_0_53.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

count_matrix <- read.table("~/count_matrix/mapkapk_h3k4me3", header = TRUE, sep = "\t")

count_matrix <- count_matrix[(count_matrix$nc_01<53) & 
                               (count_matrix$nc_02<53) & 
                               (count_matrix$kd_01<53) & 
                               (count_matrix$kd_02<53) & 
                               (count_matrix$nc_01_input<53) & 
                               (count_matrix$nc_02_input<53) & 
                               (count_matrix$kd_01_input<53) & 
                               (count_matrix$kd_02_input<53), ]

count_matrix$bed <- row.names(count_matrix)
count_matrix <- str_split_fixed(count_matrix$bed, "_", 3)
write.table(count_matrix, '~/filter_bam/mapkapk5_h3k4_53.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

rm(list = ls())

