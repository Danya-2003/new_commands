library("DESeq2")
library("edgeR")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

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

rm(list = ls())
