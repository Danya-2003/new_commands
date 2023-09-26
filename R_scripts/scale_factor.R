library("DESeq2")
library("edgeR")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

### FGD5-AS1

count_matrix <- read.table("~/count_matrix/fgd5", header = TRUE, sep = "\t")

count_matrix_1 <- count_matrix[(count_matrix$nc_01!= 0) & 
                                 (count_matrix$nc_02!= 0) & 
                                 (count_matrix$kd_01!= 0) & 
                                 (count_matrix$kd_02!= 0) & 
                                 (count_matrix$nc_01_input!= 0) & 
                                 (count_matrix$nc_02_input!= 0) & 
                                 (count_matrix$kd_01_input!= 0) & 
                                 (count_matrix$kd_02_input!= 0), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01<162) & 
                                 (count_matrix$nc_02<162) & 
                                 (count_matrix$kd_01<162) & 
                                 (count_matrix$kd_02<162) & 
                                 (count_matrix$nc_01_input<162) & 
                                 (count_matrix$nc_02_input<162) & 
                                 (count_matrix$kd_01_input<162) & 
                                 (count_matrix$kd_02_input<162), ]

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

count_matrix_1 <- count_matrix[(count_matrix$nc_01!= 0) & 
                                 (count_matrix$nc_02!= 0) & 
                                 (count_matrix$kd_01!= 0) & 
                                 (count_matrix$kd_02!= 0) & 
                                 (count_matrix$nc_01_input!= 0) & 
                                 (count_matrix$nc_02_input!= 0) & 
                                 (count_matrix$kd_01_input!= 0) & 
                                 (count_matrix$kd_02_input!= 0), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01<166) & 
                                 (count_matrix$nc_02<166) & 
                                 (count_matrix$kd_01<166) & 
                                 (count_matrix$kd_02<166) & 
                                 (count_matrix$nc_01_input<166) & 
                                 (count_matrix$nc_02_input<166) & 
                                 (count_matrix$kd_01_input<166) & 
                                 (count_matrix$kd_02_input<166), ]

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

count_matrix_1 <- count_matrix[(count_matrix$nc_01!= 0) & 
                                 (count_matrix$nc_02!= 0) & 
                                 (count_matrix$kd_01!= 0) & 
                                 (count_matrix$kd_02!= 0) & 
                                 (count_matrix$nc_01_input!= 0) & 
                                 (count_matrix$nc_02_input!= 0) & 
                                 (count_matrix$kd_01_input!= 0) & 
                                 (count_matrix$kd_02_input!= 0), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01<127) & 
                                 (count_matrix$nc_02<127) & 
                                 (count_matrix$kd_01<127) & 
                                 (count_matrix$kd_02<127) & 
                                 (count_matrix$nc_01_input<127) & 
                                 (count_matrix$nc_02_input<127) & 
                                 (count_matrix$kd_01_input<127) & 
                                 (count_matrix$kd_02_input<127), ]

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

count_matrix_1 <- count_matrix[(count_matrix$nc_01!= 0) & 
                                 (count_matrix$nc_02!= 0) & 
                                 (count_matrix$kd_01!= 0) & 
                                 (count_matrix$kd_02!= 0) & 
                                 (count_matrix$nc_01_input!= 0) & 
                                 (count_matrix$nc_02_input!= 0) & 
                                 (count_matrix$kd_01_input!= 0) & 
                                 (count_matrix$kd_02_input!= 0), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01<67) & 
                                 (count_matrix$nc_02<67) & 
                                 (count_matrix$kd_01<67) & 
                                 (count_matrix$kd_02<67) & 
                                 (count_matrix$nc_01_input<67) & 
                                 (count_matrix$nc_02_input<67) & 
                                 (count_matrix$kd_01_input<67) & 
                                 (count_matrix$kd_02_input<67), ]

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

count_matrix_1 <- count_matrix[(count_matrix$nc_01!= 0) & 
                                 (count_matrix$nc_02!= 0) & 
                                 (count_matrix$kd_01!= 0) & 
                                 (count_matrix$kd_02!= 0) & 
                                 (count_matrix$nc_01_input!= 0) & 
                                 (count_matrix$nc_02_input!= 0) & 
                                 (count_matrix$kd_01_input!= 0) & 
                                 (count_matrix$kd_02_input!= 0), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01<46) & 
                                 (count_matrix$nc_02<46) & 
                                 (count_matrix$kd_01<46) & 
                                 (count_matrix$kd_02<46) & 
                                 (count_matrix$nc_01_input<46) & 
                                 (count_matrix$nc_02_input<46) & 
                                 (count_matrix$kd_01_input<46) & 
                                 (count_matrix$kd_02_input<46), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())

### JPX

count_matrix <- read.table("~/count_matrix/jpx", header = TRUE, sep = "\t")


count_matrix_0 <- count_matrix[(count_matrix$nc_01!= 0) & 
                                 (count_matrix$nc_02!= 0) & 
                                 (count_matrix$kd_01!= 0) & 
                                 (count_matrix$kd_02!= 0) & 
                                 (count_matrix$nc_01_input!= 0) & 
                                 (count_matrix$nc_02_input!= 0) & 
                                 (count_matrix$kd_01_input!= 0) & 
                                 (count_matrix$kd_02_input!= 0), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01<53) & 
                                 (count_matrix$nc_02<53) & 
                                 (count_matrix$kd_01<53) & 
                                 (count_matrix$kd_02<53) & 
                                 (count_matrix$nc_01_input<53) & 
                                 (count_matrix$nc_02_input<53) & 
                                 (count_matrix$kd_01_input<53) & 
                                 (count_matrix$kd_02_input<53), ]

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

count_matrix_1 <- count_matrix[(count_matrix$nc_01!= 0) & 
                                 (count_matrix$nc_02!= 0) & 
                                 (count_matrix$kd_01!= 0) & 
                                 (count_matrix$kd_02!= 0) & 
                                 (count_matrix$nc_01_input!= 0) & 
                                 (count_matrix$nc_02_input!= 0) & 
                                 (count_matrix$kd_01_input!= 0) & 
                                 (count_matrix$kd_02_input!= 0), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01<227) & 
                                 (count_matrix$nc_02<227) & 
                                 (count_matrix$kd_01<227) & 
                                 (count_matrix$kd_02<227) & 
                                 (count_matrix$nc_01_input<227) & 
                                 (count_matrix$nc_02_input<227) & 
                                 (count_matrix$kd_01_input<227) & 
                                 (count_matrix$kd_02_input<227), ]

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


count_matrix_1 <- count_matrix[(count_matrix$nc_01!= 0) & 
                                 (count_matrix$nc_02!= 0) & 
                                 (count_matrix$kd_01!= 0) & 
                                 (count_matrix$kd_02!= 0) & 
                                 (count_matrix$nc_01_input!= 0) & 
                                 (count_matrix$nc_02_input!= 0) & 
                                 (count_matrix$kd_01_input!= 0) & 
                                 (count_matrix$kd_02_input!= 0), ]

count_matrix_1 <- count_matrix[(count_matrix$nc_01<53) & 
                                 (count_matrix$nc_02<53) & 
                                 (count_matrix$kd_01<53) & 
                                 (count_matrix$kd_02<53) & 
                                 (count_matrix$nc_01_input<53) & 
                                 (count_matrix$nc_02_input<53) & 
                                 (count_matrix$kd_01_input<53) & 
                                 (count_matrix$kd_02_input<53), ]

NormFactor <- calcNormFactors(object = count_matrix_1, method = "TMM")
LibSize <- colSums(count_matrix_1)
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.Reciprocal <- 1/SizeFactors

col.sum <- colSums(count_matrix_1)
check <- as.data.frame(cbind(col.sum, SizeFactors.Reciprocal))
check$ration <- check$SizeFactors.Reciprocal/check$col.sum

rm(list = ls())

