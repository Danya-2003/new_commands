fgd5 <- read.table("~/fgd5", header = TRUE, sep = "\t")
rownames(fgd5) <- fgd5$FGD5.AS1_H3K9ac
fgd5_box <- fgd5[,c(2,3,4,5)]
boxplot(fgd5_box)

emx2os <- read.table("~/emx2os", header = TRUE, sep = "\t")
rownames(emx2os) <- emx2os$EMX2OS_H3K9ac
emx2os_box <- emx2os[,c(2,3,4,5)]
boxplot(emx2os_box)

rp11 <- read.table("~/rp11", header = TRUE, sep = "\t")
rownames(rp11) <- rp11$RP11.398K22.12_H3K27ac
rp11_box <- rp11[,c(4,5)]
boxplot(rp11_box)

ctd <- read.table("~/ctd", header = TRUE, sep = "\t")
rownames(ctd) <- ctd$CTD.2587H24.5_H3K36me3
ctd_box <- ctd[,c(4,5)]
boxplot(ctd_box)

jpx <- read.table("~/jpx", header = TRUE, sep = "\t")
rownames(jpx) <- jpx$JPX_H3K9me3
jpx_box <- jpx[,c(3,4,5)]
boxplot(jpx_box)

ac <- read.table("~/ac", header = TRUE, sep = "\t")
rownames(ac) <- ac$AC005592.2_H3K9me3
ac_box <- ac[,c(4,5)]
boxplot(ac_box)

mapkapk_h3k27 <- read.table("~/mapkapk_h3k27", header = TRUE, sep = "\t")
rownames(mapkapk_h3k27) <- mapkapk_h3k27$MAPKAPK5.AS1_H3K27me3
mapkapk_h3k27_box <- mapkapk_h3k27[,c(4,5)]
boxplot(mapkapk_h3k27_box)

mapkapk_h3k4 <- read.table("~/mapkapk_h3k4", header = TRUE, sep = "\t")
rownames(mapkapk_h3k4) <- mapkapk_h3k4$MAPKAPK5.AS1_H3K4me3
mapkapk_h3k4_box <- mapkapk_h3k4[,c(3,4,5)]
boxplot(mapkapk_h3k4_box)
