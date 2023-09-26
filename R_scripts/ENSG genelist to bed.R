library(biomaRt)

biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
t2g<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position', 'strand'), mart = ensembl)


ac_high <- read.table("~/Desktop/commands/py/genes/ac_high.txt", header =FALSE, sep = '\t')
ac_low <- read.table("~/Desktop/commands/py/genes/ac_low.txt", header =FALSE, sep = '\t')
ac_zero <- read.table("~/Desktop/commands/py/genes/ac_zero.txt", header =FALSE, sep = '\t')

fgd5_high <- read.table("~/Desktop/commands/py/genes/fgd5_high.txt", header =FALSE, sep = '\t')
fgd5_low <- read.table("~/Desktop/commands/py/genes/fgd5_low.txt", header =FALSE, sep = '\t')
fgd5_zero <- read.table("~/Desktop/commands/py/genes/fgd5_zero.txt", header =FALSE, sep = '\t')

rp11_high <- read.table("~/Desktop/commands/py/genes/rp11_high.txt", header =FALSE, sep = '\t')
rp11_low <- read.table("~/Desktop/commands/py/genes/rp11_low.txt", header =FALSE, sep = '\t')
rp11_zero <- read.table("~/Desktop/commands/py/genes/rp11_zero.txt", header =FALSE, sep = '\t')

mapkapk_high <- read.table("~/Desktop/commands/py/genes/mapkapk_high.txt", header =FALSE, sep = '\t')
mapkapk_low <- read.table("~/Desktop/commands/py/genes/mapkapk_low.txt", header =FALSE, sep = '\t')
mapkapk_zero <- read.table("~/Desktop/commands/py/genes/mapkapk_zero.txt", header =FALSE, sep = '\t')

ctd_high <- read.table("~/Desktop/commands/py/genes/ctd_high.txt", header =FALSE, sep = '\t')
ctd_low <- read.table("~/Desktop/commands/py/genes/ctd_low.txt", header =FALSE, sep = '\t')
ctd_zero <- read.table("~/Desktop/commands/py/genes/ctd_zero.txt", header =FALSE, sep = '\t')

jpx_high <- read.table("~/Desktop/commands/py/genes/jpx_high.txt", header =FALSE, sep = '\t')
jpx_low <- read.table("~/Desktop/commands/py/genes/jpx_low.txt", header =FALSE, sep = '\t')
jpx_zero <- read.table("~/Desktop/commands/py/genes/jpx_zero.txt", header =FALSE, sep = '\t')

emx2os_high <- read.table("~/Desktop/commands/py/genes/emx2os_high.txt", header =FALSE, sep = '\t')
emx2os_low <- read.table("~/Desktop/commands/py/genes/emx2os_low.txt", header =FALSE, sep = '\t')
emx2os_zero <- read.table("~/Desktop/commands/py/genes/emx2os_zero.txt", header =FALSE, sep = '\t')

nc_high <- read.table("~/Desktop/commands/py/genes/nc_high.txt", header =FALSE, sep = '\t')
nc_low <- read.table("~/Desktop/commands/py/genes/nc_low.txt", header =FALSE, sep = '\t')
nc_zero <- read.table("~/Desktop/commands/py/genes/nc_zero.txt", header =FALSE, sep = '\t')


colnames(ac_high) <- c("ensembl_gene_id")
colnames(ac_low) <- c("ensembl_gene_id")
colnames(ac_zero) <- c("ensembl_gene_id")
colnames(fgd5_high) <- c("ensembl_gene_id")
colnames(fgd5_low) <- c("ensembl_gene_id")
colnames(fgd5_zero) <- c("ensembl_gene_id")
colnames(rp11_high) <- c("ensembl_gene_id")
colnames(rp11_low) <- c("ensembl_gene_id")
colnames(rp11_zero) <- c("ensembl_gene_id")
colnames(mapkapk_high) <- c("ensembl_gene_id")
colnames(mapkapk_low) <- c("ensembl_gene_id")
colnames(mapkapk_zero) <- c("ensembl_gene_id")
colnames(ctd_high) <- c("ensembl_gene_id")
colnames(ctd_low) <- c("ensembl_gene_id")
colnames(ctd_zero) <- c("ensembl_gene_id")
colnames(jpx_high) <- c("ensembl_gene_id")
colnames(jpx_low) <- c("ensembl_gene_id")
colnames(jpx_zero) <- c("ensembl_gene_id")
colnames(emx2os_high) <- c("ensembl_gene_id")
colnames(emx2os_low) <- c("ensembl_gene_id")
colnames(emx2os_zero) <- c("ensembl_gene_id")
colnames(nc_high) <- c("ensembl_gene_id")
colnames(nc_low) <- c("ensembl_gene_id")
colnames(nc_zero) <- c("ensembl_gene_id")

id_ac_high <- merge(x = ac_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ac_low <- merge(x = ac_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ac_zero <- merge(x = ac_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_fgd5_high <- merge(x = fgd5_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_fgd5_low <- merge(x = fgd5_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_fgd5_zero <- merge(x = fgd5_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_rp11_high <- merge(x = rp11_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_rp11_low <- merge(x = rp11_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_rp11_zero <- merge(x = rp11_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_mapkapk_high <- merge(x = mapkapk_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_mapkapk_low <- merge(x = mapkapk_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_mapkapk_zero <- merge(x = mapkapk_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ctd_high <- merge(x = ctd_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ctd_low <- merge(x = ctd_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ctd_zero <- merge(x = ctd_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_jpx_high <- merge(x = jpx_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_jpx_low <- merge(x = jpx_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_jpx_zero <- merge(x = jpx_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_emx2os_high <- merge(x = emx2os_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_emx2os_low <- merge(x = emx2os_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_emx2os_zero <- merge(x = emx2os_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_nc_high <- merge(x = nc_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_nc_low <- merge(x = nc_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_nc_zero <- merge(x = nc_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)

id_ac_high <- na.omit(id_ac_high)
id_ac_low <- na.omit(id_ac_low)
id_ac_zero <- na.omit(id_ac_zero)
id_fgd5_high <- na.omit(id_fgd5_high)
id_fgd5_low <- na.omit(id_fgd5_low)
id_fgd5_zero <- na.omit(id_fgd5_zero)
id_rp11_high <- na.omit(id_rp11_high)
id_rp11_low <- na.omit(id_rp11_low)
id_rp11_zero <- na.omit(id_rp11_zero)
id_mapkapk_high <- na.omit(id_mapkapk_high)
id_mapkapk_low <- na.omit(id_mapkapk_low)
id_mapkapk_zero <- na.omit(id_mapkapk_zero)
id_ctd_high <- na.omit(id_ctd_high)
id_ctd_low <- na.omit(id_ctd_low)
id_ctd_zero <- na.omit(id_ctd_zero)
id_jpx_high <- na.omit(id_jpx_high)
id_jpx_low <- na.omit(id_jpx_low)
id_jpx_zero <- na.omit(id_jpx_zero)
id_emx2os_high <- na.omit(id_emx2os_high)
id_emx2os_low <- na.omit(id_emx2os_low)
id_emx2os_zero <- na.omit(id_emx2os_zero)
id_nc_high <- na.omit(id_nc_high)
id_nc_low <- na.omit(id_nc_low)
id_nc_zero <- na.omit(id_nc_zero)

id_ac_high$chromosome_name <- paste0(x = "chr", id_ac_high$chromosome_name)
id_ac_low$chromosome_name <- paste0(x = "chr", id_ac_low$chromosome_name)
id_ac_zero$chromosome_name <- paste0(x = "chr", id_ac_zero$chromosome_name)
id_fgd5_high$chromosome_name <- paste0(x = "chr", id_fgd5_high$chromosome_name)
id_fgd5_low$chromosome_name <- paste0(x = "chr", id_fgd5_low$chromosome_name)
id_fgd5_zero$chromosome_name <- paste0(x = "chr", id_fgd5_zero$chromosome_name)
id_rp11_high$chromosome_name <- paste0(x = "chr", id_rp11_high$chromosome_name)
id_rp11_low$chromosome_name <- paste0(x = "chr", id_rp11_low$chromosome_name)
id_rp11_zero$chromosome_name <- paste0(x = "chr", id_rp11_zero$chromosome_name)
id_mapkapk_high$chromosome_name <- paste0(x = "chr", id_mapkapk_high$chromosome_name)
id_mapkapk_low$chromosome_name <- paste0(x = "chr", id_mapkapk_low$chromosome_name)
id_mapkapk_zero$chromosome_name <- paste0(x = "chr", id_mapkapk_zero$chromosome_name)
id_ctd_high$chromosome_name <- paste0(x = "chr", id_ctd_high$chromosome_name)
id_ctd_low$chromosome_name <- paste0(x = "chr", id_ctd_low$chromosome_name)
id_ctd_zero$chromosome_name <- paste0(x = "chr", id_ctd_zero$chromosome_name)
id_jpx_high$chromosome_name <- paste0(x = "chr", id_jpx_high$chromosome_name)
id_jpx_low$chromosome_name <- paste0(x = "chr", id_jpx_low$chromosome_name)
id_jpx_zero$chromosome_name <- paste0(x = "chr", id_jpx_zero$chromosome_name)
id_emx2os_high$chromosome_name <- paste0(x = "chr", id_emx2os_high$chromosome_name)
id_emx2os_low$chromosome_name <- paste0(x = "chr", id_emx2os_low$chromosome_name)
id_emx2os_zero$chromosome_name <- paste0(x = "chr", id_emx2os_zero$chromosome_name)
id_nc_high$chromosome_name <- paste0(x = "chr", id_nc_high$chromosome_name)
id_nc_low$chromosome_name <- paste0(x = "chr", id_nc_low$chromosome_name)
id_nc_zero$chromosome_name <- paste0(x = "chr", id_nc_zero$chromosome_name)

id_ac_high<-id_ac_high[,c(2,3,4,5)]
id_ac_low<-id_ac_low[,c(2,3,4,5)]
id_ac_zero<-id_ac_zero[,c(2,3,4,5)]
id_fgd5_high<-id_fgd5_high[,c(2,3,4,5)]
id_fgd5_low<-id_fgd5_low[,c(2,3,4,5)]
id_fgd5_zero<-id_fgd5_zero[,c(2,3,4,5)]
id_rp11_high<-id_rp11_high[,c(2,3,4,5)]
id_rp11_low<-id_rp11_low[,c(2,3,4,5)]
id_rp11_zero<-id_rp11_zero[,c(2,3,4,5)]
id_mapkapk_high<-id_mapkapk_high[,c(2,3,4,5)]
id_mapkapk_low<-id_mapkapk_low[,c(2,3,4,5)]
id_mapkapk_zero<-id_mapkapk_zero[,c(2,3,4,5)]
id_ctd_high<-id_ctd_high[,c(2,3,4,5)]
id_ctd_low<-id_ctd_low[,c(2,3,4,5)]
id_ctd_zero<-id_ctd_zero[,c(2,3,4,5)]
id_jpx_high<-id_jpx_high[,c(2,3,4,5)]
id_jpx_low<-id_jpx_low[,c(2,3,4,5)]
id_jpx_zero<-id_jpx_zero[,c(2,3,4,5)]
id_emx2os_high<-id_emx2os_high[,c(2,3,4,5)]
id_emx2os_low<-id_emx2os_low[,c(2,3,4,5)]
id_emx2os_zero<-id_emx2os_zero[,c(2,3,4,5)]
id_nc_high<-id_nc_high[,c(2,3,4,5)]
id_nc_low<-id_nc_low[,c(2,3,4,5)]
id_nc_zero<-id_nc_zero[,c(2,3,4,5)]

id_ac_high<-id_ac_high[order(id_ac_high$chromosome_name), ]
id_ac_low<-id_ac_low[order(id_ac_low$chromosome_name), ]
id_ac_zero<-id_ac_zero[order(id_ac_zero$chromosome_name), ]
id_fgd5_high<-id_fgd5_high[order(id_fgd5_high$chromosome_name), ]
id_fgd5_low<-id_fgd5_low[order(id_fgd5_low$chromosome_name), ]
id_fgd5_zero<-id_fgd5_zero[order(id_fgd5_zero$chromosome_name), ]
id_rp11_high<-id_rp11_high[order(id_rp11_high$chromosome_name), ]
id_rp11_low<-id_rp11_low[order(id_rp11_low$chromosome_name), ]
id_rp11_zero<-id_rp11_zero[order(id_rp11_zero$chromosome_name), ]
id_mapkapk_high<-id_mapkapk_high[order(id_mapkapk_high$chromosome_name), ]
id_mapkapk_low<-id_mapkapk_low[order(id_mapkapk_low$chromosome_name), ]
id_mapkapk_zero<-id_mapkapk_zero[order(id_mapkapk_zero$chromosome_name), ]
id_ctd_high<-id_ctd_high[order(id_ctd_high$chromosome_name), ]
id_ctd_low<-id_ctd_low[order(id_ctd_low$chromosome_name), ]
id_ctd_zero<-id_ctd_zero[order(id_ctd_zero$chromosome_name), ]
id_jpx_high<-id_jpx_high[order(id_jpx_high$chromosome_name), ]
id_jpx_low<-id_jpx_low[order(id_jpx_low$chromosome_name), ]
id_jpx_zero<-id_jpx_zero[order(id_jpx_zero$chromosome_name), ]
id_emx2os_high<-id_emx2os_high[order(id_emx2os_high$chromosome_name), ]
id_emx2os_low<-id_emx2os_low[order(id_emx2os_low$chromosome_name), ]
id_emx2os_zero<-id_emx2os_zero[order(id_emx2os_zero$chromosome_name), ]
id_nc_high<-id_nc_high[order(id_nc_high$chromosome_name), ]
id_nc_low<-id_nc_low[order(id_nc_low$chromosome_name), ]
id_nc_zero<-id_nc_zero[order(id_nc_zero$chromosome_name), ]

id_ac_high$strand <-replace(id_ac_high$strand, id_ac_high$strand==-1,"-")
id_ac_low$strand <-replace(id_ac_low$strand, id_ac_low$strand==-1,"-")
id_ac_zero$strand <-replace(id_ac_zero$strand, id_ac_zero$strand==-1,"-")
id_fgd5_high$strand <-replace(id_fgd5_high$strand, id_fgd5_high$strand==-1,"-")
id_fgd5_low$strand <-replace(id_fgd5_low$strand, id_fgd5_low$strand==-1,"-")
id_fgd5_zero$strand <-replace(id_fgd5_zero$strand, id_fgd5_zero$strand==-1,"-")
id_rp11_high$strand <-replace(id_rp11_high$strand, id_rp11_high$strand==-1,"-")
id_rp11_low$strand <-replace(id_rp11_low$strand, id_rp11_low$strand==-1,"-")
id_rp11_zero$strand <-replace(id_rp11_zero$strand, id_rp11_zero$strand==-1,"-")
id_mapkapk_high$strand <-replace(id_mapkapk_high$strand, id_mapkapk_high$strand==-1,"-")
id_mapkapk_low$strand <-replace(id_mapkapk_low$strand, id_mapkapk_low$strand==-1,"-")
id_mapkapk_zero$strand <-replace(id_mapkapk_zero$strand, id_mapkapk_zero$strand==-1,"-")
id_ctd_high$strand <-replace(id_ctd_high$strand, id_ctd_high$strand==-1,"-")
id_ctd_low$strand <-replace(id_ctd_low$strand, id_ctd_low$strand==-1,"-")
id_ctd_zero$strand <-replace(id_ctd_zero$strand, id_ctd_zero$strand==-1,"-")
id_jpx_high$strand <-replace(id_jpx_high$strand, id_jpx_high$strand==-1,"-")
id_jpx_low$strand <-replace(id_jpx_low$strand, id_jpx_low$strand==-1,"-")
id_jpx_zero$strand <-replace(id_jpx_zero$strand, id_jpx_zero$strand==-1,"-")
id_emx2os_high$strand <-replace(id_emx2os_high$strand, id_emx2os_high$strand==-1,"-")
id_emx2os_low$strand <-replace(id_emx2os_low$strand, id_emx2os_low$strand==-1,"-")
id_emx2os_zero$strand <-replace(id_emx2os_zero$strand, id_emx2os_zero$strand==-1,"-")
id_nc_high$strand <-replace(id_nc_high$strand, id_nc_high$strand==-1,"-")
id_nc_low$strand <-replace(id_nc_low$strand, id_nc_low$strand==-1,"-")
id_nc_zero$strand <-replace(id_nc_zero$strand, id_nc_zero$strand==-1,"-")


id_ac_high$strand <-replace(id_ac_high$strand, id_ac_high$strand==1,"+")
id_ac_low$strand <-replace(id_ac_low$strand, id_ac_low$strand==1,"+")
id_ac_zero$strand <-replace(id_ac_zero$strand, id_ac_zero$strand==1,"+")
id_fgd5_high$strand <-replace(id_fgd5_high$strand, id_fgd5_high$strand==1,"+")
id_fgd5_low$strand <-replace(id_fgd5_low$strand, id_fgd5_low$strand==1,"+")
id_fgd5_zero$strand <-replace(id_fgd5_zero$strand, id_fgd5_zero$strand==1,"+")
id_rp11_high$strand <-replace(id_rp11_high$strand, id_rp11_high$strand==1,"+")
id_rp11_low$strand <-replace(id_rp11_low$strand, id_rp11_low$strand==1,"+")
id_rp11_zero$strand <-replace(id_rp11_zero$strand, id_rp11_zero$strand==1,"+")
id_mapkapk_high$strand <-replace(id_mapkapk_high$strand, id_mapkapk_high$strand==1,"+")
id_mapkapk_low$strand <-replace(id_mapkapk_low$strand, id_mapkapk_low$strand==1,"+")
id_mapkapk_zero$strand <-replace(id_mapkapk_zero$strand, id_mapkapk_zero$strand==1,"+")
id_ctd_high$strand <-replace(id_ctd_high$strand, id_ctd_high$strand==1,"+")
id_ctd_low$strand <-replace(id_ctd_low$strand, id_ctd_low$strand==1,"+")
id_ctd_zero$strand <-replace(id_ctd_zero$strand, id_ctd_zero$strand==1,"+")
id_jpx_high$strand <-replace(id_jpx_high$strand, id_jpx_high$strand==1,"+")
id_jpx_low$strand <-replace(id_jpx_low$strand, id_jpx_low$strand==1,"+")
id_jpx_zero$strand <-replace(id_jpx_zero$strand, id_jpx_zero$strand==1,"+")
id_emx2os_high$strand <-replace(id_emx2os_high$strand, id_emx2os_high$strand==1,"+")
id_emx2os_low$strand <-replace(id_emx2os_low$strand, id_emx2os_low$strand==1,"+")
id_emx2os_zero$strand <-replace(id_emx2os_zero$strand, id_emx2os_zero$strand==1,"+")
id_nc_high$strand <-replace(id_nc_high$strand, id_nc_high$strand==1,"+")
id_nc_low$strand <-replace(id_nc_low$strand, id_nc_low$strand==1,"+")
id_nc_zero$strand <-replace(id_nc_zero$strand, id_nc_zero$strand==1,"+")

id_ac_high$name<-(x="0")
id_ac_low$name<-(x="0")
id_ac_zero$name<-(x="0")
id_fgd5_high$name<-(x="0")
id_fgd5_low$name<-(x="0")
id_fgd5_zero$name<-(x="0")
id_rp11_high$name<-(x="0")
id_rp11_low$name<-(x="0")
id_rp11_zero$name<-(x="0")
id_mapkapk_high$name<-(x="0")
id_mapkapk_low$name<-(x="0")
id_mapkapk_zero$name<-(x="0")
id_ctd_high$name<-(x="0")
id_ctd_low$name<-(x="0")
id_ctd_zero$name<-(x="0")
id_jpx_high$name<-(x="0")
id_jpx_low$name<-(x="0")
id_jpx_zero$name<-(x="0")
id_emx2os_high$name<-(x="0")
id_emx2os_low$name<-(x="0")
id_emx2os_zero$name<-(x="0")
id_nc_high$name<-(x="0")
id_nc_low$name<-(x="0")
id_nc_zero$name<-(x="0")

id_ac_high$phase<-(x="0")
id_ac_low$phase<-(x="0")
id_ac_zero$phase<-(x="0")
id_fgd5_high$phase<-(x="0")
id_fgd5_low$phase<-(x="0")
id_fgd5_zero$phase<-(x="0")
id_rp11_high$phase<-(x="0")
id_rp11_low$phase<-(x="0")
id_rp11_zero$phase<-(x="0")
id_mapkapk_high$phase<-(x="0")
id_mapkapk_low$phase<-(x="0")
id_mapkapk_zero$phase<-(x="0")
id_ctd_high$phase<-(x="0")
id_ctd_low$phase<-(x="0")
id_ctd_zero$phase<-(x="0")
id_jpx_high$phase<-(x="0")
id_jpx_low$phase<-(x="0")
id_jpx_zero$phase<-(x="0")
id_emx2os_high$phase<-(x="0")
id_emx2os_low$phase<-(x="0")
id_emx2os_zero$phase<-(x="0")
id_nc_high$phase<-(x="0")
id_nc_low$phase<-(x="0")
id_nc_zero$phase<-(x="0")

id_ac_high<-id_ac_high[,c(1,2,3,5,6,4)]
id_ac_low<-id_ac_low[,c(1,2,3,5,6,4)]
id_ac_zero<-id_ac_zero[,c(1,2,3,5,6,4)]
id_fgd5_high<-id_fgd5_high[,c(1,2,3,5,6,4)]
id_fgd5_low<-id_fgd5_low[,c(1,2,3,5,6,4)]
id_fgd5_zero<-id_fgd5_zero[,c(1,2,3,5,6,4)]
id_rp11_high<-id_rp11_high[,c(1,2,3,5,6,4)]
id_rp11_low<-id_rp11_low[,c(1,2,3,5,6,4)]
id_rp11_zero<-id_rp11_zero[,c(1,2,3,5,6,4)]
id_mapkapk_high<-id_mapkapk_high[,c(1,2,3,5,6,4)]
id_mapkapk_low<-id_mapkapk_low[,c(1,2,3,5,6,4)]
id_mapkapk_zero<-id_mapkapk_zero[,c(1,2,3,5,6,4)]
id_ctd_high<-id_ctd_high[,c(1,2,3,5,6,4)]
id_ctd_low<-id_ctd_low[,c(1,2,3,5,6,4)]
id_ctd_zero<-id_ctd_zero[,c(1,2,3,5,6,4)]
id_jpx_high<-id_jpx_high[,c(1,2,3,5,6,4)]
id_jpx_low<-id_jpx_low[,c(1,2,3,5,6,4)]
id_jpx_zero<-id_jpx_zero[,c(1,2,3,5,6,4)]
id_emx2os_high<-id_emx2os_high[,c(1,2,3,5,6,4)]
id_emx2os_low<-id_emx2os_low[,c(1,2,3,5,6,4)]
id_emx2os_zero<-id_emx2os_zero[,c(1,2,3,5,6,4)]
id_nc_high<-id_nc_high[,c(1,2,3,5,6,4)]
id_nc_low<-id_nc_low[,c(1,2,3,5,6,4)]
id_nc_zero<-id_nc_zero[,c(1,2,3,5,6,4)]

write.table(data.frame(id_ac_high), "~/Downloads/bed/ac_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ac_low), "~/Downloads/bed/ac_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ac_zero), "~/Downloads/bed/ac_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_fgd5_high), "~/Downloads/bed/fgd5_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_fgd5_low), "~/Downloads/bed/fgd5_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_fgd5_zero), "~/Downloads/bed/fgd5_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_rp11_high), "~/Downloads/bed/rp11_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_rp11_low), "~/Downloads/bed/rp11_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_rp11_zero), "~/Downloads/bed/rp11_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_mapkapk_high), "~/Downloads/bed/mapkapk_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_mapkapk_low), "~/Downloads/bed/mapkapk_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_mapkapk_zero), "~/Downloads/bed/mapkapk_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ctd_high), "~/Downloads/bed/ctd_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ctd_low), "~/Downloads/bed/ctd_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ctd_zero), "~/Downloads/bed/ctd_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_jpx_high), "~/Downloads/bed/jpx_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_jpx_low), "~/Downloads/bed/jpx_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_jpx_zero), "~/Downloads/bed/jpx_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_emx2os_high), "~/Downloads/bed/emx2os_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_emx2os_low), "~/Downloads/bed/emx2os_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_emx2os_zero), "~/Downloads/bed/emx2os_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_nc_high), "~/Downloads/bed/nc_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_nc_low), "~/Downloads/bed/nc_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_nc_zero), "~/Downloads/bed/nc_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)




ac_02_high <- read.table("~/Desktop/commands/py/genes/ac_02_high.txt", header = FALSE, sep = '\t')
ac_02_low <- read.table("~/Desktop/commands/py/genes/ac_02_low.txt", header = FALSE, sep = '\t')
ac_02_zero <- read.table("~/Desktop/commands/py/genes/ac_02_zero.txt", header = FALSE, sep = '\t')
ac_03_high <- read.table("~/Desktop/commands/py/genes/ac_03_high.txt", header = FALSE, sep = '\t')
ac_03_low <- read.table("~/Desktop/commands/py/genes/ac_03_low.txt", header = FALSE, sep = '\t')
ac_03_zero <- read.table("~/Desktop/commands/py/genes/ac_03_zero.txt", header = FALSE, sep = '\t')

fgd5_04_high <- read.table("~/Desktop/commands/py/genes/fgd5_04_high.txt", header = FALSE, sep = '\t')
fgd5_04_low <- read.table("~/Desktop/commands/py/genes/fgd5_04_low.txt", header = FALSE, sep = '\t')
fgd5_04_zero <- read.table("~/Desktop/commands/py/genes/fgd5_04_zero.txt", header = FALSE, sep = '\t')
fgd5_ad03_high <- read.table("~/Desktop/commands/py/genes/fgd5_ad03_high.txt", header = FALSE, sep = '\t')
fgd5_ad03_low <- read.table("~/Desktop/commands/py/genes/fgd5_ad03_low.txt", header = FALSE, sep = '\t')
fgd5_ad03_zero <- read.table("~/Desktop/commands/py/genes/fgd5_ad03_zero.txt", header = FALSE, sep = '\t')

rp11_03_high <- read.table("~/Desktop/commands/py/genes/rp11_03_high.txt", header = FALSE, sep = '\t')
rp11_03_low <- read.table("~/Desktop/commands/py/genes/rp11_03_low.txt", header = FALSE, sep = '\t')
rp11_03_zero <- read.table("~/Desktop/commands/py/genes/rp11_03_zero.txt", header = FALSE, sep = '\t')
rp11_05_high <- read.table("~/Desktop/commands/py/genes/rp11_05_high.txt", header = FALSE, sep = '\t')
rp11_05_low <- read.table("~/Desktop/commands/py/genes/rp11_05_low.txt", header = FALSE, sep = '\t')
rp11_05_zero <- read.table("~/Desktop/commands/py/genes/rp11_05_zero.txt", header = FALSE, sep = '\t')

mapkapk_06_high <- read.table("~/Desktop/commands/py/genes/mapkapk_06_high.txt", header = FALSE, sep = '\t')
mapkapk_06_low <- read.table("~/Desktop/commands/py/genes/mapkapk_06_low.txt", header = FALSE, sep = '\t')
mapkapk_06_zero <- read.table("~/Desktop/commands/py/genes/mapkapk_06_zero.txt", header = FALSE, sep = '\t')
mapkapk_ad04_high <- read.table("~/Desktop/commands/py/genes/mapkapk_ad04_high.txt", header = FALSE, sep = '\t')
mapkapk_ad04_low <- read.table("~/Desktop/commands/py/genes/mapkapk_ad04_low.txt", header = FALSE, sep = '\t')
mapkapk_ad04_zero <- read.table("~/Desktop/commands/py/genes/mapkapk_ad04_zero.txt", header = FALSE, sep = '\t')

ctd_01_high <- read.table("~/Desktop/commands/py/genes/ctd_01_high.txt", header = FALSE, sep = '\t')
ctd_01_low <- read.table("~/Desktop/commands/py/genes/ctd_01_low.txt", header = FALSE, sep = '\t')
ctd_01_zero <- read.table("~/Desktop/commands/py/genes/ctd_01_zero.txt", header = FALSE, sep = '\t')
ctd_03_high <- read.table("~/Desktop/commands/py/genes/ctd_03_high.txt", header = FALSE, sep = '\t')
ctd_03_low <- read.table("~/Desktop/commands/py/genes/ctd_03_low.txt", header = FALSE, sep = '\t')
ctd_03_zero <- read.table("~/Desktop/commands/py/genes/ctd_03_zero.txt", header = FALSE, sep = '\t')

jpx_05_high <- read.table("~/Desktop/commands/py/genes/jpx_05_high.txt", header = FALSE, sep = '\t')
jpx_05_low <- read.table("~/Desktop/commands/py/genes/jpx_05_low.txt", header = FALSE, sep = '\t')
jpx_05_zero <- read.table("~/Desktop/commands/py/genes/jpx_05_zero.txt", header = FALSE, sep = '\t')
jpx_ad04_high <- read.table("~/Desktop/commands/py/genes/jpx_ad04_high.txt", header = FALSE, sep = '\t')
jpx_ad04_low <- read.table("~/Desktop/commands/py/genes/jpx_ad04_low.txt", header = FALSE, sep = '\t')
jpx_ad04_zero <- read.table("~/Desktop/commands/py/genes/jpx_ad04_zero.txt", header = FALSE, sep = '\t')

emx2os_ad02_high <- read.table("~/Desktop/commands/py/genes/emx2os_ad02_high.txt", header = FALSE, sep = '\t')
emx2os_ad02_low <- read.table("~/Desktop/commands/py/genes/emx2os_ad02_low.txt", header = FALSE, sep = '\t')
emx2os_ad02_zero <- read.table("~/Desktop/commands/py/genes/emx2os_ad02_zero.txt", header = FALSE, sep = '\t')
emx2os_ad04_high <- read.table("~/Desktop/commands/py/genes/emx2os_ad04_high.txt", header = FALSE, sep = '\t')
emx2os_ad04_low <- read.table("~/Desktop/commands/py/genes/emx2os_ad04_low.txt", header = FALSE, sep = '\t')
emx2os_ad04_zero <- read.table("~/Desktop/commands/py/genes/emx2os_ad04_zero.txt", header = FALSE, sep = '\t')


colnames(ac_02_high) <- c("ensembl_gene_id")
colnames(ac_02_low) <- c("ensembl_gene_id")
colnames(ac_02_zero) <- c("ensembl_gene_id")
colnames(ac_03_high) <- c("ensembl_gene_id")
colnames(ac_03_low) <- c("ensembl_gene_id")
colnames(ac_03_zero) <- c("ensembl_gene_id")

colnames(fgd5_04_high) <- c("ensembl_gene_id")
colnames(fgd5_04_low) <- c("ensembl_gene_id")
colnames(fgd5_04_zero) <- c("ensembl_gene_id")
colnames(fgd5_ad03_high) <- c("ensembl_gene_id")
colnames(fgd5_ad03_low) <- c("ensembl_gene_id")
colnames(fgd5_ad03_zero) <- c("ensembl_gene_id")

colnames(rp11_03_high) <- c("ensembl_gene_id")
colnames(rp11_03_low) <- c("ensembl_gene_id")
colnames(rp11_03_zero) <- c("ensembl_gene_id")
colnames(rp11_05_high) <- c("ensembl_gene_id")
colnames(rp11_05_low) <- c("ensembl_gene_id")
colnames(rp11_05_zero) <- c("ensembl_gene_id")

colnames(mapkapk_06_high) <- c("ensembl_gene_id")
colnames(mapkapk_06_low) <- c("ensembl_gene_id")
colnames(mapkapk_06_zero) <- c("ensembl_gene_id")
colnames(mapkapk_ad04_high) <- c("ensembl_gene_id")
colnames(mapkapk_ad04_low) <- c("ensembl_gene_id")
colnames(mapkapk_ad04_zero) <- c("ensembl_gene_id")

colnames(ctd_01_high) <- c("ensembl_gene_id")
colnames(ctd_01_low) <- c("ensembl_gene_id")
colnames(ctd_01_zero) <- c("ensembl_gene_id")
colnames(ctd_03_high) <- c("ensembl_gene_id")
colnames(ctd_03_low) <- c("ensembl_gene_id")
colnames(ctd_03_zero) <- c("ensembl_gene_id")

colnames(jpx_05_high) <- c("ensembl_gene_id")
colnames(jpx_05_low) <- c("ensembl_gene_id")
colnames(jpx_05_zero) <- c("ensembl_gene_id")
colnames(jpx_ad04_high) <- c("ensembl_gene_id")
colnames(jpx_ad04_low) <- c("ensembl_gene_id")
colnames(jpx_ad04_zero) <- c("ensembl_gene_id")

colnames(emx2os_ad02_high) <- c("ensembl_gene_id")
colnames(emx2os_ad02_low) <- c("ensembl_gene_id")
colnames(emx2os_ad02_zero) <- c("ensembl_gene_id")
colnames(emx2os_ad04_high) <- c("ensembl_gene_id")
colnames(emx2os_ad04_low) <- c("ensembl_gene_id")
colnames(emx2os_ad04_zero) <- c("ensembl_gene_id")


id_ac_02_high <- merge(x = ac_02_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ac_02_low <- merge(x = ac_02_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ac_02_zero <- merge(x = ac_02_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ac_03_high <- merge(x = ac_03_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ac_03_low <- merge(x = ac_03_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ac_03_zero <- merge(x = ac_03_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)

id_fgd5_04_high <- merge(x = fgd5_04_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_fgd5_04_low <- merge(x = fgd5_04_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_fgd5_04_zero <- merge(x = fgd5_04_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_fgd5_ad03_high <- merge(x = fgd5_ad03_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_fgd5_ad03_low <- merge(x = fgd5_ad03_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_fgd5_ad03_zero <- merge(x = fgd5_ad03_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)

id_rp11_03_high <- merge(x = rp11_03_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_rp11_03_low <- merge(x = rp11_03_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_rp11_03_zero <- merge(x = rp11_03_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_rp11_05_high <- merge(x = rp11_05_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_rp11_05_low <- merge(x = rp11_05_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_rp11_05_zero <- merge(x = rp11_05_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)

id_mapkapk_06_high <- merge(x = mapkapk_06_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_mapkapk_06_low <- merge(x = mapkapk_06_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_mapkapk_06_zero <- merge(x = mapkapk_06_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_mapkapk_ad04_high <- merge(x = mapkapk_ad04_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_mapkapk_ad04_low <- merge(x = mapkapk_ad04_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_mapkapk_ad04_zero <- merge(x = mapkapk_ad04_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)

id_ctd_01_high <- merge(x = ctd_01_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ctd_01_low <- merge(x = ctd_01_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ctd_01_zero <- merge(x = ctd_01_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ctd_03_high <- merge(x = ctd_03_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ctd_03_low <- merge(x = ctd_03_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_ctd_03_zero <- merge(x = ctd_03_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)

id_jpx_05_high <- merge(x = jpx_05_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_jpx_05_low <- merge(x = jpx_05_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_jpx_05_zero <- merge(x = jpx_05_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_jpx_ad04_high <- merge(x = jpx_ad04_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_jpx_ad04_low <- merge(x = jpx_ad04_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_jpx_ad04_zero <- merge(x = jpx_ad04_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)

id_emx2os_ad02_high <- merge(x = emx2os_ad02_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_emx2os_ad02_low <- merge(x = emx2os_ad02_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_emx2os_ad02_zero <- merge(x = emx2os_ad02_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_emx2os_ad04_high <- merge(x = emx2os_ad04_high, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_emx2os_ad04_low <- merge(x = emx2os_ad04_low, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)
id_emx2os_ad04_zero <- merge(x = emx2os_ad04_zero, y = t2g, by = 'ensembl_gene_id', all.x = TRUE)


id_ac_02_high <- na.omit(id_ac_02_high)
id_ac_02_low <- na.omit(id_ac_02_low)
id_ac_02_zero <- na.omit(id_ac_02_zero)
id_ac_03_high <- na.omit(id_ac_03_high)
id_ac_03_low <- na.omit(id_ac_03_low)
id_ac_03_zero <- na.omit(id_ac_03_zero)

id_fgd5_04_high <- na.omit(id_fgd5_04_high)
id_fgd5_04_low <- na.omit(id_fgd5_04_low)
id_fgd5_04_zero <- na.omit(id_fgd5_04_zero)
id_fgd5_ad03_high <- na.omit(id_fgd5_ad03_high)
id_fgd5_ad03_low <- na.omit(id_fgd5_ad03_low)
id_fgd5_ad03_zero <- na.omit(id_fgd5_ad03_zero)

id_rp11_03_high <- na.omit(id_rp11_03_high)
id_rp11_03_low <- na.omit(id_rp11_03_low)
id_rp11_03_zero <- na.omit(id_rp11_03_zero)
id_rp11_05_high <- na.omit(id_rp11_05_high)
id_rp11_05_low <- na.omit(id_rp11_05_low)
id_rp11_05_zero <- na.omit(id_rp11_05_zero)

id_mapkapk_06_high <- na.omit(id_mapkapk_06_high)
id_mapkapk_06_low <- na.omit(id_mapkapk_06_low)
id_mapkapk_06_zero <- na.omit(id_mapkapk_06_zero)
id_mapkapk_ad04_high <- na.omit(id_mapkapk_ad04_high)
id_mapkapk_ad04_low <- na.omit(id_mapkapk_ad04_low)
id_mapkapk_ad04_zero <- na.omit(id_mapkapk_ad04_zero)

id_ctd_01_high <- na.omit(id_ctd_01_high)
id_ctd_01_low <- na.omit(id_ctd_01_low)
id_ctd_01_zero <- na.omit(id_ctd_01_zero)
id_ctd_03_high <- na.omit(id_ctd_03_high)
id_ctd_03_low <- na.omit(id_ctd_03_low)
id_ctd_03_zero <- na.omit(id_ctd_03_zero)

id_jpx_05_high <- na.omit(id_jpx_05_high)
id_jpx_05_low <- na.omit(id_jpx_05_low)
id_jpx_05_zero <- na.omit(id_jpx_05_zero)
id_jpx_ad04_high <- na.omit(id_jpx_ad04_high)
id_jpx_ad04_low <- na.omit(id_jpx_ad04_low)
id_jpx_ad04_zero <- na.omit(id_jpx_ad04_zero)

id_emx2os_ad02_high <- na.omit(id_emx2os_ad02_high)
id_emx2os_ad02_low <- na.omit(id_emx2os_ad02_low)
id_emx2os_ad02_zero <- na.omit(id_emx2os_ad02_zero)
id_emx2os_ad04_high <- na.omit(id_emx2os_ad04_high)
id_emx2os_ad04_low <- na.omit(id_emx2os_ad04_low)
id_emx2os_ad04_zero <- na.omit(id_emx2os_ad04_zero)


id_ac_02_high$chromosome_name <- paste0(x = "chr", id_ac_02_high$chromosome_name)
id_ac_02_low$chromosome_name <- paste0(x = "chr", id_ac_02_low$chromosome_name)
id_ac_02_zero$chromosome_name <- paste0(x = "chr", id_ac_02_zero$chromosome_name)
id_ac_03_high$chromosome_name <- paste0(x = "chr", id_ac_03_high$chromosome_name)
id_ac_03_low$chromosome_name <- paste0(x = "chr", id_ac_03_low$chromosome_name)
id_ac_03_zero$chromosome_name <- paste0(x = "chr", id_ac_03_zero$chromosome_name)

id_fgd5_04_high$chromosome_name <- paste0(x = "chr", id_fgd5_04_high$chromosome_name)
id_fgd5_04_low$chromosome_name <- paste0(x = "chr", id_fgd5_04_low$chromosome_name)
id_fgd5_04_zero$chromosome_name <- paste0(x = "chr", id_fgd5_04_zero$chromosome_name)
id_fgd5_ad03_high$chromosome_name <- paste0(x = "chr", id_fgd5_ad03_high$chromosome_name)
id_fgd5_ad03_low$chromosome_name <- paste0(x = "chr", id_fgd5_ad03_low$chromosome_name)
id_fgd5_ad03_zero$chromosome_name <- paste0(x = "chr", id_fgd5_ad03_zero$chromosome_name)

id_rp11_03_high$chromosome_name <- paste0(x = "chr", id_rp11_03_high$chromosome_name)
id_rp11_03_low$chromosome_name <- paste0(x = "chr", id_rp11_03_low$chromosome_name)
id_rp11_03_zero$chromosome_name <- paste0(x = "chr", id_rp11_03_zero$chromosome_name)
id_rp11_05_high$chromosome_name <- paste0(x = "chr", id_rp11_05_high$chromosome_name)
id_rp11_05_low$chromosome_name <- paste0(x = "chr", id_rp11_05_low$chromosome_name)
id_rp11_05_zero$chromosome_name <- paste0(x = "chr", id_rp11_05_zero$chromosome_name)

id_mapkapk_06_high$chromosome_name <- paste0(x = "chr", id_mapkapk_06_high$chromosome_name)
id_mapkapk_06_low$chromosome_name <- paste0(x = "chr", id_mapkapk_06_low$chromosome_name)
id_mapkapk_06_zero$chromosome_name <- paste0(x = "chr", id_mapkapk_06_zero$chromosome_name)
id_mapkapk_ad04_high$chromosome_name <- paste0(x = "chr", id_mapkapk_ad04_high$chromosome_name)
id_mapkapk_ad04_low$chromosome_name <- paste0(x = "chr", id_mapkapk_ad04_low$chromosome_name)
id_mapkapk_ad04_zero$chromosome_name <- paste0(x = "chr", id_mapkapk_ad04_zero$chromosome_name)

id_ctd_01_high$chromosome_name <- paste0(x = "chr", id_ctd_01_high$chromosome_name)
id_ctd_01_low$chromosome_name <- paste0(x = "chr", id_ctd_01_low$chromosome_name)
id_ctd_01_zero$chromosome_name <- paste0(x = "chr", id_ctd_01_zero$chromosome_name)
id_ctd_03_high$chromosome_name <- paste0(x = "chr", id_ctd_03_high$chromosome_name)
id_ctd_03_low$chromosome_name <- paste0(x = "chr", id_ctd_03_low$chromosome_name)
id_ctd_03_zero$chromosome_name <- paste0(x = "chr", id_ctd_03_zero$chromosome_name)

id_jpx_05_high$chromosome_name <- paste0(x = "chr", id_jpx_05_high$chromosome_name)
id_jpx_05_low$chromosome_name <- paste0(x = "chr", id_jpx_05_low$chromosome_name)
id_jpx_05_zero$chromosome_name <- paste0(x = "chr", id_jpx_05_zero$chromosome_name)
id_jpx_ad04_high$chromosome_name <- paste0(x = "chr", id_jpx_ad04_high$chromosome_name)
id_jpx_ad04_low$chromosome_name <- paste0(x = "chr", id_jpx_ad04_low$chromosome_name)
id_jpx_ad04_zero$chromosome_name <- paste0(x = "chr", id_jpx_ad04_zero$chromosome_name)

id_emx2os_ad02_high$chromosome_name <- paste0(x = "chr", id_emx2os_ad02_high$chromosome_name)
id_emx2os_ad02_low$chromosome_name <- paste0(x = "chr", id_emx2os_ad02_low$chromosome_name)
id_emx2os_ad02_zero$chromosome_name <- paste0(x = "chr", id_emx2os_ad02_zero$chromosome_name)
id_emx2os_ad04_high$chromosome_name <- paste0(x = "chr", id_emx2os_ad04_high$chromosome_name)
id_emx2os_ad04_low$chromosome_name <- paste0(x = "chr", id_emx2os_ad04_low$chromosome_name)
id_emx2os_ad04_zero$chromosome_name <- paste0(x = "chr", id_emx2os_ad04_zero$chromosome_name)


id_ac_02_high <- id_ac_02_high[,c(2,3,4,5)]
id_ac_02_low <- id_ac_02_low[,c(2,3,4,5)]
id_ac_02_zero <- id_ac_02_zero[,c(2,3,4,5)]
id_ac_03_high <- id_ac_03_high[,c(2,3,4,5)]
id_ac_03_low <- id_ac_03_low[,c(2,3,4,5)]
id_ac_03_zero <- id_ac_03_zero[,c(2,3,4,5)]

id_fgd5_04_high <- id_fgd5_04_high[,c(2,3,4,5)]
id_fgd5_04_low <- id_fgd5_04_low[,c(2,3,4,5)]
id_fgd5_04_zero <- id_fgd5_04_zero[,c(2,3,4,5)]
id_fgd5_ad03_high <- id_fgd5_ad03_high[,c(2,3,4,5)]
id_fgd5_ad03_low <- id_fgd5_ad03_low[,c(2,3,4,5)]
id_fgd5_ad03_zero <- id_fgd5_ad03_zero[,c(2,3,4,5)]

id_rp11_03_high <- id_rp11_03_high[,c(2,3,4,5)]
id_rp11_03_low <- id_rp11_03_low[,c(2,3,4,5)]
id_rp11_03_zero <- id_rp11_03_zero[,c(2,3,4,5)]
id_rp11_05_high <- id_rp11_05_high[,c(2,3,4,5)]
id_rp11_05_low <- id_rp11_05_low[,c(2,3,4,5)]
id_rp11_05_zero <- id_rp11_05_zero[,c(2,3,4,5)]

id_mapkapk_06_high <- id_mapkapk_06_high[,c(2,3,4,5)]
id_mapkapk_06_low <- id_mapkapk_06_low[,c(2,3,4,5)]
id_mapkapk_06_zero <- id_mapkapk_06_zero[,c(2,3,4,5)]
id_mapkapk_ad04_high <- id_mapkapk_ad04_high[,c(2,3,4,5)]
id_mapkapk_ad04_low <- id_mapkapk_ad04_low[,c(2,3,4,5)]
id_mapkapk_ad04_zero <- id_mapkapk_ad04_zero[,c(2,3,4,5)]

id_ctd_01_high <- id_ctd_01_high[,c(2,3,4,5)]
id_ctd_01_low <- id_ctd_01_low[,c(2,3,4,5)]
id_ctd_01_zero <- id_ctd_01_zero[,c(2,3,4,5)]
id_ctd_03_high <- id_ctd_03_high[,c(2,3,4,5)]
id_ctd_03_low <- id_ctd_03_low[,c(2,3,4,5)]
id_ctd_03_zero <- id_ctd_03_zero[,c(2,3,4,5)]

id_jpx_05_high <- id_jpx_05_high[,c(2,3,4,5)]
id_jpx_05_low <- id_jpx_05_low[,c(2,3,4,5)]
id_jpx_05_zero <- id_jpx_05_zero[,c(2,3,4,5)]
id_jpx_ad04_high <- id_jpx_ad04_high[,c(2,3,4,5)]
id_jpx_ad04_low <- id_jpx_ad04_low[,c(2,3,4,5)]
id_jpx_ad04_zero <- id_jpx_ad04_zero[,c(2,3,4,5)]

id_emx2os_ad02_high <- id_emx2os_ad02_high[,c(2,3,4,5)]
id_emx2os_ad02_low <- id_emx2os_ad02_low[,c(2,3,4,5)]
id_emx2os_ad02_zero <- id_emx2os_ad02_zero[,c(2,3,4,5)]
id_emx2os_ad04_high <- id_emx2os_ad04_high[,c(2,3,4,5)]
id_emx2os_ad04_low <- id_emx2os_ad04_low[,c(2,3,4,5)]
id_emx2os_ad04_zero <- id_emx2os_ad04_zero[,c(2,3,4,5)]


id_ac_02_high <- id_ac_02_high[order(id_ac_02_high$chromosome_name), ]
id_ac_02_low <- id_ac_02_low[order(id_ac_02_low$chromosome_name), ]
id_ac_02_zero <- id_ac_02_zero[order(id_ac_02_zero$chromosome_name), ]
id_ac_03_high <- id_ac_03_high[order(id_ac_03_high$chromosome_name), ]
id_ac_03_low <- id_ac_03_low[order(id_ac_03_low$chromosome_name), ]
id_ac_03_zero <- id_ac_03_zero[order(id_ac_03_zero$chromosome_name), ]

id_fgd5_04_high <- id_fgd5_04_high[order(id_fgd5_04_high$chromosome_name), ]
id_fgd5_04_low <- id_fgd5_04_low[order(id_fgd5_04_low$chromosome_name), ]
id_fgd5_04_zero <- id_fgd5_04_zero[order(id_fgd5_04_zero$chromosome_name), ]
id_fgd5_ad03_high <- id_fgd5_ad03_high[order(id_fgd5_ad03_high$chromosome_name), ]
id_fgd5_ad03_low <- id_fgd5_ad03_low[order(id_fgd5_ad03_low$chromosome_name), ]
id_fgd5_ad03_zero <- id_fgd5_ad03_zero[order(id_fgd5_ad03_zero$chromosome_name), ]

id_rp11_03_high <- id_rp11_03_high[order(id_rp11_03_high$chromosome_name), ]
id_rp11_03_low <- id_rp11_03_low[order(id_rp11_03_low$chromosome_name), ]
id_rp11_03_zero <- id_rp11_03_zero[order(id_rp11_03_zero$chromosome_name), ]
id_rp11_05_high <- id_rp11_05_high[order(id_rp11_05_high$chromosome_name), ]
id_rp11_05_low <- id_rp11_05_low[order(id_rp11_05_low$chromosome_name), ]
id_rp11_05_zero <- id_rp11_05_zero[order(id_rp11_05_zero$chromosome_name), ]

id_mapkapk_06_high <- id_mapkapk_06_high[order(id_mapkapk_06_high$chromosome_name), ]
id_mapkapk_06_low <- id_mapkapk_06_low[order(id_mapkapk_06_low$chromosome_name), ]
id_mapkapk_06_zero <- id_mapkapk_06_zero[order(id_mapkapk_06_zero$chromosome_name), ]
id_mapkapk_ad04_high <- id_mapkapk_ad04_high[order(id_mapkapk_ad04_high$chromosome_name), ]
id_mapkapk_ad04_low <- id_mapkapk_ad04_low[order(id_mapkapk_ad04_low$chromosome_name), ]
id_mapkapk_ad04_zero <- id_mapkapk_ad04_zero[order(id_mapkapk_ad04_zero$chromosome_name), ]

id_ctd_01_high <- id_ctd_01_high[order(id_ctd_01_high$chromosome_name), ]
id_ctd_01_low <- id_ctd_01_low[order(id_ctd_01_low$chromosome_name), ]
id_ctd_01_zero <- id_ctd_01_zero[order(id_ctd_01_zero$chromosome_name), ]
id_ctd_03_high <- id_ctd_03_high[order(id_ctd_03_high$chromosome_name), ]
id_ctd_03_low <- id_ctd_03_low[order(id_ctd_03_low$chromosome_name), ]
id_ctd_03_zero <- id_ctd_03_zero[order(id_ctd_03_zero$chromosome_name), ]

id_jpx_05_high <- id_jpx_05_high[order(id_jpx_05_high$chromosome_name), ]
id_jpx_05_low <- id_jpx_05_low[order(id_jpx_05_low$chromosome_name), ]
id_jpx_05_zero <- id_jpx_05_zero[order(id_jpx_05_zero$chromosome_name), ]
id_jpx_ad04_high <- id_jpx_ad04_high[order(id_jpx_ad04_high$chromosome_name), ]
id_jpx_ad04_low <- id_jpx_ad04_low[order(id_jpx_ad04_low$chromosome_name), ]
id_jpx_ad04_zero <- id_jpx_ad04_zero[order(id_jpx_ad04_zero$chromosome_name), ]

id_emx2os_ad02_high <- id_emx2os_ad02_high[order(id_emx2os_ad02_high$chromosome_name), ]
id_emx2os_ad02_low <- id_emx2os_ad02_low[order(id_emx2os_ad02_low$chromosome_name), ]
id_emx2os_ad02_zero <- id_emx2os_ad02_zero[order(id_emx2os_ad02_zero$chromosome_name), ]
id_emx2os_ad04_high <- id_emx2os_ad04_high[order(id_emx2os_ad04_high$chromosome_name), ]
id_emx2os_ad04_low <- id_emx2os_ad04_low[order(id_emx2os_ad04_low$chromosome_name), ]
id_emx2os_ad04_zero <- id_emx2os_ad04_zero[order(id_emx2os_ad04_zero$chromosome_name), ]



id_ac_02_high$strand <-replace(id_ac_02_high$strand, id_ac_02_high$strand==-1,"-")
id_ac_02_low$strand <-replace(id_ac_02_low$strand, id_ac_02_low$strand==-1,"-")
id_ac_02_zero$strand <-replace(id_ac_02_zero$strand, id_ac_02_zero$strand==-1,"-")
id_ac_03_high$strand <-replace(id_ac_03_high$strand, id_ac_03_high$strand==-1,"-")
id_ac_03_low$strand <-replace(id_ac_03_low$strand, id_ac_03_low$strand==-1,"-")
id_ac_03_zero$strand <-replace(id_ac_03_zero$strand, id_ac_03_zero$strand==-1,"-")

id_fgd5_04_high$strand <-replace(id_fgd5_04_high$strand, id_fgd5_04_high$strand==-1,"-")
id_fgd5_04_low$strand <-replace(id_fgd5_04_low$strand, id_fgd5_04_low$strand==-1,"-")
id_fgd5_04_zero$strand <-replace(id_fgd5_04_zero$strand, id_fgd5_04_zero$strand==-1,"-")
id_fgd5_ad03_high$strand <-replace(id_fgd5_ad03_high$strand, id_fgd5_ad03_high$strand==-1,"-")
id_fgd5_ad03_low$strand <-replace(id_fgd5_ad03_low$strand, id_fgd5_ad03_low$strand==-1,"-")
id_fgd5_ad03_zero$strand <-replace(id_fgd5_ad03_zero$strand, id_fgd5_ad03_zero$strand==-1,"-")

id_rp11_03_high$strand <-replace(id_rp11_03_high$strand, id_rp11_03_high$strand==-1,"-")
id_rp11_03_low$strand <-replace(id_rp11_03_low$strand, id_rp11_03_low$strand==-1,"-")
id_rp11_03_zero$strand <-replace(id_rp11_03_zero$strand, id_rp11_03_zero$strand==-1,"-")
id_rp11_05_high$strand <-replace(id_rp11_05_high$strand, id_rp11_05_high$strand==-1,"-")
id_rp11_05_low$strand <-replace(id_rp11_05_low$strand, id_rp11_05_low$strand==-1,"-")
id_rp11_05_zero$strand <-replace(id_rp11_05_zero$strand, id_rp11_05_zero$strand==-1,"-")

id_mapkapk_06_high$strand <-replace(id_mapkapk_06_high$strand, id_mapkapk_06_high$strand==-1,"-")
id_mapkapk_06_low$strand <-replace(id_mapkapk_06_low$strand, id_mapkapk_06_low$strand==-1,"-")
id_mapkapk_06_zero$strand <-replace(id_mapkapk_06_zero$strand, id_mapkapk_06_zero$strand==-1,"-")
id_mapkapk_ad04_high$strand <-replace(id_mapkapk_ad04_high$strand, id_mapkapk_ad04_high$strand==-1,"-")
id_mapkapk_ad04_low$strand <-replace(id_mapkapk_ad04_low$strand, id_mapkapk_ad04_low$strand==-1,"-")
id_mapkapk_ad04_zero$strand <-replace(id_mapkapk_ad04_zero$strand, id_mapkapk_ad04_zero$strand==-1,"-")

id_ctd_01_high$strand <-replace(id_ctd_01_high$strand, id_ctd_01_high$strand==-1,"-")
id_ctd_01_low$strand <-replace(id_ctd_01_low$strand, id_ctd_01_low$strand==-1,"-")
id_ctd_01_zero$strand <-replace(id_ctd_01_zero$strand, id_ctd_01_zero$strand==-1,"-")
id_ctd_03_high$strand <-replace(id_ctd_03_high$strand, id_ctd_03_high$strand==-1,"-")
id_ctd_03_low$strand <-replace(id_ctd_03_low$strand, id_ctd_03_low$strand==-1,"-")
id_ctd_03_zero$strand <-replace(id_ctd_03_zero$strand, id_ctd_03_zero$strand==-1,"-")

id_jpx_05_high$strand <-replace(id_jpx_05_high$strand, id_jpx_05_high$strand==-1,"-")
id_jpx_05_low$strand <-replace(id_jpx_05_low$strand, id_jpx_05_low$strand==-1,"-")
id_jpx_05_zero$strand <-replace(id_jpx_05_zero$strand, id_jpx_05_zero$strand==-1,"-")
id_jpx_ad04_high$strand <-replace(id_jpx_ad04_high$strand, id_jpx_ad04_high$strand==-1,"-")
id_jpx_ad04_low$strand <-replace(id_jpx_ad04_low$strand, id_jpx_ad04_low$strand==-1,"-")
id_jpx_ad04_zero$strand <-replace(id_jpx_ad04_zero$strand, id_jpx_ad04_zero$strand==-1,"-")

id_emx2os_ad02_high$strand <-replace(id_emx2os_ad02_high$strand, id_emx2os_ad02_high$strand==-1,"-")
id_emx2os_ad02_low$strand <-replace(id_emx2os_ad02_low$strand, id_emx2os_ad02_low$strand==-1,"-")
id_emx2os_ad02_zero$strand <-replace(id_emx2os_ad02_zero$strand, id_emx2os_ad02_zero$strand==-1,"-")
id_emx2os_ad04_high$strand <-replace(id_emx2os_ad04_high$strand, id_emx2os_ad04_high$strand==-1,"-")
id_emx2os_ad04_low$strand <-replace(id_emx2os_ad04_low$strand, id_emx2os_ad04_low$strand==-1,"-")
id_emx2os_ad04_zero$strand <-replace(id_emx2os_ad04_zero$strand, id_emx2os_ad04_zero$strand==-1,"-")


id_ac_02_high$strand <-replace(id_ac_02_high$strand, id_ac_02_high$strand==1,"+")
id_ac_02_low$strand <-replace(id_ac_02_low$strand, id_ac_02_low$strand==1,"+")
id_ac_02_zero$strand <-replace(id_ac_02_zero$strand, id_ac_02_zero$strand==1,"+")
id_ac_03_high$strand <-replace(id_ac_03_high$strand, id_ac_03_high$strand==1,"+")
id_ac_03_low$strand <-replace(id_ac_03_low$strand, id_ac_03_low$strand==1,"+")
id_ac_03_zero$strand <-replace(id_ac_03_zero$strand, id_ac_03_zero$strand==1,"+")

id_fgd5_04_high$strand <-replace(id_fgd5_04_high$strand, id_fgd5_04_high$strand==1,"+")
id_fgd5_04_low$strand <-replace(id_fgd5_04_low$strand, id_fgd5_04_low$strand==1,"+")
id_fgd5_04_zero$strand <-replace(id_fgd5_04_zero$strand, id_fgd5_04_zero$strand==1,"+")
id_fgd5_ad03_high$strand <-replace(id_fgd5_ad03_high$strand, id_fgd5_ad03_high$strand==1,"+")
id_fgd5_ad03_low$strand <-replace(id_fgd5_ad03_low$strand, id_fgd5_ad03_low$strand==1,"+")
id_fgd5_ad03_zero$strand <-replace(id_fgd5_ad03_zero$strand, id_fgd5_ad03_zero$strand==1,"+")

id_rp11_03_high$strand <-replace(id_rp11_03_high$strand, id_rp11_03_high$strand==1,"+")
id_rp11_03_low$strand <-replace(id_rp11_03_low$strand, id_rp11_03_low$strand==1,"+")
id_rp11_03_zero$strand <-replace(id_rp11_03_zero$strand, id_rp11_03_zero$strand==1,"+")
id_rp11_05_high$strand <-replace(id_rp11_05_high$strand, id_rp11_05_high$strand==1,"+")
id_rp11_05_low$strand <-replace(id_rp11_05_low$strand, id_rp11_05_low$strand==1,"+")
id_rp11_05_zero$strand <-replace(id_rp11_05_zero$strand, id_rp11_05_zero$strand==1,"+")

id_mapkapk_06_high$strand <-replace(id_mapkapk_06_high$strand, id_mapkapk_06_high$strand==1,"+")
id_mapkapk_06_low$strand <-replace(id_mapkapk_06_low$strand, id_mapkapk_06_low$strand==1,"+")
id_mapkapk_06_zero$strand <-replace(id_mapkapk_06_zero$strand, id_mapkapk_06_zero$strand==1,"+")
id_mapkapk_ad04_high$strand <-replace(id_mapkapk_ad04_high$strand, id_mapkapk_ad04_high$strand==1,"+")
id_mapkapk_ad04_low$strand <-replace(id_mapkapk_ad04_low$strand, id_mapkapk_ad04_low$strand==1,"+")
id_mapkapk_ad04_zero$strand <-replace(id_mapkapk_ad04_zero$strand, id_mapkapk_ad04_zero$strand==1,"+")

id_ctd_01_high$strand <-replace(id_ctd_01_high$strand, id_ctd_01_high$strand==1,"+")
id_ctd_01_low$strand <-replace(id_ctd_01_low$strand, id_ctd_01_low$strand==1,"+")
id_ctd_01_zero$strand <-replace(id_ctd_01_zero$strand, id_ctd_01_zero$strand==1,"+")
id_ctd_03_high$strand <-replace(id_ctd_03_high$strand, id_ctd_03_high$strand==1,"+")
id_ctd_03_low$strand <-replace(id_ctd_03_low$strand, id_ctd_03_low$strand==1,"+")
id_ctd_03_zero$strand <-replace(id_ctd_03_zero$strand, id_ctd_03_zero$strand==1,"+")

id_jpx_05_high$strand <-replace(id_jpx_05_high$strand, id_jpx_05_high$strand==1,"+")
id_jpx_05_low$strand <-replace(id_jpx_05_low$strand, id_jpx_05_low$strand==1,"+")
id_jpx_05_zero$strand <-replace(id_jpx_05_zero$strand, id_jpx_05_zero$strand==1,"+")
id_jpx_ad04_high$strand <-replace(id_jpx_ad04_high$strand, id_jpx_ad04_high$strand==1,"+")
id_jpx_ad04_low$strand <-replace(id_jpx_ad04_low$strand, id_jpx_ad04_low$strand==1,"+")
id_jpx_ad04_zero$strand <-replace(id_jpx_ad04_zero$strand, id_jpx_ad04_zero$strand==1,"+")

id_emx2os_ad02_high$strand <-replace(id_emx2os_ad02_high$strand, id_emx2os_ad02_high$strand==1,"+")
id_emx2os_ad02_low$strand <-replace(id_emx2os_ad02_low$strand, id_emx2os_ad02_low$strand==1,"+")
id_emx2os_ad02_zero$strand <-replace(id_emx2os_ad02_zero$strand, id_emx2os_ad02_zero$strand==1,"+")
id_emx2os_ad04_high$strand <-replace(id_emx2os_ad04_high$strand, id_emx2os_ad04_high$strand==1,"+")
id_emx2os_ad04_low$strand <-replace(id_emx2os_ad04_low$strand, id_emx2os_ad04_low$strand==1,"+")
id_emx2os_ad04_zero$strand <-replace(id_emx2os_ad04_zero$strand, id_emx2os_ad04_zero$strand==1,"+")



id_ac_02_high$name<-(x="0")
id_ac_02_low$name<-(x="0")
id_ac_02_zero$name<-(x="0")
id_ac_03_high$name<-(x="0")
id_ac_03_low$name<-(x="0")
id_ac_03_zero$name<-(x="0")

id_fgd5_04_high$name<-(x="0")
id_fgd5_04_low$name<-(x="0")
id_fgd5_04_zero$name<-(x="0")
id_fgd5_ad03_high$name<-(x="0")
id_fgd5_ad03_low$name<-(x="0")
id_fgd5_ad03_zero$name<-(x="0")

id_rp11_03_high$name<-(x="0")
id_rp11_03_low$name<-(x="0")
id_rp11_03_zero$name<-(x="0")
id_rp11_05_high$name<-(x="0")
id_rp11_05_low$name<-(x="0")
id_rp11_05_zero$name<-(x="0")

id_mapkapk_06_high$name<-(x="0")
id_mapkapk_06_low$name<-(x="0")
id_mapkapk_06_zero$name<-(x="0")
id_mapkapk_ad04_high$name<-(x="0")
id_mapkapk_ad04_low$name<-(x="0")
id_mapkapk_ad04_zero$name<-(x="0")

id_ctd_01_high$name<-(x="0")
id_ctd_01_low$name<-(x="0")
id_ctd_01_zero$name<-(x="0")
id_ctd_03_high$name<-(x="0")
id_ctd_03_low$name<-(x="0")
id_ctd_03_zero$name<-(x="0")

id_jpx_05_high$name<-(x="0")
id_jpx_05_low$name<-(x="0")
id_jpx_05_zero$name<-(x="0")
id_jpx_ad04_high$name<-(x="0")
id_jpx_ad04_low$name<-(x="0")
id_jpx_ad04_zero$name<-(x="0")

id_emx2os_ad02_high$name<-(x="0")
id_emx2os_ad02_low$name<-(x="0")
id_emx2os_ad02_zero$name<-(x="0")
id_emx2os_ad04_high$name<-(x="0")
id_emx2os_ad04_low$name<-(x="0")
id_emx2os_ad04_zero$name<-(x="0")


id_ac_02_high$score<-(x="0")
id_ac_02_low$score<-(x="0")
id_ac_02_zero$score<-(x="0")
id_ac_03_high$score<-(x="0")
id_ac_03_low$score<-(x="0")
id_ac_03_zero$score<-(x="0")

id_fgd5_04_high$score<-(x="0")
id_fgd5_04_low$score<-(x="0")
id_fgd5_04_zero$score<-(x="0")
id_fgd5_ad03_high$score<-(x="0")
id_fgd5_ad03_low$score<-(x="0")
id_fgd5_ad03_zero$score<-(x="0")

id_rp11_03_high$score<-(x="0")
id_rp11_03_low$score<-(x="0")
id_rp11_03_zero$score<-(x="0")
id_rp11_05_high$score<-(x="0")
id_rp11_05_low$score<-(x="0")
id_rp11_05_zero$score<-(x="0")

id_mapkapk_06_high$score<-(x="0")
id_mapkapk_06_low$score<-(x="0")
id_mapkapk_06_zero$score<-(x="0")
id_mapkapk_ad04_high$score<-(x="0")
id_mapkapk_ad04_low$score<-(x="0")
id_mapkapk_ad04_zero$score<-(x="0")

id_ctd_01_high$score<-(x="0")
id_ctd_01_low$score<-(x="0")
id_ctd_01_zero$score<-(x="0")
id_ctd_03_high$score<-(x="0")
id_ctd_03_low$score<-(x="0")
id_ctd_03_zero$score<-(x="0")

id_jpx_05_high$score<-(x="0")
id_jpx_05_low$score<-(x="0")
id_jpx_05_zero$score<-(x="0")
id_jpx_ad04_high$score<-(x="0")
id_jpx_ad04_low$score<-(x="0")
id_jpx_ad04_zero$score<-(x="0")

id_emx2os_ad02_high$score<-(x="0")
id_emx2os_ad02_low$score<-(x="0")
id_emx2os_ad02_zero$score<-(x="0")
id_emx2os_ad04_high$score<-(x="0")
id_emx2os_ad04_low$score<-(x="0")
id_emx2os_ad04_zero$score<-(x="0")



id_ac_02_high<- id_ac_02_high[,c(1,2,3,5,6,4)]
id_ac_02_low<- id_ac_02_low[,c(1,2,3,5,6,4)]
id_ac_02_zero<- id_ac_02_zero[,c(1,2,3,5,6,4)]
id_ac_03_high<- id_ac_03_high[,c(1,2,3,5,6,4)]
id_ac_03_low<- id_ac_03_low[,c(1,2,3,5,6,4)]
id_ac_03_zero<- id_ac_03_zero[,c(1,2,3,5,6,4)]

id_fgd5_04_high<- id_fgd5_04_high[,c(1,2,3,5,6,4)]
id_fgd5_04_low<- id_fgd5_04_low[,c(1,2,3,5,6,4)]
id_fgd5_04_zero<- id_fgd5_04_zero[,c(1,2,3,5,6,4)]
id_fgd5_ad03_high<- id_fgd5_ad03_high[,c(1,2,3,5,6,4)]
id_fgd5_ad03_low<- id_fgd5_ad03_low[,c(1,2,3,5,6,4)]
id_fgd5_ad03_zero<- id_fgd5_ad03_zero[,c(1,2,3,5,6,4)]

id_rp11_03_high<- id_rp11_03_high[,c(1,2,3,5,6,4)]
id_rp11_03_low<- id_rp11_03_low[,c(1,2,3,5,6,4)]
id_rp11_03_zero<- id_rp11_03_zero[,c(1,2,3,5,6,4)]
id_rp11_05_high<- id_rp11_05_high[,c(1,2,3,5,6,4)]
id_rp11_05_low<- id_rp11_05_low[,c(1,2,3,5,6,4)]
id_rp11_05_zero<- id_rp11_05_zero[,c(1,2,3,5,6,4)]

id_mapkapk_06_high<- id_mapkapk_06_high[,c(1,2,3,5,6,4)]
id_mapkapk_06_low<- id_mapkapk_06_low[,c(1,2,3,5,6,4)]
id_mapkapk_06_zero<- id_mapkapk_06_zero[,c(1,2,3,5,6,4)]
id_mapkapk_ad04_high<- id_mapkapk_ad04_high[,c(1,2,3,5,6,4)]
id_mapkapk_ad04_low<- id_mapkapk_ad04_low[,c(1,2,3,5,6,4)]
id_mapkapk_ad04_zero<- id_mapkapk_ad04_zero[,c(1,2,3,5,6,4)]

id_ctd_01_high<- id_ctd_01_high[,c(1,2,3,5,6,4)]
id_ctd_01_low<- id_ctd_01_low[,c(1,2,3,5,6,4)]
id_ctd_01_zero<- id_ctd_01_zero[,c(1,2,3,5,6,4)]
id_ctd_03_high<- id_ctd_03_high[,c(1,2,3,5,6,4)]
id_ctd_03_low<- id_ctd_03_low[,c(1,2,3,5,6,4)]
id_ctd_03_zero<- id_ctd_03_zero[,c(1,2,3,5,6,4)]

id_jpx_05_high<- id_jpx_05_high[,c(1,2,3,5,6,4)]
id_jpx_05_low<- id_jpx_05_low[,c(1,2,3,5,6,4)]
id_jpx_05_zero<- id_jpx_05_zero[,c(1,2,3,5,6,4)]
id_jpx_ad04_high<- id_jpx_ad04_high[,c(1,2,3,5,6,4)]
id_jpx_ad04_low<- id_jpx_ad04_low[,c(1,2,3,5,6,4)]
id_jpx_ad04_zero<- id_jpx_ad04_zero[,c(1,2,3,5,6,4)]

id_emx2os_ad02_high<- id_emx2os_ad02_high[,c(1,2,3,5,6,4)]
id_emx2os_ad02_low<- id_emx2os_ad02_low[,c(1,2,3,5,6,4)]
id_emx2os_ad02_zero<- id_emx2os_ad02_zero[,c(1,2,3,5,6,4)]
id_emx2os_ad04_high<- id_emx2os_ad04_high[,c(1,2,3,5,6,4)]
id_emx2os_ad04_low<- id_emx2os_ad04_low[,c(1,2,3,5,6,4)]
id_emx2os_ad04_zero<- id_emx2os_ad04_zero[,c(1,2,3,5,6,4)]


write.table(data.frame(id_ac_02_high), "~/Downloads/bed/ac_02_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ac_02_low), "~/Downloads/bed/ac_02_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ac_02_zero), "~/Downloads/bed/ac_02_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ac_03_high), "~/Downloads/bed/ac_03_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ac_03_low), "~/Downloads/bed/ac_03_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ac_03_zero), "~/Downloads/bed/ac_03_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(id_fgd5_04_high), "~/Downloads/bed/fgd5_04_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_fgd5_04_low), "~/Downloads/bed/fgd5_04_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_fgd5_04_zero), "~/Downloads/bed/fgd5_04_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_fgd5_ad03_high), "~/Downloads/bed/fgd5_ad03_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_fgd5_ad03_low), "~/Downloads/bed/fgd5_ad03_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_fgd5_ad03_zero), "~/Downloads/bed/fgd5_ad03_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(id_rp11_03_high), "~/Downloads/bed/rp11_03_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_rp11_03_low), "~/Downloads/bed/rp11_03_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_rp11_03_zero), "~/Downloads/bed/rp11_03_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_rp11_05_high), "~/Downloads/bed/rp11_05_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_rp11_05_low), "~/Downloads/bed/rp11_05_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_rp11_05_zero), "~/Downloads/bed/rp11_05_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(id_mapkapk_06_high), "~/Downloads/bed/mapkapk_06_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_mapkapk_06_low), "~/Downloads/bed/mapkapk_06_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_mapkapk_06_zero), "~/Downloads/bed/mapkapk_06_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_mapkapk_ad04_high), "~/Downloads/bed/mapkapk_ad04_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_mapkapk_ad04_low), "~/Downloads/bed/mapkapk_ad04_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_mapkapk_ad04_zero), "~/Downloads/bed/mapkapk_ad04_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(id_ctd_01_high), "~/Downloads/bed/ctd_01_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ctd_01_low), "~/Downloads/bed/ctd_01_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ctd_01_zero), "~/Downloads/bed/ctd_01_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ctd_03_high), "~/Downloads/bed/ctd_03_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ctd_03_low), "~/Downloads/bed/ctd_03_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_ctd_03_zero), "~/Downloads/bed/ctd_03_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(id_jpx_05_high), "~/Downloads/bed/jpx_05_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_jpx_05_low), "~/Downloads/bed/jpx_05_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_jpx_05_zero), "~/Downloads/bed/jpx_05_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_jpx_ad04_high), "~/Downloads/bed/jpx_ad04_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_jpx_ad04_low), "~/Downloads/bed/jpx_ad04_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_jpx_ad04_zero), "~/Downloads/bed/jpx_ad04_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(id_emx2os_ad02_high), "~/Downloads/bed/emx2os_ad02_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_emx2os_ad02_low), "~/Downloads/bed/emx2os_ad02_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_emx2os_ad02_zero), "~/Downloads/bed/emx2os_ad02_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_emx2os_ad04_high), "~/Downloads/bed/emx2os_ad04_high.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_emx2os_ad04_low), "~/Downloads/bed/emx2os_ad04_low.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(id_emx2os_ad04_zero), "~/Downloads/bed/emx2os_ad04_zero.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

