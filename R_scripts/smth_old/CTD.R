CTD_corr_peaks <- read.table("C:/Users/Daria/Documents/CTD/CTD-2587H24.5_corr_peaks.bed", header =TRUE, sep = '\t')
CTD_corr_peaks_anno <- read.table("C:/Users/Daria/Documents/CTD/CTD-2587H24.5_corr_peaks_anno.bed", header =TRUE, sep = '\t')

peaks <- read.table("C:/Users/Daria/Documents/CTD/merged_peaks_first_in_biosample.bed", header =FALSE, sep = '\t')

lnc_peak_gene_fantom <- read.table("C:/Users/Daria/Documents/CTD/lncRNA_peaks_fantom_gene_association.tsv", header =TRUE, sep = '\t')
lnc_peak_gene <- read.table("C:/Users/Daria/Documents/CTD/lncRNA_peaks_gene_association.tsv", header =TRUE, sep = '\t')


rownames(peaks) <- peaks$V4

ctd_peak_gene_fantom <- subset(lnc_peak_gene_fantom, lnc_peak_gene_fantom$lncRNA == 'ENSG00000267577')

ctd_peak_gene_fantom$chrom <- peaks[ctd_peak_gene_fantom$peak, 1]
ctd_peak_gene_fantom$schromStart <- peaks[ctd_peak_gene_fantom$peak, 2]
ctd_peak_gene_fantom$chromEnd <- peaks[ctd_peak_gene_fantom$peak, 3]



ctd_h3k36me3_fantom <- as.data.frame(ctd_peak_gene_fantom$chrom)
colnames(ctd_h3k36me3_fantom)[colnames(ctd_h3k36me3_fantom) == 'ctd_peak_gene_fantom$chrom'] <- 'chrom'
ctd_h3k36me3_fantom$chromStart <- ctd_peak_gene_fantom$schromStart
ctd_h3k36me3_fantom$chromEnd <- ctd_peak_gene_fantom$chromEnd
ctd_h3k36me3_fantom$name <- ctd_peak_gene_fantom$peak


ctd_peak_gene <- subset(lnc_peak_gene, lnc_peak_gene$lncRNA == 'ENSG00000267577')

ctd_peak_gene$chrom <- peaks[ctd_peak_gene$peak, 1]
ctd_peak_gene$chromStart <- peaks[ctd_peak_gene$peak, 2]
ctd_peak_gene$chromEnd <- peaks[ctd_peak_gene$peak, 3]

ctd_h3k36me3 <- as.data.frame(ctd_peak_gene$chrom)
colnames(ctd_h3k36me3)[colnames(ctd_h3k36me3) == 'ctd_peak_gene$chrom'] <- 'chrom'
ctd_h3k36me3$chromStart <- ctd_peak_gene$chromStart
ctd_h3k36me3$chromEnd <- ctd_peak_gene$chromEnd
ctd_h3k36me3$name <- ctd_peak_gene$peak


output1 <- "C:/Users/Daria/Documents/CTD/CTD_peaks_himorna_f.bed"
output2 <- "C:/Users/Daria/Documents/CTD/CTD_peaks_himorna.bed"


write.table(data.frame(ctd_h3k36me3_fantom), output1, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_h3k36me3), output2, sep = '\t',  quote = F,  row.names = F, col.names = F)

##############################################################################################################
#### SICER
##############################################################################################################

CTD_01_promoter<- read.table("C:/Users/Daria/Documents/CTD/sicer/CTD-2587H24.5_01.sicer.promoter.csv", header =TRUE, sep = '\t')
CTD_03_promoter<- read.table("C:/Users/Daria/Documents/CTD/sicer/CTD-2587H24.5_03.sicer.promoter.csv", header =TRUE, sep = '\t')
CTD_01_both <- read.table("C:/Users/Daria/Documents/CTD/sicer/CTD-2587H24.5_01.sicer.both.csv", header =TRUE, sep = '\t')
CTD_03_both <- read.table("C:/Users/Daria/Documents/CTD/sicer/CTD-2587H24.5_03.sicer.both.csv", header =TRUE, sep = '\t')

chrom <- paste0(x = "chr", CTD_01_promoter$seqnames)
ctd_01_pr <- as.data.frame(chrom)
ctd_01_pr$chromStart <- CTD_01_promoter$start
ctd_01_pr$chromEnd <- CTD_01_promoter$end
ctd_01_pr$name <- CTD_01_promoter$peak

chrom <- paste0(x = "chr", CTD_03_promoter$seqnames)
ctd_03_pr <- as.data.frame(chrom)
ctd_03_pr$chromStart <- CTD_03_promoter$start
ctd_03_pr$chromEnd <- CTD_03_promoter$end
ctd_03_pr$name <- CTD_03_promoter$peak

chrom <- paste0(x = "chr", CTD_01_both$seqnames)
ctd_01_b <- as.data.frame(chrom)
ctd_01_b$chromStart <- CTD_01_both$start
ctd_01_b$chromEnd <- CTD_01_both$end
ctd_01_b$name <- CTD_01_both$peak

chrom <- paste0(x = "chr", CTD_03_both$seqnames)
ctd_03_b <- as.data.frame(chrom)
ctd_03_b$chromStart <- CTD_03_both$start
ctd_03_b$chromEnd <- CTD_03_both$end
ctd_03_b$name <- CTD_03_both$peak

output1 <- "C:/Users/Daria/Documents/CTD/sicer/CTD_01_sicer_promoter.bed"
output2 <- "C:/Users/Daria/Documents/CTD/sicer/CTD_03_sicer_promoter.bed"
output3 <- "C:/Users/Daria/Documents/CTD/sicer/CTD_01_sicer_both.bed"
output4 <- "C:/Users/Daria/Documents/CTD/sicer/CTD_03_sicer_both.bed"

write.table(data.frame(ctd_01_pr), output1, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_pr), output2, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_01_b), output3, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_b), output4, sep = '\t',  quote = F,  row.names = F, col.names = F)


##############################################################################################################
#### NARROW_PEAK
##############################################################################################################


CTD_01_promoter<- read.table("C:/Users/Daria/Documents/CTD/narrow/CTD-2587H24.5_01.narrowPeak.promoter.csv", header =TRUE, sep = '\t')
CTD_03_promoter<- read.table("C:/Users/Daria/Documents/CTD/narrow/CTD-2587H24.5_03.narrowPeak.promoter.csv", header =TRUE, sep = '\t')
CTD_01_both <- read.table("C:/Users/Daria/Documents/CTD/narrow/CTD-2587H24.5_01.narrowPeak.both.csv", header =TRUE, sep = '\t')
CTD_03_both <- read.table("C:/Users/Daria/Documents/CTD/narrow/CTD-2587H24.5_03.narrowPeak.both.csv", header =TRUE, sep = '\t')

chrom <- paste0(x = "chr", CTD_01_promoter$seqnames)
ctd_01_pr <- as.data.frame(chrom)
ctd_01_pr$chromStart <- CTD_01_promoter$start
ctd_01_pr$chromEnd <- CTD_01_promoter$end
ctd_01_pr$name <- CTD_01_promoter$peak

chrom <- paste0(x = "chr", CTD_03_promoter$seqnames)
ctd_03_pr <- as.data.frame(chrom)
ctd_03_pr$chromStart <- CTD_03_promoter$start
ctd_03_pr$chromEnd <- CTD_03_promoter$end
ctd_03_pr$name <- CTD_03_promoter$peak

chrom <- paste0(x = "chr", CTD_01_both$seqnames)
ctd_01_b <- as.data.frame(chrom)
ctd_01_b$chromStart <- CTD_01_both$start
ctd_01_b$chromEnd <- CTD_01_both$end
ctd_01_b$name <- CTD_01_both$peak

chrom <- paste0(x = "chr", CTD_03_both$seqnames)
ctd_03_b <- as.data.frame(chrom)
ctd_03_b$chromStart <- CTD_03_both$start
ctd_03_b$chromEnd <- CTD_03_both$end
ctd_03_b$name <- CTD_03_both$peak

output1 <- "C:/Users/Daria/Documents/CTD/narrow/CTD_01_narrow_promoter.bed"
output2 <- "C:/Users/Daria/Documents/CTD/narrow/CTD_03_narrow_promoter.bed"
output3 <- "C:/Users/Daria/Documents/CTD/narrow/CTD_01_narrow_both.bed"
output4 <- "C:/Users/Daria/Documents/CTD/narrow/CTD_03_narrow_both.bed"

write.table(data.frame(ctd_01_pr), output1, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_pr), output2, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_01_b), output3, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_b), output4, sep = '\t',  quote = F,  row.names = F, col.names = F)


##############################################################################################################
#### BROAD_PEAK
##############################################################################################################

CTD_01_promoter<- read.table("C:/Users/Daria/Documents/CTD/broad/CTD-2587H24.5_01.broadPeak.promoter.csv", header =TRUE, sep = '\t')
CTD_03_promoter<- read.table("C:/Users/Daria/Documents/CTD/broad/CTD-2587H24.5_03.broadPeak.promoter.csv", header =TRUE, sep = '\t')
CTD_01_both <- read.table("C:/Users/Daria/Documents/CTD/broad/CTD-2587H24.5_01.broadPeak.both.csv", header =TRUE, sep = '\t')
CTD_03_both <- read.table("C:/Users/Daria/Documents/CTD/broad/CTD-2587H24.5_03.broadPeak.both.csv", header =TRUE, sep = '\t')

chrom <- paste0(x = "chr", CTD_01_promoter$seqnames)
ctd_01_pr <- as.data.frame(chrom)
ctd_01_pr$chromStart <- CTD_01_promoter$start
ctd_01_pr$chromEnd <- CTD_01_promoter$end
ctd_01_pr$name <- CTD_01_promoter$peak

chrom <- paste0(x = "chr", CTD_03_promoter$seqnames)
ctd_03_pr <- as.data.frame(chrom)
ctd_03_pr$chromStart <- CTD_03_promoter$start
ctd_03_pr$chromEnd <- CTD_03_promoter$end
ctd_03_pr$name <- CTD_03_promoter$peak

chrom <- paste0(x = "chr", CTD_01_both$seqnames)
ctd_01_b <- as.data.frame(chrom)
ctd_01_b$chromStart <- CTD_01_both$start
ctd_01_b$chromEnd <- CTD_01_both$end
ctd_01_b$name <- CTD_01_both$peak

chrom <- paste0(x = "chr", CTD_03_both$seqnames)
ctd_03_b <- as.data.frame(chrom)
ctd_03_b$chromStart <- CTD_03_both$start
ctd_03_b$chromEnd <- CTD_03_both$end
ctd_03_b$name <- CTD_03_both$peak

output1 <- "C:/Users/Daria/Documents/CTD/broad/CTD_01_broad_promoter.bed"
output2 <- "C:/Users/Daria/Documents/CTD/broad/CTD_03_broad_promoter.bed"
output3 <- "C:/Users/Daria/Documents/CTD/broad/CTD_01_broad_both.bed"
output4 <- "C:/Users/Daria/Documents/CTD/broad/CTD_03_broad_both.bed"

write.table(data.frame(ctd_01_pr), output1, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_pr), output2, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_01_b), output3, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_b), output4, sep = '\t',  quote = F,  row.names = F, col.names = F)
















