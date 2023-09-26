
ctd_01_high <- read.table("C:/Users/Daria//Downloads/expr_bed_finale/ctd_01_high.bed", header =FALSE, sep = '\t')
ctd_01_low <- read.table("C:/Users/Daria//Downloads/expr_bed_finale/ctd_01_low.bed", header =FALSE, sep = '\t')
ctd_01_zero <- read.table("C:/Users/Daria//Downloads/expr_bed_finale/ctd_01_zero.bed", header =FALSE, sep = '\t')
ctd_03_high <- read.table("C:/Users/Daria/Downloads/expr_bed_finale/ctd_03_high.bed", header =FALSE, sep = '\t')
ctd_03_low <- read.table("C:/Users/Daria/Downloads/expr_bed_finale/ctd_03_low.bed", header =FALSE, sep = '\t')
ctd_03_zero <- read.table("C:/Users/Daria/Downloads/expr_bed_finale/ctd_03_zero.bed", header =FALSE, sep = '\t')

ctd_01_high$start <- ctd_01_high$V2 - 750
ctd_01_low$start <- ctd_01_low$V2 - 750
ctd_01_zero$start <- ctd_01_zero$V2 - 750
ctd_03_high$start <- ctd_03_high$V2 - 750
ctd_03_low$start <- ctd_03_low$V2 - 750
ctd_03_zero$start <- ctd_03_zero$V2 - 750

ctd_01_high$end <- ctd_01_high$V2 + 2000
ctd_01_low$end <- ctd_01_low$V2 + 2000
ctd_01_zero$end <- ctd_01_zero$V2 + 2000
ctd_03_high$end <- ctd_03_high$V2 + 2000
ctd_03_low$end <- ctd_03_low$V2 + 2000
ctd_03_zero$end <- ctd_03_zero$V2 + 2000

ctd_01_high <- ctd_01_high[,c(1,7,8,4,5,6)]
ctd_01_low <- ctd_01_low[,c(1,7,8,4,5,6)]
ctd_01_zero <- ctd_01_zero[,c(1,7,8,4,5,6)]
ctd_03_high <- ctd_03_high[,c(1,7,8,4,5,6)]
ctd_03_low <- ctd_03_low[,c(1,7,8,4,5,6)]
ctd_03_zero <- ctd_03_zero[,c(1,7,8,4,5,6)]

ctd_01_high[ctd_01_high$start < 0,] <- NA
ctd_01_low[ctd_01_low$start < 0,] <- NA
ctd_01_zero[ctd_01_zero$start < 0,] <- NA
ctd_03_high[ctd_03_high$start < 0,] <- NA
ctd_03_low[ctd_03_low$start < 0,] <- NA
ctd_03_zero[ctd_03_zero$start < 0,] <- NA

ctd_01_high <- na.omit(ctd_01_high)
ctd_01_low <- na.omit(ctd_01_low)
ctd_01_zero <- na.omit(ctd_01_zero)
ctd_03_high <- na.omit(ctd_03_high)
ctd_03_low <- na.omit(ctd_03_low)
ctd_03_zero <- na.omit(ctd_03_zero)

write.table(data.frame(ctd_01_high), "C:/Users/Daria//Downloads/ctd_promoters/ctd_01_high_promoters.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_01_low), "C:/Users/Daria//Downloads/ctd_promoters/ctd_01_low_promoters.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_01_zero), "C:/Users/Daria//Downloads/ctd_promoters/ctd_01_zero_promoters.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_high), "C:/Users/Daria//Downloads/ctd_promoters/ctd_03_high_promoters.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_low), "C:/Users/Daria//Downloads/ctd_promoters/ctd_03_low_promoters.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_zero), "C:/Users/Daria//Downloads/ctd_promoters/ctd_03_zero_promoters.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

ctd_deg <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_deg.txt", header=FALSE, sep = "_")
ctd_deg_plus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_deg_plus.txt", header=FALSE, sep = "_")
ctd_deg_minus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_deg_minus.txt", header=FALSE, sep = "_")
ctd_undeg <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_undeg.txt", header=FALSE, sep = "_")

ctd_switch <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_switch.txt", header=FALSE, sep = "_")

ctd_switch_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_switch_01.txt", header=FALSE, sep = "_")
ctd_switch_minus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_switch_minus_01.txt", header=FALSE, sep = "_")
ctd_switch_plus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_switch_plus_01.txt", header=FALSE, sep = "_")

ctd_switch_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_switch_03.txt", header=FALSE, sep = "_")
ctd_switch_minus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_switch_minus_03.txt", header=FALSE, sep = "_")
ctd_switch_plus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_switch_plus_03.txt", header=FALSE, sep = "_")

ctd_unswitch <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_unswitch.txt", header=FALSE, sep = "_")
ctd_unswitch_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_unswitch_01.txt", header=FALSE, sep = "_")
ctd_unswitch_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_unswitch_03.txt", header=FALSE, sep = "_")

ctd_result <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result.txt", header=FALSE, sep = "_")
ctd_result_minus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result_minus.txt", header=FALSE, sep = "_")
ctd_result_plus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result_plus.txt", header=FALSE, sep = "_")
ctd_result_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result_01.txt", header=FALSE, sep = "_")
ctd_result_minus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result_minus_01.txt", header=FALSE, sep = "_")
ctd_result_plus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result_plus_01.txt", header=FALSE, sep = "_")
ctd_result_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result_03.txt", header=FALSE, sep = "_")
ctd_result_minus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result_minus_03.txt", header=FALSE, sep = "_")
ctd_result_plus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_result_plus_03.txt", header=FALSE, sep = "_")

ctd_cpm <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_cpm.txt", header=FALSE, sep = "_")
ctd_cpm_minus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_cpm_minus.txt", header=FALSE, sep = "_")
ctd_cpm_plus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_cpm_plus.txt", header=FALSE, sep = "_")
#ctd_cpm_minus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_cpm_minus_01.txt", header=FALSE, sep = "_")
#ctd_cpm_plus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_cpm_plus_01.txt", header=FALSE, sep = "_")
#ctd_cpm_minus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_cpm_minus_03.txt", header=FALSE, sep = "_")
#ctd_cpm_plus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/getfasta/ctd_cpm_plus_03.txt", header=FALSE, sep = "_")

ctd_deg$V2 <- ctd_deg$V2 - 750
ctd_deg_plus$V2 <- ctd_deg_plus$V2 - 750
ctd_deg_minus$V2 <- ctd_deg_minus$V2 - 750
ctd_undeg$V2 <- ctd_undeg$V2 - 750

ctd_switch$V2 <- ctd_switch$V2 - 750

ctd_switch_01$V2 <- ctd_switch_01$V2 - 750
ctd_switch_minus_01$V2 <- ctd_switch_minus_01$V2 - 750
ctd_switch_plus_01$V2 <- ctd_switch_plus_01$V2 - 750

ctd_switch_03$V2 <- ctd_switch_03$V2 - 750
ctd_switch_minus_03$V2 <- ctd_switch_minus_03$V2 - 750
ctd_switch_plus_03$V2 <- ctd_switch_plus_03$V2 - 750

ctd_unswitch$V2 <- ctd_unswitch$V2 - 750
ctd_unswitch_01$V2 <- ctd_unswitch_01$V2 - 750
ctd_unswitch_03$V2 <- ctd_unswitch_03$V2 - 750

ctd_result$V2 <- ctd_result$V2 - 750
ctd_result_minus$V2 <- ctd_result_minus$V2 - 750
ctd_result_plus$V2 <- ctd_result_plus$V2 - 750
ctd_result_01$V2 <- ctd_result_01$V2 - 750
ctd_result_minus_01$V2 <- ctd_result_minus_01$V2 - 750
ctd_result_plus_01$V2 <- ctd_result_plus_01$V2 - 750
ctd_result_03$V2 <- ctd_result_03$V2 - 750
ctd_result_minus_03$V2 <- ctd_result_minus_03$V2 - 750
ctd_result_plus_03$V2 <- ctd_result_plus_03$V2 - 750

ctd_cpm$V2 <- ctd_cpm$V2 - 750
ctd_cpm_minus$V2 <- ctd_cpm_minus$V2 - 750
ctd_cpm_plus$V2 <- ctd_cpm_plus$V2 - 750


mean(ctd_undeg$V3-ctd_undeg$V2)
mean(ctd_switch$V3-ctd_switch$V2)
mean(ctd_cpm$V3-ctd_cpm$V2)

ctd_deg$V3 <- ctd_deg$V3 + 0
ctd_deg_plus$V3 <- ctd_deg_plus$V3 + 0
ctd_deg_minus$V3 <- ctd_deg_minus$V3 + 0
ctd_undeg$V3 <- ctd_undeg$V3 + 0

ctd_switch$V3 <- ctd_switch$V3 + 0

ctd_switch_01$V3 <- ctd_switch_01$V3 + 0
ctd_switch_minus_01$V3 <- ctd_switch_minus_01$V3 + 0
ctd_switch_plus_01$V3 <- ctd_switch_plus_01$V3 + 0

ctd_switch_03$V3 <- ctd_switch_03$V3 + 0
ctd_switch_minus_03$V3 <- ctd_switch_minus_03$V3 + 0
ctd_switch_plus_03$V3 <- ctd_switch_plus_03$V3 + 0

ctd_unswitch$V3 <- ctd_unswitch$V3 + 0
ctd_unswitch_01$V3 <- ctd_unswitch_01$V3 + 0
ctd_unswitch_03$V3 <- ctd_unswitch_03$V3 + 0

ctd_result$V3 <- ctd_result$V3 + 0
ctd_result_minus$V3 <- ctd_result_minus$V3 + 0
ctd_result_plus$V3 <- ctd_result_plus$V3 + 0
ctd_result_01$V3 <- ctd_result_01$V3 + 0
ctd_result_minus_01$V3 <- ctd_result_minus_01$V3 + 0
ctd_result_plus_01$V3 <- ctd_result_plus_01$V3 + 0
ctd_result_03$V3 <- ctd_result_03$V3 + 0
ctd_result_minus_03$V3 <- ctd_result_minus_03$V3 + 0
ctd_result_plus_03$V3 <- ctd_result_plus_03$V3 + 0

ctd_cpm$V3 <- ctd_cpm$V3 + 0
ctd_cpm_minus$V3 <- ctd_cpm_minus$V3 + 0
ctd_cpm_plus$V3 <- ctd_cpm_plus$V3 + 0


ctd_deg$V2[ ctd_deg$V2 < 0] <- 0
ctd_deg_plus$V2[ ctd_deg_plus$V2 < 0] <- 0
ctd_deg_minus$V2[ ctd_deg_minus$V2 < 0] <- 0
ctd_undeg$V2[ ctd_undeg$V2 < 0] <- 0

ctd_switch$V2[ ctd_switch$V2 < 0] <- 0

ctd_switch_01$V2[ ctd_switch_01$V2 < 0] <- 0
ctd_switch_minus_01$V2[ ctd_switch_minus_01$V2 < 0] <- 0
ctd_switch_plus_01$V2[ ctd_switch_plus_01$V2 < 0] <- 0

ctd_switch_03$V2[ ctd_switch_03$V2 < 0] <- 0
ctd_switch_minus_03$V2[ ctd_switch_minus_03$V2 < 0] <- 0
ctd_switch_plus_03$V2[ ctd_switch_plus_03$V2 < 0] <- 0

ctd_unswitch$V2[ ctd_unswitch$V2 < 0] <- 0
ctd_unswitch_01$V2[ ctd_unswitch_01$V2 < 0] <- 0
ctd_unswitch_03$V2[ ctd_unswitch_03$V2 < 0] <- 0

ctd_result$V2[ ctd_result$V2 < 0] <- 0
ctd_result_minus$V2[ ctd_result_minus$V2 < 0] <- 0
ctd_result_plus$V2[ ctd_result_plus$V2 < 0] <- 0
ctd_result_01$V2[ ctd_result_01$V2 < 0] <- 0
ctd_result_minus_01$V2[ ctd_result_minus_01$V2 < 0] <- 0
ctd_result_plus_01$V2[ ctd_result_plus_01$V2 < 0] <- 0
ctd_result_03$V2[ ctd_result_03$V2 < 0] <- 0
ctd_result_minus_03$V2[ ctd_result_minus_03$V2 < 0] <- 0
ctd_result_plus_03$V2[ ctd_result_plus_03$V2 < 0] <- 0

ctd_cpm$V2[ ctd_cpm$V2 < 0] <- 0
ctd_cpm_minus$V2[ ctd_cpm_minus$V2 < 0] <- 0
ctd_cpm_plus$V2[ ctd_cpm_plus$V2 < 0] <- 0


ctd_deg<-ctd_deg[,c(1,2,3)]
ctd_deg_plus<-ctd_deg_plus[,c(1,2,3)]
ctd_deg_minus<-ctd_deg_minus[,c(1,2,3)]
ctd_undeg<-ctd_undeg[,c(1,2,3)]

ctd_switch<-ctd_switch[,c(1,2,3)]

ctd_switch_01<-ctd_switch_01[,c(1,2,3)]
ctd_switch_minus_01<-ctd_switch_minus_01[,c(1,2,3)]
ctd_switch_plus_01<-ctd_switch_plus_01[,c(1,2,3)]

ctd_switch_03<-ctd_switch_03[,c(1,2,3)]
ctd_switch_minus_03<-ctd_switch_minus_03[,c(1,2,3)]
ctd_switch_plus_03<-ctd_switch_plus_03[,c(1,2,3)]

ctd_unswitch<-ctd_unswitch[,c(1,2,3)]
ctd_unswitch_01<-ctd_unswitch_01[,c(1,2,3)]
ctd_unswitch_03<-ctd_unswitch_03[,c(1,2,3)]

ctd_result<-ctd_result[,c(1,2,3)]
ctd_result_minus<-ctd_result_minus[,c(1,2,3)]
ctd_result_plus<-ctd_result_plus[,c(1,2,3)]
ctd_result_01<-ctd_result_01[,c(1,2,3)]
ctd_result_minus_01<-ctd_result_minus_01[,c(1,2,3)]
ctd_result_plus_01<-ctd_result_plus_01[,c(1,2,3)]
ctd_result_03<-ctd_result_03[,c(1,2,3)]
ctd_result_minus_03<-ctd_result_minus_03[,c(1,2,3)]
ctd_result_plus_03<-ctd_result_plus_03[,c(1,2,3)]

ctd_cpm<-ctd_cpm[,c(1,2,3)]
ctd_cpm_minus<-ctd_cpm_minus[,c(1,2,3)]
ctd_cpm_plus<-ctd_cpm_plus[,c(1,2,3)]

write.table(data.frame(ctd_deg), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_deg.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_deg_plus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_deg_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_deg_minus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_deg_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_undeg), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_undeg.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_switch), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_switch.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_switch_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_switch_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_switch_minus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_switch_minus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_switch_plus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_switch_plus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_switch_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_switch_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_switch_minus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_switch_minus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_switch_plus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_switch_plus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_unswitch), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_unswitch.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_unswitch_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_unswitch_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_unswitch_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_unswitch_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_result), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_minus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_plus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_minus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result_minus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_plus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result_plus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_minus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result_minus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_plus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_result_plus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_cpm), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_cpm.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_cpm_minus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_cpm_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_cpm_plus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_cpm_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
#write.table(data.frame(ctd_cpm_minus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_cpm_minus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
#write.table(data.frame(ctd_cpm_plus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_cpm_plus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
#write.table(data.frame(ctd_cpm_minus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_cpm_minus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
#write.table(data.frame(ctd_cpm_plus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed/ctd_cpm_plus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
