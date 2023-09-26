ctd_deg <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_deg.txt", header=FALSE, sep = "_")
ctd_deg_plus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_deg_plus.txt", header=FALSE, sep = "_")
ctd_deg_minus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_deg_minus.txt", header=FALSE, sep = "_")
ctd_undeg <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_undeg.txt", header=FALSE, sep = "_")

ctd_switch <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_switch.txt", header=FALSE, sep = "_")

ctd_switch_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_switch_01.txt", header=FALSE, sep = "_")
ctd_switch_minus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_switch_minus_01.txt", header=FALSE, sep = "_")
ctd_switch_plus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_switch_plus_01.txt", header=FALSE, sep = "_")

ctd_switch_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_switch_03.txt", header=FALSE, sep = "_")
ctd_switch_minus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_switch_minus_03.txt", header=FALSE, sep = "_")
ctd_switch_plus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_switch_plus_03.txt", header=FALSE, sep = "_")

ctd_unswitch <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_unswitch.txt", header=FALSE, sep = "_")
ctd_unswitch_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_unswitch_01.txt", header=FALSE, sep = "_")
ctd_unswitch_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_unswitch_03.txt", header=FALSE, sep = "_")


#ctd_cpm <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_cpm.txt", header=FALSE, sep = "_")
#ctd_cpm_minus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_cpm_minus.txt", header=FALSE, sep = "_")
#ctd_cpm_plus <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_cpm_plus.txt", header=FALSE, sep = "_")
#ctd_cpm_minus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_cpm_minus_01.txt", header=FALSE, sep = "_")
#ctd_cpm_plus_01 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_cpm_plus_01.txt", header=FALSE, sep = "_")
#ctd_cpm_minus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_cpm_minus_03.txt", header=FALSE, sep = "_")
#ctd_cpm_plus_03 <- read.table("/home/danya-2003/Downloads/ctd/ctd_switched_promoters/txt/ctd_cpm_plus_03.txt", header=FALSE, sep = "_")

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

ctd_cpm$V2 <- ctd_cpm$V2 - 750
ctd_cpm_minus$V2 <- ctd_cpm_minus$V2 - 750
ctd_cpm_plus$V2 <- ctd_cpm_plus$V2 - 750


median(ctd_undeg$V3-ctd_undeg$V2)
median(ctd_switch$V3-ctd_switch$V2)
median(ctd_cpm$V3-ctd_cpm$V2)

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


ctd_cpm<-ctd_cpm[,c(1,2,3)]
ctd_cpm_minus<-ctd_cpm_minus[,c(1,2,3)]
ctd_cpm_plus<-ctd_cpm_plus[,c(1,2,3)]

write.table(data.frame(ctd_deg), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_deg.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_deg_plus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_deg_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_deg_minus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_deg_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_undeg), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_undeg.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_switch), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_switch.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_switch_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_switch_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_switch_minus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_switch_minus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_switch_plus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_switch_plus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_switch_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_switch_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_switch_minus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_switch_minus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_switch_plus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_switch_plus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_unswitch), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_unswitch.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_unswitch_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_unswitch_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_unswitch_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_unswitch_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)


write.table(data.frame(ctd_cpm), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_cpm.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_cpm_minus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_cpm_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_cpm_plus), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_cpm_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
#write.table(data.frame(ctd_cpm_minus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_cpm_minus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
#write.table(data.frame(ctd_cpm_plus_01), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_cpm_plus_01.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
#write.table(data.frame(ctd_cpm_minus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_cpm_minus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
#write.table(data.frame(ctd_cpm_plus_03), "/home/danya-2003/Downloads/ctd/ctd_switched_promoters/bed_new/ctd_cpm_plus_03.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
