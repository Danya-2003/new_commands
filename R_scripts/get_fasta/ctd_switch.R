ctd_result_01_minus <- read.table("~/Downloads/getfasta/ctd_result_minus_01.txt", header=FALSE, sep = "_")
ctd_result_03_minus <- read.table("~/Downloads/getfasta/ctd_result_minus_03.txt", header=FALSE, sep = "_")
ctd_result_01_plus <- read.table("~/Downloads/getfasta/ctd_result_plus_01.txt", header=FALSE, sep = "_")
ctd_result_03_plus <- read.table("~/Downloads/getfasta/ctd_result_plus_03.txt", header=FALSE, sep = "_")

ctd_cpm_01_minus <- read.table("~/Downloads/getfasta/ctd_cpm_minus_01.txt", header=FALSE, sep = "_")
ctd_cpm_03_minus <- read.table("~/Downloads/getfasta/ctd_cpm_minus_03.txt", header=FALSE, sep = "_")
ctd_cpm_01_plus <- read.table("~/Downloads/getfasta/ctd_cpm_plus_01.txt", header=FALSE, sep = "_")
ctd_cpm_03_plus <- read.table("~/Downloads/getfasta/ctd_cpm_plus_03.txt", header=FALSE, sep = "_")

ctd_cpm_minus <- read.table("~/Downloads/getfasta/ctd_cpm_minus.txt", header=FALSE, sep = "_")
ctd_cpm_plus <- read.table("~/Downloads/getfasta/ctd_cpm_plus.txt", header=FALSE, sep = "_")
ctd_result_minus <- read.table("~/Downloads/getfasta/ctd_result_minus.txt", header=FALSE, sep = "_")
ctd_result_plus <- read.table("~/Downloads/getfasta/ctd_result_plus.txt", header=FALSE, sep = "_")

write.table(data.frame(ctd_result_01_minus), "~/Downloads/ctd_switch/ctd_result_01_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_03_minus), "~/Downloads/ctd_switch/ctd_result_03_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_01_plus), "~/Downloads/ctd_switch/ctd_result_01_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_03_plus), "~/Downloads/ctd_switch/ctd_result_03_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_cpm_01_minus), "~/Downloads/ctd_switch/ctd_cpm_01_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_cpm_03_minus), "~/Downloads/ctd_switch/ctd_cpm_03_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_cpm_01_plus), "~/Downloads/ctd_switch/ctd_cpm_01_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_cpm_03_plus), "~/Downloads/ctd_switch/ctd_cpm_03_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)

write.table(data.frame(ctd_cpm_minus), "~/Downloads/ctd_switch/ctd_cpm_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_cpm_plus), "~/Downloads/ctd_switch/ctd_cpm_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_minus), "~/Downloads/ctd_switch/ctd_result_minus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_result_plus), "~/Downloads/ctd_switch/ctd_result_plus.bed", sep = '\t',  quote = F,  row.names = F, col.names = F)
