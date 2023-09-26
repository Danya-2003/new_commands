ctd_01_broad <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_01_broad.filtered", header =FALSE, sep = '\t')
ctd_01_narrow <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_01_narrow.filtered", header =FALSE, sep = '\t')
ctd_01_sicer <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_01_sicer.filtered", header =FALSE, sep = '\t')
ctd_03_broad <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_03_broad.filtered", header =FALSE, sep = '\t')
ctd_03_narrow <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_03_narrow.filtered", header =FALSE, sep = '\t')
ctd_03_sicer <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_03_sicer.filtered", header =FALSE, sep = '\t')


ctd_01_b <- as.data.frame(ctd_01_broad$V1)
ctd_01_b$V2 <- ctd_01_broad$V2
ctd_01_b$V3 <- ctd_01_broad$V3
ctd_01_b$V4 <- ctd_01_broad$V4

ctd_01_n <- as.data.frame(ctd_01_narrow$V1)
ctd_01_n$V2 <- ctd_01_narrow$V2
ctd_01_n$V3 <- ctd_01_narrow$V3
ctd_01_n$V4 <- ctd_01_narrow$V4

ctd_01_s <- as.data.frame(ctd_01_sicer$V1)
ctd_01_s$V2 <- ctd_01_sicer$V2
ctd_01_s$V3 <- ctd_01_sicer$V3
ctd_01_s$V4 <- ctd_01_sicer$V4

ctd_03_b <- as.data.frame(ctd_03_broad$V1)
ctd_03_b$V2 <- ctd_03_broad$V2
ctd_03_b$V3 <- ctd_03_broad$V3
ctd_03_b$V4 <- ctd_03_broad$V4

ctd_03_n <- as.data.frame(ctd_03_narrow$V1)
ctd_03_n$V2 <- ctd_03_narrow$V2
ctd_03_n$V3 <- ctd_03_narrow$V3
ctd_03_n$V4 <- ctd_03_narrow$V4

ctd_03_s <- as.data.frame(ctd_03_sicer$V1)
ctd_03_s$V2 <- ctd_03_sicer$V2
ctd_03_s$V3 <- ctd_03_sicer$V3
ctd_03_s$V4 <- ctd_03_sicer$V4

output1 <- ("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_01_broad.bed")
output2 <- ("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_01_narrow.bed")
output3 <- ("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_01_sicer.bed")
output4 <- ("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_03_broad.bed")
output5 <- ("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_03_narrow.bed")
output6 <- ("C:/Users/Daria/Documents/CTD/peaks/KD/ctd_03_sicer.bed")

write.table(data.frame(ctd_01_b), output1, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_01_n), output2, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_01_s), output3, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_b), output4, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_n), output5, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(ctd_03_s), output6, sep = '\t',  quote = F,  row.names = F, col.names = F)


####################################################################################################################################################


nca_1_broad <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_1_broad.filtered", header =FALSE, sep = '\t')
nca_1_narrow <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_1_narrow.filtered", header =FALSE, sep = '\t')
nca_1_sicer <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_1_sicer.filtered", header =FALSE, sep = '\t')
nca_2_broad <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_2_broad.filtered", header =FALSE, sep = '\t')
nca_2_narrow <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_2_narrow.filtered", header =FALSE, sep = '\t')
nca_2_sicer <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_2_sicer.filtered", header =FALSE, sep = '\t')

nca_1_b <- as.data.frame(nca_1_broad$V1)
nca_1_b$V2 <- nca_1_broad$V2
nca_1_b$V3 <- nca_1_broad$V3
nca_1_b$V4 <- nca_1_broad$V4

nca_1_n <- as.data.frame(nca_1_narrow$V1)
nca_1_n$V2 <- nca_1_narrow$V2
nca_1_n$V3 <- nca_1_narrow$V3
nca_1_n$V4 <- nca_1_narrow$V4

nca_1_s <- as.data.frame(nca_1_sicer$V1)
nca_1_s$V2 <- nca_1_sicer$V2
nca_1_s$V3 <- nca_1_sicer$V3
nca_1_s$V4 <- nca_1_sicer$V4

nca_2_b <- as.data.frame(nca_2_broad$V1)
nca_2_b$V2 <- nca_2_broad$V2
nca_2_b$V3 <- nca_2_broad$V3
nca_2_b$V4 <- nca_2_broad$V4

nca_2_n <- as.data.frame(nca_2_narrow$V1)
nca_2_n$V2 <- nca_2_narrow$V2
nca_2_n$V3 <- nca_2_narrow$V3
nca_2_n$V4 <- nca_2_narrow$V4

nca_2_s <- as.data.frame(nca_2_sicer$V1)
nca_2_s$V2 <- nca_2_sicer$V2
nca_2_s$V3 <- nca_2_sicer$V3
nca_2_s$V4 <- nca_2_sicer$V4

output1 <- ("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_1_broad.bed")
output2 <- ("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_1_narrow.bed")
output3 <- ("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_1_sicer.bed")
output4 <- ("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_2_broad.bed")
output5 <- ("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_2_narrow.bed")
output6 <- ("C:/Users/Daria/Documents/CTD/peaks/NCA/nca_2_sicer.bed")

write.table(data.frame(nca_1_b), output1, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(nca_1_n), output2, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(nca_1_s), output3, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(nca_2_b), output4, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(nca_2_n), output5, sep = '\t',  quote = F,  row.names = F, col.names = F)
write.table(data.frame(nca_2_s), output6, sep = '\t',  quote = F,  row.names = F, col.names = F)

###################################################################################################################################################

ctd_01_broad <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/results/f-ctd_01_broad.bed", fill = TRUE, header =FALSE)
ctd_01_narrow <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/results/f-ctd_01_narrow.bed", fill = TRUE, header =FALSE)
ctd_01_sicer <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/results/f-ctd_01_sicer.bed", fill = TRUE, header =FALSE)
ctd_03_broad <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/results/f-ctd_03_broad.bed", fill = TRUE, header =FALSE)
ctd_03_narrow <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/results/f-ctd_03_narrow.bed", fill = TRUE, header =FALSE)
ctd_03_sicer <- read.table("C:/Users/Daria/Documents/CTD/peaks/KD/results/f-ctd_03_sicer.bed", fill = TRUE, header =FALSE)

ctd_01_broad <- ctd_01_broad$V1
n <- length(ctd_01_broad)
x <- seq(3, n, by = 3)
ctd_01_b <- as.numeric(ctd_01_broad[x])
ctd_01_b <- ctd_01_b[ctd_01_b < 10000 & ctd_01_b != -1]
hist(ctd_01_b)


ctd_01_narrow <- ctd_01_narrow$V1
n <- length(ctd_01_narrow)
x <- seq(3, n, by = 3)
ctd_01_n <- as.numeric(ctd_01_narrow[x])
ctd_01_n <- ctd_01_n[ctd_01_n < 10000 & ctd_01_n != -1]
hist(ctd_01_n)

ctd_01_sicer <- ctd_01_sicer$V1
n <- length(ctd_01_sicer)
x <- seq(3, n, by = 3)
ctd_01_s <- as.numeric(ctd_01_sicer[x])
ctd_01_s <- ctd_01_s[ctd_01_s < 10000 & ctd_01_s != -1]
hist(ctd_01_s)


ctd_03_broad <- ctd_03_broad$V1
n <- length(ctd_03_broad)
x <- seq(3, n, by = 3)
ctd_03_b <- as.numeric(ctd_03_broad[x])
ctd_03_b <- ctd_03_b[ctd_03_b < 10000 & ctd_03_b != -1]
hist(ctd_03_b)

ctd_03_narrow <- ctd_03_narrow$V1
n <- length(ctd_03_narrow)
x <- seq(3, n, by = 3)
ctd_03_n <- as.numeric(ctd_03_narrow[x])
ctd_03_n <- ctd_03_n[ctd_03_n < 10000 & ctd_03_n != -1]
hist(ctd_03_n)

ctd_03_sicer <- ctd_03_sicer$V1
n <- length(ctd_03_sicer)
x <- seq(3, n, by = 3)
ctd_03_s <- as.numeric(ctd_03_sicer[x])
ctd_03_s <- ctd_03_s[ctd_03_s < 10000 & ctd_03_s != -1]
hist(ctd_03_s)

###################################################################################################################################################

nca_1_broad <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/results/f-nca_1_broad.bed", header =FALSE, fill = TRUE)
nca_1_narrow <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/results/f-nca_1_narrow.bed", header =FALSE, fill = TRUE)
nca_1_sicer <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/results/f-nca_1_sicer.bed", header =FALSE, fill = TRUE)
nca_2_broad <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/results/f-nca_2_broad.bed", header =FALSE, fill = TRUE)
nca_2_narrow <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/results/f-nca_2_narrow.bed", header =FALSE, fill = TRUE)
nca_2_sicer <- read.table("C:/Users/Daria/Documents/CTD/peaks/NCA/results/f-nca_2_sicer.bed", header =FALSE, fill = TRUE)

nca_1_broad <- nca_1_broad$V1
n <- length(nca_1_broad)
x <- seq(3, n, by = 3)
nca_1_b <- as.numeric(nca_1_broad[x])
nca_1_b <- nca_1_b[nca_1_b < 10000 & nca_1_b != -1]
hist(nca_1_b)

nca_1_narrow <- nca_1_narrow$V1
n <- length(nca_1_narrow)
x <- seq(3, n, by = 3)
nca_1_n <- as.numeric(nca_1_narrow[x])
nca_1_n <- nca_1_n[nca_1_n < 10000 & nca_1_n != -1]
hist(nca_1_n)

nca_1_sicer <- nca_1_sicer$V1
n <- length(nca_1_sicer)
x <- seq(3, n, by = 3)
nca_1_s <- as.numeric(nca_1_sicer[x])
nca_1_s <- nca_1_s[nca_1_s < 10000 & nca_1_s != -1]
hist(nca_1_s)


nca_2_broad <- nca_2_broad$V1
n <- length(nca_2_broad)
x <- seq(3, n, by = 3)
nca_2_b <- as.numeric(nca_2_broad[x])
nca_2_b <- nca_2_b[nca_2_b < 10000 & nca_2_b != -1]
hist(nca_2_b)

nca_2_narrow <- nca_2_narrow$V1
n <- length(nca_2_narrow)
x <- seq(3, n, by = 3)
nca_2_n <- as.numeric(nca_2_narrow[x])
nca_2_n <- nca_2_n[nca_2_n < 10000 & nca_2_n != -1]
hist(nca_2_n)

nca_2_sicer <- nca_2_sicer$V1
n <- length(nca_2_sicer)
x <- seq(3, n, by = 3)
nca_2_s <- as.numeric(nca_2_sicer[x])
nca_2_s <- nca_2_s[nca_2_s < 10000 & nca_2_s != -1]
hist(nca_2_s)
nca_2_s








