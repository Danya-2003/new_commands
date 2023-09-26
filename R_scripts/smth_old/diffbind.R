library(clusterProfiler)
library("org.Hs.eg.db")
library(biomaRt)

#Retrieve the ensembl info
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
dsets = listDatasets(ensembl)
bm <- useDataset("hsapiens_gene_ensembl", mart=ensembl)


### CTD-2587H24.5

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/CTD-2587H24.5_01_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

#

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/CTD-2587H24.5_03_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

### EMX2OS

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/EMX2OS_AD02_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

#

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/EMX2OS_AD04_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

### RP11-398K22.12

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/RP11-398K22.12_03_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

#

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/RP11-398K22.12_05_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

### AC005592.2

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/AC005592.2_02_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/AC005592.2_02_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/AC005592.2_02_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/AC005592.2_02_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/AC005592.2_02_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/AC005592.2_02_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/AC005592.2_02_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/AC005592.2_02_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/AC005592.2_02_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/AC005592.2_02_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

#

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/AC005592.2_03_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/AC005592.2_03_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/AC005592.2_03_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/AC005592.2_03_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/AC005592.2_03_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/AC005592.2_03_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/AC005592.2_03_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/AC005592.2_03_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/AC005592.2_03_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/AC005592.2_03_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

### FGD5-AS1

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/FGD5-AS1_04_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

#

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/FGD5-AS1_AD03_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

### JPX

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/JPX_05_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/JPX_05_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/JPX_05_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/JPX_05_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/JPX_05_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/JPX_05_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/JPX_05_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/JPX_05_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/JPX_05_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/JPX_05_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

#

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/JPX_AD04_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/JPX_AD04_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/JPX_AD04_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/JPX_AD04_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/JPX_AD04_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/JPX_AD04_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/JPX_AD04_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/JPX_AD04_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/JPX_AD04_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/JPX_AD04_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

### MAPKAPK5-AS1

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K4me3_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_06_H3K27me3_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

#

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K4me3_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.both.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.both.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.both.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.both.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.both.kegg.csv"

source("C:/ChIP-seq/GO_example.r")

data <- "C:/ChIP-seq/diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.promoter.csv"
out_mf <- "C:/ChIP-seq/GO_KEGG/GO_MF_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.promoter.mf.csv"
out_bp <- "C:/ChIP-seq/GO_KEGG/GO_BP_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.promoter.bp.csv"
out_cc <- "C:/ChIP-seq/GO_KEGG/GO_CC_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.promoter.cc.csv"
out_kegg <- "C:/ChIP-seq/GO_KEGG/KEGG_diffbind_ChIPpeakAnno/MAPKAPK5-AS1_AD04_H3K27me3_sicer.promoter.kegg.csv"

source("C:/ChIP-seq/GO_example.r")
