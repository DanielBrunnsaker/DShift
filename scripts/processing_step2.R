library(readxl)
library(gtools)
library(data.table)
library(tidyverse)
library(rstudioapi)

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

b1 = read.csv('../data/metabolomics/intermediates/ft_qlsc_b1_batched.csv')
b2 = read.csv('../data/metabolomics/intermediates/ft_qlsc_b2.csv')

b2data = b2[!grepl("QC", b2$Class),]
b1data = b1[!grepl("QC", b1$Class),]

b1data$Class[b1data$Batch %like% 1] = 'Glucose'
b1data$Class[b1data$Batch %like% 2] = 'Glucose'
b1data$Class[b1data$Batch %like% 3] = 'Ethanol'
b1data$Class[b1data$Batch %like% 4] = 'Ethanol'
b1data$Class[b1data$Batch %like% 5] = 'Glucose'
b1data$Class[b1data$Batch %like% 6] = 'Ethanol'
b1data$Class[b1data$Batch %like% 7] = 'Glucose'

b2data$Class[b2data$Batch %like% 1] = 'Glucose'
b2data$Class[b2data$Batch %like% 2] = 'Glucose'
b2data$Class[b2data$Batch %like% 3] = 'Ethanol'
b2data$Class[b2data$Batch %like% 4] = 'Ethanol'
b2data$Class[b2data$Batch %like% 5] = 'Glucose'
b2data$Class[b2data$Batch %like% 6] = 'Ethanol'
b2data$Class[b2data$Batch %like% 7] = 'Glucose'

#remove unwanted columns
b1data$X <- substring(b1data$X, 19)
b2data$X <- substring(b2data$X, 19)

rownames(b1data) = b1data$X
rownames(b2data) = b2data$X
b1data = b1data[,-c(1,3,4)]
b2data = b2data[,-c(1,3,4)]

b1data = b1data[!grepl("cst6|B1E_dld3_3|B2E_ygr_3|B2E_tda1_1|B2E_WT_2|B2E_faa1_2|B1G|B2G_WT_3|B2G_rme1_3|B2G_ygr_1|B2G_tda1_3|B2G_gal11_3|B2G_faa1_3", rownames(b1data)),]
b2data = b2data[!grepl("PCL1_1-B6E|PCL1_2-B6E|OCA1_2-B6E|RME1_1-B6E", rownames(b2data)),]
b2data = b2data[!grepl("OCA1_1-B5G|GAL11_2-B5G|MEK1_2-B5G|RTS3_3-B5G|WT_1-B5G|RME1_1-B5G", rownames(b2data)),]
b2data = b2data[!grepl("GAL11_1-B7G|DLD3_1-B7G|PCL1_1-B7G|MEK1_1-B7G", rownames(b2data)),]

#log2 before release
b1data[-1] = log2(b1data[-1])
b2data[-1] = log2(b2data[-1])

write.csv(b1data,'../data/metabolomics/intermediates/b1data.csv')
write.csv(b2data,'../data/metabolomics/intermediates/b2data.csv')

