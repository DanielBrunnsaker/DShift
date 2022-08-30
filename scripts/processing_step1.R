library(gtools)
library(NormalizeMets)
library(readxl)
library(tidyverse)
library(dplyr)
library(data.table)
library(stringr)
library(impute)
library(limma)
library(rstudioapi)

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

source("supp/processing_supp.R")

src = read_excel('../data/metabolomics/B1_filtered.xlsx')
ft = src[-c(seq(1,33))]

#create naming vector
naming_col = paste(src$`Average Rt(min)`,src$`Average Mz`,src$`Adduct type`, src$`Metabolite name`,sep = '__')
ft = data.frame(ft[!(duplicated(naming_col)),])
naming_col = naming_col[!duplicated(naming_col)]
rownames(ft) = naming_col

st = data.frame(colnames(ft))
st$Class = 0

#class
st$Class[st$colnames.ft. %like% 'tda1'] = 'TDA1'
st$Class[st$colnames.ft. %like% 'rts3'] = 'RTS3'
st$Class[st$colnames.ft. %like% 'rst3'] = 'RTS3'
st$Class[st$colnames.ft. %like% 'dld3'] = 'DLD3'
st$Class[st$colnames.ft. %like% 'pcl1'] = 'PCL1'
st$Class[st$colnames.ft. %like% 'ygr'] = 'YGR067C'
st$Class[st$colnames.ft. %like% 'rme1'] = 'RME1'
st$Class[st$colnames.ft. %like% 'oca1'] = 'OCA1'
st$Class[st$colnames.ft. %like% 'mek1'] = 'MEK1'
st$Class[st$colnames.ft. %like% 'gal11'] = 'GAL11'
st$Class[st$colnames.ft. %like% 'faa1'] = 'FAA1'
st$Class[st$colnames.ft. %like% 'cst6'] = 'CST6'
st$Class[st$colnames.ft. %like% 'WT'] = 'WT'

#batch
st$Batch = 0
st$Batch[st$colnames.ft. %like% 'B1G'] = 1
st$Batch[st$colnames.ft. %like% 'B2G'] = 2
st$Batch[st$colnames.ft. %like% 'B1E'] = 3
st$Batch[st$colnames.ft. %like% 'B2E'] = 4

#remove blanks and some qcs
ft = ft[,!grepl("blank|ev", st$colnames.ft.)]
st = st[!grepl("blank|ev", st$colnames.ft.),]

#order
st$Order = 0
st$Order = st$colnames.ft. 
st$Order = as.numeric(str_split_fixed(st$colnames.ft., "_pos_", 2)[,2])
ft = t(ft)
ft <- ft[mixedorder(st$Order,),];
st <- st[mixedorder(st$Order,),];
st$Order = seq(1,dim(st)[1])
colnames(st)[1] = 'Samples'
rownames(st) = st$Samples

#replace 0 with half of min
ft[ft == 0] <- min(ft[ft != 0])*0.5

#split up into minibatches to accomodate experimental design
ft1s = 20
ft2s = 28
ft3s = 25
ft4s = 25

ft_1 = ft[c(which(st$Batch == 1),ft1s+1),]
st_1 = st[c(which(st$Batch == 1),ft1s+1),]

ft_2 = ft[c(which(st$Batch == 2),ft2s+1),]
st_2 = st[c(which(st$Batch == 2),ft2s+1),]
st_2$Order = seq(1,dim(st_2)[1])

ft_3 = ft[c(which(st$Batch == 3)),]
st_3 = st[c(which(st$Batch == 3)),]
st_3$Order = seq(1,dim(st_3)[1])

ft_4 = ft[st$Batch == 4,]
st_4 = st[st$Batch == 4,]
st_4$Order = seq(1,dim(st_4)[1])

ft_corrected_1 = NormQcsamples(data.frame(ft_1), data.frame(st_1), method = c("rlsc"), lg = FALSE, saveoutput = FALSE, span = 0.75)
ft_corrected_2 = NormQcsamples(data.frame(ft_2), data.frame(st_2), method = c("rlsc"), lg = FALSE, saveoutput = FALSE, span = 0.75)
ft_corrected_3 = NormQcsamples(data.frame(ft_3), data.frame(st_3), method = c("rlsc"), lg = FALSE, saveoutput = FALSE, span = 0.75)
ft_corrected_4 = NormQcsamples(data.frame(ft_4), data.frame(st_4), method = c("rlsc"), lg = FALSE, saveoutput = FALSE, span = 0.75)

st[st == 0]<-'QC'
st_1[st_1 == 0]<-'QC'
st_2[st_2 == 0]<-'QC'
st_3[st_3 == 0]<-'QC'
st_4[st_4 == 0]<-'QC'

ft_annotated_1 = cbind(st_1[,c(2,3,4)],ft_corrected_1$featuredata)
ft_annotated_2 = cbind(st_2[,c(2,3,4)],ft_corrected_2$featuredata)
ft_annotated_3 = cbind(st_3[,c(2,3,4)],ft_corrected_3$featuredata)
ft_annotated_4 = cbind(st_4[,c(2,3,4)],ft_corrected_4$featuredata)

#bind them together
ft_complete = rbind(ft_annotated_1,ft_annotated_2)
ft_complete = rbind(ft_complete,ft_annotated_3)
ft_complete = rbind(ft_complete,ft_annotated_4)

write.csv(ft_complete,"../data/metabolomics/intermediates/ft_qlsc_b1_batched.csv")

##########################################

src = read_excel('../data/metabolomics/B2_filtered.xlsx')
ft = src[-c(seq(1,33))]

naming_col = paste(src$`Average Rt(min)`,src$`Average Mz`,src$`Adduct type`, src$`Metabolite name`,sep = '__')
ft = data.frame(ft[!(duplicated(naming_col)),])
naming_col = naming_col[!duplicated(naming_col)]
rownames(ft) = naming_col

st = data.frame(colnames(ft))
st$Class = 0

#class
st$Class[st$colnames.ft. %like% 'TDA1'] = 'TDA1'
st$Class[st$colnames.ft. %like% 'RTS3'] = 'RTS3'
st$Class[st$colnames.ft. %like% 'RST3'] = 'RTS3'
st$Class[st$colnames.ft. %like% 'DLD3'] = 'DLD3'
st$Class[st$colnames.ft. %like% 'PCL1'] = 'PCL1'
st$Class[st$colnames.ft. %like% 'YGR'] = 'YGR067C'
st$Class[st$colnames.ft. %like% 'RME1'] = 'RME1'
st$Class[st$colnames.ft. %like% 'OCA1'] = 'OCA1'
st$Class[st$colnames.ft. %like% 'MEK1'] = 'MEK1'
st$Class[st$colnames.ft. %like% 'GAL11'] = 'GAL11'
st$Class[st$colnames.ft. %like% 'FAA1'] = 'FAA1'
st$Class[st$colnames.ft. %like% 'CST6'] = 'CST6'
st$Class[st$colnames.ft. %like% 'WT'] = 'WT'

#batch
st$Batch = 0
st$Batch[st$colnames.ft. %like% 'B5G'] = 5
st$Batch[st$colnames.ft. %like% 'B6E'] = 6
st$Batch[st$colnames.ft. %like% 'B7G'] = 7

#remove blanks and some qcs
ft = ft[,!grepl("blank|MeOH|ev|Media|QC_MS2-|RP_POS_B6E|RP_POS_B5G|RP_POS_B7G", st$colnames.ft.)]
st = st[!grepl("blank|MeOH|ev|Media|QC_MS2-|RP_POS_B6E|RP_POS_B5G|RP_POS_B7G", st$colnames.ft.),]

#order
st$Order = 0
st$Order = as.numeric(str_sub(st$colnames.ft.,-3,-1))

ft = t(ft)
ft <- ft[mixedorder(st$Order,),];
st <- st[mixedorder(st$Order,),];
st$Order = seq(1,dim(st)[1])
colnames(st)[1] = 'Samples'
rownames(st) = st$Samples

#replace 0 with half of min
ft[ft == 0] <- min(ft[ft != 0])*0.5

ft_corrected = NormQcsamples(data.frame(ft), data.frame(st), method = c("rlsc"), lg = FALSE, saveoutput = FALSE, span = 0.75)

st[st == 0]<-'QC'
ft_annotated = cbind(st[,c(2,3,4)],ft_corrected$featuredata)
write.csv(ft_annotated,"../data/metabolomics/intermediates/ft_qlsc_b2.csv")
