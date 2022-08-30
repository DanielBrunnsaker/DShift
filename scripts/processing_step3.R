library(data.table)
library(dplyr)
library(tidyverse)
library(rstudioapi)

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

source("supp/processing_supp.R")

#load in data
b1_log = read.csv('../data/metabolomics/intermediates/b1data.csv',row.names = 1, header= TRUE)
b2_log = read.csv('../data/metabolomics/intermediates/b2data.csv',row.names = 1, header= TRUE)

src_b1 = read_excel('../data/metabolomics/B1_filtered.xlsx')
src_b2 = read_excel('../data/metabolomics/B2_filtered.xlsx')

b1_filtered = b1_log[,!grepl("Unknown", colnames(b1_log))]
b2_filtered = b2_log[,!grepl("Unknown", colnames(b2_log))]

b1_filtered = b1_log
b2_filtered = b2_log

###########################
combinations = data.frame()
#can i rewrite this a little bit?
combinations = dataset_merger(0.01, 0.75, combinations)
combinations_cleaned = remove_duplicates(combinations)
data_annotated = data_formater(combinations_cleaned, b1_filtered, b2_filtered)

#check in lcms batch
data_annotated = data_annotated %>% add_column(0, .after = "Class")
colnames(data_annotated)[2] = 'Batch'
data_annotated$Batch[grepl("B1E", rownames(data_annotated))] = 1
data_annotated$Batch[grepl("B2E", rownames(data_annotated))] = 1
data_annotated$Batch[grepl("B1G", rownames(data_annotated))] = 1
data_annotated$Batch[grepl("B2G", rownames(data_annotated))] = 1
data_annotated$Batch[grepl("B5G", rownames(data_annotated))] = 2
data_annotated$Batch[grepl("B6E", rownames(data_annotated))] = 2
data_annotated$Batch[grepl("B7G", rownames(data_annotated))] = 2

write.csv(data_annotated[,-2],'../data/metabolomics/intermediates/data_phase.csv')
write.csv(data_annotated[,-1],'../data/metabolomics/intermediates/data_batched.csv')






