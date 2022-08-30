
library(hash)
library(data.table)
library(tidyverse)
library(pheatmap)
library(viridis)
library(rstudioapi)

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

source("supp/supp.R")

strains = c('WildType','YEL071W','YGR044C','YGR067C','YGR161C','YMR291W',
            'YNL099C','YNL289W','YOL051W','YOR317W','YOR351C')

#load data
ftq_e = read.csv('../data/metabolomics/ft_ethanol.csv', row.names = 1)
ftq_g = read.csv('../data/metabolomics/ft_glucose.csv', row.names = 1)

to_remove = c('243.05', '132.10147','156.07756','160.0429',
              '160.04425','160.04459','136.06218','136.06087',
              '296.06558','376.13256','396.1585','250.17918')

#filter out unwanted duplicates
ftq_e = remove_unwanted_features_ftq(ftq_e, to_remove)
ftq_g = remove_unwanted_features_ftq(ftq_g, to_remove)

st_e = read.csv('../data/metabolomics/st_ethanol.csv')
st_g = read.csv('../data/metabolomics/st_glucose.csv')

stylized_names = c('WT','YGR', 'TDA1', 'RME1','FAA1','GAL11','RTS3','MEK1','DLD3','OCA1','PCL1')
stylized_symbolnames_e = c('WT_E','\U0394ygr067c_E', '\U0394tda1_E', '\U0394rme1_E',paste0('\U0394', 'faa1_E'),
                         '\U0394gal11_E','\U0394rts3_E','\U0394mek1_E',paste0('\U0394','dld3_E'),'\U0394oca1_E','\U0394pcl1_E')
stylized_symbolnames_g = c('WT_G','\U0394ygr067c_G', '\U0394tda1_G', '\U0394rme1_G',paste0('\U0394', 'faa1_G'),
                           '\U0394gal11_G','\U0394rts3_G','\U0394mek1_G',paste0('\U0394','dld3_G'),'\U0394oca1_G','\U0394pcl1_G')

#average replicates
df = data.frame()
for (strain_x in stylized_names){
  means = colMeans(ftq_e[which(st_e$X %like% paste0(strain_x,'|',tolower(strain_x))),])
  df = rbind(df,c(strain_x, 'Post-shift',means))
}
df$X.WT. = stylized_symbolnames_e
colnames(df)[c(1,2)] = c('Strain','Phase')
colnames(df)[seq(3,433)] = colnames(ftq_e)

for (strain_x in stylized_names){
  means = colMeans(ftq_g[which(st_g$X %like% paste0(strain_x,'|',tolower(strain_x))),])
  df = rbind(df,c(strain_x, 'Pre-shift',means))
}
df$Strain[seq(12,22)] = stylized_symbolnames_g

#load limma results
signData = read.csv('../results/phase/Phase_noIQR.csv')
signData = signData[signData$P.Value < 0.05,]
#load pls-da results
plsData = read.csv('../results/metaboanalyst/oplsda_vip.csv')
plsData = plsData[order(plsData$V1, decreasing = T),]
#filter out top metabolites
filtData = plsData[match(plsData$X, signData$X),]
filtData = filtData[filtData$V1 > 1.131,] #top 100
#merge with data
filtData = drop_na(filtData)
df_heatmap = df
matching_indices = data.frame(match(colnames(df_heatmap),filtData$X))
df_heatmap = df_heatmap[which(!(is.na(matching_indices)))]
#convert to numeric
df_heatmap <- data.frame(apply(df_heatmap, 2, function(x) as.numeric(as.character(x))))
rownames(df_heatmap) = df$Strain

#plot heatmap
anno = data.frame(df$Phase)
colnames(anno) = 'Phase'
rownames(anno) = df$Strain
anno$Strain = c('WT','\U0394ygr067c', '\U0394tda1', '\U0394rme1',paste0('\U0394', 'faa1'),'\U0394gal11',
                '\U0394rts3','\U0394mek1',paste0('\U0394','dld3'),'\U0394oca1','\U0394pcl1')
annot_colors=list(Phase=c('Pre-shift'="#FF9D93",'Post-shift'="#00CEFF"))

pheatmap(t(df_heatmap),show_rownames = F,cluster_columns = T,scale = 'row',annotation_col = anno,
         clustering_distance_cols = "correlation", clustering_distrance_rows = 'correlation',fontsize = 13, 
         color = magma(100), border_color = '#181a1a', cutree_cols = 2,annotation_colors=annot_colors)





