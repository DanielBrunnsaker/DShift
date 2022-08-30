
library(BSDA)
library(CoRegFlux)
library(dplyr)
library(ggplot2)
library(rstudioapi)

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

source("supp/supp.R")

#load the csvs
strains_sims = c('WildType', 'YGR067C','YMR291W','YGR044C','YOR317W','YOL051W',
                 'YGR161C','YOR351C','YEL071W','YNL099C','YNL289W')

glucose_indices = c(1,2,3,4,5,6,7,8,9,10,11,12,13)
ethanol_indices = c(14,15,16,17,18,19,20,21,22,23,24,25,26)
models = c('M1','M1Smart')

corr_values = data.frame()
corr_strains = data.frame()
corr_lists_g = correlation_values(glucose_indices, models, strains_sims, 'Glucose', corr_values, corr_strains)
corr_values = corr_lists_g[[1]]
corr_strains = corr_lists_g[[2]]

corr_lists_e = correlation_values(ethanol_indices, models, strains_sims, 'Ethanol', corr_values, corr_strains)
corr_values = corr_lists_e[[1]]
corr_strains = corr_lists_e[[2]]

#create cute violinplot with correlations, Figure 4A
colnames(corr_strains) = c('Strain','Phase','Rho')
corr_strains$Rho = as.numeric(corr_strains$Rho)
stylized_symbolnames = c('WT','\U0394ygr067C', '\U0394tda1', '\U0394rme1',paste0('\U0394', 'faa1'),'\U0394gal11','\U0394rts3','\U0394mek1',paste0('\U0394','dld3'),'\U0394oca1','\U0394pcl1')

p = ggplot(corr_strains, aes(x=Strain, y=Rho, fill=Phase)) +geom_boxplot(position=position_dodge(0.9), width = 0.7)+theme_bw()+theme(text = element_text(size=25),
                                                                                                                                     axis.text.x = element_text(angle=60, hjust=1))
p+scale_x_discrete(labels=c("YGR067C" = '\U0394ygr067c', "YEL071W" = paste0('\U0394','dld3'),"YGR044C" = '\U0394rme1',
                            "YGR161C" = '\U0394rts3',"YMR291W" = '\U0394tda1',"YNL099C" = '\U0394oca1',
                            "YNL289W" = '\U0394pcl1',"YOL051W" = '\U0394gal11',"YOR317W" = paste0('\U0394', 'faa1'), 
                            "YOR351C" = '\U0394mek1', "WildType" = 'WT'))+ylab('Spearman rho')+ylim(c(0.52,0.72))


#Figure 4B
colnames(corr_values) = c('Strain','Phase','Rho','Reactions')
corr_values$Reactions = as.numeric(corr_values$Reactions)
#rescale to 0-1
corr_values$Reactions = (corr_values$Reactions-0)/(1577-0)
#then heatmap of changing reactions
phmap = ggplot(corr_values, aes(Strain, Phase, fill=Reactions)) + 
  geom_tile(width=0.95, height=0.95,color='black')+scale_fill_viridis_c(option = "magma",begin = 0.05,end = 0.65)+theme_bw()+theme(text = element_text(size=25),axis.text.x = element_text(angle=45, hjust=1))

phmap+scale_x_discrete(expand=c(0,0),labels=c("YGR067C" = '\U0394ygr067c', "YEL071W" = paste0('\U0394','dld3'),"YGR044C" = '\U0394rme1',
                                              "YGR161C" = '\U0394rts3',"YMR291W" = '\U0394tda1',"YNL099C" = '\U0394oca1',
                                              "YNL289W" = '\U0394pcl1',"YOL051W" = '\U0394gal11',"YOR317W" = paste0('\U0394', 'faa1'), "YOR351C" = '\U0394mek1', "WildType" = 'WT'))+
  scale_y_discrete(expand=c(0,0),labels = c('Post-shift','Pre-shift'))
