
library(readxl)
library(ggplot2)
library(reshape2)
library(patchwork)
library(data.table)
library(tidyverse)
library(rstudioapi)

#set working dir
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

source("supp/supp.R")

#ethanol phase enrichment
pathways = read_excel('../results/strain/FELLA_pathways.xlsx', sheet = 'inverse_e')
long <- melt(pathways, id.vars = c("Pathway"))
long$value[which(long$value < (1-0.01))] = 0
x_axis = c('\U0394tda1','\U0394ygr067c','\U0394rme1',paste0('\U0394','dld3'),
           '\U0394gal11', '\U0394oca1', paste0('\U0394', 'faa1'),'\U0394rts3','\U0394pcl1','\U0394mek1')
p1 = ggplot(long, aes(x=variable, y=Pathway, size = qnorm(value), 
                      color = qnorm(value)))+geom_point()+theme_linedraw()+theme(axis.text.x = element_text(angle=90),
                                                                                text = element_text(size=15))+scale_x_discrete(labels= x_axis)+scale_y_discrete(position = "right")+ggtitle('Post-shift')

#glucose phase enrichment
pathways = read_excel('../results/phase/FELLA_pathways.xlsx', sheet = 'inverse_g')
long <- melt(pathways, id.vars = c("Pathway"))
long$value[which(long$value < (1-0.01))] = 0
p2 = ggplot(long, aes(x=variable, y=Pathway, size = qnorm(value), 
                      color = qnorm(value)))+geom_point()+theme_linedraw()+theme(axis.text.x = element_text(angle=90),
                                                                                text = element_text(size=15))+scale_x_discrete(labels= x_axis)+scale_y_discrete(position = "right")+ggtitle('Pre-shift')

#plot in the correct layout
patchwork = p2/p1
patchwork[[1]] = patchwork[[1]] + theme(axis.title.x = element_blank())
patchwork[[2]] = patchwork[[2]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
patchwork+plot_layout(guides = "collect")

