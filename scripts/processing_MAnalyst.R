
data_phase = read.csv('../data/metabolomics/intermediates/data_phase.csv')
data_batch = read.csv('../data/metabolomics/intermediates/data_batched.csv')

data_phase = add_column(data_phase,'WT',.after=0)
colnames(data_phase)[1] = 'Strain'
data_phase$Strain = 'WT'
data_phase = add_column(data_phase,'Batch',.after=1)
colnames(data_phase)[2] = 'Batch'
data_phase$Batch[data_batch$Batch == 1] = 1
data_phase$Batch[data_batch$Batch == 2] = 2

data_phase$Strain[data_phase$X %like% 'tda1|TDA1'] = 'TDA1'
data_phase$Strain[data_phase$X %like% 'rme1|RME1'] = 'RME1'
data_phase$Strain[data_phase$X %like% 'pcl1|PCL1'] = 'PCL1'
data_phase$Strain[data_phase$X %like% 'rts3|rst3|RTS3'] = 'RTS3'
data_phase$Strain[data_phase$X %like% 'dld3|DLD3'] = 'DLD3'
data_phase$Strain[data_phase$X %like% 'faa1|FAA1'] = 'FAA1'
data_phase$Strain[data_phase$X %like% 'oca1|OCA1'] = 'OCA1'
data_phase$Strain[data_phase$X %like% 'mek1|MEK1'] = 'MEK1'
data_phase$Strain[data_phase$X %like% 'ygr|YGR'] = 'YGR067C'
data_phase$Strain[data_phase$X %like% 'gal11|GAL11'] = 'GAL11'
data_phase$Strain[data_phase$X %like% 'tda1|TDA1'] = 'TDA1'
rownames(data_phase) = data_phase[,3]
data_phase = data_phase[,-3]

st = data_phase[,c(1,2,3)]
ft = data_phase[,-c(1,2,3)]

rownames(st) = rownames(ft)
write.csv(st,'../data/metabolomics/intermediates/st.csv')
write.csv(ft,'../data/metabolomics/intermediates/ft.csv')

### strainwise prep

library(preprocessCore)
ft = data.frame(t(normalize.quantiles(t(ft)))) #preprocess now to save time later
rownames(ft) = rownames(data_phase)
colnames(ft) = colnames(data_phase[,-c(1,2,3)])

#do for strains
ft_ethanol = ft[grepl("Ethanol", st$Class),]
st_ethanol = st[grepl("Ethanol", st$Class),]

for (strainname in st_ethanol$Strain){
  ft_strain = ft_ethanol[grepl(paste0(strainname,"|WT"), st_ethanol$Strain),]
  st_strain = st_ethanol[grepl(paste0(strainname,"|WT"), st_ethanol$Strain),]
  write.csv(ft_strain,paste0('../data/metabolomics/intermediates/ethanol/ft_',strainname,'.csv'))
  write.csv(st_strain,paste0('../data/metabolomics/intermediates/ethanol/st_',strainname,'.csv'))
}

#do for strains
ft_glucose = ft[grepl("Glucose", st$Class),]
st_glucose = st[grepl("Glucose", st$Class),]

for (strainname in st_glucose$Strain){
  ft_strain = ft_glucose[grepl(paste0(strainname,"|WT"), st_glucose$Strain),]
  st_strain = st_glucose[grepl(paste0(strainname,"|WT"), st_glucose$Strain),]
  write.csv(ft_strain,paste0('../data/metabolomics/intermediates/glucose/ft_',strainname,'.csv'))
  write.csv(st_strain,paste0('../data/metabolomics/intermediates/glucose/st_',strainname,'.csv'))
}
