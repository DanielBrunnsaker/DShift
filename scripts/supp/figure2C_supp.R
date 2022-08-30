
remove_unwanted_features = function(mdat, toremove){
  #remove duplicates
  for (mz in toremove){
    mdat = mdat[!mdat$X %like% mz,]
  }
  return(mdat)
}

remove_unwanted_features_ftq = function(mdat, toremove){
  #remove duplicates
  for (mz in toremove){
    mdat = mdat[!colnames(mdat) %like% mz]
  }
  return(mdat)
}

add_inchi_kegg = function(mdat, annomdat, kegg_inchi){
  
  #lazy way to add KEGG-annotations due to issues with float rounding
  
  mnames = str_split(mdat$X, '__')[[12]]
  count = 1
  for (mnames in str_split(mdat$X, '__')){
    src_c = 1
    for (src_mz in annomdat$mz){
      if (between(as.numeric(src_mz),as.numeric(mnames[2])-0.0001,as.numeric(mnames[2])+0.0001)){
        idx = src_c
      }
      src_c = src_c+1
    }
    mdat$inchikey[count] = annomdat[idx,3]
    kegg_id = kegg_inchi$KEGG[which(kegg_inchi$InChIKey %like% annomdat[idx,3])]
    if (length(kegg_id) != 0){
      mdat$KEGG[count] = kegg_id
    } else {
      mdat$KEGG[count] = NA
    }
    count = count+1
  }
  return(mdat)
}

p_threshold = function(mdat, threshold){
  tempdat = mdat[!is.na(mdat$KEGG),]
  tempdat = tempdat[!(tempdat$KEGG == 'No result'),]
  tempdat = tempdat[tempdat$P.Value < threshold,]
  return(tempdat)
}


correlation_values = function(indices, models, strains_sims, phase, corr_values, corr_strains){
  
  data(iMM904)
  #calculates spearman correlation for all combinations of simulations per phase
  m1_sims = data.frame(iMM904@react_id)
  ms_sims = data.frame(iMM904@react_id)
  
  for (strains in strains_sims){
    for (model in models){
      print(strains)
      #load simulations
      sim = read.csv(paste0('/Users/danbru/',model,'_',
                            strains,'.csv'),row.names = 1, header= TRUE) #first half is glucose, second ethanol
      sim = sim[,indices] #glucose indices
      if (model == 'M1'){
        sim_m1 = sim
      } else {
        df_ps = data.frame(rownames(m1_sims))
        df_ps$p_vals = 1
        #first check how many reaction fluxes are sign different
        corrlist = c()
        for (i in 1:1577){
          wilc = wilcox.test(as.numeric(sim[i,]), as.numeric(sim_m1[i,]), paired=TRUE) #check for sign. differing reaction fluxes
          df_ps$p_vals[i] = wilc$p.value
          
        }
        for (k in 1:13){
          for (kk in (k+1):13){}
          corrps = cor.test(as.numeric(sim[,k]), as.numeric(sim_m1[,kk]), method = "spearman")
          corrlist = c(corrlist, corrps$estimate)
          corr_strains = rbind(corr_strains, c(strains,phase, corrps$estimate))
        }
        corrlist[is.na(corrlist)] <- 1
        print(mean(corrlist,na.rm = T))
        print(length(which(df_ps$p_vals < 0.05)))
        corr_values = rbind(corr_values, c(strains,phase, mean(corrlist,na.rm = T), length(which(df_ps$p_vals < 0.05)))) #count occurences of sign. changed reactions
        
      }
    }
  }
  return(list(corr_values, corr_strains))
}



library(hash)

#define strain-stylized name hashmap
map = hash()
map[["WildType"]]= 'WT'
map[['YGR067C']] = paste0('ygr|YGR')
map[['YMR291W']] = paste0('tda1|TDA1')
map[['YGR044C']] = paste0('rme1|RME1')
map[['YOR317W']] = paste0('faa1|FAA1')
map[['YOL051W']] = paste0('gal11|GAL11')
map[['YGR161C']] = paste0('rts3|RTS3|rst3')
map[['YOR351C']] = paste0('mek1|MEK1')
map[['YEL071W']] = paste0('dld3|DLD3')
map[['YNL099C']] = paste0('oca1|OCA1')
map[['YNL289W']] = paste0('pcl1|PCL1')
