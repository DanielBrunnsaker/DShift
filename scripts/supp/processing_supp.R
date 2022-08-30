
dataset_merger = function(mz_threshold, RT_threshold, combinations){
  rownumber_b1 = 2
  for (metabolite in colnames(b1_filtered)[2:dim(b1_filtered)[2]]){
    
    temp = strsplit(metabolite,'__')
    
    RT = as.numeric(substring(temp[[1]][1], 2)) #extract the retention time
    mz = as.numeric(temp[[1]][2]) #extract the m/z
    adduct = temp[[1]][3] #extract the adduct
    name = temp[[1]][4] #extract the name
    print(mz)
    rownumber_b2 = 2
    for (metabolite_b2 in colnames(b2_filtered)[2:dim(b1_filtered)[2]]){
      
      temp_2 = strsplit(metabolite_b2,'__')
      
      RT_2 = as.numeric(substring(temp_2[[1]][1], 2)) #extract the retention time
      mz_2 = as.numeric(temp_2[[1]][2], 2) #extract the m/z
      adduct_2 = temp_2[[1]][3] #extract the adduct
      name_2 = temp_2[[1]][4] #extract the name
      
      if((between(RT_2, RT-RT_threshold, RT+RT_threshold))&(between(mz_2[[1]], mz-mz_threshold, mz+mz_threshold))&(grepl(adduct, adduct_2, fixed = T))){
        print(c(mz_2, adduct_2))
        RTmin = abs(RT_2-RT)
        mzmin = abs(mz_2-mz)
        combinations = rbind(combinations, c(rownumber_b1,rownumber_b2, mzmin, RTmin))
        print(paste('Combined',colnames(b1_filtered)[rownumber_b1],'-WITH-',colnames(b2_filtered)[rownumber_b2]))
      }
      rownumber_b2 = rownumber_b2+1
    }
    rownumber_b1 = rownumber_b1+1
  }
  colnames(combinations) = c('ordering_b1', 'ordering_b2', 'delta_mz', 'delta_RT')
  return(combinations)
}


remove_duplicates = function(combinations){
  #take out matches, order by indice from b1 first, then mz_diff, so as to remove duplicates, keeping the one with the closest m/z match
  mzsorted = combinations[order(combinations$ordering_b1, combinations$delta_mz),]
  mzsorted = mzsorted[!duplicated(mzsorted$ordering_b1),]
  #take out matches, order by indice from b2 first, then mz_diff, so as to remove duplicates, keeping the one with the closest m/z match
  mzsorted = mzsorted[order(mzsorted$ordering_b2, mzsorted$delta_mz),]
  mzsorted = mzsorted[!duplicated(mzsorted$ordering_b2),]
  return(mzsorted)
}

data_formater = function(data, data_b1, data_b2){
  #insert the classrow
  combinations = data[,c(1,2)]
  combinations <- rbind(c(1,1), combinations)
  
  b1_matched = data_b1[,combinations$ordering_b1]
  b2_matched = data_b2[,combinations$ordering_b2]
  
  namevec = colnames(b1_matched)
  colnames(b1_matched) = seq(1,length(namevec))
  colnames(b2_matched) = seq(1,length(namevec))
  
  data = rbind(b1_matched, b2_matched)
  colnames(data) = namevec
  return(data)
}