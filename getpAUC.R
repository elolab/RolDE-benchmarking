library(pROC)
getAUC<-function(results, reference.prot, spike.prot, method){
  
  #Reference proteins are all the proteins in the data.
  
  results.ref=matrix(ncol = 2, nrow = length(reference.prot))
  results.ref[,1]=reference.prot
  locs=match(results.ref[,1], results[,1])
  seq.locs=seq(1:length(locs))
  
  if(any(is.na(locs))){
    results.ref[seq.locs[-which(is.na(locs))],2]=as.numeric(as.character(results[locs[-which(is.na(locs))],2]))
  } else {
    results.ref[seq.locs,2]=as.numeric(as.character(results[locs,2]))
  }
  
  results=results.ref
  results=results[sample(seq(1:nrow(results))),] #Order should not matter, but randomize result lists just in case
  classes=c(grepl(spike.prot,as.character(results[,1]))) #UPS1 spike-proteins required in a certain format
  
  values=c(as.numeric(as.character(results[,2])))
  locs=length(which(is.na(values)))
  
  #For a fair comparison of the methods, if we have NA values in the result list, we want to replace them with some large value so that features with NA result statistics get put at the bottom of the result list.
  
  if(length(locs)>0){
    if(method=="BETR" | method=="Timecourse" | method=="Limma" | method=="OmicsLonDA"){ #Our statistic is not p-values
      max.val=max(na.omit(values))
      imp.vals=sample(x = seq(from=max.val, to=0, by=0.0000001), size = locs, replace=T) #just some artificial large values placing the features randomly at the bottom of the result list, since we have no DE information regarding them.
      values[is.na(values)] <- imp.vals
    }else if(method=="RolDE"){ #just in case, RolDE should not produce any missing values
      max.val=max(na.omit(values))
      imp.vals=sample(x = seq(from=max.val, to=max.val+(locs*1000), by=1), size = locs, replace=T) #just some artificial large values placing the features randomly at the bottom of the result list, since we have no DE information regarding them.
      values[is.na(values)] <- imp.vals
    }else{ #we have p-values as our statistic
      max.val=max(na.omit(values))
      imp.vals=sample(x = seq(from=max.val, to=1, by=0.0000001), size = locs, replace=T) #just some artificial large p-value placing the features randomly at the bottom of the result list, since we have no DE information regarding them.
      values[is.na(values)]=imp.vals
    }
  }
  roc=roc(classes, values, direction = ">")
  auc.val=round(auc(roc,partial.auc=c(1.0,0.9), partial.auc.correct=T,allow.invalid.partial.auc.correct=T)[1],digits=5)# #set for partial auc
  
  return(auc.val)
}