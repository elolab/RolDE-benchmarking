##Assumes that the samples in the given data are organized as:
#Condition 1 timepoint 1 replicate 1,
#Condition 1 timepoint 1 replicate 2,
#Condition 1 timepoint 1 replicate 3,
#Condition 1 timepoint 2 replicate 1,
#Condition 1 timepoint 2 replicate 2,
#.
#.
#.
#Condition 2 timepoint 1 replicate 1,
#Condition 2 timepoint 1 replicate 2,
#Condition 2 timepoint 1 replicate 3,
#Condition 2 timepoint 2 replicate 1,
#Condition 2 timepoint 2 replicate 2,

getDE<-function(data, sample1, sample2){
  library(timecourse)
  
  control=sample1
  case=sample2
  
  nas=apply(data, 1, function(x) length(which(is.na(x))))
  rems=which(nas==ncol(data))
  
  if(length(rems)>0){
    data.filt=data[-rems,]
  } else {
    data.filt=data
  }
 
  size <- matrix(3, nrow=nrow(data.filt), ncol=2)  #replicate, 2 conditions
  
  #As mentioned above, assumes 3 replicates/condition and the data ordered in a certain way
  data1_r1=data.filt[,seq(1,ncol(data.filt),3)]
  data1_r2=data.filt[,seq(2,ncol(data.filt),3)]
  data1_r3=data.filt[,seq(3,ncol(data.filt),3)]
  
  nr.timepoints=(ncol(data)/2)/3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
  start1=1
  end1=nr.timepoints
  start2=end1+1
  end2=nr.timepoints*2
  
  # missing time point samples are not allowed while missing replicates are allowed for a subset of genes. In the latter case, the user still need
  #to input a complete data matrix with missing replicates indicated by NAs.
  for(i in 1:nrow(data1_r1)){
    x1=as.numeric(data1_r1[i,start1:end1])
    if(any(is.na(x1))){
      x1[which(!is.na(x1))]=NA
      size[i,1]=size[i,1]-1
    }
    
    x2=as.numeric(data1_r1[i,start2:end2])
    if(any(is.na(x2))){
      x2[which(!is.na(x2))]=NA
      size[i,2]=size[i,2]-1
    }
    
    data1_r1[i,]=c(x1,x2)
  }
  
  for(i in 1:nrow(data1_r2)){
    x1=as.numeric(data1_r2[i,start1:end1])
    if(any(is.na(x1))){
      x1[which(!is.na(x1))]=NA
      size[i,1]=size[i,1]-1
    }
    
    x2=as.numeric(data1_r2[i,start2:end2])
    if(any(is.na(x2))){
      x2[which(!is.na(x2))]=NA
      size[i,2]=size[i,2]-1
    }
    
    data1_r2[i,]=c(x1,x2)
  }
  
  for(i in 1:nrow(data1_r3)){
    x1=as.numeric(data1_r3[i,start1:end1])
    if(any(is.na(x1))){
      x1[which(!is.na(x1))]=NA
      size[i,1]=size[i,1]-1
    }
    
    x2=as.numeric(data1_r3[i,start2:end2])
    if(any(is.na(x2))){
      x2[which(!is.na(x2))]=NA
      size[i,2]=size[i,2]-1
    }
    
    data1_r3[i,]=c(x1,x2)
  }
  data2=matrix(ncol = ncol(data.filt), nrow = nrow(data.filt))
  data2[,seq(1,ncol(data),3)]=data.matrix(data1_r1)
  data2[,seq(2,ncol(data),3)]=data.matrix(data1_r2)
  data2[,seq(3,ncol(data),3)]=data.matrix(data1_r3)
  colnames(data2)=colnames(data.filt)
  rownames(data2)=rownames(data.filt)
  
  nas=apply(data2, 1, function(x) length(which(is.na(x))))
  rems=which(nas==ncol(data2))
  
  if(length(rems)>0){
    data.filt=data2[-rems,]
    size=size[-rems,]
  } else {
    data.filt=data2
    size=size
  }
  
  rems1=which(size[,1]==0)
  rems2=which(size[,2]==0)
  rems=unique(c(rems1, rems2))
  
  if(length(rems)>0){
    data.filt=data.filt[-rems,]
    size=size[-rems,]
  } else {
    data.filt=data.filt
    size=size
  }
  
  data.temp=data.filt[,c(grep(paste(control), colnames(data.filt)), grep(case, colnames(data.filt)))]
  
  #reorder the data into the standard format of timecourse just in case
  
  #assumes three replicates, two conditions
  
  cntrl1=seq(1,(ncol(data)/2),3)
  cntrl2=seq(2,(ncol(data)/2),3)
  cntrl3=seq(3,(ncol(data)/2),3)
    
  case1=seq((ncol(data)/2+1),ncol(data),3)
  case2=seq((ncol(data)/2+2),ncol(data),3)
  case3=seq((ncol(data)/2+3),ncol(data),3)
  
  data.temp=data.temp[,c(cntrl1,cntrl2,cntrl3,case1,case2,case3)]
  
  trt <- rep(c(sample1, sample2),each=ncol(data.temp)/2)
  assay <- rep(c("rep1", "rep2", "rep3", "rep4", "rep5", "rep6"), each=nr.timepoints) #assumes two conditions, 3 replicates in each, non-paired, doesn't matter though if you define 3 or 6 replicates here, same result

  MB.test <- mb.long(data.temp, method="2D", times=nr.timepoints, reps=size, condition.grp=trt, rep.grp=assay)
  
  prot_positions = MB.test$pos.HotellingT2[1:length(MB.test$pos.HotellingT2)]
  protnames = rownames(data.temp)
  DE_prot = protnames[prot_positions]
  DE_test_statistic=MB.test$HotellingT2[prot_positions]
    
  res=data.frame(DE_prot, DE_test_statistic, stringsAsFactors = F)
  res[,2]=res[,2]*(-1) #switch for ROC
  colnames(res)=c("ID", "hotelling")
  res[,1]=as.character(res[,1])
  rownames(res)=res[,1]
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure
  
  return(res)
}
