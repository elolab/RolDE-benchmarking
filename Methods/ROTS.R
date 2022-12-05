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
  library(ROTS)
  
  control=sample1
  case=sample2
  
  nr.timepoints=(ncol(data)/2)/3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
  
  #loop through the timepoints, assumes that the data is ordered according to timepoint and replicates
  data1=data[,grep(control, colnames(data))]
  data2=data[,grep(case, colnames(data))]
  
  rots.list=list()
  p.val.df=data.frame(matrix(nrow = nrow(data1), ncol = nr.timepoints))
  rownames(p.val.df)=rownames(data1)
  
  #As mentioned above, assumes 3 replicates/condition and the data ordered in a certain way
  ind1=1
  ind2=3
  
  for(i in 1:nr.timepoints){
    dat1=data1[,ind1:ind2]
    dat2=data2[,ind1:ind2]
    
    temp.data=cbind(dat1, dat2)
    groups=c(1,1,1,0,0,0)
    
    #take out proteins with too many missing values for ROTS
    rems1=which(rowSums(is.na(dat1))>1)
    rems2=which(rowSums(is.na(dat2))>1)
    rems=c(rems1, rems2)
    rems=unique(rems)
    rem.names=rownames(temp.data)[rems]
    if(length(rems)>0){temp.data=temp.data[-rems,]}
    
    rots.out=ROTS(data = temp.data, groups = groups, B = 100, K = nrow(temp.data)/4, paired = F, seed = 1, progress = F, log = T)
    rots.list[[i]]=rots.out
    rots.df=data.frame(p=rots.out$pvalue)
    rownames(rots.df)=rownames(rots.out$data)
    
    if(length(rems)>0){
      rem.frame=data.frame(matrix(nrow = length(rem.names), ncol = 1))
      rownames(rem.frame)=rem.names
      colnames(rem.frame)="p"
      rots.df=rbind(rots.df, rem.frame)
    }
    p.val.df[,i]=rots.df[match(rownames(p.val.df), rownames(rots.df)),]
    
    ind1=ind1+3
    ind2=ind2+3
  }
  
  rems=which(apply(p.val.df, 1, function(x) length(which(is.na(x))))==nr.timepoints)
  rem.names=rownames(p.val.df)[rems]
  rems2=rep(NA, length(rems))
  names(rems2)=rem.names
  
  if(length(rems)>0){p.val.df=p.val.df[-rems,]}
  
  p.vals=apply(p.val.df, 1, function(x) min(x, na.rm = T)) #Summarize with minimum
  names(p.vals)=rownames(p.val.df)
  
  if(length(rems)>0){p.vals=c(p.vals, rems2)}
  res=cbind(names(p.vals), p.vals)
  res[is.nan(res)]=NA
  
  res=data.frame(res, stringsAsFactors = F)
  res[,1]=as.character(res[,1])
  res[,2]=as.numeric(res[,2])
  colnames(res)=c("ID", "pval")
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure
  
  return(res)
}