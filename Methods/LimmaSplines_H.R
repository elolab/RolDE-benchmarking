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
  library(limma)
  library(splines)
  
  control=sample1
  case=sample2
  data.temp=data[,c(grep(paste(case), colnames(data)), grep(control, colnames(data)))]
  
  nr.timepoints=(ncol(data.temp)/2)/3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
  
  temp.temp=rep(sort(rep(paste(seq(1:nr.timepoints)),3)),2) 
  X <- ns(as.numeric(temp.temp), df=(nr.timepoints-1)) 
    
  Group <- factor(c(rep(case, nr.timepoints*3), rep(control,nr.timepoints*3)))
  design <- model.matrix(~Group*X)
  fit <- lmFit(data.temp, design)
  fit <- eBayes(fit)
  
  coefs=grep("Group", colnames(fit$coefficients))
  top=topTable(fit, coef=c(coefs), number = nrow(data.temp))
  res=data.frame(cbind(rownames(top), top$P.Value), stringsAsFactors = F)
  colnames(res)=c("ID", "pval")
  res[,2]=as.numeric(res[,2])
  res[,1]=as.character(res[,1])
  rownames(res)=res[,1]
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure
  
  return(res)
}
