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
  library(betr)
  library(Biobase)
 
  control=sample1
  case=sample2
  data.temp=data[,c(grep(control, colnames(data)), grep(case, colnames(data)))]
  nr.timepoints=(ncol(data.temp)/2)/3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
  
  data.temp=na.omit(data.temp) #BETR does not tolerate missing values
  temp.set <- ExpressionSet(assayData=as.matrix(data.temp))
  	
  condition <- rep(c(control, case),each=nr.timepoints*3)
  timep <- rep(rep(seq(1:nr.timepoints),each=3),2)
  reps <-  c(rep(c(1,2,3),nr.timepoints), rep(c(4,5,6),nr.timepoints))  #assumes three non-paired replicates
  prob=betr(eset=temp.set, cond=condition,
             timepoint=timep,twoCondition = T, replicate=reps, alpha=0.05)
  res=data.frame(names(prob), prob, stringsAsFactors = F)
  colnames(res)=c("ID", "deprob")
  res[,2]=as.numeric(res[,2])*(-1) #switch for ROCs
  res[,1]=as.character(res[,1])
  rownames(res)=res[,1]
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure
  
  return(res)
}


