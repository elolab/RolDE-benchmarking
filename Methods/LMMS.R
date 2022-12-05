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
  library(lmms)
  nr.timepoints=(ncol(data)/2)/3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
  
  des_matrix=matrix(nrow = ncol(data), ncol = 4)
  des_matrix[,1]=colnames(data)
  des_matrix[,2]=c(rep(sample1, nr.timepoints*3), rep(sample2, nr.timepoints*3))
  des_matrix[,3]=c(sort(rep(seq(1:nr.timepoints),3)), sort(rep(seq(1:nr.timepoints),3)))
  des_matrix[,4]=c(rep(seq(1:3),nr.timepoints),rep(seq(from=4, to=6, by=1),nr.timepoints)) 
  
  des_matrix=data.frame(des_matrix)
  colnames(des_matrix)=c("Sample Names", "Condition", "Timepoint", "Individual")
  
  des_matrix[,1]<-as.character(des_matrix[,1])
  des_matrix[,2]<-as.factor(as.character(des_matrix[,2]))
  des_matrix[,3]<-as.numeric(as.character(des_matrix[,3]))
  des_matrix[,4]<-as.numeric(as.character(des_matrix[,4]))
  
  #Run lmms
  data=t(data)
  
  time=as.numeric(des_matrix$Timepoint)
  sampleID=paste("Ind_", as.character(des_matrix$Individual), sep="")
  group=as.character(des_matrix$Condition)
  
  lmmsDEtest <-lmmsDE(data=data,time=time,
                      sampleID=sampleID,group=group, type = "all", experiment = "all", numCores = 3)
  res=lmmsDEtest@DE
  rr=cbind(res$Group, res$Group_Time)
  rownames(rr)=res$Molecule
  pvals=apply(rr, 1, min, na.rm=T)
  pvals[is.infinite(pvals)]=NA

  res=data.frame(res$Molecule, pvals, stringsAsFactors = F)
  res[,1]=as.character(res[,1])
  res[,2]=as.numeric(res[,2])
  names(res)=c("ID", "pvals")
  rownames(res)=res[,1]
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure
  
  return(res)
}