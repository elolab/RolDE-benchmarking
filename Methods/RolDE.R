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
 library(RolDE) #v.099.4
  
  nr.timepoints=(ncol(data)/2)/3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
  
  des_matrix=matrix(nrow = ncol(data), ncol = 4)
  des_matrix[,1]=colnames(data)
  des_matrix[,2]=c(rep(sample1, nr.timepoints*3), rep(sample2, nr.timepoints*3))
  des_matrix[,3]=c(sort(rep(seq(1:nr.timepoints),3)), sort(rep(seq(1:nr.timepoints),3)))
  des_matrix[,4]=c(rep(seq(1:3),nr.timepoints),rep(seq(from=4, to=6, by=1),nr.timepoints)) 

  set.seed(1)
  res=RolDE(data = data, des_matrix = des_matrix, aligned = TRUE, n_cores = 1, sigValSampN = 0)
  res=res[[1]]
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure

  return(res)
}
