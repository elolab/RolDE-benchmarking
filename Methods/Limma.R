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
  library(Biobase)
  library(limma)

  control=sample1
  case=sample2
  data.temp=data[,c(grep(case, colnames(data)), grep(control, colnames(data)))]
  nr.timepoints=(ncol(data.temp)/2)/3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
  
  case.group=sort(rep(paste(case,"_","timepoint_",seq(1:nr.timepoints), sep = ""),3))
  control.group=sort(rep(paste(control,"_","timepoint_",seq(1:nr.timepoints), sep = ""),3))    
  
  f<-factor(c(case.group, control.group))
  design <- model.matrix(~0+f)
  colnames(design) <- levels(f)
  fit <- lmFit(data.temp, design)
  
  mc=makeContrasts(contrasts=paste(unique(case.group),"-",unique(control.group), sep = ""), levels = design)
  c.fit = contrasts.fit(fit, mc)
  eb = eBayes(c.fit)
  res=data.frame(cbind(rownames(eb$p.value), eb$F), stringsAsFactors = F)
  
  colnames(res)[1]="ID"
  colnames(res)[2]="fval"
  res[,2]=as.numeric(res[,2])*(-1) #Switch for ROCs
  res[,1]=as.character(res[,1])
  rownames(res)=res[,1]
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure
  
  return(res)
}
