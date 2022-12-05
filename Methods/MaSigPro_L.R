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
  library(maSigPro)
  control=sample1
  case=sample2
  data.temp=data[,c(grep(case, colnames(data)), grep(control, colnames(data)))]
  nr.timepoints=(ncol(data.temp)/2)/3 ##Applicable for datasets with 2 conditions and 3 replicates in each condition
  thres=floor(nr.timepoints/2)
  thres=(thres+1)*2+1
  
  #maSigPro user manual: genes with less than this number of true numerical values 
  #will be excluded from the analysis. Minimum value to estimate the model 
  #is (degree+1)xGroups+1
  #this is used here.
  
  #define the design matrix
  des=matrix(nrow=ncol(data.temp), ncol=3+length(case))
  rownames(des)=colnames(data.temp)
  colnames(des)=c("Time", "Replicate", control, case)
  des[,1]=rep(rep(c(seq(1:nr.timepoints)),each=3),(1+length(case))) #assumes 3 replicates
  des[,2]=sort(rep(seq(1:(nrow(des)/3)),3)) #Replicate measurements for timepoint
  des[,3]=0
  des[grep(control, rownames(des)),3]=1
  
  des[,4]=0
  des[grep(case, rownames(des)),4]=1
  
  degree=floor(nr.timepoints/2)
  design <- make.design.matrix(des, degree = degree) 
  fit <- p.vector(data.temp, design, Q = 1, MT.adjust = "BH", min.obs = thres)
  
  tstep <- T.fit(fit, step.method = "backward", alfa = 0.05, min.obs = thres)
  p=tstep$sol
  locs=grep(case, colnames(p)) #We only want p-values for those coefficients related to differences in longitudinal patterns between the conditions / groups. Not all; not those related to time in general over both conditions / groups. Or the p-value for the overall model fit.
  p=p[,locs]
  min.p=apply(p, 1, function(x) min(as.numeric(x), na.rm = T))
  min.p[which(is.infinite(min.p))]=NA
  
  res=data.frame(cbind(rownames(p), min.p), stringsAsFactors = F)
  res[,2]=as.numeric(res[,2])
  res[,1]=as.character(res[,1])
  colnames(res)=c("ID", "pval")
  rownames(res)=res[,1]
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure
  
  return(res)
}

