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
  library(nlme)
  
  control=sample1
  case=sample2
  
  #As mentioned above, assumes 3 replicates/condition and the data ordered in a certain way
  data1_r1=data[,seq(1,ncol(data),3)]
  data1_r2=data[,seq(2,ncol(data),3)]
  data1_r3=data[,seq(3,ncol(data),3)]
  
  nr.timepoints=(ncol(data)/2)/3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
  
  res.frame=matrix(nrow=nrow(data), ncol=2)
  rownames(res.frame)=rownames(data)

  for(l in 1:nrow(data)) {
    temp.data=matrix(nrow = nr.timepoints*6, ncol = 4)
    colnames(temp.data)=c("Intensity", "Time", "Condition", "Replicate")
    temp.data[,1]=c(as.numeric(data1_r1[l,c(grep(control, colnames(data1_r1)),grep(case, colnames(data1_r1)))]),
                    as.numeric(data1_r2[l,c(grep(control, colnames(data1_r2)),grep(case, colnames(data1_r2)))]),
                    as.numeric(data1_r3[l,c(grep(control, colnames(data1_r3)),grep(case, colnames(data1_r3)))]))
    temp.data[,2]=c(rep(seq(1:nr.timepoints), 6))
    temp.data[,3]=c(rep(c(rep(1,nr.timepoints), rep(2,nr.timepoints)),3))
    temp.data[,4]=c(rep(1,nr.timepoints), rep(2,nr.timepoints), rep(3,nr.timepoints), rep(4,nr.timepoints), rep(5,nr.timepoints), rep(6,nr.timepoints))
   
    temp.data=data.frame(temp.data, stringsAsFactors = F)
    temp.data$Condition=factor(temp.data$Condition)
    temp.data$Replicate=factor(temp.data$Replicate)
    
    degr=nr.timepoints-1
    
    lme1 = tryCatch({
      lme(Intensity ~ poly(Time, degree = get("degr"))*Condition, random = ~ 1 | Replicate,
          data = temp.data, na.action = na.omit, method = "ML")
    }, error = function(e) {tryCatch({
      lme(Intensity ~ poly(Time, degree = get("degr"))*Condition, random = ~ 1 | Replicate,
          data = temp.data, na.action = na.omit, method = "ML", control = lmeControl(opt = 'optim'))
    }, error = function(e2) {
      NULL
      })
    })
      
    
    lme2 = tryCatch({
      lme(Intensity ~ poly(Time, degree = get("degr"))*Condition, random = ~ Time | Replicate,
          data = temp.data, na.action = na.omit, method = "ML")
    }, error = function(e) {
      NULL
    })
    
    base.mod = tryCatch({
      lme(Intensity ~ 1, random = ~ 1 | Replicate, data = temp.data, na.action = na.omit, method = "ML")
    }, error = function(e) {tryCatch({
      lme(Intensity ~ 1, random = ~ 1 | Replicate, data = temp.data, na.action = na.omit, method = "ML", control = lmeControl(opt = 'optim'))
    }, error = function(e2) {
      NULL
      })
    })
    
    if(!is.null(lme2) & !is.null(lme1)){
      mod.test=anova(lme1, lme2)
      if(mod.test$`p-value`[2]>0.05){
        lme.fin=lme1
      } else {
        lme.fin=lme2
      }
    } else {
      lme.fin=lme1
    }
    
    if(!is.null(lme.fin) & !is.null(base.mod)){  
      sum=summary(lme.fin)
      sum=sum$tTable
      locs=grep("Condition", rownames(sum))
      if(length(locs)==1){
        sum=t(data.frame(sum[locs,]))
      } else {
        sum=sum[locs,]
      }
      min.p=min(as.numeric(sum[,5]), na.rm = T)
      if(is.infinite(min.p)){min.p=NA}
      
      res.frame[l,1]=rownames(res.frame)[l]    
      res.frame[l,2]=min.p
    } else {
      res.frame[l,1]=rownames(res.frame)[l]    
    }
  }
  res.frame=data.frame(res.frame, stringsAsFactors = F)
  res.frame[,1]=as.character(res.frame[,1])
  res.frame[,2]=as.numeric(res.frame[,2])
  colnames(res.frame)=c("ID", "pval")
  rownames(res.frame)=res.frame[,1]
  if(any(is.nan(res.frame[,2]))){res.frame[which(is.nan(res.frame[,2])),2]=NA} #make sure
  
  return(res.frame)
}
    
    
    
    