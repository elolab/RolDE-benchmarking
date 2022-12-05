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
  library(OmicsLonDA)
  library(SummarizedExperiment)
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
  
  des_mat=data.frame(des_matrix$Individual, des_matrix$Condition, des_matrix$Timepoint)
  rownames(des_mat)=des_matrix$`Sample Names`
  colnames(des_mat)=c("Subject", "Group", "Time")
  
  omicslonda_se_object = SummarizedExperiment(assays=list(data),
                                              colData = des_mat)
  omicslonda_se_object_adjusted = adjustBaseline(se_object = omicslonda_se_object)
  points = seq(1, nr.timepoints, length.out = nr.timepoints)
  
  stats=numeric(nrow(data))
  
  for(i in 1:nrow(data)){
    omicslonda_test_object = omicslonda_se_object_adjusted[i,]
    
    res=tryCatch({
      omicslonda(se_object = omicslonda_test_object, n.perm = 10,
                 fit.method = "ssgaussian", points = points, text = "Feature_x",
                 parall = FALSE, pvalue.threshold = 0.05, 
                 adjust.method = "BH", time.unit = "points",
                 ylabel = "Normalized Count",
                 col = c("blue", "firebrick"), prefix = tempfile())
    },error = function(e) {
      NULL
    })
    
    if(!is.null(res)){
      feature.summary = as.data.frame(do.call(cbind, res$details),
                                      stringsAsFactors = FALSE)
      if(is.data.frame(feature.summary)){
        if(nrow(feature.summary)>0){
          stats[i]=max(as.numeric(feature.summary$testStat.abs, na.rm = T))
        }else{
          stats[i]=NA
        }
      }else{
        stats[i]=NA
      }
      
    }else{
      stats[i]=NA
    }
  }
  
  stats=as.numeric(stats)
  stats=stats*(-1) #Change for ROC curves
  
  res=data.frame(rownames(data), stats, stringsAsFactors = F)
  colnames(res)=c("ID", "stat")
  res[,1]=as.character(res[,1])
  res[,2]=as.numeric(res[,2])
  rownames(res)=res[,1]
  if(any(is.nan(res[,2]))){res[which(is.nan(res[,2])),2]=NA} #make sure
  
  return(res)
}