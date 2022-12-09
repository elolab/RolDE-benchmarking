##Clear workspace
rm(list=ls())

##Make sure your working directory is set to the where the contents of the Git Repo were downloaded / cloned.

##The following R packages are required
package_list<-c("ROTS", "nlme", "Biobase", "limma", "reshape2", "splines",
                "timecourse", "maSigPro", "edge", "OmicsLonDA", "SummarizedExperiment",
                "RolDE", "pROC", "mvtnorm", "lmeSplines", "gdata", "gplots")

##Tries to install the required packages. Comment out if unnecessary.
setRepositories(ind = c(1,2,3,4,5))
for (i in package_list){install.packages(i, character.only = TRUE)}

##lmms and betr discontinued in repositories, but used versions provided here in this github folder
install.packages("betr_1.32.0.tar.gz", repos = NULL, type="source")
install.packages("lmms_1.3.3.tar.gz", repos = NULL, type="source")

##Load example UPS1 Mix data
load("UPS1_Mix_Full_ExampleData.RData")

##Load the known truly changing UPS1 spike-in proteins used in the evaluation of the results
load("UPSSpikeProteins.RData")

##Define the examined conditions
sample1="Condition1"
sample2="Condition2"

##Dataset already organized in the proper format for the methods:
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
#.
#.
#.

##Define the methods
method.list=c("ROTS", "Lme", "Pme_H", "Pme_L", "LMMS", "MaSigPro_H", "MaSigPro_L",
              "BETR", "Timecourse", "Limma", "LimmaSplines_H", "LimmaSplines_L",
              "EDGE_H", "EDGE_L", "OmicsLonDA", "RolDE")

##Examine differential expression with each method and save the results
de_results=list()

for(m in 1:length(method.list)){ #Will take a while
  source(paste("Methods/", method.list[m], ".R", sep = ""))
  set.seed(1)
  de_res=getDE(data = data, sample1 = sample1, sample2 = sample2)
  de_results[[m]]=de_res
}
names(de_results)<-method.list

##Calculate partial area under the receiver operating characteristic (ROC) (pAUC)

##Function to calculate pAUC
source("getpAUC.R")

##Calculate pAUC values
pAUC_values=data.frame(matrix(nrow = length(method.list), ncol = 2), stringsAsFactors = F)
pAUC_values[,1]=method.list
rownames(pAUC_values)=method.list

for(m in 1:length(method.list)){
  set.seed(1)
  temp_pAUC<-getAUC(results = de_results[[m]], reference.prot = rownames(data), spike.prot = ups_spike, method = names(de_results)[m])
  pAUC_values[m,2]=temp_pAUC
}

colnames(pAUC_values)=c("Method", "pAUC")

##save results
save.image("Longitudinal_Differential_Expression_Example.RData")
