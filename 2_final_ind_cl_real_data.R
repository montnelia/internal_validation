#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("destiny")


library(mclust)

library(caret)
library("destiny")
library(SC3)
library(SingleCellExperiment)
library(diffusionMap)
library(pcaMethods)
library(lle)
library(vegan)
library(RANN)
library(mclust)
library(MLmetrics) #F1_Score
library(funtimes) # purity
#install.phate(envname = "r-reticulate")
library(phateR)
library(umap)
library(Rtsne)


#############
# Read data # 
#############

#setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results/")
getwd()

data_base <- list.files("inst/data")
data_base
pos_data <- 2
#18 high 23 middle 28 weak final 
data_base[pos_data]
data <- readRDS(paste0("inst/data/", data_base[pos_data]))


###########
### ARI ### 
###########

### To Do: include ### 
### Kmeans, hclust, PHATE, tSNE, UMAPS ### 
### Add state of the art clustering algorithms ###  
pos_data <-3

accuracy_output <- list(data.frame(matrix(NA, 5, 5)))
k_true <- readRDS(paste0("inst/groups/", data_base[pos_data]))
gene_filter_data <- readRDS(paste0("inst/filtered_data/", data_base[pos_data]))
truth_simdata <- readRDS(paste0("inst/truth_simdata/", data_base[pos_data]))
ind_clustering_labelled <- readRDS(paste0("inst/ind_clustering/", k_true, "/", data_base[pos_data]))
for(m in 4:length(ind_clustering_labelled)){
  
  
  clustering_method_m <- (ind_clustering_labelled[[m]])
  accuracy_output[[1]][m-3, 1] <- adjustedRandIndex(truth_simdata, clustering_method_m)
  var1 <- as.factor(truth_simdata)
  var2 <- as.factor(clustering_method_m)
  #accuracy_output[[1]][m-3, 2] <- F1_Score(var1, var2)
  #pur <- purity(factor(truth_simdata), factor(clustering_method_m))
  #accuracy_output[[1]][m-3, 3] <- pur$pur
  
  #result <- caret::confusionMatrix(as.factor(clustering_method_m), as.factor(truth_simdata), mode = "prec_recall")
  #accuracy_output[[1]][m-3, 4] <- result$byClass["Precision"]
  #accuracy_output[[1]][m-3, 5] <- result$byClass["Recall"]
  
}
colnames(accuracy_output[[1]]) <- c("ARI", "F1", "Purity", "Precision", "Recall")
rownames(accuracy_output[[1]]) <- names(ind_clustering_labelled)[4:length(ind_clustering_labelled)]

list.files("inst/accuracy")
saveRDS(accuracy_output, paste0("inst/accuracy/", data_base[pos_data]))
#for(pos_data in c(31, 25, 19)){ # 5000
 # for(pos_data in c(17,23,29)){ # 2000
    #for(pos_data in c(17,23,29)){ #1000
  for(pos_data in c(20,27,34)){ # 5000 mit seed
  print(readRDS(paste0("inst/accuracy/", data_base[pos_data])))
}
### NaN when comparison to first clustering algorithm alternative: comparision to simulation data..... 
