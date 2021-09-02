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
library(aricode)


#############
# Read data # 
#############

#setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results/")
getwd()

data_base <- list.files("inst/data")
data_base
pos_data <- 17
#18 high 23 middle 28 weak final 
data_base[pos_data]
data <- readRDS(paste0("inst/data/", data_base[pos_data]))


###########
### ARI ### 
###########

### To Do: include ### 
### Kmeans, hclust, PHATE, tSNE, UMAPS ### 
### Add state of the art clustering algorithms ###  

t1 <- Sys.time()

pos_data <-72

accuracy_output <- list(data.frame(matrix(NA, 5, 4)))
k_true <- readRDS(paste0("inst/groups/", data_base[pos_data]))
gene_filter_data <- readRDS(paste0("inst/filtered_data/", data_base[pos_data]))
truth_simdata <- readRDS(paste0("inst/truth_simdata/", data_base[pos_data]))
ind_clustering_labelled <- readRDS(paste0("inst/ind_clustering/", k_true, "/", data_base[pos_data]))
#ind_cluster_labelled$SC3 <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/deng_sc3.rds")
#ind_cluster_labelled$adaSC3 <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/deng_adasc3.rds")




for(m in 4:length(ind_clustering_labelled)){
  
  
  clustering_method_m <- (ind_clustering_labelled[[m]])
  clustering_method_m <- (ind_clustering_labelled$kmeans)
  accuracy_output[[1]][m-3, 1] <- adjustedRandIndex(truth_simdata, clustering_method_m)
  clustering_m <- relabel_class(clustering_method_m,as.integer(as.factor(truth_simdata)), return_map = FALSE)
  pur <- purity(clustering_m, as.integer(as.factor(truth_simdata)))
  
  
  accuracy_output[[1]][m-3, 2] <- pur$pur
  accuracy_output[[1]][m-3, 3] <- NMI(clustering_m,as.integer(as.factor(truth_simdata)) )
  
  accuracy_output[[1]][m-3, 4] <- F1_Score(y_pred = as.integer(as.factor(clustering_m)), y_true = as.integer(as.factor(truth_simdata)))
  
   
  
} 
colnames(accuracy_output[[1]]) <- c("ARI", "Purity", "NMI", "F1-Score")
rownames(accuracy_output[[1]]) <- names(ind_clustering_labelled)[4:length(ind_clustering_labelled)]
t2 <- Sys.time()
difftime(t2, t1, units = "secs")

list.files("inst/accuracy")
#saveRDS(accuracy_output[[1]], paste0("inst/accuracy/", data_base[pos_data]))

xtable(accuracy_output[[1]])
