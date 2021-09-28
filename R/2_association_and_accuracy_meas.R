

###################################################
### Association measures and accuracy measures ####
###################################################
library("DescTools")
library("cola")

setwd("C:/Users/ri89por/Desktop/Github/internal_validation/")
getwd()

data_base <- list.files("data")
pos_data <-1

############################
### Association measures ###
############################

assoc_meas <- function(pos_data){
  for(k_true in 2:10){
    tmp <- data_base[pos_data]
    pos_data <-which(tmp ==data_base)
    
    ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled_neu/",k_true, "/", data_base[pos_data]))
    N <- length(ind_cluster_labelled[[4]])
    
    same_gr <- list()
    for(k in c("SC3" , "adaSC3", "phate", "tsne_kmeans", "umap_kmeans" )){
      same_grouping <- data.frame(matrix(NA,N, N))
      for(i in 1:N){
        for(j in 1:N){
          same_grouping[i,j] <- ind_cluster_labelled[[k]][i] == ind_cluster_labelled[[k]][j]
        }
      }
      diag(same_grouping) <- NA
      same_gr[[k]] <- apply(same_grouping,1,as.integer)
    }
    
    
    corr_output <- data.frame()
    for(i in c("SC3" , "adaSC3",  "phate", "tsne_kmeans", "umap_kmeans" )){
      for(j in c("SC3" , "adaSC3", "phate", "tsne_kmeans", "umap_kmeans" )){
        
        prep_method_i <- as.vector(same_gr[[i]])
        prep_method_j <- as.vector(same_gr[[j]])
        prep_method_i[-which(is.na(prep_method_i))]
        corr_output[i,j] <- cor(prep_method_i[-which(is.na(prep_method_i))], prep_method_j[-which(is.na(prep_method_j))])
      }  
    }
    saveRDS(corr_output, paste0("inst/assoc_measure_neu/corr/" ,k_true, "/", data_base[pos_data]))
    
    
    CramerV_output <- data.frame()
    for(i in c("SC3" , "adaSC3", "phate", "tsne_kmeans", "umap_kmeans" )){
      for(j in c("SC3" , "adaSC3", "phate", "tsne_kmeans", "umap_kmeans" )){
       
      
        CramerV_output[i,j] <- CramerV(table(ind_cluster_labelled[[i]], ind_cluster_labelled[[j]]))
      }  
    }
    saveRDS(CramerV_output, paste0("inst/assoc_measure_neu/cramer_v/", k_true,"/", data_base[pos_data]))
    
    
    cont_coeff_corr <- data.frame() # hier Kontingenzkoeffizient nach Pearson:
    for(i in c("SC3" , "adaSC3", "phate", "tsne_kmeans", "umap_kmeans" )){
      for(j in c("SC3" , "adaSC3", "phate", "tsne_kmeans", "umap_kmeans" )){
        cont_coeff_corr[i,j] <- ContCoef(table(ind_cluster_labelled[[i]], ind_cluster_labelled[[j]]), correct= TRUE)
      }  
    }
    saveRDS(cont_coeff_corr, paste0("inst/assoc_measure_neu/cont_coeff/", k_true, "/", data_base[pos_data]))

  }
}
assoc_meas(1)
readRDS(paste0("inst/assoc_measure_neu/corr/3/", data_base[pos_data]))
readRDS(paste0("inst/assoc_measure_neu/cramer_v/3/", data_base[pos_data]))
readRDS(paste0("inst/assoc_measure_neu/cont_coeff/3/", data_base[pos_data]))



################
### Accuracy ###
################

library(mclust)
library(MLmetrics) #F1_Score
library(funtimes) # purity
#install.phate(envname = "r-reticulate")
library(phateR)
library(umap)
library(Rtsne)
library(aricode)

getwd()


acc_meas <- function(pos_data){

    
  k_true <- readRDS(paste0("inst/groups/", data_base[pos_data]))
  ind_clustering_labelled <- readRDS(paste0("inst/ind_clustering_labeled_neu/", k_true, "/", data_base[pos_data]))
  tmp <- data_base[pos_data]
  truth_simdata <- readRDS(paste0("inst/truth_simdata/", tmp))
  
  accuracy_output <- list(data.frame(matrix(NA, 5, 4)))
  
  for(m in 4:length(ind_clustering_labelled)){
    
    
    clustering_method_m <- (ind_clustering_labelled[[m]])
    accuracy_output[[1]][m-3, 1] <- adjustedRandIndex(truth_simdata, clustering_method_m)
    clustering_m <- relabel_class(clustering_method_m,as.integer(as.factor(truth_simdata)), return_map = FALSE)
    pur <- purity(clustering_m, as.integer(as.factor(truth_simdata)))
    
    
    accuracy_output[[1]][m-3, 2] <- pur$pur
    accuracy_output[[1]][m-3, 3] <- NMI(clustering_m,as.integer(as.factor(truth_simdata)) )
    
    accuracy_output[[1]][m-3, 4] <- F1_Score(y_pred = as.integer(as.factor(clustering_m)), y_true = as.integer(as.factor(truth_simdata)))
    
    
    
  } 
  colnames(accuracy_output[[1]]) <- c("ARI", "Purity", "NMI", "F1-Score")
  rownames(accuracy_output[[1]]) <- names(ind_clustering_labelled)[4:length(ind_clustering_labelled)]
  saveRDS(accuracy_output, paste0("inst/accuracy_neu/", tmp))
}
acc_meas(1)
readRDS(paste0("inst/accuracy_neu/", tmp))



