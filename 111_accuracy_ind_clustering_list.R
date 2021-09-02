

setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

library(mclust)
library(MLmetrics) # F1_score
library(caret)
library(cola)
library(MLmetrics) #F1_Score
library(funtimes) # purity

ind_clust_0.1corr <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/ind_clustering/ind_clust_0.1corr.Rds")
ind_clust_0.4corr <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/ind_clustering/ind_clust_0.4corr.Rds")
ind_clust_0.5corr <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/ind_clustering/ind_clust_0.5corr.Rds")
ind_clust_0.9corr <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/ind_clustering/ind_clust_0.9corr.Rds")
 

cluster_relabelling <- function(ind_clust_corr, level, i_low, i_high, k_true){
  #truth_simdata <- readRDS(paste0("inst/truth_simdata/", data_base[pos_data]))
  #ind_clust_corr <- ind_clust_0.1corr
  output_relabeled <- list()
  ind_clustering_labeled <- list()
  for(i in i_low:i_high){
    #ind_clust_corr <- ind_clust_0.1corr
    ind_clustering <- ind_clust_corr[[i]][4:12]
    

  for(m in 1:length(ind_clustering)){
    ind_clustering_labeled[[m]] <- relabel_class(ind_clustering[m][[1]], ind_clustering$SC3, return_map = FALSE)
  }
  names(ind_clustering_labeled) <- names(ind_clustering)
  
  output_relabeled[[i]] <- ind_clustering_labeled
  }
  saveRDS(output_relabeled, paste0("inst/ind_clustering_labeled/",k_true, "/", level, ".Rds"))
  output_relabeled
}
ind_clust_0.1corr_relabeled <- cluster_relabelling(ind_clust_0.1corr, "weak", 3, 13, 2)
ind_clust_0.4corr_relabeled <- cluster_relabelling(ind_clust_0.4corr, "moderate", 16, 26, 2)
ind_clust_0.5corr_relabeled <- cluster_relabelling(ind_clust_0.5corr, "moderate", 28, 38, 2)
ind_clust_0.9corr_relabeled <- cluster_relabelling(ind_clust_0.9corr, "high", 40, 50, 2)


###########
### ARI ### 
###########

ind_clust_corr <- ind_clust_0.1corr_relabeled 
i_low <- 3
i_high <- 13
acc_list <- function(ind_clust_corr_label, level, i_low, i_high){
  acc <- list()
  for(i in  i_low:i_high){
    
  ind_clustering_labelled <- ind_clust_corr_label[[i]]

accuracy_output <- (data.frame(matrix(NA, 5, 5)))
k_true <- 2
truth_simdata <- factor(c(rep(1, 50), rep(2, 50)))
#ind_clustering_labelled <- readRDS(paste0("inst/ind_clustering_labeled/2/", data_base[pos_data]))
for(m in 1:length(ind_clustering_labelled)){
  
  
  clustering_method_m <- factor(ind_clustering_labelled[[m]])
  accuracy_output[m, 1] <- adjustedRandIndex(truth_simdata, clustering_method_m)
  accuracy_output[m, 2] <- F1_Score(truth_simdata, clustering_method_m)
  pur <- purity(factor(truth_simdata), factor(clustering_method_m))
  accuracy_output[m, 3] <- pur$pur
  
  result <- caret::confusionMatrix(as.factor(clustering_method_m), as.factor(truth_simdata), mode = "prec_recall")
  accuracy_output[m, 4] <- result$byClass["Precision"]
  accuracy_output[m, 5] <- result$byClass["Recall"]
  
}
colnames(accuracy_output) <- c("ARI", "F1", "Purity", "Precision", "Recall")
rownames(accuracy_output) <- names(ind_clustering_labelled)[1:length(ind_clustering_labelled)]
acc[[i-i_low + 1 ]] <- accuracy_output
  }
  saveRDS(acc, paste0("inst/accuracy/", level,".Rds"))
  print(acc)
}
acc_0.1corr <- acc_list(ind_clust_0.1corr_relabeled, "weak", 3, 13)
acc_0.4corr <- acc_list(ind_clust_0.4corr_relabeled, "moderate", 16, 26)
acc_0.5corr <- acc_list(ind_clust_0.5corr_relabeled, "moderate", 28, 38)
acc_0.9corr <- acc_list(ind_clust_0.9corr_relabeled, "high", 40, 50)

######################
### Data frame 0.1 ###
######################
acc_0.1 <- mapply(c, acc_0.1corr, acc_0.1corr, SIMPLIFY=FALSE)

SC3 <- c()
adaSC3 <- c()
hclust <- c()
kmeans <- c()
phate <- c()
tsne_kmeans <- c()
tsne_hclust_scaled <- c()
umap_kmeans <- c()
umap_hclust_scaled <- c()






for(i in 1:length(acc_0.1)){
  SC3[i] <- acc_0.1[[i]]$ARI[1]
  adaSC3[i] <- acc_0.1[[i]]$ARI[2]
  hclust[i] <- acc_0.1[[i]]$ARI[3]
  kmeans[i] <- acc_0.1[[i]]$ARI[4]
  phate[i] <- acc_0.1[[i]]$ARI[5]
  tsne_kmeans[i] <- acc_0.1[[i]]$ARI[6]
  tsne_hclust_scaled[i] <- acc_0.1[[i]]$ARI[7]
  umap_kmeans[i] <- acc_0.1[[i]]$ARI[8]
  umap_hclust_scaled[i] <- acc_0.1[[i]]$ARI[9]
}


accuracy_df_weak <- data.frame(matrix(NA, nrow = 11, ncol = 10))
accuracy_df_weak[,10] <- "weak" 
accuracy_df_weak[,1] <- SC3 
accuracy_df_weak[,2] <- adaSC3
accuracy_df_weak[,3] <- hclust
accuracy_df_weak[,4] <- kmeans
accuracy_df_weak[,5] <- phate
accuracy_df_weak[,6] <- tsne_kmeans
accuracy_df_weak[,7] <- tsne_hclust_scaled 
accuracy_df_weak[,8] <- umap_kmeans
accuracy_df_weak[,9] <- umap_hclust_scaled
colnames(accuracy_df_weak) <- c("SC3", "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", 
                                "tsne_hclust_scaled", "umap_kmeans", "umap_hclust_scaled", "dependency")




######################
### Data frame 0.5 ###
######################
acc_0.5 <- mapply(c, acc_0.5corr, acc_0.5corr, SIMPLIFY=FALSE)

SC3 <- c()
adaSC3 <- c()
hclust <- c()
kmeans <- c()
phate <- c()
tsne_kmeans <- c()
tsne_hclust_scaled <- c()
umap_kmeans <- c()
umap_hclust_scaled <- c()






for(i in 1:length(acc_0.5)){
  SC3[i] <- acc_0.5[[i]]$ARI[1]
  adaSC3[i] <- acc_0.5[[i]]$ARI[2]
  hclust[i] <- acc_0.5[[i]]$ARI[3]
  kmeans[i] <- acc_0.5[[i]]$ARI[4]
  phate[i] <- acc_0.5[[i]]$ARI[5]
  tsne_kmeans[i] <- acc_0.5[[i]]$ARI[6]
  tsne_hclust_scaled[i] <- acc_0.5[[i]]$ARI[7]
  umap_kmeans[i] <- acc_0.5[[i]]$ARI[8]
  umap_hclust_scaled[i] <- acc_0.5[[i]]$ARI[9]
}


accuracy_df_moderate <- data.frame(matrix(NA, nrow = 11, ncol = 10))
accuracy_df_moderate[,10] <- "moderate" 
accuracy_df_moderate[,1] <- SC3 
accuracy_df_moderate[,2] <- adaSC3
accuracy_df_moderate[,3] <- hclust
accuracy_df_moderate[,4] <- kmeans
accuracy_df_moderate[,5] <- phate
accuracy_df_moderate[,6] <- tsne_kmeans
accuracy_df_moderate[,7] <- tsne_hclust_scaled 
accuracy_df_moderate[,8] <- umap_kmeans
accuracy_df_moderate[,9] <- umap_hclust_scaled
colnames(accuracy_df_moderate) <- c("SC3", "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", 
                                "tsne_hclust_scaled", "umap_kmeans", "umap_hclust_scaled", "dependency")







######################
### Data frame 0.9 ###
######################
acc_0.9 <- mapply(c, acc_0.9corr, acc_0.9corr, SIMPLIFY=FALSE)

SC3 <- c()
adaSC3 <- c()
hclust <- c()
kmeans <- c()
phate <- c()
tsne_kmeans <- c()
tsne_hclust_scaled <- c()
umap_kmeans <- c()
umap_hclust_scaled <- c()






for(i in 1:length(acc_0.9)){
  SC3[i] <- acc_0.9[[i]]$ARI[1]
  adaSC3[i] <- acc_0.9[[i]]$ARI[2]
  hclust[i] <- acc_0.9[[i]]$ARI[3]
  kmeans[i] <- acc_0.9[[i]]$ARI[4]
  phate[i] <- acc_0.9[[i]]$ARI[5]
  tsne_kmeans[i] <- acc_0.9[[i]]$ARI[6]
  tsne_hclust_scaled[i] <- acc_0.9[[i]]$ARI[7]
  umap_kmeans[i] <- acc_0.9[[i]]$ARI[8]
  umap_hclust_scaled[i] <- acc_0.9[[i]]$ARI[9]
}


accuracy_df_high <- data.frame(matrix(NA, nrow = 11, ncol = 10))
accuracy_df_high[,10] <- "high" 
accuracy_df_high[,1] <- SC3 
accuracy_df_high[,2] <- adaSC3
accuracy_df_high[,3] <- hclust
accuracy_df_high[,4] <- kmeans
accuracy_df_high[,5] <- phate
accuracy_df_high[,6] <- tsne_kmeans
accuracy_df_high[,7] <- tsne_hclust_scaled 
accuracy_df_high[,8] <- umap_kmeans
accuracy_df_high[,9] <- umap_hclust_scaled
colnames(accuracy_df_high) <- c("SC3", "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", 
                                    "tsne_hclust_scaled", "umap_kmeans", "umap_hclust_scaled", "dependency")





accuracy_ARI_tot <- rbind(accuracy_df_weak, accuracy_df_moderate, accuracy_df_high)
accuracy_ARI_tot

accuracy_ARI_tot_long <- melt(accuracy_ARI_tot, id= "dependency")

ARI_plot <- ggplot(accuracy_ARI_tot_long, aes(x = variable, y = value, fill = dependency)) + geom_boxplot() +
ggtitle("ARI")

ARI_plot










ARI <- rep(NA, 11)
for(i in 1:length(acc_0.1corr)){
  i<-1
  m<-1
for(m in 1:dim(acc_0.1corr[[1]])[1]){
  assign(paste0(rownames(acc_0.1corr[[1]])[m], "_ARI_", m ), acc_0.1corr[[i]]$ARI[m])
  }
}

for(m in 1:dim(acc_0.1corr[[1]])[1]){
  umap_hclust_scaled_ARI[m] <- c(umap_hclust_scaled_ARI_1, )
}

umap_hclust_scaled_ARI

SC3_ARI <- acc_0.1corr[[1]]$ARI[1]
adaSC3_ARI <- acc_0.1corr[[1]]$ARI[1]
acc_0.1corr[[1]]$ARI[1]
acc_0.1corr[[1]]$ARI[1]
acc_0.1corr[[1]]$ARI[1]
acc_0.1corr[[1]]$ARI[1]

for(i in 1:11){
  c()
  accuracy_df$ARI[i] <- acc_dep[[i]]$ARI[1]
}
acc_0.1corr
boxplot(acc_0.1corr)





list.files("inst/accuracy")
saveRDS(accuracy_output, paste0("inst/accuracy/", data_base[pos_data]))
#for(pos_data in c(31, 25, 19)){ # 5000
 # for(pos_data in c(17,23,29)){ # 2000
    #for(pos_data in c(17,23,29)){ #1000
  for(pos_data in c(20,27,34)){ # 5000 mit seed
  print(readRDS(paste0("inst/accuracy/", data_base[pos_data])))
}
### NaN when comparison to first clustering algorithm alternative: comparision to simulation data..... 
