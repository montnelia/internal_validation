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
library(cluster)

#############
# Read data # 
#############

setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

data_base <- list.files("inst/data")
data_base
pos_data <- 2

  ######################################
  ### Preparation euclidean distance ###
  ######################################

tmp <- "biase.rds"



for(k_true in 2:10){
  for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
  #for(tmp in c("biase.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
  

  pos_data <-which(tmp ==data_base)
  D_euclidean <- readRDS(paste0("inst/euclidean/D_euclidean_", tmp))
                     
                         
  D_eucl <- readRDS(paste0("inst/eucl/D_eucl_", tmp))
  D_spear <- readRDS(paste0("inst/spear/D_spear_", tmp))
  D_pear <- readRDS(paste0("inst/pear/D_pear_", tmp))

  hierarch_cl <- hclust(D_euclidean)
  res_hclust <- cutree(hierarch_cl,k_true)
  kmeans_cl <- kmeans(D_euclidean, k_true)
  res_kmeans <- kmeans_cl$cluster
  

  
  
  
  
  
  #k_true <- readRDS(paste0("inst/groups/", data_base[pos_data]))
  data <- readRDS(paste0("inst/data/", tmp))
  data <- t(data)
  #gene_filter_data <- readRDS(paste0("inst/filtered_data/", tmp))
  
  set.seed(3)
  ### Phate ### 
  res_phate <- phateR::phate(data)
  res_cl_phate <- cluster_phate(res_phate, k= k_true)
  

  
  
  ### Umaps ### 
  res_umap <- umap(data)

  
  fit_cluster_umap_kmeans=kmeans((res_umap$layout), k_true)  
  cl_umap_kmeans = factor(fit_cluster_umap_kmeans$cluster)
  fit_cluster_hierarchical=hclust(dist(scale(res_umap$layout))) #scaling of hclust
  ## setting 3 clusters as output
  cl_umap_hclust_scaled = factor(cutree(fit_cluster_hierarchical, k=k_true))  
  
  
  
  ### tsne ###
  set.seed(9)  
  tsne_model_1 = Rtsne(as.matrix((data)), check_duplicates=FALSE, pca=TRUE, perplexity=10, theta=0.5, dims=2)
  ## getting the two dimension matrix
  d_tsne_1 = as.data.frame(tsne_model_1$Y) 
  
  
  
  
  
  ### Visualisation auslagern: mit eigener Funktion ### 
  ### Idea :::: mapping of truth and mapping of consensus values.... how often grouped together
  #d_tsne_1_plot <- cbind(d_tsne_1, group =(truth_simdata))
  #str(d_tsne_1_plot)
  #plot(d_tsne_1_plot$V1[1:(length(truth_simdata)/2)], d_tsne_1_plot$V2[1:(length(truth_simdata)/2)], col = c("red"))
  #points(d_tsne_1_plot$V1[((length(truth_simdata)/2)+1):length(truth_simdata)], d_tsne_1_plot$V2[((length(truth_simdata)/2)+1):length(truth_simdata)], col = c("blue"))
  
  
  
  #plot1 <- ggplot(d_tsne_1_plot, aes(x=V1, y=V2, color = factor(group) )) +  
   # geom_point() +
   # xlab("") + ylab("") +
   # ggtitle("t-SNE") +
   # theme_light(base_size=40) +
   # theme(axis.text.x=element_blank(),
   #       axis.text.y=element_blank()) 
  #plot1 + scale_color_brewer(palette="Dark2")
  
  
  fit_cluster_kmeans=kmeans(scale(d_tsne_1), k_true)  
  cl_tsne_kmeans = factor(fit_cluster_kmeans$cluster)
  fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1))) #scaling of hclust
  ## setting 3 clusters as output
  cl_tsne_hclust_scaled = factor(cutree(fit_cluster_hierarchical, k=k_true))  
  

  
  
  
  n_DC_min <- round(0.04*dim(D_euclidean)[1])
  n_DC_max <- ceiling(0.07*dim(D_euclidean)[1])
  n_DC <- n_DC_min:n_DC_max
  #if(length(truth_simdata) ==500){
  #  n_DC <- 20:35
  #}
  #if(length(truth_simdata) ==100){
  #  n_DC <- 4:7
  #}
  #if(length(truth_simdata) ==50){
   # n_DC <- 2:4
  #}

  
  #D_eucl <- readRDS(paste0("inst/distances/filtered/D_eucl_filt/D_eucl_filt_", data_base[pos_data]))
  #D_spear <- readRDS(paste0("inst/distances/filtered/D_spear_filt/D_spear_filt_", data_base[pos_data]))
  #D_pear <- readRDS(paste0("inst/distances/filtered/D_spear_filt/D_spear_filt_", data_base[pos_data]))
  
  
  PCA_eucl <- prcomp(D_eucl, scale = TRUE, center = TRUE)
  PCA_spear <- prcomp(D_spear, scale = TRUE, center = TRUE)
  PCA_pear <- prcomp(D_pear, scale = TRUE, center = TRUE)
  
  Lapl_eucl <- norm_laplacian(as.matrix(D_eucl)) # no scaling no centering, bzw in funktion bereits 
  Lapl_spear <- norm_laplacian(as.matrix(D_spear)) 
  Lapl_pear <- norm_laplacian(as.matrix(D_pear)) 
  
  Diffus_eucl <- as.data.frame(as.matrix(D_eucl))# no scaling no centering, bzw in funktion bereits 
  opt_sigma <- destiny::find_sigmas(Diffus_eucl, verbose =FALSE)
  Diffmap_eucl <- destiny::DiffusionMap(Diffus_eucl, n_eigs = n_DC[length(n_DC)], sigma = opt_sigma@optimal_sigma)
  opt_sigma <- destiny::find_sigmas(D_spear, verbose =FALSE)
  Diffmap_spear = destiny::DiffusionMap(D_spear, n_eigs= n_DC[length(n_DC)], sigma = opt_sigma@optimal_sigma) 
  opt_sigma <- destiny::find_sigmas(D_pear, verbose =FALSE)
  Diffmap_pear = destiny::DiffusionMap(D_pear, n_eigs= n_DC[length(n_DC)], sigma = opt_sigma@optimal_sigma) 

  
  
  
  
  ### Original version PCA and Laplace
  res_PCA_eucl <- as.data.frame(matrix(NA, length(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  res_PCA_pear <- as.data.frame(matrix(NA, length(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  res_PCA_spear <- as.data.frame(matrix(NA, length(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  
  res_Lapl_eucl <- as.data.frame(matrix(NA, length(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  res_Lapl_pear <- as.data.frame(matrix(NA, length(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  res_Lapl_spear <- as.data.frame(matrix(NA, length(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  ### Addition Diffusion Maps 
  res_diff_eucl <- as.data.frame(matrix(NA, max(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  res_diff_pear <- as.data.frame(matrix(NA, max(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  res_diff_spear <- as.data.frame(matrix(NA, max(n_DC), ncol = dim(data)[1])) # NAs entstehen wenn res 36x80
  
  
  
  
  tmp <- c()
  set.seed(1)
  
  for (d in n_DC[1]:n_DC[length(n_DC)]){ #
    res_eucl_PCA <- kmeans(PCA_eucl$rotation[,1:d], k_true, iter.max = 10^9, nstart = 1000)
    res_PCA_eucl[n_DC[length(n_DC)]-d+1, ] <- res_eucl_PCA$cluster
    res_spear_PCA <- kmeans(PCA_spear$rotation[,1:d], k_true, iter.max = 10^9, nstart = 1000)
    res_PCA_spear[n_DC[length(n_DC)]-d+1, ] <- res_spear_PCA$cluster
    res_pear_PCA <- kmeans(PCA_pear$rotation[,1:d], k_true, iter.max = 10^9, nstart = 1000)
    res_PCA_pear[n_DC[length(n_DC)]-d+1,  ] <- res_pear_PCA$cluster
    
    res_eucl_Lapl <- kmeans(eigen(Lapl_eucl)$vectors[, order(eigen(Lapl_eucl)$values)][,1:d], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_Lapl_eucl[n_DC[length(n_DC)]-d+1, ] <- res_eucl_Lapl$cluster
    res_pear_Lapl <- kmeans(eigen(Lapl_pear)$vectors[, order(eigen(Lapl_pear)$values)][,1:d], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_Lapl_pear[n_DC[length(n_DC)]-d+1, ] <- res_pear_Lapl$cluster
    res_spear_Lapl <- kmeans(eigen(Lapl_spear)$vectors[, order(eigen(Lapl_spear)$values)][,1:d], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_Lapl_spear[n_DC[length(n_DC)]-d+1, ] <- res_spear_Lapl$cluster
    
    res_eucl_diff <- kmeans(Diffmap_eucl@eigenvectors[,1:d], k_true, iter.max = 10^9, nstart = 1000)
    res_diff_eucl[n_DC[length(n_DC)]-d+1, ] <- res_eucl_diff$cluster
    res_spear_diff <- kmeans(Diffmap_spear@eigenvectors[,1:d], k_true, iter.max = 10^9, nstart = 1000)
    res_diff_spear[n_DC[length(n_DC)]-d+1, ] <- res_spear_diff$cluster
    res_pear_diff <- kmeans(Diffmap_pear@eigenvectors[,1:d], k_true, iter.max = 10^9, nstart = 1000)
    res_diff_pear[n_DC[length(n_DC)]-d+1,  ] <- res_pear_diff$cluster
    
  }
  P <- rbind(res_PCA_eucl, res_PCA_spear, res_PCA_pear)
  L <- rbind(res_Lapl_eucl, res_Lapl_spear, res_Lapl_pear)
  D_prep <- rbind(res_diff_eucl, res_diff_spear, res_diff_pear)
  na_vec <- which(is.na(D_prep[,1]))
  D <- (D_prep[-na_vec,])

  SC3 <- rbind(P, L)
  adaSC3 <- rbind(L,D)
  
  consens_mat_SC3 <- as.data.frame(matrix(NA, nrow = dim(SC3)[2], ncol = dim(SC3)[2]))
  for(k in 1:dim(SC3)[2]){
    for(j in 1:dim(SC3)[2]){
      consens_mat_SC3[k,j]<- sum(SC3[,k]== SC3[,j], na.rm = TRUE)
      consens_mat_SC3[j,k]<- sum(SC3[,k]== SC3[,j], na.rm = TRUE)
    }
  }
  
  hc <- hclust(dist(consens_mat_SC3))
  res_SC3 <- cutree(hc, k=k_true)
  
  
  
  
  consens_mat_adaSC3 <- as.data.frame(matrix(NA, nrow = dim(adaSC3)[2], ncol = dim(adaSC3)[2]))
  for(k in 1:dim(SC3)[2]){
    for(j in 1:dim(SC3)[2]){
      consens_mat_adaSC3[k,j]<- sum(adaSC3[,k]== adaSC3[,j], na.rm = TRUE)
      consens_mat_adaSC3[j,k]<- sum(adaSC3[,k]== adaSC3[,j], na.rm = TRUE)
    }
  }
  
  hc <- hclust(dist(consens_mat_adaSC3))
  res_adaSC3 <- cutree(hc, k=k_true)
  
  
  
  ss <- silhouette(res_cl_phate, dist(res_phate$embedding[,1:2]))
  sil_phate <- mean(ss[,3])
  
  ss <- silhouette(fit_cluster_kmeans$cluster, dist(scale(d_tsne_1))) # tsne 
  sil_tsne <- mean(ss[,3])
  
  ss <- silhouette(fit_cluster_umap_kmeans$cluster, dist(res_umap$layout))
  sil_umap <- mean(ss[,3])
  
  
  
    d <- n_DC[length(n_DC)]
    data_SC3 <- cbind((PCA_eucl$rotation[,1:d]),
                       (PCA_spear$rotation[,1:d]),
                       PCA_pear$rotation[,1:d],
                       
                       eigen(Lapl_eucl)$vectors[, order(eigen(Lapl_eucl)$values)][,1:d],
                       eigen(Lapl_pear)$vectors[, order(eigen(Lapl_pear)$values)][,1:d],
                       eigen(Lapl_spear)$vectors[, order(eigen(Lapl_spear)$values)][,1:d])
    
    data_adaSC3 <- cbind(Diffmap_eucl@eigenvectors[,1:d],
      Diffmap_spear@eigenvectors[,1:d],
      Diffmap_pear@eigenvectors[,1:d],
                      
                      eigen(Lapl_eucl)$vectors[, order(eigen(Lapl_eucl)$values)][,1:d],
                      eigen(Lapl_pear)$vectors[, order(eigen(Lapl_pear)$values)][,1:d],
                      eigen(Lapl_spear)$vectors[, order(eigen(Lapl_spear)$values)][,1:d])
   
    

    
    
  ss <- silhouette(res_adaSC3, dist(data_adaSC3))
  sil_adaSC3 <- mean(ss[,3])
  
  ss <- silhouette(res_SC3, dist(data_SC3))
  sil_SC3 <- mean(ss[,3])
  
  silhouette_index <- c(sil_SC3, sil_adaSC3, sil_phate, sil_tsne, sil_umap)
  names(silhouette_index) <- c("SC3", "adaSC3", "phate", "t-sne", "umap")
  
  saveRDS(silhouette_index, paste0("inst/silhouette/", k_true,"/", data_base[pos_data]))
  
  
  
  
  output_trans <- list(P = P, L = L, D= D, SC3 = as.factor(res_SC3), adaSC3 = as.factor(res_adaSC3),
                       hclust = as.factor(res_hclust), kmeans = as.factor(res_kmeans), phate = as.factor(res_cl_phate),
                       tsne_kmeans = cl_tsne_kmeans, tsne_hclust_scaled = cl_tsne_hclust_scaled,
                       umap_kmeans = cl_umap_kmeans, umap_hclust_scaled = cl_umap_hclust_scaled)
                       
                       
                       
  

  saveRDS(output_trans, paste0("inst/ind_clustering_neu/", k_true,"/", data_base[pos_data]))
  output_trans
  
  ###################
  ### Relabelling ###
  ###################
  library(cola)
  ind_clustering_labeled <- list()
  for(m in 1:3){
    
    tmp <- data.frame(matrix(NA, nrow = dim(output_trans[[m]])[1], ncol = dim(output_trans[[m]])[2]))
    
    for(i in 1:dim(output_trans[[m]])[1]){
      
      tmp[i,] <- relabel_class(output_trans[[m]][i,], output_trans$SC3, return_map = FALSE)
    }  
    ind_clustering_labeled[[m]] <- apply(tmp, 2, as.factor)
  }  
  for(m in 4:length(output_trans)){
    ind_clustering_labeled[[m]] <- relabel_class(output_trans[[m]], output_trans$SC3, return_map = FALSE)
  }
  names(ind_clustering_labeled) <- names(output_trans)
  
  
  
  
  output_relabeled <- ind_clustering_labeled
  
  #saveRDS(output_relabeled, paste0("inst/ind_clustering_labeled_neu/", k_true,"/", data_base[pos_data]))
  #output_relabeled
  
  
  
  
  }
}


### Association measures ####
library("DescTools")
library("cola")
setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

data_base <- list.files("inst/data")
data_base





for(k_true in 2:10){
  for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
    
    
    pos_data <-which(tmp ==data_base)

    ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled_neu/",k_true, "/", data_base[pos_data]))
    N <- length(ind_cluster_labelled[[4]])
    
    same_gr <- list()
    for(k in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
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
    for(i in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
      for(j in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
        #for(i in c("SC3" , "adaSC3", "kmeans", "phate")){
        # for(j in c("SC3" , "adaSC3","kmeans", "phate")){
        
        prep_method_i <- as.vector(same_gr[[i]])
        prep_method_j <- as.vector(same_gr[[j]])
        prep_method_i[-which(is.na(prep_method_i))]
        corr_output[i,j] <- cor(prep_method_i[-which(is.na(prep_method_i))], prep_method_j[-which(is.na(prep_method_j))])
      }  
    }
    saveRDS(corr_output, paste0("inst/assoc_measure_neu/" ,k_true, "/corr/", data_base[pos_data]))
  }
}

##################################
### Visualisierung Association ###
##################################
library(ggplot2)


#for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
  corr_binded <- list()
  sil_binded <- list()
  for(k_true in 2:10){
  #for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
    pos_data <-which(tmp ==data_base)
    corr_binded[[k_true-1]] <- readRDS(paste0("inst/assoc_measure_neu/" ,k_true, "/corr/", data_base[pos_data]))
    sil_binded[[k_true-1]] <- readRDS(paste0("inst/silhouette/" ,k_true, "/", data_base[pos_data]))
    
  }
  saveRDS(corr_binded, paste0("inst/assoc_binded_2_10/corr/", data_base[pos_data]))
  saveRDS(sil_binded, paste0("inst/silhouette_binded_2_10/", data_base[pos_data]))
}

for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
  pos_data <-which(tmp ==data_base)
  sil_binded_tmp <- readRDS(paste0("inst/silhouette_binded_2_10/", data_base[pos_data]))
  
  data_prep <- data.frame(matrix(NA, ncol = 2, nrow = 9*5))
  names(data_prep) <- c("clusters", "sil_value")
  group_value <- c()
  sil_value <- c()
  for(k in 1:length(sil_binded_tmp)){
    group_tmp_value <- c(rep(k+1, times =5))
    sil_tmp_value <- c(sil_binded_tmp[[k]])
    group_value <- c(group_value, group_tmp_value)
    sil_value <- c(sil_value, sil_tmp_value)
  }
  data_prep[,1] <- factor(group_value, levels = c(2:10))
  data_prep[,2] <- sil_value
  title <- sapply(strsplit(tmp, split='.', fixed=TRUE), function(x) (x[1]))
  title_ggplot <- paste0(str_to_title(title), " et al.") 
  k_true <- readRDS(paste0("inst/truth_simdata/", tmp))
  k_truth <- length(levels(factor(k_true)))
  k_truth
  p<-ggplot(data_prep, aes(x=clusters, y=sil_value)) + theme(plot.title = element_text(hjust = 0.5, size = 25))+
    geom_boxplot() + ggtitle(title_ggplot) + ylim(0,1)+ xlab("Number of clusters K") + ylab("Silhouette index") +
    geom_boxplot(data=data_prep[data_prep$clusters==k_truth,],
                 aes(x = clusters, y = sil_value),color="steelblue", fill="gray") 
  p <- p+theme(axis.text=element_text(size=20),
               axis.title=element_text(size=22))
  p
  pdf(paste0("C:/Users/ri89por/Desktop/Forschung/BIBM/Conference-LaTeX-template_10-17-19/Figuras/", title, "_number_g_sil.pdf"))
  print(p)
  dev.off()
}




for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
  pos_data <-which(tmp ==data_base)
  corr_binded_tmp <- readRDS(paste0("inst/assoc_binded_2_10/corr/", data_base[pos_data]))
  
  data_prep <- data.frame(matrix(NA, ncol = 2, nrow = 9*length(which(upper.tri(corr_binded_tmp[[1]][c(1,2,5,6,7), c(1,2,5,6,7)])))))
  names(data_prep) <- c("clusters", "corr")
  group_value <- c()
  corr_value <- c()
  for(k in 1:length(corr_binded_tmp)){
    corr_tmp <- corr_binded_tmp[[k]][c(1,2,5,6,7), c(1,2,5,6,7)][upper.tri(corr_binded_tmp[[k]][c(1,2,5,6,7), c(1,2,5,6,7)])]
    group_tmp_value <- c(rep(k+1, times =length(corr_tmp)))
    corr_tmp_value <- c(corr_tmp)
    group_value <- c(group_value, group_tmp_value)
    corr_value <- c(corr_value, corr_tmp_value)
  }
  data_prep[,1] <- factor(group_value, levels = c(2:10))
  data_prep[,2] <- corr_value
  title <- sapply(strsplit(tmp, split='.', fixed=TRUE), function(x) (x[1]))
  title_ggplot <- paste0(str_to_title(title), " et al.") 
  k_true <- readRDS(paste0("inst/truth_simdata/", tmp))
  k_truth <- length(levels(factor(k_true)))
  k_truth
  p<-ggplot(data_prep, aes(x=clusters, y=corr)) + theme(plot.title = element_text(hjust = 0.5, size = 25))+
    geom_boxplot() + ggtitle(title_ggplot) + ylim(0,1) + xlab("Number of clusters K") + ylab(expression(paste(Phi, "-Coefficient"))) +
    geom_boxplot(data=data_prep[data_prep$clusters==k_truth,],
                 aes(x = clusters, y = corr),color="steelblue", fill="gray") 
  p <- p+theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22))

  pdf(paste0("C:/Users/ri89por/Desktop/Forschung/BIBM/Conference-LaTeX-template_10-17-19/Figuras/", title, "_number_g.pdf"))
  print(p)
  dev.off()
}

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


setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

data_base <- list.files("inst/data")
data_base

for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
  
  pos_data <- which(data_base == tmp)
    k_true <- readRDS(paste0("inst/groups/", data_base[pos_data]))
  pos_data <- which(data_base == tmp)
    ind_clustering_labelled <- readRDS(paste0("inst/ind_clustering_labeled_neu/", k_true, "/", data_base[pos_data]))
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


ARI <- c()
purity <- c()
NMI <- c()
F1 <- c()

setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

data_base <- list.files("inst/data")
data_base

#tmp <- readRDS("inst/accuracy_neu/deng.rds")
#tmp <- "deng.rds"

for(tmp in c("biase.rds", "deng.rds", "goolam.rds", "yan.rds", "treutlein.rds")){
  

  acc <- readRDS(paste0("inst/accuracy_neu/", tmp))
  acc
  
  
  
  ARI_tmp <- unlist(acc[[1]][c("SC3", "adaSC3", "phate", "tsne_kmeans", "umap_kmeans"),1])
  purity_tmp <- unlist(acc[[1]][c("SC3", "adaSC3", "phate", "tsne_kmeans", "umap_kmeans"),2])
  NMI_tmp <- unlist(acc[[1]][c("SC3", "adaSC3", "phate", "tsne_kmeans", "umap_kmeans"),3])
  F1_tmp <- unlist(acc[[1]][c("SC3", "adaSC3", "phate", "tsne_kmeans", "umap_kmeans"),4])
  
  ARI <- c(ARI, ARI_tmp)
  purity <- c(purity, purity_tmp)
  NMI <- c(NMI, NMI_tmp)
  F1 <- c(F1, F1_tmp)
}

value <- c(ARI, purity, NMI, F1)
measure  <- factor((c(rep("ARI", 25), rep("Purity",25), rep("NMI", 25), rep("F1-Score", 25))))
data_base <- factor(rep(c(rep("Biase et al.", 5), rep("Deng et al.", 5), rep("Goolam et al.", 5), rep("Yan et al.", 5), rep("Treutlein et al.", 5)),4))


acc_data <- data.frame("value" = value, "measure" = measure, "data" = data_base)
acc_data

setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

data_base <- list.files("inst/data")
data_base

p1 <- ggplot(acc_data, aes(x=data, y=value, fill=measure)) + 
  geom_boxplot(fill = "gray", color = "steelblue")+ ylab("Accuracy value") + xlab("Dataset") +
  facet_wrap(~measure)  + theme(legend.position = "none")+ 
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14) )+
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+ scale_fill_discrete(name = "Accuracy measure")
p1

pdf(paste0("C:/Users/ri89por/Desktop/Forschung/BIBM/Conference-LaTeX-template_10-17-19/Figuras/acc_real.pdf"))
print(p1)
dev.off()











