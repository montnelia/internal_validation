#################
### Read Data ###
#################
library("DescTools")
setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

data_base <- list.files("inst/data")
data_base
pos_data <- 18 #high
pos_data <- 23 # middle
pos_data <- 28 # weak


############################
### Association measures ###
############################

# 1. Correlation 

### transformation needed of how often each pair is grouped together, ### 
### such that the pearson correlation can be applied as at least interval scaled is needed ####

assoc_corr <- function(pos_data, k_true){
  for(pos_data in pos_data){ # weak, middle, high
    ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled/",k_true, "/", data_base[pos_data]))
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
   saveRDS(corr_output, paste0("inst/assoc_measure/" ,k_true, "/corr/", data_base[pos_data]))
  }
}
assoc_corr(20)
assoc_corr(c(20, 27, 34),2)
assoc_corr(c(18, 23, 28),2)
assoc_corr(c(19, 25, 31),2)
assoc_corr(c(18, 23, 28),3)
assoc_corr(c(18, 23, 28),4)

readRDS(paste0("inst/assoc_measure/2/corr/", data_base[18]))
readRDS(paste0("inst/assoc_measure/3/corr/", data_base[18]))
readRDS(paste0("inst/assoc_measure/4/corr/", data_base[18]))

readRDS(paste0("inst/assoc_measure/2/corr/", data_base[23]))
readRDS(paste0("inst/assoc_measure/3/corr/", data_base[23]))
readRDS(paste0("inst/assoc_measure/4/corr/", data_base[23]))

readRDS(paste0("inst/assoc_measure/2/corr/", data_base[28]))
readRDS(paste0("inst/assoc_measure/3/corr/", data_base[28]))
readRDS(paste0("inst/assoc_measure/4/corr/", data_base[28]))




high <- readRDS(paste0("inst/assoc_measure/2/corr/", data_base[20]))
middle <- readRDS(paste0("inst/assoc_measure/2/corr/", data_base[27]))
low <- readRDS(paste0("inst/assoc_measure/2/corr/", data_base[34]))
apply(high, 1, mean)
apply(middle, 1, mean)
apply(low, 1, mean)

corr_tot <- rbind(low, middle, high)
N_meth <- length(rownames(low))

rownames(corr_tot) <- c(paste0(rownames(low), rep(" (L)", N_meth)), paste0(rownames(low), rep(" (M)", N_meth)), paste0(rownames(low), rep(" (H)", N_meth)))
xtable(corr_tot)


corr_tot_k <- cbind(rep("2", N_meth), corr_tot)
rownames(corr_tot_k) <- c(paste0(rownames(low), rep(" (L)", N_meth)), paste0(rownames(low), rep(" (M)", N_meth)), paste0(rownames(low), rep(" (H)", N_meth)))
colnames(corr_tot_k)[1] <- "k"
xtable(corr_tot_k)



#2. \Chi^2 based association measures based on the nominal level !!! 






library("DescTools")


assoc_cont_coeff <- function(pos_data){
  for(pos_data in pos_data){ # weak, middle, high
    ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled/2/", data_base[pos_data]))
    N <- length(ind_cluster_labelled[[4]])
    
    cont_coeff_corr <- data.frame() # hier Kontingenzkoeffizient nach Pearson:
    for(i in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
      for(j in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
        # for(i in c("SC3" , "adaSC3", "kmeans", "phate")){
        #  for(j in c("SC3" , "adaSC3","kmeans", "phate")){  
        cont_coeff_corr[i,j] <- ContCoef(table(ind_cluster_labelled[[i]], ind_cluster_labelled[[j]]), correct= TRUE)
      }  
    }
    saveRDS(cont_coeff_corr, paste0("inst/assoc_measure/2/cont_coeff/", data_base[pos_data]))
  }
}
assoc_cont_coeff(c(20,27,34))

high_cont <- readRDS(paste0("inst/assoc_measure/2/cont_coeff/", data_base[20]))
middle_cont <- readRDS(paste0("inst/assoc_measure/2/cont_coeff/", data_base[27]))
low_cont <- readRDS(paste0("inst/assoc_measure/2/cont_coeff/", data_base[34]))
apply(high_cont, 1, mean)
apply(middle_cont, 1, mean)
apply(low_cont, 1, mean)

cont_tot <- rbind(low_cont, middle_cont, high_cont)
N_meth <- length(rownames(low_cont))
rownames(cont_tot) <- c(paste0(rownames(low_cont), rep(" (L)", N_meth)), paste0(rownames(low_cont), rep(" (M)", N_meth)), paste0(rownames(low_cont), rep(" (H)", N_meth)))
xtable(cont_tot)





assoc_cramer <- function(pos_data){
  for(pos_data in pos_data){ # weak, middle, high
    ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled/2/", data_base[pos_data]))
    N <- length(ind_cluster_labelled[[4]])
    
    
CramerV_output <- data.frame()
for(i in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
  for(j in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
    #for(i in c("SC3" , "adaSC3", "kmeans", "phate")){
    # for(j in c("SC3" , "adaSC3","kmeans", "phate")){



    #for(i in c("SC3" , "adaSC3", "kmeans", "phate")){
     # for(j in c("SC3" , "adaSC3","kmeans", "phate")){
    CramerV_output[i,j] <- CramerV(table(ind_cluster_labelled[[i]], ind_cluster_labelled[[j]]))
  }  
}
saveRDS(CramerV_output, paste0("inst/assoc_measure/2/cramer_v/", data_base[pos_data]))
  }
}
assoc_cramer(c(20,27, 34))

high_cram <- readRDS(paste0("inst/assoc_measure/2/cramer_v/", data_base[20]))
middle_cram <- readRDS(paste0("inst/assoc_measure/2/cramer_v/", data_base[27]))
low_cram <- readRDS(paste0("inst/assoc_measure/2/cramer_v/", data_base[34]))
apply(high_cram, 1, mean)
apply(middle_cram, 1, mean)
apply(low_cram, 1, mean)

cram_tot <- rbind(low_cram, middle_cram, high_cram)
N_meth <- length(rownames(low_cram))
rownames(cram_tot) <- c(paste0(rownames(low_cram), rep(" (L)", N_meth)), paste0(rownames(low_cram), rep(" (M)", N_meth)), paste0(rownames(low_cram), rep(" (H)", N_meth)))
xtable(cram_tot)





assoc_phi <- function(pos_data){
  for(pos_data in pos_data){ # weak, middle, high
    ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled/", data_base[pos_data]))
    N <- length(ind_cluster_labelled[[4]])
    
    phi_output <- data.frame()
    for(i in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
     for(j in c("SC3" , "adaSC3", "hclust", "kmeans", "phate", "tsne_kmeans", "umap_kmeans" )){
     phi_output[i,j] <- Phi(table(ind_cluster_labelled[[i]], ind_cluster_labelled[[j]]))
      }  
    }
  saveRDS(phi_output, paste0("inst/assoc_measure/phi/", data_base[pos_data]))
  }
}
assoc_phi(c(17, 22, 27))

high_phi <- readRDS(paste0("inst/assoc_measure/phi/", data_base[17]))
middle_phi <- readRDS(paste0("inst/assoc_measure/phi/", data_base[22]))
low_phi <- readRDS(paste0("inst/assoc_measure/phi/", data_base[27]))


phi_output # im Fall der Vierfelder Tafel ist Cramer's V gleich dem Phi Koeffizienten 
#signierter Phi Koeffizient nur bei dichotomen Merkmalen möglich 
#signierter Phi Koeffizient 
#phi_s > 0: gleichsinniger / positiver Zusammenhang (Konkordanz): 1 (bzw. 2) bei X
#bewirkt eher 1 (bzw. 2) bei Y






#########################
### Accuracy measures ###
#########################

pos_data <- 20 #high
pos_data <- 27 # middle
pos_data <- 34 # weak

acc_high_corr <- as.data.frame(readRDS(paste0("inst/accuracy/", data_base[20])))
acc_moderate_corr <- as.data.frame(readRDS(paste0("inst/accuracy/", data_base[27])))
acc_low_corr <- as.data.frame(readRDS(paste0("inst/accuracy/", data_base[34])))


acc_tot <- rbind(acc_low_corr, acc_moderate_corr, acc_high_corr)
N_meth <- length(rownames(acc_low_corr))
rownames(acc_tot) <- c(paste0(rownames(acc_low_corr), rep(" (L)", N_meth)), paste0(rownames(acc_low_corr), rep(" (M)", N_meth)), paste0(rownames(acc_low_corr), rep(" (H)", N_meth)))
xtable(acc_tot)


acc_LMH <- data.frame(matrix(NA, 3*dim(acc_high_corr)[1], dim(acc_high_corr)[2]+2) )
acc_LMH[,2] <- rep(c("(L)", "(M)", "(H)"), times = dim(acc_high_corr)[1])

for(i in 1:dim(acc_high_corr)[1]){
  pos <- which(is.na(acc_LMH), arr.ind = TRUE)[1,1]
  acc_LMH[pos:(pos+2),1] <- rep(rownames(acc_low_corr)[i], 3)
  acc_LMH[pos:(pos+2),3:7] <- rbind(acc_low_corr[i,], acc_moderate_corr[i,], acc_high_corr[i,])
  }
rownames(acc_LMH) <- paste(acc_LMH[,1], acc_LMH[,2])
acc_LMH <- acc_LMH[,-c(1,2)]
colnames(acc_LMH)<- colnames(acc_low_corr)

xtable(acc_LMH)


# Vorteil Assoziationsmaße auch möglich mit k > 2 
# ordinal scaled: Rangkorrelation möglich ?

##RR oder OR etc. könnte man mit aufnehmen? Nein, da keine abhängigen und unabhängigen Variablen zur Verfügung stehen
##PRE: Goodman / Kruskals Lambda vs. Tau ...


# association measures for k = 2,3,4
#for each simulation data with ground truth  = 2,3,4 and weak, middle, high correlation
### Page 1: SImulation with two groups and association measures for 2,3,4
### Page 2: SImulation with three groups and association measures for 2,3,4 
### Page 3: Simulation with four groups and association measures 2,3,4 
#Respectively in each table: High, Middle, Low 





