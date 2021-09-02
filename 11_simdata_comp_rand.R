#####################
### Preprocessing ###
#####################


#https://fredclavel.org/2019/04/17/simulating-correlated-multivariate-data/
#Libraries we'll need:
library (MASS)
library(fitdistrplus)

goolam <- readRDS("C:/Users/ri89por/Downloads/goolam.rds")
data <- counts(goolam)
dim(data)
rowsums <- apply(data, 1, sum)
data_clean <- data[-which(rowsums==0),]
dim(data_clean)

param_nv <- data.frame(matrix(NA, ncol =2, nrow = dim(data_clean)[1]))
colnames(param_nv) <- c("mean", "sd")
for(i in 1:(dim(data_clean)[1])){
  param <- fitdist(data_clean[i,], "norm")
  param_nv[i,1] <- param$estimate[1]
  param_nv[i,2] <- param$estimate[2]
}

cov_mat <- data.frame(matrix(NA, 5000, 5000))

cov_pos <- which(upper.tri(cov_mat), arr.ind = TRUE)
cov_pos[1:20,]


for(i in 1:dim(cov_pos)[1]){
#for(i in 1:5){
  vars <- cov_pos[i,]
  cov_mat[vars[1], vars[2]] <- cov(data_clean[vars[1],], data_clean[vars[2],])
}  

#########################
### Creation of sigma ###
#########################

sigma1 <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/cov_1_5000.Rds")
param_nv <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/param_nv.Rds")

sd1 <- param_nv$sd[1:5000]


tm <- t(sigma1) # or t(m[-nrow(m),])
sigma1[lower.tri(sigma1)] <- tm[lower.tri(tm)]

str(sigma1)
dim(sigma1)[1]



for(i in 1:(dim(sigma1)[1])){
  sigma1[i,i] <- sd1[i]^2
}
sigma1
#sigma1[1:5, 1:5]
#sd1[1:5]^2

#saveRDS(sigma1, "C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/sigma_1_5000.Rds")

###################################
### Creation of simulation data ###
###################################
sigma <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/sigma_1_5000.Rds")


library(rWishart)
set.seed(123)
n_rep <- 10
n_var <- 5000
N1 <- 50 


set.seed(1234)
for(i in 1:10){
set.seed(i)  
aha <- rWishart::rWishart(1, (N1-1), newMat)
#saveRDS(aha, paste0("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/cov_wishart_rep_seed_1", i, ".Rds"))
}

set.seed(1234)
for(i in 1:10){
  set.seed(i)  
  mu1 <- sample(seq(0.15, 24.62, 0.0001), n_var, TRUE)
  #saveRDS(mu1, paste0("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/mu1_rep_seed_1", i, ".Rds"))
}




setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

cov_wishart <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/cov_wishart.Rds")
cov_wishart[1:5, 1:5,1]
mu1 <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/mu_ref.Rds")
#N1 <- 50 
#sigma1 <- 1/(N1-1)*cov_wishart
#str(sigma1)
#str(aha)
library(MASS)

create_simdata <- function(N1, N2, corr, n_var){
  
  for(i in 1:10){
  #for(i in 1:10){
    path <- "C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/param/"
    mu1 <- readRDS(paste0(path, "mu1_rep_seed_1", i, ".Rds"))
    #cov_wishart <- readRDS(paste0(path, "cov_wishart_rep_seed_1", i, ".Rds"))
    
    ### Simulation Group 1
    #sigma1 <- 1/(N1-1)*cov_wishart[,,1]
    #sigma1 <- 1/(N1-1)*cov_wishart
    
    
    
    
    #df1 <- data.frame(matrix(NA, ncol = n_var, nrow =N1))
   
    newMat_aim  <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/data/newMat.Rds")
    cov_vorlage <- cov(data)
    set.seed(123)
    sigma_wish <- rWishart::rWishart(1, (N1-1), newMat_aim)
    tmp <- sigma_wish[,,1]
    tmp[1:5, 1:5]
    tmp_div <- tmp/((N1-1))
    #tmp_div[1:5, 1:5]
    
    
    sigma1 <- cov2cor(tmp_div)
  
    sigma1 <- newMat_aim

    set.seed(12345)#123
    sd1 <- sample(seq(0.72, 32.83,0.0001), 5000 )
    for(i in 1:length(diag(sigma1))){
      diag(sigma1)[i] <- sd1[i]^2
    }
    sigma_cor <- cov2cor(sigma1)
    sigma_cor[1:5, 1:5]
    newMat_aim[1:5, 1:5]
    
    newEig <- eigen(sigma_cor)
    newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
    # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
    # eig vectors
    sigma_cor <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)
    # normalize modified matrix eqn 6 from Brissette et al 2007
    sigma_cor <- sigma_cor/sqrt(diag(sigma_cor) %*% t(diag(sigma_cor)))
    
    set.seed(123)
    
    df1 <- mvrnorm(n=N1,mu=mu1,Sigma=sigma_cor)  #simulate the data, as specified above
    
    
    ### Simulation Group 2 
    mu2 <- mu1
    sigma2 <- matrix(NA, dim(sigma1), dim(sigma1))
    
    cov_pos2 <- which(upper.tri(sigma2),        # Set arr.ind = TRUE
                      arr.ind = TRUE)
    
    for(n in 1:(dim(cov_pos2)[1])){
      vars <- cov_pos2[n,]
      sigma2[vars[1], vars[2]] <- corr*sqrt(sigma1[vars[1],vars[1]]*sigma1[vars[2], vars[2]]) 
      sigma2[vars[2], vars[1]] <- corr*sqrt(sigma1[vars[1],vars[1]]*sigma1[vars[2], vars[2]]) 
    } 
    
    diag(sigma2) <- diag(sigma1)
    df2 <- mvrnorm(n=N2,mu=mu2,Sigma=sigma2)  
    data_correlated <- rbind(df1, df2)
    #return((data_correlated))
    saveRDS(data_correlated, paste0("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/data/simdata_corr_", corr, "_rep_seed_1_test", i, ".Rds"))
    
  }

}



create_simdata(50,50, 0.1, 5000)
create_simdata(50,50, 0.4, 5000)
create_simdata(50,50, 0.9, 5000)
   




