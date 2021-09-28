
goolam <- readRDS("datasets/goolam.rds")



library("fitdistrplus")
library(SingleCellExperiment)

data <- counts(goolam)
rowsums <- apply(data, 1, sum)
rowsums_4cell <- apply(data[, which(goolam$cell_type1=="4cell")], 1, sum)
rowsums_8cell <- apply(data[, which(goolam$cell_type1=="8cell")], 1, sum)


data_clean <- data[-c(which(rowsums==0), which(rowsums_4cell==0), which(rowsums_8cell==0)),]
dim(data_clean)


n_var <- 2000
N1 <- 50
N2 <- 50
sigma_gr1 <- cov(t(data_clean[1:n_var,which(goolam$cell_type1=="4cell") ]))


for(seed in 1:50){
  set.seed(seed)
 
  mu_gr1 <- sample(seq(5, 50, 0.0001), n_var)  
  mu_gr2 <- sample(seq(0.41, 10,0.0001), n_var)
  
  df1 <- mvrnorm(N1, mu_gr1, sigma_gr1)
  
  
  
  for(corr in c(0.1, 0.5, 0.9)){
    sigma_gr2 <- data.frame(matrix(NA, n_var, n_var))
    
    cov_pos2 <- which(upper.tri(sigma_gr2), arr.ind = TRUE)
    
    
    for(i in 1:(dim(cov_pos2)[1])){
      vars <- cov_pos2[i,]
      sigma_gr2[vars[1], vars[2]] <- corr*sqrt(sigma_gr1[vars[1],vars[1]]*sigma_gr1[vars[2], vars[2]])
      sigma_gr2[vars[2], vars[1]] <- corr*sqrt(sigma_gr1[vars[1],vars[1]]*sigma_gr1[vars[2], vars[2]])
    }
    diag(sigma_gr2) <- diag(sigma_gr1)
    
    
    set.seed(seed)
    
    df2 <- mvrnorm(n=N2,mu=mu_gr2,Sigma=sigma_gr2)  
    data_corr <- rbind(df1, df2)
    saveRDS(data_corr, paste0("inst/data/data_2000_corr", corr, "_run_seed_", seed, ".Rds"))
  }
}



