setwd("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/")
getwd()

library(mclust)
library(MLmetrics) # F1_score
library(caret)
library(cola)
library(MLmetrics) #F1_Score
library(funtimes) # purity
library(reshape2)
library(aricode)



read_data <- function(i_low, i_high, level){
  assoc_weak_res<- list()
for(i in i_low:i_high){
  tmp <- readRDS(paste0("inst/assoc_measure/2/corr/simdata", level, "_", i, ".Rds"))
  assoc_weak_res[[i-i_low+1]] <- tmp
}
  assoc_weak_res
}  
assoc_weak_corr <- read_data(3,13, "weak")
assoc_moderate_corr <- read_data(28,38, "moderate")
assoc_high_corr <- read_data(40,50, "high")


read_data <- function(i_low, i_high, level){
  assoc_weak_res<- list()
  for(i in i_low:i_high){
    tmp <- readRDS(paste0("inst/assoc_measure/2/cont_coeff/simdata", level, "_", i, ".Rds"))
    assoc_weak_res[[i-i_low+1]] <- tmp
  }
  assoc_weak_res
}  
assoc_weak_cont_coeff <- read_data(3,13, "weak")
assoc_moderate_cont_coeff <- read_data(28,38, "moderate")
assoc_high_cont_coeff <- read_data(40,50, "high")


read_data <- function(i_low, i_high, level){
  assoc_weak_res<- list()
  for(i in i_low:i_high){
    tmp <- readRDS(paste0("inst/assoc_measure/2/cramer_v/simdata", level, "_", i, ".Rds"))
    assoc_weak_res[[i-i_low+1]] <- tmp
  }
  assoc_weak_res
}  
assoc_weak_cramer_v <- read_data(3,13, "weak")
assoc_moderate_cramer_v <- read_data(28,38, "moderate")
assoc_high_cramer_v <- read_data(40,50, "high")



corr_weak <- c()
corr_moderate <- c()
corr_high <- c()

cont_weak <- c()
cont_moderate <- c()
cont_high <- c()

cramer_weak <- c()
cramer_moderate <- c()
cramer_high <- c()

for(i in 1:11){
  
  prep_weak_corr <- assoc_weak_corr[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]
  prep_weak_cont <- assoc_weak_cont_coeff[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]
  prep_weak_cramer <- assoc_weak_cramer_v[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]

  
  prep_moderate_corr <- assoc_moderate_corr[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]
  prep_moderate_cont <- assoc_moderate_cont_coeff[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]
  prep_moderate_cramer <- assoc_moderate_cramer_v[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]
  
  prep_high_corr <- assoc_high_corr[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]
  prep_high_cont <- assoc_high_cont_coeff[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]
  prep_high_cramer <- assoc_high_cramer_v[[i]][c(1,2,5,6,7), c(1,2,5,6,7) ]
  
  corr_weak <- c(corr_weak, prep_weak_corr[upper.tri(prep_weak_corr)])
  corr_moderate <- c(corr_moderate, prep_moderate_corr[upper.tri(prep_moderate_corr)])
  corr_high <- c(corr_high, prep_high_corr[upper.tri(prep_high_corr)])
  
  cont_weak <- c(cont_weak, prep_weak_cont[upper.tri(prep_weak_cont)])
  cont_moderate <- c(cont_moderate, prep_moderate_cont[upper.tri(prep_moderate_cont)])
  cont_high <- c(cont_high, prep_high_cont[upper.tri(prep_high_cont)])
  
  cramer_weak <- c(cramer_weak, prep_weak_cramer[upper.tri(prep_weak_cramer)])
  cramer_moderate <- c(cramer_moderate, prep_moderate_cramer[upper.tri(prep_moderate_cramer)])
  cramer_high <- c(cramer_high, prep_high_cramer[upper.tri(prep_high_cramer)])
  
} 

correlation <- data.frame("value" = c(corr_weak, corr_moderate, corr_high),
                             "dependency" = c(rep("low", length(corr_weak)), 
                                              rep("moderate", length(corr_moderate)),
                                              rep("strong", length(corr_high))))

cont <- data.frame("value" = c(cont_weak, cont_moderate, cont_high),
                           "dependency" = c(rep("low", length(cont_weak)), 
                                            rep("moderate", length(cont_moderate)),
                                            rep("strong", length(cont_high))))

cramer <- data.frame("value" = c(cramer_weak, cramer_moderate, cramer_high),
                           "dependency" = c(rep("low", length(cramer_weak)), 
                                            rep("moderate", length(cramer_moderate)),
                                            rep("strong", length(cramer_high))))
 

ass_meas <- rbind(data.frame("ass" = correlation, "meas" = "Correlation"),
  data.frame("ass" = cont, "meas" = "Contingency \n coefficient C"),
data.frame("ass" = cramer, "meas" = "Cramér's V"))
ass_meas$meas <- factor(ass_meas$meas, levels = c("Contingency \n coefficient C", "Cramér's V", "Correlation"))
levels(ass_meas$meas)

names(ass_meas)
p<-ggplot(ass_meas, aes(x=ass.dependency, y=ass.value, fill = meas)) +
  geom_boxplot()+ theme(plot.title = element_text(hjust = 0.5)) +
ggtitle("") + xlab("Dependence structure") + ylab("Association value") + ylim(-0.05,1.05) + scale_fill_manual( labels=c("Contingency coefficient C", "Cramér's V", expression(paste(Phi, "-Coefficient"))), values = c("#1B9E77", "#7570B3", "#D95F02")) +theme(plot.title = element_text(hjust = 0.5))
  
p <- p+theme(axis.text=element_text(size=16),
             axis.title=element_text(size=18),
            legend.title=element_text(size=14),
            legend.text =element_text(size=12))+theme(legend.position = "top")  +labs(fill = "Association \n measure:") 
p

pdf(paste0("C:/Users/ri89por/Desktop/Forschung/BIBM/Conference-LaTeX-template_10-17-19/Figuras/assoc_simdata.pdf"))
print(p)
dev.off()


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


###########
### ARI ### 
###########
getwd()
readRDS("inst/accuracy/weak.Rds")
readRDS("inst/accuracy/moderate.Rds")
readRDS("inst/accuracy/high.Rds")

for(tmp in c("weak", "moderate", "high")){
  output_ari <- c()
  output_f1 <- c()
  output_pur <- c()
  output_nmi <- c()
  acc_simdata <- readRDS(paste0("inst/accuracy/", tmp, ".Rds"))
  for(k in 1:11){
    output_ari_tmp <- acc_simdata[[k]][c(1,2,5,6,8),1]
    output_ari <- c(output_ari, output_ari_tmp)
    
    output_f1_tmp <- acc_simdata[[k]][c(1,2,5,6,8),2]
    output_f1 <- c(output_f1, output_f1_tmp)
    
    output_pur_tmp <- acc_simdata[[k]][c(1,2,5,6,8),3]
    output_pur <- c(output_pur, output_pur_tmp)
    
    output_nmi_tmp <- acc_simdata[[k]][c(1,2,5,6,8),4]
    output_nmi <- c(output_nmi, output_nmi_tmp)
    
  }
  assign(paste0("acc_", tmp), cbind(output_ari, output_f1, output_pur,output_nmi))
}
acc_weak
acc_moderate
acc_high

acc_simdata <- data.frame("dependency" = rep(c(rep("low", times = 5*11), rep("moderate", times = 5*11), rep("strong", times = 5*11)), times = 4),
                          "accuracy" = c(acc_weak[,1], acc_moderate[,1], acc_high[,1],
                                          acc_weak[,2], acc_moderate[,2], acc_high[,2],
                                          acc_weak[,4], acc_moderate[,4], acc_high[,4], 
                                          acc_weak[,3], acc_moderate[,3], acc_high[,3]), 
                          "measure" = c(rep("ARI", 55*3), rep("F1-Score", 55*3), rep("NMI", 55*3), rep("Purity", 55*3)))





p1 <- ggplot(acc_simdata, aes(x=dependency, y=accuracy, fill=measure)) + 
  geom_boxplot(fill = "gray", color = "steelblue") + ylab("Accuracy value") + xlab("Dataset") +
  facet_wrap(~measure)  + theme(legend.position = "none")+ 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14) )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Dependence structure of simulation data")
#+ scale_fill_discrete(name = "Accuracy measure")
p1

pdf(paste0("C:/Users/ri89por/Desktop/Forschung/BIBM/Conference-LaTeX-template_10-17-19/Figuras/acc_simdata.pdf"))
print(p1)
dev.off()







