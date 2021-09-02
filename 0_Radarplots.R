library(fmsb)
library(tikzDevice)
library(reshape2)
library(ggplot2)
library(dplyr)

citation("fmsb")
citation("ggplot2")


citation("fmsb")
citation("tikzDevice")
citation("reshape2")
citation("ggplot2")
citation("dplyr")
citation("cola")

?citation
output_trans <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/ind_clustering_labeled_neu/3/biase.rds")


data_radarchart <- c()
for(l in 4:12){
  #prop.table(table(output_trans[[l]]))
  data_radarchart <- rbind(data_radarchart, prop.table(table(output_trans[[l]])))
}

data_radarchart_prep <- cbind(factor(c(names(output_trans)[4:12]), levels =names(output_trans)[4:12]), data_radarchart[,1:3])
rownames(data_radarchart) <- c(names(output_trans)[4:5], "hclust","k-means", "phate", "t-sne + k-means", "t-sne:hclust", "umap + k-means", "umap:hclust")

biase_truth <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/truth_simdata/biase.Rds")
ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled_neu/3/biase.rds"))
biase <- relabel_class(biase_truth, ind_cluster_labelled$SC3, return_map = FALSE)
truth <- prop.table(table(as.factor(biase)))




data_radarchart_ready <- data.frame(rbind(truth, data_radarchart[c(1,2,6,5,8),]))



Method <- factor(c(rownames(data_radarchart_ready)),
                levels = c(rownames(data_radarchart_ready)))
levels(Method)[6]<- "umap + \n k-means" 
levels(Method)[4] <- "t-SNE + \n k-means"
levels(Method)[5] <- "   phate"
levels(Method)[1] <- "truth \n"
cluster1 <- data_radarchart_ready[,1]
cluster2 <- data_radarchart_ready[,2]
cluster3 <- data_radarchart_ready[,3]

# data preparation
df = data.frame(Method = Method,
                cluster1 = cluster1,
                cluster2 = cluster2,
                cluster3 = cluster3)    
df.m <- melt(df, 
             id.vars = c("Method"), 
             measure.vars = c("cluster1", "cluster2","cluster3"),
             variable.name = "Cluster",
             value.name = "Fraction")


plot_radar <- ggplot(data=df.m,  aes(x=Method, y=Fraction, group = Cluster, colour = Cluster )) + 
  annotate("text", x = 1, y = seq(0,0.6,0.1), label = seq(0,0.6,0.1), size =6, col= "black") +
  geom_polygon(size = 3.5, alpha= 0) + 
  ggtitle("Biase et al.")  + 
  scale_y_continuous(labels = NULL) +
  scale_x_discrete() +
  scale_color_manual(values= c("#666666", "#D95F02", "#7570B3"))+
  scale_fill_manual(values= c("#666666", "#D95F02", "#7570B3"))+
  theme_light()+ xlab("Method")+ ylab("Relative Fraction \n of clustering partitions \n")+
  coord_polar(clip = "on", start = -pi/18) +theme(plot.title = element_text(hjust = 0.5)) 

add_theme <- theme(plot.title = element_text(hjust = 0.5, size=30)) + 

  theme(panel.border = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank())

plot_radar + add_theme 

final_plot <- plot_radar + add_theme   + scale_color_manual(name="Cluster", 
                                                            labels = c("1", 
                                                                       "2", 
                                                                       "3"), 
                                                            values = c("#666666", "#D95F02", "#7570B3")) +theme(legend.title = element_text(size=5), legend.text = element_text(size=5))

  
fertig_plot <- final_plot +theme(axis.text=element_text(size=25),
                                  axis.title=element_text(size=25),
                                  legend.title=element_text(size=28),
                                  legend.text =element_text(size=27))+theme(legend.position = "top")  +labs(fill = "Association \n measure:") 

pdf("C:/Users/ri89por/Desktop/Forschung/BIBM/Conference-LaTeX-template_10-17-19/Figuras/biase_radar.pdf")
fertig_plot
dev.off()








#################
### Treutlein ###
#################
output_trans <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/ind_clustering_labeled_neu/5/treutlein.rds")


data_radarchart <- c()
for(l in 4:12){
  #prop.table(table(output_trans[[l]]))
  data_radarchart <- rbind(data_radarchart, prop.table(table(output_trans[[l]])))
}

data_radarchart_prep <- cbind(factor(c(names(output_trans)[4:12]), levels =names(output_trans)[4:12]), data_radarchart[,1:3])
rownames(data_radarchart) <- c(names(output_trans)[4:5], "hclust","k-means", "phate", "t-sne + k-means", "t-sne:hclust", "umap + k-means", "umap:hclust")

treutlein_truth <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/truth_simdata/treutlein.Rds")
ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled_neu/5/treutlein.rds"))
treutlein <- relabel_class(treutlein_truth, ind_cluster_labelled$SC3, return_map = FALSE)
truth <- prop.table(table(as.factor(biase)))




data_radarchart_ready <- data.frame(rbind(truth, data_radarchart[c(1,2,6,5,8),]))



Method <- factor(c(rownames(data_radarchart_ready)),
                 levels = c(rownames(data_radarchart_ready)))
levels(Method)[6]<- "umap + \n k-means" 
levels(Method)[4] <- "t-SNE + \n k-means"
levels(Method)[5] <- "   phate"
levels(Method)[1] <- "truth \n"
cluster1 <- data_radarchart_ready[,1]
cluster2 <- data_radarchart_ready[,2]
cluster3 <- data_radarchart_ready[,3]
cluster4 <- data_radarchart_ready[,4]
cluster5 <- data_radarchart_ready[,5]

# data preparation
df = data.frame(Method = Method,
                cluster1 = cluster1,
                cluster2 = cluster2,
                cluster3 = cluster3, 
                cluster4 = cluster4, 
                cluster5 = cluster5)    
df.m <- melt(df, 
             id.vars = c("Method"), 
             measure.vars = c("cluster1", "cluster2","cluster3", "cluster4","cluster5"),
             variable.name = "Cluster",
             value.name = "Fraction")



plot_radar <- ggplot(data=df.m,  aes(x=Method, y=Fraction, group = Cluster, colour = Cluster )) + 
  annotate("text", x = 1, y = seq(0,0.6,0.1), label = seq(0,0.6,0.1), size =6, col= "black") +
  geom_polygon(size = 3.5, alpha= 0) + 
  ggtitle("Treutlein et al.")  + 
  scale_y_continuous(labels = NULL) +
  scale_x_discrete() +
  scale_color_manual(values= c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77"))+
  scale_fill_manual(values= c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77"))+
  theme_light()+ xlab("Method")+ylab("Relative Fraction \n of clustering partitions \n")+
  coord_polar(clip = "on", start = -pi/18) +theme(plot.title = element_text(hjust = 0.5)) 




final_plot <- plot_radar + add_theme   + scale_color_manual(name="Cluster", 
                                                            labels = c("1", 
                                                                       "2", 
                                                                       "3", 
                                                                       "4", "5"),
                                                            values = c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77"))+
  theme(legend.title = element_text(size=5), legend.text = element_text(size=5))


fertig_plot <- final_plot +theme(axis.text=element_text(size=25),
                                 axis.title=element_text(size=25),
                                 legend.title=element_text(size=28),
                                 legend.text =element_text(size=27))+theme(legend.position = "top")  +labs(fill = "Association \n measure:") 



pdf("C:/Users/ri89por/Desktop/Forschung/BIBM/Conference-LaTeX-template_10-17-19/Figuras/treutlein_radar.pdf")
fertig_plot
dev.off()









gt <- ggplot_gtable(ggplot_build(final_plot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)


#png("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/Paper/Org Template/Figuras/Treutlein.png", width = 480, height = 480)
#grid.draw(gt)
#dev.off()

tf <- file.path("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/Paper/2021//Radarchart_treut_truth.tex")
tikz(tf)
grid.draw(gt)
dev.off()




#################
### Yan ###
#################

output_trans <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/ind_clustering_labeled/7/yan.rds")

data_radarchart <- c()
for(l in 4:12){
  #prop.table(table(output_trans[[l]]))
  data_radarchart <- rbind(data_radarchart, prop.table(table(output_trans[[l]])))
}

data_radarchart_prep <- cbind(factor(c(names(output_trans)[4:12]), levels =names(output_trans)[4:12]), data_radarchart[,1:3])
rownames(data_radarchart) <- c(names(output_trans)[4:5], "hclust","k-means", "phate", "t-sne + k-means", "t-sne:hclust", "umap + k-means", "umap:hclust")


yan_truth <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/truth_simdata/yan.Rds")
ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled/7/yan.rds"))
yan <- relabel_class(yan_truth, ind_cluster_labelled$SC3, return_map = FALSE)
truth <- prop.table(table(as.factor(yan)))





data_radarchart_ready_update <- data.frame(rbind(data_radarchart[c(1,5,7,6,9,8),], data_radarchart[c(3,4,2),], truth))
str(data_radarchart_ready)

data_radarchart_ready <- data_radarchart_ready_update[c(7,8,2,4,6,1,9, 10),]

Method <- factor(c(rownames(data_radarchart_ready)),
                 levels = c(rownames(data_radarchart_ready)))
cluster1 <- data_radarchart_ready[,1]
cluster2 <- data_radarchart_ready[,2]
cluster3 <- data_radarchart_ready[,3]
cluster4 <- data_radarchart_ready[,4]
cluster5 <- data_radarchart_ready[,5]
cluster6 <- data_radarchart_ready[,6]
cluster7 <- data_radarchart_ready[,7]

# data preparation
df = data.frame(Method = Method,
                cluster1 = cluster1,
                cluster2 = cluster2,
                cluster3 = cluster3, 
                cluster4 = cluster4, 
                cluster5 = cluster5, 
                cluster6 = cluster6,
                cluster7 = cluster7) 
df.m <- melt(df, 
             id.vars = c("Method"), 
             measure.vars = c("cluster1", "cluster2","cluster3", "cluster4","cluster5", "cluster6", "cluster7"),
             variable.name = "Cluster",
             value.name = "Fraction")





plot_radar <- ggplot(data=df.m,  aes(x=Method, y=Fraction, group = Cluster, colour = Cluster )) + 
  annotate("text", x = 1, y = seq(0,0.6,0.1), label = seq(0,0.6,0.1), size =2.5, col= "darkgrey") +
  geom_polygon(size = 1.2, alpha= 0) + 
  ggtitle("Yan")  + 
  scale_y_continuous(labels = NULL) +
  scale_x_discrete() +
  #scale_color_manual(values=c("#7FC97F", "#BEAED4", "#FDC086", "#F0027F", "#386CB0", "#FFFF99", "#BF5B17"))+
  #scale_fill_manual(values= c("#7FC97F", "#BEAED4", "#FDC086", "#F0027F", "#386CB0", "#FFFF99", "#BF5B17"))+
  scale_color_manual(values=c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77", "#E7298A", "#66A61E"))+
  scale_fill_manual(values= c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77", "#E7298A", "#66A61E"))+
  theme_light()+ xlab("Method") + ylab("Relative Fraction of clustering partitions")+
  coord_polar(clip = "on", start = -pi/12) +theme(plot.title = element_text(hjust = 0.5)) 



library(RColorBrewer)
mypalette<-brewer.pal(9,"Accent")
mypalette<-brewer.pal(8,"Dark2")
display.brewer.pal(n = 8, name = 'Dark2')
display.brewer.pal(n = 8, name = 'Accent')
mypalette


final_plot <- plot_radar + add_theme   + scale_color_manual(name="k", 
                                                            labels = c("1", 
                                                                       "2", 
                                                                       "3", 
                                                                       "4", "5", "6", "7"),
                                                  values = c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77", "#E7298A", "#66A61E"))+
 theme(legend.title = element_text(size=6), legend.text = element_text(size=7))





#png("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/Paper/Org Template/Figuras/Yan.png", width = 480, height = 480) #700 700
#final_plot
#dev.off()



gt <- ggplot_gtable(ggplot_build(final_plot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

tf <- file.path("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/Paper/2021//Radarchart_yan_truth.tex")
tikz(tf)
grid.draw(gt)
dev.off()





#################
### Deng ###
#################

output_trans <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/ind_clustering_labeled/10/deng.rds")

data_radarchart <- c()
for(l in 4:12){
  #prop.table(table(output_trans[[l]]))
  data_radarchart <- rbind(data_radarchart, prop.table(table(output_trans[[l]])))
}

data_radarchart_prep <- cbind(factor(c(names(output_trans)[4:12]), levels =names(output_trans)[4:12]), data_radarchart[,1:3])
rownames(data_radarchart) <- c(names(output_trans)[4:5], "hclust","k-means", "phate", "t-sne + k-means", "t-sne:hclust", "umap + k-means", "umap:hclust")


deng_truth <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/truth_simdata/deng.Rds")
ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled/10/deng.rds"))
ind_cluster_labelled$SC3 <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/deng_sc3.rds")
ind_cluster_labelled$adaSC3 <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/deng_adasc3.rds")


deng <- relabel_class(deng_truth, ind_cluster_labelled$SC3, return_map = FALSE)
truth <- prop.table(table(as.factor(deng)))





data_radarchart_ready_update <- data.frame(rbind(data_radarchart[c(1,5,7,6,9,8),], data_radarchart[c(3,4,2),], truth))

str(data_radarchart_ready)

data_radarchart_ready <- data_radarchart_ready_update[c(7,8,2,4,6,1,9, 10),]

Method <- factor(c(rownames(data_radarchart_ready)),
                 levels = c(rownames(data_radarchart_ready)))
cluster1 <- data_radarchart_ready[,1]
cluster2 <- data_radarchart_ready[,2]
cluster3 <- data_radarchart_ready[,3]
cluster4 <- data_radarchart_ready[,4]
cluster5 <- data_radarchart_ready[,5]
cluster6 <- data_radarchart_ready[,6]
cluster7 <- data_radarchart_ready[,7]
cluster8 <- data_radarchart_ready[,8]
cluster9 <- data_radarchart_ready[,9]
cluster10 <- data_radarchart_ready[,10]


# data preparation
df = data.frame(Method = Method,
                cluster1 = cluster1,
                cluster2 = cluster2,
                cluster3 = cluster3, 
                cluster4 = cluster4, 
                cluster5 = cluster5, 
                cluster6 = cluster6,
                cluster7 = cluster7, 
                cluster8 = cluster8, 
                cluster9 = cluster9,
                cluster10 = cluster10) 
df.m <- melt(df, 
             id.vars = c("Method"), 
             measure.vars = c("cluster1", "cluster2","cluster3", "cluster4","cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10"),
             variable.name = "Cluster",
             value.name = "Fraction")




plot_radar <- ggplot(data=df.m,  aes(x=Method, y=Fraction, group = Cluster, colour = Cluster )) + 
  annotate("text", x = 1, y = seq(0,0.6,0.1), label = seq(0,0.6,0.1), size =2.5, col= "darkgrey") +
  geom_polygon(size = 1.2, alpha= 0) + 
  ggtitle("Deng")  + 
  scale_y_continuous(labels = NULL) +
  scale_x_discrete() +
  #scale_color_manual(values=c("#7FC97F", "#BEAED4", "#FDC086", "#F0027F", "#386CB0", "#FFFF99", "#BF5B17"))+
  #scale_fill_manual(values= c("#7FC97F", "#BEAED4", "#FDC086", "#F0027F", "#386CB0", "#FFFF99", "#BF5B17"))+
  scale_color_manual(values=c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77", "#E7298A", "#66A61E","#BEAED4","#BF5B17" , "#386CB0"))+
  scale_fill_manual(values= c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77", "#E7298A", "#66A61E","#BEAED4", "#BF5B17" ,"#386CB0"))+
  theme_light()+ xlab("Method")+ ylab("Relative Fraction of clustering partitions")+
  coord_polar(clip = "on", start = -pi/12) +theme(plot.title = element_text(hjust = 0.5)) 





final_plot <- plot_radar + add_theme   + scale_color_manual(name="k", 
                                                            labels = c("1", 
                                                                       "2", 
                                                                       "3", 
                                                                       "4", "5", "6", "7", "8", "9", "10"),
                                                            values = c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77", "#E7298A", "#66A61E","#BEAED4", "#BF5B17" ,"#386CB0"))+
  theme(legend.title = element_text(size=6), legend.text = element_text(size=7))





#png("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/Paper/Org Template/Figuras/Yan.png", width = 480, height = 480) #700 700
#final_plot
#dev.off()



gt <- ggplot_gtable(ggplot_build(final_plot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

tf <- file.path("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/Paper/2021//Radarchart_deng_truth.tex")
tikz(tf)
grid.draw(gt)
dev.off()


#################
### Goolam ###
#################

output_trans <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/ind_clustering_labeled/5/goolam.rds")

data_radarchart <- c()
for(l in 4:12){
  #prop.table(table(output_trans[[l]]))
  data_radarchart <- rbind(data_radarchart, prop.table(table(output_trans[[l]])))
}

data_radarchart_prep <- cbind(factor(c(names(output_trans)[4:12]), levels =names(output_trans)[4:12]), data_radarchart[,1:3])
rownames(data_radarchart) <- c(names(output_trans)[4:5], "hclust","k-means", "phate", "t-sne + k-means", "t-sne:hclust", "umap + k-means", "umap:hclust")


goolam_truth <- readRDS("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/R Code/Results_dep_k/inst/truth_simdata/goolam.Rds")
ind_cluster_labelled <- readRDS(paste0("inst/ind_clustering_labeled/5/goolam.rds"))
goolam <- relabel_class(goolam_truth, ind_cluster_labelled$SC3, return_map = FALSE)
truth <- prop.table(table(as.factor(goolam)))





data_radarchart_ready_update <- data.frame(rbind(data_radarchart[c(1,5,7,6,9,8),], data_radarchart[c(3,4,2),], truth))
str(data_radarchart_ready)

data_radarchart_ready <- data_radarchart_ready_update[c(7,8,2,4,6,1,9,10),]

Method <- factor(c(rownames(data_radarchart_ready)),
                 levels = c(rownames(data_radarchart_ready)))
cluster1 <- data_radarchart_ready[,1]
cluster2 <- data_radarchart_ready[,2]
cluster3 <- data_radarchart_ready[,3]
cluster4 <- data_radarchart_ready[,4]
cluster5 <- data_radarchart_ready[,5]

# data preparation
df = data.frame(Method = Method,
                cluster1 = cluster1,
                cluster2 = cluster2,
                cluster3 = cluster3, 
                cluster4 = cluster4, 
                cluster5 = cluster5)    
df.m <- melt(df, 
             id.vars = c("Method"), 
             measure.vars = c("cluster1", "cluster2","cluster3", "cluster4","cluster5"),
             variable.name = "Cluster",
             value.name = "Fraction")




plot_radar <- ggplot(data=df.m,  aes(x=Method, y=Fraction, group = Cluster, colour = Cluster )) + 
  annotate("text", x = 1, y = seq(0,0.6,0.1), label = seq(0,0.6,0.1), size =2.5, col= "darkgrey") +
  geom_polygon(size = 1.2, alpha= 0) + 
  ggtitle("Goolam")  + 
  scale_y_continuous(labels = NULL) +
  scale_x_discrete() +
  scale_color_manual(values= c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77"))+
  scale_fill_manual(values= c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77"))+
  theme_light()+ xlab("Method")+ ylab("Relative Fraction of clustering partitions")+
  coord_polar(clip = "on", start = -pi/12) +theme(plot.title = element_text(hjust = 0.5)) 


final_plot <- plot_radar + add_theme   + scale_color_manual(name="k", 
                                                            labels = c("1", 
                                                                       "2", 
                                                                       "3", "4", "5"), 
                                                            values = c("#666666", "#D95F02", "#7570B3", "#E6AB02","#1B9E77"))+
  theme(legend.title = element_text(size=6), legend.text = element_text(size=7))



plot_radar
library(grid)
gt <- ggplot_gtable(ggplot_build(final_plot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)


#png("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/Paper/Org Template/Figuras/Treutlein.png", width = 480, height = 480)
#grid.draw(gt)
#dev.off()

tf <- file.path("C:/Users/ri89por/Desktop/Forschung/NIPS 2.0/Paper/2021//Radarchart_goolam_truth.tex")
tikz(tf)
grid.draw(gt)
dev.off()








