library(tidyverse)
leukemia_big <- read.csv("http://hastie.su.domains/CASI_files/DATA/leukemia_big.csv")

n<-names(leukemia_big)
nrow(leukemia_big)
n
# number of patient
p<-ncol(leukemia_big)
#number of missing value
sum(is.na(leukemia_big))

# Let's extract the label of each patient

n<-gsub('[[:digit:]]+', '', n)
n<-gsub('[[:punct:] ]+',"",n)

# Let's anotate each patient by number 
patient <- paste0("patient", 1:p)
# rename the colunm with the patient number. 
names(leukemia_big)<-patient


x<-t(leukemia_big)

####Ordered Dissimilarity Matrix
#A visual method for assessing clustering tendency is the Ordered Dissimilarity Matrix.

library(factoextra)

# Using the euclidean distance and set scal 
res.dist <- get_dist(x, stand = TRUE, method = "manhattan")

fviz_dist(res.dist, order = TRUE,
          gradient = list(low = "red", mid = "white", high = "blue"))




##### Hiearchical clustering


#  Using this algorithm you can see the relationship of individual data points and relationships of clusters.

d=dist(x)
hc=hclust(d,method="complete")
plot(hc)
# Hierarchical clustering on x dataset
fviz_dend(hclust(dist(x)), k = 3, k_colors = "jco", as.ggplot = TRUE, show_labels = T,
          cex = 0.525)


library(pheatmap)


pheatmap(x,show_rownames=T,show_colnames=FALSE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")

annotation_col<- data.frame(LeukemiaType =n)

row.names(annotation_col) <- colnames(leukemia_big)

#annotation_col<-data.matrix(LeukemiaType=n,rownames.force=colnames(leukemia_big))



pheatmap(x,show_rownames=F,show_colnames=FALSE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")


#Where to cut the tree ?

set.seed(123567)
# Hierarchical clustering on x dataset
fviz_dend(hclust(dist(x)), k = 3, k_colors = "jco", as.ggplot = TRUE, show_labels = T,
                 cex = 0.525)
# K-means on x dataset
kclu <- kmeans(x, 3)
fviz_cluster(list(data = x, cluster = kclu$cluster), ellipse.type = "norm", geom = "point",
             stand = FALSE, palette = "jco", ggtheme = theme_classic())

# number of data points in each cluster
table(kclu$cluster)


type2kclu = data.frame(
  LeukemiaType =n,
  cluster=kclu$cluster)

table(type2kclu)

#k-medoids

kmclu=cluster::pam(x,k=3) #  cluster using k-medoids

# make a data frame with Leukemia type and cluster id
type2kmclu = data.frame(
  LeukemiaType =n,
  cluster=kmclu$cluster)

table(type2kmclu)


# Calculate distances
dists=dist(x)

# calculate MDS
mds=cmdscale(dists)

# plot the patients in the 2D space
plot(mds,pch=19,col=rainbow(5)[kclu$cluster])

# set the legend for cluster colors
legend("bottomright",
       legend=paste("clu",unique(kclu$cluster)),
       fill=rainbow(5)[unique(kclu$cluster)],
       border=NA,box.col=NA)


####How to choose "k", the number of clusters


library(cluster)
set.seed(101)
pamclu=cluster::pam(x,k=3)
plot(silhouette(pamclu),main=NULL)


Ks=sapply(2:7,
          function(i) 
            summary(silhouette(pam(x,k=i)))$avg.width)
plot(2:7,Ks,xlab="k",ylab="av. silhouette",type="b",
     pch=19)


#t seems the best value for  k is 2
fviz_dend(hclust(dist(x)), k = 3, k_colors = "jco", as.ggplot = TRUE, show_labels = T,
          cex = 0.525)

library(cluster)
set.seed(101)
# define the clustering function
pam1 <- function(x,k) 
  list(cluster = pam(x,k, cluster.only=TRUE))

# calculate the gap statistic
pam.gap= clusGap(x, FUN = pam1, K.max = 8,B=50)

# plot the gap statistic accross k values
plot(pam.gap, main = "Gap statistic for the 'Leukemia' data")

### kmeans with 2 cluster

# K-means on x dataset
kclu <- kmeans(x, 2)
fviz_cluster(list(data = x, cluster = kclu$cluster), ellipse.type = "norm", geom = "point",
             stand = FALSE, palette = "jco", ggtheme = theme_classic())

# number of data points in each cluster
type2kclu = data.frame(
  LeukemiaType =n,
  cluster=kclu$cluster)

table(type2kclu)


# Dimentional reduction

library(TeachingDemos)
library(KernSmooth)
library(MASS)
library(lattice)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(factoextra)
library(cluster)
library(fpc)
library(NbClust)
pcdemo <- prcomp(x)
pc <- pcdemo$x[, 1:2]
plot(pc)
fviz_eig(pcdemo)
est <- bkde2D(pc, bandwidth = c(0.5, 0.5), gridsize = c(20, 20))


pc %>%
  as.data.frame -> pc
p1 <- ggplot(pc, aes(x = PC1, y = PC2)) + stat_density2d(aes(fill = ..density..),
                                                         geom = "tile", contour = F)
p2 <- p1 + scale_fill_distiller(palette = "RdGy")
p3 <- p1 + scale_fill_viridis(option = "inferno")
ggpubr::ggarrange(p1, p2, p3, ncol = 3)


fviz_pca_ind(pcdemo,
             col.ind = "cos2", # Colorer par le cos2
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     
)
annotation_col$LeukemiaType<-as.factor(annotation_col$LeukemiaType)
fviz_pca_ind(pcdemo,
             col.ind = as.factor(n), # colorer par groupes
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Ellipse de concentration
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)


# calculate the PCA only for our genes and all the samples
fviz_pca_biplot(pcdemo, repel = TRUE,
                col.var = "#2E9FDF", 
                col.ind = "#696969"  
)

par(mfrow=c(1,2))
d=svd(scale(leukemia_big)) # apply SVD
assays=t(d$u) %*% scale(leukemia_big) # projection on eigenassays
plot(assays[1,],assays[2,],pch=19,
     col=as.factor(annotation_col$LeukemiaType))
#plot(d$v[,1],d$v[,2],pch=19,
#     col=annotation_col$LeukemiaType)
pr=prcomp(x,center=TRUE,scale=TRUE) # apply PCA on transposed matrix

# plot new coordinates from PCA, projections on eigenvectors
# since the matrix is transposed eigenvectors represent 
plot(pr$x[,1],pr$x[,2],col=as.factor(annotation_col$LeukemiaType))
