library(tidyverse)
library(caret)
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

# Using the euclidean distance and set scale true 
re_dist <- get_dist(x, stand = TRUE, method = "manhattan")

p<-fviz_dist(re_dist, order = TRUE,
          gradient = list(low = "red", mid = "white", high = "blue"))

ggsave(p, filename = "ODM.png")



##### Hierarchical clustering


#This algorithm is used to find the relationship of individual data points and relationships of clusters.

d<-dist(x)
hc<-hclust(d,method="complete")
hc_single <-hclust(d,method = "single")
hc_ward.D2 <-hclust(d,method = "ward.D2")
hc_ward.D<-hclust(d,method = "ward.D")
hc_average <-hclust(d,method = "average")
par(mfrow=c(3,2))
plot(hc_average)
plot(hc_single)
plot(hc_ward.D)
plot(hc_ward.D2)
plot(hc)

#Let's try to identify what is the best level to cut the dendrogram by looking 
#for the highest separation between clusters. 

# complet linkage
m <- length(hc$height)
s <- c()
for (j in m:2) s <- c(s, hc$height[j]/hc$height[j - 1])

k = (2:m)[s == max(s)]
k2 <- (2:m)[s == max(s[-(k - 1)])]
k3 <- (2:m)[s == max(s[-c(k - 1, k2 - 1)])]
k4 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1)])]
k5 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1, k4 - 1)])]

clus<- data.frame( complet=c(k,k2,k3,k4,k5))

row.names(clus) <- paste0('k',1:5)

# single linkage

m <- length(hc_single$height)
s <- c()
for (j in m:2) s <- c(s, hc_single$height[j]/hc_single$height[j - 1])

k = (2:m)[s == max(s)]
k2 <- (2:m)[s == max(s[-(k - 1)])]
k3 <- (2:m)[s == max(s[-c(k - 1, k2 - 1)])]
k4 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1)])]
k5 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1, k4 - 1)])]

clus<- cbind(clus,single=c(k,k2,k3,k4,k5))


# ward.D2 linkage

m <- length(hc_ward.D2$height)
s <- c()
for (j in m:2) s <- c(s, hc_ward.D2$height[j]/hc_ward.D2$height[j - 1])

k = (2:m)[s == max(s)]
k2 <- (2:m)[s == max(s[-(k - 1)])]
k3 <- (2:m)[s == max(s[-c(k - 1, k2 - 1)])]
k4 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1)])]
k5 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1, k4 - 1)])]

clus<- cbind(clus,ward_D2=c(k,k2,k3,k4,k5))


# ward.D linkage 

m <- length(hc_ward.D$height)
s <- c()
for (j in m:2) s <- c(s, hc_ward.D$height[j]/hc_ward.D$height[j - 1])

k = (2:m)[s == max(s)]
k2 <- (2:m)[s == max(s[-(k - 1)])]
k3 <- (2:m)[s == max(s[-c(k - 1, k2 - 1)])]
k4 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1)])]
k5 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1, k4 - 1)])]

clus<- cbind(clus,ward_D=c(k,k2,k3,k4,k5))

# average linkage 

m <- length(hc_average$height)
s <- c()
for (j in m:2) s <- c(s, hc_average$height[j]/hc_average$height[j - 1])

k = (2:m)[s == max(s)]
k2 <- (2:m)[s == max(s[-(k - 1)])]
k3 <- (2:m)[s == max(s[-c(k - 1, k2 - 1)])]
k4 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1)])]
k5 <- (2:m)[s == max(s[-c(k - 1, k2 - 1, k3 - 1, k4 - 1)])]

clus<- cbind(clus,average=c(k,k2,k3,k4,k5))
print(clus)

# Hierarchical clustering on x dataset
par(mfrow=c(3,2))
plot(hc_average)
plot(hc_single)
plot(hc_ward.D)
plot(hc)
#Where to cut the tree ? k=2 from the table clus
par(mfrow=c(1,3))
fviz_dend(hc_ward.D2, k = 2, k_colors = "jco", as.ggplot = TRUE, show_labels = T,
          cex = 0.525)
fviz_dend(hc_ward.D, k = 2, k_colors = "jco", as.ggplot = TRUE, show_labels = T,
          cex = 0.525)
fviz_dend(hc, k = 2, k_colors = "jco", as.ggplot = TRUE, show_labels = T,
          cex = 0.525)

par(mfrow=c(1,1))


## Heatmap
library(pheatmap)


p<-pheatmap(x,show_rownames=F,show_colnames=FALSE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean",annotation_col = as.factor(n))

ggsave(p, filename = "heatmap.png")

'annotation_col<- data.frame(LeukemiaType =n)

row.names(annotation_col) <- colnames(leukemia_big)
'
#annotation_col<-data.matrix(LeukemiaType=n,rownames.force=colnames(leukemia_big))




# K-means on x dataset
# here we start with the number cluster identify on Ordered Dissimilarity Matrix

set.seed(125348897)

kclu <- kmeans(x, 3)
fviz_cluster(list(data = x, cluster = kclu$cluster), ellipse.type = "norm", geom = "point",
             stand = FALSE, palette = "jco", ggtheme = theme_classic())

# number of data points in each cluster
table(kclu$cluster)


type2kclu <- data.frame(
  LeukemiaType =n,
  cluster=kclu$cluster)

table(type2kclu)

confusionMatrix(as.factor(n), as.factor(type2kclu$LeukemiaType))

###k-medoids

kmclu<-cluster::pam(x,k=3) #  cluster using k-medoids

# make a data frame with Leukemia type and cluster id
type2kmclu = data.frame(
  LeukemiaType =n,
  cluster=kmclu$cluster)

table(type2kmclu)
confusionMatrix(as.factor(n), as.factor(type2kclu$LeukemiaType))

# Calculate distances
dists=dist(x)

# calculate MDS
mds=cmdscale(dists)
par(mfrow=c(1,2))

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
Ks=sapply(2:7,
          function(i) 
            summary(silhouette(pam(x,k=i)))$avg.width)
plot(2:7,Ks,xlab="k",ylab="av. silhouette",type="b",
     pch=19)
par(mfrow=c(1,1))
plot(silhouette(pamclu),main=NULL)

par(mfrow=c(1,1))

#it seems the best value for  k is 2

set.seed(101)
# define the clustering function
pam1 <- function(x,k) 
  list(cluster = pam(x,k, cluster.only=TRUE))

# calculate the gap statistic
pam.gap= clusGap(x, FUN = pam1, K.max = 8,B=50)
par(mfrow=c(1,1))
# plot the gap statistic across k values
p<-plot(pam.gap, main = "Gap statistic for the 'Leukemia' data")
ggsave(p, filename = "gap.png")

#Elbow
fviz_nbclust(x, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")

### k-means with 2 cluster

# K-means on x dataset
set.seed(125348897)
kclu <- kmeans(x, 2)
fviz_cluster(list(data = x, cluster = kclu$cluster), ellipse.type = "norm", geom = "point",
             stand = FALSE, palette = "jco", ggtheme = theme_classic())


# number of data points in each cluster
type2kclu = data.frame(
  LeukemiaType =n,
  cluster=kclu$cluster)

table(type2kclu)


###k-medoids

kmclu<-cluster::pam(x,k=2) #  cluster using k-medoids

# make a data frame with Leukemia type and cluster id
type2kmclu = data.frame(
  LeukemiaType =n,
  cluster=kmclu$cluster)

table(type2kmclu)



#### Dimentional reduction

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

fviz_pca_ind(pcdemo,
             col.ind = as.factor(n), # colorer par groupes
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Ellipse de concentration
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)



par(mfrow=c(1,2))
d=svd(scale(leukemia_big)) # apply SVD
assays=t(d$u) %*% scale(leukemia_big) # projection on eigenarray
plot(assays[1,],assays[2,],pch=19,
     col=as.factor(n))

## kmeans on principal component 
set.seed(1258)
kc <- kmeans(pc, 2)
plot(pc,col=factor(kc$cluster))

fviz_cluster(list(data = pc, cluster = kc$cluster), ellipse.type = "norm", geom = "point",
             stand = FALSE, palette = "jco", ggtheme = theme_classic())

# number of data points in each cluster
type2kclu = data.frame(
  LeukemiaType =n,
  cluster=kc$cluster)
type2kclu
table(type2kclu)
confusionMatrix(as.factor(n), as.factor(type2kclu$LeukemiaType))

