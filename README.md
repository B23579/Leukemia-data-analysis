# Exploratory data analysis of the Leukemia gene expression

## Problem statement 
Researchers in genomics are frequently interested in understanding how their samples are related to one another and how different variables in the data set are correlated. They usually want to answer the following questions: Are our replicates similar to each other? Do the samples from the same treatment group have similar genome-wide signals? Do the patients with similar diseases have similar gene expression profiles?. From these informations, one can develop an hypothesis, identify outliers in the sample, or determine which sample groups require more data points. In this perspective, the project goal is to do an exploratory data analysis of the Leukemia gene expression to answer some scientists' questions. 


## Algorithms use to Approch the problem 
The clustering techniques(k-means, k-medoids, and hierarchical clustering) will be applied to find groups of patients that are more similar to each other than the rest of the patients using the distance or similarity metric between patients’ expression profiles. 

The dimension reduction technique (PCA, SVD) will be applied to  reduce the number of genes expression or variables to obtain  a few principal highly variables that capture most of the variations inside the data instead of thousands. This is useful for clustering, predictive modeling and 2D or 3D visualization of many genes expression.
## PCA
![fg](./PCA.png) 
