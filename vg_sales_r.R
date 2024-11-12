## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
library(readxl)
vg_sales <- read_excel("vg_sales.xlt", 
    col_types = c("numeric", "text", "text", 
        "numeric", "text", "text", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric"), n_max = 200)

## ------------------------------------------------------------------------------------------------------------------------------
str(vg_sales)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
for (i in 1:nrow(vg_sales)){
  rownames(vg_sales)[i]<- vg_sales[i,1]
}

vg_sales <- vg_sales[,-1]


## ------------------------------------------------------------------------------------------------------------------------------
vg_sales$Platform = factor(vg_sales$Platform)
vg_sales$Genre = factor(vg_sales$Genre)
vg_sales$Publisher = factor(vg_sales$Publisher)



## ------------------------------------------------------------------------------------------------------------------------------
any(is.na(vg_sales))


## ------------------------------------------------------------------------------------------------------------------------------
str(vg_sales)


## ------------------------------------------------------------------------------------------------------------------------------
vg_sales


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$Name)


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$Platform)


## ------------------------------------------------------------------------------------------------------------------------------
levels(vg_sales$Platform)


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$Year)


## ------------------------------------------------------------------------------------------------------------------------------
hist(vg_sales$Year, breaks = 30)


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$Genre)


## ------------------------------------------------------------------------------------------------------------------------------
levels(vg_sales$Genre)


## ------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
ggplot(subset(vg_sales), aes(y="", fill=Genre)) +
geom_bar(aes(x=..count../sum(..count..))) +
scale_x_continuous(labels=scales::percent) +
ylab('Percentage') +
xlab('Genre') +
coord_polar()


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$Publisher)


## ------------------------------------------------------------------------------------------------------------------------------
levels(vg_sales$Publisher)


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$NA_Sales)


## ------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
boxplot(vg_sales$NA_Sales, main="NA_Sales", col="blue")


## ------------------------------------------------------------------------------------------------------------------------------
library(gamlss)
fit.EXP <- histDist(vg_sales$NA_Sales, family=EXP, nbins = 30, main="Exponential distribution")


## ------------------------------------------------------------------------------------------------------------------------------
fit.GA <- histDist(vg_sales$NA_Sales, family=GA, nbins=30, main = "Gamma Distribution")


## ------------------------------------------------------------------------------------------------------------------------------
fit.IG <- histDist(vg_sales$NA_Sales, family=IG, nbins=30, main = "Inverse Gamma Distribution")


## ------------------------------------------------------------------------------------------------------------------------------
fit.LOGNO <- histDist(vg_sales$NA_Sales, family=LOGNO, nbins=30, main = "LogNormal Distribution")


## ------------------------------------------------------------------------------------------------------------------------------
fit.WEI <- histDist(vg_sales$NA_Sales,family=WEI, nbins=30, main = "Weibull Distribution")


## ------------------------------------------------------------------------------------------------------------------------------
data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"),
AIC=c(AIC(fit.EXP), AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO), AIC(fit.WEI)),
SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc, fit.WEI$sbc))



## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$EU_Sales)


## ------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
boxplot(vg_sales$EU_Sales, main="EU_Sales", col="#7FFFD4")


## ------------------------------------------------------------------------------------------------------------------------------
library(gamlss)
fit.EXP1 <- histDist(vg_sales$EU_Sales, family=EXP, nbins = 30, main="Exponential distribution")


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$JP_Sales)


## ------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
boxplot(vg_sales$JP_Sales, main="JP_Sales", col="red")


## ------------------------------------------------------------------------------------------------------------------------------
library(gamlss)
fit.EXP2 <- histDist(vg_sales$JP_Sales, family=EXP, nbins = 30, main="Exponential distribution")


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$Other_Sales)


## ------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
boxplot(vg_sales$Other_Sales, main="Other_Sales", col="#FFD700")


## ------------------------------------------------------------------------------------------------------------------------------
library(gamlss)
fit.EXP3 <- histDist(vg_sales$JP_Sales, family=EXP, nbins = 30, main="Exponential distribution")


## ------------------------------------------------------------------------------------------------------------------------------
summary(vg_sales$Global_Sales)


## ------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
boxplot(vg_sales$Global_Sales, main="Global_Sales", col="#D8BFD8")


## ------------------------------------------------------------------------------------------------------------------------------
library(gamlss)
fit.EXP4 <- histDist(vg_sales$Global_Sales, family=EXP, nbins = 30, main="Exponential distribution")


## ------------------------------------------------------------------------------------------------------------------------------
fit.GA4 <- histDist(vg_sales$Global_Sales, family=GA, nbins=30, main = "Gamma Distribution")


## ------------------------------------------------------------------------------------------------------------------------------
fit.IG4 <- histDist(vg_sales$Global_Sales, family=IG, nbins=30, main = "Inverse Gamma Distribution")


## ------------------------------------------------------------------------------------------------------------------------------
fit.LOGNO4 <- histDist(vg_sales$Global_Sales, family=LOGNO, nbins=30, main = "LogNormal Distribution")


## ------------------------------------------------------------------------------------------------------------------------------
fit.WEI4 <- histDist(vg_sales$Global_Sales,family=WEI, nbins=30, main = "Weibull Distribution")

## ------------------------------------------------------------------------------------------------------------------------------
data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"),
AIC=c(AIC(fit.EXP4), AIC(fit.GA4), AIC(fit.IG4), AIC(fit.LOGNO4), AIC(fit.WEI4)),
SBC=c(fit.EXP4$sbc, fit.GA4$sbc, fit.IG4$sbc, fit.LOGNO4$sbc, fit.WEI4$sbc))



## ------------------------------------------------------------------------------------------------------------------------------
vg_sales.pca <- vg_sales[1:200, 6:10]


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
for (i in 1:nrow(vg_sales)){
  rownames(vg_sales.pca)[i]<- paste(vg_sales[i,1], levels(vg_sales$Platform)[as.numeric(vg_sales[i,2])], sep = " - ") 
}



## ------------------------------------------------------------------------------------------------------------------------------
head(vg_sales.pca)


## ------------------------------------------------------------------------------------------------------------------------------
apply (vg_sales.pca, 2, var)


## ------------------------------------------------------------------------------------------------------------------------------
scaled_vg_sales.pca <- apply(vg_sales.pca, 2, scale)
head(scaled_vg_sales.pca)


## ------------------------------------------------------------------------------------------------------------------------------
std_sales.cov <- cov(scaled_vg_sales.pca) 
sales_eigen <- eigen(std_sales.cov)
str(sales_eigen)


## ------------------------------------------------------------------------------------------------------------------------------
(eg <- sales_eigen$vectors[,1:2])


## ------------------------------------------------------------------------------------------------------------------------------
eg <- -eg
row.names(eg) <- c("NA_Sales", "EU_Sales", "JP_Sales", "Other_Sales", "Global_Sales")
colnames(eg) <- c("PC1", "PC2")
eg


## ------------------------------------------------------------------------------------------------------------------------------
PC1 <- scaled_vg_sales.pca %*%  eg[,1]
PC2 <- scaled_vg_sales.pca %*%  eg[,2]

PC <- data.frame(Video_Games = row.names(vg_sales.pca), PC1, PC2)
head(PC, digits = 3)


## ------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(modelr)
ggplot(PC, aes(PC1, PC2)) +
modelr::geom_ref_line(h = 0) +
modelr::geom_ref_line(v = 0) +
geom_text(aes(label = Video_Games), size = 3) +
xlab("First Principal Component") +
ylab("Second Principal Component") +
ggtitle("Scores-1PC and 2PC")


## ------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
res.pca <- prcomp(scaled_vg_sales.pca, center = TRUE, scale. = FALSE)
biplot(res.pca, cex.axis = 0.5, scale=1)
abline(h=0)
abline(v=0)


## ------------------------------------------------------------------------------------------------------------------------------
PVE <- sales_eigen$values/sum(sales_eigen$values)
round(PVE,3)


## ------------------------------------------------------------------------------------------------------------------------------
round(sales_eigen$values,3) 


## ------------------------------------------------------------------------------------------------------------------------------
PVEplot <- qplot(c(1:5), PVE) + 
  geom_line()+
  xlab("Principal Component")+
  ylab("PVE")+
  ggtitle ("Scree Plot")+
  ylim(0,1)

## ------------------------------------------------------------------------------------------------------------------------------
cumPVE <- qplot(c(1:5), cumsum(PVE)) + 
  geom_line()+
  xlab("Principal Component") + 
  ylab(NULL)+
  ggtitle("Cumulative Scree Plot") +
  ylim(0,1)

## ------------------------------------------------------------------------------------------------------------------------------
library(gridExtra)
grid.arrange(PVEplot, cumPVE, ncol=2)


## ------------------------------------------------------------------------------------------------------------------------------
vg_sales.cl <- vg_sales[1:200, 6:10]
vg_sales.scaled <- apply(vg_sales.cl, 2, scale)
head(vg_sales.scaled)


## ------------------------------------------------------------------------------------------------------------------------------
dist.eucl <- dist(vg_sales.scaled, method = "euclidean")
round(as.matrix(dist.eucl)[1:10, 1:10], 2)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
dist.eucl.200 <- dist(vg_sales.scaled, method = "euclidean")


## ------------------------------------------------------------------------------------------------------------------------------
dist.man <- dist(vg_sales.scaled, method = "manhattan")
round(as.matrix(dist.man)[1:10, 1:10], 2)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
dist.man.200 <- dist(vg_sales.scaled, method = "manhattan")


## ------------------------------------------------------------------------------------------------------------------------------
library(cluster)

gower.dist <- daisy(vg_sales [,-1])
round(as.matrix(gower.dist)[1:10,1:10],2)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
gower.dist.200 <- daisy(vg_sales [,-1])


## ------------------------------------------------------------------------------------------------------------------------------
library(factoextra)
dist.eucl.1 <- as.matrix(dist.eucl)
dist.eucl.1 <- dist.eucl.1[1:10,1:10]
dist.eucl.1 <- dist(dist.eucl.1, method = "euclidean")
fviz_dist(dist.eucl.1)


## ------------------------------------------------------------------------------------------------------------------------------
library(factoextra)
dist.eucl.2 <- as.matrix(dist.eucl.200)
dist.eucl.2 <- dist.eucl.2
dist.eucl.2 <- dist(dist.eucl.2, method = "euclidean")
fviz_dist(dist.eucl.2)


## ------------------------------------------------------------------------------------------------------------------------------
library(factoextra)
dist.man.1 <- as.matrix(dist.man)
dist.man.1 <- dist.man.1[1:10,1:10]
dist.man.1 <- dist(dist.man.1, method = "manhattan")
fviz_dist(dist.man.1)


## ------------------------------------------------------------------------------------------------------------------------------
library(factoextra)
dist.man.2 <- as.matrix(dist.man)
dist.man.2 <- dist.man.2
dist.man.2 <- dist(dist.man.2, method = "manhattan")
fviz_dist(dist.man.2)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
library(NbClust)
res.alm <- NbClust(vg_sales.scaled, distance = "euclidean", min.nc = 2, max.nc = 10,
method = "average")


## ------------------------------------------------------------------------------------------------------------------------------
fviz_nbclust(res.alm) 



## ------------------------------------------------------------------------------------------------------------------------------
res.hc2 <- hclust(dist.eucl.200, method = "average")
fviz_dend(res.hc2, cex = 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
res.hc2 <- hclust(dist.eucl.200, method = "average")
cor(dist.eucl.200, cophenetic(res.hc2))


## ------------------------------------------------------------------------------------------------------------------------------
d.alm <- cutree(res.hc2, k = 2)
fviz_dend(res.hc2, k = 2, cex = 0.5, k_colors = c("#00AFBB", "#E7B800"), color_labels_by_k = TRUE, rect = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
table(d.alm)


## ------------------------------------------------------------------------------------------------------------------------------
pairs(vg_sales.scaled, gap=0, pch=d.alm, col=c("#00AFBB", "#E7B800"),  main="Original space\n- Average Linkage Method and Euclidean Distance, K=2"[d.alm])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data = vg_sales.scaled, cluster = d.alm), palette = c("#00AFBB", "#E7B800"), ellipse.type = "convex", main="PCs space", repel = TRUE,
show.clust.aver = FALSE, ggtheme = theme_minimal())+labs(
subtitle = "Original space\n- Average Linkage Method and Euclidean Distance, K=2", cex.sub= 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster =  d.alm),geom="point",ellipse.type="convex",palette = c("#00AFBB", "#E7B800"),repel = TRUE, show.clust.cent = FALSE)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
library(NbClust)
res.alm.md <- NbClust(vg_sales.scaled, distance = "manhattan", min.nc = 2, max.nc = 10,
method = "average")


## ------------------------------------------------------------------------------------------------------------------------------
fviz_nbclust(res.alm.md) 


## ------------------------------------------------------------------------------------------------------------------------------
res.alm.md1 <- hclust(dist.man.200, method = "average")
fviz_dend(res.alm.md1, cex = 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
res.cd.alm.md <- hclust(dist.man.200, method = "average")
cor(dist.man.200, cophenetic(res.cd.alm.md ))


## ------------------------------------------------------------------------------------------------------------------------------
alm.md <- cutree(res.cd.alm.md, k = 2)
fviz_dend(res.cd.alm.md, k = 2, cex = 0.5, k_colors = c("#00AFBB", "#E7B800"), color_labels_by_k = TRUE, rect = TRUE)

## ------------------------------------------------------------------------------------------------------------------------------
table(alm.md)


## ------------------------------------------------------------------------------------------------------------------------------
pairs(vg_sales.scaled, gap=0, pch=alm.md, col=c("#00AFBB", "#E7B800"),  main="Original space\n- Average Linkage Method and Manhattan Distance, K=2"[alm.md])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data = vg_sales.scaled, cluster = alm.md), palette = c("#00AFBB", "#E7B800"), ellipse.type = "convex", main="PCs space", repel = TRUE,
show.clust.aver = FALSE, ggtheme = theme_minimal())+labs(
subtitle = "Original space\n- Average Linkage Method and Manhattan Distance, K=2", cex.sub= 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster =  alm.md),geom="point",ellipse.type="convex",palette = c("#00AFBB", "#E7B800"),repel = TRUE, show.clust.cent = FALSE)


## ------------------------------------------------------------------------------------------------------------------------------
library(NbClust)
res.wlm1 <- NbClust(vg_sales.scaled, distance = "euclidean", min.nc = 2, max.nc = 10,
method = "ward.D2")


## ------------------------------------------------------------------------------------------------------------------------------
fviz_nbclust(res.wlm1) 


## ------------------------------------------------------------------------------------------------------------------------------
res.wlm <- hclust(d = dist.eucl.200, method = "ward.D2")
fviz_dend(res.wlm, cex = 0.5)



## ------------------------------------------------------------------------------------------------------------------------------
res.coph1 <- cophenetic(res.wlm)
cor(dist.eucl.200, res.coph1) 


## ------------------------------------------------------------------------------------------------------------------------------
d.wlm <- cutree(res.wlm, k = 4)
fviz_dend(res.wlm, k = 4, cex = 0.5, k_colors = c("#2E9FDF", "#FC4E07", "#CF3E4B", "#287233"), color_labels_by_k = TRUE, rect = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
table(d.wlm)


## ------------------------------------------------------------------------------------------------------------------------------
pairs(vg_sales.scaled, gap=0, pch=d.wlm, col=c("#2E9FDF", "#FC4E07", "#CF3E4B", "#287233"),  main="Original space\n- Ward's Linkage Method and Euclidean Distance, K=4"[d.wlm])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data = vg_sales.scaled, cluster = d.wlm), palette = c("#2E9FDF", "#FC4E07", "#CF3E4B", "#287233"), ellipse.type = "convex", main="PCs space", repel = TRUE,
show.clust.aver = FALSE, ggtheme = theme_minimal())+labs(
subtitle = "Original space\n- Ward's Linkage Method and Euclidean Distance, K=4", cex.sub= 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster =  d.wlm),geom="point",ellipse.type="convex",palette = c("#2E9FDF", "#FC4E07", "#CF3E4B", "#287233"),repel = TRUE, show.clust.cent = FALSE)


## ------------------------------------------------------------------------------------------------------------------------------
library(NbClust)
res.wlm.md <- NbClust(vg_sales.scaled, distance = "manhattan", min.nc = 2, max.nc = 10,
method = "ward.D2")


## ------------------------------------------------------------------------------------------------------------------------------
fviz_nbclust(res.wlm.md) 


## ------------------------------------------------------------------------------------------------------------------------------
res.wlm.md1 <- hclust(d = dist.man.200, method = "ward.D2")
fviz_dend(res.wlm.md1, cex = 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
res.coph2 <- cophenetic(res.wlm.md1)
cor(dist.man.200, res.coph2) 


## ------------------------------------------------------------------------------------------------------------------------------
d.wlm.md1 <- cutree(res.wlm.md1, k = 3)
fviz_dend(res.wlm.md1, k = 3, cex = 0.5, k_colors = c("#00AFBB", "#E7B800", "#FC4E07"), color_labels_by_k = TRUE, rect = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
table(d.wlm.md1)


## ------------------------------------------------------------------------------------------------------------------------------
pairs(vg_sales.scaled, gap=0, pch=d.wlm.md1, col=c("#00AFBB", "#E7B800", "#FC4E07"),  main="Original space\n- Ward's Linkage Method and the Manhattan Distance, K=3"[d.wlm.md1])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data = vg_sales.scaled, cluster = d.wlm.md1), palette = c("#00AFBB", "#E7B800", "#FC4E07"), ellipse.type = "convex", main="PCs space", repel = TRUE,
show.clust.aver = FALSE, ggtheme = theme_minimal())+labs(
subtitle = "Original space\n- Ward's Linkage Method and the Manhattan Distance, K=3", cex.sub= 0.5)

## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster =  d.wlm.md1),geom="point",ellipse.type="convex",palette = c("#2E9FDF", "#FC4E07", "#CF3E4B"),repel = TRUE, show.clust.cent = FALSE)


## ------------------------------------------------------------------------------------------------------------------------------
library(NbClust)
res.slm.ed <- NbClust(vg_sales.scaled, distance = "euclidean", min.nc = 2, max.nc = 10,
method = "single")


## ------------------------------------------------------------------------------------------------------------------------------
fviz_nbclust(res.slm.ed) 


## ------------------------------------------------------------------------------------------------------------------------------
res.slm.ed1 <- hclust(d = dist.eucl.200, method = "single")
fviz_dend(res.slm.ed1, cex = 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
res.coph.sm.ed <- cophenetic(res.slm.ed1)
cor(dist.eucl.200, res.coph.sm.ed) 


## ------------------------------------------------------------------------------------------------------------------------------
res.slm.ed2 <- cutree(res.slm.ed1, k = 2)
fviz_dend(res.slm.ed1, k = 3, cex = 0.5, k_colors = c("#FC4E07", "#E7B800"), color_labels_by_k = TRUE, rect = TRUE)



## ------------------------------------------------------------------------------------------------------------------------------

table(res.slm.ed2)



## ------------------------------------------------------------------------------------------------------------------------------
pairs(vg_sales.scaled, gap=0, pch=res.slm.ed2, col=c("#00AFBB", "#E7B800"),  main="Original space\n- Single Linkage Method and the Euclidean Distance, K=2"[res.slm.ed2])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data = vg_sales.scaled, cluster = res.slm.ed2), palette = c("#00AFBB", "#E7B800"), ellipse.type = "convex", main="PCs space", repel = TRUE,
show.clust.aver = FALSE, ggtheme = theme_minimal())+labs(
subtitle = "Original space\n- Single Linkage Method and the Euclidean Distance, K=2", cex.sub= 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster = res.slm.ed2),geom="point",ellipse.type="convex",palette = c("#2E9FDF", "#FC4E07"),repel = TRUE, show.clust.cent = FALSE)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
library(NbClust)
slm.md <- NbClust(vg_sales.scaled, distance = "manhattan", min.nc = 2, max.nc = 10,
method = "single")


## ------------------------------------------------------------------------------------------------------------------------------
fviz_nbclust(slm.md) 


## ------------------------------------------------------------------------------------------------------------------------------
slm.md1 <- hclust(dist.man.200, method = "single")
fviz_dend(slm.md1, cex = 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
slm.md2 <- hclust(dist.man.200, method = "single")
cor(dist.man.200, cophenetic(slm.md2))


## ------------------------------------------------------------------------------------------------------------------------------
slm.md3 <- cutree(slm.md2, k = 2)
fviz_dend(slm.md2, k = 2, cex = 0.5, k_colors = c("#00AFBB", "#E7B800"), color_labels_by_k = TRUE, rect = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
table(slm.md3)


## ------------------------------------------------------------------------------------------------------------------------------
pairs(vg_sales.scaled, gap=0, pch=slm.md3, col=c("#00AFBB", "#E7B800"),  main="Original space\n- Single Linkage Method and Manhattan Distance, K=2"[slm.md3])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data = vg_sales.scaled, cluster = slm.md3), palette = c("#00AFBB", "#E7B800"), ellipse.type = "convex", main="PCs space", repel = TRUE,
show.clust.aver = FALSE, ggtheme = theme_minimal())+labs(
subtitle = "Original space\n- Single Linkage Method and Manhattan Distance, K=2", cex.sub= 0.5)

## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster = slm.md3),geom="point",ellipse.type="convex",palette = c("#2E9FDF", "#FC4E07"),repel = TRUE, show.clust.cent = FALSE)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
library(NbClust)
clm.ed <- NbClust(vg_sales.scaled, distance = "euclidean", min.nc = 2, max.nc = 10,
method = "complete")


## ------------------------------------------------------------------------------------------------------------------------------
fviz_nbclust(clm.ed) 


## ------------------------------------------------------------------------------------------------------------------------------
clm.ed1 <- hclust(dist.eucl.200, method = "complete")
fviz_dend(clm.ed1, cex = 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
clm.ed2 <- hclust(dist.eucl.200, method = "complete")
cor(dist.eucl.200, cophenetic(clm.ed2))


## ------------------------------------------------------------------------------------------------------------------------------
clm.ed3 <- cutree(clm.ed2, k = 2)
fviz_dend(clm.ed2, k = 2, cex = 0.5, k_colors = c("#00AFBB", "#E7B800"), color_labels_by_k = TRUE, rect = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
table(clm.ed3)


## ------------------------------------------------------------------------------------------------------------------------------
pairs(vg_sales.scaled, gap=0, pch=clm.ed3, col=c("#00AFBB", "#E7B800"),  main="Original space\n- Complete Linkage Method and Euclidean Distance, K=2"[clm.ed3])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data = vg_sales.scaled, cluster = clm.ed3), palette = c("#00AFBB", "#E7B800"), ellipse.type = "convex", main="PCs space", repel = TRUE,
show.clust.aver = FALSE, ggtheme = theme_minimal())+labs(
subtitle = "Original space\n- Complete Linkage Method and Euclidean Distance, K=2", cex.sub= 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster = clm.ed3),geom="point",ellipse.type="convex",palette = c("#2E9FDF", "#FC4E07"),repel = TRUE, show.clust.cent = FALSE)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
library(NbClust)
clm.md <- NbClust(vg_sales.scaled, distance = "manhattan", min.nc = 2, max.nc = 10,
method = "complete")


## ------------------------------------------------------------------------------------------------------------------------------
fviz_nbclust(clm.md) 


## ------------------------------------------------------------------------------------------------------------------------------
clm.md1 <- hclust(dist.man.200, method = "complete")
fviz_dend(clm.md1, cex = 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
clm.md2 <- hclust(dist.man.200, method = "complete")
cor(dist.man.200, cophenetic(clm.md2))


## ------------------------------------------------------------------------------------------------------------------------------
clm.md3 <- cutree(clm.md2, k = 2)
fviz_dend(clm.md2, k = 2, cex = 0.5, k_colors = c("#00AFBB", "#E7B800"), color_labels_by_k = TRUE, rect = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
table(clm.md3)


## ------------------------------------------------------------------------------------------------------------------------------
pairs(vg_sales.scaled, gap=0, pch=clm.md3, col=c("#00AFBB", "#E7B800"),  main="Original space\n- Complete Linkage Method and Euclidean Distance, K=2"[clm.md3])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data = vg_sales.scaled, cluster = clm.md3), palette = c("#00AFBB", "#E7B800"), ellipse.type = "convex", main="PCs space", repel = TRUE,
show.clust.aver = FALSE, ggtheme = theme_minimal())+labs(
subtitle = "Original space\n- Complete Linkage Method and Euclidean Distance, K=2", cex.sub= 0.5)


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster = clm.md3),geom="point",ellipse.type="convex",palette = c("#2E9FDF", "#FC4E07"),repel = TRUE, show.clust.cent = FALSE)


## ------------------------------------------------------------------------------------------------------------------------------
library(factoextra)
fviz_nbclust(vg_sales.scaled, kmeans, nstart = 25, method = "wss")+
   geom_vline( xintercept = 5, linetype = 2)


## ------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
km.res <- kmeans(vg_sales.scaled, 5, nstart = 25)
print(km.res)


## ------------------------------------------------------------------------------------------------------------------------------
aggregate(vg_sales.scaled, by=list(cluster = km.res$cluster), mean)


## ------------------------------------------------------------------------------------------------------------------------------
cl.km <- km.res$cluster
pairs(vg_sales.scaled, gap=0, pch=cl.km, col = c("#2E9FDF", "#FC4E07", "#CF3E4B", "#287233", "#E7B800")[cl.km])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(km.res,
             data = vg_sales.scaled,
             palette = c("#2E9FDF", "#FC4E07", "#CF3E4B", "#287233", "#E7B800"),
             ellipse.type = "euclid",
             star.plot = TRUE,
             repel = TRUE,
             ggtheme = theme_minimal())


## ------------------------------------------------------------------------------------------------------------------------------
library(factoextra)
fviz_nbclust(vg_sales.scaled, kmeans, nstart = 25, method = "wss")+
   geom_vline( xintercept = 4, linetype = 2)

## ------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
km.res1 <- kmeans(vg_sales[,6:10], 4, nstart = 25)
print(km.res1)


## ------------------------------------------------------------------------------------------------------------------------------
aggregate(vg_sales[,6:10], by=list(cluster = km.res1$cluster), mean)


## ------------------------------------------------------------------------------------------------------------------------------
cl.km1 <- km.res1$cluster
pairs(vg_sales.scaled, gap=0, pch=cl.km1, col = c("#2E9FDF", "#FC4E07", "#CF3E4B", "#287233")[cl.km1])


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(km.res1,
             data = vg_sales.scaled,
             palette = c("#2E9FDF", "#FC4E07", "#CF3E4B", "#287233"),
             ellipse.type = "euclid",
             star.plot = TRUE,
             repel = TRUE,
             ggtheme = theme_minimal())


## ------------------------------------------------------------------------------------------------------------------------------
library(cluster)
library(factoextra)
fviz_nbclust(vg_sales.scaled, pam, method = "silhouette")+
  theme_classic()



## ------------------------------------------------------------------------------------------------------------------------------
pam.res <- pam(vg_sales.scaled, 2, metric = "euclidean", stand = FALSE)
summary(pam.res)


## ------------------------------------------------------------------------------------------------------------------------------
pam.res1 <- pam.res$cluster
pairs(vg_sales.scaled, gap=0, pch=pam.res1, col = c("#2E9FDF", "#FC4E07")[pam.res1])

## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(pam.res,
             data = vg_sales.scaled,
             palette = c("#2E9FDF", "#FC4E07"),
             ellipse.type = "euclid",
             star.plot = TRUE,
             repel = TRUE,
             ggtheme = theme_minimal())


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(pam.res, geom = "point", frame.type = "norm")


## ------------------------------------------------------------------------------------------------------------------------------
library(clustertend)
hopkins(vg_sales.scaled, n = nrow(vg_sales.scaled)-1)



## ------------------------------------------------------------------------------------------------------------------------------
s.res.alm <- eclust(vg_sales.scaled, "hclust", k = 2,
                method = "average", graph = FALSE)
silinfo <- s.res.alm$silinfo
silinfo$avg.width



## ------------------------------------------------------------------------------------------------------------------------------
fviz_silhouette(s.res.alm)


## ------------------------------------------------------------------------------------------------------------------------------
silinfo$clus.avg.widths


## ------------------------------------------------------------------------------------------------------------------------------
sil <- s.res.alm$silinfo$widths[, 1:3]
neg_sil_index<- which(sil[,'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]



## ------------------------------------------------------------------------------------------------------------------------------
s.res.wlm <- eclust(vg_sales.scaled, "hclust", k = 3,
                method = "ward", graph = FALSE)
silinfo1 <- s.res.wlm$silinfo
silinfo1$avg.width

## ------------------------------------------------------------------------------------------------------------------------------
fviz_silhouette(s.res.wlm)

## ------------------------------------------------------------------------------------------------------------------------------
silinfo1$clus.avg.widths


## ------------------------------------------------------------------------------------------------------------------------------
sil1 <- s.res.wlm$silinfo$widths[, 1:3]
neg_sil_index1<- which(sil1[,'sil_width'] < 0)
sil1[neg_sil_index1, , drop = FALSE]


## ------------------------------------------------------------------------------------------------------------------------------
s.res.slm <- eclust(vg_sales.scaled, "hclust", k = 2,
                method = "single", graph = FALSE)
silinfo2 <- s.res.slm$silinfo
silinfo2$avg.width

## ------------------------------------------------------------------------------------------------------------------------------
fviz_silhouette(s.res.slm)

## ------------------------------------------------------------------------------------------------------------------------------
silinfo2$clus.avg.widths


## ------------------------------------------------------------------------------------------------------------------------------
sil2 <- s.res.slm$silinfo$widths[, 1:3]
neg_sil_index2 <- which(sil2[,'sil_width'] < 0)
sil2[neg_sil_index2, , drop = FALSE]


## ------------------------------------------------------------------------------------------------------------------------------
s.res.clm <- eclust(vg_sales.scaled, "hclust", k = 2,
                method = "complete", graph = FALSE)
silinfo3 <- s.res.clm$silinfo
silinfo3$avg.width

## ------------------------------------------------------------------------------------------------------------------------------
fviz_silhouette(s.res.clm)


## ------------------------------------------------------------------------------------------------------------------------------
silinfo3$clus.avg.widths


## ------------------------------------------------------------------------------------------------------------------------------
sil3 <- s.res.clm$silinfo$widths[, 1:3]
neg_sil_index3<- which(sil3[,'sil_width'] < 0)
sil3[neg_sil_index3, , drop = FALSE]


## ------------------------------------------------------------------------------------------------------------------------------
s.res.km <- eclust(vg_sales.scaled, "hclust", k = 5,
                method = "kmeans", graph = FALSE)
silinfo4 <- s.res.km$silinfo
silinfo4$avg.width

## ------------------------------------------------------------------------------------------------------------------------------
fviz_silhouette(s.res.km)

## ------------------------------------------------------------------------------------------------------------------------------
silinfo4$clus.avg.widths


## ------------------------------------------------------------------------------------------------------------------------------
sil4 <- s.res.km$silinfo$widths[, 1:3]
neg_sil_index4<- which(sil4[,'sil_width'] < 0)
sil4[neg_sil_index4, , drop = FALSE]


## ------------------------------------------------------------------------------------------------------------------------------
s.res.pam <- eclust(vg_sales.scaled, "hclust", k = 2,
                method = "pam", graph = FALSE)
silinfo5 <- s.res.pam$silinfo
silinfo5$avg.width


## ------------------------------------------------------------------------------------------------------------------------------
fviz_silhouette(s.res.pam)

## ------------------------------------------------------------------------------------------------------------------------------
silinfo5 <- s.res.pam$silinfo
silinfo5$avg.width


## ------------------------------------------------------------------------------------------------------------------------------
sil5 <- s.res.pam$silinfo$widths[, 1:3]
neg_sil_index5<- which(sil5[,'sil_width'] < 0)
sil5[neg_sil_index5, , drop = FALSE]



## ------------------------------------------------------------------------------------------------------------------------------
library(fpc)
di <- cluster.stats(dist(vg_sales.scaled), s.res.alm$cluster)
di$dunn



## ------------------------------------------------------------------------------------------------------------------------------
di1 <- cluster.stats(dist(vg_sales.scaled), s.res.wlm$cluster)
di1$dunn


## ------------------------------------------------------------------------------------------------------------------------------
di2 <- cluster.stats(dist(vg_sales.scaled), s.res.slm$cluster)
di2$dunn


## ------------------------------------------------------------------------------------------------------------------------------
di3 <- cluster.stats(dist(vg_sales.scaled), s.res.clm$cluster)
di3$dunn


## ------------------------------------------------------------------------------------------------------------------------------
di4 <- cluster.stats(dist(vg_sales.scaled), s.res.km$cluster)
di4$dunn


## ------------------------------------------------------------------------------------------------------------------------------
di5 <- cluster.stats(dist(vg_sales.scaled), s.res.pam$cluster)
di5$dunn


## ------------------------------------------------------------------------------------------------------------------------------
table(vg_sales$Genre, s.res.alm$cluster)


## ------------------------------------------------------------------------------------------------------------------------------
genre <- as.numeric(vg_sales$Genre)
extern <- cluster.stats(d = dist(vg_sales), genre, s.res.alm$cluster)
extern$corrected.rand



## ------------------------------------------------------------------------------------------------------------------------------
table(vg_sales$Genre, s.res.wlm$cluster)


## ------------------------------------------------------------------------------------------------------------------------------

extern1 <- cluster.stats(d = dist(vg_sales), genre, s.res.wlm$cluster)
extern1$corrected.rand



## ------------------------------------------------------------------------------------------------------------------------------
table(vg_sales$Genre, s.res.slm$cluster)


## ------------------------------------------------------------------------------------------------------------------------------
extern2 <- cluster.stats(d = dist(vg_sales), genre, s.res.slm$cluster)
extern2$corrected.rand



## ------------------------------------------------------------------------------------------------------------------------------
table(vg_sales$Genre, s.res.clm$cluster)


## ------------------------------------------------------------------------------------------------------------------------------
extern3 <- cluster.stats(d = dist(vg_sales), genre, s.res.clm$cluster)
extern3$corrected.rand



## ------------------------------------------------------------------------------------------------------------------------------
table(vg_sales$Genre, s.res.km$cluster)


## ------------------------------------------------------------------------------------------------------------------------------

extern4 <- cluster.stats(d = dist(vg_sales), genre, s.res.km$cluster)
extern4$corrected.rand



## ------------------------------------------------------------------------------------------------------------------------------
table(vg_sales$Genre, s.res.pam$cluster)


## ------------------------------------------------------------------------------------------------------------------------------
extern5 <- cluster.stats(d = dist(vg_sales), genre, s.res.pam$cluster)
extern5$corrected.rand



## ------------------------------------------------------------------------------------------------------------------------------
library(clValid)
clmethods <- c ("hierarchical", "kmeans", "pam")
intern_vg.eucl<-clValid(vg_sales.scaled, nClust=2:5, clMethods= clmethods, metric="euclidean", validation="internal")
summary(intern_vg.eucl)


## ------------------------------------------------------------------------------------------------------------------------------
clmethods1 <- c ("hierarchical", "kmeans", "pam")
intern_vg.man<-clValid(vg_sales.scaled, nClust=2:5, clMethods= clmethods, metric="manhattan", validation="internal")
summary(intern_vg.eucl)


## ------------------------------------------------------------------------------------------------------------------------------
clmethods2 <- c("hierarchical","kmeans","pam")

stab <- clValid(vg_sales.scaled, nClust = 2:5, clMethods = clmethods, validation = "stability")

optimalScores(stab)


## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(list(data=vg_sales.scaled, cluster =  d.alm),geom="point",ellipse.type="convex",palette = c("#00AFBB", "#E7B800"),repel = TRUE, show.clust.cent = FALSE)



## ------------------------------------------------------------------------------------------------------------------------------
fviz_cluster(pam.res, geom = "point", frame.type = "norm")

