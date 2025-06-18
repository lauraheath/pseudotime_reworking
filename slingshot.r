setwd("~/pseudotime_v2/")
dir.create("~/pseudotime_v2/data_objects")


install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
BiocManager::install("tradeSeq")
BiocManager::install("clusterExperiment")
install.packages('mclust')
install.packages('psych')
library(tradeSeq)
library(monocle3)
library(dplyr)
library(magrittr)
library(clusterExperiment)
library(slingshot)
library(psych)
library(RColorBrewer)
library(uwot)
library(mclust, quietly = TRUE)
synapser::synLogin()



#use matrices/metadata from monocle3 code
#create singlecellexperiment object

Dat <- readRDS(file='data_objects/Female_RNAseq_ADgenes_datamatrix.RDS')
metadata <- read.csv(file='data_objects/Female_RNAseq_metadata.csv')
metadata$X<-NULL


Dat <- readRDS(file='data_objects/Male_RNAseq_ADgenes_datamatrix.RDS')
metadata <- read.csv(file='data_objects/Male_RNAseq_metadata.csv')
metadata$X<-NULL


sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(Dat)),
  colData = metadata
)

pca <- prcomp(t(assays(sce)$counts), scale. = FALSE)
rd1 <- pca$x[,1:2]

p <- plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
p


rd2 <- uwot::umap(t(assays(sce)$counts))
colnames(rd2) <- c('UMAP1', 'UMAP2')

p2 <- plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
p2

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)


cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1


plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2
plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')



summary(sce$slingPseudotime_1)



#library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

cl2 <- as.data.frame(cl2)
cl2$specimenID <- rownames(cl2)

rd1a <- as.data.frame(rd1)
rd1a$specimenID <- rownames(rd1a)



slingpseudo <- slingPseudotime((sce))

slingpseudo <- as.data.frame(slingpseudo)
slingpseudo$specimenID <- rownames(slingpseudo)


df <- colData(sce)
gmm <- as.data.frame(subset(df, select=c(specimenID, GMM)))





sling_meta <- left_join(metadata, slingpseudo, by='specimenID')
sling_meta <- left_join(sling_meta, cl2, by='specimenID')
sling_meta <- left_join(sling_meta, rd1a, by='specimenID')
sling_meta <- left_join(sling_meta, gmm, by='specimenID')

table(sling_meta$cl2, sling_meta$diagnosis)
table(sling_meta$GMM, sling_meta$diagnosis)
ggplot(sling_meta, aes(x=diagnosis, y=Lineage1)) + geom_boxplot()

sling_meta$cl2 <- as.factor(sling_meta$cl2)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=cl2)) + geom_point()

sling_meta$GMM <- as.factor(sling_meta$GMM)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=GMM)) + geom_point()







sce2 <- slingshot(sce, clusterLabels = 'kmeans', start.clus ='1', reducedDim = 'PCA')

#for female dataset : need to reverse the direction of slingshot (start with cluster 2)
sce3 <- slingshot(sce, clusterLabels = 'GMM', start.clus ='2', reducedDim = 'PCA')


slingpseudo <- slingPseudotime((sce2))

slingpseudo <- as.data.frame(slingpseudo)
slingpseudo$specimenID <- rownames(slingpseudo)


df <- colData(sce2)
gmm <- as.data.frame(subset(df, select=c(specimenID, GMM)))


sling_meta <- left_join(metadata, slingpseudo, by='specimenID')
sling_meta <- left_join(sling_meta, cl2, by='specimenID')
sling_meta <- left_join(sling_meta, rd1a, by='specimenID')
sling_meta <- left_join(sling_meta, gmm, by='specimenID')

table(sling_meta$cl2, sling_meta$diagnosis)
table(sling_meta$GMM, sling_meta$diagnosis)
ggplot(sling_meta, aes(x=diagnosis, y=Lineage1)) + geom_boxplot()

sling_meta$cl2 <- as.factor(sling_meta$cl2)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=cl2)) + geom_point()

sling_meta$GMM <- as.factor(sling_meta$GMM)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=GMM)) + geom_point()




sce3 <- slingshot(sce, clusterLabels = 'GMM', start.clus ='2', reducedDim = 'PCA')
slingpseudo2 <- slingPseudotime((sce3))

slingpseudo2 <- as.data.frame(slingpseudo2)
slingpseudo2$specimenID <- rownames(slingpseudo2)


df <- colData(sce3)
gmm <- as.data.frame(subset(df, select=c(specimenID, GMM)))


sling_meta2 <- left_join(metadata, slingpseudo2, by='specimenID')
sling_meta2 <- left_join(sling_meta, cl2, by='specimenID')
sling_meta2 <- left_join(sling_meta, rd1a, by='specimenID')
sling_meta2 <- left_join(sling_meta, gmm, by='specimenID')

table(sling_meta2$cl2, sling_meta$diagnosis)
table(sling_meta2$GMM, sling_meta$diagnosis)
ggplot(sling_meta2, aes(x=diagnosis, y=Lineage1)) + geom_boxplot()



metaMon <- sling_meta

metaMon$braaksc <- as.factor(metaMon$braaksc)
#tiff(file='figures/FEMALE_mono3_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='figures/MALE_mets_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(metaMon, aes(x=braaksc, y=scale(Lineage1, center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
#dev.off()

metaMon$ceradsc <- as.factor(metaMon$ceradsc)
g <- ggplot2::ggplot(metaMon, aes(x=ceradsc, y=scale(Lineage1,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g



#run logistic regression comparing pseudotime between cases and controls only
casecontrol <- subset(metaMon, metaMon$diagnosis=='AD'|metaMon$diagnosis=='CT')
casecontrol$diag2 <- ifelse(casecontrol$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ Lineage1,casecontrol,family='binomial'))
g <- ggplot(casecontrol,aes(x=diagnosis,
                            y=Lineage1,
                            color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g


#save the slingshot pseudotimes:
df <- subset(sling_meta, select=c(specimenID, Lineage1))
names(df)[names(df) == 'Lineage1'] <- 'Pseudotime_rnaseq_slingshot'
#write.csv(df, file='data_objects/Female_RNAseq_slingshot_pseudotime.csv', row.names = FALSE)
write.csv(df, file='data_objects/Male_RNAseq_slingshot_pseudotime.csv', row.names = FALSE)



#join the three pseudotimes:
mono3 <- read.csv(file='data_objects/Female_RNAseq_mono3_pseudotime.csv')
mono3 <- read.csv(file='data_objects/Male_RNAseq_mono3_pseudotime.csv')
mono3$X<-NULL
#upload monocle2 pseudotimes:
mono2Obj <- synapser::synGet('syn45272849')
mono2 <- read.csv(mono2Obj$path)
mono2$X<-NULL
names(mono2)[names(mono2) == 'Pseudotime_rnaseq'] <- 'Pseudotime_rnaseq_Mono2'
names(mono2)[names(mono2) == 'SampleID_rnaseq'] <- 'specimenID'
mono2 <- subset(mono2, select=c(specimenID, Pseudotime_rnaseq_Mono2))

metadata2 <- left_join(metadata, mono2, by='specimenID')
metadata2 <- left_join(metadata2, mono3, by='specimenID')
metadata2 <- left_join(metadata2, df, by='specimenID')


#3-way correlation plots
df <- subset(metadata2, select=c(Pseudotime_rnaseq_Mono2,Pseudotime_rnaseq_Mono3,Pseudotime_rnaseq_slingshot))
cor(df[,unlist(lapply(df, is.numeric))])


#create pairs plot
pairs.panels(df)










########### proteomics slingshot analysis #################

#use matrices/metadata from monocle3 code
#create singlecellexperiment object

Dat <- readRDS(file='data_objects/Female_prot_ADgenes_matrix.rds')
metadata <- read.csv(file='data_objects/Female_prot_ADgenes_metadata.csv')
metadata$X<-NULL


Dat <- readRDS(file='data_objects/Male_prot_ADgenes_matrix.rds')
metadata <- read.csv(file='data_objects/Male_prot_ADgenes_metadata.csv')
metadata$X<-NULL


sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(Dat)),
  colData = metadata
)

pca <- prcomp(t(assays(sce)$counts), scale. = FALSE)
rd1 <- pca$x[,1:2]

p <- plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
p

#library(uwot)
rd2 <- uwot::umap(t(assays(sce)$counts))
colnames(rd2) <- c('UMAP1', 'UMAP2')

p2 <- plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
p2

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

#library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

#library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2
plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')



summary(sce$slingPseudotime_1)



#library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

cl2 <- as.data.frame(cl2)
cl2$batch.channel <- rownames(cl2)

rd1a <- as.data.frame(rd1)
rd1a$batch.channel <- rownames(rd1a)



slingpseudo <- slingPseudotime((sce))

slingpseudo <- as.data.frame(slingpseudo)
slingpseudo$batch.channel <- rownames(slingpseudo)


df <- colData(sce)
gmm <- as.data.frame(subset(df, select=c(batch.channel, GMM)))





sling_meta <- left_join(metadata, slingpseudo, by='batch.channel')
sling_meta <- left_join(sling_meta, cl2, by='batch.channel')
sling_meta <- left_join(sling_meta, rd1a, by='batch.channel')
sling_meta <- left_join(sling_meta, gmm, by='batch.channel')

table(sling_meta$cl2, sling_meta$diagnosis)
table(sling_meta$GMM, sling_meta$diagnosis)
ggplot(sling_meta, aes(x=diagnosis, y=Lineage1)) + geom_boxplot()

sling_meta$cl2 <- as.factor(sling_meta$cl2)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=cl2)) + geom_point()

sling_meta$GMM <- as.factor(sling_meta$GMM)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=GMM)) + geom_point()










#for female dataset : need to reverse the direction of slingshot (start with cluster 2)
sce2 <- slingshot(sce, clusterLabels = 'GMM', start.clus ='2', reducedDim = 'PCA')


slingpseudo <- slingPseudotime((sce2))

slingpseudo <- as.data.frame(slingpseudo)
slingpseudo$batch.channel <- rownames(slingpseudo)


df <- colData(sce2)
gmm <- as.data.frame(subset(df, select=c(batch.channel, GMM)))


sling_meta <- left_join(metadata, slingpseudo, by='batch.channel')
sling_meta <- left_join(sling_meta, cl2, by='batch.channel')
sling_meta <- left_join(sling_meta, rd1a, by='batch.channel')
sling_meta <- left_join(sling_meta, gmm, by='batch.channel')

table(sling_meta$cl2, sling_meta$diagnosis)
table(sling_meta$GMM, sling_meta$diagnosis)
ggplot(sling_meta, aes(x=diagnosis, y=Lineage1)) + geom_boxplot()

sling_meta$cl2 <- as.factor(sling_meta$cl2)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=cl2)) + geom_point()

sling_meta$GMM <- as.factor(sling_meta$GMM)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=GMM)) + geom_point()


metaMon <- sling_meta

metaMon$braaksc <- as.factor(metaMon$braaksc)
#tiff(file='figures/FEMALE_mono3_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='figures/MALE_mets_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(metaMon, aes(x=braaksc, y=scale(Lineage1, center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
#dev.off()

metaMon$ceradsc <- as.factor(metaMon$ceradsc)
g <- ggplot2::ggplot(metaMon, aes(x=ceradsc, y=scale(Lineage1,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g



#run logistic regression comparing pseudotime between cases and controls only
casecontrol <- subset(metaMon, metaMon$diagnosis=='AD'|metaMon$diagnosis=='control')
casecontrol$diag2 <- ifelse(casecontrol$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ Lineage1,casecontrol,family='binomial'))
g <- ggplot(casecontrol,aes(x=diagnosis,
                            y=Lineage1,
                            color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g


#save the slingshot pseudotimes:
df <- subset(sling_meta, select=c(batch.channel, Lineage1))
names(df)[names(df) == 'Lineage1'] <- 'Pseudotime_prot_slingshot'
#write.csv(df, file='data_objects/Female_prot_slingshot_pseudotime.csv', row.names = FALSE)
write.csv(df, file='data_objects/Male_prot_slingshot_pseudotime.csv', row.names = FALSE)




#join the three pseudotimes:
sling <- read.csv(file='data_objects/Female_prot_slingshot_pseudotime.csv')
mono3 <- read.csv(file='data_objects/Female_prot_mono3_pseudotime.csv')


#mono3 <- read.csv(file='data_objects/Male_prot_mono3_pseudotime.csv')
mono3$X<-NULL
#upload monocle2 pseudotimes:
mono2Obj <- synapser::synGet('syn45272849')
mono2 <- read.csv(mono2Obj$path)
mono2$X<-NULL
names(mono2)[names(mono2) == 'Pseudotime_prot'] <- 'Pseudotime_prot_Mono2'
names(mono2)[names(mono2) == 'SampleID_prot'] <- 'batch.channel'
mono2 <- subset(mono2, select=c(batch.channel, Pseudotime_prot_Mono2))


metadata2 <- left_join(sling, mono2, by='batch.channel')
metadata2 <- left_join(metadata2, mono3, by='batch.channel')


#3-way correlation plots
df <- subset(metadata2, select=c(Pseudotime_prot_Mono2,Pseudotime_prot_Mono3,Pseudotime_prot_slingshot))
cor(df[,unlist(lapply(df, is.numeric))])


#create pairs plot
pairs.panels(df)








########### run slingshot on metabolomics #############
#use matrices/metadata from monocle3 code
#create singlecellexperiment object

Dat <- readRDS(file='data_objects/Female_mets_ADgenes_matrix.rds')
metadata <- read.csv(file='data_objects/Female_mets_metadata.csv')
metadata$X<-NULL


Dat <- readRDS(file='data_objects/Male_mets_ADgenes_matrix.rds')
metadata <- read.csv(file='data_objects/Male_mets_metadata.csv')
metadata$X<-NULL


#force the projid to be characters by adding a letter
colnames(Dat) <- paste(colnames(Dat), "a", sep="")
metadata$projid <- paste(metadata$projid, "a", sep="")


sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(Dat)),
  colData = metadata
)

pca <- prcomp(t(assays(sce)$counts), scale. = FALSE)
rd1 <- pca$x[,1:2]

p <- plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
p

#library(uwot)
rd2 <- uwot::umap(t(assays(sce)$counts))
colnames(rd2) <- c('UMAP1', 'UMAP2')

p2 <- plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
p2

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

#library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

#library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2
plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'kmeans', reducedDim = 'PCA')



summary(sce$slingPseudotime_1)



#library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

cl2 <- as.data.frame(cl2)
cl2$projid <- rownames(cl2)

rd1a <- as.data.frame(rd1)
rd1a$projid <- rownames(rd1a)



slingpseudo <- slingPseudotime((sce))

slingpseudo <- as.data.frame(slingpseudo)
slingpseudo$projid <- rownames(slingpseudo)


df <- colData(sce)
gmm <- as.data.frame(subset(df, select=c(projid, GMM)))

sling_meta <- left_join(metadata, slingpseudo, by='projid')
sling_meta <- left_join(sling_meta, cl2, by='projid')
sling_meta <- left_join(sling_meta, rd1a, by='projid')
sling_meta <- left_join(sling_meta, gmm, by='projid')

table(sling_meta$cl2, sling_meta$diagnosis)
table(sling_meta$GMM, sling_meta$diagnosis)
ggplot(sling_meta, aes(x=diagnosis, y=Lineage1)) + geom_boxplot()

sling_meta$cl2 <- as.factor(sling_meta$cl2)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=cl2)) + geom_point()

sling_meta$GMM <- as.factor(sling_meta$GMM)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=GMM)) + geom_point()










#for female and male dataset : need to reverse the direction of slingshot (start with cluster 2)
sce2 <- slingshot(sce, clusterLabels = 'kmeans', start.clus ='1', reducedDim = 'PCA')


slingpseudo <- slingPseudotime((sce2))

slingpseudo <- as.data.frame(slingpseudo)
slingpseudo$projid <- rownames(slingpseudo)


df <- colData(sce2)
gmm <- as.data.frame(subset(df, select=c(projid, GMM)))


sling_meta <- left_join(metadata, slingpseudo, by='projid')
sling_meta <- left_join(sling_meta, cl2, by='projid')
sling_meta <- left_join(sling_meta, rd1a, by='projid')
sling_meta <- left_join(sling_meta, gmm, by='projid')

table(sling_meta$cl2, sling_meta$diagnosis)
table(sling_meta$GMM, sling_meta$diagnosis)
ggplot(sling_meta, aes(x=diagnosis, y=Lineage1)) + geom_boxplot()

sling_meta$cl2 <- as.factor(sling_meta$cl2)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=cl2)) + geom_point()

sling_meta$GMM <- as.factor(sling_meta$GMM)
ggplot(sling_meta, aes(x=PC1, y=PC2, color=GMM)) + geom_point()


metaMon <- sling_meta

metaMon$braaksc <- as.factor(metaMon$braaksc)
#tiff(file='figures/FEMALE_mono3_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='figures/MALE_mets_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(metaMon, aes(x=braaksc, y=scale(Lineage1, center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
#dev.off()

metaMon$ceradsc <- as.factor(metaMon$ceradsc)
g <- ggplot2::ggplot(metaMon, aes(x=ceradsc, y=scale(Lineage1,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g



#run logistic regression comparing pseudotime between cases and controls only
casecontrol <- subset(metaMon, metaMon$diagnosis=='AD'|metaMon$diagnosis=='CT')
casecontrol$diag2 <- ifelse(casecontrol$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ Lineage1,casecontrol,family='binomial'))
g <- ggplot(casecontrol,aes(x=diagnosis,
                            y=Lineage1,
                            color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g


#save the slingshot pseudotimes:
df <- subset(sling_meta, select=c(projid, Lineage1))
names(df)[names(df) == 'Lineage1'] <- 'Pseudotime_mets_slingshot'

write.csv(df, file='data_objects/Female_mets_slingshot_pseudotime.csv', row.names = FALSE)
#write.csv(df, file='data_objects/Male_mets_slingshot_pseudotime.csv', row.names = FALSE)




#join the three pseudotimes:
sling <- read.csv(file='data_objects/Female_prot_slingshot_pseudotime.csv')
mono3 <- read.csv(file='data_objects/Female_prot_mono3_pseudotime.csv')


#mono3 <- read.csv(file='data_objects/Male_prot_mono3_pseudotime.csv')
mono3$X<-NULL
#upload monocle2 pseudotimes:
mono2Obj <- synapser::synGet('syn45272849')
mono2 <- read.csv(mono2Obj$path)
mono2$X<-NULL
names(mono2)[names(mono2) == 'Pseudotime_prot'] <- 'Pseudotime_prot_Mono2'
names(mono2)[names(mono2) == 'SampleID_prot'] <- 'batch.channel'
mono2 <- subset(mono2, select=c(batch.channel, Pseudotime_prot_Mono2))


metadata2 <- left_join(sling, mono2, by='batch.channel')
metadata2 <- left_join(metadata2, mono3, by='batch.channel')


#3-way correlation plots
df <- subset(metadata2, select=c(Pseudotime_prot_Mono2,Pseudotime_prot_Mono3,Pseudotime_prot_slingshot))
cor(df[,unlist(lapply(df, is.numeric))])


#create pairs plot
pairs.panels(df)

