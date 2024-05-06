install.packages('ggpubr')
install.packages('VennDiagram')
install.packages('scatterplot3d')
install.packages("plotly")
library(VennDiagram)
library(gridExtra)
library(scatterplot3d)
library(plotly)
synapser::synLogin()

#set working directory and create folders to stash data objects and figures
setwd("~/pseudotime_reworking/")
#dir.create("data_objects")
#dir.create("figures")


#upload metabolomics scripts:
# p1 <- "syn45147406"
# synGet(p1, downloadLocation = "code/")
# p2 <- "syn45147407"
# synGet(p2, downloadLocation = "code/")



#upload proteomics pseudotimes & states by sex:
p1 <- synapser::synGet('syn35317137')
female_prots <- read.csv(p1$path)

p1 <- synapser::synGet('syn35317790')
male_prots <- read.csv(p1$path)

#combine
allprots <- rbind(female_prots, male_prots)


#change column names to differentiate from RNAseq
names(allprots)[names(allprots) == 'Pseudotime'] <- 'Pseudotime_prot'
names(allprots)[names(allprots) == 'State'] <- 'State_prot'
names(allprots)[names(allprots) == 'SampleID'] <- 'SampleID_prot'
names(allprots)[names(allprots) == 'pseudotime_sc'] <- 'pseudotime_sc_prot'
allprots$msex[allprots$msex == 0] <- 'female'
allprots$msex[allprots$msex == 1] <- 'male'
allprots$diagnosis[allprots$diagnosis == 'control'] <- 'CT'
allprots$diagnosis[allprots$diagnosis == 'other'] <- 'OTHER'
allprots$batch<-NULL
allprots$apoe_genotype<-NULL
allprots$educ<-NULL
allprots$age_death<-NULL
allprots$pmi<-NULL



#upload rnaseq pseudotimes & states
p2 <- synapser::synGet('syn38354143')
female_rnaseq <- read.csv(p2$path)

p2 <- synapser::synGet('syn38348029')
male_rnaseq <- read.csv(p2$path)
#combine
allRNAseq <- rbind(female_rnaseq, male_rnaseq)

names(allRNAseq)[names(allRNAseq) == 'Pseudotime'] <- 'Pseudotime_rnaseq'
names(allRNAseq)[names(allRNAseq) == 'State'] <- 'State_rnaseq'
names(allRNAseq)[names(allRNAseq) == 'SampleID'] <- 'SampleID_rnaseq'
names(allRNAseq)[names(allRNAseq) == 'pseudotime_sc'] <- 'pseudotime_sc_rnaseq'
allRNAseq$batch<-NULL
allRNAseq$tissue<-NULL
allRNAseq$apoe4_allele<-NULL
allRNAseq$RIN<-NULL
allRNAseq$pmi<-NULL


#upload metabolomics (use supervised output)
p2 <- synapser::synGet('syn50912936')
female_mets <- read.csv(p2$path)

p2 <- synapser::synGet('syn50912976')
male_mets <- read.csv(p2$path)
#combine
allmets <- rbind(female_mets, male_mets)

names(allmets)[names(allmets) == 'Pseudotime'] <- 'Pseudotime_mets'
names(allmets)[names(allmets) == 'State'] <- 'State_mets'
names(allmets)[names(allmets) == 'SampleID'] <- 'SampleID_mets'
names(allmets)[names(allmets) == 'pseudotime_sc'] <- 'pseudotime_sc_mets'
allmets$apoe_genotype<-NULL
allmets$msex[allmets$msex == 0] <- 'female'
allmets$msex[allmets$msex == 1] <- 'male'
allmets$educ<-NULL
allmets$pmi<-NULL


#merge all three omics plus clinical metadata; keep all observations
pstimes <- dplyr::full_join(allprots, allRNAseq)
pstimes <- dplyr::full_join(pstimes, allmets)

#merge all three omics, but only those with values
pstimes2 <- dplyr::inner_join(allRNAseq, allprots)
pstimes2 <- dplyr::inner_join(pstimes2, allmets)

#create a flag for whether any given patient has a sample from each datatype 
pstimes$rna_samples = 0
pstimes$rna_samples[!is.na(pstimes$Pseudotime_rnaseq)] <- 1
pstimes$prot_samples = 0
pstimes$prot_samples[!is.na(pstimes$Pseudotime_prot)] <- 1
pstimes$met_samples = 0
pstimes$met_samples[!is.na(pstimes$Pseudotime_met)] <- 1

table(pstimes$prot_samples, pstimes$rna_samples)
table(pstimes$prot_samples, pstimes$met_samples)
table(pstimes$rna_samples, pstimes$met_samples)



#save to the workspace:
write.csv(pstimes, file="data_objects/ALL_OMICS_pseudotimes.csv")
file <- synapser::File(path='data_objects/ALL_OMICS_pseudotimes.csv', parentId='syn44292253')
file <- synapser::synStore(file)




#venn diagram of which patients overlap between rnaseq, protein, and metabolomics studies
rna_ps <- subset(pstimes, pstimes$rna_samples==1)
prot_ps <- subset(pstimes, pstimes$prot_samples==1)
met_ps <- subset(pstimes, pstimes$met_samples==1)
rna_patients <- rna_ps$projid
prot_patients <- prot_ps$projid
met_patients <- met_ps$projid

venn.diagram(
  x=list(rna_patients, prot_patients, met_patients),
  category.names=c("RNAseq IDs","Proteomics IDs", "Metabolomics IDs"),
  filename='figures/ALL_Patient_overlap_venn.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer",
  cat.pos=c(-24,18,135),
  cat.dist=c(0.05,0.05,0.05)
)




#run correlations between proteomics pseudotimes and rnaseq pseudotimes (use scaled pseudotimes)
cor(pstimes$pseudotime_sc_rnaseq, pstimes$pseudotime_sc_prot, method="pearson", use="pairwise.complete.obs")
cor(pstimes$pseudotime_sc_rnaseq, pstimes$pseudotime_sc_met, method="pearson", use="pairwise.complete.obs")
cor(pstimes$pseudotime_sc_prot, pstimes$pseudotime_sc_met, method="pearson", use="pairwise.complete.obs")

#cor(pstimes$pseudotime_sc_rnaseq, pstimes$pseudotime_sc_prot, method="spearman", use="pairwise.complete.obs")


plot(pstimes$pseudotime_sc_rnaseq, pstimes$pseudotime_sc_prot)
abline(lm(pstimes$pseudotime_sc_rnaseq~pstimes$pseudotime_sc_prot))

plot(pstimes$pseudotime_sc_rnaseq, pstimes$pseudotime_sc_met)
abline(lm(pstimes$pseudotime_sc_rnaseq~pstimes$pseudotime_sc_met))

plot(pstimes$pseudotime_sc_prot, pstimes$pseudotime_sc_met)
abline(lm(pstimes$pseudotime_sc_prot~pstimes$pseudotime_sc_met))


#restrict dataframe to patients with all three omics:
#use https://waldyrious.net/viridis-palette-generator/ to determine color palette (to match previous plots)
colors <- c("#440154FF","#22A884FF","#FDE725FF")
pstimes2$diagnosis2 <- as.factor(pstimes2$diagnosis)
colors <- colors[as.numeric(pstimes2$diagnosis2)]

scatter <- scatterplot3d(pstimes2$pseudotime_sc_mets, pstimes2$pseudotime_sc_prot, pstimes2$pseudotime_sc_rnaseq, color=colors, pch=16)
legend("right", legend = levels(pstimes2$diagnosis2), 
       col = c("#440154FF","#22A884FF","#FDE725FF"), 
       pch=16, horiz=TRUE, xpd = TRUE, inset = -.2)


#color by braak score
colors <- c("#440154", "#443983", "#31688e","#21918c","#35b779","#90d743","#fde725")
pstimes2$braaksc <- as.factor(pstimes2$braaksc)
colors <- colors[as.numeric(pstimes2$braaksc)]
scatterplot3d(pstimes2$pseudotime_sc_mets, pstimes2$pseudotime_sc_prot, pstimes2$pseudotime_sc_rnaseq, color=colors, pch=16)
legend("bottom", legend = levels(pstimes2$braaksc), 
       col = c("#440154", "#443983", "#31688e","#21918c","#35b779","#90d743","#fde725"), 
       pch=16, horiz=TRUE, xpd = TRUE, inset = -.2)


#color by cerad score
colors <- c("#fde725","#35b779","#31688e","#440154")
pstimes2$ceradsc <- as.factor(pstimes2$ceradsc)
colors <- colors[as.numeric(pstimes2$ceradsc)]
scatter <- scatterplot3d(pstimes2$pseudotime_sc_mets, pstimes2$pseudotime_sc_prot, pstimes2$pseudotime_sc_rnaseq, color=colors, pch=16)
legend("bottom", legend = levels(pstimes2$ceradsc), 
       col = c("#fde725","#35b779","#31688e","#440154"), 
       pch=16, horiz=TRUE, xpd = TRUE, inset = -.2)


#color by cogdx score
colors <- c("#440154","#414487","#2a788e","#22a884","#7ad151","#fde725")
pstimes2$cogdx <- as.factor(pstimes2$cogdx)
colors <- colors[as.numeric(pstimes2$cogdx)]
scatter <- scatterplot3d(pstimes2$pseudotime_sc_mets, pstimes2$pseudotime_sc_prot, pstimes2$pseudotime_sc_rnaseq, color=colors, pch=16)
legend("bottom", legend = levels(pstimes2$cogdx), 
       col = c("#440154","#414487","#2a788e","#22a884","#7ad151","#fde725"), 
       pch=16, horiz=TRUE, xpd = TRUE, inset = -.2)



p <- plot_ly(pstimes2, x=~pseudotime_sc_mets, y=~pseudotime_sc_prot, 
             z=~pseudotime_sc_rnaseq, color=~diagnosis) %>%
  add_markers(size=1)
print(p)


#3-way correlation plots
df <- subset(pstimes2, select=c(pseudotime_sc_prot, pseudotime_sc_rnaseq, pseudotime_sc_mets))
cor(df[,unlist(lapply(df, is.numeric))])
#load psych package
library(psych)

#create pairs plot
pairs.panels(df)

df <- subset(pstimes2, select=c(Pseudotime_prot, Pseudotime_rnaseq, Pseudotime_mets))
cor(df[,unlist(lapply(df, is.numeric))])
pairs.panels(df)

#try by sex:
females <- subset(pstimes2, pstimes2$msex=='female')
df <- subset(females, select=c(pseudotime_sc_prot, pseudotime_sc_rnaseq, pseudotime_sc_mets))
cor(df[,unlist(lapply(df, is.numeric))])
pairs.panels(df)

males <- subset(pstimes2, pstimes2$msex=='male')
df <- subset(males, select=c(pseudotime_sc_prot, pseudotime_sc_rnaseq, pseudotime_sc_mets))
cor(df[,unlist(lapply(df, is.numeric))])
pairs.panels(df)




### try a plotly scatterplot

#The following monocle objects were calculated in the scripts UPDATE_rnaseq_lineage.r 
# (in AMP-AD-2.0-transcriptomics-lineage-update repo) and prot_lineage_monocle_rerun.R 
# (in prot-lineage repo), and stored in the shared working space

#female monocle objects (P for proteomics, R for RNAseq):
p4 <- synapser::synGet('syn44293198')
MonRun_FP <- readRDS(p4$path)

p5 <- synapser::synGet('syn44293242')
MonRun_FR <- readRDS(p5$path)

p6 <- synapser::synGet('syn44293269')
MonRun_MP <- readRDS(p6$path)

p7 <- synapser::synGet('syn44293357')
MonRun_MR <- readRDS(p7$path)



#add rna-seq pseudotimes to the monocle object
rnaseq_pstime <- subset(pstimes, select=c(SampleID_rnaseq, pseudotime_sc_rnaseq, pseudotime_sc_prot, rna_samples))
rnaseq_pstime <- subset(rnaseq_pstime, rnaseq_pstime$rna_samples==1)
rnaseq_pstime$rna_samples<-NULL
names(rnaseq_pstime)[names(rnaseq_pstime) == 'SampleID_rnaseq'] <- 'specimenID'



head(pData(F_MonRunT))
tail(pData(F_MonRunT))

#reorder rnaseq_pstime data frame to match the order of the specimenIDs in the monocle CDS:
sampleIDs <- as.data.frame(pData(F_MonRunT))
rnaseq_pstime <- rnaseq_pstime[order(match(rnaseq_pstime$specimenID, sampleIDs$specimenID)),]

F_MonRunT$rnaseq_pseudotime <- rnaseq_pstime$pseudotime_sc_rnaseq
F_MonRunT$prot_pseudotime <- rnaseq_pstime$pseudotime_sc_proteomics

#tiff(file='~/prot-lineage/figures/FEMALE_tree_rnaseq_pstime.tiff',height=150,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_tree_rnaseq_pstime.tiff',height=150,width=100,units='mm',res=300)
g <- plot_cell_trajectory(F_MonRunT,color_by = "prot_pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1) + scale_color_gradient (low="blue", high="red")
g
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(pstimes, x="pseudotime_sc_rnaseq", y="pseudotime_sc_prot", 
                  add="reg.line",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_diagnosis.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_diagnosis.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(pstimes, x="pseudotime_sc_rnaseq", y="pseudotime_sc_prot", color = "diagnosis",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

pstimes_combined$braaksc <- as.factor(pstimes_combined$braaksc)
#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_braak.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(pstimes_combined, x="pseudotime_sc_proteomics", y="pseudotime_sc_rnaseq", color = "braaksc",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_cerad.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", color = "ceradsc",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_cogdx.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", color = "cogdx",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()





#compare DE genes included as feature sets in bulk rnaseq and prot


inters <- intersect(rnaseq_genesF$gene_short_name,protgenes$GeneName)
diffs <- setdiff(rnaseq_genesF$gene_short_name,protgenes$GeneName)

#rerun lineage analysis on both rna-seq and proteomics data using only the shared genes as a feature set
inters <- as.data.frame(inters)
names(inters)[names(inters) == "inters"] <- "gene_short_name"
Log2_Normalized <- readRDS(file="~/prot-lineage/data/Log2_Normalized.rds")
Meta <- readRDS(file="~/prot-lineage/data/Meta.rds")
Log2_Normalized2 <- Log2_Normalized
Log2_Normalized2$proteins <- rownames(Log2_Normalized2)
Log2_Normalized2$proteins <- gsub("\\|.*", "", Log2_Normalized2$proteins)

#function to run monocle analysis
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}

#run get_proteins_script_JG.r to get metadata and log2 normalized protein matrix
synapser::synLogin()

Dat <- Log2_Normalized2
Dat[is.na(Dat)] <- 0

#select only the rows with gene short names that match the intersecting gene short names (some are repeated peptides, same gene)
genes2<-c()
for (gene in unique(c(as.vector(inters$gene_short_name)))){
  if (gene %in% Dat$proteins){
    genes2 <- c(genes2,which(Dat$proteins==gene))
  }
}
length(genes2)
Dat2 <- Dat[genes2,]
dim(Dat2)
Dat2$proteins<-NULL


#Keeping only female data (msex==0 is female, msex==1 is male; run separately for sex-specific analysis)
In_S <- which(Meta$msex == 0)
#In_S <- which(Meta$msex == 1)
Dat2 <- Dat2[,In_S]
Meta2 <- Meta[In_S,]

gene_short_name <- rownames(Dat2)
temp <- Dat2
temp2 <- Meta2


temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc <- factor(temp2$ceradsc,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx,levels = c(1:6))

rownames(temp)<-NULL
rownames(temp2)<-NULL

MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

x <- list()
x$SampleID <- MonRun$batchChannel
x$rnaseqID <- MonRun$rnaseq_id
x$State2 <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$diagnosis <- MonRun$diagnosis
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$apoe <- MonRun$APO
x$educ   <- MonRun$educ
x$pmi <- MonRun$pmi
x$batch <- MonRun$batch
x$mmse <- MonRun$cts_mmse30_lv
x$age_death <- MonRun$age_death
x$rna_seq_sample <- MonRun$rnaseq
x$SampleID <- as.character(x$SampleID)
F_inters_prots <- as.data.frame(x)
F_inters_prots$pseudotime_sc <- scale(F_inters_prots$Pseudotime, center=F)




#### now upload rna-seq data and rerun lineage analysis on intersecting proteins


#function to run monocle analysis
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}

dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- read.delim(dlpfcCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file (rosmap covariates): syn8466814
dlpfcCovObj <- synapser::synGet('syn11024258')
covars <- read.delim(dlpfcCovObj$path,stringsAsFactors = F)

#converting ENSG to gene symbols
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           #host='useast.ensembl.org')
                           host='uswest.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}
Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}
Dat$gene_short_name <- Make.Gene.Symb(Dat$ensembl_gene_id)

#select only the rows with gene short names that match the intersecting gene short names (some are repeated peptides, same gene)
genes2<-c()
for (gene in unique(c(as.vector(inters$gene_short_name)))){
  if (gene %in% Dat$gene_short_name){
    genes2 <- c(genes2,which(Dat$gene_short_name==gene))
  }
}
length(genes2)
Dat2 <- Dat[genes2,]
dim(Dat2)


Names <- colnames(Dat2)

for (i in 1:length(Names)){
  
  Names[i] <- substring(Names[i],2)
  
}


colnames(Dat2) <- Names
cNames <- covars$SampleID
l <- length(Names)

#deleting columns not in the covariate list
temp <- rep(T,l)
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    temp[i] <- F
  }
}

In <- which(temp)
#print(temp)
Dat2 <- Dat2[,In]

#deleting extra rows in covariate list
Names <- Names[In]
l <- length(cNames)
temp <- rep(T,l)
for (i in 1:l){
  if (!(cNames[i] %in% Names)){
    temp[i] <- F
  }
}
In <- which(temp)
covars <- covars[In,]

ColNorm <- function(Dat2){
  
  M = max(colSums(Dat2))
  l <- length(colnames(Dat2))
  
  for( i in 1:l){
    
    Dat2[,i] = Dat2[,i]*(M/sum(Dat2[,i]))
    
  }
  
  return(Dat2)
}

DatNorm <- ColNorm(Dat2)

#removing bad batches
DatNorm <- DatNorm[,covars$Batch<7]
covars <- covars[covars$Batch<7,] 


#Keeping only female data 
#Sex <- 'FEMALE'
In_S <- which(covars$msex == 0)
DatNorm2 <- DatNorm[,In_S]
covars2 <- covars[In_S,]

temp <- DatNorm2
temp2 <- covars2

rosmapObj <- synapser::synGet('syn3191087')
rosmap <- data.table::fread(rosmapObj$path,data.table=F)

#add in braak score & cerad score
rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)
rosmapRNAid<-dplyr::left_join(rosmapId,rosmap)
#remove duplicate rows
rosmapRNAid <- unique(rosmapRNAid)
rosmapRNAid2 <- subset(rosmapRNAid, select=c(rnaseq_id,braaksc,ceradsc))
names(rosmapRNAid2)[names(rosmapRNAid2) == "rnaseq_id"] <- "SampleID"

temp2<-dplyr::left_join(temp2,rosmapRNAid2, by="SampleID")

temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc <- factor(temp2$ceradsc,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx, levels = c(1:6))


gene_short_name <- inters$gene_short_name

rownames(temp)<-NULL
rownames(temp2)<-NULL

MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

x <- list()
x$rnaseqID <- MonRun$SampleID
x$State_rnaseq <- MonRun$State
x$Pseudotime_rnaseq <- MonRun$Pseudotime
F_inters_rnaseq <- as.data.frame(x)
F_inters_rnaseq$pseudotime_scRNA <- scale(F_inters_rnaseq$Pseudotime, center=F)


F_intersecting_pseudotimes <- merge(F_inters_prots, F_inters_rnaseq, by="rnaseqID")

cor(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA, method="pearson")
cor(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA, method="spearman")
plot(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA)
abline(lm(F_intersecting_pseudotimes$pseudotime_sc~F_intersecting_pseudotimes$pseudotime_scRNA))
lines(lowess(corrs$rnaseq_pseudotime_sc,corrs$pseudotime_sc))

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="Pseudotime", y="rnaseq_Pseudotime", 
                  add="reg.line",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()