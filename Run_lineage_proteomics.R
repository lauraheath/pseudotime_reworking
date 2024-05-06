#first run "TMTprot_DiffExpression.R" script and save the protein matrix and metadata file in your directory
library(monocle)
install.packages('DescTools')

setwd('~/pseudotime_reworking/')

#upload batch-corrected, median-abundance-centered log2-transformed matrix:
#p <- synapser::synGet('syn21266454')
#Log2_Normalized <- read.csv(p$path)
#upload batch-corrected, median-abundance-centered log2-transformed matrix with next batch of TMT proteins (N = 604),
#saved from TMTprot_DiffExp_forAgora.R so duplicate samples removed
#p <- synapser::synGet('syn28723003')
Log2_Normalized <- read.csv(file='data_objects/Log2_Normalized.csv')

rownames(Log2_Normalized) <- Log2_Normalized$X
Log2_Normalized$X<-NULL

#winsorize to minimize extreme outliers
#for( i in 1:dim(Log2_Normalized)[1] ){
#  Log2_Normalized[i,] <- DescTools::Winsorize( as.numeric(Log2_Normalized[i,]), na.rm = TRUE ) 
#}


#p2 <- synapser::synGet('syn28723027')
Meta <- read.csv(file='data_objects/TMT_metadata.csv')


# Harmonize case-control status
Meta$diagnosis <- "other"
Meta$diagnosis[Meta$cogdx == 1 & Meta$braaksc <= 3 & Meta$ceradsc_RADCnonStd >= 3] <- "control"
Meta$diagnosis[Meta$cogdx == 4 & Meta$braaksc >= 4 & Meta$ceradsc_RADCnonStd <= 2] <- "AD"
table(Meta$diagnosis, Meta$msex)


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
                         expressionFamily=gaussianff())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, norm_method='none')
  
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

Dat <- Log2_Normalized
Dat[is.na(Dat)] <- 0

#subset by sex: msex==0 is female, msex==1 is male; run separately for sex-specific analysis
In_S <- which(Meta$msex == 0)
#In_S <- which(Meta$msex == 1)
Dat2 <- Dat[,In_S]
Meta2 <- Meta[In_S,]


#one major outlier in female data, appears to have a very different expression distribution 
#compared to other patients: remove b22.130N from Meta2 and from matrix (this patient is female, diagnosis = 'other')
Meta3 <-subset(Meta2, batch.channel!='b22.130N')
In_D <- which(Meta2$batch.channel!='b22.130N')
Dat3 <- Dat2[,In_D]

#one major outlier in male data, appears to have a very different expression distribution (similar to female outlier)
#compared to other patients: remove b52.129C from Meta2 and from matrix (this patient is male, diagnosis = 'other')
Meta3 <-subset(Meta2, batch.channel!='b52.129C')
In_D <- which(Meta2$batch.channel!='b52.129C')
Dat3 <- Dat2[,In_D]

#save the matrices for DE by state analysis:
#saveRDS(Dat3, file="~/prot-lineage/data_objects/Female_prot_matrix.rds")
saveRDS(Dat3, file="~/prot-lineage/data_objects/Male_prot_matrix.rds")



#get list of proteins that are differentially expressed between AD case & control (DE analysis run for AGORA)
p <- synapser::synGet('syn50449887')
ADgenes_prot <- read.csv(p$path)
#for sex-adjusted DE proteins:
#ADgenes_prot1 <- subset(ADgenes_prot, ADgenes_prot$comparison=='AD_DLPFC - CT_DLPFC')
#ADgenes_prot1 <- subset(ADgenes_prot1, ADgenes_prot1$PVal<0.05)
#dim(ADgenes_prot1)


#female only diff exp proteins
females <- subset(ADgenes_prot, ADgenes_prot$comparison=='AD_female_DLPFC - CT_female_DLPFC')
ADgenes_prot1 <- subset(females, females$PVal<0.05)
dim(ADgenes_prot1)

males <- subset(ADgenes_prot, ADgenes_prot$comparison=='AD_male_DLPFC - CT_male_DLPFC')
ADgenes_prot1 <- subset(males, males$PVal<0.05)
dim(ADgenes_prot1)

# males <- subset(ADgenes_prot, ADgenes_prot$comparison=='AD_male_DLPFC - CT_male_DLPFC')
# females <- subset(females, females$PVal<0.05)
# males <- subset(males, males$PVal<0.05)
# inters <- intersect(females$GeneName,males$GeneName)
#N = 446 that are in both males and females (less than 50% shared)

#inters2 <- intersect(DEFemaleGenes$hgnc_symbol,ADgenes_prot2$GeneName)

#subset the protein matrix
genes2<-c()
for (gene in unique(c(as.vector(ADgenes_prot1$UniqID)))){
  if (gene %in% rownames(Dat3)){
    genes2 <- c(genes2,which(rownames(Dat3)==gene))
  }
}
length(genes2)
Dat3 <- Dat3[genes2,]
dim(Dat3)

temp <- Dat3
temp2 <- Meta3

gene_short_name <- rownames(Dat3)

temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc_RADCnonStd <- factor(temp2$ceradsc_RADCnonStd,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx,levels = c(1:6))


rownames(temp)<-NULL
rownames(temp2)<-NULL


#Run Monocle2: (ignore warning messages that occur)
#need to reverse order for male samples
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g



#check to see if State designations make sense
g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g


table(MonRun$State, MonRun$diagnosis)

#for female tree, reorder with state 2 as the root state, which has the greatest concentration of control pts & fewer AD pts
MonRun <- orderCells(MonRun, root_state = 2)

plot_cell_trajectory(MonRun,color_by = "Pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)

table(MonRun$State)
#collapse some states due to small numbers and rearrange
plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)

#Relabel the states so they are in order along the trajectory, and collapse small states
#female samples
# MonRun$State2 <- MonRun$State
# MonRun$State2[MonRun$State == 1] <- 1
# MonRun$State2[MonRun$State == 2] <- 2
# MonRun$State2[MonRun$State == 3] <- 3
# MonRun$State2[MonRun$State == 4] <- 3
# MonRun$State2[MonRun$State == 5] <- 4
# MonRun$State2[MonRun$State == 11] <- 4
# MonRun$State2[MonRun$State == 6] <- 5
# MonRun$State2[MonRun$State == 7] <- 5
# MonRun$State2[MonRun$State == 8] <- 5
# MonRun$State2[MonRun$State == 10] <- 6
# MonRun$State2[MonRun$State == 9] <- 7

MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 2] <- 1
MonRun$State2[MonRun$State == 1] <- 2
MonRun$State2[MonRun$State == 3] <- 3
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 11] <- 4
MonRun$State2[MonRun$State == 6] <- 5
MonRun$State2[MonRun$State == 7] <- 5
MonRun$State2[MonRun$State == 8] <- 5
MonRun$State2[MonRun$State == 10] <- 6
MonRun$State2[MonRun$State == 9] <- 7

#male samples
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 1] <- 1
MonRun$State2[MonRun$State == 11] <- 2
MonRun$State2[MonRun$State == 2] <- 2
MonRun$State2[MonRun$State == 10] <- 3
MonRun$State2[MonRun$State == 3] <- 3
MonRun$State2[MonRun$State == 4] <- 4
MonRun$State2[MonRun$State == 5] <- 5
MonRun$State2[MonRun$State == 6] <- 5
MonRun$State2[MonRun$State == 7] <- 5
MonRun$State2[MonRun$State == 8] <- 6
MonRun$State2[MonRun$State == 9] <- 7



MonRun$State2 <- as.character(MonRun$State2)
table(MonRun$State2)

plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)

#save Monocle object for later
saveRDS(MonRun, file='data_objects/MonRun_prot_female.RDS')
file <- synapser::File(path='data_objects/MonRun_prot_female.RDS', parentId='syn44292253')
file <- synapser::synStore(file)
saveRDS(MonRun, file='data_objects/MonRun_prot_male.RDS')
file <- synapser::File(path='data_objects/MonRun_prot_male.RDS', parentId='syn44292253')
file <- synapser::synStore(file)

MonRun$State2 <- as.numeric(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g



###### FIGURES #########
tiff(file='~/prot-lineage/figures/FEMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()

tiff(file='~/prot-lineage/figures/FEMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)



tiff(file='~/prot-lineage/figures/FEMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()

MonRun$ceradsc_RADCnonStd <- fct_rev(MonRun$ceradsc_RADCnonStd)
tiff(file='~/prot-lineage/figures/FEMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc_RADCnonStd, y=scale(Pseudotime,center=F),fill=ceradsc_RADCnonStd)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g
dev.off()


tiff(file='~/prot-lineage/figures/FEMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Cognitive\nDiagnosis",y="Pseudotime",x="Cognitive Diagnosis")
g
dev.off()





######## for stats & other figures, create a dataframe with all relevant covariates & pseudotime & state

x <- list()
x$msex <- MonRun$msex
x$SampleID <- MonRun$batch.channel
x$projid <- MonRun$projid.ROSMAP
x$State <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$diagnosis <- MonRun$diagnosis
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc_RADCnonStd
x$cogdx <- MonRun$cogdx
x$apoe_genotype <- MonRun$apoe_genotype
x$educ   <- MonRun$educ
x$pmi <- MonRun$pmi
x$age_death <- MonRun$age_death

#rename and create a scaled pseudotime variable
prot_pstime_covars_F <- as.data.frame(x)
prot_pstime_covars_F$pseudotime_sc <- scale(prot_pstime_covars_F$Pseudotime, center=F)

#save variables file for later
write.csv(prot_pstime_covars_F, file="~/prot-lineage/data_objects/female_pseudotimes_states.csv", row.names=FALSE)
file <- synapser::File(path='~/prot-lineage/data_objects/female_pseudotimes_states.csv', parentId='syn25607662')
file <- synapser::synStore(file)

write.csv(prot_pstime_covars_F, file="~/prot-lineage/data_objects/male_pseudotimes_state.csv", row.names=FALSE)
file <- synapser::File(path='~/prot-lineage/data_objects/male_pseudotimes_state.csv', parentId='syn25607662')
file <- synapser::synStore(file)

#run female and male analyses separately
Fvariables <- prot_pstime_covars_F
#Fvariables <- prot_pstime_covars_M

#run logistic regression comparing pseudotiem between cases and controls only
casecontrolF <- subset(Fvariables, Fvariables$diagnosis=='AD'|Fvariables$diagnosis=='control')
casecontrolF$diag2 <- ifelse(casecontrolF$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,casecontrolF,family='binomial'))

tiff(file='~/prot-lineage/figures/FEMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(casecontrolF,aes(x=diagnosis,
                             y=pseudotime_sc,
                             color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g
dev.off()

#run proportional odds logistic regression for neuropath/cognitive endpoints:
braakfit <- MASS::polr(braaksc ~ pseudotime_sc,Fvariables)
ceradfit <- MASS::polr(ceradsc ~ pseudotime_sc,Fvariables)
cogdxfit <- MASS::polr(cogdx ~ pseudotime_sc,Fvariables)

cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')



#look for correlations with GWAS LOAD genes
ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")

Dat3$gene_names <- rownames(Dat3)
Dat3$gene_short_name <- gsub("\\|.*", "", Dat3$gene_names)


# dlpfcCPMObj <- synapser::synGet('syn8456638')
# Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat3)#[-120]
#sampleIds <- sampleIds[-120]
sampleIds
geneIds <- Dat3$gene_short_name
Dat3$gene_short_name<-NULL
Dat3$gene_names<-NULL
Dat3 <- t(Dat3)
colnames(Dat3) <- geneIds
Dat3 <- data.frame(Dat3,stringsAsFactors=F)
Dat3$sampleId <- sampleIds
dlpfc <- dplyr::left_join(Fvariables,Dat3,by=c('SampleID'='sampleId'))
dlpfc2 <- dlpfc[,14:1850]

corvec <- cor(dlpfc2,dlpfc$Pseudotime,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$geneid %in% ad_gwas,]$cor))

mean(corDfdlpfc$cor)
mean(corDfdlpfc[corDfdlpfc$geneid %in% ad_gwas,]$cor)


corDfdlpfc$adGwas <- corDfdlpfc$geneid %in% ad_gwas
colnames(corDfdlpfc)[3] <- 'LOADGWASGene'
corDfdlpfc$LOADGWASGene2 <- ifelse(corDfdlpfc$LOADGWASGene==FALSE, "NOT GWAS GENE", "GWAS GENE")

#tiff(file='~/prot-lineage/figures/FEMALE_loadgwas_cor.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_loadgwas_cor.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(corDfdlpfc,ggplot2::aes(x=LOADGWASGene2,y=cor,fill=LOADGWASGene2))
g <- g + ggplot2::geom_boxplot() + theme(legend.position="none")
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'LOAD GWAS STATUS',y='Correlation with pseudotime')
g
dev.off()