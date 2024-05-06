library(monocle)
setwd("~/pseudotime_reworking")

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
  
  HSMM <- orderCells(HSMM, reverse=TRUE)
  #HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}


#upload spreadsheet with sex-specific matrices and metadata in an excel workbook (provided by Richa)
#female data sheet:
metsFid <- "syn48321418"
synapser::synGet(metsFid, downloadLocation = "files/")
path <- "files/rosmap_metabolomics_processed_med_covar_corrected_females.xlsx"
sheetnames <- excel_sheets(path)
metsF <- lapply(excel_sheets(path), read_excel, path = path)
names(metsF) <- sheetnames
D1 <- metsF$assay
D1 <- as.data.frame(D1)
rownames(D1) <- D1$...1
D1$...1<-NULL
metabolite_metadata <- metsF$rowData
mets_dx <- metsF$colData


#male data sheet:
metsFid <- "syn48321423"
synapser::synGet(metsFid, downloadLocation = "files/")
path <- "files/rosmap_metabolomics_processed_med_covar_corrected_males.xlsx"
sheetnames <- excel_sheets(path)
metsM <- lapply(excel_sheets(path), read_excel, path = path)
names(metsM) <- sheetnames
D1 <- metsM$assay
D1 <- as.data.frame(D1)
rownames(D1) <- D1$...1
D1$...1<-NULL
metabolite_metadata <- metsM$rowData
mets_dx <- metsM$colData

#save the matrix to local directory:
saveRDS(D1, file='metsF_matrix_all.rds')
saveRDS(D1, file='metsM_matrix_all.rds')


#code diagnosis category (AD vs CT already coded by Risha)
mets_dx$diagnosis <- 'OTHER'
mets_dx$diagnosis[mets_dx$sage_diag==0] <- 'CT'
mets_dx$diagnosis[mets_dx$sage_diag==1] <-'AD'
table(mets_dx$diagnosis)

#upload the differentially abundant metabolites by sage diagnosis:
#female DA metabolites:
p3 <- synapser::synGet('syn50912014')
#male DA metabolites
p3 <- synapser::synGet('syn50912015')

DE_abundance <- read.csv(p3$path)
DE_abundance1 <- subset(DE_abundance, DE_abundance$p.value<0.05)
#subset the metabolomics matrix
DEmets<-c()
for (met in unique(c(as.vector(DE_abundance1$var)))){
  if (met %in% rownames(D1)){
    DEmets <- c(DEmets,which(rownames(D1)==met))
  }
}
length(DEmets)
D2 <- D1[DEmets,]
dim(D2)

#hsmm_sup <-  get_trajectories(D2)
gene_short_name <- rownames(D2)

temp <- D2
temp2 <- mets_dx


rownames(temp)<-NULL
rownames(temp2)<-NULL


#Run Monocle2: (ignore warning messages that occur)
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

#for males, state 9 has the obvious highest concentration of control samples. reroot at state 9:
MonRun <- orderCells(MonRun, root_state = 9)
plot_cell_trajectory(MonRun,color_by = "Pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1)

g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

table(MonRun$State)
plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)

#recode States so they are in order and collapse very small states
#female samples
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 2] <- 2
MonRun$State2[MonRun$State == 4] <- 2
MonRun$State2[MonRun$State == 3] <- 3
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 7] <- 4
MonRun$State2[MonRun$State == 13] <- 5
MonRun$State2[MonRun$State == 8] <- 6
MonRun$State2[MonRun$State == 9] <- 6
MonRun$State2[MonRun$State == 11] <- 6
MonRun$State2[MonRun$State == 10] <- 7
MonRun$State2[MonRun$State == 12] <- 8

#male samples
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 9] <- 1
MonRun$State2[MonRun$State == 1] <- 2
MonRun$State2[MonRun$State == 2] <- 3
MonRun$State2[MonRun$State == 8] <- 3
MonRun$State2[MonRun$State == 4] <- 4
MonRun$State2[MonRun$State == 7] <- 4
MonRun$State2[MonRun$State == 6] <- 5
MonRun$State2[MonRun$State == 5] <- 6


MonRun$State2 <- as.character(MonRun$State2)
table(MonRun$State2)
plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1)

#save Monocle object for later
saveRDS(MonRun, file='MonRun_mets_female.RDS')
file <- synapser::File(path='MonRun_mets_female.RDS', parentId='syn51204295')
file <- synapser::synStore(file)

saveRDS(MonRun, file='MonRun_mets_male.RDS')
file <- synapser::File(path='MonRun_mets_male.RDS', parentId='syn51204295')
file <- synapser::synStore(file)

#tiff(file='figures/FEMALE_mets_tree_state.tiff',height=85,width=100,units='mm',res=300)
tiff(file='figures/MALE_mets_tree_state.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()

#tiff(file='figures/FEMALE_mets_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
tiff(file='figures/MALE_mets_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
dev.off()


MonRun$braaksc <- as.factor(MonRun$braaksc)
#tiff(file='figures/FEMALE_mets_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='figures/MALE_mets_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()

MonRun$ceradsc <- as.factor(MonRun$ceradsc)
MonRun$ceradsc <- fct_rev(MonRun$cerads)
#tiff(file='figures/FEMALE_mets_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='figures/MALE_mets_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=scale(Pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g
dev.off()

MonRun$cogdx <- as.factor(MonRun$cogdx)
#tiff(file='figures/FEMALE_mets_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='figures/MALE_mets_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Cognitive\nDiagnosis",y="Pseudotime",x="Cognitive Diagnosis")
g
dev.off()


x <- list()
x$projid <- MonRun$projid
x$State <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$diagnosis <- MonRun$diagnosis
x$msex <- MonRun$msex
x$educ <- MonRun$educ
x$apoe_genotype <- MonRun$apoe_genotype
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$pmi <- MonRun$pmi


#rename and create a scaled pseudotime variable
pseudo <- as.data.frame(x)
pseudo$pseudotime_sc <- scale(pseudo$Pseudotime, center=F)

#save variables file for later
write.csv(pseudo, file="data_objects/FEMALE_mets_pseudotimes_states_covars.csv", row.names=FALSE)
file <- synapser::File(path='data_objects/FEMALE_mets_pseudotimes_states_covars.csv', parentId='syn45147359')
file <- synapser::synStore(file)

write.csv(pseudo, file="data_objects/MALE_mets_pseudotimes_states_covars.csv", row.names=FALSE)
file <- synapser::File(path='data_objects/MALE_mets_pseudotimes_states_covars.csv', parentId='syn45147359')
file <- synapser::synStore(file)

#run logistic regression comparing pseudotime between cases and controls only
casecontrol <- subset(pseudo, pseudo$diagnosis=='AD'|pseudo$diagnosis=='CT')
casecontrol$diag2 <- ifelse(casecontrol$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,casecontrol,family='binomial'))

#tiff(file='figures/FEMALE_mets_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
tiff(file='figures/MALE_mets_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
g <- ggplot(casecontrol,aes(x=diagnosis,
                            y=pseudotime_sc,
                            color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g
dev.off()