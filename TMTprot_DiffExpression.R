#get differentially abundant proteins from AMP-AD 2.0 TMT proteomics matrix

library(limma)
library(biomaRt)
library(stringr)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

## Import batch-corrected, centered, Log2-transformed protein expression matrix
#p <- synapser::synGet('syn21266454')
#NEW PHASE 1 and 2 TMT PROTEOMICS MERGED: use minimally regressed, batch- and set-corrected data, log2 transformed (dataset 2a)
p <- synapser::synGet('syn28723003')
Log2_Normalized <- read.csv(p$path)

rownames(Log2_Normalized) <- Log2_Normalized$X
Log2_Normalized$X<-NULL

#Fix protein names (some out of date)
# p1 <- synapser::synGet('syn24216770')
# correct_geneIDs <- read.csv(p1$path)
# 
# names(Log2_Normalized)[names(Log2_Normalized) == 'X'] <- 'OldPeptideID'
# Log2_Normalized <- dplyr::left_join(Log2_Normalized, correct_geneIDs, by="OldPeptideID")
# rownames(Log2_Normalized) <- Log2_Normalized$NewPeptideID
# Log2_Normalized$OldPeptideID<-NULL
# Log2_Normalized$NewPeptideID<-NULL
# Log2_Normalized$Old_Gene<-NULL
# Log2_Normalized$Old_Pep<-NULL
# Log2_Normalized$New_Gene<-NULL
# Log2_Normalized$New_Pep<-NULL
# Log2_Normalized$ENSG<-NULL

#get sample metadata
#p2 <- synapser::synGet('syn21323404')
#Meta <- read.csv(p2$path)

p2 <- synapser::synGet('syn28723027')
Meta <- read.csv(p2$path)


# Harmonize case-control status
Meta$diagnosis <- "other"
Meta$diagnosis[Meta$cogdx == 1 & Meta$braaksc <= 3 & Meta$ceradsc_RADCnonStd >= 3] <- "control"
Meta$diagnosis[Meta$cogdx == 4 & Meta$braaksc >= 4 & Meta$ceradsc_RADCnonStd <= 2] <- "AD"
table(Meta$diagnosis, Meta$msex)

#There are duplicate individualIDs in this data (2 samples run on 6 individuals). Only want one sample per patient
#look at the duplicate samples metadata
n_occur <- data.frame(table(Meta$projid.ROSMAP))
names(n_occur)[names(n_occur) == 'Var1'] <- 'projid.ROSMAP'
metadata2 <- merge(Meta, n_occur)
metadata2 <- subset(metadata2, metadata2$Freq==2)

#specimen b19.130N, b49.127C, b43.127C, b02.127C, and b21.130C are outliers in e. dammer's connectivity analysis 
#also delete specimen b57.127N
rownames(Meta) <- Meta$batch.channel
Meta3 <- Meta[!(row.names(Meta) %in% c('b19.130N','b49.127C','b43.127C','b02.127C','b21.130C', 'b57.127N')),]
length(unique(Meta3$projid.ROSMAP))





####Metadata file for AD case-control Diagnosis only for analysis
Meta_D <- Meta[ Meta$diagnosis %in% c('AD','control'), ]
#Code Sex and Diagnosis as Factors
Meta$msex <- as.factor(Meta$msex)
Meta_D$diagnosis <- as.factor(Meta_D$diagnosis)
#pare log2 matrix to AD case and controls
ADmatrix <- Log2_Normalized[Meta_D$batch.channel]

#run limma, adjust for sex and pmi 
design <- model.matrix(~ 0 + diagnosis + msex + pmi + ApoE4.Dose, data = Meta_D)
design
cont.matrix <- makeContrasts(diagnosisAD-diagnosiscontrol, levels=colnames(design))
cont.matrix

fit <- lmFit(ADmatrix, design)
fit <- contrasts.fit(fit, cont.matrix)
#do robust regression to minimize effect of outlier genes
fit <- eBayes(fit, robust=TRUE)
topTable(fit)

#extract results for all genes and get 95% confidence limits
allresults <- topTable(fit, n=Inf, confint=0.95)

allresults$Peptide <- rownames(allresults)
allresults[, 10:11] <- str_split_fixed(allresults$Peptide, fixed("|"), 2)
names(allresults)[names(allresults) == 'V10'] <- 'GeneID'
names(allresults)[names(allresults) == 'V11'] <- 'UniprotID'

mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
univals <- allresults$UniprotID
ens_trans <- getBM(c("ensembl_gene_id", "uniprot_gn_id", "hgnc_symbol"), "uniprot_gn_id", univals, mart)

#In many cases, there are multiple ensembl ids per uniprotID select the most relevant ensemblID

dupes <- ens_trans[duplicated(ens_trans$uniprot_gn_id),]
ens_trans2 <- ens_trans
#the last entry for each duplicate is the correct one (according to uniprot)--keep the last entry only for duplicate entries
ens_trans2 <- aggregate(. ~ uniprot_gn_id, data=ens_trans2, FUN = tail, 1)
#rejoin to allresults
allresults2 <- allresults
names(ens_trans2)[names(ens_trans2) == 'uniprot_gn_id'] <- 'UniprotID'
allresults2 <- dplyr::left_join(allresults2, ens_trans2)

#several uniprot ids did not map to ensemblIDs. rerun biomart on those rows only, using hgnc symbol to return ensemblID
ens_missing <- subset(allresults2, is.na(allresults2$ensembl_gene_id))
genesymbols <- ens_missing$GeneID
ens_missing2 <- getBM(c("ensembl_gene_id", "hgnc_symbol"), "hgnc_symbol", genesymbols, mart)


# still 16 duplicates. look up these genes, decide which ensembl id is appropriate based on uniprot database, ensembl database
dupes <- ens_missing2[duplicated(ens_missing2$hgnc_symbol),]

#join back to allresults2
#drop the duplicates from ens_missing
ens_missing3 <- ens_missing2 %>%
  group_by(hgnc_symbol) %>%
  filter(n() == 1) %>%
  ungroup()


names(ens_missing3)[names(ens_missing3) == 'ensembl_gene_id'] <- 'ensembl_gene_id_new'
names(ens_missing3)[names(ens_missing3) == 'hgnc_symbol'] <- 'GeneID'

allresults3 <- dplyr::left_join(allresults2, ens_missing3, by = "GeneID")
#allresults3 <- merge(allresults2, ens_missing3)
#allresults3 <- allresults2
#allresults3$ensembl_gene_id[is.na(allresults3$ensembl_gene_id)] <- ens_missing3$ensembl_gene_id[match(allresults3$hgnc_symbol,ens_missing3$hgnc_symbol)][which(is.na(allresults3$ensembl_gene_id))]

allresults3$ensembl_gene_id <- ifelse(is.na(allresults3$ensembl_gene_id), allresults3$ensembl_gene_id_new, allresults3$ensembl_gene_id)

#Some could not be reconciled with biomart. fill in appropriate ensembl id in allresults3
allresults3[ allresults2$GeneID=='C1R', ]$ensembl_gene_id <- 'ENSG00000159403'
allresults3[ allresults2$GeneID=='CACNA1C', ]$ensembl_gene_id <- 'ENSG00000151067'
allresults3[ allresults2$GeneID=='COG5', ]$ensembl_gene_id <- 'ENSG00000164597'
allresults3[ allresults2$GeneID=='EPHB6', ]$ensembl_gene_id <- 'ENSG00000106123'
allresults3[ allresults2$UniprotID=='P30456', ]$ensembl_gene_id <- 'ENSG00000206503'
allresults3[ allresults2$UniprotID=='P30455', ]$ensembl_gene_id <- 'ENSG00000206503'
allresults3[ allresults2$UniprotID=='P01891', ]$ensembl_gene_id <- 'ENSG00000206503'
allresults3[ allresults2$UniprotID=='P13746', ]$ensembl_gene_id <- 'ENSG00000206503'
allresults3[ allresults2$UniprotID=='P01892', ]$ensembl_gene_id <- 'ENSG00000206503'
allresults3[ allresults2$UniprotID=='P30460', ]$ensembl_gene_id <- 'ENSG00000234745'
allresults3[ allresults2$UniprotID=='P18463', ]$ensembl_gene_id <- 'ENSG00000234745'
allresults3[ allresults2$UniprotID=='P30483', ]$ensembl_gene_id <- 'ENSG00000234745'
allresults3[ allresults2$UniprotID=='P04229', ]$ensembl_gene_id <- 'ENSG00000196126'
allresults3[ allresults2$UniprotID=='P20039', ]$ensembl_gene_id <- 'ENSG00000196126'
allresults3[ allresults2$GeneID=='HSP90AA4P', ]$ensembl_gene_id <- 'ENSG00000205100'
allresults3[ allresults2$GeneID=='IGHG3', ]$ensembl_gene_id <- 'ENSG00000211897'
allresults3[ allresults2$GeneID=='IGKV2-40', ]$ensembl_gene_id <- 'ENSG00000273962'
allresults3[ allresults2$GeneID=='INTS3', ]$ensembl_gene_id <- 'ENSG00000143624'
allresults3[ allresults2$GeneID=='MAGI1', ]$ensembl_gene_id <- 'ENSG00000151276'
allresults3[ allresults2$GeneID=='POLR2A', ]$ensembl_gene_id <- 'ENSG00000181222'
allresults3[ allresults2$GeneID=='SAA1', ]$ensembl_gene_id <- 'ENSG00000173432'
allresults3[ allresults2$GeneID=='VPS11', ]$ensembl_gene_id <- 'ENSG00000160695'
allresults3[ allresults2$GeneID=='WASH2P', ]$ensembl_gene_id <- 'ENSG00000291134'
#allresults3[ allresults2$UniprotID=='P0DN79', ]$ensembl_gene_id <- '' #no ensembl id for this gene
allresults3[ allresults2$UniprotID=='Q5VUR7', ]$ensembl_gene_id <- 'ENSG00000276203'
allresults3[ allresults2$UniprotID=='Q12912', ]$ensembl_gene_id <- 'ENSG00000118308'
#allresults3[ allresults2$UniprotID=='P0DP04', ]$ensembl_gene_id <- '' #no ensembl id for this gene
allresults3[ allresults2$UniprotID=='Q9HB07', ]$ensembl_gene_id <- 'ENSG00000139637'
allresults3[ allresults2$UniprotID=='Q6DN03', ]$ensembl_gene_id <- 'ENSG00000261716'
allresults3[ allresults2$UniprotID=='A0A0B4J2D5', ]$ensembl_gene_id <- 'ENSG00000280071'

allresults3$hgnc_symbol<-NULL
allresults3$ensembl_gene_id_new<-NULL
sum(is.na(allresults3$ensembl_gene_id))

#save the translations to a file:
gene_uniprot_ensembl_key <- subset(allresults3, select=c(Peptide, GeneID, UniprotID, ensembl_gene_id))
write.csv(gene_uniprot_ensembl_key, file="~/prot-lineage/data_objects/gene_uniprot_ensembl_key.csv", row.names=FALSE)


#relabel columns to match lfq data for agora

names(allresults3)[names(allresults3) == 'Peptide'] <- 'UniqID'
names(allresults3)[names(allresults3) == 'GeneID'] <- 'GeneName'
names(allresults3)[names(allresults3) == 'UniprotID'] <- 'UniProtID'
names(allresults3)[names(allresults3) == 'ensembl_gene_id'] <- 'ENSG'

allresults3$Tissue <- 'DLPFC'
names(allresults3)[names(allresults3) == 'logFC'] <- 'Log2_FC'
names(allresults3)[names(allresults3) == 'CI.L'] <- 'CI_Lwr'
names(allresults3)[names(allresults3) == 'CI.R'] <- 'CI_Upr'
names(allresults3)[names(allresults3) == 'P.Value'] <- 'PVal'
names(allresults3)[names(allresults3) == 'adj.P.Val'] <- 'Cor_PVal'

allresults3$comparison <- 'AD_DLPFC - CT_DLPFC'

##reorder the columns
col_order <- c("comparison", "UniqID", "GeneName", "UniProtID", "ENSG", "Tissue", "Log2_FC", "CI_Upr", "CI_Lwr", "PVal", "Cor_PVal")
allresults3 <- allresults3[, col_order]


#save to synapse working group space:
#write.csv(allresults3, file="~/prot-lineage/data_objects/ROSMAP_DiffExp_TMTproteins_.csv", row.names=FALSE)
#file <- synapser::File(path='~/prot-lineage/data_objects/ROSMAP_DiffExp_TMTproteins.csv', parentId='syn25607662')
#file <- synapser::synStore(file)

#save data objects to local folder for downstream pseudotime analyses
#prune original log2 matrix to include only the samples in Meta3 (to avoid duplicate projids)
ADmatrix2 <- Log2_Normalized[Meta3$batch.channel]
write.csv(ADmatrix2, file="~/prot-lineage/data_objects/Log2_Normalized.csv")
write.csv(Meta3, file="~/prot-lineage/data_objects/TMT_metadata.csv")


#old DE proteins
#p <- synapser::synGet('syn35221005')
#ADgenes_prot <- read.csv(p$path)



################# run each sex separately and join all results ############

#differential protein expression in females
MetaF <- Meta_D[(Meta_D$msex == 0),]
F_matrix <- Log2_Normalized[MetaF$batch.channel]

#run limma, adjust for sex and pmi 
design <- model.matrix(~ 0 + diagnosis + pmi + ApoE4.Dose, data = MetaF)
design
cont.matrix <- makeContrasts(diagnosisAD-diagnosiscontrol, levels=colnames(design))
cont.matrix

fit <- lmFit(F_matrix, design)
fit <- contrasts.fit(fit, cont.matrix)
#do robust regression to minimize effect of outlier genes
fit <- eBayes(fit, robust=TRUE)
topTable(fit)

#extract results for all genes and get 95% confidence limits
allresults <- topTable(fit, n=Inf, confint=0.95)

#relabel columns to match lfq data for agora, and get ENSG names
genes_key <- read.csv(file='~/prot-lineage/data_objects/gene_uniprot_ensembl_key.csv')
allresults$Peptide <- row.names(allresults)
allresultsF <- left_join(allresults, genes_key)

names(allresultsF)[names(allresultsF) == 'Peptide'] <- 'UniqID'
names(allresultsF)[names(allresultsF) == 'GeneID'] <- 'GeneName'
names(allresultsF)[names(allresultsF) == 'UniprotID'] <- 'UniProtID'
names(allresultsF)[names(allresultsF) == 'ensembl_gene_id'] <- 'ENSG'


allresultsF$Tissue <- 'DLPFC'
allresultsF$comparison <- 'AD_female_DLPFC - CT_female_DLPFC'
names(allresultsF)[names(allresultsF) == 'logFC'] <- 'Log2_FC'
names(allresultsF)[names(allresultsF) == 'CI.L'] <- 'CI_Lwr'
names(allresultsF)[names(allresultsF) == 'CI.R'] <- 'CI_Upr'
names(allresultsF)[names(allresultsF) == 'P.Value'] <- 'PVal'
names(allresultsF)[names(allresultsF) == 'adj.P.Val'] <- 'Cor_PVal'

##reorder the columns
col_order <- c("comparison", "UniqID", "GeneName", "UniProtID", "ENSG", "Tissue", "Log2_FC", "CI_Upr", "CI_Lwr", "PVal", "Cor_PVal")
allresultsF <- allresultsF[, col_order]

#save to synapse
#write.csv(allresults2, file="prot-lineage/data_objects/female_DEproteins.csv", row.names=FALSE)
#save to synapse
#file <- synapser::File(path='~/prot-lineage/data_objects/ROSMAP_DiffExp_TMTproteins.csv', parentId='syn35219190')
#file <- synapser::synStore(file)



#differential protein expression in males
MetaM <- Meta_D[(Meta_D$msex == 1),]
M_matrix <- Log2_Normalized[MetaM$batch.channel]

#run limma, adjust for sex and pmi 
design <- model.matrix(~ 0 + diagnosis + pmi + ApoE4.Dose, data = MetaM)
design
cont.matrix <- makeContrasts(diagnosisAD-diagnosiscontrol, levels=colnames(design))
cont.matrix

fit <- lmFit(M_matrix, design)
fit <- contrasts.fit(fit, cont.matrix)
#do robust regression to minimize effect of outlier genes
fit <- eBayes(fit, robust=TRUE)
topTable(fit)

#extract results for all genes and get 95% confidence limits
allresults <- topTable(fit, n=Inf, confint=0.95)

#relabel columns to match lfq data for agora, and get ENSG names
#genes_key <- read.csv(file='~/prot-lineage/data_objects/gene_uniprot_ensembl_key.csv')
allresults$Peptide <- row.names(allresults)
allresultsM <- left_join(allresults, genes_key)

names(allresultsM)[names(allresultsM) == 'Peptide'] <- 'UniqID'
names(allresultsM)[names(allresultsM) == 'GeneID'] <- 'GeneName'
names(allresultsM)[names(allresultsM) == 'UniprotID'] <- 'UniProtID'
names(allresultsM)[names(allresultsM) == 'ensembl_gene_id'] <- 'ENSG'


allresultsM$Tissue <- 'DLPFC'
allresultsM$comparison <- 'AD_male_DLPFC - CT_male_DLPFC'
names(allresultsM)[names(allresultsM) == 'logFC'] <- 'Log2_FC'
names(allresultsM)[names(allresultsM) == 'CI.L'] <- 'CI_Lwr'
names(allresultsM)[names(allresultsM) == 'CI.R'] <- 'CI_Upr'
names(allresultsM)[names(allresultsM) == 'P.Value'] <- 'PVal'
names(allresultsM)[names(allresultsM) == 'adj.P.Val'] <- 'Cor_PVal'

##reorder the columns
col_order <- c("comparison", "UniqID", "GeneName", "UniProtID", "ENSG", "Tissue", "Log2_FC", "CI_Upr", "CI_Lwr", "PVal", "Cor_PVal")
allresultsM <- allresultsM[, col_order]


#join female and male dataframes
allresults4 <- rbind(allresults3, allresultsF, allresultsM)


#save to synapse
write.csv(allresults4, file="~/prot-lineage/data_objects/DiffExprProteins_bysex_and_combined.csv", row.names=FALSE)
file <- synapser::File(path='~/prot-lineage/data_objects/DiffExprProteins_bysex_and_combined.csv', parentId='syn25607662')
file <- synapser::synStore(file)