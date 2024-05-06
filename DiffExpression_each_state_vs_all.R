######### RNAseq branch-specific differential expression analysis

#need Monrun object from prot_lineage_monocle_rerun.R and gene_short_name vector
#script for males following script for females (6 states in males, 7 in females)
#upload entire matrix first to look at associations with ALL genes, not just DE genes included for trajectory model

#MonRun <- readRDS(file="data_objects/MonRun_RNAseq_female.RDS")
#upload pseudotimes to get states
p <- synapser::synGet('syn38354143')
states <- read.csv(p$path)
#upload full counts matrix with all genes
#temp <- readRDS(file="data_objects/Female_fulldatamatrix.RDS")
p1 <- synapser::synGet('syn52420490')
temp <- readRDS(p1$path)
gene_short_name <- rownames(temp)
table(states$State)
#there are 7 states in the female tree


#want to run each state against all other states
temp2 <- as.data.frame(t(temp))
temp2$SampleID <- rownames(temp2)
states2 <- subset(states, select=c(SampleID, State))
temp3 <- dplyr::left_join(temp2, states2, by='SampleID')


#STATE 1 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==1, 1, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_1 <- rep(0,length(gene_short_name))
l2$d_1 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_1[i] <- tk$s[1,4]
  l2$d_1[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state1a <- dplyr::left_join(dfa1,dfb1)


#STATE 2 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==2, 2, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$d_2 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$d_2[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state2 <- dplyr::left_join(dfa1,dfb1)




#STATE 3 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==3, 3, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_3 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_3[i] <- tk$s[1,4]
  l2$d_3[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state3 <- dplyr::left_join(dfa1,dfb1)




#STATE 4 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==4, 4, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_4 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_4[i] <- tk$s[1,4]
  l2$d_4[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state4 <- dplyr::left_join(dfa1,dfb1)




#STATE 5 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==5, 5, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_5 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_5[i] <- tk$s[1,4]
  l2$d_5[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state5 <- dplyr::left_join(dfa1,dfb1)






#STATE 6 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==6, 6, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_6 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_6[i] <- tk$s[1,4]
  l2$d_6[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state6 <- dplyr::left_join(dfa1,dfb1)






#STATE 7 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==7, 7, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_7 <- rep(0,length(gene_short_name))
l2$d_7 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_7[i] <- tk$s[1,4]
  l2$d_7[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state7 <- dplyr::left_join(dfa1,dfb1)

temp3$State <- as.character(temp3$State)
ggplot(temp3, aes(x=State, y=SERPINA5)) + geom_boxplot()



#join state-specific data frames:
allDE <- rbind(state1a, state2)
allDE <- rbind(allDE, state3)
allDE <- rbind(allDE, state4)
allDE <- rbind(allDE, state5)
allDE <- rbind(allDE, state6)
allDE <- rbind(allDE, state7)

#save file:
write.csv(allDE, file="~/redo_DE_analysis/RNAseqF_DE_states_vs_all.csv", row.names=FALSE)







#########################run Male monocle object with script below:###############################################
#MonRun <- readRDS(file="data_objects/MonRun_RNAseq_female.RDS")
p <- synapser::synGet('syn38348029')
states <- read.csv(p$path)
#temp <- readRDS(file="data_objects/Female_fulldatamatrix.RDS")
p1 <- synapser::synGet('syn52420488')
temp <- readRDS(p1$path)
gene_short_name <- rownames(temp)
table(states$State)
#there are 6 states in the male tree


#want to run each state against all other states
temp2 <- as.data.frame(t(temp))
temp2$SampleID <- rownames(temp2)
states2 <- subset(states, select=c(SampleID, State))
temp3 <- dplyr::left_join(temp2, states2, by='SampleID')


#STATE 1 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==1, 1, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_1 <- rep(0,length(gene_short_name))
l2$d_1 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_1[i] <- tk$s[1,4]
  l2$d_1[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state1a <- dplyr::left_join(dfa1,dfb1)


#STATE 2 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==2, 2, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$d_2 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$d_2[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state2 <- dplyr::left_join(dfa1,dfb1)




#STATE 3 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==3, 3, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_3 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_3[i] <- tk$s[1,4]
  l2$d_3[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state3 <- dplyr::left_join(dfa1,dfb1)




#STATE 4 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==4, 4, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_4 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_4[i] <- tk$s[1,4]
  l2$d_4[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state4 <- dplyr::left_join(dfa1,dfb1)




#STATE 5 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==5, 5, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_5 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_5[i] <- tk$s[1,4]
  l2$d_5[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state5 <- dplyr::left_join(dfa1,dfb1)






#STATE 6 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==6, 6, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_6 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_6[i] <- tk$s[1,4]
  l2$d_6[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state6 <- dplyr::left_join(dfa1,dfb1)


temp3$State <- as.character(temp3$State)
ggplot(temp3, aes(x=State, y=CRH)) + geom_boxplot()



#join state-specific data frames:
allDE <- rbind(state1a, state2)
allDE <- rbind(allDE, state3)
allDE <- rbind(allDE, state4)
allDE <- rbind(allDE, state5)
allDE <- rbind(allDE, state6)

#save file:
write.csv(allDE, file="~/redo_DE_analysis/RNAseqM_DE_states_vs_all.csv", row.names=FALSE)
























########## Proteomics branch-specific differential expression analysis

#need Monrun object from prot_lineage_monocle_rerun.R and gene_short_name vector
#script for males following script for females (6 states in males, 7 in females)
#upload entire matrix first to look at associations with ALL genes, not just DE genes included for trajectory model

#get states
p <- synapser::synGet('syn35317137')
states <- read.csv(p$path)
#temp <- readRDS(file="data_objects/Female_fulldatamatrix.RDS")
p1 <- synapser::synGet('syn52297953')
temp <- readRDS(p1$path)

temp$proteins <- rownames(temp)
temp[,417:418] <- str_split_fixed(temp$proteins, fixed("|"), 2)

names(temp)[names(temp) == 'V417'] <- 'gene_names'
temp$gene_names<-make.unique(temp$gene_names)
rownames(temp) <- temp$gene_names
temp$proteins<-NULL
temp$gene_names<-NULL
temp$V418<-NULL


gene_short_name <- rownames(temp)
table(states$State)
#there are 7 states in the female tree


#want to run each state against all other states
temp2 <- as.data.frame(t(temp))
temp2$SampleID <- rownames(temp2)
states2 <- subset(states, select=c(SampleID, State))
temp3 <- dplyr::left_join(temp2, states2, by='SampleID')


#STATE 1 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==1, 1, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_1 <- rep(0,length(gene_short_name))
l2$d_1 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_1[i] <- tk$s[1,4]
  l2$d_1[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state1a <- dplyr::left_join(dfa1,dfb1)


#STATE 2 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==2, 2, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$d_2 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$d_2[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state2 <- dplyr::left_join(dfa1,dfb1)




#STATE 3 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==3, 3, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_3 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_3[i] <- tk$s[1,4]
  l2$d_3[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state3 <- dplyr::left_join(dfa1,dfb1)




#STATE 4 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==4, 4, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_4 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_4[i] <- tk$s[1,4]
  l2$d_4[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state4 <- dplyr::left_join(dfa1,dfb1)




#STATE 5 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==5, 5, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_5 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_5[i] <- tk$s[1,4]
  l2$d_5[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state5 <- dplyr::left_join(dfa1,dfb1)






#STATE 6 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==6, 6, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_6 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_6[i] <- tk$s[1,4]
  l2$d_6[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state6 <- dplyr::left_join(dfa1,dfb1)






#STATE 7 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==7, 7, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_7 <- rep(0,length(gene_short_name))
l2$d_7 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_7[i] <- tk$s[1,4]
  l2$d_7[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state7 <- dplyr::left_join(dfa1,dfb1)



#join state-specific data frames:
allDE <- rbind(state1a, state2)
allDE <- rbind(allDE, state3)
allDE <- rbind(allDE, state4)
allDE <- rbind(allDE, state5)
allDE <- rbind(allDE, state6)
allDE <- rbind(allDE, state7)

#save file:
write.csv(allDE, file="~/redo_DE_analysis/protF_DE_states_vs_all.csv", row.names=FALSE)







#########################run Male monocle object with script below:###############################################
#get states
p <- synapser::synGet('syn35317790')
states <- read.csv(p$path)
#temp <- readRDS(file="data_objects/Female_fulldatamatrix.RDS")
p1 <- synapser::synGet('syn52297952')
temp <- readRDS(p1$path)


temp$proteins <- rownames(temp)
temp[,189:190] <- str_split_fixed(temp$proteins, fixed("|"), 2)

names(temp)[names(temp) == 'V189'] <- 'gene_names'
temp$gene_names<-make.unique(temp$gene_names)
rownames(temp) <- temp$gene_names
temp$proteins<-NULL
temp$gene_names<-NULL
temp$V190<-NULL


gene_short_name <- rownames(temp)
table(states$State)
#there are 7 states in the male tree

#want to run each state against all other states
temp2 <- as.data.frame(t(temp))
temp2$SampleID <- rownames(temp2)
states2 <- subset(states, select=c(SampleID, State))
temp3 <- dplyr::left_join(temp2, states2, by='SampleID')


#STATE 1 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==1, 1, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_1 <- rep(0,length(gene_short_name))
l2$d_1 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_1[i] <- tk$s[1,4]
  l2$d_1[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state1a <- dplyr::left_join(dfa1,dfb1)


#STATE 2 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==2, 2, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$d_2 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$d_2[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state2 <- dplyr::left_join(dfa1,dfb1)




#STATE 3 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==3, 3, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_3 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_3[i] <- tk$s[1,4]
  l2$d_3[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state3 <- dplyr::left_join(dfa1,dfb1)




#STATE 4 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==4, 4, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_4 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_4[i] <- tk$s[1,4]
  l2$d_4[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state4 <- dplyr::left_join(dfa1,dfb1)




#STATE 5 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==5, 5, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_5 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_5[i] <- tk$s[1,4]
  l2$d_5[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state5 <- dplyr::left_join(dfa1,dfb1)






#STATE 6 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==6, 6, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_6 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_6[i] <- tk$s[1,4]
  l2$d_6[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state6 <- dplyr::left_join(dfa1,dfb1)



#STATE 7 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==7, 7, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_7 <- rep(0,length(gene_short_name))
l2$d_7 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_7[i] <- tk$s[1,4]
  l2$d_7[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state7 <- dplyr::left_join(dfa1,dfb1)



#join state-specific data frames:
allDE <- rbind(state1a, state2)
allDE <- rbind(allDE, state3)
allDE <- rbind(allDE, state4)
allDE <- rbind(allDE, state5)
allDE <- rbind(allDE, state6)
allDE <- rbind(allDE, state7)

#temp3$State <- as.character(temp3$State)
#ggplot(temp3, aes(x=State, y=CRH)) + geom_boxplot()


#save file:
write.csv(allDE, file="~/redo_DE_analysis/protM_DE_states_vs_all.csv", row.names=FALSE)















############# METABOLOMICS #################

#need Monrun object from prot_lineage_monocle_rerun.R and gene_short_name vector
#script for males following script for females (6 states in males, 7 in females)
#upload entire matrix first to look at associations with ALL genes, not just DE genes included for trajectory model

#MonRun <- readRDS(file="data_objects/MonRun_RNAseq_female.RDS")
#upload pseudotimes to get states 
p <- synapser::synGet('syn50912936')
states <- read.csv(p$path)
#temp <- readRDS(file="data_objects/Female_fulldatamatrix.RDS")
p1 <- synapser::synGet('syn52297961')
temp <- readRDS(p1$path)
gene_short_name <- rownames(temp)
table(states$State)
#there are 8 states in the female tree


#want to run each state against all other states
temp2 <- as.data.frame(t(temp))
temp2$projid <- rownames(temp2)
states2 <- subset(states, select=c(projid, State))
states2$projid <- as.character(states2$projid)
temp3 <- dplyr::left_join(temp2, states2, by='projid')


#STATE 1 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==1, 1, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_1 <- rep(0,length(gene_short_name))
l2$d_1 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_1[i] <- tk$s[1,4]
  l2$d_1[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state1a <- dplyr::left_join(dfa1,dfb1)


#STATE 2 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==2, 2, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$d_2 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$d_2[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state2 <- dplyr::left_join(dfa1,dfb1)




#STATE 3 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==3, 3, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_3 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_3[i] <- tk$s[1,4]
  l2$d_3[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state3 <- dplyr::left_join(dfa1,dfb1)




#STATE 4 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==4, 4, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_4 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_4[i] <- tk$s[1,4]
  l2$d_4[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state4 <- dplyr::left_join(dfa1,dfb1)




#STATE 5 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==5, 5, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_5 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_5[i] <- tk$s[1,4]
  l2$d_5[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state5 <- dplyr::left_join(dfa1,dfb1)






#STATE 6 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==6, 6, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_6 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_6[i] <- tk$s[1,4]
  l2$d_6[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state6 <- dplyr::left_join(dfa1,dfb1)


#STATE 7 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==7, 7, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_7 <- rep(0,length(gene_short_name))
l2$d_7 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_7[i] <- tk$s[1,4]
  l2$d_7[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state7 <- dplyr::left_join(dfa1,dfb1)

temp3$State <- as.character(temp3$State)
ggplot(temp3, aes(x=State, y=SERPINA5)) + geom_boxplot()


#STATE 8 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==8, 8, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_8 <- rep(0,length(gene_short_name))
l2$d_8 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_8[i] <- tk$s[1,4]
  l2$d_8[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state8 <- dplyr::left_join(dfa1,dfb1)

temp3$State <- as.character(temp3$State)
ggplot(temp3, aes(x=State, y=SERPINA5)) + geom_boxplot()



#join state-specific data frames:
allDE <- rbind(state1a, state2)
allDE <- rbind(allDE, state3)
allDE <- rbind(allDE, state4)
allDE <- rbind(allDE, state5)
allDE <- rbind(allDE, state6)
allDE <- rbind(allDE, state7)
allDE <- rbind(allDE, state8)


#save file:
write.csv(allDE, file="FEMALE_mets_diffabundance_each_state_vs_all.csv", row.names=FALSE)
file <- synapser::File(path='FEMALE_mets_diffabundance_each_state_vs_all.csv', parentId='syn45147359')
file <- synapser::synStore(file)



####### MALES METABOLOMICS ######

#script for males following script for females (6 states in males, 8 in females)
#upload entire matrix first to look at associations with ALL genes, not just DE genes included for trajectory model

#MonRun <- readRDS(file="data_objects/MonRun_RNAseq_female.RDS")
#upload pseudotimes to get states 
p <- synapser::synGet('syn50912976')
states <- read.csv(p$path)
#temp <- readRDS(file="data_objects/Female_fulldatamatrix.RDS")
p1 <- synapser::synGet('syn52297960')
temp <- readRDS(p1$path)
gene_short_name <- rownames(temp)
table(states$State)
#there are 6 states in the male tree


#want to run each state against all other states
temp2 <- as.data.frame(t(temp))
temp2$projid <- rownames(temp2)
states2 <- subset(states, select=c(projid, State))
states2$projid <- as.character(states2$projid)
temp3 <- dplyr::left_join(temp2, states2, by='projid')


#STATE 1 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==1, 1, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_1 <- rep(0,length(gene_short_name))
l2$d_1 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_1[i] <- tk$s[1,4]
  l2$d_1[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state1a <- dplyr::left_join(dfa1,dfb1)


#STATE 2 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==2, 2, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$d_2 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$d_2[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state2 <- dplyr::left_join(dfa1,dfb1)




#STATE 3 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==3, 3, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_3 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_3[i] <- tk$s[1,4]
  l2$d_3[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state3 <- dplyr::left_join(dfa1,dfb1)




#STATE 4 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==4, 4, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_4 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_4[i] <- tk$s[1,4]
  l2$d_4[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state4 <- dplyr::left_join(dfa1,dfb1)




#STATE 5 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==5, 5, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_5 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_5[i] <- tk$s[1,4]
  l2$d_5[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state5 <- dplyr::left_join(dfa1,dfb1)






#STATE 6 VS ALL 
state1 <- states
state1$State <- ifelse(state1$State==6, 6, 0)
table(state1$State)
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_6 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(state1$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_6[i] <- tk$s[1,4]
  l2$d_6[i] <- tk$s[1,1]
}
#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))
dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)
dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])
state6 <- dplyr::left_join(dfa1,dfb1)





#join state-specific data frames:
allDE <- rbind(state1a, state2)
allDE <- rbind(allDE, state3)
allDE <- rbind(allDE, state4)
allDE <- rbind(allDE, state5)
allDE <- rbind(allDE, state6)

table(allDE$state)

#save file:
write.csv(allDE, file="MALE_mets_diffabundance_each_state_vs_all.csv", row.names=FALSE)
file <- synapser::File(path='MALE_mets_diffabundance_each_state_vs_all.csv', parentId='syn45147359')
file <- synapser::synStore(file)
