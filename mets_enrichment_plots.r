library(dplyr)
library(ggplot2)
#Make metabolomics plots--number of sub pathways in each state
  
#upload differential abundance by state data
  
#female data
p <- synapser::synGet('syn58961427')
de_metsF <- read.csv(p$path)

#upload metabolite metadata
p1 <- synapser::synGet('syn59401895')
metabolites <- read.csv(p1$path)
#rename first column
names(metabolites)[names(metabolites) == '...1'] <- 'met_name'
#change column name in DE results to match
names(de_metsF)[names(de_metsF) == 'gene_names'] <- 'met_name'
#replace 'X' metabolite sub-pathway with 'uncharacterized'
metabolites[c("SUPER_PATHWAY","SUB_PATHWAY")][is.na(metabolites[c("SUPER_PATHWAY","SUB_PATHWAY")])] <- "uncharacterized"
#join by met_name
metsF_results <- merge(de_metsF, metabolites)
metsF_results$direction <- ifelse(metsF_results$effect>0, "positive", "negative")

# metsF <- subset(metsF_results, metsF_results$pvalue<0.01)
# table(metsF$SUB_PATHWAY, metsF$state)
# #delete uncharacterized pathways
# metsF <- subset(metsF, metsF$SUB_PATHWAY!="uncharacterized")
# table(metsF$SUB_PATHWAY, metsF$state)
# 
# df <- metsF_up2 %>% count(state, SUB_PATHWAY) %>%
#   ungroup()
# 
# df2 <- metsF %>%
#   tidyr::gather(key = 'SUB_PATHWAY', value = 'direction') %>%
#   group_by(state) %>%
#   count(SUB_PATHWAY)
# 
# #need at least 3 metabolites in a sub pathway
# df <- subset(df, df$n>3)
# 
# #add back in the super pathways
# supers <- subset(metabolites, select=c(SUB_PATHWAY, SUPER_PATHWAY))
# supers <- distinct(supers)
# df <- merge(df, supers)





#create separate plots for mets going up vs mets going down, incorporate directionality into n, then rejoin data frames
metsF_up <- subset(metsF_results, metsF_results$effect>0)
metsF_down <- subset(metsF_results, metsF_results$effect<0)

metsF_up <- subset(metsF_up, metsF_up$pvalue<0.01)
table(metsF_up$SUB_PATHWAY, metsF_up$state)
#delete uncharacterized pathways
metsF_up2 <- subset(metsF_up, metsF_up$SUB_PATHWAY!="uncharacterized")
table(metsF_up2$SUB_PATHWAY, metsF_up2$state)

up <- metsF_up2 %>% count(state, SUB_PATHWAY)
#need at least 3 metabolites in a sub pathway
up <- subset(up, up$n>3)

#add back in the super pathways
supers <- subset(metabolites, select=c(SUB_PATHWAY, SUPER_PATHWAY))
supers <- distinct(supers)
up <- merge(up, supers)


up <- up %>%
  group_by(state) %>%
  arrange(desc(n), .by_group=TRUE)


metsF_down <- subset(metsF_results, metsF_results$effect<0)

metsF_down <- subset(metsF_down, metsF_down$pvalue<0.01)
table(metsF_down$SUB_PATHWAY, metsF_down$state)
#delete uncharacterized pathways
metsF_down <- subset(metsF_down, metsF_down$SUB_PATHWAY!="uncharacterized")
table(metsF_down$SUB_PATHWAY, metsF_down$state)

down <- metsF_down %>% count(state, SUB_PATHWAY)
down <- down %>%
  group_by(state) %>%
  arrange(desc(n), .by_group=TRUE)

#need at least 3 metabolites in a sub pathway
down <- subset(down, down$n>3)

#add back in the super pathways
supers <- subset(metabolites, select=c(SUB_PATHWAY, SUPER_PATHWAY))
supers <- distinct(supers)
down <- merge(down, supers)

#add a negative sign to n for the downs
down$n <- down$n*(-1)

allmets <- rbind(up, down)

#write.csv(allmets, "allmets_for_colors.csv", row.names=FALSE)
aminos <- subset(allmets, allmets$SUPER_PATHWAY=="Amino Acid")
carbs <- subset(allmets, allmets$SUPER_PATHWAY=="Carbohydrate")
cofacs <- subset(allmets, allmets$SUPER_PATHWAY=='Cofactors and Vitamins')

lipids <- subset(allmets, allmets$SUPER_PATHWAY=='Lipid')

#assign colors based on super pathway groupings to each sub pathway
allmets$color <- allmets$SUB_PATHWAY
allmets$color[allmets$color == 'Leucine, Isoleucine and Valine Metabolism'] <- 'blue4'
allmets$color[allmets$color == 'Methionine, Cysteine, SAM and Taurine Metabolism'] <- 'dodgerblue3'
allmets$color[allmets$color == 'Urea cycle; Arginine and Proline Metabolism'] <- 'deepskyblue3'
allmets$color[allmets$color == 'Lysine Metabolism'] <- 'lightslateblue'
allmets$color[allmets$color == 'Glutamate Metabolism'] <- 'mediumblue'
allmets$color[allmets$color == 'Glycine, Serine and Threonine Metabolism'] <- 'lightsteelblue2'
allmets$color[allmets$color == 'Histidine Metabolism'] <- 'blueviolet'
allmets$color[allmets$color == 'Tyrosine Metabolism'] <- 'mediumpurple4'
allmets$color[allmets$color == 'Tryptophan Metabolism'] <- 'mediumorchid'
allmets$color[allmets$color == 'Glutathione Metabolism'] <- 'plum2'

allmets$color[allmets$color == 'Pentose Metabolism'] <- 'hotpink2'
allmets$color[allmets$color == 'Fructose, Mannose and Galactose Metabolism'] <- 'red2'
allmets$color[allmets$color == 'Glycolysis, Gluconeogenesis, and Pyruvate Metabolism'] <- 'deeppink'
allmets$color[allmets$color == 'Aminosugar Metabolism'] <- 'magenta1'

allmets$color[allmets$color == 'Nicotinate and Nicotinamide Metabolism'] <- 'aquamarine'

allmets$color[allmets$color == 'TCA Cycle'] <- 'yellow'

allmets$color[allmets$color == 'Phosphatidylcholine (PC)'] <- 'tomato'
allmets$color[allmets$color == 'Lysophospholipid'] <- 'salmon1'
allmets$color[allmets$color == 'Phosphatidylethanolamine (PE)'] <- 'sandybrown'
allmets$color[allmets$color == 'Phospholipid Metabolism'] <- 'tan'
allmets$color[allmets$color == 'Phosphatidylinositol (PI)'] <- 'sienna'
allmets$color[allmets$color == 'Ceramides'] <- 'navajowhite'
allmets$color[allmets$color == 'Diacylglycerol'] <- 'rosybrown'
allmets$color[allmets$color == 'Sphingomyelins'] <- 'goldenrod1'
allmets$color[allmets$color == 'Dihydrosphingomyelins'] <- 'darkorange'
allmets$color[allmets$color == 'Long Chain Polyunsaturated Fatty Acid (n3 and n6)'] <- 'gold4'
allmets$color[allmets$color == 'Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)'] <- 'darkolivegreen'
allmets$color[allmets$color == 'Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)'] <- 'mediumseagreen'
allmets$color[allmets$color == 'Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)'] <- 'olivedrab2'
allmets$color[allmets$color == 'Fatty Acid Metabolism (Acyl Carnitine, Hydroxy)'] <- 'springgreen'
allmets$color[allmets$color == 'Fatty Acid, Monohydroxy'] <- 'seagreen4'
allmets$color[allmets$color == 'Fatty Acid, Dicarboxylate'] <- 'green3'
allmets$color[allmets$color == 'Long Chain Saturated Fatty Acid'] <- 'chartreuse'
allmets$color[allmets$color == 'Endocannabinoid'] <- 'gold1'

allmets$color[allmets$color == 'Purine Metabolism, (Hypo)Xanthine/Inosine containing'] <- 'grey40'
allmets$color[allmets$color == 'Purine Metabolism, Guanine containing'] <- 'black'
allmets$color[allmets$color == 'Purine Metabolism, Adenine containing'] <- 'grey80'

allmets$color[allmets$color == 'Gamma-glutamyl Amino Acid'] <- 'palevioletred'
allmets$color[allmets$color == 'Dipeptide'] <- 'lightpink1'

#order subpathways for legend:
allmets$SUB_PATHWAY <- factor(allmets$SUB_PATHWAY, levels=c("Leucine, Isoleucine and Valine Metabolism","Methionine, Cysteine, SAM and Taurine Metabolism",
                                                            "Urea cycle; Arginine and Proline Metabolism","Lysine Metabolism","Glutamate Metabolism","Glycine, Serine and Threonine Metabolism",
                                                            "Histidine Metabolism","Tyrosine Metabolism","Tryptophan Metabolism","Glutathione Metabolism","Pentose Metabolism",
                                                            "Fructose, Mannose and Galactose Metabolism","Glycolysis, Gluconeogenesis, and Pyruvate Metabolism","Aminosugar Metabolism",
                                                            "Nicotinate and Nicotinamide Metabolism","TCA Cycle","Phosphatidylcholine (PC)","Lysophospholipid","Phosphatidylethanolamine (PE)",
                                                            "Phospholipid Metabolism","Phosphatidylinositol (PI)","Ceramides","Diacylglycerol","Sphingomyelins","Dihydrosphingomyelins",
                                                            "Long Chain Polyunsaturated Fatty Acid (n3 and n6)","Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)","Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)",
                                                            "Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)","Fatty Acid Metabolism (Acyl Carnitine, Hydroxy)","Fatty Acid, Monohydroxy",
                                                            "Fatty Acid, Dicarboxylate","Long Chain Saturated Fatty Acid","Endocannabinoid","Purine Metabolism, (Hypo)Xanthine/Inosine containing","Purine Metabolism, Guanine containing",
                                                            "Purine Metabolism, Adenine containing","Gamma-glutamyl Amino Acid","Dipeptide"))

col <- as.character(allmets$color)
names(col) <- as.character(allmets$SUB_PATHWAY)

allmets$state <- as.factor(allmets$state)
p <- ggplot(allmets, aes(fill=SUB_PATHWAY, y=n, x=state)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=col)+
  facet_wrap(~SUPER_PATHWAY) + theme(legend.text = element_text(size=6), legend.key.size = unit(0.5, 'cm'))
p








#repeat for male metabolites
#upload differential abundance by state data

#female data
p <- synapser::synGet('syn58961717')
de_metsF <- read.csv(p$path)

#upload metabolite metadata
p1 <- synapser::synGet('syn59401895')
metabolites <- read.csv(p1$path)
#rename first column
names(metabolites)[names(metabolites) == '...1'] <- 'met_name'
#change column name in DE results to match
names(de_metsF)[names(de_metsF) == 'gene_names'] <- 'met_name'
#replace 'X' metabolite sub-pathway with 'uncharacterized'
metabolites[c("SUPER_PATHWAY","SUB_PATHWAY")][is.na(metabolites[c("SUPER_PATHWAY","SUB_PATHWAY")])] <- "uncharacterized"
#join by met_name
metsF_results <- merge(de_metsF, metabolites)
metsF_results$direction <- ifelse(metsF_results$effect>0, "positive", "negative")

# metsF <- subset(metsF_results, metsF_results$pvalue<0.01)
# table(metsF$SUB_PATHWAY, metsF$state)
# #delete uncharacterized pathways
# metsF <- subset(metsF, metsF$SUB_PATHWAY!="uncharacterized")
# table(metsF$SUB_PATHWAY, metsF$state)
# 
# df <- metsF_up2 %>% count(state, SUB_PATHWAY) %>%
#   ungroup()
# 
# df2 <- metsF %>%
#   tidyr::gather(key = 'SUB_PATHWAY', value = 'direction') %>%
#   group_by(state) %>%
#   count(SUB_PATHWAY)
# 
# #need at least 3 metabolites in a sub pathway
# df <- subset(df, df$n>3)
# 
# #add back in the super pathways
# supers <- subset(metabolites, select=c(SUB_PATHWAY, SUPER_PATHWAY))
# supers <- distinct(supers)
# df <- merge(df, supers)

#create separate dataframes for mets going up vs mets going down, incorporate directionality into n, then rejoin data frames
metsF_up <- subset(metsF_results, metsF_results$effect>0)
metsF_down <- subset(metsF_results, metsF_results$effect<0)

metsF_up <- subset(metsF_up, metsF_up$pvalue<0.01)
table(metsF_up$SUB_PATHWAY, metsF_up$state)
#delete uncharacterized pathways
metsF_up2 <- subset(metsF_up, metsF_up$SUB_PATHWAY!="uncharacterized")
table(metsF_up2$SUB_PATHWAY, metsF_up2$state)

up <- metsF_up2 %>% count(state, SUB_PATHWAY)
#need at least 3 metabolites in a sub pathway
up <- subset(up, up$n>3)

#add back in the super pathways
supers <- subset(metabolites, select=c(SUB_PATHWAY, SUPER_PATHWAY))
supers <- distinct(supers)
up <- merge(up, supers)


up <- up %>%
  group_by(state) %>%
  arrange(desc(n), .by_group=TRUE)


metsF_down <- subset(metsF_results, metsF_results$effect<0)

metsF_down <- subset(metsF_down, metsF_down$pvalue<0.01)
table(metsF_down$SUB_PATHWAY, metsF_down$state)
#delete uncharacterized pathways
metsF_down <- subset(metsF_down, metsF_down$SUB_PATHWAY!="uncharacterized")
table(metsF_down$SUB_PATHWAY, metsF_down$state)

down <- metsF_down %>% count(state, SUB_PATHWAY)
down <- down %>%
  group_by(state) %>%
  arrange(desc(n), .by_group=TRUE)

#need at least 3 metabolites in a sub pathway
down <- subset(down, down$n>3)

#add back in the super pathways
supers <- subset(metabolites, select=c(SUB_PATHWAY, SUPER_PATHWAY))
supers <- distinct(supers)
down <- merge(down, supers)

#add a negative sign to n for the downs
down$n <- down$n*(-1)

allmets <- rbind(up, down)

#write.csv(allmets, "allmets_for_colors.csv", row.names=FALSE)
aminos <- subset(allmets, allmets$SUPER_PATHWAY=="Amino Acid")
carbs <- subset(allmets, allmets$SUPER_PATHWAY=="Carbohydrate")
cofacs <- subset(allmets, allmets$SUPER_PATHWAY=='Cofactors and Vitamins')

lipids <- subset(allmets, allmets$SUPER_PATHWAY=='Lipid')

#assign colors based on super pathway groupings to each sub pathway
allmets$color <- allmets$SUB_PATHWAY
allmets$color <- as.character(allmets$color)
allmets$color[allmets$color == 'Leucine, Isoleucine and Valine Metabolism'] <- 'blue4'
allmets$color[allmets$color == 'Methionine, Cysteine, SAM and Taurine Metabolism'] <- 'dodgerblue3'
allmets$color[allmets$color == 'Urea cycle; Arginine and Proline Metabolism'] <- 'deepskyblue3'
allmets$color[allmets$color == 'Histidine Metabolism'] <- 'blueviolet'
allmets$color[allmets$color == 'Tryptophan Metabolism'] <- 'mediumorchid'
allmets$color[allmets$color == 'Tyrosine Metabolism'] <- 'mediumpurple4'
allmets$color[allmets$color == 'Glutathione Metabolism'] <- 'plum2'
allmets$color[allmets$color == 'Glutamate Metabolism'] <- 'mediumblue'
allmets$color[allmets$color == 'Glycine, Serine and Threonine Metabolism'] <- 'lightsteelblue2'
#allmets$color[allmets$color == 'Lysine Metabolism'] <- 'lightslateblue'

allmets$color[allmets$color == 'Glycolysis, Gluconeogenesis, and Pyruvate Metabolism'] <- 'deeppink'
allmets$color[allmets$color == 'Nucleotide Sugar'] <- 'mediumvioletred'
#allmets$color[allmets$color == 'Pentose Metabolism'] <- 'hotpink2'
#allmets$color[allmets$color == 'Fructose, Mannose and Galactose Metabolism'] <- 'red2'
#allmets$color[allmets$color == 'Aminosugar Metabolism'] <- 'magenta1'

allmets$color[allmets$color == 'Nicotinate and Nicotinamide Metabolism'] <- 'aquamarine'

allmets$color[allmets$color == 'TCA Cycle'] <- 'yellow'

allmets$color[allmets$color == 'Phosphatidylcholine (PC)'] <- 'tomato'
allmets$color[allmets$color == 'Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)'] <- 'darkolivegreen'
allmets$color[allmets$color == 'Phosphatidylethanolamine (PE)'] <- 'sandybrown'
allmets$color[allmets$color == 'Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)'] <- 'olivedrab2'
allmets$color[allmets$color == 'Fatty Acid, Dicarboxylate'] <- 'green3'
allmets$color[allmets$color == 'Diacylglycerol'] <- 'rosybrown'
allmets$color[allmets$color == 'Phosphatidylinositol (PI)'] <- 'sienna'
allmets$color[allmets$color == 'Sphingomyelins'] <- 'goldenrod1'
allmets$color[allmets$color == 'Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)'] <- 'mediumseagreen'
allmets$color[allmets$color == 'Ceramides'] <- 'navajowhite'
allmets$color[allmets$color == 'Hexosylceramides (HCER)'] <- 'firebrick3'
allmets$color[allmets$color == 'Phospholipid Metabolism'] <- 'tan'

#allmets$color[allmets$color == 'Lysophospholipid'] <- 'salmon1'
#allmets$color[allmets$color == 'Dihydrosphingomyelins'] <- 'darkorange'
#allmets$color[allmets$color == 'Long Chain Polyunsaturated Fatty Acid (n3 and n6)'] <- 'gold4'
#allmets$color[allmets$color == 'Fatty Acid Metabolism (Acyl Carnitine, Hydroxy)'] <- 'springgreen'
#allmets$color[allmets$color == 'Fatty Acid, Monohydroxy'] <- 'seagreen4'
#allmets$color[allmets$color == 'Long Chain Saturated Fatty Acid'] <- 'chartreuse'
#allmets$color[allmets$color == 'Endocannabinoid'] <- 'gold1'

allmets$color[allmets$color == 'Purine Metabolism, (Hypo)Xanthine/Inosine containing'] <- 'grey40'
allmets$color[allmets$color == 'Purine Metabolism, Guanine containing'] <- 'black'
allmets$color[allmets$color == 'Purine Metabolism, Adenine containing'] <- 'grey80'

allmets$color[allmets$color == 'Gamma-glutamyl Amino Acid'] <- 'palevioletred'
allmets$color[allmets$color == 'Dipeptide'] <- 'lightpink1'

#order subpathways for legend:
allmets$SUB_PATHWAY <- factor(allmets$SUB_PATHWAY, levels=c("Leucine, Isoleucine and Valine Metabolism","Methionine, Cysteine, SAM and Taurine Metabolism",
                                                            "Urea cycle; Arginine and Proline Metabolism","Glutamate Metabolism","Glycine, Serine and Threonine Metabolism",
                                                            "Tyrosine Metabolism","Tryptophan Metabolism","Glutathione Metabolism",
                                                            "Glycolysis, Gluconeogenesis, and Pyruvate Metabolism","Nucleotide Sugar",
                                                            "Nicotinate and Nicotinamide Metabolism","TCA Cycle","Phosphatidylcholine (PC)","Phosphatidylethanolamine (PE)",
                                                            "Phospholipid Metabolism","Phosphatidylinositol (PI)","Ceramides","Diacylglycerol","Sphingomyelins",
                                                            "Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)","Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)",
                                                            "Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)","Fatty Acid, Dicarboxylate", "Hexosylceramides (HCER)",
                                                            "Purine Metabolism, (Hypo)Xanthine/Inosine containing","Purine Metabolism, Guanine containing",
                                                            "Purine Metabolism, Adenine containing","Gamma-glutamyl Amino Acid","Dipeptide"))

col <- as.character(allmets$color)
names(col) <- as.character(allmets$SUB_PATHWAY)

allmets$state <- as.factor(allmets$state)
p <- ggplot(allmets, aes(fill=SUB_PATHWAY, y=n, x=state)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=col)+
  facet_wrap(~SUPER_PATHWAY) + theme(legend.text = element_text(size=6), legend.key.size = unit(0.5, 'cm'))
p
