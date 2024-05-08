setwd("/home/lheath/biodomains")

#greg cary's subdomain analysis
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("ggtree")
BiocManager::install("clusterProfiler")
install.packages("UpSetR")
install.packages("tictoc")

# packages
library(synapser)
library(fgsea)
library(org.Hs.eg.db)
library(clusterProfiler)
library(UpSetR)
library(tidyverse)
library(dplyr)
theme_set(theme_bw())


#This function quantifies the number of significantly enriched GO terms
#by AD biological domain
bd.tally <- function( enrVct, biodomDefTbl){
  bdt <- bind_cols(
    domain = unique(biodomDefTbl$Biodomain),
    n_term = map_dbl( unique(biodomDefTbl$Biodomain),
                      ~ biodomDefTbl %>% filter(Biodomain == .x) %>% 
                        dplyr::select(GOterm_Name) %>% distinct() %>% nrow()),
    n_sig_term = map_dbl( unique(biodomDefTbl$Biodomain),
                          ~ enrVct[ enrVct %in% 
                                      biodomDefTbl$GOterm_Name[
                                        biodomDefTbl$Biodomain == .x]] %>% 
                            length()) ) %>% 
    bind_rows(
      ., 
      tibble(domain = 'none', 
             n_term = setdiff(enrVct, biodomDefTbl$GOterm_Name) %>% length(), 
             n_sig_term = setdiff(enrVct, biodomDefTbl$GOterm_Name) %>% length())
    ) %>% 
    mutate(domain = fct_reorder(domain, n_sig_term, .desc = F)) %>% 
    arrange(domain)
  bdt$prop <- bdt$n_sig_term / bdt$n_term
  return(bdt)
}



# biological domain annotations
biodom <- full_join(
  
  # biodomains
  readRDS(synGet('syn25428992')$path),
  
  # domain labels
  read_csv(synGet('syn26856828')$path, col_types = cols()),  
  
  by = c('Biodomain'='domain')
) %>% mutate(Biodomain = if_else(Biodomain == 'none', NA_character_, Biodomain))

domains = biodom %>% pull(Biodomain) %>% unique() %>% sort()
subdomains <- biodom %>% select(Biodomain, subdomain_idx, Subdomain) %>% 
  distinct() %>% filter(subdomain_idx != 0)




# # # enriched biodomain terms

#First specify a list of pathways to enrich.
biodom.annotated <- biodom %>%
  filter(!is.na(n_symbol)) %>%
  select(-Biodomain) %>%
  distinct() %>%
  pull(symbol, name=GOterm_Name)
# 




################# incorporate ad vs control and state-specific biodomains analysis
#upload genes by state statistics

#AD vs control: use rnaseq harmonization project results
de_file <- synapser::synGet('syn26967458')
de1 <- read.delim(de_file$path)
de2 <- dplyr::filter(de1,Comparison=='AD_female_DLPFC - CT_female_DLPFC')
de2 <- dplyr::filter(de1,Comparison=='AD_male_DLPFC - CT_male_DLPFC')

#make short gene names unique
de2$hgnc_symbol<-make.unique(de2$hgnc_symbol)

de3 <- subset(de2, de2$P.Value<0.05)

names(de3)[names(de3) == 'hgnc_symbol'] <- 'gene_names'


#set up genelist   
gl <- de3 %>% 
  arrange(desc(logFC)) %>% 
  filter(!is.na(gene_names), !duplicated(gene_names), !is.na(logFC)) %>%
  pull(logFC, name=gene_names)



#run GSEA analysis (this can take a while, and is aided by maximizing the number of cores used via the ```nproc``` parameter)  
bd.fgsea <-  fgseaMultilevel( pathways   = biodom.annotated,
                              stats      = gl, 
                              minSize    = 1,
                              maxSize    = Inf,
                              scoreType  = 'std',
                              eps        = 0,
                              nproc      = 8 )

#tally biodomain GO terms (note: function does not account for sub-domain terms, those should be in the ```none``` category)  
bdt <- bd.tally(bd.fgsea$pathway[bd.fgsea$pval <= 0.05], biodom)

DT::datatable( arrange(bdt,desc(n_sig_term)), options = list(paging = F) ) %>% 
  DT::formatSignif(columns = 'prop', digits = 2)

#Looks like nearly all of the ```r length(sd_terms)``` sub-domains are significantly enriched in this analysis.  

#look at the top-enriched sub-domains
bd.filt.fgsea <- bd.fgsea %>% 
  filter(grepl('_SD_', pathway), ES > 0, pval < 0.05 ) %>% 
  arrange(desc(NES)) %>% distinct() %>% head(n = 42) %>% pull(pathway)

gridExtra::grid.arrange( 
  plotGseaTable(biodom.annotated[bd.filt.fgsea], gl, bd.fgsea, colwidths = c(5,2,1,0,1), render = TRUE) 
)

#annotate the enriched terms and sub-domains into biodomains for visualization  
enr <- bd.fgsea %>% 
  mutate(
    sd_term = if_else(grepl('_SD_', pathway), T, F),
    bd_abbr = if_else(sd_term, substr(pathway, 1,2), NA_character_)
  ) %>% 
  left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain, color)) %>% 
  left_join(., biodom %>% select(bd_abbr = abbr, Biodomain, color), by = 'bd_abbr') %>% 
  distinct() %>% 
  mutate(
    Biodomain = coalesce(Biodomain.x, Biodomain.y),
    color = coalesce(color.x, color.y) ) %>%
  select(-contains(c('.x', '.y'))) %>% 
  filter(Biodomain %in% c('Synapse', 'Mitochondrial Metabolism','Immune Response',
                          'Lipid Metabolism','Epigenetic')) %>% 
  mutate(
    Biodomain = fct_relevel(Biodomain, 
                            bdt %>% 
                              filter(domain %in% 
                                       c('Synapse', 'Mitochondrial Metabolism',
                                         'Immune Response','Lipid Metabolism',
                                         'Epigenetic')) %>% 
                              arrange(desc(n_sig_term)) %>% 
                              pull(domain) %>% 
                              as.character()
    )
  )

#plot the domains with sub-domains:  
tiff(file='rnaseqF_subs_ADvsCT.tiff',height=200,width=300,units='mm',res=300)
tiff(file='rnaseqM_subs_ADvsCT.tiff',height=200,width=300,units='mm',res=300)
enr %>%     
  ggplot(aes( NES, -log10(pval) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr, pval <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr, pval <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr, pval <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~Biodomain)
dev.off()


#save enr file (go terms with biodomains/subbiodomains)
enr$leadingEdge <- as.character(enr$leadingEdge)
enr <- as.data.frame(enr)
enr2 <- subset(enr, enr$pval<0.05)
#move leadingEdge column to end of dataframe
leadingEdge <- as.data.frame(enr2$leadingEdge)
#enr_female <- subset(enr2, select = -c(leadingEdge))
#enr_female <- cbind(enr_female, leadingEdge)
enr_female$sex <- 'female'


enr_male <- subset(enr2, select = -c(leadingEdge))
enr_male <- cbind(enr_male, leadingEdge)
enr_male$sex <- 'male'

enr_all <- rbind(enr_female, enr_male)
write.table(enr_all, file="RNAseq_subdomains_ADvsCT.txt", sep = "\t", row.names=FALSE)

file <- synapser::File(path='RNAseq_subdomains_ADvsCT.txt', parentId='syn53158427')
file <- synapser::synStore(file)






#####################plot subdomains state by state per data modality and by sex###############

#rnaseq female:
p <- synapser::synGet('syn53158444')
de_results <- read.csv(p$path)
#de_results <- read.csv(file='RNAseqF_DE_states_vs_all.csv')
#rnaseq male:
de_results <- read.csv(file='RNAseqM_DE_states_vs_all.csv')
#state 1:
table(de_results$state)
state1 <- subset(de_results, de_results$state==1)
state1 <- subset(de_results, de_results$state==2)
state1 <- subset(de_results, de_results$state==3)
state1 <- subset(de_results, de_results$state==4)
state1 <- subset(de_results, de_results$state==5)
state1 <- subset(de_results, de_results$state==6)
#state1 <- subset(de_results, de_results$state==7)

state1 <- subset(state1, state1$pvalue<0.05)
names(state1)[names(state1) == 'gene_names'] <- 'gene'

#set up genelist - use the Overall TRS to stratify genes  
gl <- state1 %>% 
  arrange(desc(effect)) %>% 
  filter(!is.na(gene), !duplicated(gene), !is.na(effect)) %>%
  pull(effect, name=gene)

enr = gseGO(geneList = gl, 
            ont = 'ALL', 
            OrgDb = org.Hs.eg.db, 
            keyType = 'SYMBOL',
            eps = 0,
            pvalueCutoff = 1,
            scoreType = 'std') %>% list()

#run GSEA analysis (this can take a while, and is aided by maximizing the number of cores used via the ```nproc``` parameter)  
bd.fgsea <-  fgseaMultilevel( pathways   = biodom.annotated,
                              stats      = gl, 
                              minSize    = 1,
                              maxSize    = Inf,
                              scoreType  = 'std',
                              eps        = 0,
                              nproc      = 8 )




bd.fgsea = gseGO(gl, 
            ont = 'ALL', 
            OrgDb = org.Hs.eg.db, 
            keyType = 'SYMBOL',
            eps = 0,
            pvalueCutoff = 0.1,
            scoreType = 'std') %>% list()
bdt <- bd.tally(bd.fgsea$pathway[bd.fgsea$padj <= 0.05], biodom)


results <- as.data.frame(bd.fgsea)
names(results)[names(results) == 'ID'] <- 'GO_ID'
biodom2 <- subset(biodom, select=c(Biodomain, subdomain_idx, Subdomain, TopSD, GO_ID, GOterm_Name, label, color))
results2 <- merge(results, biodom2, by="GO_ID")


left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain, color)) %>% 
  left_join(., biodom %>% select(bd_abbr = abbr, Biodomain, color), by = 'bd_abbr') %>% 


tictoc::tic()
enr.pt <- state1 %>% 
  #group_by(modality, sex, state) %>% 
  mutate( signed_P = if_else(effect < 0, log10(pvalue), -log10(pvalue) ),
          signed_P = if_else( !is.finite(signed_P), max(signed_P[is.finite(signed_P)]), signed_P)) %>% 
  summarise(
    
    # gl = effect[!duplicated(gene) & pvalue <= 0.05] %>% 
    #   setNames(., gene[!duplicated(gene) & pvalue <= 0.05]) %>% 
    
    # gl = effect[!duplicated(gene)] %>% setNames(., gene[!duplicated(gene)]) %>%
    #   sort(decreasing = T) %>% list()) %>% 
    
    gl = signed_P[!duplicated(gene)] %>% setNames(., gene[!duplicated(gene)]) %>%
      sort(decreasing = T) %>% list()) %>% 
  
  #ungroup() %>% 
  rowwise() %>% 
  mutate(
    enr = gseGO(geneList = unlist(gl), 
                ont = 'ALL', 
                OrgDb = org.Hs.eg.db, 
                keyType = 'SYMBOL',
                eps = 0,
                pvalueCutoff = 1,
                scoreType = 'std') %>% list()
  )




#tally biodomain GO terms (note: function does not account for sub-domain terms, those should be in the ```none``` category)  
bdt <- bd.tally(bd.fgsea$pathway[bd.fgsea$padj <= 0.05], biodom)

DT::datatable( arrange(bdt,desc(n_sig_term)), options = list(paging = F) ) %>% 
  DT::formatSignif(columns = 'prop', digits = 2)

#Looks like nearly all of the ```r length(sd_terms)``` sub-domains are significantly enriched in this analysis.  

#look at the top-enriched sub-domains
bd.filt.fgsea <- bd.fgsea %>% 
  filter(grepl('_SD_', pathway), ES > 0, padj < 0.05 ) %>% 
  arrange(desc(NES)) %>% distinct() %>% head(n = 42) %>% pull(pathway)

gridExtra::grid.arrange( 
  plotGseaTable(biodom.annotated[bd.filt.fgsea], gl, bd.fgsea, colwidths = c(5,2,1,0,1), render = TRUE) 
)

#annotate the enriched terms and sub-domains into biodomains for visualization  
enr2 <- enr %>% 
  mutate(
    sd_term = if_else(grepl('_SD_', pathway), T, F),
    bd_abbr = if_else(sd_term, substr(pathway, 1,2), NA_character_)
  ) %>% 
  left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain, color)) %>% 
  left_join(., biodom %>% select(bd_abbr = abbr, Biodomain, color), by = 'bd_abbr') %>% 
  distinct() %>% 
  mutate(
    Biodomain = coalesce(Biodomain.x, Biodomain.y),
    color = coalesce(color.x, color.y) ) %>%
  select(-contains(c('.x', '.y'))) %>% 
  filter(Biodomain %in% c('Synapse', 'Mitochondrial Metabolism','Immune Response',
                          'Lipid Metabolism','Epigenetic')) %>% 
  mutate(
    Biodomain = fct_relevel(Biodomain, 
                            bdt %>% 
                              filter(domain %in% 
                                       c('Synapse', 'Mitochondrial Metabolism',
                                         'Immune Response','Lipid Metabolism',
                                         'Epigenetic')) %>% 
                              arrange(desc(n_sig_term)) %>% 
                              pull(domain) %>% 
                              as.character()
    )
  )

#plot the domains with sub-domains:  
results2 %>%     
  ggplot(aes( NES, -log10(p.adjust) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(results2, p.adjust <= 0.05), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(results2, p.adjust <= 0.05), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  # ggrepel::geom_label_repel(
  #   data = subset(results2, p.adjust <= 0.05), 
  #   #aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
  #   size = 3, min.segment.length = 0
  # )+
  facet_wrap(~Biodomain)

#save enr file (go terms with biodomains/subbiodomains)
enr$leadingEdge <- as.character(enr$leadingEdge)
enr <- as.data.frame(enr)
enr2 <- subset(enr, enr$padj<0.05)
#move leadingEdge column to end of dataframe
#leadingEdge <- as.data.frame(enr2$leadingEdge)


#add states
enr2$state <- '1'
enr_state1 <- enr2
enr2$state <- '2'
enr_state2 <- enr2
enr2$state <- '3'
enr_state3 <- enr2
enr2$state <- '4'
enr_state4 <- enr2
enr2$state <- '5'
enr_state5 <- enr2
enr2$state <- '6'
enr_state6 <- enr2
enr2$state <- '7'
enr_state7 <- enr2





enr_all <- rbind(enr_state1, enr_state2)
enr_all <- rbind(enr_all, enr_state3)
enr_all <- rbind(enr_all, enr_state4)
enr_all <- rbind(enr_all, enr_state5)
enr_all <- rbind(enr_all, enr_state6)
enr_all <- rbind(enr_all, enr_state7)


leadingEdge <- as.data.frame(enr_all$leadingEdge)
enr_all <- subset(enr_all, select = -c(leadingEdge))
enr_all <- cbind(enr_all, leadingEdge)

write.table(enr_all, file="RNAseqF_subdomains_bystate.txt", sep = "\t", row.names=FALSE)
file <- synapser::File(path='RNAseqF_subdomains_bystate.txt', parentId='syn53158427')
file <- synapser::synStore(file)

write.table(enr_all, file="RNAseqM_subdomains_bystate.txt", sep = "\t", row.names=FALSE)
file <- synapser::File(path='RNAseqM_subdomains_bystate.txt', parentId='syn53158427')
file <- synapser::synStore(file)






#plot each subdomain separately for each state:

enr <- read.delim(file='RNAseqF_subdomains_bystate.txt')


enr2 <- subset(enr, enr$Biodomain=="Immune Response")
tiff(file='rnaseqF_subs_IR.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Synapse")
tiff(file='rnaseqF_subs_SY.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Lipid Metabolism")
tiff(file='rnaseqF_subs_LM.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Epigenetic")
tiff(file='rnaseqF_subs_EP.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Mitochondrial Metabolism")
tiff(file='rnaseqF_subs_MM.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


#####################plot subdomains state by state per data modality and by sex###############

#prot female:
de_results <- read.csv(file='protF_DE_states_vs_all.csv')
#prot male:
de_results <- read.csv(file='protM_DE_states_vs_all.csv')
#state 1:
table(de_results$state)
state1 <- subset(de_results, de_results$state==1)
state1 <- subset(de_results, de_results$state==2)
state1 <- subset(de_results, de_results$state==3)
state1 <- subset(de_results, de_results$state==4)
state1 <- subset(de_results, de_results$state==5)
state1 <- subset(de_results, de_results$state==6)
state1 <- subset(de_results, de_results$state==7)

state1 <- subset(state1, state1$pvalue<0.05)


#set up genelist - use the Overall TRS to stratify genes  
gl <- state1 %>% 
  arrange(desc(effect)) %>% 
  filter(!is.na(gene_names), !duplicated(gene_names), !is.na(effect)) %>%
  pull(effect, name=gene_names)



#run GSEA analysis (this can take a while, and is aided by maximizing the number of cores used via the ```nproc``` parameter)  
bd.fgsea <-  fgseaMultilevel( pathways   = biodom.annotated,
                              stats      = gl, 
                              minSize    = 1,
                              maxSize    = Inf,
                              scoreType  = 'std',
                              eps        = 0,
                              nproc      = 8 )

#tally biodomain GO terms (note: function does not account for sub-domain terms, those should be in the ```none``` category)  
bdt <- bd.tally(bd.fgsea$pathway[bd.fgsea$padj <= 0.05], biodom)

DT::datatable( arrange(bdt,desc(n_sig_term)), options = list(paging = F) ) %>% 
  DT::formatSignif(columns = 'prop', digits = 2)

#Looks like nearly all of the ```r length(sd_terms)``` sub-domains are significantly enriched in this analysis.  

#look at the top-enriched sub-domains
bd.filt.fgsea <- bd.fgsea %>% 
  filter(grepl('_SD_', pathway), ES > 0, padj < 0.05 ) %>% 
  arrange(desc(NES)) %>% distinct() %>% head(n = 42) %>% pull(pathway)

gridExtra::grid.arrange( 
  plotGseaTable(biodom.annotated[bd.filt.fgsea], gl, bd.fgsea, colwidths = c(5,2,1,0,1), render = TRUE) 
)

#annotate the enriched terms and sub-domains into biodomains for visualization  
enr <- bd.fgsea %>% 
  mutate(
    sd_term = if_else(grepl('_SD_', pathway), T, F),
    bd_abbr = if_else(sd_term, substr(pathway, 1,2), NA_character_)
  ) %>% 
  left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain, color)) %>% 
  left_join(., biodom %>% select(bd_abbr = abbr, Biodomain, color), by = 'bd_abbr') %>% 
  distinct() %>% 
  mutate(
    Biodomain = coalesce(Biodomain.x, Biodomain.y),
    color = coalesce(color.x, color.y) ) %>%
  select(-contains(c('.x', '.y'))) %>% 
  filter(Biodomain %in% c('Synapse', 'Mitochondrial Metabolism','Immune Response',
                          'Lipid Metabolism','Epigenetic')) %>% 
  mutate(
    Biodomain = fct_relevel(Biodomain, 
                            bdt %>% 
                              filter(domain %in% 
                                       c('Synapse', 'Mitochondrial Metabolism',
                                         'Immune Response','Lipid Metabolism',
                                         'Epigenetic')) %>% 
                              arrange(desc(n_sig_term)) %>% 
                              pull(domain) %>% 
                              as.character()
    )
  )

#plot the domains with sub-domains:  
enr %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0
  )+
  facet_wrap(~Biodomain)

#save enr file (go terms with biodomains/subbiodomains)
enr$leadingEdge <- as.character(enr$leadingEdge)
enr <- as.data.frame(enr)
enr2 <- subset(enr, enr$padj<0.05)
#move leadingEdge column to end of dataframe
#leadingEdge <- as.data.frame(enr2$leadingEdge)


#add states
enr2$state <- '1'
enr_state1 <- enr2
enr2$state <- '2'
enr_state2 <- enr2
enr2$state <- '3'
enr_state3 <- enr2
enr2$state <- '4'
enr_state4 <- enr2
enr2$state <- '5'
enr_state5 <- enr2
enr2$state <- '6'
enr_state6 <- enr2
enr2$state <- '7'
enr_state7 <- enr2





enr_all <- rbind(enr_state1, enr_state2)
enr_all <- rbind(enr_all, enr_state3)
enr_all <- rbind(enr_all, enr_state4)
enr_all <- rbind(enr_all, enr_state5)
enr_all <- rbind(enr_all, enr_state6)
enr_all <- rbind(enr_all, enr_state7)


leadingEdge <- as.data.frame(enr_all$leadingEdge)
enr_all <- subset(enr_all, select = -c(leadingEdge))
enr_all <- cbind(enr_all, leadingEdge)

write.table(enr_all, file="protF_subdomains_bystate.txt", sep = "\t", row.names=FALSE)
file <- synapser::File(path='protF_subdomains_bystate.txt', parentId='syn53158427')
file <- synapser::synStore(file)

write.table(enr_all, file="protM_subdomains_bystate.txt", sep = "\t", row.names=FALSE)
file <- synapser::File(path='protM_subdomains_bystate.txt', parentId='syn53158427')
file <- synapser::synStore(file)








####### male subdomains by state ###########
#####RNAseqM is wrong (i accidentally replaced it with male proteomics, need to rerun 11/14/2023)
enr <- read.delim(file='RNAseqM_subdomains_bystate.txt')


enr2 <- subset(enr, enr$Biodomain=="Immune Response")
tiff(file='rnaseqM_subs_IR.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Synapse")
tiff(file='rnaseqM_subs_SY.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Lipid Metabolism")
tiff(file='rnaseqM_subs_LM.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Epigenetic")
tiff(file='rnaseqM_subs_EP.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Mitochondrial Metabolism")
tiff(file='rnaseqM_subs_MM.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()








########################## proteomics ########################################



#upload proteins by state statistics

#AD vs control:
de_file <- synapser::synGet('syn50449887')
de1 <- read.csv(de_file$path)
de2 <- dplyr::filter(de1,comparison=='AD_female_DLPFC - CT_female_DLPFC')
de2 <- dplyr::filter(de1,comparison=='AD_male_DLPFC - CT_male_DLPFC')

#make short gene names unique
de2$GeneName<-make.unique(de2$GeneName)

de3 <- subset(de2, de2$PVal<0.05)

names(de3)[names(de3) == 'GeneName'] <- 'gene_names'





#rnaseq female:
# de_results <- read.csv(file='RNAseqF_DE_states_vs_all.csv')
# #state 1:
# state1 <- subset(de_results, de_results$state==7)

#set up genelist - use the Overall TRS to stratify genes  
gl <- de3 %>% 
  arrange(desc(Log2_FC)) %>% 
  filter(!is.na(gene_names), !duplicated(gene_names), !is.na(Log2_FC)) %>%
  pull(Log2_FC, name=gene_names)



#run GSEA analysis (this can take a while, and is aided by maximizing the number of cores used via the ```nproc``` parameter)  
bd.fgsea <-  fgseaMultilevel( pathways   = biodom.annotated,
                              stats      = gl, 
                              minSize    = 1,
                              maxSize    = Inf,
                              scoreType  = 'std',
                              eps        = 0,
                              nproc      = 8 )

#tally biodomain GO terms (note: function does not account for sub-domain terms, those should be in the ```none``` category)  
bdt <- bd.tally(bd.fgsea$pathway[bd.fgsea$pval <= 0.05], biodom)

DT::datatable( arrange(bdt,desc(n_sig_term)), options = list(paging = F) ) %>% 
  DT::formatSignif(columns = 'prop', digits = 2)

#Looks like nearly all of the ```r length(sd_terms)``` sub-domains are significantly enriched in this analysis.  

#look at the top-enriched sub-domains
bd.filt.fgsea <- bd.fgsea %>% 
  filter(grepl('_SD_', pathway), ES > 0, pval < 0.05 ) %>% 
  arrange(desc(NES)) %>% distinct() %>% head(n = 42) %>% pull(pathway)

gridExtra::grid.arrange( 
  plotGseaTable(biodom.annotated[bd.filt.fgsea], gl, bd.fgsea, colwidths = c(5,2,1,0,1), render = TRUE) 
)

#annotate the enriched terms and sub-domains into biodomains for visualization  
enr <- bd.fgsea %>% 
  mutate(
    sd_term = if_else(grepl('_SD_', pathway), T, F),
    bd_abbr = if_else(sd_term, substr(pathway, 1,2), NA_character_)
  ) %>% 
  left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain, color)) %>% 
  left_join(., biodom %>% select(bd_abbr = abbr, Biodomain, color), by = 'bd_abbr') %>% 
  distinct() %>% 
  mutate(
    Biodomain = coalesce(Biodomain.x, Biodomain.y),
    color = coalesce(color.x, color.y) ) %>%
  select(-contains(c('.x', '.y'))) %>% 
  filter(Biodomain %in% c('Synapse', 'Mitochondrial Metabolism','Immune Response',
                          'Lipid Metabolism','Epigenetic')) %>% 
  mutate(
    Biodomain = fct_relevel(Biodomain, 
                            bdt %>% 
                              filter(domain %in% 
                                       c('Synapse', 'Mitochondrial Metabolism',
                                         'Immune Response','Lipid Metabolism',
                                         'Epigenetic')) %>% 
                              arrange(desc(n_sig_term)) %>% 
                              pull(domain) %>% 
                              as.character()
    )
  )

#plot the domains with sub-domains:  
tiff(file='protF_subs_ADvsCT.tiff',height=200,width=300,units='mm',res=300)
tiff(file='protM_subs_ADvsCT.tiff',height=200,width=300,units='mm',res=300)
enr %>%     
  ggplot(aes( NES, -log10(pval) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr, pval <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr, pval <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr, pval <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~Biodomain)
dev.off()


#save enr file (go terms with biodomains/subbiodomains)
enr$leadingEdge <- as.character(enr$leadingEdge)
enr <- as.data.frame(enr)
enr2 <- subset(enr, enr$pval<0.05)
#move leadingEdge column to end of dataframe
leadingEdge <- as.data.frame(enr2$leadingEdge)
#enr_female <- subset(enr2, select = -c(leadingEdge))
#enr_female <- cbind(enr_female, leadingEdge)
enr_female$sex <- 'female'


enr_male <- subset(enr2, select = -c(leadingEdge))
enr_male <- cbind(enr_male, leadingEdge)
enr_male$sex <- 'male'

enr_all <- rbind(enr_female, enr_male)
write.table(enr_all, file="Prot_subdomains_ADvsCT.txt", sep = "\t", row.names=FALSE)

file <- synapser::File(path='Prot_subdomains_ADvsCT.txt', parentId='syn53158427')
file <- synapser::synStore(file)








enr <- read.delim(file='protF_subdomains_bystate.txt')


enr2 <- subset(enr, enr$Biodomain=="Immune Response")
tiff(file='protF_subs_IR.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Synapse")
tiff(file='protF_subs_SY.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Lipid Metabolism")
tiff(file='protF_subs_LM.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Epigenetic")
tiff(file='protF_subs_EP.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Mitochondrial Metabolism")
tiff(file='protF_subs_MM.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()



####### male subdomains by state ###########
enr <- read.delim(file='protM_subdomains_bystate.txt')


enr2 <- subset(enr, enr$Biodomain=="Immune Response")
tiff(file='protM_subs_IR.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Synapse")
tiff(file='protM_subs_SY.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Lipid Metabolism")
tiff(file='protM_subs_LM.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Epigenetic")
tiff(file='protM_subs_EP.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()


enr2 <- subset(enr, enr$Biodomain=="Mitochondrial Metabolism")
tiff(file='protM_subs_MM.tiff',height=200,width=300,units='mm',res=300)
enr2 %>%     
  ggplot(aes( NES, -log10(padj) ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(range = c(.7,5), trans = 'log10')+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == F), 
    shape = 21, alpha = .2, color = 'grey30', size = 2,
    aes(fill = color) )+
  geom_point(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    shape = 23, alpha = .5, color = 'grey30', size = 5,
    aes(fill = color) )+
  ggrepel::geom_label_repel(
    data = subset(enr2, padj <= 0.05 & sd_term == T), 
    aes(label = str_remove_all(pathway, '^.{2}_SD_') ),
    size = 3, min.segment.length = 0, max.overlaps = 100
  )+
  facet_wrap(~state)
dev.off()



