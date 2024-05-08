

#### RNAseqF by state #######
enrMerge <- read.csv(file="RNAseqF_subdomains_by_state.csv")

#change state names
enrMerge$state2 <- sub("^", "State ", enrMerge$state)

#limit to only pathways that are included as subdomains
enrMerge2 <- enrMerge[!is.na(enrMerge$Subdomain),]
table(enrMerge2$Subdomain)
table(enrMerge2$Biodomain)


Ap <- subset(enrMerge2, enrMerge2$Biodomain=='Apoptosis')
enrMerge <- Ap
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_Apsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Apoptosis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


AM <- subset(enrMerge2, enrMerge2$Biodomain=='APP Metabolism')
enrMerge <- AM
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_AMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, APP Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Au <- subset(enrMerge2, enrMerge2$Biodomain=='Autophagy')
enrMerge <- Au
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_Ausubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Autophagy Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


CC <- subset(enrMerge2, enrMerge2$Biodomain=='Cell Cycle')
enrMerge <- CC
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_CCsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Cell Cycle Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



DR <- subset(enrMerge2, enrMerge2$Biodomain=='DNA Repair')
enrMerge <- DR
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_DRsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, DNA Repair Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


En <- subset(enrMerge2, enrMerge2$Biodomain=='Endolysosome')
enrMerge <- En
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_Ensubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Endolysosome Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Ep <- subset(enrMerge2, enrMerge2$Biodomain=='Epigenetic')
enrMerge <- Ep
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_Epsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Epigenetic Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



IR <- subset(enrMerge2, enrMerge2$Biodomain=='Immune Response')
enrMerge <- IR
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_IRsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Immune Response Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


LM <- subset(enrMerge2, enrMerge2$Biodomain=='Lipid Metabolism')
enrMerge <- LM
tiff(file='subdomains_figures/rnaseqF_LMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Lipid Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


MB <- subset(enrMerge2, enrMerge2$Biodomain=='Metal Binding and Homeostasis')
enrMerge <- MB
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqF_MBsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Metal Binding & Homeostasis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()

MM <- subset(enrMerge2, enrMerge2$Biodomain=='Mitochondrial Metabolism')
enrMerge <- MM
tiff(file='subdomains_figures/rnaseqF_MMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Mitochondrial Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



My <- subset(enrMerge2, enrMerge2$Biodomain=='Myelination')
enrMerge <- My
tiff(file='subdomains_figures/rnaseqF_Mysubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Myelination Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Ox <- subset(enrMerge2, enrMerge2$Biodomain=='Oxidative Stress')
enrMerge <- Ox
tiff(file='subdomains_figures/rnaseqF_Oxsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Oxidation Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



Pr <- subset(enrMerge2, enrMerge2$Biodomain=='Proteostasis')
enrMerge <- Pr
tiff(file='subdomains_figures/rnaseqF_Prsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Proteostasis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



SS <- subset(enrMerge2, enrMerge2$Biodomain=='Structural Stabilization')
enrMerge <- SS
tiff(file='subdomains_figures/rnaseqF_SSsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Structural Stabilization Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()




Sy <- subset(enrMerge2, enrMerge2$Biodomain=='Synapse')
enrMerge <- Sy
#present this a different way:
tiff(file='subdomains_figures/rnaseqF_SYsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Synapse Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Va <- subset(enrMerge2, enrMerge2$Biodomain=='Vasculature')
enrMerge <- Va
tiff(file='subdomains_figures/rnaseqF_VAsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Transcriptomics, Vasculature Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()








####### rnaseq MALES by state
enrMergeM <- read.csv(file="RNAseqM_subdomains_by_state.csv")
#change state names
enrMergeM$state2 <- sub("^", "State ", enrMergeM$state)

#limit to only pathways that are included as subdomains
enrMerge2 <- enrMergeM[!is.na(enrMergeM$Subdomain),]
table(enrMerge2$Subdomain)
table(enrMerge2$Biodomain)


Ap <- subset(enrMerge2, enrMerge2$Biodomain=='Apoptosis')
enrMerge <- Ap
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_Apsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Apoptosis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


AM <- subset(enrMerge2, enrMerge2$Biodomain=='APP Metabolism')
enrMerge <- AM
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_AMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, APP Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Au <- subset(enrMerge2, enrMerge2$Biodomain=='Autophagy')
enrMerge <- Au
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_Ausubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Autophagy Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


CC <- subset(enrMerge2, enrMerge2$Biodomain=='Cell Cycle')
enrMerge <- CC
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_CCsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Cell Cycle Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



DR <- subset(enrMerge2, enrMerge2$Biodomain=='DNA Repair')
enrMerge <- DR
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_DRsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, DNA Repair Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


En <- subset(enrMerge2, enrMerge2$Biodomain=='Endolysosome')
enrMerge <- En
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_Ensubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Endolysosome Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Ep <- subset(enrMerge2, enrMerge2$Biodomain=='Epigenetic')
enrMerge <- Ep
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_Epsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Epigenetic Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



IR <- subset(enrMerge2, enrMerge2$Biodomain=='Immune Response')
enrMerge <- IR
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_IRsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Immune Response Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


LM <- subset(enrMerge2, enrMerge2$Biodomain=='Lipid Metabolism')
enrMerge <- LM
tiff(file='subdomains_figures/rnaseqM_LMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Lipid Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


MB <- subset(enrMerge2, enrMerge2$Biodomain=='Metal Binding and Homeostasis')
enrMerge <- MB
#chart subdomains individually:
tiff(file='subdomains_figures/rnaseqM_MBsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Metal Binding & Homeostasis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()

MM <- subset(enrMerge2, enrMerge2$Biodomain=='Mitochondrial Metabolism')
enrMerge <- MM
tiff(file='subdomains_figures/rnaseqM_MMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Mitochondrial Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



My <- subset(enrMerge2, enrMerge2$Biodomain=='Myelination')
enrMerge <- My
tiff(file='subdomains_figures/rnaseqM_Mysubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Myelination Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Ox <- subset(enrMerge2, enrMerge2$Biodomain=='Oxidative Stress')
enrMerge <- Ox
tiff(file='subdomains_figures/rnaseqM_Oxsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Oxidation Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



Pr <- subset(enrMerge2, enrMerge2$Biodomain=='Proteostasis')
enrMerge <- Pr
tiff(file='subdomains_figures/rnaseqM_Prsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Proteostasis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



SS <- subset(enrMerge2, enrMerge2$Biodomain=='Structural Stabilization')
enrMerge <- SS
tiff(file='subdomains_figures/rnaseqM_SSsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Structural Stabilization Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()




Sy <- subset(enrMerge2, enrMerge2$Biodomain=='Synapse')
enrMerge <- Sy
#present this a different way:
tiff(file='subdomains_figures/rnaseqM_SYsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Synapse Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Va <- subset(enrMerge2, enrMerge2$Biodomain=='Vasculature')
enrMerge <- Va
tiff(file='subdomains_figures/rnaseqM_VAsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Transcriptomics, Vasculature Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()





















####### proteomics ##########

####### female proteomics subdomains ########
enrMerge <- read.csv(file='protF_subdomains_by_state.csv')

#change state names
enrMerge$state2 <- sub("^", "State ", enrMerge$state)

#limit to only pathways that are included as subdomains
enrMerge2 <- enrMerge[!is.na(enrMerge$Subdomain),]
table(enrMerge2$Subdomain)
table(enrMerge2$Biodomain)


Ap <- subset(enrMerge2, enrMerge2$Biodomain=='Apoptosis')
enrMerge <- Ap
#chart subdomains individually:
tiff(file='subdomains_figures/protF_Apsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Apoptosis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()

# no APP Metabolism subdomains in female proteomics results
# AM <- subset(enrMerge2, enrMerge2$Biodomain=='APP Metabolism')
# enrMerge <- AM
# #chart subdomains individually:
# tiff(file='subdomains_figures/protF_AMsubs.tiff',height=110,width=300,units='mm',res=300)
# enrMerge %>%   
#   ggplot(aes( factor(state), NES ))+ 
#   scale_color_identity()+ scale_fill_identity()+
#   scale_size_continuous(
#     limits = -log10(enrMerge$pvalue) %>% range(),
#     range = c(.7,5))+
#   geom_jitter(
#     data = subset(enrMerge, pvalue > 0.05),
#     color = 'grey85', alpha = .2, size = .5, show.legend = F
#   )+
#   geom_jitter(
#     data = subset(enrMerge, pvalue < 0.05 ),
#     aes(color = color, size = -log10(pvalue) ), #
#     alpha = .5, show.legend = T
#   )+
#   geom_violin(
#     data = subset(enrMerge, pvalue < 0.05 & NES < 0),
#     scale='width', aes(fill = color), color = 'grey50',#, color = col
#     alpha = .3, show.legend = F
#   )+
#   geom_violin(
#     data = subset(enrMerge, pvalue < 0.05 & NES > 0),
#     scale='width', aes(fill = color), color = 'grey50',#, color = col
#     alpha = .3, show.legend = F
#   )+
#   geom_hline(yintercept = 0, lty = 3)+
#   facet_wrap(~Subdomain, ncol=5)+
#   labs(x='')+  coord_flip() + theme(legend.position = 'right') +
#   ggtitle('Female Proteomics, APP Metabolism Subdomains; term pvalue < 0.05'
#   )+ xlab("State") + scale_x_discrete(limits = rev)
# dev.off()


Au <- subset(enrMerge2, enrMerge2$Biodomain=='Autophagy')
enrMerge <- Au
#chart subdomains individually:
tiff(file='subdomains_figures/protF_Ausubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Autophagy Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


CC <- subset(enrMerge2, enrMerge2$Biodomain=='Cell Cycle')
enrMerge <- CC
#chart subdomains individually:
tiff(file='subdomains_figures/protF_CCsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Cell Cycle Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



DR <- subset(enrMerge2, enrMerge2$Biodomain=='DNA Repair')
enrMerge <- DR
#chart subdomains individually:
tiff(file='subdomains_figures/protF_DRsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, DNA Repair Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


En <- subset(enrMerge2, enrMerge2$Biodomain=='Endolysosome')
enrMerge <- En
#chart subdomains individually:
tiff(file='subdomains_figures/protF_Ensubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Endolysosome Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Ep <- subset(enrMerge2, enrMerge2$Biodomain=='Epigenetic')
enrMerge <- Ep
#chart subdomains individually:
tiff(file='subdomains_figures/protF_Epsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Epigenetic Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



IR <- subset(enrMerge2, enrMerge2$Biodomain=='Immune Response')
enrMerge <- IR
#chart subdomains individually:
tiff(file='subdomains_figures/protF_IRsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Immune Response Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


LM <- subset(enrMerge2, enrMerge2$Biodomain=='Lipid Metabolism')
enrMerge <- LM
tiff(file='subdomains_figures/protF_LMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Lipid Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


MB <- subset(enrMerge2, enrMerge2$Biodomain=='Metal Binding and Homeostasis')
enrMerge <- MB
#chart subdomains individually:
tiff(file='subdomains_figures/protF_MBsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Metal Binding & Homeostasis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()

MM <- subset(enrMerge2, enrMerge2$Biodomain=='Mitochondrial Metabolism')
enrMerge <- MM
tiff(file='subdomains_figures/protF_MMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Mitochondrial Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



My <- subset(enrMerge2, enrMerge2$Biodomain=='Myelination')
enrMerge <- My
tiff(file='subdomains_figures/protF_Mysubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Myelination Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Ox <- subset(enrMerge2, enrMerge2$Biodomain=='Oxidative Stress')
enrMerge <- Ox
tiff(file='subdomains_figures/protF_Oxsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Oxidation Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



Pr <- subset(enrMerge2, enrMerge2$Biodomain=='Proteostasis')
enrMerge <- Pr
tiff(file='subdomains_figures/protF_Prsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Proteostasis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



SS <- subset(enrMerge2, enrMerge2$Biodomain=='Structural Stabilization')
enrMerge <- SS
tiff(file='subdomains_figures/protF_SSsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Structural Stabilization Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()




Sy <- subset(enrMerge2, enrMerge2$Biodomain=='Synapse')
enrMerge <- Sy
#present this a different way:
tiff(file='subdomains_figures/protF_SYsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Synapse Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Va <- subset(enrMerge2, enrMerge2$Biodomain=='Vasculature')
enrMerge <- Va
tiff(file='subdomains_figures/protF_VAsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Female Proteomics, Vasculature Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()









####### Male proteomics subdomains ########
enrMerge <- read.csv(file='protM_subdomains_by_state.csv')

#change state names
enrMerge$state2 <- sub("^", "State ", enrMerge$state)

#limit to only pathways that are included as subdomains
enrMerge2 <- enrMerge[!is.na(enrMerge$Subdomain),]
table(enrMerge2$Subdomain)
table(enrMerge2$Biodomain)


Ap <- subset(enrMerge2, enrMerge2$Biodomain=='Apoptosis')
enrMerge <- Ap
#chart subdomains individually:
tiff(file='subdomains_figures/protM_Apsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Apoptosis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()

# no APP Metabolism subdomains in Male proteomics results
# AM <- subset(enrMerge2, enrMerge2$Biodomain=='APP Metabolism')
# enrMerge <- AM
# #chart subdomains individually:
# tiff(file='subdomains_figures/protM_AMsubs.tiff',height=110,width=300,units='mm',res=300)
# enrMerge %>%   
#   ggplot(aes( factor(state), NES ))+ 
#   scale_color_identity()+ scale_fill_identity()+
#   scale_size_continuous(
#     limits = -log10(enrMerge$pvalue) %>% range(),
#     range = c(.7,5))+
#   geom_jitter(
#     data = subset(enrMerge, pvalue > 0.05),
#     color = 'grey85', alpha = .2, size = .5, show.legend = F
#   )+
#   geom_jitter(
#     data = subset(enrMerge, pvalue < 0.05 ),
#     aes(color = color, size = -log10(pvalue) ), #
#     alpha = .5, show.legend = T
#   )+
#   geom_violin(
#     data = subset(enrMerge, pvalue < 0.05 & NES < 0),
#     scale='width', aes(fill = color), color = 'grey50',#, color = col
#     alpha = .3, show.legend = F
#   )+
#   geom_violin(
#     data = subset(enrMerge, pvalue < 0.05 & NES > 0),
#     scale='width', aes(fill = color), color = 'grey50',#, color = col
#     alpha = .3, show.legend = F
#   )+
#   geom_hline(yintercept = 0, lty = 3)+
#   facet_wrap(~Subdomain, ncol=5)+
#   labs(x='')+  coord_flip() + theme(legend.position = 'right') +
#   ggtitle('Male Proteomics, APP Metabolism Subdomains; term pvalue < 0.05'
#   )+ xlab("State") + scale_x_discrete(limits = rev)
# dev.off()


Au <- subset(enrMerge2, enrMerge2$Biodomain=='Autophagy')
enrMerge <- Au
#chart subdomains individually:
tiff(file='subdomains_figures/protM_Ausubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Autophagy Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


CC <- subset(enrMerge2, enrMerge2$Biodomain=='Cell Cycle')
enrMerge <- CC
#chart subdomains individually:
tiff(file='subdomains_figures/protM_CCsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Cell Cycle Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



DR <- subset(enrMerge2, enrMerge2$Biodomain=='DNA Repair')
# enrMerge <- DR
# #chart subdomains individually:
# tiff(file='subdomains_figures/protM_DRsubs.tiff',height=110,width=300,units='mm',res=300)
# enrMerge %>%   
#   ggplot(aes( factor(state), NES ))+ 
#   scale_color_identity()+ scale_fill_identity()+
#   scale_size_continuous(
#     limits = -log10(enrMerge$pvalue) %>% range(),
#     range = c(.7,5))+
#   geom_jitter(
#     data = subset(enrMerge, pvalue > 0.05),
#     color = 'grey85', alpha = .2, size = .5, show.legend = F
#   )+
#   geom_jitter(
#     data = subset(enrMerge, pvalue < 0.05 ),
#     aes(color = color, size = -log10(pvalue) ), #
#     alpha = .5, show.legend = T
#   )+
#   geom_violin(
#     data = subset(enrMerge, pvalue < 0.05 & NES < 0),
#     scale='width', aes(fill = color), color = 'grey50',#, color = col
#     alpha = .3, show.legend = F
#   )+
#   geom_violin(
#     data = subset(enrMerge, pvalue < 0.05 & NES > 0),
#     scale='width', aes(fill = color), color = 'grey50',#, color = col
#     alpha = .3, show.legend = F
#   )+
#   geom_hline(yintercept = 0, lty = 3)+
#   facet_wrap(~Subdomain, ncol=5)+
#   labs(x='')+  coord_flip() + theme(legend.position = 'right') +
#   ggtitle('Male Proteomics, DNA Repair Subdomains; term pvalue < 0.05'
#   )+ xlab("State") + scale_x_discrete(limits = rev)
# dev.off()


En <- subset(enrMerge2, enrMerge2$Biodomain=='Endolysosome')
enrMerge <- En
#chart subdomains individually:
tiff(file='subdomains_figures/protM_Ensubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Endolysosome Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Ep <- subset(enrMerge2, enrMerge2$Biodomain=='Epigenetic')
enrMerge <- Ep
#chart subdomains individually:
tiff(file='subdomains_figures/protM_Epsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Epigenetic Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



IR <- subset(enrMerge2, enrMerge2$Biodomain=='Immune Response')
enrMerge <- IR
#chart subdomains individually:
tiff(file='subdomains_figures/protM_IRsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Immune Response Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


LM <- subset(enrMerge2, enrMerge2$Biodomain=='Lipid Metabolism')
enrMerge <- LM
tiff(file='subdomains_figures/protM_LMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Lipid Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


MB <- subset(enrMerge2, enrMerge2$Biodomain=='Metal Binding and Homeostasis')
enrMerge <- MB
#chart subdomains individually:
tiff(file='subdomains_figures/protM_MBsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Metal Binding & Homeostasis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()

MM <- subset(enrMerge2, enrMerge2$Biodomain=='Mitochondrial Metabolism')
enrMerge <- MM
tiff(file='subdomains_figures/protM_MMsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Mitochondrial Metabolism Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



My <- subset(enrMerge2, enrMerge2$Biodomain=='Myelination')
enrMerge <- My
tiff(file='subdomains_figures/protM_Mysubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Myelination Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Ox <- subset(enrMerge2, enrMerge2$Biodomain=='Oxidative Stress')
enrMerge <- Ox
tiff(file='subdomains_figures/protM_Oxsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Oxidation Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



Pr <- subset(enrMerge2, enrMerge2$Biodomain=='Proteostasis')
enrMerge <- Pr
tiff(file='subdomains_figures/protM_Prsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Proteostasis Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()



SS <- subset(enrMerge2, enrMerge2$Biodomain=='Structural Stabilization')
enrMerge <- SS
tiff(file='subdomains_figures/protM_SSsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Structural Stabilization Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()




Sy <- subset(enrMerge2, enrMerge2$Biodomain=='Synapse')
enrMerge <- Sy
#present this a different way:
tiff(file='subdomains_figures/protM_SYsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Synapse Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()


Va <- subset(enrMerge2, enrMerge2$Biodomain=='Vasculature')
enrMerge <- Va
tiff(file='subdomains_figures/protM_VAsubs.tiff',height=110,width=300,units='mm',res=300)
enrMerge %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enrMerge$pvalue) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enrMerge, pvalue > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enrMerge, pvalue < 0.05 ),
    aes(color = color, size = -log10(pvalue) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enrMerge, pvalue < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Subdomain, ncol=5)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male Proteomics, Vasculature Subdomains; term pvalue < 0.05'
  )+ xlab("State") + scale_x_discrete(limits = rev)
dev.off()

