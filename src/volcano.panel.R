library(tidyverse)

setwd("/home/dylan/WorkForaging/Teaching/Mentorships/arianna/github_erna/squirrel_erna/")

edat <- read_tsv("results/enrsa_hap1/eRNA/Enhancer.txt")

ibadat <- read_tsv("results/enrsa_hap1/IBA_vs_SA/eRNA/Enhancer_change.txt") %>% mutate(Comparison="IBA")
entdat <- read_tsv("results/enrsa_hap1/Ent_vs_SA/eRNA/Enhancer_change.txt")  %>% mutate(Comparison="ENT")
etdat <- read_tsv("results/enrsa_hap1/ET_vs_SA/eRNA/Enhancer_change.txt")  %>% mutate(Comparison="ET")
ltdat <- read_tsv("results/enrsa_hap1/LT_vs_SA/eRNA/Enhancer_change.txt")  %>% mutate(Comparison="LT")

gdat <- rbind(rbind(ibadat, entdat), rbind(etdat, ltdat)) %>% 
  mutate(Comparison=factor(Comparison, levels=c("ENT", "ET", "LT", "IBA")))

#volcano panel
ggplot(gdat)+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), fill=-log10(padj)),shape=21, size=2)+
  scale_fill_gradientn(colors = c("black", "black", "darkred", "darkorange", "yellow"),
                       values = scales::rescale(c(0,-log10(.05), 2, 7, 15), c(0,1)),
                       oob=scales::squish)+
  geom_hline(yintercept=-log10(.05), color='red')+
  geom_vline(xintercept=1, color='red')+
  geom_vline(xintercept=-1, color='red')+
  facet_wrap(~Comparison, nrow=2)+
  theme_bw()


#filter DE data
filt_de <- gdat %>% filter(padj<.05, abs(log2FoldChange)>1)

library(GenomicRanges)
library(universalmotif)
