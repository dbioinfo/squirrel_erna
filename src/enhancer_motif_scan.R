library(tidyverse)
library(universalmotif)

setwd("~/WorkForaging/Teaching/Mentorships/arianna/github_erna/squirrel_erna/")

#1. load enhancer sequences
enh_seqs <- readRDS("data/enhancer_sequenecs.rds")

gdat <- as.data.frame(enh_seqs@ranges)

g1 <- ggplot(gdat)+
  geom_density(aes(x=width))+
  geom_vline(xintercept=5000, color='red')+
  scale_x_log10()+
  theme_bw() + 
  xlab("Enhancer Width (nt)")+
  ggtitle("Enhancer Width Distribution")
g1

ggsave("figs/enh_distribution_cutoff.png", g1)


#2. trim out long enhancers
keepers <- gdat %>% filter(width<=5000) %>% pull(names)
enh_seqs_trim <- enh_seqs[keepers]


#3. set up memes
motifs <- read_meme("raw_data/jaspar_human_motifs.meme")
motif_lookup <- data.frame(motif    = sapply(motifs, function(m) m@name),
                           alt_name = sapply(motifs, function(m) m@altname))
resscan <- scan_sequences(
  motifs         = motifs,
  sequences      = enh_seqs_trim,
  threshold      = 0.8,
  threshold.type = "logodds",
  RC             = TRUE,
  nthreads=6)
resscan <- as.data.frame(resscan) %>% left_join(motif_lookup, by='motif')

write_csv(resscan, "data/motif_hits.csv")


