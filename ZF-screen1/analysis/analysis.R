library(tidyverse)
library(SequenceTools)

dat = read_tsv('filtered_n50_brec1_znf.tsv')
ref = read_alignments('../references/references.fa')
refls = tibble(Match = names(ref), RefLength = str_length(ref))

# filter
dat = dat %>% mutate(EditDistance = as.integer(gsub('NM:i:','',EditDistance))) %>% filter(EditDistance/Length < 0.08) 

# categories
splitmatch = do.call(rbind,strsplit(gsub('-','_',dat$Match),'_'))
colnames(splitmatch) = c('Znf_pos','Linker','TS','Znf_ts','State')

ratios = bind_cols(dat,as_tibble(splitmatch)) %>% group_by(Znf_pos, Linker, Znf_ts,State) %>% count() %>% group_by(Znf_pos, Linker, Znf_ts) %>% summarise(Reads = sum(n), Percentage = n[State == 'recombined']/Reads * 100)

# Plot all together
ratios %>% mutate(its = as.integer(gsub('[A-Za-z]','',Znf_ts)), iL = as.integer(gsub('[A-Za-z]','',Linker))) %>% arrange(its,iL) %>% mutate(Znf_ts = fct_inorder(as.factor(Znf_ts)), Linker = fct_inorder(as.factor(Linker))) %>% ggplot(aes(x = Linker, y = Znf_ts, fill = Percentage)) + geom_raster() + scale_fill_viridis_c() + facet_grid(cols = vars(Znf_pos), scale = 'free', space = 'free') + theme(axis.text.x = element_text(angle = 90))

# export ratios
ratios %>% write_csv('ratios.csv')
