library(tidyverse)
library(SequenceTools)

dat1 = read_tsv('filtered_zfscreen2.tsv') %>% mutate(Ara = 1)
dat200 = read_tsv('filtered_zfscreen3.tsv') %>% mutate(Ara = 200)

# join and filter, change Brec1 to fit convention
dat = bind_rows(dat1,dat200) %>% mutate(EditDistance = as.integer(gsub('NM:i:','',EditDistance))) %>% filter(EditDistance/Length < 0.08)

# categories
splitmatch = do.call(rbind, strsplit(gsub('-','_',dat$Match),'_'))
colnames(splitmatch) = c('Znf_pos','Linker','TS','Znf_ts','State')

ratios = bind_cols(dat,as_tibble(splitmatch)) %>% mutate(Znf_orientation = ifelse(grepl('REV',Znf_ts),'Reverse','Forward'), Znf_ts = gsub('REV|Zif','',Znf_ts)) %>% group_by(Ara, Znf_pos, Linker, Znf_ts, Znf_orientation, State) %>% count() %>% group_by(Ara, Znf_pos, Linker, Znf_ts, Znf_orientation) %>% summarise(Reads = sum(n), Percentage = n[State == 'recombined']/Reads * 100) %>% ungroup() %>% mutate(its = as.integer(gsub('[A-Za-z]','',Znf_ts)), iL = as.integer(gsub('[A-Za-z]','',Linker))) %>% arrange(its,iL) %>% mutate(Znf_ts = fct_inorder(as.factor(Znf_ts)), Linker = fct_inorder(as.factor(Linker)), Percentage = round(Percentage,2))

write_csv(ratios, file = 'ratios.csv')

# Plot all together
ratios %>% filter(Znf_orientation == 'Reverse') %>% ggplot(aes(x = Linker, y = Znf_ts, fill = Percentage, label = Reads)) + geom_raster() + scale_fill_viridis_c(limits = c(0,100)) + facet_grid(cols = vars(Znf_pos), rows = vars(Znf_orientation), scale = 'free', space = 'free') + theme(axis.text.x = element_text(angle = 90))
