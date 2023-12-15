library(tidyverse)
library(SequenceTools)

revratios = read_csv(file = 'ratios.csv')
forratios = read_csv(file = '../../ZF-screen1/analysis/ratios.csv') %>% mutate(Znf_orientation = 'Forward', Ara = ifelse(Znf_pos == 'middle', 200, 1), Znf_ts = gsub('Zif','',Znf_ts))

ratios_unfilt = bind_rows(revratios,forratios) %>% arrange(its,iL)
ratios = bind_rows(revratios %>% filter(Znf_orientation == 'Reverse'),forratios) %>% arrange(as.integer(ifelse(Znf_ts == 'wt', -1, Znf_ts))) %>% mutate(Znf_ts = fct_inorder(as.factor(Znf_ts))) %>% arrange(its,iL)

# Plot all together
ratios %>% ggplot(aes(x = Linker, y = Znf_ts, fill = Percentage, label = Percentage)) + geom_raster() + scale_fill_viridis_c(limits = c(0,100)) + facet_grid(cols = vars(Znf_pos), rows = vars(Znf_orientation), scale = 'free', space = 'free') + theme(axis.text.x = element_text(angle = 90))
ggsave('Znf_all.png')

# only middle
ratios %>% filter(Ara == 200, Znf_pos == 'middle') %>% mutate(Linker1 = as.integer(gsub('GGS|xGGS.*','',Linker)), Linker2 = as.integer(gsub('GGS[0-9]+xGGS','',Linker))) %>% ggplot(aes(x = Linker2, y = Linker1, fill = Percentage, label = Percentage)) + geom_raster() + scale_fill_viridis_c(limits = c(0,100)) + facet_grid(rows = vars(Znf_orientation), cols = vars(Znf_ts), scale = 'free') + labs(y = '[GGS] in first Linker', x = '[GGS] in second Linker', fill = 'Percentage\nRecombined', subtitle = 'Distance from recombinase TS to Zif')
ggsave('Znf_middle.png')

# C and N terminal
ratios %>% filter(Ara == 1, Znf_pos != 'middle') %>% mutate(Linker = as.integer(gsub('GGS','',Linker))) %>% arrange(as.integer(gsub('[A-Za-z]','',Znf_ts)), Linker) %>% mutate(Znf_ts = fct_inorder(as.factor(gsub('Zif','',Znf_ts))), Linker = fct_inorder(as.factor(Linker))) %>% ggplot(aes(x = Linker, y = Znf_ts, fill = Percentage, label = Percentage)) + geom_raster() + scale_fill_viridis_c(limits = c(0,100)) + facet_grid(cols = vars(Znf_pos), rows = vars(Znf_orientation), scale = 'free', space = 'free') + labs(x = '[GGS] in Linker', y = 'Distance Recombinase-TS to Zif-TS')
ggsave('Znf_NC.png')
