library(tidyverse)
library(SequenceTools)

# N & C terminal parts
cterm = read_alignments('Cterm.fa')
nterm = read_alignments('Nterm.fa')
nclinkers = read_alignments('NC-linkers.fa')

cterm_linked = paste0(cterm['rec'],nclinkers,cterm['ziftots'])
names(cterm_linked) = paste0('Cterminal-',names(nclinkers))
nterm_linked = paste0(nterm['zif268'],nclinkers,nterm['rectots'])
names(nterm_linked) = paste0('Nterminal-',names(nclinkers))

# zif in middle
middle = read_alignments('middle.fa')
middle_linkers = read_alignments('middle-linkers.fa')
middle_linkers_df = do.call(rbind,strsplit(middle_linkers,','))
middle_linkers_comb = expand.grid(middle_linkers_df[,1],middle_linkers_df[,2], stringsAsFactors = F)
middle_linkers_comb_names = apply(expand.grid( names(middle_linkers) ,rownames(middle_linkers_df), stringsAsFactors = F), 1, function(x) paste0(x,collapse = 'x'))

middle_linked = paste0(middle['pre-zif'],middle_linkers_comb[,1], middle['zif268'], middle_linkers_comb[,2], middle['post-zif'])
names(middle_linked) = paste0('middle-',middle_linkers_comb_names)

# combine all together
recombinases = c(nterm_linked,cterm_linked,middle_linked)

# add ts
tsunrec = read_alignments('tspart.fa')
names(tsunrec) = paste0(names(tsunrec),'_unrecombined')
tsrec = read_alignments('tspart-recombined.fa')
names(tsrec) = paste0(names(tsrec),'_recombined')
tsseqs = c(tsrec,tsunrec)

references = apply(expand_grid(recombinases,tsseqs), 1, function(x) paste0(x,collapse=''))

#add end
ends = read_alignments('end.fa')
references = paste0(references,ends)

# add names
names(references) = apply(expand_grid(names(recombinases),names(tsseqs)), 1, function(x) paste0(x, collapse = '_'))

write_fasta(toupper(references),names(references),'../references.fa')


#### make bed files for regions of interest ####
tspos = str_length(references)-str_length(ends)
recpos = ifelse(grepl('^middle.*',names(references)),250,800)

tsbed = tibble(chrom = names(references), chromStart = tspos, chromEnd = tspos+1)
recbed = tibble(chrom = names(references), chromStart = recpos, chromEnd = recpos+1)

write_tsv(tsbed,'../ts.bed')
write_tsv(recbed,'../rec.bed')
