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

# add rev-ts
tsunrec = read_alignments('tspart_rev.fa')
names(tsunrec) = paste0(names(tsunrec),'_unrecombined')
tsrec = read_alignments('tspart_rev-recombined.fa')
names(tsrec) = paste0(names(tsrec),'_recombined')
tsseqs = c(tsrec,tsunrec)

# combine with ts
references = apply(expand_grid(recombinases,tsseqs), 1, function(x) paste0(x,collapse=''))


# nonrev controls
controlrecs = recombinases[names(recombinases) %in% c("Nterminal-GGS12","Cterminal-GGS12","middle-GGS8xGGS8")]

# add for-ts
tsunrec = read_alignments('tspart.fa')
names(tsunrec) = paste0(names(tsunrec),'_unrecombined')
tsrec = read_alignments('tspart-recombined.fa')
names(tsrec) = paste0(names(tsrec),'_recombined')
tsseqs2 = c(tsrec,tsunrec)

# combine with ts
creferences = apply(expand_grid(controlrecs,tsseqs2), 1, function(x) paste0(x,collapse=''))


#add end
ends = read_alignments('end.fa')
references = paste0(references,ends)
creferences = paste0(creferences,ends)

# add names
names(references) = apply(expand_grid(names(recombinases),names(tsseqs)), 1, function(x) paste0(x, collapse = '_'))
names(creferences) = apply(expand_grid(names(controlrecs),names(tsseqs2)), 1, function(x) paste0(x, collapse = '_'))

# all references
arefs = c(references, creferences)
narefs = names(arefs) %>% unique()
arefs = (arefs) %>% unique()
names(arefs) = narefs
#arefs %>% toupper(.) %>% {write_fasta(.,names(.),'../all_references.fa')}

refszf2 = arefs[grep('Nterminal|Cterminal|Brec1',names(arefs))] 
refszf2 %>% toupper(.) %>% {write_fasta(.,names(.),'../zfscreen2_references.fa')}
refszf3 = arefs[grep('middle|Brec1',names(arefs))] 
refszf3 %>% toupper(.) %>% {write_fasta(.,names(.),'../zfscreen3_references.fa')}



#### make bed files for regions of interest ####
tspos = str_length(refszf2)-str_length(ends)
recpos = ifelse(grepl('middle',names(refszf2)),250,800)

tsbed = tibble(chrom = names(refszf2), chromStart = tspos, chromEnd = tspos+1)
recbed = tibble(chrom = names(refszf2), chromStart = recpos, chromEnd = recpos+1)

write_tsv(tsbed,'../ts_zfscreen2.bed', col_names = F)
write_tsv(recbed,'../rec_zfscreen2.bed', col_names = F)

tspos = str_length(refszf3)-str_length(ends)
recpos = ifelse(grepl('middle',names(refszf3)),250,800)

tsbed = tibble(chrom = names(refszf3), chromStart = tspos, chromEnd = tspos+1)
recbed = tibble(chrom = names(refszf3), chromStart = recpos, chromEnd = recpos+1)

write_tsv(tsbed,'../ts_zfscreen3.bed', col_names = F)
write_tsv(recbed,'../rec_zfscreen3.bed', col_names = F)
