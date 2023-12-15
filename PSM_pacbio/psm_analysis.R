library(SequenceTools)

process_exonerate_out = function(Files){
      seqs_df = do.call(rbind,lapply(Files, function(f){
	    seqs = read_alignments(f, naming = 'headers')
	    seqs_df = do.call(rbind,strsplit(names(seqs),'\t')) %>% as_tibble()
	    colnames(seqs_df) = c('ID','length','cigar')
	    seqs_df %>% mutate(missmatch = as.integer(gsub('^.*:','',ID)), ID = gsub('[ ].*$','',ID), length = as.integer(gsub('length:','',length)), cigar = gsub('cigar:|cigar: ','',cigar), File = gsub('.fa','',basename(f)), Seq = seqs)
      }))

      examp = seqs_df$cigar
      exmatches = str_match_all(examp, '([A-Z]) ([0-9]+)' )

      indels = do.call(rbind,lapply(1:length(exmatches), function(xi){
	    x = exmatches[[xi]]
	    ins = grep('I',x[,2])
	    del = grep('D',x[,2])
	    indels = c(ins,del)
	    if(length(indels != 0)){
		  xr = x[,3]
		  xr[ins] = 0
		  #xr[del] = 1
		  xa = as.integer(x[,3])
		  xa[del] = 0
		  refpos = do.call(c,lapply(indels, function(y){
			sum(as.integer(xr[1:y]))
		  }))
		  alnpos = do.call(c,lapply(indels, function(y){
			sum(as.integer(xa[1:y]))
		  }))
		  data.frame(ReadNr = xi , Type = x[indels,2] , RefPos = refpos, AlnPos = alnpos, Length = x[indels,3])
	    } else {
		  data.frame(ReadNr = xi , Type = NA , RefPos = NA, AlnPos = NA, Length = 0)
	    }
      }))

      indels = indels %>% as_tibble() %>% mutate(Length = as.integer(ifelse(is.na(Length),0,Length)))
      mutations = inner_join(seqs_df %>% mutate(ReadNr = 1:length(seqs_df$ID)), indels, by = 'ReadNr') 

      return(mutations)
}


mutations = do.call(bind_rows,lapply(c('Brec1', 'D7L', 'D7R', 'Cre'), function(x){
      process_exonerate_out(paste0('psm/', x, '.fa'))
}))

fmut = mutations %>% filter(Type == 'D',Length == 15) %>% select(ID, File, Pos = AlnPos) %>% mutate(Pos = Pos + 1, ResPos = floor(Pos/3))
write_csv(fmut, 'psm_counts.csv')

fmut %>% group_by(File,ResPos) %>% count() %>% group_by(File) %>% mutate(Freq = n/sum(n), label = ifelse(rank(-Freq) < 11, ResPos, NA)) %>% ggplot(aes(x = ResPos, y = n, label = label)) + geom_col() + geom_label() + scale_x_continuous(breaks = seq(0,340,10)) + facet_wrap(~File)
ggsave('psm_counts.png')
