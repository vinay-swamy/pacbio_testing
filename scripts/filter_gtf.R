library(tidyverse)
# args <- c('~/NIH/','~/NIH/eyeintegration_splicing/testing/Retina_Adult.Tissue.combined.gtf',
#           '~/NIH/dev_eyeintegration_splicing/ref/gencode_comp_ano.gtf', 
#           '~/NIH/eyeintegration_splicing/testing/Retina_Adult.Tissue.tracking',
#           'study_all',
#           'Retina_Adult.Tissue', 
#           '/Volumes/data/eyeintegration_splicing/sampleTableV6.tsv', 
#           'tout.gtf')

source('~/scripts/write_gtf.R')
args <- commandArgs(trailingOnly = T)
wd <- args[1]
gtf_file <- args[2]
ref_gtf_file <- args[3]
tracking_file <- args[4]
t_tissue <- args[5]
sample_table_file <- args[6]
out_gtf_file <- args[7]
setwd(wd)
save(args, file='/tmp/args.Rdata')
nm_col <- function(col){
    col=col[col!='-']
    name=str_split(col[1], '\\.\\d|:|\\|')[[1]][2] %>% str_split('_MSTRG') %>% .[[1]] %>% .[1]
    return(name)
}




gtf <- rtracklayer::readGFF(gtf_file)
ref_gtf <- rtracklayer::readGFF(ref_gtf_file)
track_tab <- read_tsv(tracking_file, col_names=F)
names <- c('transcript_id', 'gene_id','refid','code', apply(track_tab[,-(1:4)], 2, nm_col))
colnames(track_tab) <- names
track_tab <- track_tab %>% mutate(refid=str_split(refid, '\\|') %>% sapply(function(x)x[2]))
det_df <- apply(track_tab[,-(1:4)],2, function(x) x!='-') %>% as.data.frame %>%  bind_cols(track_tab[,1:4],.)
num_det_df <-det_df %>% mutate(num_det=rowSums(det_df[,-(1:4)])) %>%   
    select(transcript_id, gene_id, refid, code,num_det) 
co <- 3
keep_codes <- c('=','+','c','k','m','n','j', 'u')
num_det_df_chess_kc <- filter(num_det_df, num_det >=co, code %in% keep_codes)
ref_gtf_tx <- ref_gtf %>% filter(type == 'transcript')%>% select(seqid, strand, start, end, refid=transcript_id)
gtf_tx <- gtf %>% filter(type == 'transcript') 
gffc_ref_absmatch <- gtf_tx %>% filter(class_code == '=') %>% inner_join(ref_gtf_tx) %>% pull(transcript_id)
det_df$code[det_df$code == '='] <- '*'
num_det_df$code[num_det_df$code == '='] <- '*'
num_det_df$code[num_det_df$transcript_id %in% gffc_ref_absmatch] <- '='
det_df$code[det_df$transcript_id %in% gffc_ref_absmatch] <- '='
icpres <- (gtf$class_code == '=') %>% replace_na(F)
matchref <- gtf$transcript_id %in% gffc_ref_absmatch
gtf[gtf == '='] <- '*'
gtf$class_code[icpres & matchref] <- '='
filt_gtf <- filter(gtf, transcript_id %in% num_det_df_chess_kc$transcript_id)    
write_gtf3(filt_gtf, out_gtf_file)


