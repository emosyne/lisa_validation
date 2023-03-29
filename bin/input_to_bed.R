#!/usr/bin/env Rscript

## PART ONE: IMPORT SNV AND GWAS FILES, MAKE BED FILES
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

args = commandArgs()

print(args)
#inputs 

print(args[9])

full_PGC_GWAS_hg19 <- fread(file = args[8])
collected_bed_files_for_enhancers <- read_lines(args[9])
clumped_GWAS_hg19 <- fread(file = args[10])



#outputs
(general_out <- "GWAS_SNPs_in_initial_bed_files.bed")





#read each SNP list and join to GWAS
totalbed = data.frame()
for(bedfile in collected_bed_files_for_enhancers){
  print(bedfile)
  (bed = fread(file = bedfile, select=c(1:4), col.names = c("seqnames","start","end","name")) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T))
  seqlevelsStyle(x = bed) <- "UCSC"
  (bed=as_tibble(bed))
  
  (totalbed = rbind(totalbed,bed))
  
}
sample_n(totalbed, 20)

#remove internal range dupls
(totalbed = makeGRangesFromDataFrame(totalbed, keep.extra.columns = T) %>% 
    GenomicRanges::reduce(drop.empty.ranges = T))

#extract GWAS SNPs wihtin ranges
print("GWAS head")
print(head(full_PGC_GWAS_hg19))
(full_PGC_GWAS_hg19 = makeGRangesFromDataFrame(full_PGC_GWAS_hg19, keep.extra.columns = T,
                                               seqnames.field = "CHR", start.field = "POS", 
                                               end.field = "POS"))
seqlevelsStyle(full_PGC_GWAS_hg19) <- "UCSC"
(full_PGC_GWAS_overlap_beds = subsetByOverlaps(x = full_PGC_GWAS_hg19, ranges = totalbed, type="any"))

(clumped_GWAS_hg19 = makeGRangesFromDataFrame(clumped_GWAS_hg19, keep.extra.columns = T,
                                               seqnames.field = "CHR", start.field = "POS", 
                                               end.field = "POS"))
seqlevelsStyle(clumped_GWAS_hg19) <- "UCSC"

(SNPs_to_extract = rbind(
  as_tibble(full_PGC_GWAS_overlap_beds),
  as_tibble(clumped_GWAS_hg19)
))


(SNPs_to_extract <- SNPs_to_extract %>% 
        group_by(SNP) %>% slice_head(n=1) %>% ungroup())

#make bed
SNPs_to_extract %>% select(seqnames,  start, SNP) %>%
  mutate(POS=start, score=".", strand=".") %>%
  relocate(seqnames, start, POS, SNP, score, strand) %>%
  fwrite(file=general_out, sep="\t", col.names=F)
