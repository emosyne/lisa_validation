#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

#set max CPU processes
nthreads = as.numeric(args[8])
setDTthreads( round(nthreads))

(ENH_list = args[9])
(ENH_bed_file = args[10])
(LOO_GWAS = args[11])
(cohort = args[12])
(EPWAS_model = args[13])
(ENH_EPwas = args[14])

#OUTPUT
TS_EPs_outfilename = paste0(cohort, "_", EPWAS_model, "_", ENH_list, "_noclump_EPWAS.tsv.gz")
residual_GWAS_compartment_outfilename = paste0(cohort, "_", EPWAS_model, "_", ENH_list, "_PGC__noclump_residual_GWAS_compartment.tsv.gz")



pDivBy = 100000




(ENH_bed = fread(file = ENH_bed_file, select=c(1:4), col.names = c("seqnames","start","end","name")) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T))
seqlevelsStyle(x = ENH_bed) <- "UCSC"

#extract GWAS SNPs within ranges
print("GWAS head")
# Read in the GWAS data
             #CHR	SNP	BP	A1	A2	FRQ_A_51393	FRQ_U_75752	INFO	OR	SE	P	ngt	Direction	HetISqt	HetDf	HetPVa	Nca	Nco	Neff
LOO_GWAS <- fread(file=LOO_GWAS, select = c("CHR", "SNP", "BP", "A1", "A2",  "P", "OR")) %>%
    dplyr::select(CHR,SNP,POS=BP,A1,A2,P,OR)
(LOO_GWAS_hg19 = makeGRangesFromDataFrame(LOO_GWAS, keep.extra.columns = T,
                                               seqnames.field = "CHR", start.field = "POS", 
                                               end.field = "POS"))
seqlevelsStyle(LOO_GWAS_hg19) <- "UCSC"

# Read in the EPWAS
(ENH_EPwas <- fread(file=ENH_EPwas, select = c("CHR", "SNP", "POS", "A1", "A2",  "P", "BETA")) %>%
    dplyr::select(CHR,SNP,POS,A1,A2,P,BETA) %>%
    mutate(OR=exp(BETA), BETA=NULL)
    ) 
(ENH_EPwas_GR = makeGRangesFromDataFrame(ENH_EPwas, keep.extra.columns = T,
                                               seqnames.field = "CHR", start.field = "POS", 
                                               end.field = "POS"))
seqlevelsStyle(ENH_EPwas_GR) <- "UCSC"

#subset GWAS by bed
################################################################################################################ CAN DIVIDE P BY VALUE TO PRIORITISE ENH SNPS ########################################################
(ENH_EPwas_GR_overlap_ENH = subsetByOverlaps(x = ENH_EPwas_GR, ranges = ENH_bed, type="any"))
(ENH_EPwas_GR_overlap_ENH = as_tibble(ENH_EPwas_GR_overlap_ENH) %>%
  dplyr::select(CHR=seqnames, SNP, POS=start, A1, A2, P, OR) %>%
  dplyr::mutate(P=P/pDivBy)
  ) 


fwrite(x= ENH_EPwas_GR_overlap_ENH, file = TS_EPs_outfilename, sep="\t")
#R.utils::gzip(TS_EPs_outfilename,destname=paste0(TS_EPs_outfilename, ".gz"))


#merge annotated_OR_E_Ps and original_base, map to LD blocks, clump, separate 2 lists
#merge annotated_OR_E_Ps and original_base
(full_GWAS_NOoverlap_beds = subsetByOverlaps(x = LOO_GWAS_hg19, ranges = ENH_bed, invert = TRUE))
(full_GWAS_NOoverlap_beds = as_tibble(full_GWAS_NOoverlap_beds) %>%
  dplyr::select(CHR=seqnames, SNP, POS=start, A1, A2, P, OR))

fwrite(x = full_GWAS_NOoverlap_beds, file = residual_GWAS_compartment_outfilename, sep="\t")