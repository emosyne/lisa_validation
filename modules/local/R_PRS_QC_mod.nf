process R_PRS_QC {
    // debug true
    container 'emosyne/r_docker:1.96'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "$cohort"
    cache "lenient"
    


    input: 
    // [xclm2, xclm2.prune.in, xclm2.het, xclm2_GWAS_QC_noclump.gz, xclm2_GWAS_QC_clump.clumped, /home/osimoe/PGC_w3_data/xclm2]
    tuple val(cohort), path (prune), path (het), path(LOO_GWAS_QC), path(LOO_GWAS_QC_clumped_SNPs), path(cohort_dir)
    
    

    output:
    tuple val(cohort), path ("*_het_valid_out_vs_LOO_GWAS*.sample"), path("*_a1_cohort_bim_vs_LOO_GWAS*"), path("*_mismatching_SNPs_vs_LOO_GWAS*"),  emit: QC_het_a1_mismatch
    tuple val(cohort), path("*_clumped_LOO_GWAS.tsv.gz"),  emit: clumped_LOO_GWAS
    // path("*fullGWAS_merged_with_EP_lists_vs_UKBB_annotated_mismatching.csv")    

    
    script:
    """
    bimfile="${cohort_dir}/imputed/hardcall_genotypes/*.bgn.bim"

    R_PRS_QC2.R ${het} \$bimfile ${LOO_GWAS_QC}  ${cohort} ${LOO_GWAS_QC_clumped_SNPs}
    
    """
}
