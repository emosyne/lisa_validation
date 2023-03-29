process R_split_lists {
    // debug true
    container 'emosyne/r_docker:1.96'
    label 'process_high'
    tag "${cohort}_${ENH_list}_${EPWAS_model}"
    cache "lenient"
    errorStrategy 'ignore'
    

    input:
    //[xs234, xs234_GWAS_QC_noclump.gz, Neural_significant_enh, 
        // xs234_DOM_Neural_significant_enh_noclump_EPWAS.tsv.gz, xs234_DOM_Neural_significant_enh_PGC__noclump_residual_GWAS_compartment.tsv.gz, xs234_Neural_significant_enh_PGC_clumped_SNPs.clumped, DOM]
    tuple val(cohort), path (LOO_GWAS_QC),  val(ENH_list), \
        path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment),  path(clumped_SNPs), val(EPWAS_model),\
        val(multiplier)
    each path(EP_ES_gene_brain_exp)


    output:
    tuple val(cohort), val(ENH_list), path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz"), val(multiplier),val(EPWAS_model),       emit: partitioned
    
    
    
    script:
    """
    R_split_lists.R "${cohort}_${ENH_list}" ${clumped_SNPs} ${noclump_residual_GWAS_compartment} ${noclump_TS_ENH_GWAS_compartment} ${EP_ES_gene_brain_exp} ${multiplier}
    
    
   
    """
}
    