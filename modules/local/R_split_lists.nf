process R_split_lists {
    // debug true
    container 'emosyne/r_docker:1.96'
    label 'process_high'
    tag "${cohort}_${ENH_list}"
    cache "lenient"
    errorStrategy 'ignore'
    

    input:
    //[celso, celso_GWAS_QC.gz, NEURAL_8k_GRB_significant_EPs, 
        // celso_NEURAL_8k_GRB_significant_EPs_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, celso_NEURAL_8k_GRB_significant_EPs_PGC__noclump_residual_GWAS_compartment.tsv.gz, celso_PGC_clumped_SNPs.clumped]
    tuple val(cohort), path (LOO_GWAS_QC),  val(ENH_list), \
        path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment),  path(clumped_SNPs), val(multiplier)
    each path(EP_ES_gene_brain_exp)


    output:
    tuple val(cohort), val(ENH_list), path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz"), val(multiplier),      emit: partitioned
    
    
    
    script:
    """
    R_split_lists.R "${cohort}_${ENH_list}" ${clumped_SNPs} ${noclump_residual_GWAS_compartment} ${noclump_TS_ENH_GWAS_compartment} ${EP_ES_gene_brain_exp} ${multiplier}
    
    
   
    """
}
    