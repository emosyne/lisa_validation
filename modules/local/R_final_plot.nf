process R_final_plot {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_high'
    tag "$cohort_ENHpart"
    cache "lenient"
    // errorStrategy 'ignore'

    input: 
    // // [xs234_Neural_significant_enh_ADD, 
        // xs234_Neural_significant_enh_0.5_ADD_clumped_EPWAS_originalOR.summary, xs234_Neural_significant_enh_0.5_ADD_mult_1_clumped_EPWAS_OR_by_measure1.summary, xs234_Neural_significant_enh_0.5_ADD_mult_1_clumped_EPWAS_OR_by_measure2.summary,
        // xs234_Neural_significant_enh_0.5_ADD_clumped_EPWAS_originalOR.prsice, xs234_Neural_significant_enh_0.5_ADD_mult_1_clumped_EPWAS_OR_by_measure1.prsice, xs234_Neural_significant_enh_0.5_ADD_mult_1_clumped_EPWAS_OR_by_measure2.prsice, 
        // xs234_Neural_significant_enh_0.5_ADD_clumped_EPWAS_originalOR.best, xs234_Neural_significant_enh_0.5_ADD_mult_1_clumped_EPWAS_OR_by_measure1.best, xs234_Neural_significant_enh_0.5_ADD_mult_1_clumped_EPWAS_OR_by_measure2.best, 
        // xs234_Neural_significant_enh_0.5_ADD_clumped_residual_GWAS_compartment.summary, xs234_Neural_significant_enh_0.5_ADD_clumped_residual_GWAS_compartment.prsice, xs234_Neural_significant_enh_0.5_ADD_clumped_residual_GWAS_compartment.best, 
        // xs234_Neural_significant_enh_0.5_ADD_clumped_LOO_GWAS.summary, xs234_Neural_significant_enh_0.5_ADD_clumped_LOO_GWAS.prsice, xs234_Neural_significant_enh_0.5_ADD_clumped_LOO_GWAS.best, 
        // Neural_significant_enh, 0.5, xs234_QC.fam, 
        // ADD, enh_ES, enh_TS_tpm]
    tuple val(cohort_ENHpart), \
        path(EPWAS_originalOR_summary), path(EPWAS_OR_by_measure1_summary), path(EPWAS_OR_by_measure2_summary),  \
        path(EPWAS_originalOR_prsice), path(EPWAS_OR_by_measure1_prsice), path(EPWAS_OR_by_measure2_prsice),  \
        path(EPWAS_originalOR_best), path(EPWAS_OR_by_measure1_best), path(EPWAS_OR_by_measure2_best),  \
        path(residual_GWAS_compartment_summary), path(residual_GWAS_compartment_prsice), path (residual_GWAS_compartment_best), \
        path(original_GWAS_summary), path(original_GWAS_prsice), path (original_GWAS_best),\
        val(ENH_list), val(CTthreshold), path(cohort_fam),\
        val(EPWAS_model), val(modif_name_1),val(modif_name_2)

    output:
    path("*/*.pdf")
    path("*/*.log")

    script:
    """
    
    R_final_plot.R $task.cpus ${ENH_list} ${cohort_fam} \
        ${EPWAS_originalOR_summary} ${EPWAS_originalOR_best}\
        ${EPWAS_OR_by_measure1_summary} ${EPWAS_OR_by_measure1_best}\
        ${EPWAS_OR_by_measure2_summary} ${EPWAS_OR_by_measure2_best}\
        ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
        ${EPWAS_originalOR_prsice} ${EPWAS_OR_by_measure1_prsice} ${EPWAS_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice}   \
        ${original_GWAS_summary} ${original_GWAS_prsice} ${original_GWAS_best}\
        ${modif_name_1} ${modif_name_2} ${CTthreshold} ${condition} ${ENHlist_thresh_model} 
    """
}
    