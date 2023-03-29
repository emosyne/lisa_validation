process R_final_plot {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_high'
    tag "$cohort_ENHpart"
    cache "lenient"
    // errorStrategy 'ignore'

    input: 
    // [celso_NEURAL_8k_GRB_significant_EPs, 
    // celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, 
    // celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, 
    // celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, 
    // celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_residual_GWAS_compartment.summary, celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_residual_GWAS_compartment.prsice, celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_residual_GWAS_compartment.best, 
    // celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_merged_GWAS.summary, celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_merged_GWAS.prsice, celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_merged_GWAS.best, 
    // celso_QC.fam, 
    // celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_LOO_GWAS.summary, celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_LOO_GWAS.prsice, celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_LOO_GWAS.best, 
    // 0.5, 
    // enh_ES, enh_TS_tpm]
    tuple val(cohort_ENHpart), \
        path(TS_ENH_GWAS_compartment_originalOR_summary), path(TS_ENH_GWAS_compartment_OR_by_measure1_summary), path(TS_ENH_GWAS_compartment_OR_by_measure2_summary),  \
        path(TS_ENH_GWAS_compartment_originalOR_prsice), path(TS_ENH_GWAS_compartment_OR_by_measure1_prsice), path(TS_ENH_GWAS_compartment_OR_by_measure2_prsice),  \
        path(TS_ENH_GWAS_compartment_originalOR_best), path(TS_ENH_GWAS_compartment_OR_by_measure1_best), path(TS_ENH_GWAS_compartment_OR_by_measure2_best),  \
        path(residual_GWAS_compartment_summary), path(residual_GWAS_compartment_prsice), path (residual_GWAS_compartment_best), \
        path(merged_GWAS_summary), path(merged_GWAS_prsice), path (merged_GWAS_best),\
        path(cohort_fam),\
        path(original_LOO_GWAS_summary), path(original_LOO_GWAS_prsice), path (original_LOO_GWAS_best),\
        val(CTthreshold),\
        val(modif_name_1),val(modif_name_2)

    output:
    path("*/*.pdf")
    path("*/*.log")

    script:
    """
    
    R_final_plot.R $task.cpus ${cohort_ENHpart} ${cohort_fam} \
        ${TS_ENH_GWAS_compartment_originalOR_summary} ${TS_ENH_GWAS_compartment_originalOR_best}\
        ${TS_ENH_GWAS_compartment_OR_by_measure1_summary} ${TS_ENH_GWAS_compartment_OR_by_measure1_best}\
        ${TS_ENH_GWAS_compartment_OR_by_measure2_summary} ${TS_ENH_GWAS_compartment_OR_by_measure2_best}\
        ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
        ${merged_GWAS_summary} ${merged_GWAS_best}\
        ${TS_ENH_GWAS_compartment_originalOR_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure1_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice} ${merged_GWAS_prsice}  \
        ${original_LOO_GWAS_summary} ${original_LOO_GWAS_prsice} ${original_LOO_GWAS_best}\
        ${modif_name_1} ${modif_name_2} ${CTthreshold}
    """
}
    