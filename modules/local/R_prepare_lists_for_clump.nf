process R_prepare_lists_for_clump {
    // debug true
    container 'emosyne/r_docker:1.96'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "${cohort}_${ENH_list}"
    // cache "lenient"
    

    input:
    // [xs234, xs234_GWAS_QC_noclump.gz, /home/osimoe/PGC_w3_data/xs234, Neural_significant_enh, /project/osimoe/.nextflow/assets/emosyne/lisa_validation/input/enh_bedfiles/Neural_significant_enh.bed, 
        //REC, /project/osimoe/.nextflow/assets/emosyne/lisa_validation/input/EPWAS/UKBB_ENH_associations_REC.tsv.gz] 
    tuple val(cohort), path (LOO_GWAS_QC), path(cohort_dir), val(ENH_list), path(ENH_bed),\
        val(EPWAS_model), path(ENH_EPwas)

    output:
    tuple val(cohort), path (LOO_GWAS_QC),  val(ENH_list), path("*_noclump_EPWAS.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), \
        val(EPWAS_model),  emit: lists_before_clump
    
    script:
    """
    R_prepare_lists_for_clump.R $task.cpus ${ENH_list} ${ENH_bed}  ${LOO_GWAS_QC} ${cohort}${EPWAS_model} ${ENH_EPwas}
    
   
    """
}
    
    