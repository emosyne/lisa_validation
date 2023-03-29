process R_prepare_lists_for_clump {
    // debug true
    container 'emosyne/r_docker:1.96'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "${cohort}_${ENH_list}"
    // cache "lenient"
    

    input:
    // [clz2a, clz2a_GWAS_QC_noclump.gz, /home/osimoe/PGC_w3_data/clz2a, 34k_neg, ./input/enh_bedfiles/34k_neg.bed]
    tuple val(cohort), path (LOO_GWAS_QC), path(cohort_dir), val(ENH_list), path(ENH_bed)


    output:
    tuple val(cohort), path (LOO_GWAS_QC),  val(ENH_list), path("*_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), emit: lists_before_clump

    
    script:
    """
    R_prepare_lists_for_clump.R $task.cpus ${ENH_list} ${ENH_bed}  ${LOO_GWAS_QC} ${cohort}
    
   
    """
}
    
    