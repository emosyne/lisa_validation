process PLINK_clump {
    // debug true
    tag "${cohort}_${ENH_list}"
    label 'process_high'
    container 'emosyne/plink2:1.23'
    // cache "lenient"

    input:
   //  tuple tuple val(cohort), path (LOO_GWAS),  path(cohort_dir), val(ENH_list), path("*_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), emit: lists_before_clump
    tuple val(cohort), path (LOO_GWAS_QC), val(ENH_list), path(PGC_noclump_TS_ENH_GWAS_compartment),  path(PGC_noclump_residual_GWAS_compartment), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)



    output:
    tuple val(cohort), path (LOO_GWAS_QC), val(ENH_list), path(PGC_noclump_TS_ENH_GWAS_compartment),  path(PGC_noclump_residual_GWAS_compartment), path("*_PGC_clumped_SNPs.clumped"), emit: clumped_SNPs_and_noclump_lists
    path("*.log")


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    plink  \\
       --clump ${PGC_noclump_TS_ENH_GWAS_compartment},${PGC_noclump_residual_GWAS_compartment} \\
       --clump-p1 1 --clump-p2 1 \\
       --clump-kb 500 --clump-r2 0.1 \\
       --bfile ${LD_ref_bed.baseName} \\
       --out ${cohort}_${ENH_list}_PGC_clumped_SNPs  \\
       --threads $task.cpus \\
       --memory $mem_mb

    
    """
}
