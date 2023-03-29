process PLINK_clump {
    // debug true
    tag "${cohort}_${ENH_list}_${EPWAS_model}"
    label 'process_high'
    container 'emosyne/plink2:1.23'
    // cache "lenient"

    input:
    // [xs234, xs234_GWAS_QC_noclump.gz, Neural_significant_enh, xs234_ADD_Neural_significant_enh_noclump_EPWAS.tsv.gz, xs234_ADD_Neural_significant_enh_PGC__noclump_residual_GWAS_compartment.tsv.gz, 
        // ADD, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    tuple val(cohort), path (LOO_GWAS_QC), val(ENH_list), path(PGC_noclump_EPWAS),  path(PGC_noclump_residual_GWAS_compartment), \
        val(EPWAS_model), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)



    output:
    tuple val(cohort), path (LOO_GWAS_QC), val(ENH_list), path(PGC_noclump_EPWAS),  path(PGC_noclump_residual_GWAS_compartment), path("*_PGC_clumped_SNPs.clumped"),  val(EPWAS_model), emit: clumped_SNPs_and_noclump_lists
    path("*.log")


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    plink  \\
       --clump ${PGC_noclump_EPWAS},${PGC_noclump_residual_GWAS_compartment} \\
       --clump-p1 1 --clump-p2 1 \\
       --clump-kb 500 --clump-r2 0.1 \\
       --bfile ${LD_ref_bed.baseName} \\
       --out ${cohort}_${ENH_list}_PGC_clumped_SNPs  \\
       --threads $task.cpus \\
       --memory $mem_mb

    
    """
}
