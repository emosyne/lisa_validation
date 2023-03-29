process bash_base_GWAS_QC {
    tag "$cohort"
    // debug true
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input: 
    tuple val(cohort), path(LOO_GWAS), path(cohort_dir), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)
    

    output:
    tuple val(cohort), path ("*_GWAS_QC_noclump.gz"),                 emit: GWAS_QC_noclump
    tuple val(cohort), path ("*_GWAS_QC_clump.clumped"),                                emit: clumped_SNPs

    // path("*.log")
    


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    #remove SNPs with INFO < 0.8 and MAF < 0.01
    zcat ${LOO_GWAS} | awk 'NR==1 || (\$6 < 0.99) && (\$8 > 0.8) {print}' | gzip > ${cohort}_GWAS_QC_noclump.gz

    plink  \\
       --clump ${cohort}_GWAS_QC_noclump.gz \\
       --clump-p1 1 --clump-p2 1 \\
       --clump-kb 500 --clump-r2 0.1 \\
       --bfile ${LD_ref_bed.baseName} \\
       --out ${cohort}_GWAS_QC_clump  \\
       --threads $task.cpus \\
       --memory $mem_mb
    """
}

