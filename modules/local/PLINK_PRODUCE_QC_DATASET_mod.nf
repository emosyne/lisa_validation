process PLINK_PRODUCE_QC_DATASET {
    // debug true
    tag "$cohort"
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"


    input:
    tuple val(cohort), path (het_valid), path (a1_bim), path(mismatch), path(cohort_dir)

    output:
    tuple val(cohort), path ("*.bed"), path ("*.bim"), path ("*.fam"), emit: target_QC
    path ("*.log")


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    bedfile="${cohort_dir}/imputed/hardcall_genotypes/*.bed"
    bimfile="${cohort_dir}/imputed/hardcall_genotypes/*.bim"
    bedfile2=`echo \$bedfile | sed 's/.bed//'`
    covariates="${cohort_dir}/prin_comp/*.mds"

    echo \$bedfile2

    plink \\
        --bfile \$bedfile2 \\
        --make-bed \\
        --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 \\
        --a1-allele ${a1_bim} \\
        --keep ${het_valid} \\
        --exclude ${mismatch} \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out ${cohort}_QC
    
    #remove text before star from IID
    mv ${cohort}_QC.fam ${cohort}_QC.fam.bak
    cat ${cohort}_QC.fam.bak | sed 's/*/1/g' > ${cohort}_QC.fam
    

    """
}
