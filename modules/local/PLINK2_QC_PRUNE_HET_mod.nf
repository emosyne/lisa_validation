process PLINK2_QC_PRUNE_HET {
    tag "$cohort"
    // debug true
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input: 
    tuple val(cohort), path(LOO_GWAS), path(cohort_dir)
    

    output:
    tuple val(cohort), path ("*.prune.in"), path ("*.het"), emit: pruned_variants_het
    // path("*.log")
    


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
      --indep-pairwise 200 50 0.25 \\
      --out ${cohort} \\
      --threads $task.cpus \\
      --memory $mem_mb

    plink \\
      --bfile \$bedfile2 \\
      --extract ${cohort}.prune.in \\
      --het \\
      --out ${cohort} \\
      --threads $task.cpus \\
      --memory $mem_mb
    """
}

