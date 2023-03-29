process PRSice_calculate_PRS_split_partitions {
    // debug true
    tag "${cohort}_${ENH_list}"
        label 'process_high'
    // label 'process_high_memory'
    // clusterOptions "--partition=shared_52c_384g"
    container 'emosyne/prsice_gwama_exec:1.0'
    cache "lenient"
    // maxForks 5
    errorStrategy 'ignore'


    input:
     // [xclo3, 18k_PsychENCODE_PFCortex, 
        // xclo3_18k_PsychENCODE_PFCortex_clumped_TS_ENH_GWAS_compartment.tsv.gz, xclo3_18k_PsychENCODE_PFCortex_clumped_residual_GWAS_compartment.tsv.gz, xclo3_18k_PsychENCODE_PFCortex_clumped_merged_GWAS.tsv.gz, 
        // xclo3_QC.bed, xclo3_QC.bim, xclo3_QC.fam, 
        // /home/osimoe/PGC_w3_data/xclo3, 
        // /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam, 
        // xclo3_clumped_LOO_GWAS.tsv.gz, 0.05]
    
        
    tuple val(cohort), val(ENH_list), \
        path(clumped_TS_ENH_GWAS_compartment), path(clumped_residual_GWAS_compartment), path(clumped_merged_GWAS), val(multiplier), \
        path(cohort_bed_QC), path(cohort_bim_QC), path(cohort_fam_QC), \
        path(cohort_dir), \
        path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam), \
        path(clumped_LOO_GWAS), val(CTthreshold)
    
    output:
    tuple val("${cohort}_${ENH_list}"), path("*_clumped_TS_ENH_GWAS_compartment_*.summary"), path("*_clumped_TS_ENH_GWAS_compartment_*.prsice"), path("*_clumped_TS_ENH_GWAS_compartment_*.best"),          emit: clumped_TS_ENH_GWAS_compartment_PRS
    tuple val("${cohort}_${ENH_list}"), path("*_clumped_residual_GWAS_compartment.summary"), path("*_clumped_residual_GWAS_compartment.prsice"), path("*_clumped_residual_GWAS_compartment.best"),          emit: clumped_residual_GWAS_compartment_PRS
    tuple val("${cohort}_${ENH_list}"), path("*_clumped_merged_GWAS.summary"), path("*_clumped_merged_GWAS.prsice"), path("*_clumped_merged_GWAS.best"),  path(cohort_fam_QC),                              emit: clumped_merged_GWAS_PRS
    tuple val("${cohort}_${ENH_list}"), path("*_clumped_LOO_GWAS.summary"), path("*_clumped_LOO_GWAS.prsice"), path("*_clumped_LOO_GWAS.best"), val(CTthreshold),                                           emit: clumped_original_LOO_GWAS_PRS
    tuple  path("*.png"), path("*.txt"), path("*.log") //figures, quantiles text and log

    script:
    def mem_Gb = (task.memory * 0.95).toGiga()
    def max_cpus = Math.round(task.cpus * 4/5)
    """
    covariates="${cohort_dir}/prin_comp/*.mds"
    #remove text before star from IID
    cat \$covariates | sed 's/*/1/g' | awk '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13 }' > covariates.pheno
    #< \$covariates cut -d'*' -f2 | awk '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13 }' > covariates.pheno
    head covariates.pheno
    echo memory: ${mem_Gb}Gb
    echo cpus: $max_cpus
    
    
    echo clumped_TS_ENH_GWAS_compartment - ORIGINAL OR
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    # ORIGINAL OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${cohort}_${ENH_list}_${CTthreshold}_clumped_TS_ENH_GWAS_compartment_originalOR

    echo clumped_TS_ENH_GWAS_compartment  - OR by measure 1
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR_by_measure1 --or \\
        --out ${cohort}_${ENH_list}_${CTthreshold}_mult_${multiplier}_clumped_TS_ENH_GWAS_compartment_OR_by_measure1
    
    echo clumped_TS_ENH_GWAS_compartment  - OR by measure 2
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR_by_measure2 --or \\
        --out ${cohort}_${ENH_list}_${CTthreshold}_mult_${multiplier}_clumped_TS_ENH_GWAS_compartment_OR_by_measure2

    echo clumped_residual_GWAS_compartment - original OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_residual_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${cohort}_${ENH_list}_${CTthreshold}_clumped_residual_GWAS_compartment
        
    
    echo clumped_merged_GWAS - original OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_merged_GWAS} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${cohort}_${ENH_list}_${CTthreshold}_mult_${multiplier}_clumped_merged_GWAS
    
    echo ORIGINAL GWAS LOO
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_LOO_GWAS} \\
        --target ${cohort_bed_QC.simpleName} \\
        --ld ${LD_ref_bed.baseName} \\
        --no-clump  --score avg \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels ${CTthreshold} --no-full --fastscore \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${cohort}_${ENH_list}_${CTthreshold}_clumped_LOO_GWAS
    """

}
        