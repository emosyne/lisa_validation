include { bash_base_GWAS_QC }           from '../modules/local/bash_base_GWAS_QC.nf'
include { PLINK2_QC_PRUNE_HET }         from '../modules/local/PLINK2_QC_PRUNE_HET_mod.nf'
include { R_PRS_QC }                    from '../modules/local/R_PRS_QC_mod.nf'
include { PLINK_PRODUCE_QC_DATASET }    from '../modules/local/PLINK_PRODUCE_QC_DATASET_mod.nf'
include { R_prepare_lists_for_clump }   from '../modules/local/R_prepare_lists_for_clump.nf'
include {PLINK_clump}                   from '../modules/local/PLINK_clump_mod.nf'
include {R_split_lists}                 from '../modules/local/R_split_lists.nf'
include {PRSice_calculate_PRS_split_partitions} from '../modules/local/PRSice_calculate_PRS_split_partitions.nf'
include {R_final_plot}                  from '../modules/local/R_final_plot.nf'


test_sample_list = ["xs234"]//"celso""xirwt","xgras","xjrsa","gawli","mcqul","xclm2","xclo3","gpc2a","xboco","xswe5","xswe6","clz2a",
LOO_PGC_GWASES = 
    Channel.fromList( test_sample_list )
    .map{[it, file("/home/osimoe/sumstats_w3/autosomes/PRS/ancestry-specific/study_scores_deduped_eur/training/daner_PGC_SCZ_w3_76_0518d_eur.no${it}.gz")]}
   

// LOO_PGC_GWASES.view()


no_PCA_sample_list = ["butr","grtr","lemu","uktr"]
Channel.fromPath( '/home/osimoe/PGC_w3_data/*', type: 'dir' )
    .map{ it -> [it.simpleName , it] }
    .branch { sample_name, sample_path ->
        no_PCA: sample_name in no_PCA_sample_list
            return tuple( tuple( sample_name, sample_path ) )
        with_PCA: true//sample_name in test_sample_list
            return tuple( tuple( sample_name, sample_path ) )
    } \
    .set { inputs }


validation_samples = inputs.with_PCA
// validation_samples.view()


enhancer_lists_bed_files = 
    Channel.from("Neural_significant_enh")
        .map { ENH_list -> ["${ENH_list}", 
            file("./input/enh_bedfiles/${ENH_list}.bed", checkIfExists: true)]
            } 
// enhancer_lists_bed_files.view()

enhancer_EPWAS_files = 
    Channel.from("REC", "DOM", "ADD")
            .map { model -> ["${model}", 
                file("./input/EPWAS/UKBB_ENH_associations_${model}.tsv.gz", checkIfExists: true)]
            } 



//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("/home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect()





workflow lisa_validation {
    // BASE = LOO GWAS
    // TARGET = cohort

    LOOGWAS_and_validationSample = LOO_PGC_GWASES.join(validation_samples)
    // LOOGWAS_and_validationSample.first().view()

    // BASE (GWAS) QC: REMOVE LOW MAF AND INFO SCORES
    //produce GWAS_QC
    bash_base_GWAS_QC (
        LOOGWAS_and_validationSample
            .combine(LD_reference)
    )

    // TARGET QC 1: PRUNE AND HETEROZIGOSITY CALCULATIONS
    // produce prune.in and het files
    PLINK2_QC_PRUNE_HET (
        LOOGWAS_and_validationSample
    )

    

    // TARGET QC 2:  remove heterogeneity outliers, produced A1 alleles, and mismatching SNPs list to be removed
    // produce QC_het_a1_mismatch, 
    R_PRS_QC ( // calculates mismatching SNPs and recodes all alleles to GWAS base
        PLINK2_QC_PRUNE_HET.out.pruned_variants_het
            .join(bash_base_GWAS_QC.out.GWAS_QC_noclump)
            .join(bash_base_GWAS_QC.out.clumped_SNPs)
            .join(validation_samples)
        
    )
    // R_PRS_QC.out.QC_het_a1_mismatch.view()
    // R_PRS_QC.out.clumped_LOO_GWAS.view()
    
    // TARGET QC 3:  
    // Remove individuals with heterozigosity F coefficients that are more than 3 standard deviation (SD) units from the mean
    // also remove mismatching SNPs
    // also standard QC --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 
    PLINK_PRODUCE_QC_DATASET (
        R_PRS_QC.out.QC_het_a1_mismatch
            .join(validation_samples)

        //tuple path ("*.bed"), path ("*.bim"), path ("*.fam"), emit: target_QC
    )

    // // PLINK_PRODUCE_QC_DATASET.out.target_QC.view()

    
    
    
    bash_base_GWAS_QC.out.GWAS_QC_noclump
        .join(validation_samples)
        .combine(enhancer_lists_bed_files)
        .combine(enhancer_EPWAS_files)
        .map { it.flatten() }
        .set{cohort_GWAS_enh_list}
    
    // cohort_GWAS_enh_list.view()

    
    // BASE subsetting
    R_prepare_lists_for_clump (
        // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
        // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
        cohort_GWAS_enh_list
    )
    
    
    // R_prepare_lists_for_clump.out.lists_before_clump
    //     .combine(LD_reference)
    //     .view()    


    PLINK_clump (
        //CLUMPING of enhancer-based SNP compartments together 
        R_prepare_lists_for_clump.out.lists_before_clump
            .combine(LD_reference)
    )
    // PLINK_clump.out.clumped_SNPs_and_noclump_lists
    //     .view()


    R_split_lists (
        // first annotate SNPs with ES of relevant E-P - for ENH SNP list
        // ##################################################### GENERATE MODIFIED ORS MULT BY ES OR EXP       ###########################################################
        // output separate lists to calculate split PRSs and also merged one
        PLINK_clump.out.clumped_SNPs_and_noclump_lists.map { [it, "1"].flatten() }, //######################## multiplier can be set here ########################
        Channel.fromPath( "../private_input_files/ES_multipliers/2023-01-18_NEURAL_ENH_EXP_significant_ES_significant_contact_EPs_gene_brain_exp_plus_100_noOverlap.csv.gz", checkIfExists: true)
    )


    R_split_lists.out.partitioned
        .combine(PLINK_PRODUCE_QC_DATASET.out.target_QC, by: [0,0])//[celso, celso_QC.bed, celso_QC.bim, celso_QC.fam]
        .combine(validation_samples, by: [0,0])
        .combine(LD_reference)
        .combine(R_PRS_QC.out.clumped_LOO_GWAS, by: [0,0])
        .map { [it, "0.5"].flatten() }         // ######################## SET CT THRESHOLD FOR PRSICE ##################
        .set{combined_splitlists_bedfile_QCeddata_LDdata_05}
    R_split_lists.out.partitioned
        .combine(PLINK_PRODUCE_QC_DATASET.out.target_QC, by: [0,0])//[celso, celso_QC.bed, celso_QC.bim, celso_QC.fam]
        .combine(validation_samples, by: [0,0])
        .combine(LD_reference)
        .combine(R_PRS_QC.out.clumped_LOO_GWAS, by: [0,0])
        .map { [it, "0.05"].flatten() }         // ######################## SET CT THRESHOLD FOR PRSICE ##################
        .set{combined_splitlists_bedfile_QCeddata_LDdata_005}
    combined_splitlists_bedfile_QCeddata_LDdata = combined_splitlists_bedfile_QCeddata_LDdata_05.mix(combined_splitlists_bedfile_QCeddata_LDdata_005)
    
    // combined_splitlists_bedfile_QCeddata_LDdata.first().view()


    PRSice_calculate_PRS_split_partitions(
        combined_splitlists_bedfile_QCeddata_LDdata
    )
    
    // ########################################### CHANGE NAMES OF MULTIPLIERS ###########################################
    PRS_results = 
        PRSice_calculate_PRS_split_partitions.out.clumped_EPWAS_PRS
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_residual_GWAS_compartment_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_original_LOO_GWAS_PRS)
            .map { [it, "enh_ES", "enh_TS_tpm"].flatten() }


    // PRS_results.view()


    R_final_plot (
        PRS_results
    )


}


// 
// srun --pty -t 1-00:00 -p shared -c1 nextflow run https://github.com/emosyne/lisa_validation -latest -r master -profile lisa -resume 
// --slurmd-debug=error