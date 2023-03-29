include { bash_base_GWAS_QC }           from '../modules/local/bash_base_GWAS_QC.nf'
include { PLINK2_QC_PRUNE_HET }         from '../modules/local/PLINK2_QC_PRUNE_HET_mod.nf'
include { R_PRS_QC }                    from '../modules/local/R_PRS_QC_mod.nf'
include { PLINK_PRODUCE_QC_DATASET }    from '../modules/local/PLINK_PRODUCE_QC_DATASET_mod.nf'
include { R_prepare_lists_for_clump }   from '../modules/local/R_prepare_lists_for_clump.nf'
include {PLINK_clump}                   from '../modules/local/PLINK_clump_mod.nf'
include {R_split_lists}                 from '../modules/local/R_split_lists.nf'
include {PRSice_calculate_PRS_split_partitions} from '../modules/local/PRSice_calculate_PRS_split_partitions.nf'
include {R_final_plot}from '../modules/local/R_final_plot.nf'

// // chain file
// hg38ToHg19_chain = Channel
//     .fromPath( "./input/chainfiles/hg38ToHg19.over.chain", checkIfExists: true)
// // SCZ wave 3 full GWAS merged to LD blocks 1000 genomes
// PGC_GWAS_GW_LD_blocks_merge_hg19 = Channel
//     .fromPath( "./input/LD/PGC_Jul22_GWAS_LD_block.tsv.gz", checkIfExists: true)

// // EP annotations
// enhancer_annotations_hg19 = Channel
//     .fromPath("./input/initial_SNP_lists/annotated_EPs.csv.gz", checkIfExists: true)
// // SCZ wave 3 full GWAS
// full_GWAS_hg19 = Channel
//     .fromPath("http://data.genereg.net/emanuele/GWAS_PGC_summary/fullGWAS_SCZ_PGC3_SCZ_wave3.european.autosome.public.v3_tidy_hg19.tsv.gz", checkIfExists: true)





// INPUTS
// // scz_clz2a_eur_sr-qc
// scz_xclm2_eur_sr-qc
// scz_xmgs2_eur_sr-qc
// scz_xclo3_eur_sr-qc
// scz_xs234_eur_sr-qc
// scz_celso_eur_sr-qc
test_sample_list = ["clz2a","xs234","celso"]//"celso""xirwt","xgras","xjrsa","gawli","mcqul","xclm2","xclo3","gpc2a","xboco","xswe5","xswe6",
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
    Channel.from(
        "PsychENCODE_PFCortex_enh", 
        "Neural_significant_enh",
        "Neural_significant_enh_GRB",
        // "NEURAL_14k_noGRB_significant_EPs",
        "Non-neural_enh",
        // "687_BRAIN_EP_eQTL",
        "Non-associated_enh"
        )
        .map { ENH_list -> ["${ENH_list}", 
            file("./input/enh_bedfiles/${ENH_list}.bed", checkIfExists: true)]
            } 
    
// enhancer_lists_bed_files.view()




//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("/home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect()

// UKBB_associations = Channel
//     .fromPath("input/UKBB_res/REC_ALL_BRAIN_EPs_UKBB_eur_associations_GWAMA_format")
//     // .fromPath( ["./input/UKBB_res/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid",
//                 // "./input/UKBB_res/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.frq"])
//     // .collect()
//     // .map{ it -> ["UKBB_eur" , it[0], it[1]] }




workflow lisa_percohort_devel {
    // BASE = LOO GWAS
    // TARGET = cohort

    input = LOO_PGC_GWASES.join(validation_samples)
    // input.first().view()

    // BASE (GWAS) QC: REMOVE LOW MAF AND INFO SCORES
    //produce GWAS_QC
    bash_base_GWAS_QC (
        input
            .combine(LD_reference)
    )

    // TARGET QC 1: PRUNE AND HETEROZIGOSITY CALCULATIONS
    // produce prune.in and het files
    PLINK2_QC_PRUNE_HET (
        input
    )

    // PLINK2_QC_PRUNE_HET.out.pruned_variants_het
    //     .join(bash_base_GWAS_QC.out.GWAS_QC_noclump)
    //     .join(bash_base_GWAS_QC.out.clumped_SNPs)
    //     .join(validation_samples)
    //     .view()
    // [xclm2, xclm2.prune.in, xclm2.het, xclm2_GWAS_QC_noclump.gz, xclm2_GWAS_QC_clump.clumped, /home/osimoe/PGC_w3_data/xclm2]
    // [xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bb/078d61c85c37e7b62f30494053b420/xs234.prune.in, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bb/078d61c85c37e7b62f30494053b420/xs234.het, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/c4/94a2c72d2ddc8f55e829406ec41004/xs234_GWAS_QC_noclump.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/c4/94a2c72d2ddc8f55e829406ec41004/xs234_GWAS_QC_clump.clumped, /home/osimoe/PGC_w3_data/xs234]
    // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/0e/ce0ddf79c4e1e43923154fa08368cc/celso.prune.in, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/0e/ce0ddf79c4e1e43923154fa08368cc/celso.het, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/38/f80cdfd18d5c37ce98e89a84ada3e8/celso_GWAS_QC_noclump.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/38/f80cdfd18d5c37ce98e89a84ada3e8/celso_GWAS_QC_clump.clumped, /home/osimoe/PGC_w3_data/celso]
    // [xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/68/4af6fe63f28e14128300284d41f8a9/xclo3.prune.in, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/68/4af6fe63f28e14128300284d41f8a9/xclo3.het, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/0f/49f152fb2410d564189ea2c68517eb/xclo3_GWAS_QC_noclump.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/0f/49f152fb2410d564189ea2c68517eb/xclo3_GWAS_QC_clump.clumped, /home/osimoe/PGC_w3_data/xclo3]
    // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ff/060e35ff0fa0dd4990a129e665205c/clz2a.prune.in, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ff/060e35ff0fa0dd4990a129e665205c/clz2a.het, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/18/c4dbd8ba39c2afa0c2ff01ebc6ea1e/clz2a_GWAS_QC_noclump.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/18/c4dbd8ba39c2afa0c2ff01ebc6ea1e/clz2a_GWAS_QC_clump.clumped, /home/osimoe/PGC_w3_data/clz2a]

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
        .map { it.flatten() }
        .set{cohort_GWAS_enh_list}
    
    // cohort_GWAS_enh_list.view()
    // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/18/c4dbd8ba39c2afa0c2ff01ebc6ea1e/clz2a_GWAS_QC_noclump.gz, /home/osimoe/PGC_w3_data/clz2a, 34k_neg, ./input/enh_bedfiles/34k_neg.bed]
    
    // BASE subsetting
    R_prepare_lists_for_clump (
        // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
        // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
        cohort_GWAS_enh_list
    )
    
    
    // R_prepare_lists_for_clump.out.lists_before_clump
    //     .combine(LD_reference)
    //     .view()
    // val(cohort), path (LOO_GWAS),  val(ENH_list), path("*_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), emit: lists_before_clump
    // [xirwt, daner_PGC_SCZ_w3_76_0518d_eur.noxirwt.gz, PsychENCODE_DER_03b_PFC_enhancers_18k, xirwt_PsychENCODE_DER_03b_PFC_enhancers_18k_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, xirwt_PsychENCODE_DER_03b_PFC_enhancers_18k_PGC__noclump_residual_GWAS_compartment.tsv.gz, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.fam]


    PLINK_clump (
        //CLUMPING of enhancer-based SNP compartments together 
        R_prepare_lists_for_clump.out.lists_before_clump
            .combine(LD_reference)
    )
    // PLINK_clump.out.clumped_SNPs_and_noclump_lists
    //     .view()
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/a1/90c446ff1933f735dc31ec989df901/celso_GWAS_QC.gz, NEURAL_8k_GRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/a1/90c446ff1933f735dc31ec989df901/celso_NEURAL_8k_GRB_significant_EPs_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/a1/90c446ff1933f735dc31ec989df901/celso_NEURAL_8k_GRB_significant_EPs_PGC__noclump_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/a1/90c446ff1933f735dc31ec989df901/celso_PGC_clumped_SNPs.clumped]
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/95/93aeae17e1e7b3b5df6b320b5c8064/celso_GWAS_QC.gz, notNeural_20k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/95/93aeae17e1e7b3b5df6b320b5c8064/celso_notNeural_20k_100flank_noInternalOverlap_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/95/93aeae17e1e7b3b5df6b320b5c8064/celso_notNeural_20k_100flank_noInternalOverlap_PGC__noclump_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/95/93aeae17e1e7b3b5df6b320b5c8064/celso_PGC_clumped_SNPs.clumped]
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/54/65433ffc2328c2dec9fe2a5be9431c/celso_GWAS_QC.gz, NEURAL_14k_noGRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/54/65433ffc2328c2dec9fe2a5be9431c/celso_NEURAL_14k_noGRB_significant_EPs_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/54/65433ffc2328c2dec9fe2a5be9431c/celso_NEURAL_14k_noGRB_significant_EPs_PGC__noclump_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/54/65433ffc2328c2dec9fe2a5be9431c/celso_PGC_clumped_SNPs.clumped]
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/80/8446a04f70e95bc51c4e693a56856b/celso_GWAS_QC.gz, sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/80/8446a04f70e95bc51c4e693a56856b/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/80/8446a04f70e95bc51c4e693a56856b/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_PGC__noclump_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/80/8446a04f70e95bc51c4e693a56856b/celso_PGC_clumped_SNPs.clumped]

    

    R_split_lists (
        // first annotate SNPs with ES of relevant E-P - for ENH SNP list
        // ##################################################### GENERATE MODIFIED ORS MULT BY ES OR EXP       ###########################################################
        // output separate lists to calculate split PRSs and also merged one
        PLINK_clump.out.clumped_SNPs_and_noclump_lists.map { [it, "1"].flatten() }, //######################## multiplier can be set here ########################
        Channel.fromPath( "./input/ES_multipliers/2023-01-18_NEURAL_ENH_EXP_significant_ES_significant_contact_EPs_gene_brain_exp_plus_100_noOverlap.csv.gz", checkIfExists: true)
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
    // [xclo3, 18k_PsychENCODE_PFCortex, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9b/d351570b0db72f3534dfedfe5020bc/xclo3_18k_PsychENCODE_PFCortex_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9b/d351570b0db72f3534dfedfe5020bc/xclo3_18k_PsychENCODE_PFCortex_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9b/d351570b0db72f3534dfedfe5020bc/xclo3_18k_PsychENCODE_PFCortex_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/5e/53ff2c688bdf17a3facfe140ae21ee/xclo3_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/5e/53ff2c688bdf17a3facfe140ae21ee/xclo3_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/5e/53ff2c688bdf17a3facfe140ae21ee/xclo3_QC.fam, /home/osimoe/PGC_w3_data/xclo3, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/c3/22b22bc4e55f5d3b80c64a6e7d9827/xclo3_clumped_LOO_GWAS.tsv.gz, 0.05]
    
    PRSice_calculate_PRS_split_partitions(
        combined_splitlists_bedfile_QCeddata_LDdata
    )
    
    // ########################################### CHANGE NAMES OF MULTIPLIERS ###########################################
    PRS_results = 
        PRSice_calculate_PRS_split_partitions.out.clumped_TS_ENH_GWAS_compartment_PRS
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_residual_GWAS_compartment_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_merged_GWAS_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_original_LOO_GWAS_PRS)
            .map { [it, "enh_ES", "enh_TS_tpm"].flatten() }


    // PRS_results.view()
    // [celso_NEURAL_8k_GRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_residual_GWAS_compartment.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_residual_GWAS_compartment.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_residual_GWAS_compartment.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_merged_GWAS.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_merged_GWAS.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_mult_1_clumped_merged_GWAS.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_QC.fam, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_LOO_GWAS.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_LOO_GWAS.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/07/db66450ec5daead00eda4decd183ed/celso_NEURAL_8k_GRB_significant_EPs_0.5_clumped_LOO_GWAS.best, 0.5, enh_ES, enh_TS_tpm]
    // [celso_34k_neg, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_residual_GWAS_compartment.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_residual_GWAS_compartment.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_residual_GWAS_compartment.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_merged_GWAS.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_merged_GWAS.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_mult_1_clumped_merged_GWAS.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_QC.fam, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_LOO_GWAS.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_LOO_GWAS.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/6f/896a4443021907c26ef3c475040b46/celso_34k_neg_0.5_clumped_LOO_GWAS.best, 0.5, enh_ES, enh_TS_tpm]
    
    R_final_plot (
        PRS_results
    )


}


// 
// srun --pty -t 1-00:00 -p shared -c1 nextflow run https://github.com/emosyne/lisa_percohort_devel -latest -r master -profile lisa -resume 
// --slurmd-debug=error