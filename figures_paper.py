#!/usr/bin/env python

#%% ENV

# imports
import os, sys, re, time
import pandas as pd
import copy as cp
from ete3 import Tree, TreeStyle
import importlib
os.environ["QT_QPA_PLATFORM"] = "offscreen" # this for ete3 to rn properly


# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    threads = 4

else:
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    threads = 48

# define dirs
CurDir = "%s/Cparapsilosis_popgen"%ParentDir
DataDir = "%s/data"%CurDir
ProcessedDataDir = "%s/processed_data"%CurDir
PlotsDir = "%s/plots_paper"%CurDir
ManualDataDir = "%s/manual_data"%CurDir
GenomesDir = "%s/genomes_and_annotations"%DataDir
genome_Cparapsilosis = "%s/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%GenomesDir
gff_Cparapsilosis = "%s/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%GenomesDir
sys.path.append(CurDir); import Cparapsilosis_popgen_functions as fun
importlib.reload(fun)

# make dirs
fun.make_folder(PlotsDir)

#%% LOAD PROCESSED DATASETS

# define the test fields of ML
fields_optimize_test = ["aligner", "CVset_type_split", "consistency_measure", "AUC_diff_threshold", "sample_set", "gene_set", "multiple_test_corr"] # test_sample_set, test_split_ID
threshold_ROC_AUC = 0.7
threshold_consistency_genes = 0.4

# gete the fileterd ML dfs
df_ML_integrated_performance_all = fun.load_object("%s/df_ML_results_all_min_nmodels=4.py"%ProcessedDataDir)
ML_df_dicts = fun.get_processed_ML_classifier_dfs(df_ML_integrated_performance_all, threshold_ROC_AUC, threshold_consistency_genes, fields_optimize_test, ProcessedDataDir, DataDir, CurDir)

# load the GWAS hits df
df_GWAS_hits_all = fun.load_object(ProcessedDataDir + "/df_GWAS_filt_all.py")
df_GWAS_hits_NR = fun.load_object(ProcessedDataDir + "/all_GWAS_hits_NR.py")

# get trees df
min_support_trees = 90 # min branch len for trees
trees_df = fun.load_object("%s/trees_data_df_all.py"%ProcessedDataDir)
if not all([all((n.support>=min_support_trees) or (t.is_root()) for n in t.traverse()) for t in trees_df.tree_corrected]):
    raise ValueError("support not matching")

# load metadata
df_metadata = fun.load_object("%s/df_metadata_all_samples.py"%ProcessedDataDir)
df_metadata_binary_phenotypes = fun.load_object("%s/df_metadata_binary_phenotypes_curated.py"%ProcessedDataDir)

# load GWAS data files
df_raw_gwas_files_all = fun.load_object("%s/df_raw_gwas_files_all.py"%ProcessedDataDir)
#df_transitions_info_GWAS = fun.get_df_transitions_info_GWAS("%s/gwas_transition_data"%ProcessedDataDir, df_raw_gwas_files_all)
df_transitions_info_GWAS = fun.get_df_transitions_info_GWAS_only_rep_samples(DataDir, df_metadata_binary_phenotypes)

# add the trees_df with clade ID data
trees_df = trees_df.apply(fun.get_trees_df_row_with_Bergin2022_clade_info, df_metadata =df_metadata[["runID", "clade_Bergin2022"]].copy(), axis=1)

# define strains to discrd
samples_discard_run1 = {"OSCAR01_L21-060", "OSCAR01_L21-067"}

# tree with clades
tree_series = trees_df[(trees_df.run_set=="exclude_triploids") & (trees_df.mode_tree=="diploid") & (trees_df.type_outgroup=="Cparapsilosis") & (trees_df.aligner=="bwa_mem")].iloc[0].copy()
triploid_strains = set(df_metadata[df_metadata.possibly_triploid==True].runID)
print("There are %i/%i strains with no clade"%(sum([v==-1 for k,v in tree_series.leaf_to_clade_Bergin2022.items() if not k in triploid_strains]), len(set(tree_series.leaf_to_clade_Bergin2022).difference(triploid_strains))))
strains_tree = fun.get_tree_ordered_clades(cp.deepcopy(tree_series.tree_corrected), pd.Series(tree_series.leaf_to_clade_Bergin2022))
strains_tree_uncorrected = fun.get_tree_ordered_clades(cp.deepcopy(tree_series.tree), pd.Series(tree_series.leaf_to_clade_Bergin2022))
                   
df_metadata_trees = df_metadata.copy()
leaf_to_clade_Bergin2022_adds = fun.get_leaf_to_clade_Bergin2022_added_clades(tree_series.leaf_to_clade_Bergin2022, strains_tree, df_metadata, samples_discard_run1)
df_metadata_trees["cladeID_based_on_Bergin2022"] = df_metadata_trees.runID.map(leaf_to_clade_Bergin2022_adds)

# load the GWAS parameter benchmarking df and filter
df_GWAS_parameter_benchmark_filt = fun.get_df_GWAS_parameter_benchmark_filtered(ProcessedDataDir)
all_gwas_phenotypes = sorted(set(df_GWAS_parameter_benchmark_filt.phenotype))

# get StrainInfo table
df_reads_info_analysis_with_validation = fun.load_object(ProcessedDataDir+"/df_reads_info_analysis_with_validation.py").reset_index(drop=True)
df_phenotypes_validation = fun.load_object(ProcessedDataDir+"/df_phenotypes_validation.py")
df_StrainInfo = fun.generate_StrainInfo_table(df_metadata_trees, df_metadata_binary_phenotypes, PlotsDir, df_phenotypes_validation, df_reads_info_analysis_with_validation)

#%% FIGs GWAS

# plot the groups in trees
fun.plot_GWAS_hits_in_trees(df_GWAS_hits_NR, df_GWAS_hits_all, PlotsDir, df_metadata_binary_phenotypes, df_metadata_trees, DataDir)

# print # nodes in FLC resistance tree
tree_flc = Tree(DataDir + "/genotype_phenotype_relationships/pheno-FLC_resistance_aligner-bwa_mem/convergence_gwas/gwas_inputs_representative_samples/corrected_rooted_tree.nw")
print("nnodes tree flc rep samples:", len({n for n in tree_flc.traverse()}))

# report on expected genes -> I decide to benchmark all azoles and CAS FKS2
fun.report_expected_genes_all_phenos(df_GWAS_parameter_benchmark_filt)

# get low confidence parms
df_lowConf_parms_GWAS = fun.get_low_confidence_parameters_GWAS(df_GWAS_parameter_benchmark_filt, PlotsDir)

# get parameter effects
fun.explore_GWAS_parameters_plots(df_GWAS_parameter_benchmark_filt, PlotsDir, fun.manual_curation_gwas_filt_parms)

# print # hits
pheno_to_nhits = df_GWAS_hits_NR[(df_GWAS_hits_NR.type_filtering=="high_confidence")][["phenotype", "unique_groupID"]].drop_duplicates().groupby("phenotype").apply(len).sort_values()
print(", ".join(["'%s' (%i hits)"%(k, v) for k,v in dict(pheno_to_nhits).items()]))
print("# GWAS phenos", len(set(df_GWAS_hits_NR[(df_GWAS_hits_NR.type_filtering=="high_confidence")].phenotype)))

# define sorted phenos
sorted_phenotypes_gwas = fun.get_sorted_phenos_by_hits_and_type(pheno_to_nhits, all_gwas_phenotypes)

# heatmap hits
fun.get_heatmap_GWAS_hits(df_GWAS_hits_NR[df_GWAS_hits_NR.type_filtering=="high_confidence"].copy(), sorted_phenotypes_gwas, PlotsDir)

# get barplots defining transitions
fun.get_stacked_barplots_types_of_pheno_transitions_by_GWAS_hits_and_expected_genes(df_GWAS_hits_NR, sorted_phenotypes_gwas, PlotsDir, DataDir, ProcessedDataDir, df_GWAS_parameter_benchmark_filt, set(all_gwas_phenotypes), short_pheno=True, type_transitions="all", show_n_genes=True)

for type_transitions in ["all", "none", "only_expected"]:
    for show_n_genes in [True, False]:    
        if type_transitions=="none" and show_n_genes is True: continue
        fun.get_stacked_barplots_types_of_pheno_transitions_by_GWAS_hits_and_expected_genes(df_GWAS_hits_NR, sorted_phenotypes_gwas, PlotsDir, DataDir, ProcessedDataDir, df_GWAS_parameter_benchmark_filt, set(all_gwas_phenotypes), short_pheno=False, type_transitions=type_transitions, show_n_genes=show_n_genes)

# print hits
fun.print_GWAS_hits_gene_info(df_GWAS_hits_NR)

#%% FIG ML

# define parms
target_phenos_ML = {'age_infant', 'VRC_resistance', 'source_blood_vs_other_sterile', 'FLC_resistance'}
ML_pheno_to_color = fun.get_ML_pheno_to_color()

# visualize top models genes variants
fun.visualize_top_ML_models_in_tree_all_genes_together(ML_df_dicts, target_phenos_ML, PlotsDir, DataDir, df_metadata_binary_phenotypes, df_metadata_trees, "representative_samples")

# plot ROC CURVES
fun.plot_ROC_curves_best_parameters(ML_df_dicts, PlotsDir, threshold_ROC_AUC, ML_pheno_to_color, add_threshold_lines=True)

# plots performance of top models
for ss in ["representative_samples"]: # all_samples
    fun.plot_performance_distribution_top_models(ML_df_dicts, PlotsDir, ML_pheno_to_color, ss)
    
# print ML things
fun.print_ML_results(ML_df_dicts["df_ML_stats"])
fields_print = ["best_model_AUC_test", "all_samples_AUC_test", "redundant_samples_AUC_test"] + ["best_model_"+f for f in fields_optimize_test]
df_ML_stats_view = ML_df_dicts["df_ML_stats"][ML_df_dicts["df_ML_stats"].n_correct_test_parameters>0][fields_print].copy()

# plot, for each pheno test AUC and model consistency
fun.plot_general_corr_train_test_split_by_pheno(ML_df_dicts["df_ML_integrated_performance_filt"], PlotsDir, threshold_ROC_AUC, df_transitions_info_GWAS, ML_df_dicts["df_ML_stats"], threshold_consistency_genes, "representative_samples", add_colors=True)

# plot affected genes heatmap
fun.ML_affected_genes_heatmap_top_models(ML_df_dicts["df_ML_results_top"], PlotsDir, DataDir)

# visualize trees of the top model
filters_dict_best = {"phenotype":target_phenos_ML, "best_model":True}
df_ML_results_best = fun.get_filtered_df_by_dict(ML_df_dicts["df_ML_results_top"], filters_dict_best)
df_ML_metrics_best = fun.get_filtered_df_by_dict(ML_df_dicts["df_ML_metrics_top"], filters_dict_best)

fun.visualize_top_ML_models_in_trees(df_ML_results_best, df_ML_metrics_best, PlotsDir, DataDir, "representative_samples", df_metadata_binary_phenotypes, df_metadata_trees)

# print # samples
#print(df_ML_metrics_best[df_ML_metrics_best.test_sample_size=="representative_samples"])

# test parameter performance
fun.plot_test_parameter_performance_top_models(ML_df_dicts, PlotsDir, "representative_samples", threshold_ROC_AUC, threshold_consistency_genes)

# train parameter performance
fun.plot_train_parameter_performance_top_models(ML_df_dicts, PlotsDir, "representative_samples", ML_pheno_to_color, DataDir)


#%% FIG ML_GWAS_overlap

# venn diagrams
fun.plot_ML_GWAS_overlap(df_GWAS_hits_NR, ML_df_dicts["df_ML_results_top"], PlotsDir, DataDir)

# heatmap
importlib.reload(fun)
azoles = ['fluconazole R', 'voriconazole R', 'posaconazole R']
all_drugs = ['fluconazole R', 'voriconazole R', 'posaconazole R', 'amphotericin B R', '5-FC R', 'chlorhexidine tol']
ages = ['age infant', 'age old']
sources_phenos = ['skin exudate', 'bronchial aspirate', 'ascitic fluid', 'blood culture', 'blood vs other sterile', 'blood', 'sterile']

for target_phenotypes in ["all", azoles, all_drugs, ages, sources_phenos]:
#for target_phenotypes in ["all"]:

    if target_phenotypes=="all": add_description_list = [True, False]
    else: add_description_list = [True]

    for add_description in add_description_list:
        if target_phenotypes=="all" and add_description is True: generate_proc_table = True
        else:  generate_proc_table = False
        fun.plot_ML_GWAS_overlap_heatmap(df_GWAS_hits_NR, ML_df_dicts["df_ML_results_top"], PlotsDir, DataDir, ProcessedDataDir, add_description=add_description, target_phenotypes=target_phenotypes, generate_proc_table=generate_proc_table)
        
#%% PRINT NUMBERS

# strain numbers
fun.print_strain_numbers(ProcessedDataDir, ManualDataDir, df_metadata, df_metadata_binary_phenotypes)

#%% FIG KnownPhenos

# plot the trees with expected genes vars grouped by drug type
fun.plot_trees_expected_genes_variants_grouped_by_type_drug(DataDir, PlotsDir, df_metadata_trees, samples_discard_run1, strains_tree_uncorrected, df_metadata_binary_phenotypes)

# plot the trees with expected genes vars
#for t_ss in ["representative_samples", "all_samples"]:    
#    fun.plot_trees_expected_genes_variants(DataDir, PlotsDir, df_metadata_trees, samples_discard_run1, strains_tree, df_metadata_binary_phenotypes, t_ss)

#%% FIG Dataset

# plot tree
fun.plot_strains_tree_with_clades("%s/Fig_Dataset_tree_clades.pdf"%PlotsDir, strains_tree, df_metadata_trees, samples_discard_run1)

# countries distribution
fun.plot_clade_distribution_per_countries(df_metadata_trees, PlotsDir, samples_discard_run1)

# seq here distribution
fun.plot_clade_distribution_per_hospitals(df_metadata_trees, PlotsDir, samples_discard_run1)

#%% FIG Phenotypes

# tree wuth phenotypes presence / absence
fun.plot_Phenotypes_presence_tree(df_metadata_trees, PlotsDir, samples_discard_run1, strains_tree, df_metadata_binary_phenotypes, df_transitions_info_GWAS)

#%% ED FIG PhenoDefinition

# get df metadata
df_metadata_filt = df_metadata[(df_metadata.possibly_triploid==False) & ~(df_metadata.runID.isin(samples_discard_run1)) & (df_metadata.source_data=="sequenced_here")].copy()

# continuous phenotypes distributions
fun.Fig_PhenoDefinition_continuous_phenotypes(df_metadata_filt, df_metadata_binary_phenotypes, PlotsDir)


# fig colony formation
fun.Fig_PhenoDefinition_colony_formation(df_metadata_filt, df_metadata_binary_phenotypes, PlotsDir)

# higher in solid
fun.fig_PhenoDefinition_higher_R_solid(df_metadata_filt, df_metadata_binary_phenotypes, PlotsDir)
                                       
# flz tolerance
fun.fig_PhenoDefinition_agar_tolerance_48h(df_metadata_filt, df_metadata_binary_phenotypes, PlotsDir)
    
# age distribution
fun.fig_PhenoDefinition_age(df_metadata_filt, df_metadata_binary_phenotypes, PlotsDir)

#%% ED FIG Recombination

# keep final data
df_metadata_filt = df_metadata[(df_metadata.possibly_triploid==False) & ~(df_metadata.runID.isin(samples_discard_run1))].copy()

# plot recombination
for zygosity in ["homo", "hetero"]:
    for aligner in ["bwa_mem"]:
        for window_size in [100000, 20000, 1000]:
            for t in {"homo":[0.05], "hetero":[0.001]}[zygosity]:
                if window_size==1000 and zygosity=="hetero": continue
                
                yfield = "fraction_diff_vars_per_kb_below_%s"%t
                fun.plot_general_stats_SNP_distance_comparisons(ProcessedDataDir, aligner, window_size, df_metadata_filt, PlotsDir, zygosity, yfield)

#%% ED FIG Tree parameters

# keep trees to plot
df_tree_stats = trees_df[(trees_df.type_outgroup=="Cparapsilosis") & (trees_df.run_set=="exclude_triploids")].sort_values(by=["type_outgroup", "run_set", "mode_tree", "aligner"]).copy()

for f in ["nsites_parsimony", "pct_collapsed_nodes", "pct_collapsed_long_nodes", "pct_politiomies_diverse"]: 
    print(df_tree_stats.sort_values(by=["mode_tree", f]).set_index("mode_tree")[f], "\n")

# plot the relevant stats
aligner_to_color = {"bwa_mem":"c", "bowtie2_local":"violet", "hisat2_no_spliced":"red"}
fun.plot_treeQC_stats(df_tree_stats, PlotsDir, aligner_to_color)

# plot tree comparisons
fun.heatmap_differences_trees(df_tree_stats, PlotsDir, aligner_to_color)

#%% ED Fig Validation Azoles
fun.plot_validation_azole_phenotypes_tree(df_StrainInfo, DataDir, PlotsDir)

                
#%% ANALYSES FOR G458S PAPER

# tree with G458S data paper
fun.plot_tree_clinical_paper(df_metadata_trees, PlotsDir, samples_discard_run1, strains_tree, df_metadata_binary_phenotypes)

# table withe used metadata, and tree
fun.generate_data_clinical_paper(df_metadata_trees, PlotsDir, samples_discard_run1, strains_tree, df_metadata_binary_phenotypes)

    
