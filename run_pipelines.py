#!/usr/bin/env python


"""
This script should be ran on Cparapsilosis_popgen_env. The perSVade parts should be run on perSVade_env (tested on latest commit by Nov 2024). It is to generate all the reads and variant calls. The reads were sequenced WITHOUT pools, in projects OSCAR01 and OSCAR02 and OSCAR03. 

Ref genome from NCBI:

    genome: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/765/GCA_000182765.2_ASM18276v2/GCA_000182765.2_ASM18276v2_genomic.fna.gz
    gff: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/765/GCA_000182765.2_ASM18276v2/GCA_000182765.2_ASM18276v2_genomic.gff.gz

Ref genome from CGD:

    genome: http://www.candidagenome.org/download/sequence/C_parapsilosis_CDC317/current/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta.gz
    gff: http://www.candidagenome.org/download/gff/C_parapsilosis_CDC317/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff
    chromosomal features: http://www.candidagenome.org/download/chromosomal_feature_files/C_parapsilosis_CDC317/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab 

Other genomes for depooling:

    - C. auris genome from Cmine: the gDNA from GCA_002759435.2 and the mtDNA from ASM716870v1. 

    - C. orthopsilosis genome from Cmine: gDNA (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/315/875/GCF_000315875.1_ASM31587v1/GCF_000315875.1_ASM31587v1_genomic.fna.gz) & mtDNA (https://www.ncbi.nlm.nih.gov/nuccore/NC_006972)

    - C. metapsilosis genome from Cmine: Assembly by Veronica (~/samba_bsc40/vdepinho/Candida_metapsilosis/reference_nanopore/CANME.fa)

    - C. albicans: haplotype A from http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_version_A22-s07-m01-r110_chromosomes.fasta.gz


I take the genome from CGD, as it is the same sequence, and also it's annotations.

Annotations:

    - I take the OBO file from https://purl.obolibrary.org/obo/go/go-basic.obo at 27/6/2024.
    - I take the Reactome pathways from https://reactome.org/download/current/ReactomePathways.txt at 28/6/2024
    - I take the Reactome pathway relations from https://reactome.org/download/current/ReactomePathwaysRelation.txt at 28/6/2024
    - GO terms from  http://www.candidagenome.org/download/go/archive/gene_association.20231003.cgd.gz (same genome version). Version C. para s01-m06-r03

I obtained the SRA data by going to the SRA run selector and getting the data for C. parapsilosis (taxID 5480). Filtering by source=DNA, type=genome, layout=paired, strategy=genome, file type=fastq. Then I downloaded the results of the RunSelector. I did this on 06/11/2023 and 05/12/2024. 

"""

##### GET ENV #######
print("defining environment...")

# imports
import os, sys, re, time
import pandas as pd
import copy as cp
import multiprocessing as multiproc
import shutil
import numpy as np
#import seaborn as sns
#import matplotlib.pyplot as plt

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    bsc40_dir =  "%s/samba_bsc40"%(os.getenv("HOME"))
    threads = 4

else:
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    bsc40_dir = "/gpfs/projects/bsc40" 
    threads = 48

# define dir
CurDir = "%s/Cparapsilosis_popgen"%ParentDir
sys.path.append(CurDir); import Cparapsilosis_popgen_functions as fun

# define extra folders
ManualDataDir = "%s/manual_data"%CurDir
DataDir = "%s/data"%CurDir
ProcessedDataDir = "%s/processed_data"%CurDir; fun.make_folder(ProcessedDataDir)
VarCallDir = "%s/varcall"%DataDir; fun.make_folder(VarCallDir)
ValidationDir = "%s/validation_data"%DataDir; fun.make_folder(ValidationDir)
ReadsDir = "%s/reads"%DataDir; fun.make_folder(ReadsDir)
PlotsDir = "%s/plots"%CurDir; fun.make_folder(PlotsDir)
GenomesDir = "%s/genomes_and_annotations"%DataDir
genome_Cparapsilosis = "%s/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%GenomesDir # genome from CGD
gff_Cparapsilosis = "%s/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%GenomesDir

# define the aligners to consider
sorted_aligners = ["bwa_mem", "bowtie2_local", "hisat2_no_spliced"]

# define some parameters
min_coverage = 50
min_pct_reads_para = 95
min_pct_covered = 90
min_coverage_pos = 12

# define whether to submit jobs
submit_jobs = True

####################

####### GET READS ####

# init the all_cmds
list_cmds = [] # It has "name_job", "cmd", "env", "threads"

# compare the genomes of CGD and NCBI (they are the same)
fun.check_that_CGD_and_NCBI_genomes_are_the_same(DataDir)

# init df paths with the sequenced genomes
print("getting sequenced genomes")
df_paths = pd.read_csv("%s/UPDATED_seq_files_03012024.tsv"%ManualDataDir, sep="\t")
df_paths = df_paths[df_paths["PROJECT NUMBER*"].isin({"OSCAR01", "OSCAR02", "OSCAR03"})]
df_paths = df_paths[["PROJECT NUMBER*", "BSC file*"]].rename(columns={"PROJECT NUMBER*":"project_number"})
df_paths["sampleID"] = df_paths["BSC file*"].apply(lambda f: "-".join(fun.get_file(f).split("_")[0].upper().split("-")[0:2]))
df_paths["readID"] = df_paths["BSC file*"].apply(lambda f: fun.get_file(f).split("_")[-2])
df_paths["raw_reads"] = bsc40_dir + "/" + df_paths["BSC file*"]
df_paths = df_paths[["project_number", "sampleID", "readID", "raw_reads"]]
if len(df_paths)!=len(df_paths[["project_number", "sampleID", "readID"]].drop_duplicates()): raise ValueError("non-unique samples")
df_paths["runID"] = df_paths.project_number + "_" + df_paths.sampleID 

# save
fun.save_df_as_tab(df_paths, "%s/sequencing_runs_OSCAR.tab"%DataDir)

# get the metadata df, with consistent sample namings
df_metadata_Oscar =  fun.get_df_metadata_final(ManualDataDir, set(df_paths.sampleID))
if sorted(df_metadata_Oscar.sampleID)!=sorted(set(df_paths.sampleID)): raise ValueError("not valid dfs")

# add tuple samples
for df in [df_paths, df_metadata_Oscar]:
    df["sampleID_tuple"] = df.sampleID.apply(fun.get_numeric_tuple_sampleID)

# check consistency in sampleIDs
if not all(df_paths[["sampleID", "sampleID_tuple"]].drop_duplicates().groupby("sampleID_tuple").apply(len)==1): raise ValueError("There should be one sID for each tuple")
if not all(df_metadata_Oscar[["sampleID", "sampleID_tuple"]].drop_duplicates().groupby("sampleID_tuple").apply(len)==1): raise ValueError("There should be one sID for each tuple")
if len(set(df_paths.sampleID))!=len(set(df_paths.sampleID_tuple)): raise ValueError("There should be one sampleID correctly formated for each sample")
if len(set(df_metadata_Oscar.sampleID))!=len(set(df_metadata_Oscar.sampleID_tuple)): raise ValueError("There should be one sampleID correctly formated for each sample")
if len(df_paths)!=(len(set(df_paths.runID))*2): raise ValueError("there should be 2 reads for each sampleID")

# print samples sequenced more than once
sID_to_nruns = df_paths[df_paths.readID=="R1"].groupby("sampleID").apply(len)
print("samples >1 runs", sID_to_nruns[sID_to_nruns>1], "# samples", len(set(df_paths.sampleID)))

# keep
df_paths = df_paths[["runID", "readID", "raw_reads"]].copy()
df_paths["source_data_type"] = "sequenced_here"

# get the SRA runs from the SRA
print("adding SRA data")
df_sra_runs = pd.read_csv("%s/SraRunTable_Cpara_WGS_06112023.csv"%ManualDataDir, sep=",")
df_sra_runs = df_sra_runs[~df_sra_runs.BioProject.isin({"PRJEB57083", "PRJEB45347", "PRJEB45642", "PRJEB60076", "PRJEB56880"})] # remove IV-evolved samples with no paper
df_sra_runs = df_sra_runs[df_sra_runs.BioProject!="PRJNA866533"] # discard a bioproject that were Crispr-Cas9 modified strains
df_sra_runs = df_sra_runs[df_sra_runs.LibrarySource=="GENOMIC"] # discard some metagenomes
df_sra_runs = df_sra_runs[df_sra_runs.apply(fun.sra_run_is_clinical_or_enviromental_WGS, axis=1)]
df_sra_runs = df_sra_runs[~df_sra_runs.Platform.isin({"BGISEQ", "DNBSEQ"})] # discard minoritary platforms, which do not give paired-end SRRs
df_sra_runs = df_sra_runs[~df_sra_runs.Run.isin(fun.wrong_srrs)].copy() # discard bad srrs
initial_sra_srrs = set(df_sra_runs.Run)
print("Total of %i SRA runs of para. Types of data"%(len(df_sra_runs)))
print(fun.Counter(df_sra_runs.apply(fun.sra_run_get_type_of_isolate, axis=1)))

# load the extra runs for validation, checking that matching bioprojects have the same data
print("loading extra SRA data")
df_sra_runs_Dec2024 = pd.read_csv("%s/SraRunTable_Cpara_WGS_05122024.csv"%ManualDataDir, sep=",")
df_sra_runs_validation = fun.get_df_sra_runs_validation(df_sra_runs_Dec2024, df_sra_runs)
df_phenotypes_validation = fun.get_df_phenotypes_validation(df_sra_runs_validation, ManualDataDir, ProcessedDataDir)
fun.save_object(df_phenotypes_validation, ProcessedDataDir+"/df_phenotypes_validation.py")
df_sra_runs_validation = df_sra_runs_validation[df_sra_runs_validation.Run.isin(set(df_phenotypes_validation.Run))].copy()
print("Total of %i SRA runs for validation"%(len(df_sra_runs_validation)))

# add the extra runs to df_sra_runs
df_sra_runs_validation_extra = df_sra_runs_validation[~df_sra_runs_validation.Run.isin(set(df_sra_runs.Run))].copy()
extra_srrs_validation = set(df_sra_runs_validation_extra.Run) # define extra SRRs
df_sra_runs = pd.concat([df_sra_runs, df_sra_runs_validation_extra])

# add homozygous C. orthopsilosis strains
homozygous_Corthopsilosis_srrs = fun.get_homozygous_Corthopsilosis_srrs(ManualDataDir)
df_sra_runs = pd.concat([df_sra_runs, pd.DataFrame({"Run":sorted(homozygous_Corthopsilosis_srrs)})])
print("Total of %i SRA runs of ortho."%(len(homozygous_Corthopsilosis_srrs)))

# check
if len(df_sra_runs)!=len(set(df_sra_runs.Run)):
    raise ValueError("should be unique Run")

# define a boolean that indicates that read info has been geenrated. This is useful to be able to remove reads
read_info_missing = fun.file_is_empty("%s/df_reads_info.py"%ProcessedDataDir)

# download SRR files (precompiled files with reads)
print("downloading SRA srr files")
sorted_srrs = sorted(set(df_sra_runs.Run))
srr_to_raw_reads = {}
inputs_srr_downloads = []
for I, srr in enumerate(sorted_srrs):
    folder_srr = "%s/%s"%(ReadsDir, srr); fun.make_folder(folder_srr)
    folder_reads = "%s/reads"%folder_srr; fun.make_folder(folder_reads)
    srr_file = "%s/%s.srr"%(folder_reads, srr)
    srr_to_raw_reads[srr] = ("%s_1.fastq.gz"%srr_file, "%s_2.fastq.gz"%srr_file)

    if not read_info_missing: continue
    if not fun.file_is_empty(srr_to_raw_reads[srr][1]): continue # if the raw reads were already generated, do not repeat
    inputs_srr_downloads.append((srr, srr_file, I+1, len(sorted_srrs), False))

if len(inputs_srr_downloads)>0:
    with multiproc.Pool(threads) as pool:
        pool.starmap(fun.download_srr_with_prefetch, inputs_srr_downloads, chunksize=1)
        pool.close()


# fastqdump the files and add to df_paths
print("getting fastqdump cmds")
for srr in sorted_srrs:

    # get reads
    folder_reads = "%s/%s/reads"%(ReadsDir, srr)
    raw_reads1, raw_reads2 = srr_to_raw_reads[srr]

    cmd_tuple = ("fastqdump", "%s --srr %s --outdir %s --threads %i --stop_after_fastqdump"%(fun.get_trimmed_reads_for_srr_py, srr, folder_reads, 56), "perSVade_env", 56)
    if fun.file_is_empty(raw_reads2) and read_info_missing: 
        if submit_jobs is True: list_cmds.append(cmd_tuple); continue
        else: fun.run_cmd_tuple(cmd_tuple)

    # check
    if any([fun.file_is_empty(f) for f in [raw_reads1, raw_reads2]]) and read_info_missing: 
        raise ValueError("raw reads can't be empty")

    # clean
    for f in os.listdir(folder_reads):
        file_path = "%s/%s"%(folder_reads, f)
        if not file_path in {raw_reads1, raw_reads2}: fun.remove_file(file_path)

    # define the source data type
    if fun.get_source_data(srr)=="SRA" and srr in extra_srrs_validation: source_data_type = "SRA_extra_validation"
    elif fun.get_source_data(srr)=="SRA" and srr in homozygous_Corthopsilosis_srrs: source_data_type = "SRA_Corthopsilosis"
    elif fun.get_source_data(srr)=="SRA" and srr in initial_sra_srrs: source_data_type = "SRA_initial"
    else: raise ValueError("invalid")

    # add to df
    df_srr =  pd.DataFrame({"raw_reads":[raw_reads1, raw_reads2], "readID":["R1", "R2"]})
    df_srr["runID"] = srr
    df_srr["source_data_type"] = source_data_type
    df_paths = pd.concat([df_paths, df_srr]).reset_index(drop=True)

# define the type of data
df_paths["source_data"] = df_paths.runID.apply(fun.get_source_data)
df_paths = df_paths.set_index(["runID", "readID"], drop=False)

# get through each runID and get read stats (mapping, coverage). This performs downsampling, trimming and mapping
if len(list_cmds)==0:
    print("adding read_info_data...")
    df_reads_info = fun.get_df_reads_info(df_paths, CurDir, list_cmds, ProcessedDataDir, ReadsDir, genome_Cparapsilosis, GenomesDir, min_coverage, sorted_aligners, VarCallDir, gff_Cparapsilosis, submit_jobs, homozygous_Corthopsilosis_srrs)

# split df_reads_info to have the used / extra validation ones
if len(list_cmds)==0:
    all_runs_before_split = sorted(set(df_reads_info.runID))
    fun.check_no_nans_in_df(df_reads_info[["source_data_type"]])
    runID_to_sourceDataType = dict(df_reads_info.set_index("runID").source_data_type)
    df_reads_info_extra_validation = df_reads_info[df_reads_info.source_data_type=="SRA_extra_validation"].copy()
    df_reads_info = df_reads_info[df_reads_info.source_data_type.isin({"sequenced_here", "SRA_initial", "SRA_Corthopsilosis"})].copy()

# run multiqc, for different projects
if len(list_cmds)==0:

    for suffix in ["SRA", "OSCAR01", "OSCAR02", "OSCAR03", "all"]:
        if suffix=="SRA": sorted_runs = sorted(set(df_reads_info[df_reads_info.source_data=="SRA"].runID))
        elif suffix=="all": sorted_runs = sorted(set(df_reads_info.runID))
        else: sorted_runs = sorted(set(df_reads_info[df_reads_info.runID.apply(lambda x: x.split("_")[0]==suffix)].runID))
        fun.run_multiqc_several_trimmed_reads("%s/multiqc_outdir_%s"%(ProcessedDataDir, suffix), sorted_runs, function_origin_zip=fun.define_origin_zip_Cparapsilosis_popgen)

###################

######## VAR CALLING ########

# infer repeats
if len(list_cmds)==0:
    print("Inferring repeats")
    simple_repeats_file = fun.get_simple_repeats_file(genome_Cparapsilosis, "%s/repeat_inference"%DataDir)

# run small variant calling
if len(list_cmds)==0:
    print("getting small variant calling")
    print("All runs:", len(all_runs_before_split))
    print("len df_reads_info", len(df_reads_info))
    print("All runs SRA_inital + sequenced here (there are two duplicates):", len(df_reads_info[df_reads_info.source_data_type!="SRA_Corthopsilosis"]))
    print("All runs extra:", len(df_reads_info_extra_validation))
    print("All C. ortho:", len(df_reads_info[df_reads_info.source_data_type=="SRA_Corthopsilosis"]))


    for aligner in sorted_aligners:
        for runID in all_runs_before_split:

            # discard some cases
            if runID in homozygous_Corthopsilosis_srrs and aligner!="bwa_mem": continue # for orthopsilosis I could only find the bwa_mem ones

            # get data
            outdir_aligner = "%s/%s/aligner-%s"%(VarCallDir, runID, aligner)
            outdir_small_vars = "%s/small_vars"%outdir_aligner
            sorted_bam = "%s/align_reads/aligned_reads.bam.sorted"%outdir_aligner
            if fun.file_is_empty(sorted_bam): raise ValueError("invalid bam")
            
            cmd_tuple = ("small_vars", "%s call_small_variants -o %s --threads 28 -r %s --min_chromosome_len 5000 -sbam %s --repeats_file skip -p 2 --callers HaplotypeCaller,bcftools,freebayes  -c %i --min_AF 0.25 --window_freebayes_bp 200000 --fractionRAM_to_dedicate 0.8 --verbose --fraction_available_mem 0.5"%(fun.perSVade_wrapper, outdir_small_vars, genome_Cparapsilosis, sorted_bam, min_coverage_pos), "perSVade_env", 28)

            if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir_small_vars): 
                if submit_jobs is True: list_cmds.append(cmd_tuple); continue
                else: fun.run_cmd_tuple(cmd_tuple)

# integrate the coverage dataset into one dataset (this is per gene)
if len(list_cmds)==0:
        
    # get raw coverage
    coverage_df_file = "%s/coverage_per_gene_all_samples.py"%ProcessedDataDir
    full_sampleID_to_cov_outdir = {"%s---%s"%(runID, aligner) : "%s/%s/aligner-%s/cov_per_gene"%(VarCallDir, runID, aligner) for runID in sorted(set(df_reads_info.runID)) for aligner in sorted_aligners if not (aligner=="hisat2_no_spliced" and runID in homozygous_Corthopsilosis_srrs)}

    fields_cov = ['chromosome', 'start', 'end', 'median_reads_per_gene', 'percentcovered_1', 'ID', 'relative_coverage']
    fun.integrate_get_cov_genes_several_samples(full_sampleID_to_cov_outdir, coverage_df_file, multiproc.cpu_count(), fields=fields_cov)

    # get the df with rel coverage
    df_cov_per_gene_file = fun.get_df_cov_per_gene_relative_gDNA(ProcessedDataDir)

# generate cov per windows
window_size_SNP_densities = 1000
window_sizes_recombination_blocks = [200000, 100000, 50000, 20000]

if len(list_cmds)==0:
    print("getting cov per window and stats")
    aligner_to_wsize_to_df_coverage_per_window_file = {}
    for aligner in sorted_aligners:
        for wsize in (window_sizes_recombination_blocks + [window_size_SNP_densities]):

            cov_file = fun.generate_df_coverage_per_window(list_cmds, aligner, wsize, VarCallDir, sorted(set(df_reads_info.runID)), genome_Cparapsilosis, simple_repeats_file, ProcessedDataDir, homozygous_Corthopsilosis_srrs)
            aligner_to_wsize_to_df_coverage_per_window_file.setdefault(aligner, {}).setdefault(wsize, cov_file)

# integrate varcalls
if len(list_cmds)==0:

    # get the integrated variants
    print("integrate and filter small vars")    
    aligner_to_df_small_vars_filt_file = {}
    for aligner in sorted_aligners:
        aligner_to_df_small_vars_filt_file[aligner] = fun.get_filtered_small_vars_df(ProcessedDataDir, VarCallDir, sorted(set(df_reads_info.runID)), aligner, simple_repeats_file, homozygous_Corthopsilosis_srrs, aligner_to_wsize_to_df_coverage_per_window_file[aligner][window_size_SNP_densities], window_size_SNP_densities)

    # get a df with SNP / INDEL densities
    print("calculate SNP densities")
    aligner_SNP_densities = "bwa_mem"
    df_variant_densities_file = fun.generate_variant_per_window_densities_df(aligner_to_df_small_vars_filt_file[aligner_SNP_densities], ProcessedDataDir, genome_Cparapsilosis, aligner_SNP_densities, window_size_SNP_densities, sorted(set(df_reads_info.runID)), simple_repeats_file, aligner_to_wsize_to_df_coverage_per_window_file[aligner_SNP_densities][window_size_SNP_densities])

    # get stats about pairwise distances between windows (recombination analysis)
    target_runs_distances = sorted(set(df_reads_info.runID).difference(homozygous_Corthopsilosis_srrs))
    df_stats_pairwise_SNP_distances = fun.get_stats_pairwise_SNP_distances_per_window(aligner_to_df_small_vars_filt_file, aligner_to_wsize_to_df_coverage_per_window_file, window_sizes_recombination_blocks + [window_size_SNP_densities], ProcessedDataDir, CurDir, target_runs_distances, genome_Cparapsilosis, simple_repeats_file, list_cmds)

#############################

########## FILTER DATASETS #########

# reports about dataset quality
print("getting reports of reads...")
if len(list_cmds)==0:

    # create df with added info about coverage and multiqc
    df_reads_info = fun.get_df_reads_info_with_multiqc_and_coverage_data(df_reads_info, ProcessedDataDir, sorted_aligners, homozygous_Corthopsilosis_srrs, df_variant_densities_file, genome_Cparapsilosis, PlotsDir, aligner_SNP_densities, window_size_SNP_densities)
    print("Initial df_reads_info has %i strains"%len(df_reads_info))

    # generate all the QC reports, getting df to be run in parallel
    df_reads_info = fun.run_QC_and_filtering_samples(df_reads_info, homozygous_Corthopsilosis_srrs, ProcessedDataDir, ReadsDir, min_coverage, min_pct_reads_para, min_pct_covered, PlotsDir, aligner_to_df_small_vars_filt_file, df_metadata_Oscar, replace=False, skip_plots=False)

    print("Final df_reads_info strain type %s "%fun.Counter(df_reads_info.source_data))
    print("Final df_reads_info has %i strains"%len(df_reads_info))
    print("# possibly triploid strains: %i"%(sum(df_reads_info.possibly_triploid)))
    print("# ortho strains: %i"%(sum(df_reads_info.species_run=="ortho")))

    # define, and save, the df to do analyses
    samples_discard_run1 = {"OSCAR01_L21-060", "OSCAR01_L21-067"}
    df_reads_info_analysis = df_reads_info[(df_reads_info.possibly_triploid==False) & (df_reads_info.species_run=="para") & ~(df_reads_info.runID.isin(samples_discard_run1))].copy()
    fun.save_object(df_reads_info_analysis, "%s/df_reads_info_filtered_afterQC_samples_analyze.py"%ProcessedDataDir)
    print("# strains analysis: %i"%(len(df_reads_info_analysis)))
    print("# strains smiley: %i"%(sum(df_reads_info_analysis.has_smiley)))

####################################

#### GET OTHER VAR CALLS ####

# SV & CNV calls, integration and annotation
if len(list_cmds)==0:
    
    # get the df of samples to run SV calling on
    df_reads_info_SVs = df_reads_info_analysis[(df_reads_info_analysis.has_smiley==False)].copy()

    # get cmds
    fun.get_cmds_perSVade_SVs_and_CNVs_calls(VarCallDir, sorted_aligners, simple_repeats_file, sorted(set(df_reads_info_SVs.runID)), genome_Cparapsilosis, ProcessedDataDir, list_cmds, DataDir, gff_Cparapsilosis)

# get the small variants integrated and annotated
if len(list_cmds)==0:
    fun.integrate_and_annotate_small_vars(VarCallDir, sorted_aligners,  sorted(set(df_reads_info_analysis.runID)), genome_Cparapsilosis, ProcessedDataDir, list_cmds, DataDir, gff_Cparapsilosis, simple_repeats_file)

##### RUN SNP TREES #####

if len(list_cmds)==0:
    print("getting tree cmds")
    trees_data_df = fun.generate_SNPs_trees(df_reads_info, ProcessedDataDir, DataDir, list_cmds, sorted_aligners, min_coverage_pos, homozygous_Corthopsilosis_srrs, VarCallDir, genome_Cparapsilosis, PlotsDir)

########################

##### GENO-PHENO RELATIONSHIPS #######

# GWAS
if len(list_cmds)==0:

    # load df metadata dfs
    df_metadata = fun.get_df_metadata_all_samples(CurDir, PlotsDir)
    df_metadata_binary_phenotypes = fun.get_df_metadata_binary_phenotypes_curated(df_metadata, ProcessedDataDir)

    # generate a file with the clade IDs 
    cladeID_info_file = fun.get_cladeID_info_file_Bergin2022(df_metadata, trees_data_df, ProcessedDataDir)

    # run interproscan
    ancestral_GWAS_dir = "%s/scripts/ancestral_GWAS/scripts"%ParentDir
    outdir_IP = "%s/interproscan_annot"%GenomesDir

    if fun.file_is_empty("%s/pipeline_finished.txt"%outdir_IP):

        CondaDir = "%s/anaconda3"%ParentDir
        lib_definition = "export LD_LIBRARY_PATH=%s/envs/InterProScan_env/lib:$LD_LIBRARY_PATH"%(CondaDir)
        interproscan_sh = "%s/software/interproscan_v5.52-86.0/interproscan-5.52-86.0/interproscan.sh"%ParentDir
        prefix_env = '%s && export INTERPROSCAN_SH="source %s/etc/profile.d/conda.sh && conda activate InterProScan_env && %s"'%(lib_definition, CondaDir, interproscan_sh)
        
        fun.run_cmd('%s && %s/run_interproscan_from_gff.py --gff %s --ref_genome %s --outdir %s --gDNA_code %i --mtDNA_code %i --threads 56 --mtDNA_chrom %s'%(prefix_env, ancestral_GWAS_dir, gff_Cparapsilosis, genome_Cparapsilosis, outdir_IP, fun.Cpara_gDNA_code, fun.Cpara_mtDNA_code, fun.mtDNA_chrom), env="ancestral_GWAS_env")

    # prepare pathway files
    outdir_pathways = "%s/pathway_annotations"%GenomesDir
    if fun.file_is_empty("%s/pipeline_finished.txt"%outdir_pathways):

        # You should have run in a terminal with internet: conda activate ancestral_GWAS_env; python -c "from ete3 import NCBITaxa; NCBITaxa()".
        # This needs to be run running '/home/mschikora/software/pathway-tools/pathway-tools -lisp -python'
        fun.run_cmd('%s/prepare_pathway_annotations.py --gff %s --outdir %s --threads 56 --interproscan_output %s/interproscan_annotation.out --reactome_pathways %s/ReactomePathways_28062024.txt --reactome_pathway_relations %s/ReactomePathwaysRelation_28062024.txt --max_fraction_genes_pathway 0.05 --reactome_pathways_keep Candida --MetaCyc_pathways_clade Ascomycota --obo_file %s/go-basic_27062024.obo --GOterms %s/gene_association.20231003.cgd.gz --GOterms_type gene_association_CGD --NCBI_taxa_db %s/NCBI_taxa_28June2024.sqlite'%(ancestral_GWAS_dir, gff_Cparapsilosis, outdir_pathways, outdir_IP, GenomesDir, GenomesDir, GenomesDir, GenomesDir, GenomesDir), env="ancestral_GWAS_env")

    # run different jobs for each phenotype
    all_phenos = sorted(set(df_metadata_binary_phenotypes.keys().difference({"sampleID", "runID"})))
    GenoPheno_dir_all = "%s/genotype_phenotype_relationships"%DataDir; fun.make_folder(GenoPheno_dir_all)
    df_raw_gwas_files_all = pd.DataFrame()
    phenos_sufficient_numbers = set()    
    phenos_sufficient_rep_transitions_some_aln = set()

    for pheno in all_phenos:
        for aligner in sorted_aligners:
            #print("GWAS", pheno, aligner)

            # debugging, discard some phenotypes
            #if pheno.startswith("source_"): continue # debug !!!
            #if pheno.startswith("age_"): continue # debug !!!!
            #if pheno!="POS_resistance": continue # debug!!!

            # dir dir
            GenoPheno_dir = "%s/pheno-%s_aligner-%s"%(GenoPheno_dir_all, pheno.replace(" ", "_"), aligner)

            # get the phenotypes df
            df_phenotypes = df_metadata_binary_phenotypes[df_metadata_binary_phenotypes[pheno].isin({0, 1})][["runID", pheno]].copy()

            # discard phenotypes with <3 for one cathegory
            n_zeros = sum(df_phenotypes[pheno]==0)
            n_ones = sum(df_phenotypes[pheno]==1)
            if n_zeros<3 or n_ones<3:
                continue
            
            # keep for recording
            phenos_sufficient_numbers.add(pheno)

            # get the tree for this subset of strains
            fun.make_folder(GenoPheno_dir)
            tree_dir = "%s/tree_subset_strains"%GenoPheno_dir
            treefile_all_samples = fun.generate_tree_one_phenotype_and_aligner(list_cmds, trees_data_df, tree_dir, aligner, set(df_phenotypes.runID), VarCallDir, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis)

            # debug
            if treefile_all_samples is None:
                continue

            # run convergence GWAS
            gwas_outdir = "%s/convergence_gwas"%GenoPheno_dir
            df_raw_gwas_files = fun.get_cmds_convergence_GWAS(gwas_outdir, list_cmds, df_phenotypes, treefile_all_samples, VarCallDir, pheno, aligner, trees_data_df, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, min_coverage_pos, simple_repeats_file, DataDir, gff_Cparapsilosis, cladeID_info_file, replace=False)

            # keep files
            if not df_raw_gwas_files is None:
                df_raw_gwas_files["aligner"] = aligner
                df_raw_gwas_files["phenotype"] = pheno
                phenos_sufficient_rep_transitions_some_aln.add(pheno)
                df_raw_gwas_files_all = pd.concat([df_raw_gwas_files_all, df_raw_gwas_files]).reset_index(drop=True)

# logs
if len(list_cmds)==0:    
    print("# initial phenotypes", (len(all_phenos)))
    print("# phenotypes with >=3 1s and >=3 0s", len(phenos_sufficient_numbers))
    print("# phenotypes with >=3 transitions in representative samples by some aligner", len(phenos_sufficient_rep_transitions_some_aln))


# ML models
if len(list_cmds)==0:    
    print("Getting ML models...")

    # get models
    min_nmodels = 4
    fun.run_ML_predictors_phenotypes(df_raw_gwas_files_all, GenoPheno_dir_all, list_cmds, DataDir, gff_Cparapsilosis, genome_Cparapsilosis, sorted_aligners, "%s/df_ML_results_all_min_nmodels=%i.py"%(ProcessedDataDir, min_nmodels), ProcessedDataDir, min_nmodels, replace=False)

# bechmark the GWAS results for some conditions to get best parameters
if len(list_cmds)==0:

    # filter out phenotypes and save
    df_raw_gwas_files_all = df_raw_gwas_files_all[df_raw_gwas_files_all.phenotype!="source_SANGRE_COMPLETA"].copy() # this is redundant with 'source_blood', and less informative

    all_phenos_gwas_raw = sorted(set(df_raw_gwas_files_all.phenotype))
    pheno_to_ntrans = pd.Series(dict(zip(all_phenos_gwas_raw, map(lambda p: fun.get_min_ntransitons_across_aligners_one_phenotype(p, DataDir), all_phenos_gwas_raw)))).sort_values()
    gwas_phenos = set(pheno_to_ntrans[pheno_to_ntrans>=3].index)
    print("GWAS Phenos\n", "\n".join(sorted(gwas_phenos)), "n=", len(gwas_phenos))

    df_raw_gwas_files_all = df_raw_gwas_files_all[(df_raw_gwas_files_all.phenotype.isin(gwas_phenos))].copy()
    fun.save_object(df_raw_gwas_files_all, "%s/df_raw_gwas_files_all.py"%ProcessedDataDir)

    # get the benchmarking df for different GWAS parameters. It contains, for each combination of aligner, ASR_method, min_support, sample_set, filter_I info on the correct genes (unique_filterID). Only for 
    GWAS_benchmarking_df_file = fun.bechmark_GWAS_parameters(df_raw_gwas_files_all, GenoPheno_dir_all, ProcessedDataDir, DataDir)

    # explore parameters -> OLD, Now done in plots of figures_paper.py
    #fun.print_stats_GWAS_benchmarking_df_all(GWAS_benchmarking_df_file, CurDir+"/plots_paper", ProcessedDataDir)
    #fun.print_stats_different_expected_gene_sets(GWAS_benchmarking_df_file, CurDir+"/plots_paper")

    # get the filteted GWAS dfs (high and low confidence)
    df_GWAS_filt_all_file = fun.get_filtered_GWAS_dfs(GWAS_benchmarking_df_file, df_raw_gwas_files_all, GenoPheno_dir_all, ProcessedDataDir, DataDir, CurDir+"/plots_paper")

    # get the NR GWAS tables
    fun.get_NR_GWAS_tables(df_GWAS_filt_all_file, CurDir+"/plots_paper", DataDir, GenoPheno_dir_all, list_cmds, ProcessedDataDir)

######################################

##### APPLY THE PIPELINE TO THE VALIDATION DATASETS #######

# get the analysis data of the extra samples
if len(list_cmds)==0:
    print("There are %i samples in df_reads_info_extra_validation"%(len(df_reads_info_extra_validation)))
    df_reads_info_analysis_extra_validation = fun.get_df_reads_info_analysis_extra_validation(ValidationDir, VarCallDir, df_reads_info_extra_validation, DataDir, CurDir, list_cmds, simple_repeats_file, homozygous_Corthopsilosis_srrs, window_size_SNP_densities, sorted_aligners, genome_Cparapsilosis, aligner_SNP_densities, ReadsDir, min_coverage, min_pct_reads_para, min_pct_covered)
    print("There are %i samples in df_reads_info_analysis_extra_validation"%(len(df_reads_info_analysis_extra_validation)))

# get the integrated varcalls
if len(list_cmds)==0:

    # define the dataframe with all data, including the extra ones
    fields_keep = ["runID", "has_smiley"]
    df_reads_info_analysis_with_validation = pd.concat([df_reads_info_analysis, df_reads_info_analysis_extra_validation])[fields_keep].copy()
    df_reads_info_analysis_with_validation["source_data_type"] = df_reads_info_analysis_with_validation.runID.map(runID_to_sourceDataType)
    fun.check_no_nans_in_df(df_reads_info_analysis_with_validation)
    if len(df_reads_info_analysis_with_validation)!=len(set(df_reads_info_analysis_with_validation.runID)): raise ValueError("should be unique")    
    print("There are %i samples in df_reads_info_analysis_with_validation"%(len(df_reads_info_analysis_with_validation)))
    fun.save_object(df_reads_info_analysis_with_validation, ProcessedDataDir+"/df_reads_info_analysis_with_validation.py")
    
    # get the variant calls
    fun.get_integrated_variant_calls_all_datasets_with_validation(df_reads_info_analysis_with_validation, ValidationDir, VarCallDir, CurDir, list_cmds, sorted_aligners, simple_repeats_file, genome_Cparapsilosis, gff_Cparapsilosis)

# get the commands related to the geno-pheno validation
if len(list_cmds)==0:
    fun.prepare_files_for_geno_pheno_validation(df_reads_info_analysis_with_validation, df_phenotypes_validation, df_metadata_binary_phenotypes, ValidationDir, CurDir, trees_data_df, VarCallDir, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, list_cmds, simple_repeats_file, gff_Cparapsilosis)

###########################################################

#### RUN JOBS #####

# run in cluster
if len(list_cmds)>0:
    print("Running %i jobs in cluster"%(len(list_cmds)))

    # get as a df
    df_cmds = pd.DataFrame(list_cmds, columns=["name_job", "cmd", "env", "threads"])

    # define jobs dir
    jobsdir_all = "%s/running_jobs"%CurDir; fun.make_folder(jobsdir_all)

    # run one command for each name job
    for name_job in sorted(set(df_cmds.name_job)):
        print(name_job)

        # get df
        df_c = df_cmds[df_cmds.name_job==name_job]
        if len(set(df_c.env))!=1: raise ValueError("env should be the same")
        if len(set(df_c.threads))!=1: raise ValueError("threads should be the same")
        env = df_c.env.iloc[0]
        threads_per_job = df_c.threads.iloc[0]

        # define dirs
        jobsdir = "%s/%s"%(jobsdir_all, name_job); fun.make_folder(jobsdir)
        stddir = "%s/STDfiles"%jobsdir; fun.delete_folder(stddir); fun.make_folder(stddir)

        # define the jobs filename
        jobs_filename = "%s/jobs.%s"%(jobsdir, name_job)
        for f in [x for x in os.listdir(jobsdir) if x.startswith(fun.get_file(jobs_filename))]: fun.remove_file("%s/%s"%(jobsdir, f))
        open(jobs_filename, "w").write("\n".join(["source /gpfs/projects/bsc40/mschikora/anaconda3/etc/profile.d/conda.sh && conda activate %s > /dev/null && %s > %s/job.%i.std 2>&1"%(env, cmd, stddir, (I + 1)) for I, cmd in enumerate(df_c.cmd)]))
       
        # submit jobs
        os.chdir(jobsdir)        
        fun.run_jobarray_file_MN5_greasy(jobs_filename, name_job, time="02:00:00", queue="gp_debug", threads_per_job=threads_per_job, nodes=32, nsubmits=1, constraint_line="") # gp_debug, "#SBATCH --constraint=highmem"
        #fun.run_jobarray_file_MN5_greasy(jobs_filename, name_job, time="08:00:00", queue="gp_bscls", threads_per_job=threads_per_job, nodes=125, nsubmits=1, constraint_line="") # gp_debug, "#SBATCH --constraint=highmem"

    print("exiting execution since some jobs were submitted...")
    sys.exit(0)

###################