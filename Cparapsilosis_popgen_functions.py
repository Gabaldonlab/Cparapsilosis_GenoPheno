#!/usr/bin/env python

### ENV ####

print("importing functions...")

import os, shutil, subprocess, sys, itertools, json
from Bio import SeqIO
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import multiprocessing as multiproc
import scipy
from ete3 import Tree, TreeStyle, NodeStyle, RectFace, TextFace, CircleFace
import copy as cp
import matplotlib.patches as mpatches
import scipy, re
from collections import Counter
from matplotlib import colors
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pandas.errors import PerformanceWarning
import sklearn.metrics as sk_metrics
import matplotlib.colors as mcolors
from matplotlib import cm as matplotlib_cm
from matplotlib_venn import venn3
from matplotlib_venn import venn2
from io import StringIO

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    

else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"

# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir # I used version 1.02.8
perSVade_wrapper = "%s/perSVade"%perSVade_dir # I used version 1.03.1
get_trimmed_reads_for_srr_py = "%s/get_trimmed_reads_for_srr.py"%perSVade_dir

# define SRRs which are wrong due to various reasons
wrong_srrs = {"ERR9707015"}

"""
ERR9707015: C. metapsilosis, need to discard.
"""

# define example srrs
# srrs_Corthopsilosis = {"SRR3547357", "SRR3547474", "SRR3547496", "SRR3547626", "SRR3547629"} # one each clade from Cmine
# srrs_Cmetapsilosis = {}

# define mtDNA
mtDNA_chrom = "mito_C_parapsilosis_CDC317"

# define parms
Cpara_gDNA_code = 12
Cpara_mtDNA_code = 4

# deinfe Nans
vcf_strings_as_NaNs = ['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'n/a', 'nan', 'null']

# GWAS filtering parms
manual_curation_gwas_filt_parms = {"min_support":70, "correction_all_groups":False, "aligner":"hisat2_no_spliced", "ASR_method":"MPPA,DOWNPASS", "sample_set":"balanced_samples", "min_epsilon":0.3, "min_nodes_GenoAndPheno":2, "type pvalues general":"only_p_GP_bonf", "pvalue_threshold":0.05} # validated

manual_curation_gwas_filt_parms_low_conf = {"min_support":10, "correction_all_groups":False, "aligner":"hisat2_no_spliced", "ASR_method":"MPPA", "sample_set":"all_samples", "min_epsilon":0.2, "min_nodes_GenoAndPheno":2, "type pvalues general":"only_p_GP<0.05", "pvalue_threshold":0.05} # validated

# define visuals
source_data_to_color = {"sequenced_here" : "black", "SRA":"gray"}
country_to_color = {'Canada': 'black', 'China': 'red', 'France': 'cyan', 'Germany': 'orange', 'Korea': 'peru', 'Kuwait': 'olive', 'Malaysia': 'lightgreen', 'Portugal': 'green', 'Puerto Rico': 'salmon', 'Spain': 'yellow', 'Turkey': 'dodgerblue', 'UK': 'fuchsia', 'USA': 'blue', 'unknown': 'gray'}

hospital_short_to_color = {'COMPLEJO ASISTENCIAL UNIVERSITARIO DE BURGOS (Burgos)': 'sienna', 
                           'H. CLINICO UNIVERS DE VALLADOLID (Valladolid)': 'black', 
                           'H. UNIVERS MARQUES DE VALDECILLA (Santander)': 'olive',
                           'H. VIRGEN DEL ROCIO (Sevilla)': 'gold',
                           'H. SAN JUAN DE DIOS DEL ALJARAFE (Sevilla)':'orange',
                           'COMPLEJO H.ARIO VIRGEN MACARENA (Sevilla)':'yellow',

                           'H. 12 DE OCTUBRE (Madrid)': 'red', 
                           'H. CENTRAL DE LA DEFENSA GOMEZ ULLA (Madrid)': 'fuchsia', 
                           'H. CLINICO SAN CARLOS (Madrid)': 'salmon',  
                           'H. DE MOSTOLES (Madrid)': 'purple', 
                           'H. PUERTA DE HIERRO (MAJADAHONDA) (Madrid)': 'orange', 

                           "H. UNIVERSITARI VALL D'HEBRON (Barcelona)": 'blue', 
                           'H. UNIVERSITARIO DE BELLVITGE (Barcelona)': 'cyan',  
                           'LABORATORIO REFERENCIA CATALUNYA (EL PRAT DE LLOBREGAT) (Prat de Llobregat)': 'teal',
                           'H. DE SANT JOAN DE DEU (Barcelona)': 'c'}

# breakpoints of R
drug_to_MIC_breakpoints = {"FLC":(2.0, 8.0) , "ITR":(0.125, 0.5), "VRC":(0.125, 0.5), "POS":(0.06, 0.25), "IVZ":(0.125, 0.5), "AMB":(0.125, 0.5), "5FC":(0.125, 0.5), "CAS":(1.0, 4.0), "MIF":(0.5, 2.0), "ANI":(1.0, 4.0)} 
disinfectant_to_breakpoints = {"ETHANOL":(0.5, 2), "HYPLOCHLORITE":(0.1, 0.5), "H202":(50, 75), "CHLORHEXIDINE":(5, 10), "ZWITTERGENT 3-14":(10, 50), "SURFANIO":(1, 5)}

############


##### FUNCTIONS ######

def sra_run_is_clinical_or_enviromental_WGS(r):

    """Returns whether the run is clinical. Strange things:
    - I could not classify most sample types: sample_type ['Cell culture', 'Whole Organism', 'cell culture', 'clinical isolate', 'clonal culture', 'isolate', 'laboratory strain', 'nan'] 
    """

    # filter out some iv evolved samples that had r.BioSampleModel == 'Microbe\\, viral or environmental'
    if r.Run in {"SRR7898455", "SRR7898456", "SRR7898458", "SRR12333222", "SRR12333223", "SRR12333224"}: return False


    # cfilter on BioSampleModel
    if r.BioSampleModel in {'MIGS.eu,MIGS/MIMS/MIMARKS.human-associated', 'Microbe\\, viral or environmental', 'Pathogen.cl', 'Pathogen.env'}: return True

    # based on source
    if r.isolation_source in {'Blood Culture', 'Clinic', 'Host blood', 'Host nose', 'Host respiratory tract', 'Human gut', 'NRRL Collection', 'Peritoneal fluid', 'Urine sample', 'bile', 'blood', 'catheter-related bloodstream infection', 'clinical isolate', 'cobweb', 'hospital', 'isolated from a patient with diarrhea in Puerto Rico by Ashford (1928)', 'soil sample', 'urine'}: return True
    
    elif r.isolation_source in {'Lab evolution experiment'}: return False

    # based on host
    if r.Host in {'"Homo sapiens"', 'Homo sapiens'}: return True

    if r.Host_disease in ['Bloodstream infection', 'Candida infection', 'Candida parapsilosis infection', 'Candidiaisis or Colonization', 'Catheter-related bloodstream infection (CRBSI)', 'Chronic kidney disease', 'Kidney graft failure', 'Maligancy of soft tissue', 'Sepsis', 'Tropical sprue', 'candidal disease', 'severe diarrhea\\, sepsis']: return True

    # based on sample_type
    if r.sample_type in {'clinical isolate', 'isolate'}: return True
    elif r.sample_type in {'laboratory strain'}: return False

    # based on different cases
    if r.env_local_scale in {'MSK hosptial', 'human-associated habitat'}: return True
    if r.env_medium in {'faeces', 'human blood', 'human stool', 'orbital tissue'}: return True
    if r.env_broad_scale in {'human gut mycobiome'}: return True
    if r.Host_age in {'1.0', '79.0'}: return True
    if r.host_disease_outcome in {"Recovery"}: return True
    if r.host_disease_stage in {'Critical, sepsis', 'Severe'}: return True
    if r.host_sex in {'male'}: return True
    if r.environmental_history in {'clinical sample'}: return True
    if r.host_associated in {'yes'}: return True
    if r.host_scientific_name in {'Homo sapiens'}: return True
    if r.Source in {'clinical'}: return True

    # should not have been here
    for k in r.keys(): 
        if str(r[k])!="nan": print(k, r[k])
    raise ValueError("Invalid r")

def sra_run_get_type_of_isolate(r):

    """
    Returns the type of isolate
    """

    # cfilter on BioSampleModel
    #if r.BioSampleModel in {'MIGS.eu,MIGS/MIMS/MIMARKS.human-associated', 'Microbe\\, viral or environmental', 'Pathogen.cl', 'Pathogen.env'}: return True

    # manual cases
    if r.Run in {"SRR24010542", "SRR16974334"}: return "human_associated"
    if r.BioProject in {"PRJNA326748"}: return "human_associated"

    # based on source (missing 'NRRL Collection')
    if r.isolation_source in {'Blood Culture', 'Clinic', 'Host blood', 'Host nose', 'Host respiratory tract', 'Human gut', 'Peritoneal fluid', 'Urine sample', 'bile', 'blood', 'catheter-related bloodstream infection', 'clinical isolate', 'hospital', 'isolated from a patient with diarrhea in Puerto Rico by Ashford (1928)', 'urine'}: return "human_associated"

    elif r.isolation_source in {'cobweb', 'soil sample'}: return "environmental"
    
    # based on host
    if r.Host in {'"Homo sapiens"', 'Homo sapiens'}: return "human_associated"

    if r.Host_disease in ['Bloodstream infection', 'Candida infection', 'Candida parapsilosis infection', 'Candidiaisis or Colonization', 'Catheter-related bloodstream infection (CRBSI)', 'Chronic kidney disease', 'Kidney graft failure', 'Maligancy of soft tissue', 'Sepsis', 'Tropical sprue', 'candidal disease', 'severe diarrhea\\, sepsis']: return "human_associated"

    # based on sample_type
    if r.sample_type in {'clinical isolate'}: return "human_associated" # isolate

    # based on different cases
    if r.env_local_scale in {'MSK hosptial', 'human-associated habitat'}: return "human_associated"
    if r.env_medium in {'faeces', 'human blood', 'human stool', 'orbital tissue'}: return "human_associated"
    if r.env_broad_scale in {'human gut mycobiome'}: return "human_associated"
    if r.Host_age in {'1.0', '79.0'}: return "human_associated"
    if r.host_disease_outcome in {"Recovery"}: return "human_associated"
    if r.host_disease_stage in {'Critical, sepsis', 'Severe'}: return "human_associated"
    if r.host_sex in {'male'}: return "human_associated"
    if r.environmental_history in {'clinical sample'}: return "human_associated"
    if r.host_associated in {'yes'}: return "human_associated"
    if r.host_scientific_name in {'Homo sapiens'}: return "human_associated"
    if r.Source in {'clinical'}: return "human_associated"

    # should not have been here
    for k in r.keys(): 
        if str(r[k])!="nan": print(k, r[k])    
    raise ValueError("Invalid r")


def remove_file(f):

    if os.path.isfile(f): 

        try: run_cmd("rm %s > /dev/null 2>&1"%f)
        except: pass

def delete_folder(f):

    if os.path.isdir(f): shutil.rmtree(f)


def make_folder(f):

    if not os.path.isdir(f): os.mkdir(f)


def delete_file_or_folder(f):

    """Takes a path and removes it"""

    if os.path.isdir(f): shutil.rmtree(f)
    if os.path.isfile(f): os.unlink(f)
    
def get_dir(filename): return "/".join(filename.split("/")[0:-1])

def get_file(filename): return filename.split("/")[-1]


def file_is_empty(path): 
    
    """ask if a file is empty or does not exist """

    if not os.path.isfile(path):
        return_val = True
    elif os.stat(path).st_size==0:
        return_val = True
    else:
        return_val = False
            
    return return_val


def run_cmd_pipelines_env(cmd, env="meTAline"):

    """This function runs a cmd with a given env"""

    # define the cmds
    CondaDir = "/gpfs/projects/bsc40/project/pipelines/anaconda3"
    SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
    cmd_prefix = "%s && conda activate %s &&"%(SOURCE_CONDA_CMD, env)

    # define the running
    cmd_to_run = "%s %s"%(cmd_prefix, cmd)

    # run
    out_stat = os.system(cmd_to_run) 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd_to_run, out_stat))




def run_cmd(cmd, env="Cparapsilosis_popgen_env"):

    """This function runs a cmd with a given env"""

    # define the cmds
    CondaDir = "/gpfs/projects/bsc40/mschikora/anaconda3"
    SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
    cmd_prefix = "%s && conda activate %s &&"%(SOURCE_CONDA_CMD, env)

    # define the running
    cmd_to_run = "%s %s"%(cmd_prefix, cmd)

    # run
    out_stat = os.system(cmd_to_run) 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd_to_run, out_stat))



def download_srr_with_prefetch(srr, SRRfile, Isrr, total_srrs, replace=False):

    """This function downloads an srr file for an srr if not already done"""

    # define the downloading dir
    downloading_dir = get_dir(SRRfile)

    # make the downloading dir
    make_folder(downloading_dir)

    # define the std files
    prefetch_std = "%s.std.txt"%SRRfile
    prefetch_std_copy = "%s.copy"%prefetch_std

    # try a couple of times
    for Itry in range(2):

        if file_is_empty(SRRfile) or replace is True:
            print("running prefetch for %s (%i/%i)"%(srr, Isrr, total_srrs))

            # remove the locks of previous runs
            for file in ["%s/%s"%(downloading_dir, f) for f in os.listdir(downloading_dir) if ".lock" in f or ".tmp." in f]: remove_file(file)

            # remove the actual srr
            remove_file(SRRfile)

            # run prefetch
            print("running prefetch. The std can be found in %s"%prefetch_std)
            prefetch_cmd = "prefetch -o %s --max-size 500G --resume no --verify yes --progress %s > %s 2>&1"%(SRRfile, srr, prefetch_std)
            try: run_cmd(prefetch_cmd, env="Cparapsilosis_popgen_env")
            except: print("prefetch did not work for %s. Command: '%s'"%(srr, prefetch_cmd))

            # test that the std of prefetch states that there are no unresolved dependencies
            std_lines = open(prefetch_std, "r").readlines()
            successful_download = any(["was downloaded successfully" in l for l in std_lines])
            no_dependencies_left = any(["has 0 unresolved dependencies" in l for l in std_lines])
            has_dependencies_line = any(["unresolved dependencies" in l for l in std_lines])

            if not successful_download or (has_dependencies_line and not no_dependencies_left): 
                run_cmd("cp %s %s"%(prefetch_std, prefetch_std_copy))
                print("prefetch did not work for %s. Test the log in %s"%(srr, prefetch_std_copy))
                remove_file(SRRfile)

    # check that the prefetch works 
    if file_is_empty(SRRfile): 
        raise ValueError("prefetch did not work for %s"%srr)

    remove_file(prefetch_std)
    remove_file(prefetch_std_copy)

    return SRRfile

def check_that_CGD_and_NCBI_genomes_are_the_same(DataDir):

    """Checks genomes"""

    genome_ncbi = "%s/genomes_and_annotations/Cparapsilosis_GenBank_GCA_000182765.2_ASM18276v2.fasta"%DataDir
    genome_CGD = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir
    mtDNA_ncbi = "HE605210.1"
    mtDNA_CGD = "mito_C_parapsilosis_CDC317"
    ncbi_chrom_to_seq = {chrom.id : chrom for chrom in SeqIO.parse(genome_ncbi, "fasta")}
    CGD_chrom_to_seq = {chrom.id : chrom for chrom in SeqIO.parse(genome_CGD, "fasta")}

    for ID_ncbi, seq_ncbi in ncbi_chrom_to_seq.items():
        if ID_ncbi==mtDNA_ncbi: ID_CGD = mtDNA_CGD
        else: ID_CGD = "Contig%s_C_parapsilosis_CDC317"%(seq_ncbi.description.split()[-1])
        seq_CGD = CGD_chrom_to_seq[ID_CGD]
        print("%s, %s. Same seq: %s"%(ID_ncbi, ID_CGD, seq_ncbi.seq.upper()==seq_CGD.seq.upper()))


def get_int_or_minus1(x):
    if str(x)=="nan": return -1
    else: return int(x)

def get_str_or_unknown(x):
    if str(x)=="nan": return "unknown"
    elif type(x)==str: return x
    else: raise ValueError("invalid %s. %s"%(x, type(x)))

def get_erg11_mut_str(x):

    if x=="Y / Y132F (heterocigota)": return "Y132F (het)"
    elif x=="Y / Y132F (heterozygous)": return "Y132F (het)"
    elif str(x)=="nan": return "WT"
    elif type(x)==str: return x
    else: raise ValueError("invalid %s"%x)

def get_float_or_minus1(x):
    if str(x)=="nan": return -1
    else: return float(x)

def get_mic_float(x):

    if str(x) in {"nan", "-1.0", "*****", "-"} or x==" ": return -1
    elif all(map(lambda c: c.isdigit() or c==".", str(x))): return float(x)
    elif str(x).startswith(">"): return float(x.split(">")[1])*2
    else: raise ValueError("invalid mic: '%s'"%x)

def covert_yes_no(x): return {"sí":"yes", "Sí":"yes", "no":"no",'nan':"unknown", "unknown":"unknown", "yes":"yes", "no":"no"}[str(x)]

def get_df_metadata_final(ManualDataDir, expected_sampleIDs):

    """Get final metadata. I am taking data from the following excels:

    - 'metadata_mail_240209_table.xlsx'. This has the tabs 'ALL DATA', 'PHENOTYPES' and 'ETEST' (which have redundant data)
    - 'metadata_mail_final_table_231121.xlsx' it has the table ‘DISINFECTANTS’ with disinfectant data that is not in the one above

    """

    # check that the sampleIDs are as expected
    expected_sampleIDs_tuples = {get_numeric_tuple_sampleID(s) for s in expected_sampleIDs}
    if len(expected_sampleIDs_tuples)!=len(expected_sampleIDs): raise ValueError("sIDs not formatted")

    # load the metadata excels
    print("loading metadata")
    rename_dict_all = {"CMI ETEST 24H":'MIC ETEST 24H', 'CMI ETEST 48H':'MIC ETEST 48H', 'COLONIA DENTRO DE HALO':'Colonies inside halo', 'CMI colonia dentro 24h':'MIC of recovered heteroresistant 24 h', 'CMI colonia dentro 48h':'MIC of recovered heteroresistant 48h', 'nuevo genotipo':"genotype", "Numero paciente duplicado":"duplicated_patientID", 'FLC (S/I/R)':"FLC_resistance", 'VOR (S/I/R)':"VRC_resistance", "ITRA (S/R)":"ITR_resistance", 'POSA (S/R)':"POS_resistance", "anfo":"AMB_MIC", "5FC":"5FC_MIC", "fluco":"FLC_MIC", "itra":"ITR_MIC", "vori":"VRC_MIC", "posa":"POS_MIC", "isa":"IVZ_MIC", "caspo":"CAS_MIC", "mica":"MIF_MIC", "anidula":"ANI_MIC", "Fecha toma":"date_collection", "Año toma":"year_collection", 'fecha nacimiento':"date_birth", 'Edad':"age" , "colección":"collection"}

    df_metadata_all = pd.read_excel("%s/metadata_mail_240209_table.xlsx"%ManualDataDir, 'ALL DATA').rename(columns=rename_dict_all)
    df_etest = pd.read_excel("%s/metadata_mail_240209_table.xlsx"%ManualDataDir, 'ETEST').rename(columns={"Strain number":"Labo"})
    df_phenos = pd.read_excel("%s/metadata_mail_240209_table.xlsx"%ManualDataDir, 'PHENOTYPES').rename(columns={"Strain number":"Labo"})
    df_disinfectants = pd.read_excel("%s/metadata_mail_final_table_231121.xlsx"%ManualDataDir, 'DISINFECTANTS').rename(columns={"Strain number":"Labo"})

    # remove some fields
    df_metadata_all = df_metadata_all.drop(['SECUENCIADA', 'Labo.2', 'Labo.1', 'nuevo genotipo.1', 'peticion'], axis=1)

    # add sample ID and general checks
    for type_df, df in [("all",df_metadata_all), ("etest",df_etest), ("pheno",df_phenos), ("dis",df_disinfectants)]:

        # checks unique
        if len(df)!=len(set(df.Labo)): raise ValueError("Labo is not unique")

        # add sampleID
        df["sampleID"] = df.Labo.apply(lambda x: "-".join(x.split("/")).upper()) # this is how the reads are called 
        if len(df)!=len(set(df["sampleID"])): raise ValueError("Labo should be unique")

        # rename the sampleIDs to match those from the cgenomics files
        sID_to_new_sID = {'L21-066' : 'L21-66', 'L21-062' : 'L21-62',  'L21-065':'L21-65'}

        def get_renamed_sampleID(r):
            if r.sampleID in sID_to_new_sID: return sID_to_new_sID[r.sampleID]
            else: return r.sampleID

        df["sampleID"] = df.apply(get_renamed_sampleID, axis=1)

        # get strange and missing samples
        strange_samples = sorted(set(df.sampleID).difference(expected_sampleIDs))
        missing_samples = sorted(expected_sampleIDs.difference(set(df.sampleID)))

        # special checks
        if type_df in {"all", "dis"} and len(strange_samples)>0: raise ValueError("there can't be strange samples in all,dis")
        if type_df=="all" and len(missing_samples)>0: raise ValueError("there can't be missing samples in all")

        if any([get_numeric_tuple_sampleID(s) in expected_sampleIDs_tuples for s in strange_samples]): # to validate that names are correct
            raise ValueError("there are some strange samples that are malformed")

        if len(df)!=len(set(df.sampleID)): raise ValueError("sampleID is not unique")

    # filter and keep fields
    df_etest = df_etest[df_etest.sampleID.isin(expected_sampleIDs)].copy()
    df_phenos = df_phenos[df_phenos.sampleID.isin(expected_sampleIDs)].copy()

    # change indices
    df_disinfectants = df_disinfectants.set_index("sampleID", drop=False)
    df_etest = df_etest.set_index("sampleID", drop=False)
    df_phenos = df_phenos.set_index("sampleID", drop=False)
    df_metadata_all = df_metadata_all.set_index("sampleID", drop=False)

    # modify pheno fields, checking that they are all the same in different dfs
    pheno_fiels = ['MIC ETEST 24H', 'MIC ETEST 48H', 'Colonies inside halo', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h', 'Agar invasiveness', 'MORPHOLOGY', 'Microfluidic behavior', 'Biofilm (High / Low / No)'] + ['% VIABILITY ETHANOL', '% VIABILITY HYPLOCHLORITE', '% VIABILITY H202', '% VIABILITY CHLORHEXIDINE', '% VIABILITY ZWITTERGENT 3-14', '% VIABILITY SURFANIO']

    for f in pheno_fiels:

        # define phenos df
        if f in {'Agar invasiveness', 'MORPHOLOGY', 'Microfluidic behavior', 'Biofilm (High / Low / No)'}: dfs_pheno = [df_metadata_all, df_etest, df_phenos]
        elif f in {'MIC ETEST 24H', 'MIC ETEST 48H', 'Colonies inside halo', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h'}: dfs_pheno = [df_metadata_all, df_etest]
        elif f.startswith("% VIABILITY"): dfs_pheno = [df_disinfectants]
        else: raise ValueError("invalid f %s"%f)

        # map each field to a function
        if f.startswith("% VIABILITY"): fn_modifier = get_float_or_minus1
        elif f in {'Agar invasiveness', 'MORPHOLOGY', 'Microfluidic behavior', 'Biofilm (High / Low / No)'}: fn_modifier = get_str_or_unknown
        elif f in {'MIC ETEST 24H', 'MIC ETEST 48H', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h'}: fn_modifier = get_mic_float
        elif f in {'Colonies inside halo'}: fn_modifier = covert_yes_no
        else: raise ValueError("invalid f %s"%f)

        # map each sample to a value
        dicts_pheno = [dict(df[f].apply(str).apply(fn_modifier)) for df in dfs_pheno]

        sampleID_to_value = {}
        for sID in expected_sampleIDs:
            all_vals = {d[sID] for d in dicts_pheno if sID in d}

            # debug
            if len(all_vals)==0: continue

            # try to remove some cases
            if len(all_vals)>1: all_vals = all_vals.difference({"unknown", -1})

            # debug
            if len(all_vals)!=1: raise ValueError("%s: %s"%(sID, all_vals))

            # keep
            sampleID_to_value[sID] = list(all_vals)[0]

        df_metadata_all[f] = df_metadata_all.sampleID.map(sampleID_to_value).map(str).map(fn_modifier)

    # add ERG11_mutation
    def get_corrected_var(x):
        if x=="Y / Y132F (heterocigota)": return "Y132F (het)"
        else: return x

    def get_erg11_mut_from_r(r):
        r = r[~pd.isna(r)]
        if len(r)==0: raise ValueError("no muts")
        all_muts = set(map(get_corrected_var, set(r.values).difference({"WT", "WT ", " "})))
        strange_muts = all_muts.difference({"Y132F", "G458S", "R398I", "Y132F (het)", "K143R"})
        if len(strange_muts)>0: raise ValueError("strange_muts: %s"%strange_muts)

        if len(all_muts)==0: return "REF"
        else: return ";".join(sorted(all_muts))

    df_metadata_all["ERG11_mutation"] = df_metadata_all[['Mutacion 132', 'Mutacion 143', 'Mutacion 458', 'Mutacion 398']].apply(get_erg11_mut_from_r, axis=1)

    # add mic fields
    for mic_f in ['AMB_MIC', '5FC_MIC', 'FLC_MIC', 'ITR_MIC', 'VRC_MIC', 'POS_MIC', 'IVZ_MIC', 'CAS_MIC', 'MIF_MIC', 'ANI_MIC']:
        df_metadata_all[mic_f] = df_metadata_all[mic_f].apply(get_mic_float)

    # change the R/S fields
    def get_dict_same_fields(x): return {v:v for v in x}
    for f in ['FLC_resistance', 'VRC_resistance', 'ITR_resistance', 'POS_resistance']:
        strange_vals = set(df_metadata_all[f]).difference({"R", "S", "I"})
        df_metadata_all[f] = df_metadata_all[f].apply(str).apply(get_str_or_unknown)
        if len(strange_vals)>0: raise ValueError("strange_vals: %s"%strange_vals)

    # check specific values
    for f, expected_vals in [('Colonies inside halo',{"yes", "no"}), ('Agar invasiveness', {"Yes", "No", "Low"}), ("MORPHOLOGY", {"Yeast", "Psedohyphae"}), ("Microfluidic behavior", {"Little clump", "Adherent", "Big Clumps"}), ("Biofilm (High / Low / No)", {"NF", "LF", "HF"})]:
        strange_vals = (set(df_metadata_all[f]).difference({"unknown"})).difference(expected_vals)
        if len(strange_vals)>0: raise ValueError("strange_vals: %s"%strange_vals)


    # change str fields
    for f in ["Labo", "genotype", "date_collection", "muestra", "date_birth", "Observaciones", "hospital", "CCAA", 'Referencia cepa Hospital']:
        df_metadata_all[f] = df_metadata_all[f].apply(str).apply(get_str_or_unknown)

    # change int fields
    for f in ['duplicated_patientID', "year_collection", "age", "collection"]:
        df_metadata_all[f] = df_metadata_all[f].apply(get_int_or_minus1)


    # remove fields
    df_metadata_all = df_metadata_all.drop(['Mutacion 132', 'Mutacion 143', 'Mutacion 458', 'Mutacion 398'], axis=1)

    # add the city
    hospital_to_city = {'COMPLEJO ASISTENCIAL UNIVERSITARIO DE BURGOS': 'Burgos', "H. UNIVERSITARI VALL D'HEBRON": 'Barcelona', 'HOSPITAL 12 DE OCTUBRE': 'Madrid', 'HOSPITAL CENTRAL DE LA DEFENSA GOMEZ ULLA': 'Madrid', 'HOSPITAL CLINICO SAN CARLOS': 'Madrid', 'HOSPITAL CLINICO UNIVERS DE VALLADOLID': 'Valladolid', 'HOSPITAL DE MOSTOLES': 'Madrid', 'HOSPITAL PUERTA DE HIERRO (MAJADAHONDA)': 'Madrid', 'HOSPITAL UNIVERS MARQUES DE VALDECILLA': 'Santander', 'HOSPITAL UNIVERSITARIO DE BELLVITGE': 'Barcelona', 'HOSPITAL VIRGEN DEL ROCIO': 'Sevilla', 'LABORATORIO REFERENCIA CATALUNYA (EL PRAT DE LLOBREGAT)': 'Prat de Llobregat', 'unknown':'unknown', 'COMPLEJO HOSPITALARIO UNIVERSITARIO DE ALBACETE':"Albacete", 'HOSPITAL DE LA LINEA':"Cadiz", 'HOSPITAL UNIVERSITARIO DEL RIO HORTEGA':"Valladolid", 'HOSPITAL DE SANT JOAN DE DEU':"Barcelona", 'HOSPITAL SAN JUAN DE DIOS DEL ALJARAFE':"Sevilla", 'COMPLEJO HOSPITALARIO VIRGEN MACARENA':"Sevilla", 'HOSPITAL DE CABUEÑES':"Gijón"}

    df_metadata_all["city_collection"] = df_metadata_all.hospital.apply(lambda x: hospital_to_city[x])

    # rename the CCAA
    renamed_ccaa = {'ANDALUCIA':"ANDALUCIA", 'ASTURIAS':"ASTURIAS", 'CANTABRIA':"CANTABRIA", 'CASTILLA Y LEON':'CASTILLA Y LEON', 'CATALUÑA':'CATALUNYA', 'COMUNIDAD DE MADRID':'COMUNIDAD DE MADRID', 'Cantabria':'CANTABRIA', 'Cataluña':'CATALUNYA', 'MADRID':'COMUNIDAD DE MADRID', 'Madrid':'COMUNIDAD DE MADRID'}
    df_metadata_all["CCAA"] = df_metadata_all["CCAA"].apply(lambda x: renamed_ccaa[x])

    # check
    check_no_nans_in_df(df_metadata_all)

    return df_metadata_all

def get_df_metadata(ManualDataDir):

    """Gets the metadata of OSCAR data"""

    this_does_not_work

    # get all dfs
    df_metadata_230626 = pd.read_excel("%s/metadata_mail_230626.xlsx"%ManualDataDir)
    df_metadata_230705_initial = pd.read_excel("%s/metadata_mail_230705.xlsx"%ManualDataDir) # from a mail sent first
    df_metadata_230705 = pd.read_excel("%s/metadata_mail_230705_resent_231123.xlsx"%ManualDataDir).rename(columns={'Mutacion ERG11':"ERG11_mutation"})
    df_metadata_231121_phenotypes = pd.read_excel("%s/metadata_mail_final_table_231121.xlsx"%ManualDataDir, sheet_name="PHENOTYPIC STUDY").rename(columns={"Strain number":"Labo"})
    df_metadata_231121_disinfectants = pd.read_excel("%s/metadata_mail_final_table_231121.xlsx"%ManualDataDir, sheet_name="DISINFECTANTS").rename(columns={"Strain number":"Labo"})

    # check differences between strains (Labo)
    if sorted(df_metadata_230626.Labo)!=sorted(df_metadata_230705.Labo): raise ValueError("different samples")
    if sorted(df_metadata_230626.Labo)!=sorted(df_metadata_230705_initial.Labo): raise ValueError("different samples")

    # check that, within 231121, all strains from disinfectants are also in the phenotypes
    strains_disinfectants_not_phenotypes = set(df_metadata_231121_disinfectants.Labo).difference(set(df_metadata_231121_phenotypes.Labo))
    if len(strains_disinfectants_not_phenotypes)>0: raise ValueError("strains_disinfectants_not_phenotypes:%s"%strains_disinfectants_not_phenotypes)
    if any(df_metadata_231121_phenotypes["STRAIN NUMBER"]!=df_metadata_231121_phenotypes.Labo): raise ValueError("some Labos are different")

    # print the numbers of strains that are in one dataframe but not the otehrs
    df_metadata_231121_not_230705 = set(df_metadata_231121_phenotypes.Labo).difference(set(df_metadata_230705.Labo))
    df_metadata_230705_not_231121 = set(df_metadata_230705.Labo).difference(set(df_metadata_231121_phenotypes.Labo))
    df_metadata_230705_and_231121 = set(df_metadata_230705.Labo).intersection(set(df_metadata_231121_phenotypes.Labo))

    print("df_metadata_231121_not_230705", len(df_metadata_231121_not_230705))
    print("df_metadata_230705_not_231121", len(df_metadata_230705_not_231121))
    print("df_metadata_230705_and_231121", len(df_metadata_230705_and_231121))

    # change specific fields
    df_metadata_230705["ERG11_mutation"] = df_metadata_230705["ERG11_mutation"].apply(get_erg11_mut_str)
    df_metadata_231121_phenotypes["ERG11_mutation"] = df_metadata_231121_phenotypes["ERG11 MUTATION"].apply(get_erg11_mut_str)

    # format the dates
    for df in [df_metadata_230626, df_metadata_230705]:
        df["date_collection"] = df["Fecha toma"].apply(str)
        df["year_collection"] = df["Año toma"].apply(get_int_or_minus1)

    df_metadata_231121_phenotypes["date_collection"] = df_metadata_231121_phenotypes['Sample isolation date'].apply(str)
    df_metadata_231121_phenotypes["year_collection"] = df_metadata_231121_phenotypes['Year'].apply(get_int_or_minus1)

    # change the fluconazole / itra / posa
    for df in [df_metadata_230626, df_metadata_230705, df_metadata_231121_phenotypes]:
        for f in ["FLC (S/I/R)", "ITRA (S/R)", "POSA (S/R)", "VOR (S/I/R)"]: df[f] = df[f].apply(get_str_or_unknown)

    # init df with Labo
    sorted_labos = sorted(set(df_metadata_230705.Labo).union(set(df_metadata_231121_phenotypes.Labo)))
    print("df_metadata has %i labos"%len(sorted_labos))    
    df_metadata = pd.DataFrame({"Labo": sorted_labos})
    df_metadata["sampleID"] = df_metadata.Labo.apply(lambda x: "-".join(x.split("/")).upper()) # this is how the reads are called 
    if len(df_metadata)!=len(set(df_metadata["sampleID"])): raise ValueError("Labo should be unique")

    # rename the sampleIDs to match those from the cgenomics files
    sID_to_new_sID = {'L21-066' : 'L21-66', 'L21-062' : 'L21-62',  'L21-065':'L21-65'} # 'L21-282',

    def get_renamed_sampleID(r):
        if r.sampleID in sID_to_new_sID: return sID_to_new_sID[r.sampleID]
        else: return r.sampleID

    df_metadata["sampleID"] = df_metadata.apply(get_renamed_sampleID, axis=1)


    # add metadata, either from 231121 (first) or 230705 (second) 
    for k in ["hospital", "year_collection", "FLC (S/I/R)", "ITRA (S/R)", "POSA (S/R)", "VOR (S/I/R)", "ERG11_mutation"]:

        # get the values
        labo_to_value_230705 = dict(df_metadata_230705.set_index("Labo")[k])
        labo_to_value_231121 = dict(df_metadata_231121_phenotypes.set_index("Labo")[k])

        # get the value in a hierarchical manner
        def get_val_for_k_one_row(l):

            if l in labo_to_value_230705 and l in labo_to_value_231121: 

                if labo_to_value_230705[l]!=labo_to_value_231121[l]: raise ValueError("for %s l=%s, 230705 is %s and 231121 is %s "%(k, l, labo_to_value_230705[l], labo_to_value_231121[l]))
                return labo_to_value_230705[l] 

            elif l in labo_to_value_230705: return labo_to_value_230705[l]

            elif l in labo_to_value_231121: return labo_to_value_231121[l]

            else: raise ValueError("invalid %s %s"%(k, l))


        df_metadata[k] = df_metadata.Labo.apply(get_val_for_k_one_row)

    # add the city
    hospital_to_city = {'COMPLEJO ASISTENCIAL UNIVERSITARIO DE BURGOS': 'Burgos', "H. UNIVERSITARI VALL D'HEBRON": 'Barcelona', 'HOSPITAL 12 DE OCTUBRE': 'Madrid', 'HOSPITAL CENTRAL DE LA DEFENSA GOMEZ ULLA': 'Madrid', 'HOSPITAL CLINICO SAN CARLOS': 'Madrid', 'HOSPITAL CLINICO UNIVERS DE VALLADOLID': 'Valladolid', 'HOSPITAL DE MOSTOLES': 'Madrid', 'HOSPITAL PUERTA DE HIERRO (MAJADAHONDA)': 'Madrid', 'HOSPITAL UNIVERS MARQUES DE VALDECILLA': 'Santander', 'HOSPITAL UNIVERSITARIO DE BELLVITGE': 'Barcelona', 'HOSPITAL VIRGEN DEL ROCIO': 'Sevilla', 'LABORATORIO REFERENCIA CATALUNYA (EL PRAT DE LLOBREGAT)': 'Prat de Llobregat', 'unknown':'unknown', 'COMPLEJO HOSPITALARIO UNIVERSITARIO DE ALBACETE':"Albacete", 'HOSPITAL DE LA LINEA':"Cadiz", 'HOSPITAL UNIVERSITARIO DEL RIO HORTEGA':"Valladolid"}

    df_metadata["city_collection"] = df_metadata.hospital.apply(lambda x: hospital_to_city[x])

    # check
    check_no_nans_in_df(df_metadata)

    return df_metadata


def get_df_metadata_first_2_runs(ManualDataDir):

    """Gets the metadata"""

    this_does_not_work

    # get all dfs
    df_metadata_230626 = pd.read_excel("%s/metadata_mail_230626.xlsx"%ManualDataDir)
    #df_metadata_230705 = pd.read_excel("%s/metadata_mail_230705.xlsx"%ManualDataDir) # from a mail sent first
    df_metadata_230705 = pd.read_excel("%s/metadata_mail_230705_resent_231123.xlsx"%ManualDataDir).rename(columns={'Mutacion ERG11':"ERG11_mutation"})

    df_metadata_231121_phenotypes = pd.read_excel("%s/metadata_mail_final_table_231121.xlsx"%ManualDataDir, sheet_name="PHENOTYPIC STUDY").rename(columns={"Strain number":"Labo"})
    df_metadata_231121_disinfectants = pd.read_excel("%s/metadata_mail_final_table_231121.xlsx"%ManualDataDir, sheet_name="DISINFECTANTS").rename(columns={"Strain number":"Labo"})

    # check differences between strains (Labo)
    if sorted(df_metadata_230626.Labo)!=sorted(df_metadata_230705.Labo): raise ValueError("different samples")

    print(df_metadata_230705, df_metadata_231121_phenotypes, df_metadata_231121_disinfectants)

    strains_disinfectants_not_phenotypes = set(df_metadata_231121_disinfectants.Labo).difference(set(df_metadata_231121_phenotypes.Labo))
    if len(strains_disinfectants_not_phenotypes)>0: raise ValueError("strains_disinfectants_not_phenotypes:%s"%strains_disinfectants_not_phenotypes)

    df_metadata_231121_not_230705 = set(df_metadata_231121_phenotypes.Labo).difference(set(df_metadata_230705.Labo))
    df_metadata_230705_not_231121 = set(df_metadata_230705.Labo).difference(set(df_metadata_231121_phenotypes.Labo))
    df_metadata_230705_and_231121 = set(df_metadata_230705.Labo).intersection(set(df_metadata_231121_phenotypes.Labo))

    print("df_metadata_231121_not_230705", len(df_metadata_231121_not_230705))
    print("df_metadata_230705_not_231121", len(df_metadata_230705_not_231121))
    print("df_metadata_230705_and_231121", len(df_metadata_230705_and_231121))

    # change specific fields
    df_metadata_230705["ERG11_mutation"] = df_metadata_230705["ERG11_mutation"].apply(get_erg11_mut_str)

    # change some fields to 
    for df in [df_metadata_230626, df_metadata_230705]:
        df["date_collection"] = df["Fecha toma"].apply(str)
        df["year_collection"] = df["Año toma"].apply(get_int_or_minus1)
        
        for f in ["FLC (S/I/R)", "ITRA (S/R)", "POSA (S/R)", "VOR (S/I/R)"]: df[f] = df[f].apply(get_str_or_unknown)

    # check

    # differences with 

    # init df with sample name, formated as in 
    df_metadata = pd.DataFrame({"Labo": sorted(set(df_metadata_230626.Labo))})
    df_metadata["sampleID"] = df_metadata.Labo.apply(lambda x: "-".join(x.split("/")).upper())
    if len(df_metadata)!=len(set(df_metadata["sampleID"])): raise ValueError("Labo should be unique")

    # rename the sampleIDs to match those from the cgenomics files
    sID_to_new_sID = {'L21-066' : 'L21-66', 'L21-062' : 'L21-62',  'L21-065':'L21-65'} # 'L21-282',

    def get_renamed_sampleID(r):
        if r.sampleID in sID_to_new_sID: return sID_to_new_sID[r.sampleID]
        else: return r.sampleID

    df_metadata["sampleID"] = df_metadata.apply(get_renamed_sampleID, axis=1)

    #print("df_metadata_230705", list(df_metadata_230705.keys()))
    #print("df_metadata_230626", list(df_metadata_230626.keys()))
    print("common", sorted(set(df_metadata_230705.keys()).intersection(set(df_metadata_230626.keys()))))
    print("unique df_metadata_230705",  sorted(set(df_metadata_230705.keys()).difference(set(df_metadata_230626.keys()))))
    print("unique df_metadata_230626",  sorted(set(df_metadata_230626.keys()).difference(set(df_metadata_230705.keys()))))

    print("\n", "values for fields unique in df_metadata_230705")
    for f in sorted(set(df_metadata_230705.keys()).difference(set(df_metadata_230626.keys()))):
        print(f, sorted(set(df_metadata_230705[f].apply(str))))
    
    print("\n", "values for fields unique in df_metadata_230626")
    for f in sorted(set(df_metadata_230626.keys()).difference(set(df_metadata_230705.keys()))):
        print(f, sorted(set(df_metadata_230626[f].apply(str))))


    print("\n", "values for common fields")
    for f in sorted(set(df_metadata_230626.keys()).intersection(set(df_metadata_230705.keys()))):
        print(f, sorted(set(df_metadata_230626[f].apply(str))), sorted(set(df_metadata_230705[f].apply(str))))

    # add metadata
    print("\n\n")
    for k in ["hospital", "year_collection", "FLC (S/I/R)", "ITRA (S/R)", "POSA (S/R)", "VOR (S/I/R)"]:
    # for k in sorted(set(df_metadata_230705.keys()).intersection(set(df_metadata_230626.keys())).difference({"Labo"})): # common fields
        print("\n", k)

        labo_to_value_230626 = df_metadata_230626.set_index("Labo")[k].sort_index()
        labo_to_value_230705 = df_metadata_230705.set_index("Labo")[k].sort_index()

        # check
        if any(labo_to_value_230626!=labo_to_value_230705): 
            different_labos = sorted(labo_to_value_230626[labo_to_value_230626!=labo_to_value_230705].index)

            if not all([pd.isna(labo_to_value_230626.loc[l]) and pd.isna(labo_to_value_230705.loc[l]) for l in different_labos]):


                #print("labo_to_value_230626", labo_to_value_230626.loc[different_labos])
                #print("labo_to_value_230705", labo_to_value_230705.loc[different_labos])
                print ("WARNING: the values for %s are not all the same. 230626 : %s, 230705: %s"%(k,  {dl : labo_to_value_230626.loc[dl]  for dl in different_labos}, {dl : labo_to_value_230705.loc[dl]  for dl in different_labos}))

        df_metadata[k] = df_metadata.Labo.apply(lambda x: labo_to_value_230705[x])


    # add metadata specific from df_metadata_230705
    for k in ["ERG11_mutation"]:
        labo_to_value_230705 = df_metadata_230705.set_index("Labo")[k].sort_index()
        df_metadata[k] = df_metadata.Labo.apply(lambda x: labo_to_value_230705[x])

    # add strains for which I don't have metadata
    df_21_282 = pd.DataFrame({"sampleID":["L21-282"], "Labo":["L21-282"]})
    df_21_282["year_collection"] = -1
    for f in ["FLC (S/I/R)", "ITRA (S/R)", "POSA (S/R)", "VOR (S/I/R)", "ERG11_mutation", "hospital"]: df_21_282[f] = "unknown"
    df_metadata = pd.concat([df_metadata, df_21_282]).reset_index(drop=True)

    # add the city
    hospital_to_city = {'COMPLEJO ASISTENCIAL UNIVERSITARIO DE BURGOS': 'Burgos', "H. UNIVERSITARI VALL D'HEBRON": 'Barcelona', 'HOSPITAL 12 DE OCTUBRE': 'Madrid', 'HOSPITAL CENTRAL DE LA DEFENSA GOMEZ ULLA': 'Madrid', 'HOSPITAL CLINICO SAN CARLOS': 'Madrid', 'HOSPITAL CLINICO UNIVERS DE VALLADOLID': 'Valladolid', 'HOSPITAL DE MOSTOLES': 'Madrid', 'HOSPITAL PUERTA DE HIERRO (MAJADAHONDA)': 'Madrid', 'HOSPITAL UNIVERS MARQUES DE VALDECILLA': 'Santander', 'HOSPITAL UNIVERSITARIO DE BELLVITGE': 'Barcelona', 'HOSPITAL VIRGEN DEL ROCIO': 'Sevilla', 'LABORATORIO REFERENCIA CATALUNYA (EL PRAT DE LLOBREGAT)': 'Prat de Llobregat', 'unknown':'unknown'}

    df_metadata["city_collection"] = df_metadata.hospital.apply(lambda x: hospital_to_city[x])

    # check
    check_no_nans_in_df(df_metadata)

    return df_metadata


def get_numeric_tuple_sampleID(sID):

    """Gets a sampleID and returns a numeric value"""

    if not sID.startswith("L"): raise ValueError("invalid %s"%sID)

    sID_split = sID.split("-")
    if len(sID_split)!=2: raise ValueError("invalid %s"%sID)

    year_ID = int(sID_split[0].lstrip("L"))
    if year_ID not in {19, 20, 21, 22, 23}: raise ValueError("invalid %s"%sID)
    sID_list = [year_ID]

    if all([x.isdigit() for x in sID_split[1]]): sID_list.append(int(sID_split[1]))
    elif not sID_split[1][-1].isdigit() and all([x.isdigit() for x in sID_split[1][0:-1]]): sID_list += [int(sID_split[1][0:-1]), sID_split[1][-1]]
    else: raise ValueError("invalid %s"%sID)

    return tuple(sID_list)
    # return (, int(sID_split[1]))

def get_data_source_BioProjects_Bergin2022(x):

    """Gets the bioprohect of Bergin 2022"""

    if x.startswith("PRJNA") or x.startswith("PRJE"): return {x.replace(" ", "")}
    elif x=="This study ": return {"PRJNA795920", "PRJNA748054"}
    else: raise ValueError("invalid %s"%x)


def get_city_and_country_collection_Bergin2022(r):

    """Adds to the df the city and country of the colleection"""

    loc_to_data = {'Belgium  ': ('Belgium', 'unknown'), 'Bergamo, Italy  ': ('Italy', 'Bergamo'), 'Durham (Duke), NC': ('USA', 'Durham'), 'Germany  ': ('Germany', 'unknown'), 'Ireland  ': ('Ireland', 'unknown'), 'Italy  ': ('Italy', 'unknown'), 'Japan  ': ('Japan', 'unknown'), 'Korea  ': ('Korea', 'unknown'), 'Kuwait ': ('Kuwait', 'unknown'), 'Leeds, UK ': ('UK', 'Leeds'), 'Lisbon, Portugal ': ('Portugal', 'Lisbon'), 'London, UK ': ('UK', 'London'), 'London, UK  ': ('UK', 'London'), 'London,UK  ': ('UK', 'London'), 'Mississippi, USA ': ('USA', 'Mississippi'), 'Nantes, France ': ('France', 'Nantes'), 'New York': ('USA', 'New York'), 'New York ': ('USA', 'New York'), 'Plymouth MN, USA  ': ('USA', 'Plymouth'), 'Porto, Portugal  ': ('Portugal', 'Porto'), 'Puerto Rico ': ('Puerto Rico', 'unknown'), 'San Jose, USA  ': ('USA', 'San Jose'), 'Tuscany  ': ('Italy', 'Tuscany'), 'UM-Madison': ('USA', 'Madison'), 'USA': ('USA', 'unknown'), 'USA  ': ('USA', 'unknown'), 'Unknown  ': ('unknown', 'unknown')}

    r["country_collection"] = loc_to_data[r.location][0]
    r["city_collection"] = loc_to_data[r.location][1]

    if not r.country_collection in {"Belgium", "USA", "Italy", "Germany", "Ireland", "Japan", "Kuwait", "Korea", "UK", "Portugal", "France", "Puerto Rico", "unknown"}: raise ValueError(r.country_collection)
    if not r.city_collection in {"unknown", "Bergamo", "Durham", "Leeds", "Lisbon", "London", "Mississippi", "Nantes", "New York", "Plymouth", "Porto", "San Jose", "Tuscany", "Madison"}: raise ValueError(r.city_collection)

    return r

def get_nan_to_unknown(x):
    if pd.isna(x): return "unknown"
    else: return x



def sra_run_return_type_sample_clinical_or_enviromental_WGS(r):

    """Returns whether the run is clinical or environmental    """

    # define some predefined experiments
    if r.name in {"SRR24010542"}: return "clinical"

    # based on source
    if r.isolation_source in {'Blood Culture', 'Clinic', 'Host blood', 'Host nose', 'Host respiratory tract', 'Human gut', 'Peritoneal fluid', 'Urine sample', 'bile', 'blood', 'catheter-related bloodstream infection', 'clinical isolate', 'hospital', 'isolated from a patient with diarrhea in Puerto Rico by Ashford (1928)', 'urine'}: return 'clinical' # 'NRRL Collection', 'cobweb'

    elif r.isolation_source in {'soil sample'}: return 'environmental'

    elif r.isolation_source in {'Lab evolution experiment'}: raise ValueError("no lab samples allowed")

    # based on host
    if r.Host in {'"Homo sapiens"', 'Homo sapiens'}: return 'clinical'

    if r.Host_disease in ['Bloodstream infection', 'Candida infection', 'Candida parapsilosis infection', 'Candidiaisis or Colonization', 'Catheter-related bloodstream infection (CRBSI)', 'Chronic kidney disease', 'Kidney graft failure', 'Maligancy of soft tissue', 'Sepsis', 'Tropical sprue', 'candidal disease', 'severe diarrhea\\, sepsis']: return 'clinical'

    # based on sample_type
    if r.sample_type in {'clinical isolate'}: return 'clinical' # isolate
    elif r.sample_type in {'Lab evolution experiment'}: raise ValueError("no lab samples allowed")

    # based on different cases
    if r.env_local_scale in {'MSK hosptial', 'human-associated habitat'}: return 'clinical'
    if r.env_medium in {'faeces', 'human blood', 'human stool', 'orbital tissue'}: return 'clinical'
    if r.env_broad_scale in {'human gut mycobiome'}: return 'clinical'
    if r.Host_age in {'1.0', '79.0'}: return 'clinical'
    if r.host_disease_outcome in {"Recovery"}: return 'clinical'
    if r.host_disease_stage in {'Critical, sepsis', 'Severe'}: return 'clinical'
    if r.host_sex in {'male'}: return 'clinical'
    if r.environmental_history in {'clinical sample'}: return 'clinical'
    if r.host_associated in {'yes'}: return 'clinical'
    if r.host_scientific_name in {'Homo sapiens'}: return 'clinical'
    if r.Source in {'clinical'}: return 'clinical'

    # cfilter on BioSampleModel
    # if r.BioSampleModel in {'Microbe\\, viral or environmental'}: return 'unknown'

    # {'MIGS.eu,MIGS/MIMS/MIMARKS.human-associated', 'Microbe\\, viral or environmental', 'Pathogen.cl', 'Pathogen.env'}: return True


    # should not have been here
    for k in r.keys(): 
        if not pd.isna(r[k]): print(k, r[k])
    raise ValueError("Invalid r: %s"%r.name)

def get_sample_name_compatible_w_Bergin(r):

    """Gets some sample names renamed so that they are compatible with Bergin2022"""

    if r.BioProject=="PRJNA661583" and r["Sample Name"].startswith("Kw") and "/" in r["Sample Name"]: return r["Sample Name"].replace("/", "-")
    
    elif r.BioProject=="PRJNA748054" and all(map(lambda x: x.isdigit(), r["Sample Name"])): return "MSK%s"%(r["Sample Name"])

    elif r.BioProject=="PRJNA795920" and "-" in r["Sample Name"]: return r["Sample Name"].replace("-", "/") 

    else: return r["Sample Name"]

    # "Kw3259/15" : "Kw3259-15", "Kw2006/15":"Kw2006-15", "Kw1590/18":"Kw1590-18", "2159":"MSK2159"

def get_tree_with_cladeID_from_supportAndBranchLength(tree, min_fraction_branchLen, min_branch_suport):

    """Takes a tree and returns it with the cladeID"""

    tree = cp.deepcopy(tree)

    # change the support of the children of the root
    for root_child in tree.get_children(): root_child.support=100

    # define the minimum branch length as the min_fraction_branchLen*(maximum distance to the root)
    max_dist_btw_leaves = max(map(lambda l: l.get_farthest_node()[1], tree.get_leaves()))
    min_branchLen = max_dist_btw_leaves*min_fraction_branchLen

    # define branches that could be clades (potential_clade)
    for n in tree.traverse():

        if n.is_leaf(): n.potential_clade = False
        elif not n.is_leaf() and n.support>=min_branch_suport and n.dist>=min_branchLen: n.potential_clade = True
        else: n.potential_clade = False

    # define the clades as those branches that could be potential clades and have no internal branch that could be it as well.
    current_clade = 1
    type_traverse = "preorder"
    for n in tree.traverse(type_traverse):

        # if you find a clade, define all the leafs of the node to that clade
        if n.potential_clade is True and all([child.potential_clade is False for child in n.traverse() if child!=n]): 

            # map all the children to this clade
            n.cladeID = current_clade
            for child in n.traverse(): child.cladeID = current_clade
            current_clade+=1

    # add missing clade IDs 
    for n in tree.traverse(type_traverse):

        # if the clade ID has not been assigned, set to nan
        if "cladeID" not in set(n.__dict__): n.cladeID = -1

    return tree


def get_clade_info_one_min_fraction_branchLen(strains_tree, min_fraction_branchLen, min_support_tree):

    """For a given tree, with min_fraction_branchLen and min_support_tree, return a dict that has clade info"""

    # log
    #print("min_fraction_branchLen=%s"%min_fraction_branchLen)

    # get the tree with the clade IDs
    tree_withCladeIDs = get_tree_with_cladeID_from_supportAndBranchLength(strains_tree, min_fraction_branchLen, min_support_tree)
    leaf_to_clade = pd.Series({l.name : l.cladeID for l in tree_withCladeIDs.get_leaves()})

    fraction_samples_clade = sum(leaf_to_clade!=-1)/len(leaf_to_clade)
    nclades = len(set(leaf_to_clade[leaf_to_clade!=-1].values))

    return {"min_fraction_branchLen":min_fraction_branchLen, "fraction_samples_clade":fraction_samples_clade, "nclades":nclades}



def get_sorted_target_phenotypes_ML(df_ML_integrated_performance, threshold_ROC_AUC, df_ML_integrated_performance_all, fields_optimize_test, df_ML_metrics_top):

    """Get the phenotypes for which we can do some predictive MLs"""

    # check
    if len(df_ML_integrated_performance)!=len(df_ML_integrated_performance[["phenotype", "test_parmID"]].drop_duplicates()):
        raise ValueError("not valid")

    # get some fields
    df_ML_integrated_performance = df_ML_integrated_performance[["ROC_AUC_test_integrated_onlyCorrectModels", "mean_min_ROC_AUC_CVset_resamples", "mean_cross_splits_consistency_jaccard_samples", "phenotype", "min_n_train_parameters", "max_n_train_parameters", "correct_model", "mean_cross_splits_consistency_jaccard_genes", "test_parmID_exceptTestSS", "best_model", "redundant_samples_n_transitions", "mean_redundant_samples_ROC_AUC", "AUC_test"] + fields_optimize_test].copy()

    # get df_ML_integrated_performance_all parsed
    df_ML_integrated_performance_all = df_ML_integrated_performance_all.set_index(["phenotype", "test_sample_set", "test_parmID_exceptTestSS"]).copy()
    df_ML_metrics_top = df_ML_metrics_top.set_index(["phenotype", "test_sample_set", "test_parmID_exceptTestSS"]).copy()

    # define the stats df
    def get_df_ML_phenotype_stats_one_pheno(df_p):

        # init with data
        data_dict = {"phenotype":df_p.iloc[0].phenotype, 
                     "min_n_train_parameters": min(df_p[df_p.min_n_train_parameters!=-1].min_n_train_parameters),
                     "max_n_train_parameters": max(df_p.max_n_train_parameters),
                     "n_test_parameters": len(df_p)}

        # add
        for f in ["ROC_AUC_test_integrated_onlyCorrectModels", "mean_min_ROC_AUC_CVset_resamples", "mean_cross_splits_consistency_jaccard_samples"]:

            for fn_obj, fn_name in [(max, "max"), ((lambda x: np.percentile(x, 95)), "pct95")]:
                data_dict["%s_%s"%(fn_name, f)] = fn_obj(df_p[f].values)

        # add the correct models
        data_dict["n_correct_test_parameters"] = sum(df_p.correct_model)
        data_dict["fraction_correct_test_parameters"] = sum(df_p.correct_model) / len(df_p)

        # add info about correct models
        if data_dict["n_correct_test_parameters"]==0:
            data_dict["correct_models_max_AUC_test"] = 0.5
            #data_dict["correct_models_max_consistency_genes"] = 0.0
            data_dict["all_samples_AUC_test"] = 0.5
            data_dict["redundant_samples_AUC_test"] = -1
            data_dict["redundant_samples_n_transitions"] = 0
            data_dict["best_model_AUC_test"] = 0.5

            #data_dict["best_model_precision"] = 0.0
            #data_dict["best_model_recall"] = 0.0
            #data_dict["best_model_specificity"] = 0.0
            #data_dict["best_model_Fvalue"] = 0.0

            for f in fields_optimize_test:
                data_dict["best_model_%s"%f] = -1

        else:
            data_dict["correct_models_max_AUC_test"] = max(df_p[df_p.correct_model].ROC_AUC_test_integrated_onlyCorrectModels)
            #data_dict["correct_models_max_consistency_genes"] = max(df_p[df_p.correct_model].mean_cross_splits_consistency_jaccard_genes)

            # get data of best model
            best_model_info = df_p[df_p.best_model].copy()
            if len(best_model_info)!=1: raise ValueError("should be one best")
            best_model_info = best_model_info.iloc[0]
            data_dict["best_model_AUC_test"] = best_model_info.AUC_test

            # add the ROC AUC of all samples
            df_test_all_samples = df_ML_integrated_performance_all.loc[[(data_dict["phenotype"], "all_samples", best_model_info.test_parmID_exceptTestSS)]].copy()
            if len(df_test_all_samples)!=1: raise ValueError("should be 1")
            all_samples_info = df_test_all_samples.iloc[0]
            data_dict["all_samples_AUC_test"] = all_samples_info.AUC_test

            # get extra metrics for the best model
            #df_metrics_best = df_ML_metrics_top.loc[[(data_dict["phenotype"], "representative_samples", best_model_info.test_parmID_exceptTestSS)]].copy()
            #if len(df_metrics_best)!=len(df_metrics_best[["eval_data", "test_split_ID", ""]])
            #print(df_metrics_best)
            #dajhddjagdajhgda

            # add the AUC on the redundant samples
            data_dict["redundant_samples_n_transitions"] = int(best_model_info.redundant_samples_n_transitions)
            if data_dict["redundant_samples_n_transitions"]>0:
                data_dict["redundant_samples_AUC_test"] = best_model_info.mean_redundant_samples_ROC_AUC

            else:
                data_dict["redundant_samples_AUC_test"] = -1

            # add the parameters
            for f in fields_optimize_test:
                data_dict["best_model_%s"%f] = best_model_info[f]

        # check
        #if set(data_dict.keys()) != {'min_n_train_parameters', 'phenotype', 'all_samples_AUC_test', 'max_n_train_parameters', 'correct_models_max_AUC_test', 'fraction_correct_test_parameters', 'pct95_mean_min_ROC_AUC_CVset_resamples', 'n_test_parameters', 'pct95_ROC_AUC_test_integrated_onlyCorrectModels', 'max_mean_cross_splits_consistency_jaccard_samples', 'pct95_mean_cross_splits_consistency_jaccard_samples', 'redundant_samples_AUC_test', 'max_ROC_AUC_test_integrated_onlyCorrectModels', 'redundant_samples_n_transitions', 'n_correct_test_parameters', 'max_mean_min_ROC_AUC_CVset_resamples'}:
        #    raise ValueError("invalid keys")

        return data_dict

    df_ML_stats = pd.DataFrame(dict(df_ML_integrated_performance.groupby("phenotype").apply(get_df_ML_phenotype_stats_one_pheno).reset_index()[0])).transpose().set_index("phenotype", drop=False)#.reset_index(drop=True)
    check_no_nans_in_df(df_ML_stats)

    # sort
    df_ML_stats = df_ML_stats.sort_values(by=["fraction_correct_test_parameters", "pct95_ROC_AUC_test_integrated_onlyCorrectModels"], ascending=False)

    # define the correct phenos (threshold based on the observation that there is a sharp transition between these)
    print(df_ML_stats[["n_correct_test_parameters", "correct_models_max_AUC_test"]]) # fraction_correct_test_parameters
    #target_phenotypes_ML = list(df_ML_stats[(df_ML_stats.fraction_correct_test_parameters>=0.01) & (df_ML_stats.pct95_ROC_AUC_test_integrated_onlyCorrectModels>=0.6)].phenotype)
    target_phenotypes_ML = list(df_ML_stats[(df_ML_stats.fraction_correct_test_parameters>0)].phenotype)

    return target_phenotypes_ML, df_ML_stats

def get_df_ML_integrated_performance_filt_one_pheno_w_top_model_info(df_pheno, ntop_models):

    """adds the ML models"""

    # checks
    if len(df_pheno)!=len(set(df_pheno.test_parmID)):
        raise ValueError("invalid parms")

    # debug
    if sum(df_pheno.correct_model)==0:
        df_pheno["top_model"] = False
        df_pheno["best_model"] = False
        df_pheno["top_model_idx"] = -1

    else:

        # define df with top AUC
        max_AUC_test = max(df_pheno[(df_pheno.correct_model==True)].AUC_test)
        df_pheno_top = df_pheno[(df_pheno.correct_model==True) & (df_pheno.AUC_test>=(max_AUC_test-0.05))].copy()

        # sort
        sorting_fields = ["mean_cross_splits_consistency_jaccard_genes", "AUC_test", "AUC_train", "consistency_measure_idx", "AUC_diff_threshold_idx"]
        ascending_fields = [False, False, False, True, True]
        df_pheno_top = df_pheno_top.sort_values(by=sorting_fields, ascending=ascending_fields)

        # filter to get the top 10 consistency
        if len(df_pheno_top)>ntop_models:
            min_consistency = df_pheno_top.iloc[(ntop_models-1)].mean_cross_splits_consistency_jaccard_genes
            df_pheno_top = df_pheno_top[df_pheno_top.mean_cross_splits_consistency_jaccard_genes>=min_consistency].copy()

        # get the top models idx
        test_parmID_to_topModel_idx = dict(zip(df_pheno_top.test_parmID, range(len(df_pheno_top))))
        for pID in set(df_pheno.test_parmID).difference(set(df_pheno_top.test_parmID)):
            test_parmID_to_topModel_idx[pID] = -1

        # add to df
        top_parms = set(df_pheno_top.test_parmID)
        df_pheno["top_model"] = df_pheno.test_parmID.isin(top_parms)
        df_pheno["top_model_idx"] = df_pheno.test_parmID.map(test_parmID_to_topModel_idx)
        check_no_nans_series(df_pheno["top_model_idx"])

        # get the best model
        best_model_parm_df = get_df_like_first_row(df_pheno_top, sorting_fields, ascending_fields)
        if len(best_model_parm_df)!=1: 
            for f in best_model_parm_df.keys(): print(f, set(best_model_parm_df[f]))
            raise ValueError("two values same case")
        best_model_parm = best_model_parm_df.test_parmID.iloc[0]
        df_pheno["best_model"] = df_pheno.test_parmID==best_model_parm

    return df_pheno

def get_df_ML_results_top_with_single_threshold_metrics_best_t(df_ML_results_top, df_ML_metrics_top):

    """gets the df_ML_results_top with the metrics of the best model"""

    # keep
    df_ML_results_top = df_ML_results_top.copy()
    df_ML_metrics_top = df_ML_metrics_top.copy()

    # checks
    unique_model_fields = ["phenotype", "test_parmID_exceptTestSS", "test_split_ID"] 
    if len(df_ML_results_top)!=len(df_ML_results_top[unique_model_fields + ["test_sample_set"]].drop_duplicates()):
        raise ValueError("not unique combinations df_ML_results_top")

    if len(df_ML_metrics_top)!=len(df_ML_metrics_top[unique_model_fields + ["test_sample_set", "eval_data", "threshold"]].drop_duplicates()):
        raise ValueError("not unique combinations df_ML_metrics_top")

    # define the best threshold for each model, based on the cross splits data or train
    for target_eval_data in ["cross_splits", "test_integrated_onlyCorrectModels"]:

        # define best set
        #df_ML_metric_setBestT = df_ML_metrics_top[(df_ML_metrics_top.optimal_F1==True) & (df_ML_metrics_top.eval_data==target_eval_data)].copy()
        #df_ML_metric_setBestT = df_ML_metrics_top[(df_ML_metrics_top.optimal_balanced_accuracy==True) & (df_ML_metrics_top.eval_data==target_eval_data)].copy()
        #df_ML_metric_setBestT = df_ML_metrics_top[(df_ML_metrics_top.optimal_accuracy==True) & (df_ML_metrics_top.eval_data==target_eval_data)].copy()
        df_ML_metric_setBestT = df_ML_metrics_top[(df_ML_metrics_top.optimal_accuracy_custom==True) & (df_ML_metrics_top.eval_data==target_eval_data) & (df_ML_metrics_top.test_sample_set=="representative_samples")].copy()

        if len(df_ML_metric_setBestT)!=len(df_ML_metric_setBestT[unique_model_fields].drop_duplicates()):
            raise ValueError("not unique combinations df_ML_metric_setBestT")

        optT_field = "optimal_threshold_"+target_eval_data
        df_ML_metric_setBestT = df_ML_metric_setBestT[unique_model_fields + ["threshold"]].copy().rename(columns={"threshold":optT_field})

        # add to dfs
        df_ML_results_top = df_ML_results_top.merge(df_ML_metric_setBestT, on=unique_model_fields, how="left", validate="many_to_one")
        #df_ML_metrics_top = df_ML_metrics_top.merge(df_ML_metric_setBestT, on=unique_model_fields, how="left", validate="many_to_many")        
        check_no_nans_in_df(df_ML_results_top[unique_model_fields + [optT_field, "test_sample_set"]])
        #check_no_nans_in_df(df_ML_metrics_top[unique_model_fields + [optT_field]])

    # prepare df_ML_metrics_subset
    accuracy_fields = ['recall', 'precision', 'fpr', 'neg_predictive_val', 'F_score_precision_recall']
    fields_metrics = unique_model_fields + ["test_sample_set", "eval_data", "threshold"] + accuracy_fields
    index_fields = unique_model_fields + ["test_sample_set", "eval_data"]

    if len(df_ML_metrics_top)!=len(df_ML_metrics_top[index_fields + ["threshold"]].drop_duplicates()):
        raise ValueError("not unique combinations df_ML_metrics_top")

    df_ML_metrics_subset = df_ML_metrics_top[(df_ML_metrics_top.eval_data.isin({"cross_splits", "test_integrated_onlyCorrectModels", "redundant_samples"}))].copy()[fields_metrics].copy().set_index(index_fields)

    # add each of the datasets
    for eval_data in ["cross_splits", "redundant_samples", "test_integrated_onlyCorrectModels"]:
        for accuracy_f in accuracy_fields:

            # define the accuracy vals
            accuracy_df = df_ML_metrics_subset[[accuracy_f, "threshold"]].copy()

            # add to df
            optimalT_field = {"cross_splits":"optimal_threshold_cross_splits", "redundant_samples":"optimal_threshold_cross_splits", "test_integrated_onlyCorrectModels":"optimal_threshold_test_integrated_onlyCorrectModels"}[eval_data]

            def get_optimal_accuracy_one_r(r):
                idx_tuple = tuple(list(r[unique_model_fields]) + [r.test_sample_set, eval_data])
                accuracy_df_filt = accuracy_df.loc[idx_tuple].copy()
                accuracy_df_filt["diff_to_t"] = (accuracy_df_filt.threshold - r[optimalT_field]).apply(abs) 
                accuracy_series = accuracy_df_filt.sort_values(by=["diff_to_t", "threshold"], ascending=True).iloc[0]

                if accuracy_series.diff_to_t>=0.02: # picking the closest threshold, always 
                    raise ValueError("threshold should be almost the same, %s"%accuracy_series)

                return accuracy_series[accuracy_f]

            df_ML_results_top["%s_optimal_%s"%(eval_data, accuracy_f)] = df_ML_results_top[unique_model_fields + ["test_sample_set", optimalT_field]].apply(get_optimal_accuracy_one_r, axis=1)

    return df_ML_results_top


def get_set_representative_samples_target_samples(treefile_all, target_samples, outdir, df_phenotypes_all):

    """gets a representative set of samples ofor """

    # get phenotypes file
    make_folder(outdir)
    phenotypes_file = "%s/subset_phenotypes.tab"%outdir
    df_phenotypes = df_phenotypes_all[df_phenotypes_all.sampleID.isin(target_samples)].copy()
    if sorted(df_phenotypes.sampleID)!=sorted(target_samples): raise ValueError("all samples should be there")
    save_df_as_tab(df_phenotypes, phenotypes_file)

    # get tree
    treefile = "%s/tree_target_samples.nw"%outdir
    tree = Tree(treefile_all); tree.prune(target_samples, preserve_branch_length=True)
    tree.write(outfile=treefile, format=2)

    # get representative samples
    ancestral_GWAS_dir = ParentDir+"/scripts/ancestral_GWAS/scripts"
    outdir_rep_samples = "%s/getting_rep_samples"%outdir; make_folder(outdir_rep_samples)

    prunned_phenotypes_file = None
    while prunned_phenotypes_file is None:

        list_cmds = []
        prunned_phenotypes_file = get_balanced_samples_for_GWAS(phenotypes_file, treefile, ancestral_GWAS_dir, 70, "DOWNPASS", outdir_rep_samples, list_cmds, nrepresentative_samples=3, extra_args="--skip_prunning_uncertain_politomies --skip_plots")[0]

        if prunned_phenotypes_file is None: run_cmd(list_cmds[0][1], env="ancestral_GWAS_env")
        else: break

    return set(get_tab_as_df_or_empty_df(prunned_phenotypes_file).sampleID)

def get_affected_geneID_one_unique_groupID(unique_gID):

    """Takes one GWAS group and adds the affected geneIDs"""

    # split
    grouping_strategy, gID = unique_gID.split("|||")
    type_collapsing = grouping_strategy.split("|")[-1]


    # get geneID
    if type_collapsing=="genes":
        return gID

    elif type_collapsing=="domains":
        return gID.split("#")[1]

    elif type_collapsing=="variant_effect":
        return gID.split("|")[0]

    else:
        raise ValueError("invalid grouping %s"%unique_gID)



def get_df_ML_results_with_validation_accuracy_metrics(df_ML_results_all, CurDir, df_ML_metrics_all):

    """add the results of the validation"""
    print("adding validation data...")

    # define dirs
    DataDir = CurDir+"/data"
    ValidationDir = DataDir+"/validation_data"
    ValidationMetricsDir = ValidationDir+"/processed_data/validation_metrics"; make_folder(ValidationMetricsDir) 
    ancestral_GWAS_dir = ParentDir+"/scripts/ancestral_GWAS/scripts"

    # keep
    df_ML_results_all = df_ML_results_all.copy()
    df_ML_metrics_all = df_ML_metrics_all.copy()

    # iterate through phenotypes
    df_ML_results_all_wMetrics = pd.DataFrame()
    df_ML_metrics_all_wMetrics = pd.DataFrame()

    for (phenotype, aligner, sample_set, test_sample_set), df_ML_results in df_ML_results_all.copy().groupby(["phenotype", "aligner", "sample_set", "test_sample_set"]):
        print(phenotype, aligner, sample_set, test_sample_set)

        # get metrics df
        metrics_filters_dict = {"phenotype":phenotype, "aligner":aligner, "sample_set":sample_set, "test_sample_set":test_sample_set}
        df_ML_metrics = get_filtered_df_by_dict(df_ML_metrics_all, metrics_filters_dict)

        # checks
        if len(df_ML_results)!=len(df_ML_results[["test_parmID_exceptTestSS", "test_split_ID"]].drop_duplicates()):
            raise ValueError("invalid parms")

        if len(df_ML_metrics)!=len(df_ML_metrics[["test_parmID_exceptTestSS", "test_split_ID", "threshold", "eval_data"]].drop_duplicates()):
            raise ValueError("invalid df_ML_metrics")

        # get validation for certain phenotypes
        if phenotype in {"FLC_resistance", "VRC_resistance"}:

            # define an outdir for the files
            outdir_models = "%s/%s_%s_%s_%s"%(ValidationMetricsDir, phenotype, aligner, sample_set, test_sample_set)
            #delete_folder(outdir_models)
            make_folder(outdir_models)

            # define validation samples
            pheno_aln_dir = "%s/genotype_phenotype_relationships/pheno-%s_aligner-%s"%(ValidationDir, phenotype, aligner)
            gwas_inputs_dir = pheno_aln_dir + "/convergence_gwas/gwas_inputs_all_samples"
            validation_treefile_all = gwas_inputs_dir + "/corrected_rooted_tree.nw"
            validation_samples = set(Tree(validation_treefile_all).get_leaf_names())
            validation_samples_seqHere = {sID for sID in validation_samples if get_source_data(sID)=="sequenced_here"}
            validation_samples_SRA = {sID for sID in validation_samples if get_source_data(sID)=="SRA"}
            df_phenotypes_validation = get_tab_as_df_or_empty_df(pheno_aln_dir + "/convergence_gwas/phenotypes_all.tab")

            # define the training samples
            if sample_set=="all_samples": 
                train_samples = validation_samples_seqHere
            
            elif sample_set=="representative_samples": 
                train_samples = set(Tree("%s/data/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas/gwas_inputs_representative_samples/corrected_rooted_tree.nw"%(CurDir, phenotype, aligner)).get_leaf_names())

            else: raise ValueError("invalid ss")

            # define the testing samples
            if test_sample_set=="all_samples": 
                test_samples = validation_samples_SRA
            
            elif test_sample_set=="representative_samples":
                test_samples = get_set_representative_samples_target_samples(validation_treefile_all, validation_samples_SRA, outdir_models+"/getting_rep_samples_test", df_phenotypes_validation.copy())

            else: raise ValueError("invalid test ss") 

            # write the samples dataset
            train_test_samples = train_samples.union(test_samples)
            df_samples = df_phenotypes_validation[df_phenotypes_validation.sampleID.isin(train_test_samples)].copy()
            df_samples["type_data"] = df_samples.sampleID.isin(train_samples).map({True:"train", False:"test"})
            check_no_nans_in_df(df_samples)
            samples_file = "%s/sample_info.tab"%outdir_models
            save_df_as_tab(df_samples, samples_file)

            # checks
            if len(train_samples.difference(validation_samples_seqHere))>0: raise ValueError("some invalid train samples")
            if len(test_samples.difference(validation_samples_SRA))>0: raise ValueError("some invalid test samples")
            if len(train_samples.intersection(test_samples))>0: raise ValueError("cannot have same train / test")
            if set(test_samples)!=set(df_samples[df_samples.type_data=="test"].sampleID): raise ValueError("invalid test")

            # define treefile with subset samples
            tree = Tree(validation_treefile_all); tree.prune(train_test_samples, preserve_branch_length=True)
            validation_treefile = "%s/treefile_prunned.nw"%outdir_models
            tree.write(outfile=validation_treefile, format=2)

            # define the model's IDs
            models_fields = ["test_parmID_exceptTestSS", "test_split_ID", "selected_predictors"] # basic fields
            models_fields += ["model_name", "class_weight", "type_phenotype", "transition_definition", 'real_type_phenotype', 'aligner', 'sample_set', 'test_sample_set'] # parameter fields
            df_models = df_ML_results[~pd.isna(df_ML_results.selected_predictors)][models_fields].copy().sort_values(by=["test_parmID_exceptTestSS", "test_split_ID"])
            df_models["model_idx"] = list(range(len(df_models)))
            check_no_nans_in_df(df_models)
            models_file = "%s/model_info.tab"%outdir_models
            save_df_as_tab(df_models, models_file)

            # define all thresholds
            all_thresholds = sorted(set(get_filtered_df_by_dict(df_ML_metrics_all, {"phenotype":phenotype, "aligner":aligner}).threshold))
            thresholds_file = "%s/thresholds.txt"%outdir_models
            open(thresholds_file, "w").write("\n".join(list(map(lambda x: "%.3f"%x, all_thresholds))))

            # define the target genes (all features are single-gene related, as they are from the subset of genes)
            all_features_set = set.union(*df_models.selected_predictors.apply(lambda x: set.union(*[set(p.split(" * ")) for p in x.split(",")])))
            target_geneIDs = ",".join(sorted({get_affected_geneID_one_unique_groupID(f) for f in all_features_set}))

            # run the models' evaluation
            variants_file = gwas_inputs_dir+"/variants_df_all.tab"
            variant_annotations_file = gwas_inputs_dir+"/variant_annotations_df_all.tab"
            pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
            CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
            genome_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir # genome from CGD
            gff_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir
            interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir

            outdir_evaluation = outdir_models+"/models_eval"
            if file_is_empty(outdir_evaluation + "/pipeline_finished.txt"):

                print("running cmd")
                run_cmd(f"{ancestral_GWAS_dir}/evaluate_ML_models.py --outdir {outdir_evaluation} --variants_file {variants_file} --variant_annotations_file {variant_annotations_file} --samples_file {samples_file} --treefile {validation_treefile} --models_file {models_file} --threads 112 --gff {gff_Cparapsilosis} --ref_genome {genome_Cparapsilosis} --target_geneIDs {target_geneIDs} --pathway_annotations_dir {pathway_annotations_dir} --CGD_chromosomal_features {CGD_chromosomal_features} --interproscan_output {interproscan_output} --gDNA_code {Cpara_gDNA_code} --mtDNA_code {Cpara_mtDNA_code} --mtDNA_chrom {mtDNA_chrom} --default_pheno_probability 0.5 --thresholds_file {thresholds_file}"
                    , env="ancestral_GWAS_env")

            # add data about optimal metric
            #prefix_str = "validation_"+type_testing
            #optimal_metrics = df_metrics[df_metrics.threshold==r.optimal_threshold].iloc[0]
            #for f in ['recall', 'precision', 'fpr', 'neg_predictive_val', 'F_score_precision_recall']:
            #    r[prefix_str+"_optimalThreshold_"+f] = optimal_metrics.loc[f]

            # load dfs
            df_validation_results = get_tab_as_df_or_empty_df(outdir_evaluation+"/validation_results.tab")
            df_validation_metrics = get_tab_as_df_or_empty_df(outdir_evaluation+"/validation_per_threshold_metrics.tab")

            if len(df_validation_metrics[["test_parmID_exceptTestSS", "test_split_ID", "threshold", "type_testing"]].drop_duplicates())!=len(df_validation_metrics):
                raise ValueError("should be unique")

            # filter
            for df in [df_models, df_validation_results, df_validation_metrics]:
                df["tuple_parm_and_sID"] = df[["test_parmID_exceptTestSS", "test_split_ID"]].apply(tuple, axis=1)

            valid_tuples = set(df_models.tuple_parm_and_sID)
            df_validation_results = get_filtered_df_by_dict(df_validation_results, {"tuple_parm_and_sID":valid_tuples}, check_all_present=True)
            df_validation_metrics = get_filtered_df_by_dict(df_validation_metrics, {"tuple_parm_and_sID":valid_tuples}, check_all_present=True)

            # add the optimal threshold's metrics to df_validation_results
            parmID_sID_to_optThreshold = dict(df_ML_results.set_index(["test_parmID_exceptTestSS", "test_split_ID"]).optimal_threshold_cross_splits)
            for type_testing in ["transitions_TrainTest", "transitions_Test"]:

                # get the optimal dataset
                df_validation_metrics_optimalT = df_validation_metrics[(df_validation_metrics.apply(lambda r: r.threshold == parmID_sID_to_optThreshold[(r.test_parmID_exceptTestSS, r.test_split_ID)], axis=1)) & (df_validation_metrics.type_testing==type_testing)].copy()
                if len(df_validation_metrics_optimalT)!=len(df_ML_results): raise ValueError("should be same len")

                # add to df_validation_results
                accuracy_fields_l = ['recall', 'precision', 'fpr', 'neg_predictive_val', 'F_score_precision_recall']
                df_validation_metrics_optimalT = df_validation_metrics_optimalT[["test_parmID_exceptTestSS", "test_split_ID"]+accuracy_fields_l].copy().rename(columns={f : "validation_%s_optimalThreshold_%s"%(type_testing, f)  for f in accuracy_fields_l})
                
                df_validation_results = df_validation_results.merge(df_validation_metrics_optimalT, on=["test_parmID_exceptTestSS", "test_split_ID"] , how="left", validate="one_to_one")
                check_no_nans_in_df(df_validation_results)

            # add to df results
            res_fields = ["validation_%s_%s"%(type_testing, f) for type_testing in ["transitions_TrainTest", "transitions_Test"] for f in ["ROC_AUC", "optimalThreshold_recall", "optimalThreshold_precision", "optimalThreshold_fpr", "optimalThreshold_neg_predictive_val", "optimalThreshold_F_score_precision_recall", "PRC_AUC"]]
            df_validation_results = df_validation_results[["test_parmID_exceptTestSS", "test_split_ID"] + res_fields].copy()
            df_ML_results = df_ML_results.merge(df_validation_results, on=["test_parmID_exceptTestSS", "test_split_ID"] , how="left", validate="one_to_one")

            # keep all metrics
            df_validation_metrics = df_validation_metrics[["threshold", 'recall', 'precision', 'fpr', 'neg_predictive_val', 'F_score_precision_recall', 'har_mean_neg_pred_val_recall', 'test_parmID_exceptTestSS', 'test_split_ID', 'type_testing']].copy().rename(columns={"type_testing":"eval_data"})

            for k,v in metrics_filters_dict.items():
                df_validation_metrics[k] = v
            df_ML_metrics = pd.concat([df_ML_metrics, df_validation_metrics])

        # keeep
        df_ML_results_all_wMetrics = pd.concat([df_ML_results_all_wMetrics, df_ML_results]).reset_index(drop=True)
        df_ML_metrics_all_wMetrics = pd.concat([df_ML_metrics_all_wMetrics, df_ML_metrics]).reset_index(drop=True)

    return df_ML_results_all_wMetrics, df_ML_metrics_all_wMetrics


def get_accuracy_based_on_numbers(r, type_accuracy):

    """Calculates the accuracy based on TP, TN, FP, FN"""

    # get the lists
    y_true = [1] * int(r.TP) + [0] * int(r.FP) + [1] * int(r.FN) + [0] * int(r.TN)
    y_pred = [1] * int(r.TP) + [1] * int(r.FP) + [0] * int(r.FN) + [0] * int(r.TN)

    if len(y_true)==0:
        return 0.0

    if len(y_true)!=len(y_pred):
        raise ValueError("should be same")

    # get accuracy
    if type_accuracy=="normal":
        return sk_metrics.accuracy_score(y_true, y_pred)

    elif type_accuracy=="balanced":
        return sk_metrics.balanced_accuracy_score(y_true, y_pred)

    else:
        raise ValueError("invalid type_accuracy")

def get_df_ML_metrics_with_optimal_accuracy_custom(df_metrics, min_recall=0.0, max_fpr=1.0, threshold_lims=(0, 1), accuracy_f="recall"):

    """Gets the custom"""

    # check
    if len(df_metrics)!=len(set(df_metrics.threshold)): raise ValueError("should be unique t")

    # filter the df
    df_metrics_t = df_metrics[(df_metrics.recall>=min_recall) & (df_metrics.fpr<=max_fpr)].copy()
    df_metrics_t = df_metrics_t[(df_metrics_t.threshold>=threshold_lims[0]) & (df_metrics_t.threshold<=threshold_lims[1])].copy()

    # define the best threshold
    if len(df_metrics_t)>0:
        best_threshold = df_metrics_t.sort_values(by=[accuracy_f, "threshold"], ascending=[False, True]).threshold.iloc[0]

    else:
        #raise ValueError("cannot be o")
        best_threshold = -1

    # add to df
    df_metrics["optimal_accuracy_custom"] = df_metrics.threshold==best_threshold
    return df_metrics.copy()


def get_df_ML_metrics_with_optimal_accuracies_fields(df_metrics, min_recall=0.0, max_fpr=1.0):

    """gets df with optimal accuracies"""

    # check
    if len(df_metrics)!=len(set(df_metrics.threshold)): raise ValueError("should be unique t")

    # add best threshold each
    for accuracy_f in ["accuracy", "balanced_accuracy"]:
        df_metrics_t = df_metrics[(df_metrics.recall>=min_recall) & (df_metrics.fpr<=max_fpr)].copy()
        df_metrics_t = df_metrics_t[(df_metrics_t.threshold>=0.1) & (df_metrics_t.threshold<=0.9)].copy()

        if len(df_metrics_t)>0:
            best_threshold = df_metrics_t.sort_values(by=[accuracy_f, "threshold"], ascending=[False, True]).threshold.iloc[0]

        else:
            best_threshold = -1
        
        df_metrics["optimal_%s"%accuracy_f] = df_metrics.threshold==best_threshold

    return df_metrics.copy()

def get_top_models_and_metrics_only_highAUC_ones(df_ML_results_top, df_ML_metrics_top):

    """gets only the top models"""

    # filter the results
    def get_filt_df_results_one_set(df):

        # checks
        if len(df)!=len(df[["test_sample_set", "test_split_ID"]].drop_duplicates()):
            raise ValueError("should be unique")

        # define the top models
        df_rep_samples = df[df.test_sample_set=="representative_samples"].copy()
        max_AUC = max(df_rep_samples.cross_splits_ROC_AUC)
        df_rep_samples = df_rep_samples[(df_rep_samples.cross_splits_ROC_AUC >= (max_AUC-0.05))].copy()

        # get the best model
        #df_rep_samples = get_df_like_first_row(df_rep_samples, ["cross_splits_consistency_jaccard_genes", "cross_splits_ROC_AUC", "ROC_AUC_test", "test_split_ID"], [False, False, False, True])
        #if len(df_rep_samples)!=1: 
        #    raise ValueError("should be 1")
        
        # define the valid split
        valid_split_IDs = set(df_rep_samples.test_split_ID)

        # filter
        return df[df.test_split_ID.isin(valid_split_IDs)].copy()

    print("filtering top models")
    df_ML_results_top = df_ML_results_top[~pd.isna(df_ML_results_top.selected_predictors)].copy().groupby(["phenotype", "test_parmID_exceptTestSS"]).apply(get_filt_df_results_one_set).reset_index(drop=True)

    # filter the metrics
    merge_fields = ["phenotype", "test_parmID_exceptTestSS", "test_split_ID", "test_sample_set"]
    df_valid_parms = df_ML_results_top[merge_fields].copy()
    if len(df_valid_parms)!=len(df_valid_parms.drop_duplicates()): raise ValueError("should be unique")
    df_valid_parms["parm_highAUC_cross_splits"] = True

    if len(df_ML_metrics_top)!=len(df_ML_metrics_top[merge_fields + ["eval_data", "threshold"]].drop_duplicates()): raise ValueError("should be unique df_ML_metrics_top")
    df_ML_metrics_top = df_ML_metrics_top.merge(df_valid_parms, on=merge_fields, validate="many_to_one", how="left")
    df_ML_metrics_top = df_ML_metrics_top[df_ML_metrics_top.parm_highAUC_cross_splits==True].copy()

    return df_ML_results_top, df_ML_metrics_top

def get_unk_to_slash(x):
    if type(x)!=str: raise ValueError("should be str")
    if x in {"unknown", "-1"}: return "-"
    else: return x
            
def get_StrainInfo_r_with_added_metadata(r, df_metadata, df_metadata_binary_phenotypes, df_phenotypes_validation, SRA_runs_w_phenotypes, sorted_binary_phenos):
    
    """Get r with all added data"""
    
    # get the corresponding data
    if r.source_data=="sequenced_here":
        metadata_r = df_metadata.loc[r.sampleID]
        binary_r = df_metadata_binary_phenotypes.loc[r.sampleID]
        validation_r = None
        
    elif r.source_data=="SRA_initial":
        metadata_r = df_metadata.loc[r.sampleID]
        binary_r = None
        if r.sampleID in SRA_runs_w_phenotypes: validation_r = df_phenotypes_validation.loc[r.sampleID]
        else: validation_r = None
    
    elif r.source_data=="SRA_extra_validation":
        metadata_r = None
        binary_r = None
        validation_r = df_phenotypes_validation.loc[r.sampleID]
        
    else:
        raise ValueError("invalid source")
        
    # add isolate info
    if not validation_r is None:
        r["validation_strain_name"] = validation_r.Isolate
        r["validation_study"] = validation_r.paper_name
        
    else:
        for f in ["validation_strain_name", "validation_study"]: r[f] = "-"
        
    # add info from metadata_r
    if not metadata_r is None:
        r["country_collection"] = get_unk_to_slash(metadata_r.country_collection)
        r["city_collection"] = get_unk_to_slash(metadata_r.city_collection)
        
        if r.source_data=="sequenced_here": r["hospital"] = get_unk_to_slash(metadata_r.hospital_short)
        elif r.source_data=="SRA_initial":  r["hospital"]  = "-"
        else: raise ValueError("invalid source")
        
        r["patient_age"] = metadata_r.age
        r["patientID"] = metadata_r.duplicated_patientID
        r["year_collection"] = metadata_r.year_collection
        
        r["microsatellite_genotype"] = get_unk_to_slash(metadata_r.genotype)
        r["clade_Bergin2022"] = get_unk_to_slash(str(metadata_r.clade_Bergin2022))
        r["clade_basedOn_Bergin2022"] = get_unk_to_slash(str(metadata_r.cladeID_based_on_Bergin2022))    
        
    else:
        for f in ["country_collection", "city_collection", "hospital"]: r[f] = "-"        
        for f in ["patient_age", "patientID", "year_collection"]: r[f] = -1
        for f in ["microsatellite_genotype", "clade_Bergin2022", "clade_basedOn_Bergin2022"]: r[f] = "-"
        
        
    # add MICs
    sorted_drugs = ["FLC", "VRC", "POS", "IVZ", "ITR"] +  ["ANI", "CAS", "MIF"] + ["AMB", "5FC"]
    mic_fields = [f+"_MIC" for f in sorted_drugs]
    
    for mic_f in  mic_fields:
        
        if r.source_data=="sequenced_here": r[mic_f] = metadata_r[mic_f]
        elif r.sampleID in SRA_runs_w_phenotypes: r[mic_f] = validation_r[mic_f]
        else: r[mic_f] = -1.0
        
    # add disinfectant tolerance
    for disinf in ["chlorhexidine", "ethanol", "H202", "hypochlorite", "surfanio", "zwittergent 3-14"]:        
        f_disinf = f"pct viability {disinf}"
        
        if disinf=="hypochlorite": target_disinf = "HYPLOCHLORITE"
        else: target_disinf = disinf.upper()
    
        if r.source_data=="sequenced_here": r[f_disinf] = metadata_r['%s VIABILITY %s'%("%", target_disinf)]
        else: r[f_disinf] = -1.0
        
    # add etest phenotypes
    if r.source_data=="sequenced_here":
        r["Etest_colonies_in_halo"] = get_unk_to_slash(metadata_r['Colonies inside halo'])
        r["Etest_FLC_MIC_24h"] = metadata_r["MIC ETEST 24H"]
        r["Etest_FLC_MIC_48h"] = metadata_r["MIC ETEST 48H"]

    else:
        r["Etest_colonies_in_halo"] = "-"
        for f in ["Etest_FLC_MIC_24h", "Etest_FLC_MIC_48h"]: r[f] = -1
        
    # add binary phenotypes
    for p in sorted_binary_phenos:
        p_field = "binaryPheno "+get_paper_dicotomic_pheno(p)
        
        if r.source_data=="sequenced_here": 
            r[p_field]  = binary_r[p]
        
        else:            
            if p in {"FLC_resistance", "VRC_resistance"}:
                
                if r.sampleID in SRA_runs_w_phenotypes: 
                    mic_f_p = p.split("_")[0] + "_MIC"
                    if validation_r[mic_f_p]>0:
                        r[p_field] = {"R":1, "I":-1, "S":0}[get_dicothomic_resistance(validation_r[mic_f_p], mic_f_p)]
                    else:
                        r[p_field] = -1
                    
                else:
                    r[p_field] = -1
            else:
                r[p_field] = -1
            
            
    
    # get the r
    return r
        
        


    
    
def generate_StrainInfo_table(df_metadata, df_metadata_binary_phenotypes, PlotsDir, df_phenotypes_validation, df_reads_info_analysis_with_validation):
    
    """Generates table with all strains"""
    
    # filter the metadata df
    samples_discard_run1 = {"OSCAR01_L21-060", "OSCAR01_L21-067"}
    df_metadata = df_metadata[(df_metadata.species_run=="para") & (df_metadata.possibly_triploid==False) & (~df_metadata.runID.isin(samples_discard_run1))].copy()
    
    # change the indices
    df_metadata = df_metadata.set_index("sampleID")
    df_metadata_binary_phenotypes = df_metadata_binary_phenotypes.copy().set_index("sampleID")
    df_phenotypes_validation = df_phenotypes_validation.copy().set_index("Run")
    
    # add
    df_phenotypes_validation["IVZ_MIC"] = -1.0

    # init df with all samples
    df_strains = df_reads_info_analysis_with_validation[["runID", "source_data_type", "has_smiley"]].copy().rename(columns={"source_data_type":"source_data"})
    df_strains["sampleID"] = df_strains.apply(get_sampleID_from_runID, axis=1)
    df_strains = df_strains.drop(["runID"], axis=1).sort_values(by=["source_data", "sampleID"], ascending=[False, True]).reset_index(drop=True)
    if len(df_strains)!=len(set(df_strains.sampleID)): raise ValueError("should be unique")
        
    # add all the necessary extra rows
    print("adding extra rows")
    SRA_runs_w_phenotypes = set(df_phenotypes_validation.index)
    sorted_binary_phenos = [x for x in list(df_metadata_binary_phenotypes.keys()) if not x in {"runID", "sampleID"}]
    df_strains = df_strains.apply(get_StrainInfo_r_with_added_metadata, df_metadata=df_metadata, df_metadata_binary_phenotypes=df_metadata_binary_phenotypes, df_phenotypes_validation=df_phenotypes_validation, SRA_runs_w_phenotypes=SRA_runs_w_phenotypes, sorted_binary_phenos=sorted_binary_phenos, axis=1)
    check_no_nans_in_df(df_strains)
    
    # change the Lab. Catalunya
    df_strains["hospital"] = df_strains.hospital.apply(lambda x: x.replace('LABORATORIO REFERENCIA CATALUNYA (EL PRAT DE LLOBREGAT) (Prat de Llobregat)', 'MYCOLOGY REFERENCE CENTRE (BARCELONA) (Barcelona)'))
            
    # write   
    df_strains.to_excel(PlotsDir+"/StrainInfo.xlsx")
    
    return df_strains
    
    
  

def get_processed_ML_classifier_dfs(df_ML_integrated_performance, threshold_ROC_AUC, threshold_consistency_genes, fields_optimize_test, ProcessedDataDir, DataDir, CurDir):

    """gets the """

    # keep
    df_ML_integrated_performance = df_ML_integrated_performance.copy()

    # change some fields
    for f, min_v in [("ROC_AUC_test_integrated_onlyCorrectModels", 0.5), 
                     ("mean_min_ROC_AUC_CVset_resamples", 0.5), 
                     ("mean_cross_splits_consistency_jaccard_samples", 0.0)]:

        df_ML_integrated_performance[f] = df_ML_integrated_performance[f].apply(lambda x: max([x, min_v]))

    # print vals
    for f in fields_optimize_test:
        print(f, sorted(set(df_ML_integrated_performance[f])))

    # add fields
    df_ML_integrated_performance["AUC_train"] = df_ML_integrated_performance.mean_min_ROC_AUC_CVset_resamples.copy()
    df_ML_integrated_performance["AUC_test"] = df_ML_integrated_performance.ROC_AUC_test_integrated_onlyCorrectModels.copy()
    df_ML_integrated_performance["AUC_train_test_diff"] = (df_ML_integrated_performance.AUC_train - df_ML_integrated_performance.AUC_test)
    df_ML_integrated_performance["type_gene_set"] = df_ML_integrated_performance.gene_set.apply(get_type_gene_set_ML)
    df_ML_integrated_performance["test_parmID_exceptTestSS"] = df_ML_integrated_performance[fields_optimize_test].apply(lambda r: "---".join(list(map(str, r))), axis=1)

    # checks
    if len(df_ML_integrated_performance)!=len(df_ML_integrated_performance[["test_parmID_exceptTestSS", "test_sample_set", "phenotype"]].drop_duplicates()):
        raise ValueError("combinations not correct")

    # add fields of idxs
    for f, sorted_vals in [("consistency_measure", ['jaccard_genes', 'jaccard_variants', 'jaccard_samples', 'none']),
                           ("AUC_diff_threshold", [0.05, 0.1, 0.15])]:

        val_to_idx = dict(zip(sorted_vals, range(len(sorted_vals))))
        df_ML_integrated_performance["%s_idx"%f] = df_ML_integrated_performance[f].map(val_to_idx)
        check_no_nans_series(df_ML_integrated_performance["%s_idx"%f])

    # get the interesting dfs
    filters_dict = {"test_sample_set":{"representative_samples"}, "multiple_test_corr":{"none"}}
    df_ML_integrated_performance_filt = get_filtered_df_by_dict(df_ML_integrated_performance, filters_dict, check_all_present=True).copy()

    # add the definition of correct model
    df_ML_integrated_performance_filt["correct_model"] =  (df_ML_integrated_performance_filt.AUC_test>=threshold_ROC_AUC) & (df_ML_integrated_performance_filt.AUC_train>=threshold_ROC_AUC) & ((df_ML_integrated_performance_filt.AUC_train_test_diff.apply(abs)<=0.2)) & (df_ML_integrated_performance_filt.mean_cross_splits_consistency_jaccard_genes>=threshold_consistency_genes)

    # add whether it is a top model
    df_ML_integrated_performance_filt = df_ML_integrated_performance_filt.groupby("phenotype").apply(get_df_ML_integrated_performance_filt_one_pheno_w_top_model_info, ntop_models=10).reset_index(drop=True)
    check_no_nans_in_df(df_ML_integrated_performance_filt, ignore_keys={"types_phenotype"})

    # load the results and metrics of top models
    df_ML_results_top_file = "%s/df_ML_results_top_models_tAUC=%s_tCons=%s.py"%(ProcessedDataDir, threshold_ROC_AUC, threshold_consistency_genes)
    df_ML_metrics_top_file = "%s/df_ML_metrics_top_models_tAUC=%s_tCons=%s.py"%(ProcessedDataDir, threshold_ROC_AUC, threshold_consistency_genes)

    if file_is_empty(df_ML_results_top_file) or file_is_empty(df_ML_metrics_top_file):
        print("Getting top models")

        # init
        df_ML_results_top = pd.DataFrame()
        df_ML_metrics_top = pd.DataFrame()

        # each phen
        for (phenotype, aligner, gene_set), df_performance in df_ML_integrated_performance_filt[df_ML_integrated_performance_filt.top_model==True].copy().groupby(["phenotype", "aligner", "gene_set"]):
            print(phenotype, aligner, gene_set)

            # load dfs
            classifier_dir = "%s/ML_classifiers_integrated_results/pheno-%s/gene_set=%s_aligner=%s/classifier_results"%(DataDir, phenotype, gene_set, aligner)
            df_ML_results = get_tab_as_df_or_empty_df(classifier_dir + "/best_models_test_performance.tab")
            df_ML_metrics = get_tab_as_df_or_empty_df(classifier_dir + "/best_models_test_threshold_metrics.tab")

            # add ID to filter
            print("adding IDs")
            for df in [df_ML_results, df_ML_metrics]:
                df["test_parmID_exceptTestSS"] = df[fields_optimize_test].apply(lambda r: "---".join(list(map(str, r))), axis=1)
                df["phenotype"] = phenotype
                df["aligner"] = aligner
                df["gene_set"] = gene_set

            # filter 
            print("filtering")
            top_test_parmID_exceptTestSS = set(df_performance.test_parmID_exceptTestSS)
            df_ML_results = df_ML_results[df_ML_results.test_parmID_exceptTestSS.isin(top_test_parmID_exceptTestSS)].copy()
            df_ML_metrics = df_ML_metrics[df_ML_metrics.test_parmID_exceptTestSS.isin(top_test_parmID_exceptTestSS)].copy()

            # check that lens make sense
            if len(df_ML_results) != (len(df_performance)*5*2):
                raise ValueError("invalid len df_ML_results")

            unique_fields_metrics = ["test_parmID_exceptTestSS", "test_sample_set", "threshold", "eval_data", "test_split_ID"]
            if len(df_ML_metrics) != len(df_ML_metrics[unique_fields_metrics].drop_duplicates()):
                raise ValueError("invalid len df_ML_metrics")

            # add extra fields to the metrics df
            print("adding accuracy")
            df_ML_metrics["accuracy"] = df_ML_metrics[["TP", "FP", "FN", "TN"]].apply(get_accuracy_based_on_numbers, type_accuracy="normal", axis=1)
            df_ML_metrics["balanced_accuracy"] = df_ML_metrics[["TP", "FP", "FN", "TN"]].apply(get_accuracy_based_on_numbers, type_accuracy="balanced", axis=1)
            check_no_nans_in_df(df_ML_metrics[["accuracy", "balanced_accuracy"]])
            len_before = len(df_ML_metrics)
            df_ML_metrics = df_ML_metrics.groupby(["test_parmID_exceptTestSS", "test_sample_set", "eval_data", "test_split_ID"]).apply(get_df_ML_metrics_with_optimal_accuracies_fields).reset_index(drop=True)
            if len_before!=len(df_ML_metrics): raise ValueError("invalid lens")

            # add
            df_ML_results_top = pd.concat([df_ML_results_top, df_ML_results]).reset_index(drop=True)
            df_ML_metrics_top = pd.concat([df_ML_metrics_top, df_ML_metrics]).reset_index(drop=True)

        # save
        print("saving...")
        save_object(df_ML_results_top, df_ML_results_top_file)
        save_object(df_ML_metrics_top, df_ML_metrics_top_file)

    df_ML_results_top = load_object(df_ML_results_top_file)
    df_ML_metrics_top = load_object(df_ML_metrics_top_file)

    # filter some
    print("filter...")
    target_phenos_manual = {'age_infant', 'FLC_resistance', 'VRC_resistance', 'source_blood_vs_other_sterile'}
    df_ML_results_top = df_ML_results_top[df_ML_results_top.phenotype.isin(target_phenos_manual)].copy()
    df_ML_metrics_top = df_ML_metrics_top[df_ML_metrics_top.phenotype.isin(target_phenos_manual)].copy()

    # add fields
    df_ML_metrics_top["specificity"] = 1 - df_ML_metrics_top.fpr

    # keep only top models that have high AUC
    df_ML_results_top, df_ML_metrics_top = get_top_models_and_metrics_only_highAUC_ones(df_ML_results_top, df_ML_metrics_top)

    # define the balanced accuracy
    print("adding accuracy")
    #custom_filers = {"min_recall":0.0, "max_fpr":1.0, "threshold_lims":(0.1, 0.9), "accuracy_f":"balanced_accuracy"}
    custom_filers = {"min_recall":0.0, "max_fpr":0.3, "threshold_lims":(0.05, 0.95), "accuracy_f":"recall"}
    #custom_filers = {"min_recall":0.0, "max_fpr":0.3, "threshold_lims":(0.1, 0.9), "accuracy_f":"balanced_accuracy"}


    df_ML_metrics_top = df_ML_metrics_top.groupby(["phenotype", "test_parmID_exceptTestSS", "test_sample_set", "eval_data", "test_split_ID"]).apply(get_df_ML_metrics_with_optimal_accuracy_custom, **custom_filers).reset_index(drop=True)

    # add metrics of the single-threshold data
    df_ML_results_top = get_df_ML_results_top_with_single_threshold_metrics_best_t(df_ML_results_top, df_ML_metrics_top)

    # add to df_ML_results_top the evaluation on independent validation datasets
    df_ML_results_top, df_ML_metrics_top = get_df_ML_results_with_validation_accuracy_metrics(df_ML_results_top, CurDir, df_ML_metrics_top)

    # add info
    pheno_to_best_modelID = dict(df_ML_integrated_performance_filt[df_ML_integrated_performance_filt.best_model==True][["phenotype", "test_parmID_exceptTestSS"]].drop_duplicates().set_index("phenotype").test_parmID_exceptTestSS)
    for df in [df_ML_results_top, df_ML_metrics_top]:
        df["best_model"] = df[["phenotype", "test_parmID_exceptTestSS"]].apply(lambda r: r.test_parmID_exceptTestSS == pheno_to_best_modelID[r.phenotype], axis=1)

    # define the phenotypes for which I can get reliable ML models
    target_phenotypes_ML, df_ML_stats = get_sorted_target_phenotypes_ML(df_ML_integrated_performance_filt, threshold_ROC_AUC, df_ML_integrated_performance, fields_optimize_test, df_ML_metrics_top)

    # get dict
    return {"target_phenotypes_ML":target_phenotypes_ML, "df_ML_stats":df_ML_stats, "df_ML_integrated_performance_filt":df_ML_integrated_performance_filt, "df_ML_integrated_performance":df_ML_integrated_performance, "df_ML_results_top":df_ML_results_top, "df_ML_metrics_top":df_ML_metrics_top}



def get_selected_predictors_demult_one_tuple_predictors(x):

    """Gets one set of predictors"""

    if len(x)==0:
        return set()

    else:
        return set.union(*[set(p.split(" * ")) for p in x])

def get_type_genes_ML_GWAS_overlaps(gID, genes_ML, genes_GWAS_high, genes_GWAS_low, geneIDs_manual):

    """Gets different gene types"""

    if gID in genes_ML and gID in genes_GWAS_high:
        return "ML & GWAS high-confidence"

    elif gID in genes_ML and gID in genes_GWAS_low:
        return "ML & GWAS low-confidence"

    elif gID in genes_ML:
        return "ML only"

    elif gID in genes_GWAS_high:
        return "GWAS high-confidence"

    elif gID in genes_GWAS_low:
        return "GWAS low-confidence"

    elif gID in geneIDs_manual:
        return "manual curation only"

    elif not any([gID in x for x in [genes_ML, genes_GWAS_high, genes_GWAS_low, geneIDs_manual]]):
        return "none"

    else: raise ValueError("invalid gene %s, %s, %s, %s"%(gID, genes_ML, genes_GWAS_high, genes_GWAS_low))


def get_inclusion_criteria_annot(gID, dict_target_genes, criteria_to_annot):

    """Gets the inclusion criteria of the genes"""

    all_inclusion_criteria = {criteria for criteria, geneIDs in dict_target_genes.items() if gID in geneIDs}

    if len(all_inclusion_criteria)>1: raise ValueError("invalid, %s, %s"%(all_inclusion_criteria, all_inclusion_criteria))
    if len(all_inclusion_criteria)==1: return criteria_to_annot[next(iter(all_inclusion_criteria))]
    if len(all_inclusion_criteria)==0: return ""

def calculate_epsilon(nodes_GenoAndPheno, nodes_withGeno, nodes_withPheno):

    """Calculates epsilon"""

    numerator_epsilon = 2*nodes_GenoAndPheno
    denominator_epsilon = nodes_withGeno + nodes_withPheno
    if denominator_epsilon==0: epsilon = 0.0
    else: epsilon = numerator_epsilon / denominator_epsilon

    return epsilon

def get_df_gwas_one_pheno_with_per_gene_stats(df_gwas):

    """Adds fields to df"""

    # checks
    for f in ["type_filtering", "phenotype", "nodes_withPheno"]:
        if len(set(df_gwas[f]))!=1: raise ValueError("should be 1 %s"%f)

    if any(df_gwas.Gene=="-"): raise ValueError("not all genes as expected")

    # map each Gene to the fraction of transitions attributed to it
    gene_to_nnodes = df_gwas.groupby("Gene").apply(lambda df_g: set.union(*df_g.nodes_transitions.apply(lambda x: set(x.split(","))))).apply(len)
    gene_to_fraction_nodes = gene_to_nnodes / df_gwas.nodes_withPheno.iloc[0]
    if any(gene_to_fraction_nodes>1): raise ValueError("cannot be 1")
    df_gwas["fraction_transitions_by_Gene"] = df_gwas.Gene.map(dict(gene_to_fraction_nodes))

    # add other stats
    for f in ["epsilon", "OR", "nodes_GenoAndPheno"]:
        gene_to_max_f =  df_gwas.groupby("Gene").apply(lambda df_g: max(df_g[f]))
        df_gwas["max_"+f] = df_gwas.Gene.map(dict(gene_to_max_f))
    
    # checks
    check_no_nans_in_df(df_gwas[["fraction_transitions_by_Gene", "max_epsilon", "max_OR", "max_nodes_GenoAndPheno"]])

    return df_gwas

def get_annot_with_fraction_transitions_by_Gene(r, phenotype, gene_to_fraction_nodes_explained_low, gene_to_fraction_nodes_explained_high):

    """Adds the fractoin nodes explained"""

    # only change those that are empty
    if r[phenotype]=="": 

        if r.geneID in gene_to_fraction_nodes_explained_high: 
            fraction_nodes = round(gene_to_fraction_nodes_explained_high[r.geneID], 1)

        elif r.geneID in gene_to_fraction_nodes_explained_low: 
            fraction_nodes = round(gene_to_fraction_nodes_explained_low[r.geneID], 1)

        else: return ""

        if fraction_nodes==1: return "1"
        else: return "." + str(fraction_nodes).split(".")[1]

    else:
        return r[phenotype]

def get_df_var_info_NR_certain_variant_groups(dict_data, outdir_NR_groups):

    """Gets a df with variant info for certain NR groups, as provided by dict_data"""

    # define general dirs
    ancestral_GWAS_dir = "%s/scripts/ancestral_GWAS/scripts"%ParentDir
    DataDir = ParentDir + "/Cparapsilosis_popgen/data"
    pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
    interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir
    CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
    genome_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir # genome from CGD
    gff_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir
    mapping_genes = "Scer_HAL9:TAC1" 
    reactome_organism_priority = "Saccharomyces_cerevisiae,Schizosaccharomyces_pombe"

    # get the NR 
    make_folder(outdir_NR_groups)
    if file_is_empty(outdir_NR_groups + "/pipeline_finished.txt"):
        print(dict_data)

        # checks
        if len(dict_data["sorted_groups"])==0: raise ValueError("should have some sorted groups")

        # check that the groups exist
        gwas_inputs_dir = "%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas/gwas_inputs_%s"%(DataDir, dict_data["phenotype"], dict_data["aligner"], dict_data["sample_set"])
        df_groups = load_bash_cmd_output_into_df("cut -f2,7 %s/all_groupings.tab"%gwas_inputs_dir)
        df_groups["unique_groupID"] = df_groups.grouping_ID + "|||" + df_groups.group

        strange_groups = set(dict_data["sorted_groups"]).difference(set(df_groups.unique_groupID))
        if len(strange_groups)>0: raise ValueError("invalid %s, sorted_groups: %s"%(strange_groups, dict_data["sorted_groups"]))

        # define a file with the input groups
        inputs_groups_file = outdir_NR_groups + "_input_groups.tab"
        df_inputs = pd.DataFrame({"unique_groupID":dict_data["sorted_groups"], "sorting_f":0})

        df_inputs["type_vars"] = df_inputs.unique_groupID.apply(lambda x: x.split("|")[0])
        df_inputs["type_genes"] = df_inputs.unique_groupID.apply(lambda x: x.split("|")[1])
        df_inputs["type_mutations"] = df_inputs.unique_groupID.apply(lambda x: x.split("|")[2])
        df_inputs["type_collapsing"] = df_inputs.unique_groupID.apply(lambda x: x.split("|")[3])

        save_df_as_tab(df_inputs, inputs_groups_file)

        # run pipeline
        cmd = f'{ancestral_GWAS_dir}/generate_NR_groups.py --outdir {outdir_NR_groups} --gwas_hits {inputs_groups_file} --gff {gff_Cparapsilosis} --ref_genome {genome_Cparapsilosis} --custom_geneID_maps {mapping_genes} --gwas_inputs_dir {gwas_inputs_dir} --pathway_annotations_dir {pathway_annotations_dir} --CGD_chromosomal_features {CGD_chromosomal_features} --interproscan_output {interproscan_output} --initial_sorting_fields sorting_f --initial_sorting_ascending False --reactome_organism_priority {reactome_organism_priority} --min_support {dict_data["min_support"]}  --ASR_method {dict_data["ASR_method"]} {dict_data["extra_args"]}'

        run_cmd(cmd, env="ancestral_GWAS_env")

    df_groups = get_tab_as_df_or_empty_df(outdir_NR_groups + "/NR_hits.tab")

    # add fields
    gwas_results_tab = "%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas/gwas_results_%s_min_sup_%i_ASR_%s/integrated_GWAS_stats.tab"%(DataDir, dict_data["phenotype"], dict_data["aligner"], dict_data["sample_set"], dict_data["min_support"], dict_data["ASR_method"])
    df_gwas_head = load_bash_cmd_output_into_df("head -n500 %s"%gwas_results_tab)
    df_gwas_head = df_gwas_head[df_gwas_head.gwas_method=="synchronous"].copy()

    if len(set(df_gwas_head.nodes_withPheno))!=1: raise ValueError("There should be the same total # transitions")
    df_groups["nodes_withPheno"] = df_gwas_head.nodes_withPheno.iloc[0]

    return df_groups

def get_command_stdout(cmd, env="Cparapsilosis_popgen_env", skip_rstrip=False):
    CondaDir = "/gpfs/projects/bsc40/mschikora/anaconda3"
    SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
    cmd_prefix = "%s && conda activate %s &&"%(SOURCE_CONDA_CMD, env)
    cmd_to_run = "%s %s"%(cmd_prefix, cmd)

    str_output = str(subprocess.check_output(cmd_to_run, shell=True)).split("'")[1]

    if skip_rstrip is True:
        return str_output

    else:
        return str_output.rstrip("\\n")

def load_bash_cmd_output_into_df(cmd, write_file=None, env="Cparapsilosis_popgen_env"):

    """loads bash cmd output into df"""

    # parse on the fly
    if write_file is None:

        # get text
        text_df = "\n".join(map(lambda l: l.replace("\\t", "\t"), get_command_stdout(cmd, env=env).split("\\n")))

        # get df
        return pd.read_csv(StringIO(text_df), sep="\t")

    else:

        # write and read
        run_cmd(cmd + " > %s"%write_file, env=env)
        return pd.read_csv(write_file, sep="\t")

def generate_GWAS_and_ML_integrated_table(table_file, df_plot, df_GWAS_hits_NR, df_ML_results_top, df_genes, df_annot, paper_pheno_to_pheno, criteria_to_annot, criteria_to_label_criteria):

    """Gets the table integrated"""

    print("Getting integrated table")

    if not file_is_empty(table_file):
        print("table already generated")
        return pd.read_excel(table_file)

    # keep dfs
    df_GWAS_hits_NR = df_GWAS_hits_NR.copy()
    df_ML_results_top = df_ML_results_top.copy()
    df_genes = df_genes.set_index("geneID").copy()
    df_annot = df_annot.set_index("geneID").copy()
    df_plot = df_plot.set_index("geneID").copy()

    # define general dirs
    DataDir = ParentDir + "/Cparapsilosis_popgen/data"
    ProcessedDataDir = ParentDir + "/Cparapsilosis_popgen/processed_data"

    # define different dfs
    df_GWAS_highc = df_GWAS_hits_NR[df_GWAS_hits_NR.type_filtering=="high_confidence"].copy()
    df_GWAS_lowc = df_GWAS_hits_NR[df_GWAS_hits_NR.type_filtering=="low_confidence"].copy()

    # map each geneID to the gene name
    gene_features_df = get_gene_features_df(DataDir).set_index("geneID", drop=False)
    geneID_to_name = dict(gene_features_df.final_name)
    tac1_ID = gene_features_df[gene_features_df.final_name=="Scer_HAL9"].iloc[0].geneID
    geneID_to_name[tac1_ID] = "TAC1"
    geneName_to_ID = {v:k for k,v in geneID_to_name.items()}

    # map each phenotype to the prioritized genes
    priotization_annot = set(criteria_to_annot.values())
    paper_phenotype_to_priorized_genes = df_annot.apply(lambda c: set(c[c.isin(priotization_annot)].index), axis=0)

    # define genes that are priorized in ergosterol-related phenotypes (azoles R, AMB R)
    ergosterol_phenotypes = {"fluconazole R", "voriconazole R", "posaconazole R", "amphotericin B R"}
    genes_ergosterol = set.union(*paper_phenotype_to_priorized_genes[list(ergosterol_phenotypes)])
    print("genes_ergosterol", {geneID_to_name[x] for x in genes_ergosterol})

    # define dicts
    annot_to_criteria = {v:k for k,v in criteria_to_annot.items()}

    """
    flz_only = []
    flz_and_vrc = []
    panazole_genes = ["ERG11", "MRR1", "TAC1", "NDT80", "Scer_TEL1", "Scer_NCP1", "CPAR2_806400"] 
    """

    # check that geneIDs are as expected
    sorted_geneIDs = sorted(set(df_plot.index))
    #sorted_geneIDs = [geneName_to_ID[x] for x in sorted_gene_names]
    #if len(set(sorted_geneIDs))!=len(sorted_geneIDs): raise ValueError("should be unique")
    #if set(sorted_geneIDs)!=set(df_plot.index): raise ValueError("some missings")

    # keep adding annotations for each combination of gene and phenotype
    integrated_table_rows = []
    for paper_phenotype in df_plot.columns:
        for geneID in sorted_geneIDs:
            gene_name = geneID_to_name[geneID]

            # get the type of hit
            type_hit = df_plot.loc[geneID, paper_phenotype]  

            # discard some cases
            if type_hit=="none": continue

            # keep only certain genes
            if not ((geneID in paper_phenotype_to_priorized_genes[paper_phenotype]) or (geneID in genes_ergosterol and paper_phenotype in ergosterol_phenotypes)): continue

            # checks
            annot_val = df_annot.loc[geneID, paper_phenotype]
            if annot_val=="": raise ValueError("annot val cannot be empty")
            if not annot_val in priotization_annot:
                if not (geneID in genes_ergosterol and paper_phenotype in ergosterol_phenotypes): raise ValueError("invalid annot")

            # define the inclusion criteria
            if annot_val in priotization_annot:
                inclusion_criteria = criteria_to_label_criteria[annot_to_criteria[annot_val]]

            elif (geneID in genes_ergosterol and paper_phenotype in ergosterol_phenotypes):
                inclusion_criteria = "ergosterol phenotypes"

            else: raise ValueError("invalid")

            # depending on the type of hit, get different GWAS NR-like dfs
            if type_hit in {'GWAS low-confidence'}:
                df_vars = df_GWAS_lowc[(df_GWAS_lowc.phenotype==paper_phenotype) & (df_GWAS_lowc.Gene==geneID)].copy()

            elif type_hit in {'GWAS high-confidence', 'ML & GWAS high-confidence'}:
                df_vars = df_GWAS_highc[(df_GWAS_highc.phenotype==paper_phenotype) & (df_GWAS_highc.Gene==geneID)].copy()

            elif type_hit in {'ML only', 'ML & GWAS low-confidence', 'manual curation only'}:

                # define a dict with the data
                if type_hit=='manual curation only':
                    dict_data = {"aligner":"bwa_mem", "sorted_groups":["all_vars|all_genes|non_syn_muts|genes|||"+geneID], "sample_set":"representative_samples", "ASR_method":"MPPA,DOWNPASS", "min_support":70, "phenotype":paper_pheno_to_pheno[paper_phenotype], "extra_args":""}

                elif type_hit in {'ML only', 'ML & GWAS low-confidence'}:

                    # filter df with the gene and pheno
                    df_ML = df_ML_results_top[(df_ML_results_top.phenotype==paper_phenotype) & (df_ML_results_top.all_affected_genes_set.apply(lambda x: geneID in x))].copy()

                    # keep some better cases
                    sorting_fields_idx = []
                    for f, vals in [("aligner", ["hisat2_no_spliced", "bowtie2_local", "bwa_mem"]), ("sample_set", ["representative_samples", "all_samples"])]:

                        val_to_I = dict(zip(vals, range(len(vals))))
                        df_ML[f+"_idx"] = df_ML[f].map(val_to_I); check_no_nans_series(df_ML[f+"_idx"])
                        sorting_fields_idx.append(f+"_idx")

                    df_ML = get_df_like_first_row(df_ML, sorting_fields_idx, [True, True])
                    for f in ["aligner", "sample_set", "ASR_method", "min_support"]:
                        if len(set(df_ML[f]))!=1: raise ValueError("invalid %s"%(df_ML[["aligner", "sample_set", "ASR_method", "min_support"]]))

                    # add to df
                    dict_data = {"sorted_groups": sorted({g for g in set.union(*df_ML.selected_predictors.apply(set)) if get_affected_geneID_one_unique_groupID(g)==geneID}), "phenotype":paper_pheno_to_pheno[paper_phenotype], "extra_args":""}
                    for f in ["aligner", "ASR_method", "min_support", "sample_set"]:
                        dict_data[f] = df_ML.iloc[0][f]

                    # change the sample set for a few cases
                    if (paper_phenotype, gene_name)==("blood vs other sterile", "RBT1"): 
                        dict_data["sample_set"] = "all_samples" # representative_samples missing some features
                    
                    if (paper_phenotype, gene_name)==("fluconazole R", "MRR1"): 
                        dict_data["sorted_groups"] = sorted(set(dict_data["sorted_groups"]).difference({'small_vars|only_protein_coding|truncating_muts|domains|||pChunk_291_300#CPAR2_807270#291_300'})) # missing in some, but redundant

                else: raise ValueError("invalid")

                # get the vars
                redundancy_reduction_dir = ProcessedDataDir + "/redundancy_reduction_integrated_tables"; make_folder(redundancy_reduction_dir)
                outdir_NR = "%s/%s_%s"%(redundancy_reduction_dir, dict_data["phenotype"], geneID)
                df_vars = get_df_var_info_NR_certain_variant_groups(dict_data, outdir_NR)

                # filter
                if any(df_vars.n_transitions_var>=1):
                    df_vars = df_vars[df_vars.n_transitions_var>=1].copy()

                else:
                    print("WARNING: For %s-%s there are no attributable transitions"%(paper_phenotype, gene_name))

            else: raise ValueError("invalid %s"%type_hit)

            # checks
            if len(df_vars)==0: raise ValueError("cannot be 0, %s - %s"%(paper_phenotype, gene_name))
            if any(df_vars.n_transitions_var==0): ValueError("all vars should be associated to some transition")
            if len(set(df_vars.nodes_withPheno))!=1: raise ValueError("should have the same nodes_withPheno")
            if len(set(df_vars.Gene))!=1: raise ValueError("should have the same nodes_withPheno")

            # keep certain rows
            set_nodes_transitions = set.union(*df_vars.nodes_transitions.apply(lambda x: set(x.split(",")))).difference({"none"})

            row_dict = {"phenotype":paper_phenotype,
                        "Gene":geneID, 
                        "final_name":gene_name,
                        "type evidence":type_hit,
                        "inclusion criterion": inclusion_criteria,
                        "variant effects": "; ".join(sorted(df_vars["variant effect"])),
                        "n_transitions_total": df_vars.nodes_withPheno.iloc[0],
                        "n_transitions_gene": len(set_nodes_transitions),
                        "nodes_transitions_gene": ",".join(sorted(set_nodes_transitions)),
                        "description": df_vars.description.iloc[0],
                        "GO_BP": df_vars["GO_BP_Gene"].iloc[0],
                        "GO_CC": df_vars["GO_CC_Gene"].iloc[0],
                        "GO_MF": df_vars["GO_MF_Gene"].iloc[0],
                        "MetaCyc": df_vars["MetaCyc_Gene"].iloc[0],
                        "Reactome": df_vars["Reactome_Gene"].iloc[0]}

            integrated_table_rows.append(row_dict)

            # checks
            if row_dict["n_transitions_gene"]>row_dict["n_transitions_total"]: raise ValueError("invalid")

    df_integrated_table = pd.DataFrame(integrated_table_rows)
    if set(df_integrated_table.Gene)!=set(sorted_geneIDs): raise ValueError("invalid IDs")

    df_integrated_table["functional effects"] = df_integrated_table["variant effects"].apply(get_functional_effects_integrated_table)
    df_integrated_table["type_vars"] = df_integrated_table["variant effects"].apply(get_set_type_vars_integrated_table).apply(lambda x: "/".join(sorted(x)))

    # add the other phenotypes
    gene_to_phenotypes = df_integrated_table.groupby("Gene").apply(lambda df: ",".join(sorted(set(df.phenotype.apply(get_short_pheno_name)))))
    df_integrated_table["phenotypes_gene"] = df_integrated_table.Gene.map(gene_to_phenotypes)
    # write
    df_integrated_table.to_excel(table_file+".tmp.xlsx")
    os.rename(table_file+".tmp.xlsx", table_file)

    return pd.read_excel(table_file)

def get_functional_effects_integrated_table(variant_effects):

    """Gets whether it has functional effects"""

    variant_effects = set(variant_effects.split("; "))

    all_functional_effects = set()
    for eff in variant_effects:

        if any({eff.startswith(x) for x in {"3' INDEL", "5' INDEL", "3' SNP", "5' SNP", "syn|c.", "5' coverageDEL", "5' coverageDUP"}}): 
            functional_eff = False
        
        elif any({eff.startswith(x) for x in {"mis|p.", "trunc coverageDEL", "DUP coverageDUP", "FS|p.", "PTC|p.", "UTR5 coverageDUP", "UTR3 coverageDUP", "trunc DEL [c.", "codSeqVar_BND DELlike c.", "elongation coverageDUP [c."}}): 
            functional_eff = True

        else: raise ValueError("invalid %s"%eff)

        all_functional_effects.add(functional_eff)

    if any(all_functional_effects): return "yes"
    else: return "no"

def get_set_type_vars_integrated_table(variant_effects):

    """get the type vars set"""

    variant_effects = set(variant_effects.split("; "))

    all_type_vars = set()
    for eff in variant_effects:

        if any({eff.startswith(x) for x in {"3' SNP", "5' SNP", "mis|p.", "syn|c.", "PTC|p."}}): 
            type_var = "SNP"

        elif any({eff.startswith(x) for x in {"3' INDEL", "5' INDEL", "FS|p."}}): 
            type_var = "INDEL"
     
        elif any({eff.startswith(x) for x in {"5' coverageDEL", "5' coverageDUP", "trunc coverageDEL", "DUP coverageDUP", "UTR5 coverageDUP", "UTR3 coverageDUP", "elongation coverageDUP [c."}}): 
            type_var = "coverageCNV"

        elif any({eff.startswith(x) for x in {"trunc DEL [c.", "codSeqVar_BND DELlike c."}}): 
            type_var = "SV"

        else: raise ValueError("invalid '%s'"%eff)

        all_type_vars.add(type_var)

    return all_type_vars

def plot_ML_GWAS_overlap_heatmap(df_GWAS_hits_NR, df_ML_results_top, PlotsDir, DataDir, ProcessedDataDir, add_description=False, target_phenotypes="all", generate_proc_table=False):

    """Gets a heatmap with all the overlaps between all genes of high confidence / low GWAS / ML"""


    # get GWAS results with epsilon_gene
    df_GWAS_hits_NR = df_GWAS_hits_NR.copy().groupby(["type_filtering", "phenotype"]).apply(get_df_gwas_one_pheno_with_per_gene_stats).reset_index(drop=True)

    # get the ML models as in the heatmap
    df_ML_results_top, df_genes = ML_affected_genes_heatmap_top_models(df_ML_results_top, PlotsDir, DataDir, only_get_dfs=True)

    # get geneID to name
    gene_features_df = get_gene_features_df(DataDir).set_index("geneID", drop=False)
    geneID_to_name = dict(gene_features_df.final_name)
    tac1_ID = gene_features_df[gene_features_df.final_name=="Scer_HAL9"].iloc[0].geneID
    geneID_to_name[tac1_ID] = "TAC1"

    gene_name_to_geneID = {gname : gID for gID,gname in geneID_to_name.items()}
    gene_name_to_geneID["Scer_HAL9"] = tac1_ID

    geneID_to_description = dict(gene_features_df.description)

    # define the gwas dfs
    df_GWAS_high = df_GWAS_hits_NR[df_GWAS_hits_NR.type_filtering=="high_confidence"].copy()
    df_GWAS_low = df_GWAS_hits_NR[df_GWAS_hits_NR.type_filtering=="low_confidence"].copy()

    # define the genes that are manually attributed to each phenotype
    phenotype_to_manual_genes = {"fluconazole R": ["ERG11", "TAC1", "NDT80", "MRR1"],
                                 "voriconazole R": ["ERG11", "TAC1", "NDT80", "MRR1"],
                                 "posaconazole R": ["ERG11", "TAC1", "NDT80", "MRR1"]}

    phenotype_to_manual_genes = {p : {gene_name_to_geneID[x] for x in g} for p,g in phenotype_to_manual_genes.items()}
    all_manually_curated_genes = set.union(*phenotype_to_manual_genes.values())

    # add the final name
    df_ML_results_top["phenotype"] = df_ML_results_top.phenotype.apply(get_paper_dicotomic_pheno)
    target_phenos_ML = sorted(set(df_ML_results_top["phenotype"]))

    # load the expected genes df
    df_expected_gene_info = get_phenotype_to_all_possible_genes_for_ML(DataDir, ParentDir+"/Cparapsilosis_popgen/manual_data", "%s/df_expected_gene_info_ML.py"%ProcessedDataDir)

    # get pheno mappings
    pheno_to_paper_pheno =  {p : get_paper_dicotomic_pheno(p) for p in set(df_expected_gene_info.phenotype)}
    paper_pheno_to_pheno = {v:k for k,v in pheno_to_paper_pheno.items()}
    df_expected_gene_info["phenotype"] = df_expected_gene_info.phenotype.apply(get_paper_dicotomic_pheno)

    sorted_gene_groups = ['known', 'GWAS_high', 'by_description', 'selection_SNP']
    df_expected_gene_info = df_expected_gene_info[df_expected_gene_info.gene_group.isin(sorted_gene_groups)].copy()
    # all groups: {'known', 'by_description', 'selection_SNP', 'selection_all', 'GWAS_high', 'GWAS_low'}

    # define the target genes
    expected_phenotypes = set(df_ML_results_top["phenotype"]).union(set(df_GWAS_high.phenotype))
    sorted_phenotypes = ['fluconazole R', 'voriconazole R', 'posaconazole R', 'amphotericin B R', '5-FC R', 'chlorhexidine tol', 'age infant', 'age old', 'skin exudate', 'bronchial aspirate', 'ascitic fluid', 'blood culture', 'blood vs other sterile', 'blood', 'sterile']
    if expected_phenotypes!=set(sorted_phenotypes): raise ValueError("should be same gene set")
    sorted_geneIDs = sorted(set(df_genes.geneID).union(set(df_GWAS_high.Gene)).union(all_manually_curated_genes))

    # init the target genes as those that are expected (the ones to be plotted)
    #target_genes_all = {gene_name_to_geneID[x] for x in set.union(*[set(gene_list) for pheno,gene_list in get_literature_expected_genes_drug_types().items()])}
    target_genes_all = set()

    # print numbers of genes
    print("There are %i GWAS high-confidence, %i ML and %i manual genes. %i in total, for %i phenotypes"%(len(set(df_GWAS_high.Gene)), len(set(df_genes.geneID)), len(all_manually_curated_genes), len(sorted_geneIDs), len(sorted_phenotypes)))

    # get a df where each column is a phenotype, and each column is a gene
    df_plot = pd.DataFrame()
    df_plot["geneID"] = sorted_geneIDs

    df_annot = pd.DataFrame()
    df_annot["geneID"] = sorted_geneIDs

    for phenotype in sorted_phenotypes:

        # add type of gene
        geneIDs_gwas_highc = set(df_GWAS_high[df_GWAS_high.phenotype==phenotype].Gene).difference({"-"})
        geneIDs_gwas_lowc = set(df_GWAS_low[df_GWAS_low.phenotype==phenotype].Gene).difference({"-"})

        if phenotype in target_phenos_ML:
            geneIDs_ML = set.union(*df_ML_results_top[df_ML_results_top.phenotype==phenotype].all_affected_genes_set).intersection(set(df_genes.geneID))

        else: geneIDs_ML = set()

        if phenotype in phenotype_to_manual_genes:
            geneIDs_manual = phenotype_to_manual_genes[phenotype]

        else: geneIDs_manual = set()

        df_plot[phenotype] = df_plot.geneID.apply(get_type_genes_ML_GWAS_overlaps, genes_ML=geneIDs_ML, genes_GWAS_high=geneIDs_gwas_highc, genes_GWAS_low=geneIDs_gwas_lowc, geneIDs_manual=geneIDs_manual)

        ##### DEFINE THE TARGET GENES #######

        # define the target genes
        df_expected_genes_pheno = df_expected_gene_info[df_expected_gene_info.phenotype==phenotype].copy()
        gene_group_to_genes = dict(df_expected_genes_pheno.groupby("gene_group").apply(lambda df_g: set(df_g.geneID)))
        if len(df_expected_genes_pheno)==0: raise ValueError("invalid %s, %s"%(phenotype, df_expected_genes_pheno))
        if not "known" in gene_group_to_genes: gene_group_to_genes["known"] = set()

        # init with the ML and manually curated genes
        known_genes_found = (gene_group_to_genes["known"]).intersection(geneIDs_ML.union(geneIDs_gwas_highc).union(geneIDs_manual))
        dict_target_genes = {"known": known_genes_found, "ML": geneIDs_ML.difference(known_genes_found)}

        # add the target genes by different expected sets of genes
        sorted_gene_groups_pheno = [gg for gg in sorted_gene_groups if gg in gene_group_to_genes and gg!="known"]
        for gene_group in sorted_gene_groups_pheno: dict_target_genes[gene_group] = set()
        for gID in geneIDs_gwas_highc:
            if gID in geneIDs_ML or gID in known_genes_found: continue

            for gene_group in sorted_gene_groups_pheno:
                if gID in gene_group_to_genes[gene_group]: 
                    dict_target_genes[gene_group].add(gID)
                    break

        # get dfs
        df_GWAS_high_pheno = df_GWAS_high[df_GWAS_high.phenotype==phenotype].copy()
        df_GWAS_low_pheno = df_GWAS_low[df_GWAS_low.phenotype==phenotype].copy()

        # add top hits if none of the others worked
        if len(set.union(*dict_target_genes.values()))==0:
            #sorting_fields = ["epsilon", "OR", "nodes_GenoAndPheno", "unique_groupID"]
            #top_hits = set(df_GWAS_high_pheno.sort_values(by=sorting_fields, ascending=[False, False, False, True])[sorting_fields].drop_duplicates().iloc[0:3].unique_groupID)
            #dict_target_genes["top 3 hits"] = set(df_GWAS_high_pheno[df_GWAS_high_pheno.unique_groupID.isin(top_hits)].Gene)
            sorting_fields = ["fraction_transitions_by_Gene", "max_epsilon", "max_OR", "max_nodes_GenoAndPheno", "Gene"]
            dict_target_genes["top 3 genes"] = set(df_GWAS_high_pheno.sort_values(by=sorting_fields, ascending=[False, False, False, False, True])[sorting_fields].drop_duplicates().iloc[0:3].Gene)
            
        # add the inclusion criteria as annotations
        criteria_to_annot = {"known":"k", "ML":"m", "GWAS_high":"g", "by_description":"d", "selection_SNP":"s", "top 3 genes":"t"}
        criteria_to_label_criteria = {"ML":"in ML features", "known":"known gene", "GWAS_high":"GWAS hit Schikora-Tamarit 2024", "by_description":"functional description", "selection_SNP":"under selection Schikora-Tamarit 2024", "top 3 genes":"top 3 genes"}
        df_annot[phenotype] = df_annot.geneID.apply(get_inclusion_criteria_annot, dict_target_genes=dict_target_genes, criteria_to_annot=criteria_to_annot)

        # fill empty annotations with the fraction_transitions_by_Gene
        gene_to_fraction_nodes_explained_low = dict(df_GWAS_low_pheno[["Gene", "fraction_transitions_by_Gene"]].drop_duplicates().set_index("Gene").fraction_transitions_by_Gene)
        gene_to_fraction_nodes_explained_high = dict(df_GWAS_high_pheno[["Gene", "fraction_transitions_by_Gene"]].drop_duplicates().set_index("Gene").fraction_transitions_by_Gene)
        df_annot[phenotype] = df_annot[[phenotype, "geneID"]].apply(get_annot_with_fraction_transitions_by_Gene, phenotype=phenotype, gene_to_fraction_nodes_explained_low=gene_to_fraction_nodes_explained_low, gene_to_fraction_nodes_explained_high=gene_to_fraction_nodes_explained_high, axis=1)

        # add
        target_genes_all.update(set.union(*dict_target_genes.values()))

        ######################################

    # filter to only keep target genes
    df_plot = df_plot[df_plot.geneID.isin(target_genes_all)].copy()
    df_annot = df_annot[df_annot.geneID.isin(target_genes_all)].copy()
    print("There are %i target genes"%(len(df_plot)))

    # generate the processed table
    if generate_proc_table is True:
        df_integrated = generate_GWAS_and_ML_integrated_table(PlotsDir+"/ML_and_GWAS_integrated_table.xlsx", df_plot, df_GWAS_hits_NR, df_ML_results_top, df_genes, df_annot, paper_pheno_to_pheno, criteria_to_annot, criteria_to_label_criteria)

        # prints
        print("There are %i functional effects"%(sum(df_integrated["functional effects"]=="yes")))
        print("There are %i functional effects with SVs / CNVs"%(sum(df_integrated["type_vars"].apply(lambda x: len(set(x.split("/")).intersection({"SV", "coverageCNV"}))>0))))

        df_integrated_CNVs = df_integrated[df_integrated["type_vars"].apply(lambda x: len(set(x.split("/")).intersection({"SV", "coverageCNV"}))>0)].copy()
        for I,r in df_integrated_CNVs.iterrows(): print(r.phenotype, r.final_name, [x for x in r["variant effects"].split("; ") if next(iter(get_set_type_vars_integrated_table(x))) in {"SV", "coverageCNV"}])

        print()
        df_nonFunctional = df_integrated[df_integrated["functional effects"]=="no"].copy()
        for I,r in df_nonFunctional.iterrows(): print(r.phenotype, r.final_name, r["variant effects"])

    # change some phenos
    if target_phenotypes!="all":
        strange_phenos = set(target_phenotypes).difference(expected_phenotypes)
        if len(strange_phenos)>0: raise ValueError("invalid %s"%strange_phenos)
        sorted_phenotypes = [p for p in sorted_phenotypes if p in target_phenotypes]

    df_plot = df_plot[["geneID"] + sorted_phenotypes].copy()
    df_annot = df_annot[["geneID"] + sorted_phenotypes].copy()

    # filter some
    df_plot = df_plot[df_plot[sorted_phenotypes].apply(lambda r: any(r!="none"), axis=1)].copy()
    df_annot = df_annot[df_annot.geneID.isin(set(df_plot.geneID))].copy()

    # change to numeric representation
    unique_values_df = set.union(*df_plot[sorted_phenotypes].apply(lambda r: set(r), axis=1))
    sorted_type_genes = ["ML & GWAS high-confidence", "ML & GWAS low-confidence", "ML only", "GWAS high-confidence", "GWAS low-confidence", "manual curation only", "none"]
    sorted_type_genes = [x for x in sorted_type_genes if x in unique_values_df]

    type_gene_to_color = {"ML & GWAS high-confidence":"purple", "ML & GWAS low-confidence":"tab:red", "ML only":"salmon", "GWAS high-confidence":"navy", "GWAS low-confidence":"c", "none":"white", "manual curation only":"silver"}
    type_genes_to_I = dict(zip(sorted_type_genes, range(len(sorted_type_genes))))

    df_plot = df_plot.set_index("geneID").map(lambda x: type_genes_to_I[x])
    df_plot.index = [geneID_to_name[gID] for gID in df_plot.index]

    df_annot = df_annot.set_index("geneID")
    df_annot.index = [geneID_to_name[gID] for gID in df_annot.index]


    # plot
    cm = sns.clustermap(df_plot, square=True, cmap=[type_gene_to_color[x] for x in sorted_type_genes], col_cluster=False, row_cluster=True, yticklabels=1, cbar_kws={"label":"", "orientation":"vertical"}, alpha=1, linewidth=0.1, linecolor="gray", annot=df_annot,fmt="") #  , vmin=0.0


    # edit
    adjust_cm_positions(cm, df_plot.copy(), hm_height_multiplier=0.02, hm_width_multiplier=0.02, cc_height_multiplier=0.017, rc_width_multiplier=0.017, idx_delimiter="-", distance_btw_boxes=0.0025, cd_height=0.03, rd_width=0.03, cbar_width=0.15, cbar_height=0.02, nrow_fields=None, ncol_fields=None) # hm_height_multiplier = 0.02

    hm_pos = cm.ax_heatmap.get_position()
    if add_description is True: x_offset_bar = 1.8
    else: x_offset_bar = 0.14

    cm.ax_cbar.set_position([hm_pos.x0 + hm_pos.width + x_offset_bar, hm_pos.y0 + hm_pos.height - 0.15, 0.02, 0.15])

    # set ticks
    yrange = cm.ax_cbar.get_ylim()[1] - cm.ax_cbar.get_ylim()[0]
    increment_y = yrange / len(sorted_type_genes)
    cm.ax_cbar.set_yticks([increment_y/2 + increment_y*y for y in type_genes_to_I.values()])
    cm.ax_cbar.set_yticklabels([ylabel for ylabel in sorted_type_genes])
    cm.ax_cbar.text(0, yrange+0.25, "type association", fontweight="bold")

    # add y ticks
    gene_names_sorted_y = [y.get_text() for y in cm.ax_heatmap.get_yticklabels()]
    if add_description is True:
        cm.ax_heatmap.set_yticklabels(["%s: %s"%(gname, geneID_to_description[gene_name_to_geneID[gname]]) for gname in gene_names_sorted_y], rotation=0)
    
    else:
        cm.ax_heatmap.set_yticklabels(gene_names_sorted_y, rotation=0)

    # add x ticks
    cm.ax_heatmap.set_xticklabels([x.get_text() for x in cm.ax_heatmap.get_xticklabels()], rotation=90)



    # add title
    cm.ax_heatmap.set_title("%i / %i prioritized genes"%(len(df_plot), len(sorted_geneIDs)), fontsize=10)

    # add legend about the inclusion criteria
    cm.ax_cbar.text(0, -2, "inclusion criteria", fontweight="bold")
    for I, (criteria, annot) in enumerate(criteria_to_annot.items()):
        cm.ax_cbar.text(0, -I-3, "%s: %s"%(annot, criteria_to_label_criteria[criteria]))


    # add legend about the fraction transitions
    cbar_pos = cm.ax_cbar.get_position()

    print(cbar_pos.width)
    cm.ax_cbar.text(0, -I-5.5, ".0-.9:", fontweight="bold")
    cm.ax_cbar.text(2.5, -I-5.5, "fraction transitions explained by gene")
    # save
    save_fig(cm, "prioritized_genes_ML_and_GWAS/heatmap_description=%s_target_phenotypes=%s.pdf"%(add_description, target_phenotypes), PlotsDir)
    plt.show()


def plot_ML_GWAS_overlap(df_GWAS_hits_NR, df_ML_results_top, PlotsDir, DataDir):


    """Plots the overlaps between GWAS and ML results"""

    # get the ML models as in the heatmap
    df_ML_results_top, df_genes = ML_affected_genes_heatmap_top_models(df_ML_results_top, PlotsDir, DataDir, only_get_dfs=True)

    # get geneID to name
    gene_features_df = get_gene_features_df(DataDir).set_index("geneID", drop=False)
    geneID_to_name = dict(gene_features_df.final_name)
    tac1_ID = gene_features_df[gene_features_df.final_name=="Scer_HAL9"].iloc[0].geneID
    geneID_to_name[tac1_ID] = "TAC1"

    # add the final name
    df_ML_results_top["phenotype"] = df_ML_results_top.phenotype.apply(get_paper_dicotomic_pheno)

    # go through each phenotype
    plots_dir = PlotsDir + "/venn_diagrams_ML_vs_GWAS"; make_folder(plots_dir)
    target_phenos_ML = sorted(set(df_ML_results_top["phenotype"]))

    for phenotype in target_phenos_ML:

        # get dfs
        df_GWAS_pheno = df_GWAS_hits_NR[df_GWAS_hits_NR.phenotype==phenotype].copy()
        df_ML_pheno = df_ML_results_top[df_ML_results_top.phenotype==phenotype].copy()

        print()
        # make venn
        geneIDs_gwas_highc = set(df_GWAS_pheno[df_GWAS_pheno.type_filtering=="high_confidence"].Gene).difference({"-"})
        geneIDs_gwas_lowc = set(df_GWAS_pheno[df_GWAS_pheno.type_filtering=="low_confidence"].Gene).difference({"-"})
        geneIDs_ML = set.union(*df_ML_pheno.all_affected_genes_set).intersection(set(df_genes.geneID))

        # check the overlaps in venn
        gene_set_to_genes = {"GWAS high-confidence": geneIDs_gwas_highc, "GWAS low-confidence": geneIDs_gwas_lowc, "ML":geneIDs_ML}
        gene_set_to_color = {"GWAS high-confidence": "navy", "GWAS low-confidence": "c", "ML":"tab:red"}
        sorted_gene_sets = [gs for gs,genes in gene_set_to_genes.items() if len(genes)>0]
        sorted_gene_sets_genes = [gene_set_to_genes[gs] for gs in sorted_gene_sets]
        sorted_colors = [gene_set_to_color[gs] for gs in sorted_gene_sets]
        labels = ["" for x in sorted_colors]

        fig = plt.figure(figsize=(4, 4))
        if len(sorted_gene_sets)==2:
            venn = venn2(sorted_gene_sets_genes, labels, sorted_colors)

        elif len(sorted_gene_sets)==3:
            venn = venn3(sorted_gene_sets_genes, labels, sorted_colors)

        else: raise ValueError("invalid n")

        plt.title(phenotype, fontsize=12)
        plt.show()

        save_fig(fig, "venn_ML_vs_GWAS_%s.pdf"%phenotype, plots_dir)

        # print overlaps between gene sets
        already_compared_gs = set()
        for gs1, genes1 in gene_set_to_genes.items():
            for gs2, genes2 in gene_set_to_genes.items():

                gs_tuple = tuple(sorted([gs1,gs2]))
                if gs_tuple in already_compared_gs or gs1==gs2: continue

                overlapping_genes = genes1.intersection(genes2)
                if len(overlapping_genes)>0: 
                    print("\n", phenotype, gs1, "overlap", gs2, ", ".join(sorted([geneID_to_name[g] for g in overlapping_genes])))

                already_compared_gs.add(gs_tuple)



def ML_affected_genes_heatmap_top_models(df_ML_results_top, PlotsDir, DataDir, only_get_dfs=False):

    """Gets the results of top models"""

    # keep
    df_ML_results_top = df_ML_results_top.copy()

    # add as set
    def get_set_sel_pred(x):
        if x=="nan": return set()
        else:
            preds = set(x.split(","))
            if not all(["|||" in x for x in preds]): raise ValueError("invalid %s"%x)
            return preds

    df_ML_results_top["selected_predictors"] = df_ML_results_top.selected_predictors.apply(str).apply(get_set_sel_pred).apply(get_selected_predictors_demult_one_tuple_predictors)

    #df_ML_results_top["all_affected_genes_set"] = df_ML_results_top.all_affected_genes_str.apply(str).apply(lambda x: set(x.split(";")))
    df_ML_results_top["all_affected_genes_set"] = df_ML_results_top.selected_predictors.apply(lambda x: set(map(get_affected_geneID_one_unique_groupID, x)))

    # keep
    df_ML_results_top = df_ML_results_top[(df_ML_results_top.test_sample_set=="representative_samples") & ~(df_ML_results_top.phenotype.isin({"MIF_resistance", "agar_invasiveness"})) & (df_ML_results_top.selected_predictors.apply(len)>0)].copy()
    if len(df_ML_results_top)!=len(df_ML_results_top[["phenotype", "test_split_ID", "test_parmID_exceptTestSS"]].drop_duplicates()):
        raise ValueError("invalid lens")

    # define sorted
    sorted_phenotypes = ['VRC_resistance', 'FLC_resistance', 'age_infant', 'source_blood_vs_other_sterile']
    if set(sorted_phenotypes)!=set(df_ML_results_top.phenotype):
        raise ValueError("invalid phenos")

    # get the genes
    all_genes_sorted = sorted(set.union(*df_ML_results_top.all_affected_genes_set))
    df_genes = pd.DataFrame({"geneID":all_genes_sorted})
    df_genes_n = pd.DataFrame({"geneID":all_genes_sorted})
    pheno_to_nmodels = {}

    for phenotype in sorted_phenotypes:
        df_ML_results_top_pheno = df_ML_results_top[df_ML_results_top.phenotype==phenotype].copy()
        list_gene_sets = list(df_ML_results_top_pheno.all_affected_genes_set)
        if any([len(x)<1 for x in list_gene_sets]): raise ValueError("should be many genes")
        df_genes_n[phenotype] = df_genes_n.geneID.apply(lambda g: sum([g in target_set for target_set in list_gene_sets])) 
        df_genes[phenotype] = df_genes_n[phenotype] / len(df_ML_results_top_pheno)
        pheno_to_nmodels[phenotype] = len(df_ML_results_top_pheno)

    # add name
    gene_features_df = get_gene_features_df(DataDir).set_index("geneID", drop=False)
    geneID_to_name = dict(gene_features_df.final_name)
    tac1_ID = gene_features_df[gene_features_df.final_name=="Scer_HAL9"].iloc[0].geneID
    geneID_to_name[tac1_ID] = "TAC1"
    df_genes["gene_name"] = df_genes.geneID.map(geneID_to_name)
    df_genes_n["gene_name"] = df_genes_n.geneID.map(geneID_to_name)

    # print genes
    target_gene_IDs = set(df_genes.geneID)

    for phenotype in sorted_phenotypes:
        if only_get_dfs is False: print("\n", phenotype)
        df_ML_results_top_pheno = df_ML_results_top[(df_ML_results_top.phenotype==phenotype)].copy()    

        for I,r in df_ML_results_top_pheno.iterrows():
            r_genes = r.all_affected_genes_set.intersection(target_gene_IDs)

            predictors_changed = []
            for p in r.selected_predictors:
                gene_p = get_affected_geneID_one_unique_groupID(p)
                predictors_changed.append(p.replace(gene_p, geneID_to_name[gene_p]))

            #print(r.test_parmID_exceptTestSS, r.test_split_ID, predictors_changed)
            if only_get_dfs is False: print(r.test_parmID_exceptTestSS, r.test_split_ID, [geneID_to_name[g] for g in sorted(r_genes)], "ngenes=%i"%len(r_genes))


    # filter genes in suff models
    df_genes = df_genes[(df_genes[sorted_phenotypes]>=0.1).apply(any, axis=1)].copy()

    if only_get_dfs is True:
        return df_ML_results_top, df_genes

    # define dfs
    df_genes = df_genes.set_index("gene_name")[sorted_phenotypes].copy()
    def get_n_models_str(x):
        if x==0: return ""
        else: return str(x)
    annot_df = df_genes_n[df_genes_n.gene_name.isin(set(df_genes.index))].set_index("gene_name")[sorted_phenotypes].loc[df_genes.index].copy().map(get_n_models_str)



    # make a clustermap
    cm = sns.clustermap(df_genes, square=True, cbar="rocket", col_cluster=True, yticklabels=1, cbar_kws={"label":"fraction top models w/ gene"}, annot=annot_df, fmt="", vmin=0.0)

    # edit
    adjust_cm_positions(cm, df_genes.copy(), hm_height_multiplier=0.02, hm_width_multiplier=0.02, cc_height_multiplier=0.017, rc_width_multiplier=0.017, idx_delimiter="-", distance_btw_boxes=0.0025, cd_height=0.03, rd_width=0.03, cbar_width=0.02, cbar_height=0.15, nrow_fields=None, ncol_fields=None)


    hm_pos = cm.ax_heatmap.get_position()
    cm.ax_cbar.set_position([hm_pos.x0 + hm_pos.width + 0.12, hm_pos.y0 + 0.1, 0.02, 0.15])


    cm.ax_heatmap.set_xticklabels([get_paper_dicotomic_pheno(x.get_text()) + " (%i)"%(pheno_to_nmodels[x.get_text()]) for x in cm.ax_heatmap.get_xticklabels()])
    cm.ax_heatmap.set_ylabel("")

    # save
    save_fig(cm, "clustermap_ML_top_genes.pdf", PlotsDir)

    

def get_sorted_phenos_by_hits_and_type(phenoPaper_to_nhits, all_gwas_phenotypes):

    """gets the sorted hits"""

    # get the nhits
    phenoPaper_to_nhits_all = {}
    for pheno in all_gwas_phenotypes:
        pheno_paper = get_paper_dicotomic_pheno(pheno)

        if pheno_paper in phenoPaper_to_nhits.index: 
            phenoPaper_to_nhits_all[pheno_paper] = phenoPaper_to_nhits[pheno_paper]

        else:
            phenoPaper_to_nhits_all[pheno_paper] = 0

    # define types of pheno
    groups_phenos = [['5-FC R', 'amphotericin B R', 'anidulafungin R', 'caspofungin R', 'fluconazole R', 'micafungin R', 'posaconazole R', 'voriconazole R'],
                     ['chlorhexidine tol', 'hypochlorite tol', 'zwittergent 3-14 tol'],
                     ['agar invasiveness'],
                     ['abscess', 'ascitic fluid', 'blood', 'blood culture', 'blood vs other sterile', 'bronchial aspirate', 'bronchoalveolar lavage', 'skin', 'skin exudate', 'sputum', 'sterile', 'urine'],
                     ['age infant', 'age old', 'age young']]

    sorted_phenos = []
    for group_phenos in groups_phenos:
        sorted_phenos += list(reversed(sorted(group_phenos, key=(lambda x: (phenoPaper_to_nhits_all[x], x) ))))

    # check
    if sorted(map(get_paper_dicotomic_pheno, all_gwas_phenotypes))!=sorted(sorted_phenos):
        raise ValueError("should be the same")

    return sorted_phenos

def get_stacked_barplots_gwas_hits(df_GWAS_hits_NR, PlotsDir, colors_dict, target_fields, val_to_bgcolor, val_to_symbol):

    """get # GWAS hits"""

    # get df
    df_GWAS_hits_NR = df_GWAS_hits_NR.copy()

    # create df
    df_plot = pd.DataFrame()
    sorted_huevalues = []
    for f in target_fields:

        val_to_nhits = df_GWAS_hits_NR[[f, "phenotype", "unique_groupID"]].drop_duplicates().groupby(f).apply(len).sort_values()
        print(val_to_nhits)

        for val in list(val_to_nhits.index):
            print(val)

            n_hits = len(df_GWAS_hits_NR[df_GWAS_hits_NR[f]==val][["phenotype", "unique_groupID"]].drop_duplicates().copy())
            df_val = pd.DataFrame({"n_hits": [n_hits]})
            df_val["grouping"] = f
            df_val["val"] = val 

            df_plot = pd.concat([df_plot, df_val]).reset_index(drop=True)

            # add hues
            if not val in sorted_huevalues:
                sorted_huevalues.append(val)

    # plot
    xfield = "grouping"
    yfield = "n_hits"
    hue = "val"
    palette_plot = {}
    sorted_xvalues = target_fields
    for f in target_fields: palette_plot = {**palette_plot, **colors_dict[f]}
    sorted_huevalues = list(df_plot["val"])

    fig = plt.figure(figsize=(1.5, 2.5))
    ax = get_stacked_barplot_noLegend(df_plot, xfield, yfield, hue, sorted_xvalues, sorted_huevalues, palette_plot, linewidth_edges=None)

    ax = fig.axes[0]
    ax.set_xticklabels(target_fields, rotation=90)
    ax.set_ylabel("# NR GWAS hits")

    sns.despine()
    save_fig(fig, "stacked_bars_GWAS_hit_types.pdf", PlotsDir)

def get_heatmap_GWAS_hits(df_GWAS_hits_NR, sorted_phenotypes_gwas, PlotsDir):

    """gets the heatmap"""

    # add group ID
    df_GWAS_hits_NR = df_GWAS_hits_NR.copy()
    df_GWAS_hits_NR["grouping_ID"] = df_GWAS_hits_NR.unique_groupID.apply(lambda x: [x.split("|||")[0].split("|")[I] for I in [0, 2, 3]]).apply(tuple)


    # get the heatmap
    df_nhits = df_GWAS_hits_NR.groupby(["phenotype", "grouping_ID"]).apply(lambda df: pd.Series({"n_hits": len(set(df.unique_groupID))})).reset_index()
    df_square = df_nhits.pivot(index="phenotype", columns="grouping_ID", values="n_hits") #.map(get_nan_to_zero, fun_app=int)
    pheno_to_nhits = dict(df_nhits.groupby("phenotype").apply(lambda df: sum(df.n_hits)))

    # add missings
    missing_phenos = sorted(set(sorted_phenotypes_gwas).difference(set(df_square.index)))
    df_square = pd.concat([df_square, pd.DataFrame(np.nan, index=missing_phenos, columns=sorted(set(df_nhits.grouping_ID)))])
    for p in missing_phenos:
        pheno_to_nhits[p] = 0

    # define the sorted groupoings
    #sorted_groupings = sorted(set(df_GWAS_hits_NR["grouping_ID"]))
    sorted_groupings = list(df_nhits.groupby("grouping_ID").apply(lambda df: sum(df.n_hits)).sort_values(ascending=False).index)

    # get df
    df_square = df_square.loc[sorted_phenotypes_gwas, sorted_groupings].copy()

    # get df annot
    def get_annot_one_cell(x):
        if x==0 or pd.isna(x): return ""
        else: return str(int(x))

    df_annot = df_square.map(get_annot_one_cell)


    # define col colors
    colors_dict = {"type_vars": {'all_vars': 'black', 'coverageCNVs': 'blue', 'small_vars': 'red', 'small_vars_and_CNVs': 'pink'},
                   "type_mutations": {'all_muts': 'black', 'non_syn_muts': 'purple', 'non_syn_non_truncating_muts': 'violet', 'truncating_muts': 'indigo'},
                   "type_collapsing": {"none":"black", "domains":"peru", "genes":"gold", "GO":"olive", "MetaCyc":"lightgreen", "Reactome":"green", "variant_effect":"gray"}}

    target_fields = ["type_vars", "type_mutations", "type_collapsing"]
    df_col_colors = pd.DataFrame({g : {target_fields[I]: colors_dict[target_fields[I]][g[I]] for I in [0, 1, 2]} for g in sorted_groupings}).transpose()

    # change things
    for df in [df_annot, df_square]:
        df.columns = ["|".join(c) for c in df.columns]

    df_col_colors.index = ["|".join(idx) for idx in df_col_colors.index]
    df_col_colors = df_col_colors.loc[df_square.columns].copy()
    sorted_groupings  = ["|".join(gID) for gID in sorted_groupings]

    # make a clustermap
    cm = sns.clustermap(df_square, col_cluster=False, row_cluster=False, row_colors=None, col_colors=df_col_colors, cmap="rocket_r", xticklabels=False, yticklabels=1, cbar_kws={"label":"# NR GWAS hits", "orientation":'vertical'},  annot=df_annot, fmt="", linewidth=.1, linecolor="gray")

    # edit
    adjust_cm_positions(cm, df_square.copy(), hm_height_multiplier=0.02, hm_width_multiplier=0.02, cc_height_multiplier=0.02, rc_width_multiplier=0.017, idx_delimiter="|", distance_btw_boxes=0.005, cd_height=0.03, rd_width=0.03, cbar_width=0.02, cbar_height=0.15, nrow_fields=None, ncol_fields=3)

    # adjusts
    cm.ax_heatmap.set_ylabel("")

    # add labels
    symbols_dict = {"type_vars": {},
                    "type_mutations": {"non_syn_muts":"n", "truncating_muts":"x"},
                    "type_collapsing": {"domains":"d", "genes":"g", "Reactome":"r"}}

    val_to_symbol = {}
    for k,cdict in symbols_dict.items(): 
        for val, symb in cdict.items(): val_to_symbol[val] = symb

    val_to_bgcolor = {"non_syn_muts":"white", "truncating_muts":"white", "domains":"black", "genes":"black", "Reactome":"white"}

    for Ic, gID in enumerate(sorted_groupings):
        for Ir in [0, 1, 2]:

            f = target_fields[Ir]
            val = gID.split("|")[Ir]
            f_dict = symbols_dict[f]

            if val in f_dict:
                cm.ax_col_colors.text(Ic+0.5, Ir+0.5, f_dict[val], ha="center", va="center", fontsize=10, color=val_to_bgcolor[val]) # ontweight="bold",

    # add the ylabels
    cm.ax_heatmap.set_yticklabels(["%s (%i)"%(pheno, pheno_to_nhits[pheno]) for pheno in sorted_phenotypes_gwas])

    # add the legend
    plt.rcParams['legend.handlelength'] = 1.2
    plt.rcParams['legend.handleheight'] = plt.rcParams['legend.handlelength'] * 1.125

    # add a legend for the row colors
    def get_lel(facecolor, label, edgecolor="gray"): return mpatches.Patch(facecolor=facecolor, edgecolor=edgecolor, label=label)

    legend_elements = []
    for f in target_fields:
        if f!="type_vars": legend_elements.append(get_lel("white", "", edgecolor="white"))
        legend_elements.append(get_lel("white", f, edgecolor="white"))
        legend_elements += [get_lel(color, val) for val, color in colors_dict[f].items()]

    legend = cm.ax_heatmap.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(0, 1))
    
    # add symbols
    ax_pos = cm.ax_heatmap.get_position()
    for Ip, patch in enumerate(legend_elements):

        x = ax_pos.x0  - 12.3
        y = ax_pos.y1 + Ip*1.134 + 0.38
        label = patch.get_label()

        if label in val_to_symbol: 
            print(label, val_to_symbol[label])
            cm.ax_heatmap.text(x, y, val_to_symbol[label], color=val_to_bgcolor[label], ha='center', va='center', fontsize=10, zorder=10)#, fontweight='bold') # val_to_bgcolor[val]

    print("hm_height", cm.ax_heatmap.get_position().height)

    # save
    save_fig(cm, "get_heatmap_GWAS_hits.pdf", PlotsDir)

    # get the stacked error bar
    get_stacked_barplots_gwas_hits(df_GWAS_hits_NR, PlotsDir, colors_dict, target_fields, val_to_bgcolor, val_to_symbol)


def get_short_pheno_name(paper_pheno):

    """gets short pheno name"""

    specific_changes_dict = {"5-FC":"5FC", "age old":"old", "age young":"young", "age infant":"infant", "agar invasiveness":"inv", "chlorhexidine tol":"chlor", "hypochlorite tol":"hypo", "bronchial aspirate":
    "b.a.", "skin exudate":"sk. ex", "blood culture":"b.c.", "blood":"blood", "skin":"skin", "bronchoalveolar lavage":"b.l.", "blood vs other sterile":"b-vs-o"}

    if paper_pheno in specific_changes_dict:
        return specific_changes_dict[paper_pheno]

    else:
        return paper_pheno[0:3]




def get_stacked_barplots_types_of_pheno_transitions_by_GWAS_hits_and_expected_genes(df_GWAS_hits_NR, sorted_phenotypes_gwas, PlotsDir, DataDir, ProcessedDataDir, df_GWAS_parameter_benchmark_filt, set_raw_phenotypes, short_pheno=True, type_transitions="all", show_n_genes=False):

    """plots stacked bars"""

    # keep high conf
    df_GWAS_hits_NR = df_GWAS_hits_NR[df_GWAS_hits_NR.type_filtering=="high_confidence"].copy()
    #if set(sorted_phenotypes_gwas)!=set(df_GWAS_hits_NR.phenotype): raise ValueError("should be the same")

    pheno_to_paperPheno = {p : get_paper_dicotomic_pheno(p) for p in set_raw_phenotypes}
    paperPheno_to_pheno = {v:k for k,v in pheno_to_paperPheno.items()}

    # keep only cases where some targeet genes where found
    expected_genes_fields = [g+"_found" for g in ["MRR1", "ERG11", "Scer_HAL9", "NDT80", "FKS2"]]
    df_GWAS_parameter_benchmark_filt = df_GWAS_parameter_benchmark_filt[df_GWAS_parameter_benchmark_filt[expected_genes_fields].apply(any, axis=1)].copy()
    phenos_w_expected_genes = {"FLC_resistance", "VRC_resistance", "POS_resistance", "CAS_resistance"}
    df_GWAS_parameter_benchmark_filt = get_filtered_df_by_dict(df_GWAS_parameter_benchmark_filt, {"phenotype": phenos_w_expected_genes})
    df_GWAS_parameter_benchmark_filt["paper_phenotype"] = df_GWAS_parameter_benchmark_filt.phenotype.map(get_paper_dicotomic_pheno)

    # add transitions as sets
    df_GWAS_hits_NR["nodes_transitions_set"] = df_GWAS_hits_NR.nodes_transitions.apply(lambda x: set(x.split(",")))

    # define the filters
    filters_dict = cp.deepcopy(manual_curation_gwas_filt_parms)

    # get the df
    df_plot_all = pd.DataFrame()
    for phenotype in sorted_phenotypes_gwas:

        # get the transitions
        outdir_run_GWAS = "%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas/gwas_results_%s_min_sup_%i_ASR_%s"%(DataDir, paperPheno_to_pheno[phenotype].replace(" ", "_"), filters_dict["aligner"], filters_dict["sample_set"], filters_dict["min_support"], filters_dict["ASR_method"])
        df_pheno_info = get_tab_as_df_or_empty_df(outdir_run_GWAS + "/phenotypes_info.tab")
        all_nodes_transitions = set(df_pheno_info[df_pheno_info.phenotype_transition==1].node_name)

        # log trans
        #print(phenotype, len(all_nodes_transitions), "transitions")

        # get the transitions explained by GWAS hits
        df_GWAS_hits_NR_pheno = df_GWAS_hits_NR[df_GWAS_hits_NR.phenotype==phenotype].copy()
        if len(df_GWAS_hits_NR_pheno)>0: 
            GWAS_nodes_transitions = set.union(*df_GWAS_hits_NR_pheno.nodes_transitions_set)

        else:
            GWAS_nodes_transitions = set()

        strange_nodes_NR = GWAS_nodes_transitions.difference(all_nodes_transitions)
        if len(strange_nodes_NR)>0: raise ValueError("invalid nodes: %s"%strange_nodes_NR)

        # get the nodes that are explained by expected genes
        if phenotype.split()[0] in {"fluconazole", "posaconazole", "voriconazole", "caspofungin"}:
            df_parameters = df_GWAS_parameter_benchmark_filt[df_GWAS_parameter_benchmark_filt.phenotype==paperPheno_to_pheno[phenotype]].copy()

            # define fields
            if phenotype.split()[0] in {"fluconazole", "posaconazole", "voriconazole"}: expected_gene_names = ["MRR1", "ERG11", "Scer_HAL9", "NDT80"]
            elif phenotype.split()[0] in {"caspofungin"}: expected_gene_names = ["FKS2"]
            else: raise ValueError("no")

            expected_genes_fields_pheno = [g+"_found" for g in expected_gene_names]
            nodes_expected_genes_fields_pheno = [g+"_affected_nodes" for g in expected_gene_names]

            # redefine to keep the expected nodes
            expectedGenes_nodes_transitions = set()
            for f in nodes_expected_genes_fields_pheno:
                df_parameters[f] = df_parameters[f].apply(lambda x: x.intersection(all_nodes_transitions))
                gene_f = f.split("_affected_nodes")[0]
                df_parameters[gene_f+"_found"] =  (df_parameters[f].apply(len)>0)
                expectedGenes_nodes_transitions.update(set.union(*df_parameters[f]))

            # keep parameters with expected nodes
            df_parameters = df_parameters[df_parameters[expected_genes_fields_pheno].apply(any, axis=1)].copy()

            # get the genes
            parm_to_genes = df_parameters[expected_genes_fields_pheno].apply(lambda r: {x.split("_found")[0] for x in r[r].index}, axis=1)
            if any(parm_to_genes.apply(len)==0): 
                #print("No expected gene")
                expectedGenes_found = set()

            else:
                expectedGenes_found = set.union(*parm_to_genes)

            print(phenotype, expectedGenes_found)

        else:
            expectedGenes_nodes_transitions = set()
            expectedGenes_found = set()

        # add to df
        n_expected_genes = len(expectedGenes_nodes_transitions)
        n_additional_GWAS = len(GWAS_nodes_transitions.difference(expectedGenes_nodes_transitions))
        n_none = len(all_nodes_transitions.difference(GWAS_nodes_transitions.union(expectedGenes_nodes_transitions)))
        df_plot = pd.DataFrame({"type transitions":["by expected genes", "by GWAS", "unexplained"], "# nodes":[n_expected_genes, n_additional_GWAS, n_none], "phenotype":phenotype})

        if show_n_genes is True:
            n_genes_GWAS = len(set(df_GWAS_hits_NR_pheno.Gene))
            df_plot["annotation_str"] = [str(len(expectedGenes_found)), str(n_genes_GWAS), ""]

        else:
            df_plot["annotation_str"] = ""

        # keep
        df_plot_all = pd.concat([df_plot_all, df_plot]).reset_index(drop=True)

    # filter to keep only some
    if type_transitions=="all":
        pass

    elif type_transitions=="none":
        df_plot_all["type transitions"] = df_plot_all["type transitions"].map({"by expected genes":"unexplained", "by GWAS":"unexplained", "unexplained":"unexplained"})

    elif type_transitions=="only_expected":
        df_plot_all["type transitions"] = df_plot_all["type transitions"].map({"by expected genes":"by expected genes", "by GWAS":"unexplained", "unexplained":"unexplained"})


    # merge the unexplained ones
    def get_series_one_type_trans_and_pheno(df): 
        pheno, type_t = df.name
        annotation_str = df.annotation_str.iloc[0] 
        if annotation_str=="0" or type_t=="unexplained": annotation_str=""
        return pd.Series({"phenotype":pheno, "type transitions": type_t, "# nodes":sum(df["# nodes"]), "annotation_str":annotation_str})

    df_plot_all = df_plot_all.groupby(["phenotype", "type transitions"]).apply(get_series_one_type_trans_and_pheno).reset_index(drop=True)

    # check
    check_no_nans_in_df(df_plot_all)

    # plot
    fig = plt.figure(figsize=(len(sorted_phenotypes_gwas)*0.28, 2.5))

    xfield = "phenotype"
    yfield = "# nodes"
    hue = "type transitions"
    palette_plot = {"by expected genes":"tab:blue", "by GWAS":"c", "unexplained":"gray"}
    sorted_xvalues = sorted_phenotypes_gwas
    sorted_huevalues = [x for x in list(palette_plot.keys()) if x in set(df_plot_all["type transitions"])]
    ax = get_stacked_barplot_noLegend(df_plot_all, xfield, yfield, hue, sorted_xvalues, sorted_huevalues, palette_plot, linewidth_edges=None, field_annot="annotation_str", hue_to_annot_color={"by expected genes":"white", "by GWAS":"black", "unexplained":"black"})

    ax = fig.axes[0]
    ax.set_yticks([0, 3, 6, 9, 12])

    if short_pheno is True:
        ax.set_xticklabels([get_short_pheno_name(p) for p in sorted_xvalues], rotation=90)

    else:
        ax.set_xticklabels(sorted_xvalues, rotation=90)

    ax.set_ylabel("# transition nodes")
    ax.axhline(3, ls="--", lw=.7, c="black")

    for I,x in enumerate(sorted_xvalues):
        if x in {"anidulafungin R", "hypochlorite tol" ,"agar invasiveness", "abscess"}:
            ax.axvline(I+0.5, ls="--", lw=.7, c="black")

    for label,x in [("resistance", 3), ("disinf.", 9), ("isolation source", 17.5), ("age", 25)]:
        ax.text(x, 14, label, fontsize=10, ha="center")


    # save
    sns.despine()
    save_fig(fig, "stacked_bars_transitions_nodes_covered_GWAS/transitions_GWAS_short_pheno=%s_type_transitions=%s_show_n_genes=%s.pdf"%(short_pheno, type_transitions, show_n_genes), PlotsDir)
    plt.show()
    
def print_ML_results(df_ML_stats):

    """print"""

    df_ML_stats = df_ML_stats.copy()

    # print # parameteres
    print("# train parameter combinations:", min(df_ML_stats.min_n_train_parameters), max(df_ML_stats.max_n_train_parameters))
    print("# test parameter combinations:", min(df_ML_stats.n_test_parameters), max(df_ML_stats.n_test_parameters))
    print(df_ML_stats.n_test_parameters)

    # print performance of best models
    print("best models info:")
    df_ML_stats_best = df_ML_stats[df_ML_stats.n_correct_test_parameters>0].copy()
    fields_print = ["best_model_AUC_test", "all_samples_AUC_test", "redundant_samples_AUC_test", "best_model_sample_set"]
    print(df_ML_stats_best[fields_print])


def get_df_metadata_all_samples_w_cladeID_definedBy_branchSupports_and_BL(df_metadata, strains_tree, min_support_tree, PlotsDir, best_min_fraction_branchLen):

    """Gets the df metadata according to the strains tree. Similar to the approach in Schikora-Tamarit 2023, Nature Communications."""

    # keep
    strains_tree = cp.deepcopy(strains_tree)

    #### INVESTIGATE BEST MIN_BL ####

    # get data
    print("Getting df_plot...")
    #linspace_list = list(np.linspace(0, 0.22, 100)) # whole range
    #linspace_list = list(np.linspace(0, 0.25, 500))
    #linspace_list = list(np.linspace(0.05, 0.09, 100))
    #linspace_list = list(np.linspace(0.055, 0.066, 100))
    linspace_list = list(np.linspace(0.055, 0.066, 30))

    inputs_fn = [(cp.deepcopy(strains_tree), min_fraction_branchLen, min_support_tree) for min_fraction_branchLen in (linspace_list + [best_min_fraction_branchLen])]

    with multiproc.Pool(multiproc.cpu_count()) as pool:
        list_dicts = pool.starmap(get_clade_info_one_min_fraction_branchLen, inputs_fn, chunksize=1)
        pool.close()

    df_plot = pd.DataFrame(list_dicts).sort_values(by="min_fraction_branchLen")

    # define the title to include the selected number of clades and the numbers of samples with some clade
    print("Getting plot")
    selected_r = df_plot[df_plot.min_fraction_branchLen==best_min_fraction_branchLen].iloc[0]
    nsamples = len(strains_tree.get_leaf_names())
    nsamples_with_clade = int(nsamples*selected_r.fraction_samples_clade)
    title = "#clades=%i; %i/%i (%.2f%s) samples with clade"%(int(selected_r.nclades), nsamples_with_clade, nsamples, selected_r.fraction_samples_clade*100, "%")

    # add plot
    ax = df_plot.plot(x="min_fraction_branchLen", y="fraction_samples_clade", legend=False, color="b")
    
    # add line with all
    plt.axhline(1.0, color="b", linestyle="--")

    # switch to ax2
    ax2 = ax.twinx()
    df_plot.plot(x="min_fraction_branchLen", y="nclades", ax=ax2, legend=False, color="r", linewidth=2)
    ax.set_title(title)

    # add best_min_fraction_branchLen
    plt.axvline(best_min_fraction_branchLen, color="k", linestyle="--")
    expected_nclades = 5
    if expected_nclades!=0: plt.axhline(expected_nclades, color="r", linestyle="--")

    # change the axis
    ax.set_ylabel("fraction samples with clade", color="b")
    ax.spines['left'].set_color('b')
    ax.tick_params(axis='y', colors='b')

    ax2.set_ylabel("total # clades", color="r")
    ax2.spines['right'].set_color('r')
    ax2.tick_params(axis='y', colors='r')

    plt.show()

    # save
    fig = ax.get_figure()
    fig.savefig("%s/choosing_thresholds_clades.pdf"%PlotsDir, bbox_inches="tight")


    ##############################

    #### DEFINING CLADES ####

    # map each leaf to the clade
    print("adding clade to df_metadata")
    strains_tree = get_tree_with_cladeID_from_supportAndBranchLength(strains_tree, best_min_fraction_branchLen, min_support_tree)
    leaf_to_clade = {l.name : l.cladeID for l in strains_tree.get_leaves()}
    df_metadata["cladeID_BL_and_support"] = df_metadata.runID.apply(lambda x: leaf_to_clade[x])

    #########################

    return df_metadata

def get_formatted_metadata_field(val, f):

    """for a given phenotypic field, get the empt"""

    # map each field to a function
    if f.startswith("% VIABILITY"): fn_modifier = get_float_or_minus1
    elif f in {'Agar invasiveness', 'MORPHOLOGY', 'Microfluidic behavior', 'Biofilm (High / Low / No)', 'genotype', 'Referencia cepa Hospital', 'date_collection', 'muestra', 'hospital', 'date_birth', 'Observaciones', 'CCAA', 'Agar invasiveness', 'MORPHOLOGY', 'Microfluidic behavior', 'Biofilm (High / Low / No)', 'ERG11_mutation', 'country_collection', 'city_collection'}: fn_modifier = get_str_or_unknown
    elif f in {'MIC ETEST 24H', 'MIC ETEST 48H', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h'}: fn_modifier = get_mic_float
    elif f in {'Colonies inside halo'}: fn_modifier = covert_yes_no
    elif f in {'duplicated_patientID', "year_collection", "age", "collection"}: fn_modifier = get_int_or_minus1
    elif f.endswith("_MIC"): fn_modifier = get_mic_float
    elif f.endswith("_resistance"): fn_modifier = get_str_or_unknown
    else: raise ValueError("invalid f %s"%f)

    return fn_modifier(val)

def visualize_continuous_phenotypes(df_metadata, PlotsDir):

    """Visualizing continuous phenos"""

    df_metadata = df_metadata.copy()
    df_metadata["xfield"] = 1

    # defin the phenotypes
    continuous_phenos =  ['FLC_MIC', 'MIC ETEST 24H', 'MIC ETEST 48H', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h', 'ITR_MIC', 'VRC_MIC', 'POS_MIC', 'IVZ_MIC', 'AMB_MIC', '5FC_MIC', 'CAS_MIC', 'MIF_MIC', 'ANI_MIC', '% VIABILITY ETHANOL', '% VIABILITY HYPLOCHLORITE', '% VIABILITY H202', '% VIABILITY CHLORHEXIDINE', '% VIABILITY ZWITTERGENT 3-14', '% VIABILITY SURFANIO']

    # init fig
    nrows = 4
    ncols = 5
    fig = plt.figure(figsize=(ncols*2, nrows*2.2))
    #fig = plt.figure(figsize=(ncols*6, nrows*6))

    if len(continuous_phenos)!=(nrows*ncols): raise ValueError("should have same dist")

    for I, pheno in enumerate(continuous_phenos):
        print(I, pheno)

        # get plot
        ax = plt.subplot(nrows, ncols, I+1)
        df_plot = df_metadata[df_metadata[pheno]!=-1].copy()

        # define hue f
        if pheno.endswith("_MIC") and pheno.split("_")[0] in {"FLC", "VRC", "ITR", "POS"}:
            hue_f = "%s_resistance"%(pheno.split("_")[0])
            palette = {"R":"red", "I":"gray", "S":"blue"}
            df_plot = df_plot.sort_values(by=[hue_f], ascending=False)
            edgecolor = [palette[x] for x in df_plot[hue_f]]

        else:
            hue_f = None
            palette = None
            edgecolor = "black"

        # get plot
        #ax = sns.swarmplot(data=df_plot, x="xfield", y=pheno, hue=hue_f, palette=palette, s=3, dodge=True)
        ax = sns.stripplot(data=df_plot, x="xfield", y=pheno, hue=hue_f, palette=palette, s=3, dodge=True, alpha=.5)#, linewidth=.5, edgecolor=edgecolor)

        # modify the face
        #for artist in ax.collections: artist.set_facecolor('none')



        if "MIC" in pheno:
            ax.set_yscale("log")


        ax.set_ylabel(pheno.replace("recovered heteroresistant", "rec. heteroR"))

        ax.set_xticks([])
        ax.set_xlabel("")
        ax.set_xticklabels([])

        if not palette is None:
            ax.get_legend().remove()

        #if I==5: break


    plt.subplots_adjust(wspace=0.8, hspace=0.1)# wspace (between cols) hspaces (btw rows)

    save_fig(fig, "distribution_continuous_phenotypes.pdf", PlotsDir)


def visualize_continuous_phenotypes_hist(df_metadata, PlotsDir):

    """Visualizing continuous phenos"""

    df_metadata = df_metadata.copy()
    df_metadata["xfield"] = 1

    # defin the phenotypes
    continuous_phenos =  ['FLC_MIC', 'MIC ETEST 24H', 'MIC ETEST 48H', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h', 'ITR_MIC', 'VRC_MIC', 'POS_MIC', 'IVZ_MIC', 'AMB_MIC', '5FC_MIC', 'CAS_MIC', 'MIF_MIC', 'ANI_MIC', '% VIABILITY ETHANOL', '% VIABILITY HYPLOCHLORITE', '% VIABILITY H202', '% VIABILITY CHLORHEXIDINE', '% VIABILITY ZWITTERGENT 3-14', '% VIABILITY SURFANIO']

    # init fig
    nrows = 4
    ncols = 5
    fig = plt.figure(figsize=(ncols*3, nrows*2))
    #fig = plt.figure(figsize=(ncols*6, nrows*6))

    if len(continuous_phenos)!=(nrows*ncols): raise ValueError("should have same dist")

    for I, pheno in enumerate(continuous_phenos):

        # get plot
        ax = plt.subplot(nrows, ncols, I+1)
        df_plot = df_metadata[df_metadata[pheno]!=-1].copy().sort_values(by=[pheno], ascending=True)

        # define hue f
        if pheno.endswith("_MIC") and pheno.split("_")[0] in {"FLC", "VRC", "ITR", "POS"}:
            hue_f = "%s_resistance"%(pheno.split("_")[0])
            palette = {"R":"red", "I":"gray", "S":"blue"}

        else:
            hue_f = None
            palette = None

        # get log
        df_plot["log %s"%pheno] = np.log10(df_plot[pheno])

        if "MIC" in pheno: 
            xfield = "log %s"%pheno

        else:
            xfield = pheno


        # print for MICs
        if "MIC" in pheno:

            if pheno in {"FLC_MIC", "VRC_MIC", "ITR_MIC", "POS_MIC"}:

                for rf in ["S", "I", "R"]:
                    if any(df_plot[hue_f]==rf):
                        print(pheno, rf, dict(Counter(df_plot[df_plot[hue_f]==rf][pheno])))
            else: 
                print(pheno, dict(Counter(df_plot[pheno])))


        # get plot
        ax = sns.histplot(data=df_plot, x=xfield, hue=hue_f, palette=palette)

        ax.set_xlabel(xfield.replace("recovered heteroresistant", "rec. heteroR"))
        ax.set_ylabel("")


        if not palette is None:
            ax.get_legend().remove()

        #if I==5: break


    plt.subplots_adjust(wspace=0.2, hspace=0.6)# wspace (between cols) hspaces (btw rows)


    save_fig(fig, "distribution_continuous_phenotypes_hist.pdf", PlotsDir)


def get_bps_for_pheno_field(pheno):


    """gets bps"""

    if pheno.endswith("_MIC"):
        bp_tuple = drug_to_MIC_breakpoints[pheno.split("_")[0]]

    elif pheno.startswith("% VIABILITY"):
        bp_tuple = disinfectant_to_breakpoints[pheno.split("% VIABILITY ")[1]]

    elif pheno in {'MIC ETEST 24H', 'MIC ETEST 48H', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h'}:
         bp_tuple = drug_to_MIC_breakpoints["FLC"]

    else:
        raise ValueError(pheno)

    return bp_tuple

def get_dicothomic_resistance(val, pheno):

    """gets the dicothomic resistance"""

    # return some cases
    if val==-1: return "unknown"

    # define resistance breakpooints
    bp_tuple = get_bps_for_pheno_field(pheno)

    # get type
    if val<=bp_tuple[0]:
        return "S"

    elif val>=bp_tuple[1]:
        return "R"

    else:
        return "I"


def visualize_continuous_phenotypes_hist_own_bps(df_metadata, PlotsDir):

    """Visualizing continuous phenos"""

    df_metadata = df_metadata.copy()
    df_metadata["xfield"] = 1

    # defin the phenotypes
    continuous_phenos =  ['FLC_MIC', 'MIC ETEST 24H', 'MIC ETEST 48H', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h', 'ITR_MIC', 'VRC_MIC', 'POS_MIC', 'IVZ_MIC', 'AMB_MIC', '5FC_MIC', 'CAS_MIC', 'MIF_MIC', 'ANI_MIC', '% VIABILITY ETHANOL', '% VIABILITY HYPLOCHLORITE', '% VIABILITY H202', '% VIABILITY CHLORHEXIDINE', '% VIABILITY ZWITTERGENT 3-14', '% VIABILITY SURFANIO']

    # init fig
    nrows = 4
    ncols = 5
    fig = plt.figure(figsize=(ncols*3.4, nrows*2.3))
    #fig = plt.figure(figsize=(ncols*6, nrows*6))

    if len(continuous_phenos)!=(nrows*ncols): raise ValueError("should have same dist")

    for I, pheno in enumerate(continuous_phenos):

        # get plot
        ax = plt.subplot(nrows, ncols, I+1)
        df_plot = df_metadata[df_metadata[pheno]!=-1].copy().sort_values(by=[pheno], ascending=True)

        # get resistance
        df_plot["resistance"] = df_plot[pheno].apply(get_dicothomic_resistance, pheno=pheno)
   
        # define xfield
        xfield = pheno

        # change some fields 
        if "MIC" in pheno and not pheno.endswith("_MIC"):
            df_plot[xfield] = df_plot[xfield].apply(lambda x: round(x, 2))


        # define the type of plot
        palette = {"R":"red", "I":"gray", "S":"blue"}
        if "MIC" in pheno:

            sorted_vals = sorted(set(df_plot[xfield]))
            df_plot["str_x"] = df_plot[xfield].apply(str)
            ax = sns.histplot(df_plot, x="str_x", hue='resistance',  multiple='stack', palette=palette, alpha=.8)
            ax.set_xticklabels(sorted_vals, rotation=90)


        else: 
            ax = sns.histplot(data=df_plot, x=xfield, hue="resistance", palette=palette, alpha=.8, bins=10)


        # get plot
        ax.set_xlabel(xfield.replace("recovered heteroresistant", "rec. heteroR"))
        ax.set_ylabel("")

        # add lines
        array_vals = np.array(sorted(set(df_plot[xfield])))


        for x in (list(get_bps_for_pheno_field(pheno)) + [np.median(array_vals)]):
        #for x in [np.median(array_vals)]:

            if "MIC" in pheno:

                nearest_val = find_nearest(array_vals, x)
                xval = [I for I,v in enumerate(array_vals) if v==nearest_val][0]

                # modify slightly
                if nearest_val>x:
                    xval = xval - 0.5

                elif nearest_val<x:
                    xval = xval + 0.5


            else: xval = x

            plt.axvline(xval, ls="--", color="k", lw=.7)

        ax.get_legend().remove()

    plt.subplots_adjust(wspace=0.2, hspace=0.8)# wspace (between cols) hspaces (btw rows)


    save_fig(fig, "distribution_continuous_phenotypes_hist.pdf", PlotsDir)


def get_df_metadata_all_samples(CurDir, PlotsDir):

    """Gets a df with all the samples to analyze. Note that Bergin2022_tableStrains.xlsx has information about many of these strains, and also clades."""

    print("getting df metadata")

    # define dirs
    ManualDataDir = "%s/manual_data"%CurDir
    ProcessedDataDir = "%s/processed_data"%CurDir
    df_metadata_file = "%s/df_metadata_all_samples.py"%ProcessedDataDir

    if file_is_empty(df_metadata_file):

        # init metadata with all the runs 
        df_metadata = load_object("%s/df_reads_info_filtered_afterQC.py"%ProcessedDataDir)[["runID", "source_data", "species_run", "sampleID", "possibly_triploid"]] # these are all the correct runs

        # keep only para
        df_metadata = df_metadata[df_metadata.species_run=="para"].copy()

        # define all the runs sequenced here
        df_read_info_bad = load_object("%s/df_reads_info_bad_runs.py"%ProcessedDataDir)
        all_sampleIDs_sequenced_here = set()
        for df in [df_metadata, df_read_info_bad]:
            all_sampleIDs_sequenced_here.update(set(df[df.source_data=="sequenced_here"].sampleID))

        # checks
        check_no_nans_in_df(df_metadata)

        # print the sequenced twice
        print(df_metadata[df_metadata.sampleID.isin({"L21-060", "L21-067"})].runID)

        # get the df metadata sequenced here
        df_metadata_Oscar = get_df_metadata_final(ManualDataDir, all_sampleIDs_sequenced_here).set_index("sampleID")
        df_metadata_Oscar["country_collection"] = "Spain"
        keys_df_Oscar = list(df_metadata_Oscar.keys())
        check_no_nans_in_df(df_metadata_Oscar)

        # load the SRA information
        df_sra = pd.read_csv("%s/SraRunTable_Cpara_WGS_06112023.csv"%ManualDataDir, sep=",")    
        df_sra = df_sra[df_sra.Run.isin(set(df_metadata.runID))].set_index("Run")
        # df_sra = df_sra[['BioProject', 'Center Name', 'Sample Name', 'Library Name', 'strain', 'geo_loc_name_country', 'geo_loc_name_country_continent', 'geo_loc_name', 'isolation_source', 'Collection_Date', 'isolate', 'Host', 'lat_lon', 'collected_by', 'Host_disease', 'Sample_name', 'external_id', 'INSDC_first_public', 'Organism_part', 'sample_type', 'env_local_scale', 'env_medium', 'Host_age', 'host_disease_outcome', 'host_disease_stage', 'host_sex', 'host_description', 'strain_name_alias', 'environmental_history', 'COUNTRY', 'host_associated', 'host_scientific_name', 'Source',  'year_isolated', 'geographic_location_(country_and/or_sea)', 'geographic_location_(region_and_locality)', 'host_tissue_sampled']]
        df_sra["type_isolate"] = df_sra.apply(sra_run_return_type_sample_clinical_or_enviromental_WGS, axis=1)
        if any(df_sra["type_isolate"]!="clinical"): raise ValueError("all isolates are expected to be clinical")
        df_sra["Sample Name"] = df_sra["Sample Name"].apply(str)
        df_sra["Sample Name"] = df_sra.apply(get_sample_name_compatible_w_Bergin, axis=1)


        # load the Bergin2022_tableStrains
        df_Bergin2022 = pd.read_excel("%s/Bergin2022_tableStrains.xlsx"%ManualDataDir)
        df_Bergin2022 = df_Bergin2022[df_Bergin2022["Species "]=="C. parapsilosis "]
        df_Bergin2022["data_source_BioProject_sets"] = df_Bergin2022['Data Source '].apply(get_data_source_BioProjects_Bergin2022)
        all_bps_Bergin2022 = set.union(*df_Bergin2022.data_source_BioProject_sets)
        df_Bergin2022 = df_Bergin2022.apply(get_city_and_country_collection_Bergin2022, axis=1)
        df_Bergin2022["strain"] = df_Bergin2022["Strain "].apply(lambda x: str(x).replace(" ", ""))
        check_no_nans_series(df_Bergin2022.strain)
        if len(df_Bergin2022)!=len(set(df_Bergin2022.strain)): raise ValueError("non unique strains")

        # define the fields of Oscar that are relevant
        pheno_fiels_Oscar =  ['FLC_resistance', 'VRC_resistance', 'ITR_resistance', 'POS_resistance', 'AMB_MIC', '5FC_MIC', 'FLC_MIC', 'ITR_MIC', 'VRC_MIC', 'POS_MIC', 'IVZ_MIC', 'CAS_MIC', 'MIF_MIC', 'ANI_MIC', 'MIC ETEST 24H', 'MIC ETEST 48H', 'Colonies inside halo', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h', 'Agar invasiveness', 'MORPHOLOGY', 'Microfluidic behavior', 'Biofilm (High / Low / No)', '% VIABILITY ETHANOL', '% VIABILITY HYPLOCHLORITE', '% VIABILITY H202', '% VIABILITY CHLORHEXIDINE', '% VIABILITY ZWITTERGENT 3-14', '% VIABILITY SURFANIO', 'ERG11_mutation']

        clinical_fields_Oscar = ['genotype', 'Referencia cepa Hospital', 'duplicated_patientID', 'date_collection', 'muestra', 'hospital', 'date_birth', 'age', 'Observaciones', 'CCAA']

        # add features to df_metadata
        def get_metadata_row(r):
            if r.source_data=="sequenced_here": 

                # info from Oscar
                for f in (pheno_fiels_Oscar + clinical_fields_Oscar + ["country_collection", "year_collection", "city_collection"]):
                    r[f] = get_formatted_metadata_field(str(df_metadata_Oscar.loc[r.sampleID, f]), f)

                # specific fields
                r["clade_Bergin2022"] = -1
                r["type_isolate"] = "clinical"

            elif r.source_data=="SRA":

                # add empty fields (no data)
                for f in (pheno_fiels_Oscar + clinical_fields_Oscar):
                    r[f] = get_formatted_metadata_field("nan", f)

                # get SRA run
                r_sra = df_sra.loc[r.runID]

                # get the type of sample
                r["type_isolate"] = r_sra.type_isolate

                # get the year_collection
                if not pd.isna(r_sra.year_isolated): r["year_collection"] = int(r_sra.year_isolated)
                elif not pd.isna(r_sra.Collection_Date) and not r_sra.Collection_Date in {"missing", "Missing", "not collected"}: r["year_collection"] = int(r_sra.Collection_Date.split("-")[0])
                else: r["year_collection"] = -1

                # for bergin2022 samples
                if r_sra.BioProject in all_bps_Bergin2022: 

                    # get df of bergin, with the geographic info
                    df = df_Bergin2022[(df_Bergin2022.data_source_BioProject_sets.apply(lambda x: r_sra.BioProject in x))]

                    sample_name_to_berginStrain = {"MSK1206":"MSK1206CP.20477.048"}
                    if r_sra["Sample Name"] in sample_name_to_berginStrain: df = df[df.strain==sample_name_to_berginStrain[r_sra["Sample Name"]]]
                    else: df = df[df.strain.isin({str(r_sra[f]) for f in ['Sample Name', 'Library Name', 'strain']})]

                    # check
                    if len(df)!=1: 
                        print(df)
                        print(sorted(set(df_Bergin2022[df_Bergin2022.data_source_BioProject_sets.apply(lambda x: r_sra.BioProject in x)].strain)))
                        print(r_sra.BioProject)
                        print(r_sra[['Sample Name', 'Library Name', 'strain']])
                        raise ValueError("df should be 1")

                    r["country_collection"] = df.iloc[0]["country_collection"]
                    r["city_collection"] = df.iloc[0]["city_collection"]
                    r["clade_Bergin2022"] = int(df.iloc[0]["Clade "])

                else:

                    # unset the clade
                    r["clade_Bergin2022"] = -1

                    if not pd.isna(r_sra.geo_loc_name): 
                        country, city = {"Germany: Gottingen":("Germany", "Gottingen"), "USA":("USA", "unknown"), "China: Wuhan": ("China", "Wuhan"), "Turkey: Izmir":("Turkey", "Izmir"), "USA: DC":("USA", "Washington DC"), "Malaysia":("Malaysia", "unknown"), "Hungary: Szeged":("Hungary", "Szeged"), 'Puerto Rico': ('Puerto Rico', "unknown"), "Spain":("Spain", "unknown"), "Canada":("Canada", "unknown"), "Portugal":("Portugal", "unknown"), 'Korea: Seoul':('Korea', 'Seoul'), 'not collected':("unknown", "unknown")}[r_sra.geo_loc_name]

                        r["country_collection"] = country
                        r["city_collection"] = city
                    
                    else:
                        r["country_collection"] = "unknown"
                        r["city_collection"] = "unknown"

            else: raise ValueError(r)

            # add the city_and_hospital
            if r.hospital!="unknown" and r.city_collection!="unknown": r["hospital_short"] = "%s (%s)"%(r.hospital.replace("HOSPITAL", "H."), r.city_collection)
            elif r.hospital=="unknown" and r.city_collection!="unknown": r["hospital_short"] = r.city_collection
            elif r.hospital=="unknown" and r.city_collection=="unknown": r["hospital_short"] = "unknown"
            else: raise ValueError("invalid r: %s"%r)

            return r

        df_metadata = df_metadata.apply(get_metadata_row, axis=1)
        check_no_nans_in_df(df_metadata)

        # change some fields
        df_metadata["muestra"] = df_metadata.muestra.apply(lambda x: "_".join(x.split()).upper().replace("Í", "I").replace("Á", "A").replace("Ó", "O"))

        # print the keys extra or missing
        print("Fields not in Oscar", set(df_metadata.keys()).difference(set(keys_df_Oscar)))
        print("Fields not in df_metadata", set(keys_df_Oscar).difference(set(df_metadata.keys())))

        save_object(df_metadata, df_metadata_file)

    return load_object(df_metadata_file)


def get_series_with_missing_runs_to_minus1(series_obj, expected_runs):

    """returns a series with missing runs set to -1"""

    missing_runs = expected_runs.difference(set(series_obj.index))
    return pd.concat([series_obj, pd.Series(-1, index=sorted(missing_runs))])


def get_df_metadata_etest_solid_different_than_liquid_24h(df_metadata):

    """gets a df of samples where the liquid and solid MICs are different"""

    return df_metadata[(df_metadata.FLC_MIC>=16) & (df_metadata.FLC_MIC<=64) & ((df_metadata["MIC ETEST 24H"]==512) | ((df_metadata["MIC ETEST 24H"]>=6) & (df_metadata["MIC ETEST 24H"]<=48)))].copy()


def get_df_metadata_etest_mismatch_24h_48h(df_metadata):

    """Gets a df where the MICs of the 12h and 24h are different"""

    field_24 = "MIC ETEST 24H"
    field_48 = "MIC ETEST 48H"

    return df_metadata[(df_metadata[field_24]>=6.0) & (df_metadata[field_24]<=16.0)].copy() # mildly resistant

def get_type_of_sample_source(sample_type, sources_1, sources_0, sources_minus1):

    """get the type of sample"""

    if sample_type in sources_1: return 1

    elif sample_type in sources_0: return 0

    elif sample_type in sources_minus1: return -1

    else:
        raise ValueError("unconsidered %s"%sample_type)

##### DEFINE PHENO MAPPINGS ######

abb_drug_to_name = {"AMB":"amphotericin B", "5FC":"5-FC", "FLC":"fluconazole", "VRC":"voriconazole", "ITR":"itraconazole", "POS":"posaconazole", "IVZ":"isavuconazole", "CAS":"caspofungin", "MIF":"micafungin", "ANI":"anidulafungin"}

disinfectant_to_name = {"ETHANOL":"ethanol", "HYPLOCHLORITE":"hypochlorite", "H202":"H2O2", "CHLORHEXIDINE":"chlorhexidine", "ZWITTERGENT 3-14":"zwittergent 3-14", "SURFANIO":"surfanio"}


dicotomic_pheno_to_name = {'etest_colony_formation_FLC_MIC_8-32': "fluconazole hetero-R", 
                           "FLC_MIC_higher_in_solid": "fluconazole agar R", 
                           "FLC_MIC_higher_48h": "fluconazole tolerance", 

                           'agar_invasiveness': "agar invasiveness",
                           "pseudohyphae_formation": "pseudohyphae",
                           "biofilm_formation": "biofilm",
                           "microfluidics_Adherent": "microfluidics adherent",
                           "microfluidics_Big_Clumps": "microfluidics big clumps",
                           "microfluidics_Little_clump": "microfluidics little clumps"}

clincial_sources_translated = {
    'source_ABSCESO': 'abscess',
    'source_ABSCESO_ABDOMINAL': 'abdominal abscess',
    'source_ASPIRADO_BRONQUIAL': 'bronchial aspirate',
    'source_ASPIRADO_TRAQUEAL': 'tracheal aspirate',
    'source_BIOPSIA_PULMONAR': 'lung biopsy',
    'source_EOF': 'oropharyngeal sputum',
    'source_ESPUTO': 'sputum',
    'source_EXUDADO_CUTANEO': 'skin exudate',
    'source_EXUDADO_DE_HERIDA': 'wound exudate',
    'source_EXUDADO_OROFARINGEO': 'oropharyngeal exudate',
    'source_EXUDADO_OTICO': 'ear exudate',
    'source_EXUDADO_RECTAL': 'rectal exudate',
    'source_EXUDADO_RESPIRATORIO_ALTO': 'upper respiratory exudate',
    'source_HC': 'blood culture',
    'source_HERIDA': 'wound',
    'source_LAVADO_BRONCOALVEOLAR': 'bronchoalveolar lavage',
    'source_LCR': 'cerebrospinal fluid',
    'source_LIQUIDO_ASCITICO': 'ascitic fluid',
    'source_LIQUIDO_BILIAR': 'bile fluid',
    'source_MUESTRA_AMBIENTAL': 'environmental sample',
    'source_ORINA': 'urine',
    'source_PIEL': 'skin surface',
    'source_RESPIRATORIAS': 'respiratory disease',
    'source_SANGRE_COMPLETA': 'whole blood',
    'source_TEJIDO_VASCULAR': 'vascular tissue',
    'source_VALVULA_CARDIACA': 'heart valve'
}

##################################



def get_paper_dicotomic_pheno(pheno):

    """Gets a dicotomic pheno and returns the paper version"""

    # get
    if pheno.endswith("_resistance"):
        return "%s R"%(abb_drug_to_name[pheno.split("_")[0]])

    elif pheno.endswith("_tolerance"):
        return "%s tol"%(disinfectant_to_name[pheno.split("_")[0]])

    elif pheno in dicotomic_pheno_to_name:
        return dicotomic_pheno_to_name[pheno]

    elif pheno in clincial_sources_translated:
        return clincial_sources_translated[pheno]

    elif pheno.startswith("source_") and not pheno in clincial_sources_translated:
        return " ".join(pheno.split("_")[1:])

    elif pheno.startswith("age_"):
        return pheno.replace("_", " ")

    elif pheno in {'Barcelona - Bellvitge', 'Barcelona - J. Deu', 'Barcelona - V. Hebron', 'Burgos', 'Lab. Catalunya', 'Madrid - 12 Oct', 'Madrid - C. Defensa', 'Madrid - Mostoles', 'Madrid - P. Hierro', 'Madrid - S. Carlos', 'Santander', 'Sevilla - S. Aljarafe', 'Sevilla - V. Macarena', 'Sevilla - V. Rocio', 'Valladolid'}:
        return pheno

    elif pheno.startswith("ERG11 - "):
        return pheno

    else:
        raise ValueError("invalid %s"%pheno)

def Fig_PhenoDefinition_continuous_phenotypes(df_metadata, df_metadata_binary_phenotypes, PlotsDir):

    """Plots the distribution of continuous phenos"""

    continuous_phenos = [ 'AMB_resistance', '5FC_resistance', 'FLC_resistance', 'ITR_resistance', 'VRC_resistance', 'POS_resistance','IVZ_resistance', 'CAS_resistance', 'MIF_resistance', 'ANI_resistance', 'ETHANOL_tolerance', 'HYPLOCHLORITE_tolerance', 'H202_tolerance', 'CHLORHEXIDINE_tolerance', 'ZWITTERGENT 3-14_tolerance', 'SURFANIO_tolerance']

    # init fig
    nrows = 4
    ncols = 5
    fig = plt.figure(figsize=(ncols*2.3, nrows*2))

    for I, pheno in enumerate(continuous_phenos):

        # define the continuous pheno field
        if pheno.endswith("_resistance"):
            continuous_f = pheno.split("_")[0] + "_MIC"
            continuous_name = "MIC " + abb_drug_to_name[pheno.split("_")[0]]

        elif pheno.endswith("_tolerance"):
            continuous_f = '% VIABILITY ' + pheno.split("_")[0]
            continuous_name = "% viability " + disinfectant_to_name[pheno.split("_")[0]]

        else:
            raise ValueError("invalid %s"%pheno)

        # get plot
        ax = plt.subplot(nrows, ncols, I+1)
        df_plot = df_metadata[df_metadata[continuous_f]!=-1].copy().sort_values(by=[continuous_f], ascending=True)[[continuous_f, "runID"]].copy()

        print(pheno, len(df_plot))

        # add resistance
        df_plot = df_plot.merge(df_metadata_binary_phenotypes[["runID", pheno]], on="runID", validate="one_to_one", how="left")
        check_no_nans_in_df(df_plot)

        # get bp tuple
        bp_tuple = get_bps_for_pheno_field(continuous_f)

        # define the type of plot
        palette = {0:"gray", -1:"white", 1:"black"}
        if "MIC" in continuous_f:

            sorted_vals = sorted(set(df_plot[continuous_f]))
            df_plot["str_x"] = df_plot[continuous_f].apply(str)

            ax = sns.histplot(df_plot, x="str_x", hue=pheno,  multiple='stack', palette=palette, alpha=.8)
            ax.set_xticklabels(sorted_vals, rotation=90)

            # get x
            for x, val in enumerate(sorted_vals):

                if val>bp_tuple[0] and val<bp_tuple[1]:
                    ax.axvline(x, c="red", lw=.8, ls="--")



        else: 
            ax = sns.histplot(data=df_plot, x=continuous_f, hue=pheno, palette=palette, alpha=.8, bins=10, multiple='stack')

            for bp in bp_tuple:
                ax.axvline(bp, c="red", lw=.8, ls="--")

        # add the labels
        ax.get_legend().remove()
        ax.set_xlabel(continuous_name)

        if I in {0, 5, 10, 15}:
            ax.set_ylabel("# isolates")

        else:
            ax.set_ylabel("")

    plt.subplots_adjust(wspace=0.2, hspace=0.85)# wspace (between cols) hspaces (btw rows)

    sns.despine()
    save_fig(fig, "Distribution_continuous_phenotypes.pdf", PlotsDir)

def Fig_PhenoDefinition_colony_formation(df_metadata, df_metadata_binary_phenotypes, PlotsDir):

    """Gets the colony formation plots"""

    # keep fields in whcich this was measured
    df_metadata = df_metadata[(df_metadata['Colonies inside halo']!="unknown")][["Colonies inside halo", "FLC_MIC", "runID"]].copy()

    # add binaries
    df_metadata = df_metadata.merge(df_metadata_binary_phenotypes[["runID", "etest_colony_formation_FLC_MIC_8-32"]], on="runID", validate="one_to_one", how="left")

    print(sum((df_metadata.FLC_MIC>=8) & (df_metadata.FLC_MIC<=32)), "strains MIC 8-32")
    print(len(df_metadata), "strains Etest")

    fig = plt.figure(figsize=(1.8, 2))
    for y in [8, 32]:
        plt.axhline(y, c="red", lw=.8, ls="--")

    palette = {0:"gray", -1:"white", 1:"black"}
    ax = sns.swarmplot(data=df_metadata, x="Colonies inside halo", y="FLC_MIC", hue="etest_colony_formation_FLC_MIC_8-32", palette=palette, linewidth=.7)
    #ax = sns.stripplot(data=df_metadata, x="Colonies inside halo", y="FLC_MIC", hue="etest_colony_formation_FLC_MIC_8-32", palette=palette, linewidth=.7, jitter=0.25)

    ax.set_yscale("log")
    ax.get_legend().remove()
    ax.set_ylabel("MIC fluconazole")
    ax.set_xlabel("Etest colonies\nwithin halo")

    sns.despine()
    save_fig(fig, "Etest_colonies_halo.pdf", PlotsDir)

def fig_PhenoDefinition_higher_R_solid(df_metadata, df_metadata_binary_phenotypes, PlotsDir):

    """Gets the higher in solid case"""

    # get good strains
    df_metadata = df_metadata[df_metadata["MIC ETEST 24H"]!=-1][["MIC ETEST 24H", "runID", "FLC_MIC", 'Colonies inside halo']].copy()

    # add phenos
    df_metadata = df_metadata.merge(df_metadata_binary_phenotypes[["runID", "FLC_MIC_higher_in_solid"]], on="runID", validate="one_to_one", how="left")

    print(len(df_metadata), "strains Etest")
    check_no_nans_in_df(df_metadata)

    fig = plt.figure(figsize=(2, 2))

    # plot
    palette = {0:"gray", -1:"white", 1:"black"}
    ax = sns.scatterplot(data=df_metadata, x="FLC_MIC", y="MIC ETEST 24H", hue="FLC_MIC_higher_in_solid", palette=palette, linewidth=.7, edgecolor="black")

    for x in [16, 64]:
        plt.axvline(x, c="red", lw=.8, ls="--", zorder=0)

    for y in [6, 48]:
        plt.axhline(y, c="red", lw=.8, ls="--", zorder=0)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.get_legend().remove()
    ax.set_xlabel("MIC fluconazole liquid")
    ax.set_ylabel("MIC fluconazole Etest 24h  ")
    sns.despine()

    save_fig(fig, "Etest_higher_solid_R.pdf", PlotsDir)


def fig_PhenoDefinition_agar_tolerance_48h(df_metadata, df_metadata_binary_phenotypes, PlotsDir):

    """Gets the higher in solid case"""

    # get good strains
    df_metadata = df_metadata[(df_metadata["MIC ETEST 24H"]!=-1) & (df_metadata["MIC ETEST 48H"]!=-1)][["MIC ETEST 24H", "MIC ETEST 48H", "runID", "FLC_MIC", 'Colonies inside halo']].copy()

    # add phenos
    df_metadata = df_metadata.merge(df_metadata_binary_phenotypes[["runID", "FLC_MIC_higher_48h"]], on="runID", validate="one_to_one", how="left")

    print(len(df_metadata), "strains Etest 24 and 48")
    check_no_nans_in_df(df_metadata)

    fig = plt.figure(figsize=(2, 2))

    # plot
    palette = {0:"gray", -1:"white", 1:"black"}
    ax = sns.scatterplot(data=df_metadata, x="MIC ETEST 24H", y="MIC ETEST 48H", hue="FLC_MIC_higher_48h", palette=palette, linewidth=.7, edgecolor="black")

    for x in [6, 16]:
        plt.axvline(x, c="red", lw=.8, ls="--", zorder=0)



    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.get_legend().remove()
    ax.set_xlabel("MIC fluconazole Etest 24h  ")
    ax.set_ylabel("MIC fluconazole Etest 48h  ")
    sns.despine()

    save_fig(fig, "Etest_48h.pdf", PlotsDir)


def fig_PhenoDefinition_age(df_metadata_filt, df_metadata_binary_phenotypes, PlotsDir):

    """Defines the ages"""

    df_metadata = df_metadata_filt[df_metadata_filt.age!=124].copy()
    print("len df ages", len(df_metadata))

    # fig
    fig = plt.figure(figsize=(2.5, 5))


    control_age_thresholds = [60, 75]
    for Iage, (type_age, age_range) in enumerate([("infant", (0,10)), ("young", (20,50)), ("old", (80, 110))]):

        ax = plt.subplot(3, 1, Iage+1)
        age_f = "age_"+type_age
        df_plot = df_metadata[["runID", "age"]].merge(df_metadata_binary_phenotypes[["runID", age_f]], on="runID", how="left", validate="one_to_one")
        
        palette = {0:"darkgray", -1:"white", 1:"black"}
        ax = sns.histplot(data=df_plot, x="age", hue=age_f, palette=palette, alpha=1, bins=30, multiple='stack')
        ax.legend().remove()

        #ax.text(3, 25, "age "+type_age)
        ax.set_title("age "+type_age, fontsize=10)


        for x in (set(control_age_thresholds).union(set(age_range))).difference({0, 110}):
            ax.axvline(x, c="red", lw=.8, ls="--")

        if Iage==2:
            ax.set_xlabel("Patient age (years)")

        else:
            ax.set_xticklabels([])
            ax.set_xlabel("")

        if Iage==1:
            ax.set_ylabel("Number of isolates")

        else:
            ax.set_ylabel("")

    sns.despine()
    plt.subplots_adjust(wspace=0, hspace=0.4)# wspace (between cols) hspaces (btw rows)

    save_fig(fig, "age_distributions.pdf", PlotsDir)



def get_df_binary_phenotypes_validation(df_phenotypes_validation):

    """gets the binary phenotypes of the validation dataset"""
    print("get binary phenotypes validation...")

    # keep
    df_phenotypes_validation = df_phenotypes_validation.copy().set_index("Run", drop=False)

    # init the df_phenos 
    df_phenos = df_phenotypes_validation[["Run"]].copy().rename(columns={"Run":"runID"})

    # add the MIC-based phenos
    MIC_fields = [k for k in df_phenotypes_validation.keys() if k.endswith("_MIC")]

    for mic_f in MIC_fields:
        drug = mic_f.split("_")[0]
        runID_to_R = df_phenotypes_validation[mic_f].apply(get_dicothomic_resistance, pheno=mic_f).map({"R":1, "S":0, "I":-1, "unknown":-1})
        df_phenos["%s_resistance"%drug] = df_phenos.runID.map(runID_to_R)
        check_no_nans_series(df_phenos["%s_resistance"%drug])

    return df_phenos


def get_df_metadata_binary_phenotypes_curated(df_metadata, ProcessedDataDir):

    """Gets the df with all the binary phenotypes 1/0 curated"""

    df_metadata_binary_phenos_file = "%s/df_metadata_binary_phenotypes_curated.py"%ProcessedDataDir
    print("getting curated metadata")

    if file_is_empty(df_metadata_binary_phenos_file):

        # filter the df
        samples_discard_run1 = {"OSCAR01_L21-060", "OSCAR01_L21-067"}
        df_metadata = df_metadata[(df_metadata.possibly_triploid==False) & ~(df_metadata.runID.isin(samples_discard_run1)) & (df_metadata.source_data=="sequenced_here")].copy()
        all_runs = set(df_metadata.runID)

        # checks
        if len(df_metadata)!=len(set(df_metadata.sampleID)): raise ValueError("sID should be unique")

        # init the df metadata 
        df_phenos = df_metadata[["sampleID", "runID"]].copy()

        # add the MIC-based phenos
        for mic_f in ['AMB_MIC', '5FC_MIC', 'FLC_MIC', 'ITR_MIC', 'VRC_MIC', 'POS_MIC', 'IVZ_MIC', 'CAS_MIC', 'MIF_MIC', 'ANI_MIC']:
            drug = mic_f.split("_")[0]
            runID_to_R = df_metadata[mic_f].apply(get_dicothomic_resistance, pheno=mic_f).map({"R":1, "S":0, "I":-1, "unknown":-1})
            df_phenos["%s_resistance"%drug] = df_phenos.runID.map(runID_to_R)
            check_no_nans_series(df_phenos["%s_resistance"%drug])

        # add the disinfectant-based phenos
        for disinf in ["ETHANOL", "HYPLOCHLORITE", "H202", "CHLORHEXIDINE", "ZWITTERGENT 3-14", "SURFANIO"]:

            pheno_f = f'% VIABILITY {disinf}'
            runID_to_R = df_metadata[pheno_f].apply(get_dicothomic_resistance, pheno=pheno_f).map({"R":1, "S":0, "I":-1, "unknown":-1})
            df_phenos["%s_tolerance"%disinf] = df_phenos.runID.map(runID_to_R)

        check_no_nans_in_df(df_phenos)

        # add whether it forms or not colonies (only samples with FLC_MIC 8-32, as in these ones it makes sense)
        runID_to_colony_formation = df_metadata[(df_metadata['Colonies inside halo']!="unknown") & (df_metadata.FLC_MIC>=8) & (df_metadata.FLC_MIC<=32)].copy()['Colonies inside halo'].map({"yes":1, "no":0})

        runID_to_colony_formation = get_series_with_missing_runs_to_minus1(runID_to_colony_formation, all_runs)
        df_phenos["etest_colony_formation_FLC_MIC_8-32"] = df_phenos.runID.map(runID_to_colony_formation)

        # add whether it has a missmatch between liquid and solid resistance (only in samples with FLC_MIC 16-64)
        bool_dict = {True:1, False:0}
        runID_to_higherMIC_solid = (get_df_metadata_etest_solid_different_than_liquid_24h(df_metadata)["MIC ETEST 24H"]==512).map(bool_dict)
        runID_to_higherMIC_solid = get_series_with_missing_runs_to_minus1(runID_to_higherMIC_solid, all_runs)
        df_phenos["FLC_MIC_higher_in_solid"] = df_phenos.runID.map(runID_to_higherMIC_solid)

        # add whether there is a missmatch between 24h and 48h
        runID_to_higher_MIC_48h = (get_df_metadata_etest_mismatch_24h_48h(df_metadata)["MIC ETEST 48H"]==512).map(bool_dict)
        runID_to_higher_MIC_48h = get_series_with_missing_runs_to_minus1(runID_to_higher_MIC_48h, all_runs)
        df_phenos["FLC_MIC_higher_48h"] = df_phenos.runID.map(runID_to_higher_MIC_48h)

        # add the age groups
        df_metadata_age = df_metadata[df_metadata.age!=124].copy()
        control_age_runs = set(df_metadata_age[(df_metadata_age.age>=60) & (df_metadata_age.age<=75)].runID)

        for type_age, age_range in [("infant", (0,10)), ("young", (20,50)), ("old", (80, 110))]:
            age_runs = set(df_metadata_age[(df_metadata_age.age>=age_range[0]) & (df_metadata_age.age<=age_range[1])].runID)

            runID_to_age_g = pd.Series(1, index=sorted(age_runs))
            runID_to_age_g = pd.concat([runID_to_age_g, pd.Series(0, index=sorted(control_age_runs))])
            runID_to_age_g = get_series_with_missing_runs_to_minus1(runID_to_age_g, all_runs)

            if len(runID_to_age_g)!=len(set(runID_to_age_g.index)): raise ValueError("index should be unique")
            df_phenos["age_%s"%type_age] = df_phenos.runID.map(runID_to_age_g)


        # add dicotomic phenos
        df_phenos["agar_invasiveness"] = df_phenos.runID.map(df_metadata["Agar invasiveness"]).map({"Yes":1, "No":0, "Low":-1, 'unknown':-1})
        df_phenos["pseudohyphae_formation"] = df_phenos.runID.map(df_metadata["MORPHOLOGY"]).map({'Psedohyphae':1, 'Yeast':0, 'unknown':-1})
        df_phenos["biofilm_formation"] = df_phenos.runID.map(df_metadata["Biofilm (High / Low / No)"]).map({'HF':1, 'NF':0, 'LF':-1, 'unknown':-1})

        # microfluidic behavior. There are three distinct cathegories 'Adherent', 'Big Clumps', 'Little clump', not ordered
        all_microfluidic_behaviors = {'Adherent', 'Big Clumps', 'Little clump'}
        runs_w_microfluidic_behaviors = set(df_metadata[df_metadata['Microfluidic behavior']!="unknown"].runID)

        for b in all_microfluidic_behaviors:

            runs_behavior = set(df_metadata[df_metadata['Microfluidic behavior']==b].runID)
            runs_not_behavior = runs_w_microfluidic_behaviors.difference(runs_behavior)
            if len(runs_behavior)==0: raise ValueError("invalid b")

            runID_to_behaviorID = pd.concat([pd.Series(1, index=sorted(runs_behavior)), pd.Series(0, index=sorted(runs_not_behavior))])
            runID_to_behaviorID = get_series_with_missing_runs_to_minus1(runID_to_behaviorID, all_runs)
            df_phenos["microfluidics_%s"%("_".join(b.split()))] = df_phenos.runID.map(runID_to_behaviorID) 

        # add the sample type, one for each type
        for sample_type in sorted(set(df_metadata.muestra)):
            df_phenos["source_%s"%sample_type] = df_phenos.runID.map(df_metadata.muestra==sample_type).map(bool_dict)
            #print(sample_type, sum(df_phenos["source_%s"%sample_type]))

        # add different types of sample types depending on source
        for field in ["source_sterile", "source_sterile_vs_environmental", "source_blood", "source_blood_vs_environmental", "source_blood_vs_other_sterile", "source_skin", "source_respiratory"]:
            print(field)

            # sterile compartments vs all the other samples
            if field=="source_sterile":
                sources_1 = {'HC', 'LCR', 'LIQUIDO_ASCITICO', 'LIQUIDO_BILIAR', 'SANGRE_COMPLETA', 'TEJIDO_VASCULAR', 'VALVULA_CARDIACA'}
                sources_0 = {'ABSCESO', 'ABSCESO_ABDOMINAL', 'ASPIRADO_BRONQUIAL', 'ASPIRADO_TRAQUEAL', 'BIOPSIA_PULMONAR', 'EOF', 'ESPUTO', 'EXUDADO_CUTANEO', 'EXUDADO_DE_HERIDA', 'EXUDADO_OROFARINGEO', 'EXUDADO_OTICO', 'EXUDADO_RECTAL', 'EXUDADO_RESPIRATORIO_ALTO', 'HERIDA', 'LAVADO_BRONCOALVEOLAR', 'MUESTRA_AMBIENTAL', 'ORINA', 'PIEL', 'RESPIRATORIAS'}
                sources_minus1 = set()

            # sterile compartments vs environmental isolates
            elif field=="source_sterile_vs_environmental":
                sources_1 = {'HC', 'LCR', 'LIQUIDO_ASCITICO', 'LIQUIDO_BILIAR', 'SANGRE_COMPLETA', 'TEJIDO_VASCULAR', 'VALVULA_CARDIACA'}
                sources_0 = {'MUESTRA_AMBIENTAL'}
                sources_minus1 = {'ABSCESO', 'ABSCESO_ABDOMINAL', 'ASPIRADO_BRONQUIAL', 'ASPIRADO_TRAQUEAL', 'BIOPSIA_PULMONAR', 'EOF', 'ESPUTO', 'EXUDADO_CUTANEO', 'EXUDADO_DE_HERIDA', 'EXUDADO_OROFARINGEO', 'EXUDADO_OTICO', 'EXUDADO_RECTAL', 'EXUDADO_RESPIRATORIO_ALTO', 'HERIDA', 'LAVADO_BRONCOALVEOLAR', 'ORINA', 'PIEL', 'RESPIRATORIAS'}

            # from blood vs all other samples
            elif field=="source_blood":
                sources_1 = {'HC', 'SANGRE_COMPLETA', 'TEJIDO_VASCULAR'}
                sources_0 = {'MUESTRA_AMBIENTAL', 'ABSCESO', 'ABSCESO_ABDOMINAL', 'ASPIRADO_BRONQUIAL', 'ASPIRADO_TRAQUEAL', 'BIOPSIA_PULMONAR', 'EOF', 'ESPUTO', 'EXUDADO_CUTANEO', 'EXUDADO_DE_HERIDA', 'EXUDADO_OROFARINGEO', 'EXUDADO_OTICO', 'EXUDADO_RECTAL', 'EXUDADO_RESPIRATORIO_ALTO', 'HERIDA', 'LAVADO_BRONCOALVEOLAR', 'LCR', 'LIQUIDO_ASCITICO', 'LIQUIDO_BILIAR', 'ORINA', 'PIEL', 'RESPIRATORIAS', 'VALVULA_CARDIACA'}
                sources_minus1 = set()

            # from blood vs environmental isolates
            elif field=="source_blood_vs_environmental":
                sources_1 = {'HC', 'SANGRE_COMPLETA', 'TEJIDO_VASCULAR'}
                sources_0 = {'MUESTRA_AMBIENTAL'}
                sources_minus1 = {'ABSCESO', 'ABSCESO_ABDOMINAL', 'ASPIRADO_BRONQUIAL', 'ASPIRADO_TRAQUEAL', 'BIOPSIA_PULMONAR', 'EOF', 'ESPUTO', 'EXUDADO_CUTANEO', 'EXUDADO_DE_HERIDA', 'EXUDADO_OROFARINGEO', 'EXUDADO_OTICO', 'EXUDADO_RECTAL', 'EXUDADO_RESPIRATORIO_ALTO', 'HERIDA', 'LAVADO_BRONCOALVEOLAR', 'LCR', 'LIQUIDO_ASCITICO', 'LIQUIDO_BILIAR', 'ORINA', 'PIEL', 'RESPIRATORIAS', 'VALVULA_CARDIACA'}

            # from blood vs other sterile compartments
            elif field=="source_blood_vs_other_sterile":
                sources_1 = {'HC', 'SANGRE_COMPLETA', 'TEJIDO_VASCULAR'}
                sources_0 = {'LCR', 'LIQUIDO_ASCITICO', 'LIQUIDO_BILIAR'}
                sources_minus1 = {'ABSCESO', 'ABSCESO_ABDOMINAL', 'ASPIRADO_BRONQUIAL', 'ASPIRADO_TRAQUEAL', 'BIOPSIA_PULMONAR', 'EOF', 'ESPUTO', 'EXUDADO_CUTANEO', 'EXUDADO_DE_HERIDA', 'EXUDADO_OROFARINGEO', 'EXUDADO_OTICO', 'EXUDADO_RECTAL', 'EXUDADO_RESPIRATORIO_ALTO', 'HERIDA', 'LAVADO_BRONCOALVEOLAR', 'MUESTRA_AMBIENTAL', 'ORINA', 'PIEL', 'RESPIRATORIAS', 'VALVULA_CARDIACA'}

            # from the skin vs non other non-blood compartments. I discard blood because it could come from skin. I discard environmental because it could be from skin
            elif field=="source_skin":
                sources_1 = {'EXUDADO_CUTANEO', 'EXUDADO_DE_HERIDA', 'HERIDA', 'PIEL'}
                sources_0 = {'ABSCESO', 'ABSCESO_ABDOMINAL', 'ASPIRADO_BRONQUIAL', 'ASPIRADO_TRAQUEAL', 'BIOPSIA_PULMONAR', 'EOF', 'ESPUTO', 'EXUDADO_OROFARINGEO', 'EXUDADO_OTICO', 'EXUDADO_RECTAL', 'EXUDADO_RESPIRATORIO_ALTO', 'LAVADO_BRONCOALVEOLAR', 'LCR', 'LIQUIDO_ASCITICO', 'LIQUIDO_BILIAR', 'ORINA', 'RESPIRATORIAS', 'VALVULA_CARDIACA'}
                sources_minus1 = {'HC', 'SANGRE_COMPLETA', 'TEJIDO_VASCULAR', 'MUESTRA_AMBIENTAL'}

            # from the respiratory system. I discard blood because it could come from respiratory
            elif field=="source_respiratory":
                sources_1 = {'ASPIRADO_BRONQUIAL', 'ASPIRADO_TRAQUEAL', 'BIOPSIA_PULMONAR', 'EOF', 'ESPUTO', 'EXUDADO_OROFARINGEO', 'EXUDADO_RESPIRATORIO_ALTO', 'LAVADO_BRONCOALVEOLAR', 'RESPIRATORIAS'}
                sources_0 = {'ABSCESO', 'ABSCESO_ABDOMINAL', 'EXUDADO_CUTANEO', 'EXUDADO_DE_HERIDA', 'EXUDADO_OTICO', 'EXUDADO_RECTAL',  'HERIDA', 'LCR', 'LIQUIDO_ASCITICO', 'LIQUIDO_BILIAR', 'MUESTRA_AMBIENTAL', 'ORINA', 'PIEL', 'VALVULA_CARDIACA'}
                sources_minus1 = {'HC', 'SANGRE_COMPLETA', 'TEJIDO_VASCULAR'}

            else: 
                raise ValueError("invalid source")

            """
            elif field=="source_":
                sources_1 = {}
                sources_0 = {}
                sources_minus1 = {}
            """

            # checks
            if len(sources_1.intersection(sources_0))>0: raise ValueError("samples in 1 and 0")
            if len(sources_1.intersection(sources_minus1))>0: raise ValueError("samples in 1 and -1")
            if len(sources_0.intersection(sources_minus1))>0: raise ValueError("samples in 0 and -1")

            sources_considered = sources_1.union(sources_0).union(sources_minus1)
            if sources_considered!=set(df_metadata.muestra): 

                print("considered, not in muestra", sources_considered.difference(set(df_metadata.muestra)))
                print(" in muestra, not considered", set(df_metadata.muestra).difference(sources_considered))
                raise ValueError("not all sources considered")

            # add to df
            df_phenos[field] = df_phenos.runID.map(df_metadata.muestra.apply(lambda x: get_type_of_sample_source(x, sources_1, sources_0, sources_minus1)))
        
        # print all the pheno numbers
        for f in df_phenos.keys():
            if f in {"runID", "sampleID"}: continue
            print(f, "%i/%i 1/0"%(sum(df_phenos[f]==1), sum(df_phenos[f]==0)))

        # checks
        check_no_nans_in_df(df_phenos)

        # checks and save
        save_object(df_phenos, df_metadata_binary_phenos_file)

    return load_object(df_metadata_binary_phenos_file)

def get_formatted_date_collection_tuple(date_col):

    """formats the date collection"""

    # split
    tuple_ints = list(map(int, re.split("\*|-", date_col.split()[0])))

    # define
    if tuple_ints[-1]>2000:
        day, month, year = tuple_ints

    elif tuple_ints[0]>2000:
        year, month, day = tuple_ints

    # checks
    if not year in {2019, 2020, 2021, 2022}: raise ValueError("invalid y: %i"%year)
    if month>12: raise ValueError("invalid month: %i"%month)
    if day>31: raise ValueError("invalid day: %i"%day)

    return (year, month, day)

def analyze_patients_changing_phenotypes(df_phenos, PlotsDir, trees_df, df_metadata_filt):

    """checks if there are some patients w/ duplicated IDs that have different phenotypes"""

    # keep
    df_phenos = df_phenos.copy()

    # add fields
    tree_series = trees_df[(trees_df.run_set=="exclude_triploids") & (trees_df.mode_tree=="diploid") & (trees_df.type_outgroup=="Cparapsilosis") & (trees_df.aligner=="bwa_mem")].iloc[0]
    strains_tree = cp.deepcopy(tree_series.tree_corrected)
    leaf_to_clade_Bergin2022 = tree_series.leaf_to_clade_Bergin2022
    df_phenos["cladeID"] = df_phenos.runID.map(leaf_to_clade_Bergin2022)
    df_phenos["duplicated_patientID"] = df_phenos.runID.map(df_metadata_filt.duplicated_patientID)
    df_phenos["date_collection"] = df_phenos.runID.map(df_metadata_filt.date_collection)
    df_phenos["source"] = df_phenos.runID.map(df_metadata_filt.muestra)


    check_no_nans_in_df(df_phenos)

    # filter for serial isolates
    df_phenos = df_phenos[df_phenos.duplicated_patientID!=-1].copy()
    pID_to_nSamples = df_phenos.groupby("duplicated_patientID").apply(len)
    serial_isolates_pIDs = set(pID_to_nSamples[pID_to_nSamples>1].index)
    df_phenos = df_phenos[df_phenos.duplicated_patientID.isin(serial_isolates_pIDs)].copy()

    # sort by date
    df_phenos["date_collection_tuple"] = df_phenos.date_collection.apply(get_formatted_date_collection_tuple)
    df_phenos = df_phenos.sort_values(by=["duplicated_patientID", "date_collection_tuple"])


    # go through serial isolates that have different phenotypes and print them
    all_phenos = sorted(set(df_phenos.keys()).difference({"cladeID", "duplicated_patientID", "runID", "sampleID", "date_collection", "date_collection_tuple", "source"}))

    for pID in sorted(set(df_phenos.duplicated_patientID)):
        df_p = df_phenos[df_phenos.duplicated_patientID==pID].set_index("runID")

        for pheno in all_phenos:
            if set(df_p[(df_p[pheno]!=-1)][pheno])=={0, 1}:

                fields = ["duplicated_patientID", "date_collection", pheno, "cladeID"]
                if pheno.startswith("source_"): fields.append("source")

                print(df_p[fields]) # runID





def plot_age_distribution(df_metadata, PlotsDir):


    """get age distribution"""

    df_metadata = df_metadata[df_metadata.age!=124].copy()

    fig = plt.figure(figsize=(5,3))
    ax = sns.histplot(df_metadata, x="age", bins=20)

    for x in [10, 20, 50, 60, 75, 80]: plt.axvline(x, ls="--", lw=.7, c="k")
    save_fig(fig, "age_distribution.pdf", PlotsDir)

def get_trees_df_row_with_Bergin2022_clade_info(r, df_metadata):

    """adds to the row info about the clade of Bergin"""

    # get the leaf to clade basic
    leaf_to_clade = get_leaf_to_clade_based_on_Bergin2022(df_metadata, r.tree_corrected)

    # add the fraction of assigned clades
    if len(leaf_to_clade)>0:
        r["fraction_leaves_w_clade_Bergin2022"] = sum([c!=-1 for c in leaf_to_clade.values()]) / len(leaf_to_clade)

    else: 
        r["fraction_leaves_w_clade_Bergin2022"] = 0.0

    # add
    r["leaf_to_clade_Bergin2022"] = leaf_to_clade
    return r

def get_leaf_to_clade_Bergin2022_added_clades(leaf_to_clade, tree, df_metadata, samples_discard_run1):

    """Gets clade names"""

    # keep
    tree = cp.deepcopy(tree)

    # add the sister clades
    leaf_to_clade_final = {}
    for leaf, clade in leaf_to_clade.items():
        if not leaf in tree.get_leaf_names(): continue

        # define new clades
        if clade==-1:

            # get the closest clade
            for n in ([(tree&leaf)] + (tree&leaf).get_ancestors()):
                clades_ancestors = list(set(map(lambda l: leaf_to_clade[l], n.get_leaf_names())).difference({-1}))

                if len(clades_ancestors)==1:
                    closest_clade = clades_ancestors[0]
                    break

                elif len(clades_ancestors)>1:
                    raise ValueError("can't be >1 clade. %s. %s"%(clades_ancestors, leaf))

            # for clade 1, there is only 1
            if closest_clade==1:
                clade_final = 1.1

            elif closest_clade==4:

                # define the clade node
                for n in (tree&leaf).get_ancestors():
                    parent_n = n.get_ancestors()[0]
                    parent_clades = set(map(lambda l: leaf_to_clade[l], parent_n.get_leaf_names()))
                    if parent_clades=={-1, 4}: break

                # define the distance to root, and based on the distance, get the clade
                distance_to_root = tree.get_distance(tree, n, topology_only=True)
                clade_final = {3.0: 4.1, 5.0: 4.2, 4.0: 4.3}[distance_to_root]


            else:
                raise ValueError("clades should be 1 and 4")


        # same clade
        else:
            clade_final = clade

        # keep
        leaf_to_clade_final[leaf] = clade_final

    return leaf_to_clade_final

def get_leaf_to_clade_based_on_Bergin2022(df_metadata, strains_tree):

    """gets a dict that maps each leaf of the tree to the clade ID"""

    strains_tree = cp.deepcopy(strains_tree)
    df_metadata = df_metadata.copy()

    # add cladeID base on clade_Bergin2022. Assumption that Bergin2020 clades are monophyletic (or else return the empty df)
    leaf_to_clade_based_on_Bergin2022 = {}
    all_clades_Bergin2022 = set(df_metadata[df_metadata.clade_Bergin2022!=-1].clade_Bergin2022)
    for clade in sorted(all_clades_Bergin2022):
        other_clades = all_clades_Bergin2022.difference({clade})

        # get ancestor
        leafs_clade = sorted(df_metadata[df_metadata.clade_Bergin2022==clade].runID)
        first_common_ancestor = strains_tree.get_common_ancestor(leafs_clade)
        leafs_ancestor = set(first_common_ancestor.get_leaf_names())

        # check that the topology is the same, or return an empty df
        if any((df_metadata.runID.isin(leafs_ancestor)) & (df_metadata.clade_Bergin2022.isin(other_clades))): 
            return {}

        # keep
        for leaf in leafs_ancestor: leaf_to_clade_based_on_Bergin2022[leaf] = clade

    missing_leafs = set(df_metadata.runID).difference(set(leaf_to_clade_based_on_Bergin2022))
    for leaf in missing_leafs: leaf_to_clade_based_on_Bergin2022[leaf] = -1

    return leaf_to_clade_based_on_Bergin2022


def run_jobarray_file_MN4_greasy(jobs_filename, name, time="12:00:00", queue="bsc_ls", threads_per_job=4, nodes=1, submit=True):

    """
    This function takes a jobs filename and creates a jobscript with args (which is a list of lines to be written to the jobs cript). It works in MN4 for greasy    

    """

    # define dirs
    outdir = get_dir(jobs_filename)
    stddir = "%s/STDfiles"%outdir; 
    delete_folder(stddir); make_folder(stddir)
    os.chdir(outdir)

    # define the std files
    greasy_logfile = "%s/%s_greasy.log"%(stddir, name)
    stderr_file = "%s/%s_stderr.txt"%(stddir, name)
    stdout_file = "%s/%s_stdout.txt"%(stddir, name)

    # define the job script
    jobs_filename_run = "%s.run"%jobs_filename

    # define the number of jobs in the job array
    njobs = int("".join([x for x in str(subprocess.check_output("wc -l %s"%jobs_filename, shell=True)).split()[0] if x.isdigit()])) + 1

    # define the requested nodes. Each node has 48 threads
    total_nthreads = njobs*threads_per_job
    if total_nthreads%48==0: max_nodes = max([int((total_nthreads)/48), 1])
    else: max_nodes = max([int((total_nthreads)/48)+1, 1])
    requested_nodes = min([nodes, max_nodes])

    # define the number of tasks
    max_ntasks = int((requested_nodes*48)/threads_per_job)
    ntasks = min([njobs, max_ntasks])


    # define the arguments
    arguments = [ "#!/bin/sh",
                  "#SBATCH --error=%s"%stderr_file,
                  "#SBATCH --output=%s"%stdout_file,
                  "#SBATCH --job-name=%s"%name, 
                  "#SBATCH --get-user-env",
                  #"#SBATCH --workdir=%s"%outdir,
                  "#SBATCH --time=%s"%time,
                  "#SBATCH --qos=%s"%queue,
                  "#SBATCH --cpus-per-task=%i"%threads_per_job,
                  "#SBATCH --ntasks=%i"%ntasks,
                  "",
                  "module load greasy",
                  "export GREASY_LOGFILE=%s;"%(greasy_logfile),
                  "rm %s;"%greasy_logfile,
                  "echo 'running pipeline';",
                  "greasy %s"%jobs_filename
                ]

    # replace by the run in MN
    arguments = [a.replace(ParentDir, "/gpfs/projects/bsc40/mschikora") for a in arguments]

    # define and write the run filename
    with open(jobs_filename_run, "w") as fd: fd.write("\n".join(arguments))
    
    # run in cluster if specified
    if submit is True: run_cmd("sbatch %s"%jobs_filename_run)


def get_downsampled_reads(input_reads, output_reads, nreads):

    """Downsamples the reads"""

    if file_is_empty(output_reads):
        print("Generating %s"%output_reads)

        # define tmp
        output_reads_uncompressed = "%s.tmp.fastq"%output_reads
        output_reads_compressed_tmp = "%s.gz"%output_reads_uncompressed
        for f in [output_reads_uncompressed, output_reads_compressed_tmp]: remove_file(f)

        # get lines
        run_cmd("zcat %s 2>/dev/null | head -n %i > %s 2>/dev/null"%(input_reads, nreads*4, output_reads_uncompressed))

        # compress and exit
        run_cmd("gzip %s"%output_reads_uncompressed)

        # check
        input_nlines = int(str(subprocess.check_output("zcat %s | wc -l"%input_reads, shell=True).split()[0]).split("'")[1])
        nlines = int(str(subprocess.check_output("zcat %s | wc -l"%output_reads_compressed_tmp, shell=True).split()[0]).split("'")[1])
        if (nreads*4)<=input_nlines and nlines!=(nreads*4): raise ValueError("you should have %i lines"%nreads)

        # keep
        os.rename(output_reads_compressed_tmp, output_reads)



def save_object(obj, filename):
    
    """ This is for saving python objects """

    filename_tmp = "%s.tmp"%filename
    remove_file(filename_tmp)
    
    with open(filename_tmp, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

    os.rename(filename_tmp, filename)


def load_object(filename):
    
    """ This is for loading python objects  in a fast way"""
    
    return pickle.load(open(filename,"rb"))


def get_read_length_from_fastqgz(file, subsample_reads=100000):

    """Takes a fastqgz file and returns the avg read length for subsample_reads reads"""

    # get the number of lines into file
    file_nreads = "%s.read_len.%ireads"%(file, subsample_reads)

    file_nreads_tmp = "%s.tmp"%file_nreads
    if file_is_empty(file_nreads):

        unpigz_stderr = "%s.generating.stderr"%file_nreads_tmp
        run_cmd("unpigz -c %s 2>%s | head -n %i | awk 'NR%s4==2' | awk '{ print length }' > %s 2>%s"%(file, unpigz_stderr, subsample_reads*4, "%", file_nreads_tmp, unpigz_stderr), env="perSVade_env")
        remove_file(unpigz_stderr)

        os.rename(file_nreads_tmp, file_nreads)

    df_nreads = pd.read_csv(file_nreads, header=None, sep="\t", names=["read_len"])
    if len(df_nreads)!=subsample_reads: raise ValueError("invalid len")

    mean_read_len = int(round(np.mean(df_nreads.read_len), 0))
    return mean_read_len

def get_source_data(runID):

    """Gets rowID and returns whether it is from SRA or from sequenced"""

    if runID.split("_")[0] in {"OSCAR01", "OSCAR02", "OSCAR03"}: return "sequenced_here"
    elif runID.startswith("SRR") or runID.startswith("ERR"): return "SRA"
    else: raise ValueError(runID)

def get_chr_to_len(genome, replace=False):

    chr_to_len_file = "%s.chr_to_len.py"%genome
    chr_to_len_file_tmp = "%s.tmp"%chr_to_len_file

    if file_is_empty(chr_to_len_file) or replace is True:

        remove_file(chr_to_len_file_tmp)

        # define chromosome_to_length for a genome
        chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(genome, "fasta")}

        # save
        save_object(chr_to_len, chr_to_len_file_tmp)
        os.rename(chr_to_len_file_tmp, chr_to_len_file)

    else: chr_to_len = load_object(chr_to_len_file)

    return chr_to_len

def plot_coverage_and_Cparapsilosis_mapping_all_samples(PlotsDir, df_reads_info, min_pct_reads_para, min_coverage, ylim=[-1,101]):

    """Plots coverage vs % mapping to Cpara"""

    fig = plt.figure(figsize=(7,7))

    #df_reads_info["log_pct_reads_Cparapsilosis"] = np.log10(df_reads_info.pct_reads_Cparapsilosis + 0.1)
    ax = sns.scatterplot(data=df_reads_info, x="expected_coverage", y="pct_reads_Cparapsilosis", hue="source_data", style="source_data", edgecolor="gray")
    plt.axvline(min_coverage, color="gray", linestyle="--", linewidth=.7)
    #plt.axhline(np.log10(min_pct_reads_para + 0.1), color="gray", linestyle="--", linewidth=.7)
    plt.axhline(min_pct_reads_para, color="gray", linestyle="--", linewidth=.7)
    ax.set_ylim(ylim)

    fig.savefig("%s/all_runs_expected_cov_vs_Cpara_mapping_ylim=%s.pdf"%(PlotsDir, ylim), bbox_inches='tight');
    plt.close(fig)



def plot_x_vs_y_df_reads_info(PlotsDir, xfield, yfield, ylim, hlines_ys, vlines_xs, df_plot, prefix, hue="soutce_data"):

    """Scatterplots"""

    # set auto ylim
    df_plot = df_plot.copy()
    if ylim is None: 
        max_y = max(df_plot[yfield])
        min_y = min(df_plot[yfield])
        range_y = max_y-min_y
        ylim = [min_y-(range_y*0.05), max_y+(range_y*0.05)]

    filename = "%s/%s_%s_vs_%s_ylim=%.2f-%.2f.pdf"%(PlotsDir, prefix, xfield, yfield, ylim[0], ylim[1])
    if file_is_empty(filename) or True:

        fig = plt.figure(figsize=(3.5,3.5))

        if hue=="source_data":
            palette_dict = {"SRA":"tab:blue", "sequenced_here":"tab:orange"}
        else:
            all_hues = sorted(set(df_plot[hue]))
            palette_dict = get_value_to_color(all_hues, palette="tab10", n=len(all_hues))[0]

        ax = sns.scatterplot(data=df_plot, x=xfield, y=yfield, hue=hue, palette=palette_dict, linewidth=.8, edgecolor=[palette_dict[x] for x in df_plot[hue]])


        # modify the face
        for artist in ax.collections: artist.set_facecolor('none')

        for x in vlines_xs: plt.axvline(x, color="gray", linestyle="--", linewidth=.7)
        for y in hlines_ys: plt.axhline(y, color="gray", linestyle="--", linewidth=.7)

        if not ylim is None: ax.set_ylim(ylim)

        ax.legend(bbox_to_anchor=(1.01,1), loc="upper left")

        print("saving", filename)
        print(df_plot)
        fig.savefig(filename, bbox_inches='tight')
        plt.close(fig)



def get_df_reads_info(df_paths, CurDir, list_cmds, ProcessedDataDir, ReadsDir, genome_Cparapsilosis, GenomesDir, min_coverage, sorted_aligners, VarCallDir, gff_Cparapsilosis, submit_jobs, homozygous_Corthopsilosis_srrs):

    """Gets a df_paths and keeps adding it with also list_cmds. It generates many read-mapping strategies"""

    df_reads_info_file = "%s/df_reads_info.py"%ProcessedDataDir
    if file_is_empty(df_reads_info_file):
        print("generating read info...")

        # init dicts
        reads_info_dict = {}

        for I, runID in enumerate(sorted(set(df_paths.runID))):

            # debug
            # if not runID in homozygous_Corthopsilosis_srrs: continue # debug

            # log
            print(I+1, runID)

            # define reads
            raw_reads1 = df_paths.loc[(runID, "R1")].raw_reads
            raw_reads2 = df_paths.loc[(runID, "R2")].raw_reads
            source_data_type = df_paths.loc[(runID, "R1")].source_data_type


            # downsample 100k reads
            folder_runID = "%s/%s"%(ReadsDir, runID); make_folder(folder_runID)
            folder_100k_reads = "%s/reads_100k"%folder_runID; make_folder(folder_100k_reads)
            raw_reads1_100k = "%s/raw_reads_R1.fastq.gz"%folder_100k_reads
            raw_reads2_100k = "%s/raw_reads_R2.fastq.gz"%folder_100k_reads

            cmd_tuple = ("reads100k", "%s/downsample_reads.py %s %s %s %s %i"%(CurDir, raw_reads1, raw_reads2, raw_reads1_100k, raw_reads2_100k, 100000), "Cparapsilosis_popgen_env", 4)

            if file_is_empty(raw_reads2_100k):
                if submit_jobs is True: list_cmds.append(cmd_tuple); continue
                else: run_cmd_tuple(cmd_tuple)

            # run with other Candidas the depooling pipeline to asses contaminations
            outdir_depool = "%s/depooling_various_Candidas"%(folder_runID)
            list_species = ["Calbicans", "Cauris", "Cmetapsilosis", "Corthopsilosis"]
            list_genomes = ["%s/%s.fasta"%(GenomesDir, spp) for spp in list_species]
            list_species.append("Cparapsilosis")
            list_genomes.append(genome_Cparapsilosis)

            read_splitting_statistics_table = "%s/read_splitting_statistics.tbl"%outdir_depool

            cmd_tuple = ("depool_100k", "%s/scripts/pooledSeq_splitting/split_pooledSeqWGS2_parallel.py --fastq1 %s --fastq2 %s --outdir %s --sampleNames %s --refGenomes %s --threads %i --nreads_per_chunk 25000"%(ParentDir, raw_reads1_100k, raw_reads2_100k, outdir_depool, ",".join(list_species), ",".join(list_genomes), 16), "pooledSeq_env", 16)

            if file_is_empty(read_splitting_statistics_table):
                if submit_jobs is True: list_cmds.append(cmd_tuple); continue
                else: run_cmd_tuple(cmd_tuple)

            for f in os.listdir(outdir_depool):
                file_path = "%s/%s"%(outdir_depool, f)
                if file_path!=read_splitting_statistics_table: remove_file(file_path) 

            # count reads per dataset
            nreads_file = "%s/nreads.txt"%folder_runID
            nreads_file_tmp = "%s.tmp"%nreads_file

            cmd_tuple = ("count_reads", "zcat %s | wc -l > %s && mv %s %s"%(raw_reads1, nreads_file_tmp, nreads_file_tmp, nreads_file), "Cparapsilosis_popgen_env", 4)

            if file_is_empty(nreads_file):
                if submit_jobs is True: list_cmds.append(cmd_tuple); continue
                else: run_cmd_tuple(cmd_tuple)

            # get the read len
            read_len1 = get_read_length_from_fastqgz(raw_reads1_100k, subsample_reads=100000)
            read_len2 = get_read_length_from_fastqgz(raw_reads2_100k, subsample_reads=100000)
            if abs(read_len1-read_len2)>=10: 
                print(runID, read_len1, read_len2)
                raise ValueError("read lens btw 1 and 2 are different")

            # calculate the expected coverage
            genome_len = sum(get_chr_to_len(genome_Cparapsilosis).values())
            npairs = int(open(nreads_file, "r").readlines()[0].strip())/4
            expected_coverage = (read_len1 * 2 * npairs)/genome_len

            # get the downsampled reads to min_coverage
            downsampling_nreads = int(round((min_coverage * genome_len) / (read_len1 * 2), 0))
            folder_min_cov = "%s/reads_%ix"%(folder_runID, min_coverage); make_folder(folder_min_cov)
            raw_reads1_min_cov = "%s/raw_reads_R1.fastq.gz"%folder_min_cov
            raw_reads2_min_cov = "%s/raw_reads_R2.fastq.gz"%folder_min_cov

            cmd_tuple = ("reads_min_cov", "%s/downsample_reads.py %s %s %s %s %i"%(CurDir, raw_reads1, raw_reads2, raw_reads1_min_cov, raw_reads2_min_cov, downsampling_nreads), "Cparapsilosis_popgen_env", 4)

            if file_is_empty(raw_reads2_min_cov):
                if submit_jobs is True: list_cmds.append(cmd_tuple); continue
                else: run_cmd_tuple(cmd_tuple)

            # get the trimmed reads
            folder_trimmed_reads = "%s_trimmed"%folder_min_cov; make_folder(folder_trimmed_reads)
            trimmed_reads1 = "%s/trimmed_reads1.fastq.gz"%folder_trimmed_reads
            trimmed_reads2 = "%s/trimmed_reads2.fastq.gz"%folder_trimmed_reads

            cmd_tuple = ("trim_reads", "%s trim_reads_and_QC -o %s --threads 28 -f1 %s -f2 %s --fractionRAM_to_dedicate 0.8 --verbose --fraction_available_mem 0.5"%(perSVade_wrapper, folder_trimmed_reads, raw_reads1_min_cov, raw_reads2_min_cov), "perSVade_env", 28)

            if file_is_empty("%s/perSVade_finished_file.txt"%folder_trimmed_reads):
                if submit_jobs is True: list_cmds.append(cmd_tuple); continue
                else: run_cmd_tuple(cmd_tuple)

            # map the reads & calculate coverage with different aligners
            print("geting alignment jobs...")
            jobs_to_run = False
            varcall_dir_runID = "%s/%s"%(VarCallDir, runID); make_folder(varcall_dir_runID)
            for aligner in sorted_aligners:

                # I saw that hisat2 was not mapping sufficiently for C. orthopsilosis
                if runID in homozygous_Corthopsilosis_srrs and aligner=="hisat2_no_spliced": continue

                # define dir
                aligner_dir = "%s/aligner-%s"%(varcall_dir_runID, aligner); make_folder(aligner_dir)

                # map reads
                align_reads_outdir = "%s/align_reads"%aligner_dir

                cmd_tuple = ("map_reads", "%s align_reads -o %s --threads 28 -f1 %s -f2 %s -r %s --min_chromosome_len 5000 --aligner %s --fractionRAM_to_dedicate 0.8 --verbose --fraction_available_mem 0.5 --skip_marking_duplicates"%(perSVade_wrapper, align_reads_outdir, trimmed_reads1, trimmed_reads2, genome_Cparapsilosis, aligner), "perSVade_env", 28)

                if file_is_empty("%s/perSVade_finished_file.txt"%align_reads_outdir):
                    if submit_jobs is True: list_cmds.append(cmd_tuple); jobs_to_run = True; continue
                    else: run_cmd_tuple(cmd_tuple)

                # calculate coverage per gene
                sorted_bam = "%s/aligned_reads.bam.sorted"%align_reads_outdir
                cov_per_gene_outdir = "%s/cov_per_gene"%aligner_dir

                cmd_tuple = ("cov_per_gene", "%s get_cov_genes -o %s --threads 28 -r %s --min_chromosome_len 5000 -gff %s -sbam %s --fractionRAM_to_dedicate 0.8 --verbose --fraction_available_mem 0.5"%(perSVade_wrapper, cov_per_gene_outdir, genome_Cparapsilosis, gff_Cparapsilosis, sorted_bam), "perSVade_env", 28)

                if file_is_empty("%s/perSVade_finished_file.txt"%cov_per_gene_outdir):
                    if submit_jobs is True: list_cmds.append(cmd_tuple); jobs_to_run = True; continue
                    else: run_cmd_tuple(cmd_tuple)

            # do not continue if there are reads to map
            if jobs_to_run is True: continue

            # keep info into reads_info_dict
            reads_info_dict[runID] = {"runID":runID, "source_data": get_source_data(runID), "npairs":npairs, "read_len1":read_len1, "read_len2":read_len2, "expected_coverage":expected_coverage, "source_data_type":source_data_type}

            type_reads_to_pct = pd.read_csv(read_splitting_statistics_table, sep="\t").set_index("type_reads").pct_total
            for spp in (list_species + ["unmapped", "multimapped"]): reads_info_dict[runID]["pct_reads_%s"%spp] = type_reads_to_pct[spp]
            
        # save at the end
        if len(list_cmds)==0:
            df_reads_info = pd.DataFrame(reads_info_dict).transpose()
            save_object(df_reads_info, df_reads_info_file)


    if len(list_cmds)==0: return load_object(df_reads_info_file)
    else: return pd.DataFrame()


def get_sampleID_from_runID(r):
    if get_source_data(r.runID)=="SRA": return r.runID
    elif get_source_data(r.runID)=="sequenced_here": return r.runID.split("_")[1]

def rsync_file(origin, dest):

    """syncs one file to the other"""

    #if file_is_empty(dest):

    dest_tmp = "%s.tmp"%dest
    run_cmd("rsync %s %s"%(origin, dest_tmp))
    os.rename(dest_tmp, dest)


def define_origin_zip_Cparapsilosis_popgen(runID, rID):

    """Takes a runID and a rID and returns the origin zip"""
        
    ReadsDir = "%s/Cparapsilosis_popgen/data/reads"%ParentDir
    min_coverage = 50 # used all over
    origin_zip = "%s/%s/reads_%ix_trimmed/trimmed_reads%i_fastqc.zip"%(ReadsDir, runID, min_coverage, rID)
    if file_is_empty(origin_zip): raise ValueError("%s should exist"%origin_zip)
    return origin_zip

def link_files_for_multiqc(Ir, rID, origin_zip, dest_folder, linked_files_dir):

    """Links one folder for multiqc"""

    print("Preparing multiqc files for run %i, R%i"%(Ir+1, rID))

    dest_zip = "%s.zip"%(dest_folder)
    rsync_file(origin_zip, dest_zip)
    dest_folder_tmp = "%s_tmp"%dest_folder; delete_folder(dest_folder_tmp)
    run_cmd("unzip %s -d %s > /dev/null 2>&1"%(dest_zip, dest_folder_tmp)); remove_file(dest_zip)
    os.rename("%s/raw_reads%i.fastq.gz.trimmed_fastqc"%(dest_folder_tmp, rID), dest_folder)
    delete_folder(dest_folder_tmp)

def run_multiqc_several_trimmed_reads(outdir, sorted_runIDs, function_origin_zip):

    """Runs multiqc into outdir for all trimmed reads"""

    # define dirs
    make_folder(outdir)
    linked_files_dir = "%s/multiqc_input_files"%outdir; make_folder(linked_files_dir) # store files
    multiqc_outdir_report = "%s/multiqc_output"%(outdir) # store output

    # if not done, run
    if file_is_empty("%s/multiqc_report.html"%multiqc_outdir_report):
        print("generating %s..."%multiqc_outdir_report)
        delete_folder(multiqc_outdir_report)

        # link all the dirs to linked_files_dir
        print("preparing files for multiqc")
        os.chdir(linked_files_dir)

        # get inputs
        inputs_fn = []
        for Ir, runID in enumerate(sorted_runIDs):
            for rID in [1, 2]:
                origin_zip = function_origin_zip(runID, rID)
                dest_folder = "%s/%s_R%i"%(linked_files_dir, runID, rID)
                if not os.path.isdir(dest_folder): inputs_fn.append((Ir, rID, origin_zip, dest_folder, linked_files_dir))

        # link files
        if len(inputs_fn)>0:
            #inputs_fn= inputs_fn[0:1] # debug
            with multiproc.Pool(multiproc.cpu_count()) as pool:
                pool.starmap(link_files_for_multiqc, inputs_fn, chunksize=1)
                pool.close()

        # run multiqc on the linked files from linked_files_dir
        print("running multiqc")
        run_cmd("multiqc %s -o %s --dirs-depth 1 -d -f -s --interactive"%(linked_files_dir, multiqc_outdir_report))
        print("multiqc finished")



def get_tab_as_df_and_add_sampleID(Is, sampleID, file, expected_rows, fields="all"):

    """returns the variant calling df as a df"""
    print("getting small vars for sample %i %s"%(Is+1, sampleID))

    # get df
    df = pd.read_csv(file, sep="\t")
    if fields!="all": df = df[fields]
    if len(df)!=expected_rows: raise ValueError("There should be one row per gene") 
    df["sampleID"] = str(sampleID)

    return df

def integrate_get_cov_genes_several_samples(sample_to_outdir, filename, threads, fields="all"):

    """Generates a call set that includes all the stacked results of get_cov_genes into filename"""

    if file_is_empty(filename):
        print("running integrate_get_cov_genes_several_samples")

        # integrate in parallel (validate that there are 5960 genes in each)
        inputs_fn = [(Is, s, "%s/genes_and_regions_coverage.tab"%o, 5960, fields) for Is, (s, o) in enumerate(sample_to_outdir.items())]
        with multiproc.Pool(threads) as pool:
            coverage_df = pd.concat(pool.starmap(get_tab_as_df_and_add_sampleID, inputs_fn, chunksize=1)) 
            pool.close()
            pool.terminate()

        # save
        print("saving")
        save_object(coverage_df, filename)

def get_cov_df_per_window_one_run(runID, cov_file):

    """gets the cov per windows with only key files"""

    print("loading run ID %s"%runID)

    df_cov = pd.read_csv(cov_file, sep="\t")[["chromosome", "start", "end", "mediancov_1", "percentcovered_1", 'fraction_repeats', 'GCcontent', 'median_mappability', 'fraction_N_bases']].copy()

    # add rel cov
    median_cov_all = np.median(df_cov[df_cov.chromosome!=mtDNA_chrom].mediancov_1.apply(float))
    if pd.isna(median_cov_all): raise ValueError("invalid cov %s"%(median_cov_all))
    
    if median_cov_all>0: df_cov["median_coverage_rel_to_gDNA"] = df_cov.mediancov_1 / median_cov_all
    else: df_cov["median_coverage_rel_to_gDNA"] = -1

    check_no_nans_series(df_cov["median_coverage_rel_to_gDNA"])

    # keep
    df_cov["runID"] = runID
    return df_cov

def get_genome_excluding_chromosomes(genome, set_chroms_exclude):

    """gets a genome exlcuding some chromosomes"""

    modified_genome = "%s.exlcuding_%s.fasta"%(genome, "_".join(sorted(set_chroms_exclude)))
    if file_is_empty(modified_genome):

        # define tmp
        modified_genome_tmp = "%s.tmp.fasta"%modified_genome

        # keep chroms
        seq_records_keep = [seq for seq in SeqIO.parse(genome, "fasta") if not seq.id in set_chroms_exclude]
        SeqIO.write(seq_records_keep, modified_genome_tmp, "fasta")
        os.rename(modified_genome_tmp, modified_genome)

    return modified_genome

def get_genomic_df_excluding_chromosomes(tab_file, set_chroms_exclude):

    """gets a genome exlcuding some chromosomes"""

    modified_tab_file = "%s.excluding_%s.%s"%(tab_file, "_".join(sorted(set_chroms_exclude)), tab_file.split(".")[-1])
    if file_is_empty(modified_tab_file):

        # load
        df = get_tab_as_df_or_empty_df(tab_file)

        # keep
        df_filt = df[~df.chromosome.apply(str).isin(set_chroms_exclude)].copy()
        save_df_as_tab(df_filt, modified_tab_file)

    return modified_tab_file

def generate_df_coverage_per_window(list_cmds, aligner, window_size, VarCallDir, sorted_runs, genome_Cparapsilosis, simple_repeats_file, ProcessedDataDir, homozygous_Corthopsilosis_srrs):

    """generates a df with the coverage per window"""

    # define final file
    coverage_df_file = "%s/coverage_per_%ibp_windows_all_samples_%s.py"%(ProcessedDataDir, window_size, aligner)

    # get mappability per pos file
    outdir_mappability = "%s/mappability_outdir"%ProcessedDataDir
    final_file_map = "%s/map_calculated.txt"%outdir_mappability

    if file_is_empty(final_file_map):
        run_cmd("%s/Cparapsilosis_popgen/get_mappability_file_whole_genome.py %s 112 %s %s"%(ParentDir, genome_Cparapsilosis, final_file_map, outdir_mappability), env="perSVade_env")

    map_per_pos_file = "%s/genome.fasta.mappability_per_position.bed"%outdir_mappability

    # define the used ref_genome, depending on window_size
    if window_size<=10000:
        ref_genome = genome_Cparapsilosis
        repeats_file = simple_repeats_file
        map_per_pos = map_per_pos_file

    else:
        ref_genome = get_genome_excluding_chromosomes(genome_Cparapsilosis, {mtDNA_chrom})
        repeats_file = get_genomic_df_excluding_chromosomes(simple_repeats_file, {mtDNA_chrom})
        map_per_pos = get_genomic_df_excluding_chromosomes(map_per_pos_file, {mtDNA_chrom})

    # generate call_CNV cmds
    list_runIDs_files = []
    for runID in sorted_runs:

        # discard some runs
        if aligner=="hisat2_no_spliced" and runID in homozygous_Corthopsilosis_srrs: continue

        # get files
        outdir_aligner = "%s/%s/aligner-%s"%(VarCallDir, runID, aligner)
        sorted_bam = "%s/align_reads/aligned_reads.bam.sorted"%outdir_aligner
        if file_is_empty(sorted_bam): raise ValueError("bam should exist")

        # get call_CVNs command
        threads_per_job = 28
        outdir_coverage = "%s/coverage_%ibp_windows"%(outdir_aligner, window_size)
        list_runIDs_files.append((runID, "%s/final_df_coverage.tab"%outdir_coverage))

        if file_is_empty("%s/perSVade_finished_file.txt"%(outdir_coverage)):

            make_folder(outdir_coverage)
            tmpdir = "%s/perSVade_tmp"%outdir_coverage
            cmd_callCNVs = "export PERSVADE_TMPDIR=%s && %s call_CNVs -o %s -r %s -sbam %s --mitochondrial_chromosome no_mitochondria --skip_CNV_calling -p 2 --window_size_CNVcalling %i --min_chromosome_len 5000 --verbose --fraction_available_mem %.2f --threads %i --fractionRAM_to_dedicate 0.8 --skip_coverage_correction --mappability_file %s --average_cov_measure median --max_fraction_N_bases 1.0 --max_fraction_repeats 1.0 --repeats_file %s --min_median_mappability 0.0 --cnv_calling_algs HMMcopy --replace"%(tmpdir, perSVade_wrapper, outdir_coverage, ref_genome, sorted_bam, window_size, threads_per_job/112, threads_per_job, map_per_pos, repeats_file)

            list_cmds.append(("call_CNVs_broad", cmd_callCNVs, "perSVade_env", threads_per_job))

    # integrate
    if len(list_cmds)==0:
        print("integrating...")

        # get the integrated df
        if file_is_empty(coverage_df_file):
            print("running coverage file")

            with multiproc.Pool(112) as pool:
                coverage_df = pd.concat(pool.starmap(get_cov_df_per_window_one_run, list_runIDs_files, chunksize=1)) 
                pool.close()
                pool.terminate()

            print("saving")
            save_object(coverage_df, coverage_df_file)

    return coverage_df_file

def check_no_nans_series(x):

    """Raise value error if nans"""

    if any(pd.isna(x)): raise ValueError("There can't be nans in series %s. %s"%(x, x[pd.isna(x)]))

def check_no_nans_in_df(df, ignore_keys=set()):

    """Returns an error if there ar nans in a df"""

    keys_consider = list(set(df.keys()).difference(ignore_keys))
    if any(pd.isna(df[keys_consider]).apply(any, axis=1)): 
        for k in keys_consider: check_no_nans_series(df[k])




def get_array_distances_previous_pos_one_chrom(df_c):

    """gets array distances"""

    array_pos = df_c.iloc[0:-1].median_coverage_rel_to_gDNA.values
    array_pos_plus1 = df_c.iloc[1:].median_coverage_rel_to_gDNA.values

    if len(array_pos)!=len(array_pos_plus1): raise ValueError("should be the same")

    return array_pos_plus1-array_pos

def get_r_df_reads_info_with_added_info(r, df_snps_all, min_rel_cov, max_rel_cov, min_n_vars):


    """takes an r and adds info about the seq run """
    print(r.runID)

    # add coverage info
    for genome in ["mtDNA", "gDNA"]:

        # get dfs
        df_snps = df_snps_all.loc[[genome]].copy()

        # add the cov distance to the previous window
        df_snps = df_snps.sort_values(by=["chromosome", "pos"])
        array_distances = np.hstack(df_snps[["chromosome", "pos", "median_coverage_rel_to_gDNA"]].groupby("chromosome").apply(get_array_distances_previous_pos_one_chrom).values.flatten())

        dist_f = "rel_cov_diff_windows_%s"%genome
        r["median_%s"%dist_f] = np.median(array_distances)
        r["mean_%s"%dist_f] = np.mean(array_distances)
        r["std_%s"%dist_f] = np.std(array_distances)
        r["mad_%s"%dist_f] = scipy.stats.median_abs_deviation(array_distances, nan_policy="raise")

        for pct in [5, 25, 75, 95]:
            r["pct%i_%s"%(pct, dist_f)] = np.percentile(array_distances, pct)

        r["IQR_%s"%dist_f] = r["pct75_%s"%dist_f] - r["pct25_%s"%dist_f]

        # add the simley face effect
        if len(set(df_snps.median_coverage_rel_to_gDNA))>1:
            r["r_distTelomere_vs_cov_%s"%genome] = scipy.stats.spearmanr(df_snps.dist_to_telomere.values, df_snps.median_coverage_rel_to_gDNA.values, nan_policy="raise")[0]

        else:
            r["r_distTelomere_vs_cov_%s"%genome] = 0.0

        # add the coverage measurement values (rel and abs)
        for cov_f, cov_name in [("median_coverage_rel_to_gDNA", "rel_median_cov"), ("median_coverage", "median_cov")]:
            cov_values = df_snps[cov_f].values

            r_field = "%s_%s"%(cov_name, genome)
            r["median_%s"%r_field] = np.median(cov_values)
            r["mean_%s"%r_field] = np.mean(cov_values)
            r["std_%s"%r_field] = np.std(cov_values)
            r["pct25_%s"%r_field] = np.percentile(cov_values, 25)
            r["pct75_%s"%r_field] = np.percentile(cov_values, 75)
            r["IQR_%s"%r_field] = r["pct75_%s"%r_field] - r["pct25_%s"%r_field]
            r["mad_%s"%r_field] = scipy.stats.median_abs_deviation(cov_values, nan_policy="raise")

        # add the mean pct covered
        r["pct_cov_%s"%(genome)] = np.mean(df_snps.percentcovered_1.values)

        # add the SNPs per kb, considering all windows, for all zygosities
        all_zygosities = ["homo", "hetero", "unknown_zygosity"]
        for zygosity in all_zygosities:
            df_snps_z = df_snps[(df_snps.zygosity==zygosity)].copy()
            if len(df_snps_z)==0: raise ValueError("invalid len")
            r["%s_SNPs_per_kb_%s"%(zygosity, genome)] = (sum(df_snps_z.nvars) / sum(df_snps_z.len_window))*1000

        # add the len of the genome
        r["whole_genome_len_%s"%genome] = sum(df_snps[(df_snps.zygosity=="homo")].len_window)
        print(r["whole_genome_len_%s"%genome])

        # filter to keep only non-aneuploid regions
        df_snps = df_snps[(df_snps.median_coverage_rel_to_gDNA>=min_rel_cov) & (df_snps.median_coverage_rel_to_gDNA<=max_rel_cov)].copy()

        # add the len of the aneuploid genome
        if len(df_snps)>0:
            all_genome_lens = set(df_snps.groupby("zygosity").apply(lambda df_z: sum(df_z.len_window)))
            if len(all_genome_lens)!=1: raise ValueError("there are >1 genome lens")
            r["non_aneuploid_genome_len_%s"%genome] = next(iter(all_genome_lens))

        else: r["non_aneuploid_genome_len_%s"%genome] = 0

        # get measurements about AF 
        for zygosity in all_zygosities:

            # get the df for the zygosity, also for sufficient coverage 
            df_snps_z = df_snps[(df_snps.zygosity==zygosity)].copy()
            nvars_z = sum(df_snps_z.nvars)
            r["%s_nSNPs_non_aneuploid_%s"%(zygosity, genome)]  = nvars_z

            prefix_f = "%s_%s_AF"%(genome, zygosity)
            if nvars_z>=min_n_vars:

                array_AFs = np.hstack(df_snps_z.AF_list.values.flatten())
                if len(array_AFs)!=nvars_z: raise ValueError("expected %i vars. this array: %s"%(nvars_z, array_AFs))

                r["%s_mean"%prefix_f] = np.mean(array_AFs)
                r["%s_std"%prefix_f] = np.std(array_AFs)
                r["%s_median"%prefix_f] = np.median(array_AFs)
                r["%s_pct25"%prefix_f] = np.percentile(array_AFs, 25)
                r["%s_pct75"%prefix_f] = np.percentile(array_AFs, 75)
                r["%s_IQR"%prefix_f] = r["%s_pct75"%prefix_f] - r["%s_pct25"%prefix_f]
                r["%s_mad"%prefix_f] = scipy.stats.median_abs_deviation(array_AFs, nan_policy="raise")

            else: 

                for suff_f in ["mean", "std", "median", "pct25", "pct75", "IQR", "mad"]:
                    r["%s_%s"%(prefix_f, suff_f)] = -1

    # checks
    check_no_nans_series(r)
    return r

def run_kraken_subset_runs(outdir_all, df_runs, type_runs, ReadsDir, only_sequenced_here=True):

    """Runs kraken for some runs"""
    print("running kraken %s"%type_runs)

    # keep
    df_runs = df_runs.copy()
    if only_sequenced_here: df_runs = df_runs[df_runs.source_data=="sequenced_here"].copy()

    # check
    if len(df_runs)==0:
        print("WARNING: Kraken not run because there are no samples...")
        return

    # run kraken
    kraken_samples = sorted(df_runs.runID)
    kraken_outdir_all = "%s/runs_%s"%(outdir_all, type_runs); make_folder(kraken_outdir_all)
    kraken_outdir_reports = "%s/report_%s"%(outdir_all, type_runs); make_folder(kraken_outdir_reports)

    for runID in kraken_samples:
        for kraken_db, kraken_db_path in [("DB-all", "/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/KRAKEN2_DB_COMPLETE"), ("DB-fungi", "/gpfs/projects/bsc40/project/pipelines/WGS/KRAKEN2_DB/Broad_fungal_microbiome_db")]:

            # define kraken outdir
            kraken_outdir_run = "%s/%s_run_%s"%(kraken_outdir_all, runID, kraken_db); make_folder(kraken_outdir_run)

            # run kraken
            reads_dir_run = "%s/%s/reads_100k"%(ReadsDir, runID)
            if file_is_empty("%s/kraken.report.txt"%kraken_outdir_run):
                print("running kraken for %s, %s"%(runID, kraken_db))

                cmd = "kraken2 --db %s --gzip-compressed --paired %s/raw_reads_R1.fastq.gz %s/raw_reads_R2.fastq.gz --output %s/kraken.txt --report %s/kraken.report.txt --threads 48 --use-names"%(kraken_db_path, reads_dir_run, reads_dir_run, kraken_outdir_run, kraken_outdir_run)  # needs --highmem
                run_cmd_pipelines_env(cmd, env="meTAline")

            # link
            soft_link_files("%s/kraken.report.txt"%kraken_outdir_run, "%s/%s_run_%s.report.txt"%(kraken_outdir_reports, runID, kraken_db))

    # run Krona to see all kraken results
    outfile_krona = "%s/%s_krona.html"%(outdir_all, type_runs)
    if file_is_empty(outfile_krona):
        remove_file(outfile_krona)

        # print("updating taxonomy") # only once, needs internet
        # run_cmd_pipelines_env("/gpfs/projects/bsc40/project/pipelines/anaconda3/envs/meTAline/opt/krona/updateTaxonomy.sh", env="meTAline")

        print("running krona")
        run_cmd_pipelines_env("ktImportTaxonomy -t 5 -m 3 -o %s %s/*.report.txt "%(outfile_krona, kraken_outdir_reports), env="meTAline")


def get_df_snps_specific_samples(aligner_to_df_small_vars_filt_file, sorted_runs, outdir, replace=False):

    """
    Gets the coverage and AFs for all SNPs with variant files
    """

    make_folder(outdir)
    df_snps_file = "%s/df_snps_%iruns_%s_to_%s.py"%(outdir, len(sorted_runs), sorted_runs[0], sorted_runs[-1])
    if file_is_empty(df_snps_file):
        print("getting %s"%df_snps_file)

        df_snps = pd.DataFrame()
        for aligner, df_small_vars_filt_file in aligner_to_df_small_vars_filt_file.items():
            print(aligner)

            # load df
            df_small_vars_filt = load_object(df_small_vars_filt_file)

            # check
            strange_runs = set(sorted_runs).difference(set(df_small_vars_filt.runID))
            if len(strange_runs)>0: raise ValueError("There are some strange runs: %s"%strange_runs)

            # keep some
            df_snps_aln = df_small_vars_filt[(df_small_vars_filt.type_var=="SNP") & (df_small_vars_filt.runID.isin(set(sorted_runs)))].copy()
            df_snps = pd.concat([df_snps, df_snps_aln]).reset_index(drop=True)

        print("saving")
        save_object(df_snps, df_snps_file)

    return load_object(df_snps_file)

def plot_distribution_AFs_specific_samples(df_snps, sorted_runs, p_dir):

    """plots dist of SNPs for specific samples"""

    print("plotting strange hetero SNPs")

    make_folder(p_dir)

    # filter
    df_plot = df_snps[(df_snps.unique_map_score>=0.25) & (df_snps.nvars_SNP_in_pos==1) & (df_snps["median_coverage_rel_to_gDNA_1.0kb_window"]<=2.5) & (df_snps["percentcovered_1_1.0kb_window"]>90)].copy().rename(columns={"aligner":"a", "runID":"r"})

    # violin plot
    for target_gts, gt_set in [("all", {"1/1", "0/1", "."}), ("non_homo", {"0/1", "."})]:

        df_plot_violin = df_plot[(df_plot["median_coverage_rel_to_gDNA_1.0kb_window"]>=0.8) & (df_plot["median_coverage_rel_to_gDNA_1.0kb_window"]<=1.2) & (df_plot.GT.isin(gt_set))] # (df_plot.GT!="1/1")
        fig = plt.figure(figsize=(13,3))
        ax = sns.violinplot(data=df_plot_violin, x="r", y="AF", hue="a", fill=False, order=sorted_runs)
        for item in ax.get_xticklabels(): item.set_rotation(90)
        ax.legend(bbox_to_anchor=(1.01,1), loc="upper left")
        for y in [0.5, 1.0]: ax.axhline(y, ls="--", c="black", lw=.7) # real hetero
        for y in [0.33, 0.66]: ax.axhline(y, ls="--", c="red", lw=.7) # triploid
        for y in [0.25, 0.75]: ax.axhline(y, ls="--", c="blue", lw=.7) # tertra
        save_fig(fig, "AF_distribution_%s.pdf"%target_gts, p_dir)

    # scatterplot
    df_plot_scatter = df_plot[(df_plot.a=="bwa_mem") & (df_plot.GT!="1/1")].copy()
    fig = sns.relplot(kind="scatter", data=df_plot_scatter, x="median_coverage_rel_to_gDNA_1.0kb_window", y="AF", row="r", alpha=.4, height=2, aspect=1.2, s=4)
    for ax in fig.axes.flat:
        for x in [0.5, 1, 1.5]: ax.axvline(x, ls="--", c="gray", lw=.7)

        for y in [0.5]: ax.axhline(y, ls="--", c="black", lw=.7) # real hetero
        for y in [0.33, 0.66]: ax.axhline(y, ls="--", c="red", lw=.7) # triploid
        for y in [0.25, 0.75]: ax.axhline(y, ls="--", c="blue", lw=.7) # tertra

    save_fig(fig, "cov_vs_AF.pdf", p_dir)

def print_phenotypes_subset_runs(df_metadata, df_runs, type_runs):

    """prints phenotypes for discarded samples"""

    print("see phenotypes for %s"%type_runs)
    df_runs = df_runs[df_runs.species_run=="para"].copy().set_index("runID")
    df_metadata = df_metadata[df_metadata.sampleID.isin(set(df_runs.sampleID))].copy().set_index("sampleID")

    pheno_fields = ['hospital', 'AMB_MIC', '5FC_MIC', 'FLC_MIC', 'VRC_MIC', 'ITR_MIC', 'POS_MIC', 'IVZ_MIC', 'CAS_MIC', 'MIF_MIC', 'ANI_MIC', 'MIC ETEST 24H', 'MIC ETEST 48H', 'Colonies inside halo', 'MIC of recovered heteroresistant 24 h', 'MIC of recovered heteroresistant 48h', 'Agar invasiveness', 'MORPHOLOGY', 'Microfluidic behavior', 'Biofilm (High / Low / No)', '% VIABILITY ETHANOL', '% VIABILITY HYPLOCHLORITE', '% VIABILITY H202', '% VIABILITY CHLORHEXIDINE', '% VIABILITY ZWITTERGENT 3-14', '% VIABILITY SURFANIO']

    for runID,r in df_runs.iterrows():
        print()
        for f in pheno_fields: 

            f_val = df_metadata.loc[r.sampleID, f]
            if f.endswith("_MIC"): pct_val = "(%.1f percentile)"%((1 - (sum(df_metadata[f]>=f_val) / len(df_metadata)))*100)
            else: pct_val = ""

            print(runID, f, "=", f_val, pct_val)

def print_strain_numbers(ProcessedDataDir, ManualDataDir, df_metadata, df_metadata_binary_phenotypes):

    """Prints various strain numbers"""

    # define ortho SRRs
    homozygous_Corthopsilosis_srrs = get_homozygous_Corthopsilosis_srrs(ManualDataDir)

    # load dfs
    df_reads_info_filt = load_object("%s/df_reads_info_filtered_afterQC.py"%ProcessedDataDir)
    df_reads_info_bad = load_object("%s/df_reads_info_bad_runs.py"%ProcessedDataDir)
    df_reads_info_analysis = load_object("%s/df_reads_info_filtered_afterQC_samples_analyze.py"%ProcessedDataDir)
    df_reads_info = load_object("%s/df_reads_info_with_multiqc_and_coverage_data_and_AFs_bwa_mem_1000bp.py"%ProcessedDataDir)

    # define samples to discard because they are redundant
    samples_discard_run1 = {"OSCAR01_L21-060", "OSCAR01_L21-067"}

    # add the spp run
    df_reads_info["species_run"] = df_reads_info.runID.isin(homozygous_Corthopsilosis_srrs).map({True:"ortho", False:"para"})
    df_reads_info["sampleID"] = df_reads_info.apply(get_sampleID_from_runID, axis=1)

    # filter to keep para
    df_reads_info_filt = df_reads_info_filt[df_reads_info_filt.species_run=="para"].copy()
    df_reads_info_bad = df_reads_info_bad[df_reads_info_bad.species_run=="para"].copy()
    df_reads_info_analysis = df_reads_info_analysis[df_reads_info_analysis.species_run=="para"].copy()
    df_reads_info = df_reads_info[df_reads_info.species_run=="para"].copy()
    df_metadata = df_metadata[df_metadata.species_run=="para"].copy()
    
    # filter to discard
    df_reads_info_filt = df_reads_info_filt[~df_reads_info_filt.runID.isin(samples_discard_run1)].copy()
    df_reads_info_bad = df_reads_info_bad[~df_reads_info_bad.runID.isin(samples_discard_run1)].copy()
    df_reads_info_analysis = df_reads_info_analysis[~df_reads_info_analysis.runID.isin(samples_discard_run1)].copy()
    df_reads_info = df_reads_info[~df_reads_info.runID.isin(samples_discard_run1)].copy()
    df_metadata = df_metadata[~df_metadata.runID.isin(samples_discard_run1)].copy()

    # filter meta
    df_metadata = df_metadata[df_metadata.possibly_triploid==False].copy()

    # checks
    if (len(df_reads_info_filt)+len(df_reads_info_bad)) != len(df_reads_info):
        raise ValueError("invalid lens")

    for df in [df_reads_info_filt, df_reads_info_bad, df_reads_info_analysis, df_reads_info, df_metadata]:
        if len(df)!=len(df.sampleID.unique()):
            raise ValueError("should be unique sampleID")

    # define the type of data
    print("all # strains sequenced here", sum((df_reads_info.source_data!="SRA")))
    print("filt # strains sequenced here", sum(df_reads_info_analysis.source_data!="SRA"))
    print("all # strains SRA", sum(df_reads_info.source_data=="SRA"))
    print("filt # strains SRA", sum(df_reads_info_analysis.source_data=="SRA"))
    print("all # strains", len(df_reads_info))
    print("filt # strains", len(df_reads_info_analysis))
    print("metadata # strains", len(df_metadata))

    # print the years
    df_metadata_seqHere = df_metadata[df_metadata.source_data!="SRA"].copy()
    df_metadata_SRA = df_metadata[df_metadata.source_data=="SRA"].copy()
    print("years collection here", sorted(set(df_metadata_seqHere.year_collection)))
    print("# hospitals here", len(set(df_metadata_seqHere.hospital)))

    print("countries all", sorted(set(df_metadata.country_collection)))
    print("countries here", sorted(set(df_metadata_seqHere.country_collection)))
    print("countries SRA", sorted(set(df_metadata_SRA.country_collection)))

    print("cities here", sorted(set(df_metadata_seqHere.city_collection)))
    print("cities SRA", sorted(set(df_metadata_SRA.city_collection)))

    # binary phenos
    print()
    print("binary phenos")
    all_fields = []

    merged_sources = ["source_sterile", "source_sterile_vs_environmental", "source_blood", "source_blood_vs_environmental", "source_blood_vs_other_sterile", "source_skin", "source_respiratory"]

    for f in df_metadata_binary_phenotypes.keys():
        if f.startswith("source_") and not f in merged_sources:
            all_fields.append(get_paper_dicotomic_pheno(f))
            print(get_paper_dicotomic_pheno(f), "("+str(sum(df_metadata_binary_phenotypes[f]==1))+")")


    print()
    for f in df_metadata_binary_phenotypes.keys():
        if f in  merged_sources:
            all_fields.append(get_paper_dicotomic_pheno(f))
            print(get_paper_dicotomic_pheno(f),  "("+str(sum(df_metadata_binary_phenotypes[f]==1))+")")

    if len(all_fields)!=len(set(all_fields)):
        raise ValueError("not all unique")

    # print # of strains where virulence was measured
    df_virulence = df_metadata_seqHere[(df_metadata_seqHere["Agar invasiveness"]!="unknown") | (df_metadata_seqHere["MORPHOLOGY"]!="unknown") | (df_metadata_seqHere["Biofilm (High / Low / No)"]!="unknown")].copy()

    print(len(df_virulence), "strains some vir. measured")
    for f in ["Agar invasiveness", "MORPHOLOGY", "Biofilm (High / Low / No)"]:
        print(f, Counter(df_virulence[f]))

    df_microfl = df_metadata_seqHere[(df_metadata_seqHere['Microfluidic behavior']!="unknown")].copy()
    print(len(df_microfl), Counter(df_microfl['Microfluidic behavior']))


def get_integrated_variant_calls_all_datasets_with_validation(df_reads_info_analysis_with_validation, ValidationDir, VarCallDir, CurDir, list_cmds, sorted_aligners, simple_repeats_file, genome_Cparapsilosis, gff_Cparapsilosis):

    """Gets the datasets for all"""

    print("Getting integrated varcalls with validation")

    # define ProcessedDataDir validation
    ProcessedDataDir_validation = "%s/processed_data"%ValidationDir

    # run SV & CNV calling and integration
    df_reads_info_SVs = df_reads_info_analysis_with_validation[(df_reads_info_analysis_with_validation.has_smiley==False)].copy()
    SV_parms_file = "%s/data/SV_parameter_optimization/optimize_parameters_representative_runs/analyze_SV_parameters/more_conservative_parameters_OSCAR02_L21-190.json"%CurDir # parameters from other samples
    integrated_SVs_CNVs_dir = ValidationDir+"/integrated_SVs_CNVs"
    get_cmds_perSVade_SVs_and_CNVs_calls(VarCallDir, sorted_aligners, simple_repeats_file, sorted(set(df_reads_info_SVs.runID)), genome_Cparapsilosis, ProcessedDataDir_validation, list_cmds, None, gff_Cparapsilosis, SV_parms_file=SV_parms_file, integrated_SVs_CNVs_dir=integrated_SVs_CNVs_dir)

    # debug
    if len(list_cmds)>0: return

    # get the small variants integrated and annotated
    integrated_small_vars_dir = ValidationDir+"/integrated_small_vars"
    integrate_and_annotate_small_vars(VarCallDir, sorted_aligners,  sorted(set(df_reads_info_analysis_with_validation.runID)), genome_Cparapsilosis, None, list_cmds, None, gff_Cparapsilosis, simple_repeats_file, integrated_small_vars_dir=integrated_small_vars_dir)


def get_df_reads_info_analysis_extra_validation(ValidationDir, VarCallDir, df_reads_info_extra_validation, DataDir, CurDir, list_cmds, simple_repeats_file, homozygous_Corthopsilosis_srrs, window_size_SNP_densities, sorted_aligners, genome_Cparapsilosis, aligner_SNP_densities, ReadsDir, min_coverage, min_pct_reads_para, min_pct_covered):

    """Runs the QC and filtering of samples only for the extra validation datasets"""

    df_reads_info_analysis_extra_validation_file = "%s/df_reads_info_analysis_extra_validation.py"%ValidationDir
    if file_is_empty(df_reads_info_analysis_extra_validation_file):
        print("Getting extra validation QC and filtering...")

        # keep
        df_reads_info_extra_validation = df_reads_info_extra_validation.copy()

        # define ProcessedDataDir validation
        ProcessedDataDir_validation = "%s/processed_data"%ValidationDir; make_folder(ProcessedDataDir_validation)
        PlotsDir_validation = "%s/plots"%ValidationDir; make_folder(PlotsDir_validation)

        # run multiQC for all of them
        sorted_runs = sorted(set(df_reads_info_extra_validation.runID))
        run_multiqc_several_trimmed_reads("%s/multiqc_outdir_extra_runs"%(ValidationDir), sorted_runs, function_origin_zip=define_origin_zip_Cparapsilosis_popgen)

        # generate cov per windows, for window_size_SNP_densities, which is necessary for QC
        aligner_to_wsize_to_df_coverage_per_window_file = {}
        for aligner in sorted_aligners:

            cov_file = generate_df_coverage_per_window(list_cmds, aligner, window_size_SNP_densities, VarCallDir, sorted_runs, genome_Cparapsilosis, simple_repeats_file, ProcessedDataDir_validation, homozygous_Corthopsilosis_srrs)
            aligner_to_wsize_to_df_coverage_per_window_file.setdefault(aligner, {}).setdefault(window_size_SNP_densities, cov_file)

        # jobs to run
        if len(list_cmds)>0: return 
        
        # get filtered vars
        aligner_to_df_small_vars_filt_file = {}
        for aligner in sorted_aligners:
            aligner_to_df_small_vars_filt_file[aligner] = get_filtered_small_vars_df(ProcessedDataDir_validation, VarCallDir, sorted_runs, aligner, simple_repeats_file, homozygous_Corthopsilosis_srrs, aligner_to_wsize_to_df_coverage_per_window_file[aligner][window_size_SNP_densities], window_size_SNP_densities)

        # get a df with SNP / INDEL densities
        df_variant_densities_file = generate_variant_per_window_densities_df(aligner_to_df_small_vars_filt_file[aligner_SNP_densities], ProcessedDataDir_validation, genome_Cparapsilosis, aligner_SNP_densities, window_size_SNP_densities, sorted_runs, simple_repeats_file, aligner_to_wsize_to_df_coverage_per_window_file[aligner_SNP_densities][window_size_SNP_densities])

        # create df with added info about coverage and multiqc
        df_reads_info_extra_validation = get_df_reads_info_with_multiqc_and_coverage_data(df_reads_info_extra_validation, ProcessedDataDir_validation, sorted_aligners, homozygous_Corthopsilosis_srrs, df_variant_densities_file, genome_Cparapsilosis, PlotsDir_validation, aligner_SNP_densities, window_size_SNP_densities, multiqc_data_file="%s/multiqc_outdir_extra_runs/multiqc_output/multiqc_data/multiqc_fastqc.txt"%ValidationDir, skip_plots=True)

        # generate all the QC reports, getting df to be run in parallel
        df_reads_info_extra_validation = run_QC_and_filtering_samples(df_reads_info_extra_validation, homozygous_Corthopsilosis_srrs, ProcessedDataDir_validation, ReadsDir, min_coverage, min_pct_reads_para, min_pct_covered, PlotsDir_validation, aligner_to_df_small_vars_filt_file, None, replace=False, skip_plots=True)

        # filter and save
        df_reads_info_analysis_extra_validation = df_reads_info_extra_validation[(df_reads_info_extra_validation.possibly_triploid==False) & (df_reads_info_extra_validation.species_run=="para")].copy()
        save_object(df_reads_info_analysis_extra_validation, df_reads_info_analysis_extra_validation_file)

    return load_object(df_reads_info_analysis_extra_validation_file)
      


def run_QC_and_filtering_samples(df_reads_info, homozygous_Corthopsilosis_srrs, ProcessedDataDir, ReadsDir, min_coverage, min_pct_reads_para, min_pct_covered, PlotsDir, aligner_to_df_small_vars_filt_file, df_metadata, replace=False, skip_plots=False):

    """runs all the QC and filtering of all samples. It returns the filtered df"""

    df_reads_info_filt_file = "%s/df_reads_info_filtered_afterQC.py"%ProcessedDataDir
    df_reads_info_bad_file = "%s/df_reads_info_bad_runs.py"%ProcessedDataDir

    if file_is_empty(df_reads_info_filt_file) or file_is_empty(df_reads_info_bad_file) or replace is True:

        # keep
        print("running QC and filt...")
        df_reads_info = df_reads_info.copy()
        total_n_runs = len(df_reads_info)
        initial_runs_all = sorted(df_reads_info.runID)

        # init bad runs
        dfs_bad_runs = []

        # print filers
        print("min_coverage", min_coverage)
        print("min_pct_reads_para", min_pct_reads_para)
        print("min_pct_covered", min_pct_covered)

        # define dirs
        kraken_dir_all = "%s/kraken_runs"%ProcessedDataDir; make_folder(kraken_dir_all)

        # add info to df
        df_reads_info["sampleID"] = df_reads_info.apply(get_sampleID_from_runID, axis=1)
        df_reads_info["species_run"] = df_reads_info.runID.isin(homozygous_Corthopsilosis_srrs).map({True:"ortho", False:"para"})
        check_no_nans_series(df_reads_info["species_run"])
        df_reads_info["fraction_aneuploid_genome"] = 1- (df_reads_info.non_aneuploid_genome_len_gDNA / df_reads_info.whole_genome_len_gDNA)

        # filter out runs with a bad fastqc (in part manually curated)
        df_reads_info["bad_fastqc"] = (df_reads_info.per_sequence_gc_content=="fail") | (df_reads_info.per_base_n_content!="pass") | (df_reads_info.adapter_content!="pass") 
        df_bad_fastqc = df_reads_info[df_reads_info.bad_fastqc].copy()
        print("runs bad fastqc:\n", df_bad_fastqc[["species_run", "per_sequence_gc_content", "per_base_n_content", "adapter_content", "mean_median_cov_gDNA", "pct_cov_gDNA"]])
        run_kraken_subset_runs(kraken_dir_all, df_bad_fastqc, "bad_fastqc_runs", ReadsDir) # clarify what is there
        
        dfs_bad_runs.append(df_reads_info[df_reads_info.bad_fastqc].copy())
        df_reads_info = df_reads_info[~df_reads_info.bad_fastqc].copy() # filter
        print("%i/%i runs remaining after filtering bad fastqc"%(len(df_reads_info), total_n_runs))

        # print the coverage of C ortho cases
        print(df_reads_info[df_reads_info.runID.isin(homozygous_Corthopsilosis_srrs)][["mean_median_cov_gDNA", "pct_cov_gDNA"]]) # all have >10 coverage, so I keep them

        # plot expected coverage (based on initial # reads) vs the fraction of para assigned
        plots_dir_QC = "%s/QC_plots"%PlotsDir; make_folder(plots_dir_QC)
        for ylim_plot in [[-1, 100], [85, 101]]: 
            if skip_plots is False: plot_x_vs_y_df_reads_info(plots_dir_QC, "expected_coverage", "pct_reads_Cparapsilosis", ylim_plot, [min_pct_reads_para], [min_coverage], df_reads_info[df_reads_info.species_run=="para"], "all_runs_coverage", hue="source_data")

        # filter samples with insufficient coverage
        df_reads_info["bad_cov"] = (df_reads_info.species_run=="para") & ((df_reads_info.expected_coverage<min_coverage) | (df_reads_info.pct_reads_Cparapsilosis<min_pct_reads_para) | (df_reads_info.pct_cov_gDNA<min_pct_covered))
        df_bad_cov = df_reads_info[df_reads_info.bad_cov].copy()
        print("%i runs bad coverage:\n"%(len(df_bad_cov)), df_bad_cov[["expected_coverage", "pct_reads_Cparapsilosis", "pct_cov_gDNA", "mean_median_cov_gDNA", "mean_rel_cov_diff_windows_gDNA", "std_rel_cov_diff_windows_gDNA", "r_distTelomere_vs_cov_gDNA"]])

        run_kraken_subset_runs(kraken_dir_all, df_bad_cov, "bad_cov_runs", ReadsDir) # clarify what is there
        print("runs bad coverage low pct cov:\n", df_bad_cov[df_bad_cov.pct_reads_Cparapsilosis<85][[x for x in df_bad_cov.keys() if x.startswith("pct_reads_") and (any(df_bad_cov[x]>=5))]])

        # check the filters with some coverage
        #print_phenotypes_subset_runs(df_metadata, df_bad_cov[(df_bad_cov.expected_coverage>=35) & (df_bad_cov.source_data=="sequenced_here") & (df_bad_cov.pct_reads_Cparapsilosis>=85)], "samples bad coverage, cov >=25x")

        # filter bad cov
        dfs_bad_runs.append(df_reads_info[df_reads_info.bad_cov].copy())
        df_reads_info = df_reads_info[~df_reads_info.bad_cov].copy() # filter
        print("%i/%i runs remaining after filtering bad cov"%(len(df_reads_info), total_n_runs))

        # define the smiley
        max_r_smiley = -0.1
        df_reads_info["has_smiley"] = (df_reads_info.r_distTelomere_vs_cov_gDNA<max_r_smiley)

        # visualize to filter by coverage dispersion, also to verify the smiley
        max_std_diff_windows = 0.35
        if skip_plots is False: 
            plot_x_vs_y_df_reads_info(plots_dir_QC, "r_distTelomere_vs_cov_gDNA", "std_rel_cov_diff_windows_gDNA", None, [max_std_diff_windows], [max_r_smiley], df_reads_info[df_reads_info.species_run=="para"], "simley_info", hue="source_data")
            plot_x_vs_y_df_reads_info(plots_dir_QC, "r_distTelomere_vs_cov_gDNA", "IQR_rel_median_cov_gDNA", None, [], [max_r_smiley], df_reads_info[df_reads_info.species_run=="para"], "simley_info", hue="source_data")

        # filter by smiley face effect and coverage dispersion (only para)
        df_reads_info["bad_cov_dispersion"] = (df_reads_info.species_run=="para") & (df_reads_info.std_rel_cov_diff_windows_gDNA>max_std_diff_windows)
        print("runs bad coverage dispersion:\n", df_reads_info[df_reads_info.bad_cov_dispersion][["pct_reads_Cparapsilosis", "pct_cov_gDNA", "mean_median_cov_gDNA" , "std_rel_cov_diff_windows_gDNA", "r_distTelomere_vs_cov_gDNA"]])
        #print_phenotypes_subset_runs(df_metadata, df_reads_info[(df_reads_info.bad_cov_dispersion) & (df_reads_info.source_data=="sequenced_here")], "samples bad coverage dispersion")

        dfs_bad_runs.append(df_reads_info[df_reads_info.bad_cov_dispersion].copy())
        df_reads_info = df_reads_info[~df_reads_info.bad_cov_dispersion].copy() # filter
        print("%i/%i runs remaining after filtering bad cov dispersion"%(len(df_reads_info), total_n_runs))

        # record the smiley of the filt dfs
        print("samples with smiley:\n", sorted(set(df_reads_info[df_reads_info.has_smiley].index)), "len=", sum(df_reads_info.has_smiley))

        # investigate data for the AF filtering
        for x,y in [("homo_SNPs_per_kb_gDNA", "hetero_SNPs_per_kb_gDNA"), ("unknown_zygosity_SNPs_per_kb_gDNA", "hetero_SNPs_per_kb_gDNA"), ("hetero_SNPs_per_kb_gDNA", "gDNA_hetero_AF_median")]:
            if skip_plots is False: plot_x_vs_y_df_reads_info(plots_dir_QC, x, y, None, [], [], df_reads_info[df_reads_info.species_run=="para"], "investigate_AFs", hue="source_data")

        df_hetero_para = df_reads_info[(df_reads_info.species_run=="para") & (df_reads_info.hetero_SNPs_per_kb_gDNA>=0.1)]
        for x,y in [("gDNA_hetero_AF_median", "gDNA_hetero_AF_mean"), ("gDNA_hetero_AF_median", "fraction_aneuploid_genome"), ("gDNA_hetero_AF_median", "gDNA_hetero_AF_IQR")]:
            if skip_plots is False: plot_x_vs_y_df_reads_info(plots_dir_QC, x, y, None, [], [], df_hetero_para, "hetero_SNPS_investigate_AFs", hue="source_data")

        # We find some samples with strange heterozygous patterns. Plots to investigate them
        possible_hetero_samples = sorted(set(df_reads_info[(df_reads_info.species_run=="para") & ((df_reads_info.hetero_SNPs_per_kb_gDNA>=0.1) | (df_reads_info.gDNA_hetero_AF_median<0.4) | (df_reads_info.gDNA_hetero_AF_median>0.6))].runID))

        example_normal_samples = list(df_reads_info[(df_reads_info.species_run=="para") & ~(df_reads_info.runID.isin(possible_hetero_samples))].sort_values(by=["hetero_SNPs_per_kb_gDNA", "gDNA_hetero_AF_median"], ascending=False).iloc[0:2].runID)

        samples_investigate_hetero = example_normal_samples + possible_hetero_samples
        df_snps_hetero_samples = get_df_snps_specific_samples(aligner_to_df_small_vars_filt_file, samples_investigate_hetero, "%s/getting_AFs_hetero_samples"%ProcessedDataDir, replace=False)
        if skip_plots is False: plot_distribution_AFs_specific_samples(df_snps_hetero_samples, samples_investigate_hetero, "%s/investigate_strange_hetero_samples"%plots_dir_QC)

        # define the types of samples with bad distribution, if there are any, based on manual curation
        if len(possible_hetero_samples)>0:

            runs_strain_mix = {"OSCAR01_L20-028", "OSCAR01_L20-032", "OSCAR01_L20-053", "OSCAR01_L20-066", "OSCAR01_L20-080", "OSCAR01_L20-085", "OSCAR01_L21-084", "OSCAR01_L21-146", "OSCAR01_L21-148", "OSCAR02_L22-121"}
            runs_possibly_triploid = {"OSCAR01_L20-027", "OSCAR01_L20-048", "OSCAR02_L21-189", "OSCAR02_L22-440"}
            if sorted(runs_strain_mix.union(runs_possibly_triploid))!=possible_hetero_samples: raise ValueError("not classified samples")
            df_reads_info["strain_mix"] = df_reads_info.runID.isin(runs_strain_mix)
            df_reads_info["possibly_triploid"] = df_reads_info.runID.isin(runs_possibly_triploid)

        else:
            df_reads_info["strain_mix"] = df_reads_info.runID.isin(possible_hetero_samples)
            df_reads_info["possibly_triploid"] = df_reads_info.runID.isin(possible_hetero_samples)

        # print the stran coverage
        print(df_reads_info[df_reads_info.runID.isin(possible_hetero_samples)][["possibly_triploid", "strain_mix", "expected_coverage"]].sort_values(by=["strain_mix", "possibly_triploid"]))
        
        # get the kraken runs
        run_kraken_subset_runs(kraken_dir_all, df_reads_info[df_reads_info.possibly_triploid], "possibly_triploid", ReadsDir) 

        # filter strain mixes
        dfs_bad_runs.append(df_reads_info[df_reads_info.strain_mix].copy())
        df_reads_info = df_reads_info[~(df_reads_info.strain_mix)].copy()
        print("%i/%i runs remaining after filtering strain mixes"%(len(df_reads_info), total_n_runs))

        # investigate the stats for triploid samples
        for x,y in [("hetero_SNPs_per_kb_gDNA", "unknown_zygosity_SNPs_per_kb_gDNA"), ("hetero_SNPs_per_kb_gDNA", "gDNA_homo_AF_median")]:
            if skip_plots is False: plot_x_vs_y_df_reads_info(plots_dir_QC, x, y, None, [], [], df_reads_info[df_reads_info.species_run=="para"], "triploid_strains", hue="possibly_triploid")

        # get bad cov df
        df_reads_info_bad_runs = pd.concat(dfs_bad_runs).reset_index(drop=True)
        if initial_runs_all!=sorted(list(df_reads_info.runID) + list(df_reads_info_bad_runs.runID)):
            raise ValueError("the runs are not as expected")

        # save
        save_object(df_reads_info, df_reads_info_filt_file)
        save_object(df_reads_info_bad_runs, df_reads_info_bad_file)

    # print bad runs
    df_reads_info_bad = load_object(df_reads_info_bad_file)
    print("Filtered out runs from OSCAR: %i"%len(df_reads_info_bad))
    print(df_reads_info_bad[df_reads_info_bad.runID.apply(lambda x: x.startswith("OSCAR"))][["runID", "bad_fastqc", "bad_cov", "bad_cov_dispersion", "strain_mix", "possibly_triploid"]].reset_index(drop=True))

    print("ortho strains bad:")
    print(df_reads_info_bad[df_reads_info_bad.runID.isin(homozygous_Corthopsilosis_srrs)][["runID", "bad_fastqc", "bad_cov", "bad_cov_dispersion", "strain_mix", "possibly_triploid"]].reset_index(drop=True))


    # get the good runs
    return load_object(df_reads_info_filt_file)

def get_df_reads_info_with_multiqc_and_coverage_data(df_reads_info, ProcessedDataDir, sorted_aligners, homozygous_Corthopsilosis_srrs, df_variant_densities_file, genome_Cparapsilosis, PlotsDir, aligner_stats, window_size_SNP_densities, multiqc_data_file="auto", skip_plots=False):

    """Gets the df_reads_info with some adds"""

    df_reads_info_file = "%s/df_reads_info_with_multiqc_and_coverage_data_and_AFs_%s_%ibp.py"%(ProcessedDataDir, aligner_stats, window_size_SNP_densities)
    if file_is_empty(df_reads_info_file):
        print("adding multiqc and coverage info...")

        # add multiqc data
        if multiqc_data_file=="auto":
            multiqc_data_file = "%s/multiqc_outdir_all/multiqc_output/multiqc_data/multiqc_fastqc.txt"%ProcessedDataDir
        df_multiqc = pd.read_csv(multiqc_data_file, sep="\t")[["Sample", "per_sequence_gc_content", "per_base_n_content", "adapter_content"]]

        df_multiqc["runID"] = df_multiqc.Sample.apply(lambda s: "_".join(s.split()[0].split("_")[0:-1]))
        df_multiqc["rID"] = df_multiqc.Sample.apply(lambda s: s.split()[0].split("_")[-1])
        if len(df_multiqc)!=(len(set(df_multiqc.runID))*2): raise ValueError("there should be one multiqc for each read")

        df_multiqc = df_multiqc.set_index("runID")
        for fastq_f in ["per_sequence_gc_content", "per_base_n_content", "adapter_content"]:

            tag_to_number = {"pass":0, "warn":1, "fail":2}
            def get_worst_fastq_tag_one_runID(runID):

                # manually curated situations
                if fastq_f=="per_sequence_gc_content" and runID=="ERR263539": return "fail" # This has a strange GC content distribution

                # get the worst tag
                all_tags = list(df_multiqc.loc[runID, fastq_f].values)
                return sorted(all_tags, key=(lambda t: tag_to_number[t]))[-1]

            df_reads_info[fastq_f] = df_reads_info.runID.apply(get_worst_fastq_tag_one_runID)

        # add coverage data and AF data
        tmpdir = "%s.generating"%df_reads_info_file
        make_folder(tmpdir)

        df_SNP_densities_file = "%s/df_SNP_densities.py"%tmpdir
        if file_is_empty(df_SNP_densities_file):

            print("loading SNP df")
            df_variant_densities = load_object(df_variant_densities_file).rename(columns={"mediancov_1":"median_coverage"})

            print("add fields")
            df_SNP_densities = df_variant_densities[(df_variant_densities.type_var=="SNP") & (df_variant_densities.aligner==aligner_stats)].copy()
            df_SNP_densities["len_window"] = df_SNP_densities.end - df_SNP_densities.start

            df_SNP_densities["genome_type"] = (df_SNP_densities.chromosome==mtDNA_chrom).map({True:"mtDNA", False:"gDNA"})
            check_no_nans_series(df_SNP_densities.genome_type)

            df_SNP_densities["chrom_len"] = df_SNP_densities.chromosome.map(get_chr_to_len(genome_Cparapsilosis))
            df_SNP_densities["pos"] = df_SNP_densities.start + ((df_SNP_densities.len_window / 2).apply(int))
            df_SNP_densities["dist_to_chr_end"] = df_SNP_densities.chrom_len - df_SNP_densities.pos
            df_SNP_densities["dist_to_telomere"] = df_SNP_densities[["dist_to_chr_end", "pos"]].apply(min, axis=1)

            print("saving")
            final_fields = ["len_window", "dist_to_telomere", "nvars", "zygosity", "genome_type", "median_coverage_rel_to_gDNA", "AF_list", "fraction_repeats", "fraction_N_bases", "median_mappability", "percentcovered_1", "median_coverage", "runID", "chromosome", "pos"]
            save_object(df_SNP_densities[final_fields], df_SNP_densities_file)

        print("loading SNP df...")
        df_SNP_densities = load_object(df_SNP_densities_file)

        # define filters to calculate AF of SNP types
        max_frac_repeats = 0.05
        max_fraction_N_bases = 0.05
        min_median_mappability = 0.25

        min_rel_cov = 0.8
        max_rel_cov = 1.2
        min_n_vars = 10

        # plot different measures for all windows
        print("plotting")
        if skip_plots is False:
            for f,tsh_val in [("fraction_repeats", max_frac_repeats), ("fraction_N_bases", max_fraction_N_bases), ("median_mappability", min_median_mappability)]:

                fig = sns.displot(data=df_SNP_densities, x=f, height=2, aspect=1.2, col="genome_type", log_scale=(False, True))
                for ax in fig.axes.flat:
                    ax.axvline(tsh_val, ls="--", c="red", lw=.7)
                save_fig(fig, "QC_windows/%s_all_windows.pdf"%f, PlotsDir)


        # filter
        df_SNP_densities = df_SNP_densities[(df_SNP_densities.fraction_repeats<=max_frac_repeats) & (df_SNP_densities.fraction_N_bases<=max_fraction_N_bases) & (df_SNP_densities.median_mappability>=min_median_mappability)].copy()

        # cov filters
        if skip_plots is False:

            df_dp_plot = df_SNP_densities[(df_SNP_densities.zygosity.isin({"hetero", "homo"})) & (df_SNP_densities.median_coverage_rel_to_gDNA>0)].copy()
            fig = sns.displot(data=df_dp_plot, x="median_coverage_rel_to_gDNA", height=2, aspect=1.2, log_scale=(True, True), col="genome_type")
            for ax in fig.axes.flat:
                for x in [min_rel_cov, max_rel_cov]: ax.axvline(x, ls="--", c="red", lw=.7)

            save_fig(fig, "distribution_rel_cov_region.pdf", PlotsDir)

        # change index
        df_SNP_densities = df_SNP_densities.set_index(["runID"])

        # prepare inputs
        print("running in parallel")
        inputs_fn = list(df_reads_info.apply(lambda r: (r.copy(), df_SNP_densities.loc[[r.runID]].set_index("genome_type").copy(), min_rel_cov, max_rel_cov, min_n_vars), axis=1).values)

        print("run")
        with multiproc.Pool(multiproc.cpu_count()) as pool:
            df_reads_info = pd.DataFrame(pool.starmap(get_r_df_reads_info_with_added_info, inputs_fn, chunksize=1)) 
            pool.close()
            pool.terminate()

        # checks and save
        check_no_nans_in_df(df_reads_info)
        delete_folder(tmpdir)
        save_object(df_reads_info, df_reads_info_file)

    return load_object(df_reads_info_file)


def save_fig(fig, file_name, PlotsDir, no_pad=False):

    """plots figure"""

    full_name = "%s/%s"%(PlotsDir, file_name)
    os.makedirs(get_dir(full_name), exist_ok=True)

    print("saving %s"%full_name)

    # get args
    kwargs = {"bbox_inches":"tight"}
    if no_pad is True:
        kwargs["pad_inches"] = 0.02

    if file_name.endswith(".png"):
        kwargs["dpi"] = 500 # dots per inch

    # save
    fig.savefig(full_name, **kwargs)

def get_correct_tree(tree, min_support=0, assume_rooted=True):

    """Gets a tree with midpointRooting"""

    # load
    if type(tree)==str: tree = Tree(tree)
    else: tree = cp.deepcopy(tree)


    # add support to children
    for n in tree.get_children(): n.support = 100

    # collapse low support bracnhes
    for n in tree.traverse():

        if n.is_leaf() is False and assume_rooted is True and len(n.get_children())!=2: 
            print(n)
            raise ValueError("each node should have 2 childs")

        # for internal nodes, if the support is below the minimum, set each of the children as children of the parent
        if n.is_leaf() is False and len(n.get_ancestors())>0 and n.support<min_support: n.delete()


    return tree

def get_set_leafs_per_node(tree, min_nleafs=0):
    internal_nodes = [n for n in tree.traverse() if n.is_leaf() is False and len(n.get_leaf_names())>=min_nleafs]
    list_nodes = list(map(lambda n: tuple(sorted(n.get_leaf_names())),  internal_nodes))
    if len(list_nodes)!=len(set(list_nodes)): raise ValueError("each node should be unique")
    return set(list_nodes)

def get_pct_different_nodes_two_trees(tree1, tree2, min_nleafs=0, common_leaves=None):

    """Prints the distance between two trees"""

    tree1 = cp.deepcopy(tree1)
    tree2 = cp.deepcopy(tree2)

    # prune tree to keep the same
    if not common_leaves is None:
        tree1.prune(common_leaves)
        tree2.prune(common_leaves)

    # get sets of the two trees
    tree1_nodes_set = get_set_leafs_per_node(tree1, min_nleafs=min_nleafs) 
    tree2_nodes_set = get_set_leafs_per_node(tree2, min_nleafs=min_nleafs) 

    # calculate distances
    fraction_nodes_different1 = len(tree1_nodes_set.difference(tree2_nodes_set)) / len(tree1_nodes_set)
    fraction_nodes_different2 = len(tree2_nodes_set.difference(tree1_nodes_set)) / len(tree2_nodes_set)

    return fraction_nodes_different1*100, fraction_nodes_different2*100


def get_n_parsimony_informative_sites_from_iqtreefile(iqtree_file):

    """Gets the numbers of parsimony informative sites"""


    log_lines = open(iqtree_file, "r").readlines()
    lines_sites = [l.strip() for l in log_lines if l.startswith("Number of parsimony informative sites:")]
    if len(lines_sites)!=1: raise ValueError("there should be 1 lines")
    return int(lines_sites[0].strip().split()[-1])


def get_tree_diploid_mode_from_outdir(outdir_tree, outdir, outgroup_resampled_trees, reference_genome):

    """Gets a tree based on the outdir_tree, into outdir"""

    # rsync tree to outdir (only important files)
    print("rsyncing %s"%outdir_tree)
    if not os.path.isdir(outdir):

        # define tmp
        outdir_tmp = "%s_tmp"%outdir; delete_folder(outdir_tmp); make_folder(outdir_tmp)
        delete_folder(outdir)

        # sync only important files
        for base_dir, dir_names, file_names in os.walk(outdir_tree):
            base_dir_out = cp.deepcopy(base_dir).replace(outdir_tree, outdir_tmp); make_folder(base_dir_out)

            # make all dirs
            for d in dir_names: make_folder("%s/%s"%(base_dir_out, d))

            # sync only necessary files
            for f in file_names: 

                # skip some files, unnecessary
                if f in {"reference_genome.fasta"}: continue

                origin_f = "%s/%s"%(base_dir, f)
                dest_f = "%s/%s"%(base_dir_out, f)

                if file_is_empty(dest_f):
                    if (not "resamplings/resample_" in base_dir) or (f in {"tree_was_generated.txt", "iqtree_unroted.treefile", "iqtree_unroted.iqtree"}): rsync_file(origin_f, dest_f)

        # keep
        os.rename(outdir_tmp, outdir)

    # run the command to get the trees
    paths_table = "%s/paths_table.tab"%outdir
    outdir_tree_generation = "%s/tree_generation"%outdir
    
    print("getting tree...")
    if file_is_empty("%s/tree_generated.txt"%outdir_tree_generation): run_cmd("%s/scripts/tree_from_SNPs/scripts/get_tree.py --paths_table %s -o %s --mode diploid --threads 48 --min_coverage_pos 12 --iqtree_mem 20G --batch_mode --reference_genome %s --mode_fasta_gen python --outgroup_resampled_trees %s"%(ParentDir, paths_table, outdir_tree_generation, reference_genome, outgroup_resampled_trees), env="tree_from_SNPs_env")

    # getting parsimony informatoive sites
    print("getting parsimony sites")
    all_parsimony_informative_sites = [get_n_parsimony_informative_sites_from_iqtreefile("%s/tree_generation_logs/resample%i.iqtree"%(outdir_tree_generation, I)) for I in range(1, 100+1)]

    return Tree("%s/rooted_tree_mode=diploid_min_cov=12.nw"%outdir_tree_generation), "-".join([str(min(all_parsimony_informative_sites)), str(max(all_parsimony_informative_sites))])

def compare_trees_different_aligners_and_get_trees_df(DataDir, min_support, PlotsDir, ProcessedDataDir, reference_genome, skip_report_outgroups=False):

    """Reports information about different trees"""

    # define the dir
    trees_dir = "%s/SNP_trees"%DataDir
    processed_trees_dir = "%s/SNP_trees"%ProcessedDataDir; make_folder(processed_trees_dir)

    # init a dict with dict info
    tress_data_keys = ["outgroup", "aligner", "mode_tree", "tree", "tree_corrected", "nsites_parsimony"]
    trees_data = [] # outgroup, aligner, mode, tree, tree_corrected,  nsites_parsimony

    # load the outgroup-based trees
    print("loading outgroup-based trees...")
    ortho_outgroups = sorted({"SRR13448246", "SRR13448248", "SRR13448249", "SRR3547629", "SRR6509076"}) 
    for outg in ortho_outgroups:

        # get the diploid mode
        tree_dip, nparsimony_sites_dip = get_tree_diploid_mode_from_outdir("%s/mode-diploid_aligner-bwa_mem_outgroup-%s"%(trees_dir, outg), "%s/mode-diploid_aligner-bwa_mem_outgroup-%s_finished_tree"%(processed_trees_dir, outg), outg, reference_genome)
        trees_data.append([outg, "bwa_mem", "diploid", cp.deepcopy(tree_dip), get_correct_tree(tree_dip, min_support=min_support), nparsimony_sites_dip])

        # get the diploid homozygous mode
        outdir_tree_gen_homo = "%s/mode-diploid_homozygous_aligner-bwa_mem_outgroup-%s/tree_generation/generate_tree_from_SNPs"%(trees_dir, outg)
        tree_dip_homo = Tree("%s/iqtree_unroted.treefile"%outdir_tree_gen_homo)
        tree_dip_homo.set_outgroup(outg)
        trees_data.append([outg, "bwa_mem", "diploid_homo", cp.deepcopy(tree_dip_homo), get_correct_tree(tree_dip_homo, min_support=min_support), str(get_n_parsimony_informative_sites_from_iqtreefile("%s/iqtree_unroted.iqtree"%outdir_tree_gen_homo))])


    trees_data_df = pd.DataFrame(trees_data, columns=tress_data_keys)

    # define the para outgroups for each tree
    print("getting outgroup info")
    outgroups_dict = {}; Id = 0
    for I,r in trees_data_df.iterrows():

        if len(r.tree_corrected.get_children())!=2: raise ValueError("should be 2 trees")
        para_tree = [n for n in r.tree_corrected.get_children() if n.name!=r.outgroup][0]
        children_leafs = [tuple(sorted(set(n.get_leaf_names()))) for n in para_tree.get_children()]
        lens_children_leafs = [len(x) for x in children_leafs]
        print("%s-%s has a para root with %i children. These are their lengths: %s. Picking shortest as para root..."%(r.mode_tree, r.outgroup, len(children_leafs), lens_children_leafs))

        for para_outg_tuple in sorted(children_leafs):
            outgroups_dict[Id] = {"para_outg" : ",".join(sorted(para_outg_tuple)), "mode" : r.mode_tree, "ortho_outg" : r.outgroup}; Id+=1

    outgroups_df = pd.DataFrame(outgroups_dict).transpose()

    para_outg_to_len = outgroups_df.groupby("para_outg").apply(lambda df_o: len(df_o[["mode", "ortho_outg"]].drop_duplicates()))
    para_outg_to_short_name = pd.Series({x : "%s-%s (%i leafs)"%(x.split(",")[0], x.split(",")[-1], len(x.split(","))) for x in para_outg_to_len.index})
    if len(para_outg_to_short_name)!=len(set(para_outg_to_short_name)): raise ValueError("non unique short names")
    para_outg_to_short_to_long = pd.Series(para_outg_to_short_name.index, index=para_outg_to_short_name.values)
    para_outg_to_len.index = [para_outg_to_short_name[x] for x in para_outg_to_len.index]
    para_outg_to_len = para_outg_to_len.sort_values()

    # I find that there are two possible outgroups that are found in all diploid (156) or all diploid_homo (18 leafs) trees. Check if there are these nodes in all the trees
    if skip_report_outgroups is False:

        interesting_para_outgs_short = sorted(para_outg_to_len[para_outg_to_len==5].index)
        for outg_para_short in interesting_para_outgs_short:
            outg_para_leafs = set(para_outg_to_short_to_long[outg_para_short].split(","))

            for aligner in ["bwa_mem", "bowtie2_local", "hisat2_no_spliced"]:

                # print diploid_homo
                tree_diploid_homo = Tree("%s/mode-diploid_homozygous_aligner-%s_outgroup-none/tree_generation/generate_tree_from_SNPs/iqtree_unroted.treefile"%(trees_dir, aligner))
                print("%s-%s. Diploid_homo. Has para outgroup node: %s"%(outg_para_short, aligner, any([set(n.get_leaf_names())==outg_para_leafs for n in tree_diploid_homo.traverse()])))

                # print diploid
                bools_node_in_tree = []
                for I in range(1, 100+1):
                    tree_diploid_I = Tree("%s/mode-diploid_aligner-%s_outgroup-none/tree_generation/tmp/generate_tree_from_SNPs_resamplingHetSNPs/resamplings/resample_%i/iqtree_unroted.treefile"%(trees_dir, aligner, I))
                    bools_node_in_tree.append(any([set(n.get_leaf_names())==outg_para_leafs for n in tree_diploid_I.traverse()]))

                print("%s-%s. Diploid. %.2fpct of resamples have outgroup node."%(outg_para_short, aligner, (sum(bools_node_in_tree)/len(bools_node_in_tree)*100)))

    # I find one outhroup for each type of tree
    print("Getting para trees rooted...")
    mode_to_para_outgroup = {mode : para_outg_to_short_to_long[short_outg] for mode, short_outg in [("diploid", "ERR4367748-SRR8984834 (156 leafs)"), ("diploid_homo", "ERR4367748-SRR7898457 (18 leafs)")]} # this is based on the observation that these are good outgroups for each mode

    for aligner in ["bwa_mem", "bowtie2_local", "hisat2_no_spliced"]:

        # get the diploid mode
        tree_dip, nparsimony_sites_dip = get_tree_diploid_mode_from_outdir("%s/mode-diploid_aligner-%s_outgroup-none"%(trees_dir, aligner), "%s/mode-diploid_aligner-%s_outgroup-none_finished_tree"%(processed_trees_dir, aligner), mode_to_para_outgroup["diploid"], reference_genome)
        trees_data.append(['none', aligner, "diploid", cp.deepcopy(tree_dip), get_correct_tree(tree_dip, min_support=min_support), nparsimony_sites_dip])

        # get the diploid_homo mode
        outdir_tree_gen_homo = "%s/mode-diploid_homozygous_aligner-%s_outgroup-none/tree_generation/generate_tree_from_SNPs"%(trees_dir, aligner)
        tree_dip_homo = Tree("%s/iqtree_unroted.treefile"%outdir_tree_gen_homo)
        tree_dip_homo.set_outgroup([n for n in tree_dip_homo.traverse() if set(n.get_leaf_names())==set(mode_to_para_outgroup["diploid_homo"].split(","))][0])
        trees_data.append(['none', aligner, "diploid_homo", cp.deepcopy(tree_dip_homo), get_correct_tree(tree_dip_homo, min_support=min_support), str(get_n_parsimony_informative_sites_from_iqtreefile("%s/iqtree_unroted.iqtree"%outdir_tree_gen_homo))])

    # get trees data
    trees_data_df = pd.DataFrame(trees_data, columns=tress_data_keys)

    # add statistical info about each tree
    print("adding info to df...")
    def get_trees_data_df_r_with_added_stats(r):

        # define para outgroup depending on the mode
        para_outgroup_r = mode_to_para_outgroup[r.mode_tree]

        # define the parapsilosis tree
        if r.outgroup=="none": 
            para_tree = cp.deepcopy(r.tree)
            para_tree_corrected = cp.deepcopy(r.tree_corrected)

        else:

            # get the para tree
            para_tree = cp.deepcopy([n for n in r.tree.get_children() if n.name!=r.outgroup][0])
            para_tree.dist = 0.001
            para_tree_corrected = cp.deepcopy([n for n in r.tree_corrected.get_children() if n.name!=r.outgroup][0])
            para_tree_corrected.dist = 0.001

            # redefine the outgroup to be the one set by C. para
            outgroup_n =  [n for n in para_tree.traverse() if set(n.get_leaf_names())==set(para_outgroup_r.split(","))][0]
            para_tree.set_outgroup(outgroup_n)
            outgroup_n_corrected =  [n for n in para_tree_corrected.traverse() if set(n.get_leaf_names())==set(para_outgroup_r.split(","))][0]
            para_tree_corrected.set_outgroup(outgroup_n_corrected)


        # define the fraction of collapsed nodes
        nnodes_total = len([n for n in para_tree.traverse() if not n.is_leaf()])
        nnodes_collapsed = nnodes_total - len([n for n in para_tree_corrected.traverse() if not n.is_leaf()])
        r["pct_collapsed_nodes"] = (nnodes_collapsed / nnodes_total)*100
        #para_tree_corrected.render("%s/example_tree.pdf"%PlotsDir)

        # define the fraction of collapsed nodes that have a long brach len
        outg_nodes = [n for n in para_tree_corrected.get_children() if set(n.get_leaf_names())==set(para_outgroup_r.split(","))]
        if len(outg_nodes)!=1: raise ValueError("there should be 1 oytg_node")
        reference_branch_len = outg_nodes[0].dist

        nnodes_collapsed_long_branches = len([n for n in para_tree.traverse() if not n.is_leaf() and n.support<min_support and n.dist>=(reference_branch_len*0.25)])
        r["pct_collapsed_long_nodes"] = (nnodes_collapsed_long_branches / nnodes_total)*100

        # define the fraction of politomies that have a lot of diversity inside the nodes
        npolitomies = len([n for n in para_tree_corrected.traverse() if not n.is_leaf() and len(n.get_children())>2])
        npolitomies_diverse = len([n for n in para_tree_corrected.traverse() if not n.is_leaf() and len(n.get_children())>2 and any([len(c.get_leaf_names())>=3 for c in n.get_children()])])
        r["pct_politiomies_diverse"] = (npolitomies_diverse / nnodes_total)*100

        return r

    trees_data_df = trees_data_df.apply(get_trees_data_df_r_with_added_stats, axis=1)
    print(trees_data_df[["outgroup", "aligner", "mode_tree", "nsites_parsimony", "pct_collapsed_nodes", "pct_collapsed_long_nodes"]])

    # compare the topology of trees
    trees_data_df["unique_treeID"] = trees_data_df.outgroup + "_" + trees_data_df.aligner + "_" + trees_data_df.mode_tree
    sorted_unique_treeIDs = sorted(set(trees_data_df[trees_data_df.outgroup=="none"].unique_treeID))
    trees_data_df = trees_data_df.set_index("unique_treeID")

    for id_1, id_2 in itertools.combinations(sorted_unique_treeIDs, 2):

        tree1 = cp.deepcopy(trees_data_df.loc[id_1, "tree_corrected"])
        tree2 = cp.deepcopy(trees_data_df.loc[id_2, "tree_corrected"])

        # print distance
        dist1, dist2 = get_pct_different_nodes_two_trees(tree1, tree2, min_nleafs=3)
        print("%s vs %s. %.2fpct & %.2fpct different nodes"%(id_1.lstrip("none_"), id_2.lstrip("none_"), dist1, dist2))

    return trees_data_df



def old_fun():
    dkaghakhd

    # load the outgroup-based trees
    trees_dir = "%s/SNP_trees"%DataDir
    f_to_tree = {}
    f_to_tree_uncorrected = {}
    for f in os.listdir(trees_dir):
        tree = Tree("%s/%s/tree_generation/generate_tree_from_SNPs/iqtree_unroted.treefile"%(trees_dir, f))

        # define the tree type 
        if f.endswith("none"): outgroup =  tree.get_midpoint_outgroup()
        else: outgroup = f.split("-")[-1]

        # get tree
        tree.set_outgroup(outgroup)

        # define the total number of nodes
        total_nnodes = len(list(tree.traverse()))

        # keep
        f_to_tree_uncorrected[f] = cp.deepcopy(tree)

        # collapse low support nodes
        tree = get_correct_tree(tree, min_support=min_support)
    
        # print stats
        nnodes = len(list(tree.traverse()))
        ncollapsed_nodes = total_nnodes - nnodes
        print("%s, collapsed nodes=%i/%i (%.2fpct)"%(f, ncollapsed_nodes, total_nnodes, (ncollapsed_nodes/total_nnodes)*100 ))

        # get tree
        # ts = TreeStyle()
        # plots_dir_trees = "%s/trees_outgroup_or_midpoint_rooted"%PlotsDir; make_folder(plots_dir_trees)
        # tree.render(file_name="%s/%s.pdf"%(plots_dir_trees, f), tree_style=ts)

        # keep
        f_to_tree[f] = tree


    # compare the tree agreement with different outgroups
    sorted_outgroups = sorted([f.split("-")[-1] for f in os.listdir(trees_dir) if not f.endswith("none")])

    for outg1, outg2 in itertools.combinations(sorted_outgroups, 2):

        tree1 = cp.deepcopy(f_to_tree["aligner-bwa_mem_outgroup-%s"%outg1])
        tree2 = cp.deepcopy(f_to_tree["aligner-bwa_mem_outgroup-%s"%outg2])

        # keep only cpara nodes
        tree1_para = [n for n in tree1.get_children() if not n.name==outg1][0]
        tree2_para = [n for n in tree2.get_children() if not n.name==outg2][0]
        if set(tree1_para.get_leaf_names())!=set(tree2_para.get_leaf_names()): raise ValueError("the trees should have the same leaf names")

        # print distance
        dist1, dist2 = get_pct_different_nodes_two_trees(tree1_para, tree2_para)
        print("%s vs %s. %.2fpct & %.2fpct different nodes"%(outg1, outg2, dist1, dist2))


    # print stats about outgroup, as defined by the C. ortho outgroups
    outg_to_small_para_outgroup = {}
    for outg in sorted_outgroups:

        tree = cp.deepcopy(f_to_tree["aligner-bwa_mem_outgroup-%s"%outg])
        tree_para = [n for n in tree.get_children() if not n.name==outg][0]
        print(outg, "nrootgroups_para", len(tree_para.get_children()), "# leafs each subgroup", [len(n.get_leaf_names()) for n in tree_para.get_children()])

        # In the previous analyses I found that all para trees had two outgroups, one with 24 except SRR13448246, which I discard
        if outg!="SRR13448246": outg_to_small_para_outgroup[outg] = [n for n in tree_para.get_children() if len(n.get_leaf_names())==24][0]

    # make the comparisons between trees with outpgroups
    for outg1, outg2 in itertools.combinations(sorted(outg_to_small_para_outgroup.keys()), 2):

        t1 = cp.deepcopy(outg_to_small_para_outgroup[outg1])
        t2 = cp.deepcopy(outg_to_small_para_outgroup[outg2])

        leafs1 = set(t1.get_leaf_names())
        leafs2 = set(t2.get_leaf_names())
        
        pct_different_leafs1 = (len(leafs1.difference(leafs2)) / len(leafs1))*100
        pct_different_leafs2 = (len(leafs2.difference(leafs1)) / len(leafs2))*100

        if leafs1!=leafs2: raise ValueError("different leafs")

        dist1, dist2 = get_pct_different_nodes_two_trees(t1, t2)
        print("%s vs %s (only 24 outgroup para). %.2fpct & %.2fpct different nodes. %.2fpct %.2fpct different leafs"%(outg1, outg2, dist1, dist2, pct_different_leafs1, pct_different_leafs2))

    # In the previous analyses I realised that the 24 nodes have always the same numbers of leafs, so I try to see if there is a tree with them
    leafs_outgroup_according_to_ortho_rooting = cp.deepcopy(leafs1)
    for aligner in ["bwa_mem", "bowtie2_local", "hisat2_no_spliced"]:

        f = "aligner-%s_outgroup-none"%aligner

        # get trees
        unrooted_tree = Tree("%s/%s/tree_generation/generate_tree_from_SNPs/iqtree_unroted.treefile"%(trees_dir, f))
        unrooted_tree = get_correct_tree(unrooted_tree, min_support=min_support, assume_rooted=False)

        rooted_tree = cp.deepcopy(f_to_tree[f])

        # potential root
        print(aligner, "unrooted tree", any([set(n.get_leaf_names())==leafs_outgroup_according_to_ortho_rooting for n in unrooted_tree.traverse()]))
        print(aligner, "rooted tree (midpoint)", any([set(n.get_leaf_names())==leafs_outgroup_according_to_ortho_rooting for n in rooted_tree.traverse()]))

    # in the previous analysis I see that this leafs_outgroup_according_to_ortho_rooting can be found in all para trees. I set it as outgroup.
    aligner_to_tree = {}
    for aligner in ["bwa_mem", "bowtie2_local", "hisat2_no_spliced"]:

        # get tree
        final_tree = Tree("%s/aligner-%s_outgroup-none/tree_generation/generate_tree_from_SNPs/iqtree_unroted.treefile"%(trees_dir, aligner))
        outgroup_n = [n for n in final_tree.traverse() if set(n.get_leaf_names())==leafs_outgroup_according_to_ortho_rooting][0]
        final_tree.set_outgroup(outgroup_n)

        # collapse low support nodes
        final_tree = get_correct_tree(final_tree, min_support=0)

        # save
        plots_dir_trees = "%s/trees_parapsilosis_rooted"%PlotsDir; make_folder(plots_dir_trees)
        
        ts = TreeStyle()
        ts.mode = 'r'
        ts.arc_start = -180 # 0 degrees = 3 o'clock
        #ts.arc_span = 180    
        ts.show_branch_support = True
        ts.show_branch_length = True

        final_tree.render(file_name="%s/%s.pdf"%(plots_dir_trees, aligner), tree_style=ts)


        aligner_to_tree[aligner] = final_tree

    # compare the tree topology btw alogners
    for aln1, aln2 in [("bwa_mem", "bwa_mem"), ("bwa_mem", "bowtie2_local"), ("bwa_mem", "hisat2_no_spliced"), ("hisat2_no_spliced", "bowtie2_local")]:
        tree1 = get_correct_tree(cp.deepcopy(aligner_to_tree[aln1]), min_support=min_support)
        tree2 = get_correct_tree(cp.deepcopy(aligner_to_tree[aln2]), min_support=min_support)

        dist1, dist2 = get_pct_different_nodes_two_trees(tree1, tree2)
        print("%s vs %s. %.2fpct & %.2fpct different nodes"%(aln1, aln2, dist1, dist2))

        #largest_polytomy = max([max([len(n.get_children()) for n in t.traverse()]) for t in [tree1, tree2]])
        #rf_distance = tree1.robinson_foulds(tree2, expand_polytomies=True, polytomy_size_limit=largest_polytomy)[0]

    return aligner_to_tree


def get_phenotype_validation_bioprojects(df_sra_runs_Dec2024):

    """Filters to get target projects"""

    # define bioprojects from the literature
    target_bps = {"PRJNA592373", "PRJNA1011319",  # McTaggart_2024
                  "PRJNA844053", # Branco_2022
                  "PRJNA996760", # Brassington_2024
                  }

    # check missing
    missing_bps = target_bps.difference(set(df_sra_runs_Dec2024.BioProject))
    #print("missing_bps", missing_bps)

    # get samples sequenced in Toni's lab -> three projects
    SRA_fields_target = ["Run", 'Center Name', 'DATASTORE provider', 'Sample Name', 'collected_by', 'INSDC_center_name', 'Submitter_Id', 'Broker_name', 'INSDC_center_alias', 'culture_collection', 'COUNTRY', 'BioProject', 'Sample_name', 'isolate', 'strain']

    target_names = ["gabaldon", "genomic regulation", "supercomputing", "biomedical research", "recerca", "barcelona", "crg", "bsc", "irb"]
    def get_is_from_lab(r):
        return any([any([n in str(r[k]).lower() for n in target_names]) for k in SRA_fields_target])

    df_sra_runs_extra_GabaldonLab = df_sra_runs_Dec2024[df_sra_runs_Dec2024.apply(get_is_from_lab, axis=1)].copy()
    target_bps.update(set(df_sra_runs_extra_GabaldonLab.BioProject))

    print(df_sra_runs_extra_GabaldonLab[["Center Name", "Sample Name", 'Sample_name', 'isolate', 'strain', "Run", "BioProject"]], "\n", sorted(set(df_sra_runs_extra_GabaldonLab.BioProject))) #  ['PRJNA516045', 'PRJNA784379', 'PRJNA921967']

    return target_bps


def fetch_sra_runs_metadata(accessions, threads, tmpdir):

    """get biosample metadata of accessions"""

    # write accesiosn to tmpdir
    make_folder(tmpdir)
    sra_runs_file = "%s/sra_runs.txt"%tmpdir
    open(sra_runs_file, "w").write("\n".join(accessions))

    # get samplename
    run_info_df_file = "%s/run_info.csv"%tmpdir
    #fields = ["Run", "SampleName", "Sample", "BioSample", "Name"]
    if file_is_empty(run_info_df_file): 
        print("getting sra data")
        remove_file(run_info_df_file)

        run_cmd("epost -db sra -input %s -format acc | efetch -db sra -format runinfo > %s"%(sra_runs_file, run_info_df_file), env="Candida_mine_env")

        #run_cmd("cat %s | join-into-groups-of 500 | xargs -n 1 sh -c 'epost -db biosample -id \"$0\" -format acc | elink -target sra | efetch -db sra -format runinfo -mode xml | xtract -pattern Row -def \"NA\" -element %s' > %s"%(biosamples_file, " ".join(fields), run_info_df_file))

    run_info_df = pd.read_csv(run_info_df_file, sep=",")

    #missing_accessions = set(accessions).difference(set(run_info_df.Run))
    #print("Missing acccessions: %s (%i)"%(missing_accessions, len(missing_accessions)))
    
    return run_info_df

def get_mic_float_SRA(x):

    if str(x) in {"nan", "-1.0", "*****", "-"} or x==" ": return -1
    elif len(str(x).split("–"))==2 or len(str(x).split("-"))==2:
        split_x = re.split("–|-", x)
        conc1 = float(split_x[0])
        conc2 = float(split_x[1])
        diff_conc = conc2 - conc1
        return (conc1 + (diff_conc/2)) 

    elif all(map(lambda c: c.isdigit() or c==".", str(x))): return float(x)
    elif str(x).startswith(">"): return float(x.split(">")[1])*2
    elif str(x).startswith("≤"): return float(x.split("≤")[1])/2
    else: raise ValueError("invalid mic: '%s'"%x)

def get_df_phenotypes_validation(df_sra_runs_validation, ManualDataDir, ProcessedDataDir):

    """loads all the phenotypes of validation"""

    # go through each paper
    info_data = [("McTaggart_2024", ["PRJNA592373", "PRJNA1011319"]),
                 ("Branco_2022", ["PRJNA844053"]),
                 ("Brassington_2024", ["PRJNA996760"]),
                 ("GomezMolero_2020", ["PRJNA516045"]),
                 ("Daneshnia_2023", ["PRJNA921967"]),
                 ("Arastehfar_2022", ["PRJNA784379"])




    ]

    df_phenotypes_all = pd.DataFrame()
    for paper_name, bioprojects in info_data:
        #print(paper_name)

        # get df
        df_bps = df_sra_runs_validation[df_sra_runs_validation.BioProject.isin(bioprojects)].copy()

        # define a dataframe with the phenotypes
        if paper_name=="McTaggart_2024":

            # load tables
            df_SRAinfo = pd.read_excel(ManualDataDir+"/McTaggart_2024_TableS1.xlsx")[["Isolate", "Run"]].copy()
            df_phenotypes = pd.read_excel(ManualDataDir+"/McTaggart_2024_TableS2.xlsx", skiprows=2).iloc[0:3].copy().drop(["Unnamed: 1"], axis=1).set_index("Isolate").transpose().reset_index()[["index", "FLU MIC"]].rename(columns={"FLU MIC":"FLC_MIC", "index":"Isolate"}).reset_index(drop=True)

            # check
            if len(df_SRAinfo)!=len(df_phenotypes):
                raise ValueError("should be same len")

        elif paper_name=="Branco_2022":

            # load tables
            df_phenotypes =  get_tab_as_df_or_empty_df(ManualDataDir+"/Branco_2022_Table1.tab")
            df_SRAinfo = df_bps[["strain", "Run"]].copy().rename(columns={"strain":"Isolate"})

        elif paper_name=="Brassington_2024":

            # load tables
            df_phenotypes =  pd.read_excel(ManualDataDir+"/Brassington_2024_MICs.xlsx")
            df_SRAinfo = df_bps[["Run", "Library Name"]].copy().rename(columns={"Library Name":"Isolate"})

            # filter one wrong
            df_phenotypes = df_phenotypes[df_phenotypes.Isolate!="NRZ-2022-0386"].copy()

            # keep some
            df_SRAinfo = df_SRAinfo[df_SRAinfo.Isolate.isin(set(df_phenotypes.Isolate))].copy()
            if set(df_SRAinfo.Isolate)!=set(df_phenotypes.Isolate):
                raise ValueError("should be the same")

        elif paper_name=="GomezMolero_2020":

            # load tables
            df_SRAinfo = df_bps[["Run", "Sample Name"]].copy().rename(columns={"Sample Name":"Isolate"})
            target_idxs = [0, 5, 6, 7, 8, 9, 10]
            df_phenotypes =  pd.read_excel(ManualDataDir+"/GomezMolero_2020_Table2.xlsx").loc[target_idxs].copy().set_index("Strain acronym").transpose().iloc[0:5].copy().reset_index(drop=True).rename(columns={"Deposition number":"Isolate"})

            # checks
            if set(df_SRAinfo.Isolate)!=set(df_phenotypes.Isolate):
                raise ValueError("should be the same")

        elif paper_name=="Daneshnia_2023":

            # load tables
            df_SRAinfo = df_bps[["Run", "Sample Name"]].copy().rename(columns={"Sample Name":"Isolate"})
            df_phenotypes =  pd.read_excel(ManualDataDir+"/Daneshnia_2023_MICs.xlsx").drop(["StrainName"], axis=1)

        elif paper_name=="Arastehfar_2022":

            # load tables
            df_SRAinfo = df_bps[["Run", "Library Name"]].copy().rename(columns={"Library Name":"Isolate"})
            df_phenotypes =  pd.read_excel(ManualDataDir+"/Arastehfar_2022_MICs.xlsx")
            
            # filter
            df_phenotypes = df_phenotypes[df_phenotypes.Isolate.isin({35, 37})].copy()
            df_phenotypes["Isolate"] = df_phenotypes.Isolate.apply(str) + "_FNZ"
            df_phenotypes = df_phenotypes[df_phenotypes.Species.apply(lambda x: "parapsilosis" in x)].drop(["Species"], axis=1).copy()

        else:
            raise ValueError("invalid paper name")

            for f in df_bps.keys():
                if not all(pd.isna(df_bps[f])):
                    print(df_bps[f])


        # merge
        df_phenotypes = df_phenotypes.merge(df_SRAinfo, on="Isolate", how="left", validate="one_to_one")

        # discard one malformatted
        wrong_SRRs = {"SRR25837524" # invalid dataset with S as MIC
                     }

        df_phenotypes = df_phenotypes[~df_phenotypes.Run.isin(wrong_SRRs)].copy()

        # get data info
        run_info_df = fetch_sra_runs_metadata(sorted(set(df_phenotypes.Run)), 4, "%s/getting_srrs_for_%s"%(ProcessedDataDir, paper_name))

        # discard some wrong studies
        wrong_studies = {"SRP233635" # C. auris genomes
                        }

        valid_srrs_by_study = set(run_info_df[~run_info_df.SRAStudy.isin(wrong_studies)].Run)
        df_phenotypes = df_phenotypes[df_phenotypes.Run.isin(valid_srrs_by_study)].copy().reset_index(drop=True)

        # change data
        for f in ["FLC_MIC", "VRC_MIC", "POS_MIC", "ITR_MIC", "MIF_MIC", "CAS_MIC", "5FC_MIC", "ANI_MIC", "AMB_MIC"]:
            if f in df_phenotypes.keys(): df_phenotypes[f] = df_phenotypes[f].apply(get_mic_float_SRA)
            else: df_phenotypes[f] = -1

        # check the exoected SRRs
        missing_SRRs = set(df_phenotypes.Run).difference(set(df_bps.Run))
        if len(missing_SRRs)>0:
            raise ValueError("missing_SRRs: %s"%missing_SRRs)

        # add
        df_phenotypes["paper_name"] = paper_name

        # add
        df_phenotypes_all = pd.concat([df_phenotypes_all, df_phenotypes]).reset_index(drop=True)
        check_no_nans_in_df(df_phenotypes_all)

    # checks
    if len(df_phenotypes_all)!=len(set(df_phenotypes_all.Run)):
        print( len(df_phenotypes_all), len(set(df_phenotypes_all.Run)))
        raise ValueError("Run should be unique")

    return df_phenotypes_all


def get_df_sra_runs_validation(df_sra_runs_Dec2024, df_sra_runs):

    """gets the SRA runs for the validation"""

    df_sra_runs_Dec2024 = df_sra_runs_Dec2024.copy()
    df_sra_runs = df_sra_runs.copy()

    # get the bioprojects
    validation_BPs = get_phenotype_validation_bioprojects(df_sra_runs_Dec2024)

    # check SRRs are the same regardless
    for bp in set(df_sra_runs.BioProject).intersection(validation_BPs):

        srrs_sra = set(df_sra_runs[df_sra_runs.BioProject==bp].Run)
        srrs_sra_extra = set(df_sra_runs_Dec2024[df_sra_runs_Dec2024.BioProject==bp].Run)
        if len(srrs_sra_extra)==0: continue
        if srrs_sra!=srrs_sra_extra:
            print("extra runs old", srrs_sra.difference(srrs_sra_extra))
            print("extra runs extra", srrs_sra_extra.difference(srrs_sra))
            raise ValueError("The runs for bp %s are not equal"%bp)

    # get the df validation
    df_sra_runs_validation = pd.DataFrame()
    for bp in validation_BPs:

        if bp in set(df_sra_runs.BioProject):
            df_add = df_sra_runs

        elif bp in set(df_sra_runs_Dec2024.BioProject):
            df_add = df_sra_runs_Dec2024

        else:
            raise ValueError("invalid bp")

        df_sra_runs_validation  = pd.concat([df_sra_runs_validation, df_add[df_add.BioProject==bp].copy()]).reset_index(drop=True)

    # checks
    if set(df_sra_runs_validation.BioProject)!=validation_BPs:
        raise ValueError("should be equal")

    if len(df_sra_runs_validation)!=len(set(df_sra_runs_validation.Run)):
        raise ValueError("run should be unique")

    return df_sra_runs_validation


def get_homozygous_Corthopsilosis_srrs(ManualDataDir):

    """Gets a set of SRRs of C. orthopsilosis homozygous strains, as defined in Del Olmo 2023. 
    
    Data: 
        - DelOlmo 2023: PRJNA767198
        - Mixao 2023: PRJNA520893
        - Riccombeni et al. 2012: PRJNA431439 (inferred by looking at SRA)
        - Schroder et al. 2016: PRJNA322245
    """
    print("Getting C. orthopsilosis SRRs")

    # get homozygote strains
    df_DelOlmo2023 = pd.read_excel("%s/DelOlmo2023_DataS4.xlsx"%ManualDataDir, skiprows=3, skipfooter=8)
    df_DelOlmo2023 = df_DelOlmo2023[df_DelOlmo2023["Heteroz Variants/kb"]<1]
    df_DelOlmo2023["Strain"] = df_DelOlmo2023.Strain.apply(lambda x: str(x).lstrip().rstrip())
    df_DelOlmo2023["Source"] = df_DelOlmo2023.Source.apply(lambda x: str(x).lstrip().rstrip())

    # load the SRA data
    df_SRA = pd.read_csv("%s/SraRunTable_Cortho_WGS_10112023.csv"%ManualDataDir, sep=",")
    df_SRA = df_SRA[df_SRA["Assay Type"]=="WGS"]

    # map each source study to the bioproject
    source_to_BioProject = {"Mixão et al. 2023":"PRJNA520893", "Riccombeni et al. 2012":"PRJNA431439", "Schroder et al. 2016":"PRJNA322245", "This study":"PRJNA767198"}

    # get each srr for each strain
    strain_to_srr = {}
    for I,r in df_DelOlmo2023.iterrows():
        df_SRA_s = df_SRA[df_SRA.BioProject==source_to_BioProject[r.Source]]
        df_SRA_s = df_SRA_s[[k for k in df_SRA_s.keys() if not all(pd.isna(df_SRA_s[k]))]]

        # get the strain
        if r.Source in {"Mixão et al. 2023", "Riccombeni et al. 2012", "This study"}: df_strain = df_SRA_s[df_SRA_s.strain.apply(lambda x: x.replace(" ", "").lstrip().rstrip())==r.Strain]
        elif r.Source=="Schroder et al. 2016": df_strain = df_SRA_s[df_SRA_s["Sample Name"]=="Sample %s"%(r.Strain.replace("s", ""))]
        else: raise ValueError("invalid")

        # debug
        if len(df_strain)!=1: 
            print("\n", df_SRA_s.strain, "\n", r.Strain)
            print(df_strain.keys())
            for k in ["Library Name", "Sample Name"]: print(df_SRA_s[k])
            raise ValueError("df_strain should be 1. %s"%df_strain)

        # keep
        strain_to_srr[str(r.Strain)] = df_strain.iloc[0].Run

    return set(strain_to_srr.values())


def run_cmd_tuple(cmd_tuple): 
    print("running job %s"%(cmd_tuple[0]))
    run_cmd(cmd_tuple[1], env=cmd_tuple[2])


def rgb_to_hex(rgb):

    # Helper function to convert colour as RGB tuple to hex string
    rgb = tuple([int(255*val) for val in rgb])
    hex_val = '#%02x%02x%02x'%(rgb[0], rgb[1], rgb[2])

    if len(hex_val)!=7: raise ValueError("%s is not valid"%hex_val)

    return hex_val



def find_nearest(a, a0):
    """Element in nd array `a` closest to the scalar value `a0`"""
    
    # Debug elements that are inf
    if a0 not in [np.inf, -np.inf]:
        a = np.array(a)
        idx = np.abs(a - a0).argmin()
        closest_in_a = a.flat[idx]
        
    elif a0==np.inf:
        closest_in_a = max(a)
        
    elif a0==-np.inf:
        closest_in_a = min(a)        

    return closest_in_a

def get_value_to_color(values, palette="mako", n=100, type_color="rgb", center=None):

    """TAkes an array and returns the color that each array has. Checj http://seaborn.pydata.org/tutorial/color_palettes.html"""

    # get the colors
    colors = sns.color_palette(palette, n)

    # change the colors
    if type_color=="rgb": colors = colors
    elif type_color=="hex": colors = [rgb_to_hex(c) for c in colors]
    else: raise ValueError("%s is not valid"%palette)

    # if they are strings
    if type(list(values)[0])==str:

        palette_dict = dict(zip(values, colors))
        value_to_color = palette_dict

    # if they are numbers
    else:

        # map eaqually distant numbers to colors
        if center==None:
            min_palette = min(values)
            max_palette = max(values)
        else: 
            max_deviation = max([abs(fn(values)-center) for fn in [min, max]])
            min_palette = center - max_deviation
            max_palette = center + max_deviation

        all_values_palette = list(np.linspace(min_palette, max_palette, n))
        palette_dict = dict(zip(all_values_palette, colors))

        # get value to color
        value_to_color = {v : palette_dict[find_nearest(all_values_palette, v)] for v in values}

    return value_to_color, palette_dict

def get_cladeID_to_color(sorted_clades):

    """cladeID to color"""

    real_clades_float = [float(c) for c in sorted_clades if c!=-1]
    real_clades_str = list(map(str, real_clades_float))

    if len(real_clades_str)<=10: palette = "tab10"
    else: palette = "tab20"
    cladeID_to_color = get_value_to_color(real_clades_str, palette=palette, n=len(real_clades_str), type_color="hex")[0]
    cladeID_to_color = {float(cID) : c for cID,c in cladeID_to_color.items()}
    cladeID_to_color[-1] = "gray"

    return cladeID_to_color


def get_tree_withcladeID(tree, leaf_to_clade):

    """Gets a tree and adds the clade"""

    tree = cp.deepcopy(tree)

    for n in tree.traverse():

        # leafs are easy
        if n.is_leaf(): n.cladeID = leaf_to_clade[n.name]

        # get the clades of the leaf names
        else:
            leafs_clades = {leaf_to_clade[l] for l in n.get_leaf_names()}
            if len(leafs_clades)==1: n.cladeID = next(iter(leafs_clades))
            else: n.cladeID = -1

    return tree

def get_tree_ordered_clades(tree, leaf_to_clade):

    """Change the tree to have the input clades"""

    # clades 1 and 2 should be swapped
    clade1_leaf = leaf_to_clade[leaf_to_clade==1].index[0]
    clade2_leaf = leaf_to_clade[leaf_to_clade==2].index[0]
    ancestor_clades_12 = tree.get_common_ancestor([clade1_leaf, clade2_leaf])

    clade2_node, clade1_node = ancestor_clades_12.get_children()
    ancestor_clades_12.children = [clade1_node, clade2_node]


    # clades 4 and 5 should be swapped
    clade4_leaf = leaf_to_clade[leaf_to_clade==4].index[0]
    clade5_leaf = leaf_to_clade[leaf_to_clade==5].index[0]
    ancestor_clades_45 = tree.get_common_ancestor([clade4_leaf, clade5_leaf])

    clade5_node, clade4_node = ancestor_clades_45.get_children()
    ancestor_clades_45.children = [clade4_node, clade5_node]

    return tree



def get_cladeID_to_color_manual(): 

    """gets clades of Bergin"""

    colors_d = {1: "tab:blue",
                1.1: "tab:cyan",
                2: "olive",
                3: "tab:brown",
                4: "tab:red",
                4.1: "tab:purple",
                4.2: "tab:pink",
                4.3: "pink",
                5: "tab:orange", 
                -1 : "gray"}

    return {clade: mcolors.to_hex(color) for clade, color in colors_d.items()}



def get_stacked_barplot_noLegend(df, xfield, yfield, hue, sorted_xvalues, sorted_huevalues, palette, linewidth_edges=None, field_annot=None, hue_to_annot_color=None):

    """It plots a stacked barplot and returns the ax with it"""

    # redefine the df
    df = cp.deepcopy(df)
    df = df[[xfield, yfield, hue, field_annot]].copy()

    # redefine
    if field_annot is None:
        df["annot"] = ""
        field_annot = "annot"

    if hue_to_annot_color is None:
        hue_to_annot_color = {h:"black" for h in sorted_huevalues}


    # reverse
    sorted_huevalues = list(reversed(sorted_huevalues))

    # reformat df to inclue all hue vals for all x vals
    all_x_and_hue = set(df.apply(lambda r: (r[xfield], r[hue]), axis=1))

    for xval in sorted_xvalues:
        for hueval in sorted_huevalues:

            # in not
            if not hueval in hue_to_annot_color: hue_to_annot_color[hueval] = "black"

            # if not existing, add a 0
            if (xval, hueval) not in all_x_and_hue:
                df_dict = {(xval, hueval): {xfield : xval, hue: hueval, yfield: 0.0, field_annot:""}}
                df = pd.concat([df, pd.DataFrame(df_dict).transpose()])

    # sort by hue and the sorted xvalues
    xval_to_sortingVal = dict(zip(sorted_xvalues, range(len(sorted_xvalues))))
    df["xfield_sortingVal"] = df[xfield].apply(lambda xval: xval_to_sortingVal[xval])
    df = df.sort_values(by=[hue, "xfield_sortingVal"])

    # define the range y
    y_range = max(df.groupby(xfield).apply(lambda df: sum(df[yfield])))

    # go through ordered hue values
    for I, hue_val in enumerate(sorted_huevalues):
        color = palette[hue_val]

        # get the df hue
        df_hue = df[df[hue]==hue_val].copy()

        # get the array of bottom
        if I==0: bottom = [0]*len(df_hue)
        else: bottom = np.add(bottom, previous_Ys).tolist()

        # create the vars
        ax = plt.bar(df_hue[xfield], df_hue[yfield], color=color, bottom=bottom, linewidth=linewidth_edges) # edgecolor='white', width=barWidth

        # add the annotations
        for I, (x,y,annot_val) in enumerate(df_hue[[xfield, yfield, field_annot]].values):
            plt.text(x, bottom[I]+y-(y_range*0.06), annot_val, c=hue_to_annot_color[hue_val], ha="center", fontsize=9)

        # keep for the next 
        previous_Ys = list(df_hue[yfield])

    return ax


def map_hospital_short_to_very_short():

    """maps hospitals"""

    return {'COMPLEJO ASISTENCIAL UNIVERSITARIO DE BURGOS (Burgos)': 'Burgos', 

            'COMPLEJO H.ARIO VIRGEN MACARENA (Sevilla)': 'Sevilla - V. Macarena', 
            'H. SAN JUAN DE DIOS DEL ALJARAFE (Sevilla)': 'Sevilla - S. Aljarafe',
            'H. VIRGEN DEL ROCIO (Sevilla)': 'Sevilla - V. Rocio',

             'H. 12 DE OCTUBRE (Madrid)': 'Madrid - 12 Oct',
             'H. CENTRAL DE LA DEFENSA GOMEZ ULLA (Madrid)': 'Madrid - C. Defensa',
             'H. CLINICO SAN CARLOS (Madrid)': 'Madrid - S. Carlos',
             'H. DE MOSTOLES (Madrid)': 'Madrid - Mostoles',
             'H. PUERTA DE HIERRO (MAJADAHONDA) (Madrid)': 'Madrid - P. Hierro',

             'H. CLINICO UNIVERS DE VALLADOLID (Valladolid)': 'Valladolid',

             'H. DE SANT JOAN DE DEU (Barcelona)': 'Barcelona - J. Deu',
             "H. UNIVERSITARI VALL D'HEBRON (Barcelona)": 'Barcelona - V. Hebron',
             'H. UNIVERSITARIO DE BELLVITGE (Barcelona)': 'Barcelona - Bellvitge',

             'H. UNIVERS MARQUES DE VALDECILLA (Santander)': 'Santander',
             'LABORATORIO REFERENCIA CATALUNYA (EL PRAT DE LLOBREGAT) (Prat de Llobregat)': 'Barcelona - MRC'}


def plot_clade_distribution_per_hospitals(df_metadata,  PlotsDir, samples_discard_run1):

    """Gets per hospitals"""

    # filter
    df_metadata = df_metadata[~(df_metadata.runID.isin(samples_discard_run1)) & (df_metadata.possibly_triploid==False) & (df_metadata.source_data=="sequenced_here")].copy()

    # clade visuals
    cladeID_to_color = get_cladeID_to_color_manual()
    sorted_clades = sorted(cladeID_to_color)

    def format_val(x):
        if float(x)==int(x): return str(int(x))
        else: return str(x)


    # add the short hospital
    df_metadata["hospital_reduced"] = df_metadata.hospital_short.map(map_hospital_short_to_very_short())
    check_no_nans_series(df_metadata["hospital_reduced"])

    for hospital in sorted(set(df_metadata["hospital_reduced"])):
        df_metadata_h = df_metadata[df_metadata.hospital_reduced==hospital].copy()

        # get a pie chart
        df_plot = pd.DataFrame({"n_strains": df_metadata_h.groupby("cladeID_based_on_Bergin2022").apply(len)}).reset_index()
        df_plot["clade_idx"] = df_plot.cladeID_based_on_Bergin2022.map(dict(zip(sorted_clades, range(len(sorted_clades)))))
        df_plot = df_plot.sort_values(by=["clade_idx"])


        fig = plt.figure(figsize=(1, 1))
        #plt.pie(df_plot.n_strains, labels=df_plot.cladeID_based_on_Bergin2022.apply(format_val), colors=[cladeID_to_color[c] for c in df_plot.cladeID_based_on_Bergin2022])
        plt.pie(df_plot.n_strains, labels=df_plot.n_strains, colors=[cladeID_to_color[c] for c in df_plot.cladeID_based_on_Bergin2022], textprops={'fontsize': 8})

        ax = fig.axes[0]
        ax.set_title(hospital.replace(" - ", "\n"), fontsize=10)

        save_fig(fig, "pie_charts_hospitals/%s.pdf"%hospital, PlotsDir, no_pad=True)



def plot_clade_distribution_per_countries(df_metadata, PlotsDir, samples_discard_run1):

    """plots tree with clades per country"""

    # filter
    df_metadata = df_metadata[~(df_metadata.runID.isin(samples_discard_run1)) & (df_metadata.possibly_triploid==False) & (df_metadata.country_collection!="unknown")].copy()

    # get countires sorted by longitude
    sorted_countries = ['Canada', 'USA',  'Puerto Rico'] + ['Portugal', 'Spain', 'UK', 'France', 'Germany'] + ['Turkey', 'Kuwait', 'Malaysia', 'China', 'Korea']

    # get theme
    cladeID_to_color = get_cladeID_to_color_manual()
    sorted_clades = sorted(cladeID_to_color)

    # get a stacked barplot per country
    country_to_n_strains = df_metadata.groupby("country_collection").apply(len)
    df_plot = pd.DataFrame({"n_strains" : df_metadata[["cladeID_based_on_Bergin2022", "runID", "country_collection"]].copy().groupby(["cladeID_based_on_Bergin2022", "country_collection"]).apply(lambda df: len(df) / country_to_n_strains[df.country_collection.iloc[0]]) * 100}).reset_index()
    df_plot = df_plot[df_plot.country_collection.isin(sorted_countries)].copy()

    # init barplot
    sns.set_theme(style="white")
    fig = plt.figure(figsize=(5, 3))
    get_stacked_barplot_noLegend(df_plot, "country_collection", "n_strains", "cladeID_based_on_Bergin2022", sorted_countries, sorted_clades, cladeID_to_color)

    ax = fig.axes[0]
    ax.set_xticklabels(["%s (%i)"%(c, country_to_n_strains[c]) for c in sorted_countries], rotation=90)
    ax.set_ylabel("% strains per clade")

    sns.despine()
    save_fig(fig, "Fig_Dataset_country_bars.pdf", PlotsDir)


    # plot the number of strains per country
    df_plot_nstrains = pd.DataFrame({"n_strains" : df_metadata.copy().groupby(["country_collection", "source_data"]).apply(len)}).reset_index()
    df_plot_nstrains["country_idx"] = df_plot_nstrains.country_collection.map(dict(zip(sorted_countries, range(len(sorted_countries)))))
    df_plot_nstrains = df_plot_nstrains.sort_values(by=["country_idx"])


    fig = plt.figure(figsize=(5, 1))

    get_stacked_barplot_noLegend(df_plot_nstrains, "country_collection", "n_strains", "source_data", sorted_countries, ["sequenced_here", "SRA"], {"sequenced_here":"black", "SRA":"gray"}, linewidth_edges=0)

    ax = fig.axes[0]


    ax.set_xticklabels(sorted_countries, rotation=90)
    ax.set_ylabel("# strains")
    ax.set_xlabel("")

    sns.despine()

    save_fig(fig, "Fig_Dataset_country_bars_nstrains.pdf", PlotsDir)


def get_sorted_binary_phenotypes(target_phenos):

    """gets the sorted phenos"""

    # define sorting
    resistance_phenos =  ['CAS_resistance', 'MIF_resistance', 'ANI_resistance', '5FC_resistance', 'AMB_resistance', 'FLC_resistance', 'VRC_resistance', 'ITR_resistance', 'POS_resistance', 'IVZ_resistance']
    disinfectants_phenos = ['ETHANOL_tolerance', 'HYPLOCHLORITE_tolerance', 'H202_tolerance', 'CHLORHEXIDINE_tolerance', 'ZWITTERGENT 3-14_tolerance', 'SURFANIO_tolerance']
    etest_phenos = ['etest_colony_formation_FLC_MIC_8-32', 'FLC_MIC_higher_in_solid', 'FLC_MIC_higher_48h']
    lab_phenos = ['agar_invasiveness', 'pseudohyphae_formation', 'biofilm_formation', 'microfluidics_Adherent', 'microfluidics_Big_Clumps', 'microfluidics_Little_clump']

    age_phenos = ['age_infant', 'age_young', 'age_old']
    general_sources = ['source_HC', 'source_SANGRE_COMPLETA', 'source_LIQUIDO_ASCITICO', 'source_PIEL', 'source_EXUDADO_CUTANEO', 'source_ABSCESO', 'source_ORINA', 'source_ASPIRADO_BRONQUIAL', 'source_ASPIRADO_TRAQUEAL', 'source_ESPUTO',  'source_LAVADO_BRONCOALVEOLAR', 'source_MUESTRA_AMBIENTAL']

    grouped_sources = ['source_sterile', 'source_sterile_vs_environmental', 'source_blood', 'source_blood_vs_environmental', 'source_blood_vs_other_sterile', 'source_skin', 'source_respiratory']

    all_sorted_phenos = resistance_phenos + etest_phenos + disinfectants_phenos + lab_phenos + age_phenos + general_sources + grouped_sources

    # get the target phenos
    pheno_to_idx = dict(zip(all_sorted_phenos, range(len(all_sorted_phenos))))

    return sorted(target_phenos, key=(lambda p: pheno_to_idx[p]))

def get_contrast_color(hex_color):

    """takes color and gets the contrast"""

    # Convert hex color to RGB
    rgb = mcolors.hex2color(hex_color)
    
    # Calculate the perceived brightness using the relative luminance formula
    # Formula: 0.299*R + 0.587*G + 0.114*B
    brightness = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    
    # Return white for dark colors, black for light colors
    if brightness<0.5:
        return "white"

    else:
        return "black"

def plot_Phenotypes_presence_tree(df_metadata, PlotsDir, samples_discard_run1, tree, df_metadata_binary_phenotypes, df_transitions_info_GWAS):

    """Plot tree with phenotypes"""

    # get metadata filtered
    df_metadata = df_metadata[~(df_metadata.runID.isin(samples_discard_run1)) & (df_metadata.possibly_triploid==False) & (df_metadata.source_data=="sequenced_here")].copy()
    df_metadata_binary_phenotypes = df_metadata_binary_phenotypes.copy().drop(["sampleID"], axis=1).set_index("runID")
    #df_blank_pheno = pd.DataFrame({"blank": {p : -1 for p in df_metadata_binary_phenotypes.columns}}).transpose()
    #df_metadata_binary_phenotypes = pd.concat([df_metadata_binary_phenotypes, df_blank_pheno])

    # define the phenotypes
    df_pheno_info = pd.DataFrame(index=df_metadata_binary_phenotypes.columns)
    df_pheno_info["n_phenos1"] = df_metadata_binary_phenotypes.apply(lambda c: sum(c==1), axis=0)
    df_pheno_info["n_phenos0"] = df_metadata_binary_phenotypes.apply(lambda c: sum(c==0), axis=0)
    df_pheno_info["n_phenos-1"] = df_metadata_binary_phenotypes.apply(lambda c: sum(c==-1), axis=0)

    target_phenotypes = set(df_pheno_info[(df_pheno_info.n_phenos1>=3) & (df_pheno_info.n_phenos0>=3)].index)
    sorted_phenotypes = get_sorted_binary_phenotypes(target_phenotypes)

    # map each run to the clade
    runID_to_clade = df_metadata.set_index("runID").cladeID_based_on_Bergin2022

    # keep consistent samples
    tree.prune(sorted(set(df_metadata.runID)), preserve_branch_length=True)

    # get the clade
    cladeID_to_color = get_cladeID_to_color_manual()
    #cladeID_to_color[100] = "red"

    runID_to_clade_Bergin2022 = df_metadata.set_index("runID").clade_Bergin2022
    tree = get_tree_withcladeID(tree, runID_to_clade)

    # init style
    ts = TreeStyle()

    # get color of branches according to clades
    for n in tree.traverse():

        # define the bwidth    
        bwidth = 4

        nst = NodeStyle()
        nst["hz_line_width"] = bwidth
        nst["vt_line_width"] = bwidth
        nst["size"] = 0

        # define the color
        linecolor = cladeID_to_color[n.cladeID]
      
        nst["bgcolor"] = "white"
        nst["hz_line_color"] = linecolor
        nst["vt_line_color"] = linecolor
        
        n.set_style(nst)


    # add sepparator nodes
    for cladeID in cladeID_to_color:
        leafs_clade = [l.name for l in tree.get_leaves() if l.cladeID==cladeID]
        if len(leafs_clade)==0: continue
        clade_node = tree.get_common_ancestor(leafs_clade)

        len_sep = len(leafs_clade) * 7 + 15
        clade_node.add_face(RectFace(0, len_sep, fgcolor=None, bgcolor="white"), position="branch-right", column=0)



    # add things to leaves
    pheno_to_col = {}
    palette_phenos = {-1:"white", 0:"darkgray", 1:"black"}
    df_clades = pd.DataFrame({"leaf_name": tree.get_leaf_names(), "clade": [l.cladeID for l in tree.get_leaves()]})

    for Il, l in enumerate(tree.get_leaves()):

        # get info
        l_series = df_metadata_binary_phenotypes.loc[l.name]

        # init with a sepparator
        col = 0
        l.add_face(RectFace(3, 7, fgcolor="white", bgcolor="white"), position="aligned", column=col); col+=1

        # add the color of the clade
        clade_color = cladeID_to_color[l.cladeID]
        col_clade_color = col
        l.add_face(RectFace(9, 6, fgcolor=None, bgcolor=clade_color), position="aligned", column=col); col+=1

        # add sepparator
        l.add_face(RectFace(4, 4, fgcolor=None, bgcolor="white"), position="aligned", column=col); col+=1

        # add each phenotype
        for pheno in sorted_phenotypes:

            pheno_val = l_series.loc[pheno]
            pheno_color = mcolors.to_hex(palette_phenos[pheno_val])
            pheno_to_col[pheno] = col
            l.add_face(RectFace(18, 6, fgcolor=None, bgcolor=pheno_color), position="aligned", column=col); col+=1

            # sep
            l.add_face(RectFace(2, 4, fgcolor=None, bgcolor="white"), position="aligned", column=col); col+=1

    # define palettes
    max_nphenos = max(df_pheno_info[["n_phenos0", "n_phenos1"]].apply(max))
    nvalues_phenos = list(range(0, len(df_metadata)+1))
    nphenos_to_color = get_value_to_color(nvalues_phenos, palette="rocket_r", n=len(nvalues_phenos), type_color="hex")[0]

    # add headers and footers
    df_transitions_info_GWAS = df_transitions_info_GWAS[(df_transitions_info_GWAS.aligner=="bwa_mem") & (df_transitions_info_GWAS.sample_set.isin({"representative_samples"})) & (df_transitions_info_GWAS.phenotype.isin(target_phenotypes))].copy()

    for pheno, col in pheno_to_col.items():

        # add header
        face_lab = TextFace("  "+get_paper_dicotomic_pheno(pheno), fsize=10, bold=False); face_lab.rotation = 270
        face_lab.vt_align = 2
        face_lab.hz_align = 2
        ts.aligned_header.add_face(face_lab, column=col) 

        # add footer with # transitions
        for target_p in [1, 0]:
            n_strains = df_pheno_info.loc[pheno, "n_phenos%i"%target_p]

            text_color = get_contrast_color(nphenos_to_color[n_strains])
            if n_strains<10:
                text_str = str(n_strains)

            else:
                text_str = ""
            face_foot = RectFace(18, 24, fgcolor="black", bgcolor=nphenos_to_color[n_strains], label={"text":text_str, "color":text_color, "fontsize":12})
            #face_foot = RectFace(18, 18, fgcolor=mcolors.to_hex(palette_phenos[target_p]), bgcolor=nphenos_to_color[n_strains])
            ts.aligned_foot.add_face(face_foot, column=col)

            # spacer
            #ts.aligned_foot.add_face(RectFace(18, 5, fgcolor=None, bgcolor="white"), column=col)

        # add the number of transitions
        n_transitions = df_transitions_info_GWAS[(df_transitions_info_GWAS.phenotype==pheno)].iloc[0].n_transitions

        text_f = TextFace(str(n_transitions), bold=False, fsize=12)
        text_f.hz_align = 1
        ts.aligned_foot.add_face(text_f, column=col)

    # add to the last col the label
    for label in ["  isolates pheno.", "  isolates no pheno.", "  transitions"]:
        text_face = TextFace(label, fsize=12, bold=False)
        text_face.hz_align = 0
        ts.aligned_foot.add_face(text_face, column=col+1)
        #ts.aligned_foot.add_face(RectFace(18, 5, fgcolor=None, bgcolor="white"), column=col_clade_color)

    # add legend
    for nstrains_float in range(0, max(nvalues_phenos)+1, 20):

        print(nstrains_float)

        if nstrains_float in {0, 60, 120, 180}:
            text_str = " "+str(int(nstrains_float))

        else:
            text_str = ""

        #print(text_str)

        #text_str = "hi"
        
        ts.legend.add_face(RectFace(30, 30, fgcolor=None, bgcolor=nphenos_to_color[int(nstrains_float)]), column=0)


        ts.legend.add_face(TextFace(text_str, fsize=12, bold=False), column=1)



    # print
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_leaf_name = False
    ts.mode = "r"
    #ts.draw_guiding_lines = 0 # 0=solid, 1=dashed, 2=dotted.
    ts.show_scale = False     # Hide scale bar if not needed

    tree.render(PlotsDir+"/phenotypes_overview.pdf", tree_style=ts)

    tree.show(tree_style=ts)


def get_one_hot_encoded_val(x, runs_1, runs_0, runs_minus1):

    """gets vals"""
    
    if x in runs_minus1: return -1
    elif x in runs_1: return 1
    elif x in runs_0: return 0
    else: raise ValueError("invalid")

def get_clinical_paper_df_metadata_binary_phenotypes_w_extra_fields(df_phenos, df_metadata):

    """get extra fields"""

    # keep
    df_phenos = df_phenos.reset_index().copy()
    df_metadata = df_metadata.copy()

    # add fields
    df_metadata["hospital_reduced"] = df_metadata.hospital_short.map(map_hospital_short_to_very_short())
    check_no_nans_series(df_metadata["hospital_reduced"])

    # add one-hot encoded things
    for f in ["hospital_reduced", "ERG11_mutation"]:

        for val in sorted(set(df_metadata[f])):

            # define vals
            if f in {"hospital_reduced"}:
                runs_1 = set(df_metadata[df_metadata[f]==val].index)
                runs_0 = set(df_metadata[df_metadata[f]!=val].index)
                runs_minus1 = set()
                suffix_f = ""

            elif f=="ERG11_mutation":
                runs_1 = set(df_metadata[df_metadata[f]==val].index)
                runs_0 = set(df_metadata[df_metadata[f]!=val].index)
                runs_minus1 = set(df_metadata[df_metadata[f]=="unknown"].index)
                suffix_f = "ERG11 - "

            else: raise ValueError("invalid")

            df_phenos[suffix_f+val] = df_phenos.runID.apply(get_one_hot_encoded_val, runs_1=runs_1, runs_0=runs_0, runs_minus1=runs_minus1)

    return df_phenos.set_index("runID")

def generate_data_clinical_paper(df_metadata_trees, PlotsDir, samples_discard_run1, strains_tree, df_metadata_binary_phenotypes):

    """gets table and tree"""

    # get data
    df_metadata_clinical = df_metadata_trees[~(df_metadata_trees.runID.isin(samples_discard_run1)) & (df_metadata_trees.possibly_triploid==False) & (df_metadata_trees.source_data=="sequenced_here")].copy()

    # get table
    df_write = df_metadata_clinical[["sampleID", "cladeID_based_on_Bergin2022", "hospital_short", "ERG11_mutation", "genotype"]].reset_index().rename(columns={"cladeID_based_on_Bergin2022":"clade", "index":"runID"}).copy()

    print(df_write)

    for drug_r in ['FLC_resistance', 'VRC_resistance', 'ITR_resistance', 'POS_resistance', 'IVZ_resistance']:
        drug = drug_r.split("_")[0]


        int_dict = {0:"S", 1:"R", -1:"I/?"}
        df_write[drug_r] = df_write.runID.apply(lambda run: "%s (%.3f)"%(int_dict[df_metadata_binary_phenotypes.loc[run, drug_r]], df_metadata_clinical.loc[run, drug+"_MIC"]))

    # write
    df_write.drop(["runID"], axis=1).sort_values(by=["clade", "hospital_short", "sampleID"]).to_excel(PlotsDir+"/clinical_paper_WGS_strains.xlsx")

    # tree
    strains_tree = cp.deepcopy(strains_tree)
    for l in strains_tree.get_leaves(): l.name = l.name.split("_")[1]
    strains_tree.write(outfile=PlotsDir+"/clinical_paper_WGS_tree.nw", format=2)

def plot_validation_azole_phenotypes_tree(df_StrainInfo, DataDir, PlotsDir):
    
    """Plots the phenotypes of the validation along the tree"""
    

    # check that for extra phenotypes all flz also include vrc
    df_extra = df_StrainInfo[(df_StrainInfo.source_data!="sequenced_here") & (df_StrainInfo.has_smiley==False) &
                             ((df_StrainInfo["binaryPheno voriconazole R"]!=-1) | (df_StrainInfo["binaryPheno fluconazole R"]!=-1))].copy()
                             # ((df_StrainInfo.FLC_MIC!=-1) | (df_StrainInfo.VRC_MIC!=-1))].copy()

    #print(df_extra[["sampleID", "FLC_MIC", "VRC_MIC", "binaryPheno voriconazole R", "binaryPheno fluconazole R"]], len(df_extra))
    if any(df_extra.FLC_MIC==-1): raise ValueError("all FLC should be measured")
            
    # load tree of the FLC resistance of the validation with the top model (using hisat2)
    strains_tree = Tree(DataDir + "/validation_data/genotype_phenotype_relationships/pheno-FLC_resistance_aligner-hisat2_no_spliced/tree_subset_strains/rooted_tree.nw")
    print("Tree has %i strains."%len(strains_tree.get_leaf_names()))
    
    # keep from extra only the clear flc resistance ones (removing one strains with intermediate FLC resistance and VRC measured)
    df_extra = df_extra[df_extra["binaryPheno fluconazole R"]!=-1].copy()
    print(df_extra[["sampleID", "FLC_MIC", "VRC_MIC", "binaryPheno voriconazole R", "binaryPheno fluconazole R"]], len(df_extra))
        
    # check that all are there
    missing_strains = set(df_extra.sampleID).difference(set(strains_tree.get_leaf_names()))    
    if len(missing_strains)>0: raise ValueError("missing strains: %s"%missing_strains)
    
    # change leaf names
    for l in strains_tree.get_leaves():
        if "OSCAR" in l.name: l.name = l.name.split("_")[-1]
        
    # keep some fields
    df_StrainInfo = df_StrainInfo[["binaryPheno voriconazole R", "binaryPheno fluconazole R", "sampleID", "source_data", "validation_study", "clade_basedOn_Bergin2022"]].copy().set_index("sampleID").rename(columns={"binaryPheno voriconazole R":"voriconazole R", "binaryPheno fluconazole R":"fluconazole R"})
    df_StrainInfo = df_StrainInfo.loc[strains_tree.get_leaf_names()].copy()
    
    print(sum(df_StrainInfo.source_data=="sequenced_here"), "strains sequenced here")
    
    # init phenotypes    
    bool_to_number = {False:0, True:1}
    df_StrainInfo["this work (Spain)"] = (df_StrainInfo.source_data=="sequenced_here").map(bool_to_number)
    df_StrainInfo["blank"] = -1
    sorted_phenotypes = ["fluconazole R", "voriconazole R", "this work (Spain)"]

    # add the study of origin
    study_to_str = {"Branco_2022":"Branco 2022 (Portugal)", "Brassington_2024":"Brassington 2024 (Germany)", "McTaggart_2024":"McTaggart 2024 (Canada)"}
    for validation_study in sorted(set(df_StrainInfo[df_StrainInfo.source_data!="sequenced_here"].validation_study)):
        df_StrainInfo[study_to_str[validation_study]] = (df_StrainInfo.validation_study==validation_study).map(bool_to_number)
        sorted_phenotypes.append(study_to_str[validation_study])
        
            
    # get the clade info
    cladeID_to_color = get_cladeID_to_color_manual()
    cladeID_to_color = {str(k):v for k,v in cladeID_to_color.items()}
    cladeID_to_color["-"] = cladeID_to_color["-1"]
    for cINT in range(1,6): cladeID_to_color["%i.0"%cINT] = cladeID_to_color[str(cINT)]
    strains_tree = get_tree_withcladeID(strains_tree, df_StrainInfo.clade_basedOn_Bergin2022)
    for n in strains_tree.traverse(): n.cladeID = str(n.cladeID)

    # init style
    ts = TreeStyle()

    # get color of branches according to clades
    for n in strains_tree.traverse():

        # define the bwidth    
        bwidth = 4

        nst = NodeStyle()
        nst["hz_line_width"] = bwidth
        nst["vt_line_width"] = bwidth
        nst["size"] = 0

        # define the color
        linecolor = cladeID_to_color[(n.cladeID)]
      
        nst["bgcolor"] = "white"
        nst["hz_line_color"] = linecolor
        nst["vt_line_color"] = linecolor
        
        n.set_style(nst)
    
    # add sepparator nodes
    for cladeID in cladeID_to_color:
        leafs_clade = [l.name for l in strains_tree.get_leaves() if (l.cladeID)==cladeID]
        if len(leafs_clade)==0: continue
        clade_node = strains_tree.get_common_ancestor(leafs_clade)

        len_sep = len(leafs_clade) * 7 + 15
        clade_node.add_face(RectFace(0, len_sep, fgcolor=None, bgcolor="white"), position="branch-right", column=0)
    
        
    # add things to leaves
    pheno_to_col = {}
    palette_phenos = {-1:"white", 0:"darkgray", 1:"black"}
    for Il, l in enumerate(strains_tree.get_leaves()):

        # get info
        l_series = df_StrainInfo.loc[l.name]

        # init with a sepparator
        col = 0
        l.add_face(RectFace(3, 7, fgcolor="white", bgcolor="white"), position="aligned", column=col); col+=1

        # add the color of the clade
        clade_color = cladeID_to_color[l.cladeID]
        col_clade_color = col
        l.add_face(RectFace(9, 6, fgcolor=None, bgcolor=clade_color), position="aligned", column=col); col+=1

        # add sepparator
        l.add_face(RectFace(4, 4, fgcolor=None, bgcolor="white"), position="aligned", column=col); col+=1

        # add each phenotype
        for pheno in sorted_phenotypes:

            pheno_val = l_series.loc[pheno]
            pheno_color = mcolors.to_hex(palette_phenos[pheno_val])
            pheno_to_col[pheno] = col
            l.add_face(RectFace(18, 6, fgcolor=None, bgcolor=pheno_color), position="aligned", column=col); col+=1

            # sep
            if pheno in {"voriconazole R"}: sep_width = 12
            else: sep_width = 2

            l.add_face(RectFace(sep_width, 4, fgcolor=None, bgcolor="white"), position="aligned", column=col); col+=1
    
    # add headers and footers
    for pheno, col in pheno_to_col.items():
        
        if pheno=="blank": continue

        # add header
        face_lab = TextFace("  "+pheno+" ", fsize=10, bold=False); face_lab.rotation = 270
        face_lab.vt_align = 0
        face_lab.hz_align = 0
        ts.aligned_foot.add_face(face_lab, column=col) 

    # print
    print("saving...")    
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_leaf_name = False
    ts.mode = "r"
    #ts.draw_guiding_lines = 0 # 0=solid, 1=dashed, 2=dotted.
    ts.show_scale = False     # Hide scale bar if not needed
    strains_tree.render(PlotsDir+"/tree_MLvalidation_azoles.pdf", tree_style=ts)
    #tree.show(tree_style=ts)
    
    
    

def plot_tree_clinical_paper(df_metadata, PlotsDir, samples_discard_run1, tree, df_metadata_binary_phenotypes):

    """Plot tree with phenotypes"""

    # get metadata filtered
    df_metadata = df_metadata[~(df_metadata.runID.isin(samples_discard_run1)) & (df_metadata.possibly_triploid==False) & (df_metadata.source_data=="sequenced_here")].copy()
    df_metadata_binary_phenotypes = df_metadata_binary_phenotypes.copy().drop(["sampleID"], axis=1).set_index("runID")
    #df_blank_pheno = pd.DataFrame({"blank": {p : -1 for p in df_metadata_binary_phenotypes.columns}}).transpose()
    #df_metadata_binary_phenotypes = pd.concat([df_metadata_binary_phenotypes, df_blank_pheno])

    # add extra binary phenotypes
    df_metadata_binary_phenotypes = get_clinical_paper_df_metadata_binary_phenotypes_w_extra_fields(df_metadata_binary_phenotypes, df_metadata)

    # define the phenotypes
    df_pheno_info = pd.DataFrame(index=df_metadata_binary_phenotypes.columns)
    df_pheno_info["n_phenos1"] = df_metadata_binary_phenotypes.apply(lambda c: sum(c==1), axis=0)
    df_pheno_info["n_phenos0"] = df_metadata_binary_phenotypes.apply(lambda c: sum(c==0), axis=0)
    df_pheno_info["n_phenos-1"] = df_metadata_binary_phenotypes.apply(lambda c: sum(c==-1), axis=0)

    target_phenotypes = set(df_pheno_info.index)
    #sorted_phenotypes = get_sorted_binary_phenotypes(target_phenotypes)
    sorted_phenotypes = ['FLC_resistance', 'VRC_resistance', 'ITR_resistance', 'POS_resistance', 'IVZ_resistance'] + ['ERG11 - Y132F', 'ERG11 - Y132F (het)', 'ERG11 - R398I;Y132F', 'ERG11 - G458S', 'ERG11 - R398I'] + ['Barcelona - Bellvitge', 'Barcelona - J. Deu', 'Barcelona - V. Hebron', 'Lab. Catalunya', 'Madrid - 12 Oct', 'Madrid - C. Defensa', 'Madrid - Mostoles', 'Madrid - P. Hierro', 'Madrid - S. Carlos', 'Sevilla - S. Aljarafe', 'Sevilla - V. Macarena', 'Sevilla - V. Rocio', 'Valladolid', 'Burgos', 'Santander'] # 'ERG11 - REF',


    # map each run to the clade
    runID_to_clade = df_metadata.set_index("runID").cladeID_based_on_Bergin2022

    # keep consistent samples
    tree.prune(sorted(set(df_metadata.runID)), preserve_branch_length=True)

    # get the clade
    cladeID_to_color = get_cladeID_to_color_manual()
    #cladeID_to_color[100] = "red"

    runID_to_clade_Bergin2022 = df_metadata.set_index("runID").clade_Bergin2022
    tree = get_tree_withcladeID(tree, runID_to_clade)

    # init style
    ts = TreeStyle()

    # get color of branches according to clades
    for n in tree.traverse():

        # define the bwidth    
        bwidth = 4

        nst = NodeStyle()
        nst["hz_line_width"] = bwidth
        nst["vt_line_width"] = bwidth
        nst["size"] = 0

        # define the color
        linecolor = cladeID_to_color[n.cladeID]
      
        nst["bgcolor"] = "white"
        nst["hz_line_color"] = linecolor
        nst["vt_line_color"] = linecolor
        
        n.set_style(nst)


    # add sepparator nodes
    for cladeID in cladeID_to_color:
        leafs_clade = [l.name for l in tree.get_leaves() if l.cladeID==cladeID]
        if len(leafs_clade)==0: continue
        clade_node = tree.get_common_ancestor(leafs_clade)

        len_sep = len(leafs_clade) * 7 + 15
        clade_node.add_face(RectFace(0, len_sep, fgcolor=None, bgcolor="white"), position="branch-right", column=0)



    # add things to leaves
    pheno_to_col = {}
    palette_phenos = {-1:"white", 0:"darkgray", 1:"black"}
    df_clades = pd.DataFrame({"leaf_name": tree.get_leaf_names(), "clade": [l.cladeID for l in tree.get_leaves()]})

    for Il, l in enumerate(tree.get_leaves()):

        # get info
        l_series = df_metadata_binary_phenotypes.loc[l.name]

        # init with a sepparator
        col = 0
        l.add_face(RectFace(3, 7, fgcolor="white", bgcolor="white"), position="aligned", column=col); col+=1

        # add the color of the clade
        clade_color = cladeID_to_color[l.cladeID]
        col_clade_color = col
        l.add_face(RectFace(9, 6, fgcolor=None, bgcolor=clade_color), position="aligned", column=col); col+=1

        # add sepparator
        l.add_face(RectFace(4, 4, fgcolor=None, bgcolor="white"), position="aligned", column=col); col+=1

        # add each phenotype
        for pheno in sorted_phenotypes:

            pheno_val = l_series.loc[pheno]
            pheno_color = mcolors.to_hex(palette_phenos[pheno_val])
            pheno_to_col[pheno] = col
            l.add_face(RectFace(18, 6, fgcolor=None, bgcolor=pheno_color), position="aligned", column=col); col+=1

            # sep
            if pheno in {"IVZ_resistance", "ERG11 - R398I"}: sep_width = 12
            else: sep_width = 2

            l.add_face(RectFace(sep_width, 4, fgcolor=None, bgcolor="white"), position="aligned", column=col); col+=1

    # define palettes
    max_nphenos = max(df_pheno_info[["n_phenos0", "n_phenos1"]].apply(max))
    nvalues_phenos = list(range(0, len(df_metadata)+1))
    nphenos_to_color = get_value_to_color(nvalues_phenos, palette="rocket_r", n=len(nvalues_phenos), type_color="hex")[0]

    # add headers and footers
    for pheno, col in pheno_to_col.items():

        # add header
        face_lab = TextFace("  "+get_paper_dicotomic_pheno(pheno), fsize=10, bold=False); face_lab.rotation = 270
        face_lab.vt_align = 2
        face_lab.hz_align = 2
        ts.aligned_header.add_face(face_lab, column=col) 

        # add footer with # transitions
        for target_p in [1, 0]:
            n_strains = df_pheno_info.loc[pheno, "n_phenos%i"%target_p]

            text_color = get_contrast_color(nphenos_to_color[n_strains])
            if n_strains<10:
                text_str = str(n_strains)

            else:
                text_str = ""
            face_foot = RectFace(18, 24, fgcolor="black", bgcolor=nphenos_to_color[n_strains], label={"text":text_str, "color":text_color, "fontsize":12})
            #face_foot = RectFace(18, 18, fgcolor=mcolors.to_hex(palette_phenos[target_p]), bgcolor=nphenos_to_color[n_strains])
            ts.aligned_foot.add_face(face_foot, column=col)

            # spacer
            #ts.aligned_foot.add_face(RectFace(18, 5, fgcolor=None, bgcolor="white"), column=col)

    # add to the last col the label
    for label in ["  # isolates yes", "  # isolates no"]:
        text_face = TextFace(label, fsize=12, bold=False)
        text_face.hz_align = 0
        ts.aligned_foot.add_face(text_face, column=col+1)
        #ts.aligned_foot.add_face(RectFace(18, 5, fgcolor=None, bgcolor="white"), column=col_clade_color)

    # add legend
    for nstrains_float in range(0, max(nvalues_phenos)+1, 20):

        print(nstrains_float)

        if nstrains_float in {0, 60, 120, 180}:
            text_str = " "+str(int(nstrains_float))

        else:
            text_str = ""

        #print(text_str)

        #text_str = "hi"
        
        ts.legend.add_face(RectFace(30, 30, fgcolor=None, bgcolor=nphenos_to_color[int(nstrains_float)]), column=0)


        ts.legend.add_face(TextFace(text_str, fsize=12, bold=False), column=1)



    # print
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_leaf_name = False
    ts.mode = "r"
    #ts.draw_guiding_lines = 0 # 0=solid, 1=dashed, 2=dotted.
    ts.show_scale = False     # Hide scale bar if not needed

    tree.render(PlotsDir+"/tree_clinical_paper.pdf", tree_style=ts)

    #tree.show(tree_style=ts)



def plot_strains_tree_with_clades(filename_tree, tree, df_metadata, samples_discard_run1):

    """Plots a tree with clades"""

    # keep
    tree = cp.deepcopy(tree)
    df_metadata = df_metadata[~(df_metadata.runID.isin(samples_discard_run1)) & (df_metadata.possibly_triploid==False)].copy()
    df_metadata_idx = df_metadata.copy().set_index("runID")

    # map each run to the clade
    runID_to_clade = df_metadata.set_index("runID").cladeID_based_on_Bergin2022

    # keep consistent samples
    tree.prune(sorted(set(df_metadata.runID)), preserve_branch_length=True)

    # get the clade
    #cladeID_to_color = get_cladeID_to_color(sorted(set(df_metadata.cladeID_based_on_Bergin2022)))
    cladeID_to_color = get_cladeID_to_color_manual()

    runID_to_clade_Bergin2022 = df_metadata.set_index("runID").clade_Bergin2022
    tree = get_tree_withcladeID(tree, runID_to_clade)

    # init style
    ts = TreeStyle()

    # get color of branches according to clades
    for n in tree.traverse():

        # define the bwidth    
        bwidth = 4

        nst = NodeStyle()
        nst["hz_line_width"] = bwidth
        nst["vt_line_width"] = bwidth
        nst["size"] = 0

        # define the color
        linecolor = cladeID_to_color[n.cladeID]
      
        nst["bgcolor"] = "white"
        nst["hz_line_color"] = linecolor
        nst["vt_line_color"] = linecolor
        
        n.set_style(nst)

    # add things to leaves
    source_data_to_color = {"sequenced_here":"black", "SRA":"gray"}
    for l in tree.get_leaves():

        # get info
        l_series = df_metadata_idx.loc[l.name]

        # add type of data
        l.add_face(RectFace(3, 5, fgcolor="white", bgcolor="white"),  column=0)
        if l_series.source_data=="sequenced_here":
            l.add_face(RectFace(4, 4, fgcolor=source_data_to_color[l_series.source_data], bgcolor=source_data_to_color[l_series.source_data]),  column=1)



    # add legend
    for cladeID in sorted(cladeID_to_color):
        if cladeID==-1: continue

        ts.legend.add_face(RectFace(20, 20, fgcolor="black", bgcolor=cladeID_to_color[cladeID]),  column=0)
        ts.legend.add_face(TextFace(" clade %s"%cladeID, bold=False, fsize=8), column=1)

    ts.legend.add_face(RectFace(20, 20, fgcolor="white", bgcolor="white"),  column=0)

    # print
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_leaf_name = False
    ts.mode = "r"
    ts.arc_start = -180 # 0 degrees = 3 o'clock
    ts.arc_span = 360  
    ts.root_opening_factor = 1
    ts.draw_guiding_lines = False # 0=solid, 1=dashed, 2=dotted.
    ts.show_scale = False     # Hide scale bar if not needed

    tree.show(tree_style=ts)
    print(filename_tree)
    tree.render(filename_tree, tree_style=ts)

def visualize_tree_clinical_info(filename_tree, tree, df_metadata, orientation_tree, show_t=False):

    """Visualizes tree with clinical info"""

    tree = cp.deepcopy(tree)
    df_metadata = cp.deepcopy(df_metadata)

    # keep
    df_metadata = df_metadata[df_metadata.runID.isin(set(tree.get_leaf_names()))].copy()
    print("There are %i/%i runs w no clade"%(sum(df_metadata.cladeID_based_on_Bergin2022==-1), len(df_metadata)))



    # get variables
    #cladeID_to_color = get_cladeID_to_color(sorted(set(df_metadata.cladeID_BL_and_support)))
    cladeID_to_color = get_cladeID_to_color(sorted(set(df_metadata.cladeID_based_on_Bergin2022)))
    runID_to_clade = df_metadata.set_index("runID").cladeID_based_on_Bergin2022
    runID_to_clade_Bergin2022 = df_metadata.set_index("runID").clade_Bergin2022
    tree = get_tree_withcladeID(tree, runID_to_clade)
    resistance_to_color = {"R":"red", "S":"blue", "I":"gray", "unknown":"white"}

    erg11_mut_to_color = {'G458S': 'orange', 'R398I': 'blue', 'R398I;Y132F': 'purple', 'REF': 'gray', 'Y132F': 'red', 'Y132F (het)': 'magenta', 'unknown': 'white'}

    # define sorted resistances:
    sorted_resistances = ['FLC_resistance', 'ITR_resistance', 'POS_resistance', 'VRC_resistance']

    # init style
    ts = TreeStyle()

    # get color of branches according to clades
    for n in tree.traverse():

        # define the bwidth    
        bwidth = 4

        nst = NodeStyle()
        nst["hz_line_width"] = bwidth
        nst["vt_line_width"] = bwidth
        nst["size"] = 0

        # define the color
        linecolor = cladeID_to_color[n.cladeID]
      
        nst["bgcolor"] = "white"
        nst["hz_line_color"] = linecolor
        nst["vt_line_color"] = linecolor
        
        n.set_style(nst)

    #country_to_color
    #source_data_to_color
    df_metadata = df_metadata.set_index("runID")

    # add things
    for l in tree.get_leaves():
        metadata_l = df_metadata.loc[l.name]

        # init the colID and rect width
        colID = 0
        rect_width = 20

        # add label with Bergin2020 clades
        #if runID_to_clade_Bergin2022[l.name]==-1: label_clade_bergin = ""
        #else: label_clade_bergin = str(runID_to_clade_Bergin2022[l.name])
        #l.add_face( TextFace(label_clade_bergin, fsize=10, bold=True, fgcolor=cladeID_to_color[l.cladeID]), position="aligned", column=colID); colID+=1

        # add clade
        l.add_face(RectFace(rect_width*2, rect_width, fgcolor="black", bgcolor=cladeID_to_color[l.cladeID]), position="aligned", column=colID); colID+=1


        # add source data
        l.add_face(RectFace(rect_width*2, rect_width, fgcolor="black", bgcolor=source_data_to_color[metadata_l.source_data]), position="aligned", column=colID); colID+=1

        # add country
        l.add_face(RectFace(rect_width*2, rect_width, fgcolor="black", bgcolor=country_to_color[metadata_l.country_collection]), position="aligned", column=colID); colID+=1

        # add hospital info
        if metadata_l.source_data=="sequenced_here" and l.name!="OSCAR01_L21-282": # there is one strain, L21-282, with no info
            color = hospital_short_to_color[metadata_l.hospital_short]
            fgcolor = "black"
            color_text = hospital_short_to_color[metadata_l.hospital_short]
            label = " %s. %s. %s"%(l.name, metadata_l.hospital_short, metadata_l.year_collection)

        else:
            color = "white"
            fgcolor = "white"
            color_text = "gray"
            if metadata_l.year_collection==-1: label_year = "unknown year"
            else: label_year = str(int(metadata_l.year_collection))

            label = " %s. %s (%s). %s"%(l.name, metadata_l.city_collection, metadata_l.country_collection, label_year)

        # hospital
        l.add_face(RectFace(rect_width*2, rect_width, fgcolor=fgcolor, bgcolor=color), position="aligned", column=colID); colID+=1

        # face white
        l.add_face(RectFace(rect_width*0.75, rect_width, fgcolor="white", bgcolor="white"), position="aligned", column=colID); colID+=1

        # add metadata of resistance
        for f in sorted_resistances: 

            if metadata_l[f]!="unknown": fgcolor_r = "black"
            else: fgcolor_r = "white"

            l.add_face(RectFace(rect_width*2, rect_width, fgcolor=fgcolor_r, bgcolor=resistance_to_color[metadata_l[f]]), position="aligned", column=colID); colID+=1

        # add the erg11 mutation
        #if metadata_l.ERG11_mutation!="unknown": fgcolor_e = "black"
        #else: fgcolor_e = "white"

        # face white
        l.add_face(RectFace(rect_width*0.75, rect_width, fgcolor="white", bgcolor="white"), position="aligned", column=colID); colID+=1

        #l.add_face(RectFace(rect_width*2, rect_width, fgcolor=fgcolor, bgcolor=erg11_mut_to_color[metadata_l.ERG11_mutation]), position="aligned", column=colID); colID+=1
        l.add_face(CircleFace(rect_width*0.8, erg11_mut_to_color[metadata_l.ERG11_mutation], style='sphere'), position="aligned", column=colID); colID+=1

        # add hospital / city
        if orientation_tree=="r": l.add_face(TextFace(label, fsize=10, bold=True, fgcolor=color_text), position="aligned", column=colID); colID+=1

        # add label with 
        #if l.cladeID==-1: label_clade = ""
        #else: label_clade = str(l.cladeID)
        #l.add_face(RectFace(rect_width, rect_width, fgcolor="black", bgcolor="white", label={"text":label_clade, "color":cladeID_to_color[l.cladeID], "size":14}), position="aligned", column=colID); colID+=1

    # add legend and footer
    for colI, (label, color_dict, type_face) in enumerate([("clade (Bergin 2022)", cladeID_to_color, "rect"), ("source", source_data_to_color, "rect"), ("country", country_to_color, "rect"), ("hospital", hospital_short_to_color, "rect")] + [("", {}, "rect")] + [(f, resistance_to_color, "rect") for f in sorted_resistances] + [("", {}, "rect"), ("ERG11", erg11_mut_to_color, "circle")]):

        # add label to plot
        face_lab = TextFace(" "+label, fsize=20, bold=True); face_lab.rotation = 90
        face_lab.vt_align = 0
        ts.aligned_foot.add_face(face_lab, column=colI) 

        # add the legend
        ts.legend.add_face(TextFace(label+"\n", bold=True, fsize=12), column=0)
        ts.legend.add_face(TextFace("\n", bold=True, fsize=12), column=1)

        for color_lab, color in color_dict.items():

            if str(color_lab)=="-1": color_lab = "unknown"

            if type_face=="rect": ts.legend.add_face(RectFace(30, 30, fgcolor='black', bgcolor=color),  column=0)
            elif type_face=="circle": ts.legend.add_face(CircleFace(20, color, style='sphere'),  column=0)

            ts.legend.add_face(TextFace(color_lab, bold=True, fsize=12), column=1)


        # add circle
        for c in [0, 1]: ts.legend.add_face(TextFace("\n", bold=True, fsize=12), column=c)




    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_leaf_name = False
    ts.mode = orientation_tree
    ts.arc_start = -180 # 0 degrees = 3 o'clock
    ts.arc_span = 360  
    ts.root_opening_factor = 1
    ts.draw_guiding_lines = 2 # 0=solid, 1=dashed, 2=dotted.




    if show_t: tree.show(tree_style=ts)
    tree.render(filename_tree, tree_style=ts)

def plot_tree_compared_with_Gilead(df_trees, df_metadata, tree_Gilead, aligner, outgroup, mode_tree, PlotsDir, min_support):

    """Compares the trees with gilead's"""

    # get the tree of gilead corrected
    samples_Gilead = set(tree_Gilead.get_leaf_names()).difference({'L19-086', 'L20-008', 'L21-164', 'L21-108', 'L21-129'}) # discard uncorrected samples
    tree_Gilead = cp.deepcopy(tree_Gilead)
    tree_Gilead.prune(samples_Gilead)

    # get tree of gilead
    tree = cp.deepcopy(df_trees[(df_trees.aligner==aligner) & (df_trees.outgroup==outgroup) & (df_trees.mode_tree==mode_tree)].iloc[0].tree_corrected)
    missing_samples = samples_Gilead.difference(df_metadata.sampleID)
    if len(missing_samples)>0: raise ValueError("samples from Gilead not in metadata: %s"%missing_samples)
    runIDs_gilead = set(df_metadata[df_metadata.sampleID.isin(samples_Gilead)].runID)
    tree.prune(runIDs_gilead)

    # set outgroup as in email's plot
    outgroup_n = tree_Gilead.get_common_ancestor("L20-080", "L21-031")
    tree_Gilead.set_outgroup(outgroup_n)

    tree_Gilead =  get_correct_tree(tree_Gilead, min_support=min_support)
    
    # change the name
    for l in tree.get_leaves(): l.name = l.name.split("_")[1]

    # print differences
    dist1, dist2 = get_pct_different_nodes_two_trees(tree, tree_Gilead, min_nleafs=3)
    print("%s-%s-%s. Difference gilead: %spct and %spct."%(outgroup, aligner, mode_tree, dist1, dist2))

    # change nodes depending on if they are on tree_Gilead
    all_nodes_Gilead = {tuple(sorted(n.get_leaf_names())) for n in tree_Gilead.traverse() if n.is_leaf() is False}

    for n in tree.traverse():

        # define the bwidth    
        bwidth = 4

        nst = NodeStyle()
        nst["hz_line_width"] = bwidth
        nst["vt_line_width"] = bwidth
        nst["size"] = 0

        # define the color
        if n.is_leaf(): linecolor = "black"
        elif tuple(sorted(n.get_leaf_names())) in all_nodes_Gilead: linecolor = "black"
        else: linecolor = "red"
      
        nst["bgcolor"] = "white"
        nst["hz_line_color"] = linecolor
        nst["vt_line_color"] = linecolor
        
        n.set_style(nst)

    ts = TreeStyle()
    ts.mode = "r"
    #tree.show(tree_style=ts)
    tree.render("%s/gilead_comparison_tree.pdf"%(PlotsDir), tree_style=ts)

    # get the gilead tree pruned
    ts = TreeStyle()
    for n in tree_Gilead.traverse():

        # define the bwidth    
        bwidth = 4

        nst = NodeStyle()
        nst["hz_line_width"] = bwidth
        nst["vt_line_width"] = bwidth
        nst["size"] = 0
      
        nst["bgcolor"] = "white"
        nst["hz_line_color"] = "black"
        nst["vt_line_color"] = "black"
        
        n.set_style(nst)


    tree_Gilead.render("%s/gilead_pruned_tree.pdf"%(PlotsDir), tree_style=ts)

def get_fullpath(x):

    """Takes a path and substitutes it bu the full path"""

    if x.startswith("/"): return x
    elif x.startswith("."): path = "/".join(x.split("/")[1:])
    else: path = x

    return os.getcwd() + "/" + path

def soft_link_files(origin, target):

    """This function takes an origin file and makes it accessible through a link (target)"""

    if file_is_empty(target):

        # rename as full paths
        origin = get_fullpath(origin)
        target = get_fullpath(target)

        # check that the origin exists
        if file_is_empty(origin): raise ValueError("The origin %s should exist"%origin)

        # remove previous lisqnk
        try: run_cmd("rm %s > /dev/null 2>&1"%target)
        except: pass

        soft_linking_std = "%s.softlinking.std"%(target)
        print("softlinking. The std is in %s"%soft_linking_std)
        run_cmd("ln -s %s %s > %s 2>&1"%(origin, target, soft_linking_std))
        remove_file(soft_linking_std)

    # check that it worked
    if file_is_empty(target): raise ValueError("The target %s should exist"%target)


def get_df_cov_per_gene_relative_gDNA_one_sample(Isample, df_cov_all):

    """Gets the cov per gene for one sample"""

    print("getting rel cov sample %i"%Isample)

    # get median cov
    median_cov_all = np.median(df_cov_all[df_cov_all.chromosome!=mtDNA_chrom].median_reads_per_gene)
    df_cov_all["relative_coverage_vs_median_gDNA"] = df_cov_all.median_reads_per_gene / median_cov_all

    return df_cov_all[['chromosome', 'start', 'end', 'median_reads_per_gene', 'percentcovered_1', 'ID', 'relative_coverage', 'sampleID', 'relative_coverage_vs_median_gDNA']]

def get_df_cov_per_gene_relative_gDNA(ProcessedDataDir):

    """Gets a coverage per gene relative to the gDNA"""

    df_cov_per_gene_file = "%s/df_cov_per_gene_relative_gDNA.py"%ProcessedDataDir
    if file_is_empty(df_cov_per_gene_file):

        # load df with all
        print("loading cov df")
        df_cov_per_gene_all = load_object("%s/coverage_per_gene_all_samples.py"%ProcessedDataDir)

        # get inputs
        print("gettin rel cov")
        inputs_fn = [(I, df) for I, (s,df) in enumerate(df_cov_per_gene_all.groupby("sampleID"))]

        with multiproc.Pool(multiproc.cpu_count()) as pool:
            df_cov_per_gene = pd.concat(pool.starmap(get_df_cov_per_gene_relative_gDNA_one_sample, inputs_fn, chunksize=1)) 
            pool.close()
            pool.terminate()

        # change fields
        print("adding fields")
        df_cov_per_gene["runID"] = df_cov_per_gene.sampleID.apply(lambda x: x.split("---")[0])
        df_cov_per_gene["aligner"] = df_cov_per_gene.sampleID.apply(lambda x: x.split("---")[1])

        # add pos
        df_cov_per_gene["pos"] = (df_cov_per_gene.start + ((df_cov_per_gene.end - df_cov_per_gene.start)/2)).apply(int)

        #save 
        print("saving")
        save_object(df_cov_per_gene, df_cov_per_gene_file)

    return df_cov_per_gene_file

def get_chrom_to_color(genome):

    """Gets chrom to color"""

    sorted_chroms = sorted(set(get_chr_to_len(genome).keys()))
    chrom_to_color = get_value_to_color(sorted_chroms, palette="tab20", n=len(sorted_chroms), type_color="hex")[0]

    return sorted_chroms, chrom_to_color


def adjust_cm_positions(cm, df_square, hm_height_multiplier=0.0002, hm_width_multiplier=0.01, cc_height_multiplier=0.017, rc_width_multiplier=0.017, idx_delimiter="-", distance_btw_boxes=0.0025, cd_height=0.07, rd_width=0.07, cbar_width=0.08, cbar_height=0.015, nrow_fields=None, ncol_fields=None):

    """ Adjusts cm positions """


    # define the hm width and height
    hm_height = len(df_square)*hm_height_multiplier
    hm_width = len(df_square.columns)*hm_width_multiplier

    # set heatmap position
    hm_pos = cm.ax_heatmap.get_position()
    hm_y0 = (hm_pos.y0+hm_pos.height)-hm_height
    cm.ax_heatmap.set_position([hm_pos.x0, hm_y0, hm_width, hm_height]); hm_pos = cm.ax_heatmap.get_position()

    # set the col colors position
    if cm.ax_col_colors is not None:
        cc_height = ncol_fields*cc_height_multiplier

        cm.ax_col_colors.set_position([hm_pos.x0, hm_pos.y0+hm_height+distance_btw_boxes, hm_width, cc_height]); cc_pos = cm.ax_col_colors.get_position()
        cd_y0 = cc_pos.y0 + cc_height + distance_btw_boxes

        top_object_pos = cm.ax_col_colors.get_position()

    else: 
        cd_y0 = hm_pos.y0 + hm_height + distance_btw_boxes
        top_object_pos = cm.ax_heatmap.get_position()


    # set the col dendrogram position
    if cm.ax_col_dendrogram is not None: cm.ax_col_dendrogram.set_position([hm_pos.x0, cd_y0, hm_width, cd_height])

    # set the row colors
    if cm.ax_row_colors is not None: 
        rc_width = nrow_fields*rc_width_multiplier

        cm.ax_row_colors.set_position([hm_pos.x0-rc_width-distance_btw_boxes, hm_pos.y0, rc_width, hm_height]); rc_pos = cm.ax_row_colors.get_position()
        rd_x0 = rc_pos.x0 - rd_width - distance_btw_boxes
        leftmost_object_pos = cm.ax_row_colors.get_position()


    else: 
        rd_x0 = hm_pos.x0 - rd_width - distance_btw_boxes
        leftmost_object_pos = cm.ax_heatmap.get_position()

    # set the row dendrogram
    if cm.ax_row_dendrogram is not None: 
        cm.ax_row_dendrogram.set_position([rd_x0, hm_pos.y0, rd_width, hm_pos.height]); rd_pos = cm.ax_row_dendrogram.get_position()


    # define the cbar position
    #cm.ax_cbar.set_position([hm_pos.x0, hm_pos.y0 - cbar_height - distance_btw_boxes, cbar_width, cbar_height]) # on the bottom
    #cm.ax_cbar.set_position([leftmost_object_pos.x0 - distance_btw_boxes*5 - cbar_width, leftmost_object_pos.y0, cbar_width, cbar_height]) # on the bottom left

    cm.ax_cbar.set_position([top_object_pos.x0, top_object_pos.y0 +  top_object_pos.height + distance_btw_boxes*5, cbar_width, cbar_height])

def heatmap_windows_all_strains(strains_tree, filename, df_all, zfield, col_f, aligner, col_colors_dict, row_colors_dicts, df_metadata, col_sorting_fields = ["chromosome", "pos"], col_color_field="chromosome", row_color_fields=["cladeID_based_on_Bergin2022", "source_data"], cmap = "vlag", vmin=0, vmax=2, square_w = 0.01, hm_height_multiplier_by_square_w=3, row_separator_field="cladeID_based_on_Bergin2022"):

    """Gets a heatmap where the rows are strains and the cols are windows"""

    ############ DATA #####

    # keep
    strains_tree = cp.deepcopy(strains_tree)
    df_metadata = cp.deepcopy(df_metadata)
    
    print("getting data")

    # get square df
    df_plot_file = "%s/.%s.data.py"%(get_dir(filename), get_file(filename))
    if file_is_empty(df_plot_file):

        # checks
        df_all = cp.deepcopy(df_all)
        df_all = df_all[df_all.aligner==aligner]
        if set(df_all.runID)!=set(strains_tree.get_leaf_names()): raise ValueError("not same leafs")

        # get data
        df_plot = df_all[df_all[zfield]!=-1].pivot(index="runID", columns=col_f, values=zfield)
        check_no_nans_in_df(df_plot)

        save_object(df_plot, df_plot_file)

    df_plot = load_object(df_plot_file) # .iloc[0:10]

    #####################

    ######## PLOT ########

    print("getting plot")

    # define the sorted columns and indices
    sorted_cols = list(df_all[col_sorting_fields + [col_f]].drop_duplicates().sort_values(by=col_sorting_fields)[col_f])
    sorted_rows = strains_tree.get_leaf_names()
    #sorted_rows = [x for x in strains_tree.get_leaf_names() if x in df_plot.index]
    df_plot = df_plot.loc[sorted_rows, sorted_cols]

    # get the col colors and df samples
    print("getting col colors all fields")
    col_f_to_color_f = dict(df_all[[col_color_field, col_f]].drop_duplicates().set_index(col_f)[col_color_field])
    col_colors_df = pd.DataFrame({col_color_field : {idx : col_colors_dict[col_f_to_color_f[idx]]  for idx in df_plot.columns}}).loc[sorted_cols]

    # get the row colors (multiple labels)
    row_colors_df = pd.DataFrame(index=df_plot.index)
    for Ir, color_dict in enumerate(row_colors_dicts):

        row_color_field = row_color_fields[Ir]
        row_f_to_color_f = dict(df_metadata[[row_color_field, "runID"]].drop_duplicates().set_index("runID")[row_color_field])

        row_colors_df[row_color_field] = [color_dict[row_f_to_color_f[r]] for r in row_colors_df.index]

    row_colors_df = row_colors_df.loc[sorted_rows]

    # clustermap
    print("get clustermap")

    cm = sns.clustermap(df_plot, col_cluster=False, row_cluster=False, row_colors=row_colors_df, col_colors=col_colors_df, cmap=cmap, vmin=vmin, vmax=vmax, xticklabels=False, cbar_kws={"label":zfield, "orientation":'horizontal'}, linewidth=0, yticklabels=1) #, linewidth=.01) # figsize=(30, 10),

    # adjust
    print("adjusting positions")
    distance_btw_boxes = 0.01
    hm_height_multiplier = square_w*hm_height_multiplier_by_square_w
    #hm_height_multiplier  = square_w

    adjust_cm_positions(cm, df_plot, hm_height_multiplier=hm_height_multiplier, hm_width_multiplier=square_w, cc_height_multiplier=hm_height_multiplier, rc_width_multiplier=hm_height_multiplier, cbar_width=hm_height_multiplier*13, cbar_height=hm_height_multiplier*2, distance_btw_boxes=distance_btw_boxes, nrow_fields=len(row_colors_df.columns), ncol_fields=len(col_colors_df.columns))

    # add lines for the chromosome (col sepparator)
    print("adding lines")
    sorted_col_colors_field_vals = sorted(set(df_all[col_sorting_fields[0]]))

    df_all_sorted_cols = df_all[col_sorting_fields + [col_f]].drop_duplicates().sort_values(by=col_sorting_fields)

    for field_v in sorted_col_colors_field_vals[1:]:
        val_x = df_all_sorted_cols[df_all_sorted_cols[col_sorting_fields[0]]==field_v].index[0]
        cm.ax_heatmap.axvline(val_x, color="black", linewidth=1.4, linestyle="--")

    # add lines for each row
    for val_y in range(0, len(df_plot)): 
        color = row_colors_df.loc[sorted_rows[val_y], row_separator_field]
        cm.ax_heatmap.axhline(val_y, color=color, linewidth=.7, linestyle="--")

    # add legend
    def get_legend_element(facecolor, label, edgecolor="gray"): return mpatches.Patch(facecolor=facecolor, edgecolor=edgecolor, label=label)

    legend_elements = []
    for f, color_dict in ([(row_color_fields[Ir], color_dict) for Ir, color_dict in enumerate(row_colors_dicts)] + [(col_color_field, col_colors_dict)]):
        # header
        legend_elements.append(get_legend_element("white", f, edgecolor="white"))

        # each of the values
        for val, color in color_dict.items(): 
            legend_elements.append(get_legend_element(color, val, edgecolor="gray"))

    hm_pos = cm.ax_heatmap.get_position()
    cm.ax_row_colors.legend(handles=legend_elements, bbox_to_anchor=(0,0), loc="lower right") # upper right position of the legend is in 0, 1 (upper leftmost corner)

    print("showing...")
    plt.show()

    print("saving fig")
    cm.savefig(filename, bbox_inches="tight")

    plt.close()

    ######################



def run_jobarray_file_MN5_greasy(jobs_filename, name, time="12:00:00", queue="gp_bscls", threads_per_job=4, nodes=1, submit=True, nsubmits=1, constraint_line=""):

    """
    This function takes a jobs filename and creates a jobscript with args (which is a list of lines to be written to the jobs cript). It works in MN4 for greasy    

    """

    # define dirs
    outdir = get_dir(jobs_filename)
    stddir = "%s/STDfiles"%outdir; 
    delete_folder(stddir); make_folder(stddir)
    os.chdir(outdir)

    # define the std files
    greasy_logfile = "%s/%s_greasy.log"%(stddir, name)
    stderr_file = "%s/%s_stderr.txt"%(stddir, name)
    stdout_file = "%s/%s_stdout.txt"%(stddir, name)

    # define the job script
    jobs_filename_run = "%s.run"%jobs_filename

    # define the number of jobs in the job array
    njobs = int("".join([x for x in str(subprocess.check_output("wc -l %s"%jobs_filename, shell=True)).split()[0] if x.isdigit()])) + 1

    # define the requested nodes. Each node has 112 threads
    total_nthreads = njobs*threads_per_job
    if total_nthreads%112==0: max_nodes = max([int((total_nthreads)/112), 1])
    else: max_nodes = max([int((total_nthreads)/112)+1, 1])
    requested_nodes = min([nodes, max_nodes])

    # define the number of tasks
    max_ntasks = int((requested_nodes*112)/threads_per_job)
    ntasks = min([njobs, max_ntasks])

    # define the arguments
    arguments = [ "#!/bin/sh",
                  "#SBATCH --error=%s"%stderr_file,
                  "#SBATCH --output=%s"%stdout_file,
                  "#SBATCH --job-name=%s"%name, 
                  "#SBATCH --get-user-env",
                  "#SBATCH --account=bsc40",
                  "#SBATCH --time=%s"%time,
                  "#SBATCH --qos=%s"%queue,
                  "#SBATCH --cpus-per-task=%i"%threads_per_job,
                  "#SBATCH --ntasks=%i"%ntasks,
                  constraint_line,
                  "",
                  "module load greasy",
                  "export GREASY_NWORKERS=%i;"%(ntasks + 1),
                  "export GREASY_LOGFILE=%s;"%(greasy_logfile),
                  "rm %s;"%greasy_logfile,
                  "rm %s/job.*.std;"%stddir,
                  "echo 'running pipeline';",
                  "greasy %s"%jobs_filename
                ]

    # replace by the run in MN
    arguments = [a.replace(ParentDir, "/gpfs/projects/bsc40/mschikora") for a in arguments]

    # define and write the run filename
    with open(jobs_filename_run, "w") as fd: fd.write("\n".join(arguments))
    
    # run in cluster if specified
    if submit is True: 
        for I in range(nsubmits): run_cmd("sbatch %s"%jobs_filename_run)

def save_df_as_tab(df, file):

    """Takes a df and saves it as tab"""

    file_tmp = "%s.tmp"%file
    df.to_csv(file_tmp, sep="\t", index=False, header=True)
    os.rename(file_tmp, file)


def generate_df_repeat_masker_file(df_repeats_file, repeat_masker_outfile):

    """Generates the repeat masker df"""

    if file_is_empty(df_repeats_file):
        print("Getting repeats df...")

        # define the header
        header = ["SW_score", "perc_div", "perc_del", "perc_ins", "chromosome", "begin_repeat", "end_repeat", "left_repeat", "strand", "repeat", "type", "position_inRepeat_begin", "position_inRepeat_end", "left_positionINrepeat",  "IDrepeat"]
        function_map_type = [int, float, float, float, str, int, int, str, str, str, str, str, str, str, int]

        # initialize header
        dict_repeat_masker = {}

        # go line by line
        for IDline, line in enumerate(open(repeat_masker_outfile, "r").readlines()):

            # debug lines
            split_line = line.strip().split()
            if len(split_line)==0 or split_line[0] in {"SW", "score"}: continue
            if split_line[-1]=="*": del split_line[-1]

            # keep
            line_content = [function_map_type[I](content) for I, content in enumerate(split_line)]
            dict_repeat_masker[IDline] = dict(zip(header, line_content))

        # get as df
        df_repeats = pd.DataFrame(dict_repeat_masker).transpose().sort_values(by=["chromosome", "begin_repeat", "end_repeat"])
        df_repeats = df_repeats.drop_duplicates(subset=[x for x in df_repeats.keys() if x not in {"IDrepeat"}], keep="first")
        df_repeats["IDrepeat"] = list(range(1, len(df_repeats)+1))
        df_repeats.index = list(range(1, len(df_repeats)+1))

        save_df_as_tab(df_repeats, df_repeats_file)

def get_simple_repeats_file(genome_Cparapsilosis, outdir):

    # get the RepeatMasker outdir

    # define the folder
    repeat_masker_outfile = "%s/%s.out"%(outdir, get_file(genome_Cparapsilosis))
    if file_is_empty(repeat_masker_outfile):

        initial_dir = os.getcwd()
        outdir_tmp = "%s_tmp"%outdir
        delete_folder(outdir_tmp); make_folder(outdir_tmp)
        os.chdir(outdir_tmp)

        run_cmd("RepeatMasker -pa %i -dir %s -poly -html -gff -noint %s"%(112, outdir_tmp, genome_Cparapsilosis), env="perSVade_env_RepeatMasker_env")

        os.chdir(initial_dir)
        os.rename(outdir_tmp, outdir)

    # generate the df repeats file formatted for perSVade
    df_repeats_file = "%s/simple_repeats_perSVade_format.tab"%outdir
    generate_df_repeat_masker_file(df_repeats_file, repeat_masker_outfile)

    return df_repeats_file

def get_df_and_header_from_vcf(vcf_file):

    """Takes a vcf file and returns the df and the header lines (as a list)"""

    # define the header
    header_lines = []
    for line in open(vcf_file, "r", encoding='utf-8', errors='ignore'):
        if line.startswith("##"): header_lines.append(line.strip())
        else: break

    # get the vcf
    df = pd.read_csv(vcf_file, skiprows=len(header_lines), sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)

    return df, header_lines


def get_df_vars_only_in_specific_regions(df_vars, positions_file, df_vars_filt_file):

    """Gets a df with variants that overlap positions_file"""

    if file_is_empty(df_vars_filt_file):

        # tmp
        tmpdir = "%s.generating"%df_vars_filt_file
        delete_folder(tmpdir)
        make_folder(tmpdir)

        # get inital vars
        initial_fields = list(df_vars.keys())

        # get positions of the variants as bed
        df_vars_pos = df_vars[["#CHROM", "POS", "chrom_pos"]].drop_duplicates()
        variants_file = "%s/variant_positions.bed"%tmpdir
        df_vars_pos["start"] = df_vars_pos.POS - 1
        df_vars_pos[["#CHROM", "start", "POS", "chrom_pos"]].sort_values(by=["#CHROM", "start", "POS"]).to_csv(variants_file, sep="\t", index=False, header=False)

        # overlap them with the positions_file (only keep positions that overlap)
        intersection_bed = "%s/intersections_bed.bed"%tmpdir
        run_cmd("bedtools intersect -a %s -b %s -wa > %s"%(variants_file, positions_file, intersection_bed), env="perSVade_env")
        chrom_pos_intersect = set(pd.read_csv(intersection_bed, sep="\t", header=None, names=["#CHROM", "start", "POS", "chrom_pos"]).chrom_pos)
        if len(chrom_pos_intersect.difference(set(df_vars.chrom_pos)))>0: raise ValueError("There are extra chrom pos")

        # keep the overlapping positions
        df_vars_filt = df_vars[df_vars.chrom_pos.isin(chrom_pos_intersect)].reset_index(drop=True)[initial_fields]

        # save
        delete_folder(tmpdir)
        save_object(df_vars_filt, df_vars_filt_file)

    return load_object(df_vars_filt_file)

def get_df_small_vars_filt_one_run(runID, aligner, filt_vcf, repeats_bed, tmpdir, median_gDNA_cov):

    """Get filt variants"""

    # log
    print("getting df for %s-%s"%(runID, aligner))

    # make
    make_folder(tmpdir)

    # load df vars
    df_vars = get_df_and_header_from_vcf(filt_vcf)[0]

    # load fields of FORMAT
    f_to_fun = {"GT":str, "AF":float, "DP":float, "AD":str}
    for I, f in enumerate(df_vars.FORMAT.iloc[0].split(":")):
        df_vars[f] = df_vars.SAMPLE.apply(lambda x: x.split(":")[I]).apply(f_to_fun[f])

    # get the initial vars
    df_vars["#CHROM"] = df_vars["#CHROM"].apply(str)
    df_vars["chrom_pos"] = df_vars["#CHROM"] + "_" + df_vars.POS.apply(str)

    # filter out simple repeats
    df_vars_repeats = get_df_vars_only_in_specific_regions(df_vars, repeats_bed, "%s/df_vars_in_repeats.py"%(tmpdir))
    df_vars = df_vars[~(df_vars.chrom_pos.isin(set(df_vars_repeats.chrom_pos)))].copy()

    # add the coverage of the region
    df_vars["DP_rel_to_gDNA"] = df_vars.DP / median_gDNA_cov

    # add other cases
    df_vars["runID"] = runID
    df_vars["aligner"] = aligner

    # add type var
    chars_ref = set.union(*df_vars.REF.apply(set))
    chars_alt = set.union(*df_vars.ALT.apply(set))
    all_chars = chars_ref.union(chars_alt)
    strange_chars = all_chars.difference({"A", "C", "G", "T"})
    if len(strange_chars)>0: raise ValueError("strange_chars: %s"%strange_chars)

    actg_c = {"A", "C", "G", "T"}
    df_vars["type_var"] = ((df_vars.REF.isin(actg_c)) & (df_vars.ALT.isin(actg_c))).map({True:"SNP", False:"INDEL"})
    check_no_nans_series(df_vars["type_var"])

    # clean
    delete_folder(tmpdir)

    return df_vars

def get_df_coverage_one_pos_and_run(runID, POS, chrom, df_cov):

    """gets a df coverage for one position"""

    df_filt = df_cov[(df_cov.start_pos<=POS) & (df_cov.end>=POS)].copy()
    if len(df_filt)!=1: raise ValueError("should be 1 run")
    df_filt["#CHROM"] = chrom
    df_filt["POS"] = POS
    df_filt["runID"] = runID

    return df_filt[["#CHROM", "POS", "median_coverage_rel_to_gDNA", "percentcovered_1", "runID"]]

def get_filtered_small_vars_df(ProcessedDataDir, VarCallDir, sorted_runs, aligner, simple_repeats_file, homozygous_Corthopsilosis_srrs, df_coverage_per_window_file, cov_windows_size):

    """Gets the filtered small variants into a single df"""

    df_small_vars_filt_file = "%s/small_vars_filt_%s_cov_windows_%ibp.py"%(ProcessedDataDir, aligner, cov_windows_size)
    if file_is_empty(df_small_vars_filt_file):

        # define tmpdir
        tmpdir = "%s.generating"%df_small_vars_filt_file
        make_folder(tmpdir)

        # load df coverage per windows and get fields
        df_coverage = load_object(df_coverage_per_window_file)
        runID_to_median_cov = df_coverage[df_coverage.chromosome!=mtDNA_chrom].groupby("runID").apply(lambda df: np.median(df.mediancov_1))

        # generate repeats bed
        df_repeats = pd.read_csv(simple_repeats_file, sep="\t")
        repeats_bed = "%s/simple_repeats.bed"%tmpdir
        df_repeats["start_0based"] = df_repeats.begin_repeat - 1
        df_repeats[["chromosome", "start_0based", "end_repeat"]].sort_values(by=["chromosome", "start_0based", "end_repeat"]).drop_duplicates().to_csv(repeats_bed, sep="\t", header=False, index=False)

        # get inputs fn
        print("preparing inputs")
        inputs_fn = []
        for runID in sorted_runs:

            # discard some cases were it was not called
            if runID in homozygous_Corthopsilosis_srrs and aligner!="bwa_mem": continue

            # get data
            outdir_aligner = "%s/%s/aligner-%s"%(VarCallDir, runID, aligner)
            filt_vcf = "%s/small_vars/variants_atLeast2PASS_ploidy2.vcf"%outdir_aligner
            inputs_fn.append((runID, aligner, filt_vcf, repeats_bed, "%s/tmp_%s_%s"%(tmpdir, runID, aligner), runID_to_median_cov[runID]))

        # running 
        print("run in parallel")
        with multiproc.Pool(multiproc.cpu_count()) as pool:
            df_small_vars_filt = pd.concat(pool.starmap(get_df_small_vars_filt_one_run, inputs_fn, chunksize=1)) 
            pool.close()
            pool.terminate()

        # add the coverage of the region info
        print("adding cov info")
        df_coverage["window_ID"] = df_coverage["chromosome"] + ":" + df_coverage.start.apply(str) + "-" + df_coverage.end.apply(str)

        # get bed of windows of coverage
        bed_fields = ["chromosome", "start", "end"]
        bed_windows = "%s/windows_coverage.bed"%tmpdir
        df_coverage[bed_fields + ["window_ID"]].drop_duplicates().sort_values(by=bed_fields).to_csv(bed_windows, sep="\t", index=False, header=False)

        # get bed of positions
        df_positions_vars = df_small_vars_filt[["#CHROM", "POS", "chrom_pos"]].drop_duplicates().rename(columns={"#CHROM":"chromosome", "POS":"end"})
        df_positions_vars["start"] = df_positions_vars.end - 1
        df_positions_vars = df_positions_vars.sort_values(by=bed_fields)
        bed_positions = "%s/positions_vars.bed"%tmpdir
        df_positions_vars[bed_fields + ["chrom_pos"]].to_csv(bed_positions, sep="\t", index=False, header=False)

        # map with bedmap
        print("running bedmap")
        bedmap_outfile = "%s/outfile_bedmap.txt"%tmpdir
        run_cmd("bedmap --delim '\t' --bp-ovr 1 --echo-map-id %s %s > %s"%(bed_positions, bed_windows, bedmap_outfile), env="perSVade_env")
        df_positions_vars["window_ID"] = list(map(lambda l: l.strip(), open(bedmap_outfile, "r")))
        strange_windows = set(df_positions_vars["window_ID"]).difference(set(df_coverage.window_ID))
        if len(strange_windows)>0: raise ValueError("strange_windows: %s"%strange_windows)

        # add to df
        print("add to df")
        df_small_vars_filt = df_small_vars_filt.merge(df_positions_vars[["chrom_pos", "window_ID"]], on="chrom_pos", how="left", validate="many_to_one")
        fields_cov = ["median_coverage_rel_to_gDNA", "percentcovered_1", "median_mappability"]
        df_small_vars_filt = df_small_vars_filt.merge(df_coverage[["window_ID", "runID"] + fields_cov], on=["window_ID", "runID"], how="left", validate="many_to_one").rename(columns={f : "%s_%skb_window"%(f, cov_windows_size/1000) for f in fields_cov})
        #df_small_vars_filt = df_small_vars_filt.drop(["window_ID"], axis=1)
        check_no_nans_in_df(df_small_vars_filt)

        # add info about biallelic vars
        print("add biallelic info")
        all_pos = set(df_small_vars_filt.chrom_pos)

        for tvar, types_var in [("INDEL", {"INDEL"}), ("all", {"SNP", "INDEL"}), ("SNP", {"SNP"})]:
            print(tvar)
            df_v = df_small_vars_filt[df_small_vars_filt.type_var.isin(types_var)]
            if len(df_v)==0: raise ValueError("should be some vars")

            pos_to_nvars = df_small_vars_filt[["chrom_pos", "ID"]].drop_duplicates().groupby("chrom_pos").apply(len)
            missing_series = pd.Series(0, index=sorted(all_pos.difference(set(pos_to_nvars.index))))
            pos_to_nvars = pd.concat([pos_to_nvars, missing_series])

            df_small_vars_filt["nvars_%s_in_pos"%tvar] = df_small_vars_filt.chrom_pos.map(pos_to_nvars)
            check_no_nans_series(df_small_vars_filt["nvars_%s_in_pos"%tvar])

        # add info about mappability
        print("adding mappability info")
        mappability_df = pd.read_csv("%s/mappability_outdir/genome.fasta.mappability_per_position.bed"%ProcessedDataDir, sep="\t").rename(columns={"chromosome":"#CHROM"})
        mappability_df["POS"] = mappability_df.position + 1
        df_small_vars_filt = df_small_vars_filt.merge(mappability_df[["#CHROM", "POS", "unique_map_score"]], on=["#CHROM", "POS"], how="left", validate="many_to_one")

        # checks
        check_no_nans_in_df(df_small_vars_filt)

        # save
        print("saving")
        delete_folder(tmpdir)
        save_object(df_small_vars_filt, df_small_vars_filt_file)


    return df_small_vars_filt_file

def index_genome(genome):

    """Takes a fasta and generates a <genome>.fai file"""

    # index the genome of interest if not already done
    if file_is_empty("%s.fai"%genome): 

        faidx_std = "%s.indexing.std"%genome
        #print("running faidx. The std is in %s"%faidx_std)
        run_cmd("samtools faidx %s > %s 2>&1"%(genome, faidx_std), env="perSVade_env")
        remove_file(faidx_std)


def get_windows_for_reference_genome(reference_genome, window_size):

    """This function generates a window file for the reference genome"""

    windows_file = "%s.windows%ibp.bed"%(reference_genome, window_size)
    if file_is_empty(windows_file):

        # index the genome
        index_genome(reference_genome)

        # get the windows
        windows_file_stderr = "%s.generating.stderr"%windows_file
        windows_file_tmp = "%s.tmp"%windows_file

        #print("genearting windows_file. The stderr is in %s"%windows_file_stderr)
        run_cmd("bedtools makewindows -g %s.fai -w %i > %s 2>%s"%(reference_genome, window_size, windows_file_tmp, windows_file_stderr), env="perSVade_env")

        remove_file(windows_file_stderr)

        # sort by chrom and pos
        colnames = ["chromosome", "start", "end"]
        df_windows = pd.read_csv(windows_file_tmp, sep="\t", header=None, names=colnames).sort_values(by=colnames)
        df_windows[colnames].to_csv(windows_file_tmp, sep="\t", header=False, index=False)

        # clean
        os.rename(windows_file_tmp, windows_file)

    return windows_file


def get_df_var_density_one_window(df_vars, data_dict, sorted_runs):

    """gets stats about var density for one window"""

    print("running on %s"%data_dict["ID_window"])

    # define the window len
    window_l = data_dict["end"] - data_dict["start"]

    # init
    df_density = pd.DataFrame()

    # go through different types of zigosities to calculate 
    strange_zigosities = set(df_vars.GT).difference({"0/1", "1/1", "."})
    if len(strange_zigosities)>0: raise ValueError("strange_zigosities: %s"%strange_zigosities)

    for zygosity, gt_val in [("hetero", "0/1"), ("homo", "1/1"), ("unknown_zygosity", ".")]:

        # get df of the zigosity
        df_zyg = df_vars[df_vars.GT==gt_val].copy()

        # init df with sorted runs
        missing_runs = sorted(set(sorted_runs).difference(set(df_zyg.runID)))
        df_d = pd.DataFrame({r : {"vars_per_kb":0.0, "mean_AF":-1, "nvars":0, "mean_DP_rel":-1, "runID":r, "AF_list":np.array([]), "set_numeric_vars":set()} for r in missing_runs}).transpose().reset_index(drop=True)

        # add the real runs
        if len(df_zyg)>0:

            def get_series_one_run_df_density(df):
                dict_r = {"runID":df.runID.iloc[0]}
                nvars = len(set(df.ID))
                dict_r["nvars"] = nvars
                dict_r["vars_per_kb"] = nvars / window_l
                dict_r["mean_AF"] = np.mean(df.AF)
                dict_r["mean_DP_rel"] = np.mean(df.DP_rel_to_gDNA)
                dict_r["AF_list"] = df.AF.values
                dict_r["set_numeric_vars"] = set(df.numeric_ID)
                return pd.Series(dict_r)

            df_d = pd.concat([df_d, df_zyg.groupby("runID").apply(get_series_one_run_df_density)]).reset_index(drop=True)

        # checks
        check_no_nans_in_df(df_d)
        if sorted(df_d.runID)!=sorted_runs: raise ValueError("runs not as expected")

        # add
        df_d["zygosity"] = zygosity
        df_density = pd.concat([df_density, df_d]).reset_index(drop=True)

    # add all the data_dict things
    for k,v in data_dict.items():
        df_density[k] = v

    return df_density



def get_df_coverage_with_fraction_repeats_per_window(df_coverage, repeats_table, tmpdir):

    """filters the repeats table"""

    # define the final fields
    final_fields = list(df_coverage.keys()) + ["fraction_repeats"]

    # define the file of cov with added repeats
    df_coverage_file = "%s/df_coverage_added_repeats.py"%tmpdir
    if file_is_empty(df_coverage_file):
        print("adding repeat content per window")

        # create tmp
        delete_folder(tmpdir)
        make_folder(tmpdir)

        # keep original len
        original_cov_len = len(df_coverage) 

        # load df of repeats
        repeats_df = pd.read_csv(repeats_table, sep="\t").rename(columns={"begin_repeat":"start", "end_repeat":"end"})
        repeats_df["start"] = repeats_df.start -1 
        repeats_df["chromosome"] = repeats_df.chromosome.apply(str)
        if any(repeats_df.start<0): raise ValueError("there can't be negative starts")
        if not all(repeats_df.start<repeats_df.end): raise ValueError("end has to be after start")

        # generate bed with repeat regions
        repeats_bed = "%s/repeats.bed"%tmpdir
        bed_fields = ["chromosome", "start", "end"]
        repeats_df[bed_fields].sort_values(by=bed_fields).to_csv(repeats_bed, sep="\t", header=False, index=False) # ["unique_repID", "rep_score"]

        # add the fraction of repeats as a bed
        coverage_bed = "%s/coverage.bed"%tmpdir
        df_coverage = df_coverage.sort_values(by=bed_fields)
        df_coverage[bed_fields].to_csv(coverage_bed, sep="\t", header=False, index=False)

        bedmap_outfile = "%s/outfile_bedmap.txt"%tmpdir
        bedmap_stderr = "%s/bedmap.stderr"%tmpdir
        run_cmd("bedmap --delim '\t' --bp-ovr 1 --bases-uniq-f %s %s > %s 2>%s"%(coverage_bed, repeats_bed, bedmap_outfile, bedmap_stderr), env="perSVade_env")

        # add the fraction 
        df_coverage["fraction_repeats"] = list(map(float, open(bedmap_outfile, "r")))

        # checks
        if any(df_coverage.fraction_repeats>1): raise ValueError("fraction can't be >1")
        if len(df_coverage)!=original_cov_len: raise ValueError("len changed")
        check_no_nans_in_df(df_coverage)

        # save
        save_object(df_coverage, df_coverage_file)

    # load
    return load_object(df_coverage_file)[final_fields]


def generate_variant_per_window_densities_df(df_small_vars_filt_file, ProcessedDataDir, genome_Cparapsilosis, aligner, window_size, sorted_runs, simple_repeats_file, df_coverage_per_window_file, excluded_chroms=set()):

    """Generates a df with variant densities for different variant types"""

    # define files
    def get_excluded_chroms_str():
        if excluded_chroms==set(): return "none"
        elif len(excluded_chroms)==1: return next(iter(excluded_chroms))
        else: return "%i_%s-%s"%(len(excluded_chroms), sorted(excluded_chroms)[0], sorted(excluded_chroms)[1])

    df_var_density_dir = "%s/small_variant_density_dfs"%ProcessedDataDir; make_folder(df_var_density_dir)
    df_var_density_file = "%s/small_variant_density_df_%s_w_%i_excluded_chroms_%s_%i_runs_%s-%s.py"%(df_var_density_dir, aligner, window_size, get_excluded_chroms_str(), len(sorted_runs), sorted_runs[0], sorted_runs[-1])

    if file_is_empty(df_var_density_file):
        print("generate var density")

        # get tmp
        tmpdir = "%s.generating"%df_var_density_file
        delete_folder(tmpdir)
        make_folder(tmpdir)

        # generate windows file
        print("getting windows...")
        windows_bed = get_windows_for_reference_genome(genome_Cparapsilosis, window_size)
        df_windows = pd.read_csv(windows_bed, sep="\t", header=None, names=["chromosome", "start", "end"])
        df_windows["ID_window"] = df_windows.chromosome + ":" + df_windows.start.apply(str) + "-" + df_windows.end.apply(str)

        # discard chroms
        df_windows = df_windows[~df_windows.chromosome.apply(str).isin(excluded_chroms)].copy()

        # add data from df_coverage_per_window
        print("getting cov per window")
        df_coverage_per_window = load_object(df_coverage_per_window_file)
        df_coverage_per_window["ID_window"] = df_coverage_per_window.chromosome + ":" + df_coverage_per_window.start.apply(str) + "-" + df_coverage_per_window.end.apply(str)
        check_no_nans_in_df(df_coverage_per_window)

        # discard chroms
        df_coverage_per_window = df_coverage_per_window[~df_coverage_per_window.chromosome.apply(str).isin(excluded_chroms)].copy()

        # keep only target runs
        df_coverage_per_window = df_coverage_per_window[df_coverage_per_window.runID.isin(set(sorted_runs))].copy()

        # check all correct windows
        missing_windows = set(df_windows.ID_window).difference(set(df_coverage_per_window.ID_window))
        if len(missing_windows)>0: raise ValueError("missing_windows: %s"%missing_windows)
        if set(df_windows.ID_window)!=set(set(df_coverage_per_window.ID_window)): raise ValueError("windows should be the same")
        if set(sorted_runs)!=set(df_coverage_per_window.runID): raise ValueError("runIDs not as expected")

        # set idx
        df_windows = df_windows.set_index("chromosome", drop=False)

        # get vars df
        print("loading df...")
        df_vars_file = "%s/df_vars.py"%tmpdir
        if file_is_empty(df_vars_file):

            # filter
            fields_vars = ['#CHROM', 'POS', 'ID', 'GT', 'AF', 'DP', 'chrom_pos', 'DP_rel_to_gDNA', 'runID', 'type_var']
            df_vars = load_object(df_small_vars_filt_file)
            df_vars = df_vars[df_vars.runID.isin(set(sorted_runs))].copy()
            df_vars = df_vars.set_index("runID", drop=False)
            df_vars = df_vars[~df_vars["#CHROM"].apply(str).isin(excluded_chroms)].copy()
            df_vars = df_vars[df_vars.aligner==aligner][fields_vars].copy()

            # add the numeric variant IDs
            sorted_vars = sorted(set(df_vars.ID))
            var_to_numericID = dict(zip(sorted_vars, range(len(sorted_vars))))
            df_vars["numeric_ID"] = df_vars.ID.map(var_to_numericID)
            check_no_nans_series(df_vars.numeric_ID)

            # save
            save_object(df_vars, df_vars_file)

        df_vars = load_object(df_vars_file).set_index("#CHROM")

        # init
        df_var_density = pd.DataFrame()
    
        # go through different types
        for type_var in ["INDEL", "SNP"]:

            # keep
            df_vars_v = df_vars[df_vars.type_var==type_var].copy()

            # prepare inputs
            inputs_fn = []
            for chrom in sorted(set(df_vars.index)):
                print("prepare inputs %s, %s"%(type_var, chrom))

                # get df of this chrom, also the windows
                df_v = df_vars_v.loc[[chrom]].copy()
                df_windows_c = df_windows.loc[[chrom]].copy()

                # get the df windows of the chrom
                inputs_fn_c = list(df_windows_c.apply(lambda r: (df_v[(df_v.POS>r.start) & (df_v.POS<=r.end)].copy(), dict(r), sorted_runs), axis=1))
                inputs_fn += inputs_fn_c

            # get in parallel the variant info tables
            print("running in parallel...")
            with multiproc.Pool(multiproc.cpu_count()) as pool:
                df_var_density_v = pd.concat(pool.starmap(get_df_var_density_one_window, inputs_fn, chunksize=1)) 
                pool.close()
                pool.terminate()

            # add the df_coverage_per_window data
            print("adding cov per window data")
            df_var_density_v = df_var_density_v.merge(df_coverage_per_window[["ID_window", "runID", 'mediancov_1', 'percentcovered_1', 'fraction_repeats', 'GCcontent', 'median_mappability', 'fraction_N_bases', 'median_coverage_rel_to_gDNA']], on=["ID_window", "runID"], validate="many_to_one", how="left")

            check_no_nans_in_df(df_var_density_v)

            # add fields
            print("keep")
            df_var_density_v["window_size"] = window_size
            df_var_density_v["aligner"] = aligner
            df_var_density_v["type_var"] = type_var
            df_var_density = pd.concat([df_var_density, df_var_density_v]).reset_index(drop=True)

        # save
        delete_folder(tmpdir)
        save_object(df_var_density, df_var_density_file)

    return df_var_density_file


def get_runIDs_one_run_set(df_reads_info, run_set):

    """gets the sorted run IDs"""

    # define the ortho runs
    runIDs_ortho = sorted(set(df_reads_info[df_reads_info.species_run=="ortho"].runID))

    # define the para runs
    df_para = df_reads_info[df_reads_info.species_run=="para"].copy()
    if run_set=="all":
        runIDs_para = sorted(set(df_para.runID))

    elif run_set=="exclude_triploids":
        runIDs_para = sorted(set(df_para[df_para.possibly_triploid==False].runID))

    else: raise ValueError("invalid")

    return runIDs_ortho, runIDs_para

def get_jaccard_index(x1, x2):

    """get jaccard dist between sets"""

    if len(x1.union(x2))==0: raise ValueError("0 sig terms")
    return len(x1.intersection(x2)) / len(x1.union(x2))


def get_outgroup_rooted_tree(tree, outgroup):

    """gets an outgroup rooted tree"""

    # definethe leafs
    outg_leafs = set(outgroup.split(","))
    if len(outgroup.split(","))!=len(outg_leafs): raise ValueError("leaf nodes should be unique")

    # check that they are all in the tree
    stange_leafs = outg_leafs.difference(set(tree.get_leaf_names()))
    if len(stange_leafs)>0: raise ValueError("stange_leafs: %s"%stange_leafs)

    # root for only one
    if len(outg_leafs)==1: outg_node = outgroup

    # root for multiple
    else:

        # map each node to the similarity to the outgroup node
        similarity_df = pd.DataFrame({I : {"nodeI": I, "node_obj": n, "similarity": get_jaccard_index(outg_leafs, set(n.get_leaf_names())), "node_len": len(n.get_leaf_names()), "node_str":",".join(sorted(n.get_leaf_names())), "support":int(n.support)} for I,n in enumerate(tree.traverse())}).transpose()

        # get the best node, closest to the outgroup
        best_node_series = similarity_df.sort_values(by=["similarity", "node_len", "support", "node_str"], ascending=[False, True, False, True]).iloc[0]
        outg_node = best_node_series.node_obj

        # checks
        if best_node_series.similarity!=1.0: 
            print("WARNING: We could not find the provided outgroup %s. Setting the following one:\n%s"%(outgroup, outg_node))

    # set outgroup
    tree.set_outgroup(outg_node)

    return tree


def add_tree_cmds_to_list(cmd_tuple, list_cmds, outdir_tree_generation, mode, threads_tree, mem_Gb, skip_before_resample_running=False):

    """adds tree generation cmd"""

    if file_is_empty("%s/tree_generated.txt"%outdir_tree_generation):

        # for mode=="diploid_homozygous", just get the cmd
        if mode=="diploid_homozygous": 
            list_cmds.append(cmd_tuple)

        # for the diploid mode, it needs to be done in different phases
        elif mode=="diploid": 

            # phase 1, generate all the resample jobs by running normally
            jobs_file_resamples = "%s/jobs_tree_generation.txt"%outdir_tree_generation
            if file_is_empty(jobs_file_resamples):
                list_cmds.append(cmd_tuple)

            else:

                # skip
                if skip_before_resample_running is True: return

                # get the cmds to run
                cmds_resmples_to_run = []
                for l in open(jobs_file_resamples, "r").readlines():

                    # get the cmd
                    cmd = l.strip().split(">")[0] 
                    outdir_cmd = cmd.split("--outdir")[1].split()[0]

                    # modify resources
                    prev_arg = ""
                    arg_list = []
                    for arg in cmd.split():
                        if prev_arg=="--threads": arg_list.append(str(threads_tree))
                        elif prev_arg=="--iqtree_mem": arg_list.append(str(mem_Gb))
                        else: arg_list.append(arg)
                        prev_arg = arg

                    cmd = " ".join(arg_list)

                    # keep
                    if file_is_empty("%s/tree_was_generated.txt"%outdir_cmd): cmds_resmples_to_run.append(cmd)

                # phase 2, run the individual jobs
                if len(cmds_resmples_to_run)>0:
                    list_cmds += [("trees_resamples", c, "tree_from_SNPs_env", threads_tree) for c in cmds_resmples_to_run]

                # phase 3, generate final tree, by running the command
                else: 
                    list_cmds.append(cmd_tuple)

        else: raise ValueError("invalid mode")


def get_cmds_tree_reconstruction_one_configuration(runIDs_para, run_set, mode, outgroup, aligner, VarCallDir, outdir_trees_all, threads, mem_Gb, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, list_cmds, threads_tree, min_support_trees, skip_before_resample_running=False):

    """gets the cmds to run a tree in a specific configuration"""

    # define samples, depending on the outgroup
    runIDs_tree = cp.deepcopy(runIDs_para)
    if outgroup in homozygous_Corthopsilosis_srrs: 
        outgroup_str = outgroup
        runIDs_tree.append(outgroup)
        type_outgroup = "Corthopsilosis"
    
    elif len(outgroup.split(","))==1: 
        outgroup_str = outgroup
        type_outgroup = "Cparapsilosis"

    elif len(outgroup.split(","))>1:
        split_outg = sorted(outgroup.split(","))
        outgroup_str = "%i_strains_%s_to_%s"%(len(split_outg), split_outg[0], split_outg[-1])
        type_outgroup = "Cparapsilosis"

    else: raise ValueError("invalid outgroup")

    # define outdirs
    outdir_tree = "%s/runs-%s_mode-%s_aligner-%s_outgroup-%s"%(outdir_trees_all, run_set, mode, aligner, outgroup_str)
    outdir_tree_generation = "%s/tree_generation"%outdir_tree
    make_folder(outdir_tree)

    # getting paths table
    df_paths_dict = {}
    for runID in runIDs_tree:
        outdir_aligner = "%s/%s/aligner-%s"%(VarCallDir, runID, aligner)
        df_paths_dict[runID] = {"sampleID": runID, "sorted_bam":"%s/align_reads/aligned_reads.bam.sorted"%outdir_aligner, "diploid_vcf":"%s/small_vars/variants_atLeast2PASS_ploidy2.vcf"%outdir_aligner}

    df_paths = pd.DataFrame(df_paths_dict).transpose()
    paths_table = "%s/paths_table.tab"%outdir_tree
    df_paths.to_csv(paths_table, sep="\t")

    # check that all files exist
    for f in ["sorted_bam", "diploid_vcf"]:
        if any(df_paths[f].apply(file_is_empty)): raise ValueError("there are empty files")

    # define the backbone cmd
    cmd_tuple = ("tree_from_SNPs", "%s/scripts/tree_from_SNPs/scripts/get_tree.py --paths_table %s -o %s --mode %s --min_coverage_pos %i --batch_mode --reference_genome %s --mode_fasta_gen python --threads %i --iqtree_mem %s --outgroup_resampled_trees %s --batch_mode --repeat_tree_integration"%(ParentDir, paths_table, outdir_tree_generation, mode, min_coverage_pos, genome_Cparapsilosis, int(threads/4), mem_Gb, outgroup), "tree_from_SNPs_env", threads)

    add_tree_cmds_to_list(cmd_tuple, list_cmds, outdir_tree_generation, mode, threads_tree, mem_Gb, skip_before_resample_running=skip_before_resample_running)

    # return empty
    if len(list_cmds)>0: return

    # get the rooted tree (depends on the type), and the parsimony informative sites
    if mode=="diploid_homozygous":  
        rooted_tree = get_outgroup_rooted_tree(Tree("%s/generate_tree_from_SNPs/iqtree_unroted.treefile"%(outdir_tree_generation)), outgroup)
        tree_file = "%s/tree_rooted_by_%s.nw"%(outdir_tree_generation, outgroup_str)
        rooted_tree.write(outfile=tree_file, format=2)
        iqtree_list = ["%s/generate_tree_from_SNPs/iqtree_unroted.iqtree"%(outdir_tree_generation)]
        rooted_trees_list = [Tree(tree_file)]

    elif mode=="diploid":
        tree_file = "%s/rooted_tree_mode=%s_min_cov=%i.nw"%(outdir_tree_generation, mode, min_coverage_pos)
        iqtree_list = ["%s/tree_generation_logs/resample%i.iqtree"%(outdir_tree_generation, I) for I in range(1, 100+1)]
        rooted_trees_list = [Tree(l.strip()) for l in open("%s/resampled_trees.txt"%outdir_tree_generation, "r").readlines()]

    all_parsimony_informative_sites = [get_n_parsimony_informative_sites_from_iqtreefile(x) for x in iqtree_list]

    # get the dict with tree information
    str_pars_sites = "-".join([str(min(all_parsimony_informative_sites)), str(max(all_parsimony_informative_sites))])
    tree_obj = Tree(tree_file)
    tree_corrected = get_correct_tree(tree_file, min_support=min_support_trees)
    tree_info_dict = {"run_set" : run_set, 
                      "mode_tree" : mode, 
                      "outgroup_str" : outgroup_str, 
                      "aligner" : aligner, 
                      "tree_file" : tree_file,
                      "tree" : cp.deepcopy(tree_obj), 
                      "tree_corrected" : cp.deepcopy(tree_corrected), 
                      "nsites_parsimony" : str_pars_sites, 
                      "type_outgroup" : type_outgroup, 
                      "Cpara_all_leaves" : set(runIDs_para)}

    # check rooted
    if len(tree_corrected.get_children())!=2: raise ValueError("should be 2 nodes")

    # identify the most frequent outgroup
    all_outgroups = [",".join(sorted([tuple(sorted(c.get_leaf_names())) for c in t.get_children()], key=len)[0]) for t in rooted_trees_list]
    outg_to_ntrees = pd.Series(Counter(all_outgroups))
    frequent_outgroup = outg_to_ntrees.sort_values(ascending=False).index[0]

    # get the fraction of roots that are like the outgroup
    for target_outgroup, target_outgroup_name in [(outgroup, "outgroup"), (frequent_outgroup, "most_frequent_outgroup")]:

        target_outgroup_leaves = set(target_outgroup.split(","))
        tree_info_dict["fraction_trees_expected_%s"%target_outgroup_name] = sum([any([set(n.get_leaf_names())==target_outgroup_leaves for n in t.get_children()]) for t in rooted_trees_list]) / len(rooted_trees_list)


    # define the fraction of collapsed nodes
    nnodes_total = len([n for n in tree_obj.traverse() if not n.is_leaf()])
    nnodes_collapsed = nnodes_total - len([n for n in tree_corrected.traverse() if not n.is_leaf()])
    tree_info_dict["pct_collapsed_nodes"] = (nnodes_collapsed / nnodes_total)*100

    # define the fraction of collapsed nodes that have a long brach len
    max_distance = max(list(map(lambda n: n.get_farthest_node()[1], tree_obj.get_leaves())))
    nnodes_collapsed_long_branches = len([n for n in tree_obj.traverse() if not n.is_leaf() and n.support<min_support_trees and n.dist>=(max_distance*0.01)])
    tree_info_dict["pct_collapsed_long_nodes"] = (nnodes_collapsed_long_branches / nnodes_total)*100

    # define the fraction of politomies that have a lot of diversity inside the nodes
    npolitomies_diverse = len([n for n in tree_corrected.traverse() if not n.is_leaf() and len(n.get_children())>2 and any([len(c.get_leaf_names())>=3 for c in n.get_children()])])
    tree_info_dict["pct_politiomies_diverse"] = (npolitomies_diverse / nnodes_total)*100

    # add things depending on the type
    if type_outgroup=="Corthopsilosis":
        tree_info_dict["Cpara_outgroup_leaves"] = get_Cpara_outgroup_leaves_Cortho_outgroups(tree_corrected, outgroup)

    elif type_outgroup=="Cparapsilosis":
        tree_info_dict["Cpara_outgroup_leaves"] = sorted([set(n.get_leaf_names()) for n in tree_corrected.get_children()], key=len)[0]

    else: raise ValueError("invalid")

    return tree_info_dict

def get_Cpara_outgroup_leaves_Cortho_outgroups(tree_corrected, outgroup):

    """gets the set of Cpara outgroup"""

    # checks
    if len(tree_corrected.get_children())!=2: raise ValueError("should be 2 trees")

    # check that one of the outgroups in the ortho one
    if sum([n.name==outgroup for n in tree_corrected.get_children()])!=1: raise ValueError("should be one case of ortho")

    # get the leafs of the para tree
    para_trees = [n for n in tree_corrected.get_children() if n.name!=outgroup]
    if len(para_trees)!=1: raise ValueError("should be only one para tree")
    para_tree = para_trees[0]
    children_leafs = [set(n.get_leaf_names()) for n in para_tree.get_children()]

    return sorted(children_leafs, key=len)[0]

def get_jaccard_index_btw_sets(x1, x2, vals_to_consider="all"):

    """get jaccard dist between sets"""

    if vals_to_consider!="all": 
        x1 = x1.intersection(vals_to_consider)
        x2 = x2.intersection(vals_to_consider)

    if len(x1.union(x2))==0: raise ValueError("0 sig terms")
    return len(x1.intersection(x2)) / len(x1.union(x2))

def plot_similarities_outgroups_different_trees(trees_df, PlotsDir):

    """Plots the outgroup similarities between leafs"""

    print("plotting differences outgroups")

    # get a df long with all comparisons
    dict_distance = {}
    recorded_tuples = set()
    for idx1 in trees_df.index: 
        for idx2 in trees_df.index:
            #if idx1==idx2: continue
            tuple_idxs = tuple(sorted([idx1, idx2]))
            if tuple_idxs in recorded_tuples: continue
            recorded_tuples.add(tuple_idxs)

            all_leaves = trees_df.loc[idx1, "Cpara_all_leaves"].intersection(trees_df.loc[idx2, "Cpara_all_leaves"])
            jaccard_idx = get_jaccard_index_btw_sets(trees_df.loc[idx1, "Cpara_outgroup_leaves"], trees_df.loc[idx2, "Cpara_outgroup_leaves"], vals_to_consider=all_leaves)
            dict_distance[tuple_idxs] = {"idx1":idx1, "idx2":idx2, "jaccard_idx":jaccard_idx}

    df_distances = pd.DataFrame(dict_distance).transpose().reset_index(drop=True)

    # dfs
    df_plot = df_distances.pivot(index="idx1", columns="idx2", values="jaccard_idx")

    outg_to_color = get_value_to_color(sorted(set(trees_df.outgroup_str)), palette="tab10")[0]
    run_set_to_color = {"all":"black", "exclude_triploids":"gray"}
    mode_to_color = {"diploid":"blue", "diploid_homozygous":"cyan"}

    row_colors_df = col_colors_df = trees_df.apply(lambda r: pd.Series({"idx":r.name, "run_set": run_set_to_color[r.run_set], "mode_tree": mode_to_color[r.mode_tree], "outgroup": outg_to_color[r.outgroup_str]}), axis=1).set_index("idx")

    # plot
    mask = df_plot.isnull()
    cm = sns.clustermap(df_plot, row_cluster=False, col_cluster=False, yticklabels=True, xticklabels=True, row_colors=row_colors_df, col_colors=row_colors_df, cbar_kws={'label': "jaccard index C. para. outg."}, cmap="rocket_r", linecolor="gray", linewidth=.05, mask=mask)

    tick_labels =  list(trees_df.apply(lambda r: ",".join(["%s=%s"%(f[0], r[f]) for f in ["run_set", "mode_tree", "outgroup_str"]]), axis=1))

    cm.ax_heatmap.set_yticklabels(tick_labels)
    cm.ax_heatmap.set_xticklabels(tick_labels, rotation=90)
    cm.ax_heatmap.set_ylabel("")
    cm.ax_heatmap.set_xlabel("")

    # save
    save_fig(cm, "comparisons_Cpara_outgroups.pdf", PlotsDir)


def plot_treeQC_stats(trees_df, PlotsDir, aligner_to_color):

    """Plots the tree stats of QC"""


    # define stats
    list_stats = ["fraction_trees_expected_outgroup", "min_nsites_parsimony", "pct_collapsed_nodes", "pct_collapsed_long_nodes", "pct_politiomies_diverse"] 

    # redefine
    trees_df = trees_df.copy()
    trees_df["min_nsites_parsimony"] = trees_df.nsites_parsimony.apply(lambda x: min(list(map(int, x.split("-")))))

    # get plot
    nrows = 1
    ncols = len(list_stats)
    fig = plt.figure(figsize=(ncols*2, nrows*3))

    for Ic, stats_f in enumerate(list_stats):

        ax = plt.subplot(nrows, ncols, Ic+1)
        ax = sns.barplot(x="mode_tree", y=stats_f, hue="aligner", palette=aligner_to_color, data=trees_df, hue_order=["bwa_mem", "bowtie2_local", "hisat2_no_spliced"])
        plt.xticks(rotation=90)

        if Ic==2:
            ax.set_xlabel("tree mode")

        else:
            ax.set_xlabel("")

        if Ic==(len(list_stats)-1):
            ax.legend(bbox_to_anchor=(1,1), markerscale=1, title="aligner")

        else:
            ax.get_legend().remove()


    sns.despine()
    plt.subplots_adjust(wspace=0.75, hspace=0.05)# wspace (between cols) hspaces (btw rows) 
    save_fig(fig, "tree_QC_metrics.pdf", PlotsDir)
 


def heatmap_differences_trees(trees_df, PlotsDir, aligner_to_color):

    """views differences between trees"""

    # keep
    trees_df = trees_df.sort_values(by=["type_outgroup", "run_set", "mode_tree", "aligner"]).copy().reset_index(drop=True)

    # get a df long with all comparisons
    print("getting the tree")
    dict_distance = {}
    recorded_tuples = set()
    for idx1 in trees_df.index: 
        for idx2 in trees_df.index:
            tuple_idxs = tuple(sorted([idx1, idx2]))
            if tuple_idxs in recorded_tuples: continue
            recorded_tuples.add(tuple_idxs)

            all_leaves = trees_df.loc[idx1, "Cpara_all_leaves"].intersection(trees_df.loc[idx2, "Cpara_all_leaves"])
            dist1, dist2 = get_pct_different_nodes_two_trees(trees_df.loc[idx1, "tree_corrected"], trees_df.loc[idx2, "tree_corrected"], common_leaves=all_leaves)
            min_similarity = 100 - max([dist1, dist2])
            dict_distance[tuple_idxs] = {"idx1":idx1, "idx2":idx2, "tree similarity": min_similarity}

            # print some cases
            df = trees_df.loc[[idx1, idx2]]
            if set(df.mode_tree)=={"diploid"} and set(df.run_set)=={"exclude_triploids"} and idx1!=idx2:
                print("diploid", trees_df.loc[idx1, "aligner"], "vs", trees_df.loc[idx2, "aligner"], "%.2f"%min_similarity)

            if set(df.mode_tree)=={"diploid_homozygous"} and set(df.run_set)=={"exclude_triploids"} and idx1!=idx2:
                print("diploid_homozygous", trees_df.loc[idx1, "aligner"], "vs", trees_df.loc[idx2, "aligner"], "%.2f"%min_similarity)

            if set(df.mode_tree)=={"diploid_homozygous", "diploid"} and set(df.run_set)=={"exclude_triploids"} and idx1!=idx2:
                print("diploid_homozygous vs diploid", trees_df.loc[idx1, "aligner"], "vs", trees_df.loc[idx2, "aligner"], "%.2f"%min_similarity)


    df_distances = pd.DataFrame(dict_distance).transpose().reset_index(drop=True)

    """
    diploid bowtie2_local vs bwa_mem 97.44
    diploid bowtie2_local vs hisat2_no_spliced 98.29
    diploid bwa_mem vs hisat2_no_spliced 98.29

    diploid_homozygous bowtie2_local vs bwa_mem 59.28
    diploid_homozygous bowtie2_local vs hisat2_no_spliced 66.84
    diploid_homozygous bwa_mem vs hisat2_no_spliced 71.04

    diploid_homozygous vs diploid bowtie2_local vs bowtie2_local 57.81 * 
    diploid_homozygous vs diploid bowtie2_local vs bwa_mem 51.13
    diploid_homozygous vs diploid bowtie2_local vs hisat2_no_spliced 57.65
    diploid_homozygous vs diploid bwa_mem vs bowtie2_local 56.77
    diploid_homozygous vs diploid bwa_mem vs bwa_mem 50.68 * 
    diploid_homozygous vs diploid bwa_mem vs hisat2_no_spliced 57.14
    diploid_homozygous vs diploid hisat2_no_spliced vs bowtie2_local 57.29
    diploid_homozygous vs diploid hisat2_no_spliced vs bwa_mem 51.13
    diploid_homozygous vs diploid hisat2_no_spliced vs hisat2_no_spliced 57.65 * 
    """

    # dfs
    print("plotting")
    df_plot = df_distances.pivot(index="idx1", columns="idx2", values="tree similarity").transpose()

    run_set_to_color = {"all":"black", "exclude_triploids":"gray"}
    mode_to_color = {"diploid":"black", "diploid_homozygous":"gray"}

    #row_colors_df = col_colors_df = trees_df.apply(lambda r: pd.Series({"idx":r.name, "run_set": run_set_to_color[r.run_set], "mode_tree": mode_to_color[r.mode_tree], "aligner": aligner_to_color[r.aligner]}), axis=1).set_index("idx")ç

    row_colors_df = col_colors_df = trees_df.apply(lambda r: pd.Series({"idx":r.name, "mode_tree": mode_to_color[r.mode_tree], "aligner": aligner_to_color[r.aligner]}), axis=1).set_index("idx")


    def get_round_or_nan(x):
        if pd.isna(x): return ""
        else: return str(int(round(x, 0)))

    annot_df = df_plot.map(get_round_or_nan)

    # plot
    #mask = df_plot.isnull()
    cm = sns.clustermap(df_plot, row_cluster=False, col_cluster=False, yticklabels=True, xticklabels=True, row_colors=row_colors_df, col_colors=row_colors_df, cbar_kws={'label': "min % equal nodes"}, cmap="rocket_r", linecolor="gray", linewidth=.01, square=True) # mask=mask

    tick_labels =  list(trees_df.apply(lambda r: ",".join(["%s=%s"%(f[0], r[f]) for f in ["run_set", "mode_tree", "aligner"]]), axis=1))

    cm.ax_heatmap.set_yticklabels([], rotation=0)
    cm.ax_heatmap.set_xticklabels([], rotation=90)
    cm.ax_heatmap.set_ylabel("")
    cm.ax_heatmap.set_xlabel("")


    adjust_cm_positions(cm, df_plot.copy(), hm_height_multiplier=0.03, hm_width_multiplier=0.03, cc_height_multiplier=0.03, rc_width_multiplier=0.03,
        distance_btw_boxes=0.01, cbar_width=0.02, cbar_height=0.15, ncol_fields=2, nrow_fields=2)

    # set colorbar
    hm_pos = cm.ax_heatmap.get_position()
    cm.ax_cbar.set_position([hm_pos.x0 + hm_pos.width + 0.02, hm_pos.y0 + hm_pos.height*0.125, hm_pos.width*0.1, hm_pos.height*0.75])


    # add legend
    def get_legend_element(facecolor, label, edgecolor="gray"): return mpatches.Patch(facecolor=facecolor, edgecolor=edgecolor, label=label)

    legend_elements = []
    #for f, color_dict in [("run_set", run_set_to_color), ("mode", mode_to_color), ("aligner", aligner_to_color)]:
    for f, color_dict in [("mode", mode_to_color), ("aligner", aligner_to_color)]:

        # header
        legend_elements.append(get_legend_element("white", f, edgecolor="white"))

        # each of the values
        for val, color in color_dict.items(): 
            legend_elements.append(get_legend_element(color, val, edgecolor="gray"))

    hm_pos = cm.ax_heatmap.get_position()
    cm.ax_row_colors.legend(handles=legend_elements, bbox_to_anchor=(0,0), loc="lower right") # upper right position of the legend is in 0, 1 (upper leftmost corner)

    # save
    save_fig(cm, "comparisons_trees.pdf", PlotsDir)



def get_cladeID_info_file_Bergin2022(df_metadata, trees_df, ProcessedDataDir):

    """generates the clade info for these trees"""

    cladeID_info_file = "%s/cladeID_Bergin2022_info.tab"%ProcessedDataDir
    if file_is_empty(cladeID_info_file):

        # add the trees_df with clade ID data
        trees_df = trees_df.copy().apply(get_trees_df_row_with_Bergin2022_clade_info, df_metadata = df_metadata[["runID", "clade_Bergin2022"]].copy(), axis=1)

        # define strains to discrd
        samples_discard_run1 = {"OSCAR01_L21-060", "OSCAR01_L21-067"}

        # tree with clades
        tree_series = trees_df[(trees_df.run_set=="exclude_triploids") & (trees_df.mode_tree=="diploid") & (trees_df.type_outgroup=="Cparapsilosis") & (trees_df.aligner=="bwa_mem")].iloc[0].copy()
        strains_tree = get_tree_ordered_clades(cp.deepcopy(tree_series.tree_corrected), pd.Series(tree_series.leaf_to_clade_Bergin2022))

        # get the metadata with trees
        df_metadata_trees = df_metadata.copy()
        leaf_to_clade_Bergin2022_adds = get_leaf_to_clade_Bergin2022_added_clades(tree_series.leaf_to_clade_Bergin2022, strains_tree, df_metadata.copy(), samples_discard_run1)
        df_metadata_trees["cladeID_based_on_Bergin2022"] = df_metadata_trees.runID.map(leaf_to_clade_Bergin2022_adds)

        # get clade
        df_cladeID = df_metadata_trees[(~pd.isna(df_metadata_trees.cladeID_based_on_Bergin2022))][["runID", "cladeID_based_on_Bergin2022"]].copy().rename(columns={"runID":"sampleID", "cladeID_based_on_Bergin2022":"cladeID"})
        df_cladeID["cladeID"] = df_cladeID.cladeID.apply(str)
        save_df_as_tab(df_cladeID, cladeID_info_file)

    return cladeID_info_file

def generate_SNPs_trees(df_reads_info, ProcessedDataDir, DataDir, list_cmds, sorted_aligners, min_coverage_pos, homozygous_Corthopsilosis_srrs, VarCallDir, genome_Cparapsilosis, PlotsDir):

    """Gets the cmds to generate SNP trees"""

    # keep
    df_reads_info = df_reads_info.copy()

    # define parms
    threads = 112
    threads_tree = 14 # took ~1h for each single trees
    mem_Gb = "20G" # only to run iqtree
    outdir_trees_all = "%s/SNP_trees"%(DataDir); make_folder(outdir_trees_all)
    min_support_trees = 90 # min branch support to define corrected trees

    ###### GENERATE TREES WITH C. ORTHO, TO DEFINE OUTGROUPS #######

    print("getting trees C. ortho outgroups")
    trees_data_df_file = "%s/trees_data_df_only_with_Cortho_outgroups.py"%ProcessedDataDir
    if file_is_empty(trees_data_df_file):

        # init a dict with dict info
        tress_data_dict = {}; Id=0

        # generate trees for different sample sets
        for run_set in ["all", "exclude_triploids"]:
            runIDs_ortho, runIDs_para = get_runIDs_one_run_set(df_reads_info, run_set)

            # for each mode, get a tree with the outgroup set for orthopsilosis
            for mode in ["diploid", "diploid_homozygous"]:
                for outgroup_runID in runIDs_ortho:
                    print(run_set, mode, outgroup_runID)

                    # get trees
                    tress_data_dict[Id] = get_cmds_tree_reconstruction_one_configuration(runIDs_para, run_set, mode, outgroup_runID, "bwa_mem", VarCallDir, outdir_trees_all, threads, mem_Gb, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, list_cmds, threads_tree, min_support_trees)
                    Id+=1
                    
        # debug if not finished
        if len(list_cmds)>0: return

        # get as df
        trees_data_df = pd.DataFrame(tress_data_dict).transpose().reset_index(drop=True)
        save_object(trees_data_df, trees_data_df_file)

    # do some plotting
    df_stats = load_object(trees_data_df_file)
    df_stats = df_stats[df_stats.run_set=="exclude_triploids"].copy()
    df_stats["n_leaves_para_outg"] = df_stats.Cpara_outgroup_leaves.apply(len)
    print(df_stats[["mode_tree", "pct_collapsed_nodes", "nsites_parsimony"]]) # Cpara_outgroup_leaves

    #for mode in sorted(set(df_stats.mode_tree)):
    #    df_m = df_stats[df_stats.mode_tree==mode].copy()
    #    print("C para outg leaves mode=%s"%mode)
    #    print("\n".join(df_m.apply(lambda r: "%s: %s"%(r.n_leaves_para_outg, ",".join(sorted(r.Cpara_outgroup_leaves))), axis=1)))

    ################################################################

    ####### RUN TREES WITH THE ORTHO-INFERRED OUTGROUPS ############

    trees_data_df_final_file = "%s/trees_data_df_all.py"%ProcessedDataDir
    if file_is_empty(trees_data_df_final_file):

        # load
        trees_data_df = load_object(trees_data_df_file)
        check_no_nans_in_df(trees_data_df)

        # plot the similarities between outgroups
        plot_similarities_outgroups_different_trees(trees_data_df, PlotsDir)
        # after looking at this I decide to use a different outgroup for each run_set and mode. However, I decide to discard run_set="all"

        # add as tuple
        trees_data_df["Cpara_outgroup_leaves_tuple"] = trees_data_df.Cpara_outgroup_leaves.apply(sorted).apply(tuple)

        # get trees with pre-defined outgroups
        print("getting final trees")
        for run_set in ["exclude_triploids", "all"]:

            # get the runs
            runIDs_ortho, runIDs_para = get_runIDs_one_run_set(df_reads_info, run_set)

            # add the sharing of the outgroup across all
            df_trees = trees_data_df[(trees_data_df.run_set==run_set) & (trees_data_df.type_outgroup=="Corthopsilosis")].copy()
            df_trees["Cpara_outgroup_fraction_shared_all"] = df_trees.Cpara_outgroup_leaves_tuple.map(df_trees.groupby("Cpara_outgroup_leaves_tuple").apply(len) / len(df_trees))
            check_no_nans_series(df_trees["Cpara_outgroup_fraction_shared_all"])

            for mode in ["diploid", "diploid_homozygous"]: # diploid_homozygous, diploid

                # add the sharing for this mode
                df_trees_m = df_trees[(df_trees.mode_tree==mode)].copy()
                df_trees_m["Cpara_outgroup_fraction_shared_mode"] = df_trees_m.Cpara_outgroup_leaves_tuple.map(df_trees_m.groupby("Cpara_outgroup_leaves_tuple").apply(len) / len(df_trees_m))
                check_no_nans_series(df_trees_m["Cpara_outgroup_fraction_shared_all"])

                # log
                fields_sharing = ["Cpara_outgroup_fraction_shared_mode", "Cpara_outgroup_fraction_shared_all", "Cpara_outgroup_leaves_tuple"]
                df_sharing = df_trees_m[fields_sharing].drop_duplicates().sort_values(by=fields_sharing, ascending=False)
                print("\n", run_set, mode, "\n", df_sharing[["Cpara_outgroup_fraction_shared_mode", "Cpara_outgroup_fraction_shared_all"]])

                # get the most shared outgroup
                best_outgroup_series = df_trees_m.sort_values(by=["Cpara_outgroup_fraction_shared_mode", "Cpara_outgroup_fraction_shared_all"], ascending=False).iloc[0]
                print("\n", run_set, mode, best_outgroup_series.Cpara_outgroup_leaves_tuple)

                # one job for each aligner
                for aligner in sorted_aligners:
                    print(run_set, mode, aligner)

                    # get trees
                    tress_data_dict = get_cmds_tree_reconstruction_one_configuration(runIDs_para, run_set, mode, ",".join(best_outgroup_series.Cpara_outgroup_leaves_tuple), aligner, VarCallDir, outdir_trees_all, threads, mem_Gb, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, list_cmds, threads_tree, min_support_trees)

                    if tress_data_dict is None: continue

                    # add to df
                    df_t = pd.DataFrame({0:tress_data_dict}).transpose()
                    trees_data_df = pd.concat([trees_data_df, df_t]).reset_index(drop=True)

        # debug if not finished
        if len(list_cmds)>0: return

        # save
        save_object(trees_data_df, trees_data_df_final_file)

    trees_data_df = load_object(trees_data_df_final_file)

    return trees_data_df

    ################################################################

def scatterplots_etest_MICs(df_etest, PlotsDir):

    """gets the scatters of the different comparsisons"""

    df_etest= df_etest.copy()

    tuples_fs = [("FLC_MIC", 'MIC ETEST 24H'), ("FLC_MIC", "MIC ETEST 48H"), ("MIC ETEST 24H", "MIC ETEST 48H"), ("FLC_MIC", 'MIC of recovered heteroresistant 24 h'), ("MIC ETEST 24H", 'MIC of recovered heteroresistant 24 h'), ("MIC ETEST 48H", 'MIC of recovered heteroresistant 48h')]

    for xfield, yfield in tuples_fs:

        fig = plt.figure(figsize=(3,3))
        log_x = "log %s"%xfield
        log_y = "log %s"%yfield

        df_etest[log_x] = df_etest[xfield].apply(np.log10)
        df_etest[log_y] = df_etest[yfield].apply(np.log10)

        ax = sns.scatterplot(df_etest, x=log_x, y=log_y, edgecolor="blue", facecolor="none", linewidth=.7, s=5)

        for val in get_bps_for_pheno_field("FLC_MIC"):

            ax.axhline(np.log10(val), ls="--", c="black", lw=.7)
            ax.axvline(np.log10(val), ls="--", c="black", lw=.7)

        ax.set_ylabel(log_y.replace("heteroresistant", "\nheteroresistant"))

        save_fig(fig, "eTest_MIC_correlations/%s_vs_%s.pdf"%(xfield, yfield), PlotsDir)

def plot_colonies_MICs_swarm(df_etest, PlotsDir):

    """plots differences in MICs between the colony-forming MICs"""

    df_etest = df_etest[df_etest["Colonies inside halo"]!="unknown"].copy()

    # ("R strains", df_etest[df_etest.FLC_MIC>=drug_to_MIC_breakpoints["FLC"][1]].copy())

    df_etest_8_32 = df_etest[(df_etest.FLC_MIC>=8) & (df_etest.FLC_MIC<=32)].copy()

    for type_df, df in [("all strains", df_etest.copy()), ("FLC MIC=8-32", df_etest_8_32)]:
        for field in ["FLC_MIC", "MIC ETEST 24H", "MIC ETEST 48H"]:

            log_field = "log %s"%field
            df[log_field] = df[field].apply(np.log10)

            df["resistance"] = df[field].apply(get_dicothomic_resistance, pheno=field)

            fig = plt.figure(figsize=(1.8,2.5))

            ax = sns.swarmplot(df, x="Colonies inside halo", y=log_field, hue="resistance", palette = {"R":"red", "I":"gray", "S":"blue"}, alpha=.9, s=4)
            ax.set_title("%s (%i)"%(type_df, len(df)), fontsize=10)
            ax.get_legend().remove()

            save_fig(fig, "Colonies_inside_halo_swarms/%s_%s.pdf"%(type_df, field), PlotsDir)


def generate_tree_one_phenotype_and_aligner(list_cmds, trees_data_df, tree_dir, aligner, target_runs, VarCallDir, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, replace=False):

    """Generates tree for one phenotype and aligner"""

    make_folder(tree_dir)
    tree_file_final = "%s/rooted_tree.nw"%tree_dir
    if file_is_empty(tree_file_final) or replace is True:

        # define the tree with all strains (with most adequate parameters)
        run_set = "exclude_triploids"
        mode_tree = "diploid"
        strains_tree_all = cp.deepcopy(trees_data_df[(trees_data_df.run_set==run_set) & (trees_data_df.mode_tree==mode_tree) & (trees_data_df.type_outgroup=="Cparapsilosis") & (trees_data_df.aligner==aligner)].tree_corrected.iloc[0])

        # define the outgroup
        strains_tree_all.prune(target_runs.intersection(set(strains_tree_all.get_leaf_names())))
        if len(strains_tree_all.get_children())!=2: raise ValueError("should be 2 children")        
        outgroup_node = sorted(strains_tree_all.get_children(), key=(lambda n: len(n.get_leaf_names())))[0]
        outgroup = ",".join(sorted(outgroup_node.get_leaf_names()))

        # get the tree and info
        threads = 112
        threads_tree = 14 # took ~1h for each single trees
        mem_Gb = "20G" # only to run iqtree
        min_support_trees = 90

        tree_data_dict = get_cmds_tree_reconstruction_one_configuration(sorted(target_runs), run_set, mode_tree, outgroup, aligner, VarCallDir, tree_dir, threads, mem_Gb, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, list_cmds, threads_tree, min_support_trees, skip_before_resample_running=False)

        # skip
        if tree_data_dict is None: 
            return

        # checks, if already generated, whether the outgroups are as expected
        if tree_data_dict["fraction_trees_expected_outgroup"]!=1: 
            print("WARNING: The outgroup is not as expected, setting another one")

            tree_obj = cp.deepcopy(tree_data_dict["tree"])
            chosen_outgroup_leaves = set(sorted([tuple(sorted(c.get_leaf_names())) for c in tree_obj.get_children()], key=len)[0])
            expected_outgroup_leaves = set(outgroup_node.get_leaf_names())
            print("extra leaves chosen:", chosen_outgroup_leaves.difference(expected_outgroup_leaves))
            print("extra leaves expected:", expected_outgroup_leaves.difference(chosen_outgroup_leaves))

        if tree_data_dict["fraction_trees_expected_outgroup"]!=1 and tree_data_dict["fraction_trees_expected_most_frequent_outgroup"]!=1:
            raise ValueError("the outgroup should be consistent across trees")

        # save
        rsync_file(tree_data_dict["tree_file"], tree_file_final)

    return tree_file_final


def get_parameters_from_json(json_file):

    """Takes a json file with the parameters and returns each important thing"""

    with open(json_file) as f: parameters_dict = json.load(f)

    # define things
    gridss_blacklisted_regions = parameters_dict["gridss_blacklisted_regions"]
    gridss_maxcoverage = parameters_dict["gridss_maxcoverage"]
    gridss_filters_dict_initial = parameters_dict["gridss_filters_dict"]
    max_rel_coverage_to_consider_del = parameters_dict["max_rel_coverage_to_consider_del"]
    min_rel_coverage_to_consider_dup = parameters_dict["min_rel_coverage_to_consider_dup"]

    # modify the filters
    def convert_tuple_to_list(x):
        if type(x)==list: return tuple(x)
        else: return x

    gridss_filters_dict = {k : convert_tuple_to_list(v) for k,v in gridss_filters_dict_initial.items()}

    return gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup


def get_SVcalling_parameters_perSVade_style(SVcalling_parameters):

    """Gets the SV calling parameters"""

    # if it is default
    if SVcalling_parameters=="default":

        gridss_filters_dict = default_filtersDict_gridss
        gridss_blacklisted_regions = default_gridss_blacklisted_regions
        gridss_maxcoverage = default_gridss_maxcoverage
        max_rel_coverage_to_consider_del = default_max_rel_coverage_to_consider_del
        min_rel_coverage_to_consider_dup = default_min_rel_coverage_to_consider_dup

    # if it is a json file
    elif not file_is_empty(SVcalling_parameters): gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup = get_parameters_from_json(SVcalling_parameters)

    else: raise ValueError("The calling parameters %s are not valid"%SVcalling_parameters)

    # check 
    if type(gridss_blacklisted_regions)!=str: raise ValueError("gridss_blacklisted_regions %s is not valid"%gridss_blacklisted_regions)
    if type(gridss_maxcoverage)!=int: raise ValueError("gridss_maxcoverage %s is not valid"%gridss_maxcoverage)
    if type(gridss_filters_dict)!=dict: raise ValueError("gridss_filters_dict %s is not valid"%gridss_filters_dict)
    if type(max_rel_coverage_to_consider_del)!=float: raise ValueError("max_rel_coverage_to_consider_del %s is not valid"%max_rel_coverage_to_consider_del)
    if type(min_rel_coverage_to_consider_dup)!=float: raise ValueError("min_rel_coverage_to_consider_dup %s is not valid"%min_rel_coverage_to_consider_dup)

    # return
    return gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup


def run_perSVade_parameter_optimization_several_samples(outdir_optmize_parameters_all, sorted_runs, ref_genome, runID_to_bam, perSVade_wrapper, simple_repeats_file, df_rep_samples_optimization, list_cmds, threads):

    """Rus parameter optimization for several samples, adding cmds """


    # get folder
    make_folder(outdir_optmize_parameters_all)
        
    # get cmds
    parms_dict = {}
    optimize_parms_dict = {}

    for runID in sorted_runs:

        # define outdir
        outdir_optmize_parms = "%s/%s"%(outdir_optmize_parameters_all, runID)
        if file_is_empty("%s/perSVade_finished_file.txt"%outdir_optmize_parms): 

            list_cmds.append(("optimizeParms", "%s optimize_parameters -o %s --fraction_available_mem 1.0 --threads %i --fractionRAM_to_dedicate 0.85 --verbose -r %s --min_chromosome_len 5000 -sbam %s -mchr %s --repeats_file %s --regions_SVsimulations random --simulation_ploidies diploid_hetero --nvars 50 --nsimulations 2 --range_filtering_benchmark theoretically_meaningful --keep_simulation_files"%(perSVade_wrapper, outdir_optmize_parms, threads, ref_genome, runID_to_bam[runID], mtDNA_chrom, simple_repeats_file), "perSVade_env", threads))
            continue
        
        # log accuracy
        df_acc = get_tab_as_df_or_empty_df("%s/optimized_parameters_accuracy.tab"%outdir_optmize_parms)
        min_f_value = min(df_acc[df_acc.SV_type=="integrated"].Fvalue)

        df_run = df_rep_samples_optimization[df_rep_samples_optimization.representative_sampleID==runID].copy()
        def get_sample_source(x):
            if x.startswith("OSCAR"): return x.split("_")[0]
            else: return "SRA"
        print("rep sample=%s"%runID, "min F=%.2f"%min_f_value, "%i samples"%(len(df_run)), ", ".join(["%s~(%.2f-%.2f)"%(f, min(df_run[f]), max(df_run[f])) for f in ["gDNA_median_coverage", "read_length", "median_insert_size", "median_insert_size_sd"]]), "run sources=%s"%(sorted(set(df_run.sampleID.apply(get_sample_source)))))

        # checks
        if min_f_value<0.85: raise ValueError("min F score in %s is %.3f"%(runID, min_f_value))

        # record data for integration
        parms_dict[runID] = {"parameterID": runID, "parameters_json": "%s/optimized_parameters.json"%(outdir_optmize_parms)}
        optimize_parms_dict[runID] = {"sampleID": runID, "optimize_parameters_outdir": outdir_optmize_parms}

    # debuug un run jobs
    if len(list_cmds)>0: return

    # run analyze_SV_parameteres to see if there are common ones
    df_SV_parameters = pd.DataFrame(parms_dict).transpose().reset_index(drop=True).sort_values(by=["parameterID"])
    df_optimize_parameters_outdirs = pd.DataFrame(optimize_parms_dict).transpose().reset_index(drop=True).sort_values(by=["sampleID"])

    SV_parameters_file = "%s/all_SV_parameters.tab"%outdir_optmize_parameters_all
    optimize_parameters_outdirs_file = "%s/all_optimize_parameters_outdirs.tab"%outdir_optmize_parameters_all
    save_df_as_tab(df_SV_parameters, SV_parameters_file)
    save_df_as_tab(df_optimize_parameters_outdirs, optimize_parameters_outdirs_file)

    outdir_analyze_SV_parameters = "%s/analyze_SV_parameters"%(outdir_optmize_parameters_all)
    if file_is_empty("%s/perSVade_finished_file.txt"%outdir_analyze_SV_parameters): 

        print("running analyze_SV_parameters")
        run_cmd("%s analyze_SV_parameters -o %s --fraction_available_mem 1.0 --threads %i --fractionRAM_to_dedicate 0.85 --verbose -r %s --min_chromosome_len 5000 -mchr %s --repeats_file %s --SV_parameters %s --optimize_parameters_outdirs %s"%(perSVade_wrapper, outdir_analyze_SV_parameters, threads, ref_genome, mtDNA_chrom, simple_repeats_file, SV_parameters_file, optimize_parameters_outdirs_file), env="perSVade_env")

def get_cmds_perSVade_SVs_and_CNVs_calls(VarCallDir, sorted_aligners, simple_repeats_file, sorted_runs, genome_Cparapsilosis, ProcessedDataDir, list_cmds, DataDir, gff_Cparapsilosis, SV_parms_file="auto", integrated_SVs_CNVs_dir="auto"):

    """Runs SV/CNV calling several parameters. Note that for SV calling it is ony possible with bwa_mem."""

    # define general things
    runID_to_aligner_to_outdir = {runID : {aligner : "%s/%s/aligner-%s"%(VarCallDir, runID, aligner) for aligner in sorted_aligners} for runID in sorted_runs}
    runID_to_aligner_to_bam = {runID : {aligner : "%s/align_reads/aligned_reads.bam.sorted"%(runID_to_aligner_to_outdir[runID][aligner]) for aligner in sorted_aligners} for runID in sorted_runs}

    ##### RUN THE PARAMETER OPT #####

    # only if needed
    if SV_parms_file=="auto":

        # run module get_stats_optimization to define the samples in which to run parameter optimization (setting --min_chromosome_len 100000 to skip mtDNA)
        outdir_parameter_optimization = "%s/SV_parameter_optimization"%DataDir; make_folder(outdir_parameter_optimization)
        file_sorted_bams = "%s/bamfiles.tab"%outdir_parameter_optimization
        df_bams = pd.DataFrame({runID : {"sorted_bam":runID_to_aligner_to_bam[runID]["bwa_mem"], "sampleID":runID}  for runID in sorted_runs}).transpose().reset_index(drop=True)[["sampleID", "sorted_bam"]]
        if any(df_bams.sorted_bam.apply(file_is_empty)): raise ValueError("empty bams")
        save_df_as_tab(df_bams, file_sorted_bams)

        outdir_select_parms = "%s/stats_optimization"%outdir_parameter_optimization
        if file_is_empty("%s/perSVade_finished_file.txt"%outdir_select_parms):
            run_cmd("%s get_stats_optimization -o %s  --verbose --fraction_available_mem 1.0 --fractionRAM_to_dedicate 0.8 -thr 112 -r %s -mchr %s --min_chromosome_len 100000 --samples_file %s --overlap_coverage 25 --overlap_insert_size 25 --overlap_insert_size_sd 25 --overlap_read_len 25 --rerun_sample_clustering --skip_cleaning_tmpdir --overlap_coverage_mtDNA 100"%(perSVade_wrapper, outdir_select_parms, genome_Cparapsilosis, mtDNA_chrom, file_sorted_bams), env="perSVade_env") # --min_Gb_per_core 4
        
        # clean
        delete_folder("%s/tmp_bams"%outdir_select_parms)

        # analyze representative samples
        df_rep_samples_optimization = get_tab_as_df_or_empty_df("%s/samples_parameter_optimization.tab"%outdir_select_parms)
        representative_runs = sorted(set(df_rep_samples_optimization.representative_sampleID))
        print("There are %i/%i runs with unique parameters. Running parameter optimization on them"%(len(representative_runs), len(df_rep_samples_optimization)))

        # run parameter optimization for each representative sample
        outdir_optmize_parameters_all = "%s/optimize_parameters_representative_runs"%outdir_parameter_optimization
        runID_to_bam_bwa_mem = {runID : aln_to_bam["bwa_mem"] for runID, aln_to_bam in runID_to_aligner_to_bam.items()}
        run_perSVade_parameter_optimization_several_samples(outdir_optmize_parameters_all, representative_runs, genome_Cparapsilosis, runID_to_bam_bwa_mem, perSVade_wrapper, simple_repeats_file, df_rep_samples_optimization, list_cmds, 112)

        # debug un run jobs
        if len(list_cmds)>0: return

        # In run_perSVade_parameter_optimization_several_samples, I found that taking the most conservative parameters was a good course, as they all work the same for all 
        SV_parms_file = "%s/analyze_SV_parameters/more_conservative_parameters_OSCAR02_L21-190.json"%outdir_optmize_parameters_all


    # check
    if file_is_empty(SV_parms_file): raise ValueError(" invalid file")

    #################################

    ######## RUN SV AND CNV CALLING ########

    # for each aligner
    if integrated_SVs_CNVs_dir=="auto":
        integrated_SVs_CNVs_dir = "%s/integrated_SVs_CNVs"%DataDir
    make_folder(integrated_SVs_CNVs_dir)

    for aligner in sorted_aligners:

        # run for each run
        for runID in  sorted_runs:

            # get the dirs
            outdir_all = runID_to_aligner_to_outdir[runID][aligner]
            sorted_bam = runID_to_aligner_to_bam[runID][aligner]
            if file_is_empty(sorted_bam): raise ValueError("file should exist")

            # define the parameters
            outdir_SV_CNV_calls = "%s/SV_CNV_calling"%outdir_all
            make_folder(outdir_SV_CNV_calls)
            if file_is_empty("%s/perSVade_finished_file.txt"%outdir_SV_CNV_calls):

                # define modules to run
                if aligner=="bwa_mem": # call SVs only possible for bwa mem
                    modules_run = "call_CNVs,call_SVs,integrate_SV_CNV_calls"
                    extra_args = "--SVcalling_parameters %s"%SV_parms_file
                else:
                    modules_run = "call_CNVs,integrate_SV_CNV_calls"
                    extra_args = ""

                # get cmd
                map_per_pos_file = "%s/mappability_outdir/genome.fasta.mappability_per_position.bed"%ProcessedDataDir
                cmd_run = "%s run_several_modules %s %s -r %s -sbam %s --repeats_file %s -p 2 --min_chromosome_len 100000 --verbose --fractionRAM_to_dedicate 0.85 --mitochondrial_chromosome %s --cnv_calling_algs HMMcopy,AneuFinder --window_size_CNVcalling 300 %s --skip_coverage_correction --mappability_file %s --average_cov_measure median --max_fraction_N_bases 0.1 --max_fraction_repeats 0.1 --min_median_mappability 0.75 --fraction_available_mem 0.25 --threads 28"%(perSVade_wrapper, modules_run, outdir_SV_CNV_calls, genome_Cparapsilosis, sorted_bam, simple_repeats_file, mtDNA_chrom, extra_args, map_per_pos_file)

                list_cmds.append(("SV_CNV_calls", cmd_run, "perSVade_env", 28))

        # debug jobs to run
        if any(map(lambda x: x[0]=="SV_CNV_calls", list_cmds)): continue

        # integrate variants using perSVade's module integrate_several_samples
        print("integrating variants...")

        # prepare input file with paths
        variant_integration_dir = "%s/SV_CNV_%s"%(integrated_SVs_CNVs_dir, aligner);
        make_folder(variant_integration_dir)
        variant_paths = "%s/paths_variant_calls.tab"%variant_integration_dir
        variant_paths_df = pd.DataFrame({runID : {"sampleID":runID, "integrate_SV_CNV_calls_outdir": "%s/SV_CNV_calling/integrate_SV_CNV_calls"%runID_to_aligner_to_outdir[runID][aligner]} for runID in sorted_runs}).transpose().reset_index(drop=True)
        variant_paths_df["sorted_bam"] = variant_paths_df.sampleID.apply(lambda s:  runID_to_aligner_to_bam[s][aligner])
        save_df_as_tab(variant_paths_df.sort_values(by=["sampleID"]), variant_paths)

        # write file with fields to keep
        fields_SV_CNVs = ['#CHROM', 'POS', 'ID', 'INFO_END', 'ALT', 'INFO_merged_relative_CN', 'INFO_RELCOVERAGE', 'INFO_SVTYPE',  'INFO_median_relative_CN_AneuFinder', 'INFO_BREAKENDIDs', 'INFO_misc', 'INFO_median_relative_CN_HMMcopy', 'REF', 'INFO_RELCOVERAGE_TO_3', 'INFO_RELCOVERAGE_TO_5', 'INFO_variantID',  'INFO_median_coverage_corrected']

        if aligner=="bwa_mem":
            fields_SV_CNVs += ['INFO_QUAL_min', 'INFO_BREAKEND_length_inexactHomology', 'INFO_QUAL', 'INFO_length_inexactHomology', 'INFO_length_event_min', 'INFO_length_inexactHomology_min', 'INFO_length_inexactHomology_mean', 'INFO_BREAKPOINTID', 'INFO_real_AF_mean', 'INFO_FILTER', 'INFO_length_microHomology_min', 'INFO_BREAKEND_length_microHomology', 'INFO_overlaps_repeats', 'INFO_any_has_poly16GC', 'INFO_length_event_max', 'INFO_length_microHomology', 'INFO_real_AF_max', 'INFO_QUAL_mean', 'INFO_BREAKEND_FILTER', 'INFO_has_poly16GC', 'INFO_BREAKEND_real_AF', 'INFO_any_overlaps_repeats', 'INFO_BREAKPOINTIDs', 'INFO_BREAKEND_coordinates', 'INFO_length_inexactHomology_max', 'INFO_best_FILTER', 'INFO_length_event_mean', 'INFO_real_AF', 'INFO_length_microHomology_max', 'INFO_QUAL_max', 'INFO_BREAKEND_has_poly16GC', 'INFO_len_inserted_sequence_mean', 'INFO_len_inserted_sequence_max', 'INFO_bpIDs', 'INFO_len_inserted_sequence_min', 'INFO_real_AF_min', 'INFO_worse_FILTER', 'INFO_length_microHomology_mean', 'INFO_BREAKEND_len_inserted_sequence', 'INFO_all_FILTERs', 'INFO_BPS_TYPE', 'INFO_BREAKEND_QUAL', 'INFO_BREAKEND_overlaps_repeats']

        if len(fields_SV_CNVs)!=len(set(fields_SV_CNVs)): raise ValueError("should be unique")
        fields_SV_CNVs_file = "%s/fields_SV_CNVs.txt"%variant_integration_dir
        open(fields_SV_CNVs_file, "w").write("\n".join(fields_SV_CNVs))

        # get cmd of integration
        outdir_integration = "%s/variant_integration"%variant_integration_dir
        if file_is_empty("%s/perSVade_finished_file.txt"%outdir_integration): 

            cmd_int = "%s integrate_several_samples --fraction_available_mem 1.0 --fractionRAM_to_dedicate 0.8 --min_chromosome_len 100000 --threads 112 --verbose --paths_table %s  -o %s --fields_SV_CNVs %s -r %s -mchr %s --repeats_file %s -p 2 --tol_bp 50 --pct_overlap 75 --CNV_overlap_only_based_on_pct"%(perSVade_wrapper, variant_paths, outdir_integration, fields_SV_CNVs_file, genome_Cparapsilosis, mtDNA_chrom, simple_repeats_file)
            list_cmds.append(("integrate_SV_CNV", cmd_int, "perSVade_env", 112))
            continue

        # annotate the variants
        outdir_SV_CNV_annot = "%s/SV_CNV_annotation"%variant_integration_dir
        if file_is_empty("%s/perSVade_finished_file.txt"%outdir_SV_CNV_annot):

            run_cmd("%s annotate_SVs -o %s -r %s -mchr %s --SV_CNV_vcf  %s -gff %s --mitochondrial_code %i --gDNA_code %i --min_chromosome_len 100000 --verbose --fraction_available_mem 0.1 --threads 112 --fractionRAM_to_dedicate 0.85"%(perSVade_wrapper, outdir_SV_CNV_annot, genome_Cparapsilosis, mtDNA_chrom, "%s/integrated_SVs_CNVs_annotation_fields.vcf"%outdir_integration, gff_Cparapsilosis, Cpara_mtDNA_code, Cpara_gDNA_code), env="perSVade_env")

    ########################################

def integrate_and_annotate_small_vars(VarCallDir, sorted_aligners, sorted_runs, genome_Cparapsilosis, ProcessedDataDir, list_cmds, DataDir, gff_Cparapsilosis, simple_repeats_file, integrated_small_vars_dir="auto"):

    """gets cmds to integrate and annotate small variants"""

    print("integrating and anotating small vars")

    # define general vars
    if integrated_small_vars_dir=="auto":
        integrated_small_vars_dir = "%s/integrated_small_vars"%DataDir
    
    make_folder(integrated_small_vars_dir)
    runID_to_aligner_to_outdir = {runID : {aligner : "%s/%s/aligner-%s"%(VarCallDir, runID, aligner) for aligner in sorted_aligners} for runID in sorted_runs}
    runID_to_aligner_to_bam = {runID : {aligner : "%s/align_reads/aligned_reads.bam.sorted"%(runID_to_aligner_to_outdir[runID][aligner]) for aligner in sorted_aligners} for runID in sorted_runs}

    # go through each aligner
    for aligner in sorted_aligners:

        # get the variant paths file
        variant_integration_dir = "%s/small_vars_%s"%(integrated_small_vars_dir, aligner)
        make_folder(variant_integration_dir)
        variant_paths = "%s/paths_variant_calls.tab"%variant_integration_dir
        variant_paths_df = pd.DataFrame({runID : {"sampleID":runID, "call_small_variants_p2_outdir": "%s/small_vars"%runID_to_aligner_to_outdir[runID][aligner]} for runID in sorted_runs}).transpose().reset_index(drop=True)
        variant_paths_df["sorted_bam"] = variant_paths_df.sampleID.apply(lambda s:  runID_to_aligner_to_bam[s][aligner])
        save_df_as_tab(variant_paths_df.sort_values(by=["sampleID"]), variant_paths)

        # define the fields to keep
        fields_small_variants = ["#Uploaded_variation", "#CHROM", "POS", "REF", "ALT", "QUAL", "CALLEDALGS", "HC_GT", "ISSNP", "NCALLED", "NPASS", "PASSALGS", "bt_GT", "common_GT", "fb_GT", "mean_AD", "mean_DP", "mean_fractionReadsCov_PASS_algs", "mean_fractionReadsCov_called_algs"]
        fields_small_variants_file = "%s/fields_small_variants.txt"%variant_integration_dir
        open(fields_small_variants_file, "w").write("\n".join(fields_small_variants))

        # get cmd of integration
        outdir_integration = "%s/variant_integration"%variant_integration_dir
        if file_is_empty("%s/perSVade_finished_file.txt"%outdir_integration): 

            cmd_int = "%s integrate_several_samples --fraction_available_mem 1.0 --fractionRAM_to_dedicate 0.8 --min_chromosome_len 100000 --threads 112 --verbose --paths_table %s  -o %s --fields_small_variants %s -r %s -mchr %s --repeats_file %s -p 2"%(perSVade_wrapper, variant_paths, outdir_integration, fields_small_variants_file, genome_Cparapsilosis, mtDNA_chrom, simple_repeats_file)
            list_cmds.append(("integrate_small_vars", cmd_int, "perSVade_env", 112))
            continue

        # annotate the variants
        print("annotating")
        outdir_small_vars_annot = "%s/small_vars_annot"%variant_integration_dir
        #delete_folder(outdir_small_vars_annot)
        if file_is_empty("%s/perSVade_finished_file.txt"%outdir_small_vars_annot):

            cmd_annot = "%s annotate_small_vars -o %s -r %s -mchr %s --merged_vcf  %s -gff %s --mitochondrial_code %i --gDNA_code %i --min_chromosome_len 100000 --verbose --fraction_available_mem 1.0 --threads 112 --fractionRAM_to_dedicate 0.85"%(perSVade_wrapper, outdir_small_vars_annot, genome_Cparapsilosis, mtDNA_chrom, "%s/integrated_small_variants_annotation_fields.vcf"%outdir_integration, gff_Cparapsilosis, Cpara_mtDNA_code, Cpara_gDNA_code)
            list_cmds.append(("annotate_small_vars", cmd_annot, "perSVade_env", 112))


def get_dfs_stats_pairwise_SNP_distances_per_window_one_run(aligner, wsize_blocks, target_runID, df_distances_per_window_file):

    """loads combined df and also adds the general stats"""
    print(aligner, wsize_blocks, target_runID)

    # get stats per window df
    df_stats_all = load_object(df_distances_per_window_file)
    df_stats_all["aligner"] = aligner
    df_stats_all["window_size"] = wsize_blocks
    df_stats_all["target_runID"] = target_runID
    df_stats_all["diff_vars_per_kb"] = (df_stats_all.nvars_different / df_stats_all.len_window)*1000

    # discard some zygosities
    df_stats_all = df_stats_all[df_stats_all.zygosity!="unknown_zygosity"].copy()

    # get integrated measures
    def get_series_stats_integrated_one_comparison(df):

        # checks
        if len(df)!=len(set(df.ID_window)): raise ValueError("wID should be unique")

        # init stats
        stats_dict  = {f : df.iloc[0][f] for f in ["origin_runID", "type_var", "zygosity"]}
        stats_dict["target_runID"] = target_runID
        stats_dict["window_size"] = wsize_blocks
        stats_dict["aligner"] = aligner

        # add general metrics
        stats_dict["whole_genome_diff_vars_per_kb"] = (sum(df.nvars_different) / sum(df.len_window))*1000

        for fn_name, fn_obj in [("std", np.std), ("mean", np.mean), ("mad", scipy.stats.median_abs_deviation), ("median", np.median)]:
            stats_dict["%s_diff_vars_per_kb"%fn_name] = fn_obj(df.diff_vars_per_kb)

        # add fraction of genome with a diff_vars_per_kb below certain thresholds
        for vars_kb_tshd in [0.001, 0.005, 0.01, 0.05, 0.1, 0.15]:
            df_low = df[df.diff_vars_per_kb<vars_kb_tshd]
            stats_dict["fraction_diff_vars_per_kb_below_%s"%vars_kb_tshd] = sum(df_low.len_window) / sum(df.len_window)

        # return
        return pd.Series(stats_dict)

    df_stats = df_stats_all.groupby(["origin_runID", "type_var", "zygosity"]).apply(get_series_stats_integrated_one_comparison).reset_index(drop=True)

    # keep only some fields
    #fields = ['type_var', 'zygosity', 'origin_runID', 'target_runID', 'chromosome', 'start', 'end', 'nvars_different', 'aligner', 'window_size', 'diff_vars_per_kb']
    # 'nvars_total', 'nvars_shared'
    #df_stats_all = df_stats_all[fields].copy()
    del df_stats_all

    return df_stats


def get_stats_pairwise_SNP_distances_per_window(aligner_to_df_small_vars_filt_file, aligner_to_wsize_to_df_coverage_per_window_file, window_sizes_list, ProcessedDataDir, CurDir, target_runs_distances, genome_Cparapsilosis, simple_repeats_file, list_cmds):

    """generates df with the distances stats"""

    # define final file
    #df_stats_pairwise_SNP_distances_per_window_file = "%s/df_stats_pairwise_SNP_distances_per_window.py"%ProcessedDataDir
    df_stats_pairwise_SNP_distances_per_window_integrated_stats_file = "%s/df_stats_pairwise_SNP_distances_per_window_integrated_stats.py"%ProcessedDataDir
    if file_is_empty(df_stats_pairwise_SNP_distances_per_window_integrated_stats_file):

        # init inputs to integrate
        print("getting jobs distances")
        inputs_integrate = []

        # get cmds
        for aligner, df_small_vars_filt_file in aligner_to_df_small_vars_filt_file.items():
            for wsize_blocks in window_sizes_list:

                # define chroms to exclude, depending on the wsize
                if wsize_blocks<10000:
                    excluded_chroms_blocks = set(get_chr_to_len(genome_Cparapsilosis)).difference({mtDNA_chrom})

                else:
                    excluded_chroms_blocks = {mtDNA_chrom}

                # get densities file
                df_var_densities_blocks_file = generate_variant_per_window_densities_df(df_small_vars_filt_file, ProcessedDataDir, genome_Cparapsilosis, aligner, wsize_blocks, target_runs_distances, simple_repeats_file, aligner_to_wsize_to_df_coverage_per_window_file[aligner][wsize_blocks], excluded_chroms=excluded_chroms_blocks)

                # add the cmds to calculate stats
                distances_dir = "%s/dfs_distances_per_window"%ProcessedDataDir
                make_folder(distances_dir)
                
                for target_runID in target_runs_distances:
                    df_distances_per_window_file = "%s/df_distances_vs_%s_%s_%ibp.py"%(distances_dir, target_runID, aligner, wsize_blocks)
                    if file_is_empty(df_distances_per_window_file):
                        list_cmds.append(("dist_per_window", "%s/get_genetic_distances_one_sample.py %s %s 28 %s"%(CurDir, df_var_densities_blocks_file, df_distances_per_window_file, target_runID), "Cparapsilosis_popgen_env", 28))

                    # keep if they are unique
                    inputs_integrate.append((aligner, wsize_blocks, target_runID, df_distances_per_window_file))

        # debug
        if len(list_cmds)>0: return

        # integrate all stats
        print("integrate...")
        with multiproc.Pool(112) as pool:
            df_stats_pairwise_SNP_distances_per_window_integrated_stats = pd.concat(pool.starmap(get_dfs_stats_pairwise_SNP_distances_per_window_one_run, inputs_integrate, chunksize=1))
            pool.close()
            pool.terminate()

        # define run combinations
        print("filtering dfs")
        unique_run_combinations = set(itertools.combinations(target_runs_distances, 2))

        for df in [df_stats_pairwise_SNP_distances_per_window_integrated_stats]:
            df["tuple_samples"] = df[["origin_runID", "target_runID"]].apply(tuple, axis=1)
            missing_tuples = unique_run_combinations.difference(set(df["tuple_samples"]))
            if len(missing_tuples)>0: raise ValueError("There are missing tuples, %i"%len(missing_tuples))

        df_stats_pairwise_SNP_distances_per_window_integrated_stats = df_stats_pairwise_SNP_distances_per_window_integrated_stats[df_stats_pairwise_SNP_distances_per_window_integrated_stats.tuple_samples.isin(unique_run_combinations)].copy()

        # drop
        print("save")
        df_stats_pairwise_SNP_distances_per_window_integrated_stats = df_stats_pairwise_SNP_distances_per_window_integrated_stats.drop(["tuple_samples"], axis=1)

        # clean
        delete_folder(distances_dir)

        # write one df for different data types
        print("save diff data types")
        for aligner in aligner_to_df_small_vars_filt_file.keys():
            for window_size in window_sizes_list:
                for zygosity in sorted(set(df_stats_pairwise_SNP_distances_per_window_integrated_stats.zygosity)):
                    for type_var in sorted(set(df_stats_pairwise_SNP_distances_per_window_integrated_stats.type_var)):

                        df_save = df_stats_pairwise_SNP_distances_per_window_integrated_stats.copy()
                        df_save = df_save[(df_save.aligner==aligner) & (df_save.window_size==window_size) & (df_save.zygosity==zygosity) & (df_save.type_var==type_var)].copy()

                        dir_dfs = "%s/dfs_stats_pairwise_SNP_distances_per_window"%ProcessedDataDir
                        make_folder(dir_dfs)
                        save_object(df_save, "%s/df_%s_%ibp_%s_%s.py"%(dir_dfs, aligner, window_size, zygosity, type_var))

        # write general df
        print("save general")
        save_object(df_stats_pairwise_SNP_distances_per_window_integrated_stats, df_stats_pairwise_SNP_distances_per_window_integrated_stats_file)


def get_df_stats_distances_one_pair_samples_and_types_vars(type_var, zygosity, origin_runID, target_runID, df_vars):

    """Gets, for one pair of runs and stats, the numbers of vars that are different or similar btw them"""

    # get dfs
    df_origin = df_vars[df_vars.runID==origin_runID].copy().sort_values(by=["chromosome", "start", "end"])
    df_target = df_vars[df_vars.runID==target_runID].copy().sort_values(by=["chromosome", "start", "end"])

    # checks
    if list(df_origin.ID_window)!=list(df_target.ID_window): raise ValueError("the windows should be the same")

    # merge data
    fields = ["chromosome", "start", "end", "ID_window", "len_window", "set_numeric_vars"]
    df_merged = df_origin[fields].merge(df_target[fields], on=["chromosome", "start", "end", "ID_window", "len_window"], suffixes=("_origin", "_target"), validate="one_to_one", how="left")
    check_no_nans_in_df(df_merged)

    # add stats
    fields_vars = ["set_numeric_vars_origin", "set_numeric_vars_target"]
    df_merged["nvars_total"] = df_merged[fields_vars].apply(lambda r: set.union(*r), axis=1).apply(len)
    df_merged["nvars_shared"] = df_merged[fields_vars].apply(lambda r: set.intersection(*r), axis=1).apply(len)
    df_merged["nvars_unique_in_target"] = df_merged.apply(lambda r: r.set_numeric_vars_target.difference(r.set_numeric_vars_origin), axis=1).apply(len)
    df_merged["nvars_unique_in_origin"] = df_merged.apply(lambda r: r.set_numeric_vars_origin.difference(r.set_numeric_vars_target), axis=1).apply(len)
    df_merged["vars_jaccard_distance"] = 1 - (df_merged.nvars_shared / df_merged.nvars_total)
    df_merged["nvars_different"] = df_merged.nvars_unique_in_target + df_merged.nvars_unique_in_origin

    df_merged["type_var"] = type_var 
    df_merged["zygosity"] = zygosity 
    df_merged["origin_runID"] = origin_runID
    df_merged["target_runID"] = target_runID


    # checks
    if origin_runID==target_runID and not all(df_merged.nvars_different==0):
        raise ValueError("not valid")

    return df_merged[["type_var", "zygosity", "origin_runID", "target_runID", "chromosome", "start", "end", "ID_window", "len_window", "nvars_total", "nvars_shared", "nvars_unique_in_target", "nvars_unique_in_origin", "nvars_different"]]



def check_that_tree_has_no_politomies(tree):

    """Raises an error if tree has politomies"""

    for n in tree.traverse():
        if n.is_leaf() is False and len(n.get_children())!=2: raise ValueError("node %s has !=2 children: %i"%(n.name, len(n.get_children())))

def get_annotationColor_on_bgcolor(bgcolor, threshold_gray=0.4):

    """Takes a bgcolor and retrives the optimum color to write an annotation on. """

    # change to RGB tuple
    if type(bgcolor)==str and not bgcolor.startswith("#"): bgcolor = colors.to_rgb(bgcolor)

    # debug
    elif type(bgcolor)==str and bgcolor.startswith("#"): bgcolor = tuple([x/255 for x in ImageColor.getcolor(bgcolor, "RGB")])

    # get the gray (a number between 0(black) and 1(white))
    gray = np.mean(bgcolor, -1)


    if gray>1: raise ValueError("gray can't be above 1. bgcolor %s is strange"%bgcolor)

    if gray<threshold_gray: return "white"
    else: return "black"


def plot_tree_AF_resistance_samples_for_GWAS(tree, sample_to_resistance, outdir_drug, min_support, title_prefix):

    """Makes a tree with the representative samples"""

    #print("getting tree with samples")

    # prune the tree
    tree = cp.deepcopy(tree)
    check_that_tree_has_no_politomies(tree)

    # get the resistance to color dict
    resistance_type_to_color = {1:"red", 0:"blue"}

    # add the resistance to each node
    for n in tree.traverse():

        if n.is_leaf(): n.resistance = sample_to_resistance[n.name]
        else:

            resistance_leaves = {sample_to_resistance[l.name] for l in n.get_leaves()}
            if len(resistance_leaves)==1: n.resistance = next(iter(resistance_leaves))
            else: n.resistance = 0


    # init counters of resistance transitions
    n_transitions = 0
    n_transitions_with_correct_support = 0

    # go through each node
    for n in tree.traverse():

        # add whether the node is a transition node (phenotype is different than the parent)
        if n.is_root() is False and n.resistance!=(n.get_ancestors()[0].resistance): n.transition = True
        else: n.transition = False

        # add counters of transitions
        if n.transition is True:
            n_transitions += 1  
            if (n.support>=min_support or n.is_leaf()) and (n.get_ancestors()[0].support>=min_support): n_transitions_with_correct_support += 1

        # node
        nst = NodeStyle()
        nst["hz_line_width"] = 8
        nst["vt_line_width"] = 8
        nst["size"] = 0

        if n.support>=min_support or n.is_leaf(): 
            color = "black"
            label_text_color = "white"
        else: 
            color = "gray"
            label_text_color = "black"

        nst["hz_line_color"] = color
        nst["vt_line_color"] = color
        n.set_style(nst)

        # add circle
        if n.is_leaf() is False: n.add_face(CircleFace(10, color, style='sphere', label=dict(text=str(n.transition)[0], font="arial", color=label_text_color, fontsize=12)), column={True:1, False:0}[n.is_leaf()])

        # add the resistance
        if n.is_leaf(): 

            # add a box with the resistance
            color_leaf =  resistance_type_to_color[sample_to_resistance[n.name]]
            label = {0:"0", 1:"1"}[sample_to_resistance[n.name]]
            square_width = 20
            n.add_face(RectFace(square_width, square_width, fgcolor="black", bgcolor=color_leaf, label={"text":label, "color":get_annotationColor_on_bgcolor(color_leaf)}), position="aligned", column=0)

            n.add_face(CircleFace(10, color_leaf, style='sphere', label=dict(text=str(n.transition)[0], font="arial", color='black', fontsize=12)), column={True:1, False:0}[n.is_leaf()])

    #print("%s-%s. There are %i/%i transitions with the current and parent nodes with supp>=%i"%(species, drug, n_transitions_with_correct_support, n_transitions, min_support))

    # init treestyle
    ts = TreeStyle()
    ts.title.add_face(TextFace("%s GWAS samples. Black nodes have support>=%i\n"%(title_prefix, min_support), bold=True, fsize=12), column=0)
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_leaf_name = True

    #import matplotlib
    #matplotlib.use('Agg')
    #os.system("export QT_QPA_PLATFORM=offscreen")
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    tree.render(file_name="%s/tree_phenotypes_min_support=%s_%s.pdf"%(outdir_drug, min_support, title_prefix), tree_style=ts)#, dpi=500,  units='mm', h=nleafs*2) #w=20


def get_int_or_zero(x):
    if str(x)=="nan": return 0
    else: return int(x)

def report_politomy_types(tree_file, min_support):

    """reports the politomy types for one tree"""

    tree =  get_correct_tree(tree_file, min_support=min_support)

    # get df
    politomyI_to_data = {}; Ip = 0
    for n in tree.traverse():
        if len(n.get_children())>2:
            politomyI_to_data[Ip] = {}

            all_nleaves = {len(c.get_leaf_names()) for c in n.get_children()}
            for nl in all_nleaves:
                politomyI_to_data[Ip]["nleaves_%i"%nl] = sum([len(c.get_leaf_names())==nl for c in n.get_children()])

            Ip+=1

    df_politomies = pd.DataFrame(politomyI_to_data).transpose().map(get_int_or_zero)

    # print
    print("politomies # nodes", ", ".join(["%s=%i"%(k, sum(df_politomies[k])) for k in sorted(df_politomies.keys(), key=(lambda x: int(x.split("_")[1])))]))

def get_balanced_samples_for_GWAS(phenotypes_file, treefile, ancestral_GWAS_dir, min_support, ASR_method, outdir, list_cmds, nrepresentative_samples=3, extra_args=""):

    """Gets the balanced set of samples to do GWAS"""

    # run
    if file_is_empty("%s/pipeline_finished.txt"%outdir):

        cmd_sample_selection = "%s/get_balanced_samples.py --phenotypes_all %s --outdir %s --tree_file %s --threads 112 --min_support %i --ASR_method %s --nrepresentative_samples %i %s"%(ancestral_GWAS_dir, phenotypes_file, outdir, treefile, min_support, ASR_method, nrepresentative_samples, extra_args)
        list_cmds.append(("GWAS_sel_samples", cmd_sample_selection, "ancestral_GWAS_env", 112))

        return (None, None)

    # get files
    phenotypes_file =  "%s/prunned_phenotypes.tab"%outdir
    df_prunning = get_tab_as_df_or_empty_df("%s/prunning_info.tab"%outdir)

    # checks
    if len(df_prunning)>1:

        # check that the refinement has increased the # transitions
        if  df_prunning[df_prunning.prunning_step==1].n_transitions.iloc[0]>df_prunning[df_prunning.prunning_step==max(df_prunning.prunning_step)].n_transitions.iloc[0]:
            print("WARNING: the prunning should not decrease the # transitions. Check %s"%outdir)

    # get
    return phenotypes_file, df_prunning

def get_tab_as_df_or_empty_df(file):

    """Gets df from file or empty df"""

    nlines = len([l for l in open(file, "r").readlines() if len(l)>1])

    if nlines==0: return pd.DataFrame()
    else: return pd.read_csv(file, sep="\t")

def get_balanced_samples_for_GWAS_redundant_samples(redundant_samples, tree_initial, outdir_gwas, df_phenotypes_all, ancestral_GWAS_dir, general_min_support, general_ASR_method, list_cmds):

    """runs the sample selection for redudnant samples"""

    # keep
    df_phenotypes_all = df_phenotypes_all.copy()
    tree_initial = cp.deepcopy(tree_initial)

    # if there are some R samples
    if len(redundant_samples)>1:
        
        # write tree
        treefile_redundant_samples = "%s/tree_redundant_samples.nw"%outdir_gwas
        tree_redundant_samples = cp.deepcopy(tree_initial); tree_redundant_samples.prune(redundant_samples, preserve_branch_length=True)
        check_that_tree_has_no_politomies(tree_redundant_samples)
        tree_redundant_samples.write(outfile=treefile_redundant_samples, format=2)

        # get samples
        phenotypes_file_redundant_samples = "%s/phenotypes_redundant_samples.tab"%outdir_gwas
        df_phenotypes_redundant = df_phenotypes_all[df_phenotypes_all.sampleID.isin(redundant_samples)][["sampleID", "phenotype"]].set_index("sampleID", drop=False).copy()
        save_df_as_tab(df_phenotypes_redundant.loc[sorted(redundant_samples)], phenotypes_file_redundant_samples)

        # get the files
        phenotypes_file_redundant, df_sample_prunning_redundant = get_balanced_samples_for_GWAS(phenotypes_file_redundant_samples, treefile_redundant_samples, ancestral_GWAS_dir, general_min_support, general_ASR_method, "%s/sample_selection_redundants"%outdir_gwas, list_cmds)
    
    else:
        phenotypes_file_redundant = None
        df_sample_prunning_redundant = pd.DataFrame({"n_transitions":[0]})

    return phenotypes_file_redundant, df_sample_prunning_redundant



def load_gff3_intoDF(gff_path, replace=False):

    """ Takes the path to a gff and loads it into a df"""

    gff_df_file = "%s.df.tab"%gff_path

    if file_is_empty(gff_df_file) or replace is True:

        # define the number of rows to skip
        gff_lines = "%s.gff_lines.gff"%gff_path
        run_cmd("egrep -v '^#' %s > %s"%(gff_path, gff_lines))

        # define the gff fields
        gff_fields = ["chromosome", "source", "feature", "start", "end", "blank1", "strand", "blank2", "annotation"]

        # load
        print("loading gff")
        gff = pd.read_csv(gff_lines, header=None, names=gff_fields, sep="\t")

        # set all possible annotations
        all_annotations = set.union(*[set([x.split("=")[0] for x in an.split(";")]) for an in gff.annotation])

        def list_to_str(list_x):

            if list_x==[]: return ""
            else: return list_x[0]
        
        for anno in all_annotations:

            # define the field 
            anno_field = "ANNOTATION_%s"%anno

            # add the normalQ annotation field
            gff[anno_field] = gff.annotation.apply(lambda x: list_to_str([a.split(anno)[1].lstrip("=") for a in x.split(";") if anno in a]))

            # add the Dbxref sepparated
            if anno=="Dbxref": 

                # get all the dbxref fields
                all_Dbxrefs = set.union(*[set([x.split(":")[0] for x in dbxref.split(",")]) for dbxref in gff[anno_field]]).difference({""})

                # go through each dbxref field and add it to the df
                for dbxref in all_Dbxrefs: 

                    gff["%s_%s"%(anno_field, dbxref)] = gff[anno_field].apply(lambda x: list_to_str([d.split(dbxref)[1].lstrip(":") for d in x.split(",") if dbxref in d]))

        # get the ID
        gff["ID"] = gff.ANNOTATION_ID
        gff["Parent"] = gff.ANNOTATION_Parent

        # change the ID so that all of the features are unique, add numbers for non-unique IDs
        gff["duplicated_ID"] = gff.duplicated(subset="ID", keep=False) # marks all the duplicated IDs with a True
        gff["numeric_idx"] = list(range(0, len(gff)))

        def getuniqueIDs_gff(row):

            """Takes a row and changes the IDs if they are not unique"""
            if row["duplicated_ID"] is False: return row["ID"]
            else: return "%s-%i"%(row["ID"], row["numeric_idx"])

        gff["ID"] = gff.apply(getuniqueIDs_gff, axis=1)

        # check that it is correct
        if len(gff)!=len(set(gff["ID"])): raise ValueError("IDs are not unique in the gff")

        # set the id as index
        gff = gff.set_index("ID", drop=False)

        # define all the IDs
        all_unique_IDs = set(gff.index)

        # add the upmost_parent (which should be the geneID)
        print("getting upmost parent")
        def get_utmost_parent(row, gff_df):

            """Takes a row and the gff_df and returns the highest parent. The index has to be the ID of the GFF"""

            # when you have found the parent it has no parent, so the ID is the upmost_parent
            if row["Parent"]=="": return row["ID"]

            # else you go to the parent
            else: 

                # define the ID of the parent
                if row["Parent"] in all_unique_IDs: parentID = row["Parent"] # normal

                elif sum(gff_df.ANNOTATION_ID==row["Parent"])>1: # split parents, take the first one
                    parentID = gff_df[gff_df.ANNOTATION_ID==row["Parent"]].ID.iloc[0]

                else: raise ValueError("%s is not a valid Parent ID"%(r["Parent"])) # errors

                # get the upmost parent from the parent
                return get_utmost_parent(gff_df.loc[parentID], gff_df)

        gff["upmost_parent"] = gff.apply(lambda row: get_utmost_parent(row, gff), axis=1)
        gff = gff.set_index("upmost_parent", drop=False)

        # add the type of upmost_parent
        df_upmost_parents = gff[gff.ID==gff.upmost_parent]

        # check that the upmost_parents are unique
        if len(df_upmost_parents)!=len(set(df_upmost_parents.ID)): raise ValueError("upmost parents are not unique")

        # map each upmost_parent to the feature
        upmost_parent_to_feature = dict(df_upmost_parents.set_index("upmost_parent")["feature"])

        # add the upmost_parent_feature to gff df
        gff["upmost_parent_feature"] = gff.upmost_parent.apply(lambda x: upmost_parent_to_feature[x])

        # write df
        gff_df_file_tmp = "%s.tmp"%gff_df_file
        gff.to_csv(gff_df_file_tmp, sep="\t", index=False, header=True)
        os.rename(gff_df_file_tmp, gff_df_file)


    # load
    gff = pd.read_csv(gff_df_file, sep="\t")

    return gff

def generate_interproscan_annotation_from_gff(gff, genome, outdir, gDNA_code, mtDNA_code, threads=4, replace=False):

    """This function generates Interproscan annotations for the gff provide"""

    ############ GENERATE A .fasta WITH THE PROTEIN-CODING TRANSCRIPTS #########    

    print("running generate_interproscan_annotation_from_gff")
    make_folder(outdir)

    # load gff
    gff_df = load_gff3_intoDF(gff)

    # keep only CDSs
    all_real_genes = set(gff_df[gff_df.feature=="gene"].upmost_parent)
    gff_df = gff_df[(gff_df.feature=="CDS") & (gff_df.upmost_parent.isin(all_real_genes))].copy()
    all_protein_coding_genes = set(gff_df.upmost_parent)

    # define the parent as the mRNA id
    gff_df["mRNA_id"] = gff_df.ANNOTATION_Parent

    # check that each gene (upmost_parent) has only one mRNA id
    upmost_parent_to_n_mRNAs = gff_df.groupby("upmost_parent").apply(lambda df_p: set(df_p.mRNA_id)).apply(len)
    if any(upmost_parent_to_n_mRNAs>1): raise ValueError("There are some upmost parents that have more than 1 mRNA:\n %s"%upmost_parent_to_n_mRNAs[upmost_parent_to_n_mRNAs>1]) 

    # check that each mRNA corresponds to one same upmost_parent
    mRNA_to_nUpmostParents = gff_df.groupby("mRNA_id").apply(lambda df_r: set(df_r.upmost_parent)).apply(len)
    if any(mRNA_to_nUpmostParents>1): raise ValueError("There are some mRNAs that have more than 1 upmost_parent:\n %s"%mRNA_to_nUpmostParents[mRNA_to_nUpmostParents>1]) 

    # add fields to the gff
    gff_df["numeric_row_idx"] = list(range(len(gff_df)))
    gff_df["annotation"] = "ID=" + gff_df.upmost_parent + "_cds" + gff_df.numeric_row_idx.apply(str) + ";Parent=" + gff_df.upmost_parent
    gff_fields = ["chromosome", "source", "feature", "start", "end", "blank1", "strand", "blank2", "annotation"]

    # write a gff with the CDSs
    cds_gff = "%s/cds_annotations.gff"%outdir
    header_gff = "##gff-version 3"
    gff_lines = gff_df[gff_fields].to_csv(sep="\t", header=False, index=False)
    open(cds_gff, "w").write(header_gff + "\n" + gff_lines)

    # get the transcripts
    spliced_transcripts = "%s.transcripts.fasta"%cds_gff
    run_cmd("gffread %s -g %s -x %s"%(cds_gff, genome, spliced_transcripts)) # -x has to be after -g

    # check that there is one transcript for each gene
    n_transcripts = int("".join([x for x in str(subprocess.check_output("egrep '^>' %s | wc -l"%spliced_transcripts, shell=True)) if x.isdigit()]))
    n_proteins = len(all_protein_coding_genes)
    if n_proteins!=n_transcripts: raise ValueError("there are %i prots and %i mRNAs"%(n_proteins, n_transcripts))

    # translate according to each code
    print("writing proteins")
    geneID_to_seq = {seq.id:seq for seq in SeqIO.parse(spliced_transcripts, "fasta")}

    mtDNA_chroms = {mtDNA_chrom}
    gDNA_chroms = set(gff_df.chromosome).difference(mtDNA_chroms)

    all_proteins = []
    for type_genome, chroms, translation_code in [("mtDNA", mtDNA_chroms, mtDNA_code), ("gDNA", gDNA_chroms, gDNA_code)]:

        # get the genes for this type of genome
        geneIDs = set(gff_df[gff_df.chromosome.isin(chroms)].upmost_parent)
        if len(geneIDs)==0: raise ValueError("There are no chroms in geneIDs")

        # init the number of strange prots
        n_strange_prots = 0

        # get the proteins
        for geneID in geneIDs:

            # define the protein seq (removing the last asterisk)
            protein_seq = str(geneID_to_seq[geneID].translate(table=translation_code).seq).rstrip("*")

            # check the sequence
            if ("*" in protein_seq) or ("X" in protein_seq): n_strange_prots+=1

            # replace the asterisks with 'X'
            protein_seq = protein_seq.replace("*", "X")

            # keep
            all_proteins.append(SeqRecord(Seq(protein_seq), id=geneID, name="", description=""))

        print("There are %i/%i proteins with bad characters in %s"%(n_strange_prots, len(geneIDs), type_genome))

    # write 
    proteins_fasta = "%s/proteins.fasta"%outdir
    SeqIO.write(all_proteins, proteins_fasta, "fasta")

    ############################################################################

    ######### RUN INTERPROSCAN TO ANNOTATE THE PROTEINS #######

    # define the outfile
    outfile_interpro = "%s/interproscan_annotation.out"%outdir
    if file_is_empty(outfile_interpro) or replace is True:
        print("running interproscan")

        # run interproscan
        outfile_interpro_tmp = "%s.tmp"%outfile_interpro
        tempdir = "%s/tmp_interproscan"%outdir; delete_folder(tempdir)
        interproscan_sh = "%s/software/interproscan_v5.52-86.0/interproscan-5.52-86.0/interproscan.sh"%ParentDir
        interproscan_libdir = "/gpfs/projects/bsc40/mschikora/anaconda3/envs/InterProScan_env/lib" # this needs to be added

        run_cmd("export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH && %s -i %s -f tsv -goterms -appl Pfam,ProSitePatterns,ProSiteProfiles,PANTHER,TIGRFAM,SFLD,SUPERFAMILY,Gene3D,Hamap,Coils,SMART,CDD,PRINTS,PIRSR,MobiDBLite,PIRSF --cpu %i -o %s --pathways --seqtype p --tempdir %s"%(interproscan_libdir, interproscan_sh, proteins_fasta, threads, outfile_interpro_tmp, tempdir), env="InterProScan_env")

        delete_folder(tempdir)
        os.rename(outfile_interpro_tmp, outfile_interpro)

    ###########################################################

def prepare_files_for_geno_pheno_validation(df_reads_info_analysis_with_validation, df_phenotypes_validation, df_metadata_binary_phenotypes, ValidationDir, CurDir, trees_data_df, VarCallDir, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, list_cmds, simple_repeats_file, gff_Cparapsilosis):

    """runs cmds for the validation of ML geno-pheno analyses"""
    print("Getting files G-P analyses")

    # keeps
    df_reads_info_analysis_with_validation = df_reads_info_analysis_with_validation.copy()
    df_phenotypes_validation = df_phenotypes_validation.copy()
    df_metadata_binary_phenotypes = df_metadata_binary_phenotypes.copy()
    trees_data_df = trees_data_df.copy()

    # define binary phenotypes for the SRA validation runs
    df_binary_phenotypes_validation = get_df_binary_phenotypes_validation(df_phenotypes_validation)

    # merge with the sequenced here samples
    validation_phenotypes = ["FLC_resistance", "VRC_resistance"]
    df_binary_phenotypes_validation = pd.concat([df_binary_phenotypes_validation, df_metadata_binary_phenotypes])[["runID"] + validation_phenotypes].copy()
    df_binary_phenotypes_validation["source_data"] = df_binary_phenotypes_validation.runID.apply(get_source_data)

    # define all valid samples for analysis. Only those with no smiley because we have to consider SVs and CNVs.
    all_analysis_runs_with_validation = set(df_reads_info_analysis_with_validation[df_reads_info_analysis_with_validation.has_smiley==False].runID)
    df_binary_phenotypes_validation = df_binary_phenotypes_validation[df_binary_phenotypes_validation.runID.isin(all_analysis_runs_with_validation)].copy()

    # define dirs
    GenoPheno_dir_all = ValidationDir+"/genotype_phenotype_relationships"; make_folder(GenoPheno_dir_all)

    # get commands for each
    for pheno in validation_phenotypes: # two phenotypes with top models and some validation data

        # get data
        df_phenotypes = df_binary_phenotypes_validation[df_binary_phenotypes_validation[pheno].isin({0, 1})][["runID", pheno, "source_data"]].copy()
        #print(df_phenotypes.groupby(["source_data", pheno]).apply(len))

        # get dir
        for aligner in ["bowtie2_local", "hisat2_no_spliced"]: # top models based on these
            print(pheno, aligner)

            # get dir
            GenoPheno_dir = "%s/pheno-%s_aligner-%s"%(GenoPheno_dir_all, pheno.replace(" ", "_"), aligner)
            make_folder(GenoPheno_dir)

            # get the tree for this subset of strains with this pheno
            tree_dir = "%s/tree_subset_strains"%GenoPheno_dir
            treefile_all_samples = generate_tree_one_phenotype_and_aligner(list_cmds, trees_data_df, tree_dir, aligner, set(df_phenotypes.runID), VarCallDir, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, replace=False)

            # debug
            if treefile_all_samples is None: continue

            # run convergence GWAS to get some files
            gwas_outdir = "%s/convergence_gwas"%GenoPheno_dir
            df_raw_gwas_files = get_cmds_convergence_GWAS(gwas_outdir, list_cmds, df_phenotypes, treefile_all_samples, VarCallDir, pheno, aligner, trees_data_df.copy(), homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, min_coverage_pos, simple_repeats_file, CurDir+"/data", gff_Cparapsilosis, None, replace=False, only_all_samples=True, skip_GWAS_run=True, DataDir_variants=ValidationDir)

def get_cmds_convergence_GWAS(outdir_gwas, list_cmds, df_phenotypes_all, treefile_all_samples, VarCallDir, phenotype, aligner, trees_data_df, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis, min_coverage_pos, simple_repeats_file, DataDir, gff_Cparapsilosis, cladeID_info_file, replace=False, only_all_samples=False, skip_GWAS_run=False, DataDir_variants="auto"):

    """runs the convergence-based GWAS for one phenotype"""

    # define a df with all GWAS files
    df_gwas_files_file = "%s/df_gwas_all_raw_results_files.py"%(outdir_gwas)    
    if file_is_empty(df_gwas_files_file) or replace is True:
        #print("working on '%s'"%(outdir_gwas.replace(ParentDir+"/Cparapsilosis_popgen/", "")))
        #print("working on '%s'"%(outdir_gwas))

        # keeps
        make_folder(outdir_gwas)
        trees_data_df = trees_data_df.copy()

        # change to sampleID
        df_phenotypes_all = df_phenotypes_all.copy().rename(columns={"runID":"sampleID", phenotype:"phenotype"})
        if len(df_phenotypes_all)!=len(set(df_phenotypes_all.sampleID)): raise ValueError("sample ID should be unique")

        # get the phenotypes file for all samples
        sorted_samples = sorted(set(df_phenotypes_all.sampleID))
        phenotypes_file_all_samples = "%s/phenotypes_all.tab"%outdir_gwas
        save_df_as_tab(df_phenotypes_all[["sampleID", "phenotype"]].set_index("sampleID", drop=False).loc[sorted_samples], phenotypes_file_all_samples)

        # define dirs
        ancestral_GWAS_dir = "%s/scripts/ancestral_GWAS/scripts"%ParentDir

        # generate files related to all samples
        if only_all_samples is False:

            # define general parameters
            general_min_support = 70
            general_ASR_method = "DOWNPASS"

            # get the samples that make for an optimal GWAS (removing redundancies)
            phenotypes_file_representative, df_sample_prunning = get_balanced_samples_for_GWAS(phenotypes_file_all_samples, treefile_all_samples, ancestral_GWAS_dir, general_min_support, general_ASR_method, "%s/sample_selection"%outdir_gwas, list_cmds)      
            if phenotypes_file_representative is None: return # run 

            # discard cases in which there are not enough transitions
            n_transitions_rep_samples = df_sample_prunning.sort_values(by=["prunning_step"]).iloc[-1].n_transitions
            if n_transitions_rep_samples<3: return

            # get a tree with the prunned samples
            representative_samples = set(get_tab_as_df_or_empty_df(phenotypes_file_representative).sampleID)
            #print("There are %i/%i representative samples"%(len(representative_samples), len(sorted_samples)))
            treefile_representative_samples = generate_tree_one_phenotype_and_aligner(list_cmds, trees_data_df, "%s/tree_representative_samples"%outdir_gwas, aligner, representative_samples, VarCallDir, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis)
            if treefile_representative_samples is None: return

            # get balanced samples and their tree, which have only 1 sample for monophyletic R/S node. This may be necessary for balancing
            phenotypes_file_balanced = get_balanced_samples_for_GWAS(phenotypes_file_all_samples, treefile_all_samples, ancestral_GWAS_dir, general_min_support, general_ASR_method, "%s/sample_selection_balanced"%outdir_gwas, list_cmds, nrepresentative_samples=1)[0]
            if phenotypes_file_balanced is None: return

            balanced_samples = set(get_tab_as_df_or_empty_df(phenotypes_file_balanced).sampleID)
            treefile_balanced_samples = generate_tree_one_phenotype_and_aligner(list_cmds, trees_data_df, "%s/tree_balanced_samples"%outdir_gwas, aligner, balanced_samples, VarCallDir, min_coverage_pos, homozygous_Corthopsilosis_srrs, genome_Cparapsilosis)
            if treefile_balanced_samples is None: return

        else:
            phenotypes_file_representative = treefile_representative_samples = None
            phenotypes_file_balanced = treefile_balanced_samples = None

        # define variants used in GWAS
        if DataDir_variants=="auto":
            DataDir_variants = DataDir

        integrated_small_vars = "%s/integrated_small_vars/small_vars_%s/variant_integration/integrated_small_variants.tab"%(DataDir_variants, aligner)
        small_vars_annot =  "%s/integrated_small_vars/small_vars_%s/small_vars_annot/annotated_variants_corrrectedGene.tab"%(DataDir_variants, aligner)
        integrated_SV_CNVs = "%s/integrated_SVs_CNVs/SV_CNV_%s/variant_integration/integrated_SVs_CNVs.tab"%(DataDir_variants, aligner)
        SV_CNVs_annot = "%s/integrated_SVs_CNVs/SV_CNV_%s/SV_CNV_annotation/annotated_variants_corrrectedGene.tab"%(DataDir_variants, aligner)

        # define parameters and annotations
        SV_CNV_filters = "--SV_min_minAF 0.1 --SV_min_maxAF 0.3 --CNV_max_relative_CN_deletion 0.5 --CNV_min_relative_CN_duplication 1.5 --CNV_max_relative_coverage_deletion 0.6 --CNV_min_relative_coverage_duplication 1.3 --CNV_min_len 600"
        aneuploidy_calling_filters = "--aneuploidies_max_fraction_N_bases 0.1 --aneuploidies_max_fraction_repeats 0.1 --aneuploidies_min_median_mappability 0.75 --min_fraction_chromosome_aneuploidy 0.5 --aneuploidies_cov_thresholds 0.6,1.4"
        interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir
        pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
        CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir

        paths_CNV_calling = "%s/paths_CNV_calling_all.tab"%outdir_gwas
        df_paths_CNV_calling = pd.DataFrame({sID : {"sampleID":sID, "call_CNVs_outdir":"%s/%s/aligner-%s/SV_CNV_calling/call_CNVs"%(VarCallDir, sID, aligner)} for sID in sorted_samples}).transpose().reset_index(drop=False)
        save_df_as_tab(df_paths_CNV_calling, paths_CNV_calling)

        # init dict with GWAS result files
        df_gwas_files = pd.DataFrame()
        some_gwas_cmd_run = False

        # log
        print("running GWAS on %s"%phenotype)

        # run analysis for each sample set
        for sample_set, phenotypes_file, treefile in [("representative_samples", phenotypes_file_representative, treefile_representative_samples), 
                                                      ("balanced_samples", phenotypes_file_balanced, treefile_balanced_samples), 
                                                      ("all_samples", phenotypes_file_all_samples, treefile_all_samples)]:

            # skip the non-all_samples, if necessary
            if only_all_samples is True and sample_set!="all_samples": continue

            # get all the inputs for GWAS prepared tree (variants, domains, genes, pathways, general ASR of phenotypes and mutations)
            outdir_gwas_inputs = "%s/gwas_inputs_%s"%(outdir_gwas, sample_set); make_folder(outdir_gwas_inputs)

            if file_is_empty("%s/pipeline_finished.txt"%outdir_gwas_inputs):
                list_cmds.append(("prepare_GWAS_inputs", 
                                  "%s/prepare_GWAS_inputs.py --outdir %s --phenotypes %s --treefile %s --threads 112 --integrated_small_vars %s --small_vars_filtering typical_diploid --repeats_file %s --small_vars_min_cov_pos %i --small_vars_annot %s --integrated_SV_CNVs %s --SV_CNVs_annot %s --gff %s --ref_genome %s %s --paths_CNV_calling %s %s --interproscan_output %s --pathway_annotations_dir %s --mtDNA_chrom %s --gDNA_code %i --mtDNA_code %i --min_nmuts_per_group 2 --min_nsamples_per_group 2 --CGD_chromosomal_features %s --reactome_organism_priority Saccharomyces_cerevisiae,Schizosaccharomyces_pombe"%(ancestral_GWAS_dir, outdir_gwas_inputs, phenotypes_file, treefile, integrated_small_vars, simple_repeats_file, min_coverage_pos, small_vars_annot, integrated_SV_CNVs, SV_CNVs_annot, gff_Cparapsilosis, genome_Cparapsilosis, SV_CNV_filters, paths_CNV_calling, aneuploidy_calling_filters, interproscan_output, pathway_annotations_dir, mtDNA_chrom, Cpara_gDNA_code, Cpara_mtDNA_code, CGD_chromosomal_features), 
                                  "ancestral_GWAS_env", 112))
                some_gwas_cmd_run = True
                continue

            # skip
            if skip_GWAS_run is True: continue

            # go through different configurations and run the GWAS
            for min_support in [10, 50, 70, 90]:
                for ASR_method in ["MPPA", "DOWNPASS", "MPPA,DOWNPASS"]:

                    # get cmd
                    outdir_gwas_results = "%s/gwas_results_%s_min_sup_%i_ASR_%s"%(outdir_gwas, sample_set, min_support, ASR_method)
                    if file_is_empty("%s/pipeline_finished.txt"%outdir_gwas_results):
                        
                        list_cmds.append(("run_GWAS", 
                                          "%s/run_GWAS.py --input_variant_groups %s/all_groupings.tab --treefile %s/corrected_rooted_tree.nw --ASR_phenotypes %s/phenotypes_ASR.out --ASR_variants_dir %s/ASR_all_variants --varID_mapping %s/variants_IDmapping.tab --varID_mapping_rep %s/mutation_to_representative_mutation.tab --resampled_phenotypes_pastml_out %s/resampled_phenotypes_pastml_df.py --threads 112 --outdir %s --min_support %i --ASR_method %s --classification_method_ASR all_ignoreNaN --all_pval_methods phenotypes --interesting_gwas_methods synchronous,phylogeny_agnostic,phyC --min_n_geno_transitions 2 --min_n_geno_gains 2 --min_n_geno 2 --cladeID_info %s --skip_QQ_plots"%(ancestral_GWAS_dir, outdir_gwas_inputs, outdir_gwas_inputs, outdir_gwas_inputs, outdir_gwas_inputs, outdir_gwas_inputs, outdir_gwas_inputs, outdir_gwas_inputs, outdir_gwas_results, min_support, ASR_method, cladeID_info_file),
                                          "ancestral_GWAS_env", 112)) # phyC,phylogeny_agnostic # RelToBranchLen

                        some_gwas_cmd_run = True
                        continue

                    # keep
                    else:
                        df_gwas_files_one = pd.DataFrame({0 : {"sample_set":sample_set, "min_support":min_support, "ASR_method":ASR_method, "raw_GWAS_stats":"%s/integrated_GWAS_stats.tab"%outdir_gwas_results, "phenotypes_file":phenotypes_file, "treefile":treefile, "outdir_gwas_inputs":outdir_gwas_inputs}}).transpose()
                        df_gwas_files = pd.concat([df_gwas_files, df_gwas_files_one]).reset_index(drop=True)

        # save or return
        if some_gwas_cmd_run is False:
            save_object(df_gwas_files, df_gwas_files_file)

        else:
            return

    # at the end, return
    return load_object(df_gwas_files_file)

def get_prioritized_feature_selection_methods():

    """Sort the feat selection methods by importance"""

    # init
    all_methods = []

    # add the forward
    for scoring in ["integrated_AUC"]:
        for tol in [0.1, 0.15]:
            all_methods.append("forward|scoring=%s|tol=%s"%(scoring, tol))

    # add with the RFECV
    all_methods += ["RFECV|scoring=integrated_AUC"]

    return all_methods

def get_prioritized_model_names():

    """Gets the list of all tested ML models, so that the preferred, are on top. This is based on the function get_all_ML_models_list()"""

    # init
    all_models = []

    # random forest
    for rf_bootstrap in ["True"]:
        all_models.append("RF|n_estimators=100|bootstrap=%s"%rf_bootstrap) 

    # gradient boost (GBoost)
    for gb_n_estimators in [100]:
        for gb_max_depth in [3, 5, 7, 9]:
            all_models.append("GBoost|loss=log_loss|n_estimators=%i|max_depth=%i"%(gb_n_estimators, gb_max_depth))

    # AdaBoost. I use SAMME.R, and not SAMME, because SAMME was giving some errors
    for adaB_n_estimators in [50]: # 200 as in Li 2021,but very slow, so I remove it
        all_models.append("AdaBoost|estimator=DecisionTreeClassifier_depth_1|n_estimators=%i|algorithm=SAMME.R"%(adaB_n_estimators))

    # tree parms (inspired by Almeida et. al. 2024 (hyperparameter tunning of random forest)
    for tree_max_depth in [3, 5, 7, 9, "none"]:
        for tree_min_samples_split in [3, 2]:
            for tree_min_samples_leaf in [2, 1]:
                for tree_criterion in ["gini"]: # ["gini", "entropy", "log_loss"]:
                    for tree_max_features in ["none"]: # "sqrt", "log2",

                        all_models.append("tree|min_samples_split=%i|min_samples_leaf=%i|criterion=%s|max_depth=%s|max_features=%s"%(tree_min_samples_split, tree_min_samples_leaf, tree_criterion, tree_max_depth, tree_max_features))

    # add SVC models
    all_models += ["SVC|kernel=linear"]
    for gamma_svc in ["auto", "scale"]:
        all_models.append("SVC|kernel=rbf|gamma=%s"%gamma_svc)

    # add KNN parms
    for n_neighbors_knn in [5, 3]:
        for weights_knn in ["uniform", "distance"]:
            all_models.append("KNN|leaf_size=30|n_neighbors=%i|weights=%s"%(n_neighbors_knn, weights_knn)) # leaf_size has impact on time

    # other s models
    all_models += ["LR|penalty=l2|solver=liblinear", "GaussianNB", "BernoulliNB"]

    return all_models

def get_df_transitions_info_GWAS(outdir, df_raw_gwas_files_all):

    """Gets the # of transitions GWAS"""
    print("get GWAS trans info...")

    df_raw_gwas_files_all = df_raw_gwas_files_all.set_index(["aligner", "min_support", "ASR_method"], drop=True).copy()
    make_folder(outdir)

    # get for each pheno
    inputs_fn = []
    for phenotype in sorted(set(df_raw_gwas_files_all.phenotype)):

        # get trans df
        df_transitions_info_GWAS_file = "%s/df_transitions_info_GWAS_%s.py"%(outdir, phenotype)
        df_raw_gwas_files_pheno = df_raw_gwas_files_all[df_raw_gwas_files_all.phenotype==phenotype].copy()
        df_transitions_info_GWAS = inputs_fn.append((df_transitions_info_GWAS_file, phenotype, df_raw_gwas_files_pheno, True))

    with multiproc.Pool(112) as pool:
        df_transitions_info_GWAS = pd.concat(pool.starmap(get_integrated_GWAS_df_stats_pheno_transitions, inputs_fn, chunksize=1))
        pool.close()
        pool.terminate()

    return df_transitions_info_GWAS

def get_df_transitions_info_GWAS_only_rep_samples(DataDir, df_metadata_binary_phenotypes, target_aligners=["bwa_mem"]):

    """gets data all phenos"""

    df_metadata_binary_phenotypes = df_metadata_binary_phenotypes.set_index("runID").copy()
    all_phenos = [k for k in df_metadata_binary_phenotypes.keys() if not k in {"sampleID", "runID"}]
    data_dict = {}; I=0

    # load
    for pheno in all_phenos:

        # discard
        n_1s = sum(df_metadata_binary_phenotypes[pheno]==1)
        n_0s = sum(df_metadata_binary_phenotypes[pheno]==0)

        if n_1s<3: continue
        if n_0s<3: continue

        # add
        for aligner in target_aligners:

            # define # transitions rep samples
            df_sample_prunning = get_tab_as_df_or_empty_df("%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas/sample_selection/prunning_info.tab"%(DataDir, pheno.replace(" ", "_"), aligner))
            n_transitions_rep_samples = df_sample_prunning.sort_values(by=["prunning_step"]).iloc[-1].n_transitions

            # add to df
            data_dict[I] = {"phenotype":pheno, "aligner":aligner, "n_transitions":n_transitions_rep_samples, "sample_set":"representative_samples"}
            I+=1

    return pd.DataFrame(data_dict).transpose()

def get_integrated_GWAS_df_stats_pheno_transitions(df_gwas_pheno_file, phenotype, df_raw_gwas_files_all, add_pheno=False):

    """Gets the df with phenotypes"""

    if file_is_empty(df_gwas_pheno_file):

        # init
        df_gwas_pheno = pd.DataFrame()

        # add
        all_sample_sets = ["all_samples", "representative_samples", "balanced_samples"]
        for aligner in ["bwa_mem", "bowtie2_local", "hisat2_no_spliced"]:

            if not (aligner, 70, "DOWNPASS") in set(df_raw_gwas_files_all.index):
                df_gwas_p = pd.DataFrame({"gwas_method":"synchronous", "nodes_withPheno":-1, "aligner":aligner, "sample_set":all_sample_sets})
                df_gwas_pheno = pd.concat([df_gwas_pheno, df_gwas_p.drop_duplicates()]).reset_index(drop=True)
                continue

            gwas_paths_df = df_raw_gwas_files_all.loc[(aligner, 70, "DOWNPASS")].copy().set_index("sample_set")
            for sample_set in all_sample_sets:
                print(aligner, sample_set)

                # load df
                df_gwas_p = get_tab_as_df_or_empty_df(gwas_paths_df.loc[sample_set].raw_GWAS_stats)[["gwas_method", "nodes_withPheno"]].copy()
                df_gwas_p["aligner"] = aligner
                df_gwas_p["sample_set"] = sample_set

                # filter
                df_gwas_p = df_gwas_p[df_gwas_p.gwas_method=="synchronous"].copy()
                if len(df_gwas_p)==0:
                    df_gwas_p = pd.DataFrame({"gwas_method":"synchronous", "nodes_withPheno":-1, "aligner":aligner, "sample_set":[sample_set]})

                if len(set(df_gwas_p.nodes_withPheno))!=1:
                    raise ValueError("there should be one pheno")

                # add
                df_gwas_pheno = pd.concat([df_gwas_pheno, df_gwas_p.drop_duplicates()]).reset_index(drop=True)

        # save
        if add_pheno is True: df_gwas_pheno["phenotype"] = phenotype
        save_object(df_gwas_pheno, df_gwas_pheno_file)
    return load_object(df_gwas_pheno_file)

def get_str_split_by_comma(x):
    if pd.isna(x): return []
    else: return str(x).split(",")

def get_formatted_str_nans(x):
    if pd.isna(x): return ""
    else: return str(x)

def load_df_from_parms_dict(parms_dict):

    """loads df with parms dict"""

    df = get_tab_as_df_or_empty_df(parms_dict["df_file"])

    for k,v in parms_dict.items():
        if k=="df_file": continue
        df[k] = v

    return df

def get_df_ML_results_formatted_from_tab_to_py(df_ML_results):

    """formatted df"""

    all_fields_lists = [eval_data+"_"+suff for suff in ["phenotypes", "probabilities", "samples", "ignored_samples"] for eval_data in ["cross_splits", "redundant_samples"]]

    for f in (["selected_predictors", "selected_predictors_redundant", "train_samples", "test_samples", "train_phenos", "test_phenos", "test_probabilities", "train_probabilities", "ignored_samples_train", "ignored_samples_test", "test_integrated_onlyCorrectModels_samples", "test_integrated_samples", "test_integrated_onlyCorrectModels_ignored_samples", "test_integrated_ignored_samples", "test_integrated_onlyCorrectModels_phenotypes", "test_integrated_phenotypes", "test_integrated_probabilities", "test_integrated_onlyCorrectModels_probabilities"] + all_fields_lists):
        df_ML_results[f] = df_ML_results[f].apply(lambda x: get_str_split_by_comma(x))

    for f in ["test_probabilities", "train_probabilities", "cross_splits_probabilities", "redundant_samples_probabilities", "test_integrated_probabilities", "test_integrated_onlyCorrectModels_probabilities"]:
        df_ML_results[f] = df_ML_results[f].apply(lambda x: list(map(float, x)))

    for f in ["test_phenos", "train_phenos", "cross_splits_phenotypes", "redundant_samples_phenotypes", "test_integrated_onlyCorrectModels_phenotypes", "test_integrated_phenotypes"]:
         df_ML_results[f] = df_ML_results[f].apply(lambda x: list(map(int, x)))

    for f in ["expected_genes_str", "all_affected_genes_str", "all_affected_vars_str", "extra_info_train", "extra_info_test", "cross_splits_extra_info"]:
        df_ML_results[f] = df_ML_results[f].apply(get_formatted_str_nans)


    return df_ML_results

def generate_redundant_samples_info_file(redundant_samples_info_file, sorted_aligners, phenotype, target_sample_sets, df_raw_gwas_files, ancestral_GWAS_dir):

    """Generates the file with red samples"""


    if file_is_empty(redundant_samples_info_file):
        print("Getting redundant samples")
        tmpdir_all = "%s.gen"%redundant_samples_info_file
        #delete_folder(tmpdir_all)
        make_folder(tmpdir_all)

        # load for different ligners and sample sets
        df_rsamples_info_all = pd.DataFrame()
        for aligner in sorted_aligners:

            # get paths
            df_paths = df_raw_gwas_files.loc[(aligner, 70, "DOWNPASS")].copy().rename(columns={"outdir_gwas_inputs":"GWAS_inputs_dir"})
            df_paths = df_paths[df_paths.sample_set.isin(target_sample_sets)].copy().set_index("sample_set")

            # load dfs
            df_phenos_all_samples = get_tab_as_df_or_empty_df(df_paths.loc["all_samples"].phenotypes_file)
            tree_all_samples = Tree(df_paths.loc["all_samples"].GWAS_inputs_dir + "/corrected_rooted_tree.nw")
            all_samples = set(df_phenos_all_samples.sampleID)

            for sample_set in target_sample_sets.difference({"all_samples"}):
                print("running on %s-%s"%(aligner, sample_set))
                tmpdir = "%s/%s_%s"%(tmpdir_all, aligner, sample_set); make_folder(tmpdir)

                # define the redundant samples
                samples_ss = set(get_tab_as_df_or_empty_df(df_paths.loc[sample_set].phenotypes_file).sampleID)
                possible_redundant_samples_ss = all_samples.difference(samples_ss)

                # get phenotypes
                phenotypes_file = "%s/phenotypes.tab"%tmpdir
                save_df_as_tab(df_phenos_all_samples[df_phenos_all_samples.sampleID.isin(possible_redundant_samples_ss)].copy(), phenotypes_file)

                # get tree
                tree_redundant_samples = cp.deepcopy(tree_all_samples)
                tree_redundant_samples.prune(possible_redundant_samples_ss, preserve_branch_length=True)
                treefile = "%s/tree.nw"%tmpdir
                tree_redundant_samples.write(outfile=treefile, format=2)

                # get cmd
                nrepresentative_samples = {"representative_samples":3, "balanced_samples":1}[sample_set]
                outdir_redundancy_red = "%s/red_red"%tmpdir; make_folder(outdir_redundancy_red)

                prunned_phenotypes_file = None
                while prunned_phenotypes_file is None:

                    list_cmds = []
                    prunned_phenotypes_file = get_balanced_samples_for_GWAS(phenotypes_file, treefile, ancestral_GWAS_dir, 70, "DOWNPASS", outdir_redundancy_red, list_cmds, nrepresentative_samples=nrepresentative_samples, extra_args="--skip_prunning_uncertain_politomies --skip_plots")[0]

                    if prunned_phenotypes_file is None: run_cmd(list_cmds[0][1], env="ancestral_GWAS_env")
                    else: break

                # add to df_rsamples_info_all
                df_rsamples_info = get_tab_as_df_or_empty_df(prunned_phenotypes_file).sort_values(by=["sampleID"])[["sampleID"]].copy()
                df_rsamples_info["sample_set"] = sample_set
                df_rsamples_info["aligner"] = aligner

                if len(df_rsamples_info)==0:
                    raise ValueError("invalid!!!")

                df_rsamples_info_all = pd.concat([df_rsamples_info_all, df_rsamples_info]).reset_index(drop=True)

        delete_folder(tmpdir_all)
        save_df_as_tab(df_rsamples_info_all, redundant_samples_info_file)

def get_all_gene_sets_combinations(phenotype):

    """gets the gene set combinations"""

    # define phenotypes
    azoles = ["FLC_resistance", "VRC_resistance", "POS_resistance", "IVZ_resistance", "ITR_resistance"]
    echinocandins = ["ANI_resistance", "CAS_resistance", "MIF_resistance"]
    virulence = ['agar_invasiveness', 'age_infant', 'age_old', 'age_young', 'biofilm_formation', 'microfluidics_Adherent', 'microfluidics_Big_Clumps', 'microfluidics_Little_clump', 'pseudohyphae_formation', 'source_ABSCESO', 'source_ABSCESO_ABDOMINAL', 'source_ASPIRADO_BRONQUIAL', 'source_ASPIRADO_TRAQUEAL', 'source_BIOPSIA_PULMONAR', 'source_EOF', 'source_ESPUTO', 'source_EXUDADO_CUTANEO', 'source_EXUDADO_DE_HERIDA', 'source_EXUDADO_OROFARINGEO', 'source_EXUDADO_OTICO', 'source_EXUDADO_RECTAL', 'source_EXUDADO_RESPIRATORIO_ALTO', 'source_HC', 'source_HERIDA', 'source_LAVADO_BRONCOALVEOLAR', 'source_LCR', 'source_LIQUIDO_ASCITICO', 'source_LIQUIDO_BILIAR', 'source_MUESTRA_AMBIENTAL', 'source_ORINA', 'source_PIEL', 'source_RESPIRATORIAS', 'source_SANGRE_COMPLETA', 'source_TEJIDO_VASCULAR', 'source_VALVULA_CARDIACA', 'source_blood', 'source_blood_vs_environmental', 'source_blood_vs_other_sterile', 'source_respiratory', 'source_skin', 'source_sterile', 'source_sterile_vs_environmental']

    if phenotype in (azoles + echinocandins):
        gene_sets = ["all",
                     "known", 
                     "known+by_description",
                     
                     "known+GWAS_high",
                     "known+GWAS_high+GWAS_low",
                     "known+by_description+GWAS_high",
                     "known+by_description+GWAS_high+GWAS_low",
                     
                     "known+selection_all",
                     "known+selection_SNP",
                     "known+by_description+selection_all",
                     "known+by_description+selection_SNP",

                     "known+by_description+GWAS_high+selection_SNP",
                     "known+by_description+GWAS_high+GWAS_low+selection_all"]

    elif phenotype in virulence:
        gene_sets = ["all",
                     "by_description",
                     "by_description+selection_all",
                     "by_description+selection_SNP"]


    return gene_sets



def get_min_ntransitons_across_aligners_one_phenotype(phenotype, DataDir):

    """gets the minimum number of transitions acorss aligners"""

    # load the stats df for this pheno
    df_stats_p = load_object("%s/simplified_GWAS_stats/%s_pheno_transitions.py"%(DataDir, phenotype))
    if len(df_stats_p)!=9: raise ValueError("not enogh")

    # get some filters
    target_sample_sets = {"all_samples", "representative_samples"}
    df_stats_p = df_stats_p[df_stats_p.sample_set.isin(target_sample_sets)].copy()

    # map each aligner to the # of transitions
    aligner_to_ntrans = df_stats_p.groupby("aligner").apply(lambda df: max(df.nodes_withPheno))

    # return the min number
    return min(aligner_to_ntrans)

def run_ML_predictors_phenotypes(df_raw_gwas_files_all, GenoPheno_dir_all, list_cmds, DataDir, gff_Cparapsilosis, genome_Cparapsilosis, sorted_aligners, df_ML_results_all_file, ProcessedDataDir, min_nmodels, replace=False):

    """Runs the ML predictors for different parameters"""
    print("Running ML predictors...")

    # define the final file, and exit if already got
    if not file_is_empty(df_ML_results_all_file) and replace is False:
        print("ML results already generated")
        return

    # keep index
    df_raw_gwas_files_all = df_raw_gwas_files_all.set_index(["aligner", "min_support", "ASR_method"], drop=True).copy()

    # define general dirs
    ancestral_GWAS_dir = "%s/scripts/ancestral_GWAS/scripts"%ParentDir
    pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
    interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir
    CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
    Integrated_ML_classifier_dirs = "%s/ML_classifiers_integrated_results"%DataDir; make_folder(Integrated_ML_classifier_dirs)

    # define the min # transitions to do predictions
    min_pheno_transitions = 5

    # get expected gene tuples
    phenotype_to_expected_genes_tuple = get_phenotype_to_expected_genes(DataDir)
    df_expected_gene_info = get_phenotype_to_all_possible_genes_for_ML(DataDir, ParentDir+"/Cparapsilosis_popgen/manual_data", "%s/df_expected_gene_info_ML.py"%ProcessedDataDir)
    df_expected_gene_info = df_expected_gene_info[df_expected_gene_info.phenotype.isin(set(df_raw_gwas_files_all.phenotype))].copy().set_index("phenotype")

    # init df
    df_ML_results_all = pd.DataFrame()
    #df_ML_metrics_all = pd.DataFrame()

    # init files remove
    files_remove = []

    # go through different phenotypes
    for phenotype in sorted(set(df_raw_gwas_files_all.phenotype)):
        print(phenotype)
        

        # discard some manually curated phenotypes, which may be redundant with others
        if phenotype=="source_SANGRE_COMPLETA": continue # this is redundant with 'source_blood', and less informative

        # debug
        #if phenotype!="source_blood_vs_other_sterile": continue # debug!!!!!
        #if phenotype in {"VRC_resistance", "source_blood", "source_sterile"}: continue # debug!!!

        # get df
        df_raw_gwas_files = df_raw_gwas_files_all[df_raw_gwas_files_all.phenotype==phenotype].copy()

        # init a df with the paths of different training sets
        list_training_data = []
        all_training_done = True

        # discard phenotypes that have not enough transitions by DOWNPASS,70
        gwas_results_dir = "%s/simplified_GWAS_stats"%DataDir; make_folder(gwas_results_dir)
        df_stats_p = get_integrated_GWAS_df_stats_pheno_transitions("%s/%s_pheno_transitions.py"%(gwas_results_dir, phenotype), phenotype, df_raw_gwas_files.copy())
        if len(df_stats_p)!=9: raise ValueError("not enogh")
        target_sample_sets = {"all_samples", "representative_samples"}
        df_stats_p = df_stats_p[df_stats_p.sample_set.isin(target_sample_sets)].copy()
        aligners_sufficient_transitions = [aln for aln in sorted_aligners if any(df_stats_p[df_stats_p.aligner==aln].nodes_withPheno>=min_pheno_transitions)]

        if len(aligners_sufficient_transitions)<3:
        #if len(aligners_sufficient_transitions)==0: -> old, with all
            #print("WARNING: For %s, there are no aligners with sufficient transitions. Skipping..."%phenotype)
            continue

        # get info
        df_expected_gene_info_pheno = df_expected_gene_info.loc[phenotype].set_index("gene_group").copy()

        # log
        #print("Getting predictors for %s, aligners: %s. transitions: %s"%(phenotype, aligners_sufficient_transitions, [max(df_stats_p[df_stats_p.aligner==aln].nodes_withPheno) for aln in aligners_sufficient_transitions]))

        # define the expected genes
        if phenotype in phenotype_to_expected_genes_tuple:
            expected_genes_str = ";".join(["%s,%s"%(x[0], x[1]) for x in phenotype_to_expected_genes_tuple[phenotype]])

        else:
            #print("WARNING: No expected genes for %s"%phenotype)
            expected_genes_str = "skip"

        # go through different parameter combinations
        for aligner, min_support, ASR_method in sorted(set(df_raw_gwas_files.index)):
            if not aligner in aligners_sufficient_transitions: continue
            #print(phenotype, aligner, min_support, ASR_method, len(list_cmds))

            # define dir
            GenoPheno_dir = "%s/pheno-%s_aligner-%s"%(GenoPheno_dir_all, phenotype.replace(" ", "_"), aligner)
            outdir_ML_predictors_all = "%s/ML_predictors"%GenoPheno_dir; make_folder(outdir_ML_predictors_all)

            # go through different definitions and generate the inputs
            transition_definition_to_inputs = {}
            for transition_definition in ['each_pheno_clade', 'two_close_balanced']:             

                # define outdirs
                outdir_ML_predictors = "%s/min_sup-%i_ASR-%s_trans-%s"%(outdir_ML_predictors_all, min_support, ASR_method, transition_definition)
                make_folder(outdir_ML_predictors)

                # get df with paths
                df_paths = df_raw_gwas_files.loc[(aligner, min_support, ASR_method)].copy().rename(columns={"outdir_gwas_inputs":"GWAS_inputs_dir"})
                df_paths = df_paths[df_paths.sample_set.isin(target_sample_sets)].copy()
                paths_file = "%s/input_paths.tab"%outdir_ML_predictors
                save_df_as_tab(df_paths[["sample_set", "phenotypes_file", "raw_GWAS_stats", "GWAS_inputs_dir"]], paths_file)

                # prepare input files
                outdir_inputs = "%s/input_files"%outdir_ML_predictors
                if file_is_empty("%s/pipeline_finished.txt"%outdir_inputs): # takes 5-6h
                    cmd =  "%s/prepare_ML_classifier_inputs.py --outdir %s --min_support %s --ASR_method %s --classification_method_ASR all_ignoreNaN --paths_file %s --threads 112 --min_pheno_transitions %i --pathway_annotations_dir %s --gff %s --interproscan_output %s --CGD_chromosomal_features %s --reactome_organism_priority Saccharomyces_cerevisiae,Schizosaccharomyces_pombe --ref_genome %s --mem_per_thread 2 --n_phenotype_resamples 10 --n_CVset_resamples 5 --min_geno_transitions_GWAS_hit 2 --min_geno_pheno_transitions 2 --transition_definition %s %s"%(ancestral_GWAS_dir, outdir_inputs, min_support, ASR_method, paths_file, min_pheno_transitions, pathway_annotations_dir, gff_Cparapsilosis, interproscan_output, CGD_chromosomal_features, genome_Cparapsilosis, transition_definition, {"each_pheno_clade":"", "two_close_balanced":"--skip_feature_selection"}[transition_definition])

                    list_cmds.append(("ML_classifier_prepare", cmd, "ancestral_GWAS_env", 112))
                    all_training_done = False
                    continue

                # keep
                transition_definition_to_inputs[transition_definition] = outdir_inputs

            # debug if not both are found
            if len(transition_definition_to_inputs)!=2: continue

            # debug if any of the outdir inputs have no results
            if any([set(os.listdir(d))=={"df_paths_checked.py", "pipeline_finished.txt"} for d in transition_definition_to_inputs.values()]): 

                if ASR_method=="DOWNPASS" and min_support==70: # should not be like this
                    raise ValueError("not possible for DOWNPASS-70")
                continue # skip

            # link the each_pheno_clade to two_close_balanced
            soft_link_files("%s/target_features_CVsets.tab"%transition_definition_to_inputs["each_pheno_clade"], "%s/target_features_CVsets.tab"%transition_definition_to_inputs["two_close_balanced"])

            # map each sample set to the default_probability_pheno
            sample_set_to_pheno_prob = {}
            for I, r_ss in df_paths.iterrows():
                df_phenos_ss = get_tab_as_df_or_empty_df(r_ss.phenotypes_file)
                sample_set_to_pheno_prob[r_ss.sample_set] = sum(df_phenos_ss.phenotype==1) / len(df_phenos_ss)

            # train ML classifiers on CVsets of different types
            for transition_definition, outdir_inputs in transition_definition_to_inputs.items():
                #print("training cmds for", aligner, min_support, ASR_method, transition_definition)

                # define dirs
                outdir_training_all = "%s/min_sup-%i_ASR-%s_trans-%s_training_CVsets"%(outdir_ML_predictors_all, min_support, ASR_method, transition_definition)
                make_folder(outdir_training_all)

                # get cmds
                for sample_set in sorted(df_paths.sample_set):
                    for type_phenotype in ["phenotype_transition", "phenotype"]:

                        # debug to keep only one pheno
                        if type_phenotype=="phenotype" and transition_definition=="two_close_balanced": continue

                        # go each split
                        for test_split_ID in range(min_pheno_transitions):

                            # go through different gene sets
                            if min_support==70 and ASR_method=="DOWNPASS": all_gene_sets = get_all_gene_sets_combinations(phenotype)
                            else: 
                                all_gene_sets = ["all"]

                            for gene_set in all_gene_sets:

                                # go through different features_max_pval_GWAS
                                if gene_set=="all":
                                    features_max_pval_GWAS_list = [0.05, 1.0]
                                    gene_set_outdir_suffix = ""
                                    target_genes_training = "all"

                                else:
                                    features_max_pval_GWAS_list = [-1]
                                    gene_set_outdir_suffix  = "_geneSet=%s"%gene_set
                                    target_genes_training = sorted(set.union(*[set(df_expected_gene_info_pheno.loc[gg].geneID) for gg in gene_set.split("+")]))

                                for features_max_pval_GWAS in features_max_pval_GWAS_list:

                                    # get the outdir
                                    outdir_training = "%s/ss=%s_sID=%i_%s_maxp=%s%s"%(outdir_training_all, sample_set, test_split_ID, type_phenotype, features_max_pval_GWAS, gene_set_outdir_suffix); make_folder(outdir_training)

                                    # get the target geneIDs for training
                                    if target_genes_training=="all":
                                        target_genes_training_arg = "all"

                                    else:
                                        target_genes_training_file = "%s/target_genes_training.txt"%(outdir_training)
                                        if file_is_empty(target_genes_training_file):
                                            open(target_genes_training_file+".tmp", "w").write("\n".join(target_genes_training))
                                            os.rename(target_genes_training_file+".tmp", target_genes_training_file)

                                        target_genes_training_arg = "%s"%target_genes_training_file

                                    # remove some files
                                    #for f in ["integrated_models_paths.tab", "pipeline_finished.txt", "integrated_models", "tmp"]:
                                    #    if os.path.isdir("%s/%s"%(outdir_training, f)) or os.path.isfile("%s/%s"%(outdir_training, f)):
                                    #        files_remove.append("%s/%s"%(outdir_training, f))
                                    #continue

                                    # get cmd
                                    if file_is_empty("%s/pipeline_finished.txt"%outdir_training):
                                        cmd = "%s/train_ML_classifier_CVset.py --outdir %s --threads 112 --sample_set %s --test_split_ID %i --type_phenotype %s --features_max_pval_GWAS %s --input_files %s --pathway_annotations_dir %s --gff %s --interproscan_output %s --CGD_chromosomal_features %s --reactome_organism_priority Saccharomyces_cerevisiae,Schizosaccharomyces_pombe --ref_genome %s --expected_genes \"%s\" --SNP_kb_thresholds 0.1,0.25,1.0 --threshold_definition_transitions 0.25,0.5,0.75 --default_probability_pheno %.2f --min_AUC_feature_sel 0.6 --min_n_features 2 --target_genes_training %s"%(ancestral_GWAS_dir, outdir_training, sample_set, test_split_ID, type_phenotype, features_max_pval_GWAS, outdir_inputs, pathway_annotations_dir, gff_Cparapsilosis, interproscan_output, CGD_chromosomal_features, genome_Cparapsilosis, expected_genes_str, sample_set_to_pheno_prob[sample_set], target_genes_training_arg)

                                        all_training_done = False
                                        list_cmds.append(("ML_training_CVset", cmd, "ancestral_GWAS_env", 112))


                                    else:
                                        # debug some cases
                                        if os.listdir(outdir_training)==["pipeline_finished.txt"]: continue

                                        # keep paths
                                        list_training_data.append({"df_file":"%s/integrated_models_paths.tab"%outdir_training, "ASR_method":ASR_method, "aligner":aligner, "min_support":min_support, "gene_set":gene_set, "transition_definition":transition_definition, "outdir_inputs":transition_definition_to_inputs[transition_definition], "type_phenotype":type_phenotype, "default_probability_pheno":sample_set_to_pheno_prob[sample_set]})

        # integrate
        if all_training_done is True:
            print("Integrating data...")

            # generate outdir
            classifier_integration_outdir_all = "%s/pheno-%s"%(Integrated_ML_classifier_dirs, phenotype); make_folder(classifier_integration_outdir_all)

            # get the redundant samples for representative samples
            redundant_samples_info_file = "%s/redundant_samples_info.tab"%classifier_integration_outdir_all
            generate_redundant_samples_info_file(redundant_samples_info_file, sorted_aligners, phenotype, target_sample_sets, df_raw_gwas_files.copy(), ancestral_GWAS_dir)

            # integrate for different combinations of gene_set and aligner, to reduce the job size
            for target_gene_set, target_aligner in itertools.product(get_all_gene_sets_combinations(phenotype), sorted_aligners):
                #print(target_gene_set, target_aligner)

                # define outdir
                classifier_integration_outdir = "%s/gene_set=%s_aligner=%s"%(classifier_integration_outdir_all, target_gene_set, target_aligner)
                make_folder(classifier_integration_outdir)

                # generate paths_training_file
                paths_training_file = "%s/paths_training.tab"%(classifier_integration_outdir)

                if file_is_empty(paths_training_file):

                    # get the target files
                    list_training_data_target = [x for x in list_training_data if x["gene_set"]==target_gene_set and x["aligner"]==target_aligner]
                    if len(list_training_data_target)==0: raise ValueError("invalid")

                    # integrate
                    with multiproc.Pool(112) as pool:
                        df_paths_training_all = pd.concat(pool.starmap(load_df_from_parms_dict, list(map(lambda x: tuple([x]), list_training_data_target)), chunksize=1))
                        pool.close()
                        pool.terminate()

                    # save df
                    check_no_nans_in_df(df_paths_training_all)
                    save_df_as_tab(df_paths_training_all, paths_training_file)

                # run integration
                outdir_ML_integration = "%s/classifier_results"%(classifier_integration_outdir)

                # debug to run again
                #for f in os.listdir(outdir_ML_integration):
                #    if not f in {"best_models_each_CVset.py", "best_models_test_performance.py", "best_models_test_threshold_metrics.tab"}:
                #        files_remove.append("%s/%s"%(outdir_ML_integration, f))
                #continue

                if file_is_empty("%s/pipeline_finished.txt"%outdir_ML_integration):
                    cmd = "%s/integrate_ML_classifier_results.py --outdir %s --threads 112 --paths_training %s --n_phenotype_resamples 10 --gff %s --pathway_annotations_dir %s --ref_genome %s --redundant_samples_info %s"%(ancestral_GWAS_dir, outdir_ML_integration, paths_training_file, gff_Cparapsilosis, pathway_annotations_dir, genome_Cparapsilosis, redundant_samples_info_file)

                    list_cmds.append(("ML_integration", cmd, "ancestral_GWAS_env", 112))
                    continue

                # get dfs
                df_ML_results = get_tab_as_df_or_empty_df("%s/best_models_test_performance_integrated_min_nmodels=%i.tab"%(outdir_ML_integration, min_nmodels))

                # check
                if len(df_ML_results) not in {288}:
                    raise ValueError("invalid len %i"%(len(df_ML_results)))

                # add
                df_ML_results["phenotype"] = phenotype
                df_ML_results_all = pd.concat([df_ML_results_all, df_ML_results]).reset_index(drop=True)
                check_no_nans_in_df(df_ML_results_all, ignore_keys={"types_phenotype"}) 

                """

                # add df
                df_ML_results = get_tab_as_df_or_empty_df("%s/best_models_test_performance.tab"%outdir_ML_integration)
                df_ML_metrics = get_tab_as_df_or_empty_df("%s/best_models_test_threshold_metrics.tab"%outdir_ML_integration)

                df_ML_results["phenotype"] = phenotype
                df_ML_metrics["phenotype"] = phenotype

                # format some things
                df_ML_results = get_df_ML_results_formatted_from_tab_to_py(df_ML_results)

                for f in ["predicted_phenotypes"]:
                    df_ML_metrics[f] = df_ML_metrics[f].apply(get_formatted_str_nans)

                # keep
                df_ML_results_all = pd.concat([df_ML_results_all, df_ML_results]).reset_index(drop=True)
                df_ML_metrics_all = pd.concat([df_ML_metrics_all, df_ML_metrics]).reset_index(drop=True)

                # checks
                check_no_nans_in_df(df_ML_results_all) 
                check_no_nans_in_df(df_ML_metrics_all) 

                """

    # remove files
    if len(files_remove)>0:
        print("Removing %i files"%len(files_remove))
        with multiproc.Pool(112) as pool:
            pool.starmap(delete_file_or_folder, list(map(lambda x: tuple([x]), files_remove)), chunksize=1)
            pool.close()
            pool.terminate()
        print("all_removed")

    # save at the end
    if len(list_cmds)==0:

        print("saving final ML results")        
        check_no_nans_in_df(df_ML_results_all, ignore_keys={"types_phenotype"})
        #check_no_nans_in_df(df_ML_metrics_all) 

        #save_object(df_ML_metrics_all, "%s.per_threshold_metrics.py"%df_ML_results_all_file)
        save_object(df_ML_results_all, df_ML_results_all_file)

def get_type_gene_set_ML(x):

    """gets type"""

    if x=="all": 
        return "all genes"

    elif x=="known":
        return "expected genes"

    else:
        return "hypothesized genes"

def get_short_modelling_method(x):

    """gets the modelling method short"""

    # parse
    mod_name = x.split("|")[0]
    parms_dict = {y.split("=")[0] : y.split("=")[1] for y in x.split("|")[1:]}

    # define relevant fields
    name_to_fs = {"GBoost":["max_depth"],
                  "KNN": ["n_neighbors", "weights"],
                  "SVC": ["kernel", "gamma"],
                  "tree": ["min_samples_split", "min_samples_leaf", "max_depth"]}

    # get
    if mod_name in {"AdaBoost", "BernoulliNB", "GaussianNB", "LR", "RF"}:
        return mod_name

    elif mod_name in name_to_fs:
        return mod_name + "|" + ";".join([(f.replace("min_samples_split", "samples/split").replace("min_samples_leaf","samples/leaf").replace("n_neighbors","neighbors"))+"="+parms_dict[f] for f in name_to_fs[mod_name]])

    else: raise ValueError("invalid")


def add_transitions_hyperparameters(r):

    """get transitions definition"""

    if r.real_type_phenotype=="phenotype":
        for f in ["transitions definition", "transitions SNP/kb threshold", "transitions for prediction", "transitions ML threshold"]:
            r[f] = "none"

    elif r.real_type_phenotype.split("|")[0]=="phenotype_from_transitions":
        r["transitions definition"] = {"each_pheno_clade":"each clade", "two_close_balanced":"two close"}[r.transition_definition]
        
        parms_pheno = {x.split(":")[0]: x.split(":")[1] for x in r.real_type_phenotype.split("|")[1:]} 
        r["transitions SNP/kb threshold"] = parms_pheno["SNPkb_threshold"]
        r["transitions for prediction"] = {"all":"all", "onlyNoTrans":"only no trans.", "onlyTrans":"only trans."}[parms_pheno["background_transitions"]]
        r["transitions ML threshold"] = parms_pheno["ML_threshold_definition"]

    else: raise ValueError("invalid")


    return r

def get_df_ML_results_top_with_affected_gene_names(df_ML_results_top, DataDir):

    """add gene names"""

    # add as set the genes
    def get_set_sel_pred(x):
        if x=="nan": return set()
        else:
            preds = set(x.split(","))
            if not all(["|||" in x for x in preds]): raise ValueError("invalid %s"%x)
            return preds

    df_ML_results_top["selected_predictors"] = df_ML_results_top.selected_predictors.apply(str).apply(get_set_sel_pred).apply(get_selected_predictors_demult_one_tuple_predictors)
    df_ML_results_top["all_affected_genes_set"] = df_ML_results_top.selected_predictors.apply(lambda x: set(map(get_affected_geneID_one_unique_groupID, x)))


    # add name
    gene_features_df = get_gene_features_df(DataDir).set_index("geneID", drop=False)
    geneID_to_name = dict(gene_features_df.final_name)
    tac1_ID = gene_features_df[gene_features_df.final_name=="Scer_HAL9"].iloc[0].geneID
    geneID_to_name[tac1_ID] = "TAC1"
    df_ML_results_top["affected_gene_names_set"] = df_ML_results_top.all_affected_genes_set.apply(lambda x: {geneID_to_name[g] for g in x})

    return df_ML_results_top

def get_nan_to_zero(x ,fun_app=int):
    if pd.isna(x): return 0
    else: return fun_app(x)

def plot_train_parameter_performance_top_models(ML_df_dicts, PlotsDir, test_sample_set, ML_pheno_to_color, DataDir):

    """get the train parms representation.
    
    feature selection:
    - GWAS min support: skip, as all the top models are based on some gene set
    - GWAS ASR: skip, as all the top models are based on some gene set
    - GWAS type filtering: skip, as all the top models are based on some gene set

    - prioritize expected genes: 'no', 'yes' (all considered)
    - FS AUC threshold: 0.05, 0.1, 0.15 (all considered)
    
    modelling:
    - modelling method: a few used:

        AdaBoost
        BernoulliNB
        GBoost|max_depth=3
        GaussianNB
        KNN|neighbors=3;weights=distance
        KNN|neighbors=3;weights=uniform
        KNN|neighbors=5;weights=uniform
        LR
        RF
        SVC|kernel=rbf;gamma=auto
        tree|samples/split=3;samples/leaf=2;max_depth=3

    - class weight: yes/no (all considered)
    - type phenotype: phenotype / pheno. transition (all considered)

    transitions hyperparameters (all)
    - "transitions definition"
    - "transitions SNP/kb threshold"
    - "transitions for prediction"
    - "transitions ML threshold"

    """

    # map phenos
    paper_pheno_to_pheno = {get_paper_dicotomic_pheno(pheno):pheno for pheno in ML_pheno_to_color.keys()}

    # get df
    df_ML_results_top = get_filtered_df_by_dict(ML_df_dicts["df_ML_results_top"], {"test_sample_set":test_sample_set})

    # add affected genes by features
    df_ML_results_top = get_df_ML_results_top_with_affected_gene_names(df_ML_results_top, DataDir)
    #all_gene_names = set.union(*df_ML_results_top.affected_gene_names_set)
    #gene_to_phenotypes = {g : set(df_ML_results_top[(df_ML_results_top.affected_gene_names_set.apply(lambda x: g in x))].phenotype) for g in all_gene_names}

    # redefine fields
    redefine_dict = {"prioritize_expected_genes": "prioritize expected genes", "CVset_type_split":"type split", "sample_set":"train samples"}
    df_ML_results_top = df_ML_results_top.rename(columns=redefine_dict)

    df_ML_results_top["feat. sel. AUC threshold"] = df_ML_results_top.feature_selection_method.apply(lambda x: float(x.split("=")[-1]))
    df_ML_results_top["type phenotype"] = df_ML_results_top.type_phenotype.apply(lambda x: x.replace("_", " "))
    df_ML_results_top["modelling method"] = df_ML_results_top.model_name.apply(get_short_modelling_method)
    df_ML_results_top["class weight"] = df_ML_results_top.class_weight.map({"balanced":"yes", "none":"no"})
    df_ML_results_top = df_ML_results_top.apply(add_transitions_hyperparameters, axis=1)
    df_ML_results_top["consistency"] = df_ML_results_top.consistency_measure.apply(lambda x: x.replace("jaccard_", ""))
    df_ML_results_top["AUC diff"] = df_ML_results_top.AUC_diff_threshold.apply(str)

    # define fields
    fields_optimize_test = ["aligner", "type split", "consistency", "AUC diff", "train samples", "gene_set"]
    fields_optimize_train = ["prioritize expected genes", "feat. sel. AUC threshold", 
                             "modelling method", "class weight", "type phenotype",
                             "transitions definition", "transitions SNP/kb threshold", "transitions for prediction", "transitions ML threshold"]

    df_ML_results_top = df_ML_results_top[["test_split_ID", "test_parmID_exceptTestSS", "phenotype", "best_model", "affected_gene_names_set", "cross_splits_ROC_AUC", "validation_transitions_TrainTest_ROC_AUC", "ROC_AUC_test_integrated_onlyCorrectModels"] + fields_optimize_test + fields_optimize_train].copy()

    # change
    df_ML_results_top["validation_transitions_TrainTest_ROC_AUC"] = df_ML_results_top["validation_transitions_TrainTest_ROC_AUC"].apply(get_nan_to_zero,fun_app=float)
    check_no_nans_in_df(df_ML_results_top)


    # sort the gene names by phenotype
    sorted_gene_names = []
    for pheno in ["age_infant", "source_blood_vs_other_sterile", "FLC_resistance", "VRC_resistance"]:

        df_ML_results_top_pheno = df_ML_results_top[df_ML_results_top.phenotype==pheno].copy()
        gen_names_pheno = set.union(*df_ML_results_top_pheno.affected_gene_names_set)
        gene_to_nmodels = pd.Series({g : sum(df_ML_results_top_pheno.affected_gene_names_set.apply(lambda x: g in x)) for g in gen_names_pheno}).sort_values(ascending=False)

        sorted_gene_names_pheno = [g for g in gene_to_nmodels.index if not g in sorted_gene_names]
        sorted_gene_names += sorted_gene_names_pheno

        # prints
        print(pheno, "ROC_AUC cross val", min(df_ML_results_top_pheno.cross_splits_ROC_AUC), "-", max(df_ML_results_top_pheno.cross_splits_ROC_AUC), "median=", "%.2f"%np.median(df_ML_results_top_pheno.cross_splits_ROC_AUC))
        print(pheno, "ROC_AUC validation", min(df_ML_results_top_pheno.validation_transitions_TrainTest_ROC_AUC), "-", max(df_ML_results_top_pheno.validation_transitions_TrainTest_ROC_AUC), "median=", "%.2f"%np.median(df_ML_results_top_pheno.validation_transitions_TrainTest_ROC_AUC))

    if len(sorted_gene_names)!=len(set(sorted_gene_names)): raise ValueError("invalid")

    # create the heatmap df with parameters
    df_square = pd.DataFrame()
    #sorted_cols = df_ML_results_top[["test_split_ID", "test_parmID_exceptTestSS"]].sort_values(on=)
    df_ML_results_top["colID"] = df_ML_results_top.phenotype.apply(get_paper_dicotomic_pheno) + "|" + df_ML_results_top.test_parmID_exceptTestSS + "|" + df_ML_results_top.test_split_ID.apply(str)
    sorted_cols = sorted(set(df_ML_results_top["colID"]))
    fields_plot = fields_optimize_test + fields_optimize_train

    for f in fields_plot:

        # get df
        df_square_f  = df_ML_results_top[["colID", f]].copy()
        df_square_f["val_field"] = 1
        df_square_f = df_square_f.pivot(columns="colID", index=f, values="val_field").map(get_nan_to_zero)

        val_to_n = df_square_f.apply(sum, axis=1).sort_values(ascending=False)
        df_square_f = df_square_f.loc[list(val_to_n.index)].copy()

        # add header
        header_df = pd.DataFrame({c : {f:-1} for c in sorted_cols})
        df_square_f = pd.concat([header_df, df_square_f])

        # keep
        df_square = pd.concat([df_square, df_square_f]) 

    # add the genes
    df_square = pd.concat([df_square, pd.DataFrame({c : {"feature genes":-1} for c in sorted_cols})])
    colID_to_genes = dict(df_ML_results_top.set_index("colID").affected_gene_names_set)

    for g in sorted_gene_names:

        row_vals =  {c : {True:1, False:0}[g in colID_to_genes[c]] for c in sorted_cols}
        row_df = pd.DataFrame({g : row_vals}).transpose()
        df_square = pd.concat([df_square, row_df]) 

    # add the metrics
    df_square = pd.concat([df_square, pd.DataFrame({c : {"performance metrics":-1} for c in sorted_cols})])
    metric_to_short_metric = {"cross_splits_ROC_AUC":"AUC cross-val", "validation_transitions_TrainTest_ROC_AUC":"AUC validation", "ROC_AUC_test_integrated_onlyCorrectModels":"AUC test"}
    short_metric_to_metric = {v:k for k,v in metric_to_short_metric.items()}
    sorted_metrics = ["ROC_AUC_test_integrated_onlyCorrectModels", "cross_splits_ROC_AUC", "validation_transitions_TrainTest_ROC_AUC"]
    sorted_short_metrics = [metric_to_short_metric[m] for m in sorted_metrics]

    df_square_metrics = pd.DataFrame(0, index=sorted_short_metrics, columns=sorted_cols)
    df_square = pd.concat([df_square, df_square_metrics]) 

    # keep
    df_square = df_square[sorted_cols].copy()

    # edit the index
    def get_edited_index(idx, r):
        if idx in (fields_plot + ["feature genes", "performance metrics"] + sorted_short_metrics): return idx
        else: return str(idx)  + " (%i)"%(sum(r==1))

    df_square.index = [get_edited_index(idx, r) for idx,r in df_square.iterrows()]

    # get the sorted rows
    sorted_rows = list(df_square.index)

    # define col colors
    df_col_colors_all = pd.DataFrame()
    sID_to_color = get_value_to_color([0, 1, 2, 3, 4], palette="Greens", type_color="hex")[0]
    pID_to_color = get_value_to_color(list(range(max(df_ML_results_top[["phenotype", "test_parmID_exceptTestSS"]].drop_duplicates().groupby("phenotype").apply(len)))), palette="tab20", type_color="hex")[0]

    for phenotype, df_ML in df_ML_results_top.groupby(["phenotype"]):

        # get redefined vals
        all_test_parmID_exceptTestSS = sorted(set(df_ML.test_parmID_exceptTestSS))
        pID_to_numeric_ID = dict(zip(all_test_parmID_exceptTestSS, range(len(all_test_parmID_exceptTestSS))))

        df_col_colors = df_ML[["phenotype", "test_parmID_exceptTestSS", "test_split_ID"]].copy()
        df_col_colors["test parms"] = df_col_colors.test_parmID_exceptTestSS.map(pID_to_numeric_ID)
        df_col_colors["model ID"] = df_col_colors.test_split_ID.copy()

        # get the col ID
        df_col_colors["colID"] = df_col_colors.phenotype.apply(get_paper_dicotomic_pheno) + "|" + df_col_colors.test_parmID_exceptTestSS + "|" + df_col_colors.test_split_ID.apply(str)

        # redefine as colors
        for f, dict_colors in [("phenotype", ML_pheno_to_color), ("test parms", pID_to_color), ("model ID", sID_to_color)]:
            df_col_colors[f] = df_col_colors[f].map(dict_colors)

        df_col_colors_all = pd.concat([df_col_colors_all, df_col_colors[["phenotype", "test parms", "model ID", "colID"]].copy()])

    df_col_colors_all = df_col_colors_all.set_index("colID").loc[sorted_cols].copy()

    # get annot with the best models
    best_colIDs = set(df_ML_results_top[df_ML_results_top.best_model==True].colID)

    df_annot = pd.DataFrame(index=sorted_rows)
    for colID in sorted_cols:

        if colID in best_colIDs:
            col_series = df_square[colID].reset_index(drop=True)
            nmeric_rows_plot = set(col_series[col_series==1].index)
            df_annot[colID] = list(map(lambda r: {True:"*", False:""}[r in nmeric_rows_plot], range(len(sorted_rows))))

        else:
            df_annot[colID] = ""

    # redefine annot with metrics
    df_ML_results_top_annot = df_ML_results_top.rename(columns=metric_to_short_metric)[sorted_short_metrics + ["colID"]].copy()
    colID_to_short_metric_to_val = df_ML_results_top_annot.set_index("colID").apply(dict, axis=1)

    for colID in sorted_cols:
        for metric in sorted_short_metrics:

            def map_values_annot(x):
                if x==0: return "0"
                elif x==1: return "1"
                else: return ("%.2f"%x).split(".")[1]

            df_annot.loc[metric, colID] = map_values_annot(colID_to_short_metric_to_val[colID][metric])

    # plot
    cm = sns.clustermap(df_square, row_cluster=False, col_cluster=False, cmap=["white", "lightgrey", "black"], yticklabels=1, linewidth=.3, xticklabels=1, col_colors=df_col_colors_all, annot=df_annot, fmt="")

    # change texts
    for text in cm.ax_heatmap.texts:
        if text.get_text() == "*": text.set_fontsize(10)
        else: text.set_fontsize(6.5)

    # adjust
    adjust_cm_positions(cm, df_square.copy(), hm_height_multiplier=0.015, hm_width_multiplier=0.012, cc_height_multiplier=0.017, distance_btw_boxes=0.0025, cbar_width=0.02, cbar_height=0.15, ncol_fields=3)

    # add lines to sepparate phenotypes
    previous_pheno = None
    for Ic, colID in enumerate(sorted_cols):
        current_pheno = paper_pheno_to_pheno[colID.split("|")[0]]
        if not previous_pheno is None and previous_pheno!=current_pheno:
            #cm.ax_heatmap.axvline(Ic-0.1, c=ML_pheno_to_color[previous_pheno], ls="--", lw=.7)
            cm.ax_heatmap.axvline(Ic-0.06, c=ML_pheno_to_color[current_pheno], ls="--", lw=.9)

        previous_pheno = current_pheno


    # edit rows
    yticks = cm.ax_heatmap.get_yticks()

    # Specify which labels to make bold
    bold_ticks = {yticks[I] for I,label in enumerate(sorted_rows) if label in (fields_plot + ["feature genes", "performance metrics"])}

    # define the gene ticks
    #gene_ticks = {yticks[I] for I,label in enumerate(sorted_rows) if label in sorted_gene_names}

    # Create font properties
    from matplotlib.font_manager import FontProperties

    # Apply font properties individually
    bold_font = FontProperties(weight='bold')

    for I, (tick, label) in enumerate(zip(cm.ax_heatmap.get_yticks(), cm.ax_heatmap.get_yticklabels())):
        if tick in bold_ticks:
            label.set_fontproperties(bold_font)  # Apply bold font

        # change colors if all from the same phenotype
        r_I = df_square.iloc[I]
        phenotypes_1 = set(map(lambda x: paper_pheno_to_pheno[x.split("|")[0]], set(r_I[r_I==1].index)))

        if len(phenotypes_1)==1: color_set = ML_pheno_to_color[next(iter(phenotypes_1))]
        else: color_set = "black"
        label.set_color(color_set)  # Change color

    # edit xticks
    cm.ax_heatmap.set_xticklabels([x.split("|")[-1] for x in sorted_cols], rotation=0)
    save_fig(cm, "test_parms_picked.pdf", PlotsDir)



    # ['ASR_method', 'AUC_diff_threshold', 'CVset_model_consistency', 'CVset_type_split', 'aligner', 'all_affected_genes_str', 'all_affected_vars_str', 'class_weight', 'consistency_measure', 'expected_genes_str', 'feature_selection_method', 'features_max_pval_GWAS', 'gene_set', 'maxT_pval_consistency_mean_ROC_AUC', 'maxT_pval_consistency_min_ROC_AUC', 'mean_ROC_AUC_CVset_resamples', 'min_ROC_AUC_CVset_resamples', 'min_support', 'model_name', 'multiple_test_corr', 'prioritize_expected_genes', 'real_type_phenotype', 'sample_set', 'selected_predictors', 'selected_predictors_redundant', 'test_split_ID', 'time_spent_seconds_CVset', 'transition_definition', 'type_phenotype', 'test_sample_set', 'train_samples', 'test_samples', 'train_phenos', 'test_phenos', 'train_probabilities', 'test_probabilities', 'ignored_samples_train', 'ignored_samples_test', 'extra_info_train', 'extra_info_test', 'ROC_AUC_train', 'ROC_AUC_test', 'ROC_AUC_all', 'parameters_ID', 'cross_splits_phenotypes', 'cross_splits_samples', 'cross_splits_probabilities', 'cross_splits_ignored_samples', 'cross_splits_extra_info', 'cross_splits_ROC_AUC', 'redundant_samples_phenotypes', 'redundant_samples_samples', 'redundant_samples_probabilities', 'redundant_samples_ignored_samples', 'redundant_samples_ROC_AUC', 'redundant_samples_n_transitions', 'cross_splits_consistency_jaccard_variants', 'cross_splits_consistency_jaccard_samples', 'cross_splits_consistency_equal_variants', 'cross_splits_consistency_equal_samples', 'ROC_AUC_test_integrated_onlyCorrectModels', 'test_integrated_onlyCorrectModels_phenotypes', 'test_integrated_onlyCorrectModels_probabilities', 'test_integrated_onlyCorrectModels_samples', 'test_integrated_onlyCorrectModels_ignored_samples', 'ROC_AUC_test_integrated', 'test_integrated_phenotypes', 'test_integrated_probabilities', 'test_integrated_samples', 'test_integrated_ignored_samples', 'nmodels_fit', 'PRC_AUC_all', 'PRC_AUC_cross_splits', 'PRC_AUC_redundant_samples', 'PRC_AUC_test', 'PRC_AUC_test_integrated', 'PRC_AUC_test_integrated_onlyCorrectModels', 'PRC_AUC_train', 'test_parmID', 'n_train_parameter_combinations', 'cross_splits_consistency_jaccard_genes', 'test_parmID_exceptTestSS', 'phenotype', 'optimal_threshold_cross_splits', 'optimal_threshold_test_integrated_onlyCorrectModels', 'cross_splits_optimal_recall', 'cross_splits_optimal_precision', 'cross_splits_optimal_fpr', 'cross_splits_optimal_neg_predictive_val', 'cross_splits_optimal_F_score_precision_recall', 'redundant_samples_optimal_recall', 'redundant_samples_optimal_precision', 'redundant_samples_optimal_fpr', 'redundant_samples_optimal_neg_predictive_val', 'redundant_samples_optimal_F_score_precision_recall', 'test_integrated_onlyCorrectModels_optimal_recall', 'test_integrated_onlyCorrectModels_optimal_precision', 'test_integrated_onlyCorrectModels_optimal_fpr', 'test_integrated_onlyCorrectModels_optimal_neg_predictive_val', 'test_integrated_onlyCorrectModels_optimal_F_score_precision_recall', 'validation_transitions_TrainTest_ROC_AUC', 'validation_transitions_TrainTest_optimalThreshold_recall', 'validation_transitions_TrainTest_optimalThreshold_precision', 'validation_transitions_TrainTest_optimalThreshold_fpr', 'validation_transitions_TrainTest_optimalThreshold_neg_predictive_val', 'validation_transitions_TrainTest_optimalThreshold_F_score_precision_recall', 'validation_transitions_TrainTest_PRC_AUC', 'validation_transitions_Test_ROC_AUC', 'validation_transitions_Test_optimalThreshold_recall', 'validation_transitions_Test_optimalThreshold_precision', 'validation_transitions_Test_optimalThreshold_fpr', 'validation_transitions_Test_optimalThreshold_neg_predictive_val', 'validation_transitions_Test_optimalThreshold_F_score_precision_recall', 'validation_transitions_Test_PRC_AUC', 'best_model']


def plot_test_parameter_performance_top_models(ML_df_dicts, PlotsDir, test_sample_set, threshold_ROC_AUC, threshold_consistency_genes):

    """Plots test parameter performances"""

    # define the test fields to check out
    fields_optimize_test = ["aligner", "CVset_type_split", "consistency_measure", "AUC_diff_threshold", "sample_set", "gene_set"]

    # get df
    df_ML_results_top = get_filtered_df_by_dict(ML_df_dicts["df_ML_results_top"], {"test_sample_set":test_sample_set})[["test_parmID_exceptTestSS", "phenotype", "best_model"] + fields_optimize_test].drop_duplicates()

    # prints
    #print(df_ML_results_top[["consistency_measure", "AUC_diff_threshold"]].drop_duplicates().sort_values(by=["consistency_measure", "AUC_diff_threshold"]))

    # checks
    if len(df_ML_results_top)!=len(df_ML_results_top[["test_parmID_exceptTestSS", "phenotype"]].drop_duplicates()):
        raise ValueError("non unique vals")

    # get the per-model performance with all necessary phenotypes
    parm_vals = {f : set(df_ML_results_top[f]) for f in fields_optimize_test}
    parm_vals["gene_set"].update({"all", "known"}) # add a control
    parm_vals["phenotype"] = set(df_ML_results_top.phenotype)
    parm_vals["test_sample_set"] = test_sample_set
    parm_vals["multiple_test_corr"] = "none" # we discarded this

    df_ML_performance = get_filtered_df_by_dict(ML_df_dicts["df_ML_integrated_performance"], parm_vals)[fields_optimize_test + ["phenotype", "AUC_test", "mean_cross_splits_consistency_jaccard_genes", "AUC_train", "AUC_train_test_diff"]].copy()

    # discard some
    df_ML_performance = df_ML_performance[~((df_ML_performance.consistency_measure=="none") & (df_ML_performance.AUC_diff_threshold!=0.05))].copy()

    # merge is_top
    df_ML_results_top["is_top"] = True
    df_ML_performance = df_ML_performance.merge(df_ML_results_top[fields_optimize_test + ["phenotype", "is_top", "best_model"]], on=fields_optimize_test +["phenotype"], validate="one_to_one", how="left")

    # redefine
    def get_minus1_to_n(x, val):
        if x==-1: return val
        else: return x

    df_ML_performance["AUC_test"] = df_ML_performance["AUC_test"].apply(get_minus1_to_n, val=0.5)
    min_consistency = min(df_ML_performance[df_ML_performance["mean_cross_splits_consistency_jaccard_genes"]!=-1].mean_cross_splits_consistency_jaccard_genes)
    max_consistency = max(df_ML_performance.mean_cross_splits_consistency_jaccard_genes)
    df_ML_performance["mean_cross_splits_consistency_jaccard_genes"] = df_ML_performance["mean_cross_splits_consistency_jaccard_genes"].apply(get_minus1_to_n, val=min_consistency)

    if len(df_ML_performance)!=len(df_ML_performance[fields_optimize_test +["phenotype"]].drop_duplicates()):
        raise ValueError("non unique vals df_ML_performance")

    # add merged fields
    df_ML_performance["type split & samples"] = "split="+df_ML_performance.CVset_type_split + ", samples=" + df_ML_performance.sample_set.apply(lambda x: x.split("_")[0])

    df_ML_performance["consistency"] = df_ML_performance.consistency_measure.apply(lambda x: x.replace("jaccard_", ""))
    df_ML_performance["AUC diff"] = df_ML_performance.AUC_diff_threshold.apply(lambda x: "."+str(x).split(".")[1])
    df_ML_performance["consistency & AUC diff"] = "consistency="+df_ML_performance.consistency + ", AUC diff=" + df_ML_performance["AUC diff"]
    df_ML_performance["aligner"] = df_ML_performance.aligner.apply(lambda x: x.split("_")[0])

    # define sorted vals
    sorted_cons_AUC = ['consistency=none, AUC diff=.05', 'consistency=genes, AUC diff=.05', 'consistency=genes, AUC diff=.1', 'consistency=genes, AUC diff=.15', 'consistency=samples, AUC diff=.05', 'consistency=samples, AUC diff=.1', 'consistency=samples, AUC diff=.15', 'consistency=variants, AUC diff=.05', 'consistency=variants, AUC diff=.1', 'consistency=variants, AUC diff=.15']

    """
    markers_dict = {
        'consistency=none, AUC diff=.05': 'o',
        
        'consistency=genes, AUC diff=.05': 'v', 
        'consistency=genes, AUC diff=.1': '^',  
        'consistency=genes, AUC diff=.15': '<', 
        
        'consistency=samples, AUC diff=.05': 's', 
        'consistency=samples, AUC diff=.1': 'p', 
        'consistency=samples, AUC diff=.15': 'h',

        'consistency=variants, AUC diff=.05': 'P',
        'consistency=variants, AUC diff=.1': 'D',  
        'consistency=variants, AUC diff=.15': 'X', 
    }
    """
    consistency_to_marker = {"none":"P", "genes":"o", "variants":"v", "samples":"s"}
    #AUC_diff_to_marker = {".05":"", ".1":"|", ".15":"--"}
    AUC_diff_to_marker = {".05":"", ".1":"", ".15":""}

    gene_set_to_color = {"all":"gray", "known":"gold", 'known+by_description':"red", "known+by_description+selection_all":"salmon", 'known+by_description+selection_SNP':"magenta", "by_description":"teal", "by_description+selection_SNP":"c", 'known+by_description+GWAS_high+selection_SNP':"black"}

    f_to_sorted_vals = {"aligner": ["bwa", "bowtie2", "hisat2"],
                        "consistency": ["none", "genes", "samples", "variants"],
                        "type split & samples": ['split=random, samples=all', 'split=random, samples=representative', 'split=balanced, samples=all', 'split=balanced, samples=representative'],
                        "consistency & AUC diff": sorted_cons_AUC,
                        "gene_set": list(gene_set_to_color.keys()),
                        "AUC diff": [".05", ".1", ".15"]}

    # checks
    for f, sorted_vals in f_to_sorted_vals.items():
        if set(df_ML_performance[f])!=set(sorted_vals): raise ValueError("should be the same")

    # add idx fields
    for f, sorted_vals in f_to_sorted_vals.items():
        df_ML_performance[f+"_idx"] = df_ML_performance[f].map(dict(zip(sorted_vals, range(len(sorted_vals)))))

    df_ML_performance = df_ML_performance.sort_values(by=["gene_set_idx", "consistency & AUC diff_idx"], ascending=[False, False])

    # add valid train
    df_ML_performance["valid_train"] = (df_ML_performance.AUC_train>=threshold_ROC_AUC) & (df_ML_performance.AUC_train_test_diff.apply(abs)<=0.2)

    # get figure each parameter changing
    get_figures_each_test_parameter_changing(df_ML_performance, PlotsDir, gene_set_to_color, f_to_sorted_vals, min_consistency, max_consistency, threshold_consistency_genes, threshold_ROC_AUC)

    ##### FIGURE #####

    # one fig each phenotype
    for phenotype in sorted(set(df_ML_performance.phenotype)):
        print(phenotype)

        # get df
        df_ML_performance_pheno = df_ML_performance[df_ML_performance.phenotype==phenotype].copy().set_index(["aligner", "type split & samples"])
        top_min_AUC = min(df_ML_performance_pheno[df_ML_performance_pheno.is_top==True].AUC_test)
        top_min_consistency = min(df_ML_performance_pheno[df_ML_performance_pheno.is_top==True].mean_cross_splits_consistency_jaccard_genes)

        # init fig
        sorted_rows = f_to_sorted_vals["aligner"]
        sorted_cols = f_to_sorted_vals["type split & samples"]
        nrows = len(sorted_rows)
        ncols = len(sorted_cols)
        fig = plt.figure(figsize=(ncols*1.3, nrows*1.3)); Ip=1
        #fig = plt.figure(figsize=(ncols*3, nrows*3)); Ip=1

        for Ir, row_val in enumerate(sorted_rows):
            for Ic, col_val in enumerate(sorted_cols):

                # get plot
                ax = plt.subplot(nrows, ncols, Ip); Ip+=1

                # get df
                df_plot = df_ML_performance_pheno.loc[(row_val, col_val)].copy().sort_values(by=["gene_set_idx", "consistency & AUC diff_idx"], ascending=[False, False])
                if len(df_plot)!=len(df_plot[["gene_set", "consistency & AUC diff"]].drop_duplicates()): raise ValueError("should be unique")

                # filter out to only get the train ones
                #df_plot = df_plot[(df_plot.valid_train==True)].copy()

                # get one plot for each type                
                ax = sns.scatterplot(data=df_plot, x="mean_cross_splits_consistency_jaccard_genes", y="AUC_test", hue="gene_set", style="consistency", s=20, edgecolor="black", linewidth=.3, palette=gene_set_to_color, markers=consistency_to_marker)

                # add text with the consistecnies
                for I,r in df_plot[df_plot["AUC diff"]!=".05"].iterrows():
                    ax.text(r.mean_cross_splits_consistency_jaccard_genes, r.AUC_test, AUC_diff_to_marker[r["AUC diff"]], fontsize=12, color=gene_set_to_color[r.gene_set], va="center", ha="center", zorder=2)

                # legend
                #if Ir==0 and Ic==(len(sorted_cols)-1): ax.legend(bbox_to_anchor=(1,1), markerscale=1.5)
                if len(df_plot)>0: ax.get_legend().remove()

                def get_legend_element(color, marker, label, markeredgecolor="black"): return plt.Line2D([0], [0], marker=marker, color=color, markerfacecolor=color, markeredgecolor=markeredgecolor, markersize=5, label=label, linewidth=0, markeredgewidth=.3)

                if Ir==0 and Ic==(len(sorted_cols)-1):

                    legend_handles = [get_legend_element("white", "o", "gene set", markeredgecolor="white")] + [get_legend_element(gene_set_to_color[gene_set], "o", gene_set) for gene_set in f_to_sorted_vals["gene_set"]] + [get_legend_element("white", "o", "\nconsistency measure", markeredgecolor="white")] + [get_legend_element("black", consistency_to_marker[cons], cons) for cons in f_to_sorted_vals["consistency"]]
                    
                    ax.legend(handles=legend_handles, title="", fontsize=10, title_fontsize=10, bbox_to_anchor=(1.15,1), markerscale=1.5)

                # markers and lims
                ax.set_xlim([min_consistency-0.05, max_consistency+0.05])
                ax.set_ylim([0.45, 1.05])

                if Ic==0 and Ir==1: ax.set_ylabel("test $AUC_{ROC}$")
                else: ax.set_ylabel("")
                if Ic!=0: ax.set_yticklabels([])

                if Ir==2 and Ic==1: ax.set_xlabel("                       m1-m5 feature consistency (genes)")
                else: ax.set_xlabel("")
                if Ir!=2: ax.set_xticklabels([])

                # add titles
                if Ir==0:
                    ax.set_title("\n".join([x.split("=")[1] + " " + x.split("=")[0] for x in col_val.split(", ")]).replace("representative", "rep.").replace("balanced", "bal.").replace("random", "rand."), fontsize=10, ha="center")

                if Ic==3:
                    ax.text(1.1, 0.75, row_val, fontsize=10, rotation=90, va="center")

                if Ir==0 and Ic==1:
                    ax.text(1.05, 1.35, "performance test parameters '"+get_paper_dicotomic_pheno(phenotype) +"'", ha="center", fontsize=10)

                # add lines
                plt.axvline(threshold_consistency_genes, lw=.8, c="gray", ls="--")
                plt.axhline(threshold_ROC_AUC, lw=.8, c="gray", ls="--")

                df_plot_top = df_plot[df_plot.is_top==True].copy()
                if len(df_plot_top)>0:

                    plt.plot([top_min_consistency, top_min_consistency], [top_min_AUC, 1.05], lw=.8, c="black", ls="--") 
                    plt.plot([top_min_consistency, max_consistency+0.05], [top_min_AUC, top_min_AUC], lw=.8, c="black", ls="--") 

                    n_gs = len(set(df_plot_top.gene_set))
                    for Igs, (gs, df_gs) in  enumerate(df_plot_top.groupby("gene_set")):
                        x_text = 0.85+Igs*0.1
                        ax.text(x_text, 0.92, str(len(df_gs)), ha="right", fontsize=10, color=gene_set_to_color[gs])
                        if (Igs+1)<n_gs: ax.text(x_text+0.04, 0.92, "/", ha="right", fontsize=10, color="gray")

        # get fig
        sns.despine()
        plt.subplots_adjust(wspace=0.05, hspace=0.05)# wspace (between cols) hspaces (btw rows) 

        save_fig(fig, "test_parameters_performance/perfromance_%s.pdf"%phenotype, PlotsDir)
        plt.show()

        ##################


def get_df_GWAS_parameter_benchmark_filtered(ProcessedDataDir):

    """filter some parms"""

    df_GWAS_parameter_benchmark_filt_file = "%s/filtered_GWAS_parameter_benchmarking_df.py"%ProcessedDataDir
    if file_is_empty(df_GWAS_parameter_benchmark_filt_file):

        # load raw df
        print("loading raw df")
        df_GWAS_parameter_benchmark = load_object(ProcessedDataDir + "/GWAS_benchmarking_df_all_some_phenos.py")
        print("\n".join([get_paper_dicotomic_pheno(p) for p in sorted(set(df_GWAS_parameter_benchmark.phenotype))]))

        # keep the lab measured phenos
        #lab_measured_phenos = {p for p in set(df_GWAS_parameter_benchmark.phenotype) if not p.split("_")[0] in {"age", "source"}}
        #df_GWAS_parameter_benchmark = df_GWAS_parameter_benchmark[df_GWAS_parameter_benchmark.phenotype.isin(lab_measured_phenos)].copy()

        # keep some fields -> the target genes are because I find them
        all_expected_genes = sorted(set(make_flat_listOflists(get_literature_expected_genes_drug_types().values())))

        print("filtering...")
        target_fields = ['pval_chi_square_phenotypes', 'pval_GenoAndPheno_phenotypes', 'pval_fisher', 'pval_chi_square_phenotypes_maxT', 'pval_epsilon_phenotypes_maxT', 'correction_method', 'correction_all_groups', 'min_epsilon', 'min_nodes_GenoAndPheno', 'filter_I', 'n_sig_variants', 'n_sig_variants_single_clade', 'n_sig_variants_no_multi_clade', 'n_sig_genes', 'n_sig_genes_multi_clade', 'n_expected_genes_found', 'any_pval', 'any_pval_correction', 'ASR_method', 'min_support', 'sample_set', 'aligner', 'phenotype', 'unique_filterID', 'set_expected_genes_found', "max_n_sig_genes", "max_n_sig_variants_no_multi_clade", "fraction_transitions_explained", "pvalue_threshold", "min_fraction_transitions_explained"] + ["%s_affected_nodes"%g for g in all_expected_genes] + ["%s_found"%g for g in all_expected_genes]

        # 'max_n_sig_genes', 'min_n_sig_genes', 'max_n_sig_genes_multi_clade', 'min_n_sig_genes_multi_clade', 'max_n_sig_variants', 'min_n_sig_variants', 'max_n_sig_variants_single_clade', 'min_n_sig_variants_single_clade', 'max_n_sig_variants_no_multi_clade', 'min_n_sig_variants_no_multi_clade']
        df_GWAS_parameter_benchmark = df_GWAS_parameter_benchmark[target_fields].reset_index(drop=True).copy()

        # rename fields
        pval_to_renamed_pval = {'pval_chi_square_phenotypes':"p_X2", 'pval_GenoAndPheno_phenotypes':"p_GP", 'pval_chi_square_phenotypes_maxT':"p_X2_maxT", 'pval_epsilon_phenotypes_maxT':"p_e_maxT"}
        df_GWAS_parameter_benchmark = df_GWAS_parameter_benchmark.rename(columns=pval_to_renamed_pval)

        # tests
        if len(df_GWAS_parameter_benchmark)!=len(df_GWAS_parameter_benchmark[["phenotype", "unique_filterID"]].drop_duplicates()):
            raise ValueError("non unique filters")
        print("There are %i parms initially"%(len(set(df_GWAS_parameter_benchmark.unique_filterID))))

        # general filters
        df_GWAS_parameter_benchmark = get_filtered_df_by_dict(df_GWAS_parameter_benchmark, {"pval_fisher":False, "min_epsilon":{0.1, 0.2, 0.3}, "correction_method":{"none", "bonferroni", "fdr_bh"}})

        # filter types of pval
        df_GWAS_parameter_benchmark["any_maxT"] = df_GWAS_parameter_benchmark[['p_X2_maxT', 'p_e_maxT']].apply(any, axis=1)
        df_GWAS_parameter_benchmark["n_pvals"] = df_GWAS_parameter_benchmark[['p_X2', 'p_GP', 'p_X2_maxT', 'p_e_maxT']].apply(sum, axis=1)
        df_GWAS_parameter_benchmark["no_pval"] = (df_GWAS_parameter_benchmark.n_pvals==0)
        df_GWAS_parameter_benchmark["p_uncorrected"] = (df_GWAS_parameter_benchmark.correction_method=="none") & (df_GWAS_parameter_benchmark[["p_X2", "p_GP"]].apply(any, axis=1)) & ~(df_GWAS_parameter_benchmark.any_maxT)

        df_GWAS_parameter_benchmark["only_p_X2_maxT"] = (df_GWAS_parameter_benchmark.n_pvals==1) & (df_GWAS_parameter_benchmark.p_X2_maxT)
        df_GWAS_parameter_benchmark["only_p_e_maxT"] = (df_GWAS_parameter_benchmark.n_pvals==1) & (df_GWAS_parameter_benchmark.p_e_maxT)
        df_GWAS_parameter_benchmark["both_maxT"] = (df_GWAS_parameter_benchmark.n_pvals==2) & (df_GWAS_parameter_benchmark[["p_X2_maxT", "p_e_maxT"]].apply(all, axis=1))

        for corr in ["bonferroni", "fdr_bh"]:
            corr_name = {"bonferroni":"bonf", "fdr_bh":"fdr"}[corr]
            df_GWAS_parameter_benchmark["only_p_X2_"+corr_name] = (df_GWAS_parameter_benchmark.correction_method==corr) & (df_GWAS_parameter_benchmark.n_pvals==1) & (df_GWAS_parameter_benchmark.p_X2)
            df_GWAS_parameter_benchmark["only_p_GP_"+corr_name] = (df_GWAS_parameter_benchmark.correction_method==corr) & (df_GWAS_parameter_benchmark.n_pvals==1) & (df_GWAS_parameter_benchmark.p_GP)
            df_GWAS_parameter_benchmark["p_X2_and_GP_"+corr_name] = (df_GWAS_parameter_benchmark.correction_method==corr) & (df_GWAS_parameter_benchmark.n_pvals==2) & (df_GWAS_parameter_benchmark[["p_X2", "p_GP"]].apply(all, axis=1))

        for pvalue_threshold in [0.05, 0.001, 0.0001]:
            idx_pval = (df_GWAS_parameter_benchmark.pvalue_threshold==pvalue_threshold) & (df_GWAS_parameter_benchmark.correction_method=="none")
            df_GWAS_parameter_benchmark["only_p_X2<%s"%pvalue_threshold] = idx_pval & (df_GWAS_parameter_benchmark.n_pvals==1) & (df_GWAS_parameter_benchmark.p_X2)
            df_GWAS_parameter_benchmark["only_p_GP<%s"%pvalue_threshold] = idx_pval & (df_GWAS_parameter_benchmark.n_pvals==1) & (df_GWAS_parameter_benchmark.p_GP)
            df_GWAS_parameter_benchmark["p_X2_and_GP<%s"%pvalue_threshold] = idx_pval & (df_GWAS_parameter_benchmark.n_pvals==2) & (df_GWAS_parameter_benchmark[["p_X2", "p_GP"]].apply(all, axis=1))

        # filter
        pval_fields_unique = ['no_pval', 'only_p_X2_maxT', 'only_p_e_maxT', 'both_maxT', 'only_p_X2_bonf', 'only_p_GP_bonf', 'p_X2_and_GP_bonf', 'only_p_X2_fdr', 'only_p_GP_fdr', 'p_X2_and_GP_fdr'] + ["%s<%s"%(f, pvalue_threshold) for pvalue_threshold in [0.05, 0.001, 0.0001] for f in ["only_p_X2", "only_p_GP", "p_X2_and_GP"]]
        df_GWAS_parameter_benchmark = df_GWAS_parameter_benchmark[df_GWAS_parameter_benchmark[pval_fields_unique].apply(any, axis=1)].copy()
        if any(df_GWAS_parameter_benchmark[pval_fields_unique].apply(sum, axis=1)!=1): raise ValueError("should only be one type")

        # add the general p value type
        print("adding type pval general")
        def get_type_pval_general(r): return r[r].index[0]
        df_GWAS_parameter_benchmark["type pvalues general"] = df_GWAS_parameter_benchmark[pval_fields_unique].copy().apply(get_type_pval_general, axis=1)

        # print # parameters
        print("There are %i total parms"%(len(set(df_GWAS_parameter_benchmark.unique_filterID))))
        print("# parms per combination", set(df_GWAS_parameter_benchmark.groupby(["phenotype", "aligner", "min_support", "ASR_method", "sample_set"]).apply(len).values))

        # save
        save_object(df_GWAS_parameter_benchmark, df_GWAS_parameter_benchmark_filt_file)

    print("load benchmarking parms")
    return load_object(df_GWAS_parameter_benchmark_filt_file)


def report_expected_genes_all_phenos(df_GWAS_parameter_benchmark):

    """report expected phenos"""

    # explore expected gene found in different phenotypes
    type_to_genes = {"azole": ["ERG11", "ERG3", "Scer_HAL9", "MRR1", "NDT80", "UPC2"],
                     "echinocandin": ["FKS1", "FKS2", "ERG3"],
                     "AMB": ["ERG2", "ERG3", "Scer_ERG5", "ERG6", "ERG25"],
                     "5FC": ["Scer_FCY1", "Scer_FCY2", "Scer_FUR1", "Scer_CHS7"]}

    resistance_phenos = [p for p in sorted(set(df_GWAS_parameter_benchmark.phenotype)) if p.endswith("_resistance")]

    for pheno in resistance_phenos:
        drug = pheno.split("_")[0]
        if drug in {"VRC", "FLC", "POS"}: target_genes = type_to_genes["azole"]
        elif drug in {"MIF", "ANI", "CAS"}: target_genes = type_to_genes["echinocandin"]
        elif drug in type_to_genes: target_genes = type_to_genes[drug]
        else: raise ValueError("invalid %s"%drug)

        df_pheno = df_GWAS_parameter_benchmark[df_GWAS_parameter_benchmark.phenotype==pheno].copy()

        for gene in target_genes:
            df_pheno_g = df_pheno[df_pheno.set_expected_genes_found.apply(lambda x: gene in x)].copy()
            if len(df_pheno_g)>0:
                print(pheno, "# parms all with", gene, len(df_pheno_g))
                print(pheno, "# parms some pval with", gene, sum(df_pheno_g.any_pval))
                print(pheno, "# parms some p corrected with", gene, sum(df_pheno_g.any_pval_correction))
                print(pheno, "# parms some p corrected, min_support>10", gene, sum((df_pheno_g.any_pval_correction) & (df_pheno_g.min_support>10)))

        print()

def get_low_confidence_parameters_GWAS(df_GWAS_parameter_benchmark, PlotsDir):

    """explore to get low conf parms"""

    # explore parameters that yield all of the expected genes
    print("exploring low confidence parms")
    pheno_to_expected_genes = {"FLC_resistance":{"MRR1", "NDT80", "ERG11"},
                               "VRC_resistance":{"NDT80", "MRR1", "ERG11"}, # to find MRR1 & ERG11 I'd need min support=10
                               "POS_resistance":{"ERG11"}, # no TAC1, which led to no parms
                               "CAS_resistance":{"FKS2"}
                               }

    # keep phenos
    df_parameters = df_GWAS_parameter_benchmark[df_GWAS_parameter_benchmark.phenotype.isin(set(pheno_to_expected_genes.keys()))].copy().reset_index(drop=True)

    # check
    strange_genes = set.union(*pheno_to_expected_genes.values()).difference(set.union(*df_parameters.set_expected_genes_found))
    if len(strange_genes)>0: raise ValueError("invalid %s"%strange_genes)

    # keep parms that have all expected genes
    df_parameters["real_expected_genes_set"] = df_parameters.phenotype.map(pheno_to_expected_genes)

    def get_real_expected_genes_found(r):
        intersection_genes = r.set_expected_genes_found.intersection(r.real_expected_genes_set)
        return intersection_genes==r.real_expected_genes_set

    df_parameters["real_expected_genes_found"] = df_parameters[["set_expected_genes_found", "real_expected_genes_set"]].copy().apply(get_real_expected_genes_found, axis=1)
    parmID_to_all_expected_genes_found = df_parameters[["unique_filterID", "real_expected_genes_found"]].copy().groupby(["unique_filterID"]).apply(lambda df: all(df.real_expected_genes_found))
    df_parameters = df_parameters[df_parameters.unique_filterID.isin(set(parmID_to_all_expected_genes_found[parmID_to_all_expected_genes_found].index))].copy()


    # pre-expected filters
    #df_parameters = df_parameters[(df_parameters.any_pval_correction) | (df_parameters.min_support>10)].copy()
    #df_parameters = df_parameters[df_parameters.max_n_sig_variants_no_multi_clade<=150].copy()

    # define parm fields
    parm_fields = ["max_n_sig_variants_no_multi_clade", "min_fraction_transitions_explained", "min_support", "type pvalues general", "aligner", "ASR_method", "sample_set", "min_epsilon", "min_nodes_GenoAndPheno", "correction_all_groups"]
    ascending_fields = [True, False, False, True, True, True, True, True, True, True]

    # additional filters
    df_parms_low_conf = df_parameters[parm_fields].drop_duplicates().sort_values(by=parm_fields, ascending=ascending_fields).reset_index(drop=True)

    # write
    df_parms_low_conf.to_excel(PlotsDir + "/GWAS_parms_yielding_all_expected_genes.xlsx")

    # add the manually curated ones
    filters_dict = {k:v for k,v in manual_curation_gwas_filt_parms_low_conf.items()}
    df_parms_low_conf_best = get_filtered_df_by_dict(df_parms_low_conf, filters_dict)
    if len(df_parms_low_conf_best)!=1: raise ValueError("should be 1")
    print("best low conf parms:\n", df_parms_low_conf_best.iloc[0])


    # return
    return df_parms_low_conf


def explore_GWAS_parameters_plots(df_GWAS_parameter_benchmark, PlotsDir, gwas_filt_parms):

    """plots different benchmarking things"""

    # define general vars
    sorted_phenotypes = ["FLC_resistance", "VRC_resistance", "POS_resistance", "CAS_resistance"]
    sorted_expected_genes_all = ["MRR1", "NDT80", "ERG11", "Scer_HAL9", "FKS2"]

    df_GWAS_parameter_benchmark = df_GWAS_parameter_benchmark[df_GWAS_parameter_benchmark.phenotype.isin(set(sorted_phenotypes))].copy()

    # add fields
    df_GWAS_parameter_benchmark["correction_scope"] = df_GWAS_parameter_benchmark.correction_all_groups.map({True:"conservative", False:"relaxed"})

    # get the filter ID
    best_filterID_df = get_filtered_df_by_dict(df_GWAS_parameter_benchmark, gwas_filt_parms)
    best_filterID_df = best_filterID_df[best_filterID_df.phenotype=="FLC_resistance"].copy()
    if len(best_filterID_df)!=1: raise ValueError("should be 1")
    best_filterID = best_filterID_df.iloc[0].unique_filterID

    ###### GRAPHICAL THINGS #######

    palette_aligner = {"bwa_mem":"tab:red", "hisat2_no_spliced":"gray", "bowtie2_local":"tab:blue",}
    palette_min_support = {10:"gray", 50:"tab:cyan", 70:"navy", 90:"tab:red"}
    palette_ASR = {"MPPA":"tab:blue", "DOWNPASS":"red", "MPPA,DOWNPASS":"black"}

    palette_ASR = {"MPPA":"tab:blue", "MPPA,DOWNPASS":"black", "DOWNPASS":"red"}


    palette_sample_set = {"all_samples":"gray", "balanced_samples":"tab:red", "representative_samples":"tab:blue"}
    palette_min_epsilon = {0.1:"tab:red", 0.2:"tab:blue", 0.3: "tab:gray"}
    palette_min_GP = {2:"gray", 3:"red"}

    palette_type_p_general = {'no_pval':"gray",
                              
                              'only_p_X2<0.05': 'tan', 'only_p_X2<0.001': 'tan', 'only_p_X2<0.0001': 'tan',
                              'only_p_GP<0.05': 'sienna', 'only_p_GP<0.001': 'sienna', 'only_p_GP<0.0001': 'sienna',
                              'p_X2_and_GP<0.05': 'black', 'p_X2_and_GP<0.001': 'black', 'p_X2_and_GP<0.0001': 'black',  

                              'only_p_X2_fdr':"olive", 'only_p_GP_fdr':"lime", 'p_X2_and_GP_fdr':"green",
                              'only_p_X2_bonf':"deepskyblue", 'only_p_GP_bonf':"royalblue", 'p_X2_and_GP_bonf':"blue",
                               'only_p_e_maxT':"tab:pink", 'only_p_X2_maxT':"m", 'both_maxT':"red"}

    palette_correction_scope = {"relaxed":"gray", "conservative":"tab:red"}

    # map palettes
    field_to_palette = {"aligner":palette_aligner,
                        "min_support":palette_min_support,
                        "ASR_method":palette_ASR,
                        "sample_set":palette_sample_set,
                        "min_epsilon":palette_min_epsilon,
                        "min_nodes_GenoAndPheno":palette_min_GP,
                        "type pvalues general":palette_type_p_general,
                        "correction_scope":palette_correction_scope}

    if len(df_GWAS_parameter_benchmark)!=len(df_GWAS_parameter_benchmark[["phenotype", 'p_X2', 'p_GP'] + list(field_to_palette.keys())].drop_duplicates()):
        raise ValueError("should be unique")

    ###############################

    ###### SEQUENTIAL FILTERING PLOTS #######

    # define function
    def run_plotting_parms(df_p, field, current_parms):
        plot_GWAS_parameter_effects_scatter(df_p.copy(), sorted_phenotypes, PlotsDir, field, field_to_palette, current_parms, best_filterID, sorted_expected_genes_all)


    # initial plot with all values, but checking effect on all parms
    print("plotting all parameters")
    df_plot = df_GWAS_parameter_benchmark.copy()
    run_plotting_parms(df_plot, 'type pvalues general', "all")

    # keep some p val correction 
    df_plot = df_plot[~(df_plot.p_uncorrected) & (df_plot.any_pval)].copy(); current_filters = "some p val correction"
    run_plotting_parms(df_plot, 'min_support', current_filters)

    # keep high support
    df_plot = df_plot[df_plot.min_support>=70].copy(); current_filters += ", min support>=70"
    run_plotting_parms(df_plot, 'ASR_method', current_filters)

    # keep non MPPA filters
    df_plot = df_plot[df_plot.ASR_method!="MPPA"].copy(); current_filters += ", DOWNPASS or MPPA-DOWNPASS ASR"
    run_plotting_parms(df_plot, 'min_epsilon', current_filters)

    # keep min_epsilon
    df_plot = df_plot[df_plot.min_epsilon==0.3].copy(); current_filters += ", min e=0.3"
    run_plotting_parms(df_plot, 'sample_set', current_filters)
    run_plotting_parms(df_plot, 'min_nodes_GenoAndPheno', current_filters)

    # keep certain values
    df_plot = df_plot[(df_plot.sample_set=="balanced_samples") & (df_plot.min_nodes_GenoAndPheno==2)].copy(); current_filters += "\nmin n=2, balanced samples"
    run_plotting_parms(df_plot, 'aligner', current_filters)

    # keep hisat2
    df_plot = df_plot[df_plot.aligner=="hisat2_no_spliced"].copy(); current_filters += ", hisat2"

    #########################################

    ##### GET FINAL FILTERS ######

    print("There are %i remaining parameters"%len(set(df_plot.unique_filterID)))

    # keep parameters that yield all expected phenos
    pheno_to_expected_genes = {"POS_resistance":{"ERG11"}, "FLC_resistance":{"ERG11", "NDT80"}}
    df_parameters = df_plot[df_plot.phenotype.isin(set(pheno_to_expected_genes))].copy()
    df_parameters["real_expected_genes_set"] = df_parameters.phenotype.map(pheno_to_expected_genes)

    def get_real_expected_genes_found(r):
        intersection_genes = r.set_expected_genes_found.intersection(r.real_expected_genes_set)
        return intersection_genes==r.real_expected_genes_set

    df_parameters["real_expected_genes_found"] = df_parameters[["set_expected_genes_found", "real_expected_genes_set"]].copy().apply(get_real_expected_genes_found, axis=1)
    parmID_to_all_expected_genes_found = df_parameters[["unique_filterID", "real_expected_genes_found"]].copy().groupby(["unique_filterID"]).apply(lambda df: all(df.real_expected_genes_found))
    interesting_parmIDs = set(parmID_to_all_expected_genes_found[parmID_to_all_expected_genes_found].index)
    df_parameters = df_parameters[df_parameters.unique_filterID.isin(interesting_parmIDs)].copy()
    if len(df_parameters)==0: raise ValueError("invalid")

    # print
    parm_fields = ["min_support", "ASR_method", "type pvalues general", "correction_scope", "sample_set", "min_epsilon", "min_nodes_GenoAndPheno", "aligner"]
    for f in parm_fields:
        print(f, sorted(set(df_parameters[f])))

    # get final parms 
    df_plot = df_plot[df_plot.unique_filterID.isin(interesting_parmIDs)].copy()

    final_fields_print = ["phenotype"] + parm_fields + ["fraction_transitions_explained", "n_sig_variants_no_multi_clade", "max_n_sig_variants_no_multi_clade", "n_sig_genes", "max_n_sig_genes"]
    df_plot[final_fields_print].sort_values(by=final_fields_print).reset_index(drop=True).to_excel(PlotsDir + "/GWAS_filters_high_confidence_final.xlsx")

    ##############################

def plot_GWAS_parameter_effects_scatter(df_plot, sorted_phenotypes, PlotsDir, hue_f, field_to_palette, current_parms,  best_filterID, sorted_expected_genes):

    """plots scatter parameter effects"""

    # get palette
    palette_dict = field_to_palette[hue_f]

    # get title
    title_str = "parameters: %s"%current_parms
    type_filt = current_parms.replace(" ", "_")

    # keep some
    sorted_expected_genes = [g for g in sorted_expected_genes  if any(df_plot[g+"_found"])]
    expected_genes_fields = [g+"_found" for g in sorted_expected_genes]
    df_plot = df_plot[["phenotype", hue_f, "n_sig_genes_multi_clade", "n_sig_variants_single_clade", "n_sig_genes", "fraction_transitions_explained", "unique_filterID", "n_sig_variants_no_multi_clade"] + expected_genes_fields].copy()
    df_plot = get_filtered_df_by_dict(df_plot, {"phenotype":set(sorted_phenotypes)})
    check_no_nans_in_df(df_plot)

    # sort bu hue
    sorted_hues = list(palette_dict.keys())
    df_plot["hue_idx"] = df_plot[hue_f].map(dict(zip(sorted_hues, range(len(sorted_hues)))))
    df_plot = df_plot.set_index("phenotype")

    # get as log scale
    df_plot["log_n_sig_variants_no_multi_clade"] = np.log10(df_plot.n_sig_variants_no_multi_clade + 1)

    # define the lims
    xfield = "n_sig_variants_no_multi_clade"; xlines_x = [10] # n_sig_variants_no_multi_clade
    yfield = "fraction_transitions_explained"; ylines_y = [0.5]

    xlims = [min(df_plot[xfield]), max(df_plot[xfield])]
    ylims = [min(df_plot[yfield]), max(df_plot[yfield])]

    range_x = xlims[1]-xlims[0]
    range_y = ylims[1]-ylims[0]

    x_add = (xlims[1]-xlims[0])*0.05
    y_add = (ylims[1]-ylims[0])*0.05

    xlims = [xlims[0]-x_add, xlims[1]+x_add]
    ylims = [ylims[0]-y_add, ylims[1]+y_add]

    # get those with some recall
    df_plot = df_plot[df_plot[yfield]>0].copy()
    sorted_hues = [h for h in sorted_hues if h in set(df_plot[hue_f])]
    sorted_phenotypes = [p for p in sorted_phenotypes if p in set(df_plot.index)]

    # define marker
    all_markers = "ovsX^<>"*10
    hue_to_marker = {}
    for I,hue_val in enumerate(sorted_hues): hue_to_marker[hue_val] = all_markers[I]

    # init fig
    ncols = len(sorted_expected_genes) + 1
    nrows = len(sorted_phenotypes)
    fig = plt.figure(figsize=(ncols*2, nrows*2)); Ip=1
    #fig = plt.figure(figsize=(ncols*5, nrows*5)); Ip=1

    # add each of data slice
    for Ir, phenotype in enumerate(sorted_phenotypes):
        df_plot_pheno = df_plot.loc[[phenotype]].copy()

        for Ic, gene in enumerate(["any"] + sorted_expected_genes):

            # get df for the pheno and  the genes
            if gene=="any": df_plot_final = df_plot_pheno.copy()
            else: df_plot_final = df_plot_pheno[df_plot_pheno[gene+"_found"]==True].copy()
            df_plot_final = df_plot_final.sort_values(by=["hue_idx"], ascending=[True]).set_index(hue_f)
            found_hues = set(df_plot_final.index)

            # get plot
            ax = plt.subplot(nrows, ncols, Ip); Ip+=1            
            for f in sorted_hues:
                if not f in found_hues: continue
                sns.scatterplot(data=df_plot_final.loc[[f]], x=xfield, y=yfield, linewidth=.6, s=25, edgecolor=palette_dict[f], facecolor="none", alpha=1, marker=hue_to_marker[f])
                #sns.scatterplot(data=df_plot_final.loc[[f]], x=xfield, y=yfield, linewidth=.5, s=23, edgecolor="red", facecolor="none", alpha=1, marker=hue_to_marker[f])

            # add an arrow to the best filters
            best_filters_df = df_plot_final[df_plot_final.unique_filterID==best_filterID].copy()
            if len(best_filters_df)==0: pass
            elif len(best_filters_df)>1: raise ValueError("should be one")
            else:

                best_filters_series = best_filters_df.iloc[0]
                x_arr = best_filters_series[xfield]
                y_arr = best_filters_series[yfield]
                color_arr = palette_dict[best_filters_series.name]

                x_start, y_start = (x_arr+range_x*0.2, y_arr)
                delta_x, delta_y = (-range_x*0.13, 0) 
                arrow_len = -delta_x

                plt.arrow(x_start, y_start, delta_x, delta_y, shape="full", width=range_y*0.03, length_includes_head=True, head_length=arrow_len/3, edgecolor=color_arr, lw=1.3, facecolor="white")
                #plt.arrow(x_arr, y_arr + range_y*0.1, 0, -range_y*0.05, width=range_x*0.023, shape="full", color=color_arr) # vertical

            # legend
            def get_legend_element(color, marker, label): return plt.Line2D([0], [0], marker=marker, markerfacecolor=color, markeredgecolor="none", markersize=5, label=label, linewidth=0, markeredgewidth=.3)

            if Ic==(ncols-1) and Ir==0:
                legend_handles = [get_legend_element("white", "o", hue_f)] + [get_legend_element(palette_dict[val], hue_to_marker[val], val) for val in reversed(sorted_hues)]
                ax.legend(handles=legend_handles, title="", fontsize=10, title_fontsize=10, bbox_to_anchor=(1.2, 1), loc="upper left", markerscale=1.5)

            # lims
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)

            # markers and lims
            ax.set_ylabel("")
            if Ic==0 and Ir=={4:1, 3:1}[nrows]:
                text_ylabel ={"n_sig_genes":"# genes w/ prot-altering hit (~recall)", "fraction_transitions_explained":"fraction pheno. transitions related to prot-altering hits (~recall)"}[yfield]
                ax.text(-range_x*0.4, {4:0, 3:0.5}[nrows], text_ylabel, rotation=90, ha="center", va="center")

            if Ic!=0: ax.set_yticklabels([])

            ax.set_xlabel("")
            if Ir==(nrows-1) and Ic=={6:3, 5:2, 4:2, 3:1}[ncols]: 
                text_xlabel = {"n_sig_variants_no_multi_clade":"# variants with non-multiclade hits (~fpr)", "log_n_sig_variants_no_multi_clade":"log # variants with non-multiclade hits (~fpr)"}[xfield]
                ax.text(0, -range_y*0.4, text_xlabel, ha="center")

            if Ir!=(nrows-1): ax.set_xticklabels([])

            # add titles
            if Ic==(ncols-1): ax.text(xlims[1] + 1, ylims[0] + (ylims[1]-ylims[0])/2, get_paper_dicotomic_pheno(phenotype), fontsize=10, va="center", ha="center", rotation=270)
            if Ir==0: 
                ax.set_title(gene.replace("Scer_HAL9", "TAC1"), fontsize=10)

            # add general title
            if Ir==0 and Ic==0:
                ax.text(0, ylims[1]+y_add*5, title_str, fontsize=10, ha="left")

            # add lines
            for y in ylines_y: plt.axhline(y, lw=.8, c="gray", ls="--")
            for x in xlines_x: plt.axvline(x, lw=.8, c="gray", ls="--")

    # get fig
    sns.despine()
    plt.subplots_adjust(wspace=0.05, hspace=0.15)# wspace (between cols) hspaces (btw rows) 
    save_fig(fig, "GWAS_parameter_effects/%s_%s.pdf"%(type_filt.replace("\n", "_"), hue_f.replace(" ", "_")), PlotsDir)
    plt.show()

def get_figures_each_test_parameter_changing(df_ML_performance, PlotsDir, gene_set_to_color, f_to_sorted_vals, min_consistency, max_consistency, threshold_consistency_genes, threshold_ROC_AUC):

    """get effects of single parameters"""

    df_ML_performance = df_ML_performance.copy()

    # definge graphics
    f_to_color_dict = {"aligner": {"bwa":"grey", "bowtie2":"tab:blue", "hisat2":"tab:green"},
                       "type split & samples": {'split=random, samples=all':"tab:red", 'split=random, samples=representative':"tab:pink", 'split=balanced, samples=all':"tab:blue", 'split=balanced, samples=representative':"tab:cyan"},
                       "consistency": {"none":"grey", "genes":"black", "variants":"tab:blue", "samples":"tab:red"},
                       "AUC diff": {".05":"grey", ".1":"olive", ".15":"black"},
                       "gene_set": gene_set_to_color}

    valid_parm_to_marker = {True:"o", False:"v"}

    # get figure
    sorted_phenos = ["FLC_resistance", "VRC_resistance", "source_blood_vs_other_sterile", "age_infant"]
    sorted_parms = ["aligner", "type split & samples", "consistency", "AUC diff", "gene_set"]

    ncols = len(sorted_phenos)
    nrows = len(sorted_parms)
    fig = plt.figure(figsize=(ncols*2, nrows*2)); Ip=1

    # one plot for each parameter
    for Ir, parameter_f in  enumerate(sorted_parms):
        for Ic, phenotype in enumerate(sorted_phenos):

            # get df
            df_plot = df_ML_performance[df_ML_performance.phenotype==phenotype].copy()
            top_min_AUC = min(df_plot[df_plot.is_top==True].AUC_test)
            top_min_consistency = min(df_plot[df_plot.is_top==True].mean_cross_splits_consistency_jaccard_genes)

            # sort by value
            df_plot = df_plot.sort_values(by=[parameter_f+"_idx"], ascending=False)

            # get plot
            ax = plt.subplot(nrows, ncols, Ip); Ip+=1
            color_palette = f_to_color_dict[parameter_f]

            for valid_t in [False, True]:
                df_p = df_plot[df_plot.valid_train==valid_t].copy()
                sns.scatterplot(data=df_p, x="mean_cross_splits_consistency_jaccard_genes", y="AUC_test", linewidth=.9, marker=valid_parm_to_marker[valid_t], s=20, edgecolor=[color_palette[x] for x in df_p[parameter_f]], facecolor="none")

            # legend
            def get_legend_element(color, marker, label): return plt.Line2D([0], [0], marker=marker, markerfacecolor=color, markeredgecolor="none", markersize=5, label=label, linewidth=0, markeredgewidth=.3)

            if Ic==3:

                legend_handles = [get_legend_element("white", "o", parameter_f)] + [get_legend_element(color_palette[val], "o", val) for val in f_to_sorted_vals[parameter_f]]

                if Ir==0: legend_handles += [get_legend_element("white", "o", "non-overfit parameters")] + [get_legend_element("black", valid_parm_to_marker[x], str(x)) for x in [True, False]]

                #+ [get_legend_element(gene_set_to_color[gene_set], "o", gene_set) for gene_set in f_to_sorted_vals["gene_set"]] + [get_legend_element("white", "o", "\nconsistency measure", markeredgecolor="white")] + [get_legend_element("black", consistency_to_marker[cons], cons) for cons in f_to_sorted_vals["consistency"]]
                
                ax.legend(handles=legend_handles, title="", fontsize=10, title_fontsize=10, bbox_to_anchor=(1.15,0.5), loc="center left", markerscale=1.5)

            # markers and lims
            ax.set_xlim([min_consistency-0.05, max_consistency+0.05])
            ax.set_ylim([0.45, 1.05])

            if Ic==0 and Ir==2: ax.set_ylabel("test $AUC_{ROC}$")
            else: ax.set_ylabel("")
            if Ic!=0: ax.set_yticklabels([])

            if Ir==4 and Ic==1: ax.set_xlabel("                           m1-m5 feature consistency (genes)")
            else: ax.set_xlabel("")
            if Ir!=4: ax.set_xticklabels([])

            # add titles
            if Ir==0: ax.set_title(get_paper_dicotomic_pheno(phenotype), fontsize=10)

            # add lines
            plt.axvline(threshold_consistency_genes, lw=.8, c="gray", ls="--")
            plt.axhline(threshold_ROC_AUC, lw=.8, c="gray", ls="--")

            plt.plot([top_min_consistency, top_min_consistency], [top_min_AUC, 1.05], lw=.8, c="black", ls="--") 
            plt.plot([top_min_consistency, max_consistency+0.05], [top_min_AUC, top_min_AUC], lw=.8, c="black", ls="--") 
            
            df_plot_top = df_plot[df_plot.is_top==True].copy()
            nparms = len(set(df_plot_top[parameter_f]))

            for Iv, (parameter_val, df_parm) in  enumerate(df_plot_top.groupby(parameter_f)):
                x_text = 0.7+Iv*0.1
                ax.text(x_text, 0.92, str(len(df_parm)), ha="left", fontsize=10, color=color_palette[parameter_val])
                if (Iv+1)<nparms: ax.text(x_text+0.06, 0.92, "/", ha="left", fontsize=10, color="black")

            # add an arrow for the best model
            df_plot_best = df_plot[df_plot.best_model==True].copy()
            if len(df_plot_best)!=1: raise ValueError("should be 1 len")
            x_arr = df_plot_best.iloc[0].mean_cross_splits_consistency_jaccard_genes
            y_arr = df_plot_best.iloc[0].AUC_test

            plt.arrow(x_arr+0.07, y_arr+0.07, -0.02, -0.02, width=0.01, shape="full", color=color_palette[df_plot_best.iloc[0][parameter_f]])

    # get fig
    sns.despine()
    plt.subplots_adjust(wspace=0.05, hspace=0.15)# wspace (between cols) hspaces (btw rows) 

    save_fig(fig, "test_parameters_performance_split_by_parm.pdf", PlotsDir)
    plt.show()


def plot_model_performance_scatters(df_ML_integrated_performance_rep, PlotsDir, target_phenotypes_ML, yfield="ROC_AUC_test_integrated_onlyCorrectModels", xfield="mean_cross_splits_consistency_jaccard_samples", row_f="type_gene_set",  hue_f="multiple_test_corr", palette={"none":"gray", "maxT_min":"red", "maxT_mean":"c"}, min_yval=0.51, ylim=[0.45, 1.05], xlim=[0.15, 1.05], hue_order=["maxT_min", "maxT_mean", "none"], xfield_offset=0, ylines_vals=[0.7], xlines_vals=[0.5]):

    """plots scatters"""

    # mean_cross_splits_consistency_jaccard_variants
    # mean_cross_splits_consistency_jaccard_genes
    # mean_cross_splits_consistency_jaccard_samples

    # redefine
    if target_phenotypes_ML is None:
        df_max_performance = df_ML_integrated_performance_rep.groupby("phenotype").apply(lambda df: max(df[yfield])).sort_values(ascending=False)
        target_phenotypes_ML = list(df_max_performance.index)

    # filter
    df_ML_integrated_performance_rep = df_ML_integrated_performance_rep[(df_ML_integrated_performance_rep[yfield]>=min_yval) & (df_ML_integrated_performance_rep.phenotype.isin(set(target_phenotypes_ML)))].copy()

    sorted_rows = sorted(set(df_ML_integrated_performance_rep[row_f]))
    ncols = len(target_phenotypes_ML)
    nrows = len(sorted_rows)
    fig = plt.figure(figsize=(ncols*2, nrows*2)); Ip=1

    for Ir, row_val in enumerate(sorted_rows):
        for Ic, phenotype in enumerate(target_phenotypes_ML):

            # get plot and add phenos
            df_plot = df_ML_integrated_performance_rep[(df_ML_integrated_performance_rep[row_f]==row_val) & (df_ML_integrated_performance_rep.phenotype==phenotype)].copy()[[hue_f, yfield, xfield, "phenotype"]].copy()

            # sort
            df_plot["hue_idx"] = df_plot[hue_f].map(dict(zip(hue_order, range(len(hue_order)))))
            df_plot = df_plot.sort_values(by=["hue_idx"], ascending=False)
            check_no_nans_in_df(df_plot)

            ax = plt.subplot(nrows, ncols, Ip); Ip+=1
            ax = sns.scatterplot(data=df_plot, x=xfield, y=yfield, hue=hue_f, palette=palette, s=20, linewidth=0.2, edgecolor="black")

            # edit plot
            if Ir==1 and Ic==0: ax.set_ylabel(yfield)
            else: ax.set_ylabel("")

            if Ir==(nrows-1) and Ic==3: ax.set_xlabel((" "*xfield_offset) + xfield)
            else: ax.set_xlabel("")

            ax.set_ylim(ylim)
            ax.set_xlim(xlim)

            if Ir==2 and Ic==0: ax.legend(bbox_to_anchor=(12,1.5), loc="upper left")
            elif len(df_plot)>0: ax.get_legend().remove()

            if Ir==0: ax.set_title(get_paper_dicotomic_pheno(phenotype), fontsize=10)

            if Ic==(ncols-1):
                ax.text(1.1, 0.75, row_val, fontsize=10, ha="left")

            if Ir!=(nrows-1):
                ax.set_xticklabels([])

            if Ic!=0:
                ax.set_yticklabels([])

            for y in ylines_vals: ax.axhline(y, ls="--", c="gray", lw=.8)
            for x in xlines_vals: ax.axvline(x, ls="--", c="gray", lw=.8)

            # add correlation
            if len(df_plot)>0:
                spearmanr, spearmanp =  scipy.stats.spearmanr(df_plot[xfield].values, df_plot[yfield].values, nan_policy="raise")
                if pd.isna(spearmanr):
                    spearmanr = 0.0
                    spearmanp = 1.0

                if spearmanp<0.05: text_suff="*"
                else: text_suff = ""
                ax.text(0.2, 0.95, "r=.%s%s"%(("%.2f"%spearmanr).split(".")[1], text_suff ), ha="left", fontsize=10)



    sns.despine()
    save_fig(fig, "violin_plots_performance/%s_vs_%s.pdf"%(xfield, yfield), PlotsDir)



def get_paper_renamed_fields(x):

    """renamed fields"""

    if x=="ROC_AUC_test_integrated_onlyCorrectModels":
        return "$AUC_{ROC}$ independent"

    else:
        raise ValueError("invalid %s"%x)

def plot_max_performance_different_test_parms(df_ML_integrated_performance_rep, fields_optimize_test, PlotsDir, threshold_ROC_AUC, target_phenotypes_ML, yfield, ylim=[0.45, 1.05], min_yval_set=0.5):

    """Plot max performance across different parameters"""

    # keep
    df_ML_integrated_performance_rep = df_ML_integrated_performance_rep[df_ML_integrated_performance_rep.phenotype.isin(set(target_phenotypes_ML))].copy()

    # change
    df_ML_integrated_performance_rep[yfield] = df_ML_integrated_performance_rep[yfield].apply(lambda y: max([y, min_yval_set]))

    # define colors
    pheno_to_color = get_value_to_color(target_phenotypes_ML, palette="tab10", type_color="hex", n=len(target_phenotypes_ML))[0]
    pheno_to_color = {get_paper_dicotomic_pheno(p):c for p,c in pheno_to_color.items()}

    # go through different fields
    for f in fields_optimize_test:

        # define the sorted fields
        if f=="aligner": sorted_vals = ["bwa_mem", "hisat2_no_spliced", "bowtie2_local"]
        else: sorted_vals = sorted(set(df_ML_integrated_performance_rep[f]))

        # init fun
        fig = plt.figure(figsize=(len(sorted_vals)*1.3 , 2))

        for I, type_max in enumerate(["pct95", "max"]):

            # get max function
            if type_max=="pct95": max_fun = (lambda x: np.percentile(x, 95))
            elif type_max=="pct75": max_fun = (lambda x: np.percentile(x, 75))
            elif type_max=="max": max_fun = max
            elif type_max=="mean": max_fun = np.mean
            else: raise ValueError("err")

            # get performance df
            max_f = "%s %s"%(type_max, yfield)
            df_max_performance = df_ML_integrated_performance_rep.groupby(["phenotype", f]).apply(lambda df: pd.Series({max_f: max_fun(df[yfield])})).reset_index()
            df_max_performance["val_idx"] = df_max_performance[f].map(dict(zip(sorted_vals, range(len(sorted_vals)))))

            df_max_performance["pheno_idx"] = df_max_performance["phenotype"].map(dict(zip(target_phenotypes_ML, range(len(target_phenotypes_ML)))))
            df_max_performance = df_max_performance.sort_values(by=["pheno_idx", "val_idx"])
            df_max_performance["phenotype"] = df_max_performance.phenotype.apply(get_paper_dicotomic_pheno)
            check_no_nans_in_df(df_max_performance)

            # get plot
            ax = plt.subplot(1, 2, I+1)
            ax = sns.lineplot(data=df_max_performance, x="val_idx", y=max_f, hue="phenotype", palette=pheno_to_color, marker="o", markeredgecolor="gray", linewidth=2)

            # lines
            for y in [0.5, 0.7, 0.9]:
                plt.axhline(y, lw=.7, c="gray", ls="--")

            # edit
            if I==1:
                ax.legend(bbox_to_anchor=(1.1,1), loc="upper left")
            else:
                ax.get_legend().remove()

            ax.set_ylim(ylim)
            ax.set_xlim([-0.2, len(sorted_vals)-0.8])
            ax.set_xticks(list(range(len(sorted_vals))))
            ax.set_xticklabels(sorted_vals, rotation=90)
            ax.set_xlabel(f)
            ax.set_title(type_max, fontsize=10)

            ax.set_yticks([0.5, 0.7, 0.9])
            if I==0: ax.set_ylabel(get_paper_renamed_fields(yfield))
            else: 
                ax.set_yticklabels([])
                ax.set_ylabel("")

        sns.despine()
        plt.subplots_adjust(wspace=0.05, hspace=0.05)# wspace (between cols) hspaces (btw rows)
        save_fig(fig, "test_parm_effect_max/%s.pdf"%f, PlotsDir)



def plot_model_performance_violins(df_ML_integrated_performance_rep, PlotsDir, yfield="ROC_AUC_test_integrated_onlyCorrectModels", row_f="type_gene_set", hue_f="multiple_test_corr", palette={"none":"gray", "maxT_min":"red", "maxT_mean":"c"}, missing_y_val=0.5, missing_hue_val="none", min_yval=0.51, ylim=[0.45, 1.05], ylines_vals=[0.7, 1.0], min_yval_set=0.5):

    """plots the model performance"""

    # define the max perf
    print("max %s"%yfield)
    df_max_performance = df_ML_integrated_performance_rep.groupby("phenotype").apply(lambda df: max(df[yfield])).sort_values(ascending=False)
    print(df_max_performance)

    # sort by max
    sorted_phenotypes = list(df_max_performance.index)
    #sorted_phenotypes = ['FLC_resistance', 'VRC_resistance', 'MIF_resistance', 'agar_invasiveness', 'age_infant', 'age_old', 'age_young', 'source_ASPIRADO_BRONQUIAL', 'source_EXUDADO_CUTANEO', 'source_ORINA', 'source_blood', 'source_blood_vs_other_sterile', 'source_sterile']


    # filter
    df_ML_integrated_performance_rep = df_ML_integrated_performance_rep[df_ML_integrated_performance_rep[yfield]>=min_yval].copy()

    # reset
    df_ML_integrated_performance_rep[yfield] = df_ML_integrated_performance_rep[yfield].apply(lambda y: max([y, min_yval_set]))

    # print phenotypes with some maxT performanc
    print("phenotypes some maxT performance")
    print(sorted(set(df_ML_integrated_performance_rep[df_ML_integrated_performance_rep.multiple_test_corr!="none"].phenotype)))


    # fig

    sorted_rows = sorted(set(df_ML_integrated_performance_rep[row_f]))
    nrows = len(sorted_rows)
    fig = plt.figure(figsize=(10, nrows*3)); Ip=1

    for Ir, row_val in enumerate(sorted_rows):

        # get plot and add phenos
        df_plot = df_ML_integrated_performance_rep[df_ML_integrated_performance_rep[row_f]==row_val].copy()[[hue_f, yfield, "phenotype"]].copy().sort_values(by=[hue_f, "phenotype"])

        missing_phenos = sorted(set(sorted_phenotypes).difference(set(df_plot.phenotype)))
        df_plot_missing = pd.DataFrame({"phenotype":missing_phenos, yfield:missing_y_val, hue_f:missing_hue_val})
        df_plot = pd.concat([df_plot, df_plot_missing])

        if set(df_plot.phenotype)!=set(sorted_phenotypes):
            raise ValueError("invalid")

        # sort
        df_plot["pheno_idx"] = df_plot.phenotype.map(dict(zip(sorted_phenotypes, range(len(sorted_phenotypes)))))
        df_plot = df_plot.sort_values(by=["pheno_idx", hue_f])
        check_no_nans_in_df(df_plot)

        ax = plt.subplot(nrows, 1, Ip); Ip+=1
        #ax = sns.violinplot(data=df_plot, x="phenotype", y=yfield, hue=hue_f, palette=palette, fill=False)
        ax = sns.boxplot(data=df_plot, x="phenotype", y=yfield, hue=hue_f, palette=palette, fill=False, flierprops={"markersize": 3})

        # add max points
        df_plot_max = df_plot.groupby(["phenotype", hue_f, "pheno_idx"]).apply(lambda df: pd.Series({"maxy": max(df[yfield])})).reset_index()
        df_plot_max = df_plot_max.sort_values(by=["pheno_idx", hue_f])
        #ax = sns.stripplot(data=df_plot_max, x="phenotype", y="maxy", hue=hue_f, palette=palette, dodge=True)

        # edit plot
        if Ir==0: ax.legend(bbox_to_anchor=(1,1), loc="upper left")
        else: ax.get_legend().remove()

        if Ir=={3:1, 1:0, 2:1}[nrows]: ax.set_ylabel(yfield)
        else: ax.set_ylabel("")

        if Ir==(len(sorted_rows)-1):
            ax.set_xticklabels([get_paper_dicotomic_pheno(x) for x in sorted_phenotypes], rotation=90)

        else:
            ax.set_xticklabels([])

        ax.set_xlabel("")

        for y in ylines_vals: ax.axhline(y, ls="--", c="gray", lw=.8)
        for x in range(0, len(sorted_phenotypes)): ax.axvline(x+0.5, ls="--", c="gray", lw=.8)

        ax.set_ylim(ylim)
        ax.set_title(row_f + "=" + row_val, fontsize=10)

    sns.despine()
    save_fig(fig, "violin_plots_performance/%s_distribution_hue=%s.pdf"%(yfield, hue_f), PlotsDir)



def run_ML_predictors_parameter_info(GenoPheno_dir_all, DataDir, df_ML_results_all_file, df_ML_parameter_info_file, list_cmds):

    """Generates parameter info file"""

    if file_is_empty(df_ML_parameter_info_file):
        print("getting parameter info")

        # define dirs
        Integrated_ML_classifier_dirs = "%s/ML_classifiers_integrated_results"%DataDir
        folders_remove = []

        # go through each phenotype
        df_ML_parameter_info_all = pd.DataFrame()
        for phenotype in os.listdir(Integrated_ML_classifier_dirs):
            print(phenotype)

            # define dirs
            classifier_integration_outdir = "%s/%s"%(Integrated_ML_classifier_dirs, phenotype)
            paths_training_file = "%s/paths_training.tab"%(classifier_integration_outdir)
            outdir_ML_integration = "%s/classifier_results"%(classifier_integration_outdir)
            ML_results_file = "%s/best_models_test_performance.tab"%outdir_ML_integration

            # define tmpdir
            outdir_pheno = "%s.gen.%s"%(df_ML_parameter_info_file, phenotype)

            # get the parameter info
            final_file = "%s/final_file.txt"%outdir_pheno
            if file_is_empty(final_file):

                args_files = [outdir_pheno, final_file, ML_results_file, paths_training_file]
                cmd = "python %s/Cparapsilosis_popgen/get_ML_parm_info.py %s"%(ParentDir, " ".join(args_files))
                list_cmds.append(("parm_info", cmd, "ancestral_GWAS_env", 112))
                continue

            df_ML_parameter_info = get_tab_as_df_or_empty_df("%s/parm_info.tab"%outdir_pheno)
            df_ML_parameter_info["phenotype"] = phenotype.split("-")[1]
            df_ML_parameter_info_all = pd.concat([df_ML_parameter_info_all, df_ML_parameter_info])
            folders_remove.append(outdir_pheno)

        # debug jobs to run
        if len(list_cmds)>0: return

        # save at the end
        print("saving")
        for f in folders_remove:
            delete_folder(f)

        save_object(df_ML_parameter_info_all, df_ML_parameter_info_file)




def get_interesting_name_FeatsDF(row):
    
    """Takes a features dataframe row that has gene_name and Scerevisiae_orthologs and returns the interesting name"""

    if not pd.isna(row["gene_name"]): return row["gene_name"]
    elif not pd.isna(row["Scerevisiae_orthologs"]): return "Scer_%s"%(row["Scerevisiae_orthologs"])
    else: return row["geneID"]

def get_gene_features_df(DataDir):

    """Gets the gene features df"""

    # load df of CGD
    CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
    nheader_lines = len(list(range(len([line for line in open(CGD_chromosomal_features, "r") if line.startswith("!")]))))
    geneFeatures_df =  pd.read_csv(CGD_chromosomal_features, sep="\t", names=["geneID", "gene_name", "aliases", "feature_type", "chromosome", "start", "end", "strand", "primary_CGDID", "secondary_CGDID", "description", "date", "sequence_coordinate_version_date", "blank1", "blank2", "date_geneName_reservation", "is_it_standardName", "Scerevisiae_orthologs"], skiprows=nheader_lines) 

    # add the name
    geneFeatures_df["final_name"] = geneFeatures_df.apply(get_interesting_name_FeatsDF, axis=1)

    # check geneID
    gff_df = load_gff3_intoDF("%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir)  
    all_protein_coding_genes = set(gff_df[gff_df.feature.isin({"CDS", "mRNA"})].upmost_parent)
    gff_df = gff_df[gff_df.feature=="gene"].copy()
    missing_genes = set(gff_df.upmost_parent).difference(geneFeatures_df.geneID)

    if len(missing_genes)>0:
        raise ValueError("some missing genes: %s"%geneFeatures_df)

    # add the prot coding
    geneFeatures_df["protein_coding_gene"] = geneFeatures_df.geneID.isin(all_protein_coding_genes)

    # print geneIDs that are not in the gff
    geneIDs_not_gff = set(geneFeatures_df.geneID).difference(set(gff_df.upmost_parent))
    #print("Gene IDs not in gff:\n", geneFeatures_df[(geneFeatures_df.geneID.isin(geneIDs_not_gff)) & (geneFeatures_df.feature_type!="pseudogene")])

    return geneFeatures_df

def get_gwas_df_filt_for_benchmarking(df_gwas, all_protein_coding_genes, GenoPheno_dir_all, dict_info):

    """Gets a gwas df that has protein altering changes for each group, ready for benchmarking"""

    # omly keep sync
    df_gwas = df_gwas[(df_gwas.gwas_method=="synchronous")].copy()

    # check general things
    if not all(df_gwas.nodes_withGeno_clearNodes>=2):
        raise ValueError("nodes_withGeno_clearNodes should be >=2")

    if not all(df_gwas.nodes_withoutGeno_clearNodes>=1):
        raise ValueError("nodes_withoutGeno_clearNodes should be >=2")

    # keep only some gropings
    df_gwas = df_gwas[(df_gwas.type_collapsing.isin({"none", "variant_effect", "domains", "genes"})) 
                    & (df_gwas.type_genes.isin({"all_genes", "only_protein_coding"}))].copy()

    df_gwas = df_gwas[(df_gwas.type_collapsing.isin({"none", "variant_effect"}))
                    | ((df_gwas.type_collapsing.isin({"genes", "domains"})) & (df_gwas.type_mutations.isin({"non_syn_non_truncating_muts", "non_syn_muts", "truncating_muts"})))].copy().reset_index(drop=True)

    # check that type_genes is only one
    for tc in set(df_gwas.type_collapsing):
        if len(set(df_gwas[df_gwas.type_collapsing==tc].type_genes))!=1:
            raise ValueError("the type_genes should be always==1 for %s"%tc)

    # if GWAS field is empty, return it
    if len(df_gwas)==0:
        df_gwas["affected_geneID_set"] = []
        df_gwas["type_convergence"] = []
        return df_gwas

    # add the affected gene for each hit (keep only prot altering changes)
    annot_df = get_tab_as_df_or_empty_df("%s/pheno-%s_aligner-%s/convergence_gwas/gwas_inputs_%s/variant_annotations_df_all.tab"%(GenoPheno_dir_all, dict_info["phenotype"].replace(" ", "_"), dict_info["aligner"], dict_info["sample_set"]))
    annot_df = annot_df[(annot_df.is_synonymous==False) & (annot_df.is_protein_coding_gene==True)].copy()
    all_prot_altering_effects = set(annot_df.Gene + "|" + annot_df.short_variant_effect)
    varID_to_affected_protAlteringGenes = dict(annot_df[["variantID_across_samples", "Gene"]].drop_duplicates().groupby("variantID_across_samples").apply(lambda df: set(df.Gene)))

    def get_affected_geneID_set_one_hit(r):
        if r.type_collapsing=="genes":
            return_val = ({r.group}).intersection(all_protein_coding_genes)

        elif r.type_collapsing=="domains":
            return_val = {r.group.split("#")[1]}

        elif r.type_collapsing=="none":
            if r.group in varID_to_affected_protAlteringGenes:
                return_val = varID_to_affected_protAlteringGenes[r.group]

            else:
                return_val = set()

        elif r.type_collapsing=="variant_effect":
            if r.group in all_prot_altering_effects:
                return_val = {r.group.split("|")[0]}

            else:
                return_val = set()

        else:
            raise ValueError("invalid grouping %s"%r)

        # check and return
        if type(return_val)!=set:
            raise ValueError("invalud return: %s"%return_val)

        return return_val

    df_gwas["affected_geneID_set"] = df_gwas[["group", "type_collapsing"]].apply(get_affected_geneID_set_one_hit, axis=1)

    # keep only hits with some gene or with some variants
    df_gwas = df_gwas[(df_gwas.affected_geneID_set.apply(len)>0) | (df_gwas.type_collapsing=="none")].copy()

    # add the type of hit according to convergence
    df_gwas["set_clades"] = df_gwas.cladeID_names_GenoAndPheno.apply(lambda x: set(x.split(",")).difference({"-1", "none"}))

    def get_type_convergence_by_clades(set_clades):
        if len(set_clades)==0: "no_or_unknown_convergence"
        elif len(set_clades)==1: "single_clade_convergence"
        elif len(set_clades)>1: "multi_clade_convergence"
        else: raise ValueError("invalid")

    df_gwas["type_convergence"] = df_gwas.set_clades.apply(get_type_convergence_by_clades)

    # checks
    strange_genes = set.union(*df_gwas.affected_geneID_set).difference(all_protein_coding_genes)
    if len(strange_genes)>0:
        raise ValueError("strange_genes: %s"%strange_genes)

    return df_gwas

def get_formatted_gene_name_with_nGenoAndPheno(df_gwas_filt, target_gene):

    """Gets, for one target gene name, the text that has the # of geno and pheno"""

    # get the gene df
    df_g = df_gwas_filt[df_gwas_filt.affected_gene_set.apply(lambda x: target_gene in x)].copy()

    if len(df_g)==0:
        raise ValueError("gene cannot be empty")

    # get the formatted gene
    return "%s=%i"%(target_gene, max(df_g.nodes_GenoAndPheno))

def get_filterd_gwas_df_one_filter(r, df_gwas):

    """gets filtered GWAS"""

    # filter based on basic things
    df_gwas_filt = df_gwas[(df_gwas.epsilon>=r.min_epsilon) & (df_gwas.nodes_GenoAndPheno>=r.min_nodes_GenoAndPheno)].copy().reset_index(drop=True)

    # filter based on pvals
    idx_correct_pvals = pd.Series([True]*len(df_gwas_filt), index=df_gwas_filt.index)
    for pval_f in r.tuple_pval_fields:
        idx_correct_pvals = (idx_correct_pvals) & (df_gwas_filt[pval_f]<r.pvalue_threshold)

    # filter
    df_gwas_filt = df_gwas_filt[idx_correct_pvals].copy()
    return df_gwas_filt

def get_gwas_filters_df_with_GWAS_behchmarking_stats(r, df_gwas, expected_genes, n_nodes_pheno_transitions_all):

    """Takes a row with all the filter info and the gwas df and returns information about the filtering"""

    # filter GWAS based on the filters
    df_gwas_filt = get_filterd_gwas_df_one_filter(r, df_gwas)

    # check
    if len(df_gwas_filt)>0:
        if not all(df_gwas_filt.nodes_GenoAndPheno>=2): raise ValueError("all should have convergence")

    # add stats
    if len(df_gwas_filt)==0:        
        r["set_sig_genes"] = set()
        r["set_expected_genes_found"] = set()
        r["n_sig_variants"] = 0

        r["set_sig_genes_multi_clade"] = set()
        r["n_sig_variants_single_clade"] = 0
        r["n_sig_variants_no_multi_clade"] = 0
        r["fraction_transitions_explained"] = 0.0

    else:
        r["set_sig_genes"] = set.union(*df_gwas_filt.affected_gene_set)
        r["set_expected_genes_found"] = set.union(*df_gwas_filt.affected_gene_set).intersection(expected_genes)
        r["n_sig_variants"] = sum(df_gwas_filt.type_collapsing=="none")

        df_multi_clade = df_gwas_filt[df_gwas_filt.type_convergence=="multi_clade_convergence"].copy()
        if len(df_multi_clade)>0: r["set_sig_genes_multi_clade"] = set.union(*df_multi_clade.affected_gene_set)
        else:  r["set_sig_genes_multi_clade"] = set()

        r["n_sig_variants_single_clade"] = sum((df_gwas_filt.type_convergence=="single_clade_convergence") & (df_gwas_filt.type_collapsing=="none"))
        r["n_sig_variants_no_multi_clade"] = sum((df_gwas_filt.type_convergence!="multi_clade_convergence") & (df_gwas_filt.type_collapsing=="none"))

        df_gwas_filt_some_gene = df_gwas_filt[df_gwas_filt.affected_gene_set.apply(len)>0].copy()
        if n_nodes_pheno_transitions_all>0 and len(df_gwas_filt_some_gene)>0:
            r["fraction_transitions_explained"] = len(set.union(*df_gwas_filt_some_gene.affected_nodes_set)) / n_nodes_pheno_transitions_all
        else:
            r["fraction_transitions_explained"] = 0.0

    # general stats
    r["set_expected_genes_found_wGenoAndPheno"] = {get_formatted_gene_name_with_nGenoAndPheno(df_gwas_filt, gene) for gene in r.set_expected_genes_found}
    r["n_sig_genes"] = len(r["set_sig_genes"])
    r["n_sig_genes_multi_clade"] = len(r["set_sig_genes_multi_clade"])
    r["n_expected_genes_found"] = len(r["set_expected_genes_found"])
    r["fraction_expected_genes_found"] = r["n_expected_genes_found"] / len(expected_genes)

    # check
    if r["fraction_transitions_explained"]>1: raise ValueError("cannot be >1")

    # add the expected genes found
    for g in sorted(expected_genes):

        # add whether it is found
        r["%s_found"%g] = (g in r["set_sig_genes"])

        # add the affected nodes by this gene
        if r["%s_found"%g] is True:
            df_g = df_gwas_filt[(df_gwas_filt.affected_gene_set.apply(lambda x: g in x))].copy()
            if len(df_g)==0: raise ValueError("can't be 0. %s, %s"%(df_gwas_filt, df_gwas_filt.keys()))
            r["%s_affected_nodes"%g] = set.union(*df_g.affected_nodes_set)

        else: 
            r["%s_affected_nodes"%g] = set()

    # add 
    r["any_pval"] = any(r[['pval_chi_square_phenotypes', 'pval_GenoAndPheno_phenotypes', 'pval_chi_square_phenotypes_maxT', 'pval_epsilon_phenotypes_maxT']])
    r["any_pval_correction"] = r['pval_chi_square_phenotypes_maxT'] or r['pval_epsilon_phenotypes_maxT'] or r['correction_method']!="none"

    return r

def get_GWAS_benchmarking_df_one_result(log_text, dict_info, geneID_to_name, expected_genes, GenoPheno_dir_all, all_protein_coding_genes, all_expected_genes, gwas_filters_df_file):

    """Runs the GWAS benchmarking for one result"""

    gwas_ID = log_text + " " + "_".join([str(dict_info[f]) for f in ["phenotype", "aligner", "ASR_method", "min_support", "sample_set"]])

    # load the GWAS results and filter them to have only protein altering variants & df_gwas
    print("loading gwas", gwas_ID)
    df_gwas = get_tab_as_df_or_empty_df(dict_info["raw_GWAS_stats"])
    df_gwas = get_gwas_df_filt_for_benchmarking(df_gwas, all_protein_coding_genes, GenoPheno_dir_all, dict_info)

    # checks
    if len(df_gwas)>0 and len(set(df_gwas.nodes_withPheno))!=1:
        raise ValueError("nodes_withPheno should always be the same")

    # add the affected gene names
    df_gwas["affected_gene_set"] = df_gwas.affected_geneID_set.apply(lambda g_ids: {geneID_to_name[g] for g in g_ids})

    # add node names set
    df_gwas["affected_nodes_set"] = df_gwas.node_names_GenoAndPheno.apply(lambda x: set(x.split(",")).difference({"none"}))
    #df_gwas["affected_nodes_set"] = df_gwas.node_names_GenoAndPheno.apply(lambda x: set(x.split(",")))

    # define all the nodes that have transitions
    df_nodes_info = get_tab_as_df_or_empty_df(get_dir(dict_info["raw_GWAS_stats"]) + "/phenotypes_info.tab")
    set_nodes_pheno_transitions = set(df_nodes_info[df_nodes_info.phenotype_transition==1].node_name)
    n_nodes_pheno_transitions_all = len(set_nodes_pheno_transitions)

    if len(df_gwas)>0:
        strange_nodes_gwas = set.union(*df_gwas["affected_nodes_set"]).difference(set_nodes_pheno_transitions)
        if len(strange_nodes_gwas)>0:
            raise ValueError("there are strange nodes in affected_nodes_set")

    # keep some fields
    df_gwas = df_gwas[["affected_gene_set", "affected_geneID_set", "epsilon", "nodes_GenoAndPheno", "nodes_withPheno", "type_collapsing", "type_convergence", "affected_nodes_set"] + [f for f in df_gwas.keys() if f.startswith("pval_")]].copy()

    # load the filters df
    gwas_filters_df = load_object(gwas_filters_df_file)

    # get the df with benchmarking info
    print("getting benchmarking", gwas_ID)
    GWAS_benchmarking_df = gwas_filters_df.apply(get_gwas_filters_df_with_GWAS_behchmarking_stats, df_gwas=df_gwas, expected_genes=expected_genes, n_nodes_pheno_transitions_all=n_nodes_pheno_transitions_all, axis=1)

    # add extra fields
    GWAS_benchmarking_df["expected_genes_str"] = ",".join(sorted(expected_genes))
    GWAS_benchmarking_df["expected_genes_found_str"] = GWAS_benchmarking_df.set_expected_genes_found.apply(lambda x: ",".join(sorted(x)))
    GWAS_benchmarking_df["expected_genes_found_wGenoAndPheno_str"] = GWAS_benchmarking_df.set_expected_genes_found_wGenoAndPheno.apply(lambda x: ",".join(sorted(x)))

    if len(df_gwas)>0:
        GWAS_benchmarking_df["nodes_withPheno"] = df_gwas.nodes_withPheno.iloc[0]

    else:
        GWAS_benchmarking_df["nodes_withPheno"] = -1

    for k,v in dict_info.items():
        if k=="raw_GWAS_stats": continue
        GWAS_benchmarking_df[k] = v

    print("returning", gwas_ID)
    GWAS_benchmarking_df = GWAS_benchmarking_df.drop(["tuple_pval_fields", "set_sig_genes", "treefile"], axis=1).copy()
    return GWAS_benchmarking_df

def get_list_pval_fields_one_combination_pvals(filters_dict, all_possible_fields):

    """Gets the list of pval fields to consider, depending on the types of pvals to consider"""

    # init
    pval_fields = []

    # add the maxT pvals
    for f in ["pval_chi_square_phenotypes_maxT", "pval_epsilon_phenotypes_maxT"]:

        # debug
        if filters_dict[f] is False: continue

        # add
        if filters_dict["correction_all_groups"] is True:
            pval_fields.append("%s_allGroups"%f)

        else:
            pval_fields.append(f)

    # add the other pval fiekds
    for f in ["pval_chi_square_phenotypes", "pval_GenoAndPheno_phenotypes", "pval_fisher"]:

        # debug
        if filters_dict[f] is False: continue

        # get the uncorrected pval
        if filters_dict["correction_method"]=="none":
            pval_fields.append(f)

        # corrected p
        else:
            
            corr_f = "%s_%s"%(f, filters_dict["correction_method"])
            if filters_dict["correction_all_groups"] is True:
                pval_fields.append("%s_allGroups"%corr_f)

            else:
                pval_fields.append(corr_f)

    # strange fields
    if len(set(pval_fields).difference(all_possible_fields))>0:
        raise ValueError("not valid fields")

    if len(set(pval_fields))!=len(pval_fields):
        raise ValueError("non unique fields")

    return pval_fields
 
def get_gwas_filters_df(ProcessedDataDir, GenoPheno_dir_all):

    """Generates the GWAS filters df"""

    # go through different combinations of filters and define a dataframe with them
    filters_df_file = "%s/GWAS_filters_df_all_filters.py"%ProcessedDataDir
    if file_is_empty(filters_df_file):
        print("getting GWAS filters df...")

        # define the possible filters
        all_possible_fields = set(get_tab_as_df_or_empty_df("%s/pheno-POS_resistance_aligner-bwa_mem/convergence_gwas/gwas_results_representative_samples_min_sup_70_ASR_DOWNPASS/integrated_GWAS_stats.tab"%(GenoPheno_dir_all)).keys())

        # init dict
        filters_dict = {}; I = 0

        # define the pval fileds
        for pval_chi_square_phenotypes in [True, False]:
            for pval_GenoAndPheno_phenotypes in [True, False]:
                for pval_fisher in [False]: # I skip fisher as it is not used by any of the NatMicro good samples
                    for pval_chi_square_phenotypes_maxT in [True, False]:
                        for pval_epsilon_phenotypes_maxT in [True, False]:

                            # require at least one p val field -> old
                            #if sum([pval_chi_square_phenotypes, pval_GenoAndPheno_phenotypes, pval_fisher, pval_chi_square_phenotypes_maxT, pval_epsilon_phenotypes_maxT])==0: continue

                            # discard combinations where there is a mix of pvalue types
                            if any([pval_chi_square_phenotypes, pval_GenoAndPheno_phenotypes]) and any([pval_chi_square_phenotypes_maxT, pval_epsilon_phenotypes_maxT]):
                                continue

                            # define the correction methods (only applied to the first three pvalues)
                            if any([pval_chi_square_phenotypes, pval_GenoAndPheno_phenotypes, pval_fisher]):
                                array_correction_methods = ["bonferroni", "fdr_bh", "none"]

                            elif any([pval_chi_square_phenotypes_maxT, pval_epsilon_phenotypes_maxT]):
                                array_correction_methods = ["none"]

                            elif sum([pval_chi_square_phenotypes, pval_GenoAndPheno_phenotypes, pval_fisher, pval_chi_square_phenotypes_maxT, pval_epsilon_phenotypes_maxT])==0:
                                array_correction_methods = ["none"]

                            else:
                                raise ValueError("some correction methods should be set")

                            # iterate through each corr method, and also other parameters
                            for correction_method in array_correction_methods:

                                # define whether to correct on all groups
                                if correction_method!="none" or any([pval_chi_square_phenotypes_maxT, pval_epsilon_phenotypes_maxT]):
                                    correction_all_groups_array = [True, False]

                                else:
                                    correction_all_groups_array = [False]

                                for correction_all_groups in correction_all_groups_array:

                                    # define the array of corrections for unfiltered p values
                                    any_uncorrected = any([pval_chi_square_phenotypes, pval_GenoAndPheno_phenotypes]) and correction_method=="none"
                                    if any_uncorrected and not any([pval_chi_square_phenotypes_maxT, pval_epsilon_phenotypes_maxT]):
                                        array_pvalue_thresholds = [0.0001, 0.001, 0.05]

                                    else:
                                        array_pvalue_thresholds = [0.05]

                                    # define the numeric filters
                                    for pvalue_threshold in array_pvalue_thresholds:
                                        for min_epsilon in [0.1, 0.2, 0.3]:
                                            for min_nodes_GenoAndPheno in [2, 3]:

                                                # keep filters
                                                filters_dict[I] = {"pval_chi_square_phenotypes":pval_chi_square_phenotypes,
                                                                   "pval_GenoAndPheno_phenotypes":pval_GenoAndPheno_phenotypes,
                                                                   "pval_fisher":pval_fisher,
                                                                   "pval_chi_square_phenotypes_maxT":pval_chi_square_phenotypes_maxT,
                                                                   "pval_epsilon_phenotypes_maxT":pval_epsilon_phenotypes_maxT,

                                                                   "correction_method":correction_method, 
                                                                   "correction_all_groups":correction_all_groups,
                                                                   "min_epsilon":min_epsilon, 
                                                                   "min_nodes_GenoAndPheno":min_nodes_GenoAndPheno,
                                                                   "pvalue_threshold":pvalue_threshold,
                                                                   "filter_I":I}
                                                # add pvals
                                                filters_dict[I]["tuple_pval_fields"] = tuple(get_list_pval_fields_one_combination_pvals(filters_dict[I], all_possible_fields))

                                                # update
                                                I+=1


        filters_df = pd.DataFrame(filters_dict).transpose()
        
        # checks
        print("checking...")
        if len(filters_df)!=len(filters_df.drop_duplicates()): raise ValueError("df is not unique")

        lens_set_pvalue_combs = set()
        for idx, df in filters_df.groupby(["min_epsilon", "min_nodes_GenoAndPheno", "pvalue_threshold"]):

            set_pvalue_combs = set(df.tuple_pval_fields)
            if len(set_pvalue_combs)!=len(df): raise ValueError("the pval combs should be unique")
            if idx[2]==0.05:
                if sum(df.tuple_pval_fields.apply(len)==0)!=1: raise ValueError("should be 1 with no pvals")
                lens_set_pvalue_combs.add(len(set_pvalue_combs))

        if len(lens_set_pvalue_combs)!=1: raise ValueError("should be only one")
        print("there are %i filters"%len(filters_df))

        # save
        save_object(filters_df, filters_df_file)

    return filters_df_file

def get_jaccard_index_sets(x1, x2):

    """get jaccard dist between sets"""

    if len(x1.union(x2))==0: return 0.0
    return len(x1.intersection(x2)) / len(x1.union(x2))

def print_stats_GWAS_benchmarking_df_all(GWAS_benchmarking_df_file, PlotsDir, ProcessedDataDir):

    """print stats of the df"""

    df_benchmarking_overall_stats_file = "%s/GWAS_benchmarking_overall_stats.xlsx"%PlotsDir
    if file_is_empty(df_benchmarking_overall_stats_file):
        print("exploring GWAS results...")

        # load
        GWAS_benchmarking_df = load_object(GWAS_benchmarking_df_file).reset_index(drop=True)
        print("in total there are %i filters"%(len(set(GWAS_benchmarking_df.unique_filterID))))

        # get the total # of filters
        n_total_filters = len(set(GWAS_benchmarking_df.unique_filterID))

        # add fields
        GWAS_benchmarking_df["gene_to_ntransitions_dict"] = GWAS_benchmarking_df.set_expected_genes_found_wGenoAndPheno.apply(lambda x: {g.split("=")[0] : {int(val) for val in g.split("=")[1].split("-")} for g in x})

        print(GWAS_benchmarking_df.groupby("phenotype").apply(len))

        comb_to_n = GWAS_benchmarking_df.groupby(["phenotype", "aligner", "min_support", "ASR_method", "sample_set"]).apply(len)
        print("N parameters per comb", set(comb_to_n))
        print(comb_to_n)

        # keep only FLC and VRC
        GWAS_benchmarking_df = GWAS_benchmarking_df[GWAS_benchmarking_df.phenotype.isin({"VRC_resistance", "FLC_resistance"})].copy()

        # print filters
        df_filters = load_object(ProcessedDataDir + "/GWAS_filters_df_all_filters.py")
        print(df_filters)

        # function to get stats
        #bool_to_str = {True:"yes", False:"no"}
        def get_stats_series_benchmarking_one_comb(df):

            # add name
            phenotype, expected_genes_found_str, min_support = df.name
            data_dict = {"phenotype":phenotype, "min_support":min_support, "genes_found":expected_genes_found_str.replace("Scer_HAL9", "TAC1")}

            # add the stats about the genes
            if expected_genes_found_str!="":

                def get_min_max_n_transitions_one_gene(g):
                    list_n_transitions = set.union(*df.gene_to_ntransitions_dict.apply(lambda x: x[g]))
                    return "%s=%i-%i"%(g, min(list_n_transitions), max(list_n_transitions)) 

                data_dict["n_Geno&Pheno"] = ",".join([get_min_max_n_transitions_one_gene(g) for g in expected_genes_found_str.split(",")])

            else:
                data_dict["n_Geno&Pheno"] = ""

            # defie all pvals
            pval_fields_all = ['pval_chi_square_phenotypes', 'pval_GenoAndPheno_phenotypes', 'pval_fisher', 'pval_chi_square_phenotypes_maxT', 'pval_epsilon_phenotypes_maxT']

            # define other stats
            data_dict["n_filters"] = len(df)
            data_dict["n_any_pval"] = sum(df[pval_fields_all].apply(any, axis=1))
            data_dict["n_maxT"] = sum(df[["pval_chi_square_phenotypes_maxT", "pval_epsilon_phenotypes_maxT"]].apply(any, axis=1))
            data_dict["n_maxT_all_groups"] = sum((df[["pval_chi_square_phenotypes_maxT", "pval_epsilon_phenotypes_maxT"]].apply(any, axis=1)) & (df["correction_all_groups"]))

            for correction_method in ["bonferroni", "fdr_bh"]:
                data_dict["n_%s"%correction_method] = sum(df.correction_method==correction_method)
                data_dict["n_%s_all_groups"%correction_method] = sum((df.correction_method==correction_method) & (df["correction_all_groups"]))

            return pd.Series(data_dict)

        # get stats
        df_stats = GWAS_benchmarking_df.sort_values(by=["phenotype", "expected_genes_found_str"]).groupby(["phenotype", "expected_genes_found_str", "min_support"]).apply(get_stats_series_benchmarking_one_comb).reset_index(drop=True)

        # add the number of transitions per pheno
        pheno_to_transitions = GWAS_benchmarking_df[(GWAS_benchmarking_df.ASR_method=="DOWNPASS") & (GWAS_benchmarking_df.aligner=="bwa_mem") & (GWAS_benchmarking_df.min_support==70) & (GWAS_benchmarking_df.filter_I==max(GWAS_benchmarking_df.filter_I)) & (GWAS_benchmarking_df.sample_set=="representative_samples")].set_index("phenotype").nodes_withPheno
        print(pheno_to_transitions)
        df_stats["n_transitions_pheno"] = df_stats.phenotype.map(pheno_to_transitions)

        # sort
        df_stats = df_stats.sort_values(by=["phenotype", "min_support", "n_filters", "genes_found"], ascending=[True, False, False, True]).reset_index(drop=True)#.drop(["type_phenotype"], axis=1)

        # modify
        for f in ["n_filters", "n_any_pval", "n_maxT", "n_maxT_all_groups", "n_bonferroni", "n_bonferroni_all_groups", "n_fdr_bh", "n_fdr_bh_all_groups"]:
            df_stats[f] = df_stats[f].apply(lambda x: "%i (%.2f%s)"%(x, (x/n_total_filters)*100, "%"))

        # save
        df_stats.to_excel(df_benchmarking_overall_stats_file)

    df_stats = pd.read_excel(df_benchmarking_overall_stats_file)
    print("\n", df_stats.drop(["Unnamed: 0"], axis=1))
    print()

def get_affected_geneID_gwas_hit_NatMicro(r):

    """gets the affected geneID"""

    # check that there is some p val
    #pval_fields_NatMicro = ['pval_chi_square_maxT', 'pval_epsilon_maxT', 'pval_chi_square_phenotypes', 'pval_GenoAndPheno_phenotypes', 'pval_fisher', 'pval_chi_square_phenotypes_bonferroni', 'pval_GenoAndPheno_phenotypes_bonferroni', 'pval_fisher_bonferroni']
    #pval_fields_NatMicro = ['pval_chi_square_phenotypes', 'pval_GenoAndPheno_phenotypes', 'pval_fisher']

    #if not any(r[pval_fields_NatMicro]<0.05):
    #    raise ValueError("some p should be <0.05")

    # get gene
    if r.type_collapsing=="none":
        return r.description.split(" (")[1].split(",")[0]

    elif r.type_collapsing=="genes":
        return r.group_name

    elif r.type_collapsing=="domains":
        return r.group_name.split("#")[1]

    else:
        raise ValueError("invalid %s"%dict(r))

def load_CGD_features_file(CGD_chromosomal_features):

    """loads"""

    nheader_lines = len(list(range(len([line for line in open(CGD_chromosomal_features, "r") if line.startswith("!")]))))
    df_CGD_features =  pd.read_csv(CGD_chromosomal_features, sep="\t", names=["geneID", "gene_name", "aliases", "feature_type", "chromosome", "start", "end", "strand", "primary_CGDID", "secondary_CGDID", "description", "date", "sequence_coordinate_version_date", "blank1", "blank2", "date_geneName_reservation", "is_it_standardName", "Scerevisiae_orthologs"], skiprows=nheader_lines)  

    return df_CGD_features

def get_filtered_NatMicro_GWAS(df_GWAS_NatMicro, df_features_NatMicro):

    """gets the filtered GWAS df"""

    # keep some fields of df_features_NatMicro
    df_features_NatMicro = df_features_NatMicro[['species', 'gff_upmost_parent', 'Scerevisiae_orthologs', 'orthofinder_orthocluster']].copy()


    # filter GWAS
    df_GWAS_NatMicro = df_GWAS_NatMicro[df_GWAS_NatMicro.type_collapsing.isin({'none', 'domains', 'genes'}) & (df_GWAS_NatMicro.drug.isin({'POS', 'VRC', 'FLC', 'MIF', 'ANI', 'ITR', 'AMB'}))][["type_collapsing", "species", "drug", "group_name", "description"]].drop_duplicates().copy()

    # add gene
    df_GWAS_NatMicro["gff_upmost_parent"] = df_GWAS_NatMicro.apply(get_affected_geneID_gwas_hit_NatMicro, axis=1)

    # check filter
    df_GWAS_NatMicro = df_GWAS_NatMicro[["species", "drug", "gff_upmost_parent"]].drop_duplicates()
    check_no_nans_in_df(df_GWAS_NatMicro)
    df_GWAS_NatMicro = df_GWAS_NatMicro.merge(df_features_NatMicro.copy(), on=["species", "gff_upmost_parent"], how="left", validate="many_to_one")
    check_no_nans_in_df(df_GWAS_NatMicro[["species", "gff_upmost_parent", "drug"]])

    return df_GWAS_NatMicro

def get_NatMicro_datasets_with_Cpara_geneID(df_GWAS_NatMicro_lowc, df_GWAS_NatMicro_highc, df_selection_NatMicro, df_features_NatMicro):

    """get mappings """

    print("add mappings")

    # add orthologs as set
    for df in [df_GWAS_NatMicro_lowc, df_GWAS_NatMicro_highc, df_selection_NatMicro, df_features_NatMicro]:
        def get_nan_or_set(x):
            if pd.isna(x): return x
            else: return set(x.split(","))
        df["Scerevisiae_orthologs"] = df["Scerevisiae_orthologs"].apply(get_nan_or_set)

    # map OGs
    print("mappings")
    df_features_NatMicro_Cpara = df_features_NatMicro[df_features_NatMicro.species=="Candida_parapsilosis"].copy()
    if not all(df_features_NatMicro_Cpara.orthofinder_orthocluster.apply(lambda x: x.startswith("OG"))):
        raise ValueError("all OGs should be such")
    
    if any(df_features_NatMicro.orthofinder_orthocluster.apply(lambda x: "," in str(x))):
        raise ValueError("should not be , in orthofinder_orthocluster")

    NatMicro_OC_to_Cpara_geneIDs = dict(df_features_NatMicro_Cpara[~pd.isna(df_features_NatMicro_Cpara.orthofinder_orthocluster)].groupby("orthofinder_orthocluster").apply(lambda df_g: set(df_g.gff_upmost_parent)))
    all_genes_OC = set(NatMicro_OC_to_Cpara_geneIDs.keys())

    # map Scer values
    df_features_NatMicro_Cpara_Scer = df_features_NatMicro_Cpara[~pd.isna(df_features_NatMicro_Cpara.Scerevisiae_orthologs)][["Scerevisiae_orthologs", "gff_upmost_parent"]].copy()
    all_Scer_orthos = sorted(set.union(*df_features_NatMicro_Cpara_Scer.Scerevisiae_orthologs))
    NatMicro_ScerOrth_to_Cpara_geneIDs = dict(zip(all_Scer_orthos, 
                                         map(lambda x: set(df_features_NatMicro_Cpara_Scer[df_features_NatMicro_Cpara_Scer.Scerevisiae_orthologs.apply(lambda ort: x in ort)].gff_upmost_parent), all_Scer_orthos)))

    all_lens = {len(x) for x in NatMicro_ScerOrth_to_Cpara_geneIDs.values()}
    if all_lens!={1}: raise ValueError("invalud")
    NatMicro_ScerOrth_to_Cpara_geneID = {k:next(iter(v)) for k,v in NatMicro_ScerOrth_to_Cpara_geneIDs.items()}
    all_NatMicro_ScerOrths_Cpara = set(NatMicro_ScerOrth_to_Cpara_geneID.keys())

    # get geneIDs
    def get_Cpara_geneIDs_based_clusters(r):

        # define genes
        if not pd.isna(r.Scerevisiae_orthologs):
            Scer_genes = {NatMicro_ScerOrth_to_Cpara_geneID[x] for x in r.Scerevisiae_orthologs if x in all_NatMicro_ScerOrths_Cpara}

        else:
            Scer_genes = set()

        if r.orthofinder_orthocluster in all_genes_OC:
            OG_genes = NatMicro_OC_to_Cpara_geneIDs[r.orthofinder_orthocluster]

        else:
            OG_genes = set()

        # get final genes
        if r.species=="Candida_parapsilosis": 
            final_genes = {r.gff_upmost_parent}

        elif len(Scer_genes)>0:
            final_genes = Scer_genes

        elif len(OG_genes)>0:
            final_genes = OG_genes

        else:
            final_genes = set()

        return {x.split("gene-")[1] for x in final_genes}

    for df in [df_GWAS_NatMicro_lowc, df_GWAS_NatMicro_highc, df_selection_NatMicro]:
        df["Cparapsilosis_geneIDs"] = df.apply(get_Cpara_geneIDs_based_clusters, axis=1)
        

    # filter
    df_GWAS_NatMicro_lowc = df_GWAS_NatMicro_lowc[df_GWAS_NatMicro_lowc.Cparapsilosis_geneIDs.apply(len)>0].copy()
    df_GWAS_NatMicro_highc = df_GWAS_NatMicro_highc[df_GWAS_NatMicro_highc.Cparapsilosis_geneIDs.apply(len)>0].copy()
    df_selection_NatMicro = df_selection_NatMicro[df_selection_NatMicro.Cparapsilosis_geneIDs.apply(len)>0].copy()

    return df_GWAS_NatMicro_lowc, df_GWAS_NatMicro_highc, df_selection_NatMicro



def get_phenotype_to_all_possible_genes_for_ML(DataDir, ManualDataDir, df_expected_gene_info_file):

    """Gets all possible genes for each phenotype"""

    print("Getting all possible Genes")
    if file_is_empty(df_expected_gene_info_file):

        # filterings
        gene_features_df = get_gene_features_df(DataDir)
        gene_features_df = gene_features_df[(gene_features_df.chromosome!=mtDNA_chrom) & (gene_features_df.feature_type.isin({'ORF|Verified', 'ORF|Uncharacterized'})) & ~(gene_features_df.geneID.apply(lambda x: x.endswith(".1")))].copy()

        # define function parsing genes
        df_description = gene_features_df[["geneID", "description", "final_name"]].copy()
        all_genes = set(df_description.final_name)
        def get_geneIDs_by_description(terms_include, gene_names_exclude, print_data=False):

            matching_idx_include = df_description.description.apply(lambda d: 
                                   any(map(lambda term_i: term_i in d.lower(), terms_include)))

            if len(gene_names_exclude.difference(all_genes))>0: raise ValueError("invalid genes")
            matching_idx_exclude = df_description.final_name.isin(gene_names_exclude)

            df_match = df_description[(matching_idx_include) & ~(matching_idx_exclude)].copy()
            if print_data is True: print("\n\n---\n", "Terms: %s, excluded genes: %s\n"%(terms_include, gene_names_exclude), "\n".join(df_match.final_name + ": "+ df_match.description), "\n---\n\n")
            return set(df_match.geneID)

        # map each gene ID to the name
        gene_features_df["NatMicro_geneID"] = "gene-" + gene_features_df.geneID
        gene_name_to_ID = dict(gene_features_df.set_index("final_name").geneID)
        geneID_to_name = dict(gene_features_df.set_index("geneID").final_name)

        # load chromosomal features
        CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
        df_CGD_features = load_CGD_features_file(CGD_chromosomal_features)
        
        # load the chrom features of nat micro
        df_CGD_features_NatMicro = load_CGD_features_file(ManualDataDir+"/NatMicro2024_C_parapsilosis_s01-m03-r45_features.tab")
        df_CGD_features_NatMicro["NatMicro_geneID"] = "gene-" + df_CGD_features_NatMicro.geneID
        df_CGD_features_NatMicro = df_CGD_features_NatMicro[(df_CGD_features_NatMicro.chromosome!="mito_C_parapsilosis_CDC317") & (df_CGD_features_NatMicro.feature_type.isin({'ORF|Verified', 'ORF|Uncharacterized'}))]

        # load datasets about selection and resistance from NatMicro Schikora-Tamarit 2024 
        print("loading dfs")
        df_selection_NatMicro = get_tab_as_df_or_empty_df(ManualDataDir+"/TableS2-Genes_under_selection.csv")
        #df_selection_all_genes_NatMicro = get_tab_as_df_or_empty_df(ManualDataDir+"/TableS2-Selection_scores_all_genes.csv")
        df_features_NatMicro = get_tab_as_df_or_empty_df(ManualDataDir+"/TableS2-Gene_features.csv")
        df_GWAS_NatMicro_lowc = get_tab_as_df_or_empty_df(ManualDataDir+"/TableS3-Low-confidence_GWAS_hits.csv")
        df_GWAS_NatMicro_highc = get_tab_as_df_or_empty_df(ManualDataDir+"/TableS3-High-confidence_GWAS_hits.csv")

        # filter features to keep strange genes out
        df_features_NatMicro = df_features_NatMicro[(df_features_NatMicro.species!="Candida_parapsilosis") | (df_features_NatMicro.gff_upmost_parent.isin(set(df_CGD_features_NatMicro.NatMicro_geneID)))].copy()

        # get filtered GWAS dfs
        print("getting filtered GWAS")
        df_GWAS_NatMicro_lowc = get_filtered_NatMicro_GWAS(df_GWAS_NatMicro_lowc, df_features_NatMicro.copy())
        df_GWAS_NatMicro_highc = get_filtered_NatMicro_GWAS(df_GWAS_NatMicro_highc, df_features_NatMicro.copy())
        df_selection_NatMicro = df_selection_NatMicro[["species", "gff_upmost_parent", "orthofinder_orthocluster", "Scerevisiae_orthologs", "type_var"]].copy().drop_duplicates()

        # check that all our genes are in nat micro genes
        all_genes_Cpara_NatMicro = set(df_features_NatMicro[df_features_NatMicro.species=="Candida_parapsilosis"].gff_upmost_parent)
        strange_genes_in_gff = set(gene_features_df["NatMicro_geneID"]).difference(all_genes_Cpara_NatMicro)
        if len(strange_genes_in_gff)>0:
            raise ValueError("all genes from gene features should be in NatMicro. extra: %s"%strange_genes_in_gff)

        # check that all genes of GWAS & selection are in gene_features_df
        strange_genes_in_features = all_genes_Cpara_NatMicro.difference(set(gene_features_df["NatMicro_geneID"]))
        if len(strange_genes_in_features)>0:
            raise ValueError("all genes from NatMicro should be in features. extra: %s"%strange_genes_in_features)
      
        # add mappings
        df_GWAS_NatMicro_lowc, df_GWAS_NatMicro_highc, df_selection_NatMicro = get_NatMicro_datasets_with_Cpara_geneID(df_GWAS_NatMicro_lowc, df_GWAS_NatMicro_highc, df_selection_NatMicro, df_features_NatMicro)

        all_genes_after_map = set.union(*[set.union(*df.Cparapsilosis_geneIDs) for df in [df_GWAS_NatMicro_lowc, df_GWAS_NatMicro_highc, df_selection_NatMicro]])
        strange_genes_after_map = all_genes_after_map.difference(set(gene_features_df["geneID"]))
        if len(strange_genes_after_map)>0:
            raise ValueError("strange_genes_after_map: %s"%strange_genes_after_map)

        # define the gene IDs under selection
        geneIDs_selection_all = set.union(*df_selection_NatMicro.Cparapsilosis_geneIDs) # type_var is DEL, if_INDEL, DUP, SNP
        geneIDs_selection_SNP = set.union(*df_selection_NatMicro[df_selection_NatMicro.type_var=="SNP"].Cparapsilosis_geneIDs)


        # define phenotypes
        azoles = ["FLC_resistance", "VRC_resistance", "POS_resistance", "IVZ_resistance", "ITR_resistance"]
        echinocandins = ["ANI_resistance", "CAS_resistance", "MIF_resistance"]
        virulence = ['agar_invasiveness', 'age_infant', 'age_old', 'age_young', 'biofilm_formation', 'microfluidics_Adherent', 'microfluidics_Big_Clumps', 'microfluidics_Little_clump', 'pseudohyphae_formation', 'source_ABSCESO', 'source_ABSCESO_ABDOMINAL', 'source_ASPIRADO_BRONQUIAL', 'source_ASPIRADO_TRAQUEAL', 'source_BIOPSIA_PULMONAR', 'source_EOF', 'source_ESPUTO', 'source_EXUDADO_CUTANEO', 'source_EXUDADO_DE_HERIDA', 'source_EXUDADO_OROFARINGEO', 'source_EXUDADO_OTICO', 'source_EXUDADO_RECTAL', 'source_EXUDADO_RESPIRATORIO_ALTO', 'source_HC', 'source_HERIDA', 'source_LAVADO_BRONCOALVEOLAR', 'source_LCR', 'source_LIQUIDO_ASCITICO', 'source_LIQUIDO_BILIAR', 'source_MUESTRA_AMBIENTAL', 'source_ORINA', 'source_PIEL', 'source_RESPIRATORIAS', 'source_SANGRE_COMPLETA', 'source_TEJIDO_VASCULAR', 'source_VALVULA_CARDIACA', 'source_blood', 'source_blood_vs_environmental', 'source_blood_vs_other_sterile', 'source_respiratory', 'source_skin', 'source_sterile', 'source_sterile_vs_environmental']
        disinfectants_phenos = ['ETHANOL_tolerance', 'HYPLOCHLORITE_tolerance', 'H202_tolerance', 'CHLORHEXIDINE_tolerance', 'ZWITTERGENT 3-14_tolerance', 'SURFANIO_tolerance']

        # define the target genes
        df_expected_gene_info = pd.DataFrame()
        for phenotype in (azoles + echinocandins + virulence + disinfectants_phenos + ["AMB_resistance", "5FC_resistance"]):

            # define the different gene names
            if phenotype in azoles:
                known_gene_names = ["ERG11", "ERG3", "Scer_HAL9", "MRR1", "NDT80", "UPC2"]
                geneIDs_by_description = get_geneIDs_by_description(["azole", "ergosterol"], {"Scer_THI6", "Scer_ADE17", "Scer_ADE13", "HIS3", "Scer_THI4", "Scer_ADE1", "Scer_HIS6", "HIS7", "Scer_ADE2"})
                NatMicro_drugs = {'POS', 'VRC', 'FLC', 'ITR'}

            elif phenotype in echinocandins:
                known_gene_names = ["FKS1", "FKS2", "ERG3"]
                geneIDs_by_description = get_geneIDs_by_description(["echinocandin", "fungin", "beta-d-glucan", "beta-glucan", "beta glucan", "drug"], {"Scer_ERG1"})
                NatMicro_drugs = {'MIF', 'ANI'}

            elif phenotype in virulence:
                known_gene_names = []
                geneIDs_by_description = get_geneIDs_by_description(["virulence", "virulen",  "biofilm", "invasive", "morphology", "invasion", "adhesion", "adhesin", "adherence", "pseudohyphae", "hypha", "filamentous", "immun", "host"], set())
                NatMicro_drugs = set()

            elif phenotype=="AMB_resistance":
                known_gene_names = ["ERG2", "ERG3", "Scer_ERG5", "ERG6", "ERG25"]
                geneIDs_by_description = get_geneIDs_by_description(["ergosterol", 'drug', 'polyene', 'amphotericin'], set()) 
                NatMicro_drugs = {"AMB"}

            elif phenotype=="5FC_resistance":
                known_gene_names = ["Scer_FCY1", "Scer_FCY2", "Scer_FUR1", "Scer_CHS7"]
                geneIDs_by_description = get_geneIDs_by_description(["flucytosine", "fluorocytosine", "5-FC", "5FC", "FCY1", "FCY2", "FUR1", "CHS7", "cytosine deaminase", "cytosine permease", "uracil phosphoribo", "uracil", "drug"], set())
                NatMicro_drugs = set()

            elif phenotype in disinfectants_phenos:
                known_gene_names = []
                geneIDs_by_description = get_geneIDs_by_description(["disinfectant", "chlorhexidine", "hypochlorite", "H2O2", "stress", "drug"], set())
                NatMicro_drugs = set()


            else:
                raise ValueError("invalid pheno")

            # define the GWAS genes
            if len(NatMicro_drugs)>0:
                geneIDs_gwas_highc = set.union(*df_GWAS_NatMicro_highc[df_GWAS_NatMicro_highc.drug.isin(NatMicro_drugs)].Cparapsilosis_geneIDs)
                geneIDs_gwas_lowc = set.union(*df_GWAS_NatMicro_lowc[df_GWAS_NatMicro_lowc.drug.isin(NatMicro_drugs)].Cparapsilosis_geneIDs)

            else:
                geneIDs_gwas_highc = set()
                geneIDs_gwas_lowc = set()

            """
            # check inclusion
            if len(geneIDs_gwas_highc.difference(geneIDs_gwas_lowc))>0:
                print(geneIDs_gwas_highc.difference(geneIDs_gwas_lowc))
                raise ValueError("all high conf should be in lowc")
            """

            # get kniwn geneIDs
            geneIDs_known = {gene_name_to_ID[gn] for gn in known_gene_names}

            # add to df
            for gene_group, geneID_set in [("known", geneIDs_known), ("by_description", geneIDs_by_description), ("GWAS_high", geneIDs_gwas_highc), ("GWAS_low", geneIDs_gwas_lowc), ("selection_all", geneIDs_selection_all), ("selection_SNP", geneIDs_selection_SNP)]:
            

                # create df
                df_info = pd.DataFrame({"geneID":sorted(geneID_set), "gene_group":gene_group, "phenotype":phenotype})
                df_info["gene_name"] = df_info.geneID.map(geneID_to_name)
                check_no_nans_in_df(df_info)

                # add
                df_expected_gene_info = pd.concat([df_expected_gene_info, df_info]).reset_index(drop=True)

        # save
        save_object(df_expected_gene_info, df_expected_gene_info_file)

    return load_object(df_expected_gene_info_file)

def get_phenotype_to_expected_genes(DataDir):

    """Maps different phenotypes to the expected genes"""

    # define phenos
    azoles = ["FLC_resistance", "VRC_resistance", "POS_resistance", "IVZ_resistance", "ITR_resistance"]
    echinocandins = ["ANI_resistance", "CAS_resistance", "MIF_resistance"]
    target_phenotypes = set(azoles).union(set(echinocandins))

    # map each phenotype to the expected genes (sorted by importance)
    phenotype_to_expected_genes = {p:[] for p in target_phenotypes}
    for p in azoles:
        for g in ["ERG11", "MRR1", "Scer_HAL9", "NDT80"]: # HAL9 is the ortholog of TAC1 in C. albicans. NDT80 and HAL9 found in EPOM experiments
            phenotype_to_expected_genes[p].append(g)

    for p in echinocandins:
        for g in ["FKS1", "FKS2"]:
            phenotype_to_expected_genes[p].append(g)

    # map each gene ID to the name
    gene_features_df = get_gene_features_df(DataDir)
    gene_name_to_ID = dict(gene_features_df.set_index("final_name").geneID)

    return {pheno : [(name, gene_name_to_ID[name]) for name in gene_names] for pheno, gene_names in phenotype_to_expected_genes.items()}

def get_literature_expected_genes_drug_types():

    """gets the expected genes"""

    type_to_genes =  {"azole": ["ERG11", "ERG3", "Scer_HAL9", "MRR1", "NDT80", "UPC2"],
                      "echinocandin": ["FKS1", "FKS2", "ERG3"],
                      "AMB": ["ERG2", "ERG3", "Scer_ERG5", "ERG6", "ERG25"],
                      "5FC": ["Scer_FCY1", "Scer_FCY2", "Scer_FUR1", "Scer_CHS7"]}

    return type_to_genes

def chunks(l, n):
    
    """Yield successive n-sized chunks from a list l"""
    
    for i in range(0, len(l), n):
        yield l[i:i + n]

def run_function_in_parallel_generate_df(min_mem_per_thread, threads, mem_per_thread, inputs_fn, function_run, check_nans=False, df_creation_fun=pd.DataFrame, chunk_mult=1):

    """Runs in parallel some jobs by chunks"""

    # define the # of jobs to run in parallel
    if mem_per_thread<min_mem_per_thread:
        njobs_parallel = int(threads / (min_mem_per_thread / mem_per_thread))

    else:
        njobs_parallel = int(threads)

    # prepare inputs
    chunk_size = njobs_parallel*chunk_mult
    df_all = pd.DataFrame()

    nchunks_total = len(list(chunks(list(range(len(inputs_fn))), chunk_size)))

    for Ichunk, inputs_fn_chunk in enumerate(chunks(inputs_fn, chunk_size)):
        print("Working on chunk %i/%i. chunk_size=%i"%(Ichunk+1, nchunks_total, chunk_size))

        with multiproc.Pool(njobs_parallel) as pool:
            df_chunk = df_creation_fun(pool.starmap(function_run, inputs_fn_chunk, chunksize=1))
            pool.close()
            pool.terminate()

        df_all = pd.concat([df_all, df_chunk]).reset_index(drop=True)

        if check_nans is True:
            check_no_nans_in_df(df_all)

    return df_all


def bechmark_GWAS_parameters(df_raw_gwas_files_all, GenoPheno_dir_all, ProcessedDataDir, DataDir):

    """Runs benchmarking of different GWAS parameters to see what is best in terms of ASR, aligner and min_support"""

    print("benchmarking GWAS parameters...")
    GWAS_benchmarking_df_all_file = "%s/GWAS_benchmarking_df_all_some_phenos.py"%ProcessedDataDir

    if file_is_empty(GWAS_benchmarking_df_all_file):

        # keep some settings
        #target_phenotypes = {p for p in sorted(set(df_raw_gwas_files_all.phenotype)) if (not p.startswith("source_")) and (not p.startswith("age_"))} # benchmark only for lab measured phenotypes
        target_phenotypes = set(df_raw_gwas_files_all.phenotype)
        #target_phenotypes = {"ZWITTERGENT 3-14_tolerance"}
        #azole_phenotypes = {"FLC_resistance", "VRC_resistance", "POS_resistance", "ITR_resistance", "IVZ_resistance"} # for these we have a set of expected genes
        df_raw_gwas_files = get_filtered_df_by_dict(df_raw_gwas_files_all, {"phenotype": target_phenotypes, "sample_set":{"representative_samples", "all_samples", "balanced_samples"}})

        # define all expected genes
        type_to_genes = get_literature_expected_genes_drug_types()

        # map each phenotype to the expected genes
        phenotype_to_expected_genes = {p: set(make_flat_listOflists(type_to_genes.values())) for p in target_phenotypes}

        # map each gene ID to the name
        gene_features_df = get_gene_features_df(DataDir)
        geneID_to_name = dict(gene_features_df.set_index("geneID").final_name)

        # checks
        if not all([x in geneID_to_name.values() for x in set.union(*phenotype_to_expected_genes.values())]):
            raise ValueError("all gene names should be in the expected genes")

        strange_expected_genes = set.union(*phenotype_to_expected_genes.values()).difference(set(geneID_to_name.values()))
        if len(strange_expected_genes)>0:
            raise ValueError("strange expected genes: %s"%strange_expected_genes)


        # define the protein coding gebes
        all_protein_coding_genes = set(gene_features_df[gene_features_df["protein_coding_gene"]].geneID)

        # defin the filters df
        gwas_filters_df_file = get_gwas_filters_df(ProcessedDataDir, GenoPheno_dir_all)
        #print(load_object(gwas_filters_df_file))

        # benchmark in parallel for each file
        print("Benchmarking each filter set in parallel for %i cases..."%len(df_raw_gwas_files))
        inputs_fn = list(df_raw_gwas_files.apply(lambda r: [dict(r), geneID_to_name, phenotype_to_expected_genes[r.phenotype], GenoPheno_dir_all, all_protein_coding_genes, set(gene_features_df.geneID), gwas_filters_df_file], axis=1).values)

        inputs_fn = [tuple(["job %i/%i"%(I+1, len(df_raw_gwas_files))] + x) for I,x in enumerate(inputs_fn)]

        min_mem_per_thread = 2
        mem_per_thread = 2
        GWAS_benchmarking_df_all = run_function_in_parallel_generate_df(min_mem_per_thread, 112, mem_per_thread, inputs_fn, get_GWAS_benchmarking_df_one_result, check_nans=True, df_creation_fun=pd.concat, chunk_mult=2)

        # add the unique filterID
        print("Finished, adding cases...")
        GWAS_benchmarking_df_all["unique_filterID"] = GWAS_benchmarking_df_all.aligner
        for f in ["ASR_method", "min_support", "sample_set", "filter_I"]:
            GWAS_benchmarking_df_all["unique_filterID"] += ("_" + GWAS_benchmarking_df_all[f].apply(str))

        # add the # sig genes across all filters
        for f in ["n_sig_genes", "n_sig_genes_multi_clade", "n_sig_variants", "n_sig_variants_single_clade", "n_sig_variants_no_multi_clade", "fraction_transitions_explained"]:
            for fun_obj, fun_name in [(max, "max"), (min, "min")]:

                filter_to_val = GWAS_benchmarking_df_all.groupby("unique_filterID").apply(lambda df: fun_obj(df[f]))
                GWAS_benchmarking_df_all["%s_%s"%(fun_name, f)] = GWAS_benchmarking_df_all.unique_filterID.map(filter_to_val)

        check_no_nans_in_df(GWAS_benchmarking_df_all)

        # checks
        if len(GWAS_benchmarking_df_all)!=len(GWAS_benchmarking_df_all[["unique_filterID", "phenotype"]].drop_duplicates()):
            raise ValueError("filtID should be unique")

        # save
        print("saving...")
        save_object(GWAS_benchmarking_df_all, GWAS_benchmarking_df_all_file)

    # return 
    return GWAS_benchmarking_df_all_file

def is_good_parmID_GWAS_filtering(df):

    """defines if a parmID is correct"""

    return all(df.selected_expected_genes_tuple == df.expected_genes_found_tuple)

def get_filtered_GWAS_dfs_one_phenotype(phenotype, raw_GWAS_stats, gwas_filters_dict):

    """gets filtered GWAS for one pheno"""
    print("running gwas for %s"%phenotype)

    # load
    df_gwas = get_tab_as_df_or_empty_df(raw_GWAS_stats)

    # keep sync
    df_gwas = df_gwas[(df_gwas.gwas_method=="synchronous")].copy()

    # check general things
    if not all(df_gwas.nodes_withGeno_clearNodes>=2):
        raise ValueError("nodes_withGeno_clearNodes should be >=2")

    if not all(df_gwas.nodes_withoutGeno_clearNodes>=1):
        raise ValueError("nodes_withoutGeno_clearNodes should be >=2")

    # add
    df_gwas["phenotype"] = phenotype

    # get the filtered GWAS
    df_gwas = get_filterd_gwas_df_one_filter(pd.Series(gwas_filters_dict), df_gwas)

    # keep
    if len(df_gwas)==0: print("WARNING: No hits for %s"%phenotype)

    return df_gwas

def print_stats_different_expected_gene_sets(GWAS_benchmarking_df_file, PlotsDir):

    """print stats about different datasets"""

    df_print_file = PlotsDir + "/parameters_per_diff_genes.xlsx"
    if file_is_empty(df_print_file):

        # load df
        print("getting the stats of different gene sets")
        GWAS_benchmarking_df = load_object(GWAS_benchmarking_df_file)
        if len(GWAS_benchmarking_df)!=len(GWAS_benchmarking_df[["phenotype", "unique_filterID"]].drop_duplicates()): raise ValueError("should be unique")    

        # keep only two phenos
        target_phenotypes = {"VRC_resistance", "FLC_resistance"}
        GWAS_benchmarking_df = GWAS_benchmarking_df[GWAS_benchmarking_df.phenotype.isin(target_phenotypes)].copy()

        # add correction
        GWAS_benchmarking_df["any_correction"] = (GWAS_benchmarking_df[['pval_chi_square_phenotypes_maxT', 'pval_epsilon_phenotypes_maxT']].apply(any, axis=1)) | (GWAS_benchmarking_df.correction_method!="none")

        # filter parms
        idx_filt = (GWAS_benchmarking_df.any_correction)
        idx_filt = (idx_filt) & (GWAS_benchmarking_df.min_epsilon>0)
        GWAS_benchmarking_df = GWAS_benchmarking_df[idx_filt].copy()

        # prints
        print("getting parms...")
        flc_expected_genes_tuples = make_flat_listOflists([list(itertools.combinations(["MRR1", "ERG11", "NDT80"], n)) for n in [1, 2, 3]]) + [tuple([])]
        vrc_expected_genes_tuples = [("NDT80",), ("NDT80", "ERG11"), tuple([])]

        print_dict_all = {} 

        for flc_expected_genes in flc_expected_genes_tuples:
            for vrc_expected_genes in vrc_expected_genes_tuples:
                print(flc_expected_genes, vrc_expected_genes)

                if not "MRR1" in flc_expected_genes: continue

                # define if the genes are as expected
                pheno_to_expected_genes = {"FLC_resistance": flc_expected_genes, "VRC_resistance": vrc_expected_genes} # in the paper we justify this choice
                GWAS_benchmarking_df["selected_expected_genes_tuple"] = GWAS_benchmarking_df.phenotype.map(pheno_to_expected_genes).apply(sorted).apply(tuple)
                GWAS_benchmarking_df["expected_genes_found_tuple"] = GWAS_benchmarking_df.set_expected_genes_found.apply(sorted).apply(tuple)
                GWAS_benchmarking_df["genes_as_expected"] = GWAS_benchmarking_df.selected_expected_genes_tuple == GWAS_benchmarking_df.expected_genes_found_tuple

                # define the valid parms working for both
                parmID_to_is_valid_parm = GWAS_benchmarking_df[["genes_as_expected", "unique_filterID"]].groupby("unique_filterID").apply(lambda df: sum(df.genes_as_expected))==len(target_phenotypes)
                GWAS_benchmarking_df["valid_parmID"] = GWAS_benchmarking_df.unique_filterID.map(parmID_to_is_valid_parm)
                GWAS_benchmarking_df_filt = GWAS_benchmarking_df[GWAS_benchmarking_df.valid_parmID].copy()

                if len(GWAS_benchmarking_df_filt)>0:

                    print_fields = ['pval_chi_square_phenotypes_maxT', 'pval_epsilon_phenotypes_maxT', 'pval_chi_square_phenotypes', 'pval_GenoAndPheno_phenotypes', 'correction_method', "correction_all_groups", 'min_epsilon', 'min_support', 'aligner', 'sample_set'] # pval_fisher, min_nodes_GenoAndPheno, ASR_method
                    
                    print_dict = {f:sorted(set(GWAS_benchmarking_df_filt[f])) for f in print_fields}
                    print_dict["n_filters"] =  len(set(GWAS_benchmarking_df_filt.unique_filterID))
                    print_dict["FLC_genes"] = flc_expected_genes
                    print_dict["VRC_genes"] = vrc_expected_genes

                    print_dict_all[(flc_expected_genes, vrc_expected_genes)] = print_dict

        # get as df
        df_print = pd.DataFrame(print_dict_all).transpose().reset_index(drop=True)
        df_print.to_excel(df_print_file)

    df_print = pd.read_excel(df_print_file)
    print("\n", df_print.drop(["Unnamed: 0"], axis=1), "\n")

def get_filtered_GWAS_dfs_old(GWAS_benchmarking_df_file, df_raw_gwas_files_all, GenoPheno_dir_all, ProcessedDataDir, DataDir, PlotsDir):

    """gets the single used GWAS dfs"""

    df_GWAS_filt_all_file = ProcessedDataDir + "/df_GWAS_filt_all.py"
    if file_is_empty(df_GWAS_filt_all_file) or True:

        think_how_to_proceed_here

        # print parameters working well on flc
        print("printing FLC parameters")
        GWAS_benchmarking_df = load_object(GWAS_benchmarking_df_file).reset_index(drop=True)
        GWAS_benchmarking_df_filt = get_filtered_df_by_dict(GWAS_benchmarking_df, {"phenotype":"FLC_resistance"})
        GWAS_benchmarking_df = GWAS_benchmarking_df[(GWAS_benchmarking_df.min_epsilon>0) & (GWAS_benchmarking_df.min_support>=70) & (GWAS_benchmarking_df.min_nodes_GenoAndPheno==2)].copy()
        pvalues_dict_filt = {"pval_chi_square_phenotypes_maxT":True, "pval_epsilon_phenotypes_maxT":False, "pval_chi_square_phenotypes":False, "pval_GenoAndPheno_phenotypes":False, "pval_fisher":False, "correction_method":"none"}
        GWAS_benchmarking_df = get_filtered_df_by_dict(GWAS_benchmarking_df, pvalues_dict_filt)
        GWAS_benchmarking_df = GWAS_benchmarking_df[GWAS_benchmarking_df.set_expected_genes_found.apply(lambda x: "MRR1" in x)].copy()

        print_fields = ["sample_set", "aligner", "ASR_method", "min_support", "min_epsilon", "min_nodes_GenoAndPheno", "correction_all_groups", "expected_genes_found_str"]
        if len(GWAS_benchmarking_df)!=len(GWAS_benchmarking_df[print_fields].drop_duplicates()): raise ValueError("invalid numbers")
        df_print = GWAS_benchmarking_df.sort_values(by=print_fields)[print_fields + ["n_sig_genes", "n_sig_variants"]].reset_index(drop=True)
        print(len(GWAS_benchmarking_df), df_print)
        df_print.to_excel(PlotsDir + "/parameters_MRR1_flc.xlsx")

        # define the best filters (manual curation)
        print("getting best filters")
        all_possible_fields_GWAS = set(get_tab_as_df_or_empty_df("%s/pheno-FLC_resistance_aligner-bwa_mem/convergence_gwas/gwas_results_representative_samples_min_sup_70_ASR_DOWNPASS/integrated_GWAS_stats.tab"%(GenoPheno_dir_all)).keys())
        pvalues_dict = { "pval_chi_square_phenotypes_maxT":True, "pval_chi_square_phenotypes":False, "pval_GenoAndPheno_phenotypes":False, "pval_fisher":False, "pval_epsilon_phenotypes_maxT":False, "correction_method":"none", "correction_all_groups":False}
        tuple_pval_fields = tuple(get_list_pval_fields_one_combination_pvals(pvalues_dict, all_possible_fields_GWAS))               
        best_filters_dict = {"tuple_pval_fields":tuple_pval_fields, "min_epsilon":0.1, "min_nodes_GenoAndPheno":2, "ASR_method":"MPPA,DOWNPASS", "min_support":70, "aligner":"bwa_mem", "sample_set":"representative_samples"}

        # print the effect of these filters
        print("printing filter effects...")
        GWAS_benchmarking_df_dict = {}
        for f in ["min_epsilon", "min_nodes_GenoAndPheno", "ASR_method", "min_support", "aligner", "sample_set"]:
            GWAS_benchmarking_df_dict[f] = best_filters_dict[f]

        for k,v in pvalues_dict.items(): GWAS_benchmarking_df_dict[k] = v
        GWAS_benchmarking_df = load_object(GWAS_benchmarking_df_file).reset_index(drop=True)
        GWAS_benchmarking_df_filt = get_filtered_df_by_dict(GWAS_benchmarking_df, GWAS_benchmarking_df_dict)
        print(GWAS_benchmarking_df_filt[["phenotype", "expected_genes_found_str", "n_sig_genes"]])

        # get the filtered results for all phenotypes
        df_raw_gwas_files_all = get_filtered_df_by_dict(df_raw_gwas_files_all, {f : best_filters_dict[f] for f in ["ASR_method", "min_support", "aligner", "sample_set"]})
        if len(df_raw_gwas_files_all)!=len(set(df_raw_gwas_files_all.phenotype)): raise ValueError("pheno should be unique")
        print("GWAS done on %i phenotypes"%(len(df_raw_gwas_files_all)))

        gwas_filters_dict = {f : best_filters_dict[f] for f in ["tuple_pval_fields", "min_epsilon", "min_nodes_GenoAndPheno"]}
        inputs_fn = list(df_raw_gwas_files_all.apply(lambda r: (r.phenotype, r.raw_GWAS_stats, gwas_filters_dict), axis=1))
        
        print("get the filtered GWAS")
        with multiproc.Pool(112) as pool:
            df_GWAS_filt_all = pd.concat(pool.starmap(get_filtered_GWAS_dfs_one_phenotype, inputs_fn, chunksize=1))
            pool.close()
            pool.terminate()

        # checks
        if len(df_GWAS_filt_all)!=len(df_GWAS_filt_all[["phenotype", "unique_groupID"]].drop_duplicates()):
            raise ValueError("non unique ")

        print("GWAS hits on %i phenotypes"%(len(set(df_GWAS_filt_all.phenotype))))
        print(df_GWAS_filt_all.groupby("phenotype").apply(len).sort_values())

        save_object(df_GWAS_filt_all, df_GWAS_filt_all_file)

    return load_object(df_GWAS_filt_all_file)



def get_filtered_GWAS_dfs(GWAS_benchmarking_df_file, df_raw_gwas_files_all, GenoPheno_dir_all, ProcessedDataDir, DataDir, PlotsDir):

    """gets the single used GWAS dfs"""

    df_GWAS_filt_all_file = ProcessedDataDir + "/df_GWAS_filt_all.py"
    if file_is_empty(df_GWAS_filt_all_file):

        # define all of the GWAS filtered fields (OLD)
        #all_possible_fields_GWAS = set(get_tab_as_df_or_empty_df("%s/pheno-FLC_resistance_aligner-bwa_mem/convergence_gwas/gwas_results_representative_samples_min_sup_70_ASR_DOWNPASS/integrated_GWAS_stats.tab"%(GenoPheno_dir_all)).keys())
        #print(sorted([f for f in all_possible_fields_GWAS if "pval" in f and "bonferroni" in f]))

        # ['pval_GenoAndPheno_phenotypes_bonferroni', 'pval_GenoAndPheno_phenotypes_bonferroni_allGroups', 'pval_OR_phenotypes_bonferroni', 'pval_OR_phenotypes_bonferroni_allGroups', 'pval_chi_square_phenotypes_bonferroni', 'pval_chi_square_phenotypes_bonferroni_allGroups', 'pval_epsilon_phenotypes_bonferroni', 'pval_epsilon_phenotypes_bonferroni_allGroups', 'pval_fisher_bonferroni', 'pval_fisher_bonferroni_allGroups']


        # get the two types of hits
        df_GWAS_filt_all = pd.DataFrame()
        for type_filtering, filt_parms in [("low_confidence", manual_curation_gwas_filt_parms_low_conf), ("high_confidence", manual_curation_gwas_filt_parms)]:
            print("getting best filters", type_filtering)

            # define the p value tuple
            if filt_parms["type pvalues general"]=="only_p_GP_bonf" and filt_parms["correction_all_groups"] is False:
                tuple_pval_fields = ("pval_GenoAndPheno_phenotypes_bonferroni",)

            elif filt_parms["type pvalues general"]=="only_p_GP<0.05" and filt_parms["correction_all_groups"] is False:
                tuple_pval_fields = ("pval_GenoAndPheno_phenotypes",)

            else: raise ValueError("invalid %s"%filt_parms)

            # define the filters to use
            filt_parms_use = {k:v for k,v in filt_parms.items() if k!="type pvalues general"}
            filt_parms_use["tuple_pval_fields"] = tuple_pval_fields

            # get the filtered results for all phenotypes
            df_raw_gwas_files_filt = get_filtered_df_by_dict(df_raw_gwas_files_all, {f : filt_parms_use[f] for f in ["ASR_method", "min_support", "aligner", "sample_set"]})
            if len(df_raw_gwas_files_filt)!=len(set(df_raw_gwas_files_filt.phenotype)): raise ValueError("pheno should be unique")
            print("GWAS done on %i phenotypes"%(len(df_raw_gwas_files_filt)))

            gwas_filters_dict = {f : filt_parms_use[f] for f in ["tuple_pval_fields", "min_epsilon", "min_nodes_GenoAndPheno", "pvalue_threshold"]}
            inputs_fn = list(df_raw_gwas_files_filt.apply(lambda r: (r.phenotype, r.raw_GWAS_stats, gwas_filters_dict), axis=1))
            
            print("get the filtered GWAS")
            with multiproc.Pool(112) as pool:
                df_GWAS_filt = pd.concat(pool.starmap(get_filtered_GWAS_dfs_one_phenotype, inputs_fn, chunksize=1))
                pool.close()
                pool.terminate()

            # checks
            if len(df_GWAS_filt)!=len(df_GWAS_filt[["phenotype", "unique_groupID"]].drop_duplicates()):
                raise ValueError("non unique ")


            print("GWAS hits on %i phenotypes"%(len(set(df_GWAS_filt.phenotype))))
            print(df_GWAS_filt.groupby("phenotype").apply(len).sort_values())

            # append
            df_GWAS_filt["type_filtering"] = type_filtering
            df_GWAS_filt_all = pd.concat([df_GWAS_filt_all, df_GWAS_filt]).reset_index(drop=True)


        save_object(df_GWAS_filt_all, df_GWAS_filt_all_file)

    return df_GWAS_filt_all_file




def get_NR_GWAS_tables(df_GWAS_filt_all_file, PlotsDir, DataDir, GenoPheno_dir_all, list_cmds, ProcessedDataDir):

    """get NR GWAS tables"""

    # load df
    df_GWAS_filt_all = load_object(df_GWAS_filt_all_file).sort_values(by=["type_filtering", "phenotype"])

    # define general dirs
    ancestral_GWAS_dir = "%s/scripts/ancestral_GWAS/scripts"%ParentDir
    pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
    interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir
    CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
    genome_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir # genome from CGD
    gff_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir
    mapping_genes = "Scer_HAL9:TAC1" 
    reactome_organism_priority = "Saccharomyces_cerevisiae,Schizosaccharomyces_pombe"

    # go through each phenotype
    df_NR_hits_all = pd.DataFrame()
    for (type_filtering, phenotype), df_GWAS_filt in df_GWAS_filt_all.groupby(["type_filtering", "phenotype"]):
        print("NR hits", phenotype, type_filtering)

        #if type_filtering!="low_confidence": continue # debug!!!!
        #if phenotype!="CAS_resistance": continue # debug
        #if phenotype!="5FC_resistance": continue # debug

        # define outdir
        outdir_all = PlotsDir + "/nr_GWAS_tables_%s"%type_filtering; make_folder(outdir_all)
        outdir_pheno = outdir_all + "/" + phenotype; make_folder(outdir_pheno)

        # get the GWAS into a file
        pval_fields = ['pval_GenoAndPheno_phenotypes', 'pval_GenoAndPheno_phenotypes_bonferroni', 'pval_GenoAndPheno_phenotypes_bonferroni_allGroups', "pval_chi_square_phenotypes_maxT", "pval_chi_square_phenotypes_maxT_allGroups"]
        #target_fields = ['nodes_withGeno', 'nodes_withoutGeno', 'nodes_unkownGeno', 'nodes_withPheno', 'nodes_withoutPheno', 'nodes_unkownPheno', 'nodes_withGeno_clearNodes', 'nodes_withoutGeno_clearNodes', 'nodes_withPheno_clearNodes', 'nodes_withoutPheno_clearNodes', 'nodes_GenoAndPheno', 'nodes_noGenoAndPheno', 'nodes_GenoAndNoPheno', 'nodes_noGenoAndNoPheno', 'group_name', 'ASR_method', 'classification_method_ASR', 'epsilon', 'OR', 'chi_square', 'grouping_ID', 'group', 'unique_groupID', 'groupName_gwas', 'linkage_groupID', 'type_vars', 'type_genes', 'type_mutations', 'type_collapsing', 'min_support', 'min_n_geno_transitions', 'min_n_geno_gains', 'min_n_geno', 'phenotype', 'node_names_GenoAndPheno', 'cladeID_names_GenoAndPheno'] + pval_fields
        target_fields = ['nodes_GenoAndPheno', 'nodes_noGenoAndPheno', 'nodes_GenoAndNoPheno', 'nodes_noGenoAndNoPheno', 'nodes_withPheno', 'epsilon', 'OR', 'chi_square', 'unique_groupID', 'type_vars', 'type_genes', 'type_mutations', 'type_collapsing', 'node_names_GenoAndPheno'] + pval_fields

        gwas_hits_file = "%s/all_gwas_hits.tab"%outdir_pheno
        if file_is_empty(gwas_hits_file) or True:
            save_df_as_tab(df_GWAS_filt[target_fields], gwas_hits_file)

        # define files
        filters_dict = {"high_confidence":manual_curation_gwas_filt_parms, "low_confidence":manual_curation_gwas_filt_parms_low_conf}[type_filtering]
        gwas_inputs_dir = "%s/pheno-%s_aligner-%s/convergence_gwas/gwas_inputs_%s"%(GenoPheno_dir_all, phenotype, filters_dict["aligner"], filters_dict["sample_set"])

        # define the sorting fields
        initial_sorting_fields = ",".join(["epsilon", "OR", "nodes_GenoAndPheno"])
        initial_sorting_ascending = ",".join(map(str, [False, False, False]))
        # "pval_GenoAndPheno_phenotypes", "pval_chi_square_phenotypes_maxT",

        # get the NR files
        outdir_NR_groups = outdir_pheno + "/NR_groups"
        if file_is_empty(outdir_NR_groups + "/pipeline_finished.txt"):
            cmd = f'{ancestral_GWAS_dir}/generate_NR_groups.py --outdir {outdir_NR_groups} --gwas_hits {gwas_hits_file} --gff {gff_Cparapsilosis} --ref_genome {genome_Cparapsilosis} --custom_geneID_maps {mapping_genes} --gwas_inputs_dir {gwas_inputs_dir} --pathway_annotations_dir {pathway_annotations_dir} --CGD_chromosomal_features {CGD_chromosomal_features} --interproscan_output {interproscan_output} --initial_sorting_fields {initial_sorting_fields} --initial_sorting_ascending {initial_sorting_ascending} --reactome_organism_priority {reactome_organism_priority} --keep_variants_w_transitions --min_support {filters_dict["min_support"]}  --ASR_method {filters_dict["ASR_method"]}'

            list_cmds.append(("gwas_NR_reduction", cmd, "ancestral_GWAS_env", 112))
            continue

        # add to df
        df_NR_hits = get_tab_as_df_or_empty_df(outdir_NR_groups+"/NR_hits.tab")
        df_NR_hits["phenotype"] = phenotype
        df_NR_hits["type_filtering"] = type_filtering

        df_NR_hits_all = pd.concat([df_NR_hits_all, df_NR_hits]).reset_index(drop=True)

    # save
    if len(list_cmds)==0:
        print("Saving NR GWAS hits...")

        # change pheno
        df_NR_hits_all["phenotype"] = df_NR_hits_all.phenotype.apply(get_paper_dicotomic_pheno)

        # save
        save_object(df_NR_hits_all, ProcessedDataDir+"/all_GWAS_hits_NR.py")
        
        for type_filtering, df_NR_hits in df_NR_hits_all.groupby("type_filtering"):

            save_df_as_tab(df_NR_hits.copy(), "%s/GWAS_hits_NR_%s.tab"%(PlotsDir, type_filtering))
            df_NR_hits.copy().to_excel("%s/GWAS_hits_NR_%s.xlsx"%(PlotsDir, type_filtering))


def print_GWAS_hits_gene_info(df_GWAS_hits_NR):

    """Prints certain GWAS hits"""

    # edit
    df_GWAS_hits_NR = df_GWAS_hits_NR.copy()
    df_highc = df_GWAS_hits_NR[df_GWAS_hits_NR.type_filtering=="high_confidence"].copy()
    df_lowc = df_GWAS_hits_NR[df_GWAS_hits_NR.type_filtering=="low_confidence"].copy()

    # prints
    gname_to_phenotypes = {"CPAR2_806400": ["fluconazole R", "voriconazole R", "posaconazole R"],
                           "Scer_NCP1": ["fluconazole R", "voriconazole R", "posaconazole R", "amphotericin B R"],
                           "Scer_TEL1": ["fluconazole R", "voriconazole R", "posaconazole R"],
                           "CDR1B": ["fluconazole R", "voriconazole R", "chlorhexidine tol"],
                           "Scer_UTP22": ["fluconazole R", "voriconazole R", "sterile"],
                           "CPAR2_204120": ["voriconazole R"],
                           "Scer_UGA2": ["fluconazole R", "voriconazole R"],
                           "CPAR2_105460": ["fluconazole R", "voriconazole R"],
                           "CPAR2_105480": ["fluconazole R", "voriconazole R"],
                           "ERG25": ["fluconazole R"],
                           "Scer_ATP10": ["fluconazole R"],
                           "Scer_ERG24": ["fluconazole R"],
                           "ERG11": ["fluconazole R", "voriconazole R", "posaconazole R"],
                           "MRR1": ["fluconazole R", "voriconazole R"],
                           "NDT80": ["fluconazole R", "voriconazole R"],
                           "Scer_CLU1": ["voriconazole R"],

                           "Scer_MDS3": ["blood vs other sterile"],
                           "ALS7": ["voriconazole R", "blood culture"],
                           "Scer_NUP1": ["voriconazole R"],
                           "Scer_NST1": ["posaconazole R"],
                           "CPAR2_806400": ["amphotericin B R", "chlorhexidine tol"],
                           "CPAR2_204210": ["amphotericin B R"],
                           "Scer_FOL1": ["amphotericin B R"],
                           "CPAR2_105750": ["5-FC R"],
                           "Scer_ENA2": ["5-FC R"],

                           "CPAR2_806390": ["5-FC R", "chlorhexidine tol"],
                           "Scer_MNR2": ["chlorhexidine tol"],
                           "Scer_RAD23": ["chlorhexidine tol"],
                           "CPAR2_101330": ["chlorhexidine tol"],
                           "Scer_ORC4": ["chlorhexidine tol"],

                           "Scer_NAT4": ["age old"],
                           "Scer_ATG20": ["age old"],
                           "CPAR2_403920": ["age old"],

                           "Scer_DFG10": ["skin exudate"],
                           "Scer_NTG2": ["bronchial aspirate"],
                           "CDR1": ["bronchial aspirate"],
                           "LPD1": ["sterile"]
                           }

    for gene_name, phenotypes in gname_to_phenotypes.items():
        if gene_name=="": continue
        print("\ngene: %s"%gene_name)

        for pheno in phenotypes:

            # get df
            df_target_high = df_highc[(df_highc.phenotype==pheno) & (df_highc.final_name==gene_name)].copy()
            df_target_low = df_lowc[(df_lowc.phenotype==pheno) & (df_lowc.final_name==gene_name)].copy()

            if len(df_target_high)>0: 
                df_target = df_target_high
                type_filt = "high_confidence"
            
            elif len(df_target_low)>0: 
                df_target = df_target_low
                type_filt = "low_confidence"
            
            else: raise ValueError("invalid %s-%s"%(gene_name, pheno))

            # print
            print("%s (%s): %s"%(pheno, type_filt, "; ".join(df_target["variant effect"])))











def get_GWAS_benchmarking_df_target_phenotypes_and_genes(GWAS_benchmarking_df, pheno_to_expected_genes):

    """Gets a GWAS filtered to give the expected genes in different phenotypes"""

    # keep the phenos
    target_phenotypes = set(pheno_to_expected_genes.keys())
    GWAS_benchmarking_df = GWAS_benchmarking_df[GWAS_benchmarking_df.phenotype.isin(target_phenotypes)].copy()

    # check all have
    if set(GWAS_benchmarking_df.phenotype)!=target_phenotypes:
        raise ValueError("all phenotypes should be there.")


    # redefine the genes
    GWAS_benchmarking_df["set_redefined_expected_genes"] = GWAS_benchmarking_df.phenotype.apply(lambda p: pheno_to_expected_genes[p])
    GWAS_benchmarking_df["all_redefined_expected_genes_found"] = (GWAS_benchmarking_df[["set_redefined_expected_genes", "set_expected_genes_found"]].apply(lambda r: r.set_redefined_expected_genes.difference(r.set_expected_genes_found), axis=1).apply(len))==0

    # define filters that work for all
    filterID_to_correctAllPhenos = GWAS_benchmarking_df[GWAS_benchmarking_df.all_redefined_expected_genes_found][["phenotype", "unique_filterID"]].drop_duplicates().groupby("unique_filterID").apply(lambda df: set(df.phenotype)==target_phenotypes)
    filterIDs_correctAllPhenos = set(filterID_to_correctAllPhenos[filterID_to_correctAllPhenos==True].index)


    # filter
    GWAS_benchmarking_df = GWAS_benchmarking_df[GWAS_benchmarking_df.unique_filterID.isin(filterIDs_correctAllPhenos)].copy()

    # log
    df_filt = GWAS_benchmarking_df[["pval_chi_square_phenotypes_maxT", "pval_epsilon_phenotypes_maxT", "unique_filterID"]].drop_duplicates()
    text_phenos = ", ".join(["%s %s"%(p, pheno_to_expected_genes[p]) for p in sorted(target_phenotypes)])
    print("%s: %i correct filters, %i maxT filters, aligners=%s, ASR=%s, min_support=%s"%(text_phenos, len(filterIDs_correctAllPhenos), sum((df_filt[["pval_chi_square_phenotypes_maxT", "pval_epsilon_phenotypes_maxT"]].apply(any, axis=1))), sorted(set(GWAS_benchmarking_df.aligner)), sorted(set(GWAS_benchmarking_df.ASR_method)), sorted(set(GWAS_benchmarking_df.min_support))))

def get_tuple_date_date_collection(x):

    """gets the tuple of date collection"""

    # '2021-03-18 00:00:00'
    if len(x.split())==2:
        year, month, day = list(map(int, x.split()[0].split("-")))

    # '03*04*2021'
    elif len(x.split())==1 and len(x.split("*"))==3:
        day, month, year = list(map(int, x.split("*")))

    else:
        raise ValueError("invalid %s"%x)

    # checks
    if not year in {2019, 2020, 2021, 2022, 2023}:
        raise ValueError("invalid year: %s"%year)

    if not month in set(range(1, 12+1)):
        raise ValueError("invalid month: %s"%month)

    if day<1 or day>31:
        raise ValueError("invalid day: %s"%day)

    return (year, month, day)

def get_days_divergence_from_tuples_dates(r):

    """Gets the days of divergence"""

    # sort the tuples
    t0_tuple, t1_tuple = sorted(r[["origin_runID_tuple_date_date_collection", "target_runID_tuple_date_date_collection"]])

    # add divergence
    time_divergence = (t1_tuple[0] - t0_tuple[0])*365
    time_divergence += (t1_tuple[1] - t0_tuple[1])*30
    time_divergence += (t1_tuple[2] - t0_tuple[2])*1

    return time_divergence

def plot_correlation_btw_time_and_genetic_divergence(df_metadata, ProcessedDataDir, aligner, PlotsDir, zygosity, trees_df):

    """Gets the correlation btw time and divergence"""

    # filter
    df_metadata = df_metadata[(df_metadata.source_data=="sequenced_here")].copy() # & (df_metadata.clade_Bergin2022!=-1)
    df_metadata["tuple_date_date_collection"] = df_metadata.date_collection.apply(get_tuple_date_date_collection)

    # add the clade 
    leaf_to_cladeID = trees_df[(trees_df.run_set=="exclude_triploids") & (trees_df.aligner==aligner) & (trees_df.type_outgroup=="Cparapsilosis") & (trees_df.mode_tree=="diploid")].iloc[0].leaf_to_clade_Bergin2022
    df_metadata["cladeID"] = df_metadata.runID.map(leaf_to_cladeID)
    check_no_nans_series(df_metadata["cladeID"])

    # filter on clade
    df_metadata = df_metadata[df_metadata.cladeID!=-1].copy()

    # checks
    if any(df_metadata.hospital=="unknown"):
        raise ValueError("unknown hospital")

    # load the divergence 
    window_size = 200000
    df_pairwise_distances_all = load_object("%s/dfs_stats_pairwise_SNP_distances_per_window/df_%s_%ibp_%s_SNP.py"%(ProcessedDataDir, aligner, window_size, zygosity)) 

    # filt divergence df
    all_runs = set(df_metadata.runID)
    df_pairwise_distances = df_pairwise_distances_all[(df_pairwise_distances_all.origin_runID.isin(all_runs)) & (df_pairwise_distances_all.target_runID.isin(all_runs)) & (df_pairwise_distances_all.type_var=="SNP")][["origin_runID", "target_runID", "whole_genome_diff_vars_per_kb"]].copy().reset_index(drop=True)

    # add fields
    for f in ["hospital", "tuple_date_date_collection", "cladeID"]:

        runID_to_fval = df_metadata.set_index("runID")[f]
        df_pairwise_distances["origin_runID_%s"%f] = df_pairwise_distances.origin_runID.map(runID_to_fval)
        df_pairwise_distances["target_runID_%s"%f] = df_pairwise_distances.target_runID.map(runID_to_fval)

    # keep those from the same hospital and clade
    df_pairwise_distances = df_pairwise_distances[(df_pairwise_distances.origin_runID_hospital == df_pairwise_distances.target_runID_hospital)].copy()
    df_pairwise_distances = df_pairwise_distances[(df_pairwise_distances.origin_runID_cladeID == df_pairwise_distances.target_runID_cladeID)].copy()

    # add the divergence
    df_pairwise_distances["days_divergence"] = df_pairwise_distances.apply(get_days_divergence_from_tuples_dates, axis=1)

    # get the spearman r
    spearmanr, spearmanp  = scipy.stats.spearmanr(df_pairwise_distances.whole_genome_diff_vars_per_kb.values, df_pairwise_distances.days_divergence.values, nan_policy="raise")


    # plot
    fig = plt.figure(figsize=(3.5, 3.5))
    ax = sns.scatterplot(data=df_pairwise_distances, x="whole_genome_diff_vars_per_kb", y="days_divergence", edgecolor="black", s=10, lw=.9)
    ax.set_title("pairwise %s SNP vs time distances, aligner=%s\nspearman r=%.2f, p=%.3f "%(zygosity, aligner, spearmanr, spearmanp), fontsize=10)

    # modify the face
    for artist in ax.collections: artist.set_facecolor('none')

    # keep
    save_fig(fig, "correlation_genetic_vs_time_diff_%s_%s.pdf"%(aligner, zygosity), PlotsDir)
    plt.show()


def plot_general_stats_SNP_distance_comparisons(ProcessedDataDir, aligner, window_size, df_metadata, PlotsDir, zygosity, yfield):

    """plots the pairwise differences stats between all samples"""

    # keep dfs
    print("loading data")
    df_metadata = df_metadata.copy()
    df_pairwise_distances_all = load_object("%s/dfs_stats_pairwise_SNP_distances_per_window/df_%s_%ibp_%s_SNP.py"%(ProcessedDataDir, aligner, window_size, zygosity)) 

    # filter runs
    print("filtering")
    relevant_runs = set(df_metadata.runID)
    df_plot = df_pairwise_distances_all[(df_pairwise_distances_all.target_runID.isin(relevant_runs)) & (df_pairwise_distances_all.origin_runID.isin(relevant_runs))].copy()

    # plot
    print("plotting")
    fig = plt.figure(figsize=(2, 2))
    ax = sns.scatterplot(data=df_plot, x="whole_genome_diff_vars_per_kb", y=yfield, edgecolor="black", s=5, lw=.9)
    if window_size==1000: suffix = " (mtDNA)"
    else: suffix = ""
    ax.set_title("pairwise %s SNP distances\naligner=%s, window=%ikb%s"%(zygosity, aligner, window_size/1000, suffix), fontsize=10)
    #ax = sns.scatterplot(data=df_plot, x=xfield, y=yfield, hue=hue, palette=palette_dict, linewidth=.8, edgecolor=[palette_dict[x] for x in df_plot[hue]])

    # threshold
    plt.axvline(float(yfield.split("_")[-1]), lw=1, ls="--", c="red")
    # ID function
    #max_x = max(df_plot.whole_genome_diff_vars_per_kb)
    #plt.plot([0, max_x], [0, max_x], lw=1, ls="--", c="red")


    # modify the face
    for artist in ax.collections: artist.set_facecolor('none')

    sns.despine()

    # keep
    save_fig(fig, "pairwise_SNP_distances/%s_%s_wsize_%i_%s.png"%(yfield, aligner, window_size, zygosity), PlotsDir)
    plt.show()

def get_ML_df_with_test_parmID(df_ML, extra_fields, test_fields):

    """adds the parameters ID"""


    # checks
    if len(df_ML)!=len(df_ML[test_fields + extra_fields].drop_duplicates()):
        raise ValueError("not expected fields")

    # add
    df_ML["test_parmID"] = df_ML[test_fields].apply(lambda r: "---".join(list(map(str, r.values))), axis=1)

    return df_ML

def get_ML_integrated_performance_df(df_ML_results_all, field_AUC_train, field_AUC_test, min_nmodels):

    """gets integrated performance across all"""

    # load df
    print("getting integrated performance df...")
    df_ML_results_all = df_ML_results_all.copy()

    # check
    group_fields = ["phenotype", "test_parmID"]
    unique_fields = group_fields + ["test_split_ID"]
    if len(df_ML_results_all)!=len(df_ML_results_all[unique_fields].drop_duplicates()):
        raise ValueError("not expected fields")

    # get df processed
    def get_series_integrated_ML_performance_info(df_all, group_fields):

        # checks
        if len(df_all)!=5: raise ValueError("invalid dfs")
        if len(set(df_all[field_AUC_test]))!=1: raise ValueError("invalid field_AUC_test")

        # init fields dict
        fields_dict = dict(df_all.iloc[0][group_fields])
        fields_dict["AUC_test"] = df_all.iloc[0][field_AUC_test] 
        fields_dict["tuple_expected_genes"] = tuple(sorted(set.union(*df_all.expected_genes_str.apply(lambda x: set(x.split(",")))).difference({""})))
        fields_dict["n_models_fit"] = sum(df_all.selected_predictors.apply(len)>0)

        # add idxs
        sorted_expected_genes = [('ERG11', 'ERG11;MRR1'), ('ERG11', 'ERG11;Scer_HAL9'), ('ERG11', 'Scer_HAL9'), ('ERG11',), ('Scer_HAL9',), ('FKS1',), tuple([])]
        exp_genes_idx_dict = dict(zip(sorted_expected_genes, range(len(sorted_expected_genes))))
        fields_dict["expected_genes_idx"] = exp_genes_idx_dict[fields_dict["tuple_expected_genes"]]
        fields_dict["all_models_found_idx"] = {True:0, False:1}[all(df_all.model_name!="none") and all(df_all.selected_predictors.apply(len)>0)]

        # filter to keep only models with predictors
        df = df_all[df_all.selected_predictors.apply(len)>0].copy()
        if len(df)>=min_nmodels:
            fields_dict["mean_AUC_train"] = np.mean(df[field_AUC_train])
            fields_dict["min_n_predictors"] = min(df.selected_predictors.apply(len))

        else:
            fields_dict["mean_AUC_train"] = 0.5
            fields_dict["min_n_predictors"] = 0

        return pd.Series(fields_dict)

    df_plot_all = df_ML_results_all.groupby(group_fields).apply(get_series_integrated_ML_performance_info, group_fields=group_fields).reset_index(drop=True)

    # edit AUC
    for f in ["AUC_test", "mean_AUC_train"]:
        df_plot_all[f] = df_plot_all[f].apply(lambda x: max([x, 0.5]))

    return df_plot_all

def get_df_like_first_row(df, sorting_fields, ascending_fields):

    """gets a df that is like the first row, according to sorting_fields"""

    # sort df
    if len(sorting_fields)==0: raise ValueError("should be some ascending fields.")
    check_no_nans_in_df(df[sorting_fields])
    df = df.sort_values(by=sorting_fields, ascending=ascending_fields).copy()

    # get df like the first row
    first_r = df.iloc[0]
    df = df[df[sorting_fields].apply(lambda r: all(map(lambda f: r[f]==first_r[f], sorting_fields)), axis=1)].copy()

    return df


def get_df_ML_results_best_each_phenotype(df_ML_performance, df_ML_results_all, threshold_ROC_AUC):

    """Gets best performance"""

    # keep only all_samples
    #df_ML_performance = df_ML_performance[df_ML_performance.sample_set=="all_samples"].copy()

    # keep
    df_ML_performance = df_ML_performance[(df_ML_performance.AUC_test>=threshold_ROC_AUC) & (df_ML_performance.mean_AUC_train>=threshold_ROC_AUC)].copy() # (df_ML_performance.all_models_found_idx==0)
    df_ML_results_all = df_ML_results_all[df_ML_results_all.phenotype.isin(set(df_ML_performance.phenotype))].copy()

    # add the parm ID
    parm_fields = ["type_phenotype", "CVset_type_split", "consistency_measure", "multiple_test_corr", "sample_set"]
    for df in [df_ML_performance, df_ML_results_all]:
        df["parm_id"] = df[parm_fields].apply(lambda r: "|".join(r), axis=1)

    if len(df_ML_performance)!=len(df_ML_performance[["parm_id", "phenotype"]].drop_duplicates()):
        raise ValueError("non unique df")

    # add fields
    df_ML_performance["consistency_measure_idx"] = df_ML_performance.consistency_measure.map({"jaccard_variants":0, "jaccard_samples":1, "none":2})
    check_no_nans_series(df_ML_performance["consistency_measure_idx"])

    df_ML_performance["sample_set_idx"] = df_ML_performance.sample_set.map({"all_samples":0, "representative_samples":1, "balanced_samples":2})
    check_no_nans_series(df_ML_performance["sample_set_idx"])

    # get best parmID
    def get_best_parmID_one_pheno(df):

        # models with high AUC
        min_AUC_select = max(df.AUC_test) - 0.05
        df_filt = df[df.AUC_test>=min_AUC_select].copy()

        # further filters
        df_filt = get_df_like_first_row(df_filt, ["expected_genes_idx", "all_models_found_idx", "sample_set_idx", "AUC_test", "mean_AUC_train", "consistency_measure_idx"], [True, True, True, False, False, True])
        if len(df_filt)!=1: raise ValueError("should be one row. %s"%list(df_filt.parm_id))
        return df_filt.iloc[0].parm_id
     
    pheno_to_best_parm = dict(df_ML_performance.groupby("phenotype").apply(get_best_parmID_one_pheno))

    # keep
    return df_ML_results_all[df_ML_results_all.apply(lambda r: r.parm_id == pheno_to_best_parm[r.phenotype], axis=1)].copy()

def general_ML_performance_plot(df_plot_all, PlotsDir, threshold_AUC=0.75):

    """Visualize the performances"""

    df_plot_all = df_plot_all.copy()

    # get plot
    sorted_rows = sorted(set(df_plot_all.type_phenotype))
    sorted_cols = sorted(set(df_plot_all.consistency_measure))
    df_plot_all = df_plot_all.set_index(["type_phenotype", "consistency_measure"])

    nrows = len(sorted_rows)
    ncols = len(sorted_cols)
    fig = plt.figure(figsize=(ncols*2, nrows*2)); I=1

    for Ir, type_phenotype in enumerate(sorted_rows):
        for Ic, consistency_measure in enumerate(sorted_cols):

            ax = plt.subplot(nrows, ncols, I); I+=1
            df_plot = df_plot_all.loc[(type_phenotype, consistency_measure)].copy().sort_values(by=["CVset_type_split", "multiple_test_corr", "phenotype"])
            ax = sns.scatterplot(data=df_plot, x="mean_AUC_train", y="AUC_test", hue="CVset_type_split", style="multiple_test_corr", palette={"random":"gray", "balanced":"salmon"}, linewidth=.7, edgecolor="black")

            #spearmanr, spearmanp  = scipy.stats.spearmanr(df_plot.mean_AUC_train.values, df_plot.AUC_test.values, nan_policy="raise")
            #ax.text(0.5, 0.97, "r=%.2f"%spearmanr, fontsize=10)

            lims = [0.45, 1.05]
            ax.set_xlim(lims)
            ax.set_ylim(lims)

            plt.plot(lims, lims, c="black", lw=.5, ls="--")
            ax.axvline(threshold_AUC, c="black", lw=.5, ls="--")
            ax.axhline(threshold_AUC, c="black", lw=.5, ls="--")

            if Ir==0 and Ic==(len(sorted_cols)-1):
                ax.legend(bbox_to_anchor=(1,1), loc="upper left")

            else:
                ax.get_legend().remove()

            if Ic==0:
                ax.set_ylabel("%s\ntest $AUC_{ROC}$"%(type_phenotype))

            else:
                ax.set_ylabel("")
                ax.set_yticklabels([])

            if Ir==1 and Ic==2:
                ax.set_xlabel("average train $AUC_{ROC}$")

            elif Ir==1:
                ax.set_xlabel("")

            else:
                ax.set_xlabel("")
                ax.set_xticklabels([])

            if Ir==0:
                ax.set_title("consistency by\n" + consistency_measure, fontsize=10)



    plt.subplots_adjust(wspace=0.05, hspace=0.05)# wspace (between cols) hspaces (btw rows)
    plt.show()

def get_ML_sig_phenotypes(df_ML_performance, threshold_AUC):

    """Get sig phenos sorted"""

    # get phenos
    phenotypes_sig = set(df_ML_performance[(df_ML_performance.AUC_test>=threshold_AUC) & (df_ML_performance.mean_AUC_train>=threshold_AUC)].phenotype)

    # sort
    df_ML_performance = df_ML_performance[df_ML_performance.phenotype.isin(phenotypes_sig)].copy()
    #return list(df_ML_performance.groupby("phenotype").apply(lambda df: sum(df.AUC_test>=threshold_AUC)).sort_values(ascending=False).index)
    return list(df_ML_performance.groupby("phenotype").apply(lambda df: max(df.AUC_test)).sort_values(ascending=False).index)

def plot_effect_sample_set(df_ML_performance, PlotsDir, phenotype_to_color, significant_ML_phenotypes, threshold_ROC_AUC):
    
    """Plots the effect of the sample set"""

    # keep
    df_plot_all = df_ML_performance.copy()

    # checks
    unique_fields = ['type_phenotype', 'consistency_measure', 'CVset_type_split', "sample_set", "multiple_test_corr", "phenotype"]
    if len(df_plot_all)!=len(df_plot_all[unique_fields].drop_duplicates()):
        raise ValueError("non unique fields")

    # add parmID
    df_plot_all["parmID"] = df_plot_all[["multiple_test_corr", "phenotype", "CVset_type_split"]].apply(lambda r: "|".join(r.values), axis=1)

    # add threshold
    high_AUC_train_f = "average train $AUC_{ROC}>=%s$"%threshold_ROC_AUC
    df_plot_all[high_AUC_train_f] = df_plot_all.mean_AUC_train>=threshold_ROC_AUC

    # get plot
    sorted_rows = sorted(set(df_plot_all.type_phenotype))
    sorted_cols = sorted(set(df_plot_all.consistency_measure))
    df_plot_all = df_plot_all.set_index(["type_phenotype", "consistency_measure"])

    nrows = len(sorted_rows)
    ncols = len(sorted_cols)
    fig = plt.figure(figsize=(ncols*2.1, nrows*2.1)); I=1

    for Ir, type_phenotype in enumerate(sorted_rows):
        for Ic, consistency_measure in enumerate(sorted_cols):

            # get plot long
            ax = plt.subplot(nrows, ncols, I); I+=1
            df_plot_long = df_plot_all.loc[(type_phenotype, consistency_measure)].copy()
            if len(df_plot_long)!=len(df_plot_long[["multiple_test_corr", "phenotype", "CVset_type_split", "sample_set"]].drop_duplicates()): raise ValueError("invalid")
            pID_to_test_AUC = df_plot_long[df_plot_long.sample_set=="all_samples"].set_index("parmID").AUC_test
            pID_to_train_AUC = df_plot_long[df_plot_long.sample_set=="all_samples"].set_index("parmID").mean_AUC_train

            # get df
            df_plot = pd.DataFrame()
            for ss in ["balanced_samples", "representative_samples"]:
                df_plot_ss = df_plot_long[df_plot_long.sample_set==ss][["AUC_test", "parmID", "mean_AUC_train", high_AUC_train_f]].copy()
                df_plot_ss["sample_set"] = ss
                df_plot_ss["AUC_test_no_reduction"] = df_plot_ss.parmID.map(pID_to_test_AUC)
                df_plot_ss["AUC_train_no_reduction"] = df_plot_ss.parmID.map(pID_to_train_AUC)

                check_no_nans_series(df_plot_ss["AUC_test_no_reduction"])
                df_plot = pd.concat([df_plot, df_plot_ss])

            # filter out cases in which a higher AUC_test with reduction makes no sense (i.e. the AUC train is low)
            #idx_bad_rows = (df_plot.mean_AUC_train<threshold_ROC_AUC) | (df_plot.AUC_train_no_reduction<threshold_ROC_AUC)
            #df_plot = df_plot[~idx_bad_rows].copy()

            # plot
            ax = sns.scatterplot(data=df_plot.sort_values(by=["sample_set", high_AUC_train_f]), x="AUC_test_no_reduction", y="AUC_test", hue="sample_set", linewidth=.7, edgecolor="black", palette={"representative_samples":"salmon", "balanced_samples":"c"}, style=high_AUC_train_f)

            lims = [0.45, 1.05]
            ax.set_xlim(lims)
            ax.set_ylim(lims)

            plt.plot(lims, lims, c="gray", lw=.5, ls="--")
            for tval in sorted(set([threshold_ROC_AUC])):
                ax.axvline(tval, c="black", lw=.5, ls="--")
                ax.axhline(tval, c="black", lw=.5, ls="--")

            if Ir==0 and Ic==(len(sorted_cols)-1):
                ax.legend(bbox_to_anchor=(1.15,1), loc="upper left")

            else:
                ax.get_legend().remove()

            if Ic==0:
                ax.set_ylabel("test $AUC_{ROC}$")

            else:
                ax.set_ylabel("")
                ax.set_yticklabels([])

            if Ic==(len(sorted_cols)-1):
                ax2 = ax.twinx()
                ax2.set_ylabel("%s"%(type_phenotype))
                ax2.set_yticks([])
                ax2.set_yticklabels([])

            if Ir==(len(sorted_rows)-1) and Ic==1:
                ax.set_xlabel("test $AUC_{ROC}$ all_samples")

            elif Ir==(len(sorted_rows)-1):
                ax.set_xlabel("")

            else:
                ax.set_xlabel("")
                ax.set_xticklabels([])

            if Ir==0:
                ax.set_title("consistency by:\n" + consistency_measure, fontsize=10)



    plt.subplots_adjust(wspace=0.05, hspace=0.05)# wspace (between cols) hspaces (btw rows)
    plt.show()
    save_fig(fig, "ML_model_results_check_sample_set.pdf", PlotsDir)





def visualize_ML_parameters_results_heatmap(df_ML_performance, PlotsDir, threshold_AUC, significant_ML_phenotypes, df_ML_results_best):

    """Perform heatmap with the parameters"""

    # filter
    df_ML_performance = df_ML_performance[df_ML_performance.phenotype.isin(significant_ML_phenotypes)].copy()

    # map each pheno to the best parm
    phenotype_to_best_parm_id = df_ML_results_best.groupby("phenotype").apply(lambda df: sorted(set(df.parm_id))[0])

    # add the row ID
    row_fields = ["type_phenotype", "CVset_type_split", "consistency_measure", "multiple_test_corr", "sample_set"]
    df_ML_performance["row_id"] = df_ML_performance[row_fields].apply(lambda r: "|".join(r), axis=1)

    if len(df_ML_performance)!=len(df_ML_performance[["row_id", "phenotype"]].drop_duplicates()):
        raise ValueError("non unique df")

    # define the sorted phenotypes
    sorted_rows = sorted(set(df_ML_performance["row_id"]))
    sorted_phenotypes = cp.deepcopy(significant_ML_phenotypes)

    # get the square
    df_plot = df_ML_performance.pivot(index="row_id", columns="phenotype", values="AUC_test")[sorted_phenotypes].loc[sorted_rows].copy()
    check_no_nans_in_df(df_plot)

    # get the df
    f_to_color_dict = {"type_phenotype": {"phenotype":"red", "phenotype_transition":"magenta"},
                       "CVset_type_split": {"balanced":"olive", "random":"gray"},
                       "consistency_measure": {"jaccard_samples":"blue", "jaccard_variants":"cyan", "none":"gray"},
                       "multiple_test_corr": {"maxT":"salmon", "none":"gray"},
                       "sample_set":{"representative_samples":"salmon", "balanced_samples":"c", "all_samples":"black"}}

    row_colors_dict = {}
    for row_val in sorted_rows:
        f_to_val = dict(zip(row_fields, row_val.split("|")))
        row_colors_dict[row_val] = {f : f_to_color_dict[f][f_to_val[f]] for f in row_fields}

    row_colors_df = pd.DataFrame(row_colors_dict).transpose()[row_fields].loc[sorted_rows]

    # get annot df
    annot_df = pd.DataFrame(index=df_plot.index)

    def get_annot_val(pID, best_pID, auc_val):
        
        if pID in best_pID: suff = "*"
        else: suff = ""

        if auc_val>=threshold_AUC: return "." + ("%.2f"%auc_val).split(".")[1] + suff
        else: return ""

    for pheno in sorted_phenotypes:
        annot_df[pheno] = list(map(lambda x: get_annot_val(x, phenotype_to_best_parm_id[pheno], df_plot.loc[x, pheno]), annot_df.index))

    # get the clustermap
    cm = sns.clustermap(df_plot, col_cluster=False, row_cluster=True, row_colors=row_colors_df, col_colors=None, cmap="Blues", vmin=0.5, vmax=1.0, xticklabels=True, cbar_kws={"label":"test $AUC_{ROC}$", "orientation":'vertical'}, yticklabels=1, linewidth=.5, square=True, annot=annot_df, fmt="", figsize=(12, 18)) #, linewidth=.01) # figsize=(30, 10),

    # annot_kws={"size": 10},
    plt.show()
    save_fig(cm, "clustermap_models_sig_phenotypes.pdf", PlotsDir)




def general_ML_performance_plot_split_old(df_plot_all, PlotsDir, phenotype_to_color, sorted_phenotypes, threshold_AUC=0.75):

    """Visualize the performances"""

    df_plot_all = df_plot_all.copy()

    # def visuals
    phenotype_to_color = cp.deepcopy(phenotype_to_color)
    missing_ps = set(sorted_phenotypes).difference(set(phenotype_to_color))
    if missing_ps!=set():
        raise ValueError("missing phenos: %s"%missing_ps)

    phenotype_to_color["other phenotypes"] = "gray"

    # add fields
    def get_pheno_val(x):
        if x in sorted_phenotypes: return x
        else: return "other phenotypes"

    df_plot_all["pheno"] = df_plot_all.phenotype.apply(get_pheno_val)
    df_plot_all["pheno_and_split"] = df_plot_all.type_phenotype + "|||" + df_plot_all.CVset_type_split

    pheno_to_max_AUC = df_plot_all.groupby("pheno").apply(lambda df: max(df.AUC_test)).sort_values(ascending=False)

    df_plot_all["pheno_idx"] = df_plot_all.pheno.map(dict(zip(pheno_to_max_AUC.index, range(len(pheno_to_max_AUC)))))

    # get plot
    sorted_rows = sorted(set(df_plot_all.pheno_and_split))
    sorted_cols = sorted(set(df_plot_all.consistency_measure))
    df_plot_all = df_plot_all.set_index(["pheno_and_split", "consistency_measure"])

    nrows = len(sorted_rows)
    ncols = len(sorted_cols)
    fig = plt.figure(figsize=(ncols*2.1, nrows*2.1)); I=1

    for Ir, pheno_and_split in enumerate(sorted_rows):
        type_phenotype, CVset_type_split = pheno_and_split.split("|||")
        for Ic, consistency_measure in enumerate(sorted_cols):

            ax = plt.subplot(nrows, ncols, I); I+=1
            df_plot = df_plot_all.loc[(pheno_and_split, consistency_measure)].copy().sort_values(by=["pheno_idx", "multiple_test_corr"])

            if len(df_plot)!=len(df_plot[["multiple_test_corr", "phenotype"]].drop_duplicates()): raise ValueError("invalid")

            ax = sns.scatterplot(data=df_plot, x="mean_AUC_train", y="AUC_test", hue="pheno", style="multiple_test_corr", linewidth=.7, edgecolor="black", palette=phenotype_to_color)

            lims = [0.45, 1.05]
            ax.set_xlim(lims)
            ax.set_ylim(lims)

            plt.plot(lims, lims, c="gray", lw=.5, ls="--")
            for tval in sorted(set([threshold_AUC, 0.9])):
                ax.axvline(tval, c="black", lw=.5, ls="--")
                ax.axhline(tval, c="black", lw=.5, ls="--")

            if Ir==0 and Ic==(len(sorted_cols)-1):
                ax.legend(bbox_to_anchor=(1.35,1), loc="upper left")

            else:
                ax.get_legend().remove()

            if Ic==0:
                ax.set_ylabel("test $AUC_{ROC}$")

            else:
                ax.set_ylabel("")
                ax.set_yticklabels([])

            if Ic==(len(sorted_cols)-1):
                ax2 = ax.twinx()
                ax2.set_ylabel("%s\n%s split"%(type_phenotype, CVset_type_split))
                ax2.set_yticks([])
                ax2.set_yticklabels([])

            if Ir==(len(sorted_rows)-1) and Ic==1:
                ax.set_xlabel("average train $AUC_{ROC}$")

            elif Ir==(len(sorted_rows)-1):
                ax.set_xlabel("")

            else:
                ax.set_xlabel("")
                ax.set_xticklabels([])

            if Ir==0:
                ax.set_title("consistency by:\n" + consistency_measure, fontsize=10)



    plt.subplots_adjust(wspace=0.05, hspace=0.05)# wspace (between cols) hspaces (btw rows)
    plt.show()
    save_fig(fig, "ML_models_results_color_models.pdf", PlotsDir)



def get_filtered_df_by_dict(df, parms_dict, check_all_present=False):

    """filter df"""

    # check
    if len(df)!=len(set(df.index)):
        raise ValueError("index should be unique")

    # init index
    idx = pd.Series(True, index=df.index)
    for k,v in parms_dict.items():

        if type(v)==set: 
            idx = (idx) & (df[k].isin(v))

        else:
            idx = (idx) & (df[k]==v)

    df_filt = df[idx].copy()

    # check all present
    if check_all_present is True:
        for k,v in parms_dict.items():

            if type(v)==set:
                expected_vals = v

            elif type(v) in {int, str, float}:
                expected_vals = {v}

            else: raise ValueError("not val types")

            if set(df_filt[k])!=expected_vals: raise ValueError("should be filts")

    return df_filt


def visualize_df_ML_parm_stats(df_ML_parm_stats, PlotsDir, max_trainTest_diff):

    """Visualize df"""

    df_ML_parm_stats = df_ML_parm_stats.copy().sort_values(by=["n_phenotypes_good_model"], ascending=[True])

    #fig = plt.figure(figsize=(2.5, 2.5))
    fig = plt.figure(figsize=(2, 2.8))

    #all_n_phenotypes_overfit = sorted(set(df_ML_parm_stats["n_phenotypes_overfit"]))
    all_n_phenotypes_good_model = sorted(set(df_ML_parm_stats["n_phenotypes_good_model"]))

    print(sum(df_ML_parm_stats.valid_parm), "valid parameters")

    #df_ML_parm_stats["n_phenotypes_overfit"] = df_ML_parm_stats.n_phenotypes_overfit + np.random.uniform(-0.1, 0.1, len(df_ML_parm_stats))
    df_ML_parm_stats["n_phenotypes_good_model_numeric"] = df_ML_parm_stats.n_phenotypes_good_model + np.random.uniform(-0.1, 0.1, len(df_ML_parm_stats))


    #ax = sns.scatterplot(data=df_ML_parm_stats, x="mean_train_minus_test", y="n_phenotypes_good_model", style=None, linewidth=.7, edgecolor="black", color="none")
    #ax = sns.scatterplot(data=df_ML_parm_stats, x="mean_train_minus_test", y="std_train_minus_test", style="n_phenotypes_good_model", linewidth=.7, edgecolor="black", color="none", hue="n_phenotypes_good_model", palette={0:"gray", 1:"red", 2:"cyan"}, hue_order=[0, 1, 2])
    #ax = sns.scatterplot(data=df_ML_parm_stats, x="mean_train_minus_test_highAUCtrain", y="n_phenotypes_good_model_numeric", style="n_phenotypes_good_model", linewidth=.7, edgecolor="black", color="none", hue="n_phenotypes_good_model", palette={0:"gray", 1:"red", 2:"cyan"}, hue_order=[0, 1, 2])
    #ax = sns.scatterplot(data=df_ML_parm_stats, x="mean_train_minus_test_highAUCtrain", y="n_phenotypes_good_model_numeric", linewidth=.7, edgecolor="black", color="none")
    #ax = sns.scatterplot(data=df_ML_parm_stats, x="mean_train_minus_test_highAUCtrain", y="n_phenotypes_good_model_numeric", linewidth=.7, edgecolor="black", color="none", hue="valid_parm", palette={True:"red", False:"gray"})
    #ax = sns.swarmplot(data=df_ML_parm_stats, x="n_phenotypes_good_model", y="mean_train_minus_test_highAUCtrain", linewidth=.7, edgecolor="black", hue="valid_parm", palette={True:"red", False:"gray"})
    ax = sns.stripplot(data=df_ML_parm_stats, x="n_phenotypes_good_model", y="mean_train_minus_test_highAUCtrain", linewidth=.7, edgecolor="black", hue="valid_parm", palette={True:"red", False:"gray"})

    # hue="valid_parm", palette={True:"red", False:"gray"}

    ax.legend(bbox_to_anchor=(1.01,1), loc="upper left", title="valid parameter")

    plt.axhline(max_trainTest_diff, c="gray", ls="--", lw=.7)

    ax.set_title("Assessing overfitting", fontsize=10)

    #for x in all_n_phenotypes_overfit:
    #    plt.axvline(x, c="gray", ls="--", lw=.7)
 
    #for y in all_n_phenotypes_good_model:
    #    plt.axhline(y, c="gray", ls="--", lw=.7)

    save_fig(fig, "ML_parameters_results.pdf", PlotsDir)


def visualize_df_ML_integrated_performance_filt_each_pheno(df_ML_integrated_performance, PlotsDir, threshold_AUC, phenotype_to_color, filters_dict, row_field=None, col_field=None, hue_field=None, symbol_field=None, add_best_test_parmID=False):

    """Visualize the performances"""

    df_plot_all = get_filtered_df_by_dict(df_ML_integrated_performance, filters_dict)

    ######## DEFINE VISUALS #######

    # keep
    df_plot_all["real_phenotype"] = df_plot_all.phenotype.copy()

    # define rows and cols
    if row_field is None:
        row_field = "row_f"
        df_plot_all[row_field] = ""
        sorted_rows = [""]

    else:
        sorted_rows = sorted(set(df_plot_all[row_field]))

    if col_field is None:
        col_field = "col_f"
        df_plot_all[col_field] = ""
        sorted_cols = [""]

    else:
        sorted_cols = sorted(set(df_plot_all[col_field]))

    # define visuals
    if hue_field=="phenotype":

        palette_dict = cp.deepcopy(phenotype_to_color)
        if "gray" in palette_dict: raise ValueError("can't be empty ones")
        palette_dict["other phenotypes"] = "gray"

        def get_modified_pheno(p):
            if p in phenotype_to_color: return p
            else: return "other phenotypes"

        df_plot_all["phenotype"] = df_plot_all.phenotype.apply(get_modified_pheno)

    else:
        raise ValueError("invalid hue_f")

    if symbol_field=="aligner":
        symbol_dict = {"bwa_mem":"o", "bowtie2_local":"X", "hisat2_no_spliced":"s"}

    elif symbol_field=="multiple_test_corr":
        symbol_dict = {"none":"o", "maxT_min":"X", "maxT_mean":"s"}

    elif symbol_field=="CVset_type_split":
        symbol_dict = {"balanced":"o", "random":"X"}

    elif symbol_field=="consistency_measure":
        symbol_dict = {"none":"o", "jaccard_variants":"X", "equal_variants":"s", "jaccard_samples":"P", "equal_samples":"^"}

    elif symbol_field is None:
        symbol_dict = None


    else:
        raise ValueError("invalid hue_f")

    ###############################

    ##### PLOT ####

    df_plot_all = df_plot_all.set_index([row_field, col_field])

    nrows = len(sorted_rows)
    ncols = len(sorted_cols)
    fig = plt.figure(figsize=(ncols*2.1, nrows*2.1)); I=1

    for Ir, row_val in enumerate(sorted_rows):
        for Ic, col_val in enumerate(sorted_cols):

            ax = plt.subplot(nrows, ncols, I); I+=1
            df_plot = df_plot_all.loc[(row_val, col_val)].copy()

            targetfields = ["real_phenotype", hue_field, symbol_field]
            if row_field==symbol_field or col_field==symbol_field: targetfields = ["real_phenotype", hue_field]
            if symbol_field is None:  targetfields = ["real_phenotype"]

            if len(df_plot)!=len(df_plot[targetfields].drop_duplicates()): 
                print(df_plot)
                raise ValueError("invalid")

            ax = sns.scatterplot(data=df_plot, x="mean_AUC_train", y="AUC_test", hue=hue_field, style=symbol_field, linewidth=.7, edgecolor="black", palette=palette_dict, markers=symbol_dict)

            # add best parmID
            if add_best_test_parmID is True:

                df_plot_best = df_plot[df_plot.best_test_parmID==True].copy()
                for Ibest, rbest in df_plot_best.iterrows():

                    color_best = palette_dict[rbest[hue_field]]

                    x = rbest.mean_AUC_train
                    y = rbest.AUC_test
                    plt.arrow(x+0.07, y+0.07, -0.02, -0.02, width=0.01, shape="full", color=color_best)


            lims = [0.45, 1.05]
            ax.set_xlim(lims)
            ax.set_ylim(lims)

            plt.plot(lims, lims, c="gray", lw=.5, ls="--")
            for tval in sorted(set([threshold_AUC])):
                ax.axvline(tval, c="black", lw=.5, ls="--")
                ax.axhline(tval, c="black", lw=.5, ls="--")

            #ax.axhline(threshold_ROC_AUC_overfit, c="black", lw=.5, ls="--")

            if Ir==0 and Ic==(len(sorted_cols)-1):
                ax.legend(bbox_to_anchor=(1.35,1), loc="upper left")

            else:
                ax.get_legend().remove()

            if Ic==0:
                ax.set_ylabel("test $AUC_{ROC}$")

            else:
                ax.set_ylabel("")
                ax.set_yticklabels([])

            if Ic==(len(sorted_cols)-1):
                ax2 = ax.twinx()

                if " & " in row_field: ylabel = "\n".join(["%s=%s"%(rf[0:1], row_val.split(" & ")[I]) for I, rf in enumerate(row_field.split(" & "))])
                elif row_field=="row_f": ylabel = ""
                else: ylabel = "%s:\n%s"%(row_field, row_val)

                ax2.set_ylabel(ylabel)
                ax2.set_yticks([])
                ax2.set_yticklabels([])

            #if Ir==(len(sorted_rows)-1) and Ic==int(np.median(range(len(sorted_cols)))):
            if Ir==(len(sorted_rows)-1):
                ax.set_xlabel("average train $AUC_{ROC}$")

            elif Ir==(len(sorted_rows)-1):
                ax.set_xlabel("")

            else:
                ax.set_xlabel("")
                ax.set_xticklabels([])

            if Ir==0:

                if Ic==int(np.median(range(len(sorted_cols)))): title_prefix = "%s\n\n"%filters_dict
                else: title_prefix = ""

                if " & " in col_field: title_str = "\n".join(["%s=%s"%(cf[0:1], col_val.split(" & ")[I]) for I, cf in enumerate(col_field.split(" & "))])
                else: title_str = "%s:\n%s"%(col_field, col_val)

                ax.set_title("%s%s"%(title_prefix, title_str), fontsize=10)



    plt.subplots_adjust(wspace=0.05, hspace=0.05)# wspace (between cols) hspaces (btw rows)
    plt.show()


    def format_val(x):
        if type(x)==set: return "_".join(sorted(map(str, x)))
        else: return x

    filters_dict_str = "_".join(["%s=%s"%(f[0:2], str(format_val(filters_dict[f]))[0:2]) for f in sorted(filters_dict)])
    filename_plot = "ML_models_integrated_per_pheno/row=%s_col=%s_hue=%s_symb=%s_%s.pdf"%(row_field, col_field, hue_field, symbol_field, filters_dict_str)
    save_fig(fig, filename_plot.replace(" & ", "_"), PlotsDir)

    ###############


def get_performance_stats_one_parm_across_phenotypes(df_parmID, threshold_ROC_AUC, fields_optimize_test):

    """Gets the info of one parm ID across all phenos"""

    if len(df_parmID)!=13: raise ValueError("should be all phenos")

    # init
    fields_dict = dict(df_parmID[fields_optimize_test].iloc[0])

    # add indices
    for f, sorted_vals in [("aligner", ["bwa_mem", "bowtie2_local", "hisat2_no_spliced"]),
                           ("CVset_type_split", ["balanced", "random"]),
                           ("consistency_measure", ["jaccard_variants", "jaccard_samples", "equal_variants", "equal_samples", "none"]),
                           ("AUC_diff_threshold", [0.05, 0.1, 0.15, 0.2]),
                           ("sample_set", ["representative_samples", "all_samples"]),
                           ("multiple_test_corr", ["maxT_min", "maxT_mean", "none"])]:

        map_dict = dict(zip(sorted_vals, range(len(sorted_vals))))
        fields_dict["%s_idx"%f] = map_dict[fields_dict[f]]

    # add samples
    #fields_dict["phenotypes_overfit"] = set(df_parmID[(df_parmID.AUC_test<threshold_ROC_AUC_overfit) & (df_parmID.mean_AUC_train>=threshold_ROC_AUC)].phenotype)
    #fields_dict["phenotypes_good_model"] = set(df_parmID[(df_parmID.AUC_test>=threshold_ROC_AUC) & (df_parmID.mean_AUC_train>=threshold_ROC_AUC)].phenotype)
    #fields_dict["phenotypes_good_model"] = set(df_parmID[(df_parmID.AUC_test>=threshold_ROC_AUC) & (df_parmID.mean_AUC_train>=threshold_ROC_AUC) & (df_parmID.all_models_found_idx==0)].phenotype)

    # define the phenotypes that have a good moel by this parm ID
    fields_dict["phenotypes_good_model"] = set(df_parmID[(df_parmID.AUC_test>=threshold_ROC_AUC) & (df_parmID.mean_AUC_train>=threshold_ROC_AUC)].phenotype)
    fields_dict["n_phenotypes_good_model"] = len(fields_dict["phenotypes_good_model"])

    # define a df with the parms to evaluate
    df_parmID_filt = df_parmID[df_parmID.mean_AUC_train>=(threshold_ROC_AUC)].copy()

    # define the mean abslute distance between each 
    for calc_df, suff in [(df_parmID, ""), (df_parmID_filt, "_highAUCtrain")]:
        all_diffs = calc_df.mean_AUC_train - calc_df.AUC_test

        for fun_name, fun_obj in [("mean", np.mean), ("max", max), ("std", np.std), ("median", np.median)]:
            if len(calc_df)==0:
                fields_dict["%s_train_minus_test%s"%(fun_name, suff)] = 0.0

            else:
                fields_dict["%s_train_minus_test%s"%(fun_name, suff)] = fun_obj(all_diffs)

    return pd.Series(fields_dict)

def make_flat_listOflists(LoL):

    return list(itertools.chain.from_iterable(LoL))


def get_harmonic_mean(x, y):

    """Gets the harmonic mean between two values"""

    if x==0 or y==0: return 0.0
    else: return (2*x*y)/(x+y)


def get_r_df_metrics_with_added_vals(r, labels, probabilities):

    """Adds metrics to vals"""

    # get predicted labels
    bool_dict = {True:1, False:0}
    predicted_labels = np.array(list(map(lambda b: bool_dict[b], probabilities>=r.threshold)))
    labels = cp.deepcopy(np.array(labels))

    # get vals
    r["TP"] = int(sum((labels==1) & (predicted_labels==1)))
    r["FP"] = int(sum((labels==0) & (predicted_labels==1)))
    r["FN"] = int(sum((labels==1) & (predicted_labels==0)))
    r["TN"] = int(sum((labels==0) & (predicted_labels==0)))

    # functions 
    def calc_recall():
        if r.TP==0 and r.FN==0: return 1.0
        else: return round(r.TP / (r.TP + r.FN), 4)

    def calc_precision():
        if r.TP==0 and r.FP==0: return 1.0
        else: return round(r.TP / (r.TP + r.FP), 4)

    def calc_specificity():
        if r.TN==0 and r.FP==0: return 0.0
        else: return round(r.TN / (r.TN + r.FP), 4)

    def calc_neg_pred_val():
        if r.TN==0 and r.FN==0: return 0.0
        else: return round(r.TN / (r.TN + r.FN), 4)

    # checks
    if round(r.recall, 4) != calc_recall():
        raise ValueError("recall not properly calculated. %s vs %s"%(r.recall, calc_recall()))

    if round(r.precision, 4) != calc_precision():
        raise ValueError("precision not properly calculated. %s vs %s"%(r.precision, calc_precision()))

    if round((1-r.fpr), 4) != calc_specificity():
        raise ValueError("spec. not properly calculated. %s vs %s"%(1-r.fpr, calc_specificity()))

    # add negative predictive value
    r["neg_predictive_val"] = calc_neg_pred_val()

    # change some stats
    if r.TP==0 and r.FN==0:
        r["recall"] = 0.0

    if r.TP==0 and r.FP==0:
        r["precision"] = 0.0

    return r


def get_df_metrics_labels_and_probabilities(labels, probabilities):

    """get metrics df"""

    # get metrics
    precision, recall, threshold_prc = sk_metrics.precision_recall_curve(labels, probabilities, drop_intermediate=False)
    fpr, tpr, threshold_roc = sk_metrics.roc_curve(labels, probabilities, drop_intermediate=False)
    threshold_prc = np.append(threshold_prc, [np.inf])

    # get as df
    df_PRC = pd.DataFrame({"precision":precision, "recall":recall, "threshold":threshold_prc})
    df_ROC = pd.DataFrame({"fpr":fpr, "tpr":tpr, "threshold":threshold_roc})
    df_merged = df_PRC.merge(df_ROC, validate="one_to_one", how="left", on="threshold")
    check_no_nans_in_df(df_merged)

    # add the calculated metrics (TP, FP, TN, FN)
    df_merged = df_merged.apply(get_r_df_metrics_with_added_vals, labels=labels, probabilities=probabilities, axis=1)

    # add the harmonic means
    df_merged["F_score_precision_recall"] = df_merged.apply(lambda r: get_harmonic_mean(r.precision, r.recall), axis=1)
    df_merged["har_mean_neg_pred_val_recall"] = df_merged.apply(lambda r: get_harmonic_mean(r.neg_predictive_val, r.recall), axis=1)

    # checks
    if not all(df_merged.tpr == df_merged.recall):
        raise ValueError("TPR and recall should be the same")

    return df_merged

def get_df_ML_metrics_best_integrated(df_ML_results_best):

    """Gets the metrics of different thresholds"""

    # keep
    df_ML_results_best = df_ML_results_best.copy()

    # prepare data to plot
    df_metrics_all = pd.DataFrame()
    for phenotype, df_pheno in df_ML_results_best.groupby("phenotype"):

        # get values
        labels = make_flat_listOflists(df_pheno.test_phenos.values)
        probabilities = make_flat_listOflists(df_pheno.test_probabilities.values)

        # get the metrics
        df_merged = get_df_metrics_labels_and_probabilities(labels, probabilities)

        # add the expected genes
        tuple_expected_genes = tuple(sorted(set(df_pheno.expected_genes_str).difference({""})))
        if not tuple_expected_genes in {tuple([]), ('ERG11',), ('Scer_HAL9',), ('ERG11', 'Scer_HAL9')}:
            raise ValueError("tuple_expected_genes: %s"%str(tuple_expected_genes))

        df_merged["tuple_expected_genes"] = [tuple_expected_genes]*len(df_merged)

        # define best threshold by F1 score
        best_threshold = df_merged.sort_values(by=["F_score_precision_recall", "threshold"], ascending=[False, True]).threshold.iloc[0]

        # define best threshold by optimizing recall and neg_predictive_val
        #best_threshold = df_merged.sort_values(by=["har_mean_neg_pred_val_recall", "threshold"], ascending=[False, True]).threshold.iloc[0]

        # add to df
        df_merged["best_threshold"] = df_merged.threshold==best_threshold

        # keep
        df_merged["phenotype"] = phenotype
        df_merged["AUC_test"] = df_pheno.ROC_AUC_test_integrated.iloc[0]
        df_metrics_all = pd.concat([df_metrics_all, df_merged]).reset_index(drop=True)

    return df_metrics_all

def get_df_ML_metrics_best_with_best_threshold(df_ML_metrics_best_input):

    """Get the best threshold"""

    df_metrics_all = pd.DataFrame()
    for (test_split_ID, phenotype), df_metrics in df_ML_metrics_best_input.copy().groupby(["test_split_ID", "phenotype"]):

        if len(df_metrics)!=len(df_metrics[["threshold", "eval_data"]].drop_duplicates()):
            raise ValueError("invalid")

        # define best threshold by F1 score on train
        df_metrics_train = df_metrics[df_metrics.eval_data=="train"].copy()
        best_threshold = df_metrics_train.sort_values(by=["F_score_precision_recall", "threshold"], ascending=[False, True]).threshold.iloc[0]

        df_metrics = df_metrics.copy()
        df_metrics["best_threshold"] = df_metrics.threshold==best_threshold

        df_metrics_all = pd.concat([df_metrics_all, df_metrics])

    return df_metrics_all

def get_ML_sample_to_pheno(df_ML_results):

    """Maps each sample to the pheno"""

    all_samples = make_flat_listOflists(df_ML_results.test_samples.values)
    all_phenos = make_flat_listOflists(df_ML_results.test_phenos.values)

    # get df
    df_phenos = pd.DataFrame({"sampleID":all_samples, "phenotype_val":all_phenos}).drop_duplicates()
    if len(df_phenos)!=len(set(df_phenos.sampleID)):
        raise ValueError("sID should be unique")

    return df_phenos

def get_predicted_phenotypes_one_threshold(prob_threshold, probabilities, samples, type_phenotype, sampleID_to_pheno):

    """Gets the predicted phenotypes for the sorted samples for one threshold"""

    # checks
    if len(samples)!=len(set(samples)):
        raise ValueError("samples should be unique")

    # get df
    sampleID_field = {"phenotype":"sampleID", "phenotype_transition":"comparisonID"}[type_phenotype]
    df_classification = pd.DataFrame({"probability":probabilities, sampleID_field:samples})
    df_classification["prediction"] = (df_classification.probability>=prob_threshold).map({True:1, False:0})
    check_no_nans_series(df_classification["prediction"])

    # each type_phenotype requires a different case
    if type_phenotype=="phenotype":
        return df_classification[["sampleID", "prediction"]].copy().rename({'prediction':'predicted_phenotype'})

    elif type_phenotype=="phenotype_transition":

        # add fields
        df_classification["target_sampleID"] = df_classification.comparisonID.apply(lambda x: x.split("_vs_")[0])
        df_classification["origin_sampleID"] = df_classification.comparisonID.apply(lambda x: x.split("_vs_")[1])
        df_classification["origin_phenotype"] = df_classification.origin_sampleID.map(sampleID_to_pheno)
        check_no_nans_series(df_classification["origin_phenotype"])

        # create a dict which, for each sampleID records the predictions for both transition and no transitions
        def get_predicted_phenotypes_proportions(df):
            prediction_val = df.prediction.iloc[0]
            o_pheno_to_t_pheno = {0: {0:0, 1:1},
                                  1: {0:1, 1:0}}[prediction_val]

            pred_pheno = df.origin_phenotype.map(o_pheno_to_t_pheno)
            return ",".join(["%.2f"%(sum(pred_pheno==p_pheno_val) / len(pred_pheno)) for p_pheno_val in [0, 1]])

        return df_classification.groupby(["prediction", "target_sampleID"]).apply(get_predicted_phenotypes_proportions).reset_index().rename(columns={0:'predicted_phenotype_proportions', 'target_sampleID':'sampleID', 'prediction':'phenotype_transition'})

def get_df_ML_metrics_best_train_and_test(df_ML_results_best, df_ML_phenos, ProcessedDataDir, replace=False):

    """Gets the metrics of different thresholds"""

    df_metrics_all_file = "%s/df_ML_classifier_best_models_all_metrics.py"%ProcessedDataDir
    df_predictions_all_file = "%s/df_ML_classifier_best_models_all_predictions.py"%ProcessedDataDir

    if file_is_empty(df_metrics_all_file) or file_is_empty(df_predictions_all_file) or replace is True:

        # keep
        df_ML_results_best = df_ML_results_best.copy()

        # prepare datasets
        df_metrics_all = pd.DataFrame()
        df_predictions_all = pd.DataFrame()

        for I,r in df_ML_results_best.iterrows():

            # different datasets
            for eval_data in ["train", "test", "all"]:

                # get the labels and probabilities
                phenos_fields = {"train":["train_phenos"], "test":["test_phenos"], "all":["train_phenos", "test_phenos"]}[eval_data]
                probs_fields = {"train":["train_probabilities"], "test":["test_probabilities"], "all":["train_probabilities", "test_probabilities"]}[eval_data]
                samples_fields = {"train":["train_samples"], "test":["test_samples"], "all":["train_samples", "test_samples"]}[eval_data]

                labels = make_flat_listOflists([r[f] for f in phenos_fields])
                probabilities = make_flat_listOflists([r[f] for f in probs_fields])
                samples = make_flat_listOflists([r[f] for f in samples_fields])

                # get metrics
                df_metrics = get_df_metrics_labels_and_probabilities(labels, probabilities)

                # keep
                sampleID_to_pheno = dict(df_ML_phenos[(df_ML_phenos.phenotype==r.phenotype)].set_index("sampleID").phenotype_val)

                # get the predictions
                for threshold in sorted(set(df_metrics.threshold)):
                    df_predictions = get_predicted_phenotypes_one_threshold(threshold, probabilities, samples, r.type_phenotype, sampleID_to_pheno)

                    df_predictions["eval_data"] = eval_data
                    df_predictions["threshold"] = threshold

                    for f in ["test_split_ID", "expected_genes_str", "phenotype"]:
                        df_predictions[f] = r[f]

                    df_predictions_all = pd.concat([df_predictions_all, df_predictions]).reset_index(drop=True)

                # add fields
                df_metrics["eval_data"] = eval_data
                for f in ["test_split_ID", "expected_genes_str", "phenotype"]:
                    df_metrics[f] = r[f]

                # add
                df_metrics_all = pd.concat([df_metrics_all, df_metrics]).reset_index(drop=True)

        # save
        save_object(df_metrics_all, df_metrics_all_file)
        save_object(df_predictions_all, df_predictions_all_file)

    return load_object(df_metrics_all_file), load_object(df_predictions_all_file)

def get_performance_curves_ML_classifiers(df_ML_metrics_best, PlotsDir, phenotype_to_color, significant_ML_phenotypes):

    """Plots """

    # keep
    df_ML_metrics_best = df_ML_metrics_best.copy()
    df_ML_metrics_best["phenotype_idx"] = df_ML_metrics_best.phenotype.map(dict(zip(significant_ML_phenotypes, range(len(significant_ML_phenotypes)))))
    check_no_nans_series(df_ML_metrics_best["phenotype_idx"])
    df_ML_metrics_best["specificity (1-fpr)"] = 1 - df_ML_metrics_best.fpr

    df_ML_metrics_best = df_ML_metrics_best.sort_values(by=["phenotype_idx", "threshold"])

    # def visuals
    phenotype_to_color = cp.deepcopy(phenotype_to_color)
    missing_ps = set(significant_ML_phenotypes).difference(set(phenotype_to_color))
    if missing_ps!=set():
        raise ValueError("missing phenos: %s"%missing_ps)

    # plot
    fig = plt.figure(figsize=(4, 8))

    #df_ML_metrics_best =df_ML_metrics_best[df_ML_metrics_best.phenotype.isin({"agar_invasiveness", "FLC_resistance"})]

    #for Ip, (xfield, yfield, title) in enumerate([("tpr", "fpr", ""), ("recall", "precision", "")]):
    for Ip, (xfield, yfield, title) in enumerate([("tpr", "specificity (1-fpr)", ""), ("recall", "precision", "")]):

        ax = plt.subplot(2, 1, Ip+1)

        ax = sns.lineplot(data=df_ML_metrics_best, x=xfield, y=yfield, hue="phenotype", lw=2, palette=phenotype_to_color, estimator=None)
        ax.legend(bbox_to_anchor=(1.01,1), loc="upper left")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("tpr (recall)")

        lims = [-0.05, 1.05]
        ax.set_xlim(lims)
        ax.set_ylim(lims)


    plt.subplots_adjust(wspace=0.05, hspace=0.25)# wspace (between cols) hspaces (btw rows)
    plt.show()
    save_fig(fig, "ML_model_performance_curves.pdf", PlotsDir)



def plot_PRC_ROC_curves_metrics(df_ML_metrics_best, PlotsDir, significant_ML_phenotypes, df_ML_results_best):

    """PRC curves all"""

    # keep
    df_ML_metrics_best = df_ML_metrics_best.copy().rename({"tpr":"recall"})
    df_ML_metrics_best["phenotype_idx"] = df_ML_metrics_best.phenotype.map(dict(zip(significant_ML_phenotypes, range(len(significant_ML_phenotypes)))))
    check_no_nans_series(df_ML_metrics_best["phenotype_idx"])
    df_ML_metrics_best["specificity (1-fpr)"] = 1 - df_ML_metrics_best.fpr

    df_ML_metrics_best = df_ML_metrics_best.sort_values(by=["phenotype_idx", "test_split_ID", "threshold"])

    # plot
    nrows = len(significant_ML_phenotypes)
    ncols = 5
    fig = plt.figure(figsize=(ncols*2.5, nrows*2.5)); Ip=1

    for Ir, phenotype in enumerate(significant_ML_phenotypes):
        print(phenotype)

        for Ic, test_split_ID in enumerate(range(5)):
            df_plot = df_ML_metrics_best[(df_ML_metrics_best.phenotype==phenotype) & (df_ML_metrics_best.test_split_ID==test_split_ID)].copy()

            # skip one
            ML_results = df_ML_results_best[(df_ML_results_best.phenotype==phenotype) & (df_ML_results_best.test_split_ID==test_split_ID)].iloc[0].copy()
            n_predictors = len(ML_results.selected_predictors)

            # get subplot
            ax = plt.subplot(nrows, ncols, Ip); Ip+=1

            if n_predictors==0: 
                ax.text(0.1, 0.5, "No model")

            else:

                # get long plot stacked
                df_plot_long = pd.DataFrame()
                palette_dict = {"precision":"red", "recall":"c", "neg_predictive_val":"blue", "specificity (1-fpr)":"black"}

                for f in ["precision", "recall", "specificity (1-fpr)"]: #,  neg_predictive_val
                    df_plot_long_f = df_plot[["threshold", f, "eval_data"]].copy().rename(columns={f:"value"})
                    df_plot_long_f["metric"] = f
                    df_plot_long = pd.concat([df_plot_long, df_plot_long_f]).reset_index(drop=True)

                if len(df_plot_long)!=len(df_plot_long[["metric", "eval_data", "threshold"]].drop_duplicates()):
                    raise ValueError("shold be unique")

                ax = sns.lineplot(data=df_plot_long.sort_values(by=["eval_data", "metric", "threshold"]), x="threshold", y="value", style="eval_data", hue="metric", palette=palette_dict, lw=2, estimator=None, markers={"all":"o", "train":"s", "test":"P"}, markeredgecolor="black", markeredgewidth=.8)

                if Ic==4 and Ir==0:
                    ax.legend(bbox_to_anchor=(1.2,1), loc="upper left")
                
                else:
                    ax.get_legend().remove()

                str_aucs = "$AUC_{tr}$=%s, $AUC_{te}$=%s"%(("%.2f"%ML_results.ROC_AUC_train), ("%.2f"%ML_results.ROC_AUC_test))
                ax.set_title("m%i (%s) %i features\n%s"%(test_split_ID, ML_results.model_name.split("|")[0], n_predictors, str_aucs), fontsize=10)

                # add threshold
                xval_threshold = df_plot[df_plot.best_threshold==True].threshold.iloc[0]
                plt.axvline(xval_threshold, c="k", ls=":", lw=.7)

            # general things
            lims = [-0.05, 1.05]
            ax.set_xlim(lims)
            ax.set_ylim(lims)

            if Ic==4:
                ax2 = ax.twinx()
                ax2.set_ylabel(phenotype)
                ax2.set_yticks([])
                ax2.set_yticklabels([])


            if Ir==(len(significant_ML_phenotypes)-1) and Ic==2:
                ax.set_xlabel("threshold")

            else:
                ax.set_xlabel("")

    plt.subplots_adjust(wspace=0.5, hspace=0.5)# wspace (between cols) hspaces (btw rows)
    plt.show()
    return
    save_fig(fig, "ML_model_performance_curves_each_model.pdf", PlotsDir)

def plot_PRC_ROC_curves_metrics_integrated(df_ML_metrics_best, PlotsDir, significant_ML_phenotypes, phenotype_to_color):

    """Plots """

    # keep
    df_ML_metrics_best = df_ML_metrics_best.copy().rename({"tpr":"recall"})
    df_ML_metrics_best["phenotype_idx"] = df_ML_metrics_best.phenotype.map(dict(zip(significant_ML_phenotypes, range(len(significant_ML_phenotypes)))))
    check_no_nans_series(df_ML_metrics_best["phenotype_idx"])
    df_ML_metrics_best["specificity (1-fpr)"] = 1 - df_ML_metrics_best.fpr

    df_ML_metrics_best = df_ML_metrics_best.sort_values(by=["phenotype_idx", "threshold"])

    # plot
    nrows = len(significant_ML_phenotypes)
    ncols = 3
    fig = plt.figure(figsize=(ncols*2.5, nrows*2.5)); Ip=1

    for Ir, phenotype in enumerate(significant_ML_phenotypes):
        df_plot = df_ML_metrics_best[df_ML_metrics_best.phenotype==phenotype].copy()

        for Ic, type_plot in enumerate(["ROC", "PRC", "threshold_vs_vals"]):

            # get subplot
            ax = plt.subplot(nrows, ncols, Ip); Ip+=1

            # plot
            if type_plot=="threshold_vs_vals":
                df_plot_long = pd.DataFrame()
                for f in ["precision", "recall", "specificity (1-fpr)", "neg_predictive_val"]:
                    df_plot_long_f = df_plot[["threshold", f]].copy().rename(columns={f:"value"})
                    df_plot_long_f["metric"] = f
                    df_plot_long = pd.concat([df_plot_long, df_plot_long_f]).reset_index(drop=True)

                ax = sns.lineplot(data=df_plot_long, x="threshold", y="value", style="metric", c=phenotype_to_color[phenotype], lw=2, estimator=None, markers="oPs^", markeredgecolor="black", markeredgewidth=.8)
                ax.legend(bbox_to_anchor=(1.2,1), loc="upper left")

                ax2 = ax.twinx()
                ax2.set_ylabel(phenotype)
                ax2.set_yticks([])
                ax2.set_yticklabels([])

                xval_threshold = df_plot[df_plot.best_threshold==True].threshold.iloc[0]

            elif type_plot=="ROC":
                ax = sns.lineplot(data=df_plot, x="fpr", y="recall", c=phenotype_to_color[phenotype], lw=2, estimator=None)
                xval_threshold = df_plot[df_plot.best_threshold==True].fpr.iloc[0]

            elif type_plot=="PRC":
                ax = sns.lineplot(data=df_plot, x="recall", y="precision", c=phenotype_to_color[phenotype], lw=2, estimator=None)
                xval_threshold = df_plot[df_plot.best_threshold==True].recall.iloc[0]


            lims = [-0.05, 1.05]
            ax.set_xlim(lims)
            ax.set_ylim(lims)
            plt.axvline(xval_threshold, c="k", ls=":", lw=.7)


    plt.subplots_adjust(wspace=0.5, hspace=0.4)# wspace (between cols) hspaces (btw rows)
    plt.show()
    save_fig(fig, "ML_model_performance_curves.pdf", PlotsDir)


def get_best_parameters_from_ML_models_one_phenotype(df_ML_results):

    """Gets the best ASR_method and min support and sample_set"""

    # checks
    df_ML_results = df_ML_results.copy()
    if len(df_ML_results)!=len(set(df_ML_results.test_split_ID)):
        raise ValueError("there should be 5 rows")

    if len(set(df_ML_results.aligner))!=1: raise ValueError("should be 1 aln")
    if len(set(df_ML_results.phenotype))!=1: raise ValueError("should be 1 phenotype")
    if len(set(df_ML_results.ROC_AUC_test_integrated))!=1: raise ValueError("should be 1 ROC_AUC_test_integrated")

    #  init with adapt
    best_model_params = {"aligner" : df_ML_results.iloc[0].aligner}

    # get a parameter data
    def get_stats_combination_parms(df, group_fields):
        fields_dict = dict(df[group_fields].iloc[0])
        fields_dict["n_splits"] = len(df)
        fields_dict["mean_min_ROC_AUC_CVset_resamples"] = np.mean(df.min_ROC_AUC_CVset_resamples)
        fields_dict["mean_CVset_model_consistency"] = np.mean(df.CVset_model_consistency)
        fields_dict["mean_ROC_AUC_test"] = np.mean(df.ROC_AUC_test)
        return pd.Series(fields_dict)

    group_fields = ["min_support", "ASR_method", "sample_set"]
    df_parameter_data = df_ML_results.groupby(group_fields).apply(get_stats_combination_parms, group_fields=group_fields).reset_index(drop=True)

    # keep most consistent one
    df_parameter_data = get_df_like_first_row(df_parameter_data, ["n_splits"], [False])

    # keep the one with high AUC
    AUC_t = max(df_parameter_data.mean_min_ROC_AUC_CVset_resamples) - 0.05
    df_parameter_data = df_parameter_data[df_parameter_data.mean_min_ROC_AUC_CVset_resamples>=AUC_t].copy()

    # sort by other criteria
    df_parameter_data["ASR_method_idx"] = df_parameter_data.ASR_method.map({"DOWNPASS":0, "MPPA,DOWNPASS":1, "MPPA":2})
    df_parameter_data = get_df_like_first_row(df_parameter_data, ["mean_CVset_model_consistency", "ASR_method_idx", "mean_min_ROC_AUC_CVset_resamples"], [False, True, False])

    if len(df_parameter_data)!=1: raise ValueError("invalid")

    for f in ["min_support", "ASR_method", "sample_set"]:
        best_model_params[f] = df_parameter_data[f].iloc[0]

    return best_model_params

def get_df_groups_only_subset(groups_file, df_groups_filt_file, set_target_groups):

    """Gets groups only of set_target_groups"""

    if file_is_empty(df_groups_filt_file):

        df_groups_all = get_tab_as_df_or_empty_df(groups_file)
        df_groups_all["unique_groupID"] = df_groups_all.grouping_ID + "|||" + df_groups_all.group

        strange_groups_input = set_target_groups.difference(set(df_groups_all["unique_groupID"]))
        if len(strange_groups_input)>0: raise ValueError("strange_groups_input: %s"%strange_groups_input)
        
        save_object(df_groups_all[df_groups_all.unique_groupID.isin(set_target_groups)], df_groups_filt_file)

    return load_object(df_groups_filt_file)

def get_simple_var_name_variantID_across_samples(x):

    """gets the simple var name"""

    if x.startswith("cluster"): 
        return x
    
    elif (x.startswith("Contig") or x.startswith("mito_")) and len(x.split("_"))==6:
        x_split = x.split("_")
        return "%s_%s_%s"%(x_split[0], x_split[-2], x_split[-1])

    else:
        raise ValueError("invalid %s"%x)

def get_sample_name_w_transition_samples(sID, extra_info, sorted_samples_phylogenetic):

    """Gets the sample ID"""

    # get name
    def get_sname_from_ID(x): return x.split("_")[1][1:]
    s_name = get_sname_from_ID(sID)

    # get transitions
    transitions_list = sorted({x for x in extra_info.split("|||")}.difference({""}).difference({"|transitions|", "transitions|", "transitions"}))

    if len(transitions_list)==0:
        return s_name

    else:

        # get the sample to the transitions
        sample_to_transitions_all = pd.Series()
        for trans_str in transitions_list:
            sample_to_trans = pd.Series({x.split(":")[0] : set(x.split(":")[1].split(",")) for x in trans_str.split("transitions|")[1].split(";")})
            sample_to_transitions_all = pd.concat([sample_to_transitions_all, sample_to_trans]).copy()

        # checks
        if len(sample_to_transitions_all)!=len(set(sample_to_transitions_all.index)):
            raise ValueError("s should be unique")

        # no trans
        if not sID in sample_to_transitions_all.index:
            return s_name

        # get the transitions
        all_transitions = sample_to_transitions_all.loc[sID]
        sorting_dict = dict(zip(sorted_samples_phylogenetic, range(len(sorted_samples_phylogenetic))))
        all_origin_samples = sorted(set(make_flat_listOflists(list(map(lambda x: x.split("_vs_"), all_transitions)))).difference({sID}), key=(lambda val: sorting_dict[val]))

        return "%s (vs %s)"%(s_name, ", ".join(list(map(get_sname_from_ID, all_origin_samples))))
 

def get_sample_to_predicted_phenotype_from_metrics_df(df_ML_metrics, best_threshold_field):

    """maps each sample to the predicted phenotypes"""

    # checks
    df_ML_metrics = df_ML_metrics.copy()
    if len(df_ML_metrics)!=len(df_ML_metrics[["threshold"]].drop_duplicates()):
        raise ValueError("should be unique")

    # get the best threshold
    df_ML_metrics_best_threshold = df_ML_metrics[(df_ML_metrics[best_threshold_field]==True)].copy()
    if len(df_ML_metrics_best_threshold)!=1:
        raise ValueError("should be 1")

    series_best_threshold = df_ML_metrics_best_threshold.iloc[0].copy()
    series_best_threshold["specificity"] = 1 - series_best_threshold.fpr

    # get the mapping
    pred_pheno_dict = {"1":1, "0":0, "?":-1}
    sampleID_to_pred_pheno = {x.split(":")[0] : pred_pheno_dict[x.split(":")[1].split(",")[1]] for x in series_best_threshold.predicted_phenotypes.split("|")}

    return sampleID_to_pred_pheno, dict(series_best_threshold[['recall', 'precision', 'specificity', 'ROC_AUC']])

def generate_df_predicted_phenotypes_from_metrics_df(predicted_phenotypes_file, df_ML_metrics, eval_data_type, best_threshold_field):


    """gets the pred phenos"""

    if file_is_empty(predicted_phenotypes_file):
        print("Getting pred phenos...")

        # load and check
        if len(df_ML_metrics)!=len(df_ML_metrics[["eval_data", "threshold"]].drop_duplicates()):
            raise ValueError("should be uique")

        # filter
        if eval_data_type=="cross_splits":
            df_ML_metrics = df_ML_metrics[df_ML_metrics.eval_data=="cross_splits"].copy()

        else:
            raise ValueError("invalid")

        print(df_ML_metrics)
        adhadgffdafd

        # define th e best threshold based on train data
        #df_ML_metrics_train = df_ML_metrics[(df_ML_metrics.eval_data=="train") & (df_ML_metrics[opt.best_threshold_field]==True)].copy()
        df_ML_metrics_train = df_ML_metrics[(df_ML_metrics.eval_data=="all") & (df_ML_metrics[opt.best_threshold_field]==True)].copy()
        if len(df_ML_metrics_train)!=1: raise ValueError("should be 1 row")
        best_threshold = df_ML_metrics_train.iloc[0].threshold

        # get the predicted phenotypes row in all
        df_ML_metrics_all = df_ML_metrics[(df_ML_metrics.eval_data=="all") & (df_ML_metrics.threshold==best_threshold)].copy()
        if len(df_ML_metrics_all)!=1: raise ValueError("should be 1 row")
        df_predicted_phenotypes = pd.DataFrame({I : {"sampleID":x.split(":")[0], "phenotype":x.split(":")[1].split(",")[0], "predicted_phenotype":x.split(":")[1].split(",")[1]} for I,x in enumerate(df_ML_metrics_all.iloc[0].predicted_phenotypes.split("|"))}).transpose()

        # change
        df_predicted_phenotypes["phenotype"] = df_predicted_phenotypes["phenotype"].apply(int)
        df_predicted_phenotypes["predicted_phenotype"] = df_predicted_phenotypes["predicted_phenotype"].map({"1":1, "0":0, "?":-1})

        # add the sample_type
        def get_set_samples_from_pred_phenos(pred_phenos):
            return set(map(lambda x: x.split(":")[0], pred_phenos.split("|")))

        train_samples = set.union(*df_ML_metrics[(df_ML_metrics.eval_data=="train")].predicted_phenotypes.apply(get_set_samples_from_pred_phenos))
        test_samples = set.union(*df_ML_metrics[(df_ML_metrics.eval_data=="test")].predicted_phenotypes.apply(get_set_samples_from_pred_phenos))

        if len(train_samples.intersection(test_samples))>0:
            raise ValueError("no intersect allowed")

        sID_to_sample_type = {}
        for s in train_samples: sID_to_sample_type[s] = "train"
        for s in test_samples: sID_to_sample_type[s] = "test"
        df_predicted_phenotypes["type_sample"] = df_predicted_phenotypes.sampleID.map(sID_to_sample_type)

        # checks
        check_no_nans_in_df(df_predicted_phenotypes)

        # get the performance metrics
        df_ML_metrics["specificity"] = 1 - df_ML_metrics.fpr
        df_ML_metrics_best = df_ML_metrics[df_ML_metrics.threshold==best_threshold][['recall', 'precision', 'specificity', "eval_data"]].copy()

        return df_predicted_phenotypes, df_ML_metrics_best


def get_df_groups_only_subset_genes(groups_file, df_groups_filt_file, set_genes):

    """Gets groups only of set_target_groups"""

    if file_is_empty(df_groups_filt_file):

        # get ID
        df_groups = get_tab_as_df_or_empty_df(groups_file)
        df_groups["unique_groupID"] = df_groups.grouping_ID + "|||" + df_groups.group

        # get genes
        df_groups = df_groups[(df_groups.type_collapsing=="genes") & (df_groups["group"].isin(set_genes)) & (df_groups.type_genes=="all_genes") & (df_groups.type_vars=="all_vars") & (df_groups.type_mutations.isin({"non_syn_muts", "all_muts"}))].copy()

        # get unique per group vals
        df_groups_unique = df_groups[["group", "type_mutations", "unique_groupID"]].drop_duplicates()

        # get one group for each gene, prioritizing
        #sorted_type_mut = ["all_muts", "non_syn_muts", "truncating_muts", "non_syn_non_truncating_muts", "syn_muts"]
        sorted_type_mut = ["non_syn_muts", "all_muts"]
        df_groups_unique["type_mutations_idx"] = df_groups_unique.type_mutations.map(dict(zip(sorted_type_mut, range(len(sorted_type_mut)))))
        check_no_nans_series(df_groups_unique["type_mutations_idx"])

        sorting_f = ["group", "type_mutations_idx"]
        if len(df_groups_unique)!=len(df_groups_unique[sorting_f].drop_duplicates()):
            raise ValueError("non unique")

        df_groups_unique = df_groups_unique.sort_values(by=sorting_f).drop_duplicates(subset=["group"], keep="first")

        if sorted(df_groups_unique["group"])!=sorted(set_genes):
            raise ValueError("should be unique fields")

        # filter
        df_groups = df_groups[df_groups.unique_groupID.isin(set(df_groups_unique.unique_groupID))].copy()

        # check
        if set(df_groups["group"])!=set_genes:
            raise ValueError("all genes should be there")

        # save
        save_object(df_groups, df_groups_filt_file)

    return load_object(df_groups_filt_file)


def get_str_cladeID_Bergin2022(x):
    if float(x)==int(x): return "clade_%i"%x
    else: return "clade_%s"%x

def get_gene_to_color():

    """maps genes to color"""

    gene_to_color = {"ERG11":"tab:blue",
                     "MRR1":"tab:cyan",
                     "Scer_HAL9":"teal",
                     "NDT80":"slateblue",
                     "UPC2":"slategray",
                     "ERG3":"darkviolet",

                     "FKS1":"tab:red",
                     "FKS2":"brown",

                     "ERG2":"plum",
                     "Scer_ERG5":"purple",
                     "ERG6":"magenta",
                     "ERG25":"crimson",

                     "Scer_FCY1":"peru",
                     "Scer_FCY2":"gray",
                     "Scer_FUR1":"darkorange",
                     "Scer_CHS7":"gold"}

    return {g : mcolors.to_hex(c) for g,c in gene_to_color.items()}


def get_df_with_all_combinations_values(input_df, target_fields):

    """Generates a df with all the possible combinations of the values of target_fields"""

    # get unique df
    input_df = input_df[target_fields].copy().drop_duplicates()

    # define unique values
    unique_values = {field: input_df[field].unique() for field in target_fields}

    # generate all possible combinations of these unique values
    all_combinations = list(itertools.product(*unique_values.values()))

    # get df
    all_combinations_df = pd.DataFrame(all_combinations, columns=target_fields).sort_values(by=target_fields).reset_index(drop=True)
    return all_combinations_df

def get_df_missing_combinations_values(input_df, target_fields):

    """gets missing combinations"""

    df_all_combinations = get_df_with_all_combinations_values(input_df, target_fields)
    merged_df = df_all_combinations.merge(input_df.copy(), on=target_fields, how='left', indicator=True, validate="one_to_one")
    missing_combinations_df = merged_df[merged_df['_merge']=='left_only'].copy().drop(columns=['_merge'])

    return missing_combinations_df

def print_n_parameters_ML_classifiers(ProcessedDataDir, ML_filters_manually_curated):

    """Prints # hyperparameters of each type"""

    # define fields
    fields_optimize_test = ["aligner", "CVset_type_split", "consistency_measure", "AUC_diff_threshold", "sample_set"] # test_split_ID
    #unique_fields_test = fields_optimize_test + ["multiple_test_corr", "test_sample_set"]
    fields_optimize_CVset = ["min_support", "ASR_method", "features_max_pval_GWAS", "transition_definition", "type_phenotype", "prioritize_expected_genes"]
    fields_model = ["real_type_phenotype", "model_name", "class_weight", "feature_selection_method"]

    # load parm info data
    df_parm_info = load_object("%s/df_ML_parameter_info.py"%ProcessedDataDir)
    print("df loaded")

    # check
    if set(df_parm_info.keys())!=set.union(*[set(fields_optimize_test), set(fields_optimize_CVset), set(fields_model), {"phenotype", "n_models", "test_split_ID"}]):
        raise ValueError("fields not as expected")

    for f in (fields_optimize_test):
        if any(df_parm_info[f]==-1): raise ValueError("no -1 expected")

    if len(df_parm_info)!=len(df_parm_info[["phenotype", "test_split_ID"] + fields_optimize_test + fields_optimize_CVset].drop_duplicates()):
        raise ValueError("non unique vals")

    check_no_nans_in_df(df_parm_info)

    # print the number of test parameters, ignoring multiple_test_corr & test_sample set
    print("# test parm combinations, skipping multiple_test_corr & test_sample: ")
    print(df_parm_info[["phenotype"] + fields_optimize_test].drop_duplicates().groupby("phenotype").apply(len))

    print("all test parm vals")
    for f in fields_optimize_test:
        print(f, sorted(set(df_parm_info[f])))

    print("all CVset parm vals")
    for f in fields_optimize_CVset:
        print(f, sorted(set(df_parm_info[f]).difference({"-1"})))

    # print # test parameters for each setting
    print("# CVset parameters")
    pheno_to_nmodels = df_parm_info.groupby("phenotype").apply(lambda df_p:
        sorted(set(df_p[df_p.n_models!=0].groupby(fields_optimize_test + ["test_split_ID"]).apply(lambda df_t: sum(df_t.n_models))))
        )
    print(pheno_to_nmodels.apply(lambda x: (x[0], x[-1])))
    print("min # parms: %i, max # parms: %i"%(min(pheno_to_nmodels.apply(min)), max(pheno_to_nmodels.apply(max))))


    # print, for each test, the feature_selection_method & class weight used
    for f in ["feature_selection_method", "class_weight"]:
        print("\nall", f)
        for I,x in df_parm_info.groupby("phenotype").apply(lambda df_p: sorted(set(df_p[f]))).items():
            print(I,x)

        print("\n", "all ", f, " ; ".join(sorted(set(df_parm_info[f]))))

    # add the split data
    for f in ["real_type_phenotype", "model_name"]:
        all_vals = sorted(set.union(*df_parm_info[f].apply(lambda x: set(x.split(",")))))
        print("vals", f, "\n", "\n".join(all_vals))

def get_ML_pheno_to_color_old():

    """Maps each pheno to color"""

    # init
    dict_colors = {}
    for pheno in ['MIF_resistance', 'age_old', 'age_young', 'source_EXUDADO_CUTANEO', 'source_ORINA', 'source_blood', 'source_sterile', "source_blood_vs_other_sterile", "source_ASPIRADO_BRONQUIAL"]:
        dict_colors[pheno] = "gray"

    # add
    dict_colors_d = {'FLC_resistance': 'tab:blue', 
                     'VRC_resistance': 'tab:cyan',
                     'agar_invasiveness': 'tab:brown',
                     'age_infant': 'peru'}

    for p,c in dict_colors_d.items():
        dict_colors[p] = c

    # resort
    sorted_phenos = ['source_EXUDADO_CUTANEO', 'source_ORINA', "source_ASPIRADO_BRONQUIAL", 'source_blood', 'source_sterile', "source_blood_vs_other_sterile"] + ["MIF_resistance", "FLC_resistance", "VRC_resistance"] + ["agar_invasiveness"] + ["age_infant", "age_young", "age_old"]

    dict_colors = {p:dict_colors[p] for p in sorted_phenos}

    return dict_colors

def get_ML_pheno_to_color():

    """Maps each pheno to color"""

    # init with no pheno colors
    dict_colors = {}
    for pheno in ['age_old', 'source_EXUDADO_CUTANEO', 'source_sterile', 'age_young', 'source_blood']:
        dict_colors[pheno] = "gray"

    # add
    dict_colors_d = {'FLC_resistance': 'tab:blue', 
                     'VRC_resistance': 'tab:cyan',
                     'MIF_resistance': 'tab:purple',

                     'source_blood_vs_other_sterile': 'tab:red',
                     'source_ORINA': 'tab:olive',

                     'source_ASPIRADO_BRONQUIAL': 'tab:orange',
                     'age_infant': 'tab:pink',
                     'agar_invasiveness': 'tab:brown'}

    for p,c in dict_colors_d.items():
        dict_colors[p] = c

    return dict_colors

def plot_trees_expected_genes_variants_grouped_by_type_drug(DataDir, PlotsDir, df_metadata, samples_discard_run1, strains_tree, df_metadata_binary_phenotypes):

    """Plots variants for different cases"""

    # keep
    strains_tree = cp.deepcopy(strains_tree)
    df_metadata_binary_phenotypes = df_metadata_binary_phenotypes.copy()

    # define dirs
    PlotsDir_trees = "%s/trees_expected_genes_grouped"%PlotsDir; make_folder(PlotsDir_trees)

    # define expected genes
    type_to_genes = {"azole": ["ERG11", "ERG3", "Scer_HAL9", "MRR1", "NDT80", "UPC2"],
                     "echinocandin": ["FKS1", "FKS2", "ERG3"],
                     "AMB": ["ERG2", "ERG3", "Scer_ERG5", "ERG6", "ERG25"],
                     "5FC": ["Scer_FCY1", "Scer_FCY2", "Scer_FUR1", "Scer_CHS7"]}

    mapping_genes = "Scer_HAL9:TAC1,Scer_ERG5:ERG5,Scer_FCY1:FCY1,Scer_FCY2:FCY2,Scer_FUR1:FUR1,Scer_CHS7:CHS7"

    # define drugs
    azoles = ['FLC_resistance', 'VRC_resistance', 'ITR_resistance', 'POS_resistance', 'IVZ_resistance']
    echinocandins = ['CAS_resistance', 'MIF_resistance', 'ANI_resistance']

    pheno_to_short_p = {"FLC_resistance":"Fl", "VRC_resistance":"Vo", "ITR_resistance":"It", "POS_resistance":"Po", "IVZ_resistance":"Is", 
                        "CAS_resistance":"Ca", "MIF_resistance":"Mi", "ANI_resistance":"An", "AMB_resistance":"Am", "5FC_resistance":"Fc"}


    # load gene feats
    gene_features_df = get_gene_features_df(DataDir)
    gname_to_gID = gene_features_df.set_index("final_name").geneID

    # define dirs
    aligner_pheno = "bwa_mem"
    source_blood_dir = "%s/genotype_phenotype_relationships/pheno-source_blood_aligner-%s/convergence_gwas"%(DataDir, aligner_pheno)
    variants_file = "%s/gwas_inputs_all_samples/variants_df_all.tab"%(source_blood_dir)
    variant_annotations_file = "%s/gwas_inputs_all_samples/variant_annotations_df_all.tab"%(source_blood_dir)

    pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
    CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
    genome_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir # genome from CGD
    gff_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir
    interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir

    all_samples = set(get_tab_as_df_or_empty_df(variants_file).sampleID)

    # define samples to discard (manually curated because they are redundant)
    type_to_samples_discard = {"azole":{"L21-253", "L21-175", "L22-443", "L21-62", "L21-074", "L21-038", "L21-170", "L20-077", "L22-224", "L22-222", "L21-279", "L22-439", "L21-249", "L22-122"},
                                "echinocandin":set(),
                                "AMB":set()}

    # keep drugs
    for drug_type, sorted_phenotypes in [("azole", azoles), ("echinocandin",echinocandins), ("AMB", ["AMB_resistance"]), ("5FC", ["5FC_resistance"])]:
        print(drug_type)


        # discard some grup types for which there are no vars
        if drug_type=="5FC": continue
            
        # define dir where to save all
        plots_dir = "%s/%s_drugs"%(PlotsDir_trees, drug_type); make_folder(plots_dir)

        # get all of the samples of this type
        all_samples_type = set()
        for phenotype in sorted_phenotypes:
            all_samples_type.update(set(get_tab_as_df_or_empty_df("%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas/gwas_inputs_representative_samples/variants_df_all.tab"%(DataDir, phenotype, aligner_pheno)).sampleID))

        samples_discard = {[x for x in all_samples_type if x.split("_")[1]==s][0] for s in type_to_samples_discard[drug_type]}
        strange_samples_discard = samples_discard.difference(all_samples_type)
        if len(strange_samples_discard)>0:
            raise ValueError("strange samples discard: %s"%strange_samples_discard)

        all_samples_type = all_samples_type.difference(samples_discard)

        extra_samples = all_samples_type.difference(all_samples)
        if len(extra_samples)>0:
            raise ValueError("some extra samples in phenotype: %s"%extra_samples)

        # write tree with these samples
        tree = cp.deepcopy(strains_tree)
        tree.prune(all_samples_type, preserve_branch_length=True)
        treefile = plots_dir + "/tree_pruned.nw"
        tree.write(outfile=treefile, format=2)

        # write the phenotypes
        phenotypes_file = plots_dir + "/phenos.tab"
        df_phenotypes = df_metadata_binary_phenotypes[["runID"] + sorted_phenotypes].rename(columns={"runID":"sampleID"}).copy()
        df_phenotypes = df_phenotypes[df_phenotypes.sampleID.isin(all_samples_type)].copy().reset_index(drop=True)
        df_phenotypes = df_phenotypes.rename(columns={pheno : get_paper_dicotomic_pheno(pheno).replace(" ", "_") for pheno in sorted_phenotypes})
        save_df_as_tab(df_phenotypes, phenotypes_file)

        # gets the pheno mapping
        df_pheno_mapping = pd.DataFrame({pheno : {"phenotype":get_paper_dicotomic_pheno(pheno).replace(" ", "_"), "short_phenotype":pheno_to_short_p[pheno]} for pheno in sorted_phenotypes}).transpose().reset_index(drop=True)
        phenotype_info_file = plots_dir + "/pheno_info.tab"
        save_df_as_tab(df_pheno_mapping, phenotype_info_file)

        # define the gene IDs
        gene_names = type_to_genes[drug_type]

        # discard some manually curated (no vars, none correlated to any rep samples phenos)
        genes_discard_dict = {"azole":{"ERG3"}, "echinocandin":{"ERG3"}, "AMB":{"ERG6", "ERG3"}, "5FC":{"Scer_FCY1", "Scer_CHS7"}}

        if drug_type in genes_discard_dict:
            gene_names = [g for g in gene_names if not g in genes_discard_dict[drug_type]]

        # get the genes
        gene_IDs = [gname_to_gID[gname] for gname in gene_names]
        gID_to_name = dict(zip(gene_IDs, gene_names))
        print(gID_to_name)
    
        # define the group info
        all_sorted_groups = ["all_vars|all_genes|non_syn_muts|genes|||%s"%gID for gID in gene_IDs]
        df_group_info = pd.DataFrame({"unique_groupID":all_sorted_groups})
        df_group_info["group_color"] = df_group_info.unique_groupID.apply(lambda x: x.split("|||")[-1]).map(gID_to_name).map(get_gene_to_color())
        group_info_file = "%s/group_info.tab"%(plots_dir)
        save_df_as_tab(df_group_info, group_info_file)

        # get the sample ID mapping
        sorted_samples = Tree(treefile).get_leaf_names()
        df_samples_IDmapping = pd.DataFrame({"sampleID":sorted_samples})
        df_samples_IDmapping["sample_name"] = df_samples_IDmapping.sampleID.apply(lambda x: x.split("_")[1]) # get sample ID
        df_samples_IDmapping["sample_label_color"] = "black"

        df_clade_info = df_metadata[["runID", "cladeID_based_on_Bergin2022"]].copy().rename(columns={"runID":"sampleID", "cladeID_based_on_Bergin2022":"clade"})
        df_samples_IDmapping = df_samples_IDmapping.merge(df_clade_info, on="sampleID", how="left", validate="one_to_one")
        df_samples_IDmapping["clade"] = df_samples_IDmapping["clade"].apply(get_str_cladeID_Bergin2022)
        sampleID_mapping_file = "%s/samples.tab"%plots_dir
        save_df_as_tab(df_samples_IDmapping, sampleID_mapping_file)

        # run
        outdir_p = "%s/tree_plot"%(plots_dir)
        if file_is_empty(outdir_p + "/pipeline_finished.txt"):

            cmd = "%s/scripts/ancestral_GWAS/scripts/plot_tree_variants_multiple_phenos.py --outdir %s --variants_file %s --variant_annotations_file %s --phenotypes_file %s --treefile %s --group_info %s --min_support 70 --ASR_method DOWNPASS --classification_method_ASR all_ignoreNaN --gff %s --ref_genome %s --pathway_annotations_dir %s --CGD_chromosomal_features %s --interproscan_output %s --sampleID_mapping %s --custom_geneID_maps %s --threads 112 --mtDNA_chrom %s --gDNA_code %i --mtDNA_code %i --target_geneIDs %s --only_vars_correlated_phenotype --phenotype_info %s"%(ParentDir, outdir_p, variants_file, variant_annotations_file, phenotypes_file, treefile, group_info_file, gff_Cparapsilosis, genome_Cparapsilosis, pathway_annotations_dir, CGD_chromosomal_features, interproscan_output, sampleID_mapping_file, mapping_genes, mtDNA_chrom, Cpara_gDNA_code, Cpara_mtDNA_code, ",".join(gene_IDs), phenotype_info_file) # title_str

            run_cmd(cmd, env="ancestral_GWAS_env")   # --only_vars_correlated_phenotype


            #"%s/prepare_GWAS_inputs.py --outdir %s --phenotypes %s --treefile %s --threads 112 --integrated_small_vars %s --small_vars_filtering typical_diploid --repeats_file %s --small_vars_min_cov_pos %i --small_vars_annot %s --integrated_SV_CNVs %s --SV_CNVs_annot %s --gff %s --ref_genome %s %s --paths_CNV_calling %s %s --interproscan_output %s --pathway_annotations_dir %s --mtDNA_chrom %s --gDNA_code %i --mtDNA_code %i --min_nmuts_per_group 2 --min_nsamples_per_group 2 --CGD_chromosomal_features %s --reactome_organism_priority Saccharomyces_cerevisiae,Schizosaccharomyces_pombe"%(ancestral_GWAS_dir, outdir_gwas_inputs, phenotypes_file, treefile, integrated_small_vars, simple_repeats_file, min_coverage_pos, small_vars_annot, integrated_SV_CNVs, SV_CNVs_annot, gff_Cparapsilosis, genome_Cparapsilosis, SV_CNV_filters, paths_CNV_calling, aneuploidy_calling_filters, interproscan_output, pathway_annotations_dir, mtDNA_chrom, Cpara_gDNA_code, Cpara_mtDNA_code, CGD_chromosomal_features), 



def plot_trees_expected_genes_variants(DataDir, PlotsDir, df_metadata, samples_discard_run1, strains_tree, df_metadata_binary_phenotypes, target_sample_set):

    """Get the plots of the trees of expected genes """

    # define dirs
    PlotsDir_trees = "%s/trees_expected_genes"%PlotsDir; make_folder(PlotsDir_trees)

    # define expected genes
    type_to_genes = {"azole": ["ERG11", "ERG3", "MRR1", "Scer_HAL9", "NDT80", "UPC2"], # UPC2 has no vars
                     "echinocandin": ["FKS1", "FKS2", "ERG3"],
                     "AMB": ["ERG2", "ERG3", "Scer_ERG5", "ERG6", "ERG25"],
                     "5FC": ["Scer_FCY1", "Scer_FCY2", "Scer_FUR1", "Scer_CHS7"]}

    mapping_genes = "Scer_HAL9:TAC1,Scer_ERG5:ERG5,Scer_FCY1:FCY1,Scer_FCY2:FCY2,Scer_FUR1:FUR1,Scer_CHS7:CHS7"
         

    # define the phenos
    azoles = ['FLC_resistance', 'VRC_resistance', 'ITR_resistance', 'POS_resistance', 'IVZ_resistance']
    echinocandins = ['CAS_resistance', 'MIF_resistance', 'ANI_resistance']
    resistance_phenos =  azoles + echinocandins + ['5FC_resistance', 'AMB_resistance']

    # load gene feats
    gene_features_df = get_gene_features_df(DataDir)
    gname_to_gID = gene_features_df.set_index("final_name").geneID

    # one plot for each pheno
    for phenotype in resistance_phenos:
        print(phenotype)
        drug_name = phenotype.split("_")[0]

        # debug !!!!!
        #if phenotype!="FLC_resistance": continue

        # define the genes
        if drug_name in type_to_genes: 
            gene_names = type_to_genes[drug_name]

        elif phenotype in azoles:
            gene_names = type_to_genes["azole"]

        elif phenotype in echinocandins:
            gene_names = type_to_genes["echinocandin"]

        else:
            raise ValueError("gene names")

        # discard some manually curated (no variants, or none correlated to pheno), in rep. samples
        genes_discard_dict = {"FLC_resistance":{"UPC2", "ERG3"}, 
                              "VRC_resistance":{"ERG3", "UPC2"},
                              "ITR_resistance":{"ERG3", "MRR1", "UPC2"},
                              "POS_resistance":{"ERG3", "UPC2"},
                              "IVZ_resistance":{"ERG3", "UPC2"},

                              "CAS_resistance":{"ERG3", "FKS2"},
                              "MIF_resistance":{"ERG3", "FKS2"},
                              "ANI_resistance":{"ERG3", "FKS2"},

                              "5FC_resistance":{"Scer_CHS7", "Scer_FUR1", "Scer_FCY1", "Scer_FCY2"},
                              "AMB_resistance":{"ERG3", "ERG6", "ERG2"}

        }

        if phenotype in genes_discard_dict:
            gene_names = [g for g in gene_names if not g in genes_discard_dict[phenotype]]

        # debug
        if len(gene_names)==0: continue


        # define the gene IDs
        gene_IDs = [gname_to_gID[gname] for gname in gene_names]
        gID_to_name = dict(zip(gene_IDs, gene_names))
        print(gID_to_name)

        # define dirs
        aligner_pheno = "bwa_mem"
        geno_pheno_dir = "%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas"%(DataDir, phenotype, aligner_pheno)
        groups_file = "%s/gwas_inputs_all_samples/all_groupings.tab"%(geno_pheno_dir)
        variants_file = "%s/gwas_inputs_all_samples/variants_df_all.tab"%(geno_pheno_dir)
        variant_annotations_file = "%s/gwas_inputs_all_samples/variant_annotations_df_all.tab"%(geno_pheno_dir)

        pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
        CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
        genome_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir # genome from CGD
        gff_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir
        interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir

        treefile = "%s/gwas_inputs_%s/corrected_rooted_tree.nw"%(geno_pheno_dir, target_sample_set)
        phenotypes_ASR = "%s/gwas_inputs_%s/phenotypes_ASR.out"%(geno_pheno_dir, target_sample_set)

        # define dir where to save all
        plots_dir = "%s/%s_%s"%(PlotsDir_trees, phenotype, target_sample_set); make_folder(plots_dir)

        #df_groups_pheno = get_df_groups_only_subset_genes(groups_file, "%s/df_groups_all.py"%plots_dir, set(gene_IDs))
        all_sorted_groups = ["all_vars|all_genes|non_syn_muts|genes|||%s"%gID for gID in gene_IDs]
        df_groups_pheno = get_df_groups_only_subset(groups_file, "%s/df_groups_all.py"%plots_dir, set(all_sorted_groups))
        sorted_vars_pheno = sorted(set(df_groups_pheno.variantID_across_samples))        
        df_vars_IDmapping = pd.DataFrame({"variantID_across_samples":sorted_vars_pheno}).reset_index().rename(columns={"index":"var_idx"})
        df_vars_IDmapping["var_name"] = df_vars_IDmapping.variantID_across_samples.apply(get_simple_var_name_variantID_across_samples)
        vars_IDmapping_file = "%s/varsID_mapping.tab"%(plots_dir)
        save_df_as_tab(df_vars_IDmapping, vars_IDmapping_file)

        # get the colors related to the predictors
        df_group_info = pd.DataFrame({"unique_groupID":sorted(set(df_groups_pheno.unique_groupID))})
        df_group_info["group_color"] = df_group_info.unique_groupID.apply(lambda x: x.split("|||")[-1]).map(gID_to_name).map(get_gene_to_color())
        group_info_file = "%s/group_info.tab"%(plots_dir)
        save_df_as_tab(df_group_info, group_info_file)

        # get the sample ID mapping
        sorted_samples = Tree(treefile).get_leaf_names()
        df_samples_IDmapping = pd.DataFrame({"sampleID":sorted_samples})
        df_samples_IDmapping["sample_name"] = df_samples_IDmapping.sampleID.apply(lambda x: x.split("_")[1]) # get sample ID
        df_samples_IDmapping["sample_label_color"] = "black"

        df_clade_info = df_metadata[["runID", "cladeID_based_on_Bergin2022"]].copy().rename(columns={"runID":"sampleID", "cladeID_based_on_Bergin2022":"clade"})
        df_samples_IDmapping = df_samples_IDmapping.merge(df_clade_info, on="sampleID", how="left", validate="one_to_one")
        df_samples_IDmapping["clade"] = df_samples_IDmapping["clade"].apply(get_str_cladeID_Bergin2022)
        sampleID_mapping_file = "%s/samples.tab"%plots_dir
        save_df_as_tab(df_samples_IDmapping, sampleID_mapping_file)

        # define a dir to save
        dir_only_plots = "%s/trees"%plots_dir; make_folder(dir_only_plots)

        # one plot for each gene
        for geneID, gene_name in gID_to_name.items():

            if sum(df_groups_pheno["group"]==geneID)==0:
                raise ValueError("group not in df")

            # groups gID
            group_gID = df_groups_pheno[df_groups_pheno["group"]==geneID].iloc[0].unique_groupID


            # plot predictors results
            outdir_p = "%s/plots_%s"%(plots_dir, gene_name)
            #delete_folder(outdir_p)
            if file_is_empty(outdir_p + "/pipeline_finished.txt"):
                print("Runing cmd...")
                cmd = "%s/scripts/ancestral_GWAS/scripts/plot_tree_variants.py --outdir %s --threads 112 --target_groups \"%s\" --min_support 70 --ASR_method DOWNPASS --variants_file %s --variant_annotations_file %s --groups_file %s --ASR_phenotypes %s --treefile %s --gff %s --ref_genome %s --pathway_annotations_dir %s --CGD_chromosomal_features %s --classification_method_ASR all_ignoreNaN --vars_IDmapping %s --custom_geneID_maps %s --interproscan_output %s --sampleID_mapping %s --only_vars_correlated_phenotype --group_info %s"%(ParentDir, outdir_p, group_gID, variants_file, variant_annotations_file, groups_file, phenotypes_ASR, treefile, gff_Cparapsilosis, genome_Cparapsilosis, pathway_annotations_dir, CGD_chromosomal_features, vars_IDmapping_file, mapping_genes, interproscan_output, sampleID_mapping_file, group_info_file) # optimal_F1, optimal_NegPredVal --representative_samples %s 

                run_cmd(cmd, env="ancestral_GWAS_env")   # --only_vars_correlated_phenotype

            # move tree plot to dir_only_plots
            rsync_file("%s/tree_variants.pdf"%outdir_p, "%s/%s_vars.pdf"%(dir_only_plots, gene_name))

def get_genes_per_groupID(gID):

    """get genes for one group"""

    group_info, group_name = gID.split("|||")
    
    # per each coll
    type_collapsing  = group_info.split("|")[-1]
    if type_collapsing=="genes":
        return {group_name}

    elif type_collapsing=="domains":
        return {group_name.split("#")[1]}

    elif type_collapsing in {"none", "GO", "Reactome", "MetaCyc"}:
        return set()

    else:
        raise ValueError(type_collapsing)

def get_list_colors_different_from_the_ones_to_avoid(set_colors_avoid, cbar, ncolors, similarity_threshold):

    """
    Gets list different to the ones to avoid
    """

    # convert to rgb
    set_colors_avoid_rgb = set(map(mcolors.hex2color, set_colors_avoid))

    # define a function that calculates color distance
    def color_distance(c1, c2):
        return np.sqrt(sum((x - y) ** 2 for x, y in zip(c1, c2)))

    # get the colors of the palette as RGB tuples
    all_colors_palette_rgb = sns.color_palette(cbar, ncolors + len(set_colors_avoid) + 2)

    # get the distinct_colors
    distinct_colors = [c for c in all_colors_palette_rgb if all([color_distance(input_rgb, c)>similarity_threshold for input_rgb in set_colors_avoid_rgb])]

    # check
    if len(distinct_colors)<ncolors:
        raise ValueError("%i expected colors, only %i found"%(ncolors, len(distinct_colors)))

    return [rgb_to_hex(x) for x in  distinct_colors[0:ncolors]]

def get_df_group_info_plotting_ML_models(set_target_groups, gene_name_to_ID, gene_to_color=None):

    """Gets group info"""

    # map gene to color
    if gene_to_color is None:
        gene_to_color = get_gene_to_color()

    geneID_to_color = {gene_name_to_ID[gname] : color for gname,color in gene_to_color.items()}

    # init df
    df_group_info = pd.DataFrame({"unique_groupID":sorted(set_target_groups)})

    # add genes
    df_group_info["set_geneIDs"] = df_group_info.unique_groupID.apply(get_genes_per_groupID)

    # init color with the ones from geneID_to_color
    def get_color_from_gID(set_geneIDs):
        gene_IDs_defined = set_geneIDs.intersection(set(geneID_to_color))
        if len(gene_IDs_defined)==0:
            return ""

        elif len(gene_IDs_defined)==1:
            return geneID_to_color[next(iter(gene_IDs_defined))]

        else:
            raise ValueError("invalid")

    df_group_info["group_color"] = df_group_info.set_geneIDs.apply(get_color_from_gID)

    # get color mapping remaining groups
    groups_no_color = list(df_group_info[df_group_info.group_color==""].unique_groupID)

    if len(groups_no_color)<8:
        cbar = "tab10"

    elif len(groups_no_color)<18:
        cbar = "tab20"

    else:
        raise ValueError("too many colors")

    set_colors_avoid = set(df_group_info[df_group_info.group_color!=""].group_color)
    similarity_threshold = 0.5 # higher, more threshold allowed
    colors_different = get_list_colors_different_from_the_ones_to_avoid(set_colors_avoid, cbar, len(groups_no_color), similarity_threshold)
    gID_to_color = dict(zip(groups_no_color, colors_different))

    def get_color_added(r):
        if r.group_color=="":
            return gID_to_color[r.unique_groupID]

        else:
            return r.group_color

    df_group_info["group_color"] = df_group_info.apply(get_color_added, axis=1)

    # return
    return df_group_info

def visualize_top_ML_models_in_tree_all_genes_together(ML_df_dicts, target_phenos_ML, PlotsDir, DataDir, df_metadata_binary_phenotypes, df_metadata, test_sample_set):

    """gets the top models' trees, only variants correlated to pheno"""

    # get filtered df
    filters_dict = {"phenotype":target_phenos_ML, "test_sample_set":test_sample_set}
    df_ML_results_all = get_filtered_df_by_dict(ML_df_dicts["df_ML_results_top"], filters_dict)
    #df_ML_metrics_all = get_filtered_df_by_dict(ML_df_dicts["df_ML_metrics_top"], filters_dict)

    # general vars
    pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
    CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
    genome_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir # genome from CGD
    gff_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir
    interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir

    # define gene_to_color_ML_genes
    gene_to_color_ML_genes = {g:c for g,c in get_gene_to_color().items() if g in {"ERG11", "MRR1"}}

    other_genes = ['CDR1B', 'Scer_ERG24', 'CPAR2_204120', 'Scer_ATP10', 'Scer_MOH1'] + ['Scer_CLU1', 'STP4'] + ['ALS11', 'ALS7', 'Scer_MDS3'] + ['RBT1', 'Scer_PMI40']
    other_colors = ["tab:purple", "tab:green", "tab:olive", "teal", "lime"] + ["tab:brown", "tab:orange"] + ["tab:red", "tab:pink", "magenta"] + ["gold", "salmon"]

    for I in range(len(other_genes)):
        gene_to_color_ML_genes[other_genes[I]] = mcolors.to_hex(other_colors[I])


    # get clade info        
    df_clade_info = df_metadata[["runID", "cladeID_based_on_Bergin2022"]].copy().rename(columns={"runID":"sampleID", "cladeID_based_on_Bergin2022":"clade"})

    # get names
    gene_features_df = get_gene_features_df(DataDir)
    gene_name_to_ID = dict(gene_features_df.set_index("final_name").geneID)

    geneID_to_name = dict(gene_features_df.set_index("geneID").final_name)
    tac1_ID = gene_features_df[gene_features_df.final_name=="Scer_HAL9"].iloc[0].geneID
    geneID_to_name[tac1_ID] = "TAC1"

    # define the mapping of genes
    mapping_genes = "Scer_HAL9:TAC1,Scer_ERG5:ERG5,Scer_FCY1:FCY1,Scer_FCY2:FCY2,Scer_FUR1:FUR1,Scer_CHS7:CHS7" 

    # add genes
    df_ML_results_all["selected_predictors"] = df_ML_results_all.selected_predictors.apply(str).apply(lambda x: set(x.split(","))).apply(get_selected_predictors_demult_one_tuple_predictors)
    df_ML_results_all["all_affected_genes_set"] = df_ML_results_all.selected_predictors.apply(lambda x: set(map(get_affected_geneID_one_unique_groupID, x)))

    # plots for each pheno / model
    for phenotype, df_ML_results_pheno in df_ML_results_all.groupby("phenotype"):

        # define the target genes for this pheno
        geneID_to_nmodels = pd.Series({gID : sum(df_ML_results_pheno.all_affected_genes_set.apply(lambda x: gID in x)) for gID in sorted(set.union(*df_ML_results_pheno.all_affected_genes_set))})
        target_geneIDs = set(geneID_to_nmodels[geneID_to_nmodels>=3].index)

        # extra filters of aligner, based on the manual curation. I checked that the main genes are there, all are the same as best pheno
        if phenotype=="VRC_resistance": aligner_pheno = "hisat2_no_spliced"
        elif phenotype=="FLC_resistance": aligner_pheno = "hisat2_no_spliced"
        elif phenotype=="age_infant": aligner_pheno = "bowtie2_local"
        elif phenotype=="source_blood_vs_other_sterile": aligner_pheno = "hisat2_no_spliced" # in best models
        else: raise ValueError("invalid pheno")
        
        aligners_keep = {aligner_pheno}
        df_ML_results = df_ML_results_pheno[df_ML_results_pheno.aligner.isin(aligners_keep)].copy()

        # redefine genes
        df_ML_results["all_affected_genes_set"] = df_ML_results["all_affected_genes_set"].apply(lambda x: x.intersection(target_geneIDs))

        # print missing genes
        genes_filt = set.union(*df_ML_results.all_affected_genes_set)
        missing_genes = target_geneIDs.difference(genes_filt)
        print(phenotype, "missing genes", [geneID_to_name[g] for g in sorted(missing_genes)])
        print(phenotype, "used genes", [geneID_to_name[g] for g in sorted(genes_filt)])

        # debug missing genes
        if phenotype=="source_blood_vs_other_sterile":
            if {geneID_to_name[g] for g in missing_genes}!={'BCR1', 'Scer_DFG16'}: raise ValueError("invalid missings")
        else:
            if len(missing_genes)!=0: raise ValueError("should not be missing genes, %s"%missing_genes)

        # redefine predictors to only have genes_filt
        df_ML_results["selected_predictors"] =  df_ML_results["selected_predictors"].apply(lambda x: {p for p in x if get_affected_geneID_one_unique_groupID(p) in genes_filt})

        # define a df with all groups
        df_groups_initial = pd.DataFrame({"unique_groupID":sorted(set.union(*df_ML_results["selected_predictors"]))})
        df_groups_initial["n_models_feature"] = df_groups_initial.unique_groupID.apply(lambda f: sum(df_ML_results["selected_predictors"].apply(lambda x: f in x)))
        df_groups_initial["gene_name"] = df_groups_initial.unique_groupID.apply(get_affected_geneID_one_unique_groupID).map(geneID_to_name)
        #print(df_groups_initial.sort_values(by=["gene_name", "n_models_feature"], ascending=False)[["gene_name", "n_models_feature"]])

        # get one group per each
        def get_most_common_group_gene(df):
            df = get_df_like_first_row(df, ["n_models_feature"], [False])
            if len(df)!=1: raise ValueError("should be sortable by # models")
            return df.iloc[0]
        df_groups = df_groups_initial.groupby("gene_name").apply(get_most_common_group_gene)

        # get the inputs for plotting
        geno_pheno_dir = "%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas"%(DataDir, phenotype, aligner_pheno)
        variants_file = "%s/gwas_inputs_all_samples/variants_df_all.tab"%(geno_pheno_dir)
        variant_annotations_file = "%s/gwas_inputs_all_samples/variant_annotations_df_all.tab"%(geno_pheno_dir)

        # define dirs
        p_dir_all = "%s/trees_ML_genes_%s"%(PlotsDir, test_sample_set); make_folder(p_dir_all)
        p_dir_pheno = "%s/%s"%(p_dir_all, phenotype); make_folder(p_dir_pheno)

        # get the tree, depending on train / test
        ss_to_treefile = {ss : "%s/gwas_inputs_%s/corrected_rooted_tree.nw"%(geno_pheno_dir, ss) for ss in ["all_samples", "representative_samples"]}
        if test_sample_set=="all_samples":
            treefile_plot = ss_to_treefile["all_samples"]

        elif test_sample_set=="representative_samples" and phenotype in {"age_infant", "source_blood_vs_other_sterile"}:
            treefile_plot = ss_to_treefile["representative_samples"]

        elif test_sample_set=="representative_samples" and phenotype in {"FLC_resistance", "VRC_resistance"}: # use specific aligner

            samples_aligner = {"FLC_resistance":"bwa_mem", "VRC_resistance":"bwa_mem"}[phenotype]
            treefile_plot = "%s/tree_plot.nw"%p_dir_pheno
            tree = Tree(ss_to_treefile["all_samples"])
            rep_samples = set(Tree("%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas/gwas_inputs_representative_samples/corrected_rooted_tree.nw"%(DataDir, phenotype, samples_aligner)).get_leaf_names()) # use the ones of bwa mem, for comparability w previous analyses
            tree.prune(rep_samples,  preserve_branch_length=True)
            tree.write(outfile=treefile_plot, format=2)


        else:
            raise ValueError("invalid combination")

        # define a df with the groups, and their info
        df_group_info = get_df_group_info_plotting_ML_models(set(df_groups.unique_groupID), gene_name_to_ID, gene_to_color=gene_to_color_ML_genes)
        df_group_info["groupID_interaction"] = "-"
        check_no_nans_in_df(df_group_info)
        group_info_file = "%s/group_info.tab"%(p_dir_pheno)
        save_df_as_tab(df_group_info[["unique_groupID", "group_color", "groupID_interaction"]], group_info_file)

        # get the sample ID mapping
        all_samples_sorted = Tree(treefile_plot).get_leaf_names()
        df_samples_IDmapping = pd.DataFrame({"sampleID":all_samples_sorted})
        df_samples_IDmapping["sample_name"] = df_samples_IDmapping.sampleID.apply(lambda x: x.split("_")[1][1:]) # get sample ID
        sname_to_sID = df_samples_IDmapping.set_index("sample_name").sampleID
        df_samples_IDmapping["sample_name"] = "L"+df_samples_IDmapping["sample_name"]
        df_samples_IDmapping = df_samples_IDmapping.merge(df_clade_info, on="sampleID", how="left", validate="one_to_one")
        df_samples_IDmapping["clade"] = df_samples_IDmapping["clade"].apply(get_str_cladeID_Bergin2022)

        if phenotype=="source_blood_vs_other_sterile":

            mapping_dict = {"whole blood":"blood", "blood culture":"blood", "cerebrospinal fluid":"cerebrospinal", "ascitic fluid":"ascitic", "bile fluid":"bile"}
            sID_to_source = df_metadata[df_metadata.runID.isin(set(df_samples_IDmapping.sampleID))].muestra.apply(lambda x: mapping_dict[clincial_sources_translated["source_"+x]])
            df_samples_IDmapping["sample_name"] = df_samples_IDmapping.apply(lambda r: "%s (%s)"%(r.sample_name, sID_to_source[r.sampleID]), axis=1)

        sampleID_mapping_file = "%s/samples.tab"%p_dir_pheno
        save_df_as_tab(df_samples_IDmapping, sampleID_mapping_file)

        # get the df with phenotypes
        df_phenotypes = df_metadata_binary_phenotypes[["runID", phenotype]].rename(columns={"runID":"sampleID"}).copy()
        df_phenotypes = df_phenotypes[df_phenotypes.sampleID.isin(all_samples_sorted)].copy().reset_index(drop=True)
        df_phenotypes = df_phenotypes.rename(columns={pheno : get_paper_dicotomic_pheno(pheno).replace(" ", "_") for pheno in [phenotype]})
        """
        if phenotype=="source_blood_vs_other_sterile":
            df_phenotypes["other_sterile_vs_blood"] = df_phenotypes.blood_vs_other_sterile.map({1:0, 0:1})
            df_phenotypes = df_phenotypes.drop(["blood_vs_other_sterile"], axis=1)

        """
        phenotypes_file = "%s/phenos.tab"%p_dir_pheno
        save_df_as_tab(df_phenotypes, phenotypes_file)

        # get a phenotype info file
        phenotype_f = sorted(set(df_phenotypes.keys()).difference({"sampleID"}))[0]
        df_pheno_mapping = pd.DataFrame({0 : {"phenotype":phenotype_f, "short_phenotype":"-", "use_for_correlations":True}}).transpose()
        phenotype_info_file = p_dir_pheno + "/pheno_info.tab"
        save_df_as_tab(df_pheno_mapping, phenotype_info_file)

        # get IDs
        target_geneIDs_str = ",".join(sorted(genes_filt))

        # define extra args
        if phenotype in {"FLC_resistance", "VRC_resistance", "source_blood_vs_other_sterile"}: extra_args = "--only_vars_correlated_phenotype"
        else: extra_args = ""
        
        # plot predictors results
        outdir_p = "%s/plot_results"%(p_dir_pheno)
        #delete_folder(outdir_p)

        if file_is_empty(outdir_p + "/pipeline_finished.txt"):
            print("Runing cmd...")
            
            cmd = "%s/scripts/ancestral_GWAS/scripts/plot_tree_variants_multiple_phenos.py --outdir %s --variants_file %s --variant_annotations_file %s --phenotypes_file %s --treefile %s --group_info %s --min_support 70 --ASR_method DOWNPASS --classification_method_ASR all_ignoreNaN --gff %s --ref_genome %s --pathway_annotations_dir %s --CGD_chromosomal_features %s --interproscan_output %s --sampleID_mapping %s --custom_geneID_maps %s --threads 112 --mtDNA_chrom %s --gDNA_code %i --mtDNA_code %i   --phenotype_info %s --target_geneIDs %s %s"%(ParentDir, outdir_p, variants_file, variant_annotations_file, phenotypes_file, treefile_plot, group_info_file, gff_Cparapsilosis, genome_Cparapsilosis, pathway_annotations_dir, CGD_chromosomal_features, interproscan_output, sampleID_mapping_file, mapping_genes, mtDNA_chrom, Cpara_gDNA_code, Cpara_mtDNA_code, phenotype_info_file, target_geneIDs_str, extra_args)

            run_cmd(cmd, env="ancestral_GWAS_env")   # --only_vars_correlated_phenotype


def plot_GWAS_hits_in_trees(df_GWAS_hits_NR, df_GWAS_hits_all, PlotsDir, df_metadata_binary_phenotypes, df_metadata, DataDir):

    """plot GWAS hits in trees"""

    # keep
    print("visualizing GWAS hits")
    df_GWAS_hits_NR = df_GWAS_hits_NR.copy()
    df_GWAS_hits_all = df_GWAS_hits_all.copy()

    # general vars
    pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
    CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
    genome_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir # genome from CGD
    gff_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir
    interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir
    mapping_genes = "Scer_HAL9:TAC1,Scer_ERG5:ERG5,Scer_FCY1:FCY1,Scer_FCY2:FCY2,Scer_FUR1:FUR1,Scer_CHS7:CHS7" 

    # get clade info        
    df_clade_info = df_metadata[["runID", "cladeID_based_on_Bergin2022"]].copy().rename(columns={"runID":"sampleID", "cladeID_based_on_Bergin2022":"clade"})

    # get names
    gene_features_df = get_gene_features_df(DataDir)
    gene_name_to_ID = dict(gene_features_df.set_index("final_name").geneID)
    geneID_to_name = dict(gene_features_df.set_index("geneID").final_name)

    # get the phenos
    sorted_phenotypes_hc = sorted(set(df_GWAS_hits_all[df_GWAS_hits_all.type_filtering=="high_confidence"].phenotype))

    # iterate through cases
    for type_filtering in ["high_confidence", "low_confidence"]:
        p_dir_all = "%s/trees_GWAS_hits_%s"%(PlotsDir, type_filtering); make_folder(p_dir_all)

        for phenotype in sorted_phenotypes_hc:
            phenotype_paper = get_paper_dicotomic_pheno(phenotype)

            # get the hits
            df_hits_NR = get_filtered_df_by_dict(df_GWAS_hits_NR, {"phenotype":phenotype_paper, "type_filtering":type_filtering})
            df_hits_all = get_filtered_df_by_dict(df_GWAS_hits_all, {"phenotype":phenotype, "type_filtering":type_filtering})

            # define the filters
            filters_dict = {"high_confidence":manual_curation_gwas_filt_parms, "low_confidence":manual_curation_gwas_filt_parms_low_conf}[type_filtering]

            # get the inputs for plotting
            geno_pheno_dir = "%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas"%(DataDir, phenotype, filters_dict["aligner"])
            variants_file = "%s/gwas_inputs_%s/variants_df_all.tab"%(geno_pheno_dir, filters_dict["sample_set"])
            variant_annotations_file = "%s/gwas_inputs_%s/variant_annotations_df_all.tab"%(geno_pheno_dir, filters_dict["sample_set"])
            treefile = "%s/gwas_inputs_%s/corrected_rooted_tree.nw"%(geno_pheno_dir, filters_dict["sample_set"])

            # plot for different group types
            for group_types_plot in ["all_SingleVariantHits"]:


                # get groups of all variants
                if group_types_plot=="all_SingleVariantHits":
                    target_groups = sorted(set(df_hits_all[df_hits_all.type_collapsing=="none"].unique_groupID))
                    #annot_df = get_tab_as_df_or_empty_df(variant_annotations_file)
                    #target_vars = {g.split("|||")[1] for g in target_groups}
                    #affected_geneIDs = set(annot_df[annot_df.variantID_across_samples.isin(set(target_vars))].Gene)
                    affected_geneIDs = None

                else: raise ValueError("invalid group")
                if len(target_groups)==0: continue
                print("plotting for", type_filtering, phenotype, group_types_plot, len(target_groups))

                # discard some cases
                if group_types_plot=="all_SingleVariantHits" and type_filtering=="low_confidence" and not phenotype in  {"AMB_resistance"}: continue # source_LIQUIDO_ASCITICO, source_blood
                if group_types_plot=="all_SingleVariantHits" and type_filtering=="high_confidence" and not len(target_groups)>=3: continue

                # define dirs
                p_dir_pheno = "%s/%s_%s"%(p_dir_all, phenotype.replace(" ", "_"), group_types_plot); make_folder(p_dir_pheno)

                # get the group info
                df_group_info = pd.DataFrame({"unique_groupID":target_groups, "groupID_interaction":"-"})
                df_group_info["group_color"] = df_group_info.unique_groupID.map(get_value_to_color(target_groups, palette="tab10", type_color="hex")[0])
                group_info_file = "%s/group_info.tab"%(p_dir_pheno)
                save_df_as_tab(df_group_info[["unique_groupID", "group_color", "groupID_interaction"]], group_info_file)

                # define affected genes
                if not affected_geneIDs is None: target_geneIDs_str = "--target_geneIDs " + ",".join(sorted(affected_geneIDs))
                else: target_geneIDs_str = " "

                # get the sample ID mapping
                all_samples_sorted = Tree(treefile).get_leaf_names()
                df_samples_IDmapping = pd.DataFrame({"sampleID":all_samples_sorted})
                df_samples_IDmapping["sample_name"] = df_samples_IDmapping.sampleID.apply(lambda x: x.split("_")[1][1:]) # get sample ID
                df_samples_IDmapping["sample_name"] = "L"+df_samples_IDmapping["sample_name"]
                df_samples_IDmapping = df_samples_IDmapping.merge(df_clade_info, on="sampleID", how="left", validate="one_to_one")
                df_samples_IDmapping["clade"] = df_samples_IDmapping["clade"].apply(get_str_cladeID_Bergin2022)
                sampleID_mapping_file = "%s/samples.tab"%p_dir_pheno
                save_df_as_tab(df_samples_IDmapping, sampleID_mapping_file)

                # get the df with phenotypes
                df_phenotypes = df_metadata_binary_phenotypes[["runID", phenotype]].rename(columns={"runID":"sampleID"}).copy()
                df_phenotypes = df_phenotypes[df_phenotypes.sampleID.isin(all_samples_sorted)].copy().reset_index(drop=True)
                df_phenotypes = df_phenotypes.rename(columns={pheno : get_paper_dicotomic_pheno(pheno).replace(" ", "_") for pheno in [phenotype]})
                check_no_nans_in_df(df_phenotypes)
                phenotypes_file = "%s/phenos.tab"%p_dir_pheno
                save_df_as_tab(df_phenotypes, phenotypes_file)

                # get a phenotype info file
                df_pheno_mapping = pd.DataFrame({pheno : {"phenotype":pheno, "short_phenotype":"-"} for pheno in list(df_phenotypes.keys()) if pheno!="sampleID"}).transpose().reset_index(drop=True)
                df_pheno_mapping["use_for_correlations"] = df_pheno_mapping.phenotype==(get_paper_dicotomic_pheno(phenotype).replace(" ", "_"))
                phenotype_info_file = p_dir_pheno + "/pheno_info.tab"
                save_df_as_tab(df_pheno_mapping, phenotype_info_file)

                # define a title
                title_str = "%s GWAS hits for '%s'. %s"%(type_filtering.replace("_", " "), get_paper_dicotomic_pheno(phenotype), group_types_plot)
 
                # plot results
                outdir_p = "%s/plot_results"%(p_dir_pheno)
                if file_is_empty(outdir_p + "/pipeline_finished.txt"):
                    print("Runing cmd...")
                    
                    cmd = "%s/scripts/ancestral_GWAS/scripts/plot_tree_variants_multiple_phenos.py --outdir %s --variants_file %s --variant_annotations_file %s --phenotypes_file %s --treefile %s --group_info %s --min_support 70 --ASR_method MPPA,DOWNPASS --classification_method_ASR all_ignoreNaN --gff %s --ref_genome %s --pathway_annotations_dir %s --CGD_chromosomal_features %s --interproscan_output %s --sampleID_mapping %s --custom_geneID_maps %s --threads 112 --mtDNA_chrom %s --gDNA_code %i --mtDNA_code %i   --phenotype_info %s --title_str \"%s\" %s --exit_before_clean --variant_correlation_only_transitions"%(ParentDir, outdir_p, variants_file, variant_annotations_file, phenotypes_file, treefile, group_info_file, gff_Cparapsilosis, genome_Cparapsilosis, pathway_annotations_dir, CGD_chromosomal_features, interproscan_output, sampleID_mapping_file, mapping_genes, mtDNA_chrom, Cpara_gDNA_code, Cpara_mtDNA_code, phenotype_info_file, title_str, target_geneIDs_str) # title_str, --only_vars_correlated_phenotype

                    run_cmd(cmd, env="ancestral_GWAS_env")   # --only_vars_correlated_phenotype

def visualize_top_ML_models_in_trees(df_ML_results_best, df_ML_metrics_best, PlotsDir, DataDir, test_sample_set, df_metadata_binary_phenotypes, df_metadata):

    """Prints one tree with all the selected predictors"""

    # keep
    df_ML_results_best = df_ML_results_best.copy()

    # general vars
    pathway_annotations_dir = "%s/genomes_and_annotations/pathway_annotations"%DataDir
    CGD_chromosomal_features = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomal_feature.tab"%DataDir
    genome_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_chromosomes.fasta"%DataDir # genome from CGD
    gff_Cparapsilosis = "%s/genomes_and_annotations/C_parapsilosis_CDC317_version_s01-m06-r03_features.gff"%DataDir
    interproscan_output = "%s/genomes_and_annotations/interproscan_annot/interproscan_annotation.out"%DataDir

    # get clade info        
    df_clade_info = df_metadata[["runID", "cladeID_based_on_Bergin2022"]].copy().rename(columns={"runID":"sampleID", "cladeID_based_on_Bergin2022":"clade"})

    # get names
    gene_features_df = get_gene_features_df(DataDir)
    gene_name_to_ID = dict(gene_features_df.set_index("final_name").geneID)
    geneID_to_name = dict(gene_features_df.set_index("geneID").final_name)

    # define the mapping of genes
    mapping_genes = "Scer_HAL9:TAC1,Scer_ERG5:ERG5,Scer_FCY1:FCY1,Scer_FCY2:FCY2,Scer_FUR1:FUR1,Scer_CHS7:CHS7" 

    # plots for each pheno / model
    for phenotype, df_ML_results_best_pheno in df_ML_results_best[df_ML_results_best.test_sample_set==test_sample_set].copy().groupby("phenotype"):

        # define general parms
        train_sample_set_pheno = df_ML_results_best_pheno.iloc[0].sample_set
        aligner_pheno = df_ML_results_best_pheno.iloc[0].aligner

        # get df with metrics
        df_ML_metrics = df_ML_metrics_best[(df_ML_metrics_best.phenotype==phenotype) & (df_ML_metrics_best.eval_data=="cross_splits") & (df_ML_metrics_best.test_sample_set==test_sample_set)].copy()

        if len(df_ML_metrics)!=len(df_ML_metrics[["test_split_ID", "threshold"]].drop_duplicates()):
            raise ValueError("not unique threshold")

        # old files
        #groups_file = "%s/gwas_inputs_all_samples/all_groupings.tab"%(geno_pheno_dir)
        #phenotypes_ASR = "%s/gwas_inputs_%s/phenotypes_ASR.out"%(geno_pheno_dir, test_sample_set)

        # get the inputs for plotting
        geno_pheno_dir = "%s/genotype_phenotype_relationships/pheno-%s_aligner-%s/convergence_gwas"%(DataDir, phenotype, aligner_pheno)
        variants_file = "%s/gwas_inputs_all_samples/variants_df_all.tab"%(geno_pheno_dir)
        variant_annotations_file = "%s/gwas_inputs_all_samples/variant_annotations_df_all.tab"%(geno_pheno_dir)

        # define dirs
        p_dir_all = "%s/trees_variants_ML_classifiers_test_%s"%(PlotsDir, test_sample_set); make_folder(p_dir_all)
        p_dir_pheno = "%s/%s"%(p_dir_all, phenotype); make_folder(p_dir_pheno)

        # get the tree, depending on train / test
        ss_to_treefile = {ss : "%s/gwas_inputs_%s/corrected_rooted_tree.nw"%(geno_pheno_dir, ss) for ss in ["all_samples", "representative_samples"]}
        if train_sample_set_pheno=="all_samples" and test_sample_set=="all_samples":
            treefile_plot = ss_to_treefile["all_samples"]

        elif train_sample_set_pheno=="representative_samples" and test_sample_set=="representative_samples":
            treefile_plot = ss_to_treefile["representative_samples"]

        elif train_sample_set_pheno=="representative_samples" and test_sample_set=="all_samples":
            treefile_plot = ss_to_treefile["all_samples"]

        elif train_sample_set_pheno=="all_samples" and test_sample_set=="representative_samples":
            treefile_plot = "%s/tree_plot.nw"%p_dir_pheno
            tree = Tree(ss_to_treefile["all_samples"])
            tree.prune(set(Tree(ss_to_treefile["representative_samples"]).get_leaf_names()), preserve_branch_length=True)
            tree.write(outfile=treefile_plot, format=2)

        else:
            raise ValueError("invalid combination")


        """
        # get the variants related to the predictors
        df_groups_pheno = get_df_groups_only_subset(groups_file, "%s/df_groups_all.py"%p_dir_pheno, set(all_selected_predictors))
        sorted_vars_pheno = sorted(set(df_groups_pheno.variantID_across_samples))        
        df_vars_IDmapping = pd.DataFrame({"variantID_across_samples":sorted_vars_pheno}).reset_index().rename(columns={"index":"var_idx"})
        df_vars_IDmapping["var_name"] = df_vars_IDmapping.variantID_across_samples.apply(get_simple_var_name_variantID_across_samples)
        vars_IDmapping_file = "%s/varsID_mapping.tab"%(p_dir_pheno)
        save_df_as_tab(df_vars_IDmapping, vars_IDmapping_file)

        # map each sample to the test splitID, and from this to color
        sID_to_testSplit = pd.concat(map(pd.Series, df_ML_results_best_pheno.apply(lambda r: {s : r.test_split_ID for s in r.test_samples}, axis=1).values))
        if len(sID_to_testSplit)!=len(subset_samples): raise ValueError("invalid")
        testSplit_to_color = get_value_to_color(list(range(5)), type_color="hex", palette="tab10")[0]

        """

        # filter df
        df_ML_results = df_ML_results_best_pheno.copy()

        # one plot for each
        for Iml, ML_results in df_ML_results.copy().iterrows():

            #print(list(ML_results.keys()))

            #if phenotype!="VRC_resistance": continue
            #if ML_results.test_split_ID!=4: continue


            # change some
            ML_results["selected_predictors"] = ML_results.selected_predictors.split(",")

            # define the expected genes
            ML_results["selected_predictors_demult"] = get_selected_predictors_demult_one_tuple_predictors(ML_results["selected_predictors"])
            all_affected_gene_IDs = sorted(set(map(get_affected_geneID_one_unique_groupID, ML_results["selected_predictors_demult"])))
            all_affected_genes = list(map(lambda g: geneID_to_name[g], all_affected_gene_IDs))
            target_geneIDs_str = ",".join(sorted(all_affected_gene_IDs))

            print(phenotype, ML_results.test_split_ID, ML_results.cross_splits_ROC_AUC, all_affected_genes, ML_results.validation_transitions_TrainTest_ROC_AUC, ML_results.model_name)
        
            # define outdir
            def convert_nan_to_0(x):
                if pd.isna(x): return 0.0
                else: return x

            p_dir_model = "%s/model_%i_AUC_%.3f_AUCval_%.3f_genes_%s_%s"%(p_dir_pheno, ML_results.test_split_ID, ML_results.cross_splits_ROC_AUC, convert_nan_to_0(ML_results.validation_transitions_TrainTest_ROC_AUC), "_".join(all_affected_genes), ML_results.model_name.split("|")[0]); make_folder(p_dir_model)


            # define the target groups
            set_target_groups = set.union(*[set(g.split(" * ")) for g in ML_results.selected_predictors])

            # check that interaction makes sense
            for gID in set_target_groups:
                if gID in ML_results.selected_predictors and any([gID in p.split(" * ") for p in ML_results.selected_predictors if " * " in p]):
                    raise ValueError("group cannot be in both")

            # define a df with the groups, and their info
            df_group_info = get_df_group_info_plotting_ML_models(set_target_groups, gene_name_to_ID)

            gID_to_interacting_gID = {g.split(" * ")[0]:g.split(" * ")[1] for g in ML_results.selected_predictors if " * " in g}
            for gID in set_target_groups.difference(set(gID_to_interacting_gID)):
                gID_to_interacting_gID[gID] = "-"

            df_group_info["groupID_interaction"] = df_group_info.unique_groupID.map(gID_to_interacting_gID)
            check_no_nans_in_df(df_group_info)
            group_info_file = "%s/group_info.tab"%(p_dir_model)
            save_df_as_tab(df_group_info[["unique_groupID", "group_color", "groupID_interaction"]], group_info_file)

            # get the sample ID mapping
            all_samples_sorted = Tree(treefile_plot).get_leaf_names()
            df_samples_IDmapping = pd.DataFrame({"sampleID":all_samples_sorted})
            df_samples_IDmapping["sample_name"] = df_samples_IDmapping.sampleID.apply(lambda x: x.split("_")[1][1:]) # get sample ID
            sname_to_sID = df_samples_IDmapping.set_index("sample_name").sampleID
            df_samples_IDmapping["sample_name"] = "L"+df_samples_IDmapping["sample_name"]
            df_samples_IDmapping = df_samples_IDmapping.merge(df_clade_info, on="sampleID", how="left", validate="one_to_one")
            df_samples_IDmapping["clade"] = df_samples_IDmapping["clade"].apply(get_str_cladeID_Bergin2022)

            # add info for transitions
            if ML_results.type_phenotype=="phenotype_transition":

                # add transitions
                df_samples_IDmapping["sample_name_w_transitions"] = df_samples_IDmapping.sampleID.apply(get_sample_name_w_transition_samples, extra_info=ML_results.cross_splits_extra_info, sorted_samples_phylogenetic=all_samples_sorted)

                def get_sampleIDs_transitions_str(x):
                    if "(" in x: return "|".join([sname_to_sID[s] for s in x.split("(vs ")[1].split(")")[0].split(", ")])
                    else: return "|"

                df_samples_IDmapping["sampleIDs_transitions"] = df_samples_IDmapping.sample_name_w_transitions.apply(get_sampleIDs_transitions_str)

                # add visuals
                sampleID_to_color = get_value_to_color(all_samples_sorted, n=len(all_samples_sorted), palette="tab20", type_color="hex")[0]
                df_samples_IDmapping["sample_label_color"] = df_samples_IDmapping.sampleID.map(sampleID_to_color)
                all_linestyles = ['-', '--', '-.', ':']*100
                df_samples_IDmapping["sample_arrows_linestyle"] = all_linestyles[0:len(df_samples_IDmapping)]

            sampleID_mapping_file = "%s/samples.tab"%p_dir_model
            save_df_as_tab(df_samples_IDmapping, sampleID_mapping_file)

            # get the df with real and predicted phenotypes
            df_ML_metrics_model = df_ML_metrics[(df_ML_metrics.test_split_ID==ML_results.test_split_ID)].copy()
            sampleID_to_predicted_phenotype, metrics_dict = get_sample_to_predicted_phenotype_from_metrics_df(df_ML_metrics_model, "optimal_accuracy_custom")
            df_phenotypes = df_metadata_binary_phenotypes[["runID", phenotype]].rename(columns={"runID":"sampleID"}).copy()
            df_phenotypes = df_phenotypes[df_phenotypes.sampleID.isin(all_samples_sorted)].copy().reset_index(drop=True)
            df_phenotypes = df_phenotypes.rename(columns={pheno : get_paper_dicotomic_pheno(pheno).replace(" ", "_") for pheno in [phenotype]})

            #df_phenotypes["predicted_" + get_paper_dicotomic_pheno(phenotype).replace(" ", "_")] = df_phenotypes.sampleID.map(sampleID_to_predicted_phenotype)
            sampleID_to_phenoProb = dict(zip(ML_results.cross_splits_samples.split(","), map(float, ML_results.cross_splits_probabilities.split(","))))
            df_phenotypes["probability_" + get_paper_dicotomic_pheno(phenotype).replace(" ", "_")] = df_phenotypes.sampleID.apply(lambda x: sampleID_to_phenoProb[x])

            check_no_nans_in_df(df_phenotypes)
            phenotypes_file = "%s/phenos.tab"%p_dir_model
            save_df_as_tab(df_phenotypes, phenotypes_file)

            # get a phenotype info file
            df_pheno_mapping = pd.DataFrame({pheno : {"phenotype":pheno, "short_phenotype":"-"} for pheno in list(df_phenotypes.keys()) if pheno!="sampleID"}).transpose().reset_index(drop=True)
            df_pheno_mapping["use_for_correlations"] = df_pheno_mapping.phenotype==(get_paper_dicotomic_pheno(phenotype).replace(" ", "_"))
            phenotype_info_file = p_dir_model + "/pheno_info.tab"
            save_df_as_tab(df_pheno_mapping, phenotype_info_file)

            # define a title
            title_str = get_paper_dicotomic_pheno(phenotype) + "; m%s; "%(ML_results.test_split_ID) + "; ".join(["%s=%.2f"%(m, metrics_dict[m]) for m in ["recall", "specificity", "precision", "ROC_AUC"]]) + "; %s"%ML_results.model_name + "; %s"%ML_results.type_phenotype

            print(title_str)

            # define args
            extra_args = ""
            
            # plot predictors results
            outdir_p = "%s/plot_results"%(p_dir_model)
            #delete_folder(outdir_p)

            if file_is_empty(outdir_p + "/pipeline_finished.txt"):
                print("Runing cmd...")
                
                #cmd = "%s/scripts/ancestral_GWAS/scripts/plot_tree_variants.py --outdir %s --threads 112 --target_groups \"%s\" --min_support 70 --ASR_method DOWNPASS --variants_file %s --variant_annotations_file %s --groups_file %s --ASR_phenotypes %s --treefile %s --gff %s --ref_genome %s --pathway_annotations_dir %s --CGD_chromosomal_features %s --classification_method_ASR all_ignoreNaN --vars_IDmapping %s --custom_geneID_maps Scer_HAL9:TAC1 --interproscan_output %s --sampleID_mapping %s --only_vars_correlated_phenotype --subset_samples %s"%(ParentDir, outdir_p, ML_results.selected_predictors, variants_file, variant_annotations_file, groups_file, phenotypes_ASR, treefile, gff_Cparapsilosis, genome_Cparapsilosis, pathway_annotations_dir, CGD_chromosomal_features, vars_IDmapping_file, interproscan_output, sampleID_mapping_file, subset_samples_file) # optimal_F1, optimal_NegPredVal --representative_samples %s 

                cmd = "%s/scripts/ancestral_GWAS/scripts/plot_tree_variants_multiple_phenos.py --outdir %s --variants_file %s --variant_annotations_file %s --phenotypes_file %s --treefile %s --group_info %s --min_support 70 --ASR_method DOWNPASS --classification_method_ASR all_ignoreNaN --gff %s --ref_genome %s --pathway_annotations_dir %s --CGD_chromosomal_features %s --interproscan_output %s --sampleID_mapping %s --custom_geneID_maps %s --threads 112 --mtDNA_chrom %s --gDNA_code %i --mtDNA_code %i   --phenotype_info %s --title_str \"%s\" %s --target_geneIDs %s"%(ParentDir, outdir_p, variants_file, variant_annotations_file, phenotypes_file, treefile_plot, group_info_file, gff_Cparapsilosis, genome_Cparapsilosis, pathway_annotations_dir, CGD_chromosomal_features, interproscan_output, sampleID_mapping_file, mapping_genes, mtDNA_chrom, Cpara_gDNA_code, Cpara_mtDNA_code, phenotype_info_file, title_str, extra_args, target_geneIDs_str) # title_str, --only_vars_correlated_phenotype

                run_cmd(cmd, env="ancestral_GWAS_env")   # --only_vars_correlated_phenotype


def print_n_transitions_redundant(df_ML_results_all):

    print("# transitions redundant:")
    df_n_transitions_redundant = get_filtered_df_by_dict(df_ML_results_all, {"sample_set":"representative_samples", "test_sample_set":"representative_samples"})
    df_n_transitions_redundant = df_n_transitions_redundant[df_n_transitions_redundant.selected_predictors.apply(len)>0].copy()
    print(df_n_transitions_redundant[["phenotype", "aligner", "redundant_samples_n_transitions"]].drop_duplicates().set_index("phenotype"))


def plot_max_AUC_different_test_sample_sets(df_ML_results_maxAUC, PlotsDir, threshold_ROC_AUC):

    """Plots, for each phenotype and test_sample_set, the max AUC"""

    # get plot
    df_ML_results_maxAUC = df_ML_results_maxAUC.copy().rename(columns={"sample_set":"train_sample_set"})

    def get_max_AUC_each_pheno(df):
        fdict = dict(df[["phenotype", "max_AUC_field"]].iloc[0])
        df_square = df.pivot(values="max_ROC_AUC", index="test_sample_set", columns="train_sample_set")
        test_ss_to_maxAUC = dict(df_square.apply(max, axis=1))
        test_ss_to_maxAUC = {"max_AUC_"+k : v for k,v in test_ss_to_maxAUC.items()}
        return pd.Series({**fdict, **test_ss_to_maxAUC})

    df_plot_all = df_ML_results_maxAUC.groupby(["phenotype", "max_AUC_field"]).apply(get_max_AUC_each_pheno).reset_index(drop=True)

    # add fields
    df_plot_all["AUC"] = df_plot_all.max_AUC_field.map({"ROC_AUC_all":"all", "ROC_AUC_test_integrated":"int", "cross_splits_ROC_AUC":"cross_splits", "ROC_AUC_test_integrated_onlyCorrectModels":"int_corr"})
    check_no_nans_in_df(df_plot_all)

    # filter some
    df_plot_all = df_plot_all[df_plot_all.AUC.isin({"cross_splits", "int", "int_corr"})]

    # plot
    g = sns.relplot(data=df_plot_all, kind="scatter", x="max_AUC_representative_samples", y="max_AUC_all_samples", col="AUC", hue="phenotype", palette="tab20", aspect=1, height=3)
    for ax in g.axes.flat:

        lims = [0.45, 1.05]
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.plot(lims, lims, c="k", lw=.7, ls="--")

        for val in [threshold_ROC_AUC, 0.9]:
            ax.axvline(val, c="gray", lw=.7, ls="--")
            ax.axhline(val, c="gray", lw=.7, ls="--")

    plt.show()  
    save_fig(g, "max_AUC_ML_classifiers_test_ss.pdf", PlotsDir)

def plot_max_AUC_different_train_test_sample_sets(df_ML_results_maxAUC, PlotsDir, threshold_ROC_AUC):

    """Plots the max AUC"""

    # parse df to make long
    df_ML_results_maxAUC = df_ML_results_maxAUC.copy().rename(columns={"sample_set":"train_sample_set"})

    # keep some
    df_ML_results_maxAUC = df_ML_results_maxAUC[df_ML_results_maxAUC.max_AUC_field.isin({"ROC_AUC_test_integrated_onlyCorrectModels"})].copy()

    # get plot
    df_plot_all = pd.DataFrame()
    all_sample_sets = ["all_samples", "representative_samples"]

    for (train_sample_set, max_AUC_field),  df_long  in df_ML_results_maxAUC.groupby(["train_sample_set", "max_AUC_field"]):

        # keep df
        df_plot = df_long[["test_sample_set", "phenotype", "max_ROC_AUC"]].pivot(values="max_ROC_AUC", index="phenotype", columns="test_sample_set").rename(columns={ss : "max_AUC_"+ss for ss in all_sample_sets}).reset_index()

        df_plot["train_sample_set"] = train_sample_set
        df_plot["max_AUC_field"] = max_AUC_field
        df_plot_all = pd.concat([df_plot_all, df_plot])

    # plot
    df_plot_all = df_plot_all.rename(columns={"train_sample_set":"train", "max_AUC_field":"AUC"})
    check_no_nans_in_df(df_plot_all)

    df_plot_all["train"] = df_plot_all.train.apply(lambda x: x[0:3])
    df_plot_all["AUC"] = df_plot_all.AUC.map({"ROC_AUC_all":"all", "ROC_AUC_test_integrated":"int", "cross_splits_ROC_AUC":"cross_splits", "ROC_AUC_test_integrated_onlyCorrectModels":"int_corr"})

    g = sns.relplot(data=df_plot_all, kind="scatter", x="max_AUC_representative_samples", y="max_AUC_all_samples", col="train", row="AUC", hue="phenotype", palette="tab20", aspect=1, height=3)

    for ax in g.axes.flat:

        lims = [0.45, 1.05]
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.plot(lims, lims, c="k", lw=.7, ls="--")

        for val in [threshold_ROC_AUC, 0.9]:
            ax.axvline(val, c="gray", lw=.7, ls="--")
            ax.axhline(val, c="gray", lw=.7, ls="--")

    plt.show()  
    save_fig(g, "max_AUC_ML_classifiers.pdf", PlotsDir)


from sklearn.metrics import roc_auc_score


def plot_max_test_AUC_one_parameter(df_ML_integrated_performance, filters_dict, PlotsDir, target_field, baseline_val, threshold_ROC_AUC):

    """Evaluate how does the max AUC compare of different parameters"""

    # filter df
    if len(filters_dict)==0:
        df_ML_performance = df_ML_integrated_performance.copy()

    else:
        df_ML_performance = get_filtered_df_by_dict(df_ML_integrated_performance, filters_dict, check_all_present=True)

    # keep some fields
    df_ML_performance = df_ML_performance[["AUC_test", "test_parmID", "phenotype", target_field]].copy()
    baseline_val = str(baseline_val)
    df_ML_performance[target_field] = df_ML_performance[target_field].apply(str)

    if len(df_ML_performance)!=len(df_ML_performance[["test_parmID", "phenotype"]].drop_duplicates()):
        raise ValueError("not unique")

    # define the baseline max AUC
    def get_maxAUC_df_per_pheno(val): 
        series_val = df_ML_performance[df_ML_performance[target_field]==val].groupby("phenotype").apply(lambda df: max(df.AUC_test))
        return pd.DataFrame({"max_AUC":series_val}).reset_index()


    df_baseline = get_maxAUC_df_per_pheno(baseline_val).rename(columns={"max_AUC":"max_AUC_%s"%(baseline_val)})

    # map each phenotype to the max AUC
    pheno_to_max_AUC = df_ML_performance.groupby("phenotype").apply(lambda df: max(df.AUC_test))

    # create plot df
    df_plot_all = pd.DataFrame()
    for other_val in sorted(set(df_ML_performance[target_field]).difference({baseline_val})):

        df_plot = get_maxAUC_df_per_pheno(other_val)
        df_plot = df_plot.merge(df_baseline, on="phenotype", validate="one_to_one")
        df_plot[target_field] = other_val

        df_plot["is_max_AUC"] = df_plot.apply(lambda r: pheno_to_max_AUC[r.phenotype]==r.max_AUC, axis=1)
        df_plot["is_max_AUC_baseline"] = df_plot.apply(lambda r: pheno_to_max_AUC[r.phenotype]==r["max_AUC_%s"%(baseline_val)], axis=1)


        df_plot_all = pd.concat([df_plot_all, df_plot])
 
        #r["max_AUC_%s"%(baseline_val)]

    # plot
    pheno_to_color = get_value_to_color(sorted(set(df_plot_all.phenotype)), palette="tab20")[0]
    g = sns.relplot(data=df_plot_all, kind="scatter", x="max_AUC_%s"%(baseline_val), y="max_AUC", col=target_field, hue="phenotype", palette=pheno_to_color, aspect=1, height=4, edgecolor="black", linewidth=.7)

    for ax in g.axes.flat:

        lims = [0.45, 1.05]
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.plot(lims, lims, c="k", lw=.7, ls="--")

        for val in [threshold_ROC_AUC, 0.9]:
            ax.axvline(val, c="gray", lw=.7, ls="--")
            ax.axhline(val, c="gray", lw=.7, ls="--")

        # add the maximum AUC as arrows
        target_field_val = ax.get_title().split()[-1]
        df_plot_maxAUC = df_plot_all[(df_plot_all[target_field]==target_field_val) & (df_plot_all.is_max_AUC==True)].copy()

        for I,r in df_plot_maxAUC.iterrows():
            ax.arrow(r["max_AUC_%s"%(baseline_val)]-0.05, r["max_AUC"], 0.015, 0, width=0.005, shape="full", color=pheno_to_color[r.phenotype])

        # add the max AUC of the baseline
        df_plot_maxAUC_base = df_plot_all[(df_plot_all[target_field]==target_field_val) & (df_plot_all.is_max_AUC_baseline==True)].copy()

        for I,r in df_plot_maxAUC_base.iterrows():
            ax.arrow(r["max_AUC_%s"%(baseline_val)], r["max_AUC"]-0.05, 0, 0.015, width=0.005, shape="full", color=pheno_to_color[r.phenotype])


    plt.show()  

    filename_plot = "%s_vs_%s_%s.pdf"%(target_field, baseline_val, "-".join(["%s=%s"%(k, filters_dict[k]) for k in sorted(filters_dict)]))
    save_fig(g, "ML_results_changing_one_parameter/"+filename_plot, PlotsDir)


def get_dict_accuracy_fields_mapping():

    """gets a dict mapping pestimate, accuracy_field to the field in the df"""

    accuracy_dict = {
                    "test": {"ROC AUC":"ROC_AUC_test_integrated_onlyCorrectModels", 
                             "recall":"test_integrated_onlyCorrectModels_optimal_recall", 
                             "fpr":"test_integrated_onlyCorrectModels_optimal_fpr", 
                             "precision":"test_integrated_onlyCorrectModels_optimal_precision", 
                             "F-value":"test_integrated_onlyCorrectModels_optimal_F_score_precision_recall"},

                    "cross-val": {"ROC AUC":"cross_splits_ROC_AUC", 
                             "recall":"cross_splits_optimal_recall", 
                             "fpr":"cross_splits_optimal_fpr", 
                             "precision":"cross_splits_optimal_precision", 
                             "F-value":"cross_splits_optimal_F_score_precision_recall"},

                    "discarded": {"ROC AUC":"redundant_samples_ROC_AUC", 
                             "recall":"redundant_samples_optimal_recall", 
                             "fpr":"redundant_samples_optimal_fpr", 
                             "precision":"redundant_samples_optimal_precision", 
                             "F-value":"redundant_samples_optimal_F_score_precision_recall"},

                    "validation": {"ROC AUC":"validation_transitions_TrainTest_ROC_AUC", 
                             "recall":"validation_transitions_TrainTest_optimalThreshold_recall", 
                             "fpr":"validation_transitions_TrainTest_optimalThreshold_fpr", 
                             "precision":"validation_transitions_TrainTest_optimalThreshold_precision", 
                             "F-value":"validation_transitions_TrainTest_optimalThreshold_F_score_precision_recall"}
                    }

    return accuracy_dict

    # [60 rows x 115 columns] ['ASR_method', 'AUC_diff_threshold', 'CVset_model_consistency', 'CVset_type_split', 'aligner', 'all_affected_genes_str', 'all_affected_vars_str', 'class_weight', 'consistency_measure', 'expected_genes_str', 'feature_selection_method', 'features_max_pval_GWAS', 'gene_set', 'maxT_pval_consistency_mean_ROC_AUC', 'maxT_pval_consistency_min_ROC_AUC', 'mean_ROC_AUC_CVset_resamples', 'min_ROC_AUC_CVset_resamples', 'min_support', 'model_name', 'multiple_test_corr', 'prioritize_expected_genes', 'real_type_phenotype', 'sample_set', 'selected_predictors', 'selected_predictors_redundant', 'test_split_ID', 'time_spent_seconds_CVset', 'transition_definition', 'type_phenotype', 'test_sample_set', 'train_samples', 'test_samples', 'train_phenos', 'test_phenos', 'train_probabilities', 'test_probabilities', 'ignored_samples_train', 'ignored_samples_test', 'extra_info_train', 'extra_info_test', 'ROC_AUC_train', 'ROC_AUC_test', 'ROC_AUC_all', 'parameters_ID', 'cross_splits_phenotypes', 'cross_splits_samples', 'cross_splits_probabilities', 'cross_splits_ignored_samples', 'cross_splits_extra_info', 'cross_splits_ROC_AUC', 'redundant_samples_phenotypes', 'redundant_samples_samples', 'redundant_samples_probabilities', 'redundant_samples_ignored_samples', 'redundant_samples_ROC_AUC', 'redundant_samples_n_transitions', 'cross_splits_consistency_jaccard_variants', 'cross_splits_consistency_jaccard_samples', 'cross_splits_consistency_equal_variants', 'cross_splits_consistency_equal_samples', 'ROC_AUC_test_integrated_onlyCorrectModels', 'test_integrated_onlyCorrectModels_phenotypes', 'test_integrated_onlyCorrectModels_probabilities', 'test_integrated_onlyCorrectModels_samples', 'test_integrated_onlyCorrectModels_ignored_samples', 'ROC_AUC_test_integrated', 'test_integrated_phenotypes', 'test_integrated_probabilities', 'test_integrated_samples', 'test_integrated_ignored_samples', 'nmodels_fit', 'PRC_AUC_all', 'PRC_AUC_cross_splits', 'PRC_AUC_redundant_samples', 'PRC_AUC_test', 'PRC_AUC_test_integrated', 'PRC_AUC_test_integrated_onlyCorrectModels', 'PRC_AUC_train', 'test_parmID', 'n_train_parameter_combinations', 'cross_splits_consistency_jaccard_genes', 'test_parmID_exceptTestSS', 'phenotype', 'optimal_threshold_cross_splits', 'optimal_threshold_test_integrated_onlyCorrectModels', 'cross_splits_optimal_recall', 'cross_splits_optimal_precision', 'cross_splits_optimal_fpr', 'cross_splits_optimal_neg_predictive_val', 'cross_splits_optimal_F_score_precision_recall', 'redundant_samples_optimal_recall', 'redundant_samples_optimal_precision', 'redundant_samples_optimal_fpr', 'redundant_samples_optimal_neg_predictive_val', 'redundant_samples_optimal_F_score_precision_recall', 'test_integrated_onlyCorrectModels_optimal_recall', 'test_integrated_onlyCorrectModels_optimal_precision', 'test_integrated_onlyCorrectModels_optimal_fpr', 'test_integrated_onlyCorrectModels_optimal_neg_predictive_val', 'test_integrated_onlyCorrectModels_optimal_F_score_precision_recall', 'validation_transitions_TrainTest_ROC_AUC', 'validation_transitions_TrainTest_optimalThreshold_recall', 'validation_transitions_TrainTest_optimalThreshold_precision', 'validation_transitions_TrainTest_optimalThreshold_fpr', 'validation_transitions_TrainTest_optimalThreshold_neg_predictive_val', 'validation_transitions_TrainTest_optimalThreshold_F_score_precision_recall', 'validation_transitions_TrainTest_PRC_AUC', 'validation_transitions_Test_ROC_AUC', 'validation_transitions_Test_optimalThreshold_recall', 'validation_transitions_Test_optimalThreshold_precision', 'validation_transitions_Test_optimalThreshold_fpr', 'validation_transitions_Test_optimalThreshold_neg_predictive_val', 'validation_transitions_Test_optimalThreshold_F_score_precision_recall', 'validation_transitions_Test_PRC_AUC', 'best_model']



def plot_performance_distribution_top_models(ML_df_dicts, PlotsDir, ML_pheno_to_color, test_sample_set):

    """plots the performance of the top models"""

    # phenos
    sorted_phenotypes = ["FLC_resistance", "VRC_resistance", "age_infant", "source_blood_vs_other_sterile"]
    #sorted_phenotypes = ["FLC_resistance", "source_blood_vs_other_sterile"]

    # keep some
    df_ML_results_top = get_filtered_df_by_dict(ML_df_dicts["df_ML_results_top"].copy(), {"test_sample_set":test_sample_set, "phenotype":set(sorted_phenotypes)})
    df_ML_results_top = df_ML_results_top[(~pd.isna(df_ML_results_top.selected_predictors))].copy()

    # keep only best model
    #df_ML_results_top = df_ML_results_top[df_ML_results_top.best_model==True].copy()

    # define metrics
    performance_metrics = ["ROC AUC", "recall", "fpr", "precision", "F-value"]
    sorted_pestimates = {"representative_samples":["test", "cross-val", "validation"], "all_samples":["test", "cross-val", "validation"]}[test_sample_set] # discarded

    # get a df for plotting
    df_plot_all = pd.DataFrame()

    # init figure
    for phenotype in sorted_phenotypes:
        df_pheno = df_ML_results_top[df_ML_results_top.phenotype==phenotype].copy().sort_values(by=["test_parmID_exceptTestSS", "test_split_ID"])

        for pestimate in sorted_pestimates:

            # discard some
            if pestimate=="validation" and not phenotype in {"FLC_resistance", "VRC_resistance"}: continue

            # checks
            if len(df_pheno)!=len(df_pheno[["test_parmID_exceptTestSS", "test_split_ID"]].drop_duplicates()):
                raise ValueError("should be unique")

            # define the target df for this estimate
            if pestimate=="test":
                df_plot_estimate = pd.DataFrame(df_pheno.groupby("test_parmID_exceptTestSS").apply(lambda df: df.iloc[0]))

            elif pestimate in {"cross-val", "validation"}:
                df_plot_estimate = df_pheno.copy()

            elif pestimate=="discarded":
                df_plot_estimate = df_pheno[(df_pheno.sample_set=="representative_samples") & (df_pheno.redundant_samples_n_transitions>1)].copy()
                if len(df_plot_estimate)==0: continue

            else:
                raise ValueError("invalid pestimate")

            # checks
            if len(df_plot_estimate)==0: raise ValueError("invalid pestimate")

            # add consistency
            for performance_metric in performance_metrics:

                accuracy_f = get_dict_accuracy_fields_mapping()[pestimate][performance_metric]
                df_plot = df_plot_estimate[[accuracy_f, "phenotype", "test_parmID_exceptTestSS", "test_split_ID", "best_model"]].copy().rename(columns={accuracy_f:"performance"})
                df_plot["testing strategy"] = pestimate
                df_plot["performance metric"] = performance_metric

                # change if fpr
                if performance_metric=="fpr":
                    df_plot["performance"] = 1 - df_plot.performance
                    df_plot["performance metric"] = "specificity"


                # keep
                df_plot_all = pd.concat([df_plot_all, df_plot]).reset_index(drop=True)
                check_no_nans_in_df(df_plot_all)

    # checks
    if len(df_plot_all)!=len(df_plot_all[["test_parmID_exceptTestSS", "test_split_ID", "phenotype", "performance metric", "testing strategy"]].drop_duplicates()):
        raise ValueError("should be unique")
    
    # plot
    plot_performance_metrics = ["ROC AUC", "recall", "specificity", "precision"]
    #plot_performance_metrics = ["ROC AUC", "recall", "specificity"]
    df_plot_all["pestimate_idx"] = df_plot_all["testing strategy"].map(dict(zip(sorted_pestimates, range(len(sorted_pestimates)))))
    check_no_nans_series(df_plot_all["pestimate_idx"])

    width_mult = {4:1.5, 3:1.2}[len(sorted_pestimates)]
    ncols = len(plot_performance_metrics)
    fig = plt.figure(figsize=(ncols*width_mult, 1.5))

    for Ip, performance_metric in enumerate(plot_performance_metrics):

        # plot
        df_plot = df_plot_all[df_plot_all["performance metric"]==performance_metric].copy().sort_values(by=["pestimate_idx", "phenotype"])
        ax = plt.subplot(1, ncols, Ip+1)

        #ax = sns.barplot(data=df_plot, x="testing strategy", y="performance", hue="phenotype", palette=ML_pheno_to_color, linewidth=.7, edgecolor="white")
        

        ax = sns.stripplot(data=df_plot[df_plot.best_model==False], x="testing strategy", y="performance", hue="phenotype", palette=ML_pheno_to_color, s=3, linewidth=.3, edgecolor="black", dodge=True, jitter=0.2, alpha=0.1)
        
        ax = sns.stripplot(data=df_plot[df_plot.best_model==True], x="testing strategy", y="performance", hue="phenotype", palette=ML_pheno_to_color, s=3, linewidth=.3, edgecolor="black", dodge=True, jitter=0.2, alpha=1)

        #ax = sns.stripplot(data=df_plot, x="testing strategy", y="performance", hue="phenotype", palette=ML_pheno_to_color, s=3, linewidth=.3, edgecolor="black", dodge=True, jitter=0.2, alpha=0.1)

        plt.axhline(0.5, c="gray", lw=.7, ls="--")
        plt.axhline(0.7, c="gray", lw=.7, ls="--")


        #ax = sns.boxplot(data=df_plot, x="testing strategy", y="performance", hue="phenotype", palette=ML_pheno_to_color)
        #ax = sns.swarmplot(data=df_plot, x="testing strategy", y="performance", hue="phenotype", palette=ML_pheno_to_color, s=3, linewidth=0, edgecolor="black", dodge=True, alpha=1)

        #for pheno in sorted_phenotypes:
        #    #ax = sns.stripplot(data=df_plot[df_plot.phenotype==pheno].copy(), x="testing strategy", y="performance", s=3, linewidth=1, edgecolor=ML_pheno_to_color[pheno], dodge=True, jitter=True, alpha=1, facecolor='none')
        #    ax = sns.swarmplot(data=df_plot[df_plot.phenotype==pheno].copy(), x="testing strategy", y="performance", s=3, linewidth=1, edgecolor=ML_pheno_to_color[pheno], dodge=True, alpha=1, facecolor='none')


        for x in range(len(performance_metrics)-2):
            ax.axvline(x+0.5, ls="--", c="gray", lw=.7)

        ax.get_legend().remove()
        ax.set_ylim([-0.05, 1.05])
        ax.set_xticklabels(sorted_pestimates, rotation=90)
        ax.set_xlabel("")
        ax.set_title(performance_metric, fontsize=10)

        if Ip!=0:
            ax.set_yticklabels([])
            ax.set_yticks([])
            ax.set_ylabel("")



    # get 
    sns.despine()
    plt.subplots_adjust(wspace=0.1, hspace=0.1)# wspace (between cols) hspaces (btw rows)
    save_fig(fig, "performance_best_trhesholds_%s.pdf"%(test_sample_set), PlotsDir)
    plt.show()     




def plot_ROC_curves_best_parameters(ML_df_dicts, PlotsDir, threshold_ROC_AUC, ML_pheno_to_color, add_threshold_lines=False):

    """plot roc curves"""

    # get the top models and metrics
    df_ML_results_top = ML_df_dicts["df_ML_results_top"].copy()
    df_ML_results_top = df_ML_results_top[~pd.isna(df_ML_results_top.selected_predictors)].copy()
    df_ML_metrics_top = ML_df_dicts["df_ML_metrics_top"].copy()

    # get the performance df
    df_ML_integrated_performance_filt = get_filtered_df_by_dict(ML_df_dicts["df_ML_integrated_performance_filt"], {"top_model":True})

    # define phenos
    sorted_phenotypes = ["FLC_resistance", "VRC_resistance", "age_infant", "source_blood_vs_other_sterile"]

    # get different curves for different performance estimates
    for type_curve, xfield, yfield in [("ROC", "fpr", "recall")]:
        for test_sample_set, performance_estimates in [("representative_samples", ["test", "cross-val", "validation"]), 
                                                       ("all_samples", ["test", "cross-val", "validation"])]:

            # init figure
            ncols = len(performance_estimates)
            pestimate_to_col = dict(zip(performance_estimates, range(1, ncols+1)))
            fig = plt.figure(figsize=(ncols*1.5, 1.5)); Ic=1

            # one set of lines each pheno
            for phenotype in sorted_phenotypes:

                # get resuts to define the best 
                df_ML_results_pheno = get_filtered_df_by_dict(df_ML_results_top, {"test_sample_set":test_sample_set, "phenotype":phenotype})
                df_ML_metrics_pheno = get_filtered_df_by_dict(df_ML_metrics_top, {"test_sample_set":test_sample_set, "phenotype":phenotype})

                # get perfromance df
                df_ML_integrated_performance_pheno = get_filtered_df_by_dict(df_ML_integrated_performance_filt, {"phenotype":phenotype})[["test_parmID_exceptTestSS", "top_model_idx"]].copy()
                if len(df_ML_integrated_performance_pheno)!=len(set(df_ML_integrated_performance_pheno.top_model_idx)):
                    raise ValueError("top_model_idx should be unique")

                # one plot for each performance_estimate
                for Ip, pestimate in enumerate(performance_estimates):

                    # discard some
                    if pestimate=="validation" and not phenotype in {"FLC_resistance", "VRC_resistance"}: continue

                    # one sunplot for each pestimate
                    current_col = Ip+1
                    ax = plt.subplot(1, ncols, current_col)

                    # define the test_parmID_exceptTestSS 
                    if pestimate in {"test", "cross-val", "validation"}:
                        df_ML_results = df_ML_results_pheno[df_ML_results_pheno.best_model==True].copy()

                    elif pestimate=="discarded":

                        # get the best top model that has representative_samples
                        df_ML_results = get_filtered_df_by_dict(df_ML_results_pheno, {"sample_set":"representative_samples"})
                        df_ML_results = df_ML_results[df_ML_results.redundant_samples_n_transitions>1].copy()

                        if len(df_ML_results)==0: 
                            print("WARNING: for %s there are no top models with ss=rep_samples"%phenotype)
                            continue

                        df_ML_results = df_ML_results.merge(df_ML_integrated_performance_pheno, on="test_parmID_exceptTestSS", how="left", validate="many_to_one")
                        df_ML_results = df_ML_results[df_ML_results.top_model_idx==min(df_ML_results.top_model_idx)].copy()

                    else: raise ValueError("invalid pestimate")

                    if len(df_ML_results)==0: raise ValueError("invalid")
                    #if len(df_ML_results)==0: continue
                    if len(set(df_ML_results.test_split_ID))!=len(df_ML_results): raise ValueError("should be unique")
                    best_test_parmID_exceptTestSS = df_ML_results.test_parmID_exceptTestSS.iloc[0]

                    # get the metrics df
                    df_ML_metrics = get_filtered_df_by_dict(df_ML_metrics_pheno, {"test_parmID_exceptTestSS":best_test_parmID_exceptTestSS})

                    if pestimate=="test":
                        first_split_ID = df_ML_metrics.test_split_ID.iloc[0]
                        df_metrics = df_ML_metrics[(df_ML_metrics.eval_data=="test_integrated_onlyCorrectModels") & (df_ML_metrics.test_split_ID==first_split_ID)].copy()
                        field_opt_threshold = "optimal_threshold_test_integrated_onlyCorrectModels"

                    elif pestimate=="cross-val":
                        df_metrics = df_ML_metrics[(df_ML_metrics.eval_data=="cross_splits")].copy()
                        field_opt_threshold = "optimal_threshold_cross_splits"

                    elif pestimate=="discarded":
                        df_metrics = df_ML_metrics[(df_ML_metrics.eval_data=="redundant_samples")].copy()
                        field_opt_threshold = "optimal_threshold_cross_splits"

                    elif pestimate=="validation":
                        df_metrics = df_ML_metrics[(df_ML_metrics.eval_data=="transitions_TrainTest")].copy() # transitions_TrainTest, transitions_Test
                        field_opt_threshold = "optimal_threshold_cross_splits"

                    else: raise ValueError("invalud pestimate")

                    if len(df_metrics)!=len(df_metrics[["test_split_ID", "threshold"]]): raise ValueError("should be unique")

                    # get plot
                    def error_bar_fun(yvals):
                        median = np.median(yvals)
                        mad = scipy.stats.median_abs_deviation(yvals)
                        return (median-mad, median+mad)
             
                    ax = sns.lineplot(data=df_metrics, x=xfield, y=yfield, color=ML_pheno_to_color[phenotype], units="test_split_ID", estimator=None)
                    #ax = sns.lineplot(data=df_metrics, x=xfield, y=yfield, color=ML_pheno_to_color[phenotype], estimator="mean", errorbar="sd")

                    # add threshold lines
                    if add_threshold_lines is True:

                        df_metrics = df_metrics.merge(df_ML_results[["test_split_ID", field_opt_threshold]], how="left", validate="many_to_one")
                        df_metrics_optimal = df_metrics[df_metrics.threshold==df_metrics[field_opt_threshold]].copy()

                        for I,r in df_metrics_optimal.iterrows():
                            plt.plot([r[xfield]], [r[yfield]], c="none", marker="o", markeredgecolor=ML_pheno_to_color[phenotype], markeredgewidth=2)

                    lims = [-0.05, 1.05]
                    ticks = [0.0, 0.5, 1.0]
                    ax.set_xlim(lims)
                    ax.set_ylim(lims)

                    if type_curve=="ROC":
                        plt.plot(lims, lims, ls="--", lw=.7, c="gray")

                    elif type_curve=="PRC":
                        plt.axvline(0.7, ls="--", lw=.7, c="gray")
                        plt.axhline(0.7, ls="--", lw=.7, c="gray")

                    ax.set_xticks(ticks)
                    ax.set_yticks(ticks)

                    if current_col==1:
                        ax.set_ylabel(yfield)

                    else:
                        ax.set_ylabel("")
                        ax.set_yticklabels([])

                    ax.set_title(pestimate, fontsize=10)

                    # add legend
                    if len(performance_estimates)==current_col and phenotype=="FLC_resistance":
                        for Ip, phenotype in enumerate(sorted_phenotypes):
                            ax.text(1.1, 0.6-Ip*0.15, get_paper_dicotomic_pheno(phenotype), color=ML_pheno_to_color[phenotype], fontsize=10, ha="left")


            sns.despine()


            plt.subplots_adjust(wspace=0.1, hspace=0.1)# wspace (between cols) hspaces (btw rows)
            save_fig(fig, "ROC_PRC_curves/%s_%s.pdf"%(test_sample_set, type_curve), PlotsDir)
            plt.show()     

def plot_ROC_PRC_curves_all(ML_df_dicts, PlotsDir, threshold_ROC_AUC, ML_pheno_to_color, df_ML_results_all, df_ML_metrics_all):

    """gets the precision and recall plots"""

    # get dfs
    df_ML_integrated_performance_filt = ML_df_dicts["df_ML_integrated_performance_filt"].copy()

    # get the performance of the best parameter
    sorted_phenotypes = ["FLC_resistance", "VRC_resistance", "age_infant", "source_blood_vs_other_sterile"]
    df_ML_integrated_performance_filt = df_ML_integrated_performance_filt[(df_ML_integrated_performance_filt.phenotype.isin(sorted_phenotypes)) & (df_ML_integrated_performance_filt.best_test_parmID==True)].copy()

    fields_parms = ["CVset_type_split", "sample_set", "AUC_diff_threshold", "consistency_measure", "aligner", "multiple_test_corr"]
    print("phenotype; ", "; ".join(fields_parms))
    for I,r in df_ML_integrated_performance_filt.iterrows():
        print("; ".join(r[["phenotype"] + fields_parms].apply(str).values))

    # define coordinates
    pheno_to_coordinates = {"FLC_resistance": (.3, -.6), "VRC_resistance": (.7, -.6),
                            "agar_invasiveness": (.3, -0.75), "age_infant": (.7, -0.75)}

    # get figures
    for type_curve, xfield, yfield in [("ROC", "fpr", "recall"), ("PRC", "recall", "precision")]:
        for test_sample_set, performance_estimates in [("representative_samples", ["independent", "cross-val", "discarded"]), 
                                                       ("all_samples", ["independent", "cross-val"])]:

            # init figure
            ncols = len(performance_estimates)
            pestimate_to_col = dict(zip(performance_estimates, range(1, ncols+1)))
            fig = plt.figure(figsize=(ncols*1.5, 1.5))

            # go through each pheno
            for phenotype in sorted_phenotypes:
                print(phenotype)

                #if phenotype!="FLC_resistance": continue # debug!!!

                # filter dfs
                filters_dict = dict(df_ML_integrated_performance_filt[df_ML_integrated_performance_filt.phenotype==phenotype].iloc[0][fields_parms])
                filters_dict["phenotype"] = phenotype
                filters_dict["test_sample_set"] = test_sample_set

                df_ML_results = get_filtered_df_by_dict(df_ML_results_all, filters_dict)
                df_ML_metrics = get_filtered_df_by_dict(df_ML_metrics_all, filters_dict)

                # checks
                if len(df_ML_results)!=5: raise ValueError("df should be 5")

                # keep fit models
                df_ML_results = df_ML_results[df_ML_results.selected_predictors.apply(len)>0].copy()
                df_ML_metrics = df_ML_metrics[df_ML_metrics.test_split_ID.isin(set(df_ML_results.test_split_ID))].copy()

                # checks
                if len(df_ML_results)<4: raise ValueError("df should be >=4")

                # plot for each performance estimate
                for pestimate in performance_estimates:
                    print(pestimate)

                    # init
                    current_col = pestimate_to_col[pestimate]
                    ax = plt.subplot(1, ncols, current_col)

                    # define a dataset that has, for different model IDs, the thresholds and the precision, recall and AUCs
                    if pestimate=="independent":
                        first_split_ID = df_ML_metrics.test_split_ID.iloc[0]
                        df_metrics = df_ML_metrics[(df_ML_metrics.eval_data=="test_integrated_onlyCorrectModels") & (df_ML_metrics.test_split_ID==first_split_ID)].copy()

                    elif pestimate=="cross-val":
                        df_metrics = df_ML_metrics[(df_ML_metrics.eval_data=="cross_splits")].copy()

                    elif pestimate=="discarded":
                        df_metrics = df_ML_metrics[(df_ML_metrics.eval_data=="redundant_samples")].copy()

                    else:
                        raise ValueError("invalud pestimate")

                    # checks
                    if len(df_metrics)!=len(df_metrics[["threshold", "test_split_ID"]].drop_duplicates()):
                        raise ValueError("test split ID")

                    df_AUC = df_metrics[["test_split_ID", "ROC_AUC", "PRC_AUC"]].drop_duplicates()
                    if len(df_AUC)!=len(set(df_metrics.test_split_ID)):
                        raise ValueError("the AUCs should be per split")

                    # define function
                    def error_bar_fun(yvals):

                        median = np.median(yvals)
                        mad = scipy.stats.median_abs_deviation(yvals)
                        return (median-mad, median+mad)
             
                    # get plots
                    #sns.lineplot(data=df_metrics, x=xfield, y=yfield, color=ML_pheno_to_color[phenotype], units="test_split_ID", estimator=None)
                    #if len(set(df_metrics.test_split_ID))==1:
                    sns.lineplot(data=df_metrics, x=xfield, y=yfield, color=ML_pheno_to_color[phenotype], units="test_split_ID", estimator=None)

                    #else:
                    #    sns.lineplot(data=df_metrics, x=xfield, y=yfield, color=ML_pheno_to_color[phenotype], estimator=np.median, errorbar=error_bar_fun)

                    # add the AUC vals
                    AUC_val = "%.2f"%(np.mean(df_AUC["%s_AUC"%type_curve]))
                    if pestimate=="discarded" and "all_samples" in {filters_dict["sample_set"], filters_dict["test_sample_set"]}:
                        AUC_val = "-"

                    ax.text(pheno_to_coordinates[phenotype][0], pheno_to_coordinates[phenotype][1], AUC_val, color=ML_pheno_to_color[phenotype], ha="center", va="center")

                    # edit
                    lims = [-0.05, 1.05]
                    ticks = [0.0, 0.5, 1.0]
                    ax.set_xlim(lims)
                    ax.set_ylim(lims)

                    if type_curve=="ROC":
                        plt.plot(lims, lims, ls="--", lw=.7, c="gray")

                    elif type_curve=="PRC":
                        plt.axvline(0.7, ls="--", lw=.7, c="gray")
                        plt.axhline(0.7, ls="--", lw=.7, c="gray")

                    ax.set_xticks(ticks)
                    ax.set_yticks(ticks)

                    if current_col==1:
                        ax.set_ylabel(yfield)

                    else:
                        ax.set_ylabel("")
                        ax.set_yticklabels([])

                    if ncols==3:
                        if current_col==2:
                            ax.set_xlabel(xfield)

                        else:
                            ax.set_xlabel("")


                    elif ncols==2:
                        ax.set_xlabel("")

                        if current_col==1 and phenotype=="FLC_resistance":
                            ax.text(1.15, -0.4, xfield, ha="center", va="center")


                    else:
                        raise ValueError("invalid")

                    ax.set_title(pestimate, fontsize=10)

            sns.despine()
            plt.subplots_adjust(wspace=0.1, hspace=0.1)# wspace (between cols) hspaces (btw rows)

            save_fig(fig, "ROC_PRC_curves/%s_%s.pdf"%(test_sample_set, type_curve), PlotsDir)
            plt.show()



def old_plot_general_corr_train_test_split_by_pheno(df_ML_integrated_performance_filt, PlotsDir, threshold_ROC_AUC, df_transitions_info_GWAS, df_ML_stats, threshold_consistency_genes):

    """Gets the phenos for each"""

    # keep
    df_ML_integrated_performance_filt = df_ML_integrated_performance_filt.copy()
    df_ML_stats = df_ML_stats.copy()
    df_transitions_info_GWAS = df_transitions_info_GWAS.copy()

    # create consistency
    min_consistency = min(df_ML_integrated_performance_filt[df_ML_integrated_performance_filt.mean_cross_splits_consistency_jaccard_genes!=-1].mean_cross_splits_consistency_jaccard_genes)
    df_ML_integrated_performance_filt["mean_cross_splits_consistency_jaccard_genes"] = df_ML_integrated_performance_filt.mean_cross_splits_consistency_jaccard_genes.apply(lambda x: max([min_consistency, x]))
    max_consistency = max(df_ML_integrated_performance_filt.mean_cross_splits_consistency_jaccard_genes)

    # add fields
    AUC_train_field = "train $AUC_{ROC}$\n(mean m1-m5)"
    df_ML_integrated_performance_filt[AUC_train_field] = df_ML_integrated_performance_filt.AUC_train.apply(lambda x: max([0.5, x]))

    correct_AUC_train_f = "train $AUC_{ROC}$>=%s &\ntrain $AUC_{ROC}$ ~ test $AUC_{ROC}$"%threshold_ROC_AUC
    df_ML_integrated_performance_filt[correct_AUC_train_f] = (df_ML_integrated_performance_filt.AUC_train>=threshold_ROC_AUC) & (df_ML_integrated_performance_filt.AUC_train_test_diff.apply(abs)<=0.2)


    # init fig
    nrows = 2
    ncols = 7
    fig = plt.figure(figsize=(ncols*1.5, nrows*1.5))

    # one plot for each pheno
    sorted_phenotypes = list(df_ML_stats.sort_values(by=["fraction_correct_test_parameters", "correct_models_max_AUC_test"]).phenotype)
    for I, (phenotype) in enumerate(sorted_phenotypes):

        # get df
        df_ML_integrated_performance_p = df_ML_integrated_performance_filt[df_ML_integrated_performance_filt.phenotype==phenotype].copy()
        df_ML_integrated_performance_p = df_ML_integrated_performance_p.sort_values(by=["correct_model", correct_AUC_train_f], ascending=True)
        #df_ML_integrated_performance_p = df_ML_integrated_performance_p[df_ML_integrated_performance_p.correct_model].copy()



        # get scatter
        ax = plt.subplot(nrows, ncols, I+1)
        #ax = sns.scatterplot(data=df_ML_integrated_performance_p, x="mean_cross_splits_consistency_jaccard_genes", y="AUC_test", edgecolor="black", linewidth=.3, hue=AUC_train_field, s=20, palette="rocket_r", style="correct_model", markers={True:"X", False:"o"}, size="correct_model", size_order=[True, False], sizes=[50, 10]) # hue="phenotype", palette=phenotype_to_color # s

        ax = sns.scatterplot(data=df_ML_integrated_performance_p, x="mean_cross_splits_consistency_jaccard_genes", y="AUC_test", edgecolor="black", linewidth=.3, hue=correct_AUC_train_f, s=20, palette={True:"red", False:"gray"}, style="correct_model", markers={True:"X", False:"o"}, size="correct_model", size_order=[True, False], sizes=[35, 10]) # hue="phenotype", palette=phenotype_to_color # s

        # add lims
        xlims = [min_consistency-0.05, max_consistency+0.05]
        ylims = [0.45, 1.05]
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        #ax.plot(lims, lims, c="k", lw=.7, ls="--", zorder=0)

        # add box with top model
        df_top_models = df_ML_integrated_performance_p[df_ML_integrated_performance_p.top_model].copy()
        if len(df_top_models)>0:

            min_consistency_top = min(df_top_models.mean_cross_splits_consistency_jaccard_genes)
            max_consistency_top = max(df_top_models.mean_cross_splits_consistency_jaccard_genes)
            min_AUC_top = min(df_top_models.AUC_test)
            max_AUC_top = max(df_top_models.AUC_test)

            plt.plot([min_consistency_top, xlims[1]], [min_AUC_top, min_AUC_top], ls="--", c="black", lw=.7)
            plt.plot([min_consistency_top, min_consistency_top], [min_AUC_top, ylims[1]], ls="--", c="black", lw=.7)

        # add number of top models
        plt.text(0.95, 0.97, "%i top"%(len(df_top_models)), ha="right", fontsize=10)


        # add thresholds    
        ax.axvline(threshold_consistency_genes, c="gray", lw=.7, ls="--", zorder=0)
        ax.axhline(threshold_ROC_AUC, c="gray", lw=.7, ls="--", zorder=0)
        
        n_transitions = df_transitions_info_GWAS[df_transitions_info_GWAS.phenotype==phenotype].n_transitions.iloc[0]
        formatted_text = get_paper_dicotomic_pheno(phenotype).replace("agar ", "").replace("blood vs other sterile", "blood-sterile").replace("bronchial aspirate", "bronch. aspir.").replace("voriconazole", "vori.").replace("fluconazole", "fluc.").replace("micafungin", "mica.")

        # add title
        #ax.set_title(formatted_text+"\n#t="+str(n_transitions), fontsize=10)
        ax.set_title(formatted_text + " (%i)"%sum(df_ML_integrated_performance_p.correct_model), fontsize=10)

        # add boxes
        #box = mpatches.Rectangle((0.43, 0.95), 1, 0.15, facecolor='white', edgecolor='white', linewidth=0.5)
        #ax.add_patch(box)
        #ax.text(0.47, 0.98, "#x=%i"%(sum(df_ML_integrated_performance_p.correct_model)), fontsize=10, color="black", ha='left', va='center')
        #ax.text(0.47, 0.98, "#t=%i"%n_transitions, fontsize=10, color="black", ha='left', va='center')

        # get ticcks
        yticks = [0.5, 0.7, 1.0]
        xticks = [0.4, 0.8]
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        # set lims
        if I in {7}:
            ax.set_ylabel("          test $AUC_{ROC}$", ha="left")
            ax.set_yticklabels(list(map(str, yticks)))

        elif I==0:
            ax.set_ylabel("")

        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        if I in {10}:
            ax.set_xlabel("m1-m5 feature consistency (genes)")

        else:
            ax.set_xlabel("")

        if I in set(range(7)):
            ax.set_xticklabels([])

        else:
            ax.set_xticklabels(list(map(str, xticks)))

        if I==12:
            ax.legend(bbox_to_anchor=(1.1,-0.5), loc="lower left", markerscale=2)

        else:
            if len(df_ML_integrated_performance_p)>0:
                ax.get_legend().remove()




    sns.despine()
    plt.subplots_adjust(wspace=0.1, hspace=0.35)# wspace (between cols) hspaces (btw rows)
    save_fig(fig, "overfitting_check_filt_parms_split.pdf", PlotsDir)



def plot_general_corr_train_test_split_by_pheno(df_ML_integrated_performance_filt, PlotsDir, threshold_ROC_AUC, df_transitions_info_GWAS, df_ML_stats, threshold_consistency_genes, test_sample_set, add_colors=True):

    """Gets the phenos for each"""

    # keep
    df_ML_integrated_performance_filt = df_ML_integrated_performance_filt.copy()
    df_ML_stats = df_ML_stats.copy()
    df_transitions_info_GWAS = df_transitions_info_GWAS.copy()

    # create consistency
    min_consistency = min(df_ML_integrated_performance_filt[df_ML_integrated_performance_filt.mean_cross_splits_consistency_jaccard_genes!=-1].mean_cross_splits_consistency_jaccard_genes)
    df_ML_integrated_performance_filt["mean_cross_splits_consistency_jaccard_genes"] = df_ML_integrated_performance_filt.mean_cross_splits_consistency_jaccard_genes.apply(lambda x: max([min_consistency, x]))
    max_consistency = max(df_ML_integrated_performance_filt.mean_cross_splits_consistency_jaccard_genes)

    # add fields
    AUC_train_field = "train $AUC_{ROC}$\n(mean m1-m5)"
    df_ML_integrated_performance_filt[AUC_train_field] = df_ML_integrated_performance_filt.AUC_train.apply(lambda x: max([0.5, x]))

    correct_AUC_train_f = "train $AUC_{ROC}$>=%s &\ntrain $AUC_{ROC}$ ~ test $AUC_{ROC}$"%threshold_ROC_AUC
    df_ML_integrated_performance_filt[correct_AUC_train_f] = (df_ML_integrated_performance_filt.AUC_train>=threshold_ROC_AUC) & (df_ML_integrated_performance_filt.AUC_train_test_diff.apply(abs)<=0.2)

    # filter to keep correct train f
    df_ML_integrated_performance_filt = df_ML_integrated_performance_filt[(df_ML_integrated_performance_filt[correct_AUC_train_f]==True)].copy()

    # add the model type
    def get_type_model_correct_or_top(r):
        if r.top_model==True: return "top"
        elif r.correct_model==True: return "correct"
        else: return "low"

    df_ML_integrated_performance_filt["type performance"] = df_ML_integrated_performance_filt[["correct_model", "top_model"]].apply(get_type_model_correct_or_top, axis=1)

    if add_colors is True:
        performanceType_to_color = {"top":"maroon", "correct":"salmon", "low":"gray"}

    else:
        performanceType_to_color = {"top":"gray", "correct":"gray", "low":"gray"}

    # init fig
    nrows = 2
    ncols = 7
    fig = plt.figure(figsize=(ncols*1.5, nrows*1.5))

    # one plot for each pheno
    sorted_phenotypes = list(df_ML_stats.sort_values(by=["fraction_correct_test_parameters", "correct_models_max_AUC_test"]).phenotype)
    for I, (phenotype) in enumerate(sorted_phenotypes):

        # get df
        df_ML_integrated_performance_p = df_ML_integrated_performance_filt[df_ML_integrated_performance_filt.phenotype==phenotype].copy()
        #df_ML_integrated_performance_p = df_ML_integrated_performance_p.sort_values(by=["correct_model", correct_AUC_train_f], ascending=True)
        #df_ML_integrated_performance_p = df_ML_integrated_performance_p[df_ML_integrated_performance_p.correct_model].copy()
        df_ML_integrated_performance_p = df_ML_integrated_performance_p.sort_values(by=["top_model", "correct_model"], ascending=True)

        # get scatter
        ax = plt.subplot(nrows, ncols, I+1)
        #ax = sns.scatterplot(data=df_ML_integrated_performance_p, x="mean_cross_splits_consistency_jaccard_genes", y="AUC_test", edgecolor="black", linewidth=.3, hue=AUC_train_field, s=20, palette="rocket_r", style="correct_model", markers={True:"X", False:"o"}, size="correct_model", size_order=[True, False], sizes=[50, 10]) # hue="phenotype", palette=phenotype_to_color # s

        ax = sns.scatterplot(data=df_ML_integrated_performance_p, x="mean_cross_splits_consistency_jaccard_genes", y="AUC_test", edgecolor="black", linewidth=.3, hue="type performance", s=20, palette=performanceType_to_color) # hue="phenotype", palette=phenotype_to_color # s

        # add lims
        xlims = [min_consistency-0.05, max_consistency+0.05]
        ylims = [0.45, 1.05]
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        #ax.plot(lims, lims, c="k", lw=.7, ls="--", zorder=0)

        # add box with top model
        df_top_models = df_ML_integrated_performance_p[df_ML_integrated_performance_p.top_model].copy()
        if len(df_top_models)>0:

            min_consistency_top = min(df_top_models.mean_cross_splits_consistency_jaccard_genes)
            max_consistency_top = max(df_top_models.mean_cross_splits_consistency_jaccard_genes)
            min_AUC_top = min(df_top_models.AUC_test)
            max_AUC_top = max(df_top_models.AUC_test)

            #plt.plot([min_consistency_top, xlims[1]], [min_AUC_top, min_AUC_top], ls="--", c="black", lw=.7)
            #plt.plot([min_consistency_top, min_consistency_top], [min_AUC_top, ylims[1]], ls="--", c="black", lw=.7)

        # add number of models
        if add_colors is True:
            plt.text(0.95, 0.97, sum(df_ML_integrated_performance_p.correct_model), ha="right", fontsize=10, color=performanceType_to_color["correct"])
            plt.text(0.95, 0.87, sum(df_ML_integrated_performance_p.top_model), ha="right", fontsize=10, color=performanceType_to_color["top"])


        # add thresholds    
        if add_colors is True:
            ax.axvline(threshold_consistency_genes, c="gray", lw=.7, ls="--", zorder=0)
            ax.axhline(threshold_ROC_AUC, c="gray", lw=.7, ls="--", zorder=0)
        
        n_transitions = df_transitions_info_GWAS[df_transitions_info_GWAS.phenotype==phenotype].n_transitions.iloc[0]
        formatted_text = get_paper_dicotomic_pheno(phenotype).replace("agar ", "").replace("blood vs other sterile", "blood-v-sterile").replace("bronchial aspirate", "bronch. aspir.").replace("voriconazole", "vori.").replace("fluconazole", "fluc.").replace("micafungin", "mica.")

        # add title
        #ax.set_title(formatted_text+"\n#t="+str(n_transitions), fontsize=10)
        #ax.set_title(formatted_text + " (%i)"%sum(df_ML_integrated_performance_p.correct_model), fontsize=10)
        ax.set_title(formatted_text, fontsize=10)

        # add boxes
        #box = mpatches.Rectangle((0.43, 0.95), 1, 0.15, facecolor='white', edgecolor='white', linewidth=0.5)
        #ax.add_patch(box)
        #ax.text(0.47, 0.98, "#x=%i"%(sum(df_ML_integrated_performance_p.correct_model)), fontsize=10, color="black", ha='left', va='center')
        #ax.text(0.47, 0.98, "#t=%i"%n_transitions, fontsize=10, color="black", ha='left', va='center')

        # get ticcks
        yticks = [0.5, 0.7, 1.0]
        xticks = [0.4, 0.8]
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        # set lims
        if I in {7}:
            ax.set_ylabel("          test $AUC_{ROC}$", ha="left")
            ax.set_yticklabels(list(map(str, yticks)))

        elif I==0:
            ax.set_ylabel("")

        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        if I in {10}:
            ax.set_xlabel("m1-m5 feature consistency (genes)")

        else:
            ax.set_xlabel("")

        if I in set(range(7)):
            ax.set_xticklabels([])

        else:
            ax.set_xticklabels(list(map(str, xticks)))

        if I==12:
            if add_colors is True: ax.legend(bbox_to_anchor=(1.1,-0.1), loc="lower left", markerscale=2, title="performance")
            else: ax.get_legend().remove()

        else:
            if len(df_ML_integrated_performance_p)>0:
                ax.get_legend().remove()




    sns.despine()
    plt.subplots_adjust(wspace=0.1, hspace=0.35)# wspace (between cols) hspaces (btw rows)
    save_fig(fig, "overfitting_check_filt_parms_split_colors=%s.pdf"%add_colors, PlotsDir)



def plot_general_corr_train_test(df_ML_integrated_performance_filt, PlotsDir, threshold_ROC_AUC, phenotype_to_color):

    """Get general corr"""

    df_ML_integrated_performance_filt = df_ML_integrated_performance_filt.copy()

    fig = plt.figure(figsize=(4, 4))
    #ax = sns.scatterplot(data=df_ML_integrated_performance_filt, x="mean_AUC_train", y="AUC_test", edgecolor=list(df_ML_integrated_performance_filt.apply(lambda r: phenotype_to_color[r.phenotype], axis=1)), linewidth=1.1, color="none") # hue="phenotype", , palette=phenotype_to_color
    ax = sns.scatterplot(data=df_ML_integrated_performance_filt, x="mean_AUC_train", y="AUC_test", edgecolor="gray", linewidth=.8, hue="phenotype", palette=phenotype_to_color) # 

    lims = [0.45, 1.05]
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.plot(lims, lims, c="k", lw=.7, ls="--")

    for val in [threshold_ROC_AUC, 0.9]:
        ax.axvline(val, c="gray", lw=.7, ls="--")
        ax.axhline(val, c="gray", lw=.7, ls="--")

    ax.legend(bbox_to_anchor=(1.2,1), loc="upper left")
    ax.set_title("checking overfitting", fontsize=10)

    sys.exit(0)

    save_fig(fig, "overfitting_check_filt_parms.pdf", PlotsDir)

def ROC_AUC_vs_max_ROC_AUC(df_ML_integrated_performance_filt, df_ML_results_maxAUC, test_sample_set, AUC_f, PlotsDir, df_ML_results_all, df_ML_metrics_all, ML_pheno_to_color):

    """Get the max vs all"""

    # get the performance on the test_sample_set of the best parameters
    target_phenotypes = {"FLC_resistance", "VRC_resistance", "agar_invasiveness", "age_infant"}
    fields_parms = ["CVset_type_split", "sample_set", "AUC_diff_threshold", "consistency_measure", "aligner", "multiple_test_corr"]

    data_dict = {}
    for phenotype in target_phenotypes:
        df_ML_performance_p = get_filtered_df_by_dict(df_ML_integrated_performance_filt, {"phenotype":phenotype, "best_test_parmID":True})
        if len(df_ML_performance_p)!=1: raise ValueError("invalid")

        filters_dict = dict(df_ML_performance_p[df_ML_performance_p.phenotype==phenotype].iloc[0][fields_parms])
        filters_dict["phenotype"] = phenotype
        filters_dict["test_sample_set"] = test_sample_set

        df_ML_results_p = get_filtered_df_by_dict(df_ML_results_all, filters_dict)
        if len(df_ML_results_p)!=5: raise ValueError("should be one")
        df_ML_results_p = df_ML_results_p[df_ML_results_p.selected_predictors.apply(len)>0].copy()

        if AUC_f=="ROC_AUC_test_integrated_onlyCorrectModels":
            data_dict[phenotype] = {"phenotype":phenotype, "ROC_AUC_bestModel":df_ML_results_p.iloc[0][AUC_f]}

        else:
            raise ValueError("invalid")

    df_performance = pd.DataFrame(data_dict).transpose().reset_index(drop=True)

    # add max AUC
    df_max = get_filtered_df_by_dict(df_ML_results_maxAUC.copy().reset_index(drop=True), {"phenotype":target_phenotypes, "test_sample_set":test_sample_set, "max_AUC_field":AUC_f})
    df_performance["ROC_AUC_max"] = df_performance.phenotype.map(df_max.groupby("phenotype").apply(lambda df_p: max(df_p.max_ROC_AUC)))

    # figure
    fig = plt.figure(figsize=(1.5,1.5))
    ax = sns.scatterplot(df_performance, x="ROC_AUC_bestModel", y="ROC_AUC_max", hue="phenotype", palette=ML_pheno_to_color, s=30, edgecolor="gray", lw=.7)

    lims = [0.45, 1.05]
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    plt.plot(lims, lims, c="gray", lw=.7, ls="--")
    ax.legend().remove()
    ticks = [0.5, 0.75, 1]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.set_title("performance %s\n"%(test_sample_set.replace("_", " ")), fontsize=10)
    ax.set_xlabel("$AUC_{ROC}$ best model\n(optimized rep. samples)")
    ax.set_ylabel("max $AUC_{ROC}$ any model")

    sns.despine()

    save_fig(fig, "%s_max_AUC_vs_opt_AUC.pdf"%test_sample_set, PlotsDir)
    

def get_processed_ML_classifier_dfs_old(df_ML_results_all, df_ML_metrics_all, threshold_ROC_AUC, min_nmodels, ProcessedDataDir, fields_optimize_test, filters_manually_curated):

    """Gets the processed dfs"""

    def get_consistent_format(x):
        if type(x) in {float, str, int}: return str(x)
        elif type(x)==set: return ",".join(map(str, sorted(x)))
        else: raise ValueError("invalud")

    dict_processed_data_file = "%s/dict_ML_dataframes_tAUC=%s_nmod=%i_filt=%s.py"%(ProcessedDataDir, threshold_ROC_AUC, min_nmodels, "-".join(["%s=%s"%(k[0:2], get_consistent_format(filters_manually_curated[k])) for k in sorted(filters_manually_curated)]))

    print("Getting processed ML dfs")
    if file_is_empty(dict_processed_data_file):

        # redefine ROC_AUC_test_integrated_onlyCorrectModels to only apply for >=min_nmodels
        print("redefine ROC_AUC_test_integrated_onlyCorrectModels")
        def get_ROC_AUC_test_integrated_onlyCorrectModels_by_nmodels_fit(r):
            if r.nmodels_fit>=min_nmodels: return r.ROC_AUC_test_integrated_onlyCorrectModels
            else: return 0.5

        df_ML_results_all["ROC_AUC_test_integrated_onlyCorrectModels"] = df_ML_results_all[["ROC_AUC_test_integrated_onlyCorrectModels", "nmodels_fit"]].apply(get_ROC_AUC_test_integrated_onlyCorrectModels_by_nmodels_fit, axis=1)

        # get the maximum possible ROC AUC for each phens
        print("Add the max possible ROC AUC...")
        def get_max_ROC_AUC_series(df, max_AUC_field):
            fields_dict = dict(df[["sample_set", "test_sample_set", "phenotype"]].iloc[0])
            fields_dict["max_ROC_AUC"] = max(df[max_AUC_field])
            return pd.Series(fields_dict)

        df_ML_results_maxAUC_all = pd.DataFrame()
        for max_AUC_field in ["ROC_AUC_test_integrated", "cross_splits_ROC_AUC", "ROC_AUC_all", "ROC_AUC_test_integrated_onlyCorrectModels"]:
            df_ML_results_maxAUC = df_ML_results_all.groupby(["sample_set", "test_sample_set", "phenotype"]).apply(get_max_ROC_AUC_series, max_AUC_field=max_AUC_field).reset_index(drop=True)
            df_ML_results_maxAUC["max_AUC_field"] = max_AUC_field
            df_ML_results_maxAUC_all = pd.concat([df_ML_results_maxAUC_all, df_ML_results_maxAUC])

        # filter df_ML_results_all
        print("filtering df for test_sample_set...")
        filter_dict = {"test_sample_set": {"representative_samples"}} 
        df_ML_results_all = get_filtered_df_by_dict(df_ML_results_all, filter_dict, check_all_present=True)

        # get df with one row per parmID & pheno
        print("adding integrated performance")
        df_ML_integrated_performance = get_ML_integrated_performance_df(df_ML_results_all, "cross_splits_ROC_AUC", "ROC_AUC_test_integrated_onlyCorrectModels", min_nmodels) 
        df_ML_integrated_performance = df_ML_integrated_performance.merge(df_ML_results_all[["test_parmID", "phenotype"] + fields_optimize_test].drop_duplicates(), on=["test_parmID", "phenotype"], how="left", validate="one_to_many")
        check_no_nans_in_df(df_ML_integrated_performance)

        # add indices
        for f, sorted_vals in [("aligner", ["bwa_mem", "bowtie2_local", "hisat2_no_spliced"]),
                               ("CVset_type_split", ["balanced", "random"]),
                               ("consistency_measure", ["jaccard_variants", "jaccard_samples", "equal_variants", "equal_samples", "none"]),
                               ("AUC_diff_threshold", [0.05, 0.1, 0.15, 0.2]),
                               ("sample_set", ["representative_samples", "all_samples"]),
                               ("multiple_test_corr", ["maxT_min", "maxT_mean", "none"])]:

            map_dict = dict(zip(sorted_vals, range(len(sorted_vals))))
            df_ML_integrated_performance["%s_idx"%f] = df_ML_integrated_performance[f].map(map_dict)

        check_no_nans_in_df(df_ML_integrated_performance)
            
        # adding combinations of fields
        for field_A, field_B in [("consistency_measure", "AUC_diff_threshold"), 
                                 ("CVset_type_split", "sample_set")]:
            df_ML_integrated_performance["%s & %s"%(field_A, field_B)] = df_ML_integrated_performance[field_A].apply(str) + " & " +  df_ML_integrated_performance[field_B].apply(str)

        # print filters initial
        print("filters initial:\n", df_ML_results_all[["test_parmID", "phenotype", "test_sample_set"]].drop_duplicates().groupby(["phenotype", "test_sample_set"]).apply(len))

        # filter the df by manually curated filters
        print("filtering manual curation...")
        df_ML_results_filt = get_filtered_df_by_dict(df_ML_results_all, filters_manually_curated, check_all_present=True)
        print("filters ater:\n", df_ML_results_filt[["test_parmID", "phenotype", "test_sample_set"]].drop_duplicates().groupby(["phenotype", "test_sample_set"]).apply(len))

        # get the manually filtered parms
        all_test_parmID = set(df_ML_results_filt.test_parmID)
        df_ML_metrics_filt = df_ML_metrics_all[df_ML_metrics_all.test_parmID.isin(all_test_parmID)].copy()
        df_ML_integrated_performance_filt = df_ML_integrated_performance[df_ML_integrated_performance.test_parmID.isin(all_test_parmID)].copy()


        # get the best parms
        print("getting best parms")    
        pheno_to_max_AUC_test = {}
        pheno_to_best_parm = {}
        pheno_to_best_AUC_test = {}

        phenotypes_ML = ["FLC_resistance", "VRC_resistance", "agar_invasiveness", "age_infant"]
        for pheno in phenotypes_ML:

            # record the max AUC across all parameters
            pheno_to_max_AUC_test[pheno] = max(df_ML_integrated_performance[df_ML_integrated_performance.phenotype==pheno].AUC_test)

            # get performance of the pheno
            df_performance_pheno = df_ML_integrated_performance_filt[(df_ML_integrated_performance_filt.phenotype==pheno) & (df_ML_integrated_performance_filt.AUC_test>=threshold_ROC_AUC) & (df_ML_integrated_performance_filt.mean_AUC_train>=threshold_ROC_AUC)].copy()

            # prints
            print(df_performance_pheno[["phenotype", "mean_AUC_train", "AUC_test", "n_models_fit"]], "\n")

            # get best pheno
            idx_sort_fields = ["sample_set_idx", "multiple_test_corr_idx", "CVset_type_split_idx", "consistency_measure_idx", "AUC_diff_threshold_idx", "aligner_idx"]
            idx_asc_fields = [True]*len(idx_sort_fields)
            df_performance_pheno_best = get_df_like_first_row(df_performance_pheno, ["AUC_test", "n_models_fit", "mean_AUC_train"] + idx_sort_fields, [False, False, False] + idx_asc_fields)
            if len(df_performance_pheno_best)!=1: raise ValueError("should be diffeent cases: %s"%df_performance_pheno_best)

            # record
            pheno_to_best_parm[pheno] = df_performance_pheno_best.iloc[0].test_parmID
            pheno_to_best_AUC_test[pheno] = df_performance_pheno_best.iloc[0].AUC_test

        print("phenotypes ML:", pheno_to_best_parm, "\n", pheno_to_max_AUC_test, "\n",  pheno_to_best_AUC_test)

        # add missings
        for p in set(df_ML_integrated_performance.phenotype).difference(set(phenotypes_ML)):
            pheno_to_best_parm[p] = None

        # add to dfs
        for df in [df_ML_integrated_performance, df_ML_integrated_performance_filt, df_ML_metrics_filt, df_ML_results_filt]:
            df["best_test_parmID"] = df[["test_parmID", "phenotype"]].apply(lambda r: r.test_parmID == pheno_to_best_parm[r.phenotype], axis=1)

        final_dict = {"df_ML_results_maxAUC": df_ML_results_maxAUC_all,
                      "df_ML_integrated_performance":df_ML_integrated_performance,
                      "df_ML_integrated_performance_filt": df_ML_integrated_performance_filt,
                      "df_ML_metrics_filt": df_ML_metrics_filt,
                      "df_ML_results_filt": df_ML_results_filt,
                      "phenotypes_ML": phenotypes_ML

                      }

        # save
        print("saving...")
        save_object(final_dict, dict_processed_data_file)

    return load_object(dict_processed_data_file)

def load_ML_results_and_metrics_all(ProcessedDataDir, fields_optimize_test):

    """loads raw dfs"""


    df_ML_results_all_file = "%s/df_ML_results_all.py"%ProcessedDataDir
    df_ML_metrics_all_file = "%s/df_ML_results_all.py.per_threshold_metrics.py"%ProcessedDataDir

    df_ML_results_all_file_processed = df_ML_results_all_file+".added_parmID.py"
    df_ML_metrics_all_file_processed = df_ML_metrics_all_file+".added_parmID.py"

    print("Getting all ML dfs...")
    if file_is_empty(df_ML_results_all_file_processed) or file_is_empty(df_ML_metrics_all_file_processed):

        # load the ML results
        df_ML_results_all = load_object(df_ML_results_all_file)
        df_ML_metrics_all = load_object(df_ML_metrics_all_file)

        # add IDs
        df_ML_results_all = get_ML_df_with_test_parmID(df_ML_results_all, ["test_split_ID", "phenotype"], fields_optimize_test)
        df_ML_metrics_all = get_ML_df_with_test_parmID(df_ML_metrics_all, ["test_split_ID", "phenotype", "eval_data", "threshold"], fields_optimize_test)

        save_object(df_ML_results_all, df_ML_results_all_file_processed)
        save_object(df_ML_metrics_all, df_ML_metrics_all_file_processed)

    return load_object(df_ML_results_all_file_processed), load_object(df_ML_metrics_all_file_processed)

def report_acuraccy_all_samples_best_models(df_ML_results_all, df_ML_metrics_all, phenotypes_ML, df_ML_integrated_performance_filt, fields_optimize_test):

    """reports accuracy on all datasets"""

    # checks
    fields_parms = [f for f in fields_optimize_test if f!="test_sample_set"]

    #df_ML_results_print

    for phenotype in phenotypes_ML:

        # get performance rep samples
        df_performance = df_ML_integrated_performance_filt[(df_ML_integrated_performance_filt.phenotype==phenotype) & (df_ML_integrated_performance_filt.best_test_parmID)].copy()
        if len(df_performance)!=1: raise ValueError("should be 1")
        series_performance_rep = df_performance.iloc[0]

        # get the performance on all samples
        print()
        for test_sample_set in ["representative_samples", "all_samples"]:

            filter_dict = {**{f : series_performance_rep[f] for f in fields_parms}, **{"test_sample_set":test_sample_set, "phenotype":phenotype}}
            df_ML_results = get_filtered_df_by_dict(df_ML_results_all, filter_dict).rename(columns={"ROC_AUC_test_integrated_onlyCorrectModels":"AUC_test", "cross_splits_ROC_AUC":"AUC_train", "expected_genes_str":"exp_genes", "redundant_samples_ROC_AUC":"indep_AUC", "redundant_samples_n_transitions":"indep_n"})
            df_ML_results["n_predictors"] = df_ML_results.selected_predictors.apply(len)
            df_ML_results["mod"] = df_ML_results.model_name.apply(lambda x: x.split("|")[0])
            df_ML_results["type_p"] = df_ML_results.type_phenotype.map({"phenotype":"pheno", "phenotype_transition":"trans"})
            df_ML_results["tdef"] = df_ML_results.transition_definition.map({"each_pheno_clade":"all", "two_close_balanced":"two"})



            df_ML_results = df_ML_results[df_ML_results.n_predictors>0].copy()

            df_ML_results.index = [""]*len(df_ML_results)
            print(phenotype, test_sample_set, "nmodels fit:", df_ML_results.nmodels_fit.iloc[0], "indep_n:", df_ML_results.indep_n.iloc[0], "AUCs:")
            for f in ["AUC_test", "AUC_train", "indep_AUC"]: df_ML_results[f] = df_ML_results[f].apply(lambda x: "." + ("%.2f"%x).split(".")[1] )
            print(df_ML_results[["AUC_test", "AUC_train", "exp_genes", "n_predictors", "indep_AUC", "type_p", "tdef", "mod"]]) # transition_definition








######################

