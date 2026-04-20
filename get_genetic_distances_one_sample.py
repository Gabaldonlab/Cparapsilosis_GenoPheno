#!/usr/bin/env python

# Calculates genetic distances per window for one sample vs all the others

# imports
import os, sys, re, time
import pandas as pd
import copy as cp
import multiprocessing as multiproc
import shutil
import numpy as np
import itertools

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if not os.path.exists(ParentDir):  ParentDir = "/gpfs/projects/bsc40/mschikora"

# define dir
CurDir = "%s/Cparapsilosis_popgen"%ParentDir
sys.path.append(CurDir); import Cparapsilosis_popgen_functions as fun

# parse args
df_var_densities_file, distances_file, threads, target_runID = sys.argv[1:]

# debug
if not fun.file_is_empty(distances_file):
	sys.exit(0)

# load df
print("loading densities")
df_var_densities_all = fun.load_object(df_var_densities_file)
print("%i windows total"%(len(set(df_var_densities_all.ID_window))))

# filter to keep high confidence regions
df_var_densities_all = df_var_densities_all[(df_var_densities_all.fraction_repeats<=0.05) & (df_var_densities_all.fraction_N_bases<=0.05) & (df_var_densities_all.median_mappability>=0.9)].copy()

# filter to keep regions with no aneuploidies in any (should have suff coverage)
windows_aneuploidies = set(df_var_densities_all[(df_var_densities_all.median_coverage_rel_to_gDNA<0.25)].ID_window)

print("%i/%i windows with aneuploidies"%(len(windows_aneuploidies), len(set(df_var_densities_all.ID_window))))
df_var_densities_all = df_var_densities_all[~df_var_densities_all.ID_window.isin(windows_aneuploidies)].copy()
all_windows = set(df_var_densities_all.ID_window)

# add general fields
df_var_densities_all["len_window"] = df_var_densities_all.end - df_var_densities_all.start

# calculate stats for each window
print("prepare the inputs")
all_runIDs = sorted(set(df_var_densities_all.runID))

inputs_fn = []
for type_var in sorted(set(df_var_densities_all.type_var)):
	for zygosity in sorted(set(df_var_densities_all.zygosity)):

		df_var_densities = df_var_densities_all[(df_var_densities_all.type_var==type_var) & (df_var_densities_all.zygosity==zygosity)].copy()
		inputs_fn += [(type_var, zygosity, origin_runID, target_runID, df_var_densities[df_var_densities.runID.isin({origin_runID, target_runID})].copy()) for origin_runID in all_runIDs]

# run
print("running in parallel")
with multiproc.Pool(int(threads)) as pool:
    df_stats_distances = pd.concat(pool.starmap(fun.get_df_stats_distances_one_pair_samples_and_types_vars, inputs_fn, chunksize=1))
    pool.close()
    pool.terminate()

# save
print("saving")
fun.save_object(df_stats_distances, distances_file)