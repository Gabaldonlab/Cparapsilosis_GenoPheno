#!/usr/bin/env python

# calculates the m,appability file per chrom and per whole genome. Should be run on perSVade_env

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
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.append(perSVade_dir); import sv_functions as fun

# parse args
ref_genome, threads, final_file, outdir = sys.argv[1:]

# get as full path ome paths
outdir = fun.get_fullpath(outdir)

# log
print("getting mappability into %s."%outdir)
fun.make_folder(outdir)

# define the final file
if not fun.file_is_empty(final_file):
    print("WARNING: already ran, exiting...")
    sys.exit(0)

# process inputs
threads = int(threads)

# make tmpdir
tmpdir = "%s/tmp"%outdir
fun.make_folder(tmpdir)

# link
ref_genome_linked = "%s/genome.fasta"%tmpdir
fun.soft_link_files(ref_genome, ref_genome_linked)

# running mappability
fun.printing_verbose_mode = True
fun.generate_genome_mappability_file(ref_genome_linked, replace=False, threads=threads)

# keep some files
for f in ["genome.fasta.mappability_per_position_per_chromosome", "genome.fasta.mappability_per_position.bed"]:
	os.rename("%s/%s"%(tmpdir, f), "%s/%s"%(outdir, f))

# clean and exit
fun.delete_folder(tmpdir)
open(final_file, "w").write("finished...")
print("finished")


