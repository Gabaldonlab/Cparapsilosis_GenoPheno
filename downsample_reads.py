#!/usr/bin/env python

# Takes some reads and donwsamples them to nreads. Run wuth Cparapsilosis_popgen_env

# env
import sys, os

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if not os.path.exists(ParentDir): ParentDir = "/gpfs/projects/bsc40/mschikora"

# define dirs
CurDir = "%s/Cparapsilosis_popgen"%ParentDir

# import functions
sys.path.append(CurDir); import Cparapsilosis_popgen_functions as fun

# parse args
input_reads1, input_reads2, output_reads1, output_reads2, nreads = sys.argv[1:]
nreads = int(nreads)

# compress
fun.get_downsampled_reads(input_reads1, output_reads1, nreads)
fun.get_downsampled_reads(input_reads2, output_reads2, nreads)

print("SUCCESS.")