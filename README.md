# Overview

Code from the study "Dissection of genotype-phenotype relationships in Candida parapsilosis uncovers drivers of clinically-relevant traits". Check the Methods of that paper to understand the code. 

# Software environment

To install the necessary dependencies, create the conda environment with:

`conda env create --file ./Cparapsilosis_popgen_env.yml --name Cparapsilosis_popgen_env`

Activate environment with:

`conda activate Cparapsilosis_popgen_env`

# Running

Once all dependencies have been installed, you may reproduce our analyses by running:

`python run_pipelines.py`

Also, you may reproduce the figures and tables of the paper with:

`python figures_paper.py`


