import os
import subprocess
import numpy as np
import pandas as pd
import glob
import sys
from multiprocessing import Pool
import multiprocessing

###########################
# This script will lift over your sumstats from hg19 to hg38
# 10.04.2023 (provided by Yasmina Mekki)
# https://github.com/hakyimlab/summary-gwas-imputation
# https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation

###########################
######### Example #########
###########################
# First create a conda virtual env:
# conda env create -f /home/gokala/programs/summary-gwas-imputation/src/conda_env.yaml
# conda activate imlabtools
# Then run this script:
# python summary_statistics_harmonization.py
###########################
    
def summary_statistics_harmonization():
    
    cmd = " ".join(['python /home/gokala/programs/summary-gwas-imputation/src/gwas_parsing.py',
                    '-gwas_file /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN.tab',
                    '-output_column_map SNP variant_id',
                    '-output_column_map A2 non_effect_allele',
                    '-output_column_map A1 effect_allele',
                    '-output_column_map est effect_size',
                    '-output_column_map Pval_Estimate pvalue',
                    '-output_column_map SE standard_error',
                    '-output_column_map BP position',
                    '-output_column_map CHR chromosome',
                    '--chromosome_format',
                    '--insert_value sample_size 1745597',
                    '-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pvalue effect_size standard_error sample_size',
                    '-liftover /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/hg19ToHg38.over.chain.gz',
                    '-output /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN_hg38.txt.gz'
                   ])

    print("command summary statistics harmonization : {}".format(cmd))
    p = subprocess.check_call(cmd, shell=True)

if __name__ == "__main__":
    
    summary_statistics_harmonization()
