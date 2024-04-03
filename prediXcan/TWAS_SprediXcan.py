import os
import subprocess
import numpy as np
import pandas as pd
import glob
import sys
from multiprocessing import Pool
import multiprocessing

###########################
# This script will run SPrediXcan on 
# dyslexia and rhythm impairment CPM
# sumstats.
# 10.04.2023 (provided by Yasmina Mekki)

###########################
######### Example #########
###########################
# conda activate imlabtools
#
# python TWAS_SprediXcan.py 
#
###########################

def perform_twas(parameters):
    
    model_db_path, covariance_path, gwas_file, snp_column, effect_allele_column, non_effect_allele_column, beta_column, pvalue_column, output_file = parameters

    cmd = " ".join(['python /home/gokala/programs/MetaXcan/software/SPrediXcan.py',
                    '--model_db_path %s' % model_db_path,
                    '--covariance %s' % covariance_path,
                    '--gwas_file %s' % gwas_file,
                    '--snp_column %s' % snp_column,
                    '--effect_allele_column %s' % effect_allele_column,
                    '--non_effect_allele_column %s' % non_effect_allele_column,
                    '--beta_column %s' % beta_column,
                    '--pvalue_column %s' % pvalue_column,
                    '--output_file %s' % output_file,
                    #'--verbosity 9',
                    #'--throw'
                   ])

    print("command Predict : {}".format(cmd))
    p = subprocess.check_call(cmd, shell=True)
    
if __name__ == "__main__":
    
    weights_path = '/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/JTI'
    model_db_paths = glob.glob(os.path.join(weights_path, '*'))
    gwas_file = '/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN_hg38.txt.gz'
    snp_column = 'variant_id'
    effect_allele_column = 'effect_allele'
    non_effect_allele_column = 'non_effect_allele'
    beta_column = 'effect_size'
    pvalue_column = 'pvalue'
    
    list_cmd = []
    for model_db_path in model_db_paths:
        
        covariance_path = os.path.join(weights_path, os.path.basename(model_db_path).replace(".db", ".txt.gz"))
        output_file = os.path.join('/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/', os.path.basename(model_db_path).replace(".db", "_assoc.txt"))
        
        list_cmd.append([model_db_path, covariance_path, gwas_file, snp_column, effect_allele_column,
                        non_effect_allele_column, beta_column, pvalue_column, output_file]) 

    pool = Pool(processes = 31)
    pool.map(perform_twas, list_cmd)
    pool.close()
    pool.join()

