#!/bin/bash
#
# Gokberk Alagoz, March 2022 - updated Nov 2022

#-----Munge Sumstats-----

# Reformat GWAS summary stats before computing LDSC intercept
# munge_sumstats.py is from github.com/bulik/ldsc
# Based on the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# After final updates on the grid, it's easier to run munge on lux13. Follow the steps 
# below and then run munge.

#module purge
#module load miniconda/3.2021.10 ldsc/v1.0.1
#conda activate ldsc

#-----Variables-----
hapmap="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist"
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data"
genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM"

# 1) Rhythym impairment

#munge_sumstats.py \
# --sumstats ${inDir}/rhythm_impairment_lambdaGCcorrected_wtotalNandOR.dat \
# --snp snp.id \
# --a1 effect_allele \
# --a2 ref_allele \
# --frq maf \
# --N-cas-col im.num.1 \
# --N-con-col im.num.0 \
# --p pvalue \
# --signed-sumstats OR,1 \
# --merge-alleles ${hapmap} \
# --out ${inDir}/munged/rhythm_impairment \
# --chunksize 500000

# 2) Dyslexia

#munge_sumstats.py \
# --sumstats ${inDir}/dyslexia.filtered.2_lambdaGCcorrected_wtotalN.dat \
# --snp rsid \
# --a1 alleleB \
# --a2 alleleA \
# --N-cas-col im.num.1 \
# --N-con-col im.num.0 \
# --p pvalue \
# --signed-sumstats OR,1 \
# --merge-alleles ${hapmap} \
# --out ${inDir}/munged/dyslexia \
# --chunksize 500000

# 3) Genomic SEM common factor results

munge_sumstats.py \
	--sumstats ${genSEM}/GenomicSEM_multivarGWAS_dys_rhyimp_lambdaGCcorrected_forMunge.tab \
	--signed-sumstats Z_Estimate,0 \
	--p Pval_Estimate \
	--frq MAF \
	--out ${inDir}/munged/GenomicSEM_multivarGWAS_dys_rhyimp \
	--merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
	--chunksize 500000

