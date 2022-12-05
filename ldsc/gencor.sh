#!/bin/bash
#
# Gokberk Alagoz, March 2022 - last update Nov 2022

#-----
# PATHS

mungedDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged"
ldscores="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/eur_w_ld_chr/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor"

#-----
#module purge
#module load miniconda/3.2021.10 ldsc/v1.0.1
#conda activate ldsc
#-----
# dys vs. rhy_imp gencor

ldsc.py \
 --rg ${mungedDir}/dyslexia.sumstats.gz,${mungedDir}/rhythm_impairment.sumstats.gz \
 --ref-ld-chr ${ldscores} \
 --w-ld-chr ${ldscores} \
 --out ${outDir}/dys_rhyimp

#-----
# self gencors for LAVA

ldsc.py \
 --rg ${mungedDir}/dyslexia.sumstats.gz,${mungedDir}/dyslexia.sumstats.gz \
 --ref-ld-chr ${ldscores} \
 --w-ld-chr ${ldscores} \
 --out ${outDir}/dys_self

#ldsc.py \
# --rg ${mungedDir}/rhythm_impairment.sumstats.gz,${mungedDir}/rhythm_impairment.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --w-ld-chr ${ldscores} \
# --out ${outDir}/rhyimp_self
