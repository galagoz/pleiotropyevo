#!/bin/bash
#
# This script runs GCTB SBayesS on the dyslexia & rhythm impairment
# meta-analysis results.

module load gctb/2.02

#gctb --sbayes S \
#     --ldm /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse \
#     --gwas-summary /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/dys_aysnch_NGWAMA_results4GCTB.ma \
#     --burn-in 2000 \
#     --out dys_asynch > test.log 2>&1

gctb --sbayes S \
     --ldm /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse \
     --gwas-summary /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythm_impairment_4gctb.ma \
     --burn-in 2000 \
     --out rhythm_impairment > rhythm_impairment.log 2>&1

gctb --sbayes S \
     --ldm /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse \
     --gwas-summary /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia.filtered.2_4gctb.ma \
     --burn-in 2000 \
     --out dyslexia > dyslexia.log 2>&1
