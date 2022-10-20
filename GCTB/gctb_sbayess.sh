#!/bin/bash
#
# This script runs GCTB SBayesS on the dyslexia & rhythm impairment
# meta-analysis results.

module load gctb/2.02
genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/GCTB/"

# Run SBayesS on GWAMA, dys and rhy_imp sumstats

#gctb --sbayes S \
#     --ldm /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse \
#     --gwas-summary /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/dys_aysnch_NGWAMA_results4GCTB.ma \
#     --burn-in 2000 \
#     --out dys_asynch > test.log 2>&1

#gctb --sbayes S \
#     --ldm /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse \
#     --gwas-summary /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythm_impairment_4gctb.ma \
#     --burn-in 2000 \
#     --out rhythm_impairment > rhythm_impairment.log 2>&1

#gctb --sbayes S \
#     --ldm /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse \
#     --gwas-summary /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia.filtered.2_4gctb.ma \
#     --burn-in 2000 \
#     --out dyslexia > dyslexia.log 2>&1

# Run SBayesS on Genomic SEM Common Pathways sumstat (dys and rhy_imp meta-analysis results)

gctb --sbayes S \
     --mldm /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
     --gwas-summary ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun_4gctb.ma \
     --burn-in 2000 \
     --out ${outDir}rhyimp_dys_v2_secondrun_genomicSEM > ${outDir}rhyimp_dys_v2_secondrun_genomicSEM.log 2>&1
