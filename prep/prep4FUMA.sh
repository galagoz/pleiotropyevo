#!/bin/bash
#
# This script will add N column to Genomic SEM meta-analysis
# output file.
#
#####

genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"

awk -F" " 'NR==FNR{a[NR]=$2;next} {print $0"\t"a[FNR]}' ${genSEM}perSNP_Ns_4genomicSEM.txt ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun.tab > ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun_withN.tab
