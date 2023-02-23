#!/bin/bash
#
# This script will convert GWAS sumstats to
# BED format.
#
# Gokberk Alagoz, 08.12.22
#
############################

module load bedtools/2.29.2

inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/"

# Genomic SEM CPM results

#awk -v OFS='\t' 'NR>1 {print "chr"$2,$3-1,$3,$12,$13,$15,$16,$18}' ${inDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr.tab > ${outDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp.bed
#sort-bed ${outDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp.bed > ${outDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp.sorted.bed

# Genomic SEM IPM results - dys

#awk -v OFS='\t' 'NR>1 {print "chr"$2,$3-1,$3,$12,$13,$15,$16,$18}' ${inDir}GenomicSEM_multivarGWAS_dys_IPM_GCcorr.tab > ${outDir}GenomicSEM_multivarGWAS_dys_IPM.bed
#sort-bed ${outDir}GenomicSEM_multivarGWAS_dys_IPM.bed > ${outDir}GenomicSEM_multivarGWAS_dys_IPM.sorted.bed

# Genomic SEM IPM results - rhyimp

awk -v OFS='\t' 'NR>1 {print "chr"$2,$3-1,$3,$12,$13,$15,$16,$18}' ${inDir}GenomicSEM_multivarGWAS_rhyimp_IPM_GCcorr.tab > ${outDir}GenomicSEM_multivarGWAS_rhyimp_IPM.bed
sort-bed ${outDir}GenomicSEM_multivarGWAS_rhyimp_IPM.bed > ${outDir}GenomicSEM_multivarGWAS_rhyimp_IPM.sorted.bed
