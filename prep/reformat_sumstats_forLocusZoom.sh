#!/bin/bash
#
# This script will reformat dyslexia, rhythm impairment and 
# Genomic SEM sumstats for LocusZoom:
# i) add a marker column
# ii) sort by chr and pos
#####

sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"
genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"

module load bedtools/2.29.2

# 1) Dyslexia
# already has a chrpos column, only needs to be sorted
#sort -k16,16 -k17,17n ${sumstats}dyslexia.filtered.2_wtotalN_forLocusZoom.txt > ${sumstats}dyslexia.filtered.2_wtotalN_forLocusZoom.sorted.txt

# 2) Rhythm impairment
#awk '{print $0"\t"$37":"$24}' ${sumstats}rhythm_impairment_raw_sumstat_wtotalN.txt > ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker.txt
#sed -Ei '1s/chr:position/chrpos/' ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker.txt
#awk -F '[\t ]+' -v OFS='\t' '{$1=$1}1' ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker.txt > ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker2.txt
#sort -k37,37 -k24,24n ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker2.txt > ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker2.sorted.txt
#echo | head -n1 ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker.txt > ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker.sorted.txt
#cat ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker2.sorted.txt >> ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker.sorted.txt
awk -F '[\t ]+' -v OFS='\t' '{$1=$1}1' ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker.sorted.txt > ${sumstats}rhythm_impairment_raw_sumstat_wtotalNandMarker.sorted.tab

# 3) Genomic SEM CPM results
#awk '{print $0"\t"$2":"$3}' ${genSEM}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN.tab > ${genSEM}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withNandMarker.tab
#sed -Ei '1s/CHR:BP/chrpos/' ${genSEM}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withNandMarker.tab
#awk -F '[\t ]+' -v OFS='\t' '{$1=$1}1' ${genSEM}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withNandMarker.tab > ${genSEM}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withNandMarker_2.tab

# 4) Genomic SEM rhyimp IPM results
#awk '{print $0"\t"$2":"$3}' ${genSEM}GenomicSEM_multivarGWAS_rhyimp_IPM_GCcorr_withN.tab > ${genSEM}GenomicSEM_multivarGWAS_rhyimp_IPM_GCcorr_withNandMarker.tab
#sed -Ei '1s/CHR:BP/chrpos/' ${genSEM}GenomicSEM_multivarGWAS_rhyimp_IPM_GCcorr_withNandMarker.tab
#awk -F '[\t ]+' -v OFS='\t' '{$1=$1}1' ${genSEM}GenomicSEM_multivarGWAS_rhyimp_IPM_GCcorr_withNandMarker.tab > ${genSEM}GenomicSEM_multivarGWAS_rhyimp_IPM_GCcorr_withNandMarker_2.tab

# 5) Genomic SEM dys IPM results
#awk '{print $0"\t"$2":"$3}' ${genSEM}GenomicSEM_multivarGWAS_dys_IPM_GCcorr_withN.tab > ${genSEM}GenomicSEM_multivarGWAS_dys_IPM_GCcorr_withNandMarker.tab
#sed -Ei '1s/CHR:BP/chrpos/' ${genSEM}GenomicSEM_multivarGWAS_dys_IPM_GCcorr_withNandMarker.tab
#awk -F '[\t ]+' -v OFS='\t' '{$1=$1}1' ${genSEM}GenomicSEM_multivarGWAS_dys_IPM_GCcorr_withNandMarker.tab > ${genSEM}GenomicSEM_multivarGWAS_dys_IPM_GCcorr_withNandMarker_2.tab
