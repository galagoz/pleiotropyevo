#!/bin/bash
#
# This script will extract CHR, POS and rsID
# columns from .clumped files.
##############################
# PATHS
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/plink/"

##############################

awk -v OFS='\t' '{print $1, $4, $3}' ${inDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN.tab.clumped > ${inDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN.clumped.leadSNP_genomicCoordinates.tab
