#!/bin/bash
#
# Gokberk Alagoz
#
# This script will reformat dys and rhy_imp sumstats for
# local genetic correlations analysis with LAVA.

#-----
# PATHS
sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"

# Dyslexia

awk 'NR>1 {print $19, $25, $24, $8, $4, $3}'  ${sumstats}dyslexia.filtered.2_lambdaGCcorrected.dat > ${sumstats}dyslexia.filtered.2_lambdaGCcorrected_forLAVA1.dat
echo -e "SNP A1 A2 N BETA P" | cat - ${sumstats}dyslexia.filtered.2_lambdaGCcorrected_forLAVA1.dat > ${sumstats}dyslexia_forLAVA.txt

chmod 777 ${sumstats}dyslexia_forLAVA.txt
rm ${sumstats}dyslexia.filtered.2_lambdaGCcorrected_forLAVA1.dat

# Rhythm impairment

#awk 'NR>1 {print $42, $39, $38, $7, $4, $3}' ${sumstats}rhythm_impairment_lambdaGCcorrected.dat > ${sumstats}rhythm_impairment_raw_sumstat_forLAVA1.txt
#echo -e "SNP A1 A2 N BETA P" | cat - ${sumstats}rhythm_impairment_raw_sumstat_forLAVA1.txt > ${sumstats}rhythm_impairment_forLAVA.txt

#chmod 777 ${sumstats}rhythm_impairment_forLAVA.txt
#rm ${sumstats}rhythm_impairment_raw_sumstat_forLAVA1.txt
