#!/bin/bash
#
# Gokberk Alagoz
#
# This script will extract and order required columns from dyslexia and rhythm
# impairment summary statictics. Then will rename field names accordingly.

#-----
# PATHS
sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"

# Dyslexia

awk 'NR>1 {print $19, $16, $17, $25, $24, $20, $7+$8, $4/$5, $3}' ${sumstats}dyslexia.filtered.2_lambdaGCcorrected.dat > ${sumstats}dyslexia_forNweighted1.txt
awk '{gsub("chr", "", $2); print}' ${sumstats}dyslexia_forNweighted1.txt > ${sumstats}dyslexia_forNweighted2.txt
echo -e "SNPID CHR BP EA OA EAF N Z P" | cat - ${sumstats}dyslexia_forNweighted2.txt > ${sumstats}dyslexia_forNweighted.txt

chmod 777 ${sumstats}dyslexia_forNweighted.txt
rm ${sumstats}dyslexia_forNweighted1.txt
rm ${sumstats}dyslexia_forNweighted2.txt

# Rhythym impairment

#awk 'NR>1 {print $42, $37, $24, $39, $38, $18, $7+$9, $4/$5, $3}' ${sumstats}rhythm_impairment_lambdaGCcorrected.dat > ${sumstats}rhythm_impairment_forNweighted1.txt
#awk '{gsub("chr", "", $2); print}' ${sumstats}rhythym_impairment_forNweighted1.txt > ${sumstats}rhythm_impairment_forNweighted2.txt
#echo -e "SNPID CHR BP EA OA EAF N Z P" | cat - ${sumstats}rhythym_impairment_forNweighted2.txt > ${sumstats}rhythm_impairment_forNweighted.txt

#chmod 777 ${sumstats}rhythm_impairment_forNweighted.txt
#rm ${sumstats}rhythm_impairment_forNweighted1.txt
#rm ${sumstats}rhythm_impairment_forNweighted2.txt
