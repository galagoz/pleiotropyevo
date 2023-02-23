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

#awk 'NR>1 {print $19, $16, $17, $25, $24, $20, $7+$8, $4, $5, $3}' ${sumstats}dyslexia.filtered.2.dat > ${sumstats}dyslexia_forgenomicSEM1.txt
#awk '{gsub("chr", "", $2); print}' ${sumstats}dyslexia_forgenomicSEM1.txt > ${sumstats}dyslexia_forgenomicSEM2.txt
#echo -e "SNPID CHR BP EA OA EAF N BETA SE P" | cat - ${sumstats}dyslexia_forgenomicSEM2.txt > ${sumstats}dyslexia_forgenomicSEM.txt

#chmod 777 ${sumstats}dyslexia_forgenomicSEM.txt
#rm ${sumstats}dyslexia_forgenomicSEM1.txt
#rm ${sumstats}dyslexia_forgenomicSEM2.txt

# Rhythym impairment

awk 'NR>1 {print $42, $37, $24, $39, $38, $18, $7+$9, $4, $5, $3}' ${sumstats}rhythm_impairment_raw_sumstat.txt > ${sumstats}rhythm_impairment_forgenomicSEM1.txt
awk '{gsub("chr", "", $2); print}' ${sumstats}rhythm_impairment_forgenomicSEM1.txt > ${sumstats}rhythm_impairment_forgenomicSEM2.txt
echo -e "SNPID CHR BP EA OA EAF N BETA SE P" | cat - ${sumstats}rhythm_impairment_forgenomicSEM2.txt > ${sumstats}rhythm_impairment_forgenomicSEM.txt

chmod 777 ${sumstats}rhythm_impairment_forgenomicSEM.txt
rm ${sumstats}rhythm_impairment_forgenomicSEM1.txt
rm ${sumstats}rhythm_impairment_forgenomicSEM2.txt
