#!/bin/bash
#
# Gokberk Alagoz - February, 2022
#
# This script will remove EAF (effect allele freq.)
# and Z columns from the Meta-analysis sumstats, 
# because LDSC munge needs a single MAF and effect column.
# Also, changes colnames to A1 to A2, so that LDSC can
# recognize these fields.

awk -i inplace '{$0=gensub(/\s*\S+/,"",7)}1' /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt

sed -Ei '1s/EA/A1/;1s/OA/A2/' /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt

awk -i inplace '{$0=gensub(/\s*\S+/,"",13)}1' /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt
