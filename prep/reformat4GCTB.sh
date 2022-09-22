#!/bin/bash

# This script will reformat MA sumstats for
# GCTB SBayesS.

# Paths
MA="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/"
rawss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"

# Meta-analysis
awk '{print $2, $5, $6, $7, $11, $12, $13, $9}' ${MA}NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt > ${MA}dys_aysnch_NGWAMA_results_4gctb.txt

# Dyslexia
awk '{print $19, $24, $25, $20, $4, $5, $3, $7+$8}' ${rawss}dyslexia.filtered.2.dat > ${rawss}dyslexia.filtered.2_tmp.dat
# atm freq column is dose.b, but GCTB needs A1 freq, thus convert to dose.a
awk '{$4 = 1-$4} 1' ${rawss}dyslexia.filtered.2_tmp.dat > ${rawss}dyslexia.filtered.2_tmp2.dat
echo -e "SNPID A1 A2 A1.freq BETA SE PVAL N" | cat - ${rawss}dyslexia.filtered.2_tmp2.dat > ${rawss}dyslexia.filtered.2_4gctb.ma
rm ${rawss}dyslexia.filtered.2_tmp.dat
rm ${rawss}dyslexia.filtered.2_tmp2.dat

# Rhythm
awk '{print $42, $38, $39, $17, $4, $5, $3, $43}' ${rawss}clap_to_beat.merged.1kgenomes_02222019_wtotalN.txt > ${rawss}clap_to_beat.merged.1kgenomes_02222019_tmp.txt
#flip effect to get rhythm impairment
awk 'NR>1 {$5 = $5*-1} 1' ${rawss}clap_to_beat.merged.1kgenomes_02222019_tmp.txt > ${rawss}rhythm_impairment_4gctb.ma
sed -Ei '1s/snp.id/SNPID/;1s/ref_allele/A1/;1s/effect_allele/A2/;1s/freq.a/A1.freq/;1s/effect/BETA/;1s/stderr/SE/;1s/pvalue/PVAL/;1s/0/N/' ${rawss}rhythm_impairment_4gctb.ma
