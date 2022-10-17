#!/bin/bash

# This script will reformat MA sumstats for
# GCTB analysis.

# Paths
MA="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/"
rawss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"
genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"

# Meta-analysis
#awk '{print $2, $5, $6, $7, $11, $12, $13, $9}' ${MA}NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt > ${MA}dys_aysnch_NGWAMA_results_4gctb.txt

# Dyslexia
#awk '{print $19, $24, $25, $20, $4, $5, $3, $7+$8}' ${rawss}dyslexia.filtered.2.dat > ${rawss}dyslexia.filtered.2_tmp.dat
# atm freq column is dose.b, but GCTB needs A1 freq, thus convert to dose.a
#awk '{$4 = 1-$4} 1' ${rawss}dyslexia.filtered.2_tmp.dat > ${rawss}dyslexia.filtered.2_tmp2.dat
#echo -e "SNPID A1 A2 A1.freq BETA SE PVAL N" | cat - ${rawss}dyslexia.filtered.2_tmp2.dat > ${rawss}dyslexia.filtered.2_4gctb.ma
#rm ${rawss}dyslexia.filtered.2_tmp.dat
#rm ${rawss}dyslexia.filtered.2_tmp2.dat

# Rhythm
#awk '{print $42, $38, $39, $17, $4, $5, $3, $43}' ${rawss}clap_to_beat.merged.1kgenomes_02222019_wtotalN.txt > ${rawss}clap_to_beat.merged.1kgenomes_02222019_tmp.txt
#flip effect to get rhythm impairment
#awk 'NR>1 {$5 = $5*-1} 1' ${rawss}clap_to_beat.merged.1kgenomes_02222019_tmp.txt > ${rawss}rhythm_impairment_4gctb.ma
#sed -Ei '1s/snp.id/SNPID/;1s/ref_allele/A1/;1s/effect_allele/A2/;1s/freq.a/A1.freq/;1s/effect/BETA/;1s/stderr/SE/;1s/pvalue/PVAL/;1s/0/N/' ${rawss}rhythm_impairment_4gctb.ma

# Reformat Genomic SEM sumstat
awk '{print $1, $5, $6, $4, $12, $13, $15}' ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp.tab > ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_4gctb.ma
#awk '{print $1, $5, $6, $4, $12, $13, $15, $22}' ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2.tab > ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2_4gctb.ma

# Replace Genomic SEM's effective pop size with per-SNP N info from GWAMA sumstats
# Take SNP order into account though!!! Genomic SEM doesn't output per SNP N....

awk 'NR==FNR{a[$1];next} $1 in a {print $1, $8}' ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp.tab ${genSEM}../MA/dys_aysnch_NGWAMA_results4GCTB.ma > ${genSEM}perSNP_Ns_4genomicSEM.txt
sed -Ei '1s/^/SNP N\n/' ${genSEM}perSNP_Ns_4genomicSEM.txt
# add perSNP Ns to gSEM_4gctb.ma
awk -F" " 'NR==FNR{a[NR]=$2;next} {print $0, a[FNR]}' ${genSEM}perSNP_Ns_4genomicSEM.txt ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_4gctb.ma > tmp && mv tmp ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_4gctb.ma
#awk -F" " 'NR==FNR{a[NR]=$2;next} {$8=a[FNR]} 1' OFS=" " ${genSEM}perSNP_Ns_4genomicSEM.txt ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2_4gctb.ma > ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2_4gctb.ma
