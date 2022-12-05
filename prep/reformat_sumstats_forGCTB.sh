#!/bin/bash
#
# This script will reformat GWAS sumstats for
# GCTB analysis.

# Paths
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"
genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"

# 1) Dyslexia
# dose.b column is allele frequency of alleleB, which is the effect allele.

#awk 'NR>1 {print $19, $24, $25, $20, $4, $5, $3, $7+$8}' ${inDir}dyslexia.filtered.2_lambdaGCcorrected.dat > ${inDir}dyslexia.filtered.2_lambdaGCcorrected_forGCTB1.dat
#echo -e "SNPID A1 A2 A1.freq BETA SE PVAL N" | cat - ${inDir}dyslexia.filtered.2_lambdaGCcorrected_forGCTB1.dat > ${inDir}dyslexia_forGCTB.ma

#chmod 777 ${inDir}dyslexia_forGCTB.ma
#rm ${inDir}dyslexia.filtered.2_lambdaGCcorrected_forGCTB1.dat

# 2) Rhythm impairment

#awk '{print $42, $38, $39, $17, $4, $5, $3, $43}' ${inDir}rhythm_impairment_lambdaGCcorrected_wtotalNandOR.dat > ${inDir}rhythm_impairment_forGCTB.ma
#sed -Ei '1s/snp.id/SNPID/;1s/ref_allele/A1/;1s/effect_allele/A2/;1s/freq.a/A1.freq/;1s/effect/BETA/;1s/stderr/SE/;1s/pvalue/PVAL/;1s/0/N/' ${inDir}rhythm_impairment_forGCTB.ma

#chmod 777 ${inDir}rhythm_impairment_forGCTB.ma

# 3) Genomic SEM common factor results

awk 'NR==FNR{a[$1];next} $2 in a {print $2, $10}' ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_lambdaGCcorrected.tab ${genSEM}../GWAMA/dyslexia_rhythmImpairment.N_weighted_GWAMA.results.txt > ${genSEM}perSNP_Ns_forgenomicSEM.txt
sed -Ei '1s/^/SNP N\n/' ${genSEM}perSNP_Ns_forgenomicSEM.txt

awk '{print $1, $5, $6, $4, $12, $13, $15}' ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_lambdaGCcorrected.tab > ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_forGCTB1.ma
awk -F" " 'NR==FNR{a[NR]=$2;next} {$8=a[FNR]} 1' OFS=" " ${genSEM}perSNP_Ns_forgenomicSEM.txt ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_forGCTB1.ma > ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_forGCTB.ma

chmod 777 ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_forGCTB.ma
rm ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_forGCTB1.ma
