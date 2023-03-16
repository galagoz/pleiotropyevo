#!/bin/bash
#
# This script will extract LDSC genetic correlation
# results from log files.
#
# Gokberk Alagoz - 09.02.2023
#
#############################

munged="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/munged/"
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/external_sumstats/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/external_sumstats/extracted/"

munged2="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/"
inDir2="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/CPM_IPM_vs_49traits/"
outDir2="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/CPM_IPM_vs_49traits/extracted/"
finalSelection="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gencor/final_selection.txt"

# Extract LDSC gencor results of all external sumstats
# vs. all external sumstats

#for i in ${munged}*gz; do
#
#	tmp1_1=$(basename "$i" .sumstats.gz)
#        tmp1_2="$(cut -d'_' -f1- <<<"${tmp1_1}")"
#
#	for j in ${munged}*gz; do
#		
#		tmp2_1=$(basename "$j" .sumstats.gz)
#                tmp2_2="$(cut -d'_' -f1- <<<"${tmp2_1}")"
#
#		echo | awk -v pheno1="$tmp1_2" -v pheno2="$tmp2_2" '{if(NR==62) print pheno1, pheno2, $3, $4, $6}' ${inDir}${tmp1_2}_${tmp2_2}.log > ${outDir}${tmp1_2}_${tmp2_2}.txt
#
#	done
#done
#
## merge all results in one file
#echo -e "T1 T2 Rg SE P" > ${outDir}all_gencors.tab
#cat ${outDir}*.txt >> ${outDir}all_gencors.tab

###################################
# Extract LDSC gencor results of Common and Independent
# factors vs. 49 selected external traits.

while read i; do
 	
	for j in ${munged2}/GenomicSEM*sumstats*; do
        
	tmp2_1=$(basename "$j" .sumstats.gz)
        tmp2_2="$(cut -d'_' -f1- <<<"${tmp2_1}")"
        echo | awk -v pheno1="$i" -v pheno2="$tmp2_2" '{if(NR==62) print pheno1, pheno2, $3, $4, $6}' ${inDir2}${i}_${tmp2_2}.log > ${outDir2}${i}_${tmp2_2}.txt

       done

done < $finalSelection

echo -e "T1 T2 Rg SE P" > ${outDir2}all_gencors.tab
cat ${outDir2}*.txt >> ${outDir2}all_gencors.tab
