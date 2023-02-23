#!/bin/bash
#
# This script will extract LDSC genetic correlation
# results from log files.
#
# Gokberk Alagoz - 09.02.2023
#
#############################
#
#munged="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/munged/for_rhythm/"
#inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/rhythm_gencors/"
#outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/rhythm_gencors/extracted/"
#
#for j in ${munged}*gz; do
#
#                tmp2_1=$(basename "$j" .sumstats.gz)
#                tmp2_2="$(cut -d'_' -f1- <<<"${tmp2_1}")"
#
#                echo | awk -v pheno2="$tmp2_2" '{if(NR==62) print "rhythm", pheno2, $3, $4, $6}' ${inDir}rhy_${tmp2_2}.log > ${outDir}rhy_${tmp2_2}.txt
#
#done
#
#echo -e "T1 T2 Rg SE P" > ${outDir}all_gencors.tab
#cat ${outDir}*.txt >> ${outDir}all_gencors.tab

munged="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/munged/"
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/external_sumstats/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/external_sumstats/extracted/"

for i in ${munged}*gz; do

	tmp1_1=$(basename "$i" .sumstats.gz)
        tmp1_2="$(cut -d'_' -f1- <<<"${tmp1_1}")"

	for j in ${munged}*gz; do
		
		tmp2_1=$(basename "$j" .sumstats.gz)
                tmp2_2="$(cut -d'_' -f1- <<<"${tmp2_1}")"

		echo | awk -v pheno1="$tmp1_2" -v pheno2="$tmp2_2" '{if(NR==62) print pheno1, pheno2, $3, $4, $6}' ${inDir}${tmp1_2}_${tmp2_2}.log > ${outDir}${tmp1_2}_${tmp2_2}.txt

	done
done

# merge all results in one file

echo -e "T1 T2 Rg SE P" > ${outDir}all_gencors_v2.tab
cat ${outDir}*.txt >> ${outDir}all_gencors_v2.tab
