#!/bin/bash
#
# Gokberk Alagoz, March 2022 - last update Feb 2023

#-----
# PATHS

mungedDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged"
ldscores="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/eur_w_ld_chr/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor"
for_rhythm="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/munged/for_rhythm/"
external="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/munged/"
finalSelection="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gencor/final_selection.txt"

#-----
module purge
module load miniconda/3.2021.10 ldsc/v1.0.1
conda activate ldsc
#-----
# dys vs. rhy_imp gencor

#ldsc.py \
# --rg ${mungedDir}/dyslexia.sumstats.gz,${mungedDir}/rhythm_impairment.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --w-ld-chr ${ldscores} \
# --out ${outDir}/dys_rhyimp
#
##-----
## self gencors for LAVA
#
#ldsc.py \
# --rg ${mungedDir}/dyslexia.sumstats.gz,${mungedDir}/dyslexia.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --w-ld-chr ${ldscores} \
# --out ${outDir}/dys_self
#
#ldsc.py \
# --rg ${mungedDir}/rhythm_impairment.sumstats.gz,${mungedDir}/rhythm_impairment.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --w-ld-chr ${ldscores} \
# --out ${outDir}/rhyimp_self

#-----
# genetic correlations between rhythm and a small selection
# of language-related traits.

#for i in ${for_rhythm}*.gz; do
#	
#	tmp_1=$(basename "$i" .sumstats.gz)
#	tmp_2="$(cut -d'.' -f1- <<<"${tmp_1}")"
#
#	ldsc.py \
#		--rg ${for_rhythm}Clapbeat_02222019_forLDSC.sumstats.gz,${i} \
#		--ref-ld-chr ${ldscores} \
#		--w-ld-chr ${ldscores} \
#		--out ${outDir}/rhythm_gencors/rhy_${tmp_2}
#
#done

#-----
# genetic correlations among all external sumstats

#for i in ${external}*.gz; do
#
#	tmp1_1=$(basename "$i" .sumstats.gz)
#	tmp1_2="$(cut -d'_' -f1- <<<"${tmp1_1}")"
#
#	for j in ${external}*.gz; do
#
#		tmp2_1=$(basename "$j" .sumstats.gz)
#		tmp2_2="$(cut -d'_' -f1- <<<"${tmp2_1}")"
#
#		ldsc.py \
#			--rg ${i},${j} \
#			--ref-ld-chr ${ldscores} \
#			--w-ld-chr ${ldscores} \
#			--out ${outDir}/external_sumstats/${tmp1_2}_${tmp2_2}
#
#	done
#done

#-----
# genetic correlations between Commmon and Independent 
# factors vs. selected 49 external traits.

while read i; do
       for j in ${mungedDir}/GenomicSEM*sumstats*; do

               tmp2_1=$(basename "$j" .sumstats.gz)
               tmp2_2="$(cut -d'_' -f1- <<<"${tmp2_1}")"

               ldsc.py \
                       --rg ${external}${i}.sumstats.gz,${j} \
                       --ref-ld-chr ${ldscores} \
                       --w-ld-chr ${ldscores} \
                       --out ${outDir}/CPM_IPM_vs_49traits/${i}_${tmp2_2}
       done
done < $finalSelection
