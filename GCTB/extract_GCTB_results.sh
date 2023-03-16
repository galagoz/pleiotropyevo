#!/bin/bash
#
# Gokberk Alagoz - 02.03.2023
#
# This script will extract GCTB results from all SBayesS
# analysis .parRes files and annotate each result's trait.
#
########################
# PATHS
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/GCTB/"

for i in ${inDir}*.parRes; do

        tmp1=$(basename "$i" .results)
        tmp_trait=$(cut -d'.' -f1- <<<"${tmp1}") #TODO this get the whole file name, not only trait name, fix it!

        echo | awk -v trait="$tmp_trait" '{if(NR==7) print $0, trait}' ${i} > ${inDir}extracted/${tmp1}_SCoefficientLine.txt

done

echo -e "Measure    Mean   SD	trait" > ${inDir}all_results.tab
cat ${inDir}extracted/*SCoefficientLine.txt >> ${inDir}all_results.tab
