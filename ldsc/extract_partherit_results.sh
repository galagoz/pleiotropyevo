#!/bin/bash
#
# Gokberk Alagoz - 02.03.2023
#
# This script will extract partherit results from all LDSC partitioned
# heritability analysis .results file and annotate each result's
# annotations and trait.
#
########################
# PATHS
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/"

for i in ${inDir}*.results; do

	tmp1=$(basename "$i" .results)
	tmp_trait=$(cut -d'.' -f1- <<<"${tmp1}") #TODO this get the whole file name, not only trait name, fix it!
	tmp_annot=$(cut -d'.' -f2- <<<"${tmp1}")

	echo | awk -v trait="$tmp_trait" -v annot="$tmp_annot" '{if(NR==2) print $0, trait, annot}' ${i} > ${inDir}extracted/${tmp1}_firstLine.txt

done

echo -e "Category        Prop._SNPs      Prop._h2        Prop._h2_std_error      Enrichment      Enrichment_std_error    Enrichment_p    Coefficient     Coefficient_std_error   Coefficient_z-score	trait	annot" > ${inDir}all_results.tab
cat ${inDir}extracted/*firstLine.txt >> ${inDir}all_results.tab
