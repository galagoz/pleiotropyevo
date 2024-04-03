#!/bin/bash
#
# This script will append all TWAS results files (from
# each tissue) into one file.
#
# Gokberk Alagoz - 10.04.2023

#############################
# PATHS
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/"

#############################

# add tissue name to each file
#for i in ${inDir}*_assoc.txt; do
#
#	base_name=$(basename "$i")
#	tissue_name="$(cut -d'_' -f2,3 <<<"$base_name")"
#
#	awk -v tis_name=$tissue_name -v OFS=',' 'NR>1 {print $0, tis_name}' ${i} > tmp && mv tmp ${i}
#	echo -e "gene,gene_name,zscore,effect_size,pvalue,var_g,pred_perf_r2,pred_perf_pval,pred_perf_qval,n_snps_used,n_snps_in_cov,n_snps_in_model,tissue"| cat - ${i} > tmp && mv tmp ${i}
#
#done

# append all results files
cat ${inDir}*.txt >> ${inDir}all_tissues.csv
