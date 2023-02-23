#!/bin/bash
#
# This script will download UKB sumstats when provided
# a list of UKB phenotype codes from Neale group's 
# data repo: UKBB GWAS Imputed v3 - File Manifest Release 20180731
#
# Gokberk Alagoz - 02.02.2023
#
#############################
phenoList="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/UKB_traits/pheno_code_list.txt"
outDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/UKB_traits/"

while read pheno; do

	echo $pheno
	wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/${pheno}.gwas.imputed_v3.both_sexes.tsv.bgz -O ${outDir}${pheno}.gwas.imputed_v3.both_sexes.tsv.gz

done < $phenoList
