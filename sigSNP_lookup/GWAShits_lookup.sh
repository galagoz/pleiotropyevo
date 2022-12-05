#!/bin/bash
#
# Gokberk Alagoz - 01.11.22
#
# This script will do the following:
# 1) Make a list of gw-sig SNPs at least in one of these summary statistics
# files: Genomic SEM CPM meta-analysis results, dyslexia, rhythm impairment.
# 2) Remove recurrent SNPs from this list.
# 3) Get Effect Sizes and SEs per SNPs, calculate Z-score per SNP.
# Then, use this output file to make a heatmap in R showing effect directions
# and Z-scores per GWAS hit for MA, Dyslexia and Rhythm Impairment.

# First download FUMA outputs for each GWAS.

##########################
# Variables

#MAsigSNPs="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/FUMA/old/FUMA_job178902/GenomicRiskLoci.txt"
MAsigSNPs="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/FUMA/FUMA_job214098_genomicSEM_CPM_Dys_RhyImp/GenomicRiskLoci.txt"
DysSNPs="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/FUMA/FUMA_job182150_Dys/GenomicRiskLoci.txt"
RhyImpSNPs="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/FUMA/FUMA_job182374_RhyImp/GenomicRiskLoci.txt"

#MA_ss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/dys_aysnch_NGWAMA_results4GCTB.ma"
MA_ss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_dys_rhyimp_finalModel_4gctb.ma"
RhyImp_ss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythm_impairment_4gctb.ma"
Dys_ss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia.filtered.2_4gctb.ma"

outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/snp_lookup/"

###########################
# Merge Genomic Risk Loci lists

head -n 1 $MAsigSNPs > ${outDir}all_GenomicRiskLoci.txt
awk FNR!=1 ${MAsigSNPs} ${DysSNPs} ${RhyImpSNPs} > ${outDir}all_GenomicRiskLoci_content.txt
cat ${outDir}all_GenomicRiskLoci_content.txt >> ${outDir}all_GenomicRiskLoci.txt 

rm ${outDir}all_GenomicRiskLoci_content.txt

awk 'NR>1 {print $3}' ${outDir}all_GenomicRiskLoci.txt > ${outDir}all_GenomicRiskLoci_SNPlist.txt # has 266 SNPs
sort ${outDir}all_GenomicRiskLoci_SNPlist.txt | uniq -u > ${outDir}all_GenomicRiskLoci_SNPlist_unique.txt # has 242 SNPs

############################
# Get CHR, BP, Effect, SE, P per SNP from each summary statistics file.

grep -iFwf ${outDir}all_GenomicRiskLoci_SNPlist_unique.txt ${MA_ss} > ${outDir}all_GenomicRiskLoci_MAallContent.txt
awk 'NR>1 {print $1, $5, $6, $7}' ${outDir}all_GenomicRiskLoci_MAallContent.txt > ${outDir}all_GenomicRiskLoci_MAcontent.txt
rm ${outDir}all_GenomicRiskLoci_MAallContent.txt
awk '{print $0"\t""MA"}' ${outDir}all_GenomicRiskLoci_MAcontent.txt > ${outDir}tmp && mv ${outDir}tmp ${outDir}all_GenomicRiskLoci_MAcontent.txt
# MA has 238 risk loci

grep -iFwf ${outDir}all_GenomicRiskLoci_SNPlist_unique.txt ${RhyImp_ss} > ${outDir}all_GenomicRiskLoci_RhyImpallContent.txt
awk 'NR>1 {print $1, $5, $6, $7}' ${outDir}all_GenomicRiskLoci_RhyImpallContent.txt > ${outDir}all_GenomicRiskLoci_RhyImpcontent.txt
rm ${outDir}all_GenomicRiskLoci_RhyImpallContent.txt
awk '{print $0"\t""RhyImp"}' ${outDir}all_GenomicRiskLoci_RhyImpcontent.txt > ${outDir}tmp && mv ${outDir}tmp ${outDir}all_GenomicRiskLoci_RhyImpcontent.txt
# RhyImp has 241 risk loci

grep -iFwf ${outDir}all_GenomicRiskLoci_SNPlist_unique.txt ${Dys_ss} > ${outDir}all_GenomicRiskLoci_DysallContent.txt
awk 'NR>1 {print $1, $5, $6, $7}' ${outDir}all_GenomicRiskLoci_DysallContent.txt > ${outDir}all_GenomicRiskLoci_Dyscontent.txt
rm ${outDir}all_GenomicRiskLoci_DysallContent.txt
awk '{print $0"\t""Dys"}' ${outDir}all_GenomicRiskLoci_Dyscontent.txt > ${outDir}tmp && mv ${outDir}tmp ${outDir}all_GenomicRiskLoci_Dyscontent.txt
# Dys has 240 risk loci

###########################
# Make sure you have same SNPs in each of 3 files you've made.

awk 'NR==FNR{a[$1];next} $1 in a {print $0}' ${outDir}all_GenomicRiskLoci_MAcontent.txt ${outDir}all_GenomicRiskLoci_RhyImpcontent.txt > ${outDir}all_GenomicRiskLoci_RhyImpcontent_subset.txt
awk 'NR==FNR{a[$1];next} $1 in a {print $0}' ${outDir}all_GenomicRiskLoci_MAcontent.txt ${outDir}all_GenomicRiskLoci_Dyscontent.txt > ${outDir}all_GenomicRiskLoci_Dyscontent_subset.txt

echo "SNP Effect SE P Trait" > ${outDir}GenomicRiskLoci_allTraits.txt
cat ${outDir}all_GenomicRiskLoci_MAcontent.txt ${outDir}all_GenomicRiskLoci_RhyImpcontent_subset.txt ${outDir}all_GenomicRiskLoci_Dyscontent_subset.txt >> ${outDir}GenomicRiskLoci_allTraits.txt
