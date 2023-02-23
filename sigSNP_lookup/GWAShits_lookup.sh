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
MAsigSNPs="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/FUMA/FUMA_genSEM_dys_rhyimp/GenomicRiskLoci.txt"
#DysSNPs="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/FUMA/FUMA_dys/GenomicRiskLoci.txt"
#RhyImpSNPs="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/FUMA/FUMA_rhyimp/GenomicRiskLoci.txt"

#MA_ss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/dys_aysnch_NGWAMA_results4GCTB.ma"
MA_ss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_dys_rhyimp_forGCTB.ma"
RhyImp_ss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythm_impairment_forGCTB.ma"
Dys_ss="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia_forGCTB.ma"

outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/snp_lookup/"

###########################
# Merge Genomic Risk Loci lists
# !!! Changed design here: instead of merging all
# risk loci, use only genSEM risk loci.

#cp $MAsigSNPs $outDir
#mv ${outDir}GenomicRiskLoci.txt ${outDir}all_GenomicRiskLoci.txt

#head -n 1 $MAsigSNPs > ${outDir}all_GenomicRiskLoci.txt
#awk FNR!=1 ${MAsigSNPs} ${DysSNPs} ${RhyImpSNPs} > ${outDir}all_GenomicRiskLoci_content.txt
#cat ${outDir}all_GenomicRiskLoci_content.txt >> ${outDir}all_GenomicRiskLoci.txt 

#rm ${outDir}all_GenomicRiskLoci_content.txt

#awk 'NR>1 {print $3}' ${outDir}all_GenomicRiskLoci.txt > ${outDir}all_GenomicRiskLoci_SNPlist.txt
#sort ${outDir}all_GenomicRiskLoci_SNPlist.txt | uniq -u > ${outDir}all_GenomicRiskLoci_SNPlist_unique.txt

############################
# Get CHR, BP, Effect, SE, P per SNP from each summary statistics file.

#grep -iFwf ${outDir}all_GenomicRiskLoci_SNPlist_unique.txt ${MA_ss} > ${outDir}all_GenomicRiskLoci_MAallContent.txt
#awk 'NR>1 {print $1, $5, $6, $7}' ${outDir}all_GenomicRiskLoci_MAallContent.txt > ${outDir}all_GenomicRiskLoci_MAcontent.txt
#rm ${outDir}all_GenomicRiskLoci_MAallContent.txt
#awk '{print $0"\t""MA"}' ${outDir}all_GenomicRiskLoci_MAcontent.txt > ${outDir}tmp && mv ${outDir}tmp ${outDir}all_GenomicRiskLoci_MAcontent.txt

#grep -iFwf ${outDir}all_GenomicRiskLoci_SNPlist_unique.txt ${RhyImp_ss} > ${outDir}all_GenomicRiskLoci_RhyImpallContent.txt
#awk 'NR>1 {print $1, $5, $6, $7}' ${outDir}all_GenomicRiskLoci_RhyImpallContent.txt > ${outDir}all_GenomicRiskLoci_RhyImpcontent.txt
#rm ${outDir}all_GenomicRiskLoci_RhyImpallContent.txt
#awk '{print $0"\t""RhyImp"}' ${outDir}all_GenomicRiskLoci_RhyImpcontent.txt > ${outDir}tmp && mv ${outDir}tmp ${outDir}all_GenomicRiskLoci_RhyImpcontent.txt

#grep -iFwf ${outDir}all_GenomicRiskLoci_SNPlist_unique.txt ${Dys_ss} > ${outDir}all_GenomicRiskLoci_DysallContent.txt
#awk 'NR>1 {print $1, $5, $6, $7}' ${outDir}all_GenomicRiskLoci_DysallContent.txt > ${outDir}all_GenomicRiskLoci_Dyscontent.txt
#rm ${outDir}all_GenomicRiskLoci_DysallContent.txt
#awk '{print $0"\t""Dys"}' ${outDir}all_GenomicRiskLoci_Dyscontent.txt > ${outDir}tmp && mv ${outDir}tmp ${outDir}all_GenomicRiskLoci_Dyscontent.txt

###########################
# Make sure you have same SNPs in each of 3 files you've made.

echo "SNP Effect SE P Trait" > ${outDir}GenomicRiskLoci_allTraits.txt
cat ${outDir}all_GenomicRiskLoci_MAcontent.txt ${outDir}all_GenomicRiskLoci_RhyImpcontent.txt ${outDir}all_GenomicRiskLoci_Dyscontent.txt >> ${outDir}GenomicRiskLoci_allTraits.txt
