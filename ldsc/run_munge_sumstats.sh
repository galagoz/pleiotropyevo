#!/bin/bash
#
# Gokberk Alagoz, March 2022 - updated Nov 2022

#-----Munge Sumstats-----

# Reformat GWAS summary stats before computing LDSC intercept
# munge_sumstats.py is from github.com/bulik/ldsc
# Based on the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# After final updates on the grid, it's easier to run munge on lux13. Follow the steps 
# below and then run munge.

module purge
module load miniconda/3.2021.10 ldsc/v1.0.1
conda activate ldsc

#-----Variables-----
hapmap="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist"
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data"
genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM"
external_pgc="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/various_traits/"
external_ukb="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/UKB_traits/"

# 1) Rhythym impairment

#munge_sumstats.py \
# --sumstats ${inDir}/rhythm_impairment_raw_sumstat_wtotalNandOR.txt \
# --snp snp.id \
# --a1 effect_allele \
# --a2 ref_allele \
# --frq maf \
# --N-cas-col im.num.1 \
# --N-con-col im.num.0 \
# --p pvalue \
# --signed-sumstats OR,1 \
# --merge-alleles ${hapmap} \
# --out ${inDir}/munged/rhythm_impairment \
# --chunksize 500000

# 2) Dyslexia

#munge_sumstats.py \
# --sumstats ${inDir}/dyslexia.filtered.2_wtotalN.txt \
# --snp rsid \
# --a1 alleleB \
# --a2 alleleA \
# --N-cas-col im.num.1 \
# --N-con-col im.num.0 \
# --p pvalue \
# --signed-sumstats OR,1 \/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/various_traits
# --merge-alleles ${hapmap} \
# --out ${inDir}/munged/dyslexia \
# --chunksize 500000

# 3) Genomic SEM common factor results

#munge_sumstats.py \
#	--sumstats ${genSEM}/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_forMunge.tab \
#	--signed-sumstats Z_Estimate,0 \
#	--p Pval_Estimate \
#	--frq MAF \
#	--out ${inDir}/munged/GenomicSEM_multivarGWAS_CPM_dys_rhyimp \
#	--merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
#	--chunksize 500000

# 4) Genomic SEM common factor results - GC corrected

#munge_sumstats.py \
#        --sumstats ${genSEM}/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_forMunge.tab \
#        --signed-sumstats Z_Estimate,0 \
#        --p Pval_Estimate \
#        --frq MAF \
#        --out ${inDir}/munged/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr \
#        --merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
#        --chunksize 500000

# 5) Genomic SEM IPM results - dyslexia

#munge_sumstats.py \
#        --sumstats ${genSEM}/GenomicSEM_multivarGWAS_dys_IPM_forMunge.tab \
#        --signed-sumstats Z_Estimate,0 \
#        --p Pval_Estimate \
#        --frq MAF \
#        --out ${inDir}/munged/GenomicSEM_multivarGWAS_dys_IPM \
#        --merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
#        --chunksize 500000
#
# 6) Genomic SEM IPM results - rhythm impairment

#munge_sumstats.py \
#        --sumstats ${genSEM}/GenomicSEM_multivarGWAS_rhyimp_IPM_forMunge.tab \
#        --signed-sumstats Z_Estimate,0 \
#        --p Pval_Estimate \
#        --frq MAF \
#        --out ${inDir}/munged/GenomicSEM_multivarGWAS_rhyimp_IPM \
#        --merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
#        --chunksize 500000
#
# 7) External summary statistics

#for i in ${external_pgc}*_forMunge.txt; do
#
#	tmp_name=$(basename "$i" _forMunge.txt)
#
#	munge_sumstats.py \
#        	--sumstats ${i} \
#        	--out ${external_pgc}/munged/${tmp_name} \
#	        --merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
#        	--chunksize 500000
#
#done

#for i in ${external_ukb}*_forMunge.txt; do
#
#        tmp_name=$(basename "$i" _forMunge.txt)
#
#        munge_sumstats.py \
#                --sumstats ${i} \
#                --out ${external_ukb}/munged/${tmp_name} \
#                --merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
#                --chunksize 500000
#
#done

while read i; do

	tmp_name=$(basename "$i" .tsv)

        munge_sumstats.py \
               --sumstats ${external_ukb}${tmp_name}_forMunge.txt \
                --out ${external_ukb}munged/${tmp_name} \
                --merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
                --chunksize 500000

done < ${external_ukb}UKBsumstats_list_diffColumnNumbers.txt

