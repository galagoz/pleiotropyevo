#$ -N munge_sumstats
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

# Gokberk Alagoz, March 2022

#-----Munge Sumstats-----

# Reformat GWAS summary stats before computing LDSC intercept
# munge_sumstats.py is from github.com/bulik/ldsc
# Based on the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

#-----Variables-----
# $input - summary statistic file
# $output - outfile_name

ldsc="/home/gokala/programs/ldsc"
hapmap="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist"
inDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM"
sumstatsList="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/sumstats_toMunge.txt"
MA="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA"

# mkdir ${inDir}/munged

#-----Rhythym

#${ldsc}/munge_sumstats.py \
#--sumstats ${inDir}/rhythym_reformatted_forNweighted6.txt \
#--signed-sumstats Z,0 \
#--out ${inDir}/munged/flippedZRhythym \
#--merge-alleles ${hapmap}

#-----Dyslexia

# ${ldsc}/munge_sumstats.py \
# --sumstats ${inDir}/dyslexia_reformatted_forModelAveraged_forMunge.tab \
# --signed-sumstats Z,0 \
# --out ${inDir}/munged/dyslexia \
# --merge-alleles ${hapmap}

#-----GWAMA

# ${ldsc}/munge_sumstats.py \
# --sumstats ${MA}/NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt \
# --frq MAF \
# --N-col N_obs \
# --out ${inDir}/munged/flippedZRhythym_dyslexia_MA \
# --merge-alleles ${hapmap}

#-----Genomic SEM

${ldsc}/munge_sumstats.py \
--sumstats ${inDir}/GenomicSEM_multivarGWAS_dys_rhyimp_v2.tab \
--frq MAF \
--signed-sumstats Z_Estimate,0 \
--out ${inDir}/munged/GenomicSEM_multivarGWAS_dys_rhyimp_v2.tab \
--merge-alleles ${hapmap}

/home/gokala/programs/ldsc/munge_sumstats.py --sumstats GenomicSEM_multivarGWAS_dys_rhyimp_v2.tab --signed-sumstats Z_Estimate,0 --p Pval_Estimate --frq MAF --out GenomicSEM_multivarGWAS_dys_rhyimp_v2 --merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist
