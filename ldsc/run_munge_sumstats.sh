#$ -N munge_sumstats
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

# Gokberk Alagoz, March 2022 - updated September 2022

#-----Munge Sumstats-----

# Reformat GWAS summary stats before computing LDSC intercept
# munge_sumstats.py is from github.com/bulik/ldsc
# Based on the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# After final updates on the grid, it's easier to run munge on lux13. Follow the steps 
# below and then run munge. ALSO, don't forget to reformat your sumstats before munging,
# get rid off all the unnecessary columns to avoid weird errors.

module purge
module load miniconda/3.2021.10 ldsc/v1.0.1
conda activate ldsc

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
munge_sumstats.py \
	--sumstats GenomicSEM_multivarGWAS_dys_rhyimp_v2_forMunge.tab \
	--signed-sumstats Z_Estimate,0 \
	--p Pval_Estimate \
	--frq MAF \
	--out GenomicSEM_multivarGWAS_dys_rhyimp_v2 \
	--merge-alleles /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/w_hm3.snplist \
	--chunksize 500000

