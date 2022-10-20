#$ -N gencor
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

# Gokberk Alagoz, March 2022

# -----
# PATHS

ldscDir="/home/gokala/programs/ldsc"
mungedDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged"
ldscores="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/eur_w_ld_chr/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor"

# -----
module purge
module load miniconda/3.2021.10 ldsc/v1.0.1
source activate ldsc
# -----
# dys rhy_imp gencor

${ldscDir}/ldsc.py \
 --rg ${mungedDir}/dyslexia.sumstats.gz,${mungedDir}/flippedZRhythym.sumstats.gz \
 --ref-ld-chr ${ldscores} \
 --w-ld-chr ${ldscores} \
 --out ${outDir}/dys_rhyimp

# -----
# self gencors for LAVA

#${ldscDir}/ldsc.py \
#--rg ${mungedDir}/dyslexia.sumstats.gz,${mungedDir}/dyslexia.sumstats.gz \
#--ref-ld-chr ${ldscores} \
#--w-ld-chr ${ldscores} \
#--out ${outDir}/dys_self

#${ldscDir}/ldsc.py \
#--rg ${mungedDir}/flippedZRhythym.sumstats.gz,${mungedDir}/flippedZRhythym.sumstats.gz \
#--ref-ld-chr ${ldscores} \
#--w-ld-chr ${ldscores} \
#--out ${outDir}/rhyimp_self
