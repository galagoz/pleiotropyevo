#$ -N ldsc_snph2_estimation
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

# Gokberk Alagoz, March 2022 - last updated April 2022

#-----
#PATHS

ldscDir="/home/gokala/programs/ldsc"
mungedDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged"
ldscores="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/eur_w_ld_chr/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor"

#-----

${ldscDir}/ldsc.py \
--h2 ${mungedDir}/flippedZRhythym_dyslexia_MA.sumstats.gz \
--ref-ld-chr ${ldscores} \
--pop-prev 0.0488 \
--samp-prev 0.065 \
--out ${outDir}/flippedZRhythym_dyslexia_popPrev0.0488_sampPrev0.065_heritability \
--w-ld-chr ${ldscores}
