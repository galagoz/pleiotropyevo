#$ -N ldsc_snph2_estimation
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

# Gokberk Alagoz, March 2022 - last updated October 2022

#-----
#PATHS

ldscDir="/home/gokala/programs/ldsc/"
mungedDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/munged/"
ldscores="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/eur_w_ld_chr/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/"

#-----
module purge
module load miniconda/3.2021.10 ldsc/v1.0.1
source activate ldsc
#-----
# Run plain
${ldscDir}ldsc.py \
--h2 ${mungedDir}GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun.sumstats.gz \
--ref-ld-chr ${ldscores} \
--out ${outDir}GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun_heritability \
--w-ld-chr ${ldscores}

#-----
# Run with liability scale
#${ldscDir}ldsc.py \
#--h2 ${mungedDir}GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun.sumstats.gz \
#--ref-ld-chr ${ldscores} \
#--pop-prev 0.0488 \
#--samp-prev 0.065 \
#--out ${outDir}GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun_popPrev0.0488_sampPrev0.065_heritability \
#--w-ld-chr ${ldscores}
