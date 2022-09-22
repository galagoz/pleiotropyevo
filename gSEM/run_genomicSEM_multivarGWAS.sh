#$ -N genomicSEM
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

# This script will run the genomicSEM_multivarGWAS.R script
# see https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects

#-----
module purge
module load R/R-4.0.3
Rscript /data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/Dyslexia/dyslexiaevol/clap-beat_dyslexia/gSEM/genomicSEM_multivarGWAS.R 
