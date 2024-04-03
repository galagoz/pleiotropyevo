#$ -N partherit
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

# Genomic SEM common factor results
/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/GenomicSEM_multivarGWAS_CPM_dys_rhyimp.sumstats.gz neuronal_promoters.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/GenomicSEM_multivarGWAS_dys_rhyimp.neuronalPromoters
