#$ -N partherit
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

# Dyslexia

#/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/dyslexia.sumstats.gz neanDepRegions_hg19.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/dyslexia.neanDepRegions

# Rhythm impairment

#/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/rhythm_impairment.sumstats.gz neanDepRegions_hg19.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/rhythm_impairment.neanDepRegions

# Genomic SEM common factor results

/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/munged/GenomicSEM_multivarGWAS_dys_rhyimp.sumstats.gz neanDepRegions_hg19.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/GenomicSEM_multivarGWAS_dys_rhyimp.neanDepRegions
