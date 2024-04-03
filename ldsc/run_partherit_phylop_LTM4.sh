#$ -N partherit
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

# Dyslexia

/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/dyslexia.sumstats.gz allchr.phyloP46way.primate_lessThanMinus4.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/dyslexia.phylop_LTM4

# Rhythm_impairment

#/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline_and_extraAnnot.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/rhythm_impairment.sumstats.gz allchr.phyloP46way.primate_lessThanMinus4.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/rhythm_impairment.phylop_LTM4 allchr.phyloP46way.primate.sorted_control.bed

# Genomic SEM common factor results

#/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline_and_extraAnnot.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/GenomicSEM_multivarGWAS_CPM_dys_rhyimp.sumstats.gz allchr.phyloP46way.primate_lessThanMinus4.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/GenomicSEM_multivarGWAS_dys_rhyimp.phylop_LTM4 allchr.phyloP46way.primate.sorted_control.bed

# Genomic SEM IPM results - dyslexia

#/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline_and_extraAnnot.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/GenomicSEM_multivarGWAS_dys_IPM.sumstats.gz allchr.phyloP46way.primate_lessThanMinus4.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/GenomicSEM_multivarGWAS_dys_IPM.phylop_LTM4 allchr.phyloP46way.primate.sorted_control.bed

# Genomic SEM IPM results - rhythm impairment

#/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline_and_extraAnnot.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/GenomicSEM_multivarGWAS_rhyimp_IPM.sumstats.gz allchr.phyloP46way.primate_lessThanMinus4.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/GenomicSEM_multivarGWAS_rhyimp_IPM.phylop_LTM4 allchr.phyloP46way.primate.sorted_control.bed
