#$ -N partherit
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/ldsc/partherit_baseline_and_extraAnnot.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/munged/GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun.sumstats.gz fetal_hge_hg19.merged.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun.fetalHGE.results E081_active_marks
