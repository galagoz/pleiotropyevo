#$ -N partherit
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

/data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/Dyslexia/dyslexiaevol/clap-beat_dyslexia/ldsc/partherit_baseline_and_extraAnnot.sh /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/flippedZRhythym.sumstats.gz fetal_hge_hg19.merged.sorted /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/partherit/flippedZRhythym.fetalHGE.results E081_active_marks
