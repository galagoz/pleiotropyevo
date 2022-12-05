#!/bin/bash
#
# This script will run LAVA for all relevant annotations.
#
# Gokberk Alagoz, 29 June 2022 - last updated Nov 2022

inputDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/LAVA/input/"
refpanel_1kgp3="/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_allchr"
bedDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/new_annotations/hg19_beds/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/LAVA/"
inputInfo="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/LAVA/input/input.info.txt"
sampleOverlap="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/LAVA/input/sample.overlap.txt"
phenos="dyslexia-rhythm_impairment"

mkdir ${outDir}/scripts

for annot in ${bedDir}*.bed; do

	tmp_annot=$(basename "$annot")
	echo $tmp_annot
	tmp_run_file="${outDir}/scripts/${tmp_annot}.sh"

	echo '#$ -N LAVA_local_rg
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

module load R/R-4.0.3

Rscript /data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/LAVA/lava_script.R' $refpanel_1kgp3 $annot $inputInfo $sampleOverlap $phenos $annot > $tmp_run_file

	chmod a+x $tmp_run_file
	echo "Created the script for cluster -> submitting ${annot} to the Grid"
	qsub -wd ${outDir}/scripts $tmp_run_file;

done
