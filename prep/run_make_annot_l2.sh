#!/bin/bash

# Make .annot and l2.ldscore files for all chromosomes and annotations.
# Run this on grid as "bash run_make_annot_l2.sh"

# Variables
inDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/beds"

# Run LD Score Estimation

mkdir $inDir/scripts

for annot in ${inDir}/*sorted.bed; do

	echo $annot
	tmp_annot=$(basename "$annot")
	echo $tmp_annot
	tmp_run_file="${inDir}/scripts/${tmp_annot}.sh"
	echo '#$ -N make_annot_and_ldscore
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

export PATH="/home/gokala/programs/bedtools2/bin:$PATH"
module purge
module load miniconda/3.2021.10 ldsc/v1.0.1 
conda activate ldsc
module purge

mkdir '$inDir'/'$tmp_annot'
bash /data/workspaces/lag/workspaces/lg-genlang/Working/SCRIPTS/23andMe/clap-beat_dyslexia/pleiotropyevo/prep/make_annot_l2.sh ' $annot > $tmp_run_file

	chmod a+x $tmp_run_file
   	echo "Created the script for cluster -> submitting ${annot} to the Grid"
   	qsub -wd "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phyloP/scripts" $tmp_run_file

done
