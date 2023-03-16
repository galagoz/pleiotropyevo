#!/bin/bash
#
# Gokberk Alagoz - 13.03.2023
#
# This script will convert .starch files to .bed format
# and merge all of chromosome in one file.
#
##############################
# PATHS
inDir_phylop="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phyloP/"
inDir_phastcons="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phastCons/"

module load bedtools/2.29.2

#for i in ${inDir_phylop}*; do
#
#	tmp_name=$(cut -d'.' -f1,2,3- <<<"$i")
#	unstarch ${i} > ${tmp_name}.bed;
#
#done

for i in ${inDir_phastcons}*; do

        tmp_name=$(cut -d'.' -f1,2,3- <<<"$i")
        unstarch ${i} > ${tmp_name}.bed;

done
