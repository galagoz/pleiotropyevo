#!/bin/bash
#
# This script will extract phyloP and phastCons
# scores of SNPs in the genSEM CPM results.
#
# Gokberk Alagoz, 07.12.22
#

module load bedtools/2.29.2

##########################
# PATHS

input="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/GenomicSEM_multivarGWAS_CPM_dys_rhyimp.sorted.bed"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/"
#declare -a arr=("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phyloP/" "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phastCons/")

#--------------------------------------------
# FUNCTIONS

# This function intersects annotation and 
# bed files using bedtools.

#mkdir ${outDir}scripts
#
#for annot_folder in "${arr[@]}"; do
#    
#    annot_name="$(cut -d'/' -f9 <<<"$annot_folder")"
#    cd $annot_folder
#
#    for i in {1..22} X; do
#
#       annot="chr${i}.${annot_name}46way.primate.starch"
#       script="${outDir}scripts/${annot_name}_${i}.sh"
#       logFile="${outDir}scripts/${annot_name}_${i}.log"
#
#       echo '#$ -N bedops_perChr
##$ -cwd
##$ -q single15.q
##$ -S /bin/bash
#
#bedops --chrom chr'${i}' --element-of 1 '${annot_folder}''${annot}' '${input}' > '${outDir}'chr'${i}''${annot_name}'.intersect.bed' > ${script}
#
#	   chmod a+x ${script}
#	   echo "Created the script for cluster ->  submitting chromosome ${i} to the Grid"
#	   qsub -o ${logFile} -j y ${script};
#    
#   done;
#    
#done

# These commands will...
# 1) Merge all phyloP or PhastCons bed files.
#	2) Sort them.
#	3) Merge these files with the genomic SEM results
#	and create a master file with all info.

# Merge all phyloP
#cat ${outDir}*phyloP* > ${outDir}all_phyloP.bed
#sort -V -k1,1 -k2,2 ${outDir}all_phyloP.bed > ${outDir}all_phyloP.sorted.bed

# Merge all phastCons
#cat ${outDir}*phastCons* > ${outDir}all_phastCons.bed
#sort -V -k1,1 -k2,2 ${outDir}all_phastCons.bed > ${outDir}all_phastCons.sorted.bed

# Merge sorted phyloP and phastCons files
#join -j 2 -o 1.1,1.2,1.3,2.5,1.5 <(sort -V -k1,1 -k2,2 ${outDir}all_phyloP.sorted.bed) <(sort -V -k1,1 -k2,2 ${outDir}all_phastCons.sorted.bed) > ${outDir}all_phastCons_and_phyloP.bed

# Add evo measures to the GenomicSEM CPM results
#join -j 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.4,2.5 <(sort -V -k1,1 -k2,2 ${outDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp.sorted.bed) <(sort -V -k1,1 -k2,2 ${outDir}all_phastCons_and_phyloP.bed) > ${outDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_forEVO.txt

# Add evo measures to the GenomicSEM IPM results - dys
#join -j 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.4,2.5 <(sort -V -k1,1 -k2,2 ${outDir}GenomicSEM_multivarGWAS_dys_IPM.sorted.bed) <(sort -V -k1,1 -k2,2 ${outDir}all_phastCons_and_phyloP.bed) > ${outDir}GenomicSEM_multivarGWAS_dys_IPM_forEVO.txt

# Add evo measures to the GenomicSEM IPM results - rhyimp
join -j 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.4,2.5 <(sort -V -k1,1 -k2,2 ${outDir}GenomicSEM_multivarGWAS_rhyimp_IPM.sorted.bed) <(sort -V -k1,1 -k2,2 ${outDir}all_phastCons_and_phyloP.bed) > ${outDir}GenomicSEM_multivarGWAS_rhyimp_IPM_forEVO.txt
