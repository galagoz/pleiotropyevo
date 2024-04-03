#!/bin/bash
#
# Gokberk Alagoz - 13.03.2023
#
# This script will convert .starch files to .bed format
# and merge all of chromosome in one file. Then, it will
# reformat all genomic coordinate files into UCSC bed
# format which is necessary for LDSC to run. This format
# looks like this:
# chr1	1	2
# chr1	10	11
# chr1	56	57
#
# If you want to make an annotation using a SNP-list,
# the actual position of the SNP should be the last column,
# and  the second column should be position-1:
# chrX	position-1	position

##############################
# PATHS
inDir_phylop="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phyloP/"
inDir_phastcons="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phastCons/"
inDir_beds="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/beds/"

module load bedtools/2.29.2

##############################
# unstarch phyloP and phastCons
# per chromosome SNP lists.

# phyloP tree model for the 10 primates subset:
# hg19, panTro2, gorGor1, ponAbe2, rheMac2, papHam1, calJac1, tarSyr1,
# micMur1, otoGar1
# (human, chimp, gorilla, orangutan, rhesus, baboon, marmoset, tarsier,
# mouse lemur, bushbaby)

#for i in ${inDir_phylop}*; do
#
#	tmp_name=$(cut -d'.' -f1,2,3- <<<"$i")
#	unstarch ${i} > ${tmp_name}.bed;
#
#done

#for i in ${inDir_phastcons}*; do
#
#        tmp_name=$(cut -d'.' -f1,2,3- <<<"$i")
#        unstarch ${i} > ${tmp_name}.bed;
#
#done

##############################
# merge chromosomes

#cat ${inDir_phylop}*.bed >> ${inDir_phylop}allchr.phyloP46way.primate.bed
#cat ${inDir_phastcons}*.bed >> ${inDir_phastcons}allchr.phastCons46way.primate.bed

##############################
# Filter allchr files based on phyloP
# and phastCons scores.

#awk '{if ($5 <= -2) print $0}' ${inDir_phylop}allchr.phyloP46way.primate.bed > ${inDir_phylop}allchr.phyloP46way.primate_lessThanMinus2.bed
#awk '{if ($5 <= -3) print $0}' ${inDir_phylop}allchr.phyloP46way.primate.bed > ${inDir_phylop}allchr.phyloP46way.primate_lessThanMinus3.bed
#awk '{if ($5 <= -4) print $0}' ${inDir_phylop}allchr.phyloP46way.primate.bed > ${inDir_phylop}allchr.phyloP46way.primate_lessThanMinus4.bed
#awk '{if ($5 <= -5) print $0}' ${inDir_phylop}allchr.phyloP46way.primate.bed > ${inDir_phylop}allchr.phyloP46way.primate_lessThanMinus5.bed

#awk '{if ($5 >= 0.25) print $0}' ${inDir_phastcons}allchr.phastCons46way.primate.bed > ${inDir_phastcons}allchr.phastCons46way.primate_moreThan0.25.bed
#awk '{if ($5 >= 0.50) print $0}' ${inDir_phastcons}allchr.phastCons46way.primate.bed > ${inDir_phastcons}allchr.phastCons46way.primate_moreThan0.50.bed
#awk '{if ($5 >= 0.75) print $0}' ${inDir_phastcons}allchr.phastCons46way.primate.bed > ${inDir_phastcons}allchr.phastCons46way.primate_moreThan0.75.bed

##############################
# Reformat .bed files and sort.

## 1) phyloP primate files.

## Remove 4th and 5th columns as they are not needed
#for i in ${inDir_phylop}*lessThanMinus*; do
#
#	awk '{$4=""; $5=""; print $0}' ${i} > ${inDir_phylop}tmp && mv ${inDir_phylop}tmp ${i};
#	
#done
#
### Sort concatenated bed file based on 1st and 2nd columns
#for j in ${inDir_phylop}*lessThanMinus*; do
#	
#	sort -k1,1 -k2,2n ${j} > ${j%.bed}.sorted.bed;
#
#done
#
## Remove X and Y chromosomes as these .bed files will
## be used by LDSC, and LDSC cannot handle X chromosome
## data.
#for file2 in ${inDir_phylop}*.sorted.bed; do
#
#      awk '$2!="X" && $2!="Y"' ${file2} > ${inDir_phylop}tmp && mv ${inDir_phylop}tmp ${file2};
#
#done
#
## Convert the white space delimited text files to tab delimited format
#for l in ${inDir_phylop}*.sorted.bed; do
#	awk -v OFS='\t' '{$1=$1; print}' ${l} > ${inDir_phylop}tmp
#	chmod 777 ${inDir_phylop}tmp
#	mv ${inDir_phylop}tmp ${l};
#
#done

## 2) phastCons primate files.

# Remove 4th and 5th columns as they are not needed
#for i in ${inDir_phastcons}*moreThan*; do
#
#       awk '{$4=""; $5=""; print $0}' ${i} > ${inDir_phastcons}tmp && mv ${inDir_phastcons}tmp ${i};
#done
#
### Sort concatenated bed file based on 1st and 2nd columns
#for j in ${inDir_phastcons}*moreThan*; do
#
#       sort -k1,1 -k2,2n ${j} > ${j%.bed}.sorted.bed;
#
#done
#
### Remove X and Y chromosomes as these .bed files will
### be used by LDSC, and LDSC cannot handle X chromosome
### data.
#for file2 in ${inDir_phastcons}*.sorted.bed; do
#
#      awk '$2!="X" && $2!="Y"' ${file2} > ${inDir_phastcons}tmp && mv ${inDir_phastcons}tmp ${file2};
#
#done
#
### Convert the white space delimited text files to tab delimited format
#for l in ${inDir_phastcons}*.sorted.bed; do
#
#       awk -v OFS='\t' '{$1=$1; print}' ${l} > ${inDir_phastcons}tmp
#       chmod 777 ${inDir_phastcons}tmp
#       mv ${inDir_phastcons}tmp ${l};
#
#done

##############################
# Make a bed file including all SNPs with a phyloP score
# NOTE: You need to specify a path other than /tmp while
# sorting these bed files, cuz they are too big for sort to
# handle using /tmp folder.

#awk '{$4=""; $5=""; print $0}' ${inDir_phylop}allchr.phyloP46way.primate.bed > ${inDir_phylop}allchr.phyloP46way.primate_control.bed
#sort -T /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phyloP -k1,1 -k2,2n ${inDir_phylop}allchr.phyloP46way.primate_control.bed > ${inDir_phylop}allchr.phyloP46way.primate.sorted_control.bed
#awk '$2!="X" && $2!="Y"' ${inDir_phylop}allchr.phyloP46way.primate.sorted_control.bed > ${inDir_phylop}tmp && mv ${inDir_phylop}tmp ${inDir_phylop}allchr.phyloP46way.primate.sorted_control.bed
#awk -v OFS='\t' '{$1=$1; print}' ${inDir_phylop}allchr.phyloP46way.primate.sorted_control.bed > ${inDir_phylop}tmp
#chmod 777 ${inDir_phylop}tmp
#mv ${inDir_phylop}tmp ${inDir_phylop}allchr.phyloP46way.primate.sorted_control.bed

# Make another bed file including all SNPs with a phastCons
# score.
#awk '{$4=""; $5=""; print $0}' ${inDir_phastcons}allchr.phastCons46way.primate.bed > ${inDir_phastcons}allchr.phastCons46way.primate_control.bed
#sort -T /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/phastCons -k1,1 -k2,2n ${inDir_phastcons}allchr.phastCons46way.primate_control.bed > ${inDir_phastcons}allchr.phastCons46way.primate.sorted_control.bed
#awk '$2!="X" && $2!="Y"' ${inDir_phastcons}allchr.phastCons46way.primate.sorted_control.bed > ${inDir_phastcons}tmp && mv ${inDir_phastcons}tmp ${inDir_phastcons}allchr.phastCons46way.primate.sorted_control.bed
#awk -v OFS='\t' '{$1=$1; print}' ${inDir_phastcons}allchr.phastCons46way.primate.sorted_control.bed > ${inDir_phastcons}tmp
#chmod 777 ${inDir_phastcons}tmp
#mv ${inDir_phastcons}tmp ${inDir_phastcons}allchr.phastCons46way.primate.sorted_control.bed
      
##############################
# Make a bed file from the HAQER csv.

# Remove the 4th column
#awk '{$4=""; print $0}' ${inDir_haqer}HAQER.csv > ${inDir_haqer}HAQER.bed

# Sort concatenated bed file based on 1st and 2nd columns
#sort -k1,1 -k2,2n ${inDir_haqer}HAQER.bed > ${inDir_haqer}HAQER.sorted.bed

# Remove X and Y chromosomes as these .bed files will
# be used by LDSC, and LDSC cannot handle X chromosome
# data.
#awk '$2!="X" && $2!="Y"' ${inDir_haqer}HAQER.sorted.bed > ${inDir_haqer}tmp && mv ${inDir_haqer}tmp ${inDir_haqer}HAQER.sorted.bed

# Convert the white space delimited text files to tab delimited format
#awk -v OFS='\t' '{$1=$1; print}' ${inDir_haqer}HAQER.sorted.bed > ${inDir_haqer}tmp
#chmod 777 ${inDir_haqer}tmp
#mv ${inDir_haqer}tmp ${inDir_haqer}HAQER.sorted.bed

##############################
# Make bed files for Nott et al. annotations.

# Sort tab files based on 1st and 2nd columns
#for i in ${inDir_beds}*.tab; do
#	
#	sort -k1,1 -k2,2n ${i} > ${i%.tab}.sorted.bed;
#
#done

# Remove X and Y chromosomes.
#for i in ${inDir_beds}*.sorted.bed; do
#
#	awk '$1!="chrX" && $1!="chrY"' ${i} > ${inDir_beds}tmp && mv ${inDir_beds}tmp ${i};
#
#done
