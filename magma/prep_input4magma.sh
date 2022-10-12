#!/bin/bash

# This script will create the necessary input files
# for MAGMA annotation, gene analysis and gene-set
# analysis.

# Gokberk Alagoz, 11.10.22
##########################

# Set Paths

genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"
geneLoc="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/"
geneLists="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/new_annotations/beds4magma/"
outDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/magma/input_files/"

##########################

# 1) Prep for magma --annotate

## a) Reformat Genomic SEM meta-analysis results file
## "The SNP location file should contain three columns: SNP ID, chromosome, and base pair position"

#awk 'OFS="\t" {print $1, $2, $3}' ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2.tab > ${outDir}snp_loc_input.tab

## b) GENELOC_FILE for hg19 is provided in MAGMA website, but you still need to reformat it.
## The gene location file must contain at least four columns, in this order: gene ID, chromosome, start site, stop site.
## NOTE: I will try to use gene symbol instead of gene IDs, cuz that's what I have in my evo annotations :)))))

#awk 'OFS="\t" {print $6, $2, $3, $4}' ${geneLoc}NCBI37.3.gene.loc > ${outDir}geneLoc_file.tab;

##########################

# 2) Prep for magma Gene analysis
# Reformat Genomic SEM meta-analysis results file
# "The p-value file must be a plain text data file with each row corresponding to a SNP. If MAGMA detects a header in the file it will look for SNP IDs and p-values in the SNP and P column respectively. Sample size can also be set by using the ncol modifier to specify a column in the p-value file that contains the sample size used per SNP."

#awk 'OFS="\t" {print $1, $7, $8}' ${genSEM}GenomicSEM_multivarGWAS_dys_rhyimp_v2_4gctb_v3.ma > ${outDir}pval_file_4magma.tab
#sed -Ei '1s/Pval_Estimate/P/' ${outDir}pval_file_4magma.tab

##########################

# 3) Prep for magma Gene-set analysis
# You have to reformat your evo gene lists in .gmt format, which lookd smt like below:
# 1st column -> gene-set name
# 2nd column -> a brief description of the annot
# 3rd column -> gene names (can be unequal lengths)

awk 'BEGIN {ORS = "\t"} {print $1}' AncientSelectiveSweeps_genes_gene_list_annot_4magma.sorted.bed > AncientSelectiveSweeps.gmt
awk 'BEGIN {ORS = "\t"} {print $1}' AMH-derivedDMR_genes_gene_list_annot_4magma.sorted.bed > AMH-derivedDMR.gmt
awk 'BEGIN {ORS = "\t"} {print $1}' HAR_genes_gene_list_annot_4magma.sorted.bed > HAR.gmt
awk 'BEGIN {ORS = "\t"} {print $1}' human-chimp_DMR_genes_gene_list_annot_4magma.sorted.bed > human-chimp_DMR.gmt

awk 'OFS="\t" {print "AncientSelectiveSweeps", "AncientSelectiveSweeps", $0}' AncientSelectiveSweeps.gmt > tmp && mv tmp AncientSelectiveSweeps.gmt
awk 'OFS="\t" {print "AMH-derivedDMR", "AMH-derivedDMR", $0}' AMH-derivedDMR.gmt > tmp && mv tmp AMH-derivedDMR.gmt
awk 'OFS="\t" {print "HAR", "HAR", $0}' HAR.gmt > tmp && mv tmp HAR.gmt
awk 'OFS="\t" {print "human-chimp_DMR", "human-chimp_DMR", $0}' human-chimp_DMR.gmt > tmp && mv tmp human-chimp_DMR.gmt

# You can tell magma what column has what info, so no need to reformat your evo annot files.
# "Gene-set definitions can also be supplied in a column-based form, by adding the col modifier 
# to  the  --set-annot  flag. This modifier  should  specify  two  values,  indicating the  index  for  the  gene  ID 
# column and gene-set name column respectively (eg. ‘--set-annot sets.col col=3,1’ will read gene - gene-
# set pairs from the sets.col file, looking for gene IDs in the third column and gene-set names in the first 
# column)."

#for i in ${geneLists}*.bed; do
#
#       tmp_list=$(basename "$i" .bed)
#       awk 'OFS="\t" {print $4, $1, $2, $3}' ${i} > ${outDir}${tmp_list}_4magma.bed;
#
#done

# Sort files based on CHR and StartPos

#for f in ${outDir}/*_4magma.bed; do
#
#       sort -V -k2,2 -k3,3 $f > ${f%.bed}.sorted.bed;
#
#done
