#!/bin/bash

# This script will 1) run snp2gene (annotation),
#		   2) run gene analysis,
#		   3) run gene-set analysis on evo. gene lists.

# Gokberk Alagoz, 11.10.22
##########################

# Set Paths

programs="/home/gokala/programs/"
inDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/magma/input_files/"
res="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/"
outDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/magma/output_files/"
genotype_f="/data/workspaces/lag/shared_spaces/Resource_DB/magma_v1.10/g1000_eur/"

# 1) Annotation

#${programs}magma --annotate --snp-loc ${inDir}snp_loc_input.tab --gene-loc ${inDir}geneLoc_file.tab --out ${outDir}dys_rhyimp_genomicSEM

# 2) Gene analysis (using SNP p-values)

${programs}magma --bfile ${genotype_f}g1000_eur --pval ${inDir}pval_file_4magma.tab ncol=N --gene-annot ${outDir}dys_rhyimp_genomicSEM.genes.annot --out ${outDir}dys_rhyimp_genomicSEM

# 3) Gene-set analysis

#tmp_annot=(basename "$annot")
#${programs}magma --gene-results ${outDir}dys_rhyimp_genomicSEM.genes.raw --set-annot ${annot} col=1, --out ${outDir}${tmp_annot%.genes_gene_list_annot_4magma.sorted.bed}'
