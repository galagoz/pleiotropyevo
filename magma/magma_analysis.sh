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
annot="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/magma/input_files/evo_annots.gmt"

# 1) Annotation

#${programs}magma --annotate --snp-loc ${inDir}snp_loc_input.tab --gene-loc ${inDir}geneLoc_file.tab --out ${outDir}dys_rhyimp_CPM_genomicSEM

${programs}magma --annotate --snp-loc ${inDir}snp_loc_input_dys_IPM.tab --gene-loc ${inDir}geneLoc_file.tab --out ${outDir}dys_IPM_genomicSEM

${programs}magma --annotate --snp-loc ${inDir}snp_loc_input_rhyimp_IPM.tab --gene-loc ${inDir}geneLoc_file.tab --out ${outDir}rhyimp_IPM_genomicSEM

# 2) Gene analysis (using SNP p-values)

#${programs}magma --bfile ${genotype_f}g1000_eur --pval ${inDir}pval_file4magma.tab ncol=N --gene-annot ${outDir}dys_rhyimp_CPM_genomicSEM.genes.annot --out ${outDir}dys_rhyimp_CPM_genomicSEM

${programs}magma --bfile ${genotype_f}g1000_eur --pval ${inDir}pval_file4magma_dys_IPM.tab ncol=N --gene-annot ${outDir}dys_IPM_genomicSEM.genes.annot --out ${outDir}dys_IPM_genomicSEM

${programs}magma --bfile ${genotype_f}g1000_eur --pval ${inDir}pval_file4magma_rhyimp_IPM.tab ncol=N --gene-annot ${outDir}rhyimp_IPM_genomicSEM.genes.annot --out ${outDir}rhyimp_IPM_genomicSEM

# 3) Gene-set analysis

#${programs}magma --gene-results ${outDir}dys_rhyimp_CPM_genomicSEM.genes.raw --set-annot ${annot} --out ${outDir}evo_annots

${programs}magma --gene-results ${outDir}dys_IPM_genomicSEM.genes.raw --set-annot ${annot} --out ${outDir}evo_annots_dys_IPM

${programs}magma --gene-results ${outDir}rhyimp_IPM_genomicSEM.genes.raw --set-annot ${annot} --out ${outDir}evo_annots_rhyimp_IPM
