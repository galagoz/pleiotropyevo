#--------------------------
#---- Clumping w/PLINK ----
#--------------------------

# clump sumstats before plotting phastCons
# vs. GWAS_pvalues.

#-----Variables-----
sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN.tab"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/plink/"
genotypeFile="/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr"
#-----
echo "Starting to clump..."
echo $sumstats
base_name=$(basename "$sumstats")

echo "#!/bin/sh
#$ -N plink_clump
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
cd '$outDir'
module load plink/1.9b6
plink --bfile '$genotypeFile' --clump '$sumstats' --clump-field 'Pval_Estimate' --clump-r2 0.6 --clump-kb 100000 --clump-p1 0.05 --out '${outDir}${base_name}'" >> ${outDir}${base_name}.sh

chmod a+x ${outDir}${base_name}.sh
qsub -o ${base_name}.out -j y "${outDir}${base_name}.sh"

echo "Done!"
