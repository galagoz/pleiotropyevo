#!/bin/bash
#
# Gokberk Alagoz, March 2022 - last updated Jan 2023
#
#-----
#PATHS
mungedDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/munged/"
ldscores="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/eur_w_ld_chr/"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/snph2_estimates/"

#-----
#module purge
#module load miniconda/3.2021.10 ldsc/v1.0.1
#conda activate ldsc
#-----
# 1) Dyslexia

## Without liability scale
#ldsc.py \
# --h2 ${mungedDir}dyslexia.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --out ${outDir}dyslexia_h2 \
# --w-ld-chr ${ldscores}
#
## With liability scale
#ldsc.py \
# --h2 ${mungedDir}dyslexia.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --pop-prev 0.05 \
# --samp-prev 0.045 \
# --out ${outDir}dyslexia_pop0.05_samp0.045_h2 \
# --w-ld-chr ${ldscores}
#
##-----
## 2) Rhythm impairment
#
## Without liability scale
#ldsc.py \
# --h2 ${mungedDir}rhythm_impairment.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --out ${outDir}rhythm_impairment_h2 \
# --w-ld-chr ${ldscores}
#
## With liability scale
#ldsc.py \
# --h2 ${mungedDir}rhythm_impairment.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --pop-prev 0.0475 \
# --samp-prev 0.085 \
# --out ${outDir}rhythm_impairment_pop0.0475_samp0.085_h2 \
# --w-ld-chr ${ldscores}
#
#-----
# 3) Genomic SEM results

# CPM results

# Without liability scale
ldsc.py \
 --h2 ${mungedDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp.sumstats.gz \
 --ref-ld-chr ${ldscores} \
 --out ${outDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_heritability \
 --w-ld-chr ${ldscores}

# With liability scale
#ldsc.py \
# --h2 ${mungedDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --pop-prev 0.0488 \
# --samp-prev 0.065 \
# --out ${outDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_popPrev0.0488_sampPrev0.065_heritability \
# --w-ld-chr ${ldscores}

# IPM results:

# Dyslexia

# Without liability scale
#ldsc.py \
# --h2 ${mungedDir}GenomicSEM_multivarGWAS_dys_IPM.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --out ${outDir}GenomicSEM_multivarGWAS_dys_IPM_heritability \
# --w-ld-chr ${ldscores}

# With liability scale
#ldsc.py \
# --h2 ${mungedDir}GenomicSEM_multivarGWAS_dys_IPM.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --pop-prev 0.0 \ # ?
# --samp-prev 0.0 \ # ?
# --out ${outDir}GenomicSEM_multivarGWAS_dys_IPM_popPrev0.0488_sampPrev0.065_heritability \
# --w-ld-chr ${ldscores}

# Rhythm impairment

# Without liability scale
ldsc.py \
 --h2 ${mungedDir}GenomicSEM_multivarGWAS_rhyimp_IPM.sumstats.gz \
 --ref-ld-chr ${ldscores} \
 --out ${outDir}GenomicSEM_multivarGWAS_rhyimp_IPM_heritability \
 --w-ld-chr ${ldscores}

# With liability scale
#ldsc.py \
# --h2 ${mungedDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp.sumstats.gz \
# --ref-ld-chr ${ldscores} \
# --pop-prev 0.0 \ # ?
# --samp-prev 0.0 \ # ?
# --out ${outDir}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_popPrev0.0488_sampPrev0.065_heritability \
# --w-ld-chr ${ldscores}
