#!/bin/bash
#
# This script will reformat dyslexia, rhythm impairment and 
# Genomic SEM common factor sumstats for FUMA.
#
# All columns need for FUMA are already there, you just need
# to add the total N column to each file.
#####

sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"
genSEM="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"

# 1) Dyslexia

# add total N column
#awk 'NR>1 {print $0"\t"$7+$8}' ${sumstats}dyslexia.filtered.2.dat > ${sumstats}tmp
#echo -e "all.data.id     src     pvalue  effect  stderr  pass    im.num.0        im.num.1        AA.0    AB.0    BB.0    AA.1    AB.1    BB.1    im.data.id      chr     position        alleles rsid    dose.b  freq.b  avg.rsqr        min.rsqr        alleleA alleleB OR      chrpos	N"| cat - ${sumstats}tmp > ${sumstats}dyslexia.filtered.2_wtotalN.txt
#rm ${sumstats}tmp

# 2) Rhythm impairment

# Convert log odss to odds ratios and rename column
#awk 'NR>1 {print $0"\t"$7+$9"\t"exp($4)}' ${sumstats}rhythm_impairment_lambdaGCcorrected.dat > ${sumstats}tmp
#echo -e "all.data.id     src     pvalue  effect  stderr  pass    im.num.0        dose.b.0        im.num.1        dose.b.1        AA.0    AB.0    BB.0    AA.1    AB.1    BB.1   freq.a   freq.b  maf     gt.data.id      im.data.id      assay.name      scaffold        position        alleles ploidy  strand  cytoband        gene.context    is.v1  is.v2    is.v3   is.v4   is.v5   h550    omni    chr     ref_allele      effect_allele   genome1k.mapped rsid.mapped.mismatch    snp.id	N	OR"| cat - ${sumstats}tmp > ${sumstats}rhythm_impairment_lambdaGCcorrected_wtotalNandOR.dat
#rm ${sumstats}tmp

#awk 'BEGIN{OFS="\t"} {print $0, exp($4)}' ${sumstats}/ > ${sumstats}/rhythm_impairment_raw_sumstat_wtotalNandOR.txt
#awk 'NR==1{$44="OR"}1' ${sumstats}/rhythym_impairment_raw_sumstat_wtotalNandOR.txt > tmp && mv tmp ${sumstats}/rhythm_impairment_raw_sumstat_wtotalNandOR.txt

# 3) Genomic SEM results

awk 'NR==FNR{a[NR]=$2;next} {print $0"\t"a[FNR]}' ${genSEM}perSNP_Ns_forgenomicSEM.txt ${genSEM}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr.tab > ${genSEM}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN.tab

chmod 777 ${genSEM}GenomicSEM_multivarGWAS_CPM_dys_rhyimp_withN.tab
