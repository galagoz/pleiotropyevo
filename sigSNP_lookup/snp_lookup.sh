!/bin/bash
#
# Gokberk Alagoz - 19.05.22
#
# This script will subset dyslexia and rhythm summary statistics using
# the list of gw-sig-SNPs in the MA results.

# First download this list from FUMA. #TODO get this list manually with PLINK

##########################
# Variables
SNPs4filtering="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/FUMA/FUMA_job178902_dys_rhyimp_MA/IndSigSNPs_rsIDs.txt"
MA_sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt"
rhy_sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/clap_to_beat.merged.1kgenomes_02222019.txt"
rhyimp_sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythym_reformatted_forNweighted5.txt"
dys_sumstats="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia.filtered.2.dat"
outDir="/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/snp_lookup"

##### P-value lookup ######
###########################
# Lookup in rhy sumstats

# head -n 1 $rhy_sumstats > ${outDir}/dysrhyimpMA_indSigSNPs_lookup_in_rhy.txt
# awk 'BEGIN{OFS="\t"} FNR==NR{a[$22]=$0; next} {print a[$1]}' $rhy_sumstats $SNPs4filtering >> ${outDir}/dysrhyimpMA_indSigSNPs_lookup_in_rhy.txt

############################
# Lookup in dys sumstats

# head -n 1 $dys_sumstats > ${outDir}/dysrhyimpMA_indSigSNPs_lookup_in_dys.txt
# awk 'BEGIN{OFS="\t"} FNR==NR{a[$19]=$0; next} {print a[$1]}' $dys_sumstats $SNPs4filtering >> ${outDir}/dysrhyimpMA_indSigSNPs_lookup_in_dys.txt

############################
###### Effect lookup #######

# Get effect sizes from the MA
# head -n 1 $MA_sumstats > ${outDir}/dysrhyimpMA_indSigSNPs_effectsize_lookup_in_MA.txt
# awk 'BEGIN{OFS="\t"} FNR==NR{a[$2]=$0; next} {print a[$1]}' $MA_sumstats $SNPs4filtering >> ${outDir}/dysrhyimpMA_indSigSNPs_effectsize_lookup_in_MA.txt

# Get Zscores from the rhy imp sumstats file (rhy sumstats with flipped effect directions)
head -n 1 $rhyimp_sumstats > ${outDir}/dysrhyimpMA_indSigSNPs_effectsize_lookup_in_rhyimp.txt
awk 'BEGIN{OFS="\t"} FNR==NR{a[$1]=$0; next} {print a[$1]}' $rhyimp_sumstats $SNPs4filtering >> ${outDir}/dysrhyimpMA_indSigSNPs_effectsize_lookup_in_rhyimp.txt
