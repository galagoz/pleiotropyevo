#!/bin/bash
#
# This script will reformat summary statistics from
# various external sources for LDSC munge.
#
# Gokberk Alagoz - 24.01.23
#
#############################################################
# 1) Reformat PGC sumstats. Extract the columns below (column 
# order doesn't matter):
#
# rsID EffAllele RefAllele EffAllFreq BETA/OR SE P Ntotal CHR BP
#
# Note: Each sumstats file comes from a different resource,
# so make sure what genome assembly each one used (most of them
# are still hg19). Also, check if a sumstats was obtained by
# a metaanalysis or not. Depending on this, use N_effective
# or N_total (you can read more about this in LDSC google group).
#############################################################

#######################
## ALZHEIMER'S
#awk -v OFS='\t' '{print $6, $4, $5, $12, $13, $14, $8, $9, $2, $3}' AD_sumstats_Jansenetal_2019sept.txt > AD_sumstats_forMunge.txt
#sed -Ei '1s/Nsum/N/' AD_sumstats_forMunge.txt
#
#######################
## ADHD
## A1 is not available, so use FREQ_U as A1
## See: https://github.com/neurogenomics/MungeSumstats/issues/105
#awk -v OFS='\t' '{print $2, $4, $5, $7, $9, $10, $11, $16+$17, $1, $3}' daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta > ADHD_sumstats_forMunge.txt
#sed -Ei '1s/FRQ_U_34194/FREQ1/;1s/0/N/' ADHD_sumstats_forMunge.txt
#
#######################
## INSOMNIA
#awk -v OFS='\t' '{print $1, $5, $6, $7, $8, $9, $10, $11, $3, $4}' Insomnia_sumstats_Jansenetal.txt > Insomnia_sumstats_forMunge.txt
#
#######################
## PARKINSON'S
## Parkinson's sumstats is a bit tricky to reformat as the file
## doesn't include rsIDs. I'll add these with bedops. This is a bit
## time consuming, but a nice exercise to learn how to get rsIDs 
## if you have only genomic positions.
## (source: https://www.biostars.org/p/312369/)
#module load bedtools/2.29.2
#
## First, download the reference genome (hg19) and convert it to bed.
## (this step takes quite some time, proceed with formatting other sumstats)
#wget -qO- ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz | gunzip -c - | vcf2bed --sort-tmpdir=${PWD} --max-mem=2G - > hg19.dbSNP150.bed
#
## Then, make a bed file from your sumstats.
#awk -v OFS='\t' '{print $1}' nallsEtAl2019_excluding23andMe_allVariants.tab > parkinson_bedToBe.txt
#awk -v OFS='\t' 'NR>1 {split($0,a,":"); print a[1],a[2]-1,a[2]}' parkinson_bedToBe.txt | sort-bed > parkinson.bed
#sed -Ei 's/^chr//g' parkinson.bed
#
## Map your positions to rsIDs.
#bedmap --echo --echo-map-id --delim '\t' parkinson.bed hg19.dbSNP150.bed > parkinson_rsIDs.bed
## Some lines have the rsID twice (e.g. rs12184279;rs12184279), remove these.
#awk -F';' '{print $1}' parkinson_rsIDs.bed > parkinson_rsIDs2.bed
#
#rm parkinson_rsIDs.bed parkinson_rsIDs2.bed
#
## Some lines don't have rsIDs, remove these.
#awk  '$4!=""' parkinson_rsIDs2.bed > parkinson_rsIDs3.bed
## Split 1st column of the original sumstat.
#awk -v OFS='\t' 'NR>1 {split($1,a,":"); print a[1],a[2],$2,$3,$4,$5,$6,$7,$8+$9}' nallsEtAl2019_excluding23andMe_allVariants.tab > Parkinson_sumstats1.txt
#sed -Ei 's/^chr//g' Parkinson_sumstats1.txt
#
## Almost ready for merging! 'Join' quickly merges two files based on
## one column from each file, so make a MARKER column in each file.
#
#awk '{print $1":"$2"\t"$0}' Parkinson_sumstats1.txt > Parkinson_sumstats_forJoin.txt
#awk '{print $1":"$3"\t"$0}' parkinson_rsIDs3.bed > parkinson_rsIDs3_forJoin.txt
#
#join -j 1 -o 2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 <(sort -V -k1 Parkinson_sumstats_forJoin.txt) <(sort -V -k1 parkinson_rsIDs3_forJoin.txt) > Parkinson_sumstats2.txt
#
#echo -e "SNP\tCHR\tBP\tA1\tA2\tFREQ\tBETA\tSE\tP\tN" | cat - Parkinson_sumstats2.txt > Parkinson_sumstats_forMunge.txt
## yay, done! clean up now.
#
#rm Parkinson_sumstats_forJoin.txt Parkinson_sumstats2.txt Parkinson_sumstats1.txt Parkinson_sumstats1.txt parkinson_rsIDs3.bed parkinson_bedToBe.txt parkinson.bed parkinson_rsIDs3_forJoin.txt
#
#####################
## SCZ
## A1 is not available, so use FREQ_U as A1
## See: https://github.com/neurogenomics/MungeSumstats/issues/105
#awk -v OFS='\t' '{print $2, $4, $5, $7, $9, $10, $11, "N", $1, $3}' daner_PGC_SCZ52_0513a.hq2 > SCZ_sumstats_forMunge1.txt
#awk -v OFS='\t' 'NR>1 {$8=150064}1' SCZ_sumstats_forMunge1.txt > SCZ_sumstats_forMunge.txt
#sed -Ei '1s/FRQ_U_46839/FREQ1/' SCZ_sumstats_forMunge.txt
#
#rm SCZ_sumstats_forMunge1.txt
#####################
## BIPOLAR
## Remove VCF header. There is no MAF, so I will just use
## FCAS. FCAS and FCON are quite similar anyway. #TODO See if this is a problem for LDSC!
#
#egrep -v "^#" pgc-bip2021-all.vcf.tsv > pgc-bip2021-all.noHeader.vcf.tsv
#
#awk -v OFS='\t' '{print $3, $4, $5, $10, $6, $7, $8, $14+$15, $1, $2}' pgc-bip2021-all.noHeader.vcf.tsv > BIP_sumstats_forMunge1.txt
#echo -e "SNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tN\tCHR\tBP" | cat - BIP_sumstats_forMunge1.txt > BIP_sumstats_forMunge.txt
#
#####################
## DEPRESSION
## This sumstats has logOdds and logOddsSEs.
## Convert them back to OR and OR_SE. It
## doesn't have CHR and BP columns either...
#
#awk -v OFS='\t' '{print $1, $2, $3, $4, exp($5), exp($6), $7, $8="N"}' PGC_UKB_depression_genome-wide.txt > DEP_sumstats_forMunge1.txt
#awk -v OFS='\t' 'NR>1 {$8=500199; print $0}' DEP_sumstats_forMunge1.txt > DEP_sumstats_forMunge2.txt
#echo -e "SNP\tA1\tA2\tFREQ\tOR\tSE\tP\tN" | cat - DEP_sumstats_forMunge2.txt > DEP_sumstats_forMunge.txt
#
#rm DEP_sumstats_forMunge1.txt DEP_sumstats_forMunge1.txt
#####################
## TOURRETTE
## No allele freq info. #TODO See if this is OK for LDSC!
#
#awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $7, $8, $9, $10=14307}' TS_Oct2018.txt > TS_sumstats_forMunge.txt
#sed -Ei '1s/14307/N/' TS_sumstats_forMunge.txt
#
#####################
## Heaviness of smoking
#
#awk -v OFS='\t' '{print $3, $5, $4, $6, $9, $10, $8, $11, $1, $2}' CigarettesPerDay.txt > HeavinessOfSmoking_forMunge.txt
#sed -Ei '1s/ALT/A1/;1s/REF/A2/;1s/AF/MAF/' HeavinessOfSmoking_forMunge.txt
#
#####################
## Smoking initiation
#
#awk -v OFS='\t' '{print $3, $5, $4, $6, $9, $10, $8, $11, $1, $2}' SmokingInitiation.txt > SmokingInitiation_forMunge.txt
#sed -Ei '1s/ALT/A1/;1s/REF/A2/;1s/AF/MAF/' SmokingInitiation_forMunge.txt
#
#####################
## Autism Spectrum Disorder
# No MAF column...

#awk -v OFS='\t' '{print $2, $4, $5, $7, $8, $9, "46350",$1, $3}' iPSYCH-PGC_ASD_Nov2017 > Autism_forMunge.txt
#sed -Ei '1s/46350/N/' Autism_forMunge.txt

#####################
## Anxiety Disorder

#egrep -v "^#" pgc-panic2019.vcf.tsv > pgc-panic2019.noHeader.vcf.tsv
#awk -v OFS='\t' '{print $3, $4, $5, $10, $6, $7, $8, $14+$15, $1, $2}' pgc-panic2019.noHeader.vcf.tsv > Anxiety_sumstats_forMunge1.txt
#echo -e "SNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tN\tCHR\tBP" | cat - Anxiety_sumstats_forMunge1.txt > Anxiety_sumstats_forMunge.txt

# rm Anxiety_sumstats_forMunge1.txt
#####################
## Processing Speed
# No MAF column...

#awk -v OFS='\t' '{print $1, $4, $5, $6, $7, $9, "332050",$2, $3}' gFactor_sumstats_noMAF.txt > gFactor_forMunge.txt
#sed -Ei '1s/Estimate/BETA/;1s/Pval_Estimate/P/;1s/332050/N/' gFactor_forMunge.txt

#####################
## Risk taking behavior

#awk -v OFS='\t' '{print $1, $4, $5, $6, $7, $8, $9, "466571",$2, $3}' RISK_GWAS_MA_UKB+replication.txt > RiskTakingBehavior_forMunge.txt
#sed -Ei '1s/MarkerName/SNP/;1s/EAF_A1/EAF/;1s/Beta/BETA/;1s/Pval/P/;1s/466571/N/' RiskTakingBehavior_forMunge.txt

########################################
# UK Biobank traits
inDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/external_sumstats/UKB_traits/"

#awk -v OFS='\t' '{print $6, $5, $4, $15, $2, $3}' ${inDir}variants.tsv > ${inDir}variants_subset.tsv

while read file; do

        tmp_name=$(basename "${file}" .tsv)
        paste ${inDir}variants_subset.tsv ${inDir}${file} > ${inDir}${tmp_name}_forMunge1.txt
        awk -v OFS='\t' '{print $1, $2, $3, $9, $14, $15, $17, $11, $5, $6}' ${inDir}${tmp_name}_forMunge1.txt > ${inDir}${tmp_name}_forMunge.txt
        sed -Ei '1s/rsid/SNP/;1s/alt/A1/;1s/ref/A2/;1s/minor_AF/MAF/;1s/beta/BETA/;1s/se/SE/;1s/pval/P/;1s/n_complete_samples/N/;1s/chr/CHR/;1s/pos/BP/' ${inDir}${tmp_name}_forMunge.txt;

done < ${inDir}UKBsumstats_list_diffColumnNumbers.txt

#rm ${inDir}*forMunge1*

# Some UKB sumstats has an extra column inbetween called
# "expected_min_category_minor_AC" (why???), so you have to
# format those ones separately.

#for file in ${inDir}*.tsv; do
#
#       tmp_name=$(basename "$file" .tsv)
#       paste ${inDir}variants_subset.tsv ${file} > ${inDir}${tmp_name}_forMunge1.txt
#       awk -v OFS='\t' '{print $1, $2, $3, $9, $15, $16, $18, $12, $5, $6}' ${inDir}${tmp_name}_forMunge1.txt > ${inDir}${tmp_name}_forMunge.txt
#       sed -Ei '1s/rsid/SNP/;1s/alt/A1/;1s/ref/A2/;1s/minor_AF/MAF/;1s/beta/BETA/;1s/se/SE/;1s/pval/P/;1s/n_complete_samples/N/;1s/chr/CHR/;1s/pos/BP/' ${inDir}${tmp_name}_forMunge.txt;
#
#done
#
#rm ${inDir}*forMunge1*
