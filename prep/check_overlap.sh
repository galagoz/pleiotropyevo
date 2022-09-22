#!/bin/bash
# Gokberk Alagoz - 11.11.21
# This script will find the overlapping SNPs in clap-beat and
# dyslexia GWAS sumstats and will create subsets with only 
# overlapping SNPs for both traits.

# Subset & sort clap-beat sumstats
awk 'NR==FNR{a[$42];next} $19 in a {print $0}' clap_to_beat.merged.1kgenomes_02222019.txt dyslexia.filtered.2.dat > dyslexia_overlap_subset.tab
awk -v OFS='\t' '{print ($42, $3)}' clap-beat_overlap_subset.tab > clap-beat_overlap_rsID_pval.tab
sort -k1 -n clap-beat_overlap_rsID_pval.tab > clap-beat_overlap_rsID_pval.sorted.tab

# Subset & sort dyslexia sumstats
awk 'NR==FNR{a[$19];next} $42 in a {print $0}' dyslexia.filtered.2.dat clap_to_beat.merged.1kgenomes_02222019.txt > clap-beat_overlap_subset.tab
awk -v OFS='\t' '{print ($19, $3)}' dyslexia_overlap_subset.tab > dyslexia_overlap_rsID_pval.tab
sort -k1 -n dyslexia_overlap_rsID_pval.tab > dyslexia_overlap_rsID_pval.sorted.tab

# Merge files
awk -v OFS='\t' '{getline f1 <"clap-beat_overlap_rsID_pval.sorted.tab" ; print f1,$2}' < dyslexia_overlap_rsID_pval.sorted.tab > clapbeat_dyslexia_pvals.tab
sed '1iSNP	beatclap	dyslexia' clapbeat_dyslexia_pvals.tab > clapbeat_dyslexia_pvals2.tab

awk -v OFS='\t' '{getline f1 <"clap-beat_overlap_rsID_pval.sorted.tab" ; print f1,$2}' < dyslexia_overlap_rsID_pval.sorted.tab > clapbeat_dyslexia_pvals.tab

# Make DFs with MARKER info to make Manhattan plots

awk -v OFS='\t' '{print ($42, $23, $24, $3)}' clap-beat_overlap_subset.tab > clap-beat_overlap_forManhattan.tab
awk -v OFS='\t' '{print ($19, $16, $17,$3)}' dyslexia_overlap_subset.tab > dyslexia_overlap_forManhattan.tab

# remove "chr"
awk '{ gsub("chr", "", $2); print }' clap-beat_overlap_forManhattan.tab > clap-beat_overlap_forManhattan2.tab
awk '{ gsub("chr", "", $2); print }' dyslexia_overlap_forManhattan.tab > dyslexia_overlap_forManhattan2.tab

# merge chr and pos as marker
awk '{print $1"\t"$2":"$3"\t"$4}' < clap-beat_overlap_forManhattan2.tab > clap-beat_overlap_forManhattan3.tab
awk '{print $1"\t"$2":"$3"\t"$4}' < dyslexia_overlap_forManhattan2.tab > dyslexia_overlap_forManhattan3.tab

######## Edit after a while
# Decided to split MARKER column again

awk 'sub(/\:/," "){$1=$1}1' OFS="\t" clap-beat_overlap_forManhattan3.tab > clap-beat_overlap_forManhattan4.tab
awk 'sub(/\:/," "){$1=$1}1' OFS="\t" dyslexia_overlap_forManhattan3.tab > dyslexia_overlap_forManhattan4.tab
