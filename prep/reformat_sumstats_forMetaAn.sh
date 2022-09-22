#!/bin/bash
#
# This script will extract and order required columns from Dyslexia and Rhythm
# GWAS summary statictics. Then will rename field names accordingly.

# Dyslexia

awk '{print $19, $16, $17, $25, $24, $20, $7+$8, $4/$5, $3}' dyslexia.filtered.2.dat > dyslexia_reformatted_forNweighted.txt
awk 'NR>1 {gsub("chr", "", $2); print}' dyslexia_reformatted_forNweighted.txt > dyslexia_reformatted_forNweighted2.txt
echo -e "SNPID CHR BP EA OA EAF N Z P" | cat - dyslexia_reformatted_forNweighted2.txt > dyslexia_reformatted_forNweighted3.txt

awk '{print $16":"$17, $19, $16, $17, $25, $24, $20, $7+$8, $4/$5, $3, $4, $5}' dyslexia.filtered.2.dat > dyslexia_reformatted_forModelAveraged.txt
awk 'NR>1 {gsub("chr", "", $1); print}' dyslexia_reformatted_forModelAveraged.txt > dyslexia_reformatted_forModelAveraged1.txt
awk 'NR>1 {gsub("chr", "", $3); print}' dyslexia_reformatted_forModelAveraged1.txt > dyslexia_reformatted_forModelAveraged2.txt
echo -e "cptid RS CHR BP A1 A2 EAF N Z PVAL Beta SE" | cat - dyslexia_reformatted_forModelAveraged2.txt > dyslexia_reformatted_forModelAveraged3.txt

chmod 777 dyslexia_reformatted_forNweighted3.txt
chmod 777 dyslexia_reformatted_forModelAveraged3.txt

# Rhythym

awk '{print $42, $37, $24, $39, $38, $18, $7+$9, $4/$5, $3}' clap_to_beat.merged.1kgenomes_02222019.txt > rhythym_reformatted_forNweighted.txt
awk 'NR>1 {gsub("chr", "", $2); print}' rhythym_reformatted_forNweighted.txt > rhythym_reformatted_forNweighted2.txt
echo -e "SNPID CHR BP EA OA EAF N Z P" | cat - rhythym_reformatted_forNweighted2.txt > rhythym_reformatted_forNweighted3.txt

awk '{print $37":"$24, $42, $37, $24, $39, $38, $18, $7+$9, $4/$5, $3, $4, $5}' clap_to_beat.merged.1kgenomes_02222019.txt > rhythym_reformatted_forModelAveraged.txt
awk 'NR>1 {gsub("chr", "", $1); print}' rhythym_reformatted_forModelAveraged.txt > rhythym_reformatted_forModelAveraged1.txt
awk 'NR>1 {gsub("chr", "", $3); print}' rhythym_reformatted_forModelAveraged1.txt > rhythym_reformatted_forModelAveraged2.txt
echo -e "cptid RS CHR BP A1 A2 EAF N Z PVAL Beta SE" | cat - rhythym_reformatted_forModelAveraged2.txt > rhythym_reformatted_forModelAveraged3.txt

chmod 777 rhythym_reformatted_forNweighted3.txt
chmod 777 rhythym_reformatted_forModelAveraged3.txt
