#!/bin/bash

# This script will generate bed files for magma gene-set analysis.
#
# Gokberk Alagoz - 04.08.2022 - last updated on 08.08.2022

##########
# PATHS

resources="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/"
annots="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/resources/new_annotations/beds4magma/"

module load bedtools/2.29.2

####################
# Convert GTF to BED

# 1) Download the GTF file for hg19 from https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/
# and gunzip.

# 2) Remove first 5 lines of the GTF file, because they are irrelevant.
# Here are those lines:
#!genome-build GRCh37.p13
#!genome-version GRCh37
#!genome-date 2009-02
#!genome-build-accession NCBI:GCA_000001405.14
#!genebuild-last-updated 2013-09

sed -i 1,5d Homo_sapiens.GRCh37.87.gtf

# 3) Extract protein coding genes from the GTF file
grep "protein_coding" Homo_sapiens.GRCh37.87.gtf > Homo_sapiens.GRCh37.87_proteinCoding.gtf
awk '$3 == "gene" {print $0}' Homo_sapiens.GRCh37.87_proteinCoding.gtf > Homo_sapiens.GRCh37.87_proteinCodingGenes.gtf

# 4) Convert gtf2bed
# Format should be like this:
# CHR BP_START BP_END GENE_NAME
awk 'OFS="\t" {print $1, $4, $5, $16}' Homo_sapiens.GRCh37.87_proteinCodingGenes.gtf > Homo_sapiens.GRCh37.87_proteinCodingGenes.bed
# clean up this file: i) remove mitochondrial, X, Y chromosome and unmapped genes, ii) delete quotation marks and semi-colons around gene names
sed -Ei 's/"//;s/;//' Homo_sapiens.GRCh37.87_proteinCodingGenes.bed
sed -Ei 's/"//' Homo_sapiens.GRCh37.87_proteinCodingGenes.bed

# Your reference bed file with all genes is ready! -> Homo_sapiens.GRCh37.87_proteinCodingGenes.bed

#################################################
# Preparing BED files for MAGMA gene-set analysis

# 0) We are interested in +-1kb of gene boundries as well, so adapt your gtf2bed file accordingly:
awk 'OFS="\t" {$2=$2-1000; $3=$3+1000} {print $0}' Homo_sapiens.GRCh37.87_proteinCodingGenes.bed > Homo_sapiens.GRCh37.87_proteinCodingGenes_extended+-1kb.bed

# 0.1) Sort all bed files
sort-bed Homo_sapiens.GRCh37.87_proteinCodingGenes_extended+-1kb.bed > tmp && mv tmp Homo_sapiens.GRCh37.87_proteinCodingGenes_extended+-1kb.bed
sort-bed ${annots}*.bed

# 1) Convert the annotation files to gene lists.
bedmap --echo --echo-map-id-uniq ${annots}AMH_derived_DMR_hg19-gokhman_et_al.sorted.bed Homo_sapiens.GRCh37.87_proteinCodingGenes_extended+-1kb.bed > ${annots}AMH-derivedDMR_genes.bed
bedmap --echo --echo-map-id-uniq ${annots}AncientSelectiveSweeps_hg19.sorted.bed Homo_sapiens.GRCh37.87_proteinCodingGenes_extended+-1kb.bed > ${annots}AncientSelectiveSweeps_genes.bed
bedmap --echo --echo-map-id-uniq ${annots}human-chimp_DMR_hg19-gokhman_et_al.sorted.bed Homo_sapiens.GRCh37.87_proteinCodingGenes_extended+-1kb.bed > ${annots}human-chimp_DMR_genes.bed
bedmap --echo --echo-map-id-uniq ${annots}HAR_hg19.sorted.bed Homo_sapiens.GRCh37.87_proteinCodingGenes_extended+-1kb.bed > ${annots}HAR_genes.bed
#awk 'OFS="\t" {print $1, $2, $3}' HAR_hg19.sorted.bed | sed -Ei 's/chr//' > HAR_hg19.new.bed

# 2) Extract gene names from *_genes.bed files, make a gene list per annotation.
for bed in ${annots}*_genes.bed; do

        tmp_base=$(basename "$bed")
        tmp_annot=$(cut -d'.' -f1 <<< $tmp_base)

	while read line; do
	
		gene=$(cut -d'|' -f2 <<< $line)
		echo $gene

	done < $bed >> ${annots}${tmp_annot}_gene_list.bed;
done

sed -Ei '/^[[:space:]]*$/d' ${annots}${tmp_annot}_gene_list.bed
for i in *_genes_gene_list.bed; do cat $i | tr ";" "\n" > tmp && mv tmp $i; done
for i in *_genes_gene_list.bed; do sort $i | uniq -u > tmp && mv tmp $i; done

# 3) Using these gene lists, extract gene coordinates from Homo_sapiens.GRCh37.87_proteinCodingGenes.bed
for i  in $annots*_genes_gene_list.bed; do awk 'NR==FNR{a[$1];next} $4 in a {print $0}' $i Homo_sapiens.GRCh37.87_proteinCodingGenes.bed > tmp && mv tmp ${i%.*}_annot.bed; done

# 4) Extract the promoter coordinates per gene from epdnew_hg38ToHg19_promoters.bed (Eukaryotic Promoter Database by SBI), add this information to the annotation file by extending each gene's coordinate over its promoter.
sed -Ei 's/chr//' epdnew_hg38ToHg19_promoters.bed
awk 'OFS="\t" {print $1, $2, $3, $4}' epdnew_hg38ToHg19_promoters.bed > tmp && mv tmp epdnew_hg38ToHg19_promoters.bed
awk '{split($4, a , "_"); $4=a[1]; print $0}' epdnew_hg38ToHg19_promoters.bed > epdnew_hg38ToHg19_promoters_v2.bed

# FINAL NOTE BEFORE THE HOLIDAY: The promoter data seems like it's not suitable for adding to bed files, because it has point-coordinates.
# So maybe use the genes-only bed files or alternatively add eQTL data to bed files. Once you're back, decide this quickly and run MAGMA.i
# -> NOTE AFTER HOLIDAY: You already included +-1kb of gene coordinates from ENSEMBL, so you probably already cover promoters,
# I don't think you need to add the promoter info separately (you can ask Else about this once she is back to make sure).

# When gene boundary+promoter annotation files are ready, feed them into MAGMA gene-set analysis.
