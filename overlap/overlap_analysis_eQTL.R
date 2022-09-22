#!/bin/Rscript
#' ---
#' title: "overlap_analysis_eQTL"
#' author: "Gokberk Alagoz"
#' date: "October 10, 2021, last updated on: 21.04.22"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors=FALSE)
library(GenomicRanges)
library(biomaRt)
library(tidyr)

#' 
#' ## R Markdown
#' 
## ----paths, echo=FALSE-------------------------------------------------------------------------------------------------------------------------------

bedfile = args[1] # "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/GRCh37/fetal_hge_hg19.merged.sorted.bed"
outDir = args[2] # "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/globularity/evolution/results"
genotypeF = "/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr"

#' 
#' ## Including Plots
#' 
## ----parse, echo=FALSE-------------------------------------------------------------------------------------------------------------------------------

# eQTL data downloaded from PsychENCODE
feqtl = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/old_results/eqtl/DER-08a_hg19_eQTL.significant.txt"

# Read .clumped SNPs
clumped = read.table("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/overlap/plink_clump/NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt.clumped", header = T)

phenoname = "dys_rhy"

annot_name = unlist(strsplit(unlist(strsplit(bedfile, "/", fixed=T))[15], ".", fixed=T))[1]

#' 
## ----globularity-------------------------------------------------------------------------------------------------------------------------------------

for (i in 1:nrow(clumped)) {
  
  # Find all SNPs in LD (r2>0.6) with each clumped SNP
  
  system(paste0("module load plink/1.9b6 \
                   plink --bfile /data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",
                clumped$SNP[i], " --out /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/overlap/tmpld_", clumped$SNP[i]))
  
  # Read in the LD calculated from plink
  
  LD = read.table(paste0("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/overlap/tmpld_", clumped$SNP[i], ".ld"), header = T)
  
  # Turn into a genomic ranges object
  
  if (i==1) {
    LDSNPs = GRanges(LD$CHR_B, IRanges(LD$BP_B, LD$BP_B), SNP = LD$SNP_B, indexSNP = LD$SNP_A)
  } else {
    LDSNPs = c(GRanges(LD$CHR_B, IRanges(LD$BP_B, LD$BP_B), SNP = LD$SNP_B, indexSNP = LD$SNP_A), LDSNPs)
  }
}

# Read in the eQTL data
eqtl = read.table(feqtl, header=TRUE)

eqtl.GR = GRanges(gsub("chr", "", eqtl$SNP_chr), 
                  IRanges(eqtl$SNP_start, eqtl$SNP_end))

mcols(eqtl.GR) = eqtl[, c(1:8, 12:15)]

# Read in the BED file containing the annotation file
annot = read.table(bedfile, header=FALSE)
annotGR = GRanges(gsub("chr", "", annot$V1), IRanges(annot$V2, annot$V3), type = annot$V4)
        
# Find overlaps
olap = findOverlaps(annotGR, LDSNPs)
dys_rhy_vars = unique(LDSNPs$indexSNP[subjectHits(olap)])
dys_rhy_Annot = LDSNPs[which(!is.na(match(LDSNPs$indexSNP, dys_rhy_vars)))]

# Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with ', annot_name, ' is: ', length(unique(dys_rhy_Annot$indexSNP)), '\n')

olap2 = findOverlaps(eqtl.GR, dys_rhy_Annot)
dys_rhy_eqtl_vars = unique(dys_rhy_Annot$indexSNP[subjectHits(olap2)])

eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])

# Remove the . annotation
eqtlgenes = sapply(eqtlgenes, function (x) {unlist(strsplit(x, ".", fixed=TRUE))[1]})

if (length(dys_rhy_vars) > 0) {
  
  m = matrix(NA, nrow = length(dys_rhy_vars), ncol = 2)
  d = as.data.frame(m)
  colnames(d) = c("olap_snps", "eqtl_snps")
  
  d$olap_snps = dys_rhy_vars
  d$eqtl_snps[match(dys_rhy_eqtl_vars, d$olap_snps)] = "Yes"
  d$olap_snps[match(dys_rhy_eqtl_vars, d$olap_snps)]
  write.csv(d, file = paste0(outDir, "/dys_rhy_", annot_name, "_olap_snps.csv"),
            row.names = FALSE, quote = FALSE)
} else {
  stop("There are not any pleiotropy-associated variants - annotation overlap")
}

if (length(dys_rhy_eqtl_vars) != 0) {
  
  cat('Number of ', annot_name, ' overlapping loci that also have an eQTL: ', length(dys_rhy_eqtl_vars), '\n')
  
  # Convert these genes to hgnc_id
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                 dataset="hsapiens_gene_ensembl",
                 host="feb2014.archive.ensembl.org")
  geneannot = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id",
                    values = eqtlgenes,
                    mart = mart)
  cat('These eQTLs impact ', length(which(geneannot == "protein_coding")), ' protein-coding eGenes\n')
  write.csv(geneannot, file = paste0(outDir, "/dys_rhy_", annot_name, "_eqtl_genes.csv"),
            row.names = FALSE, quote = FALSE)
  
} else {
  stop("There are not any pleiotropy-associated variant - eQTL overlap")
}

#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()
#knitr::purl("overlap_analysis_eQTL.Rmd", "overlap_analysis_eQTL.R", documentation = 2) # save .Rmd as an .R file to submit to Grid
