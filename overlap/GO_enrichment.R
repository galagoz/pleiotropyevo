## Run it with R version 4.0.3 and biomaRt v2.45.9.
## Otherwise the getBM function does not work.
## Run gene-ontology enrichments on variants within 
## HGEs that also impact brain structure
## mapping to genes using PsychENCODE eQTL

options(stringsAsFactors=FALSE)
library(GenomicRanges)
library(biomaRt)
#library(AnnotationHub)

##BED file containing HGE_7PCW, as well as other annotations
fannot = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/old_results/partherit/beds/new_beds2/neanDepRegions_hg19.sorted.bed"
#allAnnots_MAe3ukw3_shortnames.bed

##eQTL data downloaded from PsychENCODE
feqtl = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/old_results/eqtl/DER-08a_hg19_eQTL.significant.txt"

##SNP information corresponding to eQTL data
# fsnp = "SNP_Information_Table_with_Alleles.txt"; # it is not used in the script, so I commented it out.

##ancestry regressed gwas summary stats
fGWASsumstats = "/data/clusterfs/lag/users/gokala/enigma-evol/data/european_hemave/sumstats_txt_list.txt"

##Clumped ancestry regressed data
fclumpdir = "/data/clusterfs/lag/users/gokala/enigma-evol/eqtl/clumped_sumstats/european_hemave"

## Read sumstats paths
GWASsumstats=read.table(fGWASsumstats, header=FALSE)$V1
##Parse to get trait name
tmpname = sapply(GWASsumstats,function (x) {unlist(strsplit(x,"/",fixed=TRUE))[10]})
phenoname = paste0(sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[4]}),"_",sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[5]}))
clumpfileloc = dir(fclumpdir,pattern = ".clumped",full.names = T)

####################
##Full surface area
##Start by getting all the clumped SNPs only for full surface area
fullsurfind = which(phenoname=="surface_Full")
clump = read.table(clumpfileloc[fullsurfind],header=TRUE)
##Loop over all SNPs
for (i in 1:nrow(clump)) {
    ##Find all SNPs in LD (r2>0.6) with each clumped SNP
    system(paste0("module load plink/1.9b6 \
                  plink --bfile /data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",clump$SNP[i]," --out tmpld"))
    ##Read in the LD calculated from plink
    LD = read.table('tmpld.ld',header=TRUE)
    ##Turn into a genomic ranges object
    if (i==1) { 
       LDSNPs = GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A)
    } else {
      LDSNPs = c(GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A),LDSNPs)
    }
}

##Read in the BED file containing the fetal HGE file
annot = read.table(fannot,header=FALSE)
HGE = GRanges(gsub("chr","",annot$V1),IRanges(annot$V2,annot$V3),type=annot$V4)
##Restrict to only HGE_7PCW
##HGE = annot[which(annot$type=="HGE.7")]
##Use all HGEs for this analysis
#HGE = annot[which(annot$type=="HGE.12F" | annot$type=="HGE.12O" | annot$type=="HGE.7" | annot$type=="HGE.8")]
#HGE = reduce(HGE)
##Find overlaps
olap = findOverlaps(HGE,LDSNPs)
globalSASNPsolapHGE = unique(LDSNPs$indexSNP[subjectHits(olap)])
globalHGE = LDSNPs[which(!is.na(match(LDSNPs$indexSNP,globalSASNPsolapHGE)))]
##Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with HGE is: ',length(unique(globalHGE$indexSNP)),'\n')


##Overlap with eQTL data from psychENCODE
eqtl = read.table(feqtl,header=TRUE)
eqtl.GR = GRanges(gsub("chr","",eqtl$SNP_chr),IRanges(eqtl$SNP_start,eqtl$SNP_end))
mcols(eqtl.GR) = eqtl[,c(1:8,12:15)]
olap2 = findOverlaps(eqtl.GR,globalHGE)
cat('Number of HGE overlapping loci that also have an eQTL: ',length(unique(globalHGE$indexSNP[subjectHits(olap2)])),'\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])
##Remove the . annotation
eqtlgenes = sapply(eqtlgenes,function (x) {unlist(strsplit(x,".",fixed=TRUE))[1]})
write.table(eqtlgenes,"/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/european_hemave/fullSA_eqtlGenes_neanDep_european_hemave.txt",row.names = F,col.names = F,quote = F)
#eqtlgenes=read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/fullSA_eqtl_genes.txt",header = F)

##Convert these genes to hgnc_id
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org")
geneannot = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"),filters="ensembl_gene_id",values=eqtlgenes,mart=mart)
cat('These eQTLs impact ',length(which(geneannot=="protein_coding")),' protein-coding eGenes\n')
write.csv(geneannot,file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/european_hemave/fullSA_HGNCids_neanDep_european_hemave.csv",row.names=FALSE,quote=FALSE)
#smt_to_plot=read.csv("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/regionalSAfetalHGE.csv",header = T)

####################
##All regional surface areas
##Start by getting all the clumped SNPs only for any regional surface area
surfareaind = grep("surface",phenoname)
surfareaind = surfareaind[-7] # Find a more elegant way to exclude Full Surface area index here!
##Loop over all regions
for (j in 1:length(surfareaind)) {
    if (j==1) {
       clump = read.table(clumpfileloc[surfareaind[7]],header=TRUE)
    } else {
       if (file.exists(clumpfileloc[surfareaind[j]])) {
       	  clump = rbind(clump,read.table(clumpfileloc[surfareaind[j]],header=TRUE))
       }
    }
}

##Loop over all SNPs
for (i in 1:nrow(clump)) {
    ##Find all SNPs in LD (r2>0.6) with each clumped SNP
    system(paste0("module load plink/1.9b6 \
                  plink --bfile /data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",clump$SNP[i]," --out tmpld"))
    ##Read in the LD calculated from plink
    LD = read.table('tmpld.ld',header=TRUE)
    ##Turn into a genomic ranges object
    if (i==1) { 
       LDSNPs = GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A)
    } else {
      LDSNPs = c(GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A),LDSNPs)
    }
}

##Make sure that all index SNPs are in this list
LDSNPs = c(GRanges(clump$CHR,IRanges(clump$BP,clump$BP),SNP=clump$SNP,indexSNP=clump$SNP),LDSNPs)

save(LDSNPs,file="RegionalLDSNPs.Rdata")
##Read in the BED file containing the HGE_7PCW file
annot = read.table(fannot,header=FALSE)
HGE = GRanges(gsub("chr","",annot$V1),IRanges(annot$V2,annot$V3),type=annot$V4)
##Restrict to only HGE_7PCW
##HGE = annot[which(annot$type=="HGE.7")];
##Use all HGEs for this analysis
#HGE = annot[which(annot$type=="HGE.12F" | annot$type=="HGE.12O" | annot$type=="HGE.7" | annot$type=="HGE.8")]
#HGE = reduce(HGE);
##Find overlaps
olap = findOverlaps(HGE,LDSNPs)
regionalSASNPsolapHGE = unique(LDSNPs$indexSNP[subjectHits(olap)])
regionalHGE = LDSNPs[which(!is.na(match(LDSNPs$indexSNP,regionalSASNPsolapHGE)))]
##Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with an HGE is: ',length(unique(regionalHGE$indexSNP)),'\n')

##Overlap with eQTL data from psychENCODE
eqtl = read.table(feqtl,header=TRUE)
eqtl.GR = GRanges(gsub("chr","",eqtl$SNP_chr),IRanges(eqtl$SNP_start,eqtl$SNP_end))
mcols(eqtl.GR) = eqtl[,c(1:8,12:15)]
olap2 = findOverlaps(eqtl.GR,regionalHGE)
cat('Number of HGE overlapping loci that also have an eQTL: ',length(unique(regionalHGE$indexSNP[subjectHits(olap2)])),'\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])
##Remove the . annotation
eqtlgenes = sapply(eqtlgenes,function (x) {unlist(strsplit(x,".",fixed=TRUE))[1]})
##Convert these genes to hgnc_id
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org")
geneannot = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"),filters="ensembl_gene_id",values=eqtlgenes,mart=mart)
cat('These eQTLs impact ',length(which(geneannot=="protein_coding")),' protein-coding eGenes\n')
write.csv(geneannot,file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/european_hemave/regionalSA_HGNCids_neanDep_european_hemave.csv",row.names=FALSE,quote=FALSE)

##GO enrichment done with ggprolifer2, but can't get it to install here so doing it locally