#-----------------------------------------------------------------------------------------------
# script for the data analysis of GWAS meta-analysis results using GenomicSEM
# written by Else Eising, 6 November 2020 - last updated by Gokberk Alagoz, 05 August 2022
# See https://github.com/MichelNivard/GenomicSEM/wiki/4.-Common-Factor-GWAS and https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
# Steps
#-----------------------------------------------------------------------------------------------

# 1) Munge the sumstats
# 2) Run multivariable LDSC
# 3) Run common factor GWASs
# 4) Fit 2 factor models with genomicSEM

#-----------------------------------------------------------------------------------------------
# Load libraries and set directories
#-----------------------------------------------------------------------------------------------

library(data.table)
library(qqman)
#library(devtools)
#install_github("MichelNivard/GenomicSEM")
library(GenomicSEM)
library(lavaan)
require(stats)
require(Matrix)

inDir = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"
outDir = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"
LDSCDir = "/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/"
hm3 = "/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/w_hm3.snplist"

#----------------------------------------------------------------------------------------------------------
# 1) munge the sumstats
#----------------------------------------------------------------------------------------------------------

sumstats = c(paste0(inDir, "dyslexia_reformatted_forModelAveraged3.txt"), paste0(inDir, "rhythym_impairment_reformatted_forModelAveraged5.txt"))
trait_names = c("Dyslexia","Rhythm_impairment")
sample_sizes = c(169510, 187407) # these are effective population sizes, calculated as 4.N_cases.(1-N_cases/N_total)

munge(sumstats, hm3, trait.names = trait_names, N = sample_sizes)

#----------------------------------------------------------------------------------------------------------
# 2) Run multivariable LDSC
#----------------------------------------------------------------------------------------------------------

traits = c("Dyslexia.sumstats.gz", "Rhythm_impairment.sumstats.gz")
sample.prev = c(0.045, 0.085)
population.prev = c(0.05, 0.0475)
trait_names = c("Dyslexia", "Rhythm_impairment")
LDSCoutput_2traits = ldsc(traits, sample.prev, population.prev, LDSCDir, LDSCDir, trait_names)

##optional command to save the ldsc output in case you want to use it in a later R session. 
#save(LDSCoutput_2traits, file = paste0(outDir, "GenomicSEM_LDSCoutput_dys_rhyimp.RData"))
# to load ldsc output
load(paste0(outDir, "GenomicSEM_LDSCoutput_dys_rhyimp.RData"))

#----------------------------------------------------------------------------------------------------------
# 3) Old Common Factor Model
#----------------------------------------------------------------------------------------------------------

model = "F1 =~ a*Dyslexia + a*Rhythm_impairment
F1 ~ SNP"
#run the multivariate GWAS using parallel processing
CommonFactor2_DWLS = usermodel(LDSCoutput_2traits, estimation="DWLS", model = model)

#[1] "Running primary model"
#[1] "Calculating CFI"
#[1] "Calculating Standardized Results"
#[1] "Calculating SRMR"
#elapsed 
#0.528 
#[1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"
#[1] "Please note that when equality constraints are used in the current version of Genomic SEM that the standardized output will also impose the same constraint."
#> CommonFactor2_DWLS
#$modelfit
#chisq df p_chisq AIC CFI SRMR
#df    NA  0      NA  NA  NA   NA

#$results
#lhs op               rhs Unstand_Est         Unstand_SE STD_Genotype    STD_Genotype_SE   STD_All
#1          Dyslexia ~~          Dyslexia   0.8246293 0.0409209002331068    0.7167753 0.0411889231296829 0.7167753
#5 Rhythm_impairment ~~ Rhythm_impairment   0.2563430 0.0233969529962123    0.7167753  0.044014456010507 0.7167753
#4                F1 ~~                F1   0.1908639 0.0174246241862993    0.2832247 0.0258565530578309 0.2832247
#2                F1 =~          Dyslexia   1.0000000                       1.0000000                    0.5321886
#3                F1 =~ Rhythm_impairment   1.0000000                       1.0000000                    0.5321886

#p_value
#1 2.595675e-90
#5 6.201168e-28
#4 6.379463e-28
#2           NA
#3           NA

#----------------------------------------------------------------------------------------------------------
# 4) Prepare the summary statistics for GWAS
#----------------------------------------------------------------------------------------------------------

se.logit = c(T,T)
ref = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/reference.1000G.maf.0.005.txt"

dys_and_rhyimp_sumstats = sumstats(files = sumstats, 
                                   ref = ref, 
                                   trait.names = trait_names, 
                                   se.logit = se.logit, 
                                   N = sample_sizes)

##optional command to save the formatted sumstats in case you want to use it in a later R session. 
#save(dys_and_rhyimp_sumstats, file = paste0(outDir, "GenomicSEM_sumstats_dys_rhyimp.RData"))
# to load sumstats
load(paste0(outDir, "GenomicSEM_sumstats_dys_rhyimp.RData"))

#----------------------------------------------------------------------------------------------------------
# 5) Combine the summary statistics and LDSC output and run the user-specified GWAS - This one is a very
# memory intensive step. So run it either on the grid or limit memory use on lux13/14 while running.
# For this you should run rstudio on lux13/14 as following:
# module load R/R-4.0.3 rstudio/1.1.436
# ( ulimit -v 90000000; rstudio )
# This way, you will allocate rstudio 90GB of RAM (memory) (lux13 has 220GB) and will not interrupt others'
# jobs.
#----------------------------------------------------------------------------------------------------------

# Specify the Common Pathways Model
model = "F1 =~ 1*Dyslexia + 1*Rhythm_impairment
F1 ~ SNP
Dyslexia ~~ a*Dyslexia
Rhythm_impairment ~~ a*Rhythm_impairment"

#run the multivariate GWAS using parallel processing
CorrelatedFactors = userGWAS(covstruc = LDSCoutput_2traits,
                             SNPs = dys_and_rhyimp_sumstats,
                             estimation = "DWLS",
                             model = model,
                             printwarn = TRUE,
                             sub=c("F1~SNP"),
                             cores = 5,
                             toler = FALSE,
                             SNPSE = FALSE,
                             parallel = TRUE,
                             GC="standard",
                             MPI=FALSE,
                             smooth_check=FALSE)

##optional command to save the multivariate GWAS results in case you want to use it in a later R session.
#save(CorrelatedFactors, file = paste0(outDir, "GenomicSEM_multivarGWAS_dys_rhyimp_v2_secondrun.RData"))
# to load multivarGWAS results
load(paste0(outDir, "GenomicSEM_multivarGWAS_dys_rhyimp_v2.RData"))

#----------------------------------------------------------------------------------------------------------
# 5.5) Calculate sample size for factors. 
#----------------------------------------------------------------------------------------------------------

# Before munging Genomic SEM meta-analysed sumstats, add the N column and get rid off irrelevant columns as
# shown below.

# Calculate Expected Sample Size for Factor 1
# recommendation from GenomicSEM Wiki: restrict to MAF of <40% and >10% to have more stable estimates, but this is not required,
# so I didn't do it.
# CorrelatedFactors1<-subset(CorrelatedFactors[[1]], CorrelatedFactors[[1]]$MAF <= .4 & CorrelatedFactors[[1]]$MAF >= .1)

N_hat_F1<-mean(1/((2*CorrelatedFactors[[1]]$MAF*(1-CorrelatedFactors[[1]]$MAF))*CorrelatedFactors[[1]]$SE^2), na.rm = T)
CorrelatedFactors[[1]]$N = round(N_hat_F1)

fwrite(CorrelatedFactors[[1]][,c(1,4,5,6,14,15,22)], file = paste0(outDir,  "GenomicSEM_multivarGWAS_dys_rhyimp_v2_forMunge.tab"), 
       sep = "\t",  row.names = FALSE, col.names = TRUE)

####Make a Miamiplot w/CommonFactor and GWAMA results#######################

dys_rhy_Nweighted = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/NGWAMA_beatFlippedZ.N_weighted_GWAMA.results.txt")

#png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/Dys_Rhy_GWAMA-gSEM_Miamiplot.png", type="cairo", width=2000, height=500)
#par(mfrow=c(2,1))
#par(mar=c(0,5,5,3))
#manhattan(CorrelatedFactors[[1]], ylim = c(0,25), chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP", col = c("dodgerblue3", "grey80"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
#par(mar=c(5,5,1,3))
#manhattan(dys_rhy_Nweighted, ylim = c(25,0), chr = "CHR", bp = "BP", p = "PVAL", snp = "SNP", col = c("dodgerblue3", "grey80"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
#dev.off()

#png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/Dys_Rhyasync_gSEM_QQplot.png", type="cairo")
#qq(CorrelatedFactors[[1]]$Pval_Estimate)
#dev.off()

############################################################################

#----------------------------------------------------------------------------------------------------------
# 6) Run a follow-up model to obtain heterogeneity (Q) index.
#----------------------------------------------------------------------------------------------------------
#Step 1a: specify the Independent Pathways Model
model2 <- "F1 =~ 1*Dyslexia + 1*Rhythm_impairment 
Dyslexia + Rhythm_impairment ~ SNP
Dyslexia ~~ a*Dyslexia
Rhythm_impairment ~~ a*Rhythm_impairment"

# "F1 =~ a*Dyslexia + a*Rhythm_impairment
# Dyslexia + Rhythm_impairment ~ SNP"

#Step 2b: run the Step 2 Model using parallel processing
CorrelatedFactors2 <- userGWAS(covstruc = LDSCoutput_2traits,
                               SNPs = dys_and_rhyimp_sumstats[1000000:1000010,],
                               estimation = "DWLS", 
                               model = model2, 
                               printwarn = TRUE, 
                               sub=c("Dyslexia ~ SNP", "Rhythm_impairment ~ SNP"), 
                               cores = 5,
                               toler = FALSE, 
                               SNPSE = FALSE, 
                               parallel = TRUE)

##optional command to save the multivariate GWAS results in case you want to use it in a later R session.
#save(CorrelatedFactors2, file = paste0(outDir, "GenomicSEM_multivarGWAS_model2_dys_rhyimp_v2.RData"))
fwrite(CorrelatedFactors2[[1]], file = paste0(outDir,  "GenomicSEM_multivarGWAS_model2_dys_rhyimp_v2.tab"), 
       sep = "\t",  row.names = FALSE, col.names = TRUE)
# to load multivarGWAS results
#load(paste0(outDir, "GenomicSEM_multivarGWAS_model2_dys_rhyimp_v2.RData"))

#----------------------------------------------------------------------------------------------------------
# 7) Plot manhattan plots
#----------------------------------------------------------------------------------------------------------
dyslexia_sumstats<-fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia_reformatted_forNweighted3.txt",showProgress=F,data.table=F)
rhythm_sumstats<-fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythym_reformatted_forNweighted5.txt",showProgress=F,data.table=F)

dyslexia_sumstats$CHR[dyslexia_sumstats$CHR == "X"] = 23
rhythm_sumstats$CHR[rhythm_sumstats$CHR == "X"] = 23
dyslexia_sumstats$CHR = as.integer(dyslexia_sumstats$CHR)
rhythm_sumstats$CHR = as.integer(rhythm_sumstats$CHR)

png("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/gSEM_Dys_Rhy_Manhattans_v2.png", type="cairo", width = 900, height = 750)
par(mfrow=c(3,1))
manhattan(CorrelatedFactors[[1]], ylim = c(0,25), chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP", col = c("dodgerblue3", "grey80"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
manhattan(dyslexia_sumstats, ylim = c(0,25), chr = "CHR", bp = "BP", p = "P", snp = "SNPID", col = c("dodgerblue3", "grey80"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
manhattan(rhythm_sumstats, ylim = c(0,25), chr = "CHR", bp = "BP", p = "P", snp = "SNPID", col = c("dodgerblue3", "grey80"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
dev.off()