#-----------------------------------------------------------------------------------------------
# script for the data analysis of GWAS meta-analysis results using GenomicSEM
# written by Else Eising, 6 November 2020 - last updated by Gokberk Alagoz, 16 June 2023
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
require(GenomicSEM)

inDir = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/"
outDir = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"
LDSCDir = "/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/"
hm3 = "/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/w_hm3.snplist"

#----------------------------------------------------------------------------------------------------------
# 1) munge the sumstats
#----------------------------------------------------------------------------------------------------------
setwd("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM")
sumstats = c(paste0(inDir, "dyslexia_forgenomicSEM.txt"), paste0(inDir, "rhythm_impairment_forgenomicSEM.txt"))
trait_names = c("Dyslexia","Rhythm_impairment")
#sample_sizes = c(169510, 187407) # these are effective population sizes, calculated as 4.N_cases.(1-N_cases/N_total)

#munge(sumstats, hm3, trait.names = trait_names)

#----------------------------------------------------------------------------------------------------------
# 2) Run multivariable LDSC
#----------------------------------------------------------------------------------------------------------

traits = c("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/Dyslexia.sumstats.gz", 
           "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/Rhythm_impairment.sumstats.gz")

sample.prev = c(0.045, 0.085)
population.prev = c(0.05, 0.0475)
trait_names = c("Dyslexia", "Rhythm_impairment")

#LDSCoutput_2traits = ldsc(traits, sample.prev, population.prev, LDSCDir, LDSCDir, trait_names)

##optional command to save the ldsc output in case you want to use it in a later R session. 
#save(LDSCoutput_2traits, file = paste0(outDir, "GenomicSEM_LDSCoutput_dys_rhyimp.RData"))
# to load ldsc output
load(paste0(outDir, "GenomicSEM_LDSCoutput_dys_rhyimp.RData"))

#----------------------------------------------------------------------------------------------------------
# 3) Common Factor Model (Measurement Model)
#----------------------------------------------------------------------------------------------------------

model = "F1 =~ 1*Dyslexia + 1*Rhythm_impairment"
#run the multivariate GWAS using parallel processing
CommonFactor2_DWLS = usermodel(LDSCoutput_2traits, estimation="DWLS", model = model)

# [1] "Running primary model"
# [1] "Calculating CFI"
# [1] "Calculating Standardized Results"
# [1] "Calculating SRMR"
# elapsed 
# 7.784 
# [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"
# $modelfit
# chisq df p_chisq AIC CFI SRMR
# df    NA  0      NA  NA  NA   NA
# 
# $results
# lhs op               rhs Unstand_Est          Unstand_SE STD_Genotype    STD_Genotype_SE   STD_All      p_value
# 1          Dyslexia ~~          Dyslexia  0.11002123 0.00611191127152301    0.7175451 0.0407745118456241 0.7175451 1.909290e-72
# 5 Rhythm_impairment ~~ Rhythm_impairment  0.09762516 0.00609435901716243    0.7175451 0.0433674050247672 0.7175451 9.424537e-58
# 4                F1 ~~                F1  0.04081673 0.00387841122666485    0.2824549 0.0268389034118519 0.2824549 6.690829e-26
# 2                F1 =~          Dyslexia  1.00000000                        1.0000000                    0.5314649           NA
# 3                F1 =~ Rhythm_impairment  1.00000000                        1.0000000                    0.5314649           NA

#----------------------------------------------------------------------------------------------------------
# 4) Prepare the summary statistics for GWAS
#----------------------------------------------------------------------------------------------------------

se.logit = c(T,T)
ref = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/reference.1000G.maf.0.005.txt"

# dys_and_rhyimp_sumstats = sumstats(files = sumstats,
#                                    ref = ref,
#                                    trait.names = trait_names,
#                                    se.logit = se.logit)

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
F1 ~ SNP"

# Run the multivariate GWAS using parallel processing
CorrelatedFactors = userGWAS(covstruc = LDSCoutput_2traits,
                             SNPs = dys_and_rhyimp_sumstats[dys_and_rhyimp_sumstats$SNP=="rs17151639",],
                             estimation = "DWLS",
                             model = model,
                             printwarn = TRUE,
                             sub=c("F1 ~~ F1","Dyslexia ~~ Dyslexia","Rhythm_impairment ~~ Rhythm_impairment", "F1 ~ SNP"),
                             #cores = 5,
                             toler = FALSE,
                             SNPSE = FALSE,
                             parallel = TRUE,
                             GC="standard",
                             MPI=FALSE,
                             smooth_check=FALSE)

## optional command to save the multivariate GWAS results in case you want to use it in a later R session.
#save(CorrelatedFactors, file = paste0(outDir, "GenomicSEM_multivarGWAS_CPM_dys_rhyimp.RData"))
#fwrite(CorrelatedFactors[[1]], file = paste0(outDir,  "GenomicSEM_multivarGWAS_CPM_dys_rhyimp.tab"), 
#       sep = "\t",  row.names = FALSE, col.names = TRUE)
# to load multivarGWAS results
#load(paste0(outDir, "GenomicSEM_multivarGWAS_dys_rhyimp.RData"))

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

#fwrite(CorrelatedFactors[[1]][,c(1,4,5,6,14,15,22)], file = paste0(outDir,  "GenomicSEM_multivarGWAS_CPM_dys_rhyimp_forMunge.tab"), 
#       sep = "\t",  row.names = FALSE, col.names = TRUE)

############################################################################

#----------------------------------------------------------------------------------------------------------
# 6) Run a follow-up model to obtain heterogeneity (Q) index.
#----------------------------------------------------------------------------------------------------------
#Step 1a: specify the Independent Pathways Model
model2 <- "F1 =~ 1*Dyslexia + 1*Rhythm_impairment 
Dyslexia + Rhythm_impairment ~ SNP"

#Dyslexia ~~ a*Dyslexia
#Rhythm_impairment ~~ a*Rhythm_impairment"

#Step 2b: run the Step 2 Model using parallel processing
CorrelatedFactors2 <- userGWAS(covstruc = LDSCoutput_2traits,
                               SNPs = dys_and_rhyimp_sumstats[dys_and_rhyimp_sumstats$SNP=="rs12640404",],
                               estimation = "DWLS", 
                               model = model2, 
                               printwarn = TRUE, 
                               sub=c("Dyslexia ~~ Dyslexia","Rhythm_impairment ~~ Rhythm_impairment", "Dyslexia ~ SNP", "Rhythm_impairment ~ SNP", "F1 ~~ F1"),
                               #ub=c("Dyslexia ~ SNP", "Rhythm_impairment ~ SNP"), 
                               #cores = 5,
                               toler = FALSE, 
                               SNPSE = FALSE, 
                               parallel = TRUE)

##optional command to save the multivariate GWAS results in case you want to use it in a later R session.
#save(CorrelatedFactors2, file = paste0(outDir, "GenomicSEM_multivarGWAS_IPM_dys_rhyimp.RData"))

#fwrite(CorrelatedFactors2[[1]], file = paste0(outDir,  "GenomicSEM_multivarGWAS_dys_IPM_dys_rhyimp.tab"), 
#       sep = "\t",  row.names = FALSE, col.names = TRUE)

#fwrite(CorrelatedFactors2[[2]], file = paste0(outDir,  "GenomicSEM_multivarGWAS_rhyimp_IPM_dys_rhyimp.tab"), 
#       sep = "\t",  row.names = FALSE, col.names = TRUE)
# to load multivarGWAS results
load(paste0(outDir, "GenomicSEM_multivarGWAS_IPM_dys_rhyimp.RData"))

N_hat_F1<-mean(1/((2*CorrelatedFactors2[[1]]$MAF*(1-CorrelatedFactors2[[1]]$MAF))*CorrelatedFactors2[[1]]$SE^2), na.rm = T)
CorrelatedFactors2[[1]]$N = round(N_hat_F1)

#fwrite(CorrelatedFactors2[[1]][,c(1,4,5,6,14,15,22)], file = paste0(outDir,  "GenomicSEM_multivarGWAS_dys_IPM_forMunge.tab"), 
#       sep = "\t",  row.names = FALSE, col.names = TRUE)

N_hat_F1<-mean(1/((2*CorrelatedFactors2[[2]]$MAF*(1-CorrelatedFactors2[[2]]$MAF))*CorrelatedFactors2[[2]]$SE^2), na.rm = T)
CorrelatedFactors2[[2]]$N = round(N_hat_F1)

#fwrite(CorrelatedFactors2[[2]][,c(1,4,5,6,14,15,22)], file = paste0(outDir,  "GenomicSEM_multivarGWAS_rhyimp_IPM_forMunge.tab"), 
#       sep = "\t",  row.names = FALSE, col.names = TRUE)

#----------------------------------------------------------------------------------------------------------
# 7) Plot
#----------------------------------------------------------------------------------------------------------
# Read lambdaGC corrected CPM results

dys_rhy_genSEM_CPM_raw = fread(paste0(outDir, "GenomicSEM_multivarGWAS_CPM_dys_rhyimp.tab"))
dys_rhy_genSEM_CPM = fread(paste0(outDir, "GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr.tab"))

# QQplot and Manhattan plot of CPM results

# before GC correction
png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/genSEM_dys_rhyimp_CPM_QQplot_beforeLambdaCorrection.png", type="cairo")
qq(dys_rhy_genSEM_CPM_raw$Pval_Estimate)
dev.off()

# after GC correction
#png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/genSEM_dys_rhyimp_CPM_QQplot_afterLambdaCorrection.png", type="cairo")
qq(dys_rhy_genSEM_CPM$Pval_Estimate)
#dev.off()

png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/genSEM_dys_rhyimp_CPM_Manhattan.png", type="cairo", width=1000, height=500, units="mm", res = 300)
manhattan(dys_rhy_genSEM_CPM, ylim = c(0,15), chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP", col = c("#ED6B06", "#00786A"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
dev.off()

# miami plot v1
png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/genSEM_dys_rhyimp_Qscores_Miamiplot.png", type="cairo", width=1000, height=500, units="mm", res = 300)
par(mfrow=c(2,1))
par(mar=c(0,5,5,3))
manhattan(dys_rhy_genSEM_CPM, ylim = c(0,15), chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP", col = c("#ED6B06", "#00786A"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
par(mar=c(5,5,1,3))
manhattan(dys_rhy_genSEM_CPM, ylim = c(15,0), chr = "CHR", bp = "BP", p = "chisq_pval", snp = "SNP", col = c("#ED6B06", "#00786A"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
dev.off()

# miami plot v2
# MPI orange = #ED6B06
# MPI green = #00786A
png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/genSEM_dys_rhyimp_Qscores_Miamiplot_v2.png", type="cairo", width=1000, height=500, units="mm", res = 300)
par(mfrow=c(2,1))
par(mar=c(0,5,5,3))
manhattan(dys_rhy_genSEM_CPM, ylim = c(0,20), chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP", col = c("#ED6B06", "#00786A"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
par(mar=c(5,5,1,3))
manhattan(dys_rhy_genSEM_CPM, ylim = c(20,0), chr = "CHR", bp = "BP", p = "chisq_pval", snp = "SNP", col = c("#ED6B06", "#00786A"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
dev.off()

# Plot CPM, Q, dys and rhy manhattans at once
dyslexia_sumstats<-fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia_forgenomicSEM.txt",showProgress=F,data.table=F)
rhythm_sumstats<-fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythm_impairment_forgenomicSEM.txt",showProgress=F,data.table=F)

dyslexia_sumstats$CHR[dyslexia_sumstats$CHR == "X"] = 23
rhythm_sumstats$CHR[rhythm_sumstats$CHR == "X"] = 23
dyslexia_sumstats$CHR = as.integer(dyslexia_sumstats$CHR)
rhythm_sumstats$CHR = as.integer(rhythm_sumstats$CHR)

png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/dys_rhy_CPM_Q_Manhattans_v4.png", type="cairo", width=700, height=500, units="mm", res = 300)
par(mfrow=c(4,1))
manhattan(dyslexia_sumstats, ylim = c(0,20), chr = "CHR", bp = "BP", p = "P", snp = "SNPID", col = c("#3182bd", "#9ecae1"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
manhattan(rhythm_sumstats, ylim = c(0,20), chr = "CHR", bp = "BP", p = "P", snp = "SNPID", col = c("#31a354", "#a1d99b"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
manhattan(dys_rhy_genSEM_CPM, ylim = c(0,15), chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP", col = c("#e6550d", "#fdae6b"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
manhattan(dys_rhy_genSEM_CPM, ylim = c(0,15), chr = "CHR", bp = "BP", p = "chisq_pval", snp = "SNP", col = c("#756bb1", "#bcbddc"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
dev.off()

# miami plot dys rhy
png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/dys_rhyimp_Miamiplot_v2.png", type="cairo", width=1000, height=500, units="mm", res = 300)
par(mfrow=c(2,1))
par(mar=c(0,5,5,3))
manhattan(dyslexia_sumstats, ylim = c(0,15), chr = "CHR", bp = "BP", p = "P", snp = "SNPID", col = c("#214D9D", "#90A5CE"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
par(mar=c(5,5,1,3))
manhattan(rhythm_sumstats, ylim = c(15,0), chr = "CHR", bp = "BP", p = "P", snp = "SNPID", col = c("#EEA300", "#F7D27F"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
dev.off()

# Manhattan plots w/CommonFactor CPASSOC and GWAMA results
dys_rhy_Nweighted = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/GWAMA/dys_rhyimp.N_weighted_GWAMA.results.txt")
dys_rhy_commonSig = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/CPASSOC/SHom_Dys_RhyImp_results.txt")

png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/CPM_CPASSOC_GWAMA_Manhattans.png", type="cairo", width=700, height=500, units="mm", res = 300)
par(mfrow=c(3,1))
manhattan(dys_rhy_genSEM_CPM, ylim = c(0,15), chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP", col = c("#e6550d", "#fdae6b"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
manhattan(dys_rhy_Nweighted, ylim = c(0,25), chr = "CHR", bp = "BP", p = "PVAL", snp = "SNPID", col = c("gray20", "gray44"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
manhattan(dys_rhy_commonSig, ylim = c(0,25), chr = "CHR", bp = "BP", p = "p.shom", snp = "Row.names", col = c("gray20", "gray44"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
dev.off()

# Example SEM diagrams for significantly heterogeneous SNPs.
# Make these plot for top 2 most hetero- and most homo-effect SNPs.

dys_rhy_genSEM_CPM_Qpvalranked = dys_rhy_genSEM_CPM[order(chisq_pval),]
dys_rhy_genSEM_CPM_Qranked = dys_rhy_genSEM_CPM[order(chisq),]

# most hetero: rs17151639, rs9659996
# most homo: rs448036, rs12640404

#--------------------------------------------
library("scatterplot3d") # load

for3dplot = dys_rhy_genSEM_CPM[,c(4,12,16,18)]
for3dplot$Qsig = "grey"
for3dplot[for3dplot$chisq_pval<5e-8,]$Qsig = "green"
colors=for3dplot$Qsig
#scatterplot3d(for3dplot[,c(1,2,3)])

graphics.off()
png("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/genSEM_MAF_BETA_Q_3Dplot_angle30.png",
    width = 8, height = 6, units = "in", res = 600)
par(mai = c(0.5, 0.5, 0.5, 0.5))
s3d = scatterplot3d(x = for3dplot$MAF, y = for3dplot$est, z = for3dplot$chisq, ylim = c(-0.20,0.20),
                    cex.symbols = 1, pch = 20, color = colors, angle = 30)
dev.off()

# 1. Source the function
source('~/hubiC/Documents/R/function/addgrids3d.r')
# 2. Empty 3D scatter plot using pch=""
s3d <- scatterplot3d(dys_rhy_genSEM_CPM[,c(4,12,16)], pch = "", grid=FALSE, box=FALSE)
# 3. Add grids
addgrids3d(dys_rhy_genSEM_CPM[,c(4,12,16)], grid = c("xy", "xz", "yz"))
# 4. Add points
s3d$points3d(dys_rhy_genSEM_CPM[,c(4,12,16)], pch = 16)

#----------------------------------------------------------------------------------------------------------