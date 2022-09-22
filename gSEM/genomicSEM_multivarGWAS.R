#----------------------------------------------------------------------------------------------------------
# 5) Combine the summary statistics and LDSC output and run the user-specified GWAS - This one is a very
# memory intensive step, thus RUN THIS STEP ON THE GRID!
#----------------------------------------------------------------------------------------------------------

library(devtools)
#install_github("MichelNivard/GenomicSEM")
library(GenomicSEM)
require(stats)
require(Matrix)

outDir = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/"

load(paste0(outDir, "GenomicSEM_LDSCoutput_dys_rhyimp.RData"))
model = "F1 =~ a*Dyslexia + a*Rhythm_impairment"
#run the multivariate GWAS using parallel processing
CommonFactor2_DWLS = usermodel(LDSCoutput_2traits, estimation="DWLS", model = model)
load(paste0(outDir, "GenomicSEM_sumstats_dys_rhyimp.RData"))

#specify the model
model = "F1 =~ a*Dyslexia + a*Rhythm_impairment
F1 ~ SNP"

#run the multivariate GWAS using parallel processing
CorrelatedFactors = userGWAS(covstruc = LDSCoutput_2traits,
                             SNPs = dys_and_rhyimp_sumstats,
                             estimation = "DWLS",
                             model = model,
                             printwarn = TRUE,
                             sub=c("F1~SNP"),
                             cores = NULL,
                             toler = FALSE,
                             SNPSE = FALSE,
                             parallel = TRUE,
                             GC="standard",
                             MPI=FALSE,
                             smooth_check=FALSE)

