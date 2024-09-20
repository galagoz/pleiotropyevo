#/bin/Rscript
#
# This script will run Genomic SEM on example data
# and will generate Q indeces.

# Prep
library(GenomicSEM)
# load LDSC
load("/path/to/dir/GenomicSEM_LDSCoutput_dys_rhyimp.RData") # here set the path

# 1) Common Pathway Model
model = "F1 =~ 1*Dyslexia + 1*Rhythm_impairment
F1 ~ SNP
Dyslexia ~~ a*Dyslexia
Rhythm_impairment ~~ a*Rhythm_impairment"

dys_and_rhyimp_sumstats_firstLine = data.frame(SNP = as.character("rs79373928"),
                 CHR = as.numeric(1),
                 BP = as.numeric(801536),
                 MAF = as.numeric(0.0139165),
                 A1 = as.character("T"),
                 A2 = as.character("G"),
                 beta.Dyslexia = as.numeric(0.023621718),
                 se.Dyslexia = as.numeric(0.017975108),
                 beta.Rhythm_impairment = as.numeric(-0.0603634880),
                 se.Rhythm_impairment = as.numeric(0.027779895))

# Run the multivariate GWAS with Common Pathways Model
CorrelatedFactors = userGWAS(covstruc = LDSCoutput_2traits,
                             SNPs = dys_and_rhyimp_sumstats_firstLine,
                             estimation = "DWLS",
                             model = model,
                             printwarn = TRUE,
                             sub=c("F1~SNP"),
                             toler = FALSE,
                             SNPSE = FALSE,
                             parallel = TRUE,
                             GC="standard",
                             MPI=FALSE,
                             smooth_check=FALSE)

CorrelatedFactors

# 2) Independent Pathways Model
model2 = "F1 =~ 1*Dyslexia + 1*Rhythm_impairment 
Dyslexia + Rhythm_impairment ~ SNP
Dyslexia ~~ a*Dyslexia
Rhythm_impairment ~~ a*Rhythm_impairment"

# Run the multivariate GWAS with Independent Pathways Model
orrelatedFactors2 = userGWAS(covstruc = LDSCoutput_2traits,
                               SNPs = dys_and_rhyimp_sumstats_firstLine,
                               estimation = "DWLS", 
                               model = model2, 
                               printwarn = TRUE, 
                               sub=c("Dyslexia ~ SNP", "Rhythm_impairment ~ SNP"), 
                               toler = FALSE, 
                               SNPSE = FALSE, 
                               parallel = TRUE)

CorrelatedFactors2[[1]]
CorrelatedFactors2[[2]]

# 3) Compute Q-scores (heterogeneity index)
# Calculate the chi-square for each SNP and save it in a vector called Q_chisq
Q_chisq = CorrelatedFactors[[1]]$chisq - CorrelatedFactors2[[1]]$chisq

# Calculate the df associated with this chi-square statistic and name it df. Note that only the first run
# is used to calculate df, as the degrees of freedom will be the same for every SNP.
df = CorrelatedFactors[[1]]$chisq_df[1] - CorrelatedFactors2[[1]]$chisq_df[1]

# Calculate the p-value associated with the Q-statistic, using the relevant degrees of freedom calculated
# in the previous. This saves the p-values in a vector called Q_chisq_pval, and can easily be appended to 
# either CorrelatedFactors[[1]] or CorrelatedFactors[[2]] (i.e., using cbind). 
Q_chisq_pval = pchisq(Q_chisq, df, lower.tail = FALSE)

CorrelatedFactors2[[1]]$Q_chisq = Q_chisq
CorrelatedFactors2[[1]]$Q_chisq_pval = Q_chisq_pval
