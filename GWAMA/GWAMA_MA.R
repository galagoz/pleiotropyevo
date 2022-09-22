# This script will run Model Averaged GWAMA on dyslexia and rhythm GWAS sumstats.
#
# Gokberk Alagoz - 04.03.22
#
#-----
# Libraries
source("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/N_weighted_GWAMA.function.1_2_6.R?raw=TRUE")
library(data.table)
library(sp)
library(rater)
library(RcppArmadillo)
library(unmarked)
library(VGAM)
library(AICcmodavg)
library(metafor)
#-----
# Read sumstats that you prepared for Model Averaged GWAMA analysis

dys_MA = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia_reformatted_forModelAveraged3.txt",showProgress=F,data.table=F)
rhy_MA = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythym_reformatted_forModelAveraged3.txt",showProgress=F,data.table=F)

print("Read both sumstats succesfully.")

#-----
# Calculate Z-scores from Beta and Beta_SE for LDSC and GWAMA
# Use the formula for beta and beta_se from GWAMA wiki.

dys_MA$Z = dys_MA$Beta*(sqrt(dys_MA$N*2*dys_MA$EAF*(1-dys_MA$EAF)))
rhy_MA$Z = rhy_MA$Beta*(sqrt(rhy_MA$N*2*rhy_MA$EAF*(1-rhy_MA$EAF)))

dys_MA_b1_se1 <- dys_MA[c(2,11:12)]
rhy_MA_b2_se2 <- rhy_MA[c(2,11:12)]

M1 <- merge(dys_MA_b1_se1,rhy_MA_b2_se2,by=1)

print("Sumstats merged.")

M1 <- na.omit(M1)

A1 <- merge(dys_MA,M1, by="RS")
A1 <- A1[row.names(unique(A1["RS"])),]

B <- M1[,c(2,4)]
SE <- M1[,c(3,5)]

# run ldsc gencor, create a CTI matrix, then proceed with this script.

cov_Z = as.matrix(read.table("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gencor/CTI_matrix.txt",header=T))

# list of SNP-heritabilities
h2 <- as.vector(c(sqrt(0.1465), sqrt(0.1459))) # use the h2 estimates based on 5% prevelance for dys and 5% for atypical beat sync.

#The following for loop will perform the MA GWAMA analysis

output_SNP = as.data.frame(matrix(, nrow = 2, ncol = 4))

print("Starting to run GWAMA_MA.")

for(i in 1:nrow(B)){

  V <- diag(SE[i,]) %*% cov_Z %*% diag(SE[i,])


  yi <- as.numeric(B[i,])


  grid <- expand.grid(c(0,1),c(0,1))

  grid <- t(grid)

  m31 <- rma.mv(yi~0+h2, V=V,  method="ML")

  # 2 effects
  m32 <- rma.mv(yi~0+h2 + I(grid[,2]*h2), V=V, method="ML")
  m33 <- rma.mv(yi~0+h2 + I(grid[,3]*h2), V=V, method="ML")
  m34 <- rma.mv(yi~0+h2 + I(grid[,4]*h2), V=V, method="ML")

  Mods <- list(m31,m32,m33,m34)

  K <- c(rep(1,1),rep(2,3))

  LL <- lapply(Mods,logLik)
  LL <- unlist(LL)
  aictabCustom(LL,K,nobs=2)

  XX <- sapply(Mods,predict)

  est2 <- matrix(unlist(XX[1,1]),2,1)
  est3 <- matrix(unlist(XX[1,2]),2,3)

  est<- cbind(est2,est3)

  se2 <- matrix(unlist(XX[2,1]),2,1)
  se3 <- matrix(unlist(XX[2,2]),2,3)

  se<- cbind(se2,se3)

  y1 <- modavgCustom(LL,K,nobs=2, estimate=est[1,] ,se=se[1,],second.ord = T)
  y2 <- modavgCustom(LL,K,nobs=2, estimate=est[2,] ,se=se[2,],second.ord = T)

  output_SNP[i, ] <- c(y1$Mod.avg.est,y1$Uncond.SE,y2$Mod.avg.est,y2$Uncond.SE)


}

print("GWAMA_MA ended, making the output file.")

Output <- cbind(A1[,c(1,2,3,4)],output_SNP,pchisq((output_SNP[,1]/output_SNP[,2])^2,1,lower=F),
                pchisq((output_SNP[,3]/output_SNP[,4])^2,1,lower=F))

names(Output) <- c("RS","cptid","CHR","BP","A1","A2","Beta_Dys","SE_Dys","Beta_Rhy","SE_Rhy","P_Dys","P_Rhy")

write.table("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/MA/AVGWAS.txt" ,quote=F,row.names=F,col.names=T)

print("Output file is saved, go check your results!")

