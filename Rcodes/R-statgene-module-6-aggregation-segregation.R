############################################################
## Study Heritbality
 
############################################################
## Assume a normal additive model of Y = alpha + beta*G + e
## G is coded as 0, 1 and 2 copies of the minor allele
## Study the variance decomposition: V_Y=V_Phenotype=V_G+V_E
############################################################
# split the window to 3 rows and 2 columns
par(mfrow=c(3,2))

# parameter values
# sample size
n=1000 
# regression parameters
alpha=1
beta=2
sigma=1
# MAF
p=0.2

## G component
## obtain genotype counts assuming HWE
nG=rmultinom(1,size=n,prob=c((1-p)^2,2*p*(1-p), p^2))
print(nG)
# obtain the actual genotype vector of 0, 1 and 2
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))

## Y component
Y=alpha+beta*G+rnorm(n,mean=0,sd=sigma)

## Plot Y vs. G
plot(G,Y,xlim=c(-0.5,2),ylim=c(-2,8), cex=0.7, main=paste("beta=",beta, ", sigma=", sigma,sep=""))
title(line=0.5,paste("Y = alpha + beta G + e"),cex.main=0.8)
hist(Y,xlim=c(-2,8),nclass=20)

# change sigma value
sigma=0.01
Y=alpha+beta*G+rnorm(n,mean=0,sd=sigma)
plot(G,Y,xlim=c(-0.5,2),ylim=c(-2,8), cex=0.7, main=paste("beta=",beta, ", sigma=", sigma,sep=""))
title(line=0.5,paste("Y = alpha + beta G + e"),cex.main=0.8)
hist(Y,xlim=c(-2,8),nclass=20)

# change beta value
beta=0.02
sigma=1
Y=alpha+beta*G+rnorm(n,mean=0,sd=sigma)
plot(G,Y,xlim=c(-0.5,2),ylim=c(-2,8), cex=0.7, main=paste("beta=",beta, ", sigma=", sigma,sep=""))
title(line=0.5,paste("Y = alpha + beta G + e"),cex.main=0.8)
hist(Y,xlim=c(-2,8),nclass=20)

############################################################
## Consider the general model of Y = alpha + add*G + dom*IG1 + e
## IG1 is the indicator variable I(G=1)=1, otherwise =0
## dom=0 means additive model
## Study the variance decomposition: V_G=V_A+V_D
############################################################
# parameter values
# sample size
n=1000 
# regression parameters
alpha=1
add=2
sigma=1
# MAF
p=0.2
# dom will be specified later

## G component
## obtain genotype counts assuming HWE
nG=rmultinom(1,size=n,prob=c((1-p)^2,2*p*(1-p), p^2))
print(nG)
# obtain the actual genotype vector of 0, 1 and 2
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))
IG1=c(rep(0,nG[1]),rep(1,nG[2]),rep(0,nG[3]))

## Y component then plot assuming different dom values
par(mfrow=c(2,3))

#################################
## first the additive model
dom=0
#################################
Y=alpha+add*G+dom*IG1+rnorm(n,mean=0,sd=sigma)
plot(G,Y,xlim=c(-0.5,2),ylim=c(-2,8), cex=0.7, main=paste("add=",add, ", dom=",dom, ", sigma=", sigma,sep=""))
title(line=0.5,paste("Y = alpha + add*G + dom*IG1 + e"),cex.main=0.8)
# the mean of each genotype group
segments(-0.5+0, mean(Y[1:nG[1]]), 0, mean(Y[1:nG[1]]))
segments(-0.5+1, mean(Y[(nG[1]+1):(nG[1]+nG[2])]), 1, mean(Y[(nG[1]+1):(nG[1]+nG[2])]))
segments(-0.5+2, mean(Y[(nG[1]+nG[2]+1):(nG[1]+nG[2]+nG[3])]), 2, mean(Y[(nG[1]+nG[2]+1):(nG[1]+nG[2]+nG[3])]))
#abline(a=alpha, b=add)
# the fitted line
abline(a=lm(Y~G)$coef[1],b=lm(Y~G)$coef[2],lty=2)
#hist(Y,xlim=c(-2,8),nclass=20)
#################################
## then consider other models
dom=1
# repeat the above
dom=2
# repeat the above
dom=3
# repeat the above
dom=-1
# repeat the above
dom=-2
# repeat the above


