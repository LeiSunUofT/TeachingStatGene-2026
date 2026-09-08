########################################################
## Display a normally distributed trait Y vs. genotype G
## Assume a normal model of Y = alpha + beta*G + error
#########################################################
# parameter values: can be changed to other values
# sample size
n=1000  
# regression parameters
alpha=1
beta=1
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
plot(G,Y,xlim=c(-0.5,2),cex=0.7, main=paste("Y = alpha + beta G + e"))
title(line=0.5,paste("alpha=",alpha,", beta=",beta, ", sigma=", sigma, ";  n=",n,", MAF=",p,sep=""),cex.main=0.8)

## Add the normal density and show the mean location for each G category, both sample estimate the true mean.
dx=seq(alpha-3*sigma, alpha+2*beta+3*sigma, 0.01)

thisG=0
thismu=alpha+thisG*beta
dy=(1/sqrt(2*pi*sigma^2))*exp(-(dx-thismu)^2/(2*sigma^2))
lines(-dy+thisG,dx)
segments(-0.5+thisG,thismu,thisG,thismu)
segments(-0.5+thisG,mean(Y[1:nG[1]]),thisG,mean(Y[1:nG[1]]),lty=3)

thisG=1
thismu=alpha+thisG*beta
dy=(1/sqrt(2*pi*sigma^2))*exp(-(dx-thismu)^2/(2*sigma^2))
lines(-dy+thisG,dx)
segments(-0.5+thisG,thismu,thisG,thismu)
segments(-0.5+thisG,mean(Y[(nG[1]+1):(nG[1]+nG[2])]),thisG,mean(Y[(nG[1]+1):(nG[1]+nG[2])]),lty=3)

thisG=2 
thismu=alpha+thisG*beta
dy=(1/sqrt(2*pi*sigma^2))*exp(-(dx-thismu)^2/(2*sigma^2))
lines(-dy+thisG,dx)
segments(-0.5+thisG,thismu,thisG,thismu)
segments(-0.5+thisG,mean(Y[(nG[1]+nG[2]+1):n]),thisG,mean(Y[(nG[1]+nG[2]+1):n]),lty=3)

## add the true and fitted regression lines
abline(alpha,beta)
abline(lm(Y~G)$coef[1],lm(Y~G)$coef[2],lty=3)


###########################################################
###########################################################
## Consider an interaction model Y = alpha + gamma*G*E + error
###########################################################
###########################################################

# parameter values
# sample size
n=1000
# regression parameters
alpha=1
gamma=2
sigma=1
# MAF
p=0.2
# P(E=1)=pE
pE=0.3

# E component
nE=rbinom(1,size=n,prob=pE)
E=c(rep(0,(n-nE)),rep(1,nE))
# make sure that E is indepedent of G.
E=sample(E)
# the mean and variance of the E, 
# this will be related to the mena and variance of Y|G.
muE=pE
varE=pE*(1-pE)

## G component 
nG=rmultinom(1,size=n,prob=c((1-p)^2,2*p*(1-p), p^2))
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))

## Y component
Y=alpha+gamma*G*E+rnorm(n,mean=0,sd=sigma)

## Plot Y vs. G
plot(G,Y,xlim=c(-0.5,2),cex=0.7, main=paste("Y = alpha + gamma GxE +e"))
title(line=0.5,paste("alpha=",alpha,", gamma=",gamma, ", sigma=", sigma, ";  n=",n,", MAF=",p,", P(E=1)=",pE, sep=""),cex.main=0.8)

## Add the normal density and show the mean location for each G category,
#dx=seq(alpha-4, alpha+2*gamma*pE+5, 0.01)
dx=seq(alpha-3*sigma, alpha+2*gamma*pE+3*sqrt(4*gamma*varE+sigma^2), 0.01)

## the mean and variance of the E and this will affect the mean and variance of Y for each genotype group.
thisG=0
thismu=alpha+thisG*gamma*muE
thisvar=(thisG*gamma)^2*varE+sigma^2
dy=(1/sqrt(2*pi*thisvar))*exp(-(dx-thismu)^2/(2*thisvar))
lines(-dy+thisG,dx)
segments(-0.5+thisG,thismu,thisG,thismu)

thisG=1
thismu=alpha+thisG*gamma*muE
thisvar=(thisG*gamma)^2*varE+sigma^2
dy=(1/sqrt(2*pi*thisvar))*exp(-(dx-thismu)^2/(2*thisvar))
lines(-dy+thisG,dx)
segments(-0.5+thisG,thismu,thisG,thismu)

thisG=2
thismu=alpha+thisG*gamma*muE
thisvar=(thisG*gamma)^2*varE+sigma^2
dy=(1/sqrt(2*pi*thisvar))*exp(-(dx-thismu)^2/(2*thisvar))
lines(-dy+thisG,dx)
segments(-0.5+thisG,thismu,thisG,thismu)
###########################################################
