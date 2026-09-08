####################################################
## Graphic display of the LOD scores
####################################################

## the range of theta from 0 to 1 on a grid of 0.01
theta=seq(0,0.5,0.01)

# family 1
r=1
n=4
LOD=log10(theta^r*(1-theta)^(n-r)/0.5^n)
plot(theta,LOD,type="n",ylim=c(-1.5,1))
lines(theta,LOD,lty=2)
abline(v=r/n,lty=2)
LOD1=LOD

# family 2
r=1
n=5
LOD=log10(theta^r*(1-theta)^(n-r)/0.5^n)
lines(theta,LOD,lty=3)
abline(v=r/n,lty=3)
LOD2=LOD

# combined 
LOD=LOD1+LOD2
lines(theta,LOD,lty=1)
r=1+1
n=4+5
abline(v=r/n,lty=1)
title(line=1,"r1=1, n1=4; r2=1, n2=5")

###################################################
## What would happen if the sample size increased?
####################################################
# family 1
r=10
n=40
LOD=log10(theta^r*(1-theta)^(n-r)/0.5^n)
plot(theta,LOD,type="n",ylim=c(-1.5,7))
lines(theta,LOD,lty=2)
abline(v=r/n,lty=2)
LOD1=LOD

# family 2
r=10
n=50
LOD=log10(theta^r*(1-theta)^(n-r)/0.5^n)
lines(theta,LOD,lty=3)
abline(v=r/n,lty=3)
LOD2=LOD

# combined
LOD=LOD1+LOD2
lines(theta,LOD,lty=1)
r=10+10
n=40+50
abline(v=r/n,lty=1)
title(line=1,"r1=10, n1=40; r2=10, n2=50")

########################################################
#########################################################


####################################################
## Graphic display of the log likelihood function, complex family 1
####################################################
theta=seq(0,0.5,0.01)
r=1
s=3
logL=log(theta^r*(1-theta)^s+theta^s*(1-theta)^r)
plot(theta,logL)
logL1=logL
####################################################
## Graphic display of the LOD function, complex family 1
####################################################
theta=seq(0,0.5,0.01)
r=1
s=3
LOD=log10(theta^r*(1-theta)^s+theta^s*(1-theta)^r)-log10(0.5^(r+s)+0.5^(r+s))
plot(theta,LOD)
LOD1=LOD

####################################################
## Graphic display of the log likelihood function, complex family 2
####################################################
theta=seq(0,0.5,0.01)
r=1
s=4
logL=log(theta^r*(1-theta)^s+theta^s*(1-theta)^r)
plot(theta,logL)
logL2=logL
####################################################
## Graphic display of the LOD function, complex family 2
####################################################
theta=seq(0,0.5,0.01)
r=1
s=4
LOD=log10(theta^r*(1-theta)^s+theta^s*(1-theta)^r)-log10(0.5^(r+s)+0.5^(r+s))
plot(theta,LOD)
LOD2=LOD

####################################################
## Graphic display of the log likelihood function, combining families 1 and 2
####################################################
LOD=LOD1+LOD2
logL=logL1+logL2

plot(theta,logL,type="n",ylim=c(-9,-2))
lines(theta,logL, lty=1)
lines(theta,logL1, lty=2)
lines(theta,logL2, lty=3)

plot(theta,LOD,type="n",ylim=c(-1,0.2))
lines(theta,LOD, lty=1)
lines(theta,LOD1, lty=2)
lines(theta,LOD2, lty=3)

