#################################################################################################
## Lei Sun, University of Toronto
#################################################################################################
## Basic understanding Likelihood and
## MLE and Hypothesis Testing
#################################################################################################
## Issues considered inlcuded
## 
#################################################################################################

#################################################################################################
## Use a Bionomial coin example
#################################################################################################
## Flip a coin and determine and determine if it's a fair coin or not
## X = number of heads/'success' out of n 
## X ~ Bino(n,theta)  # assumptions?
## E(X)=n*theta and Var(X)=n*theta*(1-theta)
#################################################################################################

####################################
## Basic Likelihood and MLE
## what to show that our intuitive esimate below is MLE
####################################
n=100 # sample size
x.obs.data=57 # observed data
x.obs.data/n # the intuitive estimate

####################################
# LIKELIHOOD FUNCTION in terms of THETA:
# L(theta)=P(x, theta)=choose(n, x)*theta^x*(1-theta)^(n-x)
####################################

# probability of the data x.obs.data=57 using the bionomial
# P(x, theta)
# the general expression is choose(n, x)*theta^x*(1-theta)^(n-x)
# x=x.obs.data=57 is FIXED for a given dataset/expriment
# SO, THETA is the VARIABLE: L(theta)=f(theta)=choose(n, x)*theta^x*(1-theta)^(n-x)
# Now, we have a LIKELIHOOD FUNCTION in terms of THETA: 
# L(theta)=P(x, theta)=choose(n, x)*theta^x*(1-theta)^(n-x)

# take a look at this function:
theta=seq(0, 1, 0.001) # theta from 0 to 1 on a grid of 0.001
x=x.obs.data
L.theta=choose(n, x)*theta^x*(1-theta)^(n-x)
plot(theta,L.theta)

# MLE: the theta value that maximize the L(theta)
# But how?
# several approaches: grid search, score function, or iternative procedures.
# our intutive estimate appears to be MLE.
abline(v=x.obs.data/n)

####################################
# Before proceed further,
# let's look at variations of the Likelihood
# (Why do we typtically use the log of the kernal of the Likelihood?)
####################################
# the original one
L.theta=choose(n, x)*theta^x*(1-theta)^(n-x)
# the kernal part
L.kernal.theta=theta^x*(1-theta)^(n-x)
# based on the log
l.theta=log(choose(n, x))+x*log(theta)+(n-x)*log(1-theta)
l.kernal.theta=x*log(theta)+(n-x)*log(1-theta)

par(mfrow=c(2,2))
plot(theta,L.theta)
abline(v=x.obs.data/n)
plot(theta,L.kernal.theta)
abline(v=x.obs.data/n)
plot(theta,l.theta)
abline(v=x.obs.data/n)
plot(theta,l.kernal.theta)
abline(v=x.obs.data/n)

# 1) a theta that maximize a function also maximize the log function, 
#     so the same MLE inference based on L or l=log(L)
# 2) the coefficient choose(n, x) has nothing to do with the parameter theta of interest, 
#     so it does not affect how we find the MLE for theta.
# These two points are also clear based on the above graph.
####################################
# Which one to use?
####################################
# Let's change the data to the following
n=10000 # sample size
x.obs.data=5700 # observed data
x=x.obs.data
# and repeat the calculation and graph 
L.theta=choose(n, x)*theta^x*(1-theta)^(n-x)
# the kernal part
L.kernal.theta=theta^x*(1-theta)^(n-x)
# based on the log
l.theta=log(choose(n, x))+x*log(theta)+(n-x)*log(1-theta)
l.kernal.theta=x*log(theta)+(n-x)*log(1-theta)

par(mfrow=c(2,2))
plot(theta,L.theta)
abline(v=x.obs.data/n)
plot(theta,L.kernal.theta)
abline(v=x.obs.data/n)
plot(theta,l.theta)
abline(v=x.obs.data/n)
plot(theta,l.kernal.theta)
abline(v=x.obs.data/n)

# At least for numerical reasons,
# (overflow or underflow computing errors in calcuating choose(n, x) and theta^x) 
# we know that log of the kernal part of the likelihood is preferred:
# l.kernal.theta=x*log(theta)+(n-x)*log(1-theta)
# Another reason is more statistical inference related: easier to work out the MLE using score function.

####################################
# Now we will be focusing on using the likelihood function using the version of 
# l.kernal.theta=x*log(theta)+(n-x)*log(1-theta)
# and try find MLE using grid search, score function, or iternative procedures.
####################################

# back to the previous sample size and data
n=100 # sample size
x.obs.data=57 # observed data
x=x.obs.data
l.kernal.theta=x*log(theta)+(n-x)*log(1-theta)
par(mfrow=c(1,1))
plot(theta,l.kernal.theta)

####################################
# MLE using grid search
# The previous plotting already gives us all the grids
# only need to find out where was the maximum
####################################
print(theta[which(l.kernal.theta==max(l.kernal.theta))])
# a numerically more stable way
eps=10^(-8)
print(theta[which(abs(l.kernal.theta-max(l.kernal.theta))<eps)])

####################################
# MLE using score function, a function interms of theta: score(theta)
####################################
# Score function: also a function of theta (related to the concept of tangent)
# derivative of the likelihood function, l.kernal.theta 
f.score.theta=x/theta-(n-x)/(1-theta)  

# draw a few tangent lines at different theta values:
# for theta=0.2, the slope is the score function value at theta=0.2
# b=slope=score(theta)
# intercept can be derived as well though just a bit more involved
# score(theta)=(l.kernal.theta(theta)-intercept)/theta
# so a=intercept=l.kernal.theta(theta)-score(theta)*theta
theta.plot=0.2
theta.where=which(abs(theta-theta.plot)<eps)
abline(v=theta.plot)
abline(a=(l.kernal.theta[theta.where]-(x/theta[theta.where]-(n-x)/(1-theta[theta.where]))*theta[theta.where]),b=f.score.theta[theta.where],lwd=3)
# change it to 0.8
theta.plot=0.8
theta.where=which(abs(theta-theta.plot)<eps)
abline(v=theta.plot)
abline(a=(l.kernal.theta[theta.where]-(x/theta[theta.where]-(n-x)/(1-theta[theta.where]))*theta[theta.where]),b=f.score.theta[theta.where],lwd=3)
# change it to 0.57, our MLE
theta.plot=0.57
theta.where=which(abs(theta-theta.plot)<eps)
abline(v=theta.plot)
abline(a=(l.kernal.theta[theta.where]-(x/theta[theta.where]-(n-x)/(1-theta[theta.where]))*theta[theta.where]),b=f.score.theta[theta.where],lwd=3)
# what was the slop when theta.plot=0.57, the MLE?
print(f.score.theta[theta.where])
print(x/0.57-(n-x)/(1-0.57))
# this is the same as 
print((x/n-0.57)/(0.57*(1-0.57)/n))
# where know 0.57=x.obs/n=57/100!

# So MLE is the solution of f.score.theta==0:
# score(theta)=f.score.theta=x/theta-(n-x)/(1-theta)=0 
# Thus theta=x/n would be the theta.hat=MLE
# (in practice, also need to check it's not a local maximum or minimum.)

####################################
# MLE using iterative procedure, 
# e.g. the Newton-Raphson Method
####################################

# first plot both the likelihod function and score function
n=100 # sample size
x.obs.data=57 # observed data
x=x.obs.data
l.kernal.theta=x*log(theta)+(n-x)*log(1-theta)
f.score.theta=x/theta-(n-x)/(1-theta)

par(mfrow=c(1,2))
plot(theta,l.kernal.theta)
plot(theta,f.score.theta)

# remember the score function f.score.theta is a function of theta,
# want to find the theta value that f(theta)=0

# zoon in a bit:
l.bound<-which(abs(theta-0.15)<eps)
u.bound<-which(abs(theta-0.85)<eps)
par(mfrow=c(1,1))
plot(theta[l.bound:u.bound],f.score.theta[l.bound:u.bound])
abline(h=0)


par(mfrow=c(1,2))
plot(theta[l.bound:u.bound],l.kernal.theta[l.bound:u.bound])
plot(theta[l.bound:u.bound],f.score.theta[l.bound:u.bound])
abline(h=0)


####################################
## Basic Hypothesis Testing
## Need a null hypothesis
## H0: theta=theta0=0.5
####################################
n=100 # sample size
x.obs.data=57 # observed data
theta.0=0.5 # the null huypothesis

theta.hat.data=x.obs.data/n # MLE (why and properties?) 

####################################
## Consider MANY different tests
####################################
no.test=12 # the number of the different test statistics to be considered
T.name=c("T1=obs","T2=obs-exp","T3=(obs-exp)/sd0","T4=(obs-exp)/sd1",  "T5=theta.hat","T6=theta.hat-theta.0","T7=(theta.hat-theta.0)/sd0","T8=(theta.hat-theta.0)/sd1",  "T9=L0", "T10=L.hat-L0","T11=Lhat/L0","T12=2log(Lhat/L0)") # (T7)^2 would be the Score test, and (T8)^2 the Wald test, T12=LRT
T.obs.data=rep(0,no.test) # the vector for the different tests

x.obs=x.obs.data # rename the varialbe so we can use the following codes repeatedly
T.obs=T.obs.data
theta.hat=theta.hat.data

T.obs[1]=x.obs # use the observed count itself
T.obs[2]=x.obs-n*theta.0 # x.obs-x.expected.null
T.obs[3]=(x.obs-n*theta.0)/sqrt(n*theta.0*(1-theta.0)) # standardize by the variance |H0
T.obs[4]=(x.obs-n*theta.0)/sqrt(n*theta.hat*(1-theta.hat)) # but can also use theta.hat

T.obs[5]=theta.hat # use the estimated theta instead, and the rest is the same argument
T.obs[6]=theta.hat-theta.0 # compare the estimate with the null specification
T.obs[7]=(theta.hat-theta.0)/sqrt(theta.0*(1-theta.0)/n) # standardize by the variance |H0
T.obs[8]=(theta.hat-theta.0)/sqrt(theta.hat*(1-theta.hat)/n) # but can also use theta.hat

T.obs[9]= choose(n,x.obs)*theta.0^x.obs*(1-theta.0)^(n-x.obs) # probability/density of the data assuming H0 is true 
T.obs[10]=choose(n,x.obs)*{theta.hat^x.obs*(1-theta.hat)^(n-x.obs)} - choose(n,x.obs)*{theta.0^x.obs*(1-theta.0)^(n-x.obs)} # why not 'compare' it with the maximum possible, i.e. using theta.hat which maximized the likelihood!
T.obs[11]={theta.hat^x.obs*(1-theta.hat)^(n-x.obs)} / {theta.0^x.obs*(1-theta.0)^(n-x.obs)} # just another way to compare it; # choose(n,x.obs) not used, the concept of the kernal of the likelihood
T.obs[12]=2*{{x.obs*log(theta.hat)+(n-x.obs)*log(1-theta.hat)}-{x.obs*log(theta.0)+(n-x.obs)*log(1-theta.0)}} # yet another way to compar it using the classic LRT statistic: 2*log(Lhat/L0)=2*(log(Lhat)-log(L0))  # why don't we calculate it using 2*log({theta.hat^x.obs*(1-theta.hat)^(n-x.obs)}/{theta.0^x.obs*(1-theta.0)^(n-x.obs)})=2*log(T.obs[11])?  

T.obs.data=T.obs #save the above result using the orignal variable name
####################################
## What is our decision? 
## How big is big and how small is small? 
####################################
print(T.obs) 
print(c(min(T.obs),max(T.obs)))

###############################################################
## Use simulation study!
## For each of the n.rep independent replicates/datasets,
## if we assume data was truly generated from H0: theta=theta.0=0.5,
## and we follow the same procedure of calculating each test statistic, 
## how does the statistic look like? 
## i.e. histogram/distribution of the different test statistics
###############################################################
n.rep=10000 # the number of replicates to be done # how big it should be?
x.obs.rep=rbinom(n=n.rep, size=n, prob=theta.0) # simulate x.obs from Bino(n, theta.0) n.rep times
theta.hat.rep=x.obs.rep/n  # obtain the theta.hat estimate for each replicate

####################################
## CHECK the data simulation
####################################
par(mfrow=c(1,2))

hist(x.obs.rep) # use the count itself
abline(v=n*theta.0,col="black") # what do we expect?
abline(v=mean(x.obs.rep),col="blue") # the mean of all the replicates
abline(v=x.obs.data,col="red") # where does the data sit?  

hist(theta.hat.rep) # now use the theta estimate and the rest is the same
abline(v=theta.0,col="black") 
abline(v=mean(theta.hat.rep),col="blue") 
abline(v=theta.hat.data,col="red") 

####################################
## Calculate our test statistics for each replicate
####################################
T.obs.rep=matrix(0,nrow=n.rep, ncol=no.test) 

for(i.rep in 1:n.rep) {   # using i.rep is less prone error than i # but looping generally is not a good idea in terms of speed  
	x.obs=x.obs.rep[i.rep]  # use the previous varialbe names to use the same codes as above
	theta.hat=theta.hat.rep[i.rep]
	T.obs=T.obs.rep[i.rep,]  # replicate is by the row
	
	T.obs[1]=x.obs 
	T.obs[2]=x.obs-n*theta.0 
	T.obs[3]=(x.obs-n*theta.0)/sqrt(n*theta.0*(1-theta.0)) 
	T.obs[4]=(x.obs-n*theta.0)/sqrt(n*theta.hat*(1-theta.hat)) 

	T.obs[5]=theta.hat 
	T.obs[6]=theta.hat-theta.0 
	T.obs[7]=(theta.hat-theta.0)/sqrt(theta.0*(1-theta.0)/n) 
	T.obs[8]=(theta.hat-theta.0)/sqrt(theta.hat*(1-theta.hat)/n) 

	T.obs[9]=choose(n,x.obs)*theta.0^x.obs*(1-theta.0)^(n-x.obs) 
	T.obs[10]=choose(n,x.obs)*{theta.hat^x.obs*(1-theta.hat)^(n-x.obs)}-choose(n,x.obs)*{theta.0^x.obs*(1-theta.0)^(n-x.obs)}
	T.obs[11]={theta.hat^x.obs*(1-theta.hat)^(n-x.obs)}/{theta.0^x.obs*(1-theta.0)^(n-x.obs)}
	T.obs[12]=2*{{x.obs*log(theta.hat)+(n-x.obs)*log(1-theta.hat)}-{x.obs*log(theta.0)+(n-x.obs)*log(1-theta.0)}} 

	T.obs.rep[i.rep,]=T.obs #track the values
}
####################################
## Plot the histogram of each test statistic using the n.rep replciates
## Where does the data sit?
## Do we have the same conclusion using different test statistic?
####################################
par(mfrow=c(3,4))  
count.more.extreme=rep(0,no.test)
for(i.test in 1:no.test) {
	hist(T.obs.rep[,i.test],main=T.name[i.test],freq=F)
	abline(v=T.obs.data[i.test],col="red")
	
	count.more.extreme[i.test]=sum(T.obs.rep[,i.test]>=T.obs.data[i.test])
}
print(count.more.extreme)  
print(count.more.extreme)/n.rep #this would be the emprical p-value

####################################
## Likelihood-based tests are fundamentally two-sided
## The above T1-T8 based on the counts and proportions were one-sided
## Once we make them two-sided as below, then the values are comparable
####################################
print(2*count.more.extreme[1:8]/n.rep)
print(count.more.extreme[9:12]/n.rep)  

####################################
## "T9=L0" appears to be quite different from the others!
## When the statistic is the probablity of the data, 
## should we use T.obs.rep[,i.test]>=T.obs.data[i.test]?
## Let's count in a different way
####################################
i.test=9
print(sum(T.obs.rep[,i.test]<=T.obs.data[i.test]))
print(sum(T.obs.rep[,i.test]<=T.obs.data[i.test])/n.rep)
####################################
## In fact, this would be the 2-sided Bionomial exact test, 
## although exact tests are not neccessary given the size of the n
####################################

