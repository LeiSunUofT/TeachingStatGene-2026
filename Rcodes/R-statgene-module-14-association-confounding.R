####################################
## Spurious association CF example
####################################

## prepare data
CF.data=matrix(c(118,82,78,122), byrow=T, ncol=2)
dimnames(CF.data)=list(Y=c("Case","Control"), X=c("M","m"))
CF.data

## Pearson's chisq test
chisq.test(CF.data,correct=F)
chisq.test(CF.data,correct=F)$expected

## OR inference from 2x2 table
hat.logOR=log((CF.data[1,1]*CF.data[2,2])/(CF.data[1,2]*CF.data[2,1]))
se.logOR=sqrt(sum(1/CF.data))
z.logOR=hat.logOR/se.logOR
print(c(hat.logOR, se.logOR, z.logOR,2*pnorm(-z.logOR)))

## Logistic regression
X=c(1,0)
fit=glm(CF.data~X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

#####################################################################
## Spurious Association and Population Stratification - Binary Trait
#####################################################################

## Allele frequency of the two populations
p=0.3
q=0.5
## Proportion of population 1 in cases and controls
pi1=0.4
pi2=0.5
## Sample size of cases and controls, do not have to be equal
n1=400
n2=600

## Generate the data
# first n1 case=1 then n2 control=0
Y.phenotype=c(rep(1, n1), rep(0, n2))
# population index for each individual
# among n1 cases, n1*pi1 from population 1,  n1*(1-pi1) for population 2
# among n2 controls, n2*pi2 from population 1,  n2*(1-pi2) for population 2
Z.population=c(rep(1,n1*pi1),rep(2,n1*(1-pi1)),rep(1,n2*pi2),rep(2,n2*(1-pi2)))  
# Define genotype vector
X.genotype=rep(-1, length(Y.phenotype))
# identify population 1 samples
pop.index=which(Z.population==1)
# simulation genotype data for population 1 sample use allele freq p 
nG=rmultinom(1,size=length(pop.index),prob=c((1-p)^2,2*p*(1-p), p^2))
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))
# shuffle the genotype so that it has nothing to do with the phenotype Y
G=sample(G)
X.genotype[pop.index]=G
# do it similarly for population 2
pop.index=which(Z.population==2)
nG=rmultinom(1,size=length(pop.index),prob=c((1-q)^2,2*q*(1-q), q^2))
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))
# shuffle the genotype so that it has nothing to do with the Y
G=sample(G)
X.genotype[pop.index]=G

# check the length of the vectors
length(Y.phenotype)
length(Z.population)
length(X.genotype)

summary(glm(Y.phenotype~X.genotype,family=binomial(logit)))$coefficient
summary(glm(Y.phenotype~Z.population,family=binomial(logit)))$coefficient
summary(glm(Y.phenotype~X.genotype+Z.population,family=binomial(logit)))$coefficient


#####################################################################
## Study the power of detecting/seeing the spurious association
#####################################################################
# type 1 error rate
alpha=0.05
# number of replicates for the simulation
N.rep=100
# p.value of the (incorrect) association test that could leads to spurious association
p.value=rep(-1,N.rep)

# the paramater values
p=0.3
q=0.5
pi1=0.4
pi2=0.5
n1=4000
n2=6000

for (i in 1:N.rep) {
# simulate the data	
Y.phenotype=c(rep(1, n1), rep(0, n2))
Z.population=c(rep(1,n1*pi1),rep(2,n1*(1-pi1)),rep(1,n2*pi2),rep(2,n2*(1-pi2)))  
X.genotype=rep(-1, length(Y.phenotype))
pop.index=which(Z.population==1)
nG=rmultinom(1,size=length(pop.index),prob=c((1-p)^2,2*p*(1-p), p^2))
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))
G=sample(G)
X.genotype[pop.index]=G
pop.index=which(Z.population==2)
nG=rmultinom(1,size=length(pop.index),prob=c((1-q)^2,2*q*(1-q), q^2))
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))
G=sample(G)
X.genotype[pop.index]=G
# grab the p-value from the (incorrect) association analysis
p.value[i]=summary(glm(Y.phenotype~X.genotype,family=binomial(logit)))$coefficient[2,4]
}
# plot the distribution of the p-values
hist(p.value)
# estimate the power
sum(p.value<alpha)/N.rep

#####################################################################
## Show that there is no spurious association if we include the population index
#####################################################################
# type 1 error rate
alpha=0.05
# number of replicates for the simulation
N.rep=1000
# p.value of the (incorrect) association test that could leads to spurious association
p.value=rep(-1,N.rep)

# the paramater values
p=0.3
q=0.7
pi1=0.4
pi2=0.6
n1=400
n2=600

for (i in 1:N.rep) {
# simulate the data	
Y.phenotype=c(rep(1, n1), rep(0, n2))
Z.population=c(rep(1,n1*pi1),rep(2,n1*(1-pi1)),rep(1,n2*pi2),rep(2,n2*(1-pi2)))  
X.genotype=rep(-1, length(Y.phenotype))
pop.index=which(Z.population==1)
nG=rmultinom(1,size=length(pop.index),prob=c((1-p)^2,2*p*(1-p), p^2))
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))
G=sample(G)
X.genotype[pop.index]=G
pop.index=which(Z.population==2)
nG=rmultinom(1,size=length(pop.index),prob=c((1-q)^2,2*q*(1-q), q^2))
G=c(rep(0,nG[1]),rep(1,nG[2]),rep(2,nG[3]))
G=sample(G)
X.genotype[pop.index]=G
# grab the p-value from the (incorrect) association analysis
p.value[i]=summary(glm(Y.phenotype~X.genotype+Z.population,family=binomial(logit)))$coefficient[2,4]
}
# plot the distribution of the p-values
hist(p.value)
# estimate the power
sum(p.value<alpha)/N.rep


#####################################################################
## Spurious Association and Population Stratification - QTL
#####################################################################

# allele frequency in the two populations
p=0.3
q=0.3
# population mean of the two populations
mu1=1
mu2=1
# sample size for each poulations
n1=200
n2=400

# data for indiviudals with genotype X=2
y1=rnorm(n1*p^2,mean=mu1)
y2=rnorm(n2*q^2,mean=mu2)
y=c(y1,y2)
all.y=y
all.g=rep(2,(length(y1)+length(y2)))
#population index
all.z=c(rep(1,length(y1)),rep(2,length(y2)))

# data for indiviudals with genotype X=1
y1=rnorm(n1*2*p*(1-p),mean=mu1)
y2=rnorm(n2*2*q*(1-q),mean=mu2)
y=c(y1,y2)
all.y=c(all.y,y)
all.g=c(all.g,rep(1,(length(y1)+length(y2))))
all.z=c(all.z,rep(1,length(y1)),rep(2,length(y2)))

# data for indiviudals with genotype X=0
y1=rnorm(n1*(1-p)^2,mean=mu1)
y2=rnorm(n2*(1-q)^2,mean=mu2)
y=c(y1,y2)
all.y=c(all.y,y)
all.g=c(all.g,rep(0,(length(y1)+length(y2))))
all.z=c(all.z,rep(1,length(y1)),rep(2,length(y2)))

# visual plots
which.pop1=which(all.z==1)
plot(all.g[which.pop1],all.y[which.pop1],ylim=c(min(all.y), max(all.y)))
abline(h=mu1)
points(all.g[-which.pop1],all.y[-which.pop1],col="red")
abline(h=mu2, col="red")

# Association analysis
# ignore the population index
summary(lm(all.y~all.g))$coefficient
fit.result=summary(lm(all.y~all.g))$coefficient
abline(fit.result[1,1], fit.result[2,1],col="blue")

# include the populatin index
summary(lm(all.y~all.g+all.z))$coefficient

