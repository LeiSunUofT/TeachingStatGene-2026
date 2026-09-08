####################################################################################################
### Use simulation to understand the p-value distribution, type 1 error, power, factors influencing power 
### and the multiple hypothesis testing issue via case-control association study 
### Note that codes below are written with the purpose of clarity but not necessarily "smart"; 
### should be modified if you want to scale it up.  
####################################################################################################

##########################################
### Case-Control association study
##########################################

### sample size and phenotype 
###############################

## use the "typical" sample size
no.cases<-1000
no.controls<-1000

## y being the phenotype: the case=1 or control=0 status of the samples
y<-c(rep(0, no.controls), rep(1, no.cases))

###########################################################################################
### To Understand the classical Type 1 Error
### We start with a single **null** SNP that is NOT associated with the phenotype of interest
###########################################################################################

## the MAF of this SNP in controls
maf.controls<-0.2

## the MAF of this SNP in cases: same as controls since it's not associated
maf.cases<-maf.controls ## OR=1

  ## first get the number of samples having one of the three genotypes, separately for controls and cases, 
  ## assuming HWE, distribution is multinomial 
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(n=1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))

  ## then we create the corresponding predictor vector x using the additive coding of the genotype
  ## 0, 1, 2 copies of the minor allele
  ## NOTE: make sure the coding is consistent with the geno.counts variables
 
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))

  ## take a look at the 2 by 3 table
  xtabs(~y+x)

  ## now we do a logistic regression of y on x 
  logitfit<-glm(y~x, family=binomial(link="logit"))
  
  ## summary of the results
  summary.logitfit<-summary(logitfit)
  print(summary.logitfit)

  ## grab the p-value of testing the coefficient of x is 0
  p.value<-summary.logitfit$coefficient[2,4]
  print(p.value)

## If we were to use 0.05 as the significance threshold, i.e. alpha=0.05, 
## was this SNP significant? 

###########################################################################################
## How often do you expect to declare signficance (i.e. make a mistake) at alpha=0.05?
## The above is just 1 association study of this SNP, so to answer this question, 
## we can imagine 200 INDEPENDENT association studies have been conducted 
## (say 200 different investigators conducted INDEPENDENT association studies of this SNP of interest), 
## so on average/in the long run, how many investigators/studies will make a false claim?
###########################################################################################

no.studies<-200

## a vector to store the p.values and the test statististics from all the studies
p.value<-rep(-1, no.studies)
z.value<-rep(-1, no.studies)

for( index.studies in 1:no.studies) {
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))
  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value[index.studies]<-summary.logitfit$coefficient[2,4]
  z.value[index.studies]<-summary.logitfit$coefficient[2,3]
  print(index.studies)
}

## How often do you expect to declare signficance (i.e. make a mistake, type 1 error) at alpha=0.05?
## should be 5%
## what's the observed empirical type 1 error?

alpha<-0.05
empirical.alpha<-sum(p.value<=alpha)/no.studies
print(empirical.alpha)

## empirical type 1 error should be pretty close to the nomimal the nominal level 0.05.
## to increase the accuracy of the estimate, you can increase no.studies to say 10,000.
## Q: how many replicates are large enough?

## take a look at the histogram (empirical distribution) of the p.value and the test statisitcs
par(mfrow=c(3,2))

hist(p.value, main=" ")
par(col.main="red")
title(main=paste("type 1 error =", alpha, "; empirical = ", empirical.alpha))
abline(v=alpha, col="red")
## if we to use alpha=0.1, then
alpha<-0.1
empirical.alpha<-sum(p.value<=alpha)/no.studies
print(empirical.alpha)
par(col.main="blue")
title(main=paste("type 1 error =", alpha, "; empirical = ", empirical.alpha), col="blue", line=0.5)
abline(v=alpha, col="blue")

alpha<-0.05
hist(z.value,xlim=c(-3,7), main=" ")
abline(v=qnorm(alpha/2), col="red")
abline(v= -qnorm(alpha/2), col="red")
alpha<-0.1
abline(v=qnorm(alpha/2), col="blue")
abline(v= -qnorm(alpha/2), col="blue")

###########################################################################################
## you can derive analytically that the p.value ~ unif(0,1) if all the hypotheses are null
###########################################################################################

###########################################################################################
### To understand the classical Power
### Now with a single **alternative** SNP that is truly associated with the phenotype of interest
###########################################################################################

## phenotype
no.cases<-1000
no.controls<-1000
y<-c(rep(0, no.controls), rep(1, no.cases))

## the MAF of this SNP in controls
maf.controls<-0.2

## the MAF of this SNP in cases: a bit higher. 
## The difference will determine the power (certainly also the sample size)
maf.cases<-0.22 # OR~1.13

## then we can use the codes above and obtain the empirical power, i.e. power estimate
## you can also calculate the expected power analytically!

no.studies<-200
## a vector to store the p.values from all the studies
p.value<-rep(-1, no.studies)
z.value<-rep(-1, no.studies)

for( index.studies in 1:no.studies) {
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))
  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value[index.studies]<-summary.logitfit$coefficient[2,4]
  z.value[index.studies]<-summary.logitfit$coefficient[2,3]
  print(index.studies)
}

############################################################################################
## notice the p.value distribution now looks different from Unif(0,1)
## Q: what kind of distributuion do you expect for p.values derived under the alternative?
#############################################################################################

alpha<-0.05
empirical.power<-sum(p.value<=alpha)/no.studies
print(empirical.power)

hist(p.value, main=" ")
par(col.main="red")
title(main=paste("type 1 error =", alpha, "; power= ", empirical.power))
abline(v=alpha, col="red")

hist(z.value,xlim=c(-3,7), main="")
abline(v=qnorm(alpha/2), col="red")
abline(v= -qnorm(alpha/2), col="red")

#######################################################
### Q: what are factors that influence power?
#######################################################

#####################################################################################################
## We can artifically increase the power by allowing a larger type 1 error 
## the trade-off between type 1 error and power
## say if alpha=0.1, but the rest parameters stay the same, in fact we assume the data stay the same
####################################################################################################

alpha<-0.1

empirical.power<-sum(p.value<=alpha)/no.studies
print(empirical.power)

hist(p.value, main=" ")
par(col.main="blue")
title(main=paste("type 1 error =", alpha, "; power= ", empirical.power))
abline(v=alpha, col="blue")


hist(z.value,xlim=c(-3,7), main="")
abline(v=qnorm(alpha/2), col="blue")
abline(v= -qnorm(alpha/2), col="blue")

#########################################################################################
## However, allowing alpha=0.1 means that if the SNP is actually NOT associated, 
## now 10% instead of 5% of the studies/investigators will make mistake/false positive.
########################################################################################

#######################################################################################
## We can increase the power by increase the sample size
## say no.cases=no.controls=2000, and rest stay the same
########################################################################################

####################################################
## For comparison, first redo the 1000 sample plots
#####################################################

no.cases<-1000
no.controls<-1000
y<-c(rep(0, no.controls), rep(1, no.cases))
maf.controls<-0.2
maf.cases<-0.22
no.studies<-200
p.value<-rep(-1, no.studies)
z.value<-rep(-1, no.studies)

for( index.studies in 1:no.studies) {
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))
  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value[index.studies]<-summary.logitfit$coefficient[2,4]
  z.value[index.studies]<-summary.logitfit$coefficient[2,3]
  print(index.studies)
}

alpha<-0.05
empirical.power<-sum(p.value<=alpha)/no.studies
print(empirical.power)

par(mfrow=c(3,2))
hist(p.value, main=" ")
par(col.main="red")
title(main=paste("type 1 error =", alpha, "; power= ", empirical.power))
abline(v=alpha, col="red")

hist(z.value,xlim=c(-3,7), main=paste("sample size =", no.cases, "+", no.controls, "; maf.diff =", (maf.cases-maf.controls)))
abline(v=qnorm(alpha/2), col="red")
abline(v= -qnorm(alpha/2), col="red")

#############################################
## Now increase the sample size to 2000 each
#############################################
no.cases<-2000
no.controls<-2000
y<-c(rep(0, no.controls), rep(1, no.cases))
maf.controls<-0.2
maf.cases<-0.22
no.studies<-200
p.value<-rep(-1, no.studies)
z.value<-rep(-1, no.studies)

for( index.studies in 1:no.studies) {
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))
  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value[index.studies]<-summary.logitfit$coefficient[2,4]
  z.value[index.studies]<-summary.logitfit$coefficient[2,3]
  print(index.studies)
}


alpha<-0.05
empirical.power<-sum(p.value<=alpha)/no.studies
print(empirical.power)

hist(p.value, main=" ")
par(col.main="red")
title(main=paste("type 1 error =", alpha, "; power= ", empirical.power))
abline(v=alpha, col="red")

hist(z.value,xlim=c(-3,7), main=paste("sample size =", no.cases, "+", no.controls, "; maf.diff =", (maf.cases-maf.controls)))
abline(v=qnorm(alpha/2), col="red")
abline(v= -qnorm(alpha/2), col="red")

############################################################################################
## We can "increase" the distance between the null and alternative
## say maf.cases=0.24 instead of 0.22 (maf.control=0.2)
## Note that we cannot really increase the effect of a SNP, 
## rather, we are essentially looking at a different SNP that has a bigger genetic effect;
## or the actual causal variant rather the tagging SNP.
########################################################################################
maf.controls<-0.2
maf.cases<-0.24  # OR~1.26

no.cases<-1000
no.controls<-1000
y<-c(rep(0, no.controls), rep(1, no.cases))
no.studies<-200
p.value<-rep(-1, no.studies)
z.value<-rep(-1, no.studies)

for( index.studies in 1:no.studies) {
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))
  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value[index.studies]<-summary.logitfit$coefficient[2,4]
  z.value[index.studies]<-summary.logitfit$coefficient[2,3]
  print(index.studies)
}


alpha<-0.05
empirical.power<-sum(p.value<=alpha)/no.studies
print(empirical.power)

hist(p.value, main=" ")
par(col.main="red")
title(main=paste("type 1 error =", alpha, "; power= ", empirical.power))
abline(v=alpha, col="red")

hist(z.value,xlim=c(-3,7), main=paste("sample size =", no.cases, "+", no.controls, "; maf.diff =", (maf.cases-maf.controls)))
abline(v=qnorm(alpha/2), col="red")
abline(v= -qnorm(alpha/2), col="red")

#########################################################################################
### The above is emprical evaluation via simulations. How would you do it analtyically?
### Need to revisit some of the classical sample size and power calculation.
### Main messages for a single hypothesis test:
### 1. always start with type 1 error control
### 2. given the same type 1 error, power goes up if you increase the sample size, the 
###########################################################################################

################################################################################################################
### Something to think about: potential pitfalls of p-value
##############################################################################################################
### If maf.controls<-0.2 and maf.cases<-0.202, 
### do we ever get a significant result? If so, do you really care about this SNP?
### e.g. sample size from 1K, 10K to 100K 
### many current large meta-analyses have sample size 100K or more; what's the implication? 
### Check Yang et al. (2011). 
## Genomic inflation factors under polygenic inheritance. European Journal of Human Genetics (2011) 19, 807–812


maf.controls<-0.2
maf.cases<-0.202  # OR~1.01

for(i in 1:3) {
  no.cases<-1000*(10^i)
  no.controls<-1000*(10^i)
  y<-c(rep(0, no.controls), rep(1, no.cases))
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))
  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value<-summary.logitfit$coefficient[2,4]
 print(p.value)
}

## Notice the sampling variation and also the fact that p-value depends on the sample size
## Even maf differ by say 0.0001 (OR~1.0006), as long as you increase the sample size, p-value will be very small
################################################################################################################


###########################################################################################
### Now we move to the GWAS setting where in ONE SINGLE STUDY, 
### we are interested in many SNPs, i.e. association study of many SNPs at once. 
### Remember: whenever at least one SNP is significant, then we have a "successful" 
### (could be true positive but could be also false positive) GWAS.
###########################################################################################
### Again: we start with **1 single GWAS** of **many SNPs**, 
### and we assume all SNP are null SNPs/not associated SNP, i.e. generate data under the null model
### so, whatever we declare to be significant would be FALSE positives
###########################################################################################

no.cases<-1000
no.controls<-1000
y<-c(rep(0, no.controls), rep(1, no.cases))

########################################################################
### number of SNPs and genotypes
## for simplicity, start with 500 SNPs only. 
## We can increase number if we write our codes in a more efficient way
#########################################################################

no.null.snps<-500
no.alter.snps<-0
no.snps<-no.null.snps+no.alter.snps

## create a vector to store the p.value of the association test for each of the SNPs 
p.value<-rep(-1, no.snps)
z.value<-rep(-1, no.snps)

#################################################################################################
## MAF will be randomly drawn from unif(0.05, 0.5)
## use lower bound of 0.05 because of the concern of small sample counts for the rare homo group
## expected sample counts for that group is no.samples*MAF^2

## since we are currently working with null SNPs, 
## there is no worries about the MAF differences between cases and controls yet
## maf is a vector containing the MAF for each of SNPs
################################################################################################

maf.controls<-runif(no.snps, min=0.05, max=0.5)
maf.cases<-maf.controls

## for each of the SNPs, repeat what we did for a single SNP
## note that each SNP has its own maf indexed by index.snp
for(index.snp in 1:no.snps) {
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls[index.snp])^2, 2*maf.controls[index.snp]*(1-maf.controls[index.snp]), maf.controls[index.snp]^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases[index.snp])^2, 2*maf.cases[index.snp]*(1-maf.cases[index.snp]), maf.cases[index.snp]^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))

  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value[index.snp]<-summary.logitfit$coefficient[2,4]
  z.value[index.snp]<-summary.logitfit$coefficient[2,3]
  
  print(index.snp)
}

#####################################################################################
## take a look at the histogram (empirical distribution) of the p.values and z-values
####################################################################################
par(mfrow=c(2,2))
hist(p.value)
hist(z.value)

###############################################################################
## if we were to use 0.05 as the significance threshold, i.e. alpha=0.05, 
## how many SNPs would we EXPECT to be significant and 
## how many SNPs would we actually declare "associated" with the phenotype?
################################################################################
alpha<-0.05
no.null.snps*alpha
sum(p.value<=alpha)

########################################################################################
## so if we were to use alpha=0.05, we will always have something significant to declare!
## If we were to repeat the above for say no.studies=200 INDEPENDENT studies, 
## then every single study will declare significance, i.e. make a mistake, 
## and the actually/empirical type 1 error will be 100% instead of the nomial 5%!!!
############################################################################################

############################################################################################
### Some additional points to make
## can also use QQ-plot to compare it with the unif distribution
## better to work on the -log10 (p.value) scale: interested in the very small p-value
## expected ordered p.values of no.snps is 
## 1, no.snps-1/no.snps, no.snps-2/no.snps,...,2/no.snps, 1/no.snps based on unif(0,1) distribution
#####################################################################################################
x.name<-"Expected"
y.name<-"Observed -log10(p-value)"
main.name<-"QQ-plot of -log10(p-value)"
qqplot(-log10(seq(no.snps:1)/no.snps), sort(-log10(p.value)), xlim=c(0,6), ylim=c(0,6), xlab=x.name, ylab=y.name, main=main.name)
abline(0,1)

## on the normal Z scale
x.name<-"Expected"
y.name<-"Observed Test Statistic, Z"
main.name<-"QQ-plot of Z"
qqnorm(z.value, xlab=x.name, ylab=y.name, main=main.name)
abline(0,1)

#############################################################################################################
## outliers are SNPs of interest: p-value smaller than what's expected, evidence for association/significance
## Q: confidence bands? how to determine which ones are outliers? 
## what far the points have to be away from the line? 
## This is back to the proper control of type 1 error rate control.
##############################################################################################################

#######################################################################################
## This leads to discussions on how to better control the error 
## and consider different kinds of measures of errors as well
##########################################################################################

##################################
## Bonferroni correction
####################################
## What if we used alpha.gw = alpha/no.snpss? 
## Would this study be significant (i.e any SNP with p < alpha.gw?)
#####################################################################
alpha.gw <- alpha/no.snps
print(alpha.gw)
sum(p.value<=alpha.gw)
print(min(p.value))

################################################################################
## we now imagine 100 INDEPENDENT GWAS have been conducted 
## (say 100 different investigators conducted independent GWAS, each used a 500 chip for genotyping
## so on average/in the long run, 
## how many investigators/studies will call their GWAS significant and make a false claim?
##############################################################################################################


########################################################
## for the sake of time, we start with 10 investigators
########################################################
no.studies<-10
no.null.snps<-500
no.alter.snps<-0
no.snps<-no.null.snps+no.alter.snps

alpha<-0.05
alpha.gw <- alpha/no.snps

## count the number of significant GWAS studies at a given alpha level
count.sig.studies.alpha <- 0
count.sig.studies.alpha.gw <- 0

no.cases<-1000
no.controls<-1000
y<-c(rep(0, no.controls), rep(1, no.cases))

maf.controls<-runif(no.snps, min=0.05, max=0.5)
maf.cases<-maf.controls
p.value<-rep(-1, no.snps)

for(index.study in 1:no.studies) {

   for(index.snp in 1:no.snps) {
     geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls[index.snp])^2, 2*maf.controls[index.snp]*(1-maf.controls[index.snp]), maf.controls[index.snp]^2))
     geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases[index.snp])^2, 2*maf.cases[index.snp]*(1-maf.cases[index.snp]), maf.cases[index.snp]^2))
     x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))

    logitfit<-glm(y~x, family=binomial(link="logit"))
    summary.logitfit<-summary(logitfit)
    p.value[index.snp]<-summary.logitfit$coefficient[2,4]
   }

    print(index.study)

   if (sum(p.value<=alpha) > 0)
    count.sig.studies.alpha <- count.sig.studies.alpha + 1

    if (sum(p.value<=alpha.gw) > 0)
    count.sig.studies.alpha.gw <- count.sig.studies.alpha.gw + 1

}

print(count.sig.studies.alpha/no.studies)
print(count.sig.studies.alpha.gw/no.studies)

##########################################################################################################
## Is Bonferroni correction a good approach in this context?

## Alternative control: False Discovery Rate (FDR)
## Check Craiu RV, Sun L (2008). Choosing the lesser evil: trade-off between false discovery rate and non-discovery rate. Statistica Sinica 18:861-879. for some background reading  on this topic.
##  Multiple hypothesis testing is a cruical issue in the genome-wide studies, incuding the incoming sequencing analysis: more number of tests/snps to worry about. In somesense, this type of agnostic approach to the genome is not efficient: Sun L (2011). On the efficiency of genome-wide scans: a multiple hypothesis testing perspective. U.P.B. Sci. Bull., Series A., 73(1):19-26.

## Incorporating prior information to increase power: 

##[Theory] Sun et al. (2006). Stratified false discovery control for large-scale hypothesis testing with application to genome-wide association studies. Genetic Epidemiology 30:519-530. 
##[Theory] Yoo et al. (2010). Were genome-wide linkage studies a waste of time? Exploiting candidate regions within genome-wide association studies. Genetic Epidemiology 34:107-118. 

## [Application 1] Wright et al. (2011). Genome-wide association and linkage identify modifier loci of lung disease severity in cystic fibrosis at 11p13 and 20q13.2. Nature Genetics 43:539-548.
## [Application 2] Sun et al. (2012). Multiple apical plasma membrane constituents are associated with susceptibility to meconium ileus in individuals with cystic fibirosis. Nature Genetics 44:562-569.

## [Application 3] Andreassen et al. (2013). Improved detection of common variants associated with schizophrenia by leveraging pleiotropy with cardiovascular-disease risk factors. American journal of Human Genetics 92(2):197-209.
## [Application 4] Schork et al. (2013). All SNPs Are Not Created Equal: Genome-Wide Association Studies Reveal a Consistent Pattern of Enrichment among Functionally Annotated SNPs. PLoS Genet 9(4): e1003449. doi:10.1371/journal.pgen.1003449

##########################################################################################################

############################################################################################
## q-values and FDR
## q-values are one-to-one mapped with p-values (i.e. order of the SNPs do not change)
## But, q-values are FDR-corrected p-values 
## (in contrast to the p-values corrected for family-wise error rate, FWER) 
## Rejecting SNPs with q-values <-0.05 means, among all the significant SNPs (positives),
## the expected false positive rate is 0.05, E[FP/P]=0.05
## q-values calculation follows a recursive formula starting with the least significant one
###########################################################################################

## use the last replicate from the previous study, a GWAS of a set of 500 not associated SNPs
p.val<-p.value

  ## p.val is a vection of p-values

  ## remove missing p-values 
  total.tests<-length(p.val)
  print(c("Total number of tests:",total.tests))
  p.val<-na.omit(p.val)
  total.tests.with.p <- length(p.val)
  print(c("      after removing missing data:",total.tests.with.p))
 
  ## check if the p-values are bounded
  epsilon <- 0.0001
  which.k<-which(p.val< -epsilon | p.val > (1+epsilon))
  if (length(p.val[which.k])!=0) {
    print("Error in input p-values: < 0 or > 1")
    print("Note: missing data codes should be NA") 
    quit("yes")
  }

  ## estimation of pi.0 = 1-pi.1, 
  ## the proportion of the noise, i.e. NOT associated SNPs

  ## use the method of Storey and Tibshirani (2003) 
  v.lambda<-seq(0,.95,.01)
  v.pi.0 <-rep(-1,length(v.lambda))
  for(i in 1:(length(v.lambda)))
    v.pi.0[i]<-sum(p.val>v.lambda[i])/(length(p.val)*(1-v.lambda[i]))
  fit2<-smooth.spline(v.lambda,v.pi.0,df=3,w=(1-v.lambda))
  pi.0 <- predict(fit2,x=1)$y
  ## make sure that pi.0 not > 1
  pi.0 <- min(pi.0, 1)
  print(c("pi.0=", pi.0))

  ## Calculate q-values 
  ord.p.val<-sort(p.val)
  q.val<-rep(0,length(p.val))
  q.val[length(p.val)]<-ord.p.val[length(p.val)]*pi.0
  for (i in (length(p.val)-1):1)
  q.val[i]<-min(ord.p.val[i]*length(p.val)/i*pi.0, q.val[i+1])
  print(c("min(q-value)=", q.val[1]))

####################################################################
## In comparison: the minimum Bonferroni corrected p-value would be 
#####################################################################
# not bounded version
min(p.val)*length(p.val)
# bounded version
min(min(p.val)*length(p.val), 1)


#####################################################################################################
## Additional points
## How do you define power in multiple hypothesis testing?
## What influence power?
######################################################################################################

###############################################################
## Comparison 1: 5 out 500 are associated vs. 10 vs. 100 vs. 260
## each associated with the same small genetic effect (OR~1.13)
#################################################################

par(mfrow=c(1,1))

no.snps<-500
no.alter.snps<-5 # or 5, 10 or 100 or 260
no.null.snps<-no.snps-no.alter.snps
p.value<-rep(-1,no.snps)

## For the not associated SNPs, directly draw the p-values from Unif(0,1)
p.value[1:no.null.snps]<-runif(no.null.snps)
## For the remaining associated SNPs,
maf.controls<-0.2
maf.cases<-0.22   # check what would happen if the effect size were a bit bigger from 0.22 to 0.23 to 0.24
no.cases<-1000
no.controls<-1000
for(i in 1:no.alter.snps) {
  y<-c(rep(0, no.controls), rep(1, no.cases))
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))
  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value[no.null.snps+i]<-summary.logitfit$coefficient[2,4]
  }
  
x.name<-"Expected"
y.name<-"Observed -log10(p-value)"
main.name<-"QQ-plot of -log10(p-value)"
qqplot(-log10(seq(no.snps:1)/no.snps), sort(-log10(p.value)), xlim=c(0,6), ylim=c(0,6), xlab=x.name, ylab=y.name, main=main.name)
abline(0,1)

## FDR control
p.val<-p.value
  ## use the method of Storey and Tibshirani (2003) 
  v.lambda<-seq(0,.95,.01)
  v.pi.0 <-rep(-1,length(v.lambda))
  for(i in 1:(length(v.lambda)))
    v.pi.0[i]<-sum(p.val>v.lambda[i])/(length(p.val)*(1-v.lambda[i]))
  fit2<-smooth.spline(v.lambda,v.pi.0,df=3,w=(1-v.lambda))
  pi.0 <- predict(fit2,x=1)$y
  ## make sure that pi.0 not > 1
  pi.0 <- min(pi.0, 1)
  print(c("pi.0=", pi.0))

  ## Calculate q-values 
  ord.p.val<-sort(p.val)
  q.val<-rep(0,length(p.val))
  q.val[length(p.val)]<-ord.p.val[length(p.val)]*pi.0
  for (i in (length(p.val)-1):1)
  q.val[i]<-min(ord.p.val[i]*length(p.val)/i*pi.0, q.val[i+1])
  print(c("min(q-value)=", q.val[1]))

## if we control FDR at 0.05,how many signficant SNPs?
sum(q.val<=0.05)
## compared with FWER control with Bonferroni correction
sum(p.val<=0.05/no.snps)

##########################################################################################
## Note that FDR is less stringent than control FWER, but it's not a magic bullet
## When there is only a small proprotion of the tests/SNPs associated (e.g. 5 or 10 out 500);
## FDR and FWER give similar results; But you see the difference when the proportion gets larger
##########################################################################################
## Also note the sampling variation: if you run the program several times, you may get lucky or unlucy
## By the time no.alter.snps reaches 260, think about the interpretation of the qq-plot 
## in the context of genomic control labmda, population stratificiation, do we have any here?
###############################################################################################

##########################################################################
## Comparison 2: 10 out 500 are associated vs. 100 out 5000 are associated
## All associated have the same effect
##########################################################################

no.snps<-5000          # 500 or 5000
no.alter.snps<-100     # 10  or 100
no.null.snps<-no.snps-no.alter.snps
p.value<-rep(-1,no.snps)

## For the not associated SNPs, directly draw the p-values from Unif(0,1)
p.value[1:no.null.snps]<-runif(no.null.snps)
## For the remaining associated SNPs,
maf.controls<-0.2
maf.cases<-0.22
no.cases<-1000
no.controls<-1000
for(i in 1:no.alter.snps) {
  y<-c(rep(0, no.controls), rep(1, no.cases))
  geno.counts.controls<-rmultinom(1, size=no.controls, prob=c((1-maf.controls)^2, 2*maf.controls*(1-maf.controls), maf.controls^2))
  geno.counts.cases<-rmultinom(1, size=no.cases, prob=c((1-maf.cases)^2, 2*maf.cases*(1-maf.cases), maf.cases^2))
  x<-c(rep(0, geno.counts.controls[1]), rep(1, geno.counts.controls[2]), rep(2, geno.counts.controls[3]), rep(0, geno.counts.cases[1]), rep(1, geno.counts.cases[2]), rep(2, geno.counts.cases[3]))
  logitfit<-glm(y~x, family=binomial(link="logit"))
  summary.logitfit<-summary(logitfit)
  p.value[no.null.snps+i]<-summary.logitfit$coefficient[2,4]
  }
  
x.name<-"Expected"
y.name<-"Observed -log10(p-value)"
main.name<-"QQ-plot of -log10(p-value)"
qqplot(-log10(seq(no.snps:1)/no.snps), sort(-log10(p.value)), xlim=c(0,6), ylim=c(0,6), xlab=x.name, ylab=y.name, main=main.name)
abline(0,1)

## FDR control
p.val<-p.value
  ## use the method of Storey and Tibshirani (2003) 
  v.lambda<-seq(0,.95,.01)
  v.pi.0 <-rep(-1,length(v.lambda))
  for(i in 1:(length(v.lambda)))
    v.pi.0[i]<-sum(p.val>v.lambda[i])/(length(p.val)*(1-v.lambda[i]))
  fit2<-smooth.spline(v.lambda,v.pi.0,df=3,w=(1-v.lambda))
  pi.0 <- predict(fit2,x=1)$y
  ## make sure that pi.0 not > 1
  pi.0 <- min(pi.0, 1)
  print(c("pi.0=", pi.0))

  ## Calculate q-values 
  ord.p.val<-sort(p.val)
  q.val<-rep(0,length(p.val))
  q.val[length(p.val)]<-ord.p.val[length(p.val)]*pi.0
  for (i in (length(p.val)-1):1)
  q.val[i]<-min(ord.p.val[i]*length(p.val)/i*pi.0, q.val[i+1])
  print(c("min(q-value)=", q.val[1]))

## if we control FDR at 0.05,how many signficant SNPs?
sum(q.val<=0.05)
## compared with FWER control with Bonferroni correction
sum(p.val<=0.05/no.snps)

######################################################################
## The END
### Lei Sun, December 7, 2015 
######################################################################
