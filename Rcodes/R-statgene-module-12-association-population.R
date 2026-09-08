#####################################################################
## Exercise 1 of Chapter 7
#####################################################################

## Input the data
#################
IDDM.data=matrix(c(23, 270, 293, 1343, 3284, 4627, 1366, 3554, 4920), byrow=T, ncol=3)
IDDM.data
dimnames(IDDM.data)=list(GM=c("D","d", "Total"), IDDM=c("Yes","No","Total"))
IDDM.data

# comparing two proportions without continuity correction
# Input should be just the counts of "success" and "failure" for the two categories
##################################################################################
IDDM.data[1:2,1:2]
prop.test(IDDM.data[1:2,1:2],correct=F) 
prop.test(IDDM.data[1:2,1:2],correct=F)$statistic
prop.test(IDDM.data[1:2,1:2],correct=F)$p.value

# Call the chisq.test function for the chisq test without continuity correction
# note that the R  comparing two proportions is the same as the Pearson chisq test
##################################################################################
chisq.test(IDDM.data[1:2,1:2],correct=F)
chisq.test(IDDM.data[1:2,1:2],correct=F)$statistic
chisq.test(IDDM.data[1:2,1:2],correct=F)$p.value
# Note that the expected counts are calculated from the marginal totals
chisq.test(IDDM.data[1:2,1:2],correct=F)$observed
chisq.test(IDDM.data[1:2,1:2],correct=F)$expected
matrix(c(IDDM.data[3,1]*IDDM.data[1,3]/IDDM.data[3,3],
IDDM.data[3,2]*IDDM.data[1,3]/IDDM.data[3,3],
IDDM.data[3,1]*IDDM.data[2,3]/IDDM.data[3,3],
IDDM.data[3,2]*IDDM.data[2,3]/IDDM.data[3,3]), byrow=T, ncol=2)

# If we were to calculate directly from these expected counts
# then our test statistic and p-values are 
################################################################
c.exp=matrix(c(IDDM.data[3,1]*IDDM.data[1,3]/IDDM.data[3,3],
IDDM.data[3,2]*IDDM.data[1,3]/IDDM.data[3,3],
IDDM.data[3,1]*IDDM.data[2,3]/IDDM.data[3,3],
IDDM.data[3,2]*IDDM.data[2,3]/IDDM.data[3,3]), byrow=T, ncol=2)
c.obs=IDDM.data[1:2,1:2]
sum((c.obs-c.exp)^2/c.exp)
1-pchisq(sum((c.obs-c.exp)^2/c.exp), df=1)
# This is identicial to the R function result
chisq.test(IDDM.data[1:2,1:2],correct=F)$statistic
chisq.test(IDDM.data[1:2,1:2],correct=F)$p.value

# If we were to do a normal test
##################################
pi1=IDDM.data[1,1]/IDDM.data[1,3] 
pi2=IDDM.data[2,1]/IDDM.data[2,3] 
pi=IDDM.data[3,1]/IDDM.data[3,3] 
print(c(pi1,pi2, pi))
Z.obs=(pi1-pi2)/sqrt(pi*(1-pi)/IDDM.data[1,3]+pi*(1-pi)/IDDM.data[2,3])
# note that Z^2=chisq_1, so the test statistic here is very similar to the chisq_1 stat from the Pearson's tet
print(c(Z.obs, Z.obs^2))
# so is the p-value
2*pnorm(Z.obs)

# If we were to do a likelihood ratio test
##########################################
l1=IDDM.data[1,1]*log(pi1)+IDDM.data[1,2]*log(1-pi1)+IDDM.data[2,1]*log(pi2)+IDDM.data[2,2]*log(1-pi2)
l0=IDDM.data[1,1]*log(pi)+IDDM.data[1,2]*log(1-pi)+IDDM.data[2,1]*log(pi)+IDDM.data[2,2]*log(1-pi)
T.LRT=2*(l1-l0)
p.LRT=1-pchisq(T.LRT,df=1)
# results are similar to the other tests
print(c(T.LRT, p.LRT))
# results are also similar if we used obs and exp expression shown in the class
2*sum(c.obs*log(c.obs/c.exp))

## formally conducting a logistic regression
############################################
## prepare the data
#####################
IDDM.data=matrix(c(23, 270, 293, 1343, 3284, 4627, 1366, 3554, 4920), byrow=T, ncol=3)
dimnames(IDDM.data)=list(GM=c("D","d", "Total"), IDDM=c("Yes","No","Total"))
# We need to code the GM D and d, and we choose D=1 and d=0
GM.X=c(1,0)
# Only need the yes="success" and no="failure" counts
IDDM.YN=IDDM.data[1:2,1:2]
IDDM.YN

# call the glm function with logit link function
###############################################
fit=glm(IDDM.YN~GM.X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

## we can also just use the "success" and total counts
IDDM.YT=IDDM.data[1:2,c(1,3)]
IDDM.YT
fit=glm(IDDM.YT[,1]/IDDM.YT[,2]~GM.X, family=binomial(logit), weights=IDDM.YT[,2])
fit.summary=summary(fit)
fit.summary$coefficient

## What is the logOR, estimated SE, Z test (Z^2=Wald-test) and p-value of the logOR?
fit.summary$coefficient[2,]
## Note that the test results are very similar to those comparing two proprotions!

## Results are identical to "hand calculation"
log((IDDM.YN[1,1]*IDDM.YN[2,2])/(IDDM.YN[2,1]*IDDM.YN[1,2]))
sqrt(sum(1/IDDM.YN))

## Confidence Interval for OR: first CI for logOR then convert it back
logOR=fit.summary$coefficient[2,1]
OR=exp(logOR)
ASE=fit.summary$coefficient[2,2]
# if 95% CI
c.alpha=qnorm(1-0.05/2)
CI.logOR=c(logOR-c.alpha*ASE,logOR+c.alpha*ASE)
CI.OR=exp(CI.logOR)
print(c(logOR, CI.logOR))
print(c(OR, CI.OR))
 
############################################################
## What if we changed the coding for D and d? now D=0 and d=1
#############################################################
GM.X=c(0,1)
fit=glm(IDDM.YN~GM.X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

############################################################
## What if we changed the definition of the success and failure but keep D=1, d=0
#############################################################
IDDM.YN=IDDM.data[1:2,2:1]
GM.X=c(1,0)
fit=glm(IDDM.YN~GM.X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

############################################################
## What if we changed the definition of the success and failure and also D=0, d=1
#############################################################
IDDM.YN=IDDM.data[1:2,2:1]
GM.X=c(0,1)
fit=glm(IDDM.YN~GM.X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

## Note that the statistical inference are identical, except logOR might change the sign
## The genetic interpretation does change slightly: d or D increase or decrease the success rate, and how success is defined!

############################################################
## Finally, we showed theoretically  that OR inference is robust to retropsecitve case/control design
## So let's look at the actual data
#############################################################
## Now we assume that it was a case-control study, and 
## pi1=P(Affected/D) pi2=P(Affected/d) is no longer meaningful, 
## and we are more interested in pi1*=P(D|case) and pi2* P(D|control)
## In fact, we are swithching the defintion of response variable and predictors
###############################################################################
## Same input data, but byrow now is False, so we are doing a transpose of the previous matrix
IDDM.data=matrix(c(23, 270, 293, 1343, 3284, 4627, 1366, 3554, 4920), byrow=F, ncol=3)
dimnames(IDDM.data)=list(IDDM=c("Yes","No","Total"),GM=c("D","d", "Total"))
IDDM.data
# We now code IDDM Yes/affected/case=1 and No/unaffected/control=0 
IDDM.X=c(1,0)
# Only need the yes="success" and no="failure" counts, but "success" here = D and "failure" here=d
GM.YN=IDDM.data[1:2,1:2]
GM.YN
fit=glm(GM.YN~IDDM.X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient
## Indeed, results are identical!

## Note this switching of row and columns, it will not affected the Pearson Chisq test at all
## because the expected counts stil depends on the marginal totals which do not change.
## Another way to look at this: Pearson Chisq testing two variables X and Y are independent of each other
## But, it does not specifiy which one is the predcitor.
## The Normal and LRT tests interestingly also have identical results 
## This can be understood from the LRT expressoin which also just depends on the obs and exp and exp are calculated from the margins.
## But we emphasize that the esimation of the pi itself is not meaningfu! (unlike the OR way of comparing the two)

## Note that we can also use data.frame function to input the data
IDDM.data=data.frame(GM=c(1,0), IDDM=c(23, 1343), n=c(293,4627))
IDDM.data
IDDM.data$GM
IDDM.data$IDDM
IDDM.data$n
glm(IDDM/n ~ GM, binomial(logit), data=IDDM.data, weights=n)
## or add the success and failure counts then call the glm function
IDDM.data$YN = cbind(IDDM.data$IDDM,IDDM.data$n-IDDM.data$IDDM)
IDDM.data
IDDM.data$YN
fit.lg=glm(YN ~ GM, binomial(logit), data=IDDM.data)
  
#####################################################################
#####################################################################
## Exercise 3 of Chapter 7
#####################################################################
#####################################################################

#############################################
## We now input the data in a different form by actually creating the X and Y vectors
############################################
# case and control counts
r0=500
r1=350
r2=120

s0=521
s1=270
s2=130

# marginal and grand totals
r=r0+r1+r2
s=s0+s1+s2
n0=r0+s0
n1=r1+s1
n2=r2+s2
n=r+s
#check also n0+n1+n2
print(n0+n1+n2)

##############################
## Create X and Y vectors
###############################

# Y values
Y=c(rep(1,r),rep(0,s))
# X values
# as characeters, genotypic model
X=c(rep('aa', r0), rep('Aa', r1), rep('AA',r2),rep('aa', s0), rep('Aa', s1), rep('AA',s2))
# additive coding
Xa=c(rep(0, r0), rep(1, r1), rep(2,r2),rep(0, s0), rep(1, s1), rep(2,s2))
# dominant coding
Xd=c(rep(0, r0), rep(1, r1), rep(1,r2),rep(0, s0), rep(1, s1), rep(1,s2))
# recessive coding
Xr=c(rep(0, r0), rep(0, r1), rep(1,r2),rep(0, s0), rep(0, s1), rep(1,s2))
# genotypic model using dummy variable
X1=c(rep(0, r0), rep(1, r1), rep(0,r2),rep(0, s0), rep(1, s1), rep(0,s2))
X2=c(rep(0, r0), rep(0, r1), rep(1,r2),rep(0, s0), rep(0, s1), rep(1,s2))

# take a look at the data
data=cbind(Y,Xa, Xd, Xr, X1, X2)

> data
        Y Xa Xd Xr X1 X2  Genotype X
   [1,] 1  0  0  0  0  0  aa
   [2,] 1  0  0  0  0  0  aa 
   [3,] 1  0  0  0  0  0  aa
 ...................
 [500,] 1  0  0  0  0  0  aa
 [501,] 1  1  1  0  1  0  Aa
 [502,] 1  1  1  0  1  0  Aa 
 [503,] 1  1  1  0  1  0  Aa
 ...................
 [850,] 1  1  1  0  1  0  Aa
 [851,] 1  2  1  1  0  1  AA
 [852,] 1  2  1  1  0  1  AA
 [853,] 1  2  1  1  0  1  AA
 ...................
 [970,] 1  2  1  1  0  1  AA
 [971,] 0  0  0  0  0  0  aa 
 [972,] 0  0  0  0  0  0  aa
 [973,] 0  0  0  0  0  0  aa
 ...................
[1491,] 0  0  0  0  0  0  aa
[1492,] 0  1  1  0  1  0  Aa
[1493,] 0  1  1  0  1  0  Aa
[1494,] 0  1  1  0  1  0  Aa
 ...................
[1761,] 0  1  1  0  1  0  Aa
[1762,] 0  2  1  1  0  1  AA
[1763,] 0  2  1  1  0  1  AA
[1764,] 0  2  1  1  0  1  AA
 ...................
[1889,] 0  2  1  1  0  1  AA
[1890,] 0  2  1  1  0  1  AA
[1891,] 0  2  1  1  0  1  AA

#################################
## Analysis
##################################
###################
## additive model
#####################
fit=glm(Y~Xa, family=binomial(logit))
# summary of the regression
fit.summary=summary(fit)
print(fit.summary)
fit.summary$coef

# Obtain OR (beta for X is the log OR)
exp(fit.summary$coef[2,1])
# can also print the SE of the log OR
fit.summary$coef[2,2]

# Wald test
###########
# the test statistic is estimate of logOR/SE
T.W=(fit.summary$coef[2,1]/fit.summary$coef[2,2])^2
# this is the same as (Z^2=chisq_1)
sqrt(T.W)
fit.summary$coef[2,3]
# the corresponding p-value
1-pchisq(T.W,df=1)
# this is the same as
fit.summary$coef[2,4]

# LRT
######
# grab the logL of the full model
logLik(fit)
# the logL of the null model
fit.null=glm(Y~1, family=binomial(logit))
logLik(fit.null)
# The LRT statistic
T.LRT=2*(logLik(fit)[1]-logLik(fit.null)[1])
# The p-vaulue
1-pchisq(T.LRT,df=1)

# Wald and LRT give similar result

##################
## dominant model
#####################
fit=glm(Y~Xd, family=binomial(logit))

# summary of the regression
fit.summary=summary(fit)
print(fit.summary)
fit.summary$coef

# Obtain OR (beta for X is the log OR)
exp(fit.summary$coef[2,1])
# can also print the SE of the log OR
fit.summary$coef[2,2]

# Wald test
###########
# the test statistic is estimate of logOR/SE
T.W=(fit.summary$coef[2,1]/fit.summary$coef[2,2])^2
# this is the same as (Z^2=chisq_1)
sqrt(T.W)
fit.summary$coef[2,3]
# the corresponding p-value
1-pchisq(T.W,df=1)
# this is the same as
fit.summary$coef[2,4]

# LRT
######
# the logL of the full model
logLik(fit)
# the logL of the null model
fit.null=glm(Y~1, family=binomial(logit))
logLik(fit.null)
# The LRT statistic
T.LRT=2*(logLik(fit)[1]-logLik(fit.null)[1])
# The p-vaulue
1-pchisq(T.LRT,df=1)

# Calculate these quantities directly using the 2x2 table:
##########################################################
# OR is 
((r1+r2)*s0)/(r0*(s1+s2))
# this is the same as from the regression
exp(fit.summary$coef[2,1])

# Wald test
########### 
# logOR is
log(((r1+r2)*s0)/(r0*(s1+s2)))
# SE of the logOR is
sqrt(1/r0+1/(r1+r2)+1/s0+1/(s1+s2))
# The Wald test statistic is 
T.W=(log(((r1+r2)*s0)/(r0*(s1+s2)))/sqrt(1/r0+1/(r1+r2)+1/s0+1/(s1+s2)))^2
print(sqrt(T.W))
# The p-value is 
1-pchisq(T.W,df=1)
# these quantities are also identical to the results from the regression model
print(fit.summary$coef[2,])

#############################
## recessive model is similar procesure wise
## results are quite different though
#############################

fit=glm(Y~Xr, family=binomial(logit))

# summary of the regression
fit.summary=summary(fit)
print(fit.summary)
fit.summary$coef

# Obtain OR (beta for X is the log OR)
exp(fit.summary$coef[2,1])
# can also print the SE of the log OR
fit.summary$coef[2,2]

# Wald test
###########
# the test statistic is estimate of logOR/SE
T.W=(fit.summary$coef[2,1]/fit.summary$coef[2,2])^2
# this is the same as (Z^2=chisq_1)
sqrt(T.W)
fit.summary$coef[2,3]
# the corresponding p-value
1-pchisq(T.W,df=1)
# this is the same as
fit.summary$coef[2,4]

# LRT
######
# the logL of the full model
logLik(fit)
# the logL of the null model
fit.null=glm(Y~1, family=binomial(logit))
logLik(fit.null)
# The LRT statistic
T.LRT=2*(logLik(fit)[1]-logLik(fit.null)[1])
# The p-vaulue
1-pchisq(T.LRT,df=1)

# Calculate these quantities directly using the 2x2 table:
##########################################################
# OR is 
((r1+r2)*s0)/(r0*(s1+s2))
# this is the same as from the regression
exp(fit.summary$coef[2,1])

# Wald test
########### 
# logOR is
log((r2*(s0+s1))/((r0+r1)*s2))
# SE of the logOR is
sqrt(1/(r0+r1)+1/r2+1/(s0+s1)+1/s2)
# The Wald test statistic is 
T.W=(log((r2*(s0+s1))/((r0+r1)*s2))/sqrt(1/(r0+r1)+1/r2+1/(s0+s1)+1/s2))^2
print(sqrt(T.W))
# The p-value is 
1-pchisq(T.W,df=1)
# these quantities are also identical to the results from the regression model
print(fit.summary$coef[2,])

############################################
## Genotypic model using the dummy variables
############################################
fit=glm(Y~X1+X2, family=binomial(logit))

# summary of the regression
fit.summary=summary(fit)
print(fit.summary)
fit.summary$coef

# Obtain OR for both beta1 and beta2
exp(fit.summary$coef[2,1])
exp(fit.summary$coef[3,1])

# Again, these quantities can be directly calculated using the 2x2 table:
##########################################################
# for beta1
##############
# OR  
(r1*s0)/(r0*s1)

# Wald test 
#############
# logOR is
log((r1*s0)/(r0*s1))
# SE of the logOR is
sqrt(1/r0+1/r1+1/s0+1/s1)
# The Wald test statistic is
T.W=(log((r1*s0)/(r0*s1))/sqrt(1/r0+1/r1+1/s0+1/s1))^2
print(sqrt(T.W))
# p-value is
1-pchisq(T.W,df=1)
# these quantities are also identical to the results from the regression model
print(fit.summary$coef[2,])

# LRT
######
# the logL of the full model
logLik(fit)
# the logL of the null model (H0:beta1=0), note that here we allow beta2 != 0 
fit.null=glm(Y~X2, family=binomial(logit))
logLik(fit.null)
# The LRT statistic
T.LRT=2*(logLik(fit)[1]-logLik(fit.null)[1])
# The p-vaulue
1-pchisq(T.LRT,df=1)

# for beta2
############
# OR  
(r2*s0)/(r0*s2)

# Wald test 
############
# logOR is
log((r2*s0)/(r0*s2))
# SE of the logOR is
sqrt(1/r0+1/r2+1/s0+1/s2)
# The Wald test statistic is
T.W=(log((r2*s0)/(r0*s2))/sqrt(1/r0+1/r2+1/s0+1/s2))^2
print(sqrt(T.W))
# p-value is
1-pchisq(T.W,df=1)
# these quantities are also identical to the results from the regression model
print(fit.summary$coef[3,])

# LRT
######
# the logL of the full model
logLik(fit)
# the logL of the null model (H0: beta2=0), note that here we allow beta1 != 0 
fit.null=glm(Y~X1, family=binomial(logit))
logLik(fit.null)
# The LRT statistic
T.LRT=2*(logLik(fit)[1]-logLik(fit.null)[1])
# The p-vaulue
1-pchisq(T.LRT,df=1)


# LRT joint testing H0: beta1=0 and beta2=0
############################################
# the logL of the full model
fit=glm(Y~X1+X2, family=binomial(logit))
fit.summary=summary(fit)
logLik(fit)
# the logL of the null model 
fit.null=glm(Y~1, family=binomial(logit))
logLik(fit.null)
# The LRT statistic
T.LRT=2*(logLik(fit)[1]-logLik(fit.null)[1])
# The p-vaulue, note that it's now 2 d.f.
1-pchisq(T.LRT,df=2)

# How do we do Wald test for H0: beta1=0 and beta2=0?
# need the covariance matrix of the estimated coefficients beta1 and beta2
fit.summary$cov.unscaled
             (Intercept)           X1           X2
(Intercept)  0.003919386 -0.003919386 -0.003919386
X1          -0.003919386  0.010480049  0.003919386
X2          -0.003919386  0.003919386  0.019945022

hat.beta=matrix(nrow=2, ncol=1, c(fit.summary$coef[2,1],fit.summary$coef[3,1]))
hat.cov=matrix(nrow=2,ncol=2,c(fit.summary$cov.scaled[2,2],fit.summary$cov.scaled[2,3],fit.summary$cov.scaled[3,2],fit.summary$cov.scaled[3,3]))
hat.beta
hat.cov

T.W=t(hat.beta)%*%solve(hat.cov)%*%(hat.beta)
T.W
1-pchisq(T.W,df=2)
# Again similar results as the LRT test.

# Note that estimate of beta_1 and beta_2 are correlated with each other, so we cannot use T.W=hat beta/Var(hat beta) for testing beta 1 + T.W for testing beta 2
(fit.summary$coef[2,3])^2+(fit.summary$coef[3,3])^2
# The value is quite different from the correct Wald test statistic!

##############################################

########################################
## Genotypic model using the characeter X variable
##########################################

fit=glm(Y~X, family=binomial(logit))
summary(fit)

# This is the same as using the dummy variables!

######################################################################
## If we input the data as the table format, then results are the same
######################################################################

E3.data=matrix(c(500,521,350,270,120,130), byrow=T, ncol=2)
dimnames(E3.data)=list(X=c("dd","Dd", "DD"), Y=c("Case","Control"))
E3.data

# additive coding
X=c(0,1,2)
fit=glm(E3.data~X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

# dominant coding
X=c(0,1,1)
fit=glm(E3.data~X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

# recessive coding
X=c(0,0,1)
fit=glm(E3.data~X, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

# genotypic coding, as factor
X=c(0,1,2)
fit=glm(E3.data~as.factor(X), family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

X=c("0","1","2")
fit=glm(E3.data~as.factor(X), family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

X=c("dd","Dd","DD")
fit=glm(E3.data~as.factor(X), family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient

# genotypic coding using dummy variables
X1=c(0,1,0)
X2=c(0,0,1)
fit=glm(E3.data~X1+X2, family=binomial(logit))
fit.summary=summary(fit)
fit.summary$coefficient





