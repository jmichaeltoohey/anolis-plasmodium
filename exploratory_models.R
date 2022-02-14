source('data_carpentry.R')
library(ggplot2)
library("R2jags")

###INFECTION~X###

####SVL males####
model.log=glm(diagnosed.males$infection~diagnosed.males$SVL,family=binomial)
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:104){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
infection.svl.males <- data.frame(y = diagnosed.males$infection, n = nrow(diagnosed.males), x = diagnosed.males$SVL)
jmod=jags(model.file=model,data=infection.svl.males,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(diagnosed.males, aes(SVL, infection))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value=0.7, Std. Error >2x the estimate
#JAGS model does not converge. SD >2x the MU. 
#Positive correlation but weak. Need more sample size.


####Body condition - males####
model.log=glm(diagnosed.males$infection~diagnosed.males$sex.bodycond,family=binomial)
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:104){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
infection.bodycond.males <- data.frame(y = diagnosed.males$infection, n = nrow(diagnosed.males), x = diagnosed.males$sex.bodycond)
jmod=jags(model.file=model,data=infection.bodycond.males,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(diagnosed.males, aes(sex.bodycond, infection))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value = 0.0058**, Std. Error 1/3 the estimate
#JAGS model looks better here. SD 1/3 of Mu
#Negative correlation that is on the verge. Need more sample size.

####Body condition - females####
model.log=glm(diagnosed.females$infection~diagnosed.females$sex.bodycond,family=binomial)
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:48){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
infection.bodycond.females <- data.frame(y = diagnosed.females$infection, n = nrow(diagnosed.females), x = diagnosed.females$sex.bodycond)
jmod=jags(model.file=model,data=infection.bodycond.females,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(diagnosed.females, aes(sex.bodycond, infection))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value = 0.59, Std. Error 2x the estimate
#No correlation




####Temp Difference####
model.log=glm(diagnosed$infection~diagnosed$temp.diff,family=binomial)
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:148){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
infection.tempdiff <- na.omit(data.frame(y = diagnosed$infection, n = nrow(diagnosed), x = diagnosed$temp.diff))
jmod=jags(model.file=model,data=infection.tempdiff,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(diagnosed, aes(temp.diff, infection))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value=0.36, Std. Error ~= the estimate
#JAGS model does not converge. SD ~= the MU. 
#Negative correlation confounded by urban lizards basking on metal






#Chi-squareds (categorical data)

###HABITAT####
table.habitat.infection <- table(diagnosed$infection, diagnosed$habitat) #Make 2-way contingency table with diseased or not and deer species
chisq.test(table.habitat.infection,correct=F) #run Chi squared test
mosaicplot(table.habitat.infection) #make a mosaic plot
#Fisher’s exact test for small sample sizes
fisher.test(table.habitat.infection)

#X^2 p-value=0.51
#Fisher test p-value=0.57
#Positive correlation to urban infections, but weak. Need more sample size.

###SEX####
table.sex.infection <- table(diagnosed$infection, diagnosed$SEX) #Make 2-way contingency table with diseased or not and deer species
chisq.test(table.sex.infection,correct=F) #run Chi squared test
mosaicplot(table.sex.infection) #make a mosaic plot
#Fisher’s exact test for small sample sizes
fisher.test(table.sex.infection)

#X^2 p-value=0.69
#Fisher test p-value=0.84
#No real correlation



#

###BLOOD CHEMISTRIES
#Na####
Sodium.omit <- filter(diagnosed, !is.na(Na))
model.log=glm(Sodium.omit$Na~Sodium.omit$infection,family=Gamma(link = "inverse"))
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:63){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}
infection.Na <- data.frame(y = Sodium.omit$Na, x = Sodium.omit$infection)
jmod=jags(model.file=model,data=infection.Na,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(Sodium.omit, aes(infection, Na))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'Gamma'))

#P-value=0.5, Std. Error > the estimate
#JAGS model does not converge. SD almost 2x the Mu. 
#Positive correlation but weak. Need more sample size.


#K####
K.omit <- filter(diagnosed, !is.na(K))
model.log=glm(K.omit$K~K.omit$infection,family=Gamma(link = "inverse"))
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:58){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
infection.K <- data.frame(y = K.omit$infection, n = nrow(K.omit), x = K.omit$K)
jmod=jags(model.file=model,data=infection.K,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(K.omit, aes(K, infection))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value=0.4, Std. Error > the estimate
#JAGS model does not converge. SD > the Mu. 
#Positive correlation but weak. Need more sample size.


#Cl####
Cl.omit <- filter(diagnosed, !is.na(Cl))
model.log=glm(Cl.omit$Cl~Cl.omit$infection,family=Gamma(link = "inverse"))
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:29){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
infection.Cl <- data.frame(y = Cl.omit$infection, n = nrow(Cl.omit), x = Cl.omit$Cl)
jmod=jags(model.file=model,data=infection.Cl,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(Cl.omit, aes(Cl, infection))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value=0.2, Std. Error ~= the estimate
#JAGS has high variance. SD .7 Mu 
#Positive correlation but weak. Need more sample size.


#Glu####
Glu.omit <- filter(diagnosed, !is.na(Glu))
model.log=glm(Glu.omit$Glu~Glu.omit$infection,family=Gamma(link = "inverse"))
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:59){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
infection.Glu <- data.frame(y = Glu.omit$infection, n = nrow(Glu.omit), x = Glu.omit$Glu)
jmod=jags(model.file=model,data=infection.Glu,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(Glu.omit, aes(Glu, infection))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value=0.2
#JAGS has high variance. SD = Mu 
#Negative correlation but weak. Need more sample size.


#Hct####
Hct.omit <- filter(diagnosed, !is.na(Hct))
model.log = glm(Hct.omit$Hct ~ Hct.omit$infection, family=Gamma(link = "inverse"))
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:60){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}
infection.Hct <- data.frame(y = Hct.omit$Na, x = Hct.omit$infection)
jmod=jags(model.file=model,data=infection.Hct,n.iter=100000,
          parameters.to.save=c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(Hct.omit, aes(Hct, infection))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value = 0.3, Std. Error = the estimate
#JAGS has high variance. SD > Mu 
#Negative correlation but weak. Need more sample size.


###Multiple Regression

#Hb####
Hb.omit <- filter(diagnosed, !is.na(Hb))
model.log = glm(Hb.omit$Hct ~ Hb.omit$infection, family=Gamma(link = "inverse"))
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:60){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
infection.Hct <- data.frame(y = Hct.omit$infection, n = nrow(Hct.omit), x = Hct.omit$Hct)
jmod=jags(model.file = model, data = infection.Hct, n.iter = 100000,
          parameters.to.save = c("b0","b1"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(Hct.omit, aes(infection, Hct))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'Gamma'))

#P-value = 0.3, Std. Error = the estimate
#JAGS has high variance. SD > Mu 
#Negative correlation but weak. Need more sample size.





#Multiple regressions

##SVL+SEX####
model.log=glm(diagnosed$infection~diagnosed$SVL*diagnosed$SEX,family=binomial)
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:152){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x1[i] + b2*x2[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
}
infection.SVLplusSex <- na.omit(data.frame(y = diagnosed$infection, n = nrow(diagnosed), x1 = diagnosed$SVL, x2 = diagnosed$SEX))
jmod=jags(model.file=model,data=infection.SVLplusSex,n.iter=100000,
          parameters.to.save=c("b0","b1","b2"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(diagnosed, aes(SVL, infection, col=SEX))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value=0.36, Std. Error ~= the estimate
#JAGS model does not converge. SD ~= the MU. 
#Negative correlation confounded by urban lizards basking on metal

##SVL+Habitat####
model.log=glm(diagnosed$infection~diagnosed$SVL*diagnosed$habitat,family=binomial)
summary(model.log)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:152){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- b0 + b1*x1[i] + b2*x2[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
}
infection.SVLplusSex <- na.omit(data.frame(y = diagnosed$infection, n = nrow(diagnosed), x1 = diagnosed$SVL, x2 = diagnosed$SEX))
jmod=jags(model.file=model,data=infection.SVLplusSex,n.iter=100000,
          parameters.to.save=c("b0","b1","b2"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(diagnosed, aes(SVL, infection, col=SEX))+geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))

#P-value=0.36, Std. Error ~= the estimate
#JAGS model does not converge. SD ~= the MU. 
#Negative correlation confounded by urban lizards basking on metal
####Body condition + SEX####
model.bodycond=lm(diagnosed$bodycond~diagnosed$infection*diagnosed$SEX)
summary(model.bodycond)
#coefficient is log odds, so you need to transform
exp(coef(model.log)[2])

#bayesian parameterization of binomial regression
model <- function() {
  ## Specify likelihood
  for(i in 1:152){
    y[i] ~ dnorm(mu[i], prec)
    mu[i] <- b0 + b1*x1[i] + b2*x2[i]
  }
  ## Specify priors
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  prec ~ dgamma(1,1)
}
infection.bodycond <- data.frame(y = diagnosed$bodycond, 
                                 x1 = diagnosed$infection, x2 = diagnosed$SEX)
jmod=jags(model.file=model,data=infection.bodycond,n.iter=100000,
          parameters.to.save=c("b0","b1","b2"))
jmod
exp(jmod$BUGSoutput$median$b1)
traceplot(jmod)

ggplot(diagnosed, aes(infection, bodycond, col=SEX))+geom_point()+geom_smooth(method='lm')

#P-value = 0.0058**, Std. Error 1/3 the estimate
#JAGS model looks better here. SD 1/3 of Mu
#Negative correlation that is on the verge. Need more sample size.
#***Significant effect of infection on body condition, maybe Sex on body condition as well


#substrate
model.svl = glm(diagnosed$infection ~ diagnosed$SUBSRATE,family=binomial)
summary(model.svl)


