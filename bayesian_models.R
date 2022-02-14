source('data_carpentry.R')
library(ggplot2)
library("R2jags")
library('rjags')
library(runjags)
library("BayesPostEst")

for(i in 1:nrow(diagnosed)){
  if(diagnosed$SEX[i] == "F"){diagnosed$SEX[i] = 0}
  else if(diagnosed$SEX[i] == "M"){diagnosed$SEX[i] = 1}}

for(i in 1:nrow(diagnosed)){
  if(diagnosed$habitat[i] == "Non-Forested"){diagnosed$habitat[i] = 0}
  else if(diagnosed$habitat[i] == "Forested"){diagnosed$habitat[i] = 1}}


####BAYESIAN PARAMETERIZATIONS####
set.seed(123)
#bayesian parameterization of binomial regression w/ 1 predictor
#model_binom_1 <- function() {
#  ## Specify likelihood
#  for(i in 1:length(x1)){
#    y[i] ~ statip::dbern(p[i])
#    logit(p[i]) <- mu[i]
#    mu[i] <- alpha + b1*x1[i]
#    LogLik[i] <- log(statip::dbern(y[i], p[i]))
#  }
#  ## Specify priors
#  alpha ~ dnorm(0, 0.0001)
#  b1 ~ dnorm(0, 0.0001)
#}

model_binom_1 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    #y[i] ~ dbin(p[i], 1)
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- mu[i]
    mu[i] <- alpha + b1*x1[i]
    #LogLik[i] <- log(dbin(y[i], p[i], 1))
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  #alpha ~ dunif(-10,10)
  #b1 ~ dunif(-10,10)
}

#bayesian parameterization of binomial regression w/ 2 predictors (interaction)
#model_binom_2mult <- function() {
#  ## Specify likelihood
#  for(i in 1:length(x1)){
#    y[i] ~ dbern(p[i])
#    logit(p[i]) <- alpha + b1*x1[i] + b2*x2[i] + b3*x1[i]*x2[i]
#    LogLik[i] <- log(dbern(y[i], p[i]))
#  }
#  ## Specify priors
#  alpha ~ dnorm(0, 0.0001)
#  b1 ~ dnorm(0, 0.0001)
#  b2 ~ dnorm(0, 0.0001)
#  b3 ~ dnorm(0, 0.0001)
#}

model_binom_2mult <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + b2*x2[i] + b3*x1[i]*x2[i]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
}

#bayesian parameterization of binomial regression w/ 2 predictors (additive)
model_binom_2add <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + b2*x2[i]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
}

#bayesian parameterization of binomial regression w/ 2 predictors (interaction)
model_binom_3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    #y[i] ~ dbin(p[i], 1)
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x1[i]*x2[i] +
      b5*x1[i]*x3[i]
    #LogLik[i] <- log(dbin(y[i], p[i], 1))
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
}

#bayesian parameterization of linear regression w/ 1 predictor
model_linear_1 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#bayesian parameterization of linear regression w/ 2 predictor (interaction)
model_linear_2mult <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x1[i]*x2[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#bayesian parameterization of linear regression w/ 2 predictor (additive)
model_linear_2add <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#bayesian parameterization of linear regression w/ 3 predictor
model_linear_3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x1[i]*x2[i] + b5*x1[i]*x3[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#bayesian parameterization of linear regression w/ 3 predictor + interaction
model_linear_3.2 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x2[i]*x3[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#bayesian parameterization of linear regression w/ 3 predictor all additive
model_linear_3.3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#bayesian parameterization of linear regression w/ 4 predictors
model_linear_3.4 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i] + 
      b5*x2[i]*x3[i] + b6*x2[i]*x4[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
  b6 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#bayesian parameterization of gamma regression w/ 1 predictor
model_gamma_1 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- alpha + b1*x1[i]
    LogLik[i] <- log(dgamma(y[i], shape, shape / exp(p[i])))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}

#bayesian parameterization of gamma regression w/ 2 predictors (additive)
model_gamma_2add <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- alpha + b1*x1[i] + b2*x2[i]
    LogLik[i] <- log(dgamma(y[i], shape, shape / exp(p[i])))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}

#bayesian parameterization of gamma regression w/ 2 predictors (interaction)
model_gamma_2mult <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x1[i]*x2[i]
    LogLik[i] <- log(dgamma(y[i], shape, shape / exp(p[i])))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}

#bayesian parameterization of gamma regression w/ 3 predictors
model_gamma_3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x1[i]*x2[i] + b5*x1[i]*x3[i]
    LogLik[i] <- log(dgamma(y[i], shape, shape / exp(p[i])))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}

#bayesian parameterization of gamma regression w/ 3 predictors (interaction)
model_gamma_3.2 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x2[i]*x3[i]
    LogLik[i] <- log(dgamma(y[i], shape, shape / exp(p[i])))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}

#bayesian parameterization of gamma regression w/ 3 predictors all additive
model_gamma_3.3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i]
    LogLik[i] <- log(dgamma(y[i], shape, shape / exp(p[i])))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}

#bayesian parameterization of gamma regression w/ 4 predictors
model_gamma_3.4 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i] + 
      b5*x2[i]*x3[i] + b6*x2[i]*x4[i]
    LogLik[i] <- log(dgamma(y[i], shape, shape / exp(p[i])))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
  b6 ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}


ni <- 50000
nt <- 3
nb <- 5000
nc <- 3



####p(infection)~####


#p(infection)~1
infection.null <- data.frame(y = diagnosed$infection, n = nrow(diagnosed), x1 = 1)
mod.infection.null=jags(model.file=model_binom_1,data=infection.null,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("b0","b1","LogLik"))
mod.infection.null
exp(mod.infection.null$BUGSoutput$DIC)

#p(infection)  ~ SVL
infection.svl <- data.frame(y = diagnosed$infection, n = nrow(diagnosed), x1 = diagnosed$SVL)
mod.infection.svl=jags(model.file=model_binom_1,data=infection.svl,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","LogLik"))
#mod.infection.svl=jags(model.file=model_binom_1,data=infection.svl,
#     n.chains = nc,  n.iter = ni, n.burnin = nb,
#    parameters.to.save=c("b0","b1"))
mod.infection.svl
mcmcTab(mod.infection.svl)

#BAYESIAN REGRESSION TABLE
mcmcTab(mod.infection.svl)
#ADD proportion of positive/negative draws
mcmcTab(mod.infection.svl, Pr=TRUE)

#mcmcAveProb (calculates predicted probability at pre-defined values of one covar of interest)
#create matrix of posterior draws of coefficients to pass onto function
mcmcmat.jags <- as.matrix(coda::as.mcmc(mod.infection.svl))
#generate model matrix to pass onto function
#contains as many columns as estimated regression coefficients
#first column is vector of 1s (intercept)
#remaining columns are observed values of covariates in the model
#NOTE: order of columns in model matrix must correspond to order in posterior draws
mm <- model.matrix(infection ~ SVL,
                   data = diagnosed)
#now we generate full posterior distributions of predicted probabilities for both habitats
aveprob.inf.svl.jags <- mcmcAveProb(modelmatrix = mm,
                                   mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                                   xcol = 2,
                                   xrange = c(3.25, 6.75),
                                   link = "logit",
                                   ci = c(0.025, 0.975),
                                   fullsims = FALSE)
library("ggplot2")
library("ggridges")
ggplot(data = aveprob.inf.svl.jags, 
       aes(x = x, y = median_pp)) + 
  geom_ribbon(aes(ymin = lower_pp, ymax = upper_pp), fill = "gray") + 
  geom_line() + 
  xlab("SVL") + 
  ylab("Estimated probability of infection") + 
  ylim(0, .01) + #TOO SMALL!!!!
  labs(title = "Probability based on SVL") +
  theme_minimal()

#p(infection)  ~ Sex
infection.sex <- data.frame(y = diagnosed$infection, n = nrow(diagnosed), x1 = diagnosed$SEX)
mod.infection.sex=jags(model.file=model_binom_1,data=infection.sex,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","LogLik"))

mod.infection.sex
exp(mod.infection.sex$BUGSoutput$DIC)

#p(infection)  ~ Habitat
infection.habitat <- data.frame(y = diagnosed$infection, 
                                n = nrow(diagnosed), 
                                x1 = diagnosed$habitat)
mod.infection.habitat=jags(model.file=model_binom_1,data=infection.habitat,
                           n.chains = nc,  n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","LogLik"))
#mod.infection.habitat=jags(model.file=model_binom_1,data=infection.habitat,
#                           n.chains = nc,  n.iter = ni, n.burnin = nb,
#                           parameters.to.save=c("b0","b1"))

##jags(model.file=model_binom_1,data=infection.habitat,
#                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
 #                          parameters.to.save=c("b0","b1"))
mod.infection.habitat
exp(mod.infection.habitat$BUGSoutput$DIC)

library(BayesPostEst)
#BAYESIAN REGRESSION TABLE
mcmcTab(mod.infection.habitat)
#ADD proportion of positive/negative draws
mcmcTab(mod.infection.habitat, Pr=TRUE)

#mcmcAveProb (calculates predicted probability at pre-defined values of one covar of interest)
#create matrix of posterior draws of coefficients to pass onto function
mcmcmat.jags <- as.matrix(coda::as.mcmc(mod.infection.habitat))
#generate model matrix to pass onto function
#contains as many columns as estimated regression coefficients
#first column is vector of 1s (intercept)
#remaining columns are observed values of covariates in the model
#NOTE: order of columns in model matrix must correspond to order in posterior draws
mm <- model.matrix(infection ~ habitat,
                   data = diagnosed)
#now we generate full posterior distributions of predicted probabilities for both habitats
aveprob.forest.jags <- mcmcAveProb(modelmatrix = mm,
                                   mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                                   xcol = 2,
                                   xrange = c(0, 1),
                                   link = "logit",
                                   ci = c(0.025, 0.975),
                                   fullsims = TRUE)
library("ggplot2")
library("ggridges")
ggplot(data = aveprob.forest.jags, 
       aes(y = factor(x), x = pp)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), vline_color = "white") + 
  scale_y_discrete(labels = c("Forested", "Deforested")) + 
  ylab("") + 
  xlab("Estimated probability of volunteering") + 
  labs(title = "Probability based on average-case approach") +
  theme_minimal()

#p(infection)~(SVL*Sex)
infection.svl.sex <- data.frame(y = diagnosed$infection, n = nrow(diagnosed), 
                                x1 = diagnosed$SVL, x2 = diagnosed$SEX)
mod.infection.svl.sex=jags(model.file=model_binom_2mult,data=infection.svl.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.infection.svl.sex
exp(mod.infection.svl.sex$BUGSoutput$DIC)

#p(infection)~(SVL*Habitat)
infection.svl.habitat <- data.frame(y = diagnosed$infection, 
                                    n = nrow(diagnosed), 
                                x1 = diagnosed$SVL, x2 = diagnosed$habitat)
mod.infection.svl.habitat=jags(model.file=model_binom_2mult,data=infection.svl.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("b0","b1","b2","b3","LogLik"))
output.mod.infection.svl.habitat=jags(model.file=model_binom_2mult,data=infection.svl.habitat,
                               n.chains = nc,  n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("b0","b1","b2","b3"))
traceplot(output.mod.infection.svl.habitat)

library(BayesPostEst)
#BAYESIAN REGRESSION TABLE
mcmcTab(output.mod.infection.svl.habitat)
#ADD proportion of positive/negative draws
mcmcTab(output.mod.infection.svl.habitat, Pr=TRUE)

#mcmcAveProb (calculates predicted probability at pre-defined values of one covar of interest)
#create matrix of posterior draws of coefficients to pass onto function
mcmcmat.jags <- as.matrix(coda::as.mcmc(output.mod.infection.svl.habitat))
#generate model matrix to pass onto function
#contains as many columns as estimated regression coefficients
#first column is vector of 1s (intercept)
#remaining columns are observed values of covariates in the model
#NOTE: order of columns in model matrix must correspond to order in posterior draws
mm <- model.matrix(infection ~ SVL + habitat + SVL*habitat,
                   data = diagnosed)
#now we generate full posterior distributions of predicted probabilities for both habitats
aveprob.forest.jags <- mcmcAveProb(modelmatrix = mm,
                                   mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                                   xcol = 3,
                                   xrange = c(0, 1),
                                   link = "logit",
                                   ci = c(0.025, 0.975),
                                   fullsims = TRUE)
library("ggplot2")
library("ggridges")
ggplot(data = aveprob.forest.jags, 
       aes(y = factor(x), x = pp)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), vline_color = "white") + 
  scale_y_discrete(labels = c("Forested", "Deforested")) + 
  ylab("") + 
  xlab("Estimated probability of volunteering") + 
  labs(title = "Probability based on average-case approach") +
  theme_minimal()

#(mfrow=c(3,2))
#traceplot(mod.infection.svl.habitat)

#library(plotrix)
#b0.est=mod.infection.svl.habitat$BUGSoutput$mean$b0
#=mod.infection.svl.habitat$BUGSoutput$mean$b1
#b2.est=mod.infection.svl.habitat$BUGSoutput$mean$b2
#b3.est=mod.infection.svl.habitat$BUGSoutput$mean$b3


#parms.est=c(b0.est, b1.est, b2.est, b3.est)

#CI = apply(mod.infection.svl.habitat$BUGSoutput$sims.matrix[,-5], 2, quantile, c(0.025, 0.975))
#parms.uncertain=matrix(NA, 3, 4)
#parms.uncertain[1,]=CI[1,]
#parms.uncertain[2,]=parms.est
#parms.uncertain[3,]=CI[2,]
#rownames(parms.uncertain)=c("2.5%", "mean", "97.5%")
#colnames(parms.uncertain)=c("b0", "b1", "b2", "b3")

#parms.uncertain

#plotCI(x=1:4, parms.uncertain[2,], li=parms.uncertain[1,], ui=parms.uncertain[3,], col=1:4,
#       xlab="Parameters (b0, b1, b2, b3)", ylab="Value",
#       main="Parameter Estimates with 95% CRI")
#legend(x=0.9, y=2.2, legend=c("Intercept", "SVL", "Habitat", "SVL*Habitat"), box.lty=0, 
#       pch=c(1,1,1,1,1,1), col=1:4, lty=c(1,1,1,1,1,1), cex=0.75, bg=0)
#exp(mod.infection.svl.habitat$BUGSoutput)


#p(infection)~(Sex)+(Habitat)
infection.sex.habitat <- data.frame(y = diagnosed$infection, n = nrow(diagnosed), 
                                    x1 = diagnosed$SEX, x2 = diagnosed$habitat)
mod.infection.sex.habitat=jags(model.file=model_binom_2add,data=infection.sex.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("b0","b1","b2","LogLik"))

#mod.infection.sex.habitat=jags(model.file=model_binom_2add,data=infection.sex.habitat,n.iter=100000,
#                                   parameters.to.save=c("b0","b1","b2"))
#print(mod.infection.sex.habitat)
exp(mod.infection.sex.habitat$BUGSoutput$DIC)



#p(infection)~(SVL*Sex)+(SVL*Habitat)
infection.svl.sex.habitat <- data.frame(y = diagnosed$infection, n = nrow(diagnosed), 
                                    x1 = diagnosed$SVL, x2 = diagnosed$SEX, x3 = diagnosed$habitat)
mod.infection.svl.sex.habitat=jags(model.file=model_binom_3,data=infection.svl.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("b0","b1","b2","b3","b4","b5","LogLik"))
mod.infection.svl.sex.habitat
exp(mod.infection.svl.sex.habitat$BUGSoutput$DIC)
traceplot(mod.infection.svl.sex.habitat)


# drawing samples gives mcarrays
#update(mod.infection.svl.habitat, 100000)
#samples <- jags.samples(mod.infection.svl.habitat, c('b0', 'b1','b2','b3'), 100000)
#str(samples)

####body condition####


#body condition~1
male.bodycond.null <- data.frame(y = diagnosed.males$sex.bodycond, num = nrow(diagnosed.males),  x1 = 1)
mod.male.bodycond.null=jags(model.file=model_linear_1,data=male.bodycond.null,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("b0","b1","LogLik"))
mod.male.bodycond.null
exp(mod.male.bodycond.null$BUGSoutput$DIC)

#body condition ~ infection
male.bodycond.infection <- data.frame(y = diagnosed.males$sex.bodycond, n = nrow(diagnosed.males), x1 = diagnosed.males$infection)
mod.male.bodycond.infection=jags(model.file=model_linear_1,data=male.bodycond.infection,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","LogLik"))
mod.male.bodycond.infection
exp(mod.male.bodycond.infection$BUGSoutput$DIC)

# body condition ~ SVL
male.bodycond.SVL <- data.frame(y = diagnosed.males$sex.bodycond, n = nrow(diagnosed.males), x1 = diagnosed.males$SVL)
mod.male.bodycond.SVL=jags(model.file=model_linear_1,data=male.bodycond.SVL,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("b0","b1","LogLik"))
mod.male.bodycond.SVL
exp(mod.male.bodycond.SVL$BUGSoutput$DIC)


# body condition ~ Habitat
male.bodycond.habitat <- data.frame(y = diagnosed.males$sex.bodycond, n = nrow(diagnosed.males), x1 = diagnosed.males$habitat)
mod.male.bodycond.habitat=jags(model.file=model_linear_1,data=male.bodycond.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("b0","b1","LogLik"))
mod.male.bodycond.habitat
exp(mod.male.bodycond.habitat$BUGSoutput$DIC)


#body condition~(SVL*Habitat)
male.bodycond.SVL.habitat <- data.frame(y = diagnosed.males$sex.bodycond, n = nrow(diagnosed.males), 
                               x1 = diagnosed.males$SVL, x2 = diagnosed.males$habitat)
mod.male.bodycond.SVL.habitat=jags(model.file=model_linear_2mult,data=male.bodycond.SVL.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.male.bodycond.SVL.habitat
exp(mod.male.bodycond.SVL.habitat$BUGSoutput$DIC)


# body condition ~ INFECTION + SVL
male.bodycond.infection.svl <- data.frame(y = diagnosed.males$sex.bodycond,  n = nrow(diagnosed.males),
                                   x1 = diagnosed.males$infection, x2 = diagnosed.males$SVL)
mod.male.bodycond.infection.svl=jags(model.file=model_linear_2add,data=male.bodycond.infection.svl,
                                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                              parameters.to.save=c("b0","b1","b2","LogLik"))
mod.male.bodycond.infection.svl
exp(mod.male.bodycond.infection.svl$BUGSoutput$DIC)


# body condition ~ INFECTION + Habitat
male.bodycond.infection.habitat <- data.frame(y = diagnosed.males$sex.bodycond,  n = nrow(diagnosed.males),
                                     x1 = diagnosed.males$infection, x2 = diagnosed.males$habitat)
mod.male.bodycond.infection.habitat=jags(model.file=model_linear_2add,data=male.bodycond.infection.habitat,
                                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                parameters.to.save=c("b0","b1","b2","LogLik"))
mod.male.bodycond.infection.habitat
exp(mod.male.bodycond.infection.habitat$BUGSoutput$DIC)



#body condition~INFECTION + (SVL*Habitat)
male.bodycond.infection.svl.habitat <- data.frame(y = diagnosed.males$sex.bodycond,  n = nrow(diagnosed.males),
                                         x1 = diagnosed.males$infection, x2 = diagnosed.males$SVL, 
                                         x3 = diagnosed.males$habitat)
mod.male.bodycond.infection.svl.habitat=jags(model.file=model_linear_3,
                                    data=male.bodycond.infection.svl.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.male.bodycond.infection.svl.habitat
exp(mod.male.bodycond.infection.svl.habitat$BUGSoutput$DIC)

#ni <- 100000
#nt <- 6
#nb <- 5000
#nc <- 3
#mod.male.bodycond.infection.svl.habitat=jags(model.file=model_linear_3,
#                                             data=male.bodycond.infection.svl.habitat,
#                                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                             parameters.to.save=c("b0","b1","b2","b3","b4"))
#print(mod.male.bodycond.infection.svl.habitat)
#par(mfrow=c(3,2))
#traceplot(mod.male.bodycond.infection.svl.habitat)

#library(plotrix)
#b0.est=mod.male.bodycond.infection.svl.habitat$BUGSoutput$mean$b0
#b1.est=mod.male.bodycond.infection.svl.habitat$BUGSoutput$mean$b1
#b2.est=mod.male.bodycond.infection.svl.habitat$BUGSoutput$mean$b2
#b3.est=mod.male.bodycond.infection.svl.habitat$BUGSoutput$mean$b3
#b4.est=mod.male.bodycond.infection.svl.habitat$BUGSoutput$mean$b4


#parms.est=c(b0.est, b1.est, b2.est, b3.est, b4.est)

#CI = apply(mod.male.bodycond.infection.svl.habitat$BUGSoutput$sims.matrix[,-6], 2, quantile, c(0.025, 0.975))
#parms.uncertain=matrix(NA, 3, 5)
#parms.uncertain[1,]=CI[1,]
#parms.uncertain[2,]=parms.est
#parms.uncertain[3,]=CI[2,]
#rownames(parms.uncertain)=c("2.5%", "mean", "97.5%")
#colnames(parms.uncertain)=c("b0", "b1", "b2", "b3", "b4")

#parms.uncertain

#plotCI(x=1:5, parms.uncertain[2,], li=parms.uncertain[1,], ui=parms.uncertain[3,], col=1:5,
#       xlab="Parameters (b0, b1, b2, b3, b4)", ylab="Value",
#       main="Parameter Estimates with 95% CRI")
#legend(x=0.9, y=5.2, legend=c("Intercept", "Infection","SVL", "Habitat", "SVL*Habitat"), box.lty=0, 
#       pch=c(1,1,1,1,1,1), col=1:5, lty=c(1,1,1,1,1,1), cex=0.75, bg=0)










####Na####
Sodium.omit <- filter(diagnosed, !is.na(Na))
#Na~1
Na.null <- data.frame(y = Sodium.omit$Na, x1 = 1)
mod.Na.null=jags(model.file=model_gamma_1,data=Na.null,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
          parameters.to.save=c("b0","b1","LogLik"))
mod.Na.null
exp(mod.Na.null$BUGSoutput$DIC)

#Na ~ infection
Na.infection <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$infection)
mod.Na.infection=jags(model.file=model_gamma_1,data=Na.infection,
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("b0","b1","LogLik"))
mod.Na.infection
exp(mod.Na.infection$BUGSoutput$DIC)

# Na ~ SVL
Na.svl <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL)
mod.Na.svl=jags(model.file=model_gamma_1,data=Na.svl,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("b0","b1","LogLik"))
mod.Na.svl
exp(mod.Na.svl$BUGSoutput$DIC)

# Na ~ Sex
Na.sex <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SEX)
mod.Na.sex=jags(model.file=model_gamma_1,data=Na.sex,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("b0","b1","LogLik"))
mod.Na.sex
exp(mod.Na.sex$BUGSoutput$DIC)

# Na ~ Habitat
Na.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$habitat)
mod.Na.habitat=jags(model.file=model_gamma_1,data=Na.habitat,
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("b0","b1","LogLik"))
mod.Na.habitat
exp(mod.Na.habitat$BUGSoutput$DIC)

#Na~(SVL*Sex)
Na.svl.sex <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL, x2 = Sodium.omit$SEX)
mod.Na.svl.sex=jags(model.file=model_gamma_2mult,data=Na.svl.sex,
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Na.svl.sex
exp(mod.Na.svl.sex$BUGSoutput$DIC)

#Na~(SVL*Habitat)
Na.svl.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL, x2 = Sodium.omit$habitat)
mod.Na.svl.habitat=jags(model.file=model_gamma_2mult,data=Na.svl.habitat,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Na.svl.habitat
exp(mod.Na.svl.habitat$BUGSoutput$DIC)

#Na~(Sex)+(Habitat)
Na.sex.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SEX, x2 = Sodium.omit$habitat)
mod.Na.sex.habitat=jags(model.file=model_gamma_2add,data=Na.sex.habitat,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Na.sex.habitat
exp(mod.Na.sex.habitat$BUGSoutput$DIC)

#Na~(SVL*Sex)+(SVL*Habitat)
Na.svl.sex.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL, 
                                 x2 = Sodium.omit$SEX, x3 = Sodium.omit$habitat)
mod.Na.svl.sex.habitat=jags(model.file=model_gamma_3,data=Na.svl.sex.habitat,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                    parameters.to.save=c("b0","b1","b2","b3","b4","b5","LogLik"))
mod.Na.svl.sex.habitat
exp(mod.Na.svl.sex.habitat$BUGSoutput$DIC)

# Na ~ INFECTION + SVL
Na.infection.svl <- data.frame(y = Sodium.omit$Na, 
                               x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL)
mod.Na.infection.svl=jags(model.file=model_gamma_2add,data=Na.infection.svl,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Na.infection.svl
exp(mod.Na.infection.svl$BUGSoutput$DIC)

# Na ~ INFECTION + Sex
Na.infection.sex <- data.frame(y = Sodium.omit$Na, 
                               x1 = Sodium.omit$infection, x2 = Sodium.omit$SEX)
mod.Na.infection.sex=jags(model.file=model_gamma_2add,data=Na.infection.sex,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Na.infection.sex
exp(mod.Na.infection.sex$BUGSoutput$DIC)

# Na ~ INFECTION + Habitat
Na.infection.habitat <- data.frame(y = Sodium.omit$Na, 
                               x1 = Sodium.omit$infection, x2 = Sodium.omit$habitat)
mod.Na.infection.habitat=jags(model.file=model_gamma_2add,data=Na.infection.habitat,
                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Na.infection.habitat
exp(mod.Na.infection.habitat$BUGSoutput$DIC)

#Na~INFECTION + (SVL*Sex)
Na.infection.svl.sex <- data.frame(y = Sodium.omit$Na, 
                               x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL,
                               x3 = Sodium.omit$SEX)
mod.Na.infection.svl.sex=jags(model.file=model_gamma_3.2,data=Na.infection.svl.sex,
                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.Na.infection.svl.sex
exp(mod.Na.infection.svl.sex$BUGSoutput$DIC)

#Na~INFECTION + (SVL*Habitat)
Na.infection.svl.habitat <- data.frame(y = Sodium.omit$Na, 
                                   x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL,
                                   x3 = Sodium.omit$habitat)
mod.Na.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=Na.infection.svl.habitat,
                                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                              parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.Na.infection.svl.habitat
exp(mod.Na.infection.svl.habitat$BUGSoutput$DIC)

#Na~INFECTION + (Sex)+(Habitat)
Na.infection.sex.habitat <- data.frame(y = Sodium.omit$Na, 
                                       x1 = Sodium.omit$infection, x2 = Sodium.omit$SEX,
                                       x3 = Sodium.omit$habitat)
mod.Na.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=Na.infection.sex.habitat,
                                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                  parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Na.infection.sex.habitat
exp(mod.Na.infection.sex.habitat$BUGSoutput$DIC)

#Na~INFECTION + (SVL*Sex)+(SVL*Habitat)
Na.infection.svl.sex.habitat <- data.frame(y = Sodium.omit$Na, 
                                       x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL,
                                       x3 = Sodium.omit$SEX, x4 = Sodium.omit$habitat)
mod.Na.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=Na.infection.svl.sex.habitat,
                                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                  parameters.to.save=c("b0","b1","b2","b3",
                                                       "b4","b5","b6","LogLik"))
mod.Na.infection.svl.sex.habitat
exp(mod.Na.infection.svl.sex.habitat$BUGSoutput$DIC)


####K####
K.omit <- filter(diagnosed, !is.na(K))
#K~1
K.null <- data.frame(y = K.omit$K, x1 = 1)
mod.K.null=jags(model.file=model_gamma_1,data=K.null,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.K.null
exp(mod.K.null$BUGSoutput$DIC)

#K ~ infection
K.infection <- data.frame(y = K.omit$K, x1 = K.omit$infection)
mod.K.infection=jags(model.file=model_gamma_1,data=K.infection,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.K.infection
exp(mod.K.infection$BUGSoutput$DIC)

# K ~ SVL
K.svl <- data.frame(y = K.omit$K, x1 = K.omit$SVL)
mod.K.svl=jags(model.file=model_gamma_1,data=K.svl,
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.K.svl
exp(mod.K.svl$BUGSoutput$DIC)

# K ~ Sex
K.sex <- data.frame(y = K.omit$K, x1 = K.omit$SEX)
mod.K.sex=jags(model.file=model_gamma_1,data=K.sex,
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.K.sex
exp(mod.K.sex$BUGSoutput$DIC)

# K ~ Habitat
K.habitat <- data.frame(y = K.omit$K, x1 = K.omit$habitat)
mod.K.habitat=jags(model.file=model_gamma_1,data=K.habitat,
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                    parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.K.habitat
exp(mod.K.habitat$BUGSoutput$DIC)

#K~(SVL*Sex)
K.svl.sex <- data.frame(y = K.omit$K, x1 = K.omit$SVL, x2 = K.omit$SEX)
mod.K.svl.sex=jags(model.file=model_gamma_2mult,data=K.svl.sex,
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                    parameters.to.save=c("b0","b1","b2","b3","LogLik","LogLik"))
mod.K.svl.sex
exp(mod.K.svl.sex$BUGSoutput$DIC)

#K~(SVL*Habitat)
K.svl.habitat <- data.frame(y = K.omit$K, x1 = K.omit$SVL, x2 = K.omit$habitat)
mod.K.svl.habitat=jags(model.file=model_gamma_2mult,data=K.svl.habitat,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("b0","b1","b2","b3","LogLik","LogLik"))
mod.K.svl.habitat
exp(mod.K.svl.habitat$BUGSoutput$DIC)

#K~(Sex)+(Habitat)
K.sex.habitat <- data.frame(y = K.omit$K, x1 = K.omit$SEX, x2 = K.omit$habitat)
mod.K.sex.habitat=jags(model.file=model_gamma_2add,data=K.sex.habitat,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("b0","b1","b2","LogLik","LogLik"))
mod.K.sex.habitat
exp(mod.K.sex.habitat$BUGSoutput$DIC)

#K~(SVL*Sex)+(SVL*Habitat)
K.svl.sex.habitat <- data.frame(y = K.omit$K, x1 = K.omit$SVL, 
                                 x2 = K.omit$SEX, x3 = K.omit$habitat)
mod.K.svl.sex.habitat=jags(model.file=model_gamma_3,data=K.svl.sex.habitat,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("b0","b1","b2","b3","b4","b5","LogLik","LogLik"))
mod.K.svl.sex.habitat
exp(mod.K.svl.sex.habitat$BUGSoutput$DIC)

# K ~ INFECTION + SVL
K.infection.svl <- data.frame(y = K.omit$K, 
                               x1 = K.omit$infection, x2 = K.omit$SVL)
mod.K.infection.svl=jags(model.file=model_gamma_2add,data=K.infection.svl,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("b0","b1","b2","LogLik","LogLik"))
mod.K.infection.svl
exp(mod.K.infection.svl$BUGSoutput$DIC)

# K ~ INFECTION + Sex
K.infection.sex <- data.frame(y = K.omit$K, 
                               x1 = K.omit$infection, x2 = K.omit$SEX)
mod.K.infection.sex=jags(model.file=model_gamma_2add,data=K.infection.sex,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("b0","b1","b2","LogLik"))
mod.K.infection.sex
exp(mod.K.infection.sex$BUGSoutput$DIC)

# K ~ INFECTION + Habitat
K.infection.habitat <- data.frame(y = K.omit$K, 
                                   x1 = K.omit$infection, x2 = K.omit$habitat)
mod.K.infection.habitat=jags(model.file=model_gamma_2add,data=K.infection.habitat,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                              parameters.to.save=c("b0","b1","b2","LogLik"))
mod.K.infection.habitat
exp(mod.K.infection.habitat$BUGSoutput$DIC)

#K~INFECTION + (SVL*Sex)
K.infection.svl.sex <- data.frame(y = K.omit$K, 
                                   x1 = K.omit$infection, x2 = K.omit$SVL,
                                   x3 = K.omit$SEX)
mod.K.infection.svl.sex=jags(model.file=model_gamma_3.2,data=K.infection.svl.sex,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                              parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.K.infection.svl.sex
exp(mod.K.infection.svl.sex$BUGSoutput$DIC)

#K~INFECTION + (SVL*Habitat)
K.infection.svl.habitat <- data.frame(y = K.omit$K, 
                                       x1 = K.omit$infection, x2 = K.omit$SVL,
                                       x3 = K.omit$habitat)
mod.K.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=K.infection.svl.habitat,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                  parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.K.infection.svl.habitat
exp(mod.K.infection.svl.habitat$BUGSoutput$DIC)

#K~INFECTION + (Sex)+(Habitat)
K.infection.sex.habitat <- data.frame(y = K.omit$K, 
                                       x1 = K.omit$infection, x2 = K.omit$SEX,
                                       x3 = K.omit$habitat)
mod.K.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=K.infection.sex.habitat,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                  parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.K.infection.sex.habitat
exp(mod.K.infection.sex.habitat$BUGSoutput$DIC)

#K~INFECTION + (SVL*Sex)+(SVL*Habitat)
K.infection.svl.sex.habitat <- data.frame(y = K.omit$K, 
                                           x1 = K.omit$infection, x2 = K.omit$SVL,
                                           x3 = K.omit$SEX, x4 = K.omit$habitat)
mod.K.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=K.infection.svl.sex.habitat,
                                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                      parameters.to.save=c("b0","b1","b2","b3",
                                                           "b4","b5","b6","LogLik"))
mod.K.infection.svl.sex.habitat
exp(mod.K.infection.svl.sex.habitat$BUGSoutput$DIC)


####Cl####
Cl.omit <- filter(diagnosed, !is.na(Cl))
#Cl~1
Cl.null <- data.frame(y = Cl.omit$Cl, x1 = 1)
mod.Cl.null=jags(model.file=model_gamma_1,data=Cl.null,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("b0","b1","LogLik"))
mod.Cl.null
exp(mod.Cl.null$BUGSoutput$DIC)

#Cl ~ infection
Cl.infection <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$infection)
mod.Cl.infection=jags(model.file=model_gamma_1,data=Cl.infection,
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("b0","b1","LogLik"))
mod.Cl.infection
exp(mod.Cl.infection$BUGSoutput$DIC)

# Cl ~ SVL
Cl.svl <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL)
mod.Cl.svl=jags(model.file=model_gamma_1,data=Cl.svl,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
               parameters.to.save=c("b0","b1","LogLik"))
mod.Cl.svl
exp(mod.Cl.svl$BUGSoutput$DIC)

# Cl ~ Sex
Cl.sex <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SEX)
mod.Cl.sex=jags(model.file=model_gamma_1,data=Cl.sex,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
               parameters.to.save=c("b0","b1","LogLik"))
mod.Cl.sex
exp(mod.Cl.sex$BUGSoutput$DIC)

# Cl ~ Habitat
Cl.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$habitat)
mod.Cl.habitat=jags(model.file=model_gamma_1,data=Cl.habitat,
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                   parameters.to.save=c("b0","b1","LogLik"))
mod.Cl.habitat
exp(mod.Cl.habitat$BUGSoutput$DIC)

#Cl~(SVL*Sex)
Cl.svl.sex <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL, x2 = Cl.omit$SEX)
mod.Cl.svl.sex=jags(model.file=model_gamma_2mult,data=Cl.svl.sex,
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                   parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Cl.svl.sex
exp(mod.Cl.svl.sex$BUGSoutput$DIC)

#Cl~(SVL*Habitat)
Cl.svl.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL, x2 = Cl.omit$habitat)
mod.Cl.svl.habitat=jags(model.file=model_gamma_2mult,data=Cl.svl.habitat,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Cl.svl.habitat
exp(mod.Cl.svl.habitat$BUGSoutput$DIC)

#Cl~(Sex)+(Habitat)
Cl.sex.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SEX, x2 = Cl.omit$habitat)
mod.Cl.sex.habitat=jags(model.file=model_gamma_2add,data=Cl.sex.habitat,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Cl.sex.habitat
exp(mod.Cl.sex.habitat$BUGSoutput$DIC)

#Cl~(SVL*Sex)+(SVL*Habitat)
Cl.svl.sex.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL, 
                                x2 = Cl.omit$SEX, x3 = Cl.omit$habitat)
mod.Cl.svl.sex.habitat=jags(model.file=model_gamma_3,data=Cl.svl.sex.habitat,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2","b3","b4","b5","LogLik"))
mod.Cl.svl.sex.habitat
exp(mod.Cl.svl.sex.habitat$BUGSoutput$DIC)

# Cl ~ INFECTION + SVL
Cl.infection.svl <- data.frame(y = Cl.omit$Cl, 
                              x1 = Cl.omit$infection, x2 = Cl.omit$SVL)
mod.Cl.infection.svl=jags(model.file=model_gamma_2add,data=Cl.infection.svl,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Cl.infection.svl
exp(mod.Cl.infection.svl$BUGSoutput$DIC)

# Cl ~ INFECTION + Sex
Cl.infection.sex <- data.frame(y = Cl.omit$Cl, 
                              x1 = Cl.omit$infection, x2 = Cl.omit$SEX)
mod.Cl.infection.sex=jags(model.file=model_gamma_2add,data=Cl.infection.sex,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Cl.infection.sex
exp(mod.Cl.infection.sex$BUGSoutput$DIC)

# Cl ~ INFECTION + Habitat
Cl.infection.habitat <- data.frame(y = Cl.omit$Cl, 
                                  x1 = Cl.omit$infection, x2 = Cl.omit$habitat)
mod.Cl.infection.habitat=jags(model.file=model_gamma_2add,data=Cl.infection.habitat,
                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Cl.infection.habitat
exp(mod.Cl.infection.habitat$BUGSoutput$DIC)

#Cl~INFECTION + (SVL*Sex)
Cl.infection.svl.sex <- data.frame(y = Cl.omit$Cl, 
                                  x1 = Cl.omit$infection, x2 = Cl.omit$SVL,
                                  x3 = Cl.omit$SEX)
mod.Cl.infection.svl.sex=jags(model.file=model_gamma_3.2,data=Cl.infection.svl.sex,
                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.Cl.infection.svl.sex
exp(mod.Cl.infection.svl.sex$BUGSoutput$DIC)

#Cl~INFECTION + (SVL*Habitat)
Cl.infection.svl.habitat <- data.frame(y = Cl.omit$Cl, 
                                      x1 = Cl.omit$infection, x2 = Cl.omit$SVL,
                                      x3 = Cl.omit$habitat)
mod.Cl.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=Cl.infection.svl.habitat,
                                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.Cl.infection.svl.habitat
exp(mod.Cl.infection.svl.habitat$BUGSoutput$DIC)

#Cl~INFECTION + (Sex)+(Habitat)
Cl.infection.sex.habitat <- data.frame(y = Cl.omit$Cl, 
                                      x1 = Cl.omit$infection, x2 = Cl.omit$SEX,
                                      x3 = Cl.omit$habitat)
mod.Cl.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=Cl.infection.sex.habitat,
                                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Cl.infection.sex.habitat
exp(mod.Cl.infection.sex.habitat$BUGSoutput$DIC)

#Cl~INFECTION + (SVL*Sex)+(SVL*Habitat)
Cl.infection.svl.sex.habitat <- data.frame(y = Cl.omit$Cl, 
                                          x1 = Cl.omit$infection, x2 = Cl.omit$SVL,
                                          x3 = Cl.omit$SEX, x4 = Cl.omit$habitat)
mod.Cl.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=Cl.infection.svl.sex.habitat,
                                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                     parameters.to.save=c("b0","b1","b2","b3",
                                                          "b4","b5","b6","LogLik"))
mod.Cl.infection.svl.sex.habitat
exp(mod.Cl.infection.svl.sex.habitat$BUGSoutput$DIC)


####iCa####
iCa.omit <- filter(diagnosed, !is.na(iCa))
#iCa~1
iCa.null <- data.frame(y = iCa.omit$iCa, x1 = 1)
mod.iCa.null=jags(model.file=model_gamma_1,data=iCa.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.iCa.null
exp(mod.iCa.null$BUGSoutput$DIC)

#iCa ~ infection
iCa.infection <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$infection)
mod.iCa.infection=jags(model.file=model_gamma_1,data=iCa.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.iCa.infection
exp(mod.iCa.infection$BUGSoutput$DIC)

# iCa ~ SVL
iCa.svl <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL)
mod.iCa.svl=jags(model.file=model_gamma_1,data=iCa.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
               parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.iCa.svl
exp(mod.iCa.svl$BUGSoutput$DIC)

# iCa ~ Sex
iCa.sex <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SEX)
mod.iCa.sex=jags(model.file=model_gamma_1,data=iCa.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
               parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.iCa.sex
exp(mod.iCa.sex$BUGSoutput$DIC)

# iCa ~ Habitat
iCa.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$habitat)
mod.iCa.habitat=jags(model.file=model_gamma_1,data=iCa.habitat,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                   parameters.to.save=c("b0","b1","LogLik","LogLik"))
mod.iCa.habitat
exp(mod.iCa.habitat$BUGSoutput$DIC)

#iCa~(SVL*Sex)
iCa.svl.sex <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, x2 = iCa.omit$SEX)
mod.iCa.svl.sex=jags(model.file=model_gamma_2mult,data=iCa.svl.sex,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                   parameters.to.save=c("b0","b1","b2","b3","LogLik","LogLik"))
mod.iCa.svl.sex
exp(mod.iCa.svl.sex$BUGSoutput$DIC)

#iCa~(SVL*Habitat)
iCa.svl.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, x2 = iCa.omit$habitat)
mod.iCa.svl.habitat=jags(model.file=model_gamma_2mult,data=iCa.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","b2","b3","LogLik","LogLik"))
mod.iCa.svl.habitat
exp(mod.iCa.svl.habitat$BUGSoutput$DIC)

#iCa~(Sex)+(Habitat)
iCa.sex.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SEX, x2 = iCa.omit$habitat)
mod.iCa.sex.habitat=jags(model.file=model_gamma_2add,data=iCa.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","b2","LogLik","LogLik"))
mod.iCa.sex.habitat
exp(mod.iCa.sex.habitat$BUGSoutput$DIC)

#iCa~(SVL*Sex)+(SVL*Habitat)
iCa.svl.sex.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, 
                                x2 = iCa.omit$SEX, x3 = iCa.omit$habitat)
mod.iCa.svl.sex.habitat=jags(model.file=model_gamma_3,data=iCa.svl.sex.habitat,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2","b3","b4","b5","LogLik","LogLik"))
mod.iCa.svl.sex.habitat
exp(mod.iCa.svl.sex.habitat$BUGSoutput$DIC)

# iCa ~ INFECTION + SVL
iCa.infection.svl <- data.frame(y = iCa.omit$iCa, 
                              x1 = iCa.omit$infection, x2 = iCa.omit$SVL)
mod.iCa.infection.svl=jags(model.file=model_gamma_2add,data=iCa.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik","LogLik"))
mod.iCa.infection.svl
exp(mod.iCa.infection.svl$BUGSoutput$DIC)

# iCa ~ INFECTION + Sex
iCa.infection.sex <- data.frame(y = iCa.omit$iCa, 
                              x1 = iCa.omit$infection, x2 = iCa.omit$SEX)
mod.iCa.infection.sex=jags(model.file=model_gamma_2add,data=iCa.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik","LogLik"))
mod.iCa.infection.sex
exp(mod.iCa.infection.sex$BUGSoutput$DIC)

# iCa ~ INFECTION + Habitat
iCa.infection.habitat <- data.frame(y = iCa.omit$iCa, 
                                  x1 = iCa.omit$infection, x2 = iCa.omit$habitat)
mod.iCa.infection.habitat=jags(model.file=model_gamma_2add,data=iCa.infection.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("b0","b1","b2","LogLik","LogLik"))
mod.iCa.infection.habitat
exp(mod.iCa.infection.habitat$BUGSoutput$DIC)

#iCa~INFECTION + (SVL*Sex)
iCa.infection.svl.sex <- data.frame(y = iCa.omit$iCa, 
                                  x1 = iCa.omit$infection, x2 = iCa.omit$SVL,
                                  x3 = iCa.omit$SEX)
mod.iCa.infection.svl.sex=jags(model.file=model_gamma_3.2,data=iCa.infection.svl.sex,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("b0","b1","b2","b3","b4","LogLik","LogLik"))
mod.iCa.infection.svl.sex
exp(mod.iCa.infection.svl.sex$BUGSoutput$DIC)

#iCa~INFECTION + (SVL*Habitat)
iCa.infection.svl.habitat <- data.frame(y = iCa.omit$iCa, 
                                      x1 = iCa.omit$infection, x2 = iCa.omit$SVL,
                                      x3 = iCa.omit$habitat)
mod.iCa.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=iCa.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("b0","b1","b2","b3","b4","LogLik","LogLik"))
mod.iCa.infection.svl.habitat
exp(mod.iCa.infection.svl.habitat$BUGSoutput$DIC)

#iCa~INFECTION + (Sex)+(Habitat)
iCa.infection.sex.habitat <- data.frame(y = iCa.omit$iCa, 
                                      x1 = iCa.omit$infection, x2 = iCa.omit$SEX,
                                      x3 = iCa.omit$habitat)
mod.iCa.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=iCa.infection.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("b0","b1","b2","b3","LogLik","LogLik"))
mod.iCa.infection.sex.habitat
exp(mod.iCa.infection.sex.habitat$BUGSoutput$DIC)

#iCa~INFECTION + (SVL*Sex)+(SVL*Habitat)
iCa.infection.svl.sex.habitat <- data.frame(y = iCa.omit$iCa, 
                                          x1 = iCa.omit$infection, x2 = iCa.omit$SVL,
                                          x3 = iCa.omit$SEX, x4 = iCa.omit$habitat)
mod.iCa.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=iCa.infection.svl.sex.habitat,
                                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                     parameters.to.save=c("b0","b1","b2","b3",
                                                          "b4","b5","b6","LogLik","LogLik"))
mod.iCa.infection.svl.sex.habitat
exp(mod.iCa.infection.svl.sex.habitat$BUGSoutput$DIC)


####Hct####
Hct.omit <- filter(diagnosed, !is.na(Hct))
#Hct~1
Hct.null <- data.frame(y = Hct.omit$Hct, x1 = 0)
mod.Hct.null=jags(model.file=model_beta_1,data=Hct.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("b0","LogLik"))
mod.Hct.null
exp(mod.Hct.null$BUGSoutput$DIC)

#Hct ~ infection
Hct.infection <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$infection)
mod.Hct.infection=jags(model.file=model_beta_1,data=Hct.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("b0","b1","LogLik"))
mod.Hct.infection
exp(mod.Hct.infection$BUGSoutput$DIC)


# Hct ~ SVL
Hct.svl <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SVL)
mod.Hct.svl=jags(model.file=model_beta_1,data=Hct.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
               parameters.to.save=c("b0","b1","LogLik"))
mod.Hct.svl
exp(mod.Hct.svl$BUGSoutput$DIC)

# Hct ~ Sex
Hct.sex <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SEX)
mod.Hct.sex=jags(model.file=model_beta_1,data=Hct.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
               parameters.to.save=c("b0","b1","LogLik"))
mod.Hct.sex
exp(mod.Hct.sex$BUGSoutput$DIC)

# Hct ~ Habitat
Hct.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$habitat)
mod.Hct.habitat=jags(model.file=model_beta_1,data=Hct.habitat,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                   parameters.to.save=c("b0","b1","LogLik"))
mod.Hct.habitat
exp(mod.Hct.habitat$BUGSoutput$DIC)

#Hct~(SVL*Sex)
Hct.svl.sex <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SVL, x2 = Hct.omit$SEX)
mod.Hct.svl.sex=jags(model.file=model_beta_2mult,data=Hct.svl.sex,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                   parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Hct.svl.sex
exp(mod.Hct.svl.sex$BUGSoutput$DIC)

#Hct~(SVL*Habitat)
Hct.svl.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SVL, x2 = Hct.omit$habitat)
mod.Hct.svl.habitat=jags(model.file=model_beta_2mult,data=Hct.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Hct.svl.habitat
exp(mod.Hct.svl.habitat$BUGSoutput$DIC)

#Hct~(Sex)+(Habitat)
Hct.sex.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SEX, x2 = Hct.omit$habitat)
mod.Hct.sex.habitat=jags(model.file=model_beta_2add,data=Hct.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Hct.sex.habitat
exp(mod.Hct.sex.habitat$BUGSoutput$DIC)

#Hct~(SVL*Sex)+(SVL*Habitat)
Hct.svl.sex.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SVL, 
                                x2 = Hct.omit$SEX, x3 = Hct.omit$habitat)
mod.Hct.svl.sex.habitat=jags(model.file=model_beta_3,data=Hct.svl.sex.habitat,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2","b3","b4","b5","LogLik"))
mod.Hct.svl.sex.habitat
exp(mod.Hct.svl.sex.habitat$BUGSoutput$DIC)

# Hct ~ INFECTION + SVL
Hct.infection.svl <- data.frame(y = Hct.omit$Hct, 
                              x1 = Hct.omit$infection, x2 = Hct.omit$SVL)
mod.Hct.infection.svl=jags(model.file=model_beta_2add,data=Hct.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Hct.infection.svl
exp(mod.Hct.infection.svl$BUGSoutput$DIC)

# Hct ~ INFECTION + Sex
Hct.infection.sex <- data.frame(y = Hct.omit$Hct, 
                              x1 = Hct.omit$infection, x2 = Hct.omit$SEX)
mod.Hct.infection.sex=jags(model.file=model_beta_2add,data=Hct.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Hct.infection.sex
exp(mod.Hct.infection.sex$BUGSoutput$DIC)

# Hct ~ INFECTION + Habitat
Hct.infection.habitat <- data.frame(y = Hct.omit$Hct, 
                                  x1 = Hct.omit$infection, x2 = Hct.omit$habitat)
mod.Hct.infection.habitat=jags(model.file=model_beta_2add,data=Hct.infection.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Hct.infection.habitat
exp(mod.Hct.infection.habitat$BUGSoutput$DIC)

#Hct~INFECTION + (SVL*Sex)
Hct.infection.svl.sex <- data.frame(y = Hct.omit$Hct, 
                                  x1 = Hct.omit$infection, x2 = Hct.omit$SVL,
                                  x3 = Hct.omit$SEX)
mod.Hct.infection.svl.sex=jags(model.file=model_beta_3.2,data=Hct.infection.svl.sex,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.Hct.infection.svl.sex
exp(mod.Hct.infection.svl.sex$BUGSoutput$DIC)

#Hct~INFECTION + (SVL*Habitat)
Hct.infection.svl.habitat <- data.frame(y = Hct.omit$Hct, 
                                      x1 = Hct.omit$infection, x2 = Hct.omit$SVL,
                                      x3 = Hct.omit$habitat)
mod.Hct.infection.svl.habitat=jags(model.file=model_beta_3.2,data=Hct.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.Hct.infection.svl.habitat
exp(mod.Hct.infection.svl.habitat$BUGSoutput$DIC)

#Hct~INFECTION + (Sex)+(Habitat)
Hct.infection.sex.habitat <- data.frame(y = Hct.omit$Hct, 
                                      x1 = Hct.omit$infection, x2 = Hct.omit$SEX,
                                      x3 = Hct.omit$habitat)
mod.Hct.infection.sex.habitat=jags(model.file=model_beta_3.3,data=Hct.infection.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Hct.infection.sex.habitat
exp(mod.Hct.infection.sex.habitat$BUGSoutput$DIC)

#Hct~INFECTION + (SVL*Sex)+(SVL*Habitat)
Hct.infection.svl.sex.habitat <- data.frame(y = Hct.omit$Hct, 
                                          x1 = Hct.omit$infection, x2 = Hct.omit$SVL,
                                          x3 = Hct.omit$SEX, x4 = Hct.omit$habitat)
mod.Hct.infection.svl.sex.habitat=jags(model.file=model_beta_3.4,data=Hct.infection.svl.sex.habitat,
                                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                       parameters.to.save=c("b0","b1","b2","b3",
                                                          "b4","b5","b6","LogLik"))
mod.Hct.infection.svl.sex.habitat
exp(mod.Hct.infection.svl.sex.habitat$BUGSoutput$DIC)

####Glu####
Glu.omit <- filter(diagnosed, !is.na(Glu))
#Glu~1
Glu.null <- data.frame(y = Glu.omit$Glu, x1 = 1)
mod.Glu.null=jags(model.file=model_gamma_1,data=Glu.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("b0","b1","LogLik"))
mod.Glu.null
exp(mod.Glu.null$BUGSoutput$DIC)

#Glu ~ infection
Glu.infection <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$infection)
mod.Glu.infection=jags(model.file=model_gamma_1,data=Glu.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","LogLik"))
mod.Glu.infection
exp(mod.Glu.infection$BUGSoutput$DIC)

# Glu ~ SVL
Glu.svl <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL)
mod.Glu.svl=jags(model.file=model_gamma_1,data=Glu.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("b0","b1","LogLik"))
mod.Glu.svl
exp(mod.Glu.svl$BUGSoutput$DIC)

# Glu ~ Sex
Glu.sex <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SEX)
mod.Glu.sex=jags(model.file=model_gamma_1,data=Glu.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("b0","b1","LogLik"))
mod.Glu.sex
exp(mod.Glu.sex$BUGSoutput$DIC)

# Glu ~ Habitat
Glu.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$habitat)
mod.Glu.habitat=jags(model.file=model_gamma_1,data=Glu.habitat,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("b0","b1","LogLik"))
mod.Glu.habitat
exp(mod.Glu.habitat$BUGSoutput$DIC)

#Glu~(SVL*Sex)
Glu.svl.sex <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL, x2 = Glu.omit$SEX)
mod.Glu.svl.sex=jags(model.file=model_gamma_2mult,data=Glu.svl.sex,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Glu.svl.sex
exp(mod.Glu.svl.sex$BUGSoutput$DIC)

#Glu~(SVL*Habitat)
Glu.svl.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL, x2 = Glu.omit$habitat)
mod.Glu.svl.habitat=jags(model.file=model_gamma_2mult,data=Glu.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Glu.svl.habitat
exp(mod.Glu.svl.habitat$BUGSoutput$DIC)

#Glu~(Sex)+(Habitat)
Glu.sex.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SEX, x2 = Glu.omit$habitat)
mod.Glu.sex.habitat=jags(model.file=model_gamma_2add,data=Glu.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Glu.sex.habitat
exp(mod.Glu.sex.habitat$BUGSoutput$DIC)

#Glu~(SVL*Sex)+(SVL*Habitat)
Glu.svl.sex.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL, 
                                  x2 = Glu.omit$SEX, x3 = Glu.omit$habitat)
mod.Glu.svl.sex.habitat=jags(model.file=model_gamma_3,data=Glu.svl.sex.habitat,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("b0","b1","b2","b3","b4","b5","LogLik"))
mod.Glu.svl.sex.habitat
exp(mod.Glu.svl.sex.habitat$BUGSoutput$DIC)

# Glu ~ INFECTION + SVL
Glu.infection.svl <- data.frame(y = Glu.omit$Glu, 
                                x1 = Glu.omit$infection, x2 = Glu.omit$SVL)
mod.Glu.infection.svl=jags(model.file=model_gamma_2add,data=Glu.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Glu.infection.svl
exp(mod.Glu.infection.svl$BUGSoutput$DIC)

# Glu ~ INFECTION + Sex
Glu.infection.sex <- data.frame(y = Glu.omit$Glu, 
                                x1 = Glu.omit$infection, x2 = Glu.omit$SEX)
mod.Glu.infection.sex=jags(model.file=model_gamma_2add,data=Glu.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Glu.infection.sex
exp(mod.Glu.infection.sex$BUGSoutput$DIC)

# Glu ~ INFECTION + Habitat
Glu.infection.habitat <- data.frame(y = Glu.omit$Glu, 
                                    x1 = Glu.omit$infection, x2 = Glu.omit$habitat)
mod.Glu.infection.habitat=jags(model.file=model_gamma_2add,data=Glu.infection.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("b0","b1","b2","LogLik"))
mod.Glu.infection.habitat
exp(mod.Glu.infection.habitat$BUGSoutput$DIC)

#Glu~INFECTION + (SVL*Sex)
Glu.infection.svl.sex <- data.frame(y = Glu.omit$Glu, 
                                    x1 = Glu.omit$infection, x2 = Glu.omit$SVL,
                                    x3 = Glu.omit$SEX)
mod.Glu.infection.svl.sex=jags(model.file=model_gamma_3.2,data=Glu.infection.svl.sex,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.Glu.infection.svl.sex
exp(mod.Glu.infection.svl.sex$BUGSoutput$DIC)

#Glu~INFECTION + (SVL*Habitat)
Glu.infection.svl.habitat <- data.frame(y = Glu.omit$Glu, 
                                        x1 = Glu.omit$infection, x2 = Glu.omit$SVL,
                                        x3 = Glu.omit$habitat)
mod.Glu.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=Glu.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))
mod.Glu.infection.svl.habitat
exp(mod.Glu.infection.svl.habitat$BUGSoutput$DIC)

#Glu~INFECTION + (Sex)+(Habitat)
Glu.infection.sex.habitat <- data.frame(y = Glu.omit$Glu, 
                                        x1 = Glu.omit$infection, x2 = Glu.omit$SEX,
                                        x3 = Glu.omit$habitat)
mod.Glu.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=Glu.infection.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("b0","b1","b2","b3","LogLik"))
mod.Glu.infection.sex.habitat
exp(mod.Glu.infection.sex.habitat$BUGSoutput$DIC)

#Glu~INFECTION + (SVL*Sex)+(SVL*Habitat)
Glu.infection.svl.sex.habitat <- data.frame(y = Glu.omit$Glu, 
                                            x1 = Glu.omit$infection, x2 = Glu.omit$SVL,
                                            x3 = Glu.omit$SEX, x4 = Glu.omit$habitat)
mod.Glu.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=Glu.infection.svl.sex.habitat,
                                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                       parameters.to.save=c("b0","b1","b2","b3",
                                                            "b4","b5","b6","LogLik"))
mod.Glu.infection.svl.sex.habitat
exp(mod.Glu.infection.svl.sex.habitat$BUGSoutput$DIC)

####PC1####
PC1.omit <- filter(diagnosed, !is.na(PC1))

#PC1~1
PC1.null <- data.frame(y = PC1.omit$PC1, x1 = 0)
mod.PC1.null=jags(model.file=model_linear_1,data=PC1.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("alpha","LogLik"))

#PC1 ~ infection
PC1.infection <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$infection)
mod.PC1.infection=jags(model.file=model_linear_1,data=PC1.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

# PC1 ~ SVL
PC1.svl <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SVL)
mod.PC1.svl=jags(model.file=model_linear_1,data=PC1.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1","LogLik"))

# PC1 ~ Sex
PC1.sex <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SEX)
mod.PC1.sex=jags(model.file=model_linear_1,data=PC1.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1","LogLik"))

# PC1 ~ Habitat
PC1.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$habitat)
mod.PC1.habitat=jags(model.file=model_linear_1,data=PC1.habitat,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("alpha","b1","LogLik"))

#PC1~(SVL*Sex)
PC1.svl.sex <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SVL, x2 = PC1.omit$SEX)
mod.PC1.svl.sex=jags(model.file=model_linear_2mult,data=PC1.svl.sex,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

#PC1~(SVL*Habitat)
PC1.svl.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SVL, x2 = PC1.omit$habitat)
mod.PC1.svl.habitat=jags(model.file=model_linear_2mult,data=PC1.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

#PC1~(Sex)+(Habitat)
PC1.sex.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SEX, x2 = PC1.omit$habitat)
mod.PC1.sex.habitat=jags(model.file=model_linear_2add,data=PC1.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","LogLik"))

#PC1~(SVL*Sex)+(SVL*Habitat)
PC1.svl.sex.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SVL, 
                                  x2 = PC1.omit$SEX, x3 = PC1.omit$habitat)
mod.PC1.svl.sex.habitat=jags(model.file=model_linear_3,data=PC1.svl.sex.habitat,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("alpha","b1","b2","b3","b4","b5","LogLik"))
# PC1 ~ INFECTION + SVL
PC1.infection.svl <- data.frame(y = PC1.omit$PC1, 
                                x1 = PC1.omit$infection, x2 = PC1.omit$SVL)
mod.PC1.infection.svl=jags(model.file=model_linear_2add,data=PC1.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC1 ~ INFECTION + Sex
PC1.infection.sex <- data.frame(y = PC1.omit$PC1, 
                                x1 = PC1.omit$infection, x2 = PC1.omit$SEX)
mod.PC1.infection.sex=jags(model.file=model_linear_2add,data=PC1.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC1 ~ INFECTION + Habitat
PC1.infection.habitat <- data.frame(y = PC1.omit$PC1, 
                                    x1 = PC1.omit$infection, x2 = PC1.omit$habitat)
mod.PC1.infection.habitat=jags(model.file=model_linear_2add,data=PC1.infection.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","LogLik"))

#PC1~INFECTION + (SVL*Sex)
PC1.infection.svl.sex <- data.frame(y = PC1.omit$PC1, 
                                    x1 = PC1.omit$infection, x2 = PC1.omit$SVL,
                                    x3 = PC1.omit$SEX)
mod.PC1.infection.svl.sex=jags(model.file=model_linear_3.2,data=PC1.infection.svl.sex,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

#PC1~INFECTION + (SVL*Habitat)
PC1.infection.svl.habitat <- data.frame(y = PC1.omit$PC1, 
                                        x1 = PC1.omit$infection, x2 = PC1.omit$SVL,
                                        x3 = PC1.omit$habitat)
mod.PC1.infection.svl.habitat=jags(model.file=model_linear_3.2,data=PC1.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

#PC1~INFECTION + (Sex)+(Habitat)
PC1.infection.sex.habitat <- data.frame(y = PC1.omit$PC1, 
                                        x1 = PC1.omit$infection, x2 = PC1.omit$SEX,
                                        x3 = PC1.omit$habitat)
mod.PC1.infection.sex.habitat=jags(model.file=model_linear_3.3,data=PC1.infection.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
#PC1~INFECTION + (SVL*Sex)+(SVL*Habitat)
PC1.infection.svl.sex.habitat <- data.frame(y = PC1.omit$PC1, 
                                            x1 = PC1.omit$infection, x2 = PC1.omit$SVL,
                                            x3 = PC1.omit$SEX, x4 = PC1.omit$habitat)
mod.PC1.infection.svl.sex.habitat=jags(model.file=model_linear_3.4,data=PC1.infection.svl.sex.habitat,
                                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                       parameters.to.save=c("alpha","b1","b2","b3",
                                                            "b4","b5","b6","LogLik"))

####PC2####

PC2.omit <- filter(diagnosed, !is.na(PC2))
#PC2~1
PC2.null <- data.frame(y = PC2.omit$PC2, x1 = 0)
mod.PC2.null=jags(model.file=model_linear_1,data=PC2.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("alpha","LogLik"))

#PC2 ~ infection
PC2.infection <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$infection)
mod.PC2.infection=jags(model.file=model_linear_1,data=PC2.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

# PC2 ~ SVL
PC2.svl <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL)
mod.PC2.svl=jags(model.file=model_linear_1,data=PC2.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1","LogLik"))

# PC2 ~ Sex
PC2.sex <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SEX)
mod.PC2.sex=jags(model.file=model_linear_1,data=PC2.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1","LogLik"))

# PC2 ~ Habitat
PC2.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$habitat)
mod.PC2.habitat=jags(model.file=model_linear_1,data=PC2.habitat,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("alpha","b1","LogLik"))

#PC2~(SVL*Sex)
PC2.svl.sex <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL, x2 = PC2.omit$SEX)
mod.PC2.svl.sex=jags(model.file=model_linear_2mult,data=PC2.svl.sex,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

#PC2~(SVL*Habitat)
PC2.svl.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL, x2 = PC2.omit$habitat)
mod.PC2.svl.habitat=jags(model.file=model_linear_2mult,data=PC2.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

#PC2~(Sex)+(Habitat)
PC2.sex.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SEX, x2 = PC2.omit$habitat)
mod.PC2.sex.habitat=jags(model.file=model_linear_2add,data=PC2.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","LogLik"))

#PC2~(SVL*Sex)+(SVL*Habitat)
PC2.svl.sex.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL, 
                                  x2 = PC2.omit$SEX, x3 = PC2.omit$habitat)
mod.PC2.svl.sex.habitat=jags(model.file=model_linear_3,data=PC2.svl.sex.habitat,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("alpha","b1","b2","b3","b4","b5","LogLik"))

# PC2 ~ INFECTION + SVL
PC2.infection.svl <- data.frame(y = PC2.omit$PC2, 
                                x1 = PC2.omit$infection, x2 = PC2.omit$SVL)
mod.PC2.infection.svl=jags(model.file=model_linear_2add,data=PC2.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC2 ~ INFECTION + Sex
PC2.infection.sex <- data.frame(y = PC2.omit$PC2, 
                                x1 = PC2.omit$infection, x2 = PC2.omit$SEX)
mod.PC2.infection.sex=jags(model.file=model_linear_2add,data=PC2.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC2 ~ INFECTION + Habitat
PC2.infection.habitat <- data.frame(y = PC2.omit$PC2, 
                                    x1 = PC2.omit$infection, x2 = PC2.omit$habitat)
mod.PC2.infection.habitat=jags(model.file=model_linear_2add,data=PC2.infection.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","LogLik"))

#PC2~INFECTION + (SVL*Sex)
PC2.infection.svl.sex <- data.frame(y = PC2.omit$PC2, 
                                    x1 = PC2.omit$infection, x2 = PC2.omit$SVL,
                                    x3 = PC2.omit$SEX)
mod.PC2.infection.svl.sex=jags(model.file=model_linear_3.2,data=PC2.infection.svl.sex,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

#PC2~INFECTION + (SVL*Habitat)
PC2.infection.svl.habitat <- data.frame(y = PC2.omit$PC2, 
                                        x1 = PC2.omit$infection, x2 = PC2.omit$SVL,
                                        x3 = PC2.omit$habitat)
mod.PC2.infection.svl.habitat=jags(model.file=model_linear_3.2,data=PC2.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

#PC2~INFECTION + (Sex)+(Habitat)
PC2.infection.sex.habitat <- data.frame(y = PC2.omit$PC2, 
                                        x1 = PC2.omit$infection, x2 = PC2.omit$SEX,
                                        x3 = PC2.omit$habitat)
mod.PC2.infection.sex.habitat=jags(model.file=model_linear_3.3,data=PC2.infection.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

#PC2~INFECTION + (SVL*Sex)+(SVL*Habitat)
PC2.infection.svl.sex.habitat <- data.frame(y = PC2.omit$PC2, 
                                            x1 = PC2.omit$infection, x2 = PC2.omit$SVL,
                                            x3 = PC2.omit$SEX, x4 = PC2.omit$habitat)
mod.PC2.infection.svl.sex.habitat=jags(model.file=model_linear_3.4,data=PC2.infection.svl.sex.habitat,
                                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                       parameters.to.save=c("alpha","b1","b2","b3",
                                                            "b4","b5","b6","LogLik"))



