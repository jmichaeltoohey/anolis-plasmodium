source('data_carpentry.R')
library("R2jags")
library(BayesPostEst)
devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")
#library('rjags')

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

#bayesian parameterization of beta regression w/ 1 predictor
model_beta_1 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
}
#bayesian parameterization of beta regression w/ 2 predictors (additive)
model_beta_2add <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    mu[i] <- 1/(1+exp(-(b0 + b1*x1[i] + b2*x2[i])))
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  mu ~ dunif(min=15,max=75)
  phi ~ dgamma(0.1,0.1)
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
}

#bayesian parameterization of beta regression w/ 2 predictors (interaction)
model_beta_2mult <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    mu[i] <- 1/(1+exp(-(b0 + b1*x1[i] + b2*x2[i] + b3*x1[i]*x2[i])))
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  mu ~ dunif(min=15,max=75)
  phi ~ dgamma(0.1,0.1)
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
}

#bayesian parameterization of beta regression w/ 3 predictors
model_beta_3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    mu[i] <- 1/(1+exp(-(b0 + b1*x1[i] + b2*x2[i] + b3*x1[i]*x2[i] + b4*x1[i]*x2[i] + b5*x1[i]*x3[i])))
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  mu ~ dunif(min=15,max=75)
  phi ~ dgamma(0.1,0.1)
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
}

#bayesian parameterization of beta regression w/ 3 predictors (interaction)
model_beta_3.2 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    mu[i] <- 1/(1+exp(-(b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x2[i]*x3[i])))
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  mu ~ dunif(.15, .75)
  phi ~ dgamma(0.1,0.1)
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
}
#bayesian parameterization of beta regression w/ 3 predictors (interaction)
model_beta_3.2 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x2[i]*x3[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  mu ~ dunif(.15, .75)
  phi ~ dgamma(0.1,0.1)
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
}

#bayesian parameterization of beta regression w/ 3 predictors all additive
model_beta_3.3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    mu[i] <- 1/(1+exp(-(b0 + b1*x1[i] + b2*x2[i] + b3*x3[i])))
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  mu ~ dunif(min=15,max=75)
  phi ~ dgamma(0.1,0.1)
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
}

#bayesian parameterization of beta regression w/ 4 predictors
model_beta_3.4 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    mu[i] <- 1/(1+exp(-(b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i] + 
                          b5*x2[i]*x3[i] + b6*x2[i]*x4[i])))
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  mu ~ dunif(min=15,max=75)
  phi ~ dgamma(0.1,0.1)
  b0 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
  b6 ~ dnorm(0, 0.0001)
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

ni <- 100000
nt <- 100
nb <- 5000
nc <- 3



####p(infection)~####

#p(infection)~(svl.sex)+(SVL*Habitat)
infection.svl.sex.habitat <- data.frame(y = diagnosed$infection, n = nrow(diagnosed), 
                                        x1 = diagnosed$SVL, x2 = diagnosed$SEX, x3 = diagnosed$habitat)
mod.infection.svl.sex.habitat=jags(model.file=model_binom_3,data=infection.svl.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","b4","b5"))
mcmc.infection.svl.sex.habitat <- as.mcmc(mod.infection.svl.sex.habitat)
mat.infection.svl.sex.habitat <- as.matrix(mcmc.infection.svl.sex.habitat)
dat.infection.svl.sex.habitat <- as.data.frame(mat.infection.svl.sex.habitat)
betas.infection.svl.sex.habitat <- dat.infection.svl.sex.habitat[, grep(x = colnames(dat.infection.svl.sex.habitat), 
                                          pattern = "b", 
                                          fixed = TRUE)]
devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")
mcmctab(betas.infection.svl.sex.habitat)

write.table(mcmcTab(mod.infection.svl.sex.habitat, pars=c("alpha","b1","b2","b3","b4","b5")),file="mod.infection.svl.sex.habitat.txt", sep=",", quote=FALSE,row.names=F)

gelman.diag(as.mcmc(mod.infection.svl.sex.habitat))
traceplot(mod.infection.svl.sex.habitat)

#p(infection)~(SVL*Habitat)
nt <- 100
infection.svl.habitat <- data.frame(y = diagnosed$infection, 
                                    #n = nrow(diagnosed), 
                                    x1 = diagnosed$SVL, x2 = diagnosed$habitat)
mod.infection.svl.habitat=jags(model.file=model_binom_2mult,data=infection.svl.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","b3"))
traceplot(mod.infection.svl.habitat)
autocorr.plot(mod.infection.svl.habitat)
autocorr.diag(mod.infection.svl.habitat)
gelman.diag(as.mcmc(mod.infection.svl.habitat))
write.table(mcmcTab(mod.infection.svl.habitat, pars=c("alpha","b1","b2","b3")),file="mod.infection.svl.habitat.txt", sep=",", quote=FALSE,row.names=F)


####body condition####

devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")


#body condition ~ infection
male.bodycond.infection <- data.frame(y = diagnosed.males$sex.bodycond, x1 = diagnosed.males$infection)
mod.male.bodycond.infection=jags(model.file=model_linear_1,data=male.bodycond.infection,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("alpha","b1"))
write.table(mcmctab(mod.male.bodycond.infection, pars=c("alpha","b1")),file="mod.male.bodycond.infection.txt", sep=",", quote=FALSE,row.names=F)


autocorr.plot(mod.male.bodycond.infection)
gelman.diag(as.mcmc(mod.male.bodycond.infection))

# body condition ~ INFECTION + Habitat
male.bodycond.infection.habitat <- data.frame(y = diagnosed.males$sex.bodycond,
                                              x1 = diagnosed.males$infection, x2 = diagnosed.males$habitat)
mod.male.bodycond.infection.habitat=jags(model.file=model_linear_2add,data=male.bodycond.infection.habitat,
                                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                         parameters.to.save=c("alpha","b1","b2"))
write.table(mcmctab(mod.male.bodycond.infection.habitat, pars=c("alpha","b1","b2")),file="mod.male.bodycond.infection.habitat.txt", sep=",", quote=FALSE,row.names=F)

autocorr.plot(mod.male.bodycond.infection.habitat)
gelman.diag(as.mcmc(mod.male.bodycond.infection.habitat))

#body condition~1
male.bodycond.null <- data.frame(y = diagnosed.males$sex.bodycond, x1 = 0)
mod.male.bodycond.null=jags(model.file=model_linear_1,data=male.bodycond.null,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha"))
write.table(mcmctab(mod.male.bodycond.null, pars=c("alpha")),file="mod.male.bodycond.null.txt", sep=",", quote=FALSE,row.names=F)

autocorr.plot(mod.male.bodycond.null)
gelman.diag(as.mcmc(mod.male.bodycond.null))

# body condition ~ Habitat
male.bodycond.habitat <- data.frame(y = diagnosed.males$sex.bodycond, x1 = diagnosed.males$habitat)
mod.male.bodycond.habitat=jags(model.file=model_linear_1,data=male.bodycond.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1"))
write.table(mcmctab(mod.male.bodycond.habitat, pars=c("alpha","b1")),file="mod.male.bodycond.habitat.txt", sep=",", quote=FALSE,row.names=F)

autocorr.plot(mod.male.bodycond.habitat)
gelman.diag(as.mcmc(mod.male.bodycond.habitat))

####Hct####
Hct.omit <- filter(diagnosed, !is.na(Hct))
#Hct~INFECTION + (svl.sex)
Hct.infection.svl.sex <- data.frame(y = Hct.omit$Hct, 
                                    x1 = Hct.omit$infection, x2 = Hct.omit$SVL,
                                    x3 = Hct.omit$SEX)
mod.Hct.infection.svl.sex=jags(model.file=model_beta_3.2,data=Hct.infection.svl.sex,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("b0","b1","b2","b3","b4"))
write.table(mcmcTab(mod.Hct.infection.svl.sex, pars=c("b0","b1","b2","b3","b4")),file="mod.Hct.infection.svl.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Hct.infection.svl.sex)
gelman.diag(as.mcmc(mod.Hct.infection.svl.sex))

#Hct~1
Hct.null <- data.frame(y = Hct.omit$Hct, x1 = 0)
mod.Hct.null=jags(model.file=model_beta_1,data=Hct.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("b0"))
write.table(mcmcTab(mod.Hct.null, pars=c("b0")),file="mod.Hct.null.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Hct.null)
gelman.diag(as.mcmc(mod.Hct.null))

#Hct ~ infection
Hct.infection <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$infection)

mod.Hct.infection=jags(model.file=model_beta_1,data=Hct.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1"))
write.table(mcmcTab(mod.Hct.infection, pars=c("b0","b1")),file="mod.Hct.infection.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Hct.infection)
gelman.diag(as.mcmc(mod.Hct.infection))
traceplot(mod.Hct.infection)

#Hct~(svl.sex)
Hct.svl.sex <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SVL, x2 = Hct.omit$SEX)
mod.Hct.svl.sex=jags(model.file=model_beta_2mult,data=Hct.svl.sex,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("b0","b1","b2","b3"))
write.table(mcmcTab(mod.Hct.svl.sex, pars=c("b0","b1","b2","b3")),file="mod.Hct.svl.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Hct.svl.sex)
gelman.diag(as.mcmc(mod.Hct.svl.sex))
traceplot(mod.Hct.svl.sex)

# Hct ~ Sex
Hct.sex <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SEX)
mod.Hct.sex=jags(model.file=model_beta_1,data=Hct.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("b0","b1"))
write.table(mcmcTab(mod.Hct.sex, pars=c("b0","b1")),file="mod.Hct.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Hct.sex)
gelman.diag(as.mcmc(mod.Hct.sex))

# Hct ~ INFECTION + Sex
Hct.infection.sex <- data.frame(y = Hct.omit$Hct, 
                                x1 = Hct.omit$infection, x2 = Hct.omit$SEX)
mod.Hct.infection.sex=jags(model.file=model_beta_2add,data=Hct.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2"))
write.table(mcmcTab(mod.Hct.infection.sex, pars=c("b0","b1","b2")),file="mod.Hct.infection.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Hct.infection.sex)
gelman.diag(as.mcmc(mod.Hct.infection.sex))

# Hct ~ SVL
Hct.svl <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$svl)
mod.Hct.svl=jags(model.file=model_beta_1,data=Hct.svl,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("b0","b1"))
write.table(mcmcTab(mod.Hct.svl, pars=c("b0","b1")),file="mod.Hct.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Hct.svl)
gelman.diag(as.mcmc(mod.Hct.svl))

####K####
K.omit <- filter(diagnosed, !is.na(K))
# K ~ SVL
K.svl <- data.frame(y = K.omit$K, x1 = K.omit$SVL)
mod.K.svl=jags(model.file=model_gamma_1,data=K.svl,
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
               parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.K.svl, pars=c("alpha","b1")),file="mod.K.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.K.svl)
gelman.diag(as.mcmc(mod.K.svl))

# K ~ INFECTION + SVL
K.infection.svl <- data.frame(y = K.omit$K, 
                              x1 = K.omit$infection, x2 = K.omit$SVL)
mod.K.infection.svl=jags(model.file=model_gamma_2add,data=K.infection.svl,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2"))
write.table(mcmcTab(mod.K.infection.svl, pars=c("alpha","b1","b2")),file="mod.K.infection.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.K.infection.svl)
gelman.diag(as.mcmc(mod.K.infection.svl))
# K ~ Sex
K.sex <- data.frame(y = K.omit$K, x1 = K.omit$SEX)
mod.K.sex=jags(model.file=model_gamma_1,data=K.sex,
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
               parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.K.sex, pars=c("alpha","b1")),file="mod.K.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.K.sex)
gelman.diag(as.mcmc(mod.K.sex))

# K ~ INFECTION + Sex
K.infection.sex <- data.frame(y = K.omit$K, 
                              x1 = K.omit$infection, x2 = K.omit$SEX)
mod.K.infection.sex=jags(model.file=model_gamma_2add,data=K.infection.sex,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2"))
write.table(mcmcTab(mod.K.infection.sex, pars=c("alpha","b1","b2")),file="mod.K.infection.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.K.infection.sex)
gelman.diag(as.mcmc(mod.K.infection.sex))
traceplot(mod.K.infection.sex)


####Na####
Sodium.omit <- filter(diagnosed, !is.na(Na))
#Na~1
Na.null <- data.frame(y = Sodium.omit$Na, x1 = 0)
mod.Na.null=jags(model.file=model_gamma_1,data=Na.null,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha"))
write.table(mcmcTab(mod.Na.null, pars=c("alpha")),file="mod.Na.null.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Na.null)
gelman.diag(as.mcmc(mod.Na.null))
# Na ~ SVL
Na.svl <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL)
mod.Na.svl=jags(model.file=model_gamma_1,data=Na.svl,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Na.svl, pars=c("alpha","b1")),file="mod.Na.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Na.svl)
gelman.diag(as.mcmc(mod.Na.svl))
# Na ~ Sex
Na.sex <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SEX)
mod.Na.sex=jags(model.file=model_gamma_1,data=Na.sex,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Na.sex, pars=c("alpha","b1")),file="mod.Na.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Na.sex)
gelman.diag(as.mcmc(mod.Na.sex))
# Na ~ Habitat
Na.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$habitat)
mod.Na.habitat=jags(model.file=model_gamma_1,data=Na.habitat,
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                    parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Na.habitat, pars=c("alpha","b1")),file="mod.Na.habitat.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Na.habitat)
gelman.diag(as.mcmc(mod.Na.habitat))
#Na ~ infection
Na.infection <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$infection)
mod.Na.infection=jags(model.file=model_gamma_1,data=Na.infection,
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Na.infection, pars=c("alpha","b1")),file="mod.Na.infection.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Na.infection)
gelman.diag(as.mcmc(mod.Na.infection))

####Cl####
Cl.omit <- filter(diagnosed, !is.na(Cl))

#Cl~1
Cl.null <- data.frame(y = Cl.omit$Cl, x1 = 0)
mod.Cl.null=jags(model.file=model_gamma_1,data=Cl.null,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha"))
write.table(mcmcTab(mod.Cl.null, pars=c("alpha")),file="mod.Cl.null.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Cl.null)
gelman.diag(as.mcmc(mod.Cl.null))

# Cl ~ SVL
Cl.svl <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL)
mod.Cl.svl=jags(model.file=model_gamma_1,data=Cl.svl,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Cl.svl, pars=c("alpha","b1")),file="mod.Cl.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Cl.svl)
gelman.diag(as.mcmc(mod.Cl.svl))

# Cl ~ Habitat
Cl.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$habitat)
mod.Cl.habitat=jags(model.file=model_gamma_1,data=Cl.habitat,
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                    parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Cl.habitat, pars=c("alpha","b1")),file="mod.Cl.habitat.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Cl.habitat)
gelman.diag(as.mcmc(mod.Cl.habitat))

#Cl ~ infection
Cl.infection <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$infection)
mod.Cl.infection=jags(model.file=model_gamma_1,data=Cl.infection,
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Cl.infection, pars=c("alpha","b1")),file="mod.Cl.infection.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Cl.infection)
gelman.diag(as.mcmc(mod.Cl.infection))

# Cl ~ Sex
Cl.sex <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SEX)
mod.Cl.sex=jags(model.file=model_gamma_1,data=Cl.sex,
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Cl.sex, pars=c("alpha","b1")),file="mod.Cl.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Cl.sex)
gelman.diag(as.mcmc(mod.Cl.sex))

# Cl ~ INFECTION + SVL
Cl.infection.svl <- data.frame(y = Cl.omit$Cl, 
                               x1 = Cl.omit$infection, x2 = Cl.omit$SVL)
mod.Cl.infection.svl=jags(model.file=model_gamma_2add,data=Cl.infection.svl,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("alpha","b1","b2"))
write.table(mcmcTab(mod.Cl.infection.svl, pars=c("alpha","b1","b2")),file="mod.Cl.infection.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Cl.infection.svl)
gelman.diag(as.mcmc(mod.Cl.infection.svl))

#Cl~(Sex)+(Habitat)
Cl.sex.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SEX, x2 = Cl.omit$habitat)
mod.Cl.sex.habitat=jags(model.file=model_gamma_2add,data=Cl.sex.habitat,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("alpha","b1","b2"))
write.table(mcmcTab(mod.Cl.sex.habitat, pars=c("alpha","b1","b2")),file="mod.Cl.sex.habitat.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Cl.sex.habitat)
gelman.diag(as.mcmc(mod.Cl.sex.habitat))

# Cl ~ INFECTION + Habitat
Cl.infection.habitat <- data.frame(y = Cl.omit$Cl, 
                                   x1 = Cl.omit$infection, x2 = Cl.omit$habitat)
mod.Cl.infection.habitat=jags(model.file=model_gamma_2add,data=Cl.infection.habitat,
                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                              parameters.to.save=c("alpha","b1","b2"))
write.table(mcmcTab(mod.Cl.infection.habitat, pars=c("alpha","b1","b2")),file="mod.Cl.infection.habitat.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Cl.infection.habitat)
gelman.diag(as.mcmc(mod.Cl.infection.habitat))

####iCa####
iCa.omit <- filter(diagnosed, !is.na(iCa))

# iCa ~ INFECTION + SVL
iCa.infection.svl <- data.frame(y = iCa.omit$iCa, 
                                x1 = iCa.omit$infection, x2 = iCa.omit$SVL)
mod.iCa.infection.svl=jags(model.file=model_gamma_2add,data=iCa.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2"))
write.table(mcmcTab(mod.iCa.infection.svl, pars=c("alpha","b1","b2")),file="mod.iCa.infection.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.iCa.infection.svl)
gelman.diag(as.mcmc(mod.iCa.infection.svl))

# iCa ~ SVL
iCa.svl <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL)
mod.iCa.svl=jags(model.file=model_gamma_1,data=iCa.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.iCa.svl, pars=c("alpha","b1")),file="mod.iCa.svl.txt", sep=",", quote=FALSE,row.names=F)
#autocorr.plot(mod.iCa.svl)
#gelman.diag(as.mcmc(mod.iCa.svl))

#iCa~(SVL*Habitat)
iCa.svl.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, x2 = iCa.omit$habitat)
mod.iCa.svl.habitat=jags(model.file=model_gamma_2mult,data=iCa.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","b3"))
write.table(mcmcTab(mod.iCa.svl.habitat, pars=c("alpha","b1","b2","b3")),file="mod.iCa.svl.habitat.txt", sep=",", quote=FALSE,row.names=F)
#autocorr.plot(mod.iCa.svl.habitat)
#gelman.diag(as.mcmc(mod.iCa.svl.habitat))

#iCa~(svl.sex)
iCa.svl.sex <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, x2 = iCa.omit$SEX)
mod.iCa.svl.sex=jags(model.file=model_gamma_2mult,data=iCa.svl.sex,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("alpha","b1","b2","b3"))
write.table(mcmcTab(mod.iCa.svl.sex, pars=c("alpha","b1","b2","b3")),file="mod.iCa.svl.sex.txt", sep=",", quote=FALSE,row.names=F)
#autocorr.plot(mod.iCa.svl.sex)
#gelman.diag(as.mcmc(mod.iCa.svl.sex))

#iCa~INFECTION + (svl.sex)
iCa.infection.svl.sex <- data.frame(y = iCa.omit$iCa, 
                                    x1 = iCa.omit$infection, x2 = iCa.omit$SVL,
                                    x3 = iCa.omit$SEX)
mod.iCa.infection.svl.sex=jags(model.file=model_gamma_3.2,data=iCa.infection.svl.sex,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","b3","b4"))
write.table(mcmcTab(mod.iCa.infection.svl.sex, pars=c("alpha","b1","b2","b3","b4")),file="mod.iCa.infection.svl.sex.txt", sep=",", quote=FALSE,row.names=F)
#autocorr.plot(mod.iCa.infection.svl.sex)
#gelman.diag(as.mcmc(mod.iCa.infection.svl.sex))

#iCa~(svl.sex)+(SVL*Habitat)
iCa.svl.sex.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, 
                                  x2 = iCa.omit$SEX, x3 = iCa.omit$habitat)
mod.iCa.svl.sex.habitat=jags(model.file=model_gamma_3,data=iCa.svl.sex.habitat,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("alpha","b1","b2","b3","b4","b5"))
write.table(mcmcTab(mod.iCa.svl.sex.habitat, pars=c("alpha","b1","b2","b3","b4","b5")),file="mod.iCa.svl.sex.habitat.txt", sep=",", quote=FALSE,row.names=F)
#autocorr.plot(mod.iCa.svl.sex.habitat)
#gelman.diag(as.mcmc(mod.iCa.svl.sex.habitat))

#iCa~INFECTION + (SVL*Habitat)
iCa.infection.svl.habitat <- data.frame(y = iCa.omit$iCa, 
                                        x1 = iCa.omit$infection, x2 = iCa.omit$SVL,
                                        x3 = iCa.omit$habitat)
mod.iCa.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=iCa.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","b4"))
write.table(mcmcTab(mod.iCa.infection.svl.habitat, pars=c("alpha","b1","b2","b3","b4")),file="mod.iCa.infection.svl.habitat.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.iCa.infection.svl.habitat)
gelman.diag(as.mcmc(mod.iCa.infection.svl.habitat))


####Glu####
Glu.omit <- filter(diagnosed, !is.na(Glu))

# Glu ~ Sex
Glu.sex <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SEX)
mod.Glu.sex=jags(model.file=model_gamma_1,data=Glu.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Glu.sex, pars=c("alpha","b1")),file="mod.Glu.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Glu.sex)
gelman.diag(as.mcmc(mod.Glu.sex))

# Glu ~ SVL
Glu.svl <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL)
mod.Glu.svl=jags(model.file=model_gamma_1,data=Glu.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.Glu.svl, pars=c("alpha","b1")),file="mod.Glu.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Glu.svl)
gelman.diag(as.mcmc(mod.Glu.svl))

#Glu~(svl.sex)
Glu.svl.sex <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL, x2 = Glu.omit$SEX)
mod.Glu.svl.sex=jags(model.file=model_gamma_2mult,data=Glu.svl.sex,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","b3"))
write.table(mcmcTab(mod.Glu.svl.sex, pars=c("alpha","b1","b2","b3")),file="mod.Glu.svl.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Glu.svl.sex)
gelman.diag(as.mcmc(mod.Glu.svl.sex))

#Glu~(Sex)+(Habitat)
Glu.sex.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SEX, x2 = Glu.omit$habitat)
mod.Glu.sex.habitat=jags(model.file=model_gamma_2add,data=Glu.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2"))
write.table(mcmcTab(mod.Glu.sex.habitat, pars=c("alpha","b1","b2")),file="mod.Glu.sex.habitat.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Glu.sex.habitat)
gelman.diag(as.mcmc(mod.Glu.sex.habitat))

# Glu ~ INFECTION + Sex
Glu.infection.sex <- data.frame(y = Glu.omit$Glu, 
                                x1 = Glu.omit$infection, x2 = Glu.omit$SEX)
mod.Glu.infection.sex=jags(model.file=model_gamma_2add,data=Glu.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2"))
write.table(mcmcTab(mod.Glu.infection.sex, pars=c("alpha","b1","b2")),file="mod.Glu.infection.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.Glu.infection.sex)
gelman.diag(as.mcmc(mod.Glu.infection.sex))

####PC1####
PC1.omit <- filter(diagnosed, !is.na(PC1))

#PC1~1
PC1.null <- data.frame(y = PC1.omit$PC1, x1 = 0)
mod.PC1.null=jags(model.file=model_linear_1,data=PC1.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("alpha"))
write.table(mcmcTab(mod.PC1.null, pars=c("alpha")),file="mod.PC1.null.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.PC1.null)
gelman.diag(as.mcmc(mod.PC1.null))

#PC1 ~ infection
PC1.infection <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$infection)
mod.PC1.infection=jags(model.file=model_linear_1,data=PC1.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.PC1.infection, pars=c("alpha","b1")),file="mod.PC1.infection.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.PC1.infection)
gelman.diag(as.mcmc(mod.PC1.infection))

# PC1 ~ SVL
PC1.svl <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SVL)
mod.PC1.svl=jags(model.file=model_linear_1,data=PC1.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1"))
write.table(mcmcTab(mod.PC1.svl, pars=c("alpha","b1")),file="mod.PC1.svl.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.PC1.svl)
gelman.diag(as.mcmc(mod.PC1.svl))

####PC2####

PC2.omit <- filter(diagnosed, !is.na(PC2))

#PC2~(svl.sex)+(SVL*Habitat)
PC2.svl.sex.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL, 
                                  x2 = PC2.omit$SEX, x3 = PC2.omit$habitat)
mod.PC2.svl.sex.habitat=jags(model.file=model_linear_3,data=PC2.svl.sex.habitat,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("alpha","b1","b2","b3","b4","b5"))
write.table(mcmcTab(mod.PC2.svl.sex.habitat, pars=c("alpha","b1","b2","b3","b4","b5")),file="mod.PC2.svl.sex.habitat.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.PC2.svl.sex.habitat)
gelman.diag(as.mcmc(mod.PC2.svl.sex.habitat))
#PC2~(svl.sex)
PC2.svl.sex <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL, x2 = PC2.omit$SEX)
mod.PC2.svl.sex=jags(model.file=model_linear_2mult,data=PC2.svl.sex,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("alpha","b1","b2","b3"))
write.table(mcmcTab(mod.PC2.svl.sex, pars=c("alpha","b1","b2","b3")),file="mod.PC2.svl.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.PC2.svl.sex)
gelman.diag(as.mcmc(mod.PC2.svl.sex))

#PC2~INFECTION + (svl.sex)
PC2.infection.svl.sex <- data.frame(y = PC2.omit$PC2, 
                                    x1 = PC2.omit$infection, x2 = PC2.omit$SVL,
                                    x3 = PC2.omit$SEX)
mod.PC2.infection.svl.sex=jags(model.file=model_linear_3.2,data=PC2.infection.svl.sex,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","b3","b4"))
write.table(mcmcTab(mod.PC2.infection.svl.sex, pars=c("alpha","b1","b2","b3","b4")),file="mod.PC2.infection.svl.sex.txt", sep=",", quote=FALSE,row.names=F)
autocorr.plot(mod.PC2.infection.svl.sex)
gelman.diag(as.mcmc(mod.PC2.infection.svl.sex))
