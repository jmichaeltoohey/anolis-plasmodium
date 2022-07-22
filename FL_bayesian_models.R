source('FL_data_carpentry.R')
library("R2jags")
library(statip)

####BAYESIAN PARAMETERIZATIONS####
ni <- 100000
nt <- 100
nb <- 5000
nc <- 3
set.seed(123)
####--------------------------------BINOMIAL--------------------------------####

#---------------------------------------------------------------------NULL MODEL
model_binom_0 <- function() {                                        
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- mu[i]
    mu[i] <- alpha
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
}

#--------------------------------------------------------MODELS WITH 1 COVARIATE
model_binom_1 <- function() {                           
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- mu[i]
    mu[i] <- alpha + b1*x1[i]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
#-----------------------------------------MODELS WITH SITE AS THE ONLY COVARIATE
model_binom_1.site <- function() {       
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- mu[i]
    mu[i] <- bSite[x1[i]]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  # cat distrib
  #bSite[1] <- 0
  for(j in 1:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  #alpha ~ dnorm(0, 0.0001)
}

#-------------------------------------------INTERACTIVE MODELS WITH 2 COVARIATES
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
#----------------------------INTERACTIVE MODELS WITH 2 COVARIATES (1 BEING SITE)
model_binom_2mult.site <- function() {                                          
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + bSite[x2[i]] + bSite.Int[x2[i]]*x1[i]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int[1] <- 0
  for(j in 2:3){
    bSite.Int[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}
#----------------------------------------------ADDITIVE MODELS WITH 2 COVARIATES
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
#-------------------------------ADDITIVE MODELS WITH 2 COVARIATES (1 BEING SITE)
model_binom_2add.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + bSite[x2[i]]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
}

#-------------------------------------------INTERACTIVE MODELS WITH 3 COVARIATES
model_binom_3int <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x1[i]*x2[i] +
      b5*x1[i]*x3[i] + b6*x2[i]*x3[i]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
  b6 ~ dnorm(0, 0.0001)
}
#----------------------------INTERACTIVE MODELS WITH 3 COVARIATES (1 BEING SITE)
model_binom_3int.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    #y[i] ~ dbin(p[i], 1)
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + b2*x2[i] + bSite[x3[i]] + b4*x1[i]*x2[i] +
      bSite.Int1[x3[i]]*x1[i] + bSite.Int2[x3[i]]*x2[i]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int1[1] <- 0
  for(j in 2:3){
    bSite.Int1[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int2[1] <- 0
  for(j in 2:3){
    bSite.Int2[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
}
#----------------------------------------------ADDITIVE MODELS WITH 3 COVARIATES
model_binom_3add <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
}
#-------------------------------ADDITIVE MODELS WITH 3 COVARIATES (1 BEING SITE)
model_binom_3add.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbin(p[i], 1)
    logit(p[i]) <- alpha + b1*x1[i] + b2*x2[i] + bSite[x3[i]]
    LogLik[i] <- log(dbin(y[i], p[i], 1))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
}

####---------------------------------LINEAR---------------------------------####

#---------------------------------------------------------------------NULL MODEL
model_linear_0 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#--------------------------------------------------------MODELS WITH 1 COVARIATE
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
#-----------------------------------------MODELS WITH SITE AS THE ONLY COVARIATE
model_linear_1.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- bSite[x1[i]]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  # cat distrib
  #bSite[1] <- 0
  for(j in 1:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  #alpha ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#-------------------------------------------INTERACTIVE MODELS WITH 2 COVARIATES
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
#----------------------------INTERACTIVE MODELS WITH 2 COVARIATES (1 BEING SITE)
model_linear_2mult.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + bSite[x2[i]] + bSite.Int[x2[i]]*x1[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int[1] <- 0
  for(j in 2:3){
    bSite.Int[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}
#----------------------------------------------ADDITIVE MODELS WITH 2 COVARIATES
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
#-------------------------------ADDITIVE MODELS WITH 2 COVARIATES (1 BEING SITE)
model_linear_2add.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + bSite[x2[i]]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#-------------------------------------------INTERACTIVE MODELS WITH 3 COVARIATES
model_linear_3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x1[i]*x2[i] + b5*x1[i]*x3[i] + b6*x2[i]*x3[i]
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
#----------------------------INTERACTIVE MODELS WITH 3 COVARIATES (1 BEING SITE)
model_linear_3.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + bSite[x3[i]] + b4*x1[i]*x2[i] + 
      bSite.Int1[x3[i]]*x1[i] + bSite.Int2[x3[i]]*x2[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int1[1] <- 0
  for(j in 2:3){
    bSite.Int1[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int2[1] <- 0
  for(j in 2:3){
    bSite.Int2[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}
#--------------------------------MODELS WITH 3 COVARIATES AND ONLY 1 INTERACTION
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
#--------------MODELS WITH 3 COVARIATES AND ONLY 1 INTERACTION (BEING WITH SITE)
model_linear_3.2.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + bSite[x3[i]] + bSite.Int[x3[i]]*x2[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int[1] <- 0
  for(j in 2:3){
    bSite.Int[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}
#----------------------------------------------ADDITIVE MODELS WITH 3 COVARIATES
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
#-------------------------------ADDITIVE MODELS WITH 3 COVARIATES (1 BEING SITE)
model_linear_3.3.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + bSite[x3[i]]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

#-------------------------------ADDITIVE MODELS WITH 4 COVARIATES (1 BEING SITE)
model_linear_4 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + bSite[x4[i]]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}
#----------------------------INTERACTIVE MODELS WITH 4 COVARIATES (1 BEING SITE)
model_linear_4.2 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + b1*x1[i] + b2*x2[i] + b3*x3[i] + bSite[x4[i]] + 
      b5*x2[i]*x3[i] + bSite.Int1[x4[i]]*x2[i] + bSite.Int2[x4[i]]*x3[i]
    LogLik[i] <- log(dnorm(y[i],mu[i], tau))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int1[1] <- 0
  for(j in 2:3){
    bSite.Int1[j] ~ dnorm(0, 0.0001)
  }
  # interaction distrib
  bSite.Int2[1] <- 0
  for(j in 2:3){
    bSite.Int2[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}


####----------------------------------BETA----------------------------------####

#---------------------------------------------------------------------NULL MODEL
model_beta_0 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
}
#--------------------------------------------------------MODELS WITH 1 COVARIATE
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
#-----------------------------------------MODELS WITH SITE AS THE ONLY COVARIATE
model_beta_1.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- bSite[x1[i]]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  # cat distrib
  bSite[1] ~ dunif(log(.15), log(.75))
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
}
#----------------------------------------------ADDITIVE MODELS WITH 2 COVARIATES
model_beta_2add <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
}
#-------------------------------ADDITIVE MODELS WITH 2 COVARIATES (1 BEING SITE)
model_beta_2add.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + bSite[x2[i]]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
}
#-------------------------------------------INTERACTIVE MODELS WITH 2 COVARIATES
model_beta_2mult <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + b3*x1[i]*x2[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
}
#----------------------------INTERACTIVE MODELS WITH 2 COVARIATES (1 BEING SITE)
model_beta_2mult.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + bSite[x2[i]] + bSite.Int[x2[i]]*x1[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # Interaction distrib
  bSite.Int[1] <- 0
  for(j in 2:3){
    bSite.Int[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
}
#-------------------------------------------INTERACTIVE MODELS WITH 3 COVARIATES
model_beta_3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x1[i]*x2[i] + b5*x1[i]*x3[i] + b6*x2[i]*x3[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
  b6 ~ dnorm(0, 0.0001)
}
#----------------------------INTERACTIVE MODELS WITH 3 COVARIATES (1 BEING SITE)
model_beta_3.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + bSite[x3[i]] + b4*x1[i]*x2[i] + bSite.Int1[x3[i]]*x1[i] + bSite.Int2[x3[i]]*x2[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # Interaction distrib
  bSite.Int1[1] <- 0
  for(j in 2:3){
    bSite.Int1[j] ~ dnorm(0, 0.0001)
  }
  # Interaction distrib
  bSite.Int2[1] <- 0
  for(j in 2:3){
    bSite.Int2[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
}
#--------------------------------MODELS WITH 3 COVARIATES AND ONLY 1 INTERACTION
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
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b4 ~ dnorm(0, 0.0001)
}
#--------------MODELS WITH 3 COVARIATES AND ONLY 1 INTERACTION (BEING WITH SITE)
model_beta_3.2.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + bSite[x3[i]] + bSite.Int[x3[i]]*x2[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # Interaction distrib
  bSite.Int[1] <- 0
  for(j in 2:3){
    bSite.Int[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
}
#----------------------------------------------ADDITIVE MODELS WITH 3 COVARIATES
model_beta_3.3 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + b3*x3[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
}
#-------------------------------ADDITIVE MODELS WITH 3 COVARIATES (1 BEING SITE)
model_beta_3.3.site <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + bSite[x3[i]]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
}

#-------------------------------ADDITIVE MODELS WITH 4 COVARIATES (1 BEING SITE)
model_beta_4 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + bSite[x4[i]]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
}
#----------------------------INTERACTIVE MODELS WITH 4 COVARIATES (1 BEING SITE)
model_beta_4.2 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i]  <- (1-mu[i]) * phi
    logit(mu[i]) <- b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + bSite[x4[i]] + 
      b5*x2[i]*x3[i] + bSite.Int1[x4[i]]*x2[i] + bSite.Int2[x4[i]]*x3[i]
    LogLik[i] <- log(dbeta(y[i], alpha[i], beta[i]))
  }
  # cat distrib
  bSite[1] <- 0
  for(j in 2:3){
    bSite[j] ~ dnorm(0, 0.0001)
  }
  # Interaction distrib
  bSite.Int1[1] <- 0
  for(j in 2:3){
    bSite.Int1[j] ~ dnorm(0, 0.0001)
  }
  # Interaction distrib
  bSite.Int2[1] <- 0
  for(j in 2:3){
    bSite.Int2[j] ~ dnorm(0, 0.0001)
  }
  ## Specify priors
  phi ~ dgamma(0.1,0.1)
  b0 ~ dunif(log(.15), log(.75))
  b1 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  b5 ~ dnorm(0, 0.0001)
}

####---------------------------------GAMMA----------------------------------####

#---------------------------------------------------------------------NULL MODEL
#bayesian parameterization of gamma regression w/ 1 predictor
model_gamma_0 <- function() {
  ## Specify likelihood
  for(i in 1:length(x1)){
    y[i] ~ dgamma(shape, shape / exp(p[i]))
    p[i] <- alpha
    LogLik[i] <- log(dgamma(y[i], shape, shape / exp(p[i])))
  }
  ## Specify priors
  alpha ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 100)
}
#--------------------------------------------------------MODELS WITH 1 COVARIATE
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

#----------------------------------------------ADDITIVE MODELS WITH 2 COVARIATES
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
#-------------------------------------------INTERACTIVE MODELS WITH 2 COVARIATES
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

#-------------------------------------------INTERACTIVE MODELS WITH 3 COVARIATES
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
#--------------------------------MODELS WITH 3 COVARIATES AND ONLY 1 INTERACTION
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
#----------------------------------------------ADDITIVE MODELS WITH 3 COVARIATES
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

#-------------------------------------------INTERACTIVE MODELS WITH 4 COVARIATES
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





####MODELS####
####------------------------Probability of Infection------------------------####


# p(infection) ~ 1
infection.null <- data.frame(y = diagnosed$infection, x1 = 0)
mod.infection.null=jags(model.file=model_binom_0,data=infection.null,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("alpha","LogLik"))

# p(infection) ~ SITE
infection.site <- data.frame(y = diagnosed$infection, x1 = diagnosed$sitepairnum)
mod.infection.site=jags(model.file=model_binom_1.site,data=infection.site,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("bSite","LogLik"))

# p(infection) ~ SVL
infection.svl <- data.frame(y = diagnosed$infection, x1 = diagnosed$SVL)
mod.infection.svl=jags(model.file=model_binom_1,data=infection.svl,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

# p(infection) ~ SEX
infection.sex <- data.frame(y = diagnosed$infection, x1 = diagnosed$SEX)
mod.infection.sex=jags(model.file=model_binom_1,data=infection.sex,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

# p(infection) ~ HABITAT
infection.habitat <- data.frame(y = diagnosed$infection, 
                                x1 = diagnosed$habitat)
mod.infection.habitat=jags(model.file=model_binom_1,data=infection.habitat,
                           n.chains = nc, n.thin = nt,  n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","LogLik"))

# p(infection)~ SVL + HABITAT
infection.svl.habitat <- data.frame(y = diagnosed$infection, 
                                    x1 = diagnosed$SVL, x2 = diagnosed$habitat)
mod.infection.svl.habitat=jags(model.file=model_binom_2add,data=infection.svl.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","LogLik"))

# p(infection) ~ (SVL*HABITAT)
infection.svl.habitat <- data.frame(y = diagnosed$infection, 
                                    x1 = diagnosed$SVL, x2 = diagnosed$habitat)
mod.infection.svl..habitat=jags(model.file=model_binom_2mult,data=infection.svl.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# p(infection) ~ SEX + HABITAT
infection.sex.habitat <- data.frame(y = diagnosed$infection,  
                                    x1 = diagnosed$SEX, x2 = diagnosed$habitat)
mod.infection.sex.habitat=jags(model.file=model_binom_2add,data=infection.sex.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","LogLik"))

# p(infection) ~ (SEX*HABITAT)
infection.sex.habitat <- data.frame(y = diagnosed$infection,  
                                    x1 = diagnosed$SEX, x2 = diagnosed$habitat)
mod.infection.sex..habitat=jags(model.file=model_binom_2mult,data=infection.sex.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# p(infection) ~ SEX + SITE
infection.site.sex <- data.frame(y = diagnosed$infection, 
                                 x1 = diagnosed$SEX, x2 = diagnosed$sitepairnum)
mod.infection.site.sex=jags(model.file=model_binom_2add.site,data=infection.site.sex,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha","b1","bSite","LogLik"))

# p(infection) ~ (SEX*SITE)
infection.site.sex <- data.frame(y = diagnosed$infection, 
                                x1 = diagnosed$SEX, x2 = diagnosed$sitepairnum)
mod.infection.site..sex=jags(model.file=model_binom_2mult.site,data=infection.site.sex,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha","b1","bSite","bSite.Int","LogLik"))

# p(infection) ~ HABITAT + SITE
infection.site.habitat <- data.frame(y = diagnosed$infection, 
                                     x1 = diagnosed$habitat, x2 = diagnosed$sitepairnum)
mod.infection.site.habitat=jags(model.file=model_binom_2add.site,data=infection.site.habitat,
                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                parameters.to.save=c("alpha","b1","bSite","LogLik"))

# p(infection) ~ (HABITAT*SITE)
infection.site.habitat <- data.frame(y = diagnosed$infection, 
                                x1 = diagnosed$habitat, x2 = diagnosed$sitepairnum)
mod.infection.site..habitat=jags(model.file=model_binom_2mult.site,data=infection.svl.sex,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha","b1","bSite","bSite.Int","LogLik"))

# p(infection) ~ SEX + SITE + HABITAT
infection.sex.site.habitat <- data.frame(y = diagnosed$infection,
                                         x1 = diagnosed$SEX, x2 = diagnosed$habitat, x3 = diagnosed$sitepairnum)
mod.infection.sex.site.habitat=jags(model.file=model_binom_3add.site,data=infection.sex.site.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("alpha","b1","b2","bSite","LogLik"))

# p(infection) ~ (SEX*SITE) + (SEX*HABITAT) + (SITE*HABITAT)
infection.sex.site.habitat <- data.frame(y = diagnosed$infection,
                                         x1 = diagnosed$SEX, x2 = diagnosed$habitat, x3 = diagnosed$sitepairnum)
mod.infection.sex..site..habitat=jags(model.file=model_binom_3int.site,data=infection.sex.site.habitat,
                                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                      parameters.to.save=c("alpha","b1","b2","bSite","b4","bSite.Int1","bSite.Int2","LogLik"))

####-----------------------------Body Condition-----------------------------####

###Filter for males with complete tails
completetails <- diagnosed %>% filter(SEX == 1) %>%filter(TAIL=="C")

###Calculate body condition index
completetails$bodycond <- NA
male.mass.svl <- lm(MASS~SVL, data=completetails, na.action=na.exclude)
completetails$bodycond <- resid(male.mass.svl)


# body condition ~ 1
male.bodycond.null <- data.frame(y = completetails$bodycond,  x1 = 0)
mod.male.bodycond.null=jags(model.file=model_linear_0,data=male.bodycond.null,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha","LogLik"))

# body condition ~ INFECTION
male.bodycond.infection <- data.frame(y = completetails$bodycond, x1 = completetails$infection)
mod.male.bodycond.infection=jags(model.file=model_linear_1,data=male.bodycond.infection,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("alpha","b1","LogLik"))

# body condition ~ HABITAT
male.bodycond.habitat <- data.frame(y = completetails$bodycond, x1 = completetails$habitat)
mod.male.bodycond.habitat=jags(model.file=model_linear_1,data=male.bodycond.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","LogLik"))

# body condition ~ SITE
male.bodycond.site <- data.frame(y = completetails$bodycond, x1 = completetails$sitepairnum)
mod.male.bodycond.site=jags(model.file=model_linear_1.site,data=male.bodycond.site,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("bSite","LogLik"))

# body condition ~ HABITAT + SITE
male.bodycond.habitat.site <- data.frame(y = completetails$bodycond,
                                              x1 = completetails$habitat, x2 = completetails$sitepairnum)
mod.male.bodycond.habitat.site=jags(model.file=model_linear_2add.site,data=male.bodycond.habitat.site,
                                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                         parameters.to.save=c("alpha","b1","bSite","LogLik"))

# body condition ~ (HABITAT*SITE)
male.bodycond.habitat.site <- data.frame(y = completetails$bodycond,
                                         x1 = completetails$habitat, x2 = completetails$sitepairnum)
mod.male.bodycond.habitat..site=jags(model.file=model_linear_2mult.site,data=male.bodycond.habitat.site,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("alpha","b1","bSite","bSite.Int","LogLik"))

# body condition ~ INFECTION + HABITAT
male.bodycond.infection.habitat <- data.frame(y = completetails$bodycond,
                                              x1 = completetails$infection, x2 = completetails$habitat)
mod.male.bodycond.infection.habitat=jags(model.file=model_linear_2add,data=male.bodycond.infection.habitat,
                                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                         parameters.to.save=c("alpha","b1","b2","LogLik"))

# body condition ~ INFECTION + SITE
male.bodycond.infection.site <- data.frame(y = completetails$bodycond,
                                              x1 = completetails$infection, x2 = completetails$sitepairnum)
mod.male.bodycond.infection.site=jags(model.file=model_linear_2add.site,data=male.bodycond.infection.site,
                                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                         parameters.to.save=c("alpha","b1","bSite","LogLik"))

# body condition ~ INFECTION + SITE + HABITAT
male.bodycond.infection.habitat.site <- data.frame(y = completetails$bodycond,
                                           x1 = completetails$infection, x2 = completetails$habitat, x3 = completetails$sitepairnum)
mod.male.bodycond.infection.habitat.site=jags(model.file=model_linear_3.site,data=male.bodycond.infection.habitat.site,
                                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                      parameters.to.save=c("alpha","b1","b2","bSite","LogLik"))

# body condition ~ INFECTION + (SITE*HABITAT)
male.bodycond.infection.habitat..site <- data.frame(y = completetails$bodycond,
                                           x1 = completetails$infection, x2 = completetails$habitat, x3 = completetails$sitepairnum)
mod.male.bodycond.infection.habitat..site=jags(model.file=model_linear_3.2.site,data=male.bodycond.infection.habitat..site,
                                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                      parameters.to.save=c("alpha","b1","b2","bSite","bSite.Int","LogLik"))



####-------------------------------Hematocrit-------------------------------####
Hct.omit <- filter(diagnosed, !is.na(Hct))
# Hct ~ 1
Hct.null <- data.frame(y = Hct.omit$Hct, x1 = 0)
mod.Hct.null=jags(model.file=model_beta_1,data=Hct.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("b0","LogLik"))

# Hct ~ INFECTION
Hct.infection <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$infection)
mod.Hct.infection=jags(model.file=model_beta_1,data=Hct.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","LogLik"))

# Hct ~ SITE
Hct.site <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$sitepairnum)
mod.Hct.site=jags(model.file=model_beta_1.site,data=Hct.site,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("bSite","LogLik"))

# Hct ~ SVL
Hct.svl <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SVL)
mod.Hct.svl=jags(model.file=model_beta_1,data=Hct.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("b0","b1","LogLik"))

# Hct ~ Sex
Hct.sex <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SEX)
mod.Hct.sex=jags(model.file=model_beta_1,data=Hct.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("b0","b1","LogLik"))

# Hct ~ HABITAT
Hct.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$habitat)
mod.Hct.habitat=jags(model.file=model_beta_1,data=Hct.habitat,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("b0","b1","LogLik"))

# Hct ~ SVL + HABITAT
Hct.svl.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SVL, x2 = Hct.omit$habitat)
mod.Hct.svl.habitat=jags(model.file=model_beta_2add,data=Hct.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik"))

# Hct ~ (SVL*HABITAT)
Hct.svl.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SVL, x2 = Hct.omit$habitat)
mod.Hct.svl..habitat=jags(model.file=model_beta_2mult,data=Hct.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","b3","LogLik"))

# Hct ~ SEX + HABITAT
Hct.sex.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SEX, x2 = Hct.omit$habitat)
mod.Hct.sex.habitat=jags(model.file=model_beta_2add,data=Hct.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","LogLik"))

# Hct ~ (SEX*HABITAT)
Hct.sex.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$SEX, x2 = Hct.omit$habitat)
mod.Hct.sex..habitat=jags(model.file=model_beta_2mult,data=Hct.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","b3","LogLik"))

# Hct ~ SEX + SITE
Hct.site.sex <- data.frame(y = Hct.omit$Hct, 
                           x1 = Hct.omit$SEX, x2 = Hct.omit$sitepairnum)
mod.Hct.site.sex=jags(model.file=model_beta_2add.site,data=Hct.site.sex,
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("b0","b1","bSite","LogLik"))

# Hct ~ (SEX*SITE)
Hct.site.sex <- data.frame(y = Hct.omit$Hct, 
                                 x1 = Hct.omit$SEX, x2 = Hct.omit$sitepairnum)
mod.Hct.site..sex=jags(model.file=model_beta_2mult.site,data=Hct.site.sex,
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                             parameters.to.save=c("b0","b1","bSite","bSite.Int","LogLik"))

# Hct ~ HABITAT + SITE
Hct.site.habitat <- data.frame(y = Hct.omit$Hct, 
                               x1 = Hct.omit$habitat, x2 = Hct.omit$sitepairnum)
mod.Hct.site.habitat=jags(model.file=model_beta_2add.site,data=Hct.site.habitat,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("b0","b1","bSite","LogLik"))

# Hct ~ (HABITAT*SITE)
Hct.site.habitat <- data.frame(y = Hct.omit$Hct, 
                                     x1 = Hct.omit$habitat, x2 = Hct.omit$sitepairnum)
mod.Hct.site..habitat=jags(model.file=model_beta_2mult.site,data=Hct.site.habitat,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("b0","b1","bSite","bSite.Int","LogLik"))

# Hct ~ SEX + SITE + HABITAT
Hct.sex.site.habitat <- data.frame(y = Hct.omit$Hct,
                                         x1 = Hct.omit$SEX, x2 = Hct.omit$habitat, x3 = Hct.omit$sitepairnum)
mod.Hct.sex.site.habitat=jags(model.file=model_beta_3.3.site,data=Hct.sex.site.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("b0","b1","b2","bSite","LogLik"))

# Hct ~ (SEX*SITE) + (SEX*HABITAT) + (SITE*HABITAT)
Hct.sex.site.habitat <- data.frame(y = Hct.omit$Hct,
                                         x1 = Hct.omit$SEX, x2 = Hct.omit$habitat, x3 = Hct.omit$sitepairnum)
mod.Hct.sex..site..habitat=jags(model.file=model_beta_3.site,data=Hct.sex.site.habitat,
                                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                      parameters.to.save=c("b0","b1","b2","bSite","b4","bSite.Int1","bSite.Int2","LogLik"))


# Hct ~ INFECTION + SVL
Hct.infection.svl <- data.frame(y = Hct.omit$Hct, 
                                x1 = Hct.omit$infection, x2 = Hct.omit$SVL)
mod.Hct.infection.svl=jags(model.file=model_beta_2add,data=Hct.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2","LogLik"))

# Hct ~ INFECTION + SEX
Hct.infection.sex <- data.frame(y = Hct.omit$Hct, 
                                x1 = Hct.omit$infection, x2 = Hct.omit$SEX)
mod.Hct.infection.sex=jags(model.file=model_beta_2add,data=Hct.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2","LogLik"))

# Hct ~ INFECTION + HABITAT
Hct.infection.habitat <- data.frame(y = Hct.omit$Hct, 
                                    x1 = Hct.omit$infection, x2 = Hct.omit$habitat)
mod.Hct.infection.habitat=jags(model.file=model_beta_2add,data=Hct.infection.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("b0","b1","b2","LogLik"))

# Hct ~ INFECTION + SITE
Hct.infection.site <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$infection, x2 = Hct.omit$sitepairnum)
mod.Hct.infection.site=jags(model.file=model_beta_2add.site,data=Hct.infection.site,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("b0","b1","bSite","LogLik"))


# Hct ~ INFECTION + SVL + HABITAT
Hct.infection.svl.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$infection, x2 = Hct.omit$SVL, x3 = Hct.omit$habitat)
mod.Hct.infection.svl.habitat=jags(model.file=model_beta_3.3,data=Hct.infection.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("b0","b1","b2","b3","LogLik"))

# Hct ~ INFECTION + (SVL*HABITAT)
Hct.infection.svl.habitat <- data.frame(y = Hct.omit$Hct, 
                                        x1 = Hct.omit$infection, x2 = Hct.omit$SVL,
                                        x3 = Hct.omit$habitat)
mod.Hct.infection.svl..habitat=jags(model.file=model_beta_3.2,data=Hct.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))

# Hct ~ INFECTION + SEX + HABITAT
Hct.infection.sex.habitat <- data.frame(y = Hct.omit$Hct, 
                                        x1 = Hct.omit$infection, x2 = Hct.omit$SEX,
                                        x3 = Hct.omit$habitat)
mod.Hct.infection.sex.habitat=jags(model.file=model_beta_3.3,data=Hct.infection.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("b0","b1","b2","b3","LogLik"))

# Hct ~ INFECTION + (SEX*HABITAT)
Hct.infection.sex.habitat <- data.frame(y = Hct.omit$Hct, x1 = Hct.omit$infection, x2 = Hct.omit$SEX, x3 = Hct.omit$habitat)
mod.Hct.infection.sex..habitat=jags(model.file=model_beta_3.2,data=Hct.infection.sex.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("b0","b1","b2","b3","b4","LogLik"))

# Hct ~ INFECTION + SEX + SITE
Hct.infection.site.sex <- data.frame(y = Hct.omit$Hct, 
                                     x1 = Hct.omit$infection, x2 = Hct.omit$SEX, x3 = Hct.omit$sitepairnum)
mod.Hct.infection.site.sex=jags(model.file=model_beta_3.3.site,data=Hct.infection.site.sex,
                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                parameters.to.save=c("b0","b1","b2","bSite","LogLik"))

# Hct ~ INFECTION + (SEX*SITE)
Hct.infection.site.sex <- data.frame(y = Hct.omit$Hct, 
                                     x1 = Hct.omit$infection, x2 = Hct.omit$SEX, x3 = Hct.omit$sitepairnum)
mod.Hct.infection.site..sex=jags(model.file=model_beta_3.2.site,data=Hct.infection.site.sex,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("b0","b1","b2","bSite","bSite.Int","LogLik"))

# Hct ~ INFECTION + HABITAT + SITE
Hct.infection.site.habitat <- data.frame(y = Hct.omit$Hct, 
                                         x1 = Hct.omit$infection, x2 = Hct.omit$habitat, x3 = Hct.omit$sitepairnum)
mod.Hct.infection.site.habitat=jags(model.file=model_beta_3.3.site,data=Hct.infection.site.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("b0","b1","b2","bSite","LogLik"))

# Hct ~ INFECTION + (HABITAT*SITE)
Hct.infection.site.habitat <- data.frame(y = Hct.omit$Hct, 
                                         x1 = Hct.omit$infection, x2 = Hct.omit$habitat, x3 = Hct.omit$sitepairnum)
mod.Hct.infection.site..habitat=jags(model.file=model_beta_3.2.site,data=Hct.infection.site.habitat,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("b0","b1","b2","bSite","bSite.Int","LogLik"))

# Hct ~ INFECTION + SEX + SITE + HABITAT
Hct.infection.sex.site.habitat <- data.frame(y = Hct.omit$Hct,
                                   x1 = Hct.omit$infection, x2 = Hct.omit$SEX, x3 = Hct.omit$habitat, x4 = Hct.omit$sitepairnum)
mod.Hct.infection.sex.site.habitat=jags(model.file=model_beta_4,data=Hct.infection.sex.site.habitat,
                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                              parameters.to.save=c("b0","b1","b2","b3","bSite","LogLik"))

# Hct ~ INFECTION + (SEX*SITE) + (SEX*HABITAT) + (SITE*HABITAT)
Hct.infection.sex.site.habitat <- data.frame(y = Hct.omit$Hct,
                                   x1 = Hct.omit$infection, x2 = Hct.omit$SEX, x3 = Hct.omit$habitat, x4 = Hct.omit$sitepairnum)
mod.Hct.infection.sex..site..habitat=jags(model.file=model_beta_4.2,data=Hct.infection.sex.site.habitat,
                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                parameters.to.save=c("b0","b1","b2","b3","bSite","b5","bSite.Int1","bSite.Int2","LogLik"))

####--------------------------Principal Component 1-------------------------####
PC1.omit <- filter(diagnosed, !is.na(PC1))

# PC1 ~ 1
PC1.null <- data.frame(y = PC1.omit$PC1, x1 = 0)
mod.PC1.null=jags(model.file=model_linear_0,data=PC1.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("alpha","LogLik"))

# PC1 ~ INFECTION
PC1.infection <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$infection)
mod.PC1.infection=jags(model.file=model_linear_1,data=PC1.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

# PC1 ~ SVL
PC1.svl <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SVL)
mod.PC1.svl=jags(model.file=model_linear_1,data=PC1.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1","LogLik"))

# PC1 ~ SEX
PC1.sex <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SEX)
mod.PC1.sex=jags(model.file=model_linear_1,data=PC1.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1","LogLik"))

# PC1 ~ HABITAT
PC1.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$habitat)
mod.PC1.habitat=jags(model.file=model_linear_1,data=PC1.habitat,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("alpha","b1","LogLik"))

# PC1 ~ SITE
PC1.site <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$sitepairnum)
mod.PC1.site=jags(model.file=model_linear_1.site,data=PC1.site,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("bSite","LogLik"))

# PC1 ~ SVL + HABITAT
PC1.svl.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SVL, x2 = PC1.omit$habitat)
mod.PC1.svl.habitat=jags(model.file=model_linear_2add,data=PC1.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC1 ~ (SVL*HABITAT)
PC1.svl.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SVL, x2 = PC1.omit$habitat)
mod.PC1.svl..habitat=jags(model.file=model_linear_2mult,data=PC1.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# PC1 ~ SEX + HABITAT
PC1.sex.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SEX, x2 = PC1.omit$habitat)
mod.PC1.sex.habitat=jags(model.file=model_linear_2add,data=PC1.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC1 ~ (SEX*HABITAT)
PC1.sex.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$SEX, x2 = PC1.omit$habitat)
mod.PC1.sex..habitat=jags(model.file=model_linear_2mult,data=PC1.sex.habitat,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# PC1 ~ SEX + SITE
PC1.site.sex <- data.frame(y = PC1.omit$PC1, 
                           x1 = PC1.omit$SEX, x2 = PC1.omit$sitepairnum)
mod.PC1.site.sex=jags(model.file=model_linear_2add.site,data=PC1.site.sex,
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("alpha","b1","bSite","LogLik"))

# PC1 ~ (SEX*SITE)
PC1.site.sex <- data.frame(y = PC1.omit$PC1, 
                           x1 = PC1.omit$SEX, x2 = PC1.omit$sitepairnum)
mod.PC1.site..sex=jags(model.file=model_linear_2mult.site,data=PC1.site.sex,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","bSite","bSite.Int","LogLik"))

# PC1 ~ HABITAT + SITE
PC1.site.habitat <- data.frame(y = PC1.omit$PC1, 
                               x1 = PC1.omit$habitat, x2 = PC1.omit$sitepairnum)
mod.PC1.site.habitat=jags(model.file=model_linear_2add.site,data=PC1.site.habitat,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("alpha","b1","bSite","LogLik"))

# PC1 ~ (HABITAT*SITE)
PC1.site.habitat <- data.frame(y = PC1.omit$PC1, 
                               x1 = PC1.omit$habitat, x2 = PC1.omit$sitepairnum)
mod.PC1.site..habitat=jags(model.file=model_linear_2mult.site,data=PC1.site.habitat,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","bSite","bSite.Int","LogLik"))

# PC1 ~ SEX + SITE + HABITAT
PC1.sex.site.habitat <- data.frame(y = PC1.omit$PC1,
                                   x1 = PC1.omit$SEX, x2 = PC1.omit$habitat, x3 = PC1.omit$sitepairnum)
mod.PC1.sex.site.habitat=jags(model.file=model_linear_3.3.site,data=PC1.sex.site.habitat,
                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                              parameters.to.save=c("alpha","b1","b2","bSite","LogLik"))

# PC1 ~ (SEX*SITE) + (SEX*HABITAT) + (SITE*HABITAT)
PC1.sex.site.habitat <- data.frame(y = PC1.omit$PC1,
                                   x1 = PC1.omit$SEX, x2 = PC1.omit$habitat, x3 = PC1.omit$sitepairnum)
mod.PC1.sex..site..habitat=jags(model.file=model_linear_3.site,data=PC1.sex.site.habitat,
                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                parameters.to.save=c("alpha","b1","b2","bSite","b4","bSite.Int1","bSite.Int2","LogLik"))


# PC1 ~ INFECTION + SVL
PC1.infection.svl <- data.frame(y = PC1.omit$PC1, 
                                x1 = PC1.omit$infection, x2 = PC1.omit$SVL)
mod.PC1.infection.svl=jags(model.file=model_linear_2add,data=PC1.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC1 ~ INFECTION + SEX
PC1.infection.sex <- data.frame(y = PC1.omit$PC1, 
                                x1 = PC1.omit$infection, x2 = PC1.omit$SEX)
mod.PC1.infection.sex=jags(model.file=model_linear_2add,data=PC1.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC1 ~ INFECTION + HABITAT
PC1.infection.habitat <- data.frame(y = PC1.omit$PC1, 
                                    x1 = PC1.omit$infection, x2 = PC1.omit$habitat)
mod.PC1.infection.habitat=jags(model.file=model_linear_2add,data=PC1.infection.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC1 ~ INFECTION + SITE
PC1.infection.site <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$infection, x2 = PC1.omit$sitepairnum)
mod.PC1.infection.site=jags(model.file=model_linear_2add.site,data=PC1.infection.site,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha","b1","bSite","LogLik"))


# PC1 ~ INFECTION + SVL + HABITAT
PC1.infection.svl.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$infection, x2 = PC1.omit$SVL, x3 = PC1.omit$habitat)
mod.PC1.infection.svl.habitat=jags(model.file=model_linear_3.3,data=PC1.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# PC1 ~ INFECTION + (SVL*HABITAT)
PC1.infection.svl.habitat <- data.frame(y = PC1.omit$PC1, 
                                        x1 = PC1.omit$infection, x2 = PC1.omit$SVL,
                                        x3 = PC1.omit$habitat)
mod.PC1.infection.svl..habitat=jags(model.file=model_linear_3.2,data=PC1.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

# PC1 ~ INFECTION + SEX + HABITAT
PC1.infection.sex.habitat <- data.frame(y = PC1.omit$PC1, 
                                        x1 = PC1.omit$infection, x2 = PC1.omit$SEX,
                                        x3 = PC1.omit$habitat)
mod.PC1.infection.sex.habitat=jags(model.file=model_linear_3.3,data=PC1.infection.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# PC1 ~ INFECTION + (SEX*HABITAT)
PC1.infection.sex.habitat <- data.frame(y = PC1.omit$PC1, x1 = PC1.omit$infection, x2 = PC1.omit$SEX, x3 = PC1.omit$habitat)
mod.PC1.infection.sex..habitat=jags(model.file=model_linear_3.2,data=PC1.infection.sex.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

# PC1 ~ INFECTION + (SEX*SITE)
PC1.infection.site.sex <- data.frame(y = PC1.omit$PC1, 
                                     x1 = PC1.omit$infection, x2 = PC1.omit$SEX, x3 = PC1.omit$sitepairnum)
mod.PC1.infection.site..sex=jags(model.file=model_linear_3.2.site,data=PC1.infection.site.sex,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("alpha","b1","b2","bSite","bSite.Int","LogLik"))

# PC1~ INFECTION + SEX + SITE
PC1.infection.site.sex <- data.frame(y = PC1.omit$PC1, 
                                     x1 = PC1.omit$infection, x2 = PC1.omit$SEX, x3 = PC1.omit$sitepairnum)
mod.PC1.infection.site.sex=jags(model.file=model_linear_3.3.site,data=PC1.infection.site.sex,
                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                parameters.to.save=c("alpha","b1","b2","bSite","LogLik"))

# PC1 ~ INFECTION + (HABITAT*SITE)
PC1.infection.site.habitat <- data.frame(y = PC1.omit$PC1, 
                                         x1 = PC1.omit$infection, x2 = PC1.omit$habitat, x3 = PC1.omit$sitepairnum)
mod.PC1.infection.site..habitat=jags(model.file=model_linear_3.2.site,data=PC1.infection.site.habitat,
                                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                     parameters.to.save=c("alpha","b1","b2","bSite","bSite.Int","LogLik"))

# PC1 ~ INFECTION + HABITAT + SITE
PC1.infection.site.habitat <- data.frame(y = PC1.omit$PC1, 
                                         x1 = PC1.omit$infection, x2 = PC1.omit$habitat, x3 = PC1.omit$sitepairnum)
mod.PC1.infection.site.habitat=jags(model.file=model_linear_3.3.site,data=PC1.infection.site.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("alpha","b1","b2","bSite","LogLik"))

# PC1 ~ INFECTION + SEX + SITE + HABITAT
PC1.infection.sex.site.habitat <- data.frame(y = PC1.omit$PC1,
                                             x1 = PC1.omit$infection, x2 = PC1.omit$SEX, x3 = PC1.omit$habitat, x4 = PC1.omit$sitepairnum)
mod.PC1.infection.sex.site.habitat=jags(model.file=model_linear_4,data=PC1.infection.sex.site.habitat,
                                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                        parameters.to.save=c("alpha","b1","b2","b3","bSite","LogLik"))

# PC1 ~ INFECTION + (SEX*SITE) + (SEX*HABITAT) + (SITE*HABITAT)
PC1.infection.sex.site.habitat <- data.frame(y = PC1.omit$PC1,
                                             x1 = PC1.omit$infection, x2 = PC1.omit$SEX, x3 = PC1.omit$habitat, x4 = PC1.omit$sitepairnum)
mod.PC1.infection.sex..site..habitat=jags(model.file=model_linear_4.2,data=PC1.infection.sex.site.habitat,
                                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                          parameters.to.save=c("alpha","b1","b2","b3","bSite","b5","bSite.Int1","bSite.Int2","LogLik"))

####--------------------------Principal Component 2-------------------------####

PC2.omit <- filter(diagnosed, !is.na(PC2))
# PC2 ~ 1
PC2.null <- data.frame(y = PC2.omit$PC2, x1 = 0)
mod.PC2.null=jags(model.file=model_linear_0,data=PC2.null,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("alpha","LogLik"))

# PC2 ~ INFECTION
PC2.infection <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$infection)
mod.PC2.infection=jags(model.file=model_linear_1,data=PC2.infection,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

# PC2 ~ SVL
PC2.svl <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL)
mod.PC2.svl=jags(model.file=model_linear_1,data=PC2.svl,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1","LogLik"))

# PC2 ~ SEX
PC2.sex <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SEX)
mod.PC2.sex=jags(model.file=model_linear_1,data=PC2.sex,
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 parameters.to.save=c("alpha","b1","LogLik"))

# PC2 ~ HABITAT
PC2.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$habitat)
mod.PC2.habitat=jags(model.file=model_linear_1,data=PC2.habitat,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                     parameters.to.save=c("alpha","b1","LogLik"))


# PC2 ~ SITE
PC2.site <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$sitepairnum)
mod.PC2.site=jags(model.file=model_linear_1.site,data=PC2.site,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  parameters.to.save=c("bSite","LogLik"))

# PC2 ~ SVL + HABITAT
PC2.svl.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL, x2 = PC2.omit$habitat)
mod.PC2.svl.habitat=jags(model.file=model_linear_2add,data=PC2.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC2 ~ (SVL*HABITAT)
PC2.svl.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SVL, x2 = PC2.omit$habitat)
mod.PC2.svl..habitat=jags(model.file=model_linear_2mult,data=PC2.svl.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# PC2 ~ SEX + HABITAT
PC2.sex.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SEX, x2 = PC2.omit$habitat)
mod.PC2.sex.habitat=jags(model.file=model_linear_2add,data=PC2.sex.habitat,
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                         parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC2 ~ (SEX*HABITAT)
PC2.sex.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$SEX, x2 = PC2.omit$habitat)
mod.PC2.sex..habitat=jags(model.file=model_linear_2mult,data=PC2.sex.habitat,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# PC2 ~ SEX + SITE
PC2.site.sex <- data.frame(y = PC2.omit$PC2, 
                           x1 = PC2.omit$SEX, x2 = PC2.omit$sitepairnum)
mod.PC2.site.sex=jags(model.file=model_linear_2add.site,data=PC2.site.sex,
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parameters.to.save=c("alpha","b1","bSite","LogLik"))

# PC2 ~ (SEX*SITE)
PC2.site.sex <- data.frame(y = PC2.omit$PC2, 
                           x1 = PC2.omit$SEX, x2 = PC2.omit$sitepairnum)
mod.PC2.site..sex=jags(model.file=model_linear_2mult.site.site,data=PC2.site.sex,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","bSite","bSite.Int","LogLik"))

# PC2 ~ HABITAT + SITE
PC2.site.habitat <- data.frame(y = PC2.omit$PC2, 
                               x1 = PC2.omit$habitat, x2 = PC2.omit$sitepairnum)
mod.PC2.site.habitat=jags(model.file=model_linear_2add.site,data=PC2.site.habitat,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("alpha","b1","bSite","LogLik"))

# PC2 ~ (HABITAT*SITE)
PC2.site.habitat <- data.frame(y = PC2.omit$PC2, 
                               x1 = PC2.omit$habitat, x2 = PC2.omit$sitepairnum)
mod.PC2.site..habitat=jags(model.file=model_linear_2mult.site,data=PC2.site.habitat,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","bSite","bSite.Int","LogLik"))

# PC2 ~ SEX + SITE + HABITAT
PC2.sex.site.habitat <- data.frame(y = PC2.omit$PC2,
                                   x1 = PC2.omit$SEX, x2 = PC2.omit$habitat, x3 = PC2.omit$sitepairnum)
mod.PC2.sex.site.habitat=jags(model.file=model_linear_3.3.site,data=PC2.sex.site.habitat,
                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                              parameters.to.save=c("alpha","b1","b2","bSite","LogLik"))

# PC2 ~ (SEX*SITE) + (SEX*HABITAT) + (SITE*HABITAT)
PC2.sex.site.habitat <- data.frame(y = PC2.omit$PC2,
                                   x1 = PC2.omit$SEX, x2 = PC2.omit$habitat, x3 = PC2.omit$sitepairnum)
mod.PC2.sex..site..habitat=jags(model.file=model_linear_3.site,data=PC2.sex.site.habitat,
                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                parameters.to.save=c("alpha","b1","b2","bSite","b4","bSite.Int1","bSite.Int2","LogLik"))


# PC2 ~ INFECTION + SVL
PC2.infection.svl <- data.frame(y = PC2.omit$PC2, 
                                x1 = PC2.omit$infection, x2 = PC2.omit$SVL)
mod.PC2.infection.svl=jags(model.file=model_linear_2add,data=PC2.infection.svl,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC2 ~ INFECTION + SEX
PC2.infection.sex <- data.frame(y = PC2.omit$PC2, 
                                x1 = PC2.omit$infection, x2 = PC2.omit$SEX)
mod.PC2.infection.sex=jags(model.file=model_linear_2add,data=PC2.infection.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC2 ~ INFECTION + HABITAT
PC2.infection.habitat <- data.frame(y = PC2.omit$PC2, 
                                    x1 = PC2.omit$infection, x2 = PC2.omit$habitat)
mod.PC2.infection.habitat=jags(model.file=model_linear_2add,data=PC2.infection.habitat,
                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                               parameters.to.save=c("alpha","b1","b2","LogLik"))

# PC2 ~ INFECTION + SITE
PC2.infection.site <- data.frame(y = PC2.omit$PC2, x1 = PC1.omit$infection, x2 = PC2.omit$sitepairnum)
mod.PC2.infection.site=jags(model.file=model_linear_2add.site,data=PC2.infection.site,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha","b1","bSite","LogLik"))

# PC2 ~ INFECTION + SVL + HABITAT
PC2.infection.svl.habitat <- data.frame(y = PC2.omit$PC2, 
                                        x1 = PC2.omit$infection, x2 = PC2.omit$SVL,
                                        x3 = PC2.omit$habitat)
mod.PC2.infection.svl.habitat=jags(model.file=model_linear_3.3,data=PC2.infection.svl.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

# PC2 ~ INFECTION + (SVL*HABITAT)
PC2.infection.svl.habitat <- data.frame(y = PC2.omit$PC2, 
                                        x1 = PC2.omit$infection, x2 = PC2.omit$SVL,
                                        x3 = PC2.omit$habitat)
mod.PC2.infection.svl..habitat=jags(model.file=model_linear_3.2,data=PC2.infection.svl.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

# PC2 ~ INFECTION + SEX + HABITAT
PC2.infection.sex.habitat <- data.frame(y = PC2.omit$PC2, 
                                        x1 = PC2.omit$infection, x2 = PC2.omit$SEX,
                                        x3 = PC2.omit$habitat)
mod.PC2.infection.sex.habitat=jags(model.file=model_linear_3.3,data=PC2.infection.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","LogLik"))

# PC2 ~ INFECTION + (SEX*HABITAT)
PC2.infection.sex.habitat <- data.frame(y = PC2.omit$PC2, x1 = PC2.omit$infection, x2 = PC2.omit$SEX, x3 = PC2.omit$habitat)
mod.PC2.infection.sex..habitat=jags(model.file=model_linear_3.2,data=PC2.infection.sex.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))

# PC2 ~ INFECTION + SEX + SITE
PC2.infection.site.sex <- data.frame(y = PC2.omit$PC2, 
                                     x1 = PC2.omit$infection, x2 = PC2.omit$SEX, x3 = PC2.omit$sitepairnum)
mod.PC2.infection.site.sex=jags(model.file=model_linear_3.3.site,data=PC2.infection.site.sex,
                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                parameters.to.save=c("alpha","b1","b2","bSite","LogLik"))

# PC2 ~ INFECTION + (SEX*SITE)
PC2.infection.site.sex <- data.frame(y = PC2.omit$PC2, 
                                     x1 = PC2.omit$infection, x2 = PC2.omit$SEX, x3 = PC2.omit$sitepairnum)
mod.PC2.infection.site..sex=jags(model.file=model_linear_3.2.site,data=PC2.infection.site.sex,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("alpha","b1","b2","bSite","bSite.Int","LogLik"))

# PC2 ~ INFECTION + HABITAT + SITE
PC2.infection.site.habitat <- data.frame(y = PC2.omit$PC2, 
                                         x1 = PC2.omit$infection, x2 = PC2.omit$habitat, x3 = PC2.omit$sitepairnum)
mod.PC2.infection.site.habitat=jags(model.file=model_linear_3.3.site,data=PC2.infection.site.habitat,
                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                    parameters.to.save=c("alpha","b1","b2","bSite","LogLik"))

# PC2 ~ INFECTION + (HABITAT*SITE)
PC2.infection.site.habitat <- data.frame(y = PC2.omit$PC2, 
                                         x1 = PC2.omit$infection, x2 = PC2.omit$habitat, x3 = PC2.omit$sitepairnum)
mod.PC2.infection.site..habitat=jags(model.file=model_linear_3.2.site,data=PC2.infection.site.habitat,
                                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                     parameters.to.save=c("alpha","b1","b2","bSite","bSite.Int","LogLik"))

# PC2 ~ INFECTION + SEX + SITE + HABITAT
PC2.infection.sex.site.habitat <- data.frame(y = PC2.omit$PC2,
                                             x1 = PC2.omit$infection, x2 = PC2.omit$SEX, x3 = PC2.omit$habitat, x4 = PC2.omit$sitepairnum)
mod.PC2.infection.sex.site.habitat=jags(model.file=model_linear_4,data=PC2.infection.sex.site.habitat,
                                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                        parameters.to.save=c("alpha","b1","b2","b3","bSite","LogLik"))

# PC2 ~ INFECTION + (SEX*SITE) + (SEX*HABITAT) + (SITE*HABITAT)
PC2.infection.sex.site.habitat <- data.frame(y = PC2.omit$PC2,
                                             x1 = PC2.omit$infection, x2 = PC2.omit$SEX, x3 = PC2.omit$habitat, x4 = PC2.omit$sitepairnum)
mod.PC2.infection.sex..site..habitat=jags(model.file=model_linear_4.2,data=PC2.infection.sex.site.habitat,
                                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                          parameters.to.save=c("alpha","b1","b2","b3","bSite","b5","bSite.Int1","bSite.Int2","LogLik"))


####---------------Individual Dimensions of PCA (Preliminary)---------------####
##Preliminary analysis showed no clear effect of infection on these individual
##measurements, warranting their combination into PCA.
# ####Na####
# Sodium.omit <- filter(diagnosed, !is.na(Na))
# #Na~1
# Na.null <- data.frame(y = Sodium.omit$Na, x1 = 0)
# mod.Na.null=jags(model.file=model_gamma_0,data=Na.null,
#                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                  parameters.to.save=c("alpha","LogLik"))
# 
# #Na ~ infection
# Na.infection <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$infection)
# mod.Na.infection=jags(model.file=model_gamma_1,data=Na.infection,
#                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                       parameters.to.save=c("alpha","b1","LogLik"))
# 
# # Na ~ SVL
# Na.svl <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL)
# mod.Na.svl=jags(model.file=model_gamma_1,data=Na.svl,
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                 parameters.to.save=c("alpha","b1","LogLik"))
# # Na ~ Sex
# Na.sex <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SEX)
# mod.Na.sex=jags(model.file=model_gamma_1,data=Na.sex,
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                 parameters.to.save=c("alpha","b1","LogLik"))
# 
# # Na ~ Habitat
# Na.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$habitat)
# mod.Na.habitat=jags(model.file=model_gamma_1,data=Na.habitat,
#                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                     parameters.to.save=c("alpha","b1","LogLik"))
# 
# #Na~ (SVL*Sex)
# Na.svl.sex <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL, x2 = Sodium.omit$SEX)
# mod.Na.svl.sex=jags(model.file=model_gamma_2mult,data=Na.svl.sex,
#                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                     parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Na~ SVL + Sex
# Na.svl.sex <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL, x2 = Sodium.omit$SEX)
# mod.Na.svl.sex=jags(model.file=model_gamma_2add,data=Na.svl.sex,
#                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                     parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #Na~(SVL*Habitat)
# Na.svl.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL, x2 = Sodium.omit$habitat)
# mod.Na.svl.habitat=jags(model.file=model_gamma_2mult,data=Na.svl.habitat,
#                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                         parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Na~ SVL + Habitat
# Na.svl.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL, x2 = Sodium.omit$habitat)
# mod.Na.svl.habitat=jags(model.file=model_gamma_2add,data=Na.svl.habitat,
#                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                         parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #Na~(Sex)+(Habitat)
# Na.sex.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SEX, x2 = Sodium.omit$habitat)
# mod.Na.sex.habitat=jags(model.file=model_gamma_2add,data=Na.sex.habitat,
#                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                         parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #Na~(SVL*Sex)+(SVL*Habitat)
# Na.svl.sex.habitat <- data.frame(y = Sodium.omit$Na, x1 = Sodium.omit$SVL, 
#                                  x2 = Sodium.omit$SEX, x3 = Sodium.omit$habitat)
# mod.Na.svl.sex.habitat=jags(model.file=model_gamma_3,data=Na.svl.sex.habitat,
#                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                             parameters.to.save=c("alpha","b1","b2","b3","b4","b5","LogLik"))
# 
# # Na ~ INFECTION + SVL
# Na.infection.svl <- data.frame(y = Sodium.omit$Na, 
#                                x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL)
# mod.Na.infection.svl=jags(model.file=model_gamma_2add,data=Na.infection.svl,
#                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                           parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # Na ~ INFECTION + Sex
# Na.infection.sex <- data.frame(y = Sodium.omit$Na, 
#                                x1 = Sodium.omit$infection, x2 = Sodium.omit$SEX)
# mod.Na.infection.sex=jags(model.file=model_gamma_2add,data=Na.infection.sex,
#                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                           parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # Na ~ INFECTION + Habitat
# Na.infection.habitat <- data.frame(y = Sodium.omit$Na, 
#                                    x1 = Sodium.omit$infection, x2 = Sodium.omit$habitat)
# mod.Na.infection.habitat=jags(model.file=model_gamma_2add,data=Na.infection.habitat,
#                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                               parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #Na~INFECTION + (SVL*Sex)
# Na.infection.svl.sex <- data.frame(y = Sodium.omit$Na, 
#                                    x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL,
#                                    x3 = Sodium.omit$SEX)
# mod.Na.infection.svl.sex=jags(model.file=model_gamma_3.2,data=Na.infection.svl.sex,
#                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                               parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #Na~INFECTION + (SVL*Sex)
# Na.infection.svl.sex <- data.frame(y = Sodium.omit$Na, 
#                                    x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL,
#                                    x3 = Sodium.omit$SEX)
# mod.Na.infection.svl.sex=jags(model.file=model_gamma_3.2,data=Na.infection.svl.sex,
#                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                               parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #Na~INFECTION + (SVL*Habitat)
# Na.infection.svl.habitat <- data.frame(y = Sodium.omit$Na, 
#                                        x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL,
#                                        x3 = Sodium.omit$habitat)
# mod.Na.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=Na.infection.svl.habitat,
#                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                   parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #Na~INFECTION + (SVL*Habitat)
# Na.infection.svl.habitat <- data.frame(y = Sodium.omit$Na, 
#                                        x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL,
#                                        x3 = Sodium.omit$habitat)
# mod.Na.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=Na.infection.svl.habitat,
#                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                   parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #Na~INFECTION + (Sex)+(Habitat)
# Na.infection.sex.habitat <- data.frame(y = Sodium.omit$Na, 
#                                        x1 = Sodium.omit$infection, x2 = Sodium.omit$SEX,
#                                        x3 = Sodium.omit$habitat)
# mod.Na.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=Na.infection.sex.habitat,
#                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                   parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Na~INFECTION + (SVL*Sex)+(SVL*Habitat)
# Na.infection.svl.sex.habitat <- data.frame(y = Sodium.omit$Na, 
#                                            x1 = Sodium.omit$infection, x2 = Sodium.omit$SVL,
#                                            x3 = Sodium.omit$SEX, x4 = Sodium.omit$habitat)
# mod.Na.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=Na.infection.svl.sex.habitat,
#                                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                       parameters.to.save=c("alpha","b1","b2","b3",
#                                                            "b4","b5","b6","LogLik"))
# 
# 
# ####K####
# K.omit <- filter(diagnosed, !is.na(K))
# #K~1
# K.null <- data.frame(y = K.omit$K, x1 = 0)
# mod.K.null=jags(model.file=model_gamma_0,data=K.null,
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                 parameters.to.save=c("alpha","LogLik"))
# 
# #K ~ infection
# K.infection <- data.frame(y = K.omit$K, x1 = K.omit$infection)
# mod.K.infection=jags(model.file=model_gamma_1,data=K.infection,
#                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                      parameters.to.save=c("alpha","b1","LogLik"))
# 
# # K ~ SVL
# K.svl <- data.frame(y = K.omit$K, x1 = K.omit$SVL)
# mod.K.svl=jags(model.file=model_gamma_1,data=K.svl,
#                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                parameters.to.save=c("alpha","b1","LogLik"))
# 
# # K ~ Sex
# K.sex <- data.frame(y = K.omit$K, x1 = K.omit$SEX)
# mod.K.sex=jags(model.file=model_gamma_1,data=K.sex,
#                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                parameters.to.save=c("alpha","b1","LogLik"))
# 
# # K ~ Habitat
# K.habitat <- data.frame(y = K.omit$K, x1 = K.omit$habitat)
# mod.K.habitat=jags(model.file=model_gamma_1,data=K.habitat,
#                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                    parameters.to.save=c("alpha","b1","LogLik"))
# 
# #K~(SVL*Sex)
# K.svl.sex <- data.frame(y = K.omit$K, x1 = K.omit$SVL, x2 = K.omit$SEX)
# mod.K.svl.sex=jags(model.file=model_gamma_2mult,data=K.svl.sex,
#                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                    parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #K~(SVL*Habitat)
# K.svl.habitat <- data.frame(y = K.omit$K, x1 = K.omit$SVL, x2 = K.omit$habitat)
# mod.K.svl.habitat=jags(model.file=model_gamma_2mult,data=K.svl.habitat,
#                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                        parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #K~(Sex)+(Habitat)
# K.sex.habitat <- data.frame(y = K.omit$K, x1 = K.omit$SEX, x2 = K.omit$habitat)
# mod.K.sex.habitat=jags(model.file=model_gamma_2add,data=K.sex.habitat,
#                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                        parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #K~(SVL*Sex)+(SVL*Habitat)
# K.svl.sex.habitat <- data.frame(y = K.omit$K, x1 = K.omit$SVL, 
#                                 x2 = K.omit$SEX, x3 = K.omit$habitat)
# mod.K.svl.sex.habitat=jags(model.file=model_gamma_3,data=K.svl.sex.habitat,
#                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                            parameters.to.save=c("alpha","b1","b2","b3","b4","b5","LogLik"))
# 
# # K ~ INFECTION + SVL
# K.infection.svl <- data.frame(y = K.omit$K, 
#                               x1 = K.omit$infection, x2 = K.omit$SVL)
# mod.K.infection.svl=jags(model.file=model_gamma_2add,data=K.infection.svl,
#                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                          parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # K ~ INFECTION + Sex
# K.infection.sex <- data.frame(y = K.omit$K, 
#                               x1 = K.omit$infection, x2 = K.omit$SEX)
# mod.K.infection.sex=jags(model.file=model_gamma_2add,data=K.infection.sex,
#                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                          parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # K ~ INFECTION + Habitat
# K.infection.habitat <- data.frame(y = K.omit$K, 
#                                   x1 = K.omit$infection, x2 = K.omit$habitat)
# mod.K.infection.habitat=jags(model.file=model_gamma_2add,data=K.infection.habitat,
#                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                              parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #K~INFECTION + (SVL*Sex)
# K.infection.svl.sex <- data.frame(y = K.omit$K, 
#                                   x1 = K.omit$infection, x2 = K.omit$SVL,
#                                   x3 = K.omit$SEX)
# mod.K.infection.svl.sex=jags(model.file=model_gamma_3.2,data=K.infection.svl.sex,
#                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                              parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #K~INFECTION + (SVL*Habitat)
# K.infection.svl.habitat <- data.frame(y = K.omit$K, 
#                                       x1 = K.omit$infection, x2 = K.omit$SVL,
#                                       x3 = K.omit$habitat)
# mod.K.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=K.infection.svl.habitat,
#                                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                  parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #K~INFECTION + (Sex)+(Habitat)
# K.infection.sex.habitat <- data.frame(y = K.omit$K, 
#                                       x1 = K.omit$infection, x2 = K.omit$SEX,
#                                       x3 = K.omit$habitat)
# mod.K.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=K.infection.sex.habitat,
#                                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                  parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #K~INFECTION + (SVL*Sex)+(SVL*Habitat)
# K.infection.svl.sex.habitat <- data.frame(y = K.omit$K, 
#                                           x1 = K.omit$infection, x2 = K.omit$SVL,
#                                           x3 = K.omit$SEX, x4 = K.omit$habitat)
# mod.K.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=K.infection.svl.sex.habitat,
#                                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                      parameters.to.save=c("alpha","b1","b2","b3",
#                                                           "b4","b5","b6","LogLik"))
# ####Cl####
# Cl.omit <- filter(diagnosed, !is.na(Cl))
# #Cl~1
# Cl.null <- data.frame(y = Cl.omit$Cl, x1 = 0)
# mod.Cl.null=jags(model.file=model_gamma_0,data=Cl.null,
#                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                  parameters.to.save=c("alpha","LogLik"))
# 
# #Cl ~ infection
# Cl.infection <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$infection)
# mod.Cl.infection=jags(model.file=model_gamma_1,data=Cl.infection,
#                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                       parameters.to.save=c("alpha","b1","LogLik"))
# 
# # Cl ~ SVL
# Cl.svl <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL)
# mod.Cl.svl=jags(model.file=model_gamma_1,data=Cl.svl,
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                 parameters.to.save=c("alpha","b1","LogLik"))
# 
# # Cl ~ Sex
# Cl.sex <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SEX)
# mod.Cl.sex=jags(model.file=model_gamma_1,data=Cl.sex,
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                 parameters.to.save=c("alpha","b1","LogLik"))
# 
# # Cl ~ Habitat
# Cl.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$habitat)
# mod.Cl.habitat=jags(model.file=model_gamma_1,data=Cl.habitat,
#                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                     parameters.to.save=c("alpha","b1","LogLik"))
# 
# #Cl~(SVL*Sex)
# Cl.svl.sex <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL, x2 = Cl.omit$SEX)
# mod.Cl.svl.sex=jags(model.file=model_gamma_2mult,data=Cl.svl.sex,
#                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                     parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Cl~(SVL*Habitat)
# Cl.svl.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL, x2 = Cl.omit$habitat)
# mod.Cl.svl.habitat=jags(model.file=model_gamma_2mult,data=Cl.svl.habitat,
#                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                         parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Cl~(Sex)+(Habitat)
# Cl.sex.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SEX, x2 = Cl.omit$habitat)
# mod.Cl.sex.habitat=jags(model.file=model_gamma_2add,data=Cl.sex.habitat,
#                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                         parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #Cl~(SVL*Sex)+(SVL*Habitat)
# Cl.svl.sex.habitat <- data.frame(y = Cl.omit$Cl, x1 = Cl.omit$SVL, 
#                                  x2 = Cl.omit$SEX, x3 = Cl.omit$habitat)
# mod.Cl.svl.sex.habitat=jags(model.file=model_gamma_3,data=Cl.svl.sex.habitat,
#                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                             parameters.to.save=c("alpha","b1","b2","b3","b4","b5","LogLik"))
# 
# # Cl ~ INFECTION + SVL
# Cl.infection.svl <- data.frame(y = Cl.omit$Cl, 
#                                x1 = Cl.omit$infection, x2 = Cl.omit$SVL)
# mod.Cl.infection.svl=jags(model.file=model_gamma_2add,data=Cl.infection.svl,
#                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                           parameters.to.save=c("alpha","b1","b2","LogLik"))
# # Cl ~ INFECTION + Sex
# Cl.infection.sex <- data.frame(y = Cl.omit$Cl, 
#                                x1 = Cl.omit$infection, x2 = Cl.omit$SEX)
# mod.Cl.infection.sex=jags(model.file=model_gamma_2add,data=Cl.infection.sex,
#                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                           parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # Cl ~ INFECTION + Habitat
# Cl.infection.habitat <- data.frame(y = Cl.omit$Cl, 
#                                    x1 = Cl.omit$infection, x2 = Cl.omit$habitat)
# mod.Cl.infection.habitat=jags(model.file=model_gamma_2add,data=Cl.infection.habitat,
#                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                               parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #Cl~INFECTION + (SVL*Sex)
# Cl.infection.svl.sex <- data.frame(y = Cl.omit$Cl, 
#                                    x1 = Cl.omit$infection, x2 = Cl.omit$SVL,
#                                    x3 = Cl.omit$SEX)
# mod.Cl.infection.svl.sex=jags(model.file=model_gamma_3.2,data=Cl.infection.svl.sex,
#                               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                               parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #Cl~INFECTION + (SVL*Habitat)
# Cl.infection.svl.habitat <- data.frame(y = Cl.omit$Cl, 
#                                        x1 = Cl.omit$infection, x2 = Cl.omit$SVL,
#                                        x3 = Cl.omit$habitat)
# mod.Cl.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=Cl.infection.svl.habitat,
#                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                   parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #Cl~INFECTION + (Sex)+(Habitat)
# Cl.infection.sex.habitat <- data.frame(y = Cl.omit$Cl, 
#                                        x1 = Cl.omit$infection, x2 = Cl.omit$SEX,
#                                        x3 = Cl.omit$habitat)
# mod.Cl.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=Cl.infection.sex.habitat,
#                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                   parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Cl~INFECTION + (SVL*Sex)+(SVL*Habitat)
# Cl.infection.svl.sex.habitat <- data.frame(y = Cl.omit$Cl, 
#                                            x1 = Cl.omit$infection, x2 = Cl.omit$SVL,
#                                            x3 = Cl.omit$SEX, x4 = Cl.omit$habitat)
# mod.Cl.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=Cl.infection.svl.sex.habitat,
#                                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                       parameters.to.save=c("alpha","b1","b2","b3",
#                                                            "b4","b5","b6","LogLik"))
# 
# 
# ####iCa####
# iCa.omit <- filter(diagnosed, !is.na(iCa))
# #iCa~1
# iCa.null <- data.frame(y = iCa.omit$iCa, x1 = 0)
# mod.iCa.null=jags(model.file=model_gamma_0,data=iCa.null,
#                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                   parameters.to.save=c("alpha","LogLik"))
# 
# #iCa ~ infection
# iCa.infection <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$infection)
# mod.iCa.infection=jags(model.file=model_gamma_1,data=iCa.infection,
#                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                        parameters.to.save=c("alpha","b1","LogLik"))
# 
# # iCa ~ SVL
# iCa.svl <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL)
# mod.iCa.svl=jags(model.file=model_gamma_1,data=iCa.svl,
#                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                  parameters.to.save=c("alpha","b1","LogLik"))
# 
# # iCa ~ Sex
# iCa.sex <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SEX)
# mod.iCa.sex=jags(model.file=model_gamma_1,data=iCa.sex,
#                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                  parameters.to.save=c("alpha","b1","LogLik"))
# 
# # iCa ~ Habitat
# iCa.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$habitat)
# mod.iCa.habitat=jags(model.file=model_gamma_1,data=iCa.habitat,
#                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                      parameters.to.save=c("alpha","b1","LogLik"))
# 
# #iCa~(SVL*Sex)
# iCa.svl.sex <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, x2 = iCa.omit$SEX)
# mod.iCa.svl.sex=jags(model.file=model_gamma_2mult,data=iCa.svl.sex,
#                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                      parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #iCa~(SVL*Habitat)
# iCa.svl.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, x2 = iCa.omit$habitat)
# mod.iCa.svl.habitat=jags(model.file=model_gamma_2mult,data=iCa.svl.habitat,
#                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                          parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #iCa~(Sex)+(Habitat)
# iCa.sex.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SEX, x2 = iCa.omit$habitat)
# mod.iCa.sex.habitat=jags(model.file=model_gamma_2add,data=iCa.sex.habitat,
#                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                          parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #iCa~(SVL*Sex)+(SVL*Habitat)
# iCa.svl.sex.habitat <- data.frame(y = iCa.omit$iCa, x1 = iCa.omit$SVL, 
#                                   x2 = iCa.omit$SEX, x3 = iCa.omit$habitat)
# mod.iCa.svl.sex.habitat=jags(model.file=model_gamma_3,data=iCa.svl.sex.habitat,
#                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                              parameters.to.save=c("alpha","b1","b2","b3","b4","b5","LogLik"))
# 
# # iCa ~ INFECTION + SVL
# iCa.infection.svl <- data.frame(y = iCa.omit$iCa, 
#                                 x1 = iCa.omit$infection, x2 = iCa.omit$SVL)
# mod.iCa.infection.svl=jags(model.file=model_gamma_2add,data=iCa.infection.svl,
#                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                            parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # iCa ~ INFECTION + Sex
# iCa.infection.sex <- data.frame(y = iCa.omit$iCa, 
#                                 x1 = iCa.omit$infection, x2 = iCa.omit$SEX)
# mod.iCa.infection.sex=jags(model.file=model_gamma_2add,data=iCa.infection.sex,
#                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                            parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # iCa ~ INFECTION + Habitat
# iCa.infection.habitat <- data.frame(y = iCa.omit$iCa, 
#                                     x1 = iCa.omit$infection, x2 = iCa.omit$habitat)
# mod.iCa.infection.habitat=jags(model.file=model_gamma_2add,data=iCa.infection.habitat,
#                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #iCa~INFECTION + (SVL*Sex)
# iCa.infection.svl.sex <- data.frame(y = iCa.omit$iCa, 
#                                     x1 = iCa.omit$infection, x2 = iCa.omit$SVL,
#                                     x3 = iCa.omit$SEX)
# mod.iCa.infection.svl.sex=jags(model.file=model_gamma_3.2,data=iCa.infection.svl.sex,
#                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #iCa~INFECTION + (SVL*Habitat)
# iCa.infection.svl.habitat <- data.frame(y = iCa.omit$iCa, 
#                                         x1 = iCa.omit$infection, x2 = iCa.omit$SVL,
#                                         x3 = iCa.omit$habitat)
# mod.iCa.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=iCa.infection.svl.habitat,
#                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                    parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #iCa~INFECTION + (Sex)+(Habitat)
# iCa.infection.sex.habitat <- data.frame(y = iCa.omit$iCa, 
#                                         x1 = iCa.omit$infection, x2 = iCa.omit$SEX,
#                                         x3 = iCa.omit$habitat)
# mod.iCa.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=iCa.infection.sex.habitat,
#                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                    parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #iCa~INFECTION + (SVL*Sex)+(SVL*Habitat)
# iCa.infection.svl.sex.habitat <- data.frame(y = iCa.omit$iCa, 
#                                             x1 = iCa.omit$infection, x2 = iCa.omit$SVL,
#                                             x3 = iCa.omit$SEX, x4 = iCa.omit$habitat)
# mod.iCa.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=iCa.infection.svl.sex.habitat,
#                                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                        parameters.to.save=c("alpha","b1","b2","b3",
#                                                             "b4","b5","b6","LogLik"))
# 
# ####Glu####
# Glu.omit <- filter(diagnosed, !is.na(Glu))
# #Glu~1
# Glu.null <- data.frame(y = Glu.omit$Glu, x1 = 0)
# mod.Glu.null=jags(model.file=model_gamma_0,data=Glu.null,
#                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                   parameters.to.save=c("alpha","LogLik"))
# 
# #Glu ~ infection
# Glu.infection <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$infection)
# mod.Glu.infection=jags(model.file=model_gamma_1,data=Glu.infection,
#                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                        parameters.to.save=c("alpha","b1","LogLik"))
# 
# # Glu ~ SVL
# Glu.svl <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL)
# mod.Glu.svl=jags(model.file=model_gamma_1,data=Glu.svl,
#                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                  parameters.to.save=c("alpha","b1","LogLik"))
# 
# # Glu ~ Sex
# Glu.sex <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SEX)
# mod.Glu.sex=jags(model.file=model_gamma_1,data=Glu.sex,
#                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                  parameters.to.save=c("alpha","b1","LogLik"))
# 
# # Glu ~ Habitat
# Glu.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$habitat)
# mod.Glu.habitat=jags(model.file=model_gamma_1,data=Glu.habitat,
#                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                      parameters.to.save=c("alpha","b1","LogLik"))
# 
# #Glu~(SVL*Sex)
# Glu.svl.sex <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL, x2 = Glu.omit$SEX)
# mod.Glu.svl.sex=jags(model.file=model_gamma_2mult,data=Glu.svl.sex,
#                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                      parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Glu~(SVL*Habitat)
# Glu.svl.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL, x2 = Glu.omit$habitat)
# mod.Glu.svl.habitat=jags(model.file=model_gamma_2mult,data=Glu.svl.habitat,
#                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                          parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Glu~(Sex)+(Habitat)
# Glu.sex.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SEX, x2 = Glu.omit$habitat)
# mod.Glu.sex.habitat=jags(model.file=model_gamma_2add,data=Glu.sex.habitat,
#                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                          parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #Glu~(SVL*Sex)+(SVL*Habitat)
# Glu.svl.sex.habitat <- data.frame(y = Glu.omit$Glu, x1 = Glu.omit$SVL, 
#                                   x2 = Glu.omit$SEX, x3 = Glu.omit$habitat)
# mod.Glu.svl.sex.habitat=jags(model.file=model_gamma_3,data=Glu.svl.sex.habitat,
#                              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                              parameters.to.save=c("alpha","b1","b2","b3","b4","b5","LogLik"))
# 
# # Glu ~ INFECTION + SVL
# Glu.infection.svl <- data.frame(y = Glu.omit$Glu, 
#                                 x1 = Glu.omit$infection, x2 = Glu.omit$SVL)
# mod.Glu.infection.svl=jags(model.file=model_gamma_2add,data=Glu.infection.svl,
#                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                            parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # Glu ~ INFECTION + Sex
# Glu.infection.sex <- data.frame(y = Glu.omit$Glu, 
#                                 x1 = Glu.omit$infection, x2 = Glu.omit$SEX)
# mod.Glu.infection.sex=jags(model.file=model_gamma_2add,data=Glu.infection.sex,
#                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                            parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# # Glu ~ INFECTION + Habitat
# Glu.infection.habitat <- data.frame(y = Glu.omit$Glu, 
#                                     x1 = Glu.omit$infection, x2 = Glu.omit$habitat)
# mod.Glu.infection.habitat=jags(model.file=model_gamma_2add,data=Glu.infection.habitat,
#                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                parameters.to.save=c("alpha","b1","b2","LogLik"))
# 
# #Glu~INFECTION + (SVL*Sex)
# Glu.infection.svl.sex <- data.frame(y = Glu.omit$Glu, 
#                                     x1 = Glu.omit$infection, x2 = Glu.omit$SVL,
#                                     x3 = Glu.omit$SEX)
# mod.Glu.infection.svl.sex=jags(model.file=model_gamma_3.2,data=Glu.infection.svl.sex,
#                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #Glu~INFECTION + (SVL*Habitat)
# Glu.infection.svl.habitat <- data.frame(y = Glu.omit$Glu, 
#                                         x1 = Glu.omit$infection, x2 = Glu.omit$SVL,
#                                         x3 = Glu.omit$habitat)
# mod.Glu.infection.svl.habitat=jags(model.file=model_gamma_3.2,data=Glu.infection.svl.habitat,
#                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                    parameters.to.save=c("alpha","b1","b2","b3","b4","LogLik"))
# 
# #Glu~INFECTION + (Sex)+(Habitat)
# Glu.infection.sex.habitat <- data.frame(y = Glu.omit$Glu, 
#                                         x1 = Glu.omit$infection, x2 = Glu.omit$SEX,
#                                         x3 = Glu.omit$habitat)
# mod.Glu.infection.sex.habitat=jags(model.file=model_gamma_3.3,data=Glu.infection.sex.habitat,
#                                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                    parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
# 
# #Glu~INFECTION + (SVL*Sex)+(SVL*Habitat)
# Glu.infection.svl.sex.habitat <- data.frame(y = Glu.omit$Glu, 
#                                             x1 = Glu.omit$infection, x2 = Glu.omit$SVL,
#                                             x3 = Glu.omit$SEX, x4 = Glu.omit$habitat)
# mod.Glu.infection.svl.sex.habitat=jags(model.file=model_gamma_3.4,data=Glu.infection.svl.sex.habitat,
#                                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                                        parameters.to.save=c("alpha","b1","b2","b3",
#                                                             "b4","b5","b6","LogLik"))