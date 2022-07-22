data <- read.csv("datasheets/Spatial Heterogeneity.csv",header=TRUE)

library("R2jags")
library(runjags)
library(dplyr)

####DATA CARPENTRY####
#only use individuals for which we have all the data (no missing info)
data <- filter(data, !is.na(SVL))

#filter our the stratalus
data <- filter(data, Spp == "A. gundlachi")

#filter our non-diagnosed
data <- filter(data, !is.na(Infeccion))

#change sex from character to binary
for(i in 1:nrow(data)){
  if(data$Sex[i] == "Female"){data$Sex[i] = 0}
  else if(data$Sex[i] == "Male"){data$Sex[i] = 1}}
data$Sex <- as.numeric(data$Sex)

prevalence.females <- (data %>% filter(Sex == '0') %>% filter(!is.na(Infeccion), Infeccion == 1) %>% nrow())/
  (data %>% filter(Sex == '0') %>% nrow())#.23 female infection rate
prevalence.males <- (data %>% filter(Sex == '1') %>% filter(!is.na(Infeccion), Infeccion == 1) %>% nrow())/
  (data %>% filter(Sex == '1') %>% nrow()) #.44 male infection rate

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





####MODELS####
#--------------------------Probability of Infection-------------------------####

# p(infection) ~ 1
PR.infection.null <- data.frame(y = d$y, x1 = 0)
mod.PR.infection.null=jags(model.file=model_binom_1,data=PR.infection.null,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","LogLik"))

# p(infection) ~ SVL
PR.infection.svl <- data.frame(y = d$y, x1 = d$svl)
mod.PR.infection.svl=jags(model.file=model_binom_1,data=PR.infection.svl,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("alpha","b1","LogLik"))

# p(infection) ~ Sex
PR.infection.sex <- data.frame(y = d$y, x1 = d$sex)
mod.PR.infection.sex=jags(model.file=model_binom_1,data=PR.infection.sex,
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          parameters.to.save=c("alpha","b1","LogLik"))

#-------------------------------Spatial Models------------------------------####
#The package glmmfields can be used to fit a variety of spatial GLM. 
#The method is described in detail in this paper:
# https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2403 
#We are going to use some data we collected in Puerto Rico on Anolis gundlachi to describe the method.

#time = 1 for simplicity, pt = id for each point, lon+lat = coordinates,
#y = response variable (could be any distribution); covariates x1. x2. x3
d <- data.frame(time=rep(1,dim(data)[1]),pt=1:dim(data)[1],lon=data$xcoord,
                lat=data$ycoord,station_id=1:dim(data)[1],y=data$Infeccion,
                svl=data$SVL,sex=data$Sex)

d <- d[complete.cases(d), ]

#MODEL BUILDING#
#The models are Bayesian GLM that use a MVT spatial process. 
#These are fitted in rstan. This is an alternative to JAGS that also estimate
#posterior probabilities using MCMC.
#Access at the package documentation to see how we set up vague priors.
#We are using the argument save_log_lik = TRUE to save the log likelihood 
#so we can estimate WAIC and other measures for model selection.

library(glmmfields)
# p(infection) ~ Space
m_spatial.0 <- glmmfields(y ~ 1,
                          data = d, family = binomial (link = "logit"),
                          lat = "lat", lon = "lon", nknots = 20,
                          iter = 100000, chains = 4,
                          prior_intercept = student_t(3, 0, 10),
                          prior_beta = student_t(3, 0, 3),
                          prior_sigma = half_t(3, 0, 3),
                          prior_gp_theta = half_t(3, 0, 10),
                          prior_gp_sigma = half_t(3, 0, 3),
                          seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)

# p(infection) ~ Space + SVL
m_spatial.svl <- glmmfields(y ~ svl,
                            data = d, family = binomial (link = "logit"),
                            lat = "lat", lon = "lon", nknots = 20,
                            iter = 100000, chains = 4,
                            prior_intercept = student_t(3, 0, 10),
                            prior_beta = student_t(3, 0, 3),
                            prior_sigma = half_t(3, 0, 3),
                            prior_gp_theta = half_t(3, 0, 10),
                            prior_gp_sigma = half_t(3, 0, 3),
                            seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)

# p(infection) ~ Space + SEX
m_spatial.sex <- glmmfields(y ~ sex,
                            data = d, family = binomial (link = "logit"),
                            lat = "lat", lon = "lon", nknots = 20, 
                            iter = 100000, chains = 4,
                            prior_intercept = student_t(3, 0, 10),
                            prior_beta = student_t(3, 0, 3),
                            prior_sigma = half_t(3, 0, 3),
                            prior_gp_theta = half_t(3, 0, 10),
                            prior_gp_sigma = half_t(3, 0, 3),
                            seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



####Model Selection####
library(AICcmodavg)
library(loo)

#Extract log likelihoods
ll.PR.spatial.null <- extract_log_lik(m_spatial.0$model)
ll.PR.spatial.svl <- extract_log_lik(m_spatial.svl$model)
ll.PR.spatial.sex <- extract_log_lik(m_spatial.sex$model)
ll.PR.infection.null <- mod.PR.infection.null$BUGSoutput$sims.list$LogLik
ll.PR.infection.svl <- mod.PR.infection.svl$BUGSoutput$sims.list$LogLik
ll.PR.infection.sex <- mod.PR.infection.sex$BUGSoutput$sims.list$LogLik

#Extract wAIC for each model
waic.PR.spatial.null <- waic(ll.PR.spatial.null)
waic.PR.spatial.svl <- waic(ll.PR.spatial.svl)
waic.PR.spatial.sex <- waic(ll.PR.spatial.sex)
waic.PR.infection.null <- waic(ll.PR.infection.null)
waic.PR.infection.svl <- waic(ll.PR.infection.svl)
waic.PR.infection.sex <- waic(ll.PR.infection.sex)

#Order wAIC in list to compare
PR.infection.list <- list(waic.PR.infection.null, 
                          waic.PR.infection.svl, 
                          waic.PR.infection.sex)
PR.infection.compare = loo_compare(PR.infection.list)

PR.space.list <- list(waic.PR.spatial.0,
                      waic.PR.spatial.svl,
                      waic.PR.spatial.sex)
PR.space.compare = loo_compare(PR.space.list)

PR.full.list <- list(waic.PR.infection.null, 
                     waic.PR.infection.svl, 
                     waic.PR.infection.sex,
                     waic.PR.spatial.0,
                     waic.PR.spatial.svl,
                     waic.PR.spatial.sex)
PR.full.compare = loo_compare(PR.full.list)

#Given names for clarity
#Done manually by comparing numbered output of loo_compare() to the order of models in the list
PR.infection_names <- c("p(infection) ~ SVL",
                        "p(infection) ~ 1",
                        "p(infection) ~ Sex")
rownames(PR.infection.compare) <- PR.infection_names
#print(PR.infection.compare, simplify=FALSE)

PR.space_names <- c("p(infection) ~ Space+SVL",
                    "p(infection) ~ Space+Sex",
                    "p(infection) ~ Space")
rownames(PR.space.compare) <- PR.space_names
#print(PR.space.compare, simplify=FALSE)

PR.full_names <- c("p(infection) ~ Space+SVL",
                   "p(infection) ~ SVL",
                   "p(infection) ~ Sex",
                   "p(infection) ~ Space+Sex",
                   "p(infection) ~ Space",
                   "p(infection) ~ 1")
rownames(PR.full.compare) <- PR.full_names
#print(PR.full.compare, simplify=FALSE)

#Calculate wAIC weights for each model in set
weights.PR <- data.frame(length=6)
sums=0
for(i in 1:length(PR.full.list)){
  sums <- sums+ exp(0.5*(PR.full.compare[1,7]-PR.full.compare[i,7]))
}
for(i in 1:length(PR.full.list)){
  weights.PR[i] <- exp(0.5*(PR.full.compare[1,7]-PR.full.compare[i,7]))/sums
}

#create and export model selection table
PR_SPACE_modsel_table <- as.data.frame(PR.full.compare)[7:8]
PR_SPACE_modsel_table$weight <- t(weights.PR)
#print(PR_SPACE_modsel_table)
write.table(PR_SPACE_modsel_table, file="PR_SPACE_modsel_table.txt", sep=",", quote=FALSE,row.names=T)


#-------------------------------Body Condition------------------------------####

###Filter for males with complete tails
PR.data.males <- filter(data, Sex == 1)
PR.completetails <- PR.data.males %>%filter(Tail=="Complete")

###Calculate body condition index
PR.male.mass.svl <- lm(Weight~SVL, data=PR.completetails, na.action=na.exclude)
PR.completetails$bodycond <- resid(PR.male.mass.svl)

# body condition ~ 1
PR.male.bodycond.null <- data.frame(y = PR.completetails$bodycond,   x1 = 0)
PR.mod.male.bodycond.null=jags(model.file=model_linear_0,data=PR.male.bodycond.null,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha","LogLik"))

# body condition ~ infection
PR.male.bodycond.infection <- data.frame(y = PR.completetails$bodycond,  x1 = PR.completetails$Infeccion)
PR.mod.male.bodycond.infection=jags(model.file=model_linear_1,data=PR.male.bodycond.infection,
                                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                 parameters.to.save=c("alpha","b1","LogLik"))

#Model Selection
loglik.PR.male.bodycond.null <- PR.mod.male.bodycond.null$BUGSoutput$sims.list$LogLik
loglik.PR.male.bodycond.infection <- PR.mod.male.bodycond.infection$BUGSoutput$sims.list$LogLik

waic.PR.male.bodycond.null <- waic(loglik.PR.male.bodycond.null)
waic.PR.male.bodycond.infection <- waic(loglik.PR.male.bodycond.infection)

PR.male.bodycond.list <- list(waic.PR.male.bodycond.null, waic.PR.male.bodycond.infection)
PR.male.bodycond.compare = loo_compare(PR.male.bodycond.list)

PR.male.bodycond_names <- c("M Body Condition ~ Infection",
                        "M Body Condition ~ 1")

rownames(PR.male.bodycond.compare) <- PR.male.bodycond_names
#print(PR.male.bodycond.compare, simplify=FALSE)

weights.PR.bodycond <- data.frame(length=8)
sums=0
for(i in 1:length(PR.male.bodycond.list)){
  sums <- sums+ exp(0.5*(PR.male.bodycond.compare[1,7]-PR.male.bodycond.compare[i,7]))
}
for(i in 1:length(PR.male.bodycond.list)){
  weights.PR.bodycond[i] <- exp(0.5*(PR.male.bodycond.compare[1,7]-PR.male.bodycond.compare[i,7]))/sums
}

PR.bodycond_modsel_table <- as.data.frame(PR.male.bodycond.compare)[7:8]
PR.bodycond_modsel_table$weight <- t(weights.PR.bodycond)
print(PR.bodycond_modsel_table)
write.table(PR.bodycond_modsel_table, file="PR.bodycond_mod.sel.txt", sep=",", quote=FALSE,row.names=T)

write.table(mcmcTab(PR.mod.male.bodycond.null, pars=c("alpha")),file="PR.mod.male.bodycond.null.txt", sep=",", quote=FALSE,row.names=F)
write.table(mcmcTab(PR.mod.male.bodycond.infection, pars=c("alpha","b1")),file="PR.mod.male.bodycond.infection.txt", sep=",", quote=FALSE,row.names=F)