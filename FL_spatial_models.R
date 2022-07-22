#The package glmmfields can be used to fit a variety of spatial GLM. 
#The method is described in detail in this paper:
# https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2403 

library(glmmfields)
library(loo)
library(R2jags)
set.seed(123)


#DATA FORMATING####
source("FL_data_carpentry.R")
EnvF <- filter(diagnosed, SITE == "Env. F")

#rotate by 45 degrees
degree <- pi * 45 / 180
l <- sqrt(EnvF$x_coord^2 + EnvF$y_coord^2)
teta <- atan(EnvF$y_coord / EnvF$x_coord)
EnvF$x_coord <- l * cos(teta - degree)
EnvF$y_coord <- l * sin(teta - degree)

#time = 1 for simplicity, pt = id for each point, lon+lat = coordinates,
#y = response variable (could be any distribution); covariates x1. x2. x3
d.EnvF <- data.frame(time=rep(1,dim(EnvF)[1]),pt=1:dim(EnvF)[1],lon=-EnvF$x_coord,
                     lat=EnvF$y_coord,station_id=1:dim(EnvF)[1],y=EnvF$infection,
                     svl=EnvF$SVL,sex=EnvF$SEX)

#only use individuals for which we have all the data (no missing info)
d.EnvF <- d.EnvF[complete.cases(d.EnvF), ] %>% na.omit()

#convert to 0-100 axes
d.EnvF$lon = (((d.EnvF$lon - min(d.EnvF$lon)) * (100 - 0)) / (max(d.EnvF$lon) - min(d.EnvF$lon))) + 0
d.EnvF$lat = (((d.EnvF$lat - min(d.EnvF$lat)) * (100 - 0)) / (max(d.EnvF$lat) - min(d.EnvF$lat))) + 0

#NON-SPATIAL MODELS####
EnvF <- filter(EnvF, !is.na(x_coord))
# p(infection) ~ 1
EnvF.infection.null <- data.frame(y = EnvF$infection, x1 = 0)
mod.EnvF.infection.null=jags(model.file=model_binom_0,data=EnvF.infection.null,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("alpha","LogLik"))

# p(infection) ~ SVL
EnvF.infection.svl <- data.frame(y = EnvF$infection, x1 = EnvF$SVL)
mod.EnvF.infection.svl=jags(model.file=model_binom_1,data=EnvF.infection.svl,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

# p(infection) ~ Sex
EnvF.infection.sex <- data.frame(y = EnvF$infection, x1 = EnvF$SEX)
mod.EnvF.infection.sex=jags(model.file=model_binom_1,data=EnvF.infection.sex,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))




#SPATIAL MODELS#####
#The models are Bayesian GLM that use a MVT spatial process. 
#These are fitted in rstan. This is an alternative to JAGS that also estimate
#posterior probabilities using MCMC.
#Access at the package documentation to see how we set up vague priors.
#We are using the argument save_log_lik = TRUE to save the log likelihood 
#so we can estimate WAIC and other measures for model selection.

# p(infection) ~ Space
m.EnvF_spatial.0 <- glmmfields(y ~ 1,
                               data = d.EnvF, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 20,
                               iter = 100000, chains = 4,
                               prior_intercept = student_t(3, 0, 10),
                               prior_beta = student_t(3, 0, 3),
                               prior_sigma = half_t(3, 0, 3),
                               prior_gp_theta = half_t(3, 0, 10),
                               prior_gp_sigma = half_t(3, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)

# p(infection) ~ Space + SVL
m.EnvF_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.EnvF, family = binomial (link = "logit"),
                                 time = NULL,lat = "lat", lon = "lon", nknots = 20,
                                 iter = 100000, chains = 4,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)

# p(infection) ~ Space + SEX
m.EnvF_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.EnvF, family = binomial (link = "logit"),
                                 time = NULL,lat = "lat", lon = "lon", nknots = 20,
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
ll.EnvF.spatial.null <- extract_log_lik(m.EnvF_spatial.0$model)
ll.EnvF.spatial.svl <- extract_log_lik(m.EnvF_spatial.svl$model)
ll.EnvF.spatial.sex <- extract_log_lik(m.EnvF_spatial.sex$model)
ll.EnvF.infection.null <- mod.EnvF.infection.null$BUGSoutput$sims.list$LogLik
ll.EnvF.infection.svl <- mod.EnvF.infection.svl$BUGSoutput$sims.list$LogLik
ll.EnvF.infection.sex <- mod.EnvF.infection.sex$BUGSoutput$sims.list$LogLik

#Extract wAIC for each model
waic.EnvF.spatial.0 <- waic(ll.EnvF.spatial.null)
waic.EnvF.spatial.svl <- waic(ll.EnvF.spatial.svl)
waic.EnvF.spatial.sex <- waic(ll.EnvF.spatial.sex)
waic.EnvF.infection.null <- waic(ll.EnvF.infection.null)
waic.EnvF.infection.svl <- waic(ll.EnvF.infection.svl)
waic.EnvF.infection.sex <- waic(ll.EnvF.infection.sex)

#Order wAIC in list to compare
EnvF.infection.list <- list(waic.EnvF.infection.null, waic.EnvF.infection.svl, waic.EnvF.infection.sex)
EnvF.infection.compare = loo_compare(EnvF.infection.list)

EnvF.space.list <- list(waic.EnvF.spatial.0,
                      waic.EnvF.spatial.svl,
                      waic.EnvF.spatial.sex)
EnvF.space.compare = loo_compare(EnvF.space.list)

EnvF.full.list <- c(waic.EnvF.infection.null, 
                    waic.EnvF.infection.svl, 
                    waic.EnvF.infection.sex,
                    waic.EnvF.spatial.0,
                    waic.EnvF.spatial.svl,
                    waic.EnvF.spatial.sex)
EnvF.full.compare = loo_compare(EnvF.full.list)

#Given names for clarity
#Done manually by comparing numbered output of loo_compare() to the order of models in the list
EnvF.infection_names <- c("p(infection) ~ SVL","p(infection) ~ 1",
                        "p(infection) ~ Sex")
rownames(EnvF.infection.compare) <- EnvF.infection_names
#print(EnvF.infection.compare, simplify=FALSE)

EnvF.space_names <- c("p(infection) ~ Space",
                      "p(infection) ~ Space+Sex",
                      "p(infection) ~ Space+SVL")
rownames(EnvF.space.compare) <- EnvF.space_names
#print(EnvF.space.compare, simplify=FALSE)

EnvF.full_names <- c("p(infection) ~ Space",
                     "p(infection) ~ Space+Sex",
                     "p(infection) ~ Space+SVL",
                     "p(infection) ~ SVL",
                     "p(infection) ~ 1",
                     "p(infection) ~ Sex")
rownames(EnvF.full.compare) <- EnvF.full_names
#print(EnvF.full.compare, simplify=FALSE)

#Calculate wAIC weights for each model in set
weights.EnvF <- data.frame(length=8)
sums=0
for(i in 1:length(EnvF.full.list)){
  sums <- sums+ exp(0.5*(EnvF.full.compare[1,7]-EnvF.full.compare[i,7]))
}
for(i in 1:length(EnvF.full.list)){
  weights.EnvF[i] <- exp(0.5*(EnvF.full.compare[1,7]-EnvF.full.compare[i,7]))/sums
}
weights.EnvF

#create and export model selection table
FL_SPACE_modsel_table <- as.data.frame(EnvF.full.compare)[7:8]
FL_SPACE_modsel_table$weight <- t(weights.EnvF)
print(FL_SPACE_modsel_table)
write.table(FL_SPACE_modsel_table, file="FL_SPACE_modsel_table.txt", sep=",", quote=FALSE,row.names=T)
