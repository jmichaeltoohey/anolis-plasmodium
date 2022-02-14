#The package glmmfields can be used to fit a variety of spatial GLM. 
#The method is described in detail in this paper:
# https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2403 
#We are going to use some data we collected in Puerto Rico on Anolis gundlachi to describe the method.
library(ggplot2)
library(glmmfields)
library(loo)
library(R2jags)
set.seed(123)


#DATA FORMATING####
source("data_carpentry.R")
EnvF <- filter(diagnosed, SITE == "Env. F")

#rotate by 45 degrees
degree <- pi * 45 / 180
l <- sqrt(EnvF$lon^2 + EnvF$lat^2)
teta <- atan(EnvF$lat / EnvF$lon)
EnvF$lon <- l * cos(teta - degree)
EnvF$lat <- l * sin(teta - degree)

#time = 1 for simplicity, pt = id for each point, lon+lat = coordinates,
#y = response variable (could be any distribution); covariates x1. x2. x3
d.EnvF <- data.frame(time=rep(1,dim(EnvF)[1]),pt=1:dim(EnvF)[1],lon=-EnvF$lon,
                     lat=EnvF$lat,station_id=1:dim(EnvF)[1],y=EnvF$infection,
                     svl=EnvF$SVL,sex=EnvF$SEX)

#only use individuals for which we have all the data (no missing info)
d.EnvF <- d.EnvF[complete.cases(d.EnvF), ] %>% na.omit()

#convert to 0-100 axes
d.EnvF$lon = (((d.EnvF$lon - min(d.EnvF$lon)) * (100 - 0)) / (max(d.EnvF$lon) - min(d.EnvF$lon))) + 0
d.EnvF$lat = (((d.EnvF$lat - min(d.EnvF$lat)) * (100 - 0)) / (max(d.EnvF$lat) - min(d.EnvF$lat))) + 0



#visualize raw data spatially
ggplot(d.EnvF, aes(lon, lat, color = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)

#NON-SPATIAL MODELS####
EnvF <- filter(EnvF, !is.na(lon))
#p(infection)~1
EnvF.infection.null <- data.frame(y = EnvF$infection, n = nrow(EnvF), x1 = 0)
mod.EnvF.infection.null=jags(model.file=model_binom_1,data=EnvF.infection.null,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("alpha","LogLik"))

#p(infection)  ~ SVL
EnvF.infection.svl <- data.frame(y = EnvF$infection, n = nrow(EnvF), x1 = EnvF$SVL)
mod.EnvF.infection.svl=jags(model.file=model_binom_1,data=EnvF.infection.svl,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

#p(infection)  ~ Sex
EnvF.infection.sex <- data.frame(y = EnvF$infection, n = nrow(EnvF), x1 = EnvF$SEX)
mod.EnvF.infection.sex=jags(model.file=model_binom_1,data=EnvF.infection.sex,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))

#p(infection)~(SVL*Sex)
EnvF.infection.svl.sex <- data.frame(y = EnvF$infection, n = nrow(EnvF), 
                                x1 = EnvF$SVL, x2 = EnvF$SEX)
mod.EnvF.infection.svl.sex=jags(model.file=model_binom_2mult,data=EnvF.infection.svl.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","b3","LogLik"))


#MODEL BUILDING####
#The models are Bayesian GLM that use a MVT spatial process. 
#These are fitted in rstan. This is an alternative to JAGS that also estimate
#posterior probabilities using MCMC. Here we will use 3000 iterations
#and 1 chain for demonstration purposes but the documentation suggests 
#that we should use > 3000 iterations and at least 4 chains. 
#I would look at the package documentation to see how we set up vague priors. 
#We are using the argument save_log_lik = TRUE to save the log likelihood 
#so we can estimate WAIC and other measures for model selection.

#Let’s start with an intercept only model y∼1. 
#We are using the family binomial and the link logit because our response is infection. 
#For positive and continuous responses you can use a gamma family and a log link.
m.EnvF_spatial.0 <- glmmfields(y ~ 1,
                               data = d.EnvF, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 100000, chains = 3,
                               prior_intercept = student_t(3, 0, 10),
                               prior_beta = student_t(3, 0, 3),
                               prior_sigma = half_t(3, 0, 3),
                               prior_gp_theta = half_t(3, 0, 10),
                               prior_gp_sigma = half_t(3, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvF_spatial.0
#gp_sigma and gp_theta are parameters for the spatial process while B[1] is the intercept.
#n_eff is the effective sample size and Rhat is a measure of model convergence
# Usually we want this measure to be close to 1.00.

#Now we can fit a model with svl as a covariate and another with svl+sex.
m.EnvF_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.EnvF, family = binomial (link = "logit"),
                                 time = NULL,lat = "lat", lon = "lon", nknots = 6,
                                 iter = 21000, chains = 3,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvF_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.EnvF, family = binomial (link = "logit"),
                                 time = NULL,lat = "lat", lon = "lon", nknots = 6,
                                 iter = 100000, chains = 3,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvF_spatial.svl.sex <- glmmfields(y ~ svl*sex,
                                     data = d.EnvF, family = binomial (link = "logit"),
                                     time = NULL,lat = "lat", lon = "lon", nknots = 6, 
                                     iter = 100000, chains = 3,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)

#MODEL SELECTION####
#There are a few tools available to conduct model selection in a Bayesian framework using rstan. 
#There are relatively easy functions to estimate WAIC and LOOIC (leave-one-out information criterion). 
#This is a really good reference for model selection in a Bayesian framework:
# https://www.stat.colostate.edu/~hooten/papers/pdf/Hooten_Hobbs_EcolMono_2015.pdf.

#WAIC can be estimated using functions from the loo package. 
#First, we have to extract the log-likelihoods from the model ouputs.

#If you want to do LOOIC the glmmfields package has its own convenient function.
c(sum(EnvF$infection),nrow(EnvF))
sum(EnvF$SEX)

ll.EnvF.0 <- extract_log_lik(m.EnvF_spatial.0$model) #null
ll.EnvF.svl <- extract_log_lik(m.EnvF_spatial.svl$model) #svl
ll.EnvF.sex <- extract_log_lik(m.EnvF_spatial.sex$model) #svl+sex
ll.EnvF.svl.sex <- extract_log_lik(m.EnvF_spatial.svl.sex$model) #svl+sex

waic.EnvF.spatial.0 <- waic(ll.EnvF.0)
waic.EnvF.spatial.svl <- waic(ll.EnvF.svl)
waic.EnvF.spatial.sex <- waic(ll.EnvF.sex)
waic.EnvF.spatial.svl.sex <- waic(ll.EnvF.svl.sex)

#loo(m.EnvF_spatial.0) #null
#loo(m.EnvF_spatial.svl) #svl
#loo(m.EnvF_spatial.sex) #svl
#loo(m.EnvF_spatial.svl.sex) #svl+sex

ll.EnvF.infection.null <- mod.EnvF.infection.null$BUGSoutput$sims.list$LogLik
ll.EnvF.infection.svl <- mod.EnvF.infection.svl$BUGSoutput$sims.list$LogLik
ll.EnvF.infection.sex <- mod.EnvF.infection.sex$BUGSoutput$sims.list$LogLik
ll.EnvF.infection.svl.sex <- mod.EnvF.infection.svl.sex$BUGSoutput$sims.list$LogLik

waic.EnvF.infection.null <- waic(ll.EnvF.infection.null)
waic.EnvF.infection.svl <- waic(ll.EnvF.infection.svl)
waic.EnvF.infection.sex <- waic(ll.EnvF.infection.sex)
waic.EnvF.infection.svl.sex <- waic(ll.EnvF.infection.svl.sex)

EnvF.infection.list <- list(waic.EnvF.infection.null, waic.EnvF.infection.svl, waic.EnvF.infection.sex,
                          waic.EnvF.infection.svl.sex)
EnvF.infection.compare = loo_compare(EnvF.infection.list)
EnvF.infection.compare

EnvF.space.list <- list(waic.EnvF.spatial.0,
                      waic.EnvF.spatial.svl,
                      waic.EnvF.spatial.sex,
                      waic.EnvF.spatial.svl.sex)
EnvF.space.compare = loo_compare(EnvF.space.list)
EnvF.space.compare

EnvF.full.list <- c(EnvF.infection.list,EnvF.space.list)
EnvF.full.compare = loo_compare(EnvF.full.list)
EnvF.full.compare

EnvF.infection_names <- c("p(infection) ~ SVL","p(infection) ~ 1",
                        "p(infection) ~ Sex",
                        "p(infection) ~ (SVL*Sex)")

rownames(EnvF.infection.compare) <- EnvF.infection_names
print(EnvF.infection.compare, simplify=FALSE)

EnvF.space_names <- c("p(infection) ~ Space+Sex",
                    "p(infection) ~ Space",
                    "p(infection) ~ Space+SVL",
                    "p(infection) ~ Space+(SVL*Sex)")

rownames(EnvF.space.compare) <- EnvF.space_names
print(EnvF.space.compare, simplify=FALSE)

EnvF.full_names <- c(EnvF.space_names,EnvF.infection_names)
rownames(EnvF.full.compare) <- EnvF.full_names
print(EnvF.full.compare, simplify=FALSE)

weights.EnvF <- data.frame(length=8)
sums=0
for(i in 1:length(EnvF.full.list)){
  sums <- sums+ exp(0.5*(EnvF.full.compare[1,7]-EnvF.full.compare[i,7]))
}
for(i in 1:length(EnvF.full.list)){
  weights.EnvF[i] <- exp(0.5*(EnvF.full.compare[1,7]-EnvF.full.compare[i,7]))/sums
}
weights.EnvF

FL_SPACE_modsel_table <- as.data.frame(EnvF.full.compare)[7:8]
FL_SPACE_modsel_table$weight <- t(weights.EnvF)
print(FL_SPACE_modsel_table)
write.table(FL_SPACE_modsel_table, file="FL_SPACE_modsel_table.txt", sep=",", quote=FALSE,row.names=T)


weights.best.EnvF <- data.frame(length=2) #three strong models
sums=0
for(i in 1:2){
  sums <- sums+ exp(0.5*(EnvF.full.compare[1,7]-EnvF.full.compare[i,7]))
}
for(i in 1:2){
  weights.best.EnvF[i] <- exp(0.5*(EnvF.full.compare[1,7]-EnvF.full.compare[i,7]))/sums
}
weights.best.EnvF

library(BayesPostEst)
mcmcTab(m.EnvF_spatial.sex, Pr = TRUE)
mcmcTab(m.EnvF_spatial.0, Pr = TRUE)

#MAKING PREDICTIVE PLOTS####
#We will use the function predict in the same way we use it to make predictions
#from non-spatial linear models or GLM. 
#This is going to give us a probability of infection for each sample with their credible intervals. 
#We can also calculate the highest posterior density for each parameter.
p <- predict(m.EnvF_spatial.0, type = "response")
head(p)
head(tidy(m.EnvF_spatial.0, conf.int = TRUE, conf.method = "HPDinterval"))

#We can also use the function expand.grid to make a surface that includes
#the range of x and y coordiantes to make spatial predictions. 
#To do this we can add the covariate svl and specific predictions for each spatial
#point in the grid. Note that we are doing it in this example for mean value of svl. 
#You can also try it for large and small lizards to see if the predictions change.
pred_grid_FL <- expand.grid(lat = seq(min(d.EnvF$lat), max(d.EnvF$lat), length.out = 50),
                         lon = seq(min(d.EnvF$lon), max(d.EnvF$lon), length.out = 50))
pred_grid$svl <- max(d.EnvF$svl)
pred_grid$sex <- max(d.EnvF$sex)
#pred_grid$sex1 <- pred_grid$svl*pred_grid$sex
pred_grid$prediction <- predict(m.EnvF_spatial.sex, newdata = pred_grid,
                                type = "response")$estimate
head(pred_grid)

#Now we are ready to make a pretty plot of these spatial predictions.
pred_grid$sex <- max(d.EnvF$sex)
pred_grid_FL$prediction <- predict(m.EnvF_spatial.sex, newdata = pred_grid,
                                type = "response")$estimate
ggplot(pred_grid_FL, aes(lon, lat, fill = prediction)) +
  geom_raster() +
  viridis::scale_fill_viridis()
ggplot(d.EnvF, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)

library(BayesPostEst)
devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")
write.table(mcmctab(m.EnvF_spatial.sex),file="m.EnvF_spatial.sex.txt", sep=",", quote=FALSE,row.names=F)
write.table(mcmctab(m.EnvF_spatial.0, pars=c("gp_sigma","gp_theta","B[1]","lp__")),file="m.EnvF_spatial.0.txt", sep=",", quote=FALSE,row.names=F)
write.table(mcmctab(m.EnvF_spatial.svl.sex, pars=c("gp_sigma","gp_theta","B[1]","B[2]","B[3]","B[4]","lp__")),file="m.EnvF_spatial.svl.sex.txt", sep=",", quote=FALSE,row.names=F)

print(m.EnvF_spatial.sex)
print(m.EnvF_spatial.0)
print(m.EnvF_spatial.svl.sex)