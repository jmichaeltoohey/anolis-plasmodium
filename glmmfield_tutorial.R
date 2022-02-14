#The package glmmfields can be used to fit a variety of spatial GLM. 
#The method is described in detail in this paper:
# https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2403 
#We are going to use some data we collected in Puerto Rico on Anolis gundlachi to describe the method.

#DATA FORMATING####
library(dplyr)
data <- read.csv("datasheets/Spatial Heterogeneity.csv",header=TRUE)

#only use individuals for which we have all the data (no missing info)
data <- filter(data, !is.na(SVL))

#filter our the stratalus
data <- filter(data, Spp == "A. gundlachi")

#center SVL data to combat collinearity with Sex
data$SVL <- data$SVL-mean(data$SVL)

#time = 1 for simplicity, pt = id for each point, lon+lat = coordinates,
#y = response variable (could be any distribution); covariates x1. x2. x3
d <- data.frame(time=rep(1,dim(data)[1]),pt=1:dim(data)[1],lon=data$xcoord,
                lat=data$ycoord,station_id=1:dim(data)[1],y=data$Infeccion,
                svl=data$SVL,sex=data$Sex)

d <- d[complete.cases(d), ]

#visualize raw data spatially
library(ggplot2)
ggplot(d, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)
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
library(glmmfields)
m_spatial.0 <- glmmfields(y ~ 1,
                          data = d, family = binomial (link = "logit"),
                          lat = "lat", lon = "lon", nknots = 6,
                          iter = 100000, chains = 3,
                          prior_intercept = student_t(3, 0, 10),
                          prior_beta = student_t(3, 0, 3),
                          prior_sigma = half_t(3, 0, 3),
                          prior_gp_theta = half_t(3, 0, 10),
                          prior_gp_sigma = half_t(3, 0, 3),
                          seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m_spatial.0
#gp_sigma and gp_theta are parameters for the spatial process while B[1] is the intercept.
#n_eff is the effective sample size and Rhat is a measure of model convergence
# Usually we want this measure to be close to 1.00.

#Now we can fit a model with svl as a covariate and another with svl+sex.
m_spatial.svl <- glmmfields(y ~ svl,
                        data = d, family = binomial (link = "logit"),
                        lat = "lat", lon = "lon", nknots = 6,
                        iter = 100000, chains = 3,
                        prior_intercept = student_t(3, 0, 10),
                        prior_beta = student_t(3, 0, 3),
                        prior_sigma = half_t(3, 0, 3),
                        prior_gp_theta = half_t(3, 0, 10),
                        prior_gp_sigma = half_t(3, 0, 3),
                        seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)

m_spatial.sex <- glmmfields(y ~ sex,
                        data = d, family = binomial (link = "logit"),
                        lat = "lat", lon = "lon", nknots = 6, 
                        iter = 100000, chains = 3,
                        prior_intercept = student_t(3, 0, 10),
                        prior_beta = student_t(3, 0, 3),
                        prior_sigma = half_t(3, 0, 3),
                        prior_gp_theta = half_t(3, 0, 10),
                        prior_gp_sigma = half_t(3, 0, 3),
                        seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m_spatial.svl.sex <- glmmfields(y ~ svl*sex,
                        data = d, family = binomial (link = "logit"),
                        lat = "lat", lon = "lon", nknots = 6, 
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
library(loo)
ll.PR.0 <- extract_log_lik(m_spatial.0$model) #space
ll.PR.svl <- extract_log_lik(m_spatial.svl$model) #svl
ll.PR.sex <- extract_log_lik(m_spatial.sex$model) #sex
ll.PR.svl.sex <- extract_log_lik(m_spatial.svl.sex$model) #svl+sex

waic.PR.spatial.0 <- waic(ll.PR.0)
waic.PR.spatial.svl <- waic(ll.PR.svl)
waic.PR.spatial.sex <- waic(ll.PR.sex)
waic.PR.spatial.svl.sex <- waic(ll.PR.svl.sex)

loo.PR.spatial.null <- loo(ll.PR.0)
loo.PR.spatial.svl <- loo(ll.PR.svl)
loo.PR.spatial.sex <- loo(ll.PR.sex)
loo.PR.spatial.svl.sex <- loo(ll.PR.svl.sex)

ll.PR.infection.null <- mod.PR.infection.null$BUGSoutput$sims.list$LogLik
ll.PR.infection.svl <- mod.PR.infection.svl$BUGSoutput$sims.list$LogLik
ll.PR.infection.sex <- mod.PR.infection.sex$BUGSoutput$sims.list$LogLik
ll.PR.infection.svl.sex <- mod.PR.infection.svl.sex$BUGSoutput$sims.list$LogLik

waic.PR.infection.null <- waic(ll.PR.infection.null)
waic.PR.infection.svl <- waic(ll.PR.infection.svl)
waic.PR.infection.sex <- waic(ll.PR.infection.sex)
waic.PR.infection.svl.sex <- waic(ll.PR.infection.svl.sex)

loo.PR.infection.null <- loo(ll.PR.infection.null)
loo.PR.infection.svl <- loo(ll.PR.infection.svl)
loo.PR.infection.sex <- loo(ll.PR.infection.sex)
loo.PR.infection.svl.sex <- loo(ll.PR.infection.svl.sex)

PR.infection.list <- list(loo.PR.infection.null, loo.PR.infection.svl, loo.PR.infection.sex,
                          loo.PR.infection.svl.sex)
PR.infection.list <- list(waic.PR.infection.null, waic.PR.infection.svl, waic.PR.infection.sex,
                          waic.PR.infection.svl.sex)
PR.infection.compare = loo_compare(PR.infection.list)
PR.infection.compare

PR.space.list <- list(loo.PR.spatial.null,
                       loo.PR.spatial.svl,
                       loo.PR.spatial.sex,
                       loo.PR.spatial.svl.sex)
PR.space.list <- list(waic.PR.spatial.0,
                      waic.PR.spatial.svl,
                      waic.PR.spatial.sex,
                      waic.PR.spatial.svl.sex)
PR.space.compare = loo_compare(PR.space.list)
PR.space.compare

PR.full.list <- c(PR.infection.list,PR.space.list)
PR.full.compare = loo_compare(PR.full.list)
PR.full.compare

PR.infection_names <- c("p(infection) ~ SVL",
                        "p(infection) ~ Sex","p(infection) ~ (SVL*Sex)","p(infection) ~ 1")

rownames(PR.infection.compare) <- PR.infection_names
print(PR.infection.compare, simplify=FALSE)

PR.space_names <- c("p(infection) ~ Space+SVL",
                     "p(infection) ~ Space+(SVL*sex)",
                     "p(infection) ~ Space+Sex",
                     "p(infection) ~ Space")

rownames(PR.space.compare) <- PR.space_names
print(PR.space.compare, simplify=FALSE)

PR.full_names <- c("p(infection) ~ Space+SVL","p(infection) ~ Space+(SVL*Sex)",
                   "p(infection) ~ SVL","p(infection) ~ Space+Sex",
                   "p(infection) ~ Sex","p(infection) ~ (SVL*Sex)",
                   "p(infection) ~ Space","p(infection) ~ 1")
rownames(PR.full.compare) <- PR.full_names
print(PR.full.compare, simplify=FALSE)

weights.PR <- data.frame(length=8)
sums=0
for(i in 1:length(PR.full.list)){
  sums <- sums+ exp(0.5*(PR.full.compare[1,7]-PR.full.compare[i,7]))
}
for(i in 1:length(PR.full.list)){
  weights.PR[i] <- exp(0.5*(PR.full.compare[1,7]-PR.full.compare[i,7]))/sums
}
weights.PR

PR_SPACE_modsel_table <- as.data.frame(PR.full.compare)[7:8]
PR_SPACE_modsel_table$weight <- t(weights.PR)
print(PR_SPACE_modsel_table)
write.table(PR_SPACE_modsel_table, file="PR_SPACE_modsel_table.txt", sep=",", quote=FALSE,row.names=T)

print(m_spatial.svl)


#If you want to do LOOIC the glmmfields package has its own convenient function.
loo(m_spatial.0) #null
loo(m_spatial.svl) #svl
loo(m_spatial.sex) #svl+sex
loo(m_spatial.svl.sex) #svl+sex
#MAKING PREDICTIVE PLOTS####
#We will use the function predict in the same way we use it to make predictions
#from non-spatial linear models or GLM. 
#This is going to give us a probability of infection for each sample with their credible intervals. 
#We can also calculate the highest posterior density for each parameter.
p <- predict(m_spatial, type = "response")
head(p)
head(tidy(m_spatial, conf.int = TRUE, conf.method = "HPDinterval"))

#We can also use the function expand.grid to make a surface that includes
#the range of x and y coordiantes to make spatial predictions. 
#To do this we can add the covariate svl and specific predictions for each spatial
#point in the grid. Note that we are doing it in this example for mean value of svl. 
#You can also try it for large and small lizards to see if the predictions change.
pred_grid <- expand.grid(lat = seq(min(d$lat), max(d$lat), length.out = 30),
                         lon = seq(min(d$lon), max(d$lon), length.out = 30))
pred_grid$svl <- mean(d$svl)
pred_grid$prediction <- predict(m_spatial.svl, newdata = pred_grid,
                                type = "response")$estimate
head(pred_grid)

#Now we are ready to make a pretty plot of these spatial predictions.
ggplot(pred_grid, aes(lon, lat, fill = prediction)) +
  geom_raster() +
  viridis::scale_fill_viridis()+ggtitle("Spatial predictions p(Infection) of average-length 
                                        A. gundlachi at a site in Puerto Rico")

pred_grid$svl <- min(d$svl)
pred_grid$prediction <- predict(m_spatial, newdata = pred_grid,
                                type = "response")$estimate
head(pred_grid)

#Now we are ready to make a pretty plot of these spatial predictions.
ggplot(pred_grid, aes(lon, lat, fill = prediction)) +
  geom_raster() +
  viridis::scale_fill_viridis()+ggtitle("Spatial predictions p(Infection) of short-length 
                                        A. gundlachi at a site in Puerto Rico")

pred_grid$svl <- max(d$svl)
pred_grid$prediction <- predict(m_spatial, newdata = pred_grid,
                                type = "response")$estimate
head(pred_grid)

#Now we are ready to make a pretty plot of these spatial predictions.
ggplot(pred_grid, aes(lon, lat, fill = prediction)) +
  geom_raster() +
  viridis::scale_fill_viridis()+ggtitle("Spatial predictions p(Infection) of long-length 
                                        A. gundlachi at a site in Puerto Rico")


