#The package glmmfields can be used to fit a variety of spatial GLM. 
#The method is described in detail in this paper:
# https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2403 
#We are going to use some data we collected in Puerto Rico on Anolis gundlachi to describe the method.
library(ggplot2)
library(glmmfields)
library(loo)

#DATA FORMATING####
source("data_carpentry.R")
LawU <- filter(diagnosed, SITE == "Law-U")
LawF <- filter(diagnosed, SITE == "Law-F")
Museum <- filter(diagnosed, SITE == "Museum")
McCarty <- filter(diagnosed, SITE == "McCarty")
Ficke <- filter(diagnosed, SITE == "Ficke")
Alice <- filter(diagnosed, SITE == "Alice")
Entom <- filter(diagnosed, SITE == "Entom")
NATL <- filter(diagnosed, SITE == "NATL")
EnvU <- filter(diagnosed, SITE == "Env. U")
EnvF <- filter(diagnosed, SITE == "Env. F")
BTPU <- filter(diagnosed, SITE == "BTP U")
BTPF <- filter(diagnosed, SITE == "BTP F")
KewU <- filter(diagnosed, SITE == "Kew U")
KewF <- filter(diagnosed, SITE == "Kew F")

Gaines <- filter(diagnosed, y_coord > 29)
Orl <- filter(diagnosed, y_coord < 29)

#time = 1 for simplicity, pt = id for each point, lon+lat = coordinates,
#y = response variable (could be any distribution); covariates x1. x2. x3
d.LawU <- data.frame(time=rep(1,dim(LawU)[1]),pt=1:dim(LawU)[1],lon=-LawU$x_coord,
                lat=LawU$y_coord,station_id=1:dim(LawU)[1],y=LawU$infection,
                svl=LawU$SVL,sex=LawU$SEX)
#only use individuals for which we have all the data (no missing info)
d.LawU <- d.LawU[complete.cases(d.LawU), ] %>% na.omit()

#visualize raw data spatially
ggplot(d.LawU, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.LawF <- data.frame(time=rep(1,dim(LawF)[1]),pt=1:dim(LawF)[1],lon=-LawF$x_coord,
                     lat=LawF$y_coord,station_id=1:dim(LawF)[1],y=LawF$infection,
                     svl=LawF$SVL,sex=LawF$SEX)
d.LawF <- d.LawF[complete.cases(d.LawF), ] %>% na.omit()
ggplot(d.LawF, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.Museum <- data.frame(time=rep(1,dim(Museum)[1]),pt=1:dim(Museum)[1],lon=-Museum$x_coord,
                     lat=Museum$y_coord,station_id=1:dim(Museum)[1],y=Museum$infection,
                     svl=Museum$SVL,sex=Museum$SEX)
d.Museum <- d.Museum[complete.cases(d.Museum), ] %>% na.omit()
ggplot(d.Museum, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.McCarty <- data.frame(time=rep(1,dim(McCarty)[1]),pt=1:dim(McCarty)[1],lon=-McCarty$x_coord,
                     lat=McCarty$y_coord,station_id=1:dim(McCarty)[1],y=McCarty$infection,
                     svl=McCarty$SVL,sex=McCarty$SEX)
d.McCarty <- d.McCarty[complete.cases(d.McCarty), ] %>% na.omit()
ggplot(d.McCarty, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.Ficke <- data.frame(time=rep(1,dim(Ficke)[1]),pt=1:dim(Ficke)[1],lon=-Ficke$x_coord,
                     lat=Ficke$y_coord,station_id=1:dim(Ficke)[1],y=Ficke$infection,
                     svl=Ficke$SVL,sex=Ficke$SEX)
d.Ficke <- d.Ficke[complete.cases(d.Ficke), ] %>% na.omit()
ggplot(d.Ficke, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.Alice <- data.frame(time=rep(1,dim(Alice)[1]),pt=1:dim(Alice)[1],lon=-Alice$x_coord,
                     lat=Alice$y_coord,station_id=1:dim(Alice)[1],y=Alice$infection,
                     svl=Alice$SVL,sex=Alice$SEX)
d.Alice <- d.Alice[complete.cases(d.Alice), ] %>% na.omit()
ggplot(d.Alice, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.Entom <- data.frame(time=rep(1,dim(Entom)[1]),pt=1:dim(Entom)[1],lon=-Entom$x_coord,
                     lat=Entom$y_coord,station_id=1:dim(Entom)[1],y=Entom$infection,
                     svl=Entom$SVL,sex=Entom$SEX)
d.Entom <- d.Entom[complete.cases(d.Entom), ] %>% na.omit()
ggplot(d.Entom, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.NATL <- data.frame(time=rep(1,dim(NATL)[1]),pt=1:dim(NATL)[1],lon=-NATL$x_coord,
                     lat=NATL$y_coord,station_id=1:dim(NATL)[1],y=NATL$infection,
                     svl=NATL$SVL,sex=NATL$SEX)
d.NATL <- d.NATL[complete.cases(d.NATL), ] %>% na.omit()
ggplot(d.NATL, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.EnvU <- data.frame(time=rep(1,dim(EnvU)[1]),pt=1:dim(EnvU)[1],lon=-EnvU$x_coord,
                     lat=EnvU$y_coord,station_id=1:dim(EnvU)[1],y=EnvU$infection,
                     svl=EnvU$SVL,sex=EnvU$SEX)
d.EnvU <- d.EnvU[complete.cases(d.EnvU), ] %>% na.omit()
ggplot(d.EnvU, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.EnvF <- data.frame(time=rep(1,dim(EnvF)[1]),pt=1:dim(EnvF)[1],lon=-EnvF$x_coord,
                     lat=EnvF$y_coord,station_id=1:dim(EnvF)[1],y=EnvF$infection,
                     svl=EnvF$SVL,sex=EnvF$SEX)
d.EnvF <- d.EnvF[complete.cases(d.EnvF), ] %>% na.omit()
ggplot(d.EnvF, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.BTPU <- data.frame(time=rep(1,dim(BTPU)[1]),pt=1:dim(BTPU)[1],lon=-BTPU$x_coord,
                     lat=BTPU$y_coord,station_id=1:dim(BTPU)[1],y=BTPU$infection,
                     svl=BTPU$SVL,sex=BTPU$SEX)
d.BTPU <- d.BTPU[complete.cases(d.BTPU), ] %>% na.omit()
ggplot(d.BTPU, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.BTPF <- data.frame(time=rep(1,dim(BTPF)[1]),pt=1:dim(BTPF)[1],lon=-BTPF$x_coord,
                     lat=BTPF$y_coord,station_id=1:dim(BTPF)[1],y=BTPF$infection,
                     svl=BTPF$SVL,sex=BTPF$SEX)
d.BTPF <- d.BTPF[complete.cases(d.BTPF), ] %>% na.omit()
ggplot(d.BTPF, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.KewU <- data.frame(time=rep(1,dim(KewU)[1]),pt=1:dim(KewU)[1],lon=-KewU$x_coord,
                     lat=KewU$y_coord,station_id=1:dim(KewU)[1],y=KewU$infection,
                     svl=KewU$SVL,sex=KewU$SEX)
d.KewU <- d.KewU[complete.cases(d.KewU), ] %>% na.omit()
ggplot(d.KewU, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)


d.KewF <- data.frame(time=rep(1,dim(KewF)[1]),pt=1:dim(KewF)[1],lon=-KewF$x_coord,
                     lat=KewF$y_coord,station_id=1:dim(KewF)[1],y=KewF$infection,
                     svl=KewF$SVL,sex=KewF$SEX)
d.KewF <- d.KewF[complete.cases(d.KewF), ] %>% na.omit()
ggplot(d.KewF, aes(lon, lat, colour = y)) +
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
m.LawU_spatial.0 <- glmmfields(y ~ 1,
                          data = d.LawU, family = binomial (link = "logit"),
                          time = NULL, lat = "lat", lon = "lon", nknots = 4,
                          iter = 3000, chains = 1,
                          prior_intercept = student_t(100, 0, 10),
                          prior_beta = student_t(100, 0, 3),
                          prior_sigma = half_t(100, 0, 3),
                          prior_gp_theta = half_t(100, 0, 10),
                          prior_gp_sigma = half_t(100, 0, 3),
                          seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.LawU_spatial.0
#gp_sigma and gp_theta are parameters for the spatial process while B[1] is the intercept.
#n_eff is the effective sample size and Rhat is a measure of model convergence
# Usually we want this measure to be close to 1.00.

#Now we can fit a model with svl as a covariate and another with svl+sex.
m.LawU_spatial.svl <- glmmfields(y ~ svl,
                        data = d.LawU, family = binomial (link = "logit"),
                        lat = "lat", lon = "lon", nknots = 4,
                        iter = 3000, chains = 1,
                        prior_intercept = student_t(3, 0, 10),
                        prior_beta = student_t(3, 0, 3),
                        prior_sigma = half_t(3, 0, 3),
                        prior_gp_theta = half_t(3, 0, 10),
                        prior_gp_sigma = half_t(3, 0, 3),
                        seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.LawU_spatial.sex <- glmmfields(y ~ sex,
                            data = d.LawU, family = binomial (link = "logit"),
                            lat = "lat", lon = "lon", nknots = 4,
                            iter = 3000, chains = 1,
                            prior_intercept = student_t(3, 0, 10),
                            prior_beta = student_t(3, 0, 3),
                            prior_sigma = half_t(3, 0, 3),
                            prior_gp_theta = half_t(3, 0, 10),
                            prior_gp_sigma = half_t(3, 0, 3),
                            seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.LawU_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                         data = d.LawU, family = binomial (link = "logit"),
                         lat = "lat", lon = "lon", nknots = 4, iter = 3000, chains = 1,
                         prior_intercept = student_t(3, 0, 10),
                         prior_beta = student_t(3, 0, 3),
                         prior_sigma = half_t(3, 0, 3),
                         prior_gp_theta = half_t(3, 0, 10),
                         prior_gp_sigma = half_t(3, 0, 3),
                         seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.LawF_spatial.0 <- glmmfields(y ~ 1,
                               data = d.LawF, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.LawF_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.LawF, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.LawF_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.LawF, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.LawF_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.LawF, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.Museum_spatial.0 <- glmmfields(y ~ 1,
                               data = d.Museum, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Museum_spatial.sex <- glmmfields(y ~ svl,
                                 data = d.Museum, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Museum_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.Museum, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Museum_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.Museum, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.McCarty_spatial.0 <- glmmfields(y ~ 1,
                               data = d.McCarty, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.McCarty_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.McCarty, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.McCarty_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.McCarty, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.McCarty_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.McCarty, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.Ficke_spatial.0 <- glmmfields(y ~ 1,
                               data = d.Ficke, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Ficke_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.Ficke, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Ficke_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.Ficke, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Ficke_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.Ficke, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.Alice_spatial.0 <- glmmfields(y ~ 1,
                               data = d.Alice, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Alice_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.Alice, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Alice_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.Alice, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Alice_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.Alice, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.Entom_spatial.0 <- glmmfields(y ~ 1,
                               data = d.Entom, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Entom_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.Entom, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Entom_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.Entom, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.Entom_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.Entom, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.NATL_spatial.0 <- glmmfields(y ~ 1,
                               data = d.NATL, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.NATL_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.NATL, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.NATL_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.NATL, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.NATL_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.NATL, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.EnvU_spatial.0 <- glmmfields(y ~ 1,
                               data = d.EnvU, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvU_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.EnvU, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvU_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.EnvU, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvU_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.EnvU, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.EnvF_spatial.0 <- glmmfields(y ~ 1,
                               data = d.EnvF, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvF_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.EnvF, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvF_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.EnvF, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.EnvF_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.EnvF, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.BTPU_spatial.0 <- glmmfields(y ~ 1,
                               data = d.BTPU, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.BTPU_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.BTPU, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.BTPU_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.BTPU, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.BTPU_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.BTPU, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.BTPF_spatial.0 <- glmmfields(y ~ 1,
                               data = d.BTPF, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.BTPF_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.BTPF, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.BTPF_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.BTPF, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.BTPF_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.BTPF, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.KewU_spatial.0 <- glmmfields(y ~ 1,
                               data = d.KewU, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.KewU_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.KewU, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.KewU_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.KewU, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.KewU_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.KewU, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
                                     prior_intercept = student_t(3, 0, 10),
                                     prior_beta = student_t(3, 0, 3),
                                     prior_sigma = half_t(3, 0, 3),
                                     prior_gp_theta = half_t(3, 0, 10),
                                     prior_gp_sigma = half_t(3, 0, 3),
                                     seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)



m.KewF_spatial.0 <- glmmfields(y ~ 1,
                               data = d.KewF, family = binomial (link = "logit"),
                               time = NULL, lat = "lat", lon = "lon", nknots = 6,
                               iter = 3000, chains = 1,
                               prior_intercept = student_t(100, 0, 10),
                               prior_beta = student_t(100, 0, 3),
                               prior_sigma = half_t(100, 0, 3),
                               prior_gp_theta = half_t(100, 0, 10),
                               prior_gp_sigma = half_t(100, 0, 3),
                               seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.KewF_spatial.svl <- glmmfields(y ~ svl,
                                 data = d.KewF, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.KewF_spatial.sex <- glmmfields(y ~ sex,
                                 data = d.KewF, family = binomial (link = "logit"),
                                 lat = "lat", lon = "lon", nknots = 6,
                                 iter = 3000, chains = 1,
                                 prior_intercept = student_t(3, 0, 10),
                                 prior_beta = student_t(3, 0, 3),
                                 prior_sigma = half_t(3, 0, 3),
                                 prior_gp_theta = half_t(3, 0, 10),
                                 prior_gp_sigma = half_t(3, 0, 3),
                                 seed = 123,save_log_lik = TRUE # passed to rstan::sampling()
)
m.KewF_spatial.svl.sex <- glmmfields(y ~ svl+sex,
                                     data = d.KewF, family = binomial (link = "logit"),
                                     lat = "lat", lon = "lon", nknots = 6, iter = 3000, chains = 1,
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

#infection by site
c(sum(LawU$infection),nrow(LawU))
c(sum(LawF$infection),nrow(LawF))
c(sum(Museum$infection),nrow(Museum))
c(sum(McCarty$infection),nrow(McCarty))
c(sum(Ficke$infection),nrow(Ficke))
c(sum(Alice$infection),nrow(Alice))
c(sum(Entom$infection),nrow(Entom))
c(sum(NATL$infection),nrow(NATL))
c(sum(EnvU$infection),nrow(EnvU))
c(sum(EnvF$infection),nrow(EnvF))
c(sum(BTPU$infection),nrow(BTPU))
c(sum(BTPF$infection),nrow(BTPF))
c(sum(KewU$infection),nrow(KewU))
c(sum(KewF$infection),nrow(KewF))



ll.LawU.0 <- extract_log_lik(m.LawU_spatial.0$model) #null
ll.LawU.svl <- extract_log_lik(m.LawU_spatial.svl$model) #svl
ll.LawU.sex <- extract_log_lik(m.LawU_spatial.sex$model) #svl+sex
ll.LawU.svl.sex <- extract_log_lik(m.LawU_spatial.svl.sex$model) #svl+sex

waic(ll.LawU.0)
waic(ll.LawU.svl)
waic(ll.LawU.sex)
waic(ll.LawU.svl.sex)

loo(m.LawU_spatial.0) #null
loo(m.LawU_spatial.svl) #svl
loo(m.LawU_spatial.sex) #svl
loo(m.LawU_spatial.svl.sex) #svl+sex




ll.LawF.0 <- extract_log_lik(m.LawF_spatial.0$model) #null
ll.LawF.svl <- extract_log_lik(m.LawF_spatial.svl$model) #svl
ll.LawF.sex <- extract_log_lik(m.LawF_spatial.sex$model) #svl+sex
ll.LawF.svl.sex <- extract_log_lik(m.LawF_spatial.svl.sex$model) #svl+sex

waic(ll.LawF.0)
waic(ll.LawF.svl)
waic(ll.LawF.sex)
waic(ll.LawF.svl.sex)

loo(m.LawF_spatial.0) #null
loo(m.LawF_spatial.svl) #svl
loo(m.LawF_spatial.sex) #svl
loo(m.LawF_spatial.svl.sex) #svl+sex




ll.Museum.0 <- extract_log_lik(m.Museum_spatial.0$model) #null
ll.Museum.svl <- extract_log_lik(m.Museum_spatial.svl$model) #svl
ll.Museum.sex <- extract_log_lik(m.Museum_spatial.sex$model) #svl+sex
ll.Museum.svl.sex <- extract_log_lik(m.Museum_spatial.svl.sex$model) #svl+sex

waic(ll.Museum.0)
waic(ll.Museum.svl)
waic(ll.Museum.sex)
waic(ll.Museum.svl.sex)

loo(m.Museum_spatial.0) #null
loo(m.Museum_spatial.svl) #svl
loo(m.Museum_spatial.sex) #svl
loo(m.Museum_spatial.svl.sex) #svl+sex




ll.McCarty.0 <- extract_log_lik(m.McCarty_spatial.0$model) #null
ll.McCarty.svl <- extract_log_lik(m.McCarty_spatial.svl$model) #svl
ll.McCarty.sex <- extract_log_lik(m.McCarty_spatial.sex$model) #svl+sex
ll.McCarty.svl.sex <- extract_log_lik(m.McCarty_spatial.svl.sex$model) #svl+sex

waic(ll.McCarty.0)
waic(ll.McCarty.svl)
waic(ll.McCarty.sex)
waic(ll.McCarty.svl.sex)

loo(m.McCarty_spatial.0) #null
loo(m.McCarty_spatial.svl) #svl
loo(m.McCarty_spatial.sex) #svl
loo(m.McCarty_spatial.svl.sex) #svl+sex




ll.Ficke.0 <- extract_log_lik(m.Ficke_spatial.0$model) #null
ll.Ficke.svl <- extract_log_lik(m.Ficke_spatial.svl$model) #svl
ll.Ficke.sex <- extract_log_lik(m.Ficke_spatial.sex$model) #svl+sex
ll.Ficke.svl.sex <- extract_log_lik(m.Ficke_spatial.svl.sex$model) #svl+sex

waic(ll.Ficke.0)
waic(ll.Ficke.svl)
waic(ll.Ficke.sex)
waic(ll.Ficke.svl.sex)

loo(m.Ficke_spatial.0) #null
loo(m.Ficke_spatial.svl) #svl
loo(m.Ficke_spatial.sex) #svl
loo(m.Ficke_spatial.svl.sex) #svl+sex




ll.Alice.0 <- extract_log_lik(m.Alice_spatial.0$model) #null
ll.Alice.svl <- extract_log_lik(m.Alice_spatial.svl$model) #svl
ll.Alice.sex <- extract_log_lik(m.Alice_spatial.sex$model) #svl+sex
ll.Alice.svl.sex <- extract_log_lik(m.Alice_spatial.svl.sex$model) #svl+sex

waic(ll.Alice.0)
waic(ll.Alice.svl)
waic(ll.Alice.sex)
waic(ll.Alice.svl.sex)

loo(m.Alice_spatial.0) #null
loo(m.Alice_spatial.svl) #svl
loo(m.Alice_spatial.sex) #svl
loo(m.Alice_spatial.svl.sex) #svl+sex




ll.Entom.0 <- extract_log_lik(m.Entom_spatial.0$model) #null
ll.Entom.svl <- extract_log_lik(m.Entom_spatial.svl$model) #svl
ll.Entom.sex <- extract_log_lik(m.Entom_spatial.sex$model) #svl+sex
ll.Entom.svl.sex <- extract_log_lik(m.Entom_spatial.svl.sex$model) #svl+sex

waic(ll.Entom.0)
waic(ll.Entom.svl)
waic(ll.Entom.sex)
waic(ll.Entom.svl.sex)

loo(m.Entom_spatial.0) #null
loo(m.Entom_spatial.svl) #svl
loo(m.Entom_spatial.sex) #svl
loo(m.Entom_spatial.svl.sex) #svl+sex




ll.NATL.0 <- extract_log_lik(m.NATL_spatial.0$model) #null
ll.NATL.svl <- extract_log_lik(m.NATL_spatial.svl$model) #svl
ll.NATL.sex <- extract_log_lik(m.NATL_spatial.sex$model) #svl+sex
ll.NATL.svl.sex <- extract_log_lik(m.NATL_spatial.svl.sex$model) #svl+sex

waic(ll.NATL.0)
waic(ll.NATL.svl)
waic(ll.NATL.sex)
waic(ll.NATL.svl.sex)

loo(m.NATL_spatial.0) #null
loo(m.NATL_spatial.svl) #svl
loo(m.NATL_spatial.sex) #svl
loo(m.NATL_spatial.svl.sex) #svl+sex




ll.EnvU.0 <- extract_log_lik(m.EnvU_spatial.0$model) #null
ll.EnvU.svl <- extract_log_lik(m.EnvU_spatial.svl$model) #svl
ll.EnvU.sex <- extract_log_lik(m.EnvU_spatial.sex$model) #svl+sex
ll.EnvU.svl.sex <- extract_log_lik(m.EnvU_spatial.svl.sex$model) #svl+sex

waic(ll.EnvU.0)
waic(ll.EnvU.svl)
waic(ll.EnvU.sex)
waic(ll.EnvU.svl.sex)

loo(m.EnvU_spatial.0) #null
loo(m.EnvU_spatial.svl) #svl
loo(m.EnvU_spatial.sex) #svl
loo(m.EnvU_spatial.svl.sex) #svl+sex




ll.EnvF.0 <- extract_log_lik(m.EnvF_spatial.0$model) #null
ll.EnvF.svl <- extract_log_lik(m.EnvF_spatial.svl$model) #svl
ll.EnvF.sex <- extract_log_lik(m.EnvF_spatial.sex$model) #svl+sex
ll.EnvF.svl.sex <- extract_log_lik(m.EnvF_spatial.svl.sex$model) #svl+sex

waic(ll.EnvF.0)
waic(ll.EnvF.svl)
waic(ll.EnvF.sex)
waic(ll.EnvF.svl.sex)

loo(m.EnvF_spatial.0) #null
loo(m.EnvF_spatial.svl) #svl
loo(m.EnvF_spatial.sex) #svl
loo(m.EnvF_spatial.svl.sex) #svl+sex




ll.BTPU.0 <- extract_log_lik(m.BTPU_spatial.0$model) #null
ll.BTPU.svl <- extract_log_lik(m.BTPU_spatial.svl$model) #svl
ll.BTPU.sex <- extract_log_lik(m.BTPU_spatial.sex$model) #svl+sex
ll.BTPU.svl.sex <- extract_log_lik(m.BTPU_spatial.svl.sex$model) #svl+sex

waic(ll.BTPU.0)
waic(ll.BTPU.svl)
waic(ll.BTPU.sex)
waic(ll.BTPU.svl.sex)

loo(m.BTPU_spatial.0) #null
loo(m.BTPU_spatial.svl) #svl
loo(m.BTPU_spatial.sex) #svl
loo(m.BTPU_spatial.svl.sex) #svl+sex




ll.BTPF.0 <- extract_log_lik(m.BTPF_spatial.0$model) #null
ll.BTPF.svl <- extract_log_lik(m.BTPF_spatial.svl$model) #svl
ll.BTPF.sex <- extract_log_lik(m.BTPF_spatial.sex$model) #svl+sex
ll.BTPF.svl.sex <- extract_log_lik(m.BTPF_spatial.svl.sex$model) #svl+sex

waic(ll.BTPF.0)
waic(ll.BTPF.svl)
waic(ll.BTPF.sex)
waic(ll.BTPF.svl.sex)

loo(m.BTPF_spatial.0) #null
loo(m.BTPF_spatial.svl) #svl
loo(m.BTPF_spatial.sex) #svl
loo(m.BTPF_spatial.svl.sex) #svl+sex




ll.KewU.0 <- extract_log_lik(m.KewU_spatial.0$model) #null
ll.KewU.svl <- extract_log_lik(m.KewU_spatial.svl$model) #svl
ll.KewU.sex <- extract_log_lik(m.KewU_spatial.sex$model) #svl+sex
ll.KewU.svl.sex <- extract_log_lik(m.KewU_spatial.svl.sex$model) #svl+sex

waic(ll.KewU.0)
waic(ll.KewU.svl)
waic(ll.KewU.sex)
waic(ll.KewU.svl.sex)

loo(m.KewU_spatial.0) #null
loo(m.KewU_spatial.svl) #svl
loo(m.KewU_spatial.sex) #svl
loo(m.KewU_spatial.svl.sex) #svl+sex




ll.KewF.0 <- extract_log_lik(m.KewF_spatial.0$model) #null
ll.KewF.svl <- extract_log_lik(m.KewF_spatial.svl$model) #svl
ll.KewF.sex <- extract_log_lik(m.KewF_spatial.sex$model) #svl+sex
ll.KewF.svl.sex <- extract_log_lik(m.KewF_spatial.svl.sex$model) #svl+sex

waic(ll.KewF.0)
waic(ll.KewF.svl)
waic(ll.KewF.sex)
waic(ll.KewF.svl.sex)

loo(m.KewF_spatial.0) #null
loo(m.KewF_spatial.svl) #svl
loo(m.KewF_spatial.sex) #svl
loo(m.KewF_spatial.svl.sex) #svl+sex






#MAKING PREDICTIVE PLOTS####
#We will use the function predict in the same way we use it to make predictions
#from non-spatial linear models or GLM. 
#This is going to give us a probability of infection for each sample with their credible intervals. 
#We can also calculate the highest posterior density for each parameter.
p <- predict(m.EnvF_spatial.svl, type = "response")
head(p)
head(tidy(m.EnvF_spatial.svl, conf.int = TRUE, conf.method = "HPDinterval"))

#We can also use the function expand.grid to make a surface that includes
#the range of x and y coordiantes to make spatial predictions. 
#To do this we can add the covariate svl and specific predictions for each spatial
#point in the grid. Note that we are doing it in this example for mean value of svl. 
#You can also try it for large and small lizards to see if the predictions change.
pred_grid <- expand.grid(lat = seq(min(d.EnvF$lat), max(d.EnvF$lat), length.out = 50),
                         lon = seq(min(d.EnvF$lon), max(d.EnvF$lon), length.out = 50))
pred_grid$svl <- mean(d.EnvF$svl)
#pred_grid$sex <- mean(as.numeric(d.LawF$sex))
#pred_grid$sex1 <- pred_grid$svl*pred_grid$sex
pred_grid$prediction <- predict(m.EnvF_spatial.svl, newdata = pred_grid,
                                type = "response")$estimate
head(pred_grid)

#Now we are ready to make a pretty plot of these spatial predictions.
ggplot(pred_grid, aes(lon, lat, fill = prediction)) +
  geom_raster() +
  viridis::scale_fill_viridis()


