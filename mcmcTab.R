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

## ------------------------------------------------------------------------
library("R2jags")
set.seed(123)
infection.svl.sex.habitat <- data.frame(y = diagnosed$infection, n = nrow(diagnosed), 
                                        x1 = diagnosed$SVL, x2 = diagnosed$SEX, x3 = diagnosed$habitat)
mod.infection.svl.sex.habitat=jags(model.file=model_binom_3,data=infection.svl.sex.habitat,
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                                   parameters.to.save=c("alpha","b1","b2","b3","b4","b5"))


## ------------------------------------------------------------------------
mod.infection.svl.sex.habitat.mcmc <- as.mcmc(mod.infection.svl.sex.habitat)


## ------------------------------------------------------------------------
mod.infection.svl.sex.habitat.mat <- as.matrix(mod.infection.svl.sex.habitat.mcmc)
mod.infection.svl.sex.habitat.dat <- as.data.frame(mod.infection.svl.sex.habitat.mat)


## ------------------------------------------------------------------------
mod.infection.svl.sex.habitat.betas <- mod.infection.svl.sex.habitat.dat[, grep(x = colnames(mod.infection.svl.sex.habitat.dat), 
                                          pattern = "b", 
                                          fixed = TRUE)]


## ------------------------------------------------------------------------
head(mod.infection.svl.sex.habitat.betas)


## ------------------------------------------------------------------------
devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")
mcmctab(mod.infection.svl.sex.habitat)
mod.infection.svl.sex.habitat.tab <- mcmctab(mod.infection.svl.sex.habitat)

## ---- fig.width = 5, fig.height = 2--------------------------------------
library("ggmcmc")
infection.svl.sex.habitat.gg <- ggs(as.mcmc(mod.infection.svl.sex.habitat))
ggs_caterpillar(D = infection.svl.sex.habitat.gg, family = "b")


## ---- fig.width = 5, fig.height = 2--------------------------------------
ggs_caterpillar(D = infection.svl.sex.habitat.gg, family = "b") + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = c("SVL", "Sex",
                              "Habitat", "SVL*Sex",
                              "SVL*Habitat")) #WRONG 

## ------------------------------------------------------------------------
infection.svl.sex.habitat.posterior.coefs <- select(mod.infection.svl.sex.habitat.dat, `b1`, `b2`, `b3`,`b4`,`b5`)
names(infection.svl.sex.habitat.posterior.coefs) <- c("SVL", "Sex",
                                   "Habitat", "SVL*Sex",
                                   "SVL*Habitat")

## ------------------------------------------------------------------------
infection.svl.sex.habitat.posterior.coefs.long <- gather(infection.svl.sex.habitat.posterior.coefs)
head(infection.svl.sex.habitat.posterior.coefs.long)


## ---- fig.width = 5, fig.height = 3--------------------------------------
library("ggridges")
ggplot(data = infection.svl.sex.habitat.posterior.coefs.long, aes(x = value, y = key)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), 
                      alpha = 0.7) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() + xlab("Posterior estimate") + ylab("") 


## ---- fig.width = 5, fig.height = 2--------------------------------------
infection.svl.sex.habitat.posterior.coefs.sum <- summarize(group_by(infection.svl.sex.habitat.posterior.coefs.long, key),
                                        median_coef = median(value),
                                        lower_coef = quantile(value, probs = c(0.025)),
                                        upper_coef = quantile(value, probs = c(0.975)))
ggplot(data = infection.svl.sex.habitat.posterior.coefs.sum, aes(x = median_coef, y = key)) + 
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) + 
  geom_point() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  xlab("Posterior estimate") + ylab("") + 
  theme_bw()


