p.infection.sim <- data.frame(svl = rep(seq(from=min(diagnosed$SVL),
                                                    to=max(diagnosed$SVL),
                                                    by=0.05),4))
p.infection.sim$habitat <- 0
p.infection.sim$habitat[1:(nrow(p.infection.sim)/2)] <- 1
p.infection.sim$sex <- 0
p.infection.sim$sex[c(1:(nrow(p.infection.sim)/4),(nrow(p.infection.sim)/2+1):(3*nrow(p.infection.sim)/4))] <- 1
p.infection.sim <- as.matrix(p.infection.sim)

mod.infection.svl.sex.habitat.mat
Xb <- t(p.infection.sim %*% t(mod.infection.svl.sex.habitat.mat))

devtools::source_url("https://github.com/jkarreth/JKmisc/raw/master/MCMC_observed_probs.R")
model.matrix(infection.svl_X_sex__t__svl_X_habitat)

library(broom)
tidy_coefs <- tidy(mod.infection.svl.sex.habitat)

library(arm)
bayes_logistic_model_function <- 
  function(SVL, SEX, habitat) {
    invlogit(mod.infection.svl.sex.habitat.tab$Median[1]
                  + mod.infection.svl.sex.habitat.tab$Median[2]*SVL
                  + mod.infection.svl.sex.habitat.tab$Median[3]*SEX
                  + mod.infection.svl.sex.habitat.tab$Median[4]*habitat
                  + mod.infection.svl.sex.habitat.tab$Median[5]*SVL*SEX
                  + mod.infection.svl.sex.habitat.tab$Median[6]*SVL*habitat)
  }

ggplot(diagnosed, aes(x = SVL, y = infection)) +
  geom_jitter(width = 0, height = .05, alpha = 0.5) +
  stat_function(fun = bayes_logistic_model_function(diagnosed$SVL,diagnosed$SEX,diagnosed$habitat), 
                color = "blue", size = 1, data = diagnosed)
