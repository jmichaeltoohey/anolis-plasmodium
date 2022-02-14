#source("bayesian_models.R")
library(AICcmodavg)
library(loo)
library(knitr)
library(BayesPostEst)
devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")


####DIC model selection table####
#p(infection)
modelsel_infection <- list(mod.infection.null, mod.infection.svl, mod.infection.sex,
                           mod.infection.habitat, mod.infection.svl.sex,
                           mod.infection.svl.habitat, mod.infection.sex.habitat,
                           mod.infection.svl.sex.habitat)
infection_names <- c("p(infection) ~ 1","p(infection) ~ SVL", "p(infection) ~ Sex",
                     "p(infection) ~ Habitat", "p(infection) ~ (SVL*Sex)", "p(infection) ~ (SVL*Habitat)",
                     "p(infection) ~ (Sex)+(Habitat)", "p(infection) ~ (SVL*Sex)+(SVL*Habitat)")
dictab(cand.set = modelsel_infection, modnames = infection_names)




####wAIC model selection table####
#p(infection)####
loglik.infection.null <- mod.infection.null$BUGSoutput$sims.list$LogLik
loglik.infection.svl <- mod.infection.svl$BUGSoutput$sims.list$LogLik
loglik.infection.sex <- mod.infection.sex$BUGSoutput$sims.list$LogLik
loglik.infection.habitat <- mod.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.infection.svl.sex <- mod.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.infection.svl.habitat <- mod.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.infection.sex.habitat <- mod.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.infection.svl.sex.habitat <- mod.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.infection.null <- waic(loglik.infection.null)
waic.infection.svl <- waic(loglik.infection.svl)
waic.infection.sex <- waic(loglik.infection.sex)
waic.infection.habitat <- waic(loglik.infection.habitat)
waic.infection.svl.sex <- waic(loglik.infection.svl.sex)
waic.infection.svl.habitat <- waic(loglik.infection.svl.habitat)
waic.infection.sex.habitat <- waic(loglik.infection.sex.habitat)
waic.infection.svl.sex.habitat <- waic(loglik.infection.svl.sex.habitat)

infection.list <- list(waic.infection.null, waic.infection.svl, waic.infection.sex,
                       waic.infection.habitat, waic.infection.svl.sex, waic.infection.svl.habitat,
                       waic.infection.sex.habitat, waic.infection.svl.sex.habitat)

infection.compare = loo_compare(infection.list)
infection.compare

infection_names <- c("p(infection) ~ (SVL*Habitat)", "p(infection) ~ (SVL*Sex)+(SVL*Habitat)",  
                     "p(infection) ~ (Sex)+(Habitat)","p(infection) ~ Habitat","p(infection) ~ SVL",  
                     "p(infection) ~ (SVL*Sex)", "p(infection) ~ Sex","p(infection) ~ 1")

rownames(infection.compare) <- infection_names
print(infection.compare, simplify=FALSE)

weights.infection <- data.frame(length=8)
sums=0
for(i in 1:length(infection.list)){
  sums <- sums+ exp(0.5*(infection.compare[1,7]-infection.compare[i,7]))
}
for(i in 1:length(infection.list)){
  weights.infection[i] <- exp(0.5*(infection.compare[1,7]-infection.compare[i,7]))/sums
}

infection_modsel_table <- as.data.frame(infection.compare)[7:8]
infection_modsel_table$weight <- t(weights.infection)
print(infection_modsel_table)
write.table(infection_modsel_table, file="infection_mod.sel.txt", sep=",", quote=FALSE,row.names=T)


weights.infection <- data.frame(length=2) #three strong models
sums=0
for(i in 1:2){
  sums <- sums+ exp(0.5*(infection.compare[1,7]-infection.compare[i,7]))
}
for(i in 1:2){
  weights.infection[i] <- exp(0.5*(infection.compare[1,7]-infection.compare[i,7]))/sums
}
weights.infection

#body condition####

loglik.M.bodycond.null <- mod.male.bodycond.null$BUGSoutput$sims.list$LogLik
loglik.M.bodycond.habitat <- mod.male.bodycond.habitat$BUGSoutput$sims.list$LogLik
loglik.M.bodycond.infection <- mod.male.bodycond.infection$BUGSoutput$sims.list$LogLik
loglik.M.bodycond.infection.habitat <- mod.male.bodycond.infection.habitat$BUGSoutput$sims.list$LogLik

waic.M.bodycond.null <- waic(loglik.M.bodycond.null)
waic.M.bodycond.habitat <- waic(loglik.M.bodycond.habitat)
waic.M.bodycond.infection <- waic(loglik.M.bodycond.infection)
waic.M.bodycond.infection.habitat <- waic(loglik.M.bodycond.infection.habitat)

M.bodycond.list <- list(waic.M.bodycond.null,
                       waic.M.bodycond.habitat,
                       waic.M.bodycond.infection,
                       waic.M.bodycond.infection.habitat)

M.bodycond.compare = loo_compare(M.bodycond.list)
M_bodycond_names <- c("M Body Condition ~ Infection","M Body Condition ~ Infection+(Habitat)",
                      "M Body Condition ~ 1","M Body Condition ~ Habitat")

rownames(M.bodycond.compare) <- M_bodycond_names
print(M.bodycond.compare, simplify=FALSE)

weights.M.bodycond <- data.frame(length=8)
sums=0
for(i in 1:length(M.bodycond.list)){
  sums <- sums+ exp(0.5*(M.bodycond.compare[1,7]-M.bodycond.compare[i,7]))
}
for(i in 1:length(M.bodycond.list)){
  weights.M.bodycond[i] <- exp(0.5*(M.bodycond.compare[1,7]-M.bodycond.compare[i,7]))/sums
}
weights.M.bodycond

M.bodycond_modsel_table <- as.data.frame(M.bodycond.compare)[7:8]
M.bodycond_modsel_table$weight <- t(weights.M.bodycond)
print(M.bodycond_modsel_table)
write.table(M.bodycond_modsel_table, file="M.bodycond_mod.sel.txt", sep=",", quote=FALSE,row.names=T)


weights.M.bodycond <- data.frame(length=4)
sums=0
for(i in 1:4){
  sums <- sums+ exp(0.5*(M.bodycond.compare[1,7]-M.bodycond.compare[i,7]))
}
for(i in 1:4){
  weights.M.bodycond[i] <- exp(0.5*(M.bodycond.compare[1,7]-M.bodycond.compare[i,7]))/sums
}
weights.M.bodycond


#Hct####
loglik.Hct.null <- mod.Hct.null$BUGSoutput$sims.list$LogLik
loglik.Hct.svl <- mod.Hct.svl$BUGSoutput$sims.list$LogLik
loglik.Hct.sex <- mod.Hct.sex$BUGSoutput$sims.list$LogLik
loglik.Hct.habitat <- mod.Hct.habitat$BUGSoutput$sims.list$LogLik
loglik.Hct.svl.sex <- mod.Hct.svl.sex$BUGSoutput$sims.list$LogLik
loglik.Hct.svl.habitat <- mod.Hct.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.Hct.sex.habitat <- mod.Hct.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Hct.svl.sex.habitat <- mod.Hct.svl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Hct.infection <- mod.Hct.infection$BUGSoutput$sims.list$LogLik
loglik.Hct.infection.svl <- mod.Hct.infection.svl$BUGSoutput$sims.list$LogLik
loglik.Hct.infection.sex <- mod.Hct.infection.sex$BUGSoutput$sims.list$LogLik
loglik.Hct.infection.habitat <- mod.Hct.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.Hct.infection.svl.sex <- mod.Hct.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.Hct.infection.svl.habitat <- mod.Hct.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.Hct.infection.sex.habitat <- mod.Hct.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Hct.infection.svl.sex.habitat <- mod.Hct.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.Hct.null <- waic(loglik.Hct.null)
waic.Hct.svl <- waic(loglik.Hct.svl)
waic.Hct.sex <- waic(loglik.Hct.sex)
waic.Hct.habitat <- waic(loglik.Hct.habitat)
waic.Hct.svl.sex <- waic(loglik.Hct.svl.sex)
waic.Hct.svl.habitat <- waic(loglik.Hct.svl.habitat)
waic.Hct.sex.habitat <- waic(loglik.Hct.sex.habitat)
waic.Hct.svl.sex.habitat <- waic(loglik.Hct.svl.sex.habitat)
waic.Hct.infection <- waic(loglik.Hct.infection)
waic.Hct.infection.svl <- waic(loglik.Hct.infection.svl)
waic.Hct.infection.sex <- waic(loglik.Hct.infection.sex)
waic.Hct.infection.habitat <- waic(loglik.Hct.infection.habitat)
waic.Hct.infection.svl.sex <- waic(loglik.Hct.infection.svl.sex)
waic.Hct.infection.svl.habitat <- waic(loglik.Hct.infection.svl.habitat)
waic.Hct.infection.sex.habitat <- waic(loglik.Hct.infection.sex.habitat)
waic.Hct.infection.svl.sex.habitat <- waic(loglik.Hct.infection.svl.sex.habitat)

Hct.list <- list(waic.Hct.null, waic.Hct.svl, waic.Hct.sex,
                      waic.Hct.habitat, waic.Hct.svl.sex, waic.Hct.svl.habitat,
                      waic.Hct.sex.habitat, waic.Hct.svl.sex.habitat,
                      waic.Hct.infection, waic.Hct.infection.svl, waic.Hct.infection.sex,
                      waic.Hct.infection.habitat, waic.Hct.infection.svl.sex, waic.Hct.infection.svl.habitat,
                      waic.Hct.infection.sex.habitat, waic.Hct.infection.svl.sex.habitat)

#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

Hct.compare <- loo_compare(Hct.list)
Hct_names <- c("Hct ~ Infection + (SVL*Sex)","Hct ~ 1","Hct ~ Infection",
               "Hct ~ (SVL*Sex)","Hct ~ Sex","Hct ~ Infection + Sex","Hct ~ SVL","Hct ~ Habitat",
               "Hct ~ Infection + Habitat","Hct ~ Infection + SVL",
               "Hct ~ Infection + (SVL*Sex) + (SVL*Habitat)",
               "Hct ~ (SVL*Sex) + (SVL*Habitat)",
               "Hct ~ Sex + Habitat","Hct ~ Infection + (Sex) + (Habitat)",
               "Hct ~ (SVL*Habitat)",
               "Hct ~ Infection + (SVL*Habitat)")

rownames(Hct.compare) <- Hct_names
print(Hct.compare, simplify=FALSE)

weights.Hct <- data.frame(length=8)
sums=0
for(i in 1:length(Hct.list)){
  sums <- sums+ exp(0.5*(Hct.compare[1,7]-Hct.compare[i,7]))
}
for(i in 1:length(Hct.list)){
  weights.Hct[i] <- exp(0.5*(Hct.compare[1,7]-Hct.compare[i,7]))/sums
}

Hct_modsel_table <- as.data.frame(Hct.compare)[7:8]
Hct_modsel_table$weight <- t(weights.Hct)
print(Hct_modsel_table)
write.table(Hct_modsel_table, file="Hct_mod.sel.txt", sep=",", quote=FALSE,row.names=T)


weights.Hct <- data.frame(length=5)
sums=0
for(i in 1:5){
  sums <- sums+ exp(0.5*(Hct.compare[1,7]-Hct.compare[i,7]))
}
for(i in 1:5){
  weights.Hct[i] <- exp(0.5*(Hct.compare[1,7]-Hct.compare[i,7]))/sums
}
weights.Hct

#Hct_output <- data.frame(intercept=double(),infection=double(),SVL=double(),sex=double(),deviance=double())
#as.data.frame(Hct_output)
#Hct_output$intercept <- weights.Hct[1]*mod.Hct.infection$BUGSoutput$sims.list$b0+
#  weights.Hct[2]*mod.Hct.infection.svl$BUGSoutput$sims.list$b0+
#  weights.Hct[3]*mod.Hct.infection.svl.sex$BUGSoutput$sims.list$b0
#Hct_pred
Hct.b0 <- weights.Hct[1]*mod.Hct.infection$BUGSoutput$sims.list$b0+
  weights.Hct[2]*mod.Hct.infection.svl$BUGSoutput$sims.list$b0+
  weights.Hct[3]*mod.Hct.infection.svl.sex$BUGSoutput$sims.list$b0
Hct.b1 <- weights.Hct[1]*mod.Hct.infection$BUGSoutput$sims.list$b1+
  weights.Hct[2]*mod.Hct.infection.svl$BUGSoutput$sims.list$b1+
  weights.Hct[3]*mod.Hct.infection.svl.sex$BUGSoutput$sims.list$b1
Hct.b0 <- rnorm(215, 2.80, 0.27)
Hct.b1 <- rnorm(215, -0.11, 0.04)
Hct.b2 <- rnorm(215, 0.14, 0.06)
Hct.b3 <- rnorm(215, 0.51, 0.24)
Hct.b4 <- rnorm(215, -0.11, 0.05)

mymodel.mcmc <- as.mcmc(mod.Hct.infection)


#Na####
loglik.Na.null <- mod.Na.null$BUGSoutput$sims.list$LogLik
loglik.Na.svl <- mod.Na.svl$BUGSoutput$sims.list$LogLik
loglik.Na.sex <- mod.Na.sex$BUGSoutput$sims.list$LogLik
loglik.Na.habitat <- mod.Na.habitat$BUGSoutput$sims.list$LogLik
loglik.Na.svl.sex <- mod.Na.svl.sex$BUGSoutput$sims.list$LogLik
loglik.Na.svl.habitat <- mod.Na.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.Na.sex.habitat <- mod.Na.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Na.svl.sex.habitat <- mod.Na.svl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Na.infection <- mod.Na.infection$BUGSoutput$sims.list$LogLik
loglik.Na.infection.svl <- mod.Na.infection.svl$BUGSoutput$sims.list$LogLik
loglik.Na.infection.sex <- mod.Na.infection.sex$BUGSoutput$sims.list$LogLik
loglik.Na.infection.habitat <- mod.Na.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.Na.infection.svl.sex <- mod.Na.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.Na.infection.svl.habitat <- mod.Na.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.Na.infection.sex.habitat <- mod.Na.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Na.infection.svl.sex.habitat <- mod.Na.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.Na.null <- waic(loglik.Na.null)
waic.Na.svl <- waic(loglik.Na.svl)
waic.Na.sex <- waic(loglik.Na.sex)
waic.Na.habitat <- waic(loglik.Na.habitat)
waic.Na.svl.sex <- waic(loglik.Na.svl.sex)
waic.Na.svl.habitat <- waic(loglik.Na.svl.habitat)
waic.Na.sex.habitat <- waic(loglik.Na.sex.habitat)
waic.Na.svl.sex.habitat <- waic(loglik.Na.svl.sex.habitat)
waic.Na.infection <- waic(loglik.Na.infection)
waic.Na.infection.svl <- waic(loglik.Na.infection.svl)
waic.Na.infection.sex <- waic(loglik.Na.infection.sex)
waic.Na.infection.habitat <- waic(loglik.Na.infection.habitat)
waic.Na.infection.svl.sex <- waic(loglik.Na.infection.svl.sex)
waic.Na.infection.svl.habitat <- waic(loglik.Na.infection.svl.habitat)
waic.Na.infection.sex.habitat <- waic(loglik.Na.infection.sex.habitat)
waic.Na.infection.svl.sex.habitat <- waic(loglik.Na.infection.svl.sex.habitat)

Na.list <- list(waic.Na.null, waic.Na.svl, waic.Na.sex,
                 waic.Na.habitat, waic.Na.svl.sex, waic.Na.svl.habitat,
                 waic.Na.sex.habitat, waic.Na.svl.sex.habitat,
                 waic.Na.infection, waic.Na.infection.svl, waic.Na.infection.sex,
                 waic.Na.infection.habitat, waic.Na.infection.svl.sex, waic.Na.infection.svl.habitat,
                 waic.Na.infection.sex.habitat, waic.Na.infection.svl.sex.habitat)

#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

Na.compare <- loo_compare(Na.list)
Na_names <- c("Na ~ 1", "Na ~ SVL", "Na ~ Sex", "Na ~ Habitat", "Na ~ Infection", 
              "Na ~ Infection + SVL","Na ~ Infection + Sex", "Na ~ (Sex) + (Habitat)",
              "Na ~ Infection + Habitat", "Na ~ (SVL*Habitat)", "Na ~ (SVL*Sex)",
              "Na ~ Infection + (Sex) + (Habitat)", "Na ~ Infection + (SVL*Habitat)",   
              "Na ~ Infection + (SVL*Sex)",
              "Na ~ (SVL*Sex) + (SVL*Habitat)",
              "Na ~ Infection + (SVL*Sex) + (SVL*Habitat)")
              
#library(textreg)
#mcmcReg(mod = list(mod.K.infection, mod.K.infection.sex.habitat),
#        pars = list(c('alpha', 'b'),c('alpha', 'b')),
#        pointest = "median",
#        coefnames = list(
#          c("Intercept", "Infection"),
#          c("Intercept", "Infection", "Sex", "Habitat")),
#        ci = 0.95,
##        caption = "K Models",
#        caption.above = TRU, doctype = F)
#mcmcReg(mod.K.infection,)

rownames(Na.compare) <- Na_names

weights.Na <- data.frame(length=8)
sums=0
for(i in 1:length(Na.list)){
  sums <- sums+ exp(0.5*(Na.compare[1,7]-Na.compare[i,7]))
}
for(i in 1:length(Na.list)){
  weights.Na[i] <- exp(0.5*(Na.compare[1,7]-Na.compare[i,7]))/sums
}
weights.Na

Na_modsel_table <- as.data.frame(Na.compare)[7:8]
Na_modsel_table$weight <- t(weights.Na)
print(Na_modsel_table)
write.table(Na_modsel_table, file="Na_mod.sel.txt", sep=",", quote=FALSE,row.names=T)

weights.Na <- data.frame(length=5)
sums=0
for(i in 1:5){
  sums <- sums+ exp(0.5*(Na.compare[1,7]-Na.compare[i,7]))
}
for(i in 1:5){
  weights.Na[i] <- exp(0.5*(Na.compare[1,7]-Na.compare[i,7]))/sums
}
weights.Na

#K####
loglik.K.null <- mod.K.null$BUGSoutput$sims.list$LogLik
loglik.K.svl <- mod.K.svl$BUGSoutput$sims.list$LogLik
loglik.K.sex <- mod.K.sex$BUGSoutput$sims.list$LogLik
loglik.K.habitat <- mod.K.habitat$BUGSoutput$sims.list$LogLik
loglik.K.svl.sex <- mod.K.svl.sex$BUGSoutput$sims.list$LogLik
loglik.K.svl.habitat <- mod.K.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.K.sex.habitat <- mod.K.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.K.svl.sex.habitat <- mod.K.svl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.K.infection <- mod.K.infection$BUGSoutput$sims.list$LogLik
loglik.K.infection.svl <- mod.K.infection.svl$BUGSoutput$sims.list$LogLik
loglik.K.infection.sex <- mod.K.infection.sex$BUGSoutput$sims.list$LogLik
loglik.K.infection.habitat <- mod.K.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.K.infection.svl.sex <- mod.K.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.K.infection.svl.habitat <- mod.K.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.K.infection.sex.habitat <- mod.K.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.K.infection.svl.sex.habitat <- mod.K.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.K.null <- waic(loglik.K.null)
waic.K.svl <- waic(loglik.K.svl)
waic.K.sex <- waic(loglik.K.sex)
waic.K.habitat <- waic(loglik.K.habitat)
waic.K.svl.sex <- waic(loglik.K.svl.sex)
waic.K.svl.habitat <- waic(loglik.K.svl.habitat)
waic.K.sex.habitat <- waic(loglik.K.sex.habitat)
waic.K.svl.sex.habitat <- waic(loglik.K.svl.sex.habitat)
waic.K.infection <- waic(loglik.K.infection)
waic.K.infection.svl <- waic(loglik.K.infection.svl)
waic.K.infection.sex <- waic(loglik.K.infection.sex)
waic.K.infection.habitat <- waic(loglik.K.infection.habitat)
waic.K.infection.svl.sex <- waic(loglik.K.infection.svl.sex)
waic.K.infection.svl.habitat <- waic(loglik.K.infection.svl.habitat)
waic.K.infection.sex.habitat <- waic(loglik.K.infection.sex.habitat)
waic.K.infection.svl.sex.habitat <- waic(loglik.K.infection.svl.sex.habitat)

K.list <- list(waic.K.null, waic.K.svl, waic.K.sex,
                 waic.K.habitat, waic.K.svl.sex, waic.K.svl.habitat,
                 waic.K.sex.habitat, waic.K.svl.sex.habitat,
                 waic.K.infection, waic.K.infection.svl, waic.K.infection.sex,
                 waic.K.infection.habitat, waic.K.infection.svl.sex, waic.K.infection.svl.habitat,
                 waic.K.infection.sex.habitat, waic.K.infection.svl.sex.habitat)

#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

K.compare <- loo_compare(K.list)
K_names <- c("K ~ SVL", "K ~ Infection + (SVL)", "K ~ Sex", "K ~ Infection + (Sex)",
             "K ~ (SVL*Sex)", "K ~ (Sex) + (Habitat)", "K ~ (SVL*Habitat)", 
             "K ~ Infection + (Sex) + (Habitat)", "K ~ Infection + (SVL*Sex)",
             "K ~ Infection + (SVL*Habitat)", "K ~ 1", "K ~ Infection",
             "K ~ (SVL*Sex) + (SVL*Habitat)", "K ~ Habitat",
             "K ~ Infection + (SVL*Sex) + (SVL*Habitat)","K ~ Infection + (Habitat)")

rownames(K.compare) <- K_names



weights.K <- data.frame(length=16)
sums=0
for(i in 1:length(K.list)){
  sums <- sums+ exp(0.5*(K.compare[1,7]-K.compare[i,7]))
}
for(i in 1:length(K.list)){
  weights.K[i] <- exp(0.5*(K.compare[1,7]-K.compare[i,7]))/sums
}
weights.K

K_modsel_table <- as.data.frame(K.compare)[7:8]
K_modsel_table$weight <- t(weights.K)
print(K_modsel_table)
write.table(K_modsel_table, file="K_mod.sel.txt", sep=",", quote=FALSE,row.names=T)

weights.K <- data.frame(length=4)
sums=0
for(i in 1:4){
  sums <- sums+ exp(0.5*(K.compare[1,7]-K.compare[i,7]))
}
for(i in 1:4){
  weights.K[i] <- exp(0.5*(K.compare[1,7]-K.compare[i,7]))/sums
}
weights.K

#iCa####
loglik.iCa.null <- mod.iCa.null$BUGSoutput$sims.list$LogLik
loglik.iCa.svl <- mod.iCa.svl$BUGSoutput$sims.list$LogLik
loglik.iCa.sex <- mod.iCa.sex$BUGSoutput$sims.list$LogLik
loglik.iCa.habitat <- mod.iCa.habitat$BUGSoutput$sims.list$LogLik
loglik.iCa.svl.sex <- mod.iCa.svl.sex$BUGSoutput$sims.list$LogLik
loglik.iCa.svl.habitat <- mod.iCa.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.iCa.sex.habitat <- mod.iCa.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.iCa.svl.sex.habitat <- mod.iCa.svl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.iCa.infection <- mod.iCa.infection$BUGSoutput$sims.list$LogLik
loglik.iCa.infection.svl <- mod.iCa.infection.svl$BUGSoutput$sims.list$LogLik
loglik.iCa.infection.sex <- mod.iCa.infection.sex$BUGSoutput$sims.list$LogLik
loglik.iCa.infection.habitat <- mod.iCa.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.iCa.infection.svl.sex <- mod.iCa.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.iCa.infection.svl.habitat <- mod.iCa.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.iCa.infection.sex.habitat <- mod.iCa.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.iCa.infection.svl.sex.habitat <- mod.iCa.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.iCa.null <- waic(loglik.iCa.null)
waic.iCa.svl <- waic(loglik.iCa.svl)
waic.iCa.sex <- waic(loglik.iCa.sex)
waic.iCa.habitat <- waic(loglik.iCa.habitat)
waic.iCa.svl.sex <- waic(loglik.iCa.svl.sex)
waic.iCa.svl.habitat <- waic(loglik.iCa.svl.habitat)
waic.iCa.sex.habitat <- waic(loglik.iCa.sex.habitat)
waic.iCa.svl.sex.habitat <- waic(loglik.iCa.svl.sex.habitat)
waic.iCa.infection <- waic(loglik.iCa.infection)
waic.iCa.infection.svl <- waic(loglik.iCa.infection.svl)
waic.iCa.infection.sex <- waic(loglik.iCa.infection.sex)
waic.iCa.infection.habitat <- waic(loglik.iCa.infection.habitat)
waic.iCa.infection.svl.sex <- waic(loglik.iCa.infection.svl.sex)
waic.iCa.infection.svl.habitat <- waic(loglik.iCa.infection.svl.habitat)
waic.iCa.infection.sex.habitat <- waic(loglik.iCa.infection.sex.habitat)
waic.iCa.infection.svl.sex.habitat <- waic(loglik.iCa.infection.svl.sex.habitat)

iCa.list <- list(waic.iCa.null, waic.iCa.svl, waic.iCa.sex,
                 waic.iCa.habitat, waic.iCa.svl.sex, waic.iCa.svl.habitat,
                 waic.iCa.sex.habitat, waic.iCa.svl.sex.habitat,
                 waic.iCa.infection, waic.iCa.infection.svl, waic.iCa.infection.sex,
                 waic.iCa.infection.habitat, waic.iCa.infection.svl.sex, waic.iCa.infection.svl.habitat,
                 waic.iCa.infection.sex.habitat, waic.iCa.infection.svl.sex.habitat)

#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

iCa.compare <- loo_compare(iCa.list)
iCa_names <- c("iCa ~ Infection + (SVL)","iCa ~ SVL","iCa ~ (SVL*Habitat)",
               "iCa ~ (SVL*Sex)","iCa ~ Infection + (SVL*Sex)",
               "iCa ~ (SVL*Sex) + (SVL*Habitat)","iCa ~ Infection + (SVL*Habitat)",
               "iCa ~ Infection + (SVL*Sex) + (SVL*Habitat)","iCa ~ (Sex) + (Habitat)",
               "iCa ~ Sex","iCa ~ Infection + (Sex)","iCa ~ Infection + (Sex) + (Habitat)",
               "iCa ~ 1","iCa ~ Habitat","iCa ~ Infection",
               "iCa ~ Infection + (Habitat)")

rownames(iCa.compare) <- iCa_names
print(iCa.compare, simplify=FALSE)

weights.iCa <- data.frame(length=8)
sums=0
for(i in 1:length(iCa.list)){
  sums <- sums+ exp(0.5*(iCa.compare[1,7]-iCa.compare[i,7]))
}
for(i in 1:length(iCa.list)){
  weights.iCa[i] <- exp(0.5*(iCa.compare[1,7]-iCa.compare[i,7]))/sums
}
weights.iCa

iCa_modsel_table <- as.data.frame(iCa.compare)[7:8]
iCa_modsel_table$weight <- t(weights.iCa)
print(iCa_modsel_table)
write.table(iCa_modsel_table, file="iCa_mod.sel.txt", sep=",", quote=FALSE,row.names=T)


weights.iCa <- data.frame(length=7)
sums=0
for(i in 1:7){
  sums <- sums+ exp(0.5*(iCa.compare[1,7]-iCa.compare[i,7]))
}
for(i in 1:7){
  weights.iCa[i] <- exp(0.5*(iCa.compare[1,7]-iCa.compare[i,7]))/sums
}
weights.iCa

#Cl####
loglik.Cl.null <- mod.Cl.null$BUGSoutput$sims.list$LogLik
loglik.Cl.svl <- mod.Cl.svl$BUGSoutput$sims.list$LogLik
loglik.Cl.sex <- mod.Cl.sex$BUGSoutput$sims.list$LogLik
loglik.Cl.habitat <- mod.Cl.habitat$BUGSoutput$sims.list$LogLik
loglik.Cl.svl.sex <- mod.Cl.svl.sex$BUGSoutput$sims.list$LogLik
loglik.Cl.svl.habitat <- mod.Cl.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.Cl.sex.habitat <- mod.Cl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Cl.svl.sex.habitat <- mod.Cl.svl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Cl.infection <- mod.Cl.infection$BUGSoutput$sims.list$LogLik
loglik.Cl.infection.svl <- mod.Cl.infection.svl$BUGSoutput$sims.list$LogLik
loglik.Cl.infection.sex <- mod.Cl.infection.sex$BUGSoutput$sims.list$LogLik
loglik.Cl.infection.habitat <- mod.Cl.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.Cl.infection.svl.sex <- mod.Cl.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.Cl.infection.svl.habitat <- mod.Cl.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.Cl.infection.sex.habitat <- mod.Cl.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Cl.infection.svl.sex.habitat <- mod.Cl.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.Cl.null <- waic(loglik.Cl.null)
waic.Cl.svl <- waic(loglik.Cl.svl)
waic.Cl.sex <- waic(loglik.Cl.sex)
waic.Cl.habitat <- waic(loglik.Cl.habitat)
waic.Cl.svl.sex <- waic(loglik.Cl.svl.sex)
waic.Cl.svl.habitat <- waic(loglik.Cl.svl.habitat)
waic.Cl.sex.habitat <- waic(loglik.Cl.sex.habitat)
waic.Cl.svl.sex.habitat <- waic(loglik.Cl.svl.sex.habitat)
waic.Cl.infection <- waic(loglik.Cl.infection)
waic.Cl.infection.svl <- waic(loglik.Cl.infection.svl)
waic.Cl.infection.sex <- waic(loglik.Cl.infection.sex)
waic.Cl.infection.habitat <- waic(loglik.Cl.infection.habitat)
waic.Cl.infection.svl.sex <- waic(loglik.Cl.infection.svl.sex)
waic.Cl.infection.svl.habitat <- waic(loglik.Cl.infection.svl.habitat)
waic.Cl.infection.sex.habitat <- waic(loglik.Cl.infection.sex.habitat)
waic.Cl.infection.svl.sex.habitat <- waic(loglik.Cl.infection.svl.sex.habitat)

Cl.list <- list(waic.Cl.null, waic.Cl.svl, waic.Cl.sex,
                 waic.Cl.habitat, waic.Cl.svl.sex, waic.Cl.svl.habitat,
                 waic.Cl.sex.habitat, waic.Cl.svl.sex.habitat,
                 waic.Cl.infection, waic.Cl.infection.svl, waic.Cl.infection.sex,
                 waic.Cl.infection.habitat, waic.Cl.infection.svl.sex, waic.Cl.infection.svl.habitat,
                 waic.Cl.infection.sex.habitat, waic.Cl.infection.svl.sex.habitat)

#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

Cl.compare <- loo_compare(Cl.list)
Cl_names <- c("Cl ~ 1","Cl ~ SVL","Cl ~ Habitat","Cl ~ Infection","Cl ~ Sex",
              "Cl ~ Infection + (SVL)","Cl ~ (Sex) + (Habitat)","Cl ~ Infection + (Habitat)",
              "Cl ~ Infection + (Sex)","Cl ~ (SVL*Sex)","Cl ~ (SVL*Habitat)",
              "Cl ~ Infection + (Sex) + (Habitat)","Cl ~ Infection + (SVL*Sex)",
              "Cl ~ Infection + (SVL*Habitat)","Cl ~ (SVL*Sex) + (SVL*Habitat)",
              "Cl ~ Infection + (SVL*Sex) + (SVL*Habitat)")

rownames(Cl.compare) <- Cl_names
print(Cl.compare, simplify=FALSE)

weights.Cl <- data.frame(length=8)
sums=0
for(i in 1:length(Cl.list)){
  sums <- sums+ exp(0.5*(Cl.compare[1,7]-Cl.compare[i,7]))
}
for(i in 1:length(Cl.list)){
  weights.Cl[i] <- exp(0.5*(Cl.compare[1,7]-Cl.compare[i,7]))/sums
}
weights.Cl

Cl_modsel_table <- as.data.frame(Cl.compare)[7:8]
Cl_modsel_table$weight <- t(weights.Cl)
print(Cl_modsel_table)
write.table(Cl_modsel_table, file="Cl_mod.sel.txt", sep=",", quote=FALSE,row.names=T)

weights.Cl <- data.frame(length=16)
sums=0
for(i in 1:8){
  sums <- sums+ exp(0.5*(Cl.compare[1,7]-Cl.compare[i,7]))
}
for(i in 1:8){
  weights.Cl[i] <- exp(0.5*(Cl.compare[1,7]-Cl.compare[i,7]))/sums
}
weights.Cl


#Glu####
loglik.Glu.null <- mod.Glu.null$BUGSoutput$sims.list$LogLik
loglik.Glu.svl <- mod.Glu.svl$BUGSoutput$sims.list$LogLik
loglik.Glu.sex <- mod.Glu.sex$BUGSoutput$sims.list$LogLik
loglik.Glu.habitat <- mod.Glu.habitat$BUGSoutput$sims.list$LogLik
loglik.Glu.svl.sex <- mod.Glu.svl.sex$BUGSoutput$sims.list$LogLik
loglik.Glu.svl.habitat <- mod.Glu.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.Glu.sex.habitat <- mod.Glu.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Glu.svl.sex.habitat <- mod.Glu.svl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Glu.infection <- mod.Glu.infection$BUGSoutput$sims.list$LogLik
loglik.Glu.infection.svl <- mod.Glu.infection.svl$BUGSoutput$sims.list$LogLik
loglik.Glu.infection.sex <- mod.Glu.infection.sex$BUGSoutput$sims.list$LogLik
loglik.Glu.infection.habitat <- mod.Glu.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.Glu.infection.svl.sex <- mod.Glu.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.Glu.infection.svl.habitat <- mod.Glu.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.Glu.infection.sex.habitat <- mod.Glu.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.Glu.infection.svl.sex.habitat <- mod.Glu.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.Glu.null <- waic(loglik.Glu.null)
waic.Glu.svl <- waic(loglik.Glu.svl)
waic.Glu.sex <- waic(loglik.Glu.sex)
waic.Glu.habitat <- waic(loglik.Glu.habitat)
waic.Glu.svl.sex <- waic(loglik.Glu.svl.sex)
waic.Glu.svl.habitat <- waic(loglik.Glu.svl.habitat)
waic.Glu.sex.habitat <- waic(loglik.Glu.sex.habitat)
waic.Glu.svl.sex.habitat <- waic(loglik.Glu.svl.sex.habitat)
waic.Glu.infection <- waic(loglik.Glu.infection)
waic.Glu.infection.svl <- waic(loglik.Glu.infection.svl)
waic.Glu.infection.sex <- waic(loglik.Glu.infection.sex)
waic.Glu.infection.habitat <- waic(loglik.Glu.infection.habitat)
waic.Glu.infection.svl.sex <- waic(loglik.Glu.infection.svl.sex)
waic.Glu.infection.svl.habitat <- waic(loglik.Glu.infection.svl.habitat)
waic.Glu.infection.sex.habitat <- waic(loglik.Glu.infection.sex.habitat)
waic.Glu.infection.svl.sex.habitat <- waic(loglik.Glu.infection.svl.sex.habitat)

Glu.list <- list(waic.Glu.null, waic.Glu.svl, waic.Glu.sex,
                 waic.Glu.habitat, waic.Glu.svl.sex, waic.Glu.svl.habitat,
                 waic.Glu.sex.habitat, waic.Glu.svl.sex.habitat,
                 waic.Glu.infection, waic.Glu.infection.svl, waic.Glu.infection.sex,
                 waic.Glu.infection.habitat, waic.Glu.infection.svl.sex, waic.Glu.infection.svl.habitat,
                 waic.Glu.infection.sex.habitat, waic.Glu.infection.svl.sex.habitat)

#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

Glu.compare <- loo_compare(Glu.list)
Glu_names <- c("Glu ~ Sex","Glu ~ SVL","Glu ~ (SVL*Sex)","Glu ~ (Sex) + (Habitat)",
               "Glu ~ Infection + (Sex)",
               "Glu ~ Infection + (Sex) + (Habitat)","Glu ~ 1",
               "Glu ~ Infection + (SVL*Sex)","Glu ~ Infection + (SVL)",
               "Glu ~ (SVL*Habitat)","Glu ~ Habitat",
               "Glu ~ (SVL*Sex) + (SVL*Habitat)","Glu ~ Infection",
               "Glu ~ Infection + (SVL*Habitat)",
               "Glu ~ Infection + (Habitat)",
               "Glu ~ Infection + (SVL*Sex) + (SVL*Habitat)")

rownames(Glu.compare) <- Glu_names
print(Glu.compare, simplify=FALSE)

weights.Glu <- data.frame(length=16)
sums=0
for(i in 1:length(Glu.list)){
  sums <- sums+ exp(0.5*(Glu.compare[1,7]-Glu.compare[i,7]))
}
for(i in 1:length(Glu.list)){
  weights.Glu[i] <- exp(0.5*(Glu.compare[1,7]-Glu.compare[i,7]))/sums
}
weights.Glu

Glu_modsel_table <- as.data.frame(Glu.compare)[7:8]
Glu_modsel_table$weight <- t(weights.Glu)
print(Glu_modsel_table)
write.table(Glu_modsel_table, file="Glu_mod.sel.txt", sep=",", quote=FALSE,row.names=T)


weights.Glu <- data.frame(length=3)
sums=0
for(i in 1:3){
  sums <- sums+ exp(0.5*(Glu.compare[1,7]-Glu.compare[i,7]))
}
for(i in 1:3){
  weights.Glu[i] <- exp(0.5*(Glu.compare[1,7]-Glu.compare[i,7]))/sums
}
weights.Glu

#PC1####
loglik.PC1.null <- mod.PC1.null$BUGSoutput$sims.list$LogLik
loglik.PC1.svl <- mod.PC1.svl$BUGSoutput$sims.list$LogLik
loglik.PC1.sex <- mod.PC1.sex$BUGSoutput$sims.list$LogLik
loglik.PC1.habitat <- mod.PC1.habitat$BUGSoutput$sims.list$LogLik
loglik.PC1.svl.sex <- mod.PC1.svl.sex$BUGSoutput$sims.list$LogLik
loglik.PC1.svl.habitat <- mod.PC1.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.PC1.sex.habitat <- mod.PC1.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.PC1.svl.sex.habitat <- mod.PC1.svl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.PC1.infection <- mod.PC1.infection$BUGSoutput$sims.list$LogLik
loglik.PC1.infection.svl <- mod.PC1.infection.svl$BUGSoutput$sims.list$LogLik
loglik.PC1.infection.sex <- mod.PC1.infection.sex$BUGSoutput$sims.list$LogLik
loglik.PC1.infection.habitat <- mod.PC1.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.PC1.infection.svl.sex <- mod.PC1.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.PC1.infection.svl.habitat <- mod.PC1.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.PC1.infection.sex.habitat <- mod.PC1.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.PC1.infection.svl.sex.habitat <- mod.PC1.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.PC1.null <- waic(loglik.PC1.null)
waic.PC1.svl <- waic(loglik.PC1.svl)
waic.PC1.sex <- waic(loglik.PC1.sex)
waic.PC1.habitat <- waic(loglik.PC1.habitat)
waic.PC1.svl.sex <- waic(loglik.PC1.svl.sex)
waic.PC1.svl.habitat <- waic(loglik.PC1.svl.habitat)
waic.PC1.sex.habitat <- waic(loglik.PC1.sex.habitat)
waic.PC1.svl.sex.habitat <- waic(loglik.PC1.svl.sex.habitat)
waic.PC1.infection <- waic(loglik.PC1.infection)
waic.PC1.infection.svl <- waic(loglik.PC1.infection.svl)
waic.PC1.infection.sex <- waic(loglik.PC1.infection.sex)
waic.PC1.infection.habitat <- waic(loglik.PC1.infection.habitat)
waic.PC1.infection.svl.sex <- waic(loglik.PC1.infection.svl.sex)
waic.PC1.infection.svl.habitat <- waic(loglik.PC1.infection.svl.habitat)
waic.PC1.infection.sex.habitat <- waic(loglik.PC1.infection.sex.habitat)
waic.PC1.infection.svl.sex.habitat <- waic(loglik.PC1.infection.svl.sex.habitat)

PC1.list <- list(waic.PC1.null, waic.PC1.svl, waic.PC1.sex,
               waic.PC1.habitat, waic.PC1.svl.sex, waic.PC1.svl.habitat,
               waic.PC1.sex.habitat, waic.PC1.svl.sex.habitat,
               waic.PC1.infection, waic.PC1.infection.svl, waic.PC1.infection.sex,
               waic.PC1.infection.habitat, waic.PC1.infection.svl.sex, waic.PC1.infection.svl.habitat,
               waic.PC1.infection.sex.habitat, waic.PC1.infection.svl.sex.habitat)

#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

PC1.compare <- loo_compare(PC1.list)
PC1_names <- c("PC1 ~ 1", "PC1 ~ Infection", "PC1 ~ SVL", "PC1 ~ Sex", "PC1 ~ Habitat",
               "PC1 ~ Infection + SVL", "PC1 ~ Infection + Habitat", "PC1 ~ Infection + Sex",
               "PC1 ~ (SVL*Sex)", "PC1 ~ (Sex) + (Habitat)", "PC1 ~ (SVL*Habitat)",
               "PC1 ~ Infection + (Sex) + (Habitat)", "PC1 ~ Infection + (SVL*Habitat)",
               "PC1 ~ Infection + (SVL*Sex)", "PC1 ~ (SVL*Sex) + (SVL*Habitat)",
               "PC1 ~ Infection + (SVL*Sex) + (SVL*Habitat)")

rownames(PC1.compare) <- PC1_names
print(PC1.compare, simplify=FALSE)

weights.PC1 <- data.frame(length=16)
sums=0
for(i in 1:length(PC1.list)){
  sums <- sums + exp(0.5*(PC1.compare[1,7]-PC1.compare[i,7]))
}
for(i in 1:length(PC1.list)){
  weights.PC1[i] <- exp(0.5*(PC1.compare[1,7]-PC1.compare[i,7]))/sums
}

PC1_modsel_table <- as.data.frame(PC1.compare)[7:8]
PC1_modsel_table$weight <- t(weights.PC1)
print(PC1_modsel_table)
write.table(PC1_modsel_table, file="PC1_mod.sel.txt", sep=",", quote=FALSE,row.names=T)


#PC2####
loglik.PC2.null <- mod.PC2.null$BUGSoutput$sims.list$LogLik
loglik.PC2.svl <- mod.PC2.svl$BUGSoutput$sims.list$LogLik
loglik.PC2.sex <- mod.PC2.sex$BUGSoutput$sims.list$LogLik
loglik.PC2.habitat <- mod.PC2.habitat$BUGSoutput$sims.list$LogLik
loglik.PC2.svl.sex <- mod.PC2.svl.sex$BUGSoutput$sims.list$LogLik
loglik.PC2.svl.habitat <- mod.PC2.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.PC2.sex.habitat <- mod.PC2.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.PC2.svl.sex.habitat <- mod.PC2.svl.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.PC2.infection <- mod.PC2.infection$BUGSoutput$sims.list$LogLik
loglik.PC2.infection.svl <- mod.PC2.infection.svl$BUGSoutput$sims.list$LogLik
loglik.PC2.infection.sex <- mod.PC2.infection.sex$BUGSoutput$sims.list$LogLik
loglik.PC2.infection.habitat <- mod.PC2.infection.habitat$BUGSoutput$sims.list$LogLik
loglik.PC2.infection.svl.sex <- mod.PC2.infection.svl.sex$BUGSoutput$sims.list$LogLik
loglik.PC2.infection.svl.habitat <- mod.PC2.infection.svl.habitat$BUGSoutput$sims.list$LogLik
loglik.PC2.infection.sex.habitat <- mod.PC2.infection.sex.habitat$BUGSoutput$sims.list$LogLik
loglik.PC2.infection.svl.sex.habitat <- mod.PC2.infection.svl.sex.habitat$BUGSoutput$sims.list$LogLik

waic.PC2.null <- waic(loglik.PC2.null)
waic.PC2.svl <- waic(loglik.PC2.svl)
waic.PC2.sex <- waic(loglik.PC2.sex)
waic.PC2.habitat <- waic(loglik.PC2.habitat)
waic.PC2.svl.sex <- waic(loglik.PC2.svl.sex)
waic.PC2.svl.habitat <- waic(loglik.PC2.svl.habitat)
waic.PC2.sex.habitat <- waic(loglik.PC2.sex.habitat)
waic.PC2.svl.sex.habitat <- waic(loglik.PC2.svl.sex.habitat)
waic.PC2.infection <- waic(loglik.PC2.infection)
waic.PC2.infection.svl <- waic(loglik.PC2.infection.svl)
waic.PC2.infection.sex <- waic(loglik.PC2.infection.sex)
waic.PC2.infection.habitat <- waic(loglik.PC2.infection.habitat)
waic.PC2.infection.svl.sex <- waic(loglik.PC2.infection.svl.sex)
waic.PC2.infection.svl.habitat <- waic(loglik.PC2.infection.svl.habitat)
waic.PC2.infection.sex.habitat <- waic(loglik.PC2.infection.sex.habitat)
waic.PC2.infection.svl.sex.habitat <- waic(loglik.PC2.infection.svl.sex.habitat)

PC2.list <- list(waic.PC2.null, waic.PC2.svl, waic.PC2.sex,
               waic.PC2.habitat, waic.PC2.svl.sex, waic.PC2.svl.habitat,
               waic.PC2.sex.habitat, waic.PC2.svl.sex.habitat,
               waic.PC2.infection, waic.PC2.infection.svl, waic.PC2.infection.sex,
               waic.PC2.infection.habitat, waic.PC2.infection.svl.sex, waic.PC2.infection.svl.habitat,
               waic.PC2.infection.sex.habitat, waic.PC2.infection.svl.sex.habitat)

#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

PC2.compare <- loo_compare(PC2.list)
PC2_names <- c("PC2 ~ (SVL*Sex) + (SVL*Habitat)", "PC2 ~ (SVL*Sex)", "PC2 ~ Infection + (SVL*Sex)", 
             "PC2 ~ Infection + (SVL*Sex) + (SVL*Habitat)", "PC2 ~ Infection", 
             "PC2 ~ Habitat", "PC2 ~ 1", "PC2 ~ Infection + Habitat",
             "PC2 ~ Infection + Sex", "PC2 ~ Sex + Habitat", "PC2 ~ Infection + SVL", 
             "PC2 ~ Infection + (Sex) + (Habitat)", "PC2 ~ Sex", "PC2 ~ SVL",
             "PC2 ~ (SVL*Habitat)", "PC2 ~ Infection + (SVL*Habitat)")

rownames(PC2.compare) <- PC2_names
print(PC2.compare, simplify=FALSE)

weights.PC2 <- data.frame(length=16)
sums=0
for(i in 1:length(PC2.list)){
  sums <- sums + exp(0.5*(PC2.compare[1,7]-PC2.compare[i,7]))
}
for(i in 1:length(PC1.list)){
  weights.PC2[i] <- exp(0.5*(PC2.compare[1,7]-PC2.compare[i,7]))/sums
}

PC2_modsel_table <- as.data.frame(PC2.compare)[7:8]
PC2_modsel_table$weight <- t(weights.PC2)
print(PC2_modsel_table)
write.table(PC2_modsel_table, file="PC2_mod.sel.txt", sep=",", quote=FALSE,row.names=T)


#PRINTING####
print(infection.compare, simplify=FALSE)
print(M.bodycond.compare, simplify=FALSE)
print(Hct.compare, simplify=FALSE)
print(Na.compare, simplify=FALSE)
print(K.compare, simplify=FALSE)
print(iCa.compare, simplify=FALSE)


library(knitr); knit('anole-thesis-main/Tables.RmdTables.Rmd')
library(BayesPostEst)
mcmcReg(mod.Hct.infection, format = 'html', doctype = F)
mod.Hct.infection




