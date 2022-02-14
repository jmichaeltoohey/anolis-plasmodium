data <- read.csv("datasheets/Spatial Heterogeneity.csv",header=TRUE)

library("R2jags")
library(runjags)
library(dplyr)

data <- filter(data, !is.na(Infeccion))
for(i in 1:nrow(data)){
  if(data$Sex[i] == "Female"){data$Sex[i] = 0}
  else if(data$Sex[i] == "Male"){data$Sex[i] = 1}}
data$Sex <- as.numeric(data$Sex)
####p(infection)~####
sum(data$Infeccion)/nrow(data)
data <- filter(data, !is.na(xcoord))


#p(infection)~1
PR.infection.null <- data.frame(y = data$Infeccion, n = nrow(data), x1 = 1)
mod.PR.infection.null=jags(model.file=model_binom_1,data=PR.infection.null,
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        parameters.to.save=c("alpha","LogLik"))
#mod.PR.infection.null

#p(infection)  ~ SVL
PR.infection.svl <- data.frame(y = data$Infeccion, n = nrow(data), x1 = data$SVL)
mod.PR.infection.svl=jags(model.file=model_binom_1,data=PR.infection.svl,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))
#mod.PR.infection.svl

#p(infection)  ~ Sex
PR.infection.sex <- data.frame(y = data$Infeccion, n = nrow(data), x1 = data$Sex)
mod.PR.infection.sex=jags(model.file=model_binom_1,data=PR.infection.sex,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parameters.to.save=c("alpha","b1","LogLik"))
#mod.PR.infection.sex


#p(infection)~(SVL*Sex)
PR.infection.svl.sex <- data.frame(y = data$Infeccion, n = nrow(data), 
                                x1 = data$SVL, x2 = data$Sex)
mod.PR.infection.svl.sex=jags(model.file=model_binom_2mult,data=PR.infection.svl.sex,
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                           parameters.to.save=c("alpha","b1","b2","b3","LogLik"))
#mod.PR.infection.svl.sex


#Model Selection####
library(AICcmodavg)
library(loo)

loglik.PR.infection.null <- mod.PR.infection.null$BUGSoutput$sims.list$LogLik
loglik.PR.infection.svl <- mod.PR.infection.svl$BUGSoutput$sims.list$LogLik
loglik.PR.infection.sex <- mod.PR.infection.sex$BUGSoutput$sims.list$LogLik
loglik.PR.infection.svl.sex <- mod.PR.infection.svl.sex$BUGSoutput$sims.list$LogLik

waic.PR.infection.null <- waic(loglik.PR.infection.null)
waic.PR.infection.svl <- waic(loglik.PR.infection.svl)
waic.PR.infection.sex <- waic(loglik.PR.infection.sex)
waic.PR.infection.svl.sex <- waic(loglik.PR.infection.svl.sex)

loo.PR.infection.null <- loo(loglik.PR.infection.null)
loo.PR.infection.svl <- loo(loglik.PR.infection.svl)
loo.PR.infection.sex <- loo(loglik.PR.infection.sex)
loo.PR.infection.svl.sex <- loo(loglik.PR.infection.svl.sex)

PR.infection.list <- list(waic.PR.infection.null, waic.PR.infection.svl, waic.PR.infection.sex,
                       waic.PR.infection.svl.sex)

#PR.infection.loo.list <- list(loo.PR.infection.null, loo.PR.infection.svl, loo.PR.infection.sex,
#                           loo.PR.infection.svl.sex)


#waic.infection.null$estimates["waic",]
#waic.infection.svl$estimates["waic",]

PR.infection.compare = loo_compare(PR.infection.list)
#PR.infection.loo.compare = loo_compare(PR.infection.loo.list)
#infection_names <- c("p(infection) ~ (SVL*Habitat)","p(infection) ~ Sex", "p(infection) ~ SVL",
#                     "p(infection) ~ (SVL*Sex)", "p(infection) ~ 1", "p(infection) ~ (SVL*Sex)+(SVL*Habitat)",
#                     "p(infection) ~ Habitat", "p(infection) ~ (Sex)+(Habitat)")
PR.infection_names <- c("p(infection) ~ SVL",
                     "p(infection) ~ Sex","p(infection) ~ (SVL*Sex)","p(infection) ~ 1")

rownames(PR.infection.compare) <- PR.infection_names
print(PR.infection.compare, simplify=FALSE)

jags(model.file=model_binom_1,data=PR.infection.svl,
     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
     parameters.to.save=c("b0","b1"))

t.sex <- t.test(data$Infeccion~data$Sex)
biserial.cor(x=data$SVL,y=data$Sex)

#print(PR.infection.loo.compare, simplify=FALSE)

#weights.infection <- data.frame(length=8)
#sums=0
#for(i in 1:length(infection.list)){
#  sums <- sums+ exp(0.5*(infection.compare[1,7]-infection.compare[i,7]))
#}
#for(i in 1:length(infection.list)){
#  weights.infection[i] <- exp(0.5*(infection.compare[1,7]-infection.compare[i,7]))/sums
#}

PR.weights.infection <- data.frame(length=2) #three strong models
sums=0
for(i in 1:2){
  sums <- sums+ exp(0.5*(PR.infection.compare[1,7]-PR.infection.compare[i,7]))
}
for(i in 1:2){
  PR.weights.infection[i] <- exp(0.5*(PR.infection.compare[1,7]-PR.infection.compare[i,7]))/sums
}

#body condition####
PR.data.males <- filter(data, Sex == 1)
PR.male.mass.svl <- lm(Weight~SVL, data=PR.data.males, na.action=na.exclude)
PR.data.males$bodycond <- resid(PR.male.mass.svl)

#body condition~1
PR.male.bodycond.null <- data.frame(y = PR.data.males$bodycond,   x1 = 0)
PR.mod.male.bodycond.null=jags(model.file=model_linear_1,data=PR.male.bodycond.null,
                            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                            parameters.to.save=c("alpha","LogLik"))

#body condition ~ infection
PR.male.bodycond.infection <- data.frame(y = PR.data.males$bodycond,  x1 = PR.data.males$Infeccion)
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
print(PR.male.bodycond.compare, simplify=FALSE)

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


#body condition~1
PR.bodycond.null = lm(PR.data.males$bodycond ~ 1)
summary(PR.bodycond.null)

BF.PR.bodycond4.1 <- BF(PR.bodycond.null)
summary(BF.PR.bodycond4.1)

#body condition ~ infection
PR.bodycond.infection = lm(PR.data.males$bodycond ~ PR.data.males$Infeccion)
summary(PR.bodycond.infection)

BF.PR.bodycond1.1 <- BF(PR.bodycond.infection,
                     hypothesis = "Infeccion=0")
summary(BF.PR.bodycond1.1)
