library(MBA)
library(coda)
library(spBayes)
library(dplyr)
library(tidyverse)

set.seed(1)
source("data_carpentry.R")

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p)))
}

####infection####

##Isolate coordinates by site
LawU <- filter(diagnosed, SITE == "Law-U") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.LawU <- nrow(LawU) ##n.LawUumber of location.LawUs

phi <- 3/.5
sigma.sq <- 2

R.LawU <- sigma.sq*exp(-phi*as.matrix(dist(LawU)))
w.LawU <- rmvn(1, rep(0,n.LawU), R.LawU)

x.LawU <- as.matrix(rep(1,n.LawU))
beta <- 0.1
p.LawU <- 1/(1+exp(-(x.LawU%*%beta+w.LawU)))

weights.LawU <- rep(1, n.LawU)
weights.LawU[LawU[,1]>mean(LawU[,1])] <- 1

LawF <- filter(diagnosed, SITE == "Law-F") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.LawF <- nrow(LawF) ##n.LawFumber of location.LawFs

phi <- 3/.5
sigma.sq <- 2

R.LawF <- sigma.sq*exp(-phi*as.matrix(dist(LawF)))
w.LawF <- rmvn(1, rep(0,n.LawF), R.LawF)

x.LawF <- as.matrix(rep(1,n.LawF))
beta <- 0.1
p.LawF <- 1/(1+exp(-(x.LawF%*%beta+w.LawF)))

weights.LawF <- rep(1, n.LawF)
weights.LawF[LawF[,1]>mean(LawF[,1])] <- 1

Museum <- filter(diagnosed, SITE == "Museum") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.Museum <- nrow(Museum) ##n.Museumumber of location.Museums

phi <- 3/.5
sigma.sq <- 2

R.Museum <- sigma.sq*exp(-phi*as.matrix(dist(Museum)))
w.Museum <- rmvn(1, rep(0,n.Museum), R.Museum)

x.Museum <- as.matrix(rep(1,n.Museum))
beta <- 0.1
p.Museum <- 1/(1+exp(-(x.Museum%*%beta+w.Museum)))

weights.Museum <- rep(1, n.Museum)
weights.Museum[Museum[,1]>mean(Museum[,1])] <- 1

McCarty <- filter(diagnosed, SITE == "McCarty") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.McCarty <- nrow(McCarty) ##n.McCartyumber of location.McCartys

phi <- 3/.5
sigma.sq <- 2

R.McCarty <- sigma.sq*exp(-phi*as.matrix(dist(McCarty)))
w.McCarty <- rmvn(1, rep(0,n.McCarty), R.McCarty)

x.McCarty <- as.matrix(rep(1,n.McCarty))
beta <- 0.1
p.McCarty <- 1/(1+exp(-(x.McCarty%*%beta+w.McCarty)))

weights.McCarty <- rep(1, n.McCarty)
weights.McCarty[McCarty[,1]>mean(McCarty[,1])] <- 1

Ficke <- filter(diagnosed, SITE == "Ficke") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.Ficke <- nrow(Ficke) ##n.Fickeumber of location.Fickes

phi <- 3/.5
sigma.sq <- 2

R.Ficke <- sigma.sq*exp(-phi*as.matrix(dist(Ficke)))
w.Ficke <- rmvn(1, rep(0,n.Ficke), R.Ficke)

x.Ficke <- as.matrix(rep(1,n.Ficke))
beta <- 0.1
p.Ficke <- 1/(1+exp(-(x.Ficke%*%beta+w.Ficke)))

weights.Ficke <- rep(1, n.Ficke)
weights.Ficke[Ficke[,1]>mean(Ficke[,1])] <- 1

Alice <- filter(diagnosed, SITE == "Alice") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.Alice <- nrow(Alice) ##n.Aliceumber of location.Alices

phi <- 3/.5
sigma.sq <- 2

R.Alice <- sigma.sq*exp(-phi*as.matrix(dist(Alice)))
w.Alice <- rmvn(1, rep(0,n.Alice), R.Alice)

x.Alice <- as.matrix(rep(1,n.Alice))
beta <- 0.1
p.Alice <- 1/(1+exp(-(x.Alice%*%beta+w.Alice)))

weights.Alice <- rep(1, n.Alice)
weights.Alice[Alice[,1]>mean(Alice[,1])] <- 1

Entom <- filter(diagnosed, SITE == "Entom") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.Entom <- nrow(Entom) ##n.Entomumber of location.Entoms

phi <- 3/.5
sigma.sq <- 2

R.Entom <- sigma.sq*exp(-phi*as.matrix(dist(Entom)))
w.Entom <- rmvn(1, rep(0,n.Entom), R.Entom)

x.Entom <- as.matrix(rep(1,n.Entom))
beta <- 0.1
p.Entom <- 1/(1+exp(-(x.Entom%*%beta+w.Entom)))

weights.Entom <- rep(1, n.Entom)
weights.Entom[Entom[,1]>mean(Entom[,1])] <- 1

NATL <- filter(diagnosed, SITE == "NATL") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.NATL <- nrow(NATL) ##n.NATLumber of location.NATLs

phi <- 3/.5
sigma.sq <- 2

R.NATL <- sigma.sq*exp(-phi*as.matrix(dist(NATL)))
w.NATL <- rmvn(1, rep(0,n.NATL), R.NATL)

x.NATL <- as.matrix(rep(1,n.NATL))
beta <- 0.1
p.NATL <- 1/(1+exp(-(x.NATL%*%beta+w.NATL)))

weights.NATL <- rep(1, n.NATL)
weights.NATL[NATL[,1]>mean(NATL[,1])] <- 1

EnvU <- filter(diagnosed, SITE == "Env. U") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.EnvU <- nrow(EnvU) ##n.EnvUumber of location.EnvUs

phi <- 3/.5
sigma.sq <- 2

R.EnvU <- sigma.sq*exp(-phi*as.matrix(dist(EnvU)))
w.EnvU <- rmvn(1, rep(0,n.EnvU), R.EnvU)

x.EnvU <- as.matrix(rep(1,n.EnvU))
beta <- 0.1
p.EnvU <- 1/(1+exp(-(x.EnvU%*%beta+w.EnvU)))

weights.EnvU <- rep(1, n.EnvU)
weights.EnvU[EnvU[,1]>mean(EnvU[,1])] <- 1

EnvF <- filter(diagnosed, SITE == "Env. F") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.EnvF <- nrow(EnvF) ##n.EnvFumber of location.EnvFs

phi <- 3/.5
sigma.sq <- 2

R.EnvF <- sigma.sq*exp(-phi*as.matrix(dist(EnvF)))
w.EnvF <- rmvn(1, rep(0,n.EnvF), R.EnvF)

x.EnvF <- as.matrix(rep(1,n.EnvF))
beta <- 0.1
p.EnvF <- 1/(1+exp(-(x.EnvF%*%beta+w.EnvF)))

weights.EnvF <- rep(1, n.EnvF)
weights.EnvF[EnvF[,1]>mean(EnvF[,1])] <- 1

btpU <- filter(diagnosed, SITE == "BTP U") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.btpU <- nrow(btpU) ##n.btpUumber of location.btpUs

phi <- 3/.5
sigma.sq <- 2

R.btpU <- sigma.sq*exp(-phi*as.matrix(dist(btpU)))
w.btpU <- rmvn(1, rep(0,n.btpU), R.btpU)

x.btpU <- as.matrix(rep(1,n.btpU))
beta <- 0.1
p.btpU <- 1/(1+exp(-(x.btpU%*%beta+w.btpU)))

weights.btpU <- rep(1, n.btpU)
weights.btpU[btpU[,1]>mean(btpU[,1])] <- 1

btpF <- filter(diagnosed, SITE == "BTP F") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.btpF <- nrow(btpF) ##n.btpFumber of location.btpFs

phi <- 3/.5
sigma.sq <- 2

R.btpF <- sigma.sq*exp(-phi*as.matrix(dist(btpF)))
w.btpF <- rmvn(1, rep(0,n.btpF), R.btpF)

x.btpF <- as.matrix(rep(1,n.btpF))
beta <- 0.1
p.btpF <- 1/(1+exp(-(x.btpF%*%beta+w.btpF)))

weights.btpF <- rep(1, n.btpF)
weights.btpF[btpF[,1]>mean(btpF[,1])] <- 1

KewU <- filter(diagnosed, SITE == "Kew U") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.KewU <- nrow(KewU) ##n.KewUumber of location.KewUs

phi <- 3/.5
sigma.sq <- 2

R.KewU <- sigma.sq*exp(-phi*as.matrix(dist(KewU)))
w.KewU <- rmvn(1, rep(0,n.KewU), R.KewU)

x.KewU <- as.matrix(rep(1,n.KewU))
beta <- 0.1
p.KewU <- 1/(1+exp(-(x.KewU%*%beta+w.KewU)))

weights.KewU <- rep(1, n.KewU)
weights.KewU[KewU[,1]>mean(KewU[,1])] <- 1

KewF <- filter(diagnosed, SITE == "Kew F") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
n.KewF <- nrow(KewF) ##n.KewFumber of location.KewFs

phi <- 3/.5
sigma.sq <- 2

R.KewF <- sigma.sq*exp(-phi*as.matrix(dist(KewF)))
w.KewF <- rmvn(1, rep(0,n.KewF), R.KewF)

x.KewF <- as.matrix(rep(1,n.KewF))
beta <- 0.1
p.KewF <- 1/(1+exp(-(x.KewF%*%beta+w.KewF)))

weights.KewF <- rep(1, n.KewF)
weights.KewF[KewF[,1]>mean(KewF[,1])] <- 1




y.LawU <- filter(diagnosed, SITE == "Law-U")$infection
y.LawF <- filter(diagnosed, SITE == "Law-F")$infection
y.Museum <- filter(diagnosed, SITE == "Museum")$infection
y.McCarty <- filter(diagnosed, SITE == "McCarty")$infection
y.Ficke <- filter(diagnosed, SITE == "Ficke")$infection
y.Alice <- filter(diagnosed, SITE == "Alice")$infection
y.Entom <- filter(diagnosed, SITE == "Entom")$infection
y.NATL <- filter(diagnosed, SITE == "NATL")$infection
y.EnvU <- filter(diagnosed, SITE == "Env. U")$infection
y.EnvF <- filter(diagnosed, SITE == "Env. F")$infection
y.btpU <- filter(diagnosed, SITE == "BTP U")$infection
y.btpF <- filter(diagnosed, SITE == "BTP F")$infection
y.KewU <- filter(diagnosed, SITE == "Kew U")$infection
y.KewF <- filter(diagnosed, SITE == "Kew F")$infection


sex.LawU <- as.numeric(filter(diagnosed, SITE == "Law-U")$SEX)
sex.LawF <- as.numeric(filter(diagnosed, SITE == "Law-F")$SEX)
sex.Museum <- as.numeric(filter(diagnosed, SITE == "Museum")$SEX)
sex.McCarty <- as.numeric(filter(diagnosed, SITE == "McCarty")$SEX)
sex.Ficke <- as.numeric(filter(diagnosed, SITE == "Ficke")$SEX)
sex.Alice <- as.numeric(filter(diagnosed, SITE == "Alice")$SEX)
sex.Entom <- as.numeric(filter(diagnosed, SITE == "Entom")$SEX)
sex.NATL <- as.numeric(filter(diagnosed, SITE == "NATL")$SEX)
sex.EnvU <- as.numeric(filter(diagnosed, SITE == "Env. U")$SEX)
sex.EnvF <- as.numeric(filter(diagnosed, SITE == "Env. F")$SEX)
sex.btpU <- as.numeric(filter(diagnosed, SITE == "BTP U")$SEX)
sex.btpF <- as.numeric(filter(diagnosed, SITE == "BTP F")$SEX)
sex.KewU <- as.numeric(filter(diagnosed, SITE == "Kew U")$SEX)
sex.KewF <- as.numeric(filter(diagnosed, SITE == "Kew F")$SEX)


svl.LawU <- filter(diagnosed, SITE == "Law-U")$SVL
svl.LawF <- filter(diagnosed, SITE == "Law-F")$SVL
svl.Museum <- filter(diagnosed, SITE == "Museum")$SVL
svl.McCarty <- filter(diagnosed, SITE == "McCarty")$SVL
svl.Ficke <- filter(diagnosed, SITE == "Ficke")$SVL
svl.Alice <- filter(diagnosed, SITE == "Alice")$SVL
svl.Entom <- filter(diagnosed, SITE == "Entom")$SVL
svl.NATL <- filter(diagnosed, SITE == "NATL")$SVL
svl.EnvU <- filter(diagnosed, SITE == "Env. U")$SVL
svl.EnvF <- filter(diagnosed, SITE == "Env. F")$SVL
svl.btpU <- filter(diagnosed, SITE == "BTP U")$SVL
svl.btpF <- filter(diagnosed, SITE == "BTP F")$SVL
svl.KewU <- filter(diagnosed, SITE == "Kew U")$SVL
svl.KewF <- filter(diagnosed, SITE == "Kew F")$SVL



####p(infection)~structure####
####Law U
fit <- glm((y.LawU/weights.LawU)~x.LawU-1, weights=weights.LawU, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.LawU~1, family="binomial", coords=LawU, weights=weights.LawU, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.LawU%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.LawU, size=weights.LawU, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(LawU,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Law U: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(LawU, label=paste("(",y.LawU,")",sep=""))

surf <- mba.surf(cbind(LawU,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Law U: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(LawU, label=paste("(",weights.LawU,")",sep=""))

####Law F
fit <- glm((y.LawF/weights.LawF)~x.LawF-1, weights=weights.LawF, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.LawF~1, family="binomial", coords=LawF, weights=weights.LawF, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.LawF%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.LawF, size=weights.LawF, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(LawF,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Law F: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(LawF, label=paste("(",y.LawF,")",sep=""))

surf <- mba.surf(cbind(LawF,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Law F: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(LawF, label=paste("(",weights.LawF,")",sep=""))

####Museum
fit <- glm((y.Museum/weights.Museum)~x.Museum-1, weights=weights.Museum, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.Museum~1, family="binomial", coords=Museum, weights=weights.Museum, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.Museum%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.Museum, size=weights.Museum, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(Museum,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Museum: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(Museum, label=paste("(",y.Museum,")",sep=""))

surf <- mba.surf(cbind(Museum,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Museum: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(Museum, label=paste("(",weights.Museum,")",sep=""))

####McCarty
fit <- glm((y.McCarty/weights.McCarty)~x.McCarty-1, weights=weights.McCarty, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.McCarty~1, family="binomial", coords=McCarty, weights=weights.McCarty, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.McCarty%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.McCarty, size=weights.McCarty, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(McCarty,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="McCarty: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(McCarty, label=paste("(",y.McCarty,")",sep=""))

surf <- mba.surf(cbind(McCarty,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="McCarty: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(McCarty, label=paste("(",weights.McCarty,")",sep=""))

####Ficke
fit <- glm((y.Ficke/weights.Ficke)~x.Ficke-1, weights=weights.Ficke, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.Ficke~1, family="binomial", coords=Ficke, weights=weights.Ficke, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.Ficke%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.Ficke, size=weights.Ficke, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(Ficke,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Ficke: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(Ficke, label=paste("(",y.Ficke,")",sep=""))

surf <- mba.surf(cbind(Ficke,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Ficke: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(Ficke, label=paste("(",weights.Ficke,")",sep=""))

####Alice
fit <- glm((y.Alice/weights.Alice)~x.Alice-1, weights=weights.Alice, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.Alice~1, family="binomial", coords=Alice, weights=weights.Alice, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.Alice%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.Alice, size=weights.Alice, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(Alice,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Alice: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(Alice, label=paste("(",y.Alice,")",sep=""))

surf <- mba.surf(cbind(Alice,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Alice: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(Alice, label=paste("(",weights.Alice,")",sep=""))

####Entom
fit <- glm((y.Entom/weights.Entom)~x.Entom-1, weights=weights.Entom, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.Entom~1, family="binomial", coords=Entom, weights=weights.Entom, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.Entom%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.Entom, size=weights.Entom, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(Entom,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Entom: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(Entom, label=paste("(",y.Entom,")",sep=""))

surf <- mba.surf(cbind(Entom,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Entom: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(Entom, label=paste("(",weights.Entom,")",sep=""))

####NATL
fit <- glm((y.NATL/weights.NATL)~x.NATL-1, weights=weights.NATL, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.NATL~1, family="binomial", coords=NATL, weights=weights.NATL, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.NATL%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.NATL, size=weights.NATL, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(NATL,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="NATL: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(NATL, label=paste("(",y.NATL,")",sep=""))

surf <- mba.surf(cbind(NATL,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="NATL: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(NATL, label=paste("(",weights.NATL,")",sep=""))

####Env. U
fit <- glm((y.EnvU/weights.EnvU)~x.EnvU-1, weights=weights.EnvU, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 100
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.EnvU~1, family="binomial", coords=EnvU, weights=weights.EnvU, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.EnvU%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.EnvU, size=weights.EnvU, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(EnvU,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Env U: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(EnvU, label=paste("(",y.EnvU,")",sep=""))

surf <- mba.surf(cbind(EnvU,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Env U: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(EnvU, label=paste("(",weights.EnvU,")",sep=""))

####Env. F
fit <- glm((y.EnvF/weights.EnvF)~x.EnvF-1, weights=weights.EnvF, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.EnvF~1, family="binomial", coords=EnvF, weights=weights.EnvF, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.EnvF%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.EnvF, size=weights.EnvF, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(EnvF,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Env F: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(EnvF, label=paste("(",y.EnvF,")",sep=""))

surf <- mba.surf(cbind(EnvF,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Env F: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(EnvF, label=paste("(",weights.EnvF,")",sep=""))

####BTP U
fit <- glm((y.btpU/weights.btpU)~x.btpU-1, weights=weights.btpU, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.btpU~1, family="binomial", coords=btpU, weights=weights.btpU, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.btpU%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.btpU, size=weights.btpU, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(btpU,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="BTP U: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(btpU, label=paste("(",y.btpU,")",sep=""))

surf <- mba.surf(cbind(btpU,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="BTP U: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(btpU, label=paste("(",weights.btpU,")",sep=""))

####BTP F
fit <- glm((y.btpF/weights.btpF)~x.btpF-1, weights=weights.btpF, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.btpF~1, family="binomial", coords=btpF, weights=weights.btpF, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.btpF%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.btpF, size=weights.btpF, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(btpF,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="BTP F: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(btpF, label=paste("(",y.btpF,")",sep=""))

surf <- mba.surf(cbind(btpF,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="BTP F: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(btpF, label=paste("(",weights.btpF,")",sep=""))

####Kew U
fit <- glm((y.KewU/weights.KewU)~x.KewU-1, weights=weights.KewU, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.KewU~1, family="binomial", coords=KewU, weights=weights.KewU, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.KewU%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.KewU, size=weights.KewU, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(KewU,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Kew U: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(KewU, label=paste("(",y.KewU,")",sep=""))

surf <- mba.surf(cbind(KewU,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Kew U: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(KewU, label=paste("(",weights.KewU,")",sep=""))

####Kew F
fit <- glm((y.KewF/weights.KewF)~x.KewF-1, weights=weights.KewF, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.KewF~1, family="binomial", coords=KewF, weights=weights.KewF, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.KewF%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.KewF, size=weights.KewF, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(KewF,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Kew F: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(KewF, label=paste("(",y.KewF,")",sep=""))

surf <- mba.surf(cbind(KewF,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Kew F: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(KewF, label=paste("(",weights.KewF,")",sep=""))

####p(infection)~Structure+SVL####

####Law U
fit <- glm((y.LawU/weights.LawU)~x.LawU-1, weights=weights.LawU, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.LawU~1, family="binomial", coords=LawU, weights=weights.LawU, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.LawU%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.LawU, size=weights.LawU, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(LawU,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Law U: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(LawU, label=paste("(",y.LawU,")",sep=""))

surf <- mba.surf(cbind(LawU,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Law U: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(LawU, label=paste("(",weights.LawU,")",sep=""))

####Law F
fit <- glm((y.LawF/weights.LawF)~x.LawF-1, weights=weights.LawF, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.LawF~1, family="binomial", coords=LawF, weights=weights.LawF, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.LawF%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.LawF, size=weights.LawF, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(LawF,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Law F: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(LawF, label=paste("(",y.LawF,")",sep=""))

surf <- mba.surf(cbind(LawF,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Law F: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(LawF, label=paste("(",weights.LawF,")",sep=""))

####Museum
fit <- glm((y.Museum/weights.Museum)~x.Museum-1, weights=weights.Museum, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.Museum~1, family="binomial", coords=Museum, weights=weights.Museum, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.Museum%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.Museum, size=weights.Museum, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(Museum,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Museum: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(Museum, label=paste("(",y.Museum,")",sep=""))

surf <- mba.surf(cbind(Museum,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Museum: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(Museum, label=paste("(",weights.Museum,")",sep=""))

####McCarty
fit <- glm((y.McCarty/weights.McCarty)~x.McCarty-1, weights=weights.McCarty, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.McCarty~1, family="binomial", coords=McCarty, weights=weights.McCarty, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.McCarty%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.McCarty, size=weights.McCarty, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(McCarty,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="McCarty: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(McCarty, label=paste("(",y.McCarty,")",sep=""))

surf <- mba.surf(cbind(McCarty,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="McCarty: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(McCarty, label=paste("(",weights.McCarty,")",sep=""))

####Ficke
fit <- glm((y.Ficke/weights.Ficke)~x.Ficke-1, weights=weights.Ficke, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.Ficke~1, family="binomial", coords=Ficke, weights=weights.Ficke, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.Ficke%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.Ficke, size=weights.Ficke, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(Ficke,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Ficke: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(Ficke, label=paste("(",y.Ficke,")",sep=""))

surf <- mba.surf(cbind(Ficke,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Ficke: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(Ficke, label=paste("(",weights.Ficke,")",sep=""))

####Alice
fit <- glm((y.Alice/weights.Alice)~x.Alice-1, weights=weights.Alice, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.Alice~1, family="binomial", coords=Alice, weights=weights.Alice, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.Alice%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.Alice, size=weights.Alice, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(Alice,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Alice: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(Alice, label=paste("(",y.Alice,")",sep=""))

surf <- mba.surf(cbind(Alice,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Alice: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(Alice, label=paste("(",weights.Alice,")",sep=""))

####Entom
fit <- glm((y.Entom/weights.Entom)~x.Entom-1, weights=weights.Entom, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.Entom~1, family="binomial", coords=Entom, weights=weights.Entom, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.Entom%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.Entom, size=weights.Entom, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(Entom,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Entom: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(Entom, label=paste("(",y.Entom,")",sep=""))

surf <- mba.surf(cbind(Entom,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Entom: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(Entom, label=paste("(",weights.Entom,")",sep=""))

####NATL
fit <- glm((y.NATL/weights.NATL)~x.NATL-1, weights=weights.NATL, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.NATL~1, family="binomial", coords=NATL, weights=weights.NATL, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.NATL%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.NATL, size=weights.NATL, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(NATL,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="NATL: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(NATL, label=paste("(",y.NATL,")",sep=""))

surf <- mba.surf(cbind(NATL,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="NATL: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(NATL, label=paste("(",weights.NATL,")",sep=""))

####Env. U
fit <- glm((y.EnvU/weights.EnvU)~x.EnvU-1, weights=weights.EnvU, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.EnvU~1, family="binomial", coords=EnvU, weights=weights.EnvU, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.EnvU%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.EnvU, size=weights.EnvU, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(EnvU,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Env U: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(EnvU, label=paste("(",y.EnvU,")",sep=""))

surf <- mba.surf(cbind(EnvU,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Env U: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(EnvU, label=paste("(",weights.EnvU,")",sep=""))

####Env. F
fit <- glm((y.EnvF/weights.EnvF)~x.EnvF-1, weights=weights.EnvF, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.EnvF~1, family="binomial", coords=EnvF, weights=weights.EnvF, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.EnvF%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.EnvF, size=weights.EnvF, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(EnvF,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Env F: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(EnvF, label=paste("(",y.EnvF,")",sep=""))

surf <- mba.surf(cbind(EnvF,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Env F: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(EnvF, label=paste("(",weights.EnvF,")",sep=""))

####BTP U
fit <- glm((y.btpU/weights.btpU)~x.btpU-1, weights=weights.btpU, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.btpU~1, family="binomial", coords=btpU, weights=weights.btpU, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.btpU%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.btpU, size=weights.btpU, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(btpU,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="BTP U: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(btpU, label=paste("(",y.btpU,")",sep=""))

surf <- mba.surf(cbind(btpU,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="BTP U: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(btpU, label=paste("(",weights.btpU,")",sep=""))

####BTP F
fit <- glm((y.btpF/weights.btpF)~x.btpF-1, weights=weights.btpF, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.btpF~1, family="binomial", coords=btpF, weights=weights.btpF, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.btpF%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.btpF, size=weights.btpF, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(btpF,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="BTP F: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(btpF, label=paste("(",y.btpF,")",sep=""))

surf <- mba.surf(cbind(btpF,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="BTP F: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(btpF, label=paste("(",weights.btpF,")",sep=""))

####Kew U
fit <- glm((y.KewU/weights.KewU)~x.KewU-1, weights=weights.KewU, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.KewU~1, family="binomial", coords=KewU, weights=weights.KewU, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.KewU%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.KewU, size=weights.KewU, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(KewU,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Kew U: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(KewU, label=paste("(",y.KewU,")",sep=""))

surf <- mba.surf(cbind(KewU,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Kew U: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(KewU, label=paste("(",weights.KewU,")",sep=""))

####Kew F
fit <- glm((y.KewF/weights.KewF)~x.KewF-1, weights=weights.KewF, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y.KewF~1, family="binomial", coords=KewF, weights=weights.KewF, 
             starting=list("beta"=beta.starting, "phi"=3/0.5,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=1, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(3/0.75, 3/0.25), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x.KewF%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n.KewF, size=weights.KewF, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(KewF,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Kew F: Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(KewF, label=paste("(",y.KewF,")",sep=""))

surf <- mba.surf(cbind(KewF,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Kew F: Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(KewF, label=paste("(",weights.KewF,")",sep=""))