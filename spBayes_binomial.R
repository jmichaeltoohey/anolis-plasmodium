library(MBA)
library(coda)
library(spBayes)
library(dplyr)
library(tidyverse)
library(AICcmodavg)
library(loo)

set.seed(1)
source("data_carpentry.R")

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p)))
}



##Generate binary data
#coords <- as.matrix(expand.grid(seq(0,100,length.out=8), seq(0,100,length.out=8)))
coords <- filter(diagnosed, SITE == "Kew F") %>% select(x_coord, y_coord) %>% 
  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>% 
  na.omit() %>% as.matrix
#coords <- diagnosed %>% select(x_coord, y_coord) %>%
#  summarize(x_coord = -1*x_coord/180, y_coord = y_coord/90) %>%
#  na.omit() %>% as.matrix
n <- nrow(coords) ##number of locations

phi <- 3/50
sigma.sq <- 2

R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- rmvn(1, rep(0,n), R)

x <- as.matrix(rep(1,n))
beta <- 0.1
p <- 1/(1+exp(-(x%*%beta+w)))

weights <- rep(1, n)
weights[coords[,1]>mean(coords[,1])] <- 1

#y <- rbinom(n, size=weights, prob=p)
y <- filter(diagnosed, SITE == "Kew F")$infection
#y <- diagnosed %>% drop_na(x_coord)
#y <- y$infection

x2 <- as.numeric(filter(diagnosed, SITE == "Kew F")$SVL)



##Collect samples
fit <- glm((y/weights)~(x-1), weights=weights, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y~1, family="binomial", coords=coords, weights=weights, 
             #starting=list("beta"=beta1.starting, "beta2"=beta2.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             #tuning=list("beta"=beta1.tuning, "beta2"=beta2.tuning, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
             starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta1.Normal"=list(0,10), "beta2.Normal"=list(0,10),  "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)
m.1$p.beta.theta.samples


burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n, size=weights, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(coords, label=paste("(",y,")",sep=""))

surf <- mba.surf(cbind(coords,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(coords, label=paste("(",weights,")",sep=""))

spDiag(m.1)


#group_by(diagnosed, SITE) %>% summarize(sum_infected=length(SITE))
#group_by(diagnosed, SITE) %>% summarize(percent_infected=sum(infection)/length(SITE))
#group_by(diagnosed, SITE)

