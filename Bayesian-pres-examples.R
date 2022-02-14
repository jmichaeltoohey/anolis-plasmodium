### BAYESIAN INFERENCE ###

###Simulate seed data####

#simulate some poisson data... say, Average Number of Seeds in a Quadrant
set.seed(3) #So that everybody has the same random numbers
data=rpois(50, lambda=6.4) #Average of 6.4 seeds per quadrant
hist(data)

#This will be our Prior knowledge, which we will update with new data 
#How do you choose a Prior based on previous knowledge?...

###Priors p(PARAMETERS)####

#Say, last year we measured 10.5 +/- 0.5 seeds on average
#Range was from 0-15

#How do we choose a Prior distribution to incorporate this information?
#Average density of seeds: continuous and >= 0... GAMMA
#mu (average) = 10.5, sd = 0.5
#Gamma uses shape (mu^2/sd^2) and scale (mu/sd^2)

#The function will take the average and sd and estimate the probability of all values of interest
theta=seq(0,15,.01)
mu.prior=10.5
sigma.prior=0.5
Prior=function(theta, mu=mu.prior, sigma=sigma.prior){
  dgamma(theta,mu^2/sigma^2,mu/sigma^2)
}
#Make a figure to visualize how the distribution looks
plot(theta,Prior(theta), typ="l", ylab=expression(paste("P(",theta,")")),
     xlab=expression(theta),main="Prior", xlim=c(5,15))

###Likelihood p(DATA|PARAMETERS)####

#Remember: Likelihood is probability of the data (y) given the parameters (theta)

#What distribution will we use? 
#We know it is discrete, positive data (number of seeds in quadrant)
#Poisson sounds right... are variance and mean equal:
var(data) #4.057551
mean(data) #5.94
#close enough!

###Build likelihood function

#Remember, probability mass function for Poisson is given by dpois()
#To calculate the likelihood of this data, we just need the product for all the values in the data

Like=function(theta,y=data){
  L=rep(0,length(theta))
  for(i in 1:length(theta)) L[i] = prod(dpois(y,theta[i], log=FALSE))
  return(L)
} 

#Note: takes in theta as input, which are all potential values for the parameter.
#For loop is calculating the likelihood, or probability of data for all potential values of theta

#Visualize likelihood
plot(theta,Like(theta), typ="l", xlab=expression(theta),
     ylab=expression(paste("P(y | ",theta,")")), main="Likelihood",xlim=c(5,15))
#The mean of this distribution is much less than that of the Priors

###Joint Distribution####

#Multiply the two probability distributions
Joint=function(theta) {
  Like(theta)*Prior(theta)
}
# function simply takes all the potential values of the parameter, theta and runs them through
# the two functions that we already built Like and Prior and multiplies them

#Visualize
plot(theta,Joint(theta),typ="l",
     ylab=expression(paste("P(y | ",theta,") x P(",theta,")" )),
     main="Joint", xlab=expression(theta), xlim=c(5,15) )
#Mean is inbetween likelihood and prior distribution, but probability is very small

###Marginal Probability####

#Simply the area under the curve of the joint distribution
#Just integrate it - it will be just a scalar number

Py=integrate(Joint,0,30,abs.tol=1e-100)
Py #Note that Py is a very small scalar

#integrate() just takes as input the function that we want to integrate, the range of
#values to evaluate the integral and the absolute tolerance, which is simply the amount of error that we can
#tolerate

###Posterior p(PARAMETERS|DATA)####

P.theta = Joint(theta)/Py$value
plot(theta,P.theta,typ="l", xlab=expression(theta),
     ylab=expression(paste("P( ",theta," | y)")),xlim=c(5,15), main="Posterior")

#Plot all probability distributions in one figure
like.v=Like(theta)
like.scaled = like.v/max(like.v)*max(P.theta)
#pdf("Scaled_overlay.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(theta,like.scaled,typ="l",col="red",xlim=c(5,15),
     xlab=expression(theta),ylab= "Probability density",
     main="Scaled Overlay")
lines(theta,Prior(theta),col="blue")
lines(theta,P.theta,col="black", lwd=4, typ="l", lty="dashed")
legend(11,1, legend= c("Scaled Likelihood", "Prior","Integrated posterior"),
       cex=.8,lwd=2, bty="n", col=c("red","blue","black"),
       lty=c("solid","solid","dashed"))

###JAGS####

#Simulate new data
N <- 1000
x <- rnorm(N,2,10)
epsilon <- rnorm(N, 0, 3)
y <- x + epsilon

#Use linear model to recuperate parameters simulated
lm(y~x)

#NOTE: JAGS USES PRECISION INSTEAD OF VARIANCE, 
#so, linear models look like:
# y ~ N(mu, 1/sd^2) [NOT] N(mu, sd^2)
# mu = beta1 + beta2*x
# beta1 and beta2 are parameters, x and y are data
# priors will be vague, so beta_i ~ N(0, 0.00001)

library("rjags")
model_string <- "model{
# Prior for beta
for(j in 1:2){
beta[j] ~ dnorm(0,0.0001)
}
# Likelihood
for(i in 1:n){
Y[i] ~ dnorm(mu[i],inv.var)
mu[i] <- beta[1] + beta[2]*x[i]
}
# Prior for the inverse variance (precision)
inv.var ~ dunif(0, 100)
sigma <- 1/sqrt(inv.var)
}"

#rjags is used for R to communicate with JAGS

#We store the model into a variable we call model_string - lets break it down

#Priors
#Vague, so the distribution should be flat
#Large variances mean very small precisions (we use 0.0001 to give us a large variance)

#Likelihood
#Just the state space notation we have in previous equation
#we had set the distribution of the data as normal with deterministic function that describes a linear model



#Now we build the model using jags.model()
#10000 iterations of burn-in (first draws of the MCMC where the algorithm hasnt converged yet)
model <- jags.model(textConnection(model_string),
                    data = list(Y=y,n=N,x=x))
update(model, 10000, progress.bar="none") # Burn-in for 10000 samples


#Run the model using coda.samples()
#one chain of 20000
samp <- coda.samples(model,
                     variable.names=c("beta","sigma"),
                     n.iter=20000, progress.bar="none")
summary(samp)
plot(samp)
