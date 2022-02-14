#### BAYESIAN FUNCTIONS ####
#For Visualizing the Output of Limited Dependent variables 

#Quick look at the Distribution of data
library(ggplot2)
glance <- function(data){
  require(ggplot2)
  plot <- ggplot(data = melt(data), aes(x = factor(value))) + 
    geom_histogram(position = "identity") + 
    facet_wrap(~ variable, scales = "free", as.table = TRUE, nrow = 3) + xlab("") + theme_bw()
  return(plot)
}


#Function to draw out the parts of the MCMC output that we care about. 
msum <- function(x){
  require(coda)
  f <- summary(as.mcmc(x))
  out <- cbind(f$statistics[,c(1,2)],f$quantiles[,c(1,5)])
  return(out)
}


# Predicted Probabilities
#Using Observed Values Approach (Hanmer Kalkan 2013)
pprobs <- function(mcmc.output=NULL,model.data=NULL,params=NULL,
                   manipulation.var=NULL,dv=NULL,link="logit",
                   type="smoother",points=F){
  require(stringr)
  require(ggplot2)
  require(coda)
  require(reshape2)
  require(dplyr)
  dd <- data.frame(model.data)
  dd[,dv] <- NULL #remove DV
  dd$N <- NULL
  dd <- cbind(data.frame(constant=1),dd)
  value.range <- dd[,str_detect(colnames(dd),manipulation.var)] 
  value.range <- sort(unique(value.range))
  value.range <- rev(value.range)
  
  run <- function(value.range){
    dd[,str_detect(colnames(dd),manipulation.var)] <- value.range
    dd.mat <- as.matrix(dd)
    if(class(mcmc.output)=="rjags"){
      c <- as.matrix(as.mcmc(mcmc.output))
    }
    if(class(mcmc.output)=="bugs"){
      c <- mcmc.output$sims.matrix
    }
    params <- c("alpha","b")
    coefs = NULL
    for(i in 1:length(params)){
      coefs <- cbind(coefs,c[,str_detect(colnames(c),pattern=params[i])])
    }
    Xb <- t(dd.mat %*% t(coefs))
    if(link=="logit"){
      pprobs <- exp(Xb) / (1 + exp(Xb))
    }
    if(link=="probit"){
      pprobs <- pnorm(Xb)
    }
    pp <- melt(pprobs, varnames = c("Iteration", manipulation.var))
    pp.sum <- summarise(pp, mean.pp = mean(value), 
                        lower.pp = quantile(value, probs = c(0.05)), 
                        upper.pp = quantile(value, probs = c(0.95)))
  }
  gathered.pp.sum <- plyr::ldply(value.range,function(x) run(x))
  gathered.pp.sum$value.range <- value.range
  if(type=="smoother"){
    plot <- ggplot(data = gathered.pp.sum, aes(x = value.range, y = mean.pp)) + 
      geom_ribbon(aes(ymin = lower.pp, ymax = upper.pp), alpha = 0.2) + 
      geom_line() + theme_bw() + ylab(paste0("Pr(",dv,")")) + xlab(manipulation.var) 
    if(points){
      plot <- plot + geom_point(x=value.range,size=3)
    }
    plot <- plot + scale_x_continuous(breaks=value.range,minor_breaks=value.range)
  }
  if(type=="discrete"){
    plot <-  ggplot(data = gathered.pp.sum, aes(x = value.range, y = mean.pp)) + 
      geom_pointrange(aes(ymin = lower.pp, ymax = upper.pp)) + theme_bw() + 
      ylab(paste0("Pr(",dv,")")) + xlab(manipulation.var) 
    plot <- plot + scale_x_continuous(breaks=value.range,minor_breaks=value.range)
  }
  return(plot)
}


#Marginal Effects Plot
#Using Observed Values Approach (Hanmer Kalkan 2013)
marginFX <- function(mcmc.output=NULL,model.data=NULL,params=NULL,
                     manipulation.var=NULL,dv=NULL,link="logit",
                     type="smoother",points=F){
  require(data.table)
  require(stringr)
  require(ggplot2)
  require(coda)
  require(reshape2)
  require(dplyr)
  dd <- data.frame(model.data)
  dd[,dv] <- NULL #remove DV
  dd$N <- NULL
  dd <- cbind(data.frame(constant=1),dd)
  value.range <- dd[,str_detect(colnames(dd),manipulation.var)] 
  value.range <- sort(unique(value.range))
  value.range <- rev(value.range)
  
  #Creating a baseline
  dd[,str_detect(colnames(dd),manipulation.var)] <- min(value.range)
  dd.mat <- as.matrix(dd)
  if(class(mcmc.output)=="rjags"){
    c <- as.matrix(as.mcmc(mcmc.output))
  }
  if(class(mcmc.output)=="bugs"){
    c <- mcmc.output$sims.matrix
  }
  params <- c("alpha","b")
  coefs = NULL
  for(i in 1:length(params)){
    coefs <- cbind(coefs,c[,str_detect(colnames(c),pattern=params[i])])
  }
  Xb <- t(dd.mat %*% t(coefs))
  if(link=="logit"){
    pprobs <- exp(Xb) / (1 + exp(Xb))
  }
  if(link=="probit"){
    pprobs <- pnorm(Xb)
  }
  pp <- melt(pprobs, varnames = c("Iteration", manipulation.var))
  pp.baseline <- data.frame(Baseline = pp, Iteration = 1:nrow(pp))
  pp.baseline <- data.table(pp.baseline,key="Iteration")
  
  gathered.sum <- NULL
  run <- function(value.range){
    dd[,str_detect(colnames(dd),manipulation.var)] <- value.range
    dd.mat <- as.matrix(dd)
    if(class(mcmc.output)=="rjags"){
      c <- as.matrix(as.mcmc(mcmc.output))
    }
    if(class(mcmc.output)=="bugs"){
      c <- mcmc.output$sims.matrix
    }
    params <- c("alpha","b")
    coefs = NULL
    for(i in 1:length(params)){
      coefs <- cbind(coefs,c[,str_detect(colnames(c),pattern=params[i])])
    }
    Xb <- t(dd.mat %*% t(coefs))
    if(link=="logit"){
      pprobs <- exp(Xb) / (1 + exp(Xb))
    }
    if(link=="probit"){
      pprobs <- pnorm(Xb)
    }
    pp <- melt(pprobs, varnames = c("Iteration", manipulation.var))
    pp <- data.table(pp,key = "Iteration")
    together <-merge(pp,pp.baseline)
    together$me <- with(together, value - Baseline.value)
    together.me.sum <- summarise(together, mean.me = mean(me), 
                                 lower.me = quantile(me, probs = c(0.05)), upper.me = quantile(me, probs = c(0.95)))
    return(together.me.sum)
  }
  gathered.sum <- plyr::ldply(value.range,function(x) run(x))
  gathered.sum$value.range <- value.range
  if(type=="smoother"){
    plot <- ggplot(data = gathered.sum, aes(x = value.range, y = mean.me)) + 
      geom_ribbon(aes(ymin = lower.me, ymax = upper.me), alpha = 0.2) + 
      geom_line() + theme_bw() + ylab(paste0("Change in Pr(",dv,") from base probability")) + xlab(manipulation.var) 
    if(points){
      plot <- plot + geom_point(x=value.range,size=3)
    }
    plot <- plot + scale_x_continuous(breaks=value.range,minor_breaks=value.range)
  }
  if(type=="discrete"){
    plot <-  ggplot(data = gathered.sum, aes(x = value.range, y = mean.me)) + 
      geom_pointrange(aes(ymin = lower.me, ymax = upper.me)) + theme_bw() + 
      ylab(paste0("Change in Pr(",dv,") from base probability")) + xlab(manipulation.var) 
    plot <- plot + scale_x_continuous(breaks=value.range,minor_breaks=value.range)
  }
  return(plot)
}



### Function to measure the distribution of the discrete effect
#Using Observed Values Approach (Hanmer Kalkan 2013)
diffFX <- function(mcmc.output=NULL,model.data=NULL,params=NULL,
                   from=NULL,to=NULL,manipulation.var=NULL,dv=NULL,link="logit",type="distribution",add=F){
  require(stringr)
  require(ggplot2)
  require(coda)
  require(reshape2)
  require(dplyr)
  dd <- data.frame(model.data)
  dd[,dv] <- NULL #remove DV
  dd$N <- NULL
  dd <- cbind(data.frame(constant=1),dd)
  value.range <- dd[,str_detect(colnames(dd),manipulation.var)] 
  value.range <- sort(unique(value.range))
  value.range <- rev(value.range)
  
  run <- function(value.range){
    dd[,str_detect(colnames(dd),manipulation.var)] <- value.range
    dd.mat <- as.matrix(dd)
    if(class(mcmc.output)=="rjags"){
      c <- as.matrix(as.mcmc(mcmc.output))
    }
    if(class(mcmc.output)=="bugs"){
      c <- mcmc.output$sims.matrix
    }
    params <- c("alpha","b")
    coefs = NULL
    for(i in 1:length(params)){
      coefs <- cbind(coefs,c[,str_detect(colnames(c),pattern=params[i])])
    }
    Xb <- t(dd.mat %*% t(coefs))
    if(link=="logit"){
      pprobs <- exp(Xb) / (1 + exp(Xb))
    }
    if(link=="probit"){
      pprobs <- pnorm(Xb)
    }
    pp <- melt(pprobs, varnames = c("Iteration", manipulation.var))
    return(pp$value)
  }
  difference <- run(to) - run(from)
  difference <- data.frame(difference)
  effect <- summarise(difference, mean.pp = mean(difference), 
                      lower.pp = quantile(difference, probs = c(0.05)), 
                      upper.pp = quantile(difference, probs = c(0.95)))
  difference$ci <- ifelse(difference >= effect$lower.pp & difference <= effect$upper.pp,1,0)
  if(type=="distribution"){ 
    plot <- ggplot(data=difference,add=add) +
      geom_density(aes(x = difference,y=..density..),fill="skyblue",alpha=.3) + theme_bw() + 
      xlab(paste0("Discrete effect in Pr(",dv,") \n moving from ",from," to ",to," on ",manipulation.var)) + 
      geom_vline(xintercept=0,colour="red", linetype = "longdash")
    return(plot)
  }
  if(type=="ci.distribution"){
    hdrcde::hdr.den(difference$difference,main="",xlab=paste0("Discrete effect in Pr(",dv,") \n moving from ",from," to ",to," on ",manipulation.var)) 
    abline(v=0,col='red',lty=4,lwd=2)
  }  
}

message("\n Purpose: Functions for Visualizing the Baysian Output of Limited Dependent Variable Models \n \n  Author: Eric Dunford \n Affiliation: University of Maryland, College Park \n Date Created: July 2015")
