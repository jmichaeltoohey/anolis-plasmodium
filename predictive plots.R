library(brms)
library(arm)
library(BayesPostEst)
devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")
#infection~ habitat####
meanfunc.inf.habitat <- function(habitat){
  return(exp(-6.7090+0.6300*habitat)/(1+exp(-6.7090+0.6300*habitat)))
}
lowerfunc.inf.habitat <- function(habitat){
  return(exp(-7.1130+0.1720*habitat)/(1+exp(-7.1130+0.1720*habitat)))
}
upperfunc.inf.habitat <- function(habitat){
  return(exp(-6.3560+1.1290*habitat)/(1+exp(-6.3560+1.1290*habitat)))
}
twentyfive.inf.habitat <- function(habitat){
  return(exp(-6.8390+0.4740*habitat)/(1+exp(-7.1130+0.1720*habitat)))
}
seventyfive.inf.habitat <- function(habitat){
  return(exp(-6.5910+0.7760*habitat)/(1+exp(-6.5910+0.7760*habitat)))
}

inf.hab.predict <- data.frame(
  x=c("Forested","Deforested"), y2.5=c(lowerfunc.inf.habitat(1),lowerfunc.inf.habitat(0)),
      y25=c(twentyfive.inf.habitat(1),twentyfive.inf.habitat(0)),
      y50=c(meanfunc.inf.habitat(1),meanfunc.inf.habitat(0)), 
      y75=c(seventyfive.inf.habitat(1),seventyfive.inf.habitat(0)),
      y97.5=c(upperfunc.inf.habitat(1),upperfunc.inf.habitat(0))
)

library(ggplot2)
plot.inf.hab <- ggplot(data=inf.hab.predict, aes(x))+
  geom_boxplot(
    aes(ymin = y2.5, lower = y25, middle = y50, upper = y75, ymax = y97.5), 
    stat = "identity"
  )+xlab("")+ylab("p(Infection)")

#infection~ svl sex habitat####

mod.infection.svl.sex.habitat.tab <- mcmctab(mod.infection.svl.sex.habitat)
meanfunc.inf.svl.sex.habitat <- function(SVL, SEX, habitat){
  return(invlogit(mod.infection.svl.sex.habitat.tab$Median[1]
                  + mod.infection.svl.sex.habitat.tab$Median[2]*SVL
                  + mod.infection.svl.sex.habitat.tab$Median[3]*SEX
                  + mod.infection.svl.sex.habitat.tab$Median[4]*habitat
                  + mod.infection.svl.sex.habitat.tab$Median[5]*SVL*SEX
                  + mod.infection.svl.sex.habitat.tab$Median[6]*SVL*habitat))
}
lowerfunc.inf.svl.sex.habitat <- function(SVL, SEX, habitat){
  return(invlogit(mod.infection.svl.sex.habitat.tab$Lower[1]
                  + mod.infection.svl.sex.habitat.tab$Lower[2]*SVL
                  + mod.infection.svl.sex.habitat.tab$Lower[3]*SEX
                  + mod.infection.svl.sex.habitat.tab$Lower[4]*habitat
                  + mod.infection.svl.sex.habitat.tab$Lower[5]*SVL*SEX
                  + mod.infection.svl.sex.habitat.tab$Lower[6]*SVL*habitat))
}
upperfunc.inf.svl.sex.habitat <- function(SVL, SEX, habitat){
  return(invlogit(mod.infection.svl.sex.habitat.tab$Upper[1]
                  + mod.infection.svl.sex.habitat.tab$Upper[2]*SVL
                  + mod.infection.svl.sex.habitat.tab$Upper[3]*SEX
                  + mod.infection.svl.sex.habitat.tab$Upper[4]*habitat
                  + mod.infection.svl.sex.habitat.tab$Upper[5]*SVL*SEX
                  + mod.infection.svl.sex.habitat.tab$Upper[6]*SVL*habitat))
}

p.infection.sim <- data.frame(svl = rep(seq(from=min(diagnosed$SVL),
                                            to=max(diagnosed$SVL),
                                            by=0.05),4))
p.infection.sim$habitat <- 0
p.infection.sim$habitat[1:(nrow(p.infection.sim)/2)] <- 1
p.infection.sim$sex <- 0
p.infection.sim$sex[c(1:(nrow(p.infection.sim)/4),(nrow(p.infection.sim)/2+1):(3*nrow(p.infection.sim)/4))] <- 1
#p.infection.sim <- as.matrix(p.infection.sim)

library(dplyr)
p.infection.sim$svl.hab.mean <- meanfunc.inf.svl.sex.habitat(p.infection.sim$svl,
                                                             p.infection.sim$sex,
                                                                 p.infection.sim$habitat)
p.infection.sim$svl.hab.lower <- lowerfunc.inf.svl.sex.habitat(p.infection.sim$svl,
                                                               p.infection.sim$sex,
                                                                   p.infection.sim$habitat)
p.infection.sim$svl.hab.upper <- upperfunc.inf.svl.sex.habitat(p.infection.sim$svl,
                                                               p.infection.sim$sex,
                                                                   p.infection.sim$habitat)

plot.inf.hab.svl <- ggplot(data=p.infection.sim, aes(x=svl,y=svl.hab.mean,group=sex))+
  geom_line(aes(color=habitat))+
  geom_line(aes(y=svl.hab.lower, color=habitat))+
  geom_line(aes(y=svl.hab.upper, color=habitat))+
  xlab("SVL")+ylab("p(Infection)")

#infection~ habitat+svl####
mod.infection.svl.habitat.tab <- mcmctab(mod.infection.svl.habitat)
meanfunc.inf.svl.habitat <- function(SVL, habitat){
  return(invlogit(mod.infection.svl.habitat.tab$Median[1]
                  + mod.infection.svl.habitat.tab$Median[2]*SVL
                  + mod.infection.svl.habitat.tab$Median[3]*habitat
                  + mod.infection.svl.habitat.tab$Median[4]*SVL*habitat))
}
lowerfunc.inf.svl.habitat <- function(SVL, habitat){
  return(invlogit(mod.infection.svl.habitat.tab$Lower[1]
                  + mod.infection.svl.habitat.tab$Lower[2]*SVL
                  + mod.infection.svl.habitat.tab$Lower[3]*habitat
                  + mod.infection.svl.habitat.tab$Lower[4]*SVL*habitat))
}
upperfunc.inf.svl.habitat <- function(SVL, habitat){
  return(invlogit(mod.infection.svl.habitat.tab$Upper[1]
                  + mod.infection.svl.habitat.tab$Upper[2]*SVL
                  + mod.infection.svl.habitat.tab$Upper[3]*habitat
                  + mod.infection.svl.habitat.tab$Upper[4]*SVL*habitat))
}

p.infection.sim <- data.frame(svl = rep(seq(from=min(diagnosed$SVL),
                                            to=max(diagnosed$SVL),
                                            by=0.5),2))
p.infection.sim$habitat <- 0
p.infection.sim$habitat[1:(nrow(p.infection.sim)/2)] <- 1
#p.infection.sim <- as.matrix(p.infection.sim)

library(dplyr)
p.infection.sim$svl.hab.mean <- meanfunc.inf.svl.habitat(p.infection.sim$svl,
                                                             p.infection.sim$habitat)
p.infection.sim$svl.hab.lower <- lowerfunc.inf.svl.habitat(p.infection.sim$svl,
                                                               p.infection.sim$habitat)
p.infection.sim$svl.hab.upper <- upperfunc.inf.svl.habitat(p.infection.sim$svl,
                                                               p.infection.sim$habitat)

plot.inf.hab.svl <- ggplot(data=p.infection.sim, aes(x=svl,y=svl.hab.mean))+
  geom_point(aes(color=habitat))+
  geom_point(aes(y=svl.hab.lower, color=habitat))+
  geom_point(aes(y=svl.hab.upper, color=habitat))+
  xlab("SVL")+ylab("p(Infection)")

meanfunc.inf.svl.habitat <- function(svl,habitat){
  return(exp(-2.7118+0.246*svl-1.533*habitat+0.351*svl*habitat)/
           (1+exp(-2.7118+0.246*svl-1.533*habitat+0.351*svl*habitat)))
}

lowerfunc.inf.svl.habitat <- function(svl,habitat){
  return(exp(-5.6288511-0.189*svl-4.949*habitat-0.179*svl*habitat)/
           (1+exp(-5.6288511-0.189*svl-4.949*habitat-0.179*svl*habitat)))
}

upperfunc.inf.svl.habitat <- function(svl,habitat){
  return(exp(0.003656582+0.672*svl+1.507*habitat+0.935*svl*habitat)/
           (1+exp(0.003656582+0.672*svl+1.507*habitat+0.935*svl*habitat)))
}

p.infection.predictions <- data.frame(svl = rep(seq(from=min(diagnosed$SVL),
                                                    to=max(diagnosed$SVL),
                                                    by=0.05),2))
p.infection.predictions$habitat <- 0
p.infection.predictions$habitat[1:(nrow(p.infection.predictions)/2)] <- 1

library(dplyr)
p.infection.predictions$svl.hab.mean <- meanfunc.inf.svl.habitat(p.infection.predictions$svl,
                                                                 p.infection.predictions$habitat)
p.infection.predictions$svl.hab.lower <- lowerfunc.inf.svl.habitat(p.infection.predictions$svl,
                                                                   p.infection.predictions$habitat)
p.infection.predictions$svl.hab.upper <- upperfunc.inf.svl.habitat(p.infection.predictions$svl,
                                                                   p.infection.predictions$habitat)

plot.inf.hab.svl <- ggplot(data=p.infection.predictions, aes(x=svl,y=svl.hab.mean,group=habitat))+
  geom_line(aes(color=habitat))+
  geom_line(aes(y=svl.hab.lower, color=habitat),linetype=2)+
  geom_line(aes(y=svl.hab.upper, color=habitat),linetype=2)+
  xlab("SVL")+ylab("p(Infection)")

infection.svl_X_habitat
confint(infection.svl_X_habitat)

#infection~ svl####
meanfunc.inf.svl <- function(svl){
  return(exp(0.3830*svl)/(1+exp(0.3830*svl)))
}
lowerfunc.inf.svl <- function(svl){
  return(exp(0.0400*svl)/(1+exp(0.0400*svl)))
}
upperfunc.inf.svl <- function(svl){
  return(exp(0.7320*svl)/(1+exp(0.7320*svl)))
}

p.infection.predictions$svl.mean <- meanfunc.inf.svl(p.infection.predictions$svl)
p.infection.predictions$svl.lower <- lowerfunc.inf.svl(p.infection.predictions$svl)
p.infection.predictions$svl.upper <- upperfunc.inf.svl(p.infection.predictions$svl)

plot.inf.svl <- ggplot(data=p.infection.predictions, aes(x=svl,y=svl.mean))+
  geom_line()+
  geom_ribbon(aes(ymin = svl.lower, ymax = svl.upper), fill = "blue", alpha = 0.2)+
  geom_line(aes(y = svl.lower), color = "blue", alpha = 0.8, size = 0.5) + 
  geom_line(aes(y = svl.upper), color = "blue", alpha = 0.8, size = 0.5)+
  #geom_line(aes(y=svl.lower),linetype=2)+
  #geom_line(aes(y=svl.upper),linetype=2)+
  xlab("SVL")+ylab("p(Infection)")

plot_grid(plot.inf.hab,plot.inf.hab.svl,plot.inf.svl,labels="AUTO",nrow=1)

#infection~ avg model####
meanfunc.inf <- function(svl,habitat){
  return(exp(-7.4433+0.1412*svl+0.1887*habitat+0.0540*svl*habitat)/
           (1+exp(-7.4433+0.1412*svl+0.1887*habitat+0.0540*svl*habitat)))
}
lowerfunc.inf <- function(svl,habitat){
  return(exp(-8.8450-0.0663*svl-1.2960*habitat-0.1386*svl*habitat)/
           (1+exp(-8.8450-0.0663*svl-1.2960*habitat-0.1386*svl*habitat)))
}
upperfunc.inf <- function(svl,habitat){
  return(exp(-6.0828+0.3484*svl+1.5555*habitat+0.1281*svl*habitat)/
           (1+exp(-6.0828+0.3484*svl+1.5555*habitat+0.1281*svl*habitat)))
}

p.infection.predictions$avg.mean <- meanfunc.inf(p.infection.predictions$svl,
                                             p.infection.predictions$habitat)
p.infection.predictions$avg.lower <- lowerfunc.inf(p.infection.predictions$svl,
                                              p.infection.predictions$habitat)
p.infection.predictions$avg.upper <- upperfunc.inf(p.infection.predictions$svl,
                                               p.infection.predictions$habitat)


p.infection.predictions$habitat <- as.factor(p.infection.predictions$habitat)
plot.inf.avg <- ggplot(data=p.infection.predictions, aes(x=svl,y=avg.mean,group=as.factor(habitat)))+
  geom_line(aes(color=habitat))+
  geom_line(aes(y=avg.lower, color=habitat),linetype=2)+
  geom_line(aes(y=avg.upper, color=habitat),linetype=2)+
  xlab("SVL")+ylab("p(Infection)")


#hemoglobin####
hemoglobin.predictions <- data.frame(svl = rep(seq(from=min(diagnosed$SVL),
                                               to=max(diagnosed$SVL),
                                               by=0.05),4))
hemoglobin.predictions$sex <- 0
hemoglobin.predictions$sex[1:(nrow(hemoglobin.predictions)/2)] <- 1
hemoglobin.predictions$infection <- 0
hemoglobin.predictions$infection[1:(nrow(hemoglobin.predictions)/4)] <- 1
hemoglobin.predictions$infection[(nrow(hemoglobin.predictions)/2+1):(3*nrow(hemoglobin.predictions)/4)] <- 1

meanfunc.hemo <- function(svl,sex,infection){
  return(2.79784-0.11127*infection+0.13652*svl+0.51471*sex-0.10830*svl*sex)
}
twentyfive.hemo <- function(svl,sex,infection){
  return(2.60799-0.14055*infection+0.09890*svl+0.35682*sex-0.13937*svl*sex)
}
lowerfunc.hemo <- function(svl,sex,infection){
  return(2.24981-0.19612*infection+0.03078*svl+0.03819*sex-0.20121*svl*sex)
}
seventyfive.hemo <- function(svl,sex,infection){
  return(2.98968-0.08270*infection+0.17300*svl+0.67004*sex-0.07638*svl*sex)
}
upperfunc.hemo <- function(svl,sex,infection){
  return(3.34544-0.02627*infection+0.24254*svl+0.97841*sex-0.01397*svl*sex)
}

hemoglobin.predictions$avg.mean <- meanfunc.hemo(hemoglobin.predictions$svl,hemoglobin.predictions$sex,hemoglobin.predictions$infection)
hemoglobin.predictions$avg.lower <- lowerfunc.hemo(hemoglobin.predictions$svl,hemoglobin.predictions$sex,hemoglobin.predictions$infection)
hemoglobin.predictions$avg.upper <- upperfunc.hemo(hemoglobin.predictions$svl,hemoglobin.predictions$sex,hemoglobin.predictions$infection)
hemoglobin.predictions$avg.twentyfive <- twentyfive.hemo(hemoglobin.predictions$svl,hemoglobin.predictions$sex,hemoglobin.predictions$infection)
hemoglobin.predictions$avg.seventyfive <- seventyfive.hemo(hemoglobin.predictions$svl,hemoglobin.predictions$sex,hemoglobin.predictions$infection)

hemo.M <- data.frame()
hemo.F <- data.frame()
hemo.M <- hemoglobin.predictions[1:(nrow(hemoglobin.predictions)/2),]
hemo.F <- hemoglobin.predictions[(nrow(hemoglobin.predictions)/2+1):(nrow(hemoglobin.predictions)),]
hemo.F <- hemo.F %>% filter(svl <= 5.45,svl >= 4.2)

hemo.M$infection <- as.factor(hemo.M$infection)
hemo.M.plot <- ggplot(data=hemo.M, aes(x=svl,y=avg.mean,group=infection))+
  geom_line(aes(color=infection),size=1)+
  geom_line(aes(y=avg.lower, color=infection),linetype=2,size=1)+
  geom_line(aes(y=avg.upper, color=infection),linetype=2,size=1)+
  xlab("SVL")+ylab("Hemoglobin (%PCV)")+ggtitle("Males")+
  theme(legend.position = "none")

hemo.F$infection <- as.factor(hemo.F$infection)
hemo.F.plot <- ggplot(data=hemo.F, aes(x=svl,y=avg.mean,group=infection))+
  geom_line(aes(color=infection),size=1)+
  geom_line(aes(y=avg.lower, color=infection),linetype=2,size=1)+
  geom_line(aes(y=avg.upper, color=infection),linetype=2,size=1)+
  xlab("SVL")+ylab("Hemoglobin (%PCV)")+ggtitle("Females")

library(cowplot)
library(dplyr)
plot_grid(hemo.M.plot,hemo.F.plot,labels="AUTO",nrow=1)
max(diagnosed.females$SVL)
min(diagnosed.females$SVL)
max(diagnosed.males$SVL)
min(diagnosed.males$SVL)
