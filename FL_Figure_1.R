source("FL_data_carpentry.R")

p.infection.sim <- data.frame(svl = rep(seq(from=min(diagnosed$SVL),
                                            to=max(diagnosed$SVL),
                                            by=0.5),2))
p.infection.sim$habitat <- 0
p.infection.sim$habitat[1:(nrow(p.infection.sim)/2)] <- 1

inf.mod.mcmc <- as.mcmc(mod.infection.svl..habitat)
inf.mod.mcmc.combined <- as.mcmc(rbind(inf.mod.mcmc[[1]],
                                       inf.mod.mcmc[[2]],
                                       inf.mod.mcmc[[3]]))
pred.mean.mean <- 1/(1+exp(-(mean(inf.mod.mcmc.combined[,"alpha"])+
                               mean(inf.mod.mcmc.combined[,"b1"])*p.infection.sim$svl+
                               mean(inf.mod.mcmc.combined[,"b2"])*p.infection.sim$habitat+
                               mean(inf.mod.mcmc.combined[,"b3"])*p.infection.sim$svl*p.infection.sim$habitat)))
                           
pred.mean.dist <- matrix(NA, nrow=nrow(inf.mod.mcmc.combined),
                         ncol=length(p.infection.sim$svl))
for (i in 1:nrow(pred.mean.dist)){
  pred.mean.dist[i,] <- 1/(1+exp(-(inf.mod.mcmc.combined[i,"alpha"]+
                                     inf.mod.mcmc.combined[i,"b1"]*p.infection.sim$svl+
                                     inf.mod.mcmc.combined[i,"b2"]*p.infection.sim$habitat+
                                     inf.mod.mcmc.combined[i,"b3"]*p.infection.sim$svl*p.infection.sim$habitat)))
}
credible_lower <- apply(pred.mean.dist, MARGIN = 2, quantile, prob = 0.025)
credible_upper <- apply(pred.mean.dist, MARGIN = 2, quantile, prob = 0.975)


inf.predictions <- data.frame(mean=0,upper=0,lower=0)
inf.predictions <- data.frame(pred.mean.mean,credible_upper,credible_lower,
                              p.infection.sim$svl,p.infection.sim$habitat)
colnames(inf.predictions) <- c('infection','upper','lower','svl','habitat')

for(i in 1:nrow(inf.predictions)){
  if(inf.predictions$habitat[i]==1){inf.predictions$habitat[i] <- 'Forested'}
  else{inf.predictions$habitat[i] <- 'Deforested'}
}

ggplot(data=diagnosed, aes(y=infection,x=SVL))+theme_bw()+
  labs(y="p(Infection)",x="SVL (mm)")+
  geom_jitter(aes(color=habitat,shape=SEX),height=0.03, size=1)+
  geom_ribbon(data=inf.predictions[72:142,],aes(ymax=upper,ymin=lower,x=svl,fill=habitat,color=habitat),alpha=0.3,linetype=1)+
  geom_line(data=inf.predictions[72:142,],aes(y=infection,x=svl,color=habitat),size=1,linetype=1)+ 
  geom_ribbon(data=inf.predictions[1:71,],aes(ymax=upper,ymin=lower,x=svl,fill=habitat,color=habitat),alpha=0.3,linetype=5)+
  geom_line(data=inf.predictions[1:71,],aes(y=infection,x=svl,color=habitat),size=1,linetype=5)+
  geom_vline(xintercept=55, linetype="dotted", size=.75)+
  scale_shape_manual(name='Sex',values=c('F'=2,'M'=19), labels=c('Female','Male'))+
  scale_color_manual(name='Habitat',values=c('Forested'="gold3",'Deforested'='cornflowerblue'))+
  scale_fill_manual(name='Habitat',values=c('Forested'="gold3",'Deforested'='cornflowerblue'))+theme(panel.grid.major = element_blank(), 
                                                                                                       panel.grid.minor = element_blank())


  