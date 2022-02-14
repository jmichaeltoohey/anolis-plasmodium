

p.infection.sim <- data.frame(svl = rep(seq(from=min(diagnosed$SVL),
                                            to=max(diagnosed$SVL),
                                            by=0.5),2))
p.infection.sim$habitat <- 0
p.infection.sim$habitat[1:(nrow(p.infection.sim)/2)] <- 1
#svl.seq <- seq(min(diagnosed$SVL), max(diagnosed$SVL),by=0.5)
inf.mod.mcmc <- as.mcmc(mod.infection.svl.habitat)
plot(inf.mod.mcmc)
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
plot(infection~SVL*habitat,data=diagnosed,family=binomial)


inf.predictions <- data.frame(mean=0,upper=0,lower=0)
inf.predictions <- data.frame(pred.mean.mean,credible_upper,credible_lower,
                              p.infection.sim$svl,p.infection.sim$habitat)
colnames(inf.predictions) <- c('infection','upper','lower','svl','habitat')

###habitat type (Forested or Non-Forested)
data['habitat'] <- NA #add empty column
for(i in 1:nrow(data)){
  if(data$SITE[i]=='NATL' | data$SITE[i]=='McCarty' | data$SITE[i]=='Law-F' | data$SITE[i]=='Kew F' | 
     data$SITE[i]=='Env. F' | data$SITE[i]=='BTP F' | data$SITE[i]=='Alice'){data$habitat[i] <- 'Forested'}
  else{data$habitat[i] <- 'Non-Forested'}
}
#### Filter data #
diagnosed <- data[136:486,] %>% filter(!is.na(infection))
diagnosed$SVL <- diagnosed$SVL*10
diagnosed$SVL <- diagnosed$SVL-mean(diagnosed$SVL)

for(i in 1:nrow(inf.predictions)){
  if(inf.predictions$habitat[i]==1){inf.predictions$habitat[i] <- 'Forested'}
  else{inf.predictions$habitat[i] <- 'Non-Forested'}
}

ggplot(data=diagnosed, aes(y=infection,x=SVL))+theme_bw()+
  labs(y="p(Infection)",x="SVL (mm)")+
  geom_jitter(aes(color=habitat),height=0.03)+
  geom_ribbon(data=inf.predictions[72:142,],aes(ymax=upper,ymin=lower,x=svl,fill=habitat,color=habitat),alpha=0.3,linetype=1)+
  geom_line(data=inf.predictions[72:142,],aes(y=infection,x=svl,color=habitat),size=1,linetype=1)+ 
  geom_ribbon(data=inf.predictions[1:71,],aes(ymax=upper,ymin=lower,x=svl,fill=habitat,color=habitat),alpha=0.3,linetype=5)+
  geom_line(data=inf.predictions[1:71,],aes(y=infection,x=svl,color=habitat),size=1,linetype=5)+
  scale_color_manual(name='Habitat',values=c('Forested'="forestgreen",'Non-Forested'='cornflowerblue'))+
  scale_fill_manual(name='Habitat',values=c('Forested'="forestgreen",'Non-Forested'='cornflowerblue'))+guides(fill=FALSE)
                     #,values=c('Forested'="darkgreen",'Non-Forested'='darkgrey')


  