library("rstanarm")
library(rstanarm)
library(dplyr)
library(ggeffects)
library(ggiraphExtra)
library(ggiraph)
library(ggplot2)

#get data as factors####
data <- as.data.frame(read.csv("datasheets/Thesis Datasheet_ThresholdsAsNA.csv"), header=T, stringsAsFactors=F)

#### Calculate/Aggregate Columns #

data$SVL <- data$SVL*10

###habitat type (Forested or Non-Forested)
data['habitat'] <- NA #add empty column
for(i in 1:nrow(data)){
  if(data$SITE[i]=='NATL' | data$SITE[i]=='McCarty' | data$SITE[i]=='Law-F' | data$SITE[i]=='Kew F' | 
     data$SITE[i]=='Env. F' | data$SITE[i]=='BTP F' | data$SITE[i]=='Alice'){data$habitat[i] <- 'Forested'}
  else{data$habitat[i] <- 'Non-Forested'}
}

###temperatures
data$temp.diff <- data$TEMP.LIZ-data$TEMP.SUB
#ggplot(data, aes(x=temp.diff))+geom_histogram() #normal distribution

temp.ratio <- data$TEMP.LIZ/data$TEMP.SUB
#hist(temp.ratio)

###body condition column
#mass / svl
data <- mutate(data, mass.over.svl=MASS/SVL)
#hist(data$mass.over.svl) #BIMODAL... by sex?
#ggplot(data, aes(x=mass.over.svl, color=SEX))+geom_histogram() #

#residuals
mass.svl <- lm(MASS~SVL, data=data, na.action=na.exclude)
bodycond <- resid(mass.svl)
data$bodycond <- bodycond
#sex specific residuals
data$sex.bodycond <- NA
female.mass.svl <- lm(MASS~SVL, data=filter(data, SEX == 'F'), na.action=na.exclude)
female.bodycond <- resid(female.mass.svl)
male.mass.svl <- lm(MASS~SVL, data=filter(data, SEX == 'M'), na.action=na.exclude)
male.bodycond <- resid(male.mass.svl)
f = 1
m = 1
for(i in 1:nrow(data)){
  if(data$SEX[i] == 'F'){
    data$sex.bodycond[i] = female.bodycond[f]
    f = f + 1
  }
  else if(data$SEX[i] == 'M'){
    data$sex.bodycond[i] = male.bodycond[m]
    m = m + 1
  }
}


#### Filter data #
diagnosed <- data[136:486,] %>% filter(!is.na(infection))


#infection status ####
inf.svl.hab.rstanarm <- stan_glm(infection ~ SVL*habitat, 
                         data = diagnosed, family = binomial(link = "logit"),
                         prior = normal(0, 3),
                         prior_intercept = normal(0, 3),
                         chains = 3, 
                         iter = 10000,
                         seed = 123)

inf.Predict <- ggPredict(inf.svl.hab.rstanarm, se=T)+
  #theme(legend.position="None")+
  scale_color_manual(name='Habitat',values=c('Forested'="darkgreen",'Non-Forested'='darkgrey'))+
  guides(fill=FALSE)+
  labs(y="p(Infection)",x="SVL (mm)")
inf.Predict[["data"]][["ymax"]] <- inf.Predict[["data"]][["ymax"]]+
  0.96*inf.Predict[["data"]][["se.fit"]]
inf.Predict[["data"]][["ymin"]] <- inf.Predict[["data"]][["ymin"]]-
  0.96*inf.Predict[["data"]][["se.fit"]]
inf.Predict

#+theme(legend.postition="None")
#+scale_color_manual_interactive(name="HABITAT")
scale_fill
summary(inf.svl.hab.rstanarm)
print(inf.svl.hab.rstanarm)
