library(dplyr)

data <- as.data.frame(read.csv("datasheets/Thesis Datasheet_ThresholdsAsNA.csv"), header=T, stringsAsFactors=F)

#### Calculate/Aggregate Columns ####


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



#infection status
diagnosed <- data[1:135,] %>% filter(!is.na(infection))
infected <- diagnosed %>% filter(infection == 1)
non.infected <- diagnosed %>% filter(infection == 0)





#### Filter data ####
#infection status
diagnosed <- data[136:486,] %>% filter(!is.na(infection))
infected <- data[136:486,] %>% filter(infection == 1)
non.infected <- data[136:486,] %>% filter(infection == 0)
#sex
males <- data[136:486,] %>% filter(SEX == 'M')
females <- data[136:486,] %>% filter(SEX == 'F')
#habitat type
urban <- data[136:486,] %>% filter(habitat == 'Non-Forested')
forest <- data[136:486,] %>% filter(habitat == 'Forested')

#### INFECTION ####
prop.infected <- nrow(infected)/nrow(diagnosed) #0.36 infection rate
##   SEX    ##
diagnosed.females <- females %>% filter(!is.na(infection))
diagnosed.males <- males %>% filter(!is.na(infection))

infected.females <- females %>% filter(!is.na(infection), infection == 1)
infected.males <- males %>% filter(!is.na(infection), infection == 1)
non.infected.males <- males %>% filter(!is.na(infection), infection == 0)

prevalence.females <- nrow(infected.females)/nrow(diagnosed.females) #.25 female infection rate
prevalence.males <- nrow(infected.males)/nrow(diagnosed.males) #.36 male infection rate

##   HABITAT  ##
diagnosed.urban <- urban %>% filter(!is.na(infection))
diagnosed.forest <- forest %>% filter(!is.na(infection))

infected.urban <- urban %>% filter(!is.na(infection), infection == 1)
infected.forest <- forest %>% filter(!is.na(infection), infection == 1)

prevalence.urban <- nrow(infected.urban)/nrow(diagnosed.urban) #0.26 urban infection rate
prevalence.forest <- nrow(infected.forest)/nrow(diagnosed.forest) #0.49 forest infection rate

diagnosed.U.F <- diagnosed.urban %>%filter(SEX=="F")
infected.U.F <- infected.urban %>%filter(SEX=="F")
diagnosed.U.M <- diagnosed.urban %>%filter(SEX=="M")
infected.U.M <- infected.urban %>%filter(SEX=="M")
diagnosed.Fo.F <- diagnosed.forest %>%filter(SEX=="F")
infected.Fo.F <- infected.forest %>%filter(SEX=="F")
diagnosed.Fo.M <- diagnosed.forest %>%filter(SEX=="M")
infected.Fo.M <- infected.forest %>%filter(SEX=="M")
prev.U.F <- nrow(infected.U.F)/nrow(diagnosed.U.F)
prev.U.M <- nrow(infected.U.M)/nrow(diagnosed.U.M)
prev.Fo.F <- nrow(infected.Fo.F)/nrow(diagnosed.Fo.F)
prev.Fo.M <- nrow(infected.Fo.M)/nrow(diagnosed.Fo.M)



#### Change characters to binary####
for(i in 1:nrow(diagnosed)){
  if(diagnosed$SEX[i] == "F"){diagnosed$SEX[i] = 0}
  else if(diagnosed$SEX[i] == "M"){diagnosed$SEX[i] = 1}}

for(i in 1:nrow(diagnosed)){
  if(diagnosed$habitat[i] == "Non-Forested"){diagnosed$habitat[i] = 0}
  else if(diagnosed$habitat[i] == "Forested"){diagnosed$habitat[i] = 1}}
diagnosed$habitat <- as.numeric(diagnosed$habitat)
diagnosed$SEX <- as.numeric(diagnosed$SEX)


for(i in 1:nrow(diagnosed.males)){
  if(diagnosed.males$habitat[i] == "Non-Forested"){diagnosed.males$habitat[i] = 0}
  else if(diagnosed.males$habitat[i] == "Forested"){diagnosed.males$habitat[i] = 1}}

for(i in 1:nrow(diagnosed.females)){
  if(diagnosed.females$habitat[i] == "Non-Forested"){diagnosed.females$habitat[i] = 0}
  else if(diagnosed.females$habitat[i] == "Forested"){diagnosed.females$habitat[i] = 1}}





diagnosed$SVL <- diagnosed$SVL*10
diagnosed$SVL <- diagnosed$SVL-mean(diagnosed$SVL)
diagnosed$Hct <- diagnosed$Hct/100