library(dplyr)

data <- as.data.frame(read.csv("datasheets/Thesis Datasheet_ThresholdsAsNA.csv"), header=T, stringsAsFactors=F)

#### Calculate/Aggregate Columns ####

###Create column that groups an individual as an inhabitant of either a 
###Forested or Deforested habitat
data['habitat'] <- NA #add empty column
for(i in 1:nrow(data)){
  if(data$SITE[i]=='NATL' | data$SITE[i]=='McCarty' | #UNUSED PRELIMINARY SITES
     data$SITE[i]=='Law-F' | data$SITE[i]=='Alice' | #UNUSED PRELIMINARY SITES
     data$SITE[i]=='Kew F' | data$SITE[i]=='Env. F' | data$SITE[i]=='BTP F')
    {data$habitat[i] <- 'Forested'}
  else{data$habitat[i] <- 'Deforested'}
}

###Create columns that denote which site an individual was captured in
###in order to control for the paired site selection design
###(2 plots per site; 1 Forested, 1 Deforested)
###One column holds name of site
###One column holds numerical to be used as factor
data['sitepair'] <- NA #add empty column
for(i in 1:nrow(data)){
  #------------------------UNUSED PRELIMINARY DATA-----------------------------#
  #if(data$SITE[i]=='NATL' | data$SITE[i]=='Entom'){
  #  data$sitepair[i] <- 'NATL' 
  #  data$sitepairnum[i] <- 4}
  #else if(data$SITE[i]=='McCarty' | data$SITE[i]=='Museum'){
  #  data$sitepair[i] <- 'WEC'
  #  data$sitepairnum[i] <- 5}
  #else if(data$SITE[i]=='Law-F' | data$SITE[i]=='Law-U'){
  #  data$sitepair[i] <- 'Law'
  #  data$sitepairnum[i] <- 6}
  #else if(data$SITE[i]=='Alice' | data$SITE[i]=='Ficke'){
  #  data$sitepair[i] <- 'Lake Alice'
  #  data$sitepairnum[i] <- 7}
  #----------------------------------------------------------------------------#
  if(data$SITE[i]=='Kew F' | data$SITE[i]=='Kew U'){
    data$sitepair[i] <- 'Kewannee'
    data$sitepairnum[i] <- 3}
  else if(data$SITE[i]=='Env. F' | data$SITE[i]=='Env. U'){
    data$sitepair[i] <- 'Env Center'
    data$sitepairnum[i] <- 1}
  else if(data$SITE[i]=='BTP F' | data$SITE[i]=='BTP U'){
    data$sitepair[i] <- 'Big Tree Park'
    data$sitepairnum[i] <- 2}
}

#### Filter data ####

###----------------------UNUSED PRELIMINARY DATA-------------------------------#
#infection status
#diagnosed.GNV <- data[1:135,] %>% filter(!is.na(infection))
#infected.GNV <- diagnosed.GNV %>% filter(infection == 1)
#non.infected.GNV <- diagnosed.GNV %>% filter(infection == 0)
#------------------------------------------------------------------------------#

###INFECTION STATUS
diagnosed <- data[136:486,] %>% filter(!is.na(infection))
infected <- diagnosed %>% filter(infection == 1)
non.infected <- diagnosed %>% filter(infection == 0)
#Calculate overall prevalence
prop.infected <- nrow(infected)/nrow(diagnosed) #0.33 infection rate

###SEX
males <- diagnosed %>% filter(SEX == 'M')
females <- diagnosed %>% filter(SEX == 'F')

###HABITAT
urban <- diagnosed %>% filter(habitat == 'Deforested')
forest <- diagnosed %>% filter(habitat == 'Forested')

##INFECTION BY SEX
diagnosed.females <- females %>% filter(!is.na(infection))
diagnosed.males <- males %>% filter(!is.na(infection))

infected.females <- females %>% filter(!is.na(infection), infection == 1)
infected.males <- males %>% filter(!is.na(infection), infection == 1)

prevalence.females <- nrow(infected.females)/nrow(diagnosed.females) #.25 female infection rate
prevalence.males <- nrow(infected.males)/nrow(diagnosed.males) #.36 male infection rate

##INFECTION BY HABITAT
diagnosed.urban <- urban %>% filter(!is.na(infection))
diagnosed.forest <- forest %>% filter(!is.na(infection))

infected.urban <- urban %>% filter(!is.na(infection), infection == 1)
infected.forest <- forest %>% filter(!is.na(infection), infection == 1)

prevalence.urban <- nrow(infected.urban)/nrow(diagnosed.urban) #0.26 DEFORESTED infection rate
prevalence.forest <- nrow(infected.forest)/nrow(diagnosed.forest) #0.43 FORESTED infection rate

##INFECTION BY SEX AND HABITAT
diagnosed.U.F <- diagnosed.urban %>%filter(SEX=="F")
infected.U.F <- infected.urban %>%filter(SEX=="F")
diagnosed.U.M <- diagnosed.urban %>%filter(SEX=="M")
infected.U.M <- infected.urban %>%filter(SEX=="M")
diagnosed.Fo.F <- diagnosed.forest %>%filter(SEX=="F")
infected.Fo.F <- infected.forest %>%filter(SEX=="F")
diagnosed.Fo.M <- diagnosed.forest %>%filter(SEX=="M")
infected.Fo.M <- infected.forest %>%filter(SEX=="M")
prev.U.F <- nrow(infected.U.F)/nrow(diagnosed.U.F) #0.25
prev.U.M <- nrow(infected.U.M)/nrow(diagnosed.U.M) #0.27
prev.Fo.F <- nrow(infected.Fo.F)/nrow(diagnosed.Fo.F) #0.25
prev.Fo.M <- nrow(infected.Fo.M)/nrow(diagnosed.Fo.M) #0.48



####CHANGE CHARACTERS TO BINARY####
for(i in 1:nrow(diagnosed)){
  if(diagnosed$SEX[i] == "F"){diagnosed$SEX[i] = 0}
  else if(diagnosed$SEX[i] == "M"){diagnosed$SEX[i] = 1}}

for(i in 1:nrow(diagnosed)){
  if(diagnosed$habitat[i] == "Deforested"){diagnosed$habitat[i] = 0}
  else if(diagnosed$habitat[i] == "Forested"){diagnosed$habitat[i] = 1}}
diagnosed$habitat <- as.numeric(diagnosed$habitat)
diagnosed$SEX <- as.numeric(diagnosed$SEX)


for(i in 1:nrow(diagnosed.males)){
  if(diagnosed.males$habitat[i] == "Deforested"){diagnosed.males$habitat[i] = 0}
  else if(diagnosed.males$habitat[i] == "Forested"){diagnosed.males$habitat[i] = 1}}

for(i in 1:nrow(diagnosed.females)){
  if(diagnosed.females$habitat[i] == "Deforested"){diagnosed.females$habitat[i] = 0}
  else if(diagnosed.females$habitat[i] == "Forested"){diagnosed.females$habitat[i] = 1}}

####CONVERT UNITS####
diagnosed$SVL <- diagnosed$SVL*10  #SVL from cm to mm
diagnosed$Hct <- diagnosed$Hct/100 #Hct from percentage to decimal