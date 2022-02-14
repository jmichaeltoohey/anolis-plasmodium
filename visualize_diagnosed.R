source('data_carpentry.R')

## SIZE ##
#SVL
ggplot(diagnosed, aes(as.factor(infection), SVL))+geom_boxplot()          #no pattern
#MASS
ggplot(diagnosed, aes(as.factor(infection), MASS))+geom_boxplot()         #no pattern
#body condition
ggplot(diagnosed, aes(as.factor(infection), bodycond))+geom_boxplot()     #a bit lower in infecteds
ggplot(diagnosed, aes(as.factor(infection), sex.bodycond))+geom_boxplot() #a bit lower in infecteds

## Temperature Difference ##
ggplot(diagnosed, aes(as.factor(infection), y=temp.diff))+geom_boxplot()       #no pattern
#outliers (lizards basking on extremely hot surfaces) seem to be infected

## Blood Chemistry ##
#Na
ggplot(diagnosed, aes(as.factor(infection), y=Na))+geom_boxplot()              #no pattern
#K
ggplot(diagnosed, aes(as.factor(infection), y=K))+geom_boxplot()               #no pattern
#Cl
ggplot(diagnosed, aes(as.factor(infection), y=Cl))+geom_boxplot()              #slightly higher in infected
#TCO2
ggplot(diagnosed, aes(as.factor(infection), y=TCO2))+geom_boxplot()            #no pattern
#Crea
ggplot(diagnosed, aes(as.factor(infection), y=Crea))+geom_boxplot()            #no pattern
#negligible sample size
#Glu
ggplot(diagnosed, aes(as.factor(infection), y=Glu))+geom_boxplot()             #slightly higher in noninfected
#iCa
ggplot(diagnosed, aes(as.factor(infection), y=iCa))+geom_boxplot()             #no pattern
#AnGap
ggplot(diagnosed, aes(as.factor(infection), y=AnGap))+geom_boxplot()           #no pattern
#Hct
ggplot(diagnosed, aes(as.factor(infection), y=Hct))+geom_boxplot()             #no pattern
#Hb
ggplot(diagnosed, aes(as.factor(infection), y=Hb))+geom_boxplot()              #no pattern
