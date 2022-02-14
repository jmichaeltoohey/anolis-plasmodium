source('data_carpentry.R')

####Visualize####

## Habitat Type ##
#SVL
ggplot(data, aes(habitat, y=SVL))+geom_boxplot()             #slighlty longer in forested
#MASS
ggplot(data, aes(habitat, y=MASS))+geom_boxplot()            #slightly heavier in forested
#body condition
ggplot(data, aes(habitat, y=bodycond))+geom_boxplot()        #no pattern
ggplot(data, aes(habitat, y=sex.bodycond))+geom_boxplot()        #no pattern
#temp.diff
ggplot(data, aes(habitat, y=temp.diff))+geom_boxplot()       #no pattern
#extreme differences in urban where lizards bask on extremely hot surfaces
#Na
ggplot(data, aes(habitat, y=Na))+geom_boxplot()              #no pattern
#K
ggplot(data, aes(habitat, y=K))+geom_boxplot()               #higher in urban
#Cl
ggplot(data, aes(habitat, y=Cl))+geom_boxplot()              #no pattern
#TCO2
ggplot(data, aes(habitat, y=TCO2))+geom_boxplot()            #slightly higher in urban
#Crea
ggplot(data, aes(habitat, y=Crea))+geom_boxplot()            #higher in urban
#negligible sample size
#Glu
ggplot(data, aes(habitat, y=Glu))+geom_boxplot()             #no pattern
#iCa
ggplot(data, aes(habitat, y=iCa))+geom_boxplot()             #slightly higher in forest
#AnGap
ggplot(data, aes(habitat, y=AnGap))+geom_boxplot()           #no pattern
#Hct
ggplot(data, aes(habitat, y=Hct))+geom_boxplot()             #no pattern
#Hb
ggplot(data, aes(habitat, y=Hb))+geom_boxplot()              #no pattern

##  SVL  ##
#mass by sex
ggplot(data, aes(x=SVL, y=MASS, color=SEX))+geom_point()+geom_smooth(method='lm')
#two different positive patterns
#temp.diff
ggplot(data, aes(SVL, y=temp.diff))+geom_point()       #no pattern
#Na
ggplot(data, aes(SVL, y=Na))+geom_point()              #no pattern
#K
ggplot(data, aes(SVL, y=K))+geom_point()               #no pattern
#Cl
ggplot(data, aes(SVL, y=Cl))+geom_point()              #no pattern
#TCO2
ggplot(data, aes(SVL, y=TCO2))+geom_point()            #weak positive trend
#Crea
ggplot(data, aes(SVL, y=Crea))+geom_point()            #no pattern
#Glu
ggplot(data, aes(SVL, y=Glu))+geom_point()             #no pattern
#iCa
ggplot(data, aes(SVL, y=iCa))+geom_point()             #slight negative pattern
#AnGap
ggplot(data, aes(SVL, y=AnGap))+geom_point()           #slight positive pattern
#Hct
ggplot(data, aes(SVL, y=Hct))+geom_point()             #no pattern
#Hb
ggplot(data, aes(SVL, y=Hb))+geom_point()              #no pattern
#habitat
ggplot(data, aes(SVL, y=habitat))+geom_boxplot()


## SEX ##
#SVL
ggplot(data, aes(SEX, y=SVL))+geom_boxplot()             #greater for males
#MASS
ggplot(data, aes(SEX, y=MASS))+geom_boxplot()            #greater for males
#body condition
ggplot(data, aes(SEX, y=bodycond))+geom_boxplot()        #no pattern
#temp.diff
ggplot(data, aes(SEX, y=temp.diff))+geom_boxplot()       #no pattern
#Na
ggplot(data, aes(SEX, y=Na))+geom_boxplot()              #no pattern
#K
ggplot(data, aes(SEX, y=K))+geom_boxplot()               #no pattern
#Cl
ggplot(data, aes(SEX, y=Cl))+geom_boxplot()              #no pattern
#TCO2
ggplot(data, aes(SEX, y=TCO2))+geom_boxplot()            #greater for males
#Crea
ggplot(data, aes(SEX, y=Crea))+geom_boxplot()            #negligible sample size
#Glu
ggplot(data, aes(SEX, y=Glu))+geom_boxplot()             #slightly greater in males
#iCa
ggplot(data, aes(SEX, y=iCa))+geom_boxplot()             #greater for females
#AnGap
ggplot(data, aes(SEX, y=AnGap))+geom_boxplot()           #negligible sample size
#Hct
ggplot(data, aes(SEX, y=Hct))+geom_boxplot()             #no pattern
#Hb
ggplot(data, aes(SEX, y=Hb))+geom_boxplot()              #no pattern

##  SVL  ##
#mass by sex
ggplot(data, aes(x=SVL, y=MASS, color=SEX))+geom_point()+geom_smooth(method='lm')
#two different positive patterns
#temp.diff
ggplot(data, aes(SVL, y=temp.diff))+geom_point()       #no pattern
#Na
ggplot(data, aes(SVL, y=Na))+geom_point()              #no pattern
#K
ggplot(data, aes(SVL, y=K))+geom_point()               #no pattern
#Cl
ggplot(data, aes(SVL, y=Cl))+geom_point()              #no pattern
#TCO2
ggplot(data, aes(SVL, y=TCO2))+geom_point()            #no pattern
#Crea
ggplot(data, aes(SVL, y=Crea))+geom_point()            #no pattern
#Glu
ggplot(data, aes(SVL, y=Glu))+geom_point()             #no pattern
#iCa
ggplot(data, aes(SVL, y=iCa))+geom_point()             #slight negative pattern
#AnGap
ggplot(data, aes(SVL, y=AnGap))+geom_point()           #slight positive pattern
#Hct
ggplot(data, aes(SVL, y=Hct))+geom_point()             #no pattern
#Hb
ggplot(data, aes(SVL, y=Hb))+geom_point()              #no pattern
#habitat
ggplot(data, aes(SVL, y=habitat))+geom_boxplot()



## Body Condition ##
ggplot(data, aes(SEX, bodycond))+geom_boxplot(notch=TRUE)
#Na
ggplot(data, aes(sex.bodycond, y=Na))+geom_point()       #no pattern
#K
ggplot(data, aes(sex.bodycond, y=K))+geom_point()        #no pattern
#Cl
ggplot(data, aes(sex.bodycond, y=Cl))+geom_point()       #no pattern
#TCO2
ggplot(data, aes(sex.bodycond, y=TCO2))+geom_point()     #no pattern
#Crea
ggplot(data, aes(sex.bodycond, y=Crea))+geom_point()     #no pattern
#Glu
ggplot(data, aes(sex.bodycond, y=Glu))+geom_point()      #no pattern
#iCa
ggplot(data, aes(sex.bodycond, y=iCa))+geom_point()      #slight negative pattern
#AnGap
ggplot(data, aes(sex.bodycond, y=AnGap))+geom_point()    #slight positive pattern
#Hct
ggplot(data, aes(sex.bodycond, y=Hct))+geom_point()      #no pattern
#Hb
ggplot(data, aes(sex.bodycond, y=Hb))+geom_point()       #no pattern



## Temperatures ##
ggplot(data, aes(x=data$TEMP.SUB, y=data$TEMP.LIZ))+geom_point()+geom_smooth(method='lm')



















hist(data$TEMP.LIZ)
hist(data$TEMP.SUB)
mean(na.omit(data$TEMP.LIZ))
mean(na.omit(data$TEMP.SUB))
hist(data$SVL)
hist(data$MASS)
hist(na.omit(is.numeric(data$Na)))

library(ggplot2)
ggplot(data, aes(SVL, MASS, color=SEX)) + geom_point() + geom_smooth(method=lm)
ggplot()