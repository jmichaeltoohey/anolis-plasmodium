source('data_carpentry.R')

library(BayesFactor)
library(bain)
library(BFpack)


####p(infection)~####

#p(infection)~(SVL*Habitat)
infection.svl_X_habitat = glm(data=diagnosed, infection ~ SVL * habitat,family=binomial)
summary(infection.svl_X_habitat)

ct <- lmtest::coeftest(infection.svl_X_habitat)
BF.inf2.1 <- BF(infection.svl_X_habitat, Sigma = diag(ct[,2]^2), n = nrow(infection.svl_X_habitat))
summary(BF.inf2.1)
print(BF.inf2.1)
BF.inf2.2 <- BF(infection.svl_X_habitat,
                hypothesis = "habitat>0")
summary(BF.inf2.2)
BF.inf2.3 <- BF(infection.svl_X_habitat,
                hypothesis = "SVL>0")
summary(BF.inf2.3)
BF.inf2.4 <- BF(infection.svl_X_habitat,
                hypothesis = "SVL:habitat>0")
summary(BF.inf2.4)
BF.inf2.4 <- BF(infection.svl_X_habitat,
                hypothesis = "SVL:habitat+SVL>0")
summary(BF.inf2.4)



####body condition####
#body condition~1
bodycond.null = lm(diagnosed.males$sex.bodycond ~ 1)
summary(bodycond.null)

BF.bodycond4.1 <- BF(bodycond.null)
summary(BF.bodycond4.1)

#body condition ~ infection
bodycond.infection = lm(diagnosed.males$sex.bodycond ~ diagnosed.males$infection)
summary(bodycond.infection)

BF.bodycond1.1 <- BF(bodycond.infection,
                hypothesis = "infection<0")
summary(BF.bodycond1.1)

# body condition ~ Habitat
bodycond.habitat = lm(diagnosed.males$sex.bodycond ~ diagnosed.males$habitat)
summary(bodycond.habitat)

BF.bodycond3.1 <- BF(bodycond.habitat,
                hypothesis = "habitat<0")
summary(BF.bodycond3.1)

# body condition ~ INFECTION + Habitat
bodycond.infection.habitat = lm(diagnosed.males$sex.bodycond ~ diagnosed.males$infection + diagnosed.males$habitat)
summary(bodycond.infection.habitat)

BF.bodycond2.1 <- BF(bodycond.infection.habitat,
                hypothesis = "infection<0")
summary(BF.bodycond2.1)
BF.bodycond2.2 <- BF(bodycond.infection.habitat,
                     hypothesis = "habitat<0")
summary(BF.bodycond2.2)

#PR BODYCOND
bodycond.PR.infection = lm(PR.data.males$bodycond ~ PR.data.males$Infeccion)
summary(bodycond.PR.infection)

BF.bodycond2.1 <- BF(bodycond.PR.infection,
                     hypothesis = "Infeccion>0")
summary(BF.bodycond2.1)

####Na####
#Na~1
Na.null = glm(data=diagnosed, Na ~ 1, family=Gamma(link = "inverse"))
summary(Na.null)

BF.Na1.1 <- BF(Na.null,
                hypothesis = "SVL=SEX=habitat=SVL:SEX=SVL:habitat=0")
summary(BF.Na1.1)


#Na ~ infection
Na.infection = glm(data=diagnosed, Na ~ infection, family=Gamma(link = "inverse"))
summary(Na.infection)

BF.Na4.1 <- BF(Na.infection,
                hypothesis = "infection=0")
summary(BF.Na4.1)


# Na ~ SVL
Na.svl = glm(data=diagnosed, Na ~ SVL, family=Gamma(link = "inverse"))
summary(Na.svl)

BF.Na2.1 <- BF(Na.svl,
                hypothesis = "SVL=0")
summary(BF.Na2.1)


# Na ~ Sex
Na.sex = glm(data=diagnosed, Na ~ SEX, family=Gamma(link = "inverse"))
summary(Na.sex)

BF.Na3.1 <- BF(Na.sex,
                hypothesis = "SEX=0")
summary(BF.Na3.1)


# Na ~ Habitat
Na.habitat = glm(data=diagnosed, Na ~ habitat, family=Gamma(link = "inverse"))
summary(Na.habitat)

BF.Na5.1 <- BF(Na.habitat,
                hypothesis = "habitat=0")
summary(BF.Na5.1)



####K####
# K ~ SVL
K.svl = glm(data=diagnosed, K ~ SVL, family=Gamma(link = "inverse"))
summary(K.svl)

BF.K1.1 <- BF(K.svl,
                hypothesis = "SVL=0")
summary(BF.K1.1)


# K ~ Sex
K.sex = glm(data=diagnosed, K ~ SEX, family=Gamma(link = "inverse"))
summary(K.sex)

BF.K3.1 <- BF(K.sex,
                hypothesis = "SEX=0")
summary(BF.K3.1)


# K ~ INFECTION + SVL
K.infection.svl = glm(data=diagnosed, K ~ infection + SVL, family=Gamma(link = "inverse"))
summary(K.infection.svl)

BF.K2.1 <- BF(K.infection.svl,
                hypothesis = "SVL=infection=0")
summary(BF.K2.1)


# K ~ INFECTION + Sex
K.infection.sex = glm(data=diagnosed, K ~ infection + SEX, family=Gamma(link = "inverse"))
summary(K.infection.sex)

BF.K4.1 <- BF(K.infection.sex,
                hypothesis = "SEX=infection=0")
summary(BF.K4.1)

####Cl####
#Cl~1
Cl.null = glm(data=diagnosed, Cl ~ 1, family=Gamma(link = "inverse"))
summary(Cl.null)

BF.Cl1.1 <- BF(Cl.null,
                hypothesis = "SVL=SEX=habitat=SVL:SEX=SVL:habitat=0")
summary(BF.Cl1.1)

#Cl ~ infection
Cl.infection = glm(data=diagnosed, Cl ~ infection, family=Gamma(link = "inverse"))
summary(Cl.infection)

BF.Cl5.1 <- BF(Cl.infection,
                hypothesis = "infection=0")
summary(BF.Cl5.1)


# Cl ~ SVL
Cl.svl = glm(data=diagnosed, Cl ~ SVL, family=Gamma(link = "inverse"))
summary(Cl.svl)

BF.Cl2.1 <- BF(Cl.svl,
                hypothesis = "SVL=0")
summary(BF.Cl2.1)


# Cl ~ Sex
Cl.sex = glm(data=diagnosed, Cl ~ SEX, family=Gamma(link = "inverse"))
summary(Cl.sex)

BF.Cl4.1 <- BF(Cl.sex,
                hypothesis = "SEX=0")
summary(BF.Cl4.1)


# Cl ~ Habitat
Cl.habitat = glm(data=diagnosed, Cl ~ habitat, family=Gamma(link = "inverse"))
summary(Cl.habitat)

BF.Cl3.1 <- BF(Cl.habitat,
                hypothesis = "habitat=0")
summary(BF.Cl3.1)


#Cl~(Sex)+(Habitat)
Cl.sex_t_habitat = glm(data=diagnosed, Cl ~ SEX + habitat, family=Gamma(link = "inverse"))
summary(Cl.sex_t_habitat)

BF.Cl7.1 <- BF(Cl.sex_t_habitat,
                hypothesis = "SEX=habitat=0")
summary(BF.Cl7.1)


# Cl ~ INFECTION + SVL
Cl.infection.svl = glm(data=diagnosed, Cl ~ infection + SVL, family=Gamma(link = "inverse"))
summary(Cl.infection.svl)

BF.Cl6.1 <- BF(Cl.infection.svl,
                hypothesis = "SVL=infection=0")
summary(BF.Cl6.1)


# Cl ~ INFECTION + Habitat
Cl.infection.habitat = glm(data=diagnosed, Cl ~ infection + habitat, family=Gamma(link = "inverse"))
summary(Cl.infection.habitat)

BF.Cl8.1 <- BF(Cl.infection.habitat,
                hypothesis = "infection=habitat=0")
summary(BF.Cl8.1)



####iCa####

# iCa ~ SVL
iCa.svl = glm(data=diagnosed, iCa ~ SVL, family=Gamma(link = "inverse"))
summary(iCa.svl)

BF.iCa2.1 <- BF(iCa.svl,
                hypothesis = "SVL=0")
summary(BF.iCa2.1)


#iCa~(SVL*Sex)
iCa.svl_sex = glm(data=diagnosed, iCa ~ SVL * SEX, family=Gamma(link = "inverse"))
summary(iCa.svl_sex)

BF.iCa1.1 <- BF(iCa.svl_sex,
                hypothesis = "SVL=SEX=SVL:SEXt=0")
summary(BF.iCa1.1)


#iCa~(SVL*Habitat)
iCa.svl_X_habitat = glm(data=diagnosed, iCa ~ SVL * habitat, family=Gamma(link = "inverse"))
summary(iCa.svl_X_habitat)

BF.iCa4.1 <- BF(iCa.svl_X_habitat,
                hypothesis = "SVL=habitat=SVL:habitat=0")
summary(BF.iCa4.1)


#iCa~(SVL*Sex)+(SVL*Habitat)
iCa.svl_X_sex__t__svl_X_habitat = 
  glm(data=diagnosed, iCa ~ (SVL * SEX)+(SVL * habitat)
      , family=Gamma(link = "inverse"))
summary(iCa.svl_X_sex__t__svl_X_habitat)

BF.iCa5.1 <- BF(iCa.svl_X_sex__t__svl_X_habitat,
                hypothesis = "SVL=SEX=habitat=SVL:SEX=SVL:habitat=0")
summary(BF.iCa5.1)


# iCa ~ INFECTION + SVL
iCa.infection.svl = glm(data=diagnosed, iCa ~ infection + SVL, family=Gamma(link = "inverse"))
summary(iCa.infection.svl)

BF.iCa3.1 <- BF(iCa.infection.svl,
                hypothesis = "SVL=infection=0")
summary(BF.iCa3.1)


#iCa~INFECTION + (SVL*Sex)
iCa.infection.svl_sex = glm(data=diagnosed, iCa ~ infection + SVL * SEX, family=Gamma(link = "inverse"))
summary(iCa.infection.svl_sex)

BF.iCa6.1 <- BF(iCa.infection.svl_sex,
                hypothesis = "infection=SVL=SEX=SVL:SEX=0")
summary(BF.iCa6.1)


#iCa~INFECTION + (SVL*Habitat)
iCa.infection.svl_X_habitat = glm(data=diagnosed, iCa ~ infection + SVL * habitat, family=Gamma(link = "inverse"))
summary(iCa.infection.svl_X_habitat)

BF.iCa7.1 <- BF(iCa.infection.svl_X_habitat,
                hypothesis = "infection=SVL=habitat=SVL:habitat=0")
summary(BF.iCa7.1)



####Hct####

#Hct~1
Hct.null = glm(data=diagnosed, Hct ~ 1, family=Gamma(link = "inverse"))
summary(Hct.null)

BF.Hct2.1 <- BF(Hct.null,
                hypothesis = "SVL=SEX=habitat=SVL:SEX=SVL:habitat=0")
summary(BF.Hct2.1)


#Hct ~ infection
Hct.infection = glm(data=Hct.omit, Hct ~ infection, family = quasi(link = "logit", 
                                                                   variance = "mu(1-mu)"))
summary(Hct.infection)

BF.Hct3.1 <- BF(Hct.infection, #prior.hyp = mu ~ dunif(min=.15,max=.75),
                hypothesis = "infection<0")
summary(BF.Hct3.1)


# Hct ~ Sex
Hct.sex = glm(data=diagnosed, Hct ~ SEX, family = quasi(link = "logit", 
                                                        variance = "mu(1-mu)"))
summary(Hct.sex)

BF.Hct4.1 <- BF(Hct.sex,
                hypothesis = "SEX=0")
summary(BF.Hct4.1)


# Hct ~ INFECTION + Sex
Hct.infection.sex = glm(data=Hct.omit, Hct ~ infection + SEX, family = quasi(link = "logit", 
                                                                              variance = "mu(1-mu)"))
summary(Hct.infection.sex)

BF.Hct5.1 <- BF(Hct.infection.sex,
                hypothesis = "infection<0")
summary(BF.Hct5.1)


#Hct~INFECTION + (SVL*Sex)
Hct.infection.svl_sex = glm(data=Hct.omit, Hct ~ infection + SVL * SEX, family = quasi(link = "logit", 
                                                                                        variance = "mu(1-mu)"))
summary(Hct.infection.svl_sex)

BF.Hct1.1 <- BF(Hct.infection.svl_sex,
                hypothesis = "infection<0")
summary(BF.Hct1.1)
BF.Hct1.2 <- BF(Hct.infection.svl_sex,
                hypothesis = "SEX<0")
summary(BF.Hct1.2)
BF.Hct1.3 <- BF(Hct.infection.svl_sex,
                hypothesis = "SVL>0")
summary(BF.Hct1.3)
BF.Hct1.4 <- BF(Hct.infection.svl_sex,
                hypothesis = "SVL+SVL:SEX=0")
summary(BF.Hct1.4)

#Hct~INFECTION + (SVL*Sex)
Hct.svl_sex = glm(data=Hct.omit, Hct ~ SVL * SEX, family = quasi(link = "logit", 
                                                                                       variance = "mu(1-mu)"))
summary(Hct.svl_sex)

BF.Hct5.1 <- BF(Hct.svl_sex,
                hypothesis = "SEX<0")
summary(BF.Hct5.1)
BF.Hct5.2 <- BF(Hct.svl_sex,
                hypothesis = "SVL>0")
summary(BF.Hct5.2)
BF.Hct5.3 <- BF(Hct.svl_sex,
                hypothesis = "infection<0")
summary(BF.Hct5.3)



####Glu####

# Glu ~ SVL
Glu.svl = glm(data=diagnosed, Glu ~ SVL, family=Gamma(link = "inverse"))
summary(Glu.svl)

BF.Glu1.1 <- BF(Glu.svl,
                hypothesis = "SVL=0")
summary(BF.Glu1.1)


# Glu ~ Sex
Glu.sex = glm(data=diagnosed, Glu ~ SEX, family=Gamma(link = "inverse"))
summary(Glu.sex)

BF.Glu1.1 <- BF(Glu.sex,
                hypothesis = "SEX=0")
summary(BF.Glu1.1)


#Glu~(SVL*Sex)
Glu.svl_sex = glm(data=diagnosed, Glu ~ SVL * SEX, family=Gamma(link = "inverse"))
summary(Glu.svl_sex)

BF.Glu1.1 <- BF(Glu.svl_sex,
                hypothesis = "SVL=SEX=SVL:SEX=0")
summary(BF.Glu1.1)


#Glu~(Sex)+(Habitat)
Glu.sex_t_habitat = glm(data=diagnosed, Glu ~ SEX + habitat, family=Gamma(link = "inverse"))
summary(Glu.sex_t_habitat)

BF.Glu1.1 <- BF(Glu.sex_t_habitat,
                hypothesis = "SEX=habitat=0")
summary(BF.Glu1.1)


# Glu ~ INFECTION + Sex
Glu.infection.sex = glm(data=diagnosed, Glu ~ infection + SEX, family=Gamma(link = "inverse"))
summary(Glu.infection.sex)

BF.Glu1.1 <- BF(Glu.infection.sex,
                hypothesis = "infection=SEX=0")
summary(BF.Glu1.1)



#PC1 ~ infection
PC1.infection = lm(PC1.omit$PC1 ~ PC1.omit$infection)
summary(PC1.infection)

BF.PC11.1 <- BF(PC1.infection,
                     hypothesis = "infection>0")
summary(BF.PC11.1)

#PC1 ~ svl
PC1.SVL = lm(PC1.omit$PC1 ~ PC1.omit$SVL)
summary(PC1.SVL)

BF.PC12.1 <- BF(PC1.SVL,
                hypothesis = "SVL>0")
summary(BF.PC12.1)



#PC2 ~ infection
PC2.infection = lm(PC2.omit$PC2 ~ (PC2.omit$SVL*PC2.omit$SEX)+(PC2.omit$SVL*PC2.omit$habitat))
summary(PC2.infection)

BF.PC21.1 <- BF(PC2.infection,
                hypothesis = "SVL<0")
summary(BF.PC21.1)
BF.PC21.2 <- BF(PC2.infection,
                hypothesis = "SEX>0")
summary(BF.PC21.2)
BF.PC21.3 <- BF(PC2.infection,
                hypothesis = "habitat<0")
summary(BF.PC21.3)
BF.PC21.4 <- BF(PC2.infection,
                hypothesis = "SVL+SVL:SEX>0")
summary(BF.PC21.4)
BF.PC21.5 <- BF(PC2.infection,
                hypothesis = "SVL+SVL:habitat>0")
summary(BF.PC21.5)

#PC2 ~ svl
PC2.SVL = lm(PC2.omit$PC2 ~ PC2.omit$SVL)
summary(PC2.SVL)

BF.PC22.1 <- BF(PC2.SVL,
                hypothesis = "SVL>0")
summary(BF.PC22.1)
#PC2 ~ infection
PC2.infection = lm(PC2.omit$PC2 ~ PC2.omit$infection)
summary(PC2.infection)

BF.PC21.1 <- BF(PC2.infection,
                hypothesis = "infection>0")
summary(BF.PC21.1)

#PC2 ~ svl
PC2.SVL = lm(PC2.omit$PC2 ~ PC2.omit$SVL)
summary(PC2.SVL)

BF.PC22.1 <- BF(PC2.SVL,
                hypothesis = "SVL>0")
summary(BF.PC22.1)
