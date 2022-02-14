library(knitr)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(mgcv)
library(tidymv)
devtools::install_github("cardiomoon/ggiraphExtra")
require(ggiraph)
require(ggiraphExtra)
require(plyr)


#p(infection)~(SVL*Habitat)
diagnosed$habitat <- as.factor(diagnosed$habitat)
infection.svl_X_habitat = glm(infection ~ SVL * habitat,family=binomial("logit"), data=diagnosed)
predict.inf.svl.hab <- ggPredict(infection.svl_X_habitat,se=TRUE,interactive=TRUE)

infection.habitat = glm(infection ~ habitat,family=binomial, data=diagnosed)
predict.inf.hab <- ggPredict(infection.habitat,se=TRUE,interactive=TRUE,ci.lvl = 0.95)

infection.svl_X_sex__t__svl_X_habitat = 
  glm(infection ~ (SVL * habitat)+(SVL * SEX),
      family=binomial, data=diagnosed)
predict.inf.svl.hab <- ggPredict(infection.svl_X_sex__t__svl_X_habitat,se=TRUE,interactive=TRUE)

#M Body Condition
M.bodycond.infection.svl.habitat = lm(sex.bodycond ~ SVL * habitat + infection, data=diagnosed.males)
ggPredict(M.bodycond.infection.svl.habitat,se=TRUE,interactive=TRUE)

#M Body Condition
K.infection.svl = lm(K ~ SVL * SEX + infection, data=diagnosed)
ggPredict(K.infection.svl,se=TRUE,interactive=TRUE)

diagnosed$habitat <- as.factor(diagnosed$habitat)
infection.svl_X_habitat = glm(infection ~ SVL * habitat,family=binomial, data=diagnosed)
predict.inf.svl.hab <- ggPredict(infection.svl_X_habitat,se=TRUE,interactive=TRUE)


diagnosed$habitat <- as.numeric(diagnosed$habitat)
diagnosed$SEX <- as.numeric(diagnosed$SEX)
infection <- data.frame()
for(i in 1:nrow(diagnosed)){
  infection[i] <- 1/(1+exp(-(-7.44334 + 0.14124*diagnosed$SVL[i] + 0.1887*diagnosed$habitat[i] + 
    0.054*diagnosed$SVL[i]*diagnosed$habitat[i])))
}
infection <- -7.44334 + 0.14124*diagnosed$SVL
infection <- infection + 0.1887*diagnosed$habitat
infection <- infection + 0.054*diagnosed$SVL*diagnosed$habitat
 0.054*diagnosed$SVL*diagnosed$habitat

ggplot(data=diagnosed, aes(x=SVL,y=infection,col=habitat)) + 
  geom_point(shape = 1) +
  geom_abline(intercept = -7.44, slope = .4, col = "red") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)




1/(1+exp(-(-7.44334)))
1/(1+exp(-(.141)))
1/(1+exp(-(.1887)))
1/(1+exp(-(.054)))



