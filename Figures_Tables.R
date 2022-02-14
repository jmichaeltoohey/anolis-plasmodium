library(knitr)
library(ggplot2)
infection.compare = loo_compare(infection.list)
#rownames(infection.compare) <- c("p(infection) ~ 1","p(infection) ~ SVL", "p(infection) ~ Sex",
#                     "p(infection) ~ Habitat", "p(infection) ~ (SVL*Sex)", "p(infection) ~ (SVL*Habitat)",
#                     "p(infection) ~ (Sex)+(Habitat)", "p(infection) ~ (SVL*Sex)+(SVL*Habitat)")
kable(infection.compare, caption="p(Infection)")

p <- predict(mod.infection.sex, type = "response")
predict(mod.infection.sex, newdata=pred.inf.sex)
ggplot(diagnosed, aes(x=SEX, y=infection))+geom_boxplot()

table.sex.infection <- table(diagnosed$infection, diagnosed$SEX) 
mosaicplot(table.sex.infection, color=TRUE, xlab="Infection Status", ylab="Sex") #make a mosaic plot

bodycond.compare = loo_compare(bodycond.list)
kable(bodycond.compare, caption="Body Condition")

par(mfrow=c(1,2))
plot.bodycond.sex <- ggplot(diagnosed, aes(x=bodycond,y=infection,color=SEX))+
  geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))+
  labs(title="Body Condition by Infection Status and Sex",
       x ="Body Condition", y = "Infection Status")+
  scale_color_manual(labels = c("F", "M"), values = c("violet", "dodgerblue4"))

plot.bodycond.habitat <- ggplot(diagnosed, aes(x=bodycond,y=infection,color=habitat))+
  geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))+
  labs(title="Body Condition by Infection Status and Habitat Type",
       x ="Body Condition", y = "Infection Status")+
  scale_color_manual(labels = c("Non-Forested", "Forested"), values = c("lightsteelblue3", "darkgreen"))
plot.bodycond.sex + coord_flip()
plot.bodycond.habitat + coord_flip()

temp.diff.compare = loo_compare(temp.diff.list)
kable(temp.diff.compare, caption="Temperature Difference")

Hct.compare = loo_compare(Hct.list)
kable(Hct.compare, caption="Hematocrit Count")

ggplot(diagnosed, aes(x=SVL, y=Na, col=habitat))+geom_point()+
  geom_smooth(method = "glm", method.args = list(family = "Gamma"))+
  scale_color_manual(labels = c("Non-Forested", "Forested"), 
                     values = c("lightsteelblue3", "darkgreen"))

ggplot(diagnosed, aes(x=SVL, y=Na, col=SEX))+geom_point()+
  geom_smooth(method = "glm", method.args = list(family = "Gamma"))+
  scale_color_manual(labels = c("F", "M"), 
                     values = c("violet", "dodgerblue4"))

ggplot(diagnosed, aes(x=Hct,y=infection,color=SEX))+geom_point()+
  geom_smooth(method='glm', method.args = list(family = 'binomial'))+
  labs(title="Hematocrit Count by Infection Status and Sex", x ="Hct", y = "Infection Status")+
  scale_color_manual(labels = c("F", "M"), values = c("violet", "dodgerblue4"))


library(knitr)
library(ggplot2)
infection.compare = loo_compare(infection.list)
#rownames(infection.compare) <- c("p(infection) ~ 1","p(infection) ~ SVL", "p(infection) ~ Sex",
#                     "p(infection) ~ Habitat", "p(infection) ~ (SVL*Sex)", "p(infection) ~ (SVL*Habitat)",
#                     "p(infection) ~ (Sex)+(Habitat)", "p(infection) ~ (SVL*Sex)+(SVL*Habitat)")
kable(infection.compare, caption="p(Infection)")

table.sex.infection <- table(diagnosed$infection, diagnosed$SEX) 
mosaicplot(table.sex.infection, color=TRUE, xlab="Infection Status", ylab="Sex") #make a mosaic plot

bodycond.compare = loo_compare(bodycond.list)
kable(bodycond.compare, caption="Body Condition")

par(mfrow=c(1,2))
plot.bodycond.sex <- ggplot(diagnosed, aes(x=bodycond,y=infection,color=SEX))+
  geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))+
  labs(title="Body Condition by Infection Status and Sex",
       x ="Body Condition", y = "Infection Status")+
  scale_color_manual(labels = c("F", "M"), values = c("violet", "dodgerblue4"))
plot.bodycond.habitat <- ggplot(diagnosed, aes(x=bodycond,y=infection,color=habitat))+
  geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))+
  labs(title="Body Condition by Infection Status and Habitat Type",
      x ="Body Condition", y = "Infection Status")+
  scale_color_manual(labels = c("Non-Forested", "Forested"), values = c("lightsteelblue3", "darkgreen"))
plot.bodycond.sex + coord_flip()
plot.bodycond.habitat + coord_flip()

temp.diff.compare = loo_compare(temp.diff.list)
kable(temp.diff.compare, caption="Temperature Difference")

Hct.compare = loo_compare(Hct.list)
kable(Hct.compare, caption="Hematocrit Count")

ggplot(diagnosed, aes(x=SVL, y=Na, col=habitat))+
  geom_point()+geom_smooth(method = "glm", method.args = list(family = "Gamma"))+
  scale_color_manual(labels = c("Non-Forested", "Forested"), values = c("lightsteelblue3", "darkgreen"))

ggplot(diagnosed, aes(x=SVL, y=Na, col=SEX))+geom_point()+
  geom_smooth(method = "glm", method.args = list(family = "Gamma"))+
  scale_color_manual(labels = c("F", "M"), values = c("violet", "dodgerblue4"))

ggplot(diagnosed, aes(x=Hct,y=infection,color=SEX))+
  geom_point()+geom_smooth(method='glm', method.args = list(family = 'binomial'))+
  labs(title="Hematocrit Count by Infection Status and Sex",
       x ="Hct", y = "Infection Status")+
  scale_color_manual(labels = c("F", "M"), values = c("violet", "dodgerblue4"))