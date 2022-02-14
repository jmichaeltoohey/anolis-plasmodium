predicted.inf <- data.frame(intercept=rnorm(300,-7.44,.74), svl=rnorm(300,.14,.11),
                            habitat=rnorm(300,.19,.74),svl.habitat=rnorm(300,.05,.11))
predict.function.inf <- function(svl,hab){
  rnorm(300,-7.44,.74)+rnorm(300,.14,.11)*svl+rnorm(300,.19,.74)*hab+
    rnorm(300,.05,.11)*svl*hab
}
ggplot(data=diagnosed,aes(y=infection,x=SVL,col=habitat))+
  geom_point()+
  geom_line()

predict.function.inf(5.5,1)




library(broom)
library(tibble)
library(dotwhisker)
library(dplyr)

#infection####
table.infection <- as.data.frame(read.csv("datasheets/Parameter Estimates.csv"), header=T, stringsAsFactors=F)
table.infection <- select(table.infection, !X.1)
colnames(table.infection)[1] <- "term"
colnames(table.infection)[7] <- "model"
colnames(table.infection)[2] <- "estimate"
colnames(table.infection)[4] <- "conf.low"
colnames(table.infection)[6] <- "conf.high"
colnames(table.infection)[3] <- "std.error"

df.infection <- as_tibble(table.infection)%>% 
  filter(term != "deviance",term!="intercept")

dwplot(df.infection,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

dwplot(df.infection,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2), # plot line at zero _behind_ coefs
       whisker_args = list(aes(linetype = model))) +
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("Predicting Probability of infection") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.007, 0.01),
        legend.justification = c(0, 0),
        legend.background = element_rect(colour="grey80"),
        legend.title.align = .5)+
  scale_colour_manual(labels = c("~Habitat", "~SVL*Habitat",
                                 "~SVL", "Model Average"),
                      values = c("darkblue", "green",
                                 "orange", "red"))

df.modavg.infection <- as_tibble(table.infection)%>% 
  filter(term != "deviance",term!="intercept",model=="Model Average")
dwplot(df.modavg.infection,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

#bodycond####
table.bodycond <- as.data.frame(read.csv("datasheets/Parameter Estimates_bodycond.csv"), header=T, stringsAsFactors=F)
table.bodycond <- select(table.bodycond, !X.1)
colnames(table.bodycond)[1] <- "term"
colnames(table.bodycond)[7] <- "model"
colnames(table.bodycond)[2] <- "estimate"
colnames(table.bodycond)[4] <- "conf.low"
colnames(table.bodycond)[6] <- "conf.high"
colnames(table.bodycond)[3] <- "std.error"

df.bodycond <- as_tibble(table.bodycond)%>% 
  filter(term != "deviance")

dwplot(df.bodycond,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

dwplot(df.bodycond,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2), # plot line at zero _behind_ coefs
       whisker_args = list(aes(linetype = model))) +
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("Predicting Probability of bodycond") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.007, 0.01),
        legend.justification = c(0, 0),
        legend.background = element_rect(colour="grey80"),
        legend.title.align = .5)+
  scale_colour_manual(labels = c("~Habitat", "~SVL*Habitat",
                                 "~SVL", "Model Average"),
                      values = c("darkblue", "green",
                                 "orange", "red"))

df.modavg.bodycond <- as_tibble(table.bodycond)%>% 
  filter(term != "deviance",term!="intercept",model=="Model Average")
dwplot(df.modavg.bodycond,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

#Hemoglobin####
table.hemoglobin <- as.data.frame(read.csv("datasheets/Parameter Estimates_hemoglobin.csv"), header=T, stringsAsFactors=F)
table.hemoglobin <- select(table.hemoglobin, !X.1)
colnames(table.hemoglobin)[1] <- "term"
colnames(table.hemoglobin)[7] <- "model"
colnames(table.hemoglobin)[2] <- "estimate"
colnames(table.hemoglobin)[4] <- "conf.low"
colnames(table.hemoglobin)[6] <- "conf.high"
colnames(table.hemoglobin)[3] <- "std.error"

df.hemoglobin <- as_tibble(table.hemoglobin)%>% 
  filter(term != "deviance",term!="intercept")

dwplot(df.hemoglobin,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

dwplot(df.hemoglobin,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2), # plot line at zero _behind_ coefs
       whisker_args = list(aes(linetype = model))) +
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("Predicting Hemoglobin Levels") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.45, 0.55),
        legend.justification = c(0, 0),
        legend.background = element_rect(colour="grey80"),
        legend.title.align = .5)+
  scale_colour_manual(labels = c("~infection", "~infection+SVL",
                                 "~infection+(SVL*Sex)", "Model Average"),
                      values = c("darkblue", "green",
                                 "orange", "red"))

df.modavg.hemoglobin <- as_tibble(table.hemoglobin)%>% 
  filter(term != "deviance",term!="intercept",model=="Model Average")
dwplot(df.modavg.hemoglobin,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

#Sodium####
table.Sodium <- as.data.frame(read.csv("datasheets/Parameter Estimates_sodium.csv"), header=T, stringsAsFactors=F)
table.Sodium <- select(table.Sodium, !X.1)
colnames(table.Sodium)[1] <- "term"
colnames(table.Sodium)[7] <- "model"
colnames(table.Sodium)[2] <- "estimate"
colnames(table.Sodium)[4] <- "conf.low"
colnames(table.Sodium)[6] <- "conf.high"
colnames(table.Sodium)[3] <- "std.error"

df.Sodium <- as_tibble(table.Sodium)%>% 
  filter(term != "deviance",term!="intercept")

dwplot(df.Sodium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

dwplot(df.Sodium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2), # plot line at zero _behind_ coefs
       whisker_args = list(aes(linetype = model))) +
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("Predicting Sodium Levels") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.45, 0.55),
        legend.justification = c(0, 0),
        legend.background = element_rect(colour="grey80"),
        legend.title.align = .5)+
  scale_colour_manual(labels = c("Null","~infection", "~infection+SVL",
                                 "~infection+(SVL*Sex)", "Model Average"),
                      values = c("magenta","darkblue", "green",
                                 "orange", "red"))

df.modavg.sodium <- as_tibble(table.Sodium)%>% 
  filter(term != "deviance",term!="intercept",model=="Model Average")
dwplot(df.modavg.sodium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

#Potassium####
table.Potassium <- as.data.frame(read.csv("datasheets/Parameter Estimates_potassium.csv"), header=T, stringsAsFactors=F)
table.Potassium <- select(table.Potassium, !X.1)
colnames(table.Potassium)[1] <- "term"
colnames(table.Potassium)[7] <- "model"
colnames(table.Potassium)[2] <- "estimate"
colnames(table.Potassium)[4] <- "conf.low"
colnames(table.Potassium)[6] <- "conf.high"
colnames(table.Potassium)[3] <- "std.error"

df.Potassium <- as_tibble(table.Potassium)%>% 
  filter(term != "deviance",term!="intercept",model!="~1")

dwplot(df.Potassium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

dwplot(df.Potassium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2), # plot line at zero _behind_ coefs
       whisker_args = list(aes(linetype = model))) +
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("Predicting Potassium Levels") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.007, 0.01),
        legend.justification = c(0, 0),
        legend.background = element_rect(colour="grey80"),
        legend.title.align = .5)+
  scale_colour_manual(labels = c("~SVL","~Sex","~SVL*Sex", "~infection+SVL",
                                 "~infection+Sex", "Model Average"),
                      values = c("magenta","darkblue", "green","turquoise",
                                 "orange", "red"))

df.modavg.Potassium <- as_tibble(table.Potassium)%>% 
  filter(term != "deviance",term!="intercept",model=="Model Average")
dwplot(df.modavg.Potassium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))


#Calcium####
table.Calcium <- as.data.frame(read.csv("datasheets/Parameter Estimates_Calcium.csv"), header=T, stringsAsFactors=F)
table.Calcium <- select(table.Calcium, !X.1)
colnames(table.Calcium)[1] <- "term"
colnames(table.Calcium)[7] <- "model"
colnames(table.Calcium)[2] <- "estimate"
colnames(table.Calcium)[4] <- "conf.low"
colnames(table.Calcium)[6] <- "conf.high"
colnames(table.Calcium)[3] <- "std.error"
df.Calcium <- as_tibble(table.Calcium)%>% 
  filter(term != "deviance",term!="intercept")

dwplot(df.Calcium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

dwplot(df.Calcium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2), # plot line at zero _behind_ coefs
       whisker_args = list(aes(linetype = model))) +
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("Predicting Calcium Levels") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.55, 0.45),
        legend.justification = c(0, 0),
        legend.background = element_rect(colour="grey80"),
        legend.title.align = .5)+
  scale_colour_manual(labels = c("~SVL*Habitat","~(SVL*Sex)+(SVL*Habitat)", "~infection+(SVL*Habitat)",
                                 "~SVL", "Model Average"),
                      values = c("magenta","darkblue", "green",
                                 "orange", "red"))

df.modavg.Calcium <- as_tibble(table.Calcium)%>% 
  filter(term != "deviance",term!="intercept",model=="Model Average")
dwplot(df.modavg.Calcium,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))


#Glucose####
table.Glucose <- as.data.frame(read.csv("datasheets/Parameter Estimates_Glucose.csv"), header=T, stringsAsFactors=F)
table.Glucose <- select(table.Glucose, !X.1)
colnames(table.Glucose)[1] <- "term"
colnames(table.Glucose)[7] <- "model"
colnames(table.Glucose)[2] <- "estimate"
colnames(table.Glucose)[4] <- "conf.low"
colnames(table.Glucose)[6] <- "conf.high"
colnames(table.Glucose)[3] <- "std.error"

df.Glucose <- as_tibble(table.Glucose)%>% 
  filter(term != "deviance", term!="intercept")

dwplot(df.Glucose,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))

dwplot(df.Glucose,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2), # plot line at zero _behind_ coefs
       whisker_args = list(aes(linetype = model))) +
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  ggtitle("Predicting Probability of Glucose") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.007, 0.01),
        legend.justification = c(0, 0),
        legend.background = element_rect(colour="grey80"),
        legend.title.align = .5)+
  scale_colour_manual(labels = c("~Habitat", "~SVL*Habitat",
                                 "~SVL", "Model Average"),
                      values = c("darkblue", "green",
                                 "orange", "red"))

df.modavg.Glucose <- as_tibble(table.Glucose)%>% 
  filter(term != "deviance",term!="intercept",model=="Model Average")
dwplot(df.modavg.Glucose,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))
