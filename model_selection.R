source("all models.R")
library(MuMIn)

####Create Model Selection Tables####
modelsel.infection <- model.sel(infection.null, infection.svl, infection.sex, infection.habitat, infection.svl_sex,
                                infection.svl_X_habitat, infection.svl_X_sex__t__svl_X_habitat)

modelsel.bodycond <- model.sel(bodycond.null, bodycond.infection, bodycond.svl,  bodycond.habitat, 
                              bodycond.svl_X_habitat,
                              bodycond.infection.svl, bodycond.infection.habitat,
                              bodycond.infection.svl_X_habitat)

modelsel.temp.diff <- model.sel(temp.diff.null, temp.diff.infection, temp.diff.svl, temp.diff.sex, temp.diff.habitat, 
                                temp.diff.svl_sex, temp.diff.svl_X_habitat, temp.diff.sex_t_habitat,
                                temp.diff.svl_X_sex__t__svl_X_habitat,
                                temp.diff.infection.svl, temp.diff.infection.sex, temp.diff.infection.habitat,
                                temp.diff.infection.svl_sex, temp.diff.infection.svl_X_habitat, temp.diff.infection.sex_t_habitat,
                                temp.diff.infection.svl_X_sex__t__svl_X_habitat)

modelsel.Na <- model.sel(Na.null, Na.infection, Na.svl, Na.sex, Na.habitat, 
                         Na.svl_sex, Na.svl_X_habitat, Na.sex_t_habitat,
                         Na.svl_X_sex__t__svl_X_habitat,
                         Na.infection.svl, Na.infection.sex, Na.infection.habitat,
                         Na.infection.svl_sex, Na.infection.svl_X_habitat, Na.infection.sex_t_habitat,
                         Na.infection.svl_X_sex__t__svl_X_habitat)

modelsel.K <- model.sel(K.null, K.infection, K.svl, K.sex, K.habitat, 
                        K.svl_sex, K.svl_X_habitat, K.sex_t_habitat,
                        K.svl_X_sex__t__svl_X_habitat,
                        K.infection.svl, K.infection.sex, K.infection.habitat,
                        K.infection.svl_sex, K.infection.svl_X_habitat, K.infection.sex_t_habitat,
                        K.infection.svl_X_sex__t__svl_X_habitat)

modelsel.Cl <- model.sel(Cl.null, Cl.infection, Cl.svl, Cl.sex, Cl.habitat, 
                         Cl.svl_sex, Cl.svl_X_habitat, Cl.sex_t_habitat,
                         Cl.svl_X_sex__t__svl_X_habitat,
                         Cl.infection.svl, Cl.infection.sex, Cl.infection.habitat,
                         Cl.infection.svl_sex, Cl.infection.svl_X_habitat, Cl.infection.sex_t_habitat,
                         Cl.infection.svl_X_sex__t__svl_X_habitat)

modelsel.iCa <- model.sel(iCa.null, iCa.infection, iCa.svl, iCa.sex, iCa.habitat, 
                          iCa.svl_sex, iCa.svl_X_habitat, iCa.sex_t_habitat,
                          iCa.svl_X_sex__t__svl_X_habitat,
                          iCa.infection.svl, iCa.infection.sex, iCa.infection.habitat,
                          iCa.infection.svl_sex, iCa.infection.svl_X_habitat, iCa.infection.sex_t_habitat,
                          iCa.infection.svl_X_sex__t__svl_X_habitat)

modelsel.Hct <- model.sel(Hct.null, Hct.infection, Hct.svl, Hct.sex, Hct.habitat, 
                          Hct.svl_sex, Hct.svl_X_habitat, Hct.sex_t_habitat,
                          Hct.svl_X_sex__t__svl_X_habitat,
                          Hct.infection.svl, Hct.infection.sex, Hct.infection.habitat,
                          Hct.infection.svl_sex, Hct.infection.svl_X_habitat, Hct.infection.sex_t_habitat,
                          Hct.infection.svl_X_sex__t__svl_X_habitat)

modelsel.Glu <- model.sel(Glu.null, Glu.infection, Glu.svl, Glu.sex, Glu.habitat, 
                          Glu.svl_sex, Glu.svl_X_habitat, Glu.sex_t_habitat,
                          Glu.svl_X_sex__t__svl_X_habitat,
                          Glu.infection.svl, Glu.infection.sex, Glu.infection.habitat,
                          Glu.infection.svl_sex, Glu.infection.svl_X_habitat, Glu.infection.sex_t_habitat,
                          Glu.infection.svl_X_sex__t__svl_X_habitat)

####View Tables####
modelsel.infection
modelsel.bodycond
modelsel.temp.diff
modelsel.Na
modelsel.K
modelsel.Cl
modelsel.iCa
modelsel.Hct
modelsel.Glu

#Model Averaging
#infection####
#models with delta.aicc < 2
infection.modavg <- model.avg(modelsel.infection, subset = delta < 2)
summary(infection.modavg)
#or as a 95% confidence set:
avgmod.95p <- model.avg(modelsel.infection, cumsum(weight) <= .95)
confint(avgmod.95p)
new.infection.data <- predict(infection.modavg, full=TRUE)

#DIAGNOSTICS
#plot regressors pairwise
diagnosed %>%
  select(infection, SVL, SEX, habitat) %>%
  plot()
#check regressor correlations
diagnosed %>%
  select(infection, SVL, SEX, habitat) %>%
  cor() %>%
  round(3)
#test for collinearity
infection.additive <- glm(diagnosed$infection ~ diagnosed$SVL+diagnosed$habitat,
                          family=binomial)
vif(infection.additive)
#residual plots
scatter.smooth(fitted(infection.modavg), residuals(infection.modavg, type = "pearson"),
               mgp = c(2.2, 1, 0),
               ylab = "Residuals (Pearson)",
               xlab = "Predicted Infection")
title("Residual plot", line = 0.7)
abline(h = 0, col="blue")
