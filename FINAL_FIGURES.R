#FIGURES
library(dotwhisker)
library(glmmfields)
library(cowplot)
#p(Infection) coefficient dot plot - all models
dwplot(df.infection,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2), # plot line at zero _behind_ coefs
       whisker_args = list(aes(color = model))) +
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

#Hemoglobin coefficient dot plot - all models
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

#Spatial plots - FL
pred_grid_FL <- expand.grid(lat = seq(min(d.EnvF$lat), max(d.EnvF$lat), length.out = 50),
                         lon = seq(min(d.EnvF$lon), max(d.EnvF$lon), length.out = 50))
pred_grid_FL$sex <- 1
pred_grid_FL$prediction <- predict(m.EnvF_spatial.sex, newdata = pred_grid_FL,
                                type = "response")$estimate

#rotate by 45 degrees
degree <- pi * 45 / 180
streams <- as.data.frame(read.csv(file = "datasheets/STREAMS.csv", fileEncoding = 'UTF-8-BOM'))
l <- sqrt(streams$lat^2 + streams$lon^2)
teta <- atan(streams$lon / streams$lat)
streams$lat <- l * cos(teta - degree)
streams$lon <- l * sin(teta - degree)
#convert to 0-100 axes
streams$lon = (((streams$LON - min(streams$LON)) * (100 - 0)) / (max(streams$LON) - min(streams$LON))) + 0
streams$lat = (((streams$LAT - min(streams$LAT)) * (100 - 0)) / (max(streams$LAT) - min(streams$LAT))) + 0

spatial_model_FL <- ggplot(pred_grid_FL, aes(lon, lat, fill = prediction)) +
  geom_raster() +
  viridis::scale_fill_viridis("p(Infection)")+
  #ggtitle("Spatial predictions of p(Infection) - Florida")+
  xlab('X Coordinate')+ylab('Y Coordinate')+ggtitle("Predicted")+
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(face="bold"))
library(ggnewscale)
spatial_data_FL <- ggplot(d.EnvF, aes(lon, lat)) +
  scale_fill_manual(name='',values=c("Streambed"='dodgerblue3'))+
  geom_point(data=streams, aes(lat, lon,fill="Streambed"),shape=22,size=2)+
  geom_segment(data=streams, aes(x=lat,y=lon,
                                 xend=c(tail(lat, n=-1),NA),yend=c(tail(lon,n=-1),NA)),
               size=3,color='dodgerblue3')+
  new_scale_colour()+
  viridis::scale_colour_viridis(discrete=TRUE) +
  geom_point(aes(colour = as.factor(y)),size = 4)+
  geom_point(shape=1,size=4,colour="black")+
  #ggtitle("Observed infection status - Florida")+
  #ggtitle("Florida")+
  panel_border(color = "black", size = 1, linetype = 1, remove = FALSE)+
  xlab('X Coordinate')+ylab('Y Coordinate')+ggtitle("Observed")+
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(face="bold"))

plot_grid(spatial_data_FL,spatial_model_FL,labels="AUTO",nrow=1)

#Spatial plots - PR
pred_grid_PR <- expand.grid(lat = seq(min(d$lat), max(d$lat), length.out = 50),
                         lon = seq(min(d$lon), max(d$lon), length.out = 50))
pred_grid_PR$svl <- mean(d$svl)
pred_grid_PR$prediction <- predict(m_spatial.svl, newdata = pred_grid_PR,
                                type = "response")$estimate
spatial_model_PR <- ggplot(pred_grid_PR, aes(lon, lat, fill = prediction)) +
  geom_raster() +
  viridis::scale_fill_viridis("p(Infection)")+
  #ggtitle("Spatial predictions of p(Infection) - Puerto Rico")+
  xlab('X Coordinate')+ylab('Y Coordinate')+
  theme(legend.position="none",axis.text.x=element_blank(),
                                                   axis.text.y=element_blank(),
                                                   panel.background = element_blank())
spatial_data_PR <- ggplot(d, aes(lon, lat, colour = y)) +
  viridis::scale_colour_viridis() +
  geom_point(aes(colour = y),size = 4)+
  geom_point(shape=1,size=4,colour="black")+
  #ggtitle("Observed infection status - Puerto Rico")+
  #ggtitle("Florida")+
  panel_border(color = "black", size = 1, linetype = 1, remove = FALSE)+
  xlab('X Coordinate')+ylab('Y Coordinate')+ theme(legend.position="none",axis.text.x=element_blank(),
                               axis.text.y=element_blank(),
                               panel.background = element_blank())
spatial_PR_combined <- plot_grid(spatial_data_PR,spatial_model_PR,labels="AUTO",nrow=1)
# extract a legend that is laid out horizontally
legend_a <- get_legend(
  spatial_data_FL + 
    guides(color = guide_legend("Infection Status",ncol = 1)) +
    theme(legend.position = "bottom")
)
legend_b <- get_legend(
  spatial_data_PR + 
    guides(color = guide_colorbar("p(Infection)",nrow = 1)) +
    theme(legend.position = "bottom")
          #legend.box.margin = margin(0, 0, 0, 0)) 
          #legend.box.margin = margin(0, 0, 0, 6))
    #theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
plot_grid(spatial_PR_combined, legend_b, 
          #ncol = 1, 
          rel_widths = c(3, .4))
spatial_four_combined <- plot_grid(spatial_data_FL,spatial_model_FL,
          spatial_data_PR,spatial_model_PR,
          labels=c("FL","","PR",""),
          label_x=0,label_y=0,
          hjust = -0.5, vjust = -0.5,nrow=2)
spatial_four_combined <- plot_grid(spatial_four_combined,legend_b,ncol=1,rel_heights = c(1,.1))
spatial_four_combined

spatial_data_combined <- plot_grid(spatial_data_FL,spatial_data_PR,legend_a,
                                   labels=c("FL","PR",""),label_x=0,label_y=0,
                                   hjust = -0.5, vjust = -0.5,nrow=3,
                                   rel_heights = c(1,1,.2))
spatial_model_combined <- plot_grid(spatial_model_FL,spatial_model_PR,legend_b,
                                   labels=c("","",""),label_x=0,label_y=0,
                                   hjust = -0.5, vjust = -0.5,nrow=3,
                                   rel_heights = c(1,1,.2))
spatial_four_combined <- plot_grid(spatial_data_combined,spatial_model_combined,
                                   ncol=2)

#Observed data visualization
library(ggplot2)
df.prevalence.grouped <- data.frame(habitat=c("Forested","Non-Forested"),
                            sex=rep(c("M","F"),each=2),
                            prevalence=c(prev.Fo.M,prev.Fo.F,
                                         prev.U.M,prev.U.F))
t.habitat <- t.test(diagnosed$infection~diagnosed$habitat)
t.sex <- t.test(diagnosed$infection~diagnosed$SEX)
t.svl <- t.test(diagnosed$infection,diagnosed$SVL)
t.sex.urban <- t.test(diagnosed.urban$infection~diagnosed.urban$SEX)
t.sex.forest <- t.test(diagnosed.forest$infection~diagnosed.forest$SEX)
#habitat grouped by sex
ggplot(df.prevalence.grouped, aes(x=habitat,y=prevalence,fill=sex))+
  geom_bar(stat="identity",color="black", position="dodge")+
  scale_fill_manual(values=c("darkgreen","darkgrey"))+ 
  #theme(legend.position="none")+
  ylim(0,1)+ 
  # Between urban sexes
  geom_signif(
    y_position = c(.65, .475), xmin = c(0.8, 1.8), xmax = c(1.2, 2.2),
    annotation = c("**", "NS")
  ) +
  geom_signif(
    comparisons = list(c("Forested", "Non-Forested")),
    y_position = .9, vjust = 0,
    annotation = "p = 0.0005621"
  )
#habitat####
df.prevalence.habitat <- data.frame(habitat=c("Forested","Non-Forested"),
                                    prevalence=c(prevalence.forest,
                                                 prevalence.urban))
barplot.habitat <- ggplot(df.prevalence.habitat, aes(x=habitat,y=prevalence,fill=habitat))+
  geom_bar(stat="identity",color="black",size=1)+
  scale_fill_manual(values=c("darkgreen","darkgrey"))+ 
  #theme(legend.position="none")+
  ylim(0,1)+ 
  # Between urban sexes
  geom_signif(
    comparisons = list(c("Forested", "Non-Forested")),
    y_position = .75, vjust = 0, tip_length=c(1.175,2.175),
    annotation = "p = 0.0005621",
    size=1
  )+theme_cowplot()+ theme(legend.position="none",
                                 axis.line.x = element_line(size = 1),
                                 axis.line.y = element_line(size = 1))+
  xlab("Habitat")+ylab("Prevalence")
#sex####
df.prevalence.sex <- data.frame(sex=c("M","F"),
                                    prevalence=c(prevalence.males,
                                                 prevalence.females))
barplot.sex <- ggplot(df.prevalence.sex, aes(x=sex,y=prevalence,fill=sex))+
  geom_bar(stat="identity",color="black",size=1)+
  scale_fill_manual(values=c("hotpink","dodgerblue"))+ 
  #theme(legend.position="none")+
  ylim(0,1)+ 
  # Between urban sexes
  geom_signif(
    comparisons = list(c("M", "F")),
    y_position = .75, vjust = 0, tip_length=c(3.45,2.45),
    annotation = "p = 0.04266",
    size=1
  )+theme_cowplot()+ theme(legend.position="none",
                           axis.line.x = element_line(size = 1),
                           axis.line.y = element_line(size = 1))+
  xlab("Sex")+ylab("Prevalence")
barplot.sex

#svl####
t.svl <- t.test(diagnosed$infection,diagnosed$SVL)
t.svl.m <- t.test(diagnosed.males$infection,diagnosed.males$SVL)
t.svl.f <- t.test(diagnosed.females$infection,diagnosed.females$SVL)
ggplot(diagnosed,aes(y=SVL,x=as.factor(infection)))+
  geom_boxplot(alpha=.5,position="dodge")+
  #facet_wrap(~SEX)+
  #theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c('No','Yes'))+
  xlab("Infection Status")+theme_cowplot(12)+
  ylim(c(3.5,8))+
  theme(strip.background = element_blank(), strip.text.x = element_blank())
#0.01149
infection.svl = glm(infection ~ SVL,family=binomial, data=diagnosed)
#ggPredict(infection.svl,se=TRUE,interactive=TRUE)

plot.svl <- ggplot(diagnosed,aes(x=SVL,y=infection))+geom_jitter(height=.075)+
  geom_smooth( method="glm",
               method.args=list(family="binomial"))+theme_cowplot()+
  theme(axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1))+
  scale_y_continuous(breaks=c(0,.25,.5,1)) + ylab("Infection")

raw.data.plots <- plot_grid(barplot.habitat,barplot.sex,plot.svl,labels="AUTO",nrow=1)
ggsave("raw_data_plots.jpg",raw.data.plots,width=6.5,dpi=320)

install.packages('cowplot')
library(cowplot)
install.packages('devtools')
library(devtools)
install_github('azzoam/ggBrackets')
library(ggplot2)
library(ggBrackets)
library(ggPredict)
devtools::install_github("cardiomoon/ggiraphExtra")
require(ggiraph)
require(ggiraphExtra)
require(plyr)
library(ggpmisc)

library(ltm)
biserial.cor(x=diagnosed$SVL,y=diagnosed$habitat)
biserial.cor(x=diagnosed$SVL,y=diagnosed$SEX)
