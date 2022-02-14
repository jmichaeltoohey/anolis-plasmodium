#FIGURES
library(glmmfields)
library(cowplot)


####Spatial plots - FL####
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

####Spatial plots - PR####
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

####Combine FL and PR spatial plots for final figure####
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

