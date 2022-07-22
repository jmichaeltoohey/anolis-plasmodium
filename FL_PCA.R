library(ggplot2)
devtools::install_local(path = "ggbiplot-master.zip")
library('ggbiplot')
library(factoextra)

source("FL_data_carpentry.R")

virul <- na.omit(select(diagnosed, infection, Na, K, Glu, iCa, SVL, SEX))

virul.pca <- prcomp(virul[,c(2:5)], center = TRUE,scale. = TRUE)
summary(virul.pca)

# Results for Variables
virul.var <- get_pca_var(virul.pca)
virul.var$coord          # Coordinates
virul.var$contrib        # Contributions to the PCs
virul.var$cos2           # Quality of representation 
# Results for individuals
virul.ind <- get_pca_ind(virul.pca)
virul.ind$coord          # Coordinates
virul.ind$contrib        # Contributions to the PCs
virul.ind$cos2           # Quality of representation 
data$PC1 <- NA
data$PC2 <- NA
PC1 <- virul.ind$coord[,1]
index <- as.numeric(names(PC1))
for(i in 1:length(index)){
  data$PC1[index[i]] <- virul.ind$coord[i,1]
  data$PC2[index[i]] <- virul.ind$coord[i,2]
}

#FIGURE 3
pca.plot <- ggbiplot(virul.pca, 
                     groups=factor(virul$infection),alpha=0)+
  labs(y="PC2 (30.1% explained variance)",x="PC1 (31.8% explained variance)")+
  scale_color_manual(name='Sex',values=c('0'="goldenrod3",'1'='blue3'),
                     labels=c("Female","Male"))+
  scale_fill_manual(name='Sex',values=c('0'="gold",'1'='dodgerblue2'),
                    labels=c("Female","Male"))+
  scale_shape_manual(name="Infection Status",values=c(1,21),
                     labels=c('Non-Infected','Infected'))+
  geom_point(aes(shape=groups,size=virul$SVL,
                 fill=factor(virul$SEX),color=factor(virul$SEX)))+ scale_size_continuous(name="SVL",limits = c(4, 7)) +
  theme(panel.background = element_rect(fill = "white",colour = "black"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())
#pca.plot

####Other visualizations####

#scree plot (percentage of explained variances by principle component)
fviz_eig(virul.pca)

#Graph of individuals
fviz_pca_ind(virul.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#Graph of variables
fviz_pca_var(virul.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#Biplot of individuals and variables
fviz_pca_biplot(virul.pca, repel = TRUE,
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                alpha.ind = 0
)
#Grouped by infection
ggbiplot(virul.pca, choices=c(3,4), 
         groups=factor(virul$infection),ellipse=TRUE)