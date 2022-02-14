###PCA ANALYSIS
library(devtools)
devtools::install_local(path = "ggbiplot-master.zip")
library('ggbiplot')
library(factoextra)
library(ggplot2)

source("data_carpentry.R")

#select just virulence data columns
library(dplyr)
virul <- na.omit(select(diagnosed, infection, Na, K, Glu, iCa, SVL, SEX))
#virul <- na.omit(select(diagnosed, infection, K, bodycond))

virul.pca <- prcomp(virul[,c(2:5)], center = TRUE,scale. = TRUE)

summary(virul.pca)


#plotting
fviz_eig(virul.pca)

fviz_pca_ind(virul.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(virul.var)
fviz_pca_var(virul.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(virul.pca, repel = TRUE,
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                alpha.ind = 0
)+
  scale_color_manual(name='Sex',values=c('0'="deeppink3",'1'='cornflowerblue'),
                     labels=c("Female","Male"))+
  scale_fill_manual(name='Sex',values=c('0'="deeppink3",'1'='cornflowerblue'),
                    labels=c("Female","Male"))+
  scale_shape_manual(values=c(1,21))+
  geom_point(aes(shape=factor(virul$infection),size=virul$SVL,
                 fill=factor(virul$SEX),color=factor(virul$SEX)))

ggbiplot(virul.pca)
ggbiplot(virul.pca, choices=c(3,4), 
         groups=factor(virul$infection),ellipse=TRUE)

diagnosed <- data[136:486,] %>% filter(!is.na(infection))
for(i in 1:nrow(diagnosed)){
  if(diagnosed$SEX[i] == "F"){diagnosed$SEX[i] = 0}
  else if(diagnosed$SEX[i] == "M"){diagnosed$SEX[i] = 1}}
virul <- na.omit(select(diagnosed, infection, Na, K, Glu, iCa, SVL, SEX))
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

#Access to PCA results
library(factoextra)
# Eigenvalues
eig.val <- get_eigenvalue(virul.pca)
eig.val

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


#Theory behind PCA results
# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- virul.pca$rotation
sdev <- virul.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord[, 1:5])

# Compute Cos2
#::::::::::::::::::::::::::::::::::::::::
var.cos2 <- var.coord^2
head(var.cos2[, 1:5])

# Compute contributions
#::::::::::::::::::::::::::::::::::::::::
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:5])
