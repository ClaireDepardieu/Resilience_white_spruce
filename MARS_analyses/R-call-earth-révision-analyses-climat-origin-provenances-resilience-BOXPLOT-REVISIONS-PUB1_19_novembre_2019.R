#############################################################################################
## Introduction  
#############################################################################################
#Code for 
#"Adaptive genetic variation to drought in a widely distributed conifer 
#suggests a potential for increasing forest resilience in a drying climate""
#
#Authors: Claire Depardieu, Martin P. Girardin, Simon Nadeau, Patrick Lenz, Jean Bousquet, Nathalie Isabel
#
#Journal: New Phytologist
#Article acceptance date: 29 February 2020
#
#Author for correspondence: Claire Depardieu, claire.depardieu@canada.ca or calima45@hotmail.fr 

##################################################################################
### Downloading libraries
##################################################################################
library(earth) 
require(akima)
library(spm)
library(ggplot2)
library(gridExtra) #Arrangement des graphiques sur un meme panel

#Preliminary notes on Multivariate Adaptive Regression Splines 
# MARS build linear relationship between predictor and target by segmenting predictor variables. Then a non-linear relationship can be identified by integrating all segments. 

###set R working directory once for the whole script
setwd("C:/Users/Claire/Desktop/TRAVAIL_R_DENDRO/Data_PUB_Mastigouche")
#Data_observed_values_CORRECTED_COLONY_moy_PROV.csv
Regression_data=read.table (file.choose(), header=TRUE, row.names=1, sep=";", dec=".")
head(Regression_data)
names(Regression_data)
unique(Regression_data$Corrected_provenance)
data.ready=Regression_data

#Update du 23 octobre: CORRECTED COLONY 
#J'importe un dataset qui contient les coef de correlation pour la matrice treeclim BAI versus Summer_SMI
#Data.frame.COEF.CORRECTED.COLONY
COEF_data=read.table (file.choose(), header=TRUE, sep="\t", dec=".")
head(COEF_data)
names(COEF_data)
Regression_data=COEF_data

#Pour les nouvelles donnees cliamtiques incluant database QC:
#Fichier Climatic_dataset_Mastigouche_prov_origin_revised_Martin
#Climatic_dataset_Mastigouche_prov_origin_revised_Martin
#Climat_prov_origin_PUB2_revision_elevation en date du 11 décembre
Climate_data=read.table (file.choose(), header=TRUE, row.names=1, sep=";", dec=".")
head(Climate_data)
unique(Climate_data$Nb_proveNAnce)
#Okay pour les données climatiques j'ai pas changé de numéro de prov j'ai donc juste a renommer ma colonne et that's okay
colnames(Climate_data)[colnames(Climate_data)=="Nb_proveNAnce"] <- "Corrected_provenance"
names(Climate_data)

data.ready = merge(Climate_data, Regression_data, by="Corrected_provenance")
head(data.ready)
dim(data.ready)
unique(data.ready$Corrected_provenance)

names(data.ready)
#[40] "COEF_july_radial_lumen_diameter"     
#[41] "COEF_july_avg_lumen_diameter"        
#[42] "COEF_july_wood_density"              
#[43] "COEF_july_radial_CWR"                
#[44] "COEF_BAI_july_SMI"                   
#[45] "COEF_BAI_august_SMI"                 
#[46] "COEF_BAI_september_SMI"              
#[47] "COEF_BAI_july_prec"                  
#[48] "COEF_BAI_july_Tmax"                  
#[49] "COEF_BAI_moy_curr"   

#[6] "mean_Mean_Tmean_40"                 
#[7] "mean_Total_Precipitation_40"        
#[8] "mean_Frost_Days_40"                 
#[9] "mean_DryDay_40"                     
#[10] "SMImean_summer_40"             

#Avant de runner tes modeles tu sélectionnes juste les colonnes qui vont te servir....
names(data.ready)
str(data.ready)
#DATA.EARTH=data.ready[,c(23,7,6,9,10)]
DATA.EARTH=data.ready[,c(19,10,6,7,9)]
names(DATA.EARTH)

##Si transformation des données avant application du modele....
sqrt_COEF_BAI_july_SMI<- sqrt(DATA.EARTH$COEF_BAI_july_SMI)
DATA.EARTH <- cbind(DATA.EARTH, sqrt_COEF_BAI_july_SMI)
names(DATA.EARTH)
DATA.EARTH=DATA.EARTH[,c(1,2,3,5)]
names(DATA.EARTH)

#Test du MARS
a <- earth(Resilience_2_years_2002 ~ ., data = DATA.EARTH, nfold=10, penalty = 2, thresh = 0.001, minspan = 0, endspan = 0, newvar.penalty = 0, fast.k = 20, fast.beta = 1)
plotmo(a) 
plot(a)
summary(a, digits = 5, style = "pmax")

plot(a,which=1)

evimp(a,trim=FALSE) #Displaying variable importance even for variables not included in the model

predicted = predict(a)

class(predicted)
predicted=as.data.frame(predicted)
colnames(predicted)[1] <- "predicted.values"
predicted

predicted.dataframe=cbind(DATA.EARTH,predicted)
head(predicted.dataframe)

#Plot de tes observes versus predit
plot(predicted.dataframe$Resilience_2_years_2002,predicted.dataframe$predicted.values,xlab="Observed values",ylab="Predicted values")
abline(0,1,lwd=2,col="red")

#GRAPHIQUE PUBLICATION - SHAPEE
#Plot of the model  using ggplot2
library(ggplot2)

#http://www.sthda.com/english/wiki/ggplot2-point-shapes
#Selection des données pour faire ton graphique
names(data.ready)
DATA.EARTH=data.ready[,c(10,17,6,15)]
names(DATA.EARTH)
Group_provenance <- as.factor(DATA.EARTH$Group_provenance)

df<-data.frame("x"=DATA.EARTH$SMImean_summer_40, "y"=DATA.EARTH$Estimate_Resilience_2002_Corrected_data, "group"=DATA.EARTH$Group_provenance)
p <-ggplot(df, aes(x, y,group=Group_provenance)) 
p <- p + geom_point(aes(shape=Group_provenance, color=Group_provenance, size=Group_provenance))+
  scale_shape_manual(values=c(17,16))+
  scale_color_manual(values=c('#706c66','#140e01'))+ 
  scale_size_manual(values=c(2.5,3))
p <- p + theme_bw() 
p <- p + stat_function(fun = function(x) 0.6364+0.017909*(96.285- x),xlim=c(88.8, 96.285), colour="black",size=0.8)
p <- p + stat_function(fun = function(x) 0.6364,xlim=c(96.285, 97.7), colour="black",size=0.8)
p

#b2aea6

#GRAPHIQUES PUBLICATION - SHAPPEE 2 predicteurs
#Graph 1: avec les valeurs observees--------------------------
#Pour la résilience 2002
names(DATA.EARTH)
resolution <- 0.004 # you can increase the resolution by decreasing this number (warning: the resulting dataframe size increase very quickly)
GRAPH <- interp(y=DATA.EARTH$mean_Mean_Tmean_40, x=DATA.EARTH$SMImean_summer_40, z=DATA.EARTH$Resilience_2_years_2002, 
            yo=seq(min(DATA.EARTH$mean_Mean_Tmean_40),max(DATA.EARTH$mean_Mean_Tmean_40),by=resolution), 
            xo=seq(min(DATA.EARTH$SMImean_summer_40),max(DATA.EARTH$SMImean_summer_40),by=resolution), duplicate="mean")

image(GRAPH) #you can of course modify the color palette and the color categories. See ?image for more explanation

filled.contour(GRAPH, xlim=c(88.5,98),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow","red")))
#Alternatif
#filled.contour(GRAPH,xlim=c(88.5,98),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow","red")),xlab = "Summer SMI",ylab = "MAT", key.title = title(main = "", cex.main = 1), plot.axes=c({points(95.37988,5.5967742, pch=15, cex=1.4)},{points(95.42732,6.4806452, pch=15, cex=1.4)},{points(89.07403,3.2064516, pch=16, cex=1.4)},{points(97.11463,-0.2903226, pch=16, cex=1.4)},{points(97.63837,1.9354839, pch=16, cex=1.4)}))
#filled.contour(GRAPH,xlim=c(89,98),ylim=c(-1,7),color.palette=heat.colors,xlab = "Summer SMI",ylab = "MAT", main = "Resilience_observed",key.title = title(main = "", cex.main = 1), plot.axes=c({points(95.15737,6.4258065, pch=15, cex=0.7)},{points(94.70443,5.5483871, pch=15, cex=0.7)},{points(87.18954,3.1000000, pch=16, cex=0.7)},{points(96.82255,-0.2354839, pch=16, cex=0.7)}))

                                                                                                                                                                                                                                         # PROVENANCE 25                                     # PROVENANCE 43                   # PROVENANCE 10 ici                                                         # PROVENANCE 10
#filled.contour(GRAPH,xlim=c(88.5,98),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow","red")), plot.axes=c({points(95.42732,6.4806452, pch=15, cex=0.5)},{points(95.37988,5.5967742, pch=15, cex=0.5)},{points(89.07403,3.2064516, pch=16, cex=0.5)},{points(97.11463,-0.2903226, pch=16, cex=0.5)},{points(97.63837,1.9354839, pch=16, cex=0.5)}))

# Sans les points extremes montrés                                                                                                                                                                                                                                           # PROVENANCE 25                                     # PROVENANCE 43                   # PROVENANCE 10 ici                                                         # PROVENANCE 10
filled.contour(GRAPH,xlim=c(88.5,98),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow","red")))


#pETIT CODE POUR EXTRAIRE LES COORDONNEES DES POINTS QUE TU VEUX AJOUTER SUR TON GRAPH
names(DATA.EARTH)
vecteur.points= DATA.EARTH[c(10,34,37,25,43,7),c(1,2,3)]
vecteur.points

min(DATA.EARTH$Resilience_2_years_2002)#Prov 25 min de résilience 
max(DATA.EARTH$Resilience_2_years_2002)#Prov 10 max de résilience - okay

#Graph 2: avec les valeurs predites du modele earth--------------------------

resolution <- 0.004 # you can increase the resolution by decreasing this number (warning: the resulting dataframe size increase very quickly)
GRAPH2 <- interp(y=predicted.dataframe$mean_Mean_Tmean_40, x=predicted.dataframe$SMImean_summer_40, z=predicted.dataframe$predicted.values, 
                yo=seq(min(predicted.dataframe$mean_Mean_Tmean_40),max(predicted.dataframe$mean_Mean_Tmean_40),by=resolution), 
                xo=seq(min(predicted.dataframe$SMImean_summer_40),max(predicted.dataframe$SMImean_summer_40),by=resolution), duplicate="mean")
image(GRAPH2) #you can of course modify the color palette and the color categories. See ?image for more explanation
#Alternatif

#filled.contour(GRAPH2, xlim=c(88.5,98),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow","red")))
#Alternatif
filled.contour(GRAPH2,xlim=c(88.5,98),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow","red")),xlab = "Summer SMI",ylab = "MAT", key.title = title(main = "", cex.main = 1), plot.axes=c({points(95.37988,5.5967742, pch=15, cex=1.4)},{points(95.42732,6.4806452, pch=15, cex=1.4)},{points(89.07403,3.2064516, pch=16, cex=1.4)},{points(97.11463,-0.2903226, pch=16, cex=1.4)},{points(97.63837,1.9354839, pch=16, cex=1.4)}))


# FAIRE DES BOXPLOTS 
names(data.ready)

# Boxplot 1: avec tous les sensitivity traits 
Coef=data.ready[,c(40,42,43,44,45,46,47,48,49)]
names(Coef)
names(Coef)[1]<-"CS_LDr_SMIJul"
names(Coef)[2]<-"CS_WD_SMIJul"
names(Coef)[3]<-"CS_CWRr_SMIJul"
names(Coef)[4]<-"CS_BAI_SMIJul"
names(Coef)[5]<-"CS_BAI_SMIAug"
names(Coef)[6]<-"CS_BAI_SMISep"
names(Coef)[7]<-"CS_BAI_TAPJul"
names(Coef)[8]<-"CS_BAI_TmaxJul"
names(Coef)[9]<-"CS_BAI_SMImoyen"
names(Coef)

s <- apply(Coef,2,sd)
mn <- colMeans(Coef)
ci1 <- mn - qnorm(0.95)*s
ci2 <- mn + qnorm(0.95)*s
minm <- apply(Coef, 2, min)
maxm <- apply(Coef, 2, max)

bp <- boxplot(Coef, plot=FALSE)
bp

bp$stats <- matrix(c(minm, ci1, mn, ci2, maxm), nrow= 5, byrow=TRUE)
bxp(bp)

#Faire un boxplot juste pour une variable....
boxplot(Coef$CS_LDr_SMIJul, horizontal=TRUE,main="Climate sensitivity", col="white")

# Boxplot 2: avec tous les sensitivity traits
names(data.ready)
Coef=data.ready[,c(1,40,42,43,44,45,46,47,48,49)]
names(Coef)
names(Coef)[2]<-"CS_LDr_SMIJul"
names(Coef)[3]<-"CS_WD_SMIJul"
names(Coef)[4]<-"CS_CWRr_SMIJul"
names(Coef)[5]<-"CS_BAI_SMIJul"
names(Coef)[6]<-"CS_BAI_SMIAug"
names(Coef)[7]<-"CS_BAI_SMISep"
names(Coef)[8]<-"CS_BAI_TAPJul"
names(Coef)[9]<-"CS_BAI_TmaxJul"
names(Coef)[10]<-"CS_BAI_SMImoyen"
names(Coef)

# Boxplot 3: avec BAI les sensitivity traits
names(data.ready)
Coef=data.ready[,c(1,44,45,46,49)]
names(Coef)
names(Coef)[2]<-"CS_BAI_SMIJul"
names(Coef)[3]<-"CS_BAI_SMIAug"
names(Coef)[4]<-"CS_BAI_SMISep"
names(Coef)[5]<-"CS_BAI_SMImoyen"
names(Coef)



library(reshape)
library(ggplot2)
data <- melt(Coef,id="Nb_proveNAnce")
head(data)

#Plot du graphique formatté pour la publication 
graphique <- ggplot(data, aes(x=variable, y=value)) + geom_boxplot(notch=TRUE)
P <- graphique+ stat_summary(fun.y=mean, geom="point", shape=20, size=3, colour="red")
#Rotate the box plot and add white background
P + theme_classic()
P+ coord_flip() + theme_classic()
#Définir les couleurs du boxplot
#graphique+ scale_fill_manual(values=c("blue","blue","blue","blue","blue","blue","green","red","yellow"))
# Ici mets des codes bizarres
  
graphique+ stat_summary(colour="red", size= 1.5, fun.args= list(conf.int=.95))


#Je sauvegarde un fichier 
png(file= "Boxplot_Coef_BAI_july_prec")

#Je fais mon graphique 
boxplot(DATA.EARTH$COEF_BAI_july_prec,
        xlab= "CS_BAI_TAP",
        ylab= "Correlation coefficients",
        xlim=0,
        ylim=0,
        notch=TRUE)
,
        varwidth= TRUE,
        col= "green")
#Je sauves mon fichier 
dev.off()


##################################################################################
### VIOLIN PLOTS pour h, PVE et PGE
##################################################################################
# Data pour violin plots
data.PVE=cbind(trait.type,trait,h.moyen,pve.moyen,pge.moyen,ngamma.moyen)
data.PVE
str(data.PVE)
data.PVE$trait <- as.factor(data.PVE$trait)

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Basic violin plot pour PVE
PVE <- ggplot(DATA.EARTH$COEF_BAI_july_prec, aes(x=trait, y=pve.moyen, fill=trait.type)) + 
  geom_violin(trim=FALSE) +
  #stat_summary(fun.data="mean_sdl",geom="crossbar",width=0.1) + 
  stat_summary(fun.data=data_summary) +
  labs(x="Traits", y = "PVE") +  # title="Plot of PVE",
  theme_classic() +
  scale_fill_manual(values=c("#b7b3b3", "#e5dede", "#898686", "#605c5c"))























#Alternatif
filled.contour(GRAPH2,xlim=c(98.2,99.7),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow", "red")),xlab = "Spring SMI",ylab = "MAT", main = "H1997_predicted",key.title = title(main = "", cex.main = 1))

#Ajouter des points sur le graphique 
points()

#Shapping des graphs pour publication
GRAPH.growth.observed <- interp(y=DATA.EARTH$mean_Mean_Tmean_40, x=DATA.EARTH$SMImean_summer_40, z=DATA.EARTH$Estimate_mean_sqrt_BAI_1995_2005, 
                                    yo=seq(min(DATA.EARTH$mean_Mean_Tmean_40),max(DATA.EARTH$mean_Mean_Tmean_40),by=resolution), 
                                    xo=seq(min(DATA.EARTH$SMImean_summer_40),max(DATA.EARTH$SMImean_summer_40),by=resolution), duplicate="mean")
image(GRAPH.growth.observed) #you can of course modify the color palette and the color categories. See ?image for more explanation
GRAPH.growth.pub=filled.contour(GRAPH.growth.observed,xlim=c(89,98),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow", "red")),xlab = "SMI_summer (%)",ylab = "MAT (°C)", key.title = title(main = "", cex.main = 1))
GRAPH.growth.pub

GRAPH.resilience.observed <- interp(y=DATA.EARTH$mean_Mean_Tmean_40, x=DATA.EARTH$SMImean_summer_40, z=DATA.EARTH$Estimate_mean_resilience_BAI_2002, 
                                    yo=seq(min(DATA.EARTH$mean_Mean_Tmean_40),max(DATA.EARTH$mean_Mean_Tmean_40),by=resolution), 
                                    xo=seq(min(DATA.EARTH$SMImean_summer_40),max(DATA.EARTH$SMImean_summer_40),by=resolution), duplicate="mean")
image(GRAPH.resilience.observed) #you can of course modify the color palette and the color categories. See ?image for more explanation
GRAPH.resilience.pub=filled.contour(GRAPH.resilience.observed,xlim=c(89,98),ylim=c(-1,7),color.palette=colorRampPalette(c("blue","yellow", "red")),xlab = "SMI_summer (%)",ylab = "MAT (°C)", key.title = title(main = "", cex.main = 1))
GRAPH.resilience.pub

