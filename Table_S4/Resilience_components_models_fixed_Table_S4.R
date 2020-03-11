###############################
## Introduction
###############################
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

###############################
## Downloading libraries  
###############################
library(asreml)
library(tidyverse)

#############################################################################################
## Pour mes données de résilience ici - PAS DE STRUCTURE REPETEE DANS LE TEMPS ICI
#############################################################################################
#Import the file Data.resilience
Data_resilience=read.table (file.choose(), header=TRUE, row.names=1, sep=";", dec=".")
head(Data_resilience)
names(Data_resilience)
dim(Data_resilience)

#Definition des variables facteur
Data_resilience$Provenance =as.factor(Data_resilience$Provenance)
Data_resilience$Family =as.factor(Data_resilience$Family)
Data_resilience$Bloc =as.factor(Data_resilience$Bloc)

#######################################################################################
### Resilience components in 2002 Testing effets of provenance, family and tree size 
#######################################################################################
names(Data_resilience)
#2002 resilience components...
##Recovery_2002: growth recovery calculated in 2002
##Resilience_2002: growth resilience calculated in 2002
##Resistance_2002: growth resistance calculated in 2002
##Relative_resilience_2002: relative resilience calculated in 2002
##Sum_BAI_2002: tree size in 2002
##Provenance: provenance effect (random)
##Family: family effect (random)

##Example with Recovery_2002: see above 
model_full<-asreml(fixed=Recovery_2002~Bloc+Sum_BAI_2002,
                    random=~Provenance+Provenance:Family+Family:Bloc,
                    data=Data_resilience, maxiter = 1000, workspace = "10000mb")
wald(model_full) #testing the significance of fixed effects

model_without_provenance<-asreml(Recovery_2002~Bloc+Sum_BAI_2002,
                    random=~Provenance:Family+Family:Bloc,
                    data=Data_resilience, maxiter = 1000, workspace = "10000mb")

wald(model_without_provenance)  #testing the significance of fixed effects

#Significance of the provenance effect?
lrt.asreml(model_without_provenance,model_full, boundary=TRUE) #REML likelihood ratio test

model_without_fam_prov<-asreml(fixed=Recovery_2002~Bloc+Sum_BAI_2002,
                                 random=~Provenance+Family:Bloc,
                                 data=Data_resilience, maxiter = 1000, workspace = "10000mb")

lrt.asreml(model_without_fam_prov,model_full, boundary=TRUE)

model_without_fam_bloc<-asreml(fixed=Recovery_2002~Bloc+Sum_BAI_2002,
                               random=~Provenance+Provenance:Family,
                               data=Data_resilience, maxiter = 1000, workspace = "10000mb")

lrt.asreml(model_without_fam_bloc,model_full, boundary=TRUE)

#Checking normality of residuals...
##Is the response variable a reasonably linear function of the fitted values?
plot(fitted(model_full), Data_resilience$Recovery_2002)
abline(0,1)

##Are the errors reasonably close to normally distributed?
plot(fitted(model_full), resid(model_full))
abline(0,0)
qqnorm(resid(model_full))
qqline(resid(model_full))
resid = residuals(model_full)
hist(resid)
plot(resid)


#######################################################################################
### Mean resilience components (mean for 1997, 2002 and 2005 drought events)
###Testing effets of provenance, family and tree size 
#######################################################################################
names(Data_resilience)
#2002 resilience components...
##mean_recovery: growth recovery calculated in 2002
##mean_resilience: growth resilience calculated in 2002
##mean_resistance: growth resistance calculated in 2002
##mean_rel.resil: relative resilience calculated in 2002
##Avg_Sum_BAI_3_droughts: tree size
##Provenance: provenance effect (random)
##Family: family effect (random)

##Example with mean_recovery: see above 
model_full<-asreml(fixed=mean_recovery~Bloc+Avg_Sum_BAI_3_droughts,
                   random=~Provenance+Provenance:Family+Family:Bloc,
                   data=Data_resilience, maxiter = 1000, workspace = "10000mb")
wald(model_full) #testing the significance of fixed effects

model_without_provenance<-asreml(mean_recovery~Bloc+Avg_Sum_BAI_3_droughts,
                                 random=~Provenance:Family+Family:Bloc,
                                 data=Data_resilience, maxiter = 1000, workspace = "10000mb")

wald(model_without_provenance)  #testing the significance of fixed effects

#Significance of the provenance effect?
lrt.asreml(model_without_provenance,model_full, boundary=TRUE) #REML likelihood ratio test

model_without_fam_prov<-asreml(fixed=mean_recovery~Bloc+Avg_Sum_BAI_3_droughts,
                               random=~Provenance+Family:Bloc,
                               data=Data_resilience, maxiter = 1000, workspace = "10000mb")

lrt.asreml(model_without_fam_prov,model_full, boundary=TRUE)

model_without_fam_bloc<-asreml(fixed=mean_recovery~Bloc+Avg_Sum_BAI_3_droughts,
                               random=~Provenance+Provenance:Family,
                               data=Data_resilience, maxiter = 1000, workspace = "10000mb")

lrt.asreml(model_without_fam_bloc,model_full, boundary=TRUE)

#Checking normality of residuals...
##Is the response variable a reasonably linear function of the fitted values?
plot(fitted(model_full), Data_resilience$mean_recovery)
abline(0,1)

##Are the errors reasonably close to normally distributed?
plot(fitted(model_full), resid(model_full))
abline(0,0)
qqnorm(resid(model_full))
qqline(resid(model_full))
resid = residuals(model_full)
hist(resid)
plot(resid)

