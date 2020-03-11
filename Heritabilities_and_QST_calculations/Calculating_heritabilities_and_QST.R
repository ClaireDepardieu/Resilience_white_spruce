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

#############################################################################################
## Downloading libraries  
#############################################################################################
library(asreml)
library(asremlPlus)
library(tidyverse)

#############################################################################################
## First dataset: Wood anatomy 
## Repeated Measures Analysis  
## Results are presented in Table 1 (see Depardieu et al., 2020)
#############################################################################################

#Importing data 
#Choose the file "Data_phenotype_wood_anatomy.csv"
data.phenotype.WA=read.table (file.choose(), header=TRUE, row.names=1, sep=";", dec=".")
summary(data.phenotype.WA)
names(data.phenotype.WA)

##Description of the variables----------------------------------------------------
#Description of the 7 Variables tested 
#Below: abbreviations of the variables are indicated in parenthesis - Please refer to Table 1 in Depardieu et al., 2020 for details 

#bai (BAI) : Basal area increment, annual-based radial growth of trees
#avg_double_cell_wall_thickness (CWT): Double cell wall thickness. Average values of annual rings.
#avg_density (WD): Wood density. Average values of annual rings.
#lumen_diameter_radial (LDr): Radial lumen diameter. Average values of annual rings.
#avg_lumen_diameter (LD): Lumen diameter.  
#CWR_radial (CWRr): 
#avg_CWR_100 (CWR): 


#Setting variables as factors or continous variables 
data.phenotype.WA$Family<-as.factor(data.phenotype.WA$Family)
data.phenotype.WA$Provenance<-as.factor(data.phenotype.WA$Provenance)
data.phenotype.WA$yearf<-as.factor(data.phenotype.WA$year)# Year as a factor
data.phenotype.WA$yearc<-data.phenotype.WA$year # Year as a variate
data.phenotype.WA$tree_id<-as.factor(data.phenotype.WA$tree_id)
data.phenotype.WA$Bloc<-as.factor(data.phenotype.WA$Bloc)

#Order dataframe
data.phenotype.WA = data.phenotype.WA[order(data.phenotype.WA$tree_id, data.phenotype.WA$year),]
str(data.phenotype.WA)

#If applicable, transformation of the variable 
## A log transformation was applied to bai, lumen_diameter_radial, 
##In this script: example with lumen_diameter_radial

log_lumen_diameter_radial<- log(data.phenotype.WA$lumen_diameter_radial)
#log_avg_CWR_100 <- log(data.phenotype.WA$avg_CWR_100)
#log_bai <- log(data.phenotype.WA$bai)

data.phenotype.WA <- cbind(data.phenotype.WA,log_lumen_diameter_radial)
names(data.phenotype.WA)


#Selection of the final model and calculation of heritability and QST estimates 

model_FINAL_1<-asreml(fixed=log_lumen_diameter_radial~yearf + yearf:Bloc,
                                                 random=~Provenance+Provenance:Family + yearf:Family:Bloc + Provenance:yearf + Provenance:Family:yearf,
                                                 residual=~tree_id:ar1h(yearf),data=data.phenotype.WA, maxiter = 100, workspace = "1000mb")
summary(model_FINAL_1)$varcomp

model_FINAL_2<-asreml(fixed=log_lumen_diameter_radial~yearf + yearf:Bloc,
                         random=~Provenance+Provenance:Family + diag(yearf):Family:Bloc + Provenance:yearf + Provenance:Family:yearf,
                         residual=~tree_id:ar1h(yearf),data=data.phenotype.WA, maxiter = 100, workspace = "1000mb")
summary(model_FINAL_2)$varcomp

lrt.asreml(model_FINAL_1, model_FINAL_2, boundary = TRUE)

#Running the final model and calculation of heritability and QST estimates
model_FINAL<-asreml(fixed=log_lumen_diameter_radial~yearf + yearf:Bloc,
                      random=~Provenance+Provenance:Family + diag(yearf):Family:Bloc + Provenance:yearf + Provenance:Family:yearf,
                      residual=~tree_id:ar1h(yearf),data=data.phenotype.WA, maxiter = 100, workspace = "1000mb")

##Testing the significance of the fixed effects of the final model
wald(model_FINAL)

##Extraction of variance components
summary(model_FINAL)$varcomp

##Calculating heritability estimate
(h2 = vpredict(model_FINAL, h2~(3.88694433*V4)/(V1+V2+V3+V4+V5+V6+(V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19)/11)))

##Calculating QST estimate
(QST = vpredict(model_FINAL, QST~(V1)/(V1+3.88694433*2*V4)))

#Checking the distribution of the residuals 
##(1) Is the response variable a reasonably linear function of the fitted values?
plot(fitted(model_FINAL),data.phenotype.WA$log_lumen_diameter_radial)
abline(0,1)

##(2) Are the errors reasonably close to normally distributed?
plot(fitted(model_FINAL), resid(model_FINAL))
abline(0,0)
qqnorm(resid(model_FINAL))
qqline(resid(model_FINAL))
resid = residuals(model_FINAL)
hist(resid)
plot(resid)


#############################################################################################
## Second dataset: Growth resilience
## Growth resilience indices calculated in 2002
## Results are presented in Table 1 (see Depardieu et al., 2020)
#############################################################################################
#Importing data 
#Choose the file "Data_resilience.csv"
data.resilience=read.table (file.choose(), header=TRUE, row.names=1, sep=";", dec=".")
summary(data.resilience)
names(data.resilience)

##Description of the variables----------------------------------------------------
#Description of the 4 Variables tested 
#Below: abbreviations of the variables are indicated in parenthesis - Please refer to Table 1 in Depardieu et al., 2020 for details 

#Recovery_2002 (Rc2002) : Growth recovery
#Relative_resilience_2002 (Rr2002): Growth relative resilience 
#Resilience_2002 (Rl2002): Growth resilience
#Resistance_2002 (Rs2002): Growth resistance

#Setting variables as factors 
data.resilience$Family<-as.factor(data.resilience$Family)
data.resilience$Provenance<-as.factor(data.resilience$Provenance)
data.resilience$Bloc<-as.factor(data.resilience$Bloc)

model_full<-asreml(fixed=Resilience_2002~Bloc,
                random=~Family:Bloc+Provenance:Bloc,
                data=data.resilience, maxiter = 10000)

wald(model_full)

summary(model_full)$varcomp


#Variance components
(v <- summary(model_full)$varcomp) #Corrected_prov*bloc est pas important tres petit z ratio
lrt(model_4,model_3)

