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
library("tidyverse")  
library("hierfstat")
library(adegenet)
library(dplyr)

##################################################################################
### Analysis: basic statistics for individual SNPs (6,3)
##################################################################################
#Importing the data 
file1 <- "C:/Users/...working_directory.../Data.genotype.csv"
data.genotype <- read_delim(file1, delim = ";") 
fix(data.genotype)
dim(data.genotype)
#Description of the dataset imported:
#First column of the dataset: Provenance
#Other columns: 6,386 validated SNPs
#1481 rows: 1481 trees

# Format the dataframe... 
class(data.genotype) <- "data.frame"

#Ici je demande des stats basiques 
BASIC.STATS <- basic.stats(data.genotype, diploid=TRUE)
BASIC.STATS$overall

#Results obtained ---------------------
#Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis 
#0.3040  0.2922  0.3055  0.0133  0.3059  0.0137  0.0437  0.0447 -0.0405 
#Dest 
#0.0193 

basic.data.perloc=as.data.frame(BASIC.STATS2$perloc)
write.table(basic.data.perloc,"Basic_statistics.txt")
