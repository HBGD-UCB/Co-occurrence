





rm(list=ls())
library(tidyverse)
library(reshape2)

#merge outcomes with covariates

setwd("U:/ucb-superlearner/Co-occurrence/")

#load covariates
cov<-readRDS("U:/ucb-superlearner/stunting rallies/FINAL_clean_covariates.rds")

#load outcomes
load("co_prev.RData")
load("co_cuminc.rdata")
load("co_cuminc_nobirth.rdata")
load("pers_co.rdata")
load("co_rec.rdata")


#convert subjid to character for the merge with covariate dataset
cov$subjid <- as.character(cov$subjid)
prev$subjid <- as.character(prev$subjid)
cuminc$subjid <- as.character(cuminc$subjid)
cuminc_nobirth$subjid <- as.character(cuminc_nobirth$subjid)



#------------------------------------
# Create cumulative incidence dataset
#------------------------------------

#merge in covariates
cuminc <- cuminc %>% subset(., select = -c(tr))
d <- left_join(cuminc, cov, by=c("studyid", "subjid", "country"))
head(d)



#Vector of outcome names
Y<-c("ever_coed")

#Vector of risk factor names
A<-c( "sex",               "mage",          "mhtcm",         "mwtkg",        
      "mbmi",          "single",        "fage",          "fhtcm",       
      "nrooms",      "nchldlt5",    "nhh",              
      "hhwealth_quart", "brthmon", "parity",   "meducyrs", 
      "feducyrs", "hfoodsec")  



#Vector of covariate names
W<-c("")

#Subgroup variable
V <- c("agecat")

#clusterid ID variable
id <- c("id")


save(d, Y, A,V, id,  file="co_cuminc_rf.Rdata")


#------------------------------------
# Create cumulative incidence dataset
# - no birth incidence
#------------------------------------

#merge in covariates
cuminc_nobirth <- cuminc_nobirth %>% subset(., select = -c(tr))
cuminc_nobirth <- bind_rows(cuminc_nobirth, cuminc[cuminc$agecat=="6-24 months",])

d <- left_join(cuminc_nobirth, cov, by=c("studyid", "subjid", "country"))
head(d)



#Vector of outcome names
Y<-c("ever_coed")

#Vector of risk factor names
A<-c( "gagebrth", "birthwt",    
      "birthlen", "vagbrth", "hdlvry",    
      "enwast", "anywast06", "pers_wast", 
      "trth2o", "cleanck", "impfloor",  
      "impsan", "safeh20",
      "perdiar6", "perdiar24", 
      "predfeed3", "exclfeed3", "predfeed6", "exclfeed6", "predfeed36", "exclfeed36",
      "predexfd6", "earlybf", "month")

#Vector of covariate names
W<-c("")

#Subgroup variable
V <- c("agecat")

#clusterid ID variable
id <- c("id")


save(d, Y, A,V, id,  file="co_cuminc_nobirth_rf.Rdata")


#------------------------------------
# Create prevalence dataset
#------------------------------------


#merge in covariates
d <- left_join(prev, cov, by=c("studyid", "subjid", "country"))
head(d)


#Vector of outcome names
Y<-c("wasted","swasted")


A<-c( "sex",              "gagebrth",      "birthwt",      
      "birthlen",      "vagbrth",       "hdlvry",        "mage",          "mhtcm",         "mwtkg",        
      "mbmi",          "single",        "fage",          "fhtcm",         "nrooms",        "nhh",           "nchldlt5",     
      "hhwealth_quart", "month", "brthmon", "parity",   "meducyrs", 
      "feducyrs", "hfoodsec",  
      "enwast", "anywast06", "pers_wast", 
      "trth2o", "cleanck", "impfloor",  "impsan", "safeh20",
      "perdiar6", "perdiar24", "predexfd6", 
      "predfeed3", "exclfeed3", "predfeed6", "exclfeed6", "predfeed36", "exclfeed36",
      "earlybf")  



save(d, Y, A,V, id,  file="co_prev_rf.Rdata")

