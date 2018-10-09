
#-----------------------------------
# 
#-----------------------------------
rm(list=ls())
library(tidyverse)



load("U:/Data/Wasting/Wasting_inc_data.RData")
load("U:/Data/Wasting/Stunting_inc_data.RData")


#--------------------------------------
# Calculate prevalence of
# Wasting in certain age ranges for the
# risk factor analysis
#--------------------------------------

# define age windows
#Cut into 3 month quartiles
d <- d %>% mutate(agemonth = agedays/30.4167)
dstunt <- dstunt %>% mutate(agemonth = agedays/30.4167)

d$agecat <- cut(d$agemonth, breaks=c(0,3,6,9,12,15,18,21,24), labels=c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "12-15 months", "15-18 months", "18-21 months", "21-24 months"), include.lowest=T, right=T)
dstunt$agecat <- cut(dstunt$agemonth, breaks=c(0,3,6,9,12,15,18,21,24), labels=c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "12-15 months", "15-18 months", "18-21 months", "21-24 months"), include.lowest=T, right=T)

d$agecat <-factor(d$agecat, labels=c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "12-15 months", "15-18 months", "18-21 months", "21-24 months"))
dstunt$agecat <-factor(dstunt$agecat, labels=c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "12-15 months", "15-18 months", "18-21 months", "21-24 months"))

#--------------------------------------
# Calculate cumulative incidence of
# wasting in certain age ranges for the
# analysis of at-birth wasting and stunting
#--------------------------------------


#Mark any wasting in age ranges
wast_ci = d %>% ungroup() %>% 
  filter(!is.na(agecat)) %>%
  arrange(studyid,country,subjid, agedays) %>% 
  group_by(studyid,country,subjid, agecat) %>%
  mutate(ever_wasted=1*(sum(wast_inc, na.rm=T)>0), Nobs=n()) %>% slice(1) %>%
  mutate(N=n()) %>%
  ungroup()   
table(wast_ci$agecat, wast_ci$ever_wasted)


#Mark any stunting in age ranges
stunt_ci = dstunt %>% ungroup() %>% 
  filter(!is.na(agecat)) %>%
  arrange(studyid,country,subjid, agedays) %>% 
  group_by(studyid,country,subjid) %>%
  mutate(cuminc=cumsum(stunt_inc), cumcuminc=cumsum(cuminc)) %>%
  filter(cumcuminc<2) %>%
  group_by(studyid,country,subjid, agecat) %>%
  mutate(ever_stunted=1*(sum(stunt_inc, na.rm=T)>0),  Nobs=n()) %>% slice(1) %>%
  mutate(N=n()) %>%
  ungroup()   
table(stunt_ci$agecat, stunt_ci$ever_stunted)


wast_ci <- wast_ci %>% subset(., select = c(subjid, studyid, country, agecat, ever_wasted))
stunt_ci <- stunt_ci %>% subset(., select = c(subjid, studyid, country, agecat, ever_stunted))

d <- merge(wast_ci, stunt_ci, by=c("subjid", "studyid", "country", "agecat"), all = T)

d <- d %>% group_by(studyid,country,subjid) %>% 
  mutate(lag_ever_wasted=lag(ever_wasted),
         lag_ever_stunted=lag(ever_stunted))

table(d$ever_stunted)
table(d$ever_wasted)

table(d$agecat, d$ever_stunted)
table(d$agecat, d$ever_wasted)

table(d$lag_ever_wasted, d$ever_stunted, d$agecat)
table(d$lag_ever_stunted, d$ever_wasted, d$agecat)

df <- d


#--------------------------------
# Merge stunting datasets and covariates
#--------------------------------

setwd("U:/ucb-superlearner/Stunting rallies/")

#load covariates
cov<-readRDS("FINAL_clean_covariates.rds")

#Merge in covariates
d <- left_join(df, cov, by=c("studyid", "subjid", "country"))


#Treatment name
A<-c("lag_ever_stunted", "lag_ever_wasted")

#Vector of outcome names
Y<-c("ever_stunted", "ever_wasted")


#Vector of covariate names
W=c("arm","sex", "W_mage", "W_fage", "meducyrs", "feducyrs", "hhwealth_quart", "hfoodsec",
    "vagbrth","hdlvry",
    "single",
    "W_nrooms","W_nhh","W_nchldlt5",
    "month","brthmon","W_parity",
    "trth2o","cleanck","impfloor","impsan","safeh20")

#Subgroup variable
V <- c("agecat")

#clusterid ID variable
id <- c("id")

#Merge in ID variable from stunting 


save(d, Y, A, W, V, id, file="U:/UCB-SuperLearner/Wasting rallies/stuntwastCI_lag_rf.Rdata")




