


#-----------------------------------
# Co-occurence analysis
# Objective 1a
# Calculate cumulative incidence 
# (ever stunted and wasted) at
# 3, 6, 12, 18, and 24 mo of age

# Cumulative incidence pooled using random effects
# Wasting and stunting must occur simultaneously
#-----------------------------------
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(metafor)
theme_set(theme_bw())

# load base functions
source("U:/Scripts/Stunting/1-outcomes/0_st_basefunctions.R")

load("U:/Data/Co-occurrence/co-occurrence_data.RData")



#Subset analysis to monthly studies
d <- d %>% filter(measurefreq!="yearly")

# define age windows
d = d %>% 
  mutate(agecat=
           ifelse(agedays>0 & agedays<=6*30.4167,"0-6 months",
                  ifelse(agedays>6*30.4167 & agedays<=24*30.4167,"6-24 months",NA))) %>%
  mutate(agecat=factor(agecat,levels=c("0-6 months","6-24 months")))




#--------------------------------------------------------------
# identify children ever stunted + wasted at the same time -
# Within specefic age ranges
#--------------------------------------------------------------


#calculate any stunting from 0-6
co_ci_0_6 = d %>% ungroup() %>% 
  arrange(studyid,country,subjid, agedays) %>% 
  filter(agecat=="0-6 months") %>%
  group_by(studyid,country,subjid) %>%
  mutate( minhaz=min(haz), anystunt=ifelse(minhaz< -2,1,0), 
         minwhz=min(whz), anywast=ifelse(minwhz< -2,1,0), 
         co_occurence = as.numeric(whz < -2 & haz < -2), ever_co=as.numeric(sum(co_occurence)>0),
         anystunt_wast = as.numeric(anystunt==1 & anywast==1),
         anystunt_wast_noco = as.numeric(anystunt_wast==1 & ever_co==0),
         stunt = as.numeric(haz < -2),
         wast = as.numeric(whz < -2),
         co = as.numeric(whz < -2 & haz < -2),
         Nobs=n()) %>% 
  #Mark incidence onset of wasting and stunting
  group_by(studyid,country,subjid, stunt) %>%
  mutate(stuntinc = ifelse(stunt==1 & agedays==first(agedays),1,0)) %>%
  group_by(studyid,country,subjid, wast) %>%
  mutate(wastinc = ifelse(wast==1 & agedays==first(agedays),1,0)) %>%
  group_by(studyid,country,subjid, co) %>%
  mutate(co_inc = ifelse(co==1 & agedays==first(agedays),1,0)) %>%
  #Mark if wasting or stunting occurence happened first
  group_by(studyid,country,subjid, stuntinc) %>%
  mutate(minage_stunt_onset=min(agedays)) %>%
  group_by(studyid,country,subjid, wastinc) %>%
  mutate(minage_wast_onset=min(agedays)) %>% 
  group_by(studyid,country,subjid) %>%
  mutate(minage_stunt_onset=ifelse(stuntinc==1, minage_stunt_onset, 99999),
         minage_wast_onset=ifelse(wastinc==1, minage_wast_onset, 99999),
         minage_stunt_onset=min(minage_stunt_onset, na.rm=T),
         minage_wast_onset=min(minage_wast_onset, na.rm=T),
         which_first=ifelse(minage_wast_onset==99999 & minage_stunt_onset==99999,"",
                            ifelse(minage_wast_onset==minage_stunt_onset, "Co-occurence",
                                   ifelse(minage_wast_onset<minage_stunt_onset,"Wasting","Stunting"))),
         which_first=ifelse(anystunt_wast==1,which_first,""),
         wastfirst=as.numeric(which_first=="Wasting" & anystunt_wast_noco==1),
         stuntfirst=as.numeric(which_first=="Stunting" & anystunt_wast_noco==1)
  )%>%
  slice(1) %>%
  mutate(N=n()) %>%
  subset(., select= -c(agedays, haz, whz, measurefreq, measid, co_occurence, stuntinc, wastinc, minage_stunt_onset, minage_wast_onset)) %>%
  ungroup() %>% as.data.frame()
head(as.data.frame(co_ci_0_6))

table(co_ci_0_6$which_first)
table(co_ci_0_6$anystunt_wast)
table(co_ci_0_6$ever_co)


#calculate any co-occurence from 6-24 months

 co_ci_6_24 <- d %>% ungroup() %>% 
    arrange(studyid,country,subjid, agedays) %>% 
    group_by(studyid,country,subjid) %>%
    mutate(laghaz=lag(haz), lagwhz=lag(whz)) %>%
   filter(agecat!="0-6 months") %>%
    #mark if children started stunted and/or wasted. They need to recover to be included in the at-risk pool
    mutate(start_stunt= as.numeric(first(laghaz) < -2), cumsum_notstunted=cumsum(as.numeric(haz >= -2 & laghaz >= -2)), anystuntrecovery=max(cumsum_notstunted)>0) %>%
    mutate(start_wast= as.numeric(first(lagwhz) < -2), cumsum_notwasted=cumsum(as.numeric(whz >= -2 & lagwhz >= -2)), anywastrecovery=max(cumsum_notwasted)>0) %>%
    #drop children never at risk (start stunted and never recovered or started wasted and never recovered) and drop obs prior to recovery
    filter((anystuntrecovery & cumsum_notstunted!=0 & start_stunt==1 & anywastrecovery & cumsum_notwasted!=0 & start_wast==1) |
             (anystuntrecovery & cumsum_notstunted!=0 & start_stunt==1 & start_wast==0) |
             (anywastrecovery & cumsum_notwasted!=0 & start_wast==1 & start_stunt==0) |
             (start_stunt==0 & start_wast==0)) %>%
    subset(., select = -c(start_stunt, cumsum_notstunted, anystuntrecovery, start_wast, cumsum_notwasted, anywastrecovery)) %>%
    mutate(minhaz=min(haz), anystunt=ifelse(minhaz< -2,1,0), 
           minwhz=min(whz), anywast=ifelse(minwhz< -2,1,0), 
           co_occurence = as.numeric(whz < -2 & haz < -2), ever_co=as.numeric(sum(co_occurence)>0),
           anystunt_wast = as.numeric(anystunt==1 & anywast==1),
           anystunt_wast_noco = as.numeric(anystunt_wast==1 & ever_co==0),
           stunt = as.numeric(haz < -2),
           wast = as.numeric(whz < -2),
           co = as.numeric(whz < -2 & haz < -2),
           Nobs=n()) %>% 
    #Mark incidence onset of wasting and stunting
    group_by(studyid,country,subjid, stunt) %>%
    mutate(stuntinc = ifelse(stunt==1 & agedays==first(agedays),1,0)) %>%
    group_by(studyid,country,subjid, wast) %>%
    mutate(wastinc = ifelse(wast==1 & agedays==first(agedays),1,0)) %>%
    group_by(studyid,country,subjid, co) %>%
    mutate(co_inc = ifelse(co==1 & agedays==first(agedays),1,0)) %>%
    #Mark if wasting or stunting occurence happened first
    group_by(studyid,country,subjid, stuntinc) %>%
    mutate(minage_stunt_onset=min(agedays)) %>%
    group_by(studyid,country,subjid, wastinc) %>%
    mutate(minage_wast_onset=min(agedays)) %>% 
    group_by(studyid,country,subjid) %>%
    mutate(minage_stunt_onset=ifelse(stuntinc==1, minage_stunt_onset, 99999),
           minage_wast_onset=ifelse(wastinc==1, minage_wast_onset, 99999),
           minage_stunt_onset=min(minage_stunt_onset, na.rm=T),
           minage_wast_onset=min(minage_wast_onset, na.rm=T),
           which_first=ifelse(minage_wast_onset==99999 & minage_stunt_onset==99999,"",
                              ifelse(minage_wast_onset==minage_stunt_onset, "Co-occurence",
                                     ifelse(minage_wast_onset<minage_stunt_onset,"Wasting","Stunting"))),
           which_first=ifelse(anystunt_wast==1,which_first,""),
           wastfirst=as.numeric(which_first=="Wasting" & anystunt_wast_noco==1),
           stuntfirst=as.numeric(which_first=="Stunting" & anystunt_wast_noco==1)
    ) %>%
    slice(1) %>%
    mutate(N=n()) %>%
    subset(., select= -c(agedays, haz, whz, laghaz, lagwhz, measurefreq, measid,  co_occurence, stuntinc, wastinc, minage_stunt_onset, minage_wast_onset)) %>%
    ungroup() %>% as.data.frame()



#Combine age categories
co_ci <- rbind(co_ci_0_6, co_ci_6_24)

co_ci %>% group_by(studyid, country, agecat) %>% summarize(n(), sum(co_inc), mean(co_inc))

co_ci %>% group_by(agecat) %>% summarize(n(), sum(anystunt_wast), mean(anystunt_wast))


save(co_ci, file="U:/ucb-superlearner/Co-occurrence/co_CI_rf_outcomes.RData")
















