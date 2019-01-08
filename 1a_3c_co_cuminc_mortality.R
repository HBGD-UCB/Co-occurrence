



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

load("U:/Data/Co-occurrence/int_co-occurrence_data.RData")




#--------------------------------------------------------------
# identify children ever stunted + wasted at the same time -
# Cumulative incidence from birth to 6 months
#--------------------------------------------------------------


CI_06 <- d %>%
  group_by(studyid,country,subjid) %>%
  arrange(studyid,subjid) %>%
  filter(agedays<=6*30.4167) %>%
  group_by(studyid,country,subjid) %>%
  mutate(co_occurence = as.numeric(whz < (-2) & haz < (-2)), 
         ever_co06=as.numeric(sum(co_occurence)>0)) %>% #Co-occurence is stunting and wasting simultaneously
  slice(1) %>% subset(., select = c(studyid, subjid, ever_co06))
table(CI_06$ever_co06)

#--------------------------------------------------------------
# identify children ever stunted + wasted at the same time -
# Cumulative incidence from birth to 24 months
#--------------------------------------------------------------


CI_024 <- d %>%
  group_by(studyid,country,subjid) %>%
  arrange(studyid,subjid) %>%
  filter(agedays<=24*30.4167) %>%
  group_by(studyid,country,subjid) %>%
    mutate(co_occurence = as.numeric(whz < (-2) & haz < (-2)), 
           ever_co024=as.numeric(sum(co_occurence)>0)) %>% #Co-occurence is stunting and wasting simultaneously
  slice(1) %>% subset(., select = c(studyid, subjid, ever_co024))
table(CI_024$ever_co024)
  
#--------------------------------------------------------------
# identify children ever stunted + wasted at the same time -
# Prevalence at 18 months
#--------------------------------------------------------------


#  get observations closest to 18 months and calculate co-occurence
co_prev18 <- d %>% 
  arrange(studyid,subjid,agedays) %>%
  filter(agedays>17*30.4167 & agedays<19*30.4167) %>%
  group_by(studyid,country,subjid) %>%
  filter(!is.na(haz) & !is.na(whz)) %>%
  filter(abs(agedays-18*30.4167)==min(abs(agedays-18*30.4167))) %>%
  mutate(co_occurence = as.numeric(whz < (-2) & haz < (-2)))
table(co_prev18$co_occurence)





#-----------------------------------
# save outcome datasets
#-----------------------------------
save(CI_06, CI_024, co_prev18, file="U:/ucb-superlearner/Co-occurrence/co_prev_mortality.RData")














