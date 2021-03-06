


#-----------------------------------
# Co-occurence analysis
# Additional risk factor analysis
# Quartile of WHZ preceeding stunting
#-----------------------------------

rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(metafor)
theme_set(theme_bw())

# load base functions
source("U:/Scripts/Stunting/2-analyses/0_st_basefunctions.R")

load("U:/Data/Co-occurrence/co-occurrence_data.RData")



#Subset analysis to monthly studies
d <- d %>% filter(measurefreq=="monthly")


# define age windows
d <- d %>% filter(agedays <= 25*30.4167)

#Cut into 3 month quartiles
d <- d %>% mutate(agemonth = agedays/30.4167)

d$agecat <- cut(d$agemonth, 
                breaks=c(0,3,6,9,12,15,18,21,24),
                labels=c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "12-15 months", "15-18 months", "18-21 months", "21-24 months"),
                include.lowest=T, right=T)
d$agecat <-factor(d$agecat, labels=c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "12-15 months", "15-18 months", "18-21 months", "21-24 months"))
table(d$agecat)







#Quartile WHZ and HAZ
d_whzcat <- d %>% group_by(subjid, agecat) %>% mutate(meanWHZ = mean(whz)) %>% slice(1) %>% arrange(subjid,agecat, agedays) %>%
  group_by(agecat) %>% mutate(agelevel=as.numeric(agecat),WHZ_quart = (ntile(meanWHZ, 4))) %>% 
  subset(., select = -c(haz, whz, agedays, measurefreq)) %>% ungroup() %>%
  mutate(lag_agelevel= lag(agelevel), lag_WHZ_quart=lag(WHZ_quart)) %>%
  filter(!is.na(lag_WHZ_quart)) %>%
  filter(agelevel==lag_agelevel+1) %>%
  mutate(lag_WHZ_quart=factor(lag_WHZ_quart))


table(d_whzcat$agecat, d_whzcat$WHZ_quart)


d_hazcat <- d %>% group_by(subjid, agecat) %>% mutate(meanHAZ = mean(haz)) %>% slice(1) %>% arrange(subjid,agecat, agedays) %>%
  group_by(agecat) %>% mutate(agelevel=as.numeric(agecat),HAZ_quart = (ntile(meanHAZ, 4))) %>% 
  subset(., select = -c(haz, whz, agedays, measurefreq)) %>% ungroup() %>%
  mutate(lag_agelevel= lag(agelevel), lag_HAZ_quart=lag(HAZ_quart)) %>%
  filter(!is.na(lag_HAZ_quart)) %>%
  filter(agelevel==lag_agelevel+1) %>%
  mutate(lag_HAZ_quart=factor(lag_HAZ_quart))


#mark measure frequencies
# d$measurefreq <- NA
# 
# d$measurefreq[d$studyid %in% c(
#   "ki0047075b-MAL-ED",   
#   "ki1000108-CMC-V-BCS-2002",              
#   "ki1000108-IRC",               
#   "ki1000109-EE",           
#   "ki1000109-ResPak",  
#   "ki1017093b-PROVIDE",  
#   "ki1066203-TanzaniaChild2",           
#   "ki1101329-Keneba",  
#   "ki1112895-Guatemala BSC",       
#   "ki1113344-GMS-Nepal",             
#   "ki1114097-CONTENT"
# )] <- "monthly"



df <- left_join(d, d_whzcat, by=c("subjid", "agecat"))
df <- df %>% filter(!is.na(lag_WHZ_quart))
df <- droplevels(df)

df_HAZ <- left_join(d, d_hazcat, by=c("subjid", "agecat"))
df_HAZ <- df_HAZ %>% filter(!is.na(lag_HAZ_quart))
df_HAZ <- droplevels(df_HAZ)


#Splines of HAZ by prior WHZ category
#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728", 
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

plotdf <- df %>% filter(measurefreq=="monthly")
plotdf_HAZ <- df_HAZ %>% filter(measurefreq=="monthly")

p<-ggplot(plotdf, aes(x=agedays, y=haz, group=lag_WHZ_quart, color=lag_WHZ_quart)) + geom_smooth(method = 'loess') +
  facet_wrap(~agecat, scales="free_x", nrow=1) +
  scale_color_manual(values=tableau10, name = "Quartile of WHZ in\nthe prior 3 months")+
  xlab("Child age in days") + ylab("HAZ") + 
  ggtitle("Spline curves of HAZ over 3-month age ranges\nstratified by quartile of WHZ in prior 3-month range.") +
  theme(strip.background = element_blank(),
        legend.position="bottom",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 0, hjust = 1)) 

save(plotdf, plotdf_HAZ, file="U:/Data/Stunting/HAZ_by_WHZ.RData")
ggsave(p, file="U:/Figures/Stunting Webinar/HAZcurves_by_WHZ.png", width=5.7, height=4.6)



#Calculate CI across 3 month age ranges.

#calculate any stunting from 3-24
stunt_ci <- NULL

agecat_vec <- c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "12-15 months", "15-18 months", "18-21 months", "21-24 months")

for(i in 2:length(agecat_vec)){
  current_agecat<-agecat_vec[i]
  lag_agecat <-agecat_vec[i-1]
  
  temp = d %>% ungroup() %>% 
    group_by(studyid,country,subjid) %>%
    arrange(studyid,country,subjid, agedays) %>% 
    mutate(stunt=1*haz< (-2), numstunt=cumsum(stunt), 
           anystunt_prior = 1*(agecat==lag_agecat & numstunt>0)) %>%
    filter(agecat==lag_agecat | agecat==current_agecat) %>%
    mutate(anystunt_prior = anystunt_prior[1]) %>% 
    filter(agecat==current_agecat & !is.na(agecat) & anystunt_prior==0) %>%
    mutate( ever_stunted=ifelse(min(haz)< -2,1,0), Nobs=n()) %>% slice(1) %>%
    mutate(N=n()) %>%
    ungroup() %>%
    select(studyid,subjid, country,tr,agedays,haz, measurefreq, measid, agecat, ever_stunted,Nobs, N)
  
  stunt_ci <- rbind(stunt_ci, temp)
  
}

table(stunt_ci$agecat, stunt_ci$ever_stunted)

table(stunt_ci, stunt_ci$ever_stunted)



#--------------------------------
# Prevalence outcome 
#--------------------------------



# define age windows
dprev = d %>% 
  arrange(studyid,subjid,agedays) %>%
  mutate(agecat=ifelse(agedays==1,"Birth",
                       ifelse(agedays>2*30.4167 & agedays<4*30.4167,"3 months",
                              ifelse(agedays>5*30.4167 & agedays<7*30.4167,"6 months",
                                     ifelse(agedays>8*30.4167 & agedays<10*30.4167,"9 months",
                                            ifelse(agedays>11*30.4167 & agedays<13*30.4167,"12 months",
                                                   ifelse(agedays>14*30.4167 & agedays<16*30.4167,"15 months",
                                                          ifelse(agedays>17*30.4167 & agedays<19*30.4167,"18 months",
                                                                 ifelse(agedays>20*30.4167 & agedays<22*30.4167,"21 months",
                                                          ifelse(agedays>23*30.4167& agedays<25*30.4167,"24 months","")))))))))) %>%
  mutate(agecat=factor(agecat,levels=c("Birth","3 months","6 months","9 months",
                                       "12 months","15 months","18 months","21 months","24 months")),
         agelevel=(as.numeric(agecat)-1)*3) %>%
  filter(!is.na(agecat))

table(dprev$agecat)

#  Get the observation closest to prevalence times
dprev <- dprev %>%
  filter(!is.na(agecat)) %>%
  group_by(studyid,country,subjid,agecat) %>%
  filter(!is.na(haz)) %>%
  filter(abs(agedays-agelevel*30.4167)==min(abs(agedays-agelevel*30.4167))) %>%
  mutate(N=n())



#Calc outcomes
dprev$stunted=ifelse(dprev$haz< (-2), 1,0)

table(dprev$agecat, dprev$stunted)

dprev <- dprev %>% filter(agecat!="Birth")

#--------------------------------
# Merge stunting outcomes and WHZ quartiles.
#--------------------------------

stuntCI_whz <- left_join(stunt_ci, d_whzcat, by=c("studyid","country","subjid", "agecat"))

table(stuntCI_whz$agecat, stuntCI_whz$lag_WHZ_quart, stuntCI_whz$ever_stunted)

table(stuntCI_whz$agecat, stuntCI_whz$studyid, stuntCI_whz$ever_stunted)

table(d_whzcat$agecat)
dprev_whzcat <- d_whzcat 
dprev_whzcat$agecat <- factor(paste0(sapply(strsplit(as.character(dprev_whzcat$agecat), "-", fixed=T), `[`, 2)," months"))
table(dprev_whzcat$agecat)


stuntprev_whz <- left_join(dprev, dprev_whzcat, by=c("studyid","country","subjid", "agecat"))
table(stuntprev_whz$agecat, stuntprev_whz$lag_WHZ_quart, stuntprev_whz$stunted)



#--------------------------------
# Merge stunting datasets and covariates
#--------------------------------

setwd("U:/ucb-superlearner/Stunting rallies/")

#load covariates
cov<-readRDS("FINAL_clean_covariates.rds")

#Merge in covariates
d <- left_join(stuntprev_whz, cov, by=c("studyid", "subjid", "country"))


#Treatment name
A<-c("lag_WHZ_quart")

#Vector of outcome names
Y<-c("stunted")


#Vector of covariate names
W=c("arm","sex", "W_mage", "W_fage", "meducyrs", "feducyrs", "hhwealth_quart", "hfoodsec",
    #"W_mhtcm","W_mwtkg","W_bmi", "W_fhtcm",
    "vagbrth","hdlvry",
    "single",
    "W_nrooms","W_nhh","W_nchldlt5",
    "month","brthmon","W_parity",
    "trth2o","cleanck","impfloor","impsan","safeh20")

#Subgroup variable
V <- c("agecat")

#clusterid ID variable
id <- c("id")

save(d, Y, A, W, V, id, file="stuntprev_whz_rf.Rdata")



#CI outcomes

d <- left_join(stuntCI_whz, cov, by=c("studyid", "subjid", "country"))
Y<-c("ever_stunted")
save(d, Y, A, W, V, id, file="stuntCI_whz_rf.Rdata")