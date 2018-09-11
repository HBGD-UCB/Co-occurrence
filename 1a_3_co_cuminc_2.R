


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
source("U:/Scripts/Stunting/2-analyses/0_st_basefunctions.R")

load("U:/Data/Co-occurrence/co-occurrence_data.RData")



#Subset analysis to monthly studies
d <- d %>% filter(measurefreq=="monthly")

#TEMP
#d <- d %>% filter(studyid=="ki1066203-TanzaniaChild2")



# define age windows
d = d %>% 
  #filter(agedays>1) %>%
  mutate(agecat=
           ifelse(agedays>0 & agedays<=6*30.4167,"6 months",
                  ifelse(agedays>6*30.4167 & agedays<=12*30.4167,"12 months",
                         ifelse(agedays>12*30.4167 & agedays<=18*30.4167,"18 months",
                                ifelse(agedays>12*30.4167& agedays<=24*30.4167,"24 months",""))))) %>%
  mutate(agecat=factor(agecat,levels=c("6 months","12 months","18 months","24 months")))



#--------------------------------------------------------------
# identify children ever stunted + wasted at the same time -
#--------------------------------------------------------------


#calculate any stunting from 0-6
co_ci_0_6 = d %>% ungroup() %>% 
  arrange(studyid,country,subjid, agedays) %>% 
  filter(agecat=="6 months") %>%
  group_by(studyid,country,subjid) %>%
  mutate(agecat="0-6 months", 
         minhaz=min(haz), anystunt=ifelse(minhaz< -2,1,0), 
         minwhz=min(whz), anywast=ifelse(minwhz< -2,1,0), 
         co_occurence = as.numeric(whz < -2 & haz < -2), ever_co=as.numeric(sum(co_occurence)>0),
         anystunt_wast = as.numeric(anystunt==1 & anywast==1),
         anystunt_wast_noco = as.numeric(anystunt_wast==1 & ever_co==0),
         stunt = as.numeric(haz < -2),
         wast = as.numeric(whz < -2),
         co = as.numeric(whz < -2 & haz < -2),
         Nobs=n()) %>% 
  #Mark incidence onset of wasting and stunting
  group_by(studyid,country,subjid, stunt) %>% arrange(studyid,country,subjid,agedays) %>%
  mutate(stuntinc = ifelse(stunt==1 & agedays==first(agedays),1,0)) %>%
  group_by(studyid,country,subjid, wast) %>% arrange(studyid,country,subjid,agedays) %>%
  mutate(wastinc = ifelse(wast==1 & agedays==first(agedays),1,0)) %>%
  #Mark if wasting or stunting occurence happened first
  group_by(studyid,country,subjid, stuntinc) %>% arrange(studyid,country,subjid,agedays) %>%
  mutate(minage_stunt_onset=min(agedays)) %>%
  group_by(studyid,country,subjid, wastinc) %>% arrange(studyid,country,subjid,agedays) %>%
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
         wastfirst=as.numeric(which_first=="Wasting" & anystunt_wast==1),
         stuntfirst=as.numeric(which_first=="Stunting" & anystunt_wast==1),
         co_first=as.numeric(which_first=="Co-occurence" & anystunt_wast==1)
  )%>%
  slice(1) %>%
  mutate(N=n()) %>%
  subset(., select= -c(agedays, haz, whz, measurefreq, measid, co_occurence, stuntinc, wastinc, minage_stunt_onset, minage_wast_onset)) %>%
  ungroup() %>% as.data.frame()
head(as.data.frame(co_ci_0_6))

table(co_ci_0_6$which_first)
table(co_ci_0_6$anystunt_wast)
table(co_ci_0_6$ever_co)


sum(co_ci_0_6$anystunt_wast)
sum(co_ci_0_6$co_first + co_ci_0_6$stuntfirst + co_ci_0_6$wastfirst)


#Function for CI within an age range, dropping kids who start stunted/wasted and who do not recover
incprop <- function(d, Agecat, RangeStart){
  
  df = d %>% ungroup() %>% 
    arrange(studyid,country,subjid, agedays) %>% 
    mutate(lagHAZ=lag(haz), lagWHZ=lag(whz)) %>%
    filter(agecat==Agecat) %>%
    group_by(studyid,country,subjid) %>%  arrange(studyid,country,subjid, agedays) %>% 
    #mark if children started stunted and/or wasted. They need to recover to be included in the at-risk pool
    mutate(start_stunt= as.numeric(first(lagHAZ) < -2)) %>%
    mutate(start_wast= as.numeric(first(lagWHZ) < -2), cumsum_notwasted=cumsum(as.numeric(whz >= -2 & lag(whz) >= -2)), anywastrecovery=max(cumsum_notwasted)>0) %>%
    #drop children never at risk (start stunted or started wasted and never recovered) and drop obs prior to wasting recovery
    filter((anywastrecovery & start_wast==1 & start_stunt==0) |
             (start_stunt==0 & start_wast==0)) %>%
    subset(., select = -c(start_stunt, start_wast, cumsum_notwasted, anywastrecovery)) %>%
    mutate(agecat=paste0(RangeStart,"-",Agecat), 
           minhaz=min(haz), anystunt=ifelse(minhaz< -2,1,0), 
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
           wastfirst=as.numeric(which_first=="Wasting" & anystunt_wast==1),
           stuntfirst=as.numeric(which_first=="Stunting" & anystunt_wast==1),
           co_first=as.numeric(which_first=="Co-occurence" & anystunt_wast==1)
    ) %>%
    slice(1) %>%
    mutate(N=n()) %>%
    subset(., select= -c(agedays, haz, whz, measurefreq, measid,  co_occurence, stuntinc, wastinc, minage_stunt_onset, minage_wast_onset)) %>%
    ungroup() %>% as.data.frame()
  
  return(df)
}

#calculate any stunting from 6-12
co_ci_6_12 <- incprop(d, Agecat="12 months", RangeStart="6")
co_ci_12_18 <- incprop(d, Agecat="18 months", RangeStart="12")
co_ci_18_24 <- incprop(d, Agecat="24 months", RangeStart="18")

#Combine age categories
co_ci <- bind_rows(co_ci_0_6, co_ci_6_12, co_ci_12_18, co_ci_18_24)

co_ci %>% group_by(agecat) %>% summarize(n(), sum(anystunt_wast), mean(anystunt_wast))

save(co_ci,file="U:/Data/Co-occurrence/co_ci.RData")


table(co_ci$which_first)
table(co_ci$anystunt_wast)

sum(co_ci$anystunt_wast)
sum(co_ci$co_first + co_ci$stuntfirst + co_ci$wastfirst)





#Need to go from long from to summing the N cases and the total N children
res_co_ci <- co_ci %>% group_by(studyid, country, agecat) %>%
  summarize(Ncases_co=sum(anystunt_wast),
            Ncases_co_first=sum(co_first==1),
            Ncases_stunt_first=sum(stuntfirst==1), 
            Ncases_wast_first=sum(wastfirst==1), 
            Nchildren=sum(N)) %>% 
  filter(Nchildren>1) %>% 
  as.data.frame()  

head(res_co_ci)



#Fit RMA functions, saving outcome name 
fit.rma <- function(data,age,ni,xi,measure,nlab){
  data=filter(data,agecat==age)
  if(measure!="IR"){
    fit<-rma(ni=data[[ni]], xi=data[[xi]], 
             method="REML", measure=measure)
  }else{
    fit<-rma(ti=data[[ni]], xi=data[[xi]], 
             method="REML", measure=measure)
  }
  out=data %>%
    ungroup() %>%
    summarise(nstudies=length(unique(studyid)),
              nmeas=sum(data[[ni]][agecat==age])) %>%
    mutate(agecat=age,est=fit$beta, se=fit$se, lb=fit$ci.lb, ub=fit$ci.ub,
           outcome= xi,
           nmeas.f=paste0("N=",format(sum(data[[ni]]),big.mark=",",scientific=FALSE),
                          " ",nlab),
           nstudy.f=paste0("N=",nstudies," studies"))
  return(out)
}         

# estimate random effects, and format results
ci.res1=lapply(as.list(unique(co_ci$agecat)),function(x) 
  fit.rma(data=res_co_ci, ni="Nchildren", xi="Ncases_co",age=x,measure="PR",nlab="children"))      
ci.res2=lapply(as.list(unique(co_ci$agecat)),function(x) 
  fit.rma(data=res_co_ci, ni="Ncases_co", xi="Ncases_stunt_first",age=x,measure="PR",nlab="children"))      
ci.res3=lapply(as.list(unique(co_ci$agecat)),function(x) 
  fit.rma(data=res_co_ci, ni="Ncases_co", xi="Ncases_wast_first",age=x,measure="PR",nlab="children"))      
ci.res4=lapply(as.list(unique(co_ci$agecat)),function(x) 
  fit.rma(data=res_co_ci, ni="Ncases_co", xi="Ncases_co_first",age=x,measure="PR",nlab="children"))      



ci.res <- as.data.frame(do.call(rbind, c(ci.res1, ci.res2, ci.res3, ci.res4)))
#TEMP
ci.res <- as.data.frame(do.call(rbind, c(ci.res1)))


ci.res$est <- ci.res$est*100
ci.res$lb <- ci.res$lb*100
ci.res$ub <- ci.res$ub*100


unique(ci.res$outcome)

ci.res$OutcomeName <- NA
ci.res$OutcomeName[ci.res$outcome=="Ncases_co"] <- "Co-occurent stunting and wasting"
ci.res$OutcomeName[ci.res$outcome=="Ncases_stunt_first"] <- "Stunting onset first"
ci.res$OutcomeName[ci.res$outcome=="Ncases_wast_first"] <- "Wasting onset first"
ci.res$OutcomeName[ci.res$outcome=="Ncases_co_first"] <- "Co-occurent onset"



ci.res$OutcomeName <- factor(ci.res$OutcomeName, levels=c("Co-occurent stunting and wasting", "Stunting onset first", "Wasting onset first", "Co-occurent onset"))
ci.res$agecat <- factor(ci.res$agecat, levels=c("0-6 months", "6-12 months", "12-18 months", "18-24 months"))
ci.res$est <- as.numeric(ci.res$est)
ci.res <- ci.res %>% arrange(OutcomeName, agecat)



save(ci.res, file="U:/Data/Co-occurrence/pooled_CI_res.RData")


pdf("U:/Figures/co-occurance-ci-pool.pdf",width=8,height=8,onefile=TRUE)
ggplot(ci.res[ci.res$outcome=="Ncases_co",],aes(y=est,x=agecat))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_color_manual(values=tableau10)+xlab("Age category")+
  ylab("Cumulative incidence (95% CI)")+
  ggtitle("Pooled cumulative incidence of stunting and wasting co-occurance")
dev.off()

pdf("U:/Figures/co-occurance-ci-pool_which_first.pdf",width=10,height=8,onefile=TRUE)
ggplot(ci.res[ci.res$outcome!="Ncases_co",],aes(y=est,x=agecat))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_color_manual(values=tableau10)+xlab("Age category")+
  ylab("Proportion of cumulative incidence (95% CI)")+
  facet_wrap(~OutcomeName, scales="fixed", ncol=3) +
  ggtitle("Among children who become wasted and stunted in an age range,\n which condition occurs first?")
dev.off()




#Plot stacked bar chart
#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728", 
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

#Strip grant identifier and add country
co_ci$studyid <- gsub("^k.*?-" , "", co_ci$studyid)
co_ci$studyid <- paste0(co_ci$studyid, ", ", paste0(substring(as.character(co_ci$country),1,1), tolower(substring(as.character(co_ci$country),2))))
co_ci$studyid <- gsub("Tanzania, united republic of", "Tanzania", co_ci$studyid)
co_ci$studyid <- gsub("africa", "Africa", co_ci$studyid)

library(scales)




head(co_ci)
table(co_ci$which_first)

plotdf <- co_ci %>% filter(which_first!="")
plotdf$agecat <- factor(plotdf$agecat, levels=c("0-6 months", "6-12 months", "12-18 months", "18-24 months"))
plotdf <- plotdf %>% arrange(agecat)

table(plotdf$agecat, plotdf$which_first)

ggplot(plotdf, aes(agecat, ..count..)) + geom_bar(aes(fill = which_first), position = "fill") +
  scale_fill_manual(values=tableau10) +
  #facet_wrap(~studyid) +
  xlab("Age") + ylab("Proportion") +
  scale_y_continuous(labels = percent_format()) +
  theme(strip.background = element_blank(),
        legend.position="bottom",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        strip.text.y = element_text(angle = 0)) +
  ggtitle("Proportion of stunting-first, wasting first, and co-occurent wasting and stunting\namong children who became both wasted and stunted in the age range")








