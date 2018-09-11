#-----------------------------------
# Co-occurence analysis
# Objective 1a
# Calculate cumulative incidence (ever stunted) at
# 3, 6, 12, 18, and 24 mo of age

# Cumulative incidence pooled using random effects
# Wasting and stunting can occur at any time in the age range
# They can occur at different times in the age range.
#-----------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(metafor)
theme_set(theme_bw())

# load base functions
source("U:/Scripts/Stunting/2-analyses/0_st_basefunctions.R")

load("U:/Data/Co-occurrence/co-occurrence_data.RData")

# define age windows
d = d %>% 
  #filter(agedays>1) %>%
  mutate(agecat=ifelse(agedays<=3*30.4167,"3 months",
                       ifelse(agedays>3*30.4167 & agedays<=6*30.4167,"6 months",
                              ifelse(agedays>6*30.4167 & agedays<=12*30.4167,"12 months",
                                     ifelse(agedays>12*30.4167 & agedays<=18*30.4167,"18 months",
                                            ifelse(agedays>12*30.4167& agedays<=24*30.4167,"24 months","")))))) %>%
  mutate(agecat=factor(agecat,levels=c("3 months","6 months","12 months","18 months","24 months")))

# check age categories
d %>%
  group_by(agecat) %>%
  summarise(n=sum(!is.na(agedays)),
            min=min(agedays/30.4167),
            mean=mean(agedays/30.4167),
            max=max(agedays/30.4167))

# identify children ever stunted + wasted at the same time
evs = d %>%
  filter(!is.na(agecat)) %>%
  group_by(studyid,country,subjid) %>%
  arrange(studyid,subjid) %>%
  #create variable with minhaz by age category, cumulatively
  mutate(minhaz=ifelse(agecat=="3 months",min(haz[agecat=="3 months"]),
                       ifelse(agecat=="6 months",min(haz[agecat=="3 months" | agecat=="6 months"]),
                              ifelse(agecat=="12 months",min(haz[agecat=="3 months" | agecat=="6 months"|agecat=="12 months"]),
                                     ifelse(agecat=="18 months",min(haz[agecat=="3 months" | agecat=="6 months"|agecat=="12 months"|agecat=="18 months"]),
                                            min(haz)))))) %>%
  mutate(minwhz=ifelse(agecat=="3 months",min(whz[agecat=="3 months"]),
                       ifelse(agecat=="6 months",min(whz[agecat=="3 months" | agecat=="6 months"]),
                              ifelse(agecat=="12 months",min(whz[agecat=="3 months" | agecat=="6 months"|agecat=="12 months"]),
                                     ifelse(agecat=="18 months",min(whz[agecat=="3 months" | agecat=="6 months"|agecat=="12 months"|agecat=="18 months"]),
                                            min(whz)))))) %>%
  # create indicator for whether the child was ever stunted and wasted
  # by age category
  group_by(studyid,country,agecat,subjid) %>%
  summarise(minhaz=min(minhaz), minwhz=min(minwhz)) %>%
  mutate(ever_stunted=ifelse(minhaz< -2,1,0), ever_wasted=ifelse(minwhz< -2,1,0), ever_co=as.numeric(ever_co==1 & ever_wasted==1))



# count incident cases per study by age
# exclude time points if number of measurements per age
# in a study is <50  
cuminc.data= evs%>%
  group_by(studyid,country,agecat) %>%
  summarise(
    nchild=length(unique(subjid)),
    nstudy=length(unique(studyid)),
    ncases=sum(ever_co),
    N=sum(length(ever_co))) %>%
  filter(N>=50)

cuminc.data

# cohort specific results
ci.cohort=lapply(list("3 months","6 months","12 months","18 months","24 months"),function(x) 
  fit.escalc(data=cuminc.data,ni="N", xi="ncases",age=x,meas="PR"))
ci.cohort=as.data.frame(do.call(rbind, ci.cohort))
ci.cohort=cohort.format(ci.cohort,y=ci.cohort$yi,
                        lab=  c("3 months","6 months","12 months","18 months","24 months"))

# estimate random effects, format results
ci.res=lapply(list("3 months","6 months","12 months","18 months","24 months"),function(x)
  fit.rma(data=cuminc.data,ni="N", xi="ncases",age=x,measure="PR",nlab=" measurements"))
ci.res=as.data.frame(do.call(rbind, ci.res))
ci.res[,4]=as.numeric(ci.res[,4])
ci.res = ci.res %>%
  mutate(est=est*100, lb=lb*100, ub=ub*100)
ci.res$agecat.f=as.factor(ifelse(ci.res$agecat=="3 months","0-3 months",
                                 ifelse(ci.res$agecat=="6 months","4-6 months",
                                        ifelse(ci.res$agecat=="12 months","7-12 months",
                                               ifelse(ci.res$agecat=="18 months","13-18 months","19-24 months")))))
ci.res$agecat.f=factor(ci.res$agecat.f,levels=c("0-3 months","4-6 months",
                                                "7-12 months","13-18 months","19-24 months"))
ci.res$ptest.f=sprintf("%0.0f",ci.res$est)

ci.res


# plot cohort incidence

pdf("U:/Figures/stunting-cuminc-africa.pdf",width=11,height=5,onefile=TRUE)
ggplot(ci.cohort[ci.cohort$region=="Africa",],
       aes(y=y,x=age.f))+
  geom_point(size=2)+facet_wrap(~cohort)+
  geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
                 size=2,alpha=0.3) +
  scale_y_continuous(limits=c(0,100))+
  xlab("Age category")+
  ylab("Cumulative incidence (95% CI)")+
  ggtitle("Cohort-specific cumulative incidence of stunting and wasting co-occurence - Africa")
dev.off()

pdf("U:/Figures/co-occurence-cuminc-latamer-eur.pdf",width=11,height=5,onefile=TRUE)
ggplot(ci.cohort[ci.cohort$region=="Latin America"|
                   ci.cohort$region=="Europe",],
       aes(y=y,x=age.f))+
  geom_point(size=2)+facet_wrap(~cohort)+
  geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
                 size=2,alpha=0.3) +
  scale_y_continuous(limits=c(0,100))+
  xlab("Age category")+
  ylab("Cumulative incidence (95% CI)")+
  ggtitle("Cohort-specific cumulative incidence of stunting and wasting co-occurence - Latin America & Europe")
dev.off()

pdf("U:/Figures/co-occurence-cuminc-asia.pdf",width=17,height=7,onefile=TRUE)
ggplot(ci.cohort[ci.cohort$region=="Asia",],
       aes(y=y,x=age.f))+
  geom_point(size=2)+facet_wrap(~cohort)+
  geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
                 size=2,alpha=0.3) +
  scale_y_continuous(limits=c(0,100))+
  xlab("Age category")+
  ylab("Cumulative incidence (95% CI)")+
  ggtitle("Cohort-specific cumulative incidence of stunting and wasting co-occurence - Asia")
dev.off()

# plot pooled cumulative incidence
pdf("U:/Figures/co-occurence-cuminc-pool.pdf",width=9,height=3.5,onefile=TRUE)
ggplot(ci.res,aes(y=est,x=agecat.f))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_y_continuous(limits=c(0,80))+
  xlab("Age category")+
  ylab("Cumulative incidence per 100 children (95% CI)")+
  annotate("text",x=ci.res$agecat.f,y=5,label=ci.res$nmeas.f,size=3)+
  annotate("text",x=ci.res$agecat.f,y=1,label=ci.res$nstudy.f,size=3)+
  annotate("text",label=ci.res$ptest.f,x=ci.res$agecat.f,
           y=ci.res$est,hjust=-0.75,size=3)+
  ggtitle("Pooled cumulative incidence of stunting and wasting co-occurence")
dev.off()

#--------------------------------------
# Sensitivity analysis excluding birth from CI
# only among birth cohort studies
#--------------------------------------
# birth cohorts
bc= d %>% group_by(studyid,country) %>% 
  summarise(minage=min(agedays)) %>%
  # filter(minage==1) %>%
  mutate(birthcohort=ifelse(minage==1,1,0)) %>%
  select(-minage)

# identify ever stunted children
# excluding birth in birth cohort studies
evs.sens.nobirth = left_join(d, bc, by=c("studyid","country")) %>%
  filter(birthcohort==1 & !is.na(agecat) & agedays>1) %>%
  group_by(studyid,country,subjid) %>%
  arrange(studyid,subjid) %>%
  #create variable with minhaz by age category, cumulatively
  mutate(minhaz=ifelse(agecat=="3 months",min(haz[agecat=="3 months"]),
                       ifelse(agecat=="6 months",min(haz[agecat=="3 months" | agecat=="6 months"]),
                              ifelse(agecat=="12 months",min(haz[agecat=="3 months" | agecat=="6 months"|agecat=="12 months"]),
                                     ifelse(agecat=="18 months",min(haz[agecat=="3 months" | agecat=="6 months"|agecat=="12 months"|agecat=="18 months"]),
                                            min(haz)))))) %>%
  # create indicator for whether the child was ever stunted
  # by age category
  group_by(studyid,country,agecat,subjid) %>%
  summarise(minhaz=min(minhaz), minwhz=min(minwhz)) %>%
  mutate(ever_stunted=ifelse(minhaz< -2,1,0), ever_wasted=ifelse(minwhz< -2,1,0), ever_co=as.numeric(ever_co==1 & ever_wasted==1))



# count incident cases per study by age
# exclude time points if number of measurements per age
# in a study is <50  
cuminc.data.nobirth= evs.sens.nobirth%>%
  group_by(studyid,country,agecat) %>%
  summarise(
    nchild=length(unique(subjid)),
    nstudy=length(unique(studyid)),
    ncases=sum(ever_co),
    N=sum(length(ever_co))) %>%
  filter(N>=50)

cuminc.data.nobirth

# estimate random effects, format results
ci.res.nobirth=lapply(list("3 months","6 months","12 months","18 months","24 months"),function(x)
  fit.rma(data=cuminc.data.nobirth,ni="N", xi="ncases",age=x,measure="PR",nlab=" measurements"))
ci.res.nobirth=as.data.frame(do.call(rbind, ci.res.nobirth))
ci.res.nobirth[,4]=as.numeric(ci.res.nobirth[,4])
ci.res.nobirth = ci.res.nobirth %>%
  mutate(est=est*100, lb=lb*100, ub=ub*100)
ci.res.nobirth$agecat.f=as.factor(ifelse(ci.res.nobirth$agecat=="3 months","0-3 months",
                                         ifelse(ci.res.nobirth$agecat=="6 months","4-6 months",
                                                ifelse(ci.res.nobirth$agecat=="12 months","7-12 months",
                                                       ifelse(ci.res.nobirth$agecat=="18 months","13-18 months","19-24 months")))))
ci.res.nobirth$agecat.f=factor(ci.res.nobirth$agecat.f,levels=c("0-3 months","4-6 months",
                                                                "7-12 months","13-18 months","19-24 months"))
ci.res.nobirth$ptest.f=sprintf("%0.0f",ci.res.nobirth$est)

ci.res.nobirth


# plot pooled cumulative incidence
pdf("U:/Figures/co-occurence-cuminc-pool-bc-nobirth.pdf",width=9,height=3.5,onefile=TRUE)
ggplot(ci.res.nobirth,aes(y=est,x=agecat.f))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_y_continuous(limits=c(0,80))+
  xlab("Age category")+
  ylab("Cumulative incidence per 100 children (95% CI)")+
  annotate("text",x=ci.res.nobirth$agecat.f,y=5,label=ci.res.nobirth$nmeas.f,size=3)+
  annotate("text",x=ci.res.nobirth$agecat.f,y=1,label=ci.res.nobirth$nstudy.f,size=3)+
  annotate("text",label=ci.res.nobirth$ptest.f,x=ci.res.nobirth$agecat.f,
           y=ci.res.nobirth$est,hjust=-0.75,size=3)+
  ggtitle("Pooled cumulative incidence of stunting and wasting co-occurence - birth cohorts only - CI excludes birth")
dev.off()


# identify ever stunted children
# including birth in birth cohort studies
evs.sens.birth = left_join(d, bc, by=c("studyid","country")) %>%
  filter(birthcohort==1 & !is.na(agecat)) %>%
  group_by(studyid,country,subjid) %>%
  arrange(studyid,subjid) %>%
  #create variable with minhaz by age category, cumulatively
  mutate(minhaz=ifelse(agecat=="3 months",min(haz[agecat=="3 months"]),
                       ifelse(agecat=="6 months",min(haz[agecat=="3 months" | agecat=="6 months"]),
                              ifelse(agecat=="12 months",min(haz[agecat=="3 months" | agecat=="6 months"|agecat=="12 months"]),
                                     ifelse(agecat=="18 months",min(haz[agecat=="3 months" | agecat=="6 months"|agecat=="12 months"|agecat=="18 months"]),
                                            min(haz)))))) %>%
  # create indicator for whether the child was ever stunted
  # by age category
  group_by(studyid,country,agecat,subjid) %>%
  summarise(minhaz=min(minhaz), minwhz=min(minwhz)) %>%
  mutate(ever_stunted=ifelse(minhaz< -2,1,0), ever_wasted=ifelse(minwhz< -2,1,0), ever_co=as.numeric(ever_co==1 & ever_wasted==1))


# count incident cases per study by age
# exclude time points if number of measurements per age
# in a study is <50  
cuminc.data.birth= evs.sens.birth%>%
  group_by(studyid,country,agecat) %>%
  summarise(
    nchild=length(unique(subjid)),
    nstudy=length(unique(studyid)),
    ncases=sum(ever_co),
    N=sum(length(ever_co))) %>%
  filter(N>=50)

cuminc.data.birth

# estimate random effects, format results
ci.res.birth=lapply(list("3 months","6 months","12 months","18 months","24 months"),function(x)
  fit.rma(data=cuminc.data.birth,ni="N", xi="ncases",age=x,measure="PR",nlab=" measurements"))
ci.res.birth=as.data.frame(do.call(rbind, ci.res.birth))
ci.res.birth[,4]=as.numeric(ci.res.birth[,4])
ci.res.birth = ci.res.birth %>%
  mutate(est=est*100, lb=lb*100, ub=ub*100)
ci.res.birth$agecat.f=as.factor(ifelse(ci.res.birth$agecat=="3 months","0-3 months",
                                       ifelse(ci.res.birth$agecat=="6 months","4-6 months",
                                              ifelse(ci.res.birth$agecat=="12 months","7-12 months",
                                                     ifelse(ci.res.birth$agecat=="18 months","13-18 months","19-24 months")))))
ci.res.birth$agecat.f=factor(ci.res.birth$agecat.f,levels=c("0-3 months","4-6 months",
                                                            "7-12 months","13-18 months","19-24 months"))
ci.res.birth$ptest.f=sprintf("%0.0f",ci.res.birth$est)

ci.res.birth

# plot pooled cumulative incidence
pdf("U:/Figures/co-occurence-cuminc-pool-bc-birth.pdf",width=9,height=3.5,onefile=TRUE)
ggplot(ci.res.birth,aes(y=est,x=agecat.f))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_y_continuous(limits=c(0,80))+
  xlab("Age category")+
  ylab("Cumulative incidence per 100 children (95% CI)")+
  annotate("text",x=ci.res.birth$agecat.f,y=5,label=ci.res.birth$nmeas.f,size=3)+
  annotate("text",x=ci.res.birth$agecat.f,y=1,label=ci.res.birth$nstudy.f,size=3)+
  annotate("text",label=ci.res.birth$ptest.f,x=ci.res.birth$agecat.f,
           y=ci.res.birth$est,hjust=-0.75,size=3)+
  ggtitle("Pooled cumulative incidence of stunting and wasting co-occurence - birth cohorts only - CI includes birth")
dev.off()






#--------------------------------------
# Calculate cumulative incidence of
# co-occurence in certain age ranges only
# including children who had not yet
# become stunted in prior age ranges
#--------------------------------------

#NOTE: NEED TO FIX FOR CO-OCCURANCE OUTCOME

# stunt_ci = d %>% ungroup() %>%
#   #filter(!is.na(agecat) ) %>%
#   filter(!is.na(agecat) & agedays>1) %>%
#   group_by(studyid,country,subjid) %>%
#   arrange(studyid,country,subjid,agedays) %>%
#   mutate(stunt = as.numeric(haz < -2), stunt_nmeas=cumsum(stunt), stunt_onset= as.numeric(stunt==1 & stunt_nmeas==1)) %>%
#   mutate(wast = as.numeric(whz < -2), wast_nmeas=cumsum(wast), wast_onset= as.numeric(wast==1 & wast_nmeas==1)) %>%
#   mutate()
#   filter(stunt_nmeas<2 & wast_nmeas<2) %>%
#   ungroup() %>% group_by(studyid,country, agecat) %>% mutate(N=n()) %>%
#   ungroup() %>% group_by(studyid,country,subjid, agecat) %>% arrange(desc(stunt_onset)) %>% slice(1) %>% 
#   ungroup() 
# 
# 
# group_by(studyid,country,agecat,subjid) %>%
#   summarise(minhaz=min(minhaz), minwhz=min(minwhz)) %>%
#   mutate(ever_stunted=ifelse(minhaz< -2,1,0), ever_wasted=ifelse(minwhz< -2,1,0), ever_co=as.numeric(ever_co==1 & ever_wasted==1))
# 
# 
# 
# # count incident cases per study by age
# # exclude time points if number of measurements per age
# # in a study is <50  
# cuminc.data.agerange <- stunt_ci%>%
#   group_by(studyid,country,agecat) %>%
#   summarise(
#     nchild=length(unique(subjid)),
#     nstudy=length(unique(studyid)),
#     ncases=sum(stunt_onset),
#     N=sum(length(stunt_onset))) %>%
#   filter(N>=50)
# 
# 
# 
# # manually calculate incident cases, person-time at risk at each time point
# stunt_ci %>% group_by(studyid,country,agecat) %>% filter(N>=50) %>%
#   group_by(agecat) %>%
#   summarise(inc.case=sum(stunt_onset),N= n())
# 
# 
# #cuminc.data %>% group_by(agecat) %>% summarise(mean(ncases/nchild))
# cuminc.data.agerange %>% group_by(agecat) %>% summarize(sum(ncases), sum(N), one=mean(ncases/N), two=sum(ncases)/sum(N))
# 
# 
# # estimate random effects, format results
# gc()
# ci.res.agerange=lapply(list("3 months","6 months","12 months","18 months","24 months"),function(x)
#   fit.rma(data=cuminc.data.agerange,ni="N", xi="ncases",age=x,measure="PR",nlab=" at-risk"))
# ci.res.agerange=as.data.frame(do.call(rbind, ci.res.agerange))
# ci.res.agerange[,4]=as.numeric(ci.res.agerange[,4])
# ci.res.agerange = ci.res.agerange %>%
#   mutate(est=est*100, lb=lb*100, ub=ub*100)
# ci.res.agerange$agecat.f=as.factor(ifelse(ci.res.agerange$agecat=="3 months","2 days-3 months",
#                                           ifelse(ci.res.agerange$agecat=="6 months","3-6 months",
#                                                  ifelse(ci.res.agerange$agecat=="12 months","6-12 months",
#                                                         ifelse(ci.res.agerange$agecat=="18 months","12-18 months","18-24 months")))))
# ci.res.agerange$agecat.f=factor(ci.res.agerange$agecat.f,levels=c("2 days-3 months","3-6 months",
#                                                                   "6-12 months","12-18 months","18-24 months"))
# ci.res.agerange$ptest.f=sprintf("%0.0f",ci.res.agerange$est)
# 
# ci.res.agerange
# 
# # plot pooled cumulative incidence
# pdf("U:/Figures/co-occurence-incprop-pool.pdf",width=9,height=3.5,onefile=TRUE)
# ggplot(ci.res.agerange,aes(y=est,x=agecat.f))+
#   geom_point(size=3)+
#   geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
#   scale_y_continuous(limits=c(0,40))+
#   xlab("Age category")+
#   ylab("Incidence proportion (95% CI)")+
#   annotate("text",x=ci.res.agerange$agecat.f,y=5,label=ci.res.agerange$nmeas.f,size=3)+
#   annotate("text",x=ci.res.agerange$agecat.f,y=1,label=ci.res.agerange$nstudy.f,size=3)+
#   annotate("text",label=ci.res.agerange$ptest.f,x=ci.res.agerange$agecat.f,
#            y=ci.res.agerange$est,hjust=-0.75,size=3)+
#   ggtitle("Pooled, age-stratified incidence proportion of stunting and wasting co-occurence")
# dev.off()


# #--------------------------------------
# # cohort specific incidence proportion
# #--------------------------------------
# 
# ip.cohort=lapply(list("3 months","6 months","12 months","18 months","24 months"),function(x) 
#   fit.escalc(data=cuminc.data.agerange,ni="N", xi="ncases",age=x,meas="PR"))
# ip.cohort=as.data.frame(do.call(rbind, ip.cohort))
# ip.cohort=cohort.format(ip.cohort,y=ip.cohort$yi,
#                         lab=  c("3 months","6 months","12 months","18 months","24 months"))
# 
# 
# pdf("U:/Figures/co-occurence-incprop-africa.pdf",width=11,height=5,onefile=TRUE)
# ggplot(ip.cohort[ip.cohort$region=="Africa",],
#        aes(y=y,x=age.f))+
#   geom_point(size=2)+facet_wrap(~cohort)+
#   geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
#                  size=2,alpha=0.3) +
#   scale_y_continuous(limits=c(0,100))+
#   xlab("Age category")+
#   ylab("Percent stunted (95% CI)")+
#   ggtitle("Cohort-specific incidence proportion of stunting and wasting co-occurence within age ranges - Africa")
# 
# dev.off()
# 
# pdf("U:/Figures/co-occurence-incprop-latamer-eur.pdf",width=11,height=5,onefile=TRUE)
# ggplot(ip.cohort[ip.cohort$region=="Latin America"|
#                    ip.cohort$region=="Europe",],
#        aes(y=y,x=age.f))+
#   geom_point(size=2)+facet_wrap(~cohort)+
#   geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
#                  size=2,alpha=0.3) +
#   scale_y_continuous(limits=c(0,100))+
#   xlab("Age category")+
#   ylab("Percent stunted (95% CI)")+
#   ggtitle("Cohort-specific incidence proportion of stunting and wasting co-occurence within age ranges - Latin America & Europe")
# dev.off()
# 
# pdf("U:/Figures/co-occurence-incprop-asia.pdf",width=17,height=7,onefile=TRUE)
# ggplot(ip.cohort[ip.cohort$region=="Asia",],
#        aes(y=y,x=age.f))+
#   geom_point(size=2)+facet_wrap(~cohort)+
#   geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
#                  size=2,alpha=0.3) +
#   scale_y_continuous(limits=c(0,100))+
#   xlab("Age category")+
#   ylab("Percent stunted (95% CI)")+
#   ggtitle("Cohort-specific incidence proportion of stunting and wasting co-occurence within age ranges - Asia")
# dev.off()





#--------------------------------------
# export data 
#--------------------------------------

cuminc=evs %>% select(studyid,country,subjid,agecat,ever_co) 

save(cuminc,file="U:/Data/Co-occurence/st_cuminc.RData")




