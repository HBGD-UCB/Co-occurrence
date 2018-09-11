



#-----------------------------------
# Growth velocity preceding and after
# wasting and stunting
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
d = d %>% filter(agedays <= 24*30.4167) %>%
  mutate(agecat=
           ifelse(agedays>0 & agedays<=6*30.4167,"6 months",
                  ifelse(agedays>6*30.4167 & agedays<=12*30.4167,"12 months",
                         ifelse(agedays>12*30.4167 & agedays<=18*30.4167,"18 months",
                                ifelse(agedays>12*30.4167& agedays<=24*30.4167,"24 months",""))))) %>%
  mutate(agecat=factor(agecat,levels=c("6 months","12 months","18 months","24 months")))



#TEMP
#d <- d %>% filter(studyid=="ki1066203-TanzaniaChild2")



#--------------------------------------------------------------
# identify children ever stunted + wasted at the same time -
# Within specefic age ranges
#--------------------------------------------------------------



#calculate the time of incidence of wasting and stunting from birth to 6 months
inc_0_6 = d %>% ungroup() %>% 
  arrange(studyid,country,subjid, agedays) %>% 
  filter(agecat=="6 months") %>%
  group_by(studyid,country,subjid) %>%
  arrange(studyid,country,subjid, agedays) %>% 
  mutate(agecat="0-6 months", 
         stunt = as.numeric(haz < -2),
         wast = as.numeric(whz < -2),
         co = as.numeric(whz < -2 & haz < -2)) %>%
  group_by(studyid,country,subjid, stunt) %>%
        mutate(stuntinc = ifelse(stunt==1 & agedays==first(agedays),1,0)) %>%
  group_by(studyid,country,subjid, wast) %>%
        mutate(wastinc = ifelse(wast==1 & agedays==first(agedays),1,0)) %>%
  group_by(studyid,country,subjid, co) %>%
        mutate(co_inc = ifelse(co==1 & agedays==first(agedays),1,0)) %>%
  subset(., select= -c(measurefreq, measid)) %>%
  ungroup() %>% as.data.frame()

inc_0_6 %>% group_by(subjid) %>% summarize(sum(stuntinc))


#Function for CI within an age range, dropping kids who start stunted/wasted and who do not recover
inc_onset <- function(d, Agecat, RangeStart){
  
  df = d %>% ungroup() %>% 
    arrange(studyid,country,subjid, agedays) %>% 
    filter(agecat==Agecat) %>%
    group_by(studyid,country,subjid) %>%
    #mark if children started stunted and/or wasted. They need to recover to be included in the at-risk pool
    mutate(start_stunt= as.numeric(first(haz) < -2), cumsum_notstunted=cumsum(as.numeric(haz >= -2)), anystuntrecovery=max(cumsum_notstunted)>0) %>%
    mutate(start_wast= as.numeric(first(whz) < -2), cumsum_notwasted=cumsum(as.numeric(whz >= -2)), anywastrecovery=max(cumsum_notwasted)>0) %>%
    #drop children never at risk (start stunted and never recovered or started wasted and never recovered) and drop obs prior to recovery
    filter((anystuntrecovery & cumsum_notstunted!=0 & start_stunt==1 & anywastrecovery & cumsum_notwasted!=0 & start_wast==1) | 
             (anystuntrecovery & cumsum_notstunted!=0 & start_stunt==1 & start_wast==0) | 
             (anywastrecovery & cumsum_notwasted!=0 & start_wast==1 & start_stunt==0) | 
             (start_stunt==0 & start_wast==0)) %>% 
    subset(., select = -c(start_stunt, cumsum_notstunted, anystuntrecovery, start_wast, cumsum_notwasted, anywastrecovery)) %>%
    arrange(studyid,country,subjid, agedays) %>% 
    mutate(agecat=paste0(RangeStart,"-",Agecat), 
           stunt = as.numeric(haz < -2),
           wast = as.numeric(whz < -2),
           co = as.numeric(whz < -2 & haz < -2)) %>%
      group_by(studyid,country,subjid, stunt) %>%
      mutate(stuntinc = ifelse(stunt==1 & agedays==first(agedays),1,0)) %>%
      group_by(studyid,country,subjid, wast) %>%
      mutate(wastinc = ifelse(wast==1 & agedays==first(agedays),1,0)) %>%
      group_by(studyid,country,subjid, co) %>%
      mutate(co_inc = ifelse(co==1 & agedays==first(agedays),1,0)) %>%
    subset(., select= -c(measurefreq, measid)) %>%
    ungroup() %>% as.data.frame()
  
  return(df)
}

#calculate any stunting from 6-12
inc_6_12 <- inc_onset(d, Agecat="12 months", RangeStart="6")
inc_12_18 <- inc_onset(d, Agecat="18 months", RangeStart="12")
inc_18_24 <- inc_onset(d, Agecat="24 months", RangeStart="18")

#Combine age categories
incidence <- rbind(inc_0_6, inc_6_12, inc_12_18, inc_18_24)


#Calculate change in Z-score before and after incidence
df <- incidence %>% group_by(studyid, subjid, country, agecat) %>%
  mutate(N=n()) %>% #filter(N>4) %>%
  arrange(agedays) %>%
  mutate(delta_haz = ifelse(sum(wastinc)==0, last(haz) - first(haz), NA),
         delta_whz = ifelse(sum(stuntinc)==0, last(whz) - first(whz), NA),
         stunt_episode=cumsum(stuntinc), wast_episode=cumsum(wastinc)) %>%
  group_by(studyid, subjid, country, agecat, stunt_episode) %>% arrange(agedays) %>%
  mutate(pre_delta_whz = ifelse(stunt_episode==1, last(whz) - first(whz), NA),
         post_delta_whz = ifelse(stunt_episode==0, last(whz) - first(whz), NA),
         N=n(),
         pre_delta_whz = ifelse(N<2, NA, pre_delta_whz),
         post_delta_whz = ifelse(N<2, NA, post_delta_whz)) %>%
  group_by(studyid, subjid, country, agecat, wast_episode) %>% arrange(agedays) %>%
  mutate(pre_delta_haz = ifelse(wast_episode==1, last(haz) - first(haz), NA),
         post_delta_haz = ifelse(wast_episode==0, last(haz) - first(haz), NA),
         N=n(),
         pre_delta_haz = ifelse(N<2, NA, pre_delta_haz),
         post_delta_haz = ifelse(N<2, NA, post_delta_haz)) %>%
  arrange(studyid, country, subjid) %>%
  as.data.frame()

head(df)

df2 <- df %>% filter(!is.na(pre_delta_whz) & !is.na(post_delta_whz)) 



resdf <- df %>% 
  group_by(studyid, country, agecat) %>%
  summarize(var_delta_haz=(var(delta_haz, na.rm=T)), 
            var_delta_whz=(var(delta_whz, na.rm=T)), 
            var_pre_delta_haz=(var(pre_delta_haz, na.rm=T)), 
            var_pre_delta_whz=(var(pre_delta_whz, na.rm=T)), 
            var_post_delta_haz=(var(post_delta_haz, na.rm=T)), 
            var_post_delta_whz=(var(post_delta_whz, na.rm=T)), 
            delta_haz=mean(delta_haz, na.rm=T), 
            delta_whz=mean(delta_whz, na.rm=T), 
            pre_delta_haz=mean(pre_delta_haz, na.rm=T), 
            pre_delta_whz=mean(pre_delta_whz, na.rm=T), 
            post_delta_haz=mean(post_delta_haz, na.rm=T), 
            post_delta_whz=mean(post_delta_whz, na.rm=T), 
            n=n()) %>% as.data.frame()
head(resdf)




# random effects function, save results nicely
fit.cont.rma=function(data,age,yi,vi,ni,nlab){
  data=filter(data,agecat==age)
  
  data <- data[data[[vi]]!=0,]
  fit <- NULL
  try(fit <- rma(yi=data[[yi]], vi=data[[vi]], method="REML"))
  if(is.null(fit)) try(fit <- rma(yi=data[[yi]], vi=data[[vi]], method="ML"))
  if(is.null(fit)) try(fit <- rma(yi=data[[yi]], vi=data[[vi]], method="HE"))
  
  out = data %>%
    ungroup() %>%
    summarise(nstudies=length(unique(studyid)),
              nmeas=sum(data[[ni]][agecat==age])) %>%
    mutate(agecat=age,est=fit$beta, se=fit$se, lb=fit$ci.lb, ub=fit$ci.ub,
           outcome= yi,
           nmeas.f=paste0("N=",format(sum(data[[ni]]),big.mark=",",scientific=FALSE),
                          " ",nlab),
           nstudy.f=paste0("N=",nstudies," studies"))
  return(out)
}



ci.res1=lapply(as.list(unique(resdf$agecat)),function(x) 
  fit.cont.rma(data=resdf, ni="n", yi="delta_haz", vi="var_delta_haz", age=x,nlab="children"))      
ci.res2=lapply(as.list(unique(resdf$agecat)),function(x) 
  fit.cont.rma(data=resdf, ni="n", yi="delta_whz", vi="var_delta_whz", age=x,nlab="children"))    
ci.res3=lapply(as.list(unique(resdf$agecat)),function(x) 
  fit.cont.rma(data=resdf, ni="n", yi="pre_delta_haz", vi="var_pre_delta_haz",age=x,nlab="children")) 
ci.res4=lapply(as.list(unique(resdf$agecat)),function(x) 
  fit.cont.rma(data=resdf, ni="n", yi="pre_delta_whz", vi="var_pre_delta_whz",age=x,nlab="children")) 
ci.res5=lapply(as.list(unique(resdf$agecat)),function(x) 
  fit.cont.rma(data=resdf, ni="n", yi="post_delta_haz", vi="var_post_delta_haz",age=x,nlab="children"))       
ci.res6=lapply(as.list(unique(resdf$agecat)),function(x) 
  fit.cont.rma(data=resdf, ni="n", yi="post_delta_whz", vi="var_post_delta_whz",age=x,nlab="children")) 


ci.res <- as.data.frame(do.call(rbind, c(ci.res1, ci.res2, ci.res3, ci.res4, ci.res5, ci.res6)))


unique(ci.res$outcome)
ci.res$OutcomeName <- ci.res %>% select(outcome) %>%
  mutate(
    OutcomeName = case_when(
      outcome=="delta_haz" ~ "Change in HAZ\nwithout wasting",
      outcome=="delta_whz" ~ "Change in WHZ\nwithout stunting",
      outcome=="pre_delta_haz" ~ "Change in HAZ\nafter wasting onset",
      outcome=="pre_delta_whz" ~ "Change in WHZ\nafter stunting onset",
      outcome=="post_delta_haz" ~ "Change in HAZ\nbefore wasting onset",
      outcome=="post_delta_whz" ~ "Change in WHZ\nbefore stunting onset")) %>% 
  select(OutcomeName)
ci.res$OutcomeName <- as.character(ci.res$OutcomeName[,1])


save(ci.res, file="U:/Data/Co-occurrence/pooled_vel_res.RData")


pdf("U:/Figures/co-occurance-vel-pool.pdf",width=10,height=4,onefile=TRUE)
ggplot(ci.res,aes(y=est,x=agecat))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_color_manual(values=tableau10)+xlab("Age category")+
  ylab("Change in Z-score (95% CI)")+
  #scale_y_continuous(limits=c(-4,100))+
  # annotate("text",x=ci.res$agecat,y=0,label=ci.res$nmeas.f,size=3)+
  # annotate("text",x=ci.res$agecat,y=-3,label=ci.res$nstudy.f,size=3)+
  # annotate("text",label=ci.res$ptest.f,x=ci.res$agecat,
  #          y=ci.res$est,hjust=-0.75,size=3)+
  facet_wrap(~OutcomeName, scales="free_y") +
  ggtitle("Pooled change in Z-scores before and after wasting and stunting onset -  age stratified")
dev.off()






#No age stratification


#calculate the time of incidence of wasting and stunting from birth to 6 months
inc_0_24 = d %>% ungroup() %>% 
  arrange(studyid,country,subjid, agedays) %>% 
  group_by(studyid,country,subjid) %>%
  mutate(agecat="0-24 months", 
         stunt = as.numeric(haz < -2),
         wast = as.numeric(whz < -2),
         co = as.numeric(whz < -2 & haz < -2)) %>%
  group_by(studyid,country,subjid, stunt) %>%
  mutate(stuntinc = ifelse(stunt==1 & agedays==first(agedays),1,0)) %>%
  group_by(studyid,country,subjid, wast) %>%
  mutate(wastinc = ifelse(wast==1 & agedays==first(agedays),1,0)) %>%
  group_by(studyid,country,subjid, co) %>%
  mutate(co_inc = ifelse(co==1 & agedays==first(agedays),1,0)) %>%
  subset(., select= -c(measurefreq, measid)) %>%
  ungroup() %>% as.data.frame()



#Calculate change in Z-score before and after incidence
df0_24 <- inc_0_24 %>% group_by(studyid, subjid, country, agecat) %>%
  mutate(N=n()) %>% filter(N>4) %>%
  arrange(agedays) %>%
  mutate(delta_haz = ifelse(sum(wastinc)==0, last(haz) - first(haz), NA),
         delta_whz = ifelse(sum(stuntinc)==0, last(whz) - first(whz), NA),
         stunt_episode=cumsum(stuntinc), wast_episode=cumsum(wastinc)) %>%
  group_by(studyid, subjid, country, agecat, stunt_episode) %>% arrange(agedays) %>%
  mutate(pre_delta_whz = ifelse(stunt_episode==1, last(whz) - first(whz), NA),
         post_delta_whz = ifelse(stunt_episode==0, last(whz) - first(whz), NA),
         N=n(),
         pre_delta_whz = ifelse(N<2, NA, pre_delta_whz),
         post_delta_whz = ifelse(N<2, NA, post_delta_whz)) %>%
  group_by(studyid, subjid, country, agecat, wast_episode) %>% arrange(agedays) %>%
  mutate(pre_delta_haz = ifelse(wast_episode==1, last(haz) - first(haz), NA),
         post_delta_haz = ifelse(wast_episode==0, last(haz) - first(haz), NA),
         N=n(),
         pre_delta_haz = ifelse(N<2, NA, pre_delta_haz),
         post_delta_haz = ifelse(N<2, NA, post_delta_haz)) %>%
  arrange(studyid, country, subjid) %>%
  as.data.frame()

head(df0_24)

#Grab one observation per child
resdf_024 <- df0_24 %>% 
  group_by(studyid, country, subjid) %>%
  summarize(var_delta_haz=(var(delta_haz, na.rm=T)), 
            var_delta_whz=(var(delta_whz, na.rm=T)), 
            var_pre_delta_haz=(var(pre_delta_haz, na.rm=T)), 
            var_pre_delta_whz=(var(pre_delta_whz, na.rm=T)), 
            var_post_delta_haz=(var(post_delta_haz, na.rm=T)), 
            var_post_delta_whz=(var(post_delta_whz, na.rm=T)), 
            delta_haz=mean(delta_haz, na.rm=T), 
            delta_whz=mean(delta_whz, na.rm=T), 
            pre_delta_haz=mean(pre_delta_haz, na.rm=T), 
            pre_delta_whz=mean(pre_delta_whz, na.rm=T), 
            post_delta_haz=mean(post_delta_haz, na.rm=T), 
            post_delta_whz=mean(post_delta_whz, na.rm=T), 
            n=n()) %>% as.data.frame()
head(resdf_024)

resdf_024$agecat <- "0-24 months"


df2 <- resdf_024 %>% filter(!is.na(pre_delta_whz) & !is.na(post_delta_whz) | !is.na(delta_whz))

summary(df$delta_whz)
summary(df$pre_delta_whz)
summary(df$post_delta_whz)

#Get cohort-specific estimates

#Get pooled estimates

ci.res1=lapply(as.list(unique(resdf_024$agecat)),function(x) 
  fit.cont.rma(data=resdf_024, ni="n", yi="delta_haz", vi="var_delta_haz", age=x,nlab="children"))      
ci.res2=lapply(as.list(unique(resdf_024$agecat)),function(x) 
  fit.cont.rma(data=resdf_024, ni="n", yi="delta_whz", vi="var_delta_whz", age=x,nlab="children"))    
ci.res3=lapply(as.list(unique(resdf_024$agecat)),function(x) 
  fit.cont.rma(data=resdf_024, ni="n", yi="pre_delta_haz", vi="var_pre_delta_haz",age=x,nlab="children")) 
ci.res4=lapply(as.list(unique(resdf_024$agecat)),function(x) 
  fit.cont.rma(data=resdf_024, ni="n", yi="pre_delta_whz", vi="var_pre_delta_whz",age=x,nlab="children")) 
ci.res5=lapply(as.list(unique(resdf_024$agecat)),function(x) 
  fit.cont.rma(data=resdf_024, ni="n", yi="post_delta_haz", vi="var_post_delta_haz",age=x,nlab="children"))       
ci.res6=lapply(as.list(unique(resdf_024$agecat)),function(x) 
  fit.cont.rma(data=resdf_024, ni="n", yi="post_delta_whz", vi="var_post_delta_whz",age=x,nlab="children")) 


ci.res0_24 <- as.data.frame(do.call(rbind, c(ci.res1, ci.res2, ci.res3, ci.res4, ci.res5, ci.res6)))



pdf("U:/Figures/co-occurance-vel-pool_0_24.pdf",width=10,height=4,onefile=TRUE)
ggplot(ci.res0_24,aes(y=est,x=outcome))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_color_manual(values=tableau10)+xlab("Age category")+
  ylab("Change in Z-score (95% CI)")+
  #scale_y_continuous(limits=c(-4,100))+
  # annotate("text",x=ci.res$agecat,y=0,label=ci.res$nmeas.f,size=3)+
  # annotate("text",x=ci.res$agecat,y=-3,label=ci.res$nstudy.f,size=3)+
  # annotate("text",label=ci.res$ptest.f,x=ci.res$agecat,
  #          y=ci.res$est,hjust=-0.75,size=3)+
  ggtitle("Pooled change in Z-score before and after wasting and stunting onset")
 dev.off()