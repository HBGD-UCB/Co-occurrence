



#-----------------------------------
# Growth velocity preceding and after
# wasting and stunting
#-----------------------------------

rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(metafor)
library(broom)
library(purrr)
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





#calculate the time of incidence of wasting and stunting from birth to 24 months
stunting_df = d %>% ungroup() %>%  filter(!is.na(agecat)) %>% #Need to update so that we can calculate post-WHZ if the stunting occurs right before 24 months.
  arrange(studyid,country,subjid, agedays) %>% 
  group_by(studyid,country,subjid) %>%
  mutate(measid= row_number(),
         N_measures = max(measid),
         agecat="0-24 months", 
         stunt = as.numeric(haz < -2))

inc_0_24 <- stunting_df %>%
  group_by(studyid,country,subjid, stunt) %>%
  mutate(stuntinc = ifelse(stunt==1 & agedays==first(agedays),1,0))

inc_0_24 <- inc_0_24 %>%
  group_by(studyid,country,subjid) %>% arrange(studyid, country, subjid,agedays) %>%
  mutate(inc_age=agedays * (stuntinc==1),
         inc_measid=measid * (stuntinc==1)) %>%
  group_by(studyid,country,subjid) %>% arrange(studyid, country, subjid,agedays) %>%
  mutate(everstunt= as.numeric(sum(stuntinc)>0),
         stunt_episode=as.numeric(cumsum(stuntinc)>0)) 


#Fill in inc_age in non-controls
inc_0_24$inc_age[inc_0_24$stuntinc==0] <- NA
inc_0_24$inc_measid[inc_0_24$stuntinc==0] <- NA
inc_0_24 <- inc_0_24 %>%
  group_by(studyid,country,subjid) %>% arrange(studyid, country, subjid, -(inc_age)) %>%
  mutate(inc_age = inc_age[1], inc_measid = inc_measid[1])
table(inc_0_24$inc_measid, inc_0_24$stuntinc)


inc_0_24 <- inc_0_24 %>%
  filter((inc_measid >= 4 & inc_measid < N_measures-2) | everstunt==0)%>% # Filter out observations without enouth measurements before and after
  ungroup() %>% as.data.frame()
  
table(inc_0_24$stunt_episode, inc_0_24$everstunt)
table(inc_0_24$measid, inc_0_24$everstunt)
table(inc_0_24$stunt_episode, inc_0_24$inc_measid)


  #Split into stunted and non-stunted children

df_stunt <- inc_0_24 %>% filter(everstunt==1) %>% group_by(studyid,country) %>% mutate(matchid =  as.numeric(factor(subjid)))

df_nostunt <- inc_0_24 %>% filter(everstunt==0) %>% group_by(studyid,country) %>% arrange(agedays) %>% 
  mutate(matchid = as.numeric(factor(subjid)))  

#Create data frames before and after onset of stunting
pre_stunt <- df_stunt %>% filter(stunt_episode==0) %>% group_by(studyid,country,subjid) %>% mutate(stuntdf=1, max_age=max(agedays), N=n()) %>% filter(N>2)
post_stunt <- df_stunt %>% #filter(stunt_episode==1 & stuntinc==0 ) %>% #Drop this code so that the incident episode is kept in the post
  filter(stunt_episode==1) %>% group_by(studyid,country,subjid) %>%  arrange(studyid, country, subjid, agedays) %>%
  mutate(stuntdf=1, min_age=min(agedays), N=n()) %>% filter(N>2) 

pre_stunt_id <- pre_stunt %>% subset(., select=c(studyid, country, matchid, stuntdf, max_age)) %>% 
                             group_by(studyid, country, matchid) %>% slice(1)
post_stunt_id <- post_stunt %>% subset(., select=c(studyid, country, matchid, stuntdf, min_age)) %>% 
                             group_by(studyid, country, matchid) %>% slice(1)

pre_nostunt <- left_join(df_nostunt, pre_stunt_id, by=c("studyid", "country", "matchid")) %>%
                filter(stuntdf==1 & agedays <= max_age) %>%
                group_by(studyid, country, matchid) %>% arrange(studyid, country, matchid, agedays) %>%
                filter(agedays==max(agedays) | agedays==nth(agedays,-2) | agedays==nth(agedays,-3))
post_nostunt <- left_join(df_nostunt, post_stunt_id, by=c("studyid", "country", "matchid")) %>%
                filter(stuntdf==1 & agedays >= min_age)  %>%
                group_by(studyid, country, matchid) %>% arrange(studyid, country, matchid, agedays) %>%
                filter(agedays==min(agedays) | agedays==nth(agedays,2) | agedays==nth(agedays,3))

# ggplot(pre_nostunt, aes(x=agedays, y=whz)) + geom_smooth()
# ggplot(post_nostunt, aes(x=agedays, y=whz)) + geom_smooth()

dim(pre_stunt)
dim(post_stunt)
dim(pre_nostunt)
dim(post_nostunt)


pre_stunt$group <- "Pre stunting onset"
post_stunt$group <- "Post stunting onset"
pre_nostunt$group <- "Pre stunting onset - matched control"
post_nostunt$group <- "Post stunting onset - matched control"

pre_stunt$pre <- 1
post_stunt$pre <- 0
pre_nostunt$pre <- 1
post_nostunt$pre <- 0

pre_stunt$Control <- "Stunted"
post_stunt$Control <- "Stunted"
pre_nostunt$Control <- "Control"
post_nostunt$Control <- "Control"


nostunt <- bind_rows(pre_nostunt, post_nostunt)

#Fill in control midpoint
nostunt <- nostunt %>% group_by(studyid, country, matchid )  %>% arrange(studyid, country, matchid ,agedays) %>%
  mutate(match_stuntinc= 1* (pre==0 & lag(pre)==1),
         inc_age=agedays * match_stuntinc,
         inc_measid=measid * match_stuntinc)
nostunt$inc_age[nostunt$match_stuntinc==0] <- NA
nostunt$inc_measid[nostunt$match_stuntinc==0] <- NA
nostunt <- nostunt %>%
  group_by(studyid,country,subjid) %>% arrange(studyid, country, subjid, -(inc_age)) %>%
  mutate(inc_age = inc_age[1], inc_measid = inc_measid[1]) %>%
  filter(!is.na(inc_age))


#Create analysis dataset
df <- bind_rows(pre_stunt, post_stunt,nostunt)

df <- df %>% mutate(cohort=paste0(studyid," ",country)) %>% 
        group_by(studyid, country, matchid, group)  %>%
        arrange(agedays)


#Check matching
summary(df$inc_age)
summary(df$inc_age[df$Control=="Control"])
summary(df$inc_age[df$Control=="Stunted"])

summary(df$agedays[df$Control=="Control"])
summary(df$agedays[df$Control=="Stunted"])

dim(df[df$Control=="Control",])
dim(df[df$Control=="Stunted",])

df2 <- df %>% mutate(minage=ifelse(pre==1, 
                                   nth(agedays,-3),
                                   first(agedays)),
                     maxage=ifelse(pre==1, 
                                   last(agedays),
                                   nth(agedays,3))) %>%
              filter(agedays >= minage & agedays <= maxage) %>%
              mutate(agedays2=agedays-min(agedays), 
                     agedays2=ifelse(pre==1, agedays2-max(agedays2) +100, agedays2),
                     agedays3= agedays-inc_age, # Center age at the incident age
                     N=n()) 


df2$group <- factor(df2$group, levels=c("Pre stunting onset","Post stunting onset","Pre stunting onset - matched control", "Post stunting onset - matched control"))
df2 <- df2 %>% arrange(group)

df2 %>% group_by(group) %>% summarize(mean(whz))

df2$Control <- factor(df2$Control, levels=c("Stunted","Control"))

# p <- ggplot(df2, aes(x=agedays2, y=whz, group=cohort)) + geom_smooth(method="lm",alpha=0.1) + facet_wrap(~group, scales="fixed")
# p

p2 <- ggplot(df2[df2$agedays3>=-100 & df2$agedays3<=100,], aes(x=agedays3, y=whz)) + #geom_point(alpha=0.1)+
                  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) + 
                  #geom_smooth(method = 'loess') +
                  #geom_point(alpha=0.1) +
                  #coord_cartesian(xlim=-c(-100,100)) +
                  geom_vline(xintercept=0) +
                  facet_wrap(~Control, scales="fixed") +
                    xlab("Days around stunting incidence (x=0) and around matched control points") + ylab("WHZ") + 
                  ggtitle("WHZ 100 days before and after stunting onset, \nand in matched controls")
p2


# p2 <- ggplot(df2, aes(x=agedays2, y=whz)) + geom_point(alpha=0.1)+
#    geom_line(aes(group=subjid, alpha=0.05)) + geom_smooth() +
#     facet_wrap(~group, scales="free_y") +
#   xlab("Days") + ylab("WHZ") + ggtitle("WHZ 100 days before and after stunting onset, \nand in matched controls")
# p2



#Get mean before and after stunting onset

meanWHZ_df <- df2 %>% group_by(studyid, country, matchid, group) %>% arrange(agedays) %>%
                     #mutate(delta_whz = last(whz) - first(whz)) %>% slice(1)
                     mutate(meanwhz = mean(whz)) %>% slice(1)


#Grab one observation per child
resdf_024 <- meanWHZ_df %>% 
  group_by(studyid, country, group) %>%
  summarize(var_meanwhz=var(meanwhz, na.rm=T), 
            meanwhz=mean(meanwhz, na.rm=T), 
            n=n()) %>% as.data.frame()
head(resdf_024)

resdf_024$agecat <-resdf_024$group



# random effects function, save results nicely
fit.cont.rma=function(data,age,yi,vi,ni,nlab, method="REML"){
  data=filter(data,agecat==age)
  
  data <- data[data[[vi]]!=0,]
  fit <- NULL
  if(method=="FE"){
    fit <- rma(yi=data[[yi]], vi=data[[vi]], method="FE")
  }else{
  try(fit <- rma(yi=data[[yi]], vi=data[[vi]], method="REML"))
  if(is.null(fit)) try(fit <- rma(yi=data[[yi]], vi=data[[vi]], method="ML"))
  if(is.null(fit)) try(fit <- rma(yi=data[[yi]], vi=data[[vi]], method="HE"))
  }
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



#Get pooled estimates

resdf_024 %>% group_by(group) %>% summarize(mean(meanwhz))

resdf_024 <- as.data.frame(resdf_024)

ci.res=lapply(as.list(unique(resdf_024$group)),function(x) 
  fit.cont.rma(data=resdf_024, ni="n", yi="meanwhz", vi="var_meanwhz", age=x, nlab="children"))  



ci.res0_24 <- as.data.frame(do.call(rbind, c(ci.res)))


#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728", 
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

ci.res0_24$agecat <- factor(ci.res0_24$agecat, levels=c("Pre stunting onset", "Pre stunting onset - matched control", "Post stunting onset", "Post stunting onset - matched control"))


ci.res0_24$Control <- c("Stunted", "Stunted", "Control", "Control")
ci.res0_24$Control <- factor(ci.res0_24$Control, levels=unique(ci.res0_24$Control))


pdf("U:/Figures/co-occurance-meanWHZ-pool_0_24.pdf",width=10,height=4,onefile=TRUE)
ggplot(ci.res0_24,aes(y=est,x=agecat, color=Control, fill=Control))+
  scale_fill_manual(values=tableau10) +
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_color_manual(values=tableau10)+xlab("Matched stunted and control children,\nbefore and after stunting onset")+
  ylab("Mean Z-score (95% CI)")+
  ggtitle("Mean Z-score before and after wasting and stunting onset")
dev.off()







# Change in WHZ


# 
# #Get slope of model
# 
# change_df <- df2 %>%
#   group_by(studyid, country, matchid, group) %>% arrange(agedays) %>%
#   nest() %>%
#   mutate(model = map(data, ~lm(whz ~ 1 + agedays, data = .x) %>% #Should I fit models by cohort, with interaction term by group, and extract coefficients?
#                        tidy)) %>%
#   unnest(model) %>%
#   filter(term == 'agedays') %>%
#   subset(., select = -c(std.error, statistic, p.value)) %>%
#   rename(delta_whz=estimate)
# 
# 
# 
# 
# #Grab one observation per child
# resdf_024 <- change_df %>% 
#   group_by(studyid, country, group) %>%
#   summarize(var_delta_whz=var(delta_whz, na.rm=T), 
#             delta_whz=mean(delta_whz, na.rm=T), 
#             n=n()) %>% as.data.frame()
# head(resdf_024)
# 
# resdf_024$agecat <-resdf_024$group
# 
# 
# 
# 
# 
# #Get pooled estimates
# 
# resdf_024 %>% group_by(group) %>% summarize(mean(delta_whz))
# 
# resdf_024 <- as.data.frame(resdf_024)
# 
# ci.res=lapply(as.list(unique(resdf_024$group)),function(x) 
#   fit.cont.rma(data=resdf_024, ni="n", yi="delta_whz", vi="var_delta_whz", age=x,nlab="children"))  
# 
# ci.res_FE=lapply(as.list(unique(resdf_024$group)),function(x) 
#   fit.cont.rma(data=resdf_024, ni="n", yi="delta_whz", vi="var_delta_whz", age=x,nlab="children", method="FE"))  
# 
# 
# 
# ci.res0_24 <- as.data.frame(do.call(rbind, c(ci.res)))
# ci.res0_24_FE <- as.data.frame(do.call(rbind, c(ci.res_FE)))
# 
# #hbgdki pallet
# tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728", 
#                "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
# 
# ci.res0_24$agecat <- factor(ci.res0_24$agecat, levels=c("Pre stunting onset", "Pre stunting onset - matched control", "Post stunting onset", "Post stunting onset - matched control"))
# 
# 
# ci.res0_24$Control <- c("Stunted", "Stunted", "Control", "Control")
# ci.res0_24$Control <- factor(ci.res0_24$Control, levels=unique(ci.res0_24$Control))
# 
# 
# pdf("U:/Figures/co-occurance-vel-pool_0_24.pdf",width=10,height=4,onefile=TRUE)
# ggplot(ci.res0_24,aes(y=est,x=agecat, color=Control, fill=Control))+
#   scale_fill_manual(values=tableau10) +
#   geom_point(size=3)+
#   geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
#   scale_color_manual(values=tableau10)+xlab("Matched stunted and control children,\nbefore and after stunting onset")+
#   ylab("Change in Z-score (95% CI)")+
#   ggtitle("Pooled change in Z-score before and after wasting and stunting onset")
# dev.off()
