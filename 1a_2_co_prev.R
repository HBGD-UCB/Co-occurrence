#-----------------------------------
# Co-occurence analysis
# Objective 1a
# Calculate point prevalence at
# Birth, 3, 6, 12, 18, and 24 mo of age

# Prevalence pooled using random effects
#-----------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(binom)
library(metafor)
theme_set(theme_bw())

#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728", 
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

# load base functions
source("U:/Scripts/Stunting/2-analyses/0_st_basefunctions.R")

load("U:/Data/Co-occurrence/co-occurrence_data.RData")

# define age windows
d = d %>% 
  arrange(studyid,subjid,agedays) %>%
  mutate(agecat=ifelse(agedays==1,"Birth",
     ifelse(agedays>2*30.4167 & agedays<4*30.4167,"3 months",
      ifelse(agedays>5*30.4167 & agedays<7*30.4167,"6 months",
             ifelse(agedays>8*30.4167 & agedays<10*30.4167,"9 months",
                  ifelse(agedays>11*30.4167 & agedays<13*30.4167,"12 months",
              ifelse(agedays>17*30.4167 & agedays<19*30.4167,"18 months",
                     ifelse(agedays>23*30.4167& agedays<25*30.4167,"24 months","")))))))) %>%
    mutate(agecat=factor(agecat,levels=c("Birth","3 months","6 months","9 months",
                                         "12 months","18 months","24 months")),
           agelevel=(as.numeric(agecat)-1)*3) %>%
    filter(!is.na(agecat))

# check age categories
d %>%
  group_by(agecat) %>%
  summarise(n=sum(!is.na(agedays)),
            min=min(agedays/30.4167),
            mean=mean(agedays/30.4167),
            max=max(agedays/30.4167))

#  Get the observation closest to prevalence times
dmn <- d %>%
  filter(!is.na(agecat)) %>%
  group_by(studyid,country,subjid,agecat) %>%
  filter(!is.na(haz) & !is.na(whz)) %>%
  filter(abs(agedays-agelevel*30.4167)==min(abs(agedays-agelevel*30.4167))) %>%
  mutate(N=n())

mean(dmn$N)


#Calc outcomes

dmn$cowaststunt=ifelse(dmn$haz< -2 & dmn$whz < -2, 1,0)
dmn$onlywast=ifelse(dmn$haz>= -2 & dmn$whz < -2, 1,0)
dmn$onlystunt=ifelse(dmn$haz< -2 & dmn$whz >= -2, 1,0)
dmn$nowaststunt=ifelse(dmn$haz>= -2 & dmn$whz >= -2, 1,0)
dmn$status <- NA
dmn$status[dmn$cowaststunt==1] <- "Co-occurence"
dmn$status[dmn$onlywast==1] <- "Wasting only"
dmn$status[dmn$onlystunt==1] <- "Stunting only"
dmn$status[dmn$nowaststunt==1] <- "Healthy"
dmn$status <- factor(dmn$status, levels=rev(c("Co-occurence", "Wasting only", "Stunting only", "Healthy")))
dmn < dmn %>% arrange(status)

table(dmn$agecat, dmn$status)





# count measurements per study by age
# exclude time points if number of measurements per age
# in a study is <50
prev.data = dmn %>%
  filter(!is.na(agecat)) %>%
  group_by(studyid,country,agecat) %>%
  summarise(nmeas=n(),
            prev=mean(cowaststunt),
            nxprev=sum(cowaststunt==1)) %>%
  filter(nmeas>=50) 
  
# cohort specific results
prev.cohort=lapply(list("Birth","3 months","6 months","9 months","12 months","18 months","24 months"),function(x) 
  fit.escalc(data=prev.data,ni="nmeas", xi="nxprev",age=x,meas="PR"))
prev.cohort=as.data.frame(do.call(rbind, prev.cohort))
prev.cohort=cohort.format(prev.cohort,y=prev.cohort$yi,
                         lab=  c("Birth","3m","6m","9m","12m","18m","24m"))

# estimate random effects, format results
prev.res=lapply(list("Birth","3 months","6 months","9 months","12 months","18 months","24 months"),function(x) 
  fit.rma(data=prev.data,ni="nmeas", xi="nxprev",age=x,measure="PR",nlab="children"))
prev.res=as.data.frame(do.call(rbind, prev.res))
prev.res[,4]=as.numeric(prev.res[,4])
                prev.res = prev.res %>%
  mutate(est=est*100,lb=lb*100,ub=ub*100)
prev.res$agecat=factor(prev.res$agecat,levels=c("Birth","3 months","6 months","9 months","12 months","18 months","24 months"))
prev.res$ptest.f=sprintf("%0.0f",prev.res$est)

# plot cohort prevalence
pdf("U:/Figures/co-occurance-ptprev-africa.pdf",width=11,height=5,onefile=TRUE)
p_africa <- ggplot(prev.cohort[prev.cohort$region=="Africa",],
       aes(y=y,x=age.f))+
  geom_point(size=2)+facet_wrap(~cohort)+
  geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
                 size=2,alpha=0.3) +
  # scale_y_continuous(limits=c(0,96))+
  xlab("Age in months")+
  ylab("Point prevalence (95% CI)")+
  ggtitle("Cohort-specific point prevalence of stunting and wasting co-occurance - Africa")

print(p_africa)
dev.off()

pdf("U:/Figures/co-occurance-ptprev-latamer-eur.pdf",width=8,height=5,onefile=TRUE)
p_latamer <- ggplot(prev.cohort[prev.cohort$region=="Latin America"|
                     prev.cohort$region=="Europe",],
       aes(y=y,x=age.f))+
  geom_point(size=2)+facet_wrap(~cohort)+
  geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
                 size=2,alpha=0.3) +
  # scale_y_continuous(limits=c(0,96))+
  xlab("Age in months")+
  ylab("Point prevalence (95% CI)")+
  ggtitle("Cohort-specific point prevalence of stunting and wasting co-occurance - Latin America & Europe")
dev.off()

pdf("U:/Figures/co-occurance-ptprev-asia.pdf",width=15,height=7,onefile=TRUE)
p_asia <- ggplot(prev.cohort[prev.cohort$region=="Asia",],
       aes(y=y,x=age.f))+
  geom_point(size=2)+facet_wrap(~cohort)+
  geom_linerange(aes(ymin=ci.lb,ymax=ci.ub),
                 size=2,alpha=0.3) +
  # scale_y_continuous(limits=c(0,96))+
  xlab("Age in months")+
  ylab("Point prevalence (95% CI)")+
  ggtitle("Cohort-specific point prevalence of stunting and wasting co-occurance - Asia")

print(p_asia)
dev.off()

# plot pooled prevalence
pdf("U:/Figures/co-occurance-ptprev-pool.pdf",width=10,height=4,onefile=TRUE)
p_pool <- ggplot(prev.res,aes(y=est,x=agecat))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.05) +
  scale_color_manual(values=tableau10)+xlab("Age category")+
  ylab("Point prevalence (95% CI)")+
  scale_y_continuous(limits=c(0,10))+
  annotate("text",x=prev.res$agecat,y=0,label=prev.res$nmeas.f,size=3)+
  annotate("text",x=prev.res$agecat,y=-3,label=prev.res$nstudy.f,size=3)+
  annotate("text",label=prev.res$ptest.f,x=prev.res$agecat,
           y=prev.res$est,hjust=-0.75,size=3)+
  ggtitle("Pooled point prevalence of stunting and wasting co-occurance")

print(p_pool)
dev.off()


# export
prev = dmn %>% 
  select(studyid,subjid,country,agecat,
         cowaststunt, cowaststunt)

save(prev,file="U:/Data/Co-occurrence/co_prev.RData")


#Relative proportion bar charts (Stacked bar chart)

#Strip grant identifier and add country
dmn$studyid <- gsub("^k.*?-" , "", dmn$studyid)
dmn$studyid <- paste0(dmn$studyid, ", ", paste0(substring(as.character(dmn$country),1,1), tolower(substring(as.character(dmn$country),2))))
dmn$studyid <- gsub("Tanzania, united republic of", "Tanzania", dmn$studyid)
dmn$studyid <- gsub("africa", "Africa", dmn$studyid)

library(scales)

barplot <- ggplot(dmn, aes(agecat, ..count..)) + geom_bar(aes(fill = status), position = "fill") +
  scale_fill_manual(values=tableau10) +
  xlab("Age") + ylab("Proportion") +
  scale_y_continuous(labels = percent_format()) + 
  ggtitle("Comparing distribution of co-occurrence, wasting only, stunting only, and healthy children at each age")


barplot_cohort_strat <- ggplot(dmn, aes(agecat, ..count..)) + geom_bar(aes(fill = status), position = "fill") +
  scale_fill_manual(values=tableau10) +
  facet_wrap(~studyid) +
  xlab("Age") + ylab("Proportion") +
  scale_y_continuous(labels = percent_format()) +
  theme(strip.background = element_blank(),
        legend.position="bottom",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        strip.text.y = element_text(angle = 0))


barplot_NH <- ggplot(dmn[dmn$status!="Healthy",], aes(agecat, ..count..)) + geom_bar(aes(fill = status), position = "fill") +
  scale_fill_manual(values=tableau10[-1]) +
  xlab("Age") + ylab("Proportion") +
  scale_y_continuous(labels = percent_format()) + 
  ggtitle("Comparing distribution of co-occurrence, wasting only, and stunting only\namong children wasted or stunted at each age")

#Make plot of the proportion stunted among wasted and not wasted, and vice-versa
d_wast <- dmn %>% filter(whz < -2) %>% mutate(status = ifelse(status=="Co-occurence", "Stunted", "Not stunted"))
d_stunt <- dmn %>% filter(haz < -2) %>% mutate(status = ifelse(status=="Co-occurence", "Wasted", "Not wasted"))

barplot_wast <- ggplot(d_wast, aes(agecat, ..count..)) + geom_bar(aes(fill = status), position = "fill") +
  scale_fill_manual(values=tableau10) +
  xlab("Age") + ylab("Proportion") +
  scale_y_continuous(labels = percent_format()) + 
  ggtitle("Proportion of stunted children among children who are wasted")

barplot_stunt <- ggplot(d_stunt, aes(agecat, ..count..)) + geom_bar(aes(fill = status), position = "fill") +
  scale_fill_manual(values=tableau10[-2]) +
  xlab("Age") + ylab("Proportion") +
  scale_y_continuous(labels = percent_format()) + 
  ggtitle("Proportion of wasted children among children who are stunted")


#---------------------------------------------------------------------
# Presentation sized plots.
#---------------------------------------------------------------------

#Combined pooled prev and bar chart
multiplot(p_pool,barplot_NH)


#Region stratified prevalence



#-----------------------------------
# save plot dataframe
#-----------------------------------
save(prev.res, dmn, file="U:/Data/Stunting/co_prev.RData")
