

#-----------------------------------
# Stunting analysis
# Objective 1a
# Import data, subset to relevant variables
#-----------------------------------
rm(list=ls())
library(dplyr)
library(ggplot2)
library(data.table)

#--------------------------------------------
# Read in .csv file and save as an .rds file
#--------------------------------------------

# setwd("U:/data/")
#  d<-fread("U:/data/FINAL/UCB Rally7/Main/adam/FINAL.csv", header = T)
d<-fread("U:/data/Stunting/Full-compiled-data/FINAL.csv", header = T)



#--------------------------------------------
# Subset to  just identifying and haz data
#--------------------------------------------

#change names to lower case
colnames(d) <- tolower(colnames(d))
d<-d %>% subset(., select=c(studyid, subjid, country, tr, agedays, haz, whz))

#--------------------------------------------
#Check for duplicate agedays
#--------------------------------------------
dup_age <- d %>% group_by(studyid, subjid, agedays) %>%
  summarize(N=n())
mean(dup_age$N)

#--------------------------------------------
# Mark monthly, quarterly, and yearly-measured 
# studies
#--------------------------------------------
#Drop studies Vishak added to data product that don't meet inclusion criteria
d <- d %>% filter(studyid!="ki1000301-DIVIDS" & studyid!="ki1055867-WomenFirst" & studyid!="ki1135782-INCAP")

#mark measure frequencies
d$measurefreq <- NA

d$measurefreq[d$studyid %in% c(
  "ki0047075b-MAL-ED",   
  "ki1000108-CMC-V-BCS-2002",              
  "ki1000108-IRC",               
  "ki1000109-EE",           
  "ki1000109-ResPak",  
  "ki1017093b-PROVIDE",  
  "ki1066203-TanzaniaChild2",           
  "ki1101329-Keneba",  
  "ki1112895-Guatemala BSC",       
  "ki1113344-GMS-Nepal",             
  "ki1114097-CONTENT"
)] <- "monthly"

d$measurefreq[d$studyid %in% c(
  "ki1112895-iLiNS-Zinc",  
  "kiGH5241-JiVitA-3",          
  "kiGH5241-JiVitA-4", 
  "ki1148112-LCNI-5",          
  "ki1017093-NIH-Birth",
  "ki1017093c-NIH-Crypto",   
  "ki1119695-PROBIT",         
  "ki1000304b-SAS-CompFeed",   
  "ki1000304b-SAS-FoodSuppl",   
  "ki1126311-ZVITAMBO",   
  "ki1114097-CMIN",                 
  "ki1135781-COHORTS"
)] <- "quarterly"

d$measurefreq[d$studyid %in% c(
  "ki1000110-WASH-Bangladesh",       
  "ki1000111-WASH-Kenya",  
  "ki1148112-iLiNS-DOSE",     
  "ki1148112-iLiNS-DYAD-M", 
  "ki1033518-iLiNS-DYAD-G",
  "ki1000125-AgaKhanUniv",           
  "ki1112895-Burkina Faso Zn",    
  "ki1000304-VITAMIN-A",  
  "ki1000304-Vitamin-B12",
  "ki1000107-Serrinha-VitA",   
  "ki1000304-EU",        
  "ki1000304-ZnMort"
)] <- "yearly"



#Mark COHORTS and CMIN cohorts with different measurement frequency than quarterly
d$measurefreq[d$studyid=="ki1114097-CMIN" & d$country=="BANGLADESH"] <- "monthly"
d$measurefreq[d$studyid=="ki1114097-CMIN" & d$country=="PERU"] <- "monthly"
d<- d[!(d$studyid=="ki1135781-COHORTS" & d$country=="BRAZIL"),] #Drop because yearly but not an RCT
d<- d[!(d$studyid=="ki1135781-COHORTS" & d$country=="SOUTH AFRICA"),] #Drop because yearly but not an RCT


#--------------------------------------------
# Save just identifying and haz data
#--------------------------------------------
save(d, file="U:/data/compiled_co-occurence_dataset.RData")




#--------------------------------------------
# Subset to relevant variables
#--------------------------------------------
colnames(d)=tolower(colnames(d))
d <- d %>% select(studyid, subjid, country, tr, agedays, haz, whz, measurefreq)

nrow(d)




#--------------------------------------------
# drop unrealistic HAZ
#--------------------------------------------
nrow(d)
d = filter(d,haz >= -6 & haz <=6)
d = filter(d,whz >= -6 & whz <=6)
nrow(d)

#--------------------------------------------
# order data, create measurement id
#--------------------------------------------
d <- d %>% 
  arrange(studyid,subjid,agedays) %>%
  group_by(studyid,subjid) %>%
  arrange(studyid,subjid,agedays) %>%
  # create id for measurement within person
  mutate(measid=seq_along(subjid)) 

# count number of studies
length(names(table(d$studyid)))

# table of studies
table(d$studyid)
table(d$studyid,d$country)

#convert subjid to character for the merge with covariate dataset
d$subjid <- as.character(d$subjid)


#--------------------------------------------
# Save intervention effects dataset
#--------------------------------------------

save(d,file="U:/Data/Co-occurrence/int_co-occurrence_data.RData")

#--------------------------------------------
# Subset to and save risk factor data
#--------------------------------------------

#Keep monthly and quarterly studies
d <- d %>% filter(measurefreq!="yearly")

save(d,file="U:/Data/Co-occurrence/rf_co-occurrence_data.RData")


#--------------------------------------------
# Subset to and save descriptive epi data
#--------------------------------------------

#--------------------------------------------
# drop trial arms with intervention impact on HAZ
# -either based on published literature or analysis
# of effects on CI of stunting by 24months of age
#--------------------------------------------
d=d[-which(d$studyid=="kiGH5241-JiVitA-4" & d$tr!="Control"),]
d=d[-which(d$studyid=="ki1119695-PROBIT" & d$tr!="Control"),]
d=d[-which(d$studyid=="ki1000304b-SAS-FoodSuppl" & d$tr!="Control"),]
d=d[-which(d$studyid=="ki1112895-iLiNS-Zinc" & d$tr!="Control"),]
d=d[-which(d$studyid=="ki1000304b-SAS-CompFeed" & d$tr!="Control"),]
d=d[-which(d$studyid=="kiGH5241-JiVitA-3" & d$tr!="Control"),]
d=d[-which(d$studyid=="ki1135781-COHORTS" & d$tr=="Other"),]


save(d,file="U:/Data/Co-occurrence/co-occurrence_data.RData")


#--------------------------------------------
# plot HAZ by agedays for included studies
#--------------------------------------------

plotdf <- d %>% subset(., select=c(studyid,country,agedays,haz,whz))
plotdf <- melt(plotdf, id.vars=1:3)

pdf("U:/Figures/co-occurrence-scatter-1a.pdf",width=15,height=15,onefile=TRUE)
ggplot(plotdf[plotdf$agedays<=365*2,],aes(x=agedays,y=value, color=variable))+geom_point(alpha=0.05)+geom_smooth()+
  facet_wrap(~studyid+country)+geom_hline(yintercept=-2,linetype="dashed",col="red")
dev.off()

summary(d$whz[d$studyid=="ki1126311-ZVITAMBO"])

