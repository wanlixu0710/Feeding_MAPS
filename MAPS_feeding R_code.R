setwd("~/R_directory/MAPS/feeding paper/MAPS feeding/")
rm(list = ls())
library(plyr)
library(dplyr) #dplyr has to load after plyr, otherwise the summarize function will take plyrs 

### descriptive stats about the stool sample and counts
totalreads <- read.csv(file="otu_counts.csv", sep=",",header=T) #read in the total reads
stool <- read.csv(file="MAPS_stool_MCB_processing set up021616.csv",sep=",", header=T)
stool <- stool[,(1:5)]
totalreads2 <- inner_join(totalreads, stool)
setdiff(totalreads$Fecal_Sample_ID_OLD, totalreads2$Fecal_Sample_ID_OLD) # check how many are not matching

totalreads3 <- totalreads2%>%filter((as.numeric(as.character(totalreads2$Subject_ID))!=33)&(as.numeric(as.character(totalreads2$Subject_ID))!=34)&(as.numeric(as.character(totalreads2$Subject_ID))!=39)) #remove 33, 34, 39 which only have very few samples

reads30 <- totalreads3%>%filter(as.numeric(as.character(totalreads3$PNA))<31) # yeild total 410
summary(reads30$count)
sd(reads30$count)

# retrive the samples without sequencing
sent <- read.csv(file="sent till 30.csv",sep=",", header=T)
Fecal_Sample_ID_OLD <- setdiff(sent$Fecal_Sample_ID_OLD, totalreads3$Fecal_Sample_ID_OLD)
count <- rep(0, time=16)
noDNA <- data.frame(Fecal_Sample_ID_OLD,count)
noDNA <- noDNA[-(noDNA$Fecal_Sample_ID_OLD%in%totalreads$Fecal_Sample_ID_OLD),]
noDNA <- inner_join(noDNA,stool)
noDNA$sample <- noDNA$Fecal_Sample_ID_OLD
noDNA <- filter(noDNA, noDNA, as.numeric(as.character(PNA))<31)
meconium <- filter(noDNA, as.numeric(as.character(PNA))<7) # I count any stool that is less than 7 days old meconium, so we have 9 meconium

totalsample <- rbind(reads30,noDNA)

count <- totalsample %>% group_by(Subject_ID) %>% tally () #average  stool collection
summary(count)
sd(count$n)

# total senquencing
sum(reads30$count) #25527333
lessthan10thousands <- filter(reads30, count<10000)

### demographic data
demo <- read.csv(file="MAPS demographic data_3_11_16.csv",sep=",", header=T)
str(demo)
demo2 <- demo%>% filter(demo$Subject_ID %in% reads30$Subject_ID)
table(demo2$Gender) #17/33
table(demo2$Baby_Race) #26/32
table(demo2$Baby_hispanic) #22/33
table(demo2$Delivery) #20/33
summary(demo2$BIrth_GA)
sd(demo2$BIrth_GA)
summary(demo2$Birth_weight)
sd(demo2$Birth_weight)

### descriptives for feeding
feeding <- read.csv(file="feeding_first_50.csv", sep=",", header=T)
firstfeeding <- demo2%>%select(Subject_ID, DOB, First_enteral_feeding)%>%mutate(fed=as.Date(First_enteral_feeding, "%m/%d/%Y")-as.Date(DOB, "%m/%d/%Y")+1)
summary(as.numeric(firstfeeding$fed))
sd((as.numeric(firstfeeding$fed)))

MBM <- feeding %>% filter(feeding$Subject_ID%in% reads30$Subject_ID & feeding$PNA<31) %>% mutate(mbm=MBM/(MBM+DBM+Formula+NPO)) %>% select(Subject_ID,PNA,MBM,mbm) #percentage of the MBM
summary(MBM$mbm) # 63.5% of the total feedings 

# MOM for each infant varies from 0 to 28 days (17.7 ± 7.9 days)
MBM1 <- MBM %>% filter(MBM$MBM==1) %>% select(Subject_ID, PNA, MBM) %>% distinct () %>% group_by (Subject_ID) %>% tally()
setdiff(MBM$Subject_ID, MBM1$Subject_ID)
Subject_ID <- 7
n <- 0
MBM2 <- data.frame(Subject_ID,n)
MBM2 <- rbind(MBM1, MBM2)
summary(MBM2)
sd(MBM2$n)

# NPO
str(feeding$NPO)
NPO <- feeding %>% filter(feeding$Subject_ID%in% reads30$Subject_ID & feeding$PNA<31) 

NPO <- NPO%>% filter(NPO$NPO==1) %>% select(Subject_ID, PNA) %>% distinct () 
NPO.first <- read.csv(file="NPO.before.first.feed.csv", sep=",", header=T)
NPO <- rbind(NPO,NPO.first)

NPO <- NPO%>% distinct()%>% group_by (Subject_ID) %>% tally()
Subject_ID <- setdiff(firstfeeding$Subject_ID,NPO$Subject_ID)
n <- rep(0, times=7)
NPO.first <- data.frame (Subject_ID, n)
NPO <- rbind(NPO, NPO.first)
summary(NPO)
sd(NPO$n)

