############### Descriptive calc ########################
setwd("~/R_directory/MAPS/feeding paper/MAPS feeding/")
rm(list = ls())
library(plyr)
library(dplyr) #dplyr has to load after plyr, otherwise the summarize function will take plyrs 

### descriptive stats about the stool sample and counts
#totalreads <- read.csv(file="otu_table_stats.csv", sep=",",header=T) #read in the total reads
stool <- read.csv(file="Updated_ MAPS_stool_MCB_processing set up04_11_16.csv",sep=",", header=T)
stool$PNA <- as.numeric(as.character(stool$PNA))
stool <- stool[,(1:5)]
#totalreads2 <- inner_join(totalreads, stool)
#setdiff(totalreads$Fecal_Sample_ID_OLD, totalreads2$Fecal_Sample_ID_OLD)
#write.csv(totalreads2, file="otu_stats_upload.csv")# get rid of the duplicates for each day and just save the one with the highest sequence
totalreads <- read.csv(file="otu_stats_upload.csv", stringsAsFactors = F)
totalreads2 <- totalreads%>%filter((as.numeric(as.character(totalreads$Subject_ID))!=33)&(as.numeric(as.character(totalreads$Subject_ID))!=34)&(as.numeric(as.character(totalreads$Subject_ID))!=39)) #remove 33, 34, 39 which only have very few samples

reads30 <- totalreads2%>%filter(as.numeric(as.character(totalreads2$PNA))<31)

#write.csv(reads30, file="reads30.csv")#calculate the interval 
reads30 <- read.csv(file="reads30.csv") #yeilds 406

summary(reads30$count)
sd(reads30$count)

a <- reads30 %>% count(Subject_ID)
summary(a$n)
sd(a$n)
reads30$interval[reads30$interval<1] <- NA
summary(reads30$interval)
sd(reads30$interval, na.rm=T)

# retrive the samples without sequencing
sent <- read.csv(file="sent till 30.csv",sep=",", header=T,stringsAsFactors = F)
sent <- left_join(sent, stool)
sent <- sent%>%filter(PNA<31)
Fecal_Sample_ID <- setdiff(sent$Fecal_Sample_ID, reads30$Fecal_Sample_ID);Fecal_Sample_ID
Daily_Code <- setdiff(sent$Daily_Code,reads30$Daily_Code);Daily_Code# only the first 30 babies have the problem of no reading
setdiff(reads30$Daily_Code, sent$Daily_Code)

count <- rep(0, time=13)
noDNA <- data.frame(Daily_Code,count)
noDNA$sample <- noDNA$Fecal_Sample_ID_OLD
noDNA <- left_join(noDNA, stool[,c("Daily_Code", "Subject_ID", "PNA")])
noDNA <- noDNA[!duplicated(noDNA),]
noDNA <- filter(noDNA,  as.numeric(as.character(PNA))<31)
meconium <- filter(noDNA, as.numeric(as.character(PNA))<7) # I count any stool that is less than 7 days old meconium, so we have 9 meconium
reads30.a <- reads30 %>% select(Daily_Code,Subject_ID, PNA) %>% distinct() 
totalsample <- rbind(reads30[, c("Daily_Code","count","Subject_ID","PNA")],noDNA) # totally 419 samples

count <- totalsample %>% group_by(Subject_ID) %>% tally () #average  stool collection
summary(count)
sd(count$n)

# calculate the interval of sample collection
#write.csv(totalsample,file="test.csv") #and then make a new variable as interval
totalsample <- read.csv(file="test.csv")
totalsample$interval[totalsample$interval<1] <- NA
summary(totalsample$interval)
sd(totalsample$interval, na.rm=T)

# total senquencing
sum(reads30$count) #25384479
lessthan10thousands <- filter(reads30, count<10000)

# get rid of the sequences that is less than 10000
reads30 <- reads30 %>% filter(count>=10000)
reads30 <- reads30[reads30$Fecal_Sample_ID_OLD!="Baby.021.E1",]

### demographic data
demo <- read.csv(file="MAPS demographic data_3_11_16.csv",sep=",", header=T)
str(demo)
demo2 <- demo%>% filter(demo$Subject_ID %in% reads30$Subject_ID)
table(demo2$Gender) #17/33
table(demo2$Baby_Race) #26/33
table(demo2$Baby_hispanic) #22/33
table(demo2$Delivery) #20/33
table(demo2$PROM)
table(demo2$Twins)
table(demo2$Resuscitation_1.yes_2.no)
summary(demo2$BIrth_GA)
sd(demo2$BIrth_GA)
summary(demo2$Birth_weight)
sd(demo2$Birth_weight)
summary(demo2$Birth_length)
sd(demo2$Birth_length)
summary(demo2$Birth_head_circumference)
sd(demo2$Birth_head_circumference, na.rm=T)
summary(demo2$SNAPEII)
sd(demo2$SNAPEII, na.rm=T )
summary(demo2$Mother_age)
sd(demo2$Mother_age)
# antibiotic use
anti <- read.csv(file="Daily_antibiotic.40.50.csv",stringsAsFactors = F)
anti.30 <- anti %>% filter(anti$Subject_ID %in% reads30$Subject_ID) %>% filter(PNA<=30)
setdiff( reads30$Subject_ID, anti.30$Subject_ID)
count <- anti.30%>%count(Subject_ID)#26 out of 33 infants have antibiotic use, but baby #40 doesn't have antibiotics during the first 3 days. there for 25 have antibiotics during the first 3 days.
25/33

### descriptives for feeding
feeding <- read.csv(file="feeding_first_50.csv", sep=",", header=T)
firstfeeding <- demo2%>%select(Subject_ID, DOB, First_enteral_feeding)%>%mutate(fed=as.Date(First_enteral_feeding, "%m/%d/%Y")-as.Date(DOB, "%m/%d/%Y")+1)
summary(as.numeric(firstfeeding$fed))
sd((as.numeric(firstfeeding$fed)))

MBM <- feeding %>% filter(feeding$Subject_ID %in% reads30$Subject_ID & feeding$PNA<31) %>% mutate(mbm=MBM/(MBM+DBM+Formula+NPO)) %>% select(Subject_ID,PNA,MBM,mbm) #percentage of the MBM
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


############ Data analysis ###################
rm(list = ls())
library(vegan)
library(ecodist)
library(ggplot2)
library(plyr)
library(scatterplot3d)
library(lattice)
library(dplyr) #dplyr has to load after plyr, otherwise the summarize function will take plyrs 
library(reshape2)
library(RSvgDevice)
library(reshape)
library(indicspecies)
library(RColorBrewer)
library(heatmap3)
library(devEMF)
library(lme4)
library(lmerTest)

# read in data
map <- read.csv(file="reads30.csv", sep=",", header=T, stringsAsFactors = F)
map <- map[,-c(8)]
map <- map %>% filter(count>=10000)
map <- map[map$Fecal_Sample_ID_OLD!="Baby.021.E1",] # Baby021.E1 is contaminated with rare bacterias and oceanospirillales
stool <- read.csv(file="Updated_ MAPS_stool_MCB_processing set up04_11_16.csv", sep=",", header=T)
stool <- stool[, (1:5)]
a.div <- read.csv(file="alpha_even.csv", sep=",", header=T, stringsAsFactors = F)
map$PNA <- as.numeric(as.character(map$PNA))
stool$PNA <- as.numeric(as.character(stool$PNA))
setdiff( map$sample,a.div$sample)
setdiff(a.div$sample, map$sample)
a.div <- inner_join(a.div, map, by="sample")

demo <- read.csv(file="MAPS demographic data_3_11_16.csv",sep=",", header=T)
nnns <- read.csv(file="NNNS_first 50.csv")
nnns <- nnns[,c(1,3)]
demo <-inner_join(demo, nnns) 
a.div <- left_join(a.div, demo)

# new coding based on 50% or more for each 30 days: categorized to 3 groups
feeding <- read.csv(file="feeding_first_50.csv", sep=",", header=T)
MBM <- subset(feeding, MBM==1)
DBM <- subset(feeding, DBM==1)
Formula <- subset(feeding, Formula==1)
MBM$FT <- MBM$MBM
DBM$FT <- DBM$DBM+1
Formula$FT <- Formula$Formula+2
feeding <- bind_rows(MBM, DBM, Formula)
feeding$PNA <- as.numeric(feeding$PNA)
#check the first days of enteral feeding
#minfeeding <- feeding%>% select(Subject_ID, PNA, FT)%>%group_by(Subject_ID) %>% summarize(minfedpna=min(PNA))
#a <- a.div%>%select(Subject_ID, BIrth_GA)
#a <- a[!duplicated(a),]
#minfeeding <- inner_join(minfeeding,a)
#summary(minfeeding)
#sd(minfeeding$minfedpna)
feeding$day[feeding$PNA < 11] <- "1st 10 days"
feeding$day[feeding$PNA> 10 &feeding$PNA <= 20] <- "2nd 10 days"
feeding$day[feeding$PNA> 20 &feeding$PNA <= 30] <- "3rd 10 days"

#count <- feeding %>% filter(!is.na(day)) %>% filter(!is.na(Formula)) %>% group_by(Subject_ID,day)%>%summarise(summbm=sum(MBM, na.rm=TRUE), sumdbm=sum(DBM, na.rm=TRUE), sumformula=sum(Formula, na.rm=TRUE))
#count <- count %>% mutate(pmbm=summbm/(summbm+sumdbm+sumformula),pdbm=sumdbm/(summbm+sumdbm+sumformula),pform=sumformula/(summbm+sumdbm+sumformula))

# for 70 percent cut for each ten days
#count$ten.days.fed[count$pmbm >=0.7 & count$pdbm<0.7 & count$pform<0.7] <- "4"
#count$ten.days.fed[count$pmbm <0.7 & count$pdbm>=0.7 & count$pform<0.7] <- "5"
#count$ten.days.fed[count$pmbm >=0.7 & count$pdbm<0.7 & count$pform>=0.7] <- "6"
#count$ten.days.fed[count$pmbm >=0.7] <- "1"
#3count$ten.days.fed[count$pdbm >=0.5 &count$pmbm <0.7] <- "2"
#count$ten.days.fed[count$pform >0.7] <- "3"
#c <- select(count, Subject_ID, day, ten.days.fed)
#a.div$PNA <- as.numeric(a.div$PNA)
#a.div$day[a.div$PNA < 11] <- "1st 10 days"
#a.div$day[a.div$PNA> 10 &a.div$PNA <= 20] <- "2nd 10 days"
#a.div$day[a.div$PNA> 20 &a.div$PNA <= 30] <- "3rd 10 days"
#a.div <- left_join(a.div, c)

#manually change 21	2nd 10 days to 4, 23	1st 10 days to 4 and 25 2nd 10 days to 5
# write.csv(count, file="count.for.upload.csv")
#### 70% for 6 groups and for each 10 days
count <- read.csv(file="count.for.upload.csv", sep=",", header=T)

c <- select(count, Subject_ID, day, ten.days.fed)

a.div2 <- a.div 
a.div2$PNA <- as.numeric(as.character(a.div$PNA))

a.div2$day[a.div2$PNA < 11] <- "1st 10 days"
a.div2$day[a.div2$PNA> 10 &a.div2$PNA <= 20] <- "2nd 10 days"
a.div2$day[a.div2$PNA> 20 &a.div2$PNA <= 30] <- "3rd 10 days"
a.div.30 <- filter(a.div2, PNA<31)
a.div.30 <- inner_join(a.div.30, c)

a.div.30$subject.recode <- revalue(as.character(a.div.30$Subject_ID),  c("1"="A", "2"="B","3"="C","4"="D","5"="E","6"="F", "7"="G","8"="H","9"="I","10"="J","11"="K","12"="L", "13"="M","14"="N","15"="O","16"="P","18"="Q","19"="R","20"="S","21"="T","22"="U","23"="V","24"="W","25"="X", "26"="Y","27"="Z","28"="ZA","29"="ZB","30"="ZC", "35"="ZD","36"="ZE","40"="ZF", "41"="ZG"))

a.div.30$mo <- revalue(as.character(a.div.30$ten.days.fed), c("1"="1.MOM", "2"="4. HDM", "3"="5.Formula", "4"='2.MOM+HDM', "5"="3.MOM+Formula", "6"="6.HDM+Formula"))


a.div.30 <- a.div.30 %>% filter(!is.na(a.div.30$ten.days.fed))
SPG1 <-ggplot (a.div.30, aes(x=subject.recode,y=PNA,colour=simpson))+geom_point(size=5)+scale_colour_gradient(low="green", high="blue")+labs(x="Infant", y="Postnatal age (day)", colour="Simpson",title="Gini-Simpson Index for feeding types")+facet_grid(day~mo, scales = "free",shrink=T, space = "free")+scale_y_continuous ( breaks=1:30)+ theme_bw() #theme(axis.title.y = element_text(vjust=-5))+#scale_x_continuous (breaks=1:29)
SPG1

### plot simpson diversity
SPG2 <-ggplot (a.div.30, aes(x=PNA,y=simpson,colour=as.factor(mo)))+geom_point()+stat_smooth(colour="darkblue", alpha = 0.1,se = FALSE) +labs(x="Postnatal Day (PNA)", y="Gini-Simpson Diversity Index", colour="Feeding Group", title="Gini-Simpson Index for feeding types")+facet_grid(mo ~ day, scales = "free_x",shrink=T, space = "free_x")+scale_x_continuous (breaks=1:30)+theme(axis.text.x = element_text(angle = 45)) #theme(axis.title.y = element_text(vjust=-5))
SPG2

SPG3 <-ggplot (a.div.30, aes(x=PNA,y=simpson,colour=mo)) + geom_point(aes(color= mo, shape= mo, alpha=0.7)) + stat_summary(aes(colour= mo ,shape= mo ,group=mo), fun.y=mean, geom="line", size=1.5) + labs(x="Postnatal Day (PNA)", y="Gini-Simpson Diversity Index", colour="Feeding Group", title="Gini-Simpson Index for feeding types")+ facet_grid(~day, scales = "free_x",shrink=T, space = "free_x")+  scale_x_continuous (breaks=1:30)+ scale_color_manual(values=c("#004E00", "#33FF00", "#FF9966", "#3399FF", "#FF004C","#BCBDDC")) + theme(axis.text.x = element_text(angle = 45)) #theme(axis.title.y = element_text(vjust=-5))
SPG3


a.div.30$Feeding_type <- as.factor(a.div.30$mo)
a.div.30$meansimp <- ave(a.div.30$simpson, by=list (a.div.30$Subject_ID,a.div.30$day), FUN= mean)
fig3 <- ggplot(a.div.30,  aes( y=simpson, x=mo))+ #order=- makes the stack the same order as the legend
  geom_boxplot(alpha=0.5, aes(fill = Feeding_type, shape=factor(Feeding_type)))+
  scale_x_discrete()+
  labs(x="Feeding types", y="Gini-Simpson Diversity Index")+
  facet_grid(~day,shrink=T, space = "free")+
  theme(axis.text.x = element_text(angle = 45), plot.background= element_rect(fill=NULL, colour = NULL),legend.position = "none")
fig3



### TAXA
L4 <- read.csv(file="even_table_L4.csv", header=T, sep=",")
L4 <- inner_join(map,L4)
L4 <- L4[, -(1:5)]
#a <-a.div.30[,c("Fecal_Sample_ID","infant.fed.70","PNA")]
a <-a.div.30[,c("Fecal_Sample_ID","ten.days.fed","day","PNA")]
L4 <- inner_join(a, L4)
str(L4)
order.color<-c("Other" = "#ededed", 
               "Fusobacteriales"= "ivory",
               "Actinomycetales"="#DEEBF7",
               "Bifidobacteriales" = "#9ECAE1", 
               "Bacteroidales" = "#4292C6",
               "Bacillales" = "#08519C", 
               "Gemellales" = "#E5F5E0", 
               "Lactobacillales" = "#A1D99B", 
               "Clostridiales" = "#41AB5D", 
               "Fusobacteriales" = "#006D2C", 
               "Caulobacterales" = "#FFF7BC", 
               "Rhizobiales" = "#FEE391", 
               "Rhodobacterales" = "#FEC44F", 
               "Burkholderiales" = "#EFEDF5", 
               "Campylobacterales" ="#BCBDDC" , 
               "Aeromonadales" = "#807DBA", 
               "Alteromonadales" = "#54278F",
               "Enterobacteriales" = "#FCBBA1", 
               "Oceanospirillales" = "#FB6A4A", 
               "Pasteurellales" = "black", 
               "Pseudomonadales" = "#CB181D", 
               "Xanthomonadales" = "#67000D")

# rename the legend
L4 <- rename(L4, c("k__Bacteria.Other.Other.Other"="Other", 
                   "k__Bacteria.p__Fusobacteria.c__Fusobacteria.o__Fusobacteriales"="Fusobacteriales",
                   "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales"="Actinomycetales",
                   "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Bifidobacteriales"="Bifidobacteriales",
                   "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales" = "Bacteroidales",
                   "k__Bacteria.p__Firmicutes.Other.Other"="Other",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.Other"="Other",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales"="Bacillales",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.o__Gemellales"="Gemellales",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales"="Lactobacillales",
                   "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales"="Clostridiales",
                   "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales" = "Fusobacteriales", 
                   "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales" = "Caulobacterales", 
                   "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales" = "Rhizobiales", 
                   "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodobacterales" = "Rhodobacterales", 
                   "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales" = "Burkholderiales", 
                   "k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales" ="Campylobacterales" , 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.Other" =  "Other", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Aeromonadales" = "Aeromonadales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales" = "Alteromonadales",
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales" = "Enterobacteriales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales" = "Oceanospirillales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales" = "Pasteurellales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales" = "Pseudomonadales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales" = "Xanthomonadales"))

# melt the data
#order<- melt(L4, id.vars=c("Fecal_Sample_ID","infant.fed.70","PNA","day","f"))
order<- melt(L4, id.vars=c("Fecal_Sample_ID","ten.days.fed","PNA","day"))

# ggplot aggregation makes a box for each sample=horizontal lines in bars
#order2 <- aggregate(value~variable+Fecal_Sample_ID+day+f, data=order, FUN=mean) 
str(order)
order$mo <- revalue(as.character(order$ten.days.fed), c("1"="1.MOM", "2"="4. HDM", "3"="5.Formula", "4"='2.MOM+HDM', "5"="3.MOM+Formula", "6"="6.HDM+Formula"))
order2 <- aggregate(value~variable+Fecal_Sample_ID+day+mo, data=order, FUN=mean) 

order3 <- order2[order(order2$variable, decreasing = TRUE),]
feeding.plot3 <- ggplot(order3,  aes( y=as.numeric(value), x=factor(day), color=NULL, fill=factor(variable), order=-as.numeric(variable)))+ #order=- makes the stack the same order as the legend
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=order.color, name="Key")+
  scale_x_discrete()+
  labs(x="Time intervals within feeding type", y="Order-level of bacteria (Relative abundance %)")+ 
  facet_grid(~mo,shrink=T, space = "free")+
  theme(axis.text.x = element_text(angle = 45))#+
#theme_bw()
feeding.plot3

### bifido boxplot
L4.2 <- inner_join(L4[,c("Fecal_Sample_ID","Bifidobacteriales","Lactobacillales")], a.div.30[,c("Fecal_Sample_ID", "mo","day","Feeding_type")])
bifido.plot <- ggplot(L4.2,  aes( y=as.numeric(Bifidobacteriales), x=factor(day)))+ geom_boxplot()+labs(x= "Feeding types", y="Relative abundance (%) of Bifidobacteriales")+facet_grid(~mo,shrink=T, space = "free")+theme(axis.text.x = element_text(angle = 45));bifido.plot

### lacto boxplot
lacto.plot <- ggplot(L4.2,  aes( y=as.numeric(Lactobacillales), x=factor(day)))+ geom_boxplot()+labs(x="Feeding types", y="Relative abundance (%) of Lactobacillales")+facet_grid(~mo,shrink=T, space = "free")+theme(axis.text.x = element_text(angle = 45));lacto.plot

##############fig 2##########typical individual plot##########
L4 <- read.csv(file="even_table_L4.csv", header=T, sep=",")
L4 <- inner_join(map,L4)
L4 <- L4[, -(1:2)]
a.div.30 <- subset(a.div, as.numeric(as.character(PNA))<31)
a <-a.div.30[,c("Fecal_Sample_ID","PNA","Subject_ID")]
a.indi <- filter(a, Subject_ID==12|Subject_ID==21| Subject_ID==7|Subject_ID==24)
L4 <- inner_join(a.indi, L4)

str(L4)
order.color<-c("Other" = "#ededed", "Fusobacteriales"= "ivory","Actinomycetales"="#DEEBF7", "Bifidobacteriales" = "#9ECAE1", "Bacteroidales" = "#4292C6", "Bacillales" = "#08519C",  "Gemellales" = "#E5F5E0",  "Lactobacillales" = "#A1D99B",  "Clostridiales" = "#41AB5D",  "Fusobacteriales" = "#006D2C", "Caulobacterales" = "#FFF7BC", "Rhizobiales" = "#FEE391", "Rhodobacterales" = "#FEC44F",  "Burkholderiales" = "#EFEDF5", "Campylobacterales" ="#BCBDDC" ,  "Aeromonadales" = "#807DBA",   "Alteromonadales" = "#54278F", "Enterobacteriales" = "#FCBBA1", "Oceanospirillales" = "#FB6A4A", "Pasteurellales" = "black", "Pseudomonadales" = "#CB181D", "Xanthomonadales" = "#67000D")

# rename the legend
L4 <- rename(L4, c("k__Bacteria.Other.Other.Other"="Other", 
                   "k__Bacteria.p__Fusobacteria.c__Fusobacteria.o__Fusobacteriales"="Fusobacteriales",
                   "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales"="Actinomycetales",
                   "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Bifidobacteriales"="Bifidobacteriales",
                   "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales" = "Bacteroidales",
                   "k__Bacteria.p__Firmicutes.Other.Other"="Other",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.Other"="Other",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales"="Bacillales",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.o__Gemellales"="Gemellales",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales"="Lactobacillales",
                   "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales"="Clostridiales",
                   "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales" = "Fusobacteriales", 
                   "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales" = "Caulobacterales", 
                   "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales" = "Rhizobiales", 
                   "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodobacterales" = "Rhodobacterales", 
                   "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales" = "Burkholderiales", 
                   "k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales" ="Campylobacterales" , 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.Other" =  "Other", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Aeromonadales" = "Aeromonadales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales" = "Alteromonadales",
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales" = "Enterobacteriales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales" = "Oceanospirillales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales" = "Pasteurellales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales" = "Pseudomonadales", 
                   "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales" = "Xanthomonadales"))

# melt the data
L4$infant.rename <- revalue(as.character(L4$Subject_ID),  c("12"="A. MOM", "21"="B. MOM+HDM to MOM+HDM+Formula to HDM+Formula","7"="C. HDM","24"="D. HDM to HDM+Formula"))
ord<- melt(L4, id.vars=c("Fecal_Sample_ID","Subject_ID","PNA","infant.rename", "Daily_Code","Fecal_Sample_ID_OLD"))

ord2 <- ord[order(ord$variable, decreasing = TRUE),]
str(ord)
ord2$PNA <- as.numeric(as.character(ord2$PNA))
fig2 <- ggplot(ord2,  aes( y=as.numeric(value), x=factor(PNA), color=NULL, fill=factor(variable), order=-as.numeric(variable)))+ #order=- makes the stack the same order as the legend
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=order.color, name="Key")+
  scale_x_discrete()+
  xlab("Postnatal age (day)")+ ylab("Order-level of bacteria (Relative abundance %)")+
  facet_wrap(~infant.rename,scales = "free_x",shrink=T,  ncol=2)+
  theme_bw()
fig2


### ANOVA : one way anova , didn't work
aov(simpson~mo, data=a.div.30)

############----indicator species---genus level-##########
# the whole feeding type
otu <- read.csv(file="even_table_L6.csv",  header=T, sep=",")
otu <- inner_join(map,otu)
otu <- otu[,-c(1:5)]
# recode feeding type to one variable
#FT <- a.div.30[, c("Fecal_Sample_ID", "infant.fed.70","PNA")]

FT <- select(a.div.30, Fecal_Sample_ID, mo, ten.days.fed, subject.recode, Gender)
FT$ten.days.fed <- as.character(FT$ten.days.fed)
FT <- inner_join(FT,otu)
# delete the missing values
FT<-FT[complete.cases(FT), ]
FT <- FT[!duplicated(FT),]
otu.f<-FT[,(7:68)]
ft<- FT[,(1:6)]
# assign Fecal_Sample_ID as the otu row name
rownames(otu.f) <- ft[,c("Fecal_Sample_ID")]
# indicator species
FT.ind  <-  multipatt(otu.f, ft$mo,
                  control = how(nperm=350))
levels(ft$FT.ind)
summary(FT.ind)




################################ beta-diversity #######################################

br_dist <- vegdist(otu.f, method="bray"); br_dist

#######-----------NMS -------------
# nms <- metaMDS(otu.f, dist="bray", k=2, trymax=250, wascores=TRUE, trymin=50) # 2D, stress is too high > 2.5.

# bray-curtis
nms <- metaMDS(otu.f, dist="bray", k=3, trymax=250, wascores=TRUE, trymin=50)
stressplot(nms)
stat <- data.frame(nms$points)
stat$Fecal_Sample_ID <- rownames(stat)
stat.d <- inner_join(stat, ft)

stat.d$color <- revalue(as.character(stat.d$mo), c("1.MOM"="red", "4. HDM"="purple", "5.Formula"="yellow", "2.MOM+HDM"="green", "3.MOM+Formula"="blue", "6.HDM+Formula"="black"))
stat.d$shape <- revalue(as.character(stat.d$mo), c("1.MOM"="15", "4. HDM"="17", "5.Formula"="16", "2.MOM+HDM"="12", "3.MOM+Formula"="14", "6.HDM+Formula"="25"))
scatterplot3d(stat.d[,1:3], color=stat.d$color,  pch=as.numeric(stat.d$shape), grid=TRUE, box=FALSE)
legend("topright", legend = levels(factor(stat.d$mo)), pch =  c (15,  12, 14, 17, 16, 25), col=c("red","green", "blue", "purple", "yellow", "black"), cex = 0.75)

# Jaccard
nms <- metaMDS(otu.f, dist="jaccard", k=3, trymax=800, wascores=TRUE, trymin=50)
stressplot(nms)
stat <- data.frame(nms$points)
stat$Fecal_Sample_ID <- rownames(stat)
stat.d <- inner_join(stat, ft)

stat.d$color <- revalue(as.character(stat.d$mo), c("1.MOM"="red", "4. HDM"="purple", "5.Formula"="yellow", "2.MOM+HDM"="green", "3.MOM+Formula"="blue", "6.HDM+Formula"="black"))
scatterplot3d(stat.d[,1:3], color=stat.d$color, main="1.MOM=red,4. HDM=purple, 5.Formula=yellow, 2.MOM+HDM=green, 3.MOM+Formula=blue, 6.HDM+Formula=black", pch=16, grid=TRUE, box=FALSE)



###-------------permanova bray curtis----------------------------------------------
otu <- read.csv(file="even_table_L6.csv",  header=T, sep=",")
otu <- inner_join(map,otu)
otu <- otu[,-c(1:5)]


a.div.30$subject.recode <- revalue(as.character(a.div.30$Subject_ID),  c("1"="A", "2"="B","3"="C","4"="D","5"="E","6"="F", "7"="G","8"="H","9"="I","10"="J","11"="K","12"="L", "13"="M","14"="N","15"="O","16"="P","18"="Q","19"="R","20"="S","21"="T","22"="U","23"="V","24"="W","25"="X", "26"="Y","27"="Z","28"="ZA","29"="ZB","30"="ZC", "35"="ZD","36"="ZE","40"="ZF", "41"="ZG"))

a.div.30$mo <- revalue(as.character(a.div.30$ten.days.fed), c("1"="1.MOM", "2"="4. HDM", "3"="5.Formula", "4"='2.MOM+HDM', "5"="3.MOM+Formula", "6"="6.HDM+Formula"))

a.div.30 <- a.div.30 %>% filter(!is.na(a.div.30$ten.days.fed))
# recode feeding type to one variable
#FT <- a.div.30[, c("Fecal_Sample_ID", "infant.fed.70","PNA")]
#niss <- read.csv(file="niss_first_50.csv",sep=",",header=T)
#weight <- niss[, c("Subject_ID","PNA", "weight")]
anti <- read.csv(file="Daily_antibiotic.40.50.csv",stringsAsFactors = F)
anti.30 <- anti %>% filter(anti$Subject_ID %in% a.div.30$Subject_ID) %>% filter(PNA<=30)
anti.count <- anti.30%>%count(Subject_ID)%>%mutate(Numbers.of.antibiotic=n)

#anti <- read.csv(file="antibiotic.base.csv", sep=",",header=T)
#anti2 <- read.csv(file="antibiotic.average.number.csv", sep=",", header=T)
#anti2 <-anti2[,-1] 
#a.div.30 <- left_join(a.div.30, weight)
a.div.30 <- left_join(a.div.30,anti.count)
#a.div.30 <- left_join(a.div.30,anti2)
FT2 <- select(a.div.30, Subject_ID, Fecal_Sample_ID, mo, PNA, Numbers.of.antibiotic, Gender, BIrth_GA, PROM )
FT2$Numbers.of.antibiotic[is.na(FT2$Numbers.of.antibiotic)] <- 0
FT2 <- inner_join(FT2,otu)

otu.f<-FT2[,(9:70)]
ft<- FT2[,(1:8)]
## Calculation of bray curtis dissimilarity
br_dist <- vegdist(otu.f, method="bray")
str(br_dist)

# permanova for feeding
br_perm <- adonis(br_dist ~ft$PNA+ ft$Numbers.of.antibiotic+ ft$BIrth_GA+ ft$PROM+ ft$mo*ft$Gender, na.rm=TRUE, sqrt.dist=TRUE, permutations=999,strata=factor(ft$Subject_ID))
br_perm


## Calculation of jaccard dissimilarity
jac_dist <- vegdist(otu.f, method="jaccard")
str(jac_dist)

# permanova for feeding
adonis(jac_dist ~ft$PNA+ ft$Numbers.of.antibiotic+ ft$BIrth_GA+ ft$PROM+ ft$mo*ft$Gender, na.rm=TRUE, permutations=999,strata=factor(ft$Subject_ID))
betadisper(jac_dist, ft$mo)



### pairwise permanova
pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

PW.Adonis <- pairwise.adonis(otu.f,ft$mo,sim.method="bray",p.adjust.m = "bonferroni")
PW.Adonis

pairwise.adonis(otu.f,ft$mo,sim.method="jaccard",p.adjust.m = "bonferroni")


### check the number of participants for each category that are used in the data analysis
c <- a.div.30 %>% select(Subject_ID, day, mo) %>%unique()
a <- c %>% filter(Subject_ID %in% a.div.30$Subject_ID) %>% filter( day=="1st 10 days") 
a %>% group_by(mo) %>% tally()
a <- c %>% filter(Subject_ID %in% a.div.30$Subject_ID) %>% filter( day=="2nd 10 days")
a %>% group_by(mo) %>% tally()
a <- c %>% filter(Subject_ID %in% a.div.30$Subject_ID) %>% filter( day=="3rd 10 days")
a %>% group_by(mo) %>% tally()

### check the number of participants for each category that of all of the feedings for the 33 babies
c <- count %>% select(Subject_ID, day, ten.days.fed) %>%unique()
c$mo <- revalue(as.character(c$ten.days.fed), c("1"="1.MOM", "2"="4. HDM", "3"="5.Formula", "4"='2.MOM+HDM', "5"="3.MOM+Formula", "6"="6.HDM+Formula"))
a <- c %>% filter(Subject_ID %in% a.div.30$Subject_ID) %>% filter( day=="1st 10 days")
a %>% group_by(mo) %>% tally()
a <- c %>% filter(Subject_ID %in% a.div.30$Subject_ID) %>% filter( day=="2nd 10 days")
a %>% group_by(mo) %>% tally()
a <- c %>% filter(Subject_ID %in% a.div.30$Subject_ID) %>% filter( day=="3rd 10 days") 
a %>% group_by(mo) %>% tally()
setdiff (a.div.30$Subject_ID, a$Subject_ID)



#####data prepare for sas mixed effect
a.div <- read.csv(file="alpha_even.csv", sep=",", header=T)
a.div <- select(a.div, sample,simpson)
a.div <- inner_join(a.div, map)
a.div <- inner_join(a.div, stool)
a.div$PNA <- as.numeric(as.character(a.div$PNA))
a.div$day[a.div$PNA < 11] <- "1st 10 days"
a.div$day[a.div$PNA> 10 &a.div$PNA <= 20] <- "2nd 10 days"
a.div$day[a.div$PNA> 20 &a.div$PNA <= 30] <- "3rd 10 days"
a.div2 <- aggregate(a.div, by=list(a.div$Subject_ID, a.div$day),FUN=mean, na.rm=TRUE)
a.div2$day <- a.div2$Group.2
a.div2$Subject_ID <- a.div2$Group.1
a.div3 <- select(a.div2, Subject_ID, day, simpson)
demo <- read.csv(file="demo_first 50.csv", sep=",", header=T)
nnns <- read.csv(file="NNNS_first 50.csv", sep=",", header=T)
nnns <- nnns[,c(1,3)]
demo <-inner_join(demo, nnns) 


feeding <- read.csv(file="feeding_first_50.csv", sep=",", header=T)
feeding$day[feeding$PNA < 11] <- "1st 10 days"
feeding$day[feeding$PNA> 10 &feeding$PNA <= 20] <- "2nd 10 days"
feeding$day[feeding$PNA> 20 &feeding$PNA <= 30] <- "3rd 10 days"

feeding2 <- aggregate(feeding, by=list(feeding$Subject_ID, feeding$day), FUN=mean, na.rm=TRUE)

feeding2$day <- feeding2$Group.2
feeding2$Subject_ID <- feeding2$Group.1
feeding3 <- feeding2[,-(1:2)]

niss <- read.csv(file="niss_first_50.csv",sep=",", header=T)
niss$day[niss$PNA < 11] <- "1st 10 days"
niss$day[niss$PNA> 10 &niss$PNA <= 20] <- "2nd 10 days"
niss$day[niss$PNA> 20 &niss$PNA <= 30] <- "3rd 10 days"
#niss2 <- niss%>% group_by(Subject_ID, day)%>%summarise (meanweight=mean(weight) )%>% select( Subject_ID, day, meanweight) %>% filter (!is.na(day))
niss2<- select(niss, Subject_ID, day, weight)
niss2 <- aggregate(niss$weight, by=list(niss$Subject_ID, niss$day), FUN=mean, na.rm=TRUE)
niss2$day <- niss2$Group.2
niss2$Subject_ID <- niss2$Group.1
niss2$weight <- niss2$x
niss3 <- niss2[, (4:6)]
sas <- left_join(a.div3, niss3)
sas <- left_join(sas, feeding3)
sas <- left_join(sas, demo)
sas <- left_join(sas, nnns)

med <- read.csv(file="Daily_antibiotic.csv",sep=",", header=T)
med1 <- select(med, Subject_ID, day, antibiotic.use, Numbers.of.antibiotic)
med1$num.anti <- (ave(med1$Numbers.of.antibiotic, med1$Subject_ID, med1$day, na.rm=TRUE, FUN=sum))/7
med2 <- select(med1, Subject_ID, day, num.anti)
med2 <- med2[!duplicated(med2),]

anti <- read.csv(file="antibiotic.base.csv",sep=",", header=T)
sas <- left_join(sas, med2)
sas <- left_join(sas, anti)
write.csv(sas, file="feedingdata_for_sas_3_24_16.csv")



