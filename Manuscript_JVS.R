### DESKTOP - Set Working Directory ###
setwd("C:/Users/alyoung6/Box/THESIS DATA/Manuscript")

### LAPTOP - Set Working Directory ###
setwd("C:/Users/alyou/Box/THESIS DATA/Manuscript")

### LOAD PACKAGES ###
library(vegetarian)
library(MASS)
library(BMA)
library(gdata)
library(leaps)
library(lmSupport)
library(gtools)
library(gridExtra)
library(vegan)
library(reshape2)
library(gtable)
library(relaimpo)
library(doBy)
library(reldist)
library(tidyverse)
library(codyn)
library(ggplot2)
library(grid)
library(cowplot)
library(vegan)
library(devtools)
library(ggfortify)
library(data.table)
library(GGally)
library(gridGraphics)
library(lme4)
library(lmerTest)
library(sjPlot)
library(splines2)
library(AER)
library(xlsx)
library(MuMIn)
library(factoextra)
library(nlme)
library(AICcmodavg)
library(dplyr)
library(ggpubr)
library(lavaan)
library(readr)
theme_set(theme_bw(14))

### READ IN DATA FILES ###

# species composition #
SpComp<-read.csv("SpeciesCompositionData.csv", header=T)%>%
  mutate(SpNum=paste("sp",SpeciesID, sep=""))%>%
  select(-SpeciesID)

# species unique IDs #
SpList<-read.csv("SpeciesUniqueID.csv", header=T)%>%
  mutate(SpNum=paste("sp",SpeciesID, sep=""))

# merge datasets #
SpComp2<- merge(SpComp, SpList, by=c("SpNum"))%>%
  mutate(Site_Plot=paste(Site,Plot,sep="_"))%>%
  select(Site,Plot,Site_Plot,SpNum,Cover,SpeciesID,SpeciesName,FunctionalGroup)

# shrub data #
ShrubData<-read.csv("ShrubNumbers.csv", header=T)%>%
  select(Site,Plot,Stems)

# Site metadata, biomass x 10, light data and calculations #
SiteData<-read_csv("SiteData.csv",col_types=cols(ANPPpyd=col_number(),ANPPpine = col_number()))%>%
  group_by(Site,Plot)%>%
  mutate(Open=(OpenSkyrep1+Openskyrep2)/2)%>%
  mutate(Above=(Lightrep1ab+Lightrep2ab+Lightrep3ab)/3)%>%
  mutate(Below=(Lightrep1be+Lightrep2be+Lightrep3be)/3)%>%
  mutate(SoilSurfaceLight=(Below/Open)*100)%>%
  mutate(UnderstoryLight=(Above/Open)*100)%>%
  group_by(Site,Plot)%>%
  mutate(anppgrass=(ANPPgrass*10))%>%
  mutate(anppwoody=(ANPPwoody*10))%>%
  mutate(anppforb=(ANPPforb*10))%>%
  mutate(anpppyd=(ANPPpyd*10))%>%
  mutate(litter=(ANPPpine*10))%>%
  mutate(plotanpp=sum(anppgrass,anppforb,anppwoody,na.rm=T))%>%
  group_by(Site)%>%
  mutate(siteanpp=sum(anppgrass,anppforb,anppwoody,na.rm=T))%>%
  mutate(meansiteanpp=mean(plotanpp,na.rm=T))

            # # # #
# nutrient data #
Nutrients<-read.csv("SoilNutrients.csv",header=T)%>%
  select(Site,Plot,CEC, BaseSat.,pH,P,K,Ca,Mg,Mn,Zn,Cu,S,Na)
####################################################
# PCA with CATION stuff included data #
PCA<-prcomp(Nutrients[,3:14],scale=TRUE)
PCA
summary(PCA)

###########################################

axes <- predict(PCA, newdata = Nutrients)
head(axes, 4)
  
# put PCA axes with site and plot #   
pca<-cbind(Nutrients,axes)%>%
  select(Site,Plot,PC1,PC2)

# find contributions of nutrients to PCA axes #
var <- get_pca_var(PCA)
var
head(var$contrib, 12)
            # # # #

# calculate richness #
Richness<-SpComp2%>%
  group_by(Site,Plot)%>%
  mutate(plotrich=(count=uniqueN(SpNum)))%>%
  group_by(Site)%>%
  mutate(siterich=(count=uniqueN(SpNum)))%>%
  select(-FunctionalGroup)

# calcualte richness of just graminoids and forbs #
GramForbRichness<-SpComp2%>%
  group_by(Site,Plot)%>%
  filter(FunctionalGroup %in% c("Graminoid","Forb"))%>%
  mutate(gramforbrich=(count=uniqueN(SpNum)))%>%
  group_by(Site)%>%
  mutate(meangramforbrich=(count=uniqueN(SpNum)))%>%
  select(Site,Plot,Site_Plot,gramforbrich,meangramforbrich)
# Remove Duplicate Rows #
GramForbRichness<-subset(GramForbRichness,!duplicated(subset(GramForbRichness,select=c(Site_Plot))))

# calculate relative wiregrass cover #
Wiregrass<-SpComp2%>%
  group_by(Site,Plot)%>%
  mutate(totcov=sum(Cover))%>%
  filter(SpeciesName %in% c("Aristida stricta"))%>%
  mutate(wirecov=sum(Cover))%>%
  mutate(wirerelcov=(wirecov/totcov)*100)%>%
  group_by(Site)%>%
  mutate(meanwirerelcov=mean(wirerelcov))%>%
  select(-SpNum,-Cover,-SpeciesID,-SpeciesName,-totcov,-FunctionalGroup,-Site_Plot)

###################################################################################
 
 ### MERGE THE DATA FILES ###
a<-merge(Richness,SiteData,by=c("Site","Plot"),all=T)
b<-merge(a,Wiregrass,by=c("Site","Plot"),all=T)
c<-merge(b,ShrubData,by=c("Site","Plot"),all=T)
THEDATA<-merge(c,pca,by=c("Site","Plot"),all=T)%>%
  select(Site,Plot,Site_Plot,plotrich,plotanpp,litter,SoilSurfaceLight,UnderstoryLight,GSprecip,
        FRI,TimeSinceLastBurn,Stems,TreeNum,PC1,siterich,wirerelcov,meansiteanpp,meanwirerelcov)

### Remove Duplicate Rows From Merging ###
THEDATA<-subset(THEDATA,!duplicated(subset(THEDATA,select=c(Site_Plot))))
THEDATA<-subset(THEDATA,Site_Plot!="CNF1_20")
### Make NAs a Zero ###
THEDATA[is.na(THEDATA)]<-0
### Keep NAs For Data That Was Not Collected ###
THEDATA[,4:8][THEDATA[,4:8]==0]<-NA
### Retain Zeros For Sites That ANPP Was Collected But There Was Zero ###
THEDATA[THEDATA$Site_Plot=="CFP3_17", "plotanpp"] <- 0
THEDATA[THEDATA$Site_Plot=="WEWO_20", "plotanpp"] <- 0
    #######################

write.csv(THEDATA,"C:/Users/alyoung6/Box/THESIS DATA/Manuscript\\SEMData.csv",row.names = T)

########################################################################################

### FOR THE POLYNOMIAL GRAPHS ###

forbgram<-merge(THEDATA,GramForbRichness,by=c("Site","Plot","Site_Plot"),all=T)%>%
  select(Site,Plot,Site_Plot,wirerelcov,plotanpp,gramforbrich)

### Make NAs a Zero ###
forbgram[is.na(forbgram)]<-0

### Keep NAs For Data That Was Not Collected ###
forbgram[,5][forbgram[,5]==0]<-NA

### Retain Zeros For Sites That ANPP Was Collected But There Was Zero ###
forbgram[forbgram$Site_Plot=="CFP3_17", "plotanpp"] <- 0
forbgram[forbgram$Site_Plot=="WEWO_20", "plotanpp"] <- 0

########################################################################
 
### FIGURE 2 - Natural Variation ###
relwire<-ggplot(THEDATA, aes(x=wirerelcov)) + 
  geom_histogram(color="black",fill="deepskyblue2",binwidth=6)+
  xlab("Plot Relative Wiregrass Cover")+theme_classic(base_size=15)
relwire
freqrich<-ggplot(THEDATA, aes(x=plotrich)) + 
  geom_histogram(color="black",fill="deepskyblue2",binwidth=1)+
  xlab("Plot Richness")+theme_classic(base_size=15)
freqrich
meanrelwire<-ggplot(THEDATA, aes(x=meanwirerelcov)) + 
  geom_histogram(color="black", fill="deepskyblue2",binwidth=10)+
  xlab("Mean Site Relative Wiregrass Cover")+theme_classic(base_size=15)+
  xlim(13,100)
meanrelwire
freqtotrich<-ggplot(THEDATA, aes(x=siterich)) + 
  geom_histogram(color="black", fill="deepskyblue2",binwidth=5)+
  xlab("Total Site Richness")+theme_classic(base_size=15)
freqtotrich
plotanpp<-ggplot(THEDATA,aes(x=plotanpp)) + 
  geom_histogram(color="black", fill="deepskyblue2",binwidth=12)+
  xlab(bquote('Plot Aboveground Biomass '~(g/m^2)))+theme_classic(base_size=15)+
  xlim(0,500)
plotanpp
Meananpp<-ggplot(THEDATA,aes(x=meansiteanpp)) + 
  geom_histogram(color="black", fill="deepskyblue2",binwidth=20)+
  xlab(bquote('Mean Site Aboveground Biomass '~(g/m^2)))+theme_classic(base_size=15)+
  xlim(50,300)
Meananpp

 ### JOIN THE BAR GRAPHS INTO ONE ###
# Figure 2 #
ggarrange(plotanpp,relwire,freqrich,Meananpp,meanrelwire,freqtotrich,
          labels=c("a","b","c","d","e","f"),
          font.label=list(size=20,color="black"),ncol=3,nrow=2)

#####################################################################################

### Polynomial Stats at small scales ###

# for Table S7 #

# Richness ~ Aboveground Biomass #
bt<-unique(forbgram$Site)
for (i in 1:length(bt)){
  subset<-forbgram %>%
    filter(Site==bt[i])
  fit<-lm(gramforbrich~plotanpp+I(plotanpp^2),data=subset)
  print(summary(fit))
}

# Richness ~ Wiregrass #
bt<-unique(forbgram$Site)
for (i in 1:length(bt)){
  subset<-forbgram %>%
    filter(Site==bt[i])
  fit<-lm(gramforbrich~wirerelcov+I(wirerelcov^2),data=subset)
  print(summary(fit))
}

# Aboveground Biomass ~ Wiregrass #
bt<-unique(forbgram$Site)
for (i in 1:length(bt)){
  subset<-forbgram %>%
    filter(Site==bt[i])
  fit<-lm(plotanpp~wirerelcov+I(wirerelcov^2),data=subset)
  print(summary(fit))
}

#################################################################

### FIGURE 4 ###

# Figure 4a - Understory Richness ~ Aboveground Biomass #
polyANPP<-ggplot(forbgram,aes(x=plotanpp,y=gramforbrich))+
  geom_point(position=position_jitter(width=0.05, height=0.05),
  alpha=0.2, color="black")+geom_smooth(method=lm,formula=y~poly(x,2),se=FALSE,
  mapping=aes(group=as.character(Site)),method.args = list(family=gaussian),
  fill=NA, lwd=0.8,alpha=0.1, color="gray67")+
  geom_smooth(method=lm, formula=y~poly(x,2),mapping=aes(group=1),
  method.args=list(family=poisson),color="palevioletred3", lwd=2) +
  theme_bw()+theme(axis.line=element_line(color='black',size=0.1),
          axis.text=element_text(size=18,color="black"),
          axis.title=element_text(size=20),
          axis.ticks.length=unit(0.3,"cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+ 
  ylim(0,15)+
  ylab("Understory Richness") +
  xlab(bquote('Mean Aboveground Biomass '~(g/m^2)))
polyANPP

# mixed effect regression #
forbgram[is.na(forbgram)]<-0
polyanpp<-lmer(gramforbrich~poly(plotanpp,2)+(1|Site),data=forbgram)
summary(polyanpp)
r.squaredGLMM(polyanpp)
  ###

# Figure 4b - Understory Richness ~ Relative Wiregrass Cover #
polyWiregrass<-ggplot(forbgram,aes(x=wirerelcov,y=gramforbrich))+
  geom_point(position=position_jitter(width=0.05, height=0.05),
  alpha=0.2,color="black")+geom_smooth(method=lm,formula=y~poly(x,2),
  se=FALSE,mapping=aes(group=as.character(Site)),
  method.args = list(family=gaussian),fill=NA,lwd=0.8,alpha=0.1,color="gray67")+
  geom_smooth(method=lm, formula=y~poly(x,2),mapping=aes(group=1),
  method.args=list(family=poisson),color="palegreen4",lwd=2)+
  theme_bw()+theme(axis.line=element_line(color='black',size=0.1),
          axis.text=element_text(size=18,color="black"),
          axis.title.x=element_text(size=20),
          axis.title.y=element_text(size=20),
          axis.ticks.length=unit(0.3,"cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+ 
  ylim(0,15)+ 
  ylab("Understory Richness")+
  xlab("Relative Wiregrass Cover")
polyWiregrass

# mixed effect regression #
polywire<-lmer(gramforbrich~poly(wirerelcov,2)+(1|Site),data=forbgram)
summary(polywire)
r.squaredGLMM(polywire)
  ###

# Figure 4c - Aboveground Biomass ~ Relative Wiregrass Cover #
polybiomass<-ggplot(forbgram,aes(x=wirerelcov,y=plotanpp))+
  geom_point(position=position_jitter(width=0.05, height=0.05),
  alpha=0.2, color="black")+geom_smooth(method=lm,formula=y~poly(x,2),se=FALSE,
  mapping=aes(group=as.character(Site)),
  method.args = list(family=gaussian),fill=NA, lwd=0.8,alpha=0.1, color="gray67") +
  geom_smooth(method=lm, formula=y~poly(x,2),mapping=aes(group=1),
  method.args=list(family=poisson),color="skyblue4", lwd=2) +
  theme_bw()+theme(axis.line=element_line(color='black',size=0.1),
          axis.text=element_text(size=18,color="black"),
          axis.title=element_text(size=20),
          axis.ticks.length=unit(0.3,"cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+ 
  ylab(bquote('Mean Aboveground Biomass '~(g/m^2)))+
  xlab("Relative Wiregrass Cover")
polybiomass

# mixed effect regression #
polybio<-lmer(plotanpp~poly(wirerelcov,2)+(1|Site),data=forbgram)
summary(polybio)
r.squaredGLMM(polybio)
  ###

### JOIN THE GRAPHS INTO ONE ###
# Figure 4 #
ggarrange(polyANPP,polyWiregrass,polybiomass,labels=c("a","b","c"),
          font.label=list(size=20,color="black"),ncol=3,nrow=1)
#1600 x 700 

#######################################################

### FIGURE S1: PCA axes 1 and 2 ###
g<-autoplot(PCA,data=Nutrients,scale=0,
         colour="Site",loadings=TRUE,loadings.colour="black",size=3,
         loadings.label=TRUE,loadings.label.colour="black",loadings.label.size=6)

arrow_ends <- layer_data(g, 2)[,c(2,4)]


# if want frame around loadings for each site #
autoplot(PCA,data=Nutrients,scale=0,
         colour="Site",loadings=TRUE,loadings.colour="black",size=3,
         loadings.label=TRUE,loadings.label.colour="black",loadings.label.size=5,
         loadings.label.vjust = 1.5,frame=T,frame.colour = 'Site') +
  geom_point(data = arrow_ends, aes(xend, yend), size = 2) + 
  theme(plot.background=element_blank(),
        panel.background=element_rect(fill='transparent',color='black',size=1),
        legend.key=element_blank())

#1000 x 700


### FIGURES S2-S4: Polynomial Relationships within each site ###

# Figure S2 - Richness ~ Aboveground Biomass #
polyS2<-ggplot(forbgram,aes(x=plotanpp,y=gramforbrich,color=Site))+
  geom_point(position=position_jitter(width=0.05, height=0.05),
        alpha=0.2, color="black")+geom_smooth(method=lm,formula=y~poly(x,2),
        se=FALSE,method.args = list(family=gaussian),fill=NA,lwd=1.5,alpha=0.1,
        linetype="dashed")+theme(strip.text = element_text(size = 12),
        axis.title=element_text(size=20))+
  xlab(bquote('Mean Aboveground Biomass '~(g/m^2)))+
  ylab("Understory Richness")+
  facet_wrap(~Site, scales="free")
polyS2
#1000 x 700


# Figure S3 - Richness ~ Relative Wiregrass Cover #
polyS3<-ggplot(forbgram,aes(x=wirerelcov,y=gramforbrich,color=Site))+
  geom_point(position=position_jitter(width=0.05, height=0.05),
    alpha=0.2, color="black")+geom_smooth(method=lm,formula=y~poly(x,2),
    se=FALSE,method.args = list(family=gaussian),fill=NA, lwd=1.5,alpha=0.1,
    linetype="dashed")+theme(strip.text = element_text(size = 12),
    axis.title=element_text(size=20))+
  ylab("Understory Richness")+
  xlab("Relative Wiregrass Cover")+
  facet_wrap(~Site, scales="free")
polyS3
#1000 x 700

# Figure S4 - Aboveground Biomass ~ Relative Wiregrass Cover #
polyS4<-ggplot(forbgram,aes(x=wirerelcov,y=plotanpp,color=Site))+
  geom_point(position=position_jitter(width=0.05, height=0.05),
    alpha=0.2, color="black")+geom_smooth(method=lm,formula=y~poly(x,2),
    se=FALSE,method.args = list(family=gaussian),fill=NA, lwd=1.5,alpha=0.1,
    linetype="dashed")+theme(strip.text = element_text(size = 12),
    axis.title=element_text(size=20))+
  ylab(bquote('Mean Aboveground Biomass '~(g/m^2)))+
  xlab("Relative Wiregrass Cover")+
  facet_wrap(~Site, scales="free")
polyS4
#1000 x 700
