#Stream metabolizer applied to WCC40 SCAN sensor DO data from 2018 and 2019

install_github('streampulse/StreamPULSE', dependencies=TRUE)
install.packages('streamMetabolizer', dependencies=TRUE,
                 repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))

#Get ready to run StreamMetabolizer
library(streamMetabolizer)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
library(rstan)
library(readr)
library(unitted)
library(lubridate)
library(data.table)

setwd("S:/DOviedo/Manuscripts/Metabolism")

#Compiling DO and Temp data from scan files
DO =list.files(pattern="*.csv", full.names=TRUE)
DO = ldply(DO, read_csv)
DOplot<-ggplot(DO, aes(x=DateTime, y=DOOxi_cln)) + geom_point(size =1)
DOTEMP<-DO%>%select(DateTime,DOOxi_cln,TempOxiC_cln)

names(DOTEMP)[1] <- "date"
names(DOTEMP)[2] <-"DO.obs"
names(DOTEMP)[3] <-"temp.water"

DOTEMP$date<-(as.POSIXct(DOTEMP$dateEST, format="%m/%d/%Y %H:%M", tz = "EST"))

#depth light and pressure data

PQPAR<-read.csv("flowlightpress.csv", header = TRUE)

PQPAR$date<-(as.POSIXct(PQPAR$date, format="%m/%d/%Y %H:%M", tz = "EST"))


#Merging DLP and DO and T

data<-merge(x=DOTEMP,y=PQPAR,by="date", all.y = TRUE)


#to calculate DOsat


data$pressure<-data$pressure*1000*0.09804
names(data)[5] <- "pressure.air"
data$DO.sat<-calc_DO_sat(temp.water=data$temp.water, pressure.air=data$pressure.air, salinity.water = 0, model = "garcia-benson",)

data$X<-NULL

#to calculate depth from discharge

data$discharge<-data$discharge*0.028316846592 
data$depth<-calc_depth(Q = u(data$discharge, "m^3 s^-1"), c = u(0.409, "m"), f = u(0.294, ""))


#Prepare data set

mydata<-data
posix.time.localtz <- as.POSIXct(mydata$date, format="%m-%d-%Y %H:%M:%S", tz='Etc/GMT+5') # creates a vector with date as Posixct
lubridate::tz(posix.time.localtz) #confirm that time zone is Etc/GMT+5"
posix.time.solar <- streamMetabolizer::calc_solar_time(posix.time.localtz, longitude=-75.783109) #longitude at WCC, converts UTC to solar
mydata$solar.time<-posix.time.solar # adds solar.time column to dataset
lubridate::tz(posix.time.solar)#confirm that time zone is UTC

# to remove extra columms and NA lines
mydata$date<-NULL
mydata$pressure.air<-NULL
mydatacomplete<-mydata[complete.cases(mydata),]

#smoothing mydata to mydata.s

mydatacomplete$index <- 1:nrow(mydatacomplete)  # create index variable
mydatasmoothing<-loess(DO.obs~index, data=mydatacomplete, span=0.01)
mydata.s.pred<-predict(mydatasmoothing)
mydata.s<-mydatacomplete
mydata.s$DO.obs<-mydata.s.pred
mydata.s$index<-NULL

#modeling mydata

plot(x=mydatacomplete$solar.time,y =mydatacomplete$DO.obs)
mydatanoq<-mydata
mydatanoq$discharge<-NULL
mydatanoq<-mydata[complete.cases(mydatanoq),]


bayes_mydatanoq <- mm_name(type='bayes', pool_K600='none', err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_specs <- specs(bayes_mydatanoq)

mm_bayes_mydatanoq <- metab(bayes_specs, data=mydatanoq)
plot_DO_preds(mm_bayes_mydata)

#mle model for mydata

mle_WCC17 <- mm_name(type='mle')
mle_specs <- specs(mle_WCC17)
mm_mle_mydata.s <- metab(mle_specs, data=mydata.s, info=c(site='WCC', source='DOV'))

#examine results
plot_DO_preds(mm_mle_mydata.s)
plot_metab_preds(mm_mle_mydata.s)
mle_mydata.s_daily<-get_params(mm_mle_mydata.s)



#subset for testing

test<-subset(data, date>"2018-08-27"&date<"2018-09-05")
posix.time.localtz <- as.POSIXct(test$date, format="%m-%d-%Y %H:%M:%S", tz='Etc/GMT+5') # creates a vector with date as Posixct
lubridate::tz(posix.time.localtz) #confirm that time zone is Etc/GMT+5"
posix.time.solar <- streamMetabolizer::calc_solar_time(posix.time.localtz, longitude=-75.783109) #longitude at WCC, converts UTC to solar
test$solar.time<-posix.time.solar # adds solar.time column to dataset
lubridate::tz(posix.time.solar)#confirm that time zone is UTC


#bayesian modeling

test$date<-NULL
test$pressure.air<-NULL

plot(x=test$solar.time,y =test$temp.water)

bayes_test <- mm_name(type='bayes', pool_K600='normal', err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_specs <- specs(bayes_test)

mm_bayes_test <- metab(bayes_specs, data=test)
plot_DO_preds(mm_bayes_set1)

get_params(mm2,  date_start = NA,  date_end = NA,  uncertainty = "ci",  messages = TRUE,  attach.units = TRUE)


#extras

#ts <- seq.POSIXt(as.POSIXct("2018-01-01 00:00:00","%Y-%m-%d %H:%M:%S", tz = "EST"), as.POSIXct("2019-12-31 23:55:00","%Y-%m-%d %H:%M:%S", tz = "EST"), by="5 min", tz = "EST")
#ts <- seq.POSIXt(as.POSIXlt("2018-01-01 00:00:00", tz = "EST"), as.POSIXlt("2019-12-31 23:55:00", tz = "EST"), by="5 min", tz = "EST")
#ts <- format.POSIXct(ts,"%Y-%m-%d %H:%M:%S")
#df <- as.data.frame(ts, is.character(ts))
#names(df)[1] <-"date"
#DOTEMPwithmissing <- full_join(df,DOTEMP)



mydatacomplete$ydate=date(mydatacomplete$solar.time)
test=data.table(mydatacomplete)
test[,mean(discharge),date]

