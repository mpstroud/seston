install.packages("library")
install.packages("streamMetabolizer", dependencies=TRUE, 
                 repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))
update.packages(oldPkgs=c("streamMetabolizer","unitted"), dependencies=TRUE,
                repos=c("https://owi.usgs.gov/R", "https://cran.rstudio.com"))
devtools::find_rtools()
install.packages("Rtools")
library(streamMetabolizer)
install.packages('dplyr', repos="http://cran.cnr.berkeley.edu/")
library(dplyr)
install.packages('tidyr', repos="http://cran.cnr.berkeley.edu/")
library(tidyr)
install.packages('ggplot2', repos="http://cran.cnr.berkeley.edu/")
library(ggplot2)


#import and inspect data
#time should be in Unix units (output from miniDOT)
dat <- read.csv("MC full met input.csv",head=T,sep=",") # load csv data
dim(dat)
dat[c(1,48,96,240,288),] # some example rows
head(dat)
dat[c(1,2,3,4,5,6,7,8,9,10,11,12),]
dat[c(22,2330),]

#convert miniDOT unix time to solar time
dat$solar.time<-as.POSIXct(dat$solar.time, origin="1970-01-01")
dat$solar.time <- lubridate::force_tz(dat$solar.time, "UTC")
head(dat)

#view time series of data
#DO percent saturation and concentration
dat %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(DO.obs='DO concentration\n(mg/L)', DO.sat='Do saturation(%)')
dat %>% unitted::v() %>%
  select(solar.time, DO.obs, DO.sat) %>%
  gather(type, value, DO.obs, DO.sat) %>%
  mutate(type=ordered(type, levels=c('DO.obs','DO.sat')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line(color='black') + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('black')


#depth,water temperature, and light
labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
dat %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

#water temperature and light
labels <- c( temp.water='water temp\n(deg C)', light='PAR\n(umol m-2 s-1)')
dat %>% unitted::v() %>%
  select(solar.time, temp.water, light) %>%
  gather(type, value, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line(color='black') + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('black')

#DO concetration and saturation, water temperature, and light
labels <- c(DO.obs='DO concentration\n(mg/L)', DO.sat='Do saturation(%)', temp.water='water temp\n(deg C)', light='PAR\n(umol m-2 s-1)')
dat %>% unitted::v() %>%
  select(DO.obs,DO.sat,solar.time, temp.water, light) %>%
  gather(type, value, DO.obs,DO.sat,temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('DO.obs','DO.sat','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line(color='black') + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('black')

#########model 1: Bayesian
bayes_name <- mm_name(type='bayes', pool_K600='none', err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_name
bayes_specs <- specs(bayes_name)
bayes_specs
#specifications of model (i.e. priors)
bayes_specs <- specs(bayes_name, burnin_steps=100, saved_steps=200, n_cores=1,K600_daily_meanlog = 10)
#run model
mm <- metab(bayes_specs, data=dat)
mm
#model parameters
output<-predict_metab(mm)
output
output[c(11:25),]
plot_metab_preds(mm)
#check for collinearity between ER and K600
out<-get_params(mm)
par(mar = c(5,5,1,1))
plot<-plot(output2$K600.daily,output2$ER.daily,xlab =expression("Daily"~K[600]~(d^{-1})),
           pch=16,cex=1.5,col='black',
           ylab=expression("Daily ER" ~ (g ~ O[2]~m^{-2} ~d^{-1})),
           cex.lab=1.2)
output




######model 2:Maximum likelihood estimation
mle_name<- mm_name(type='mle')
mle_name
mle_specs <- specs(mle_name)
mle_specs
mle_specs <- specs(mle_name, init.GPP.daily=15, init.ER.daily=-12, init.K600.daily=3)
mle_specs <- specs(mle_name, day_start=5, day_end=29, init.K600.daily=100,init.GPP.daily=15)
mle_specs
mm <- metab(mle_specs, data=dat, info=c(site='Miller Creek', source='Kim Bray'))
mm <- metab(mle_specs, data=dat)
#extract specific pieces of info
get_info(mm)
head(get_data(mm))
get_fitting_time(mm) #time for model to run
get_specs(mm) #model specs

#predict and plot metabolism estimates from MLE model
predict_metab(mm)
plot_metab_preds(mm)
#predict and plot DO predictions from MLE model AND observed values
plot_DO_preds(mm)
head(predict_DO(mm))


######### Model 3:Classic: linear GPP, constant ER (also the default)
mm_classic <- 
  mm_name('mle', GPP_fun='linlight', ER_fun='constant') %>% 
  specs() %>%
  metab(dat)
mm_classic

#inspect fitted daily parameters
get_params(mm)
#DO predictions from model
predict_DO(mm) %>% head()
plot_DO_preds(mm)
mcmc <- get_mcmc(mm)

########Model 4:The Saturator (GPP saturating with light, constant ER)
mm_saturator <- 
  mm_name('mle', GPP_fun='satlight', ER_fun='constant') %>% 
  specs() %>%
  metab(dat)
mm_saturator

get_params(mm_saturator) %>% select(date, warnings, errors)
predict_metab(mm_saturator) %>% select(date, warnings, errors)
predict_DO(mm_saturator) %>% head
plot_DO_preds(mm_saturator)

mm_saturator2 <- 
  mm_name('mle', GPP_fun='satlight', ER_fun='constant') %>% 
  specs() %>%
  metab(dat, data_daily=select(get_params(mm_saturator), date, init.Pmax=Pmax, init.alpha=alpha))
get_params(mm_saturator2)
