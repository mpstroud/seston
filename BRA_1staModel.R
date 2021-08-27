#install and load packages
library(rstan)
library(shinystan)
install.packages('dplyr', repos="http://cran.cnr.berkeley.edu/")
library(dplyr)
install.packages('tidyr', repos="http://cran.cnr.berkeley.edu/")
library(tidyr)
library(gdata)
library(bayesplot)



###################model 1:Bormann and Webster 1999##############

options(mc.cores = parallel::detectCores())
#this code allows stan to count the number of cores and run 4 chains on each core. 
##load data

data<-read.csv("/Users/mpeipoch/Dropbox/BRA_BWmodel.csv")
head(data)
datain <- list(N = length(data$conc), 
               conc =data$conc, 
               H = data$H, 
               t= data$t)
datain

write("//Stan model for a simple exponential model of dChl/dt based on Borman an Webster 1999
      
      data{
      int < lower = 1 > N; // Sample size
      vector[N] H; // Predictor
      vector[N] conc; // Outcome
}

parameters {
 real alpha; // Intercept
 real beta; // Slope (regression coefficients)
 real < lower = 0 > sigma; // Error SD
}

model {
 conc ~ normal(alpha + H * beta , sigma);
}

generated quantities {
} // The posterior predictive distribution",

"BWmodel.stan")

stanc("BWmodel.stan") #check if the model is written
stan_model1 <- "BWmodel.stan" #save the file path for the model

fit <- stan(file = stan_model1, data = datain, warmup = 500, iter = 1000, chains = 4, cores = 2, thin = 1)


fit


#We can also look at the full posterior of our parameters by extracting them from the model object. 

posterior <- extract(fit)
str(posterior)






















##model 1
sink("BWmodel.stan")

cat("
    
    data {
    int <lower = 1> N;
    vector[N] H;//channel depth  //channel depth is a vector of length N
    vector[N] conc;//chl conc. //nitrate conc. is a vector of length N
    }
    
    parameters {
     real intercept;                // intercept aka. initial conc
     real U;                // growth
     real W;                // deposition velocity
     real<lower = 0> sigma; //  standard deviation
    }
    
    model {
    for (i in 1:N){
    conc[i] ~ normal(intercept*(U-(W/H[i])), sigma); // likelihood
    }
    U~normal(0,10); //priors
    W~normal(0,5);
    intercept~normal(0,0.29);
    sigma~normal(0,1);
    
    }
    "
    ,fill=TRUE)
sink()

#run the MCMC

fit <- stan("BWmodel.stan", data = datain,  iter = 1000, chains = 4,control = list(max_treedepth = 15))
load("fit.RData")
print(fit)
traceplot(fit,pars="U")
traceplot(fit)



########################model 2:hierarchical daily nitrate uptake#########################
data2<-read.csv("eco models 2.csv")
head(data2)
datain2 <- list(N = length(data2$conc), 
                conc =data2$conc, 
                Q = data2$Q, 
                A = data2$A, 
                t= data2$t,
                day=data2$day,
                nday=27)
datain2


sink("nmodel2.stan")

cat("
    
    data {
    int <lower = 1> N;          //sample size is an integer with a lower bound of 1
    int <lower = 1> nday;       //# days is an integer with a lower bound of 1
    int day[N];                // day is a vector of length N
    vector[N] Q;//discharge     //discharge is a vector of length N
    vector[N] A;//channel area  //channel area is a vector of length N
    vector[N] conc;//nitrate conc. //nitrate conc. is a vector of length N
    }
    
    parameters {
    real<lower = 0> intercept;         // predawn max peak is a vector of length nday
    vector[nday] U;                // daily uptake is a vector of length nday
    real<lower = 0> sigma;        //standard deviation is a real number with a lower bound of zero
    real<lower = 0> Umean;  //mean of all daily uptake fluxes is a real number with a lower bound of zero  
    real<lower = 0> Usd;    // sd of all daily uptake fluxes is a real number with a lower bound of zero
   }
    
    model {
    for (i in 1:N){
    conc[i]~ normal(intercept-((A[i]*U[day[i]])/Q[i]), sigma); // likelihood
    }
    for (j in 1:nday){
    U[j]~normal (Umean,Usd);
    }
    U~normal(0,150); //priors
    intercept~normal(0,0.29);
    sigma~normal(0,50);
    Umean~normal(0,150);
    Usd~normal(0,50);

    }
    "
    ,fill=TRUE)
sink()

#run the MCMC

fit2 <- stan("nmodel2.stan", data = datain2,  iter = 1000, chains = 4,control = list(max_treedepth = 15))
load("fit.RData")
print(fit2)
traceplot(fit2,pars="U")
traceplot(fit2,pars="intercept")

traceplot(fit2,pars="Umean")
traceplot(fit2,pars="Usd")



###model 2
data2<-read.csv("eco models 2.csv")
head(data2)
datain2 <- list(N = length(data2$conc), 
                conc =data2$conc, 
                Q = data2$Q, 
                A = data2$A, 
                t= data2$t,
                day=data2$day,
                ndays=data2$ndays)
datain2




#############################model 3: Daily nitrate uptake with time series#################

sink("model2.stan")

cat("
    
    data {
    int <lower = 1> N;      //N (sample size) is an integer with a lower bound of 1
    int <lower = 1> day;  //day (# of days) is an integer with a lower bound of 1
    vector[N] Q;           //Q (discharge) is a vector of length N
    vector[N] A;          //A (channel area) is a vector of length N
    vector[N] conc;       // conc (of nitrate) is a vector of length N
    }
    
    parameters {
    real intercept;     // intercept in one day (aka max predawn peak of nitrate) is a real number
    real<lower = 0> sigma; //  standard deviation of uptake within a day is an integer with a lower bound of zero
    vector[day] U;             //U (daily nitrate uptake) is a vector of length day
    vector [day] theta;     // theta is a vector of length day
    vector [day] w;         // w (random variable) is a vector of length day
    real b;                 // b (a coefficient) is a real number
    real<lower = 0> sigma_day; // S.D. of daily uptake over month is an integer with a lower bound of zero
    }
    
     model {
    for (i in 1:N){
    conc[i]~ normal(intercept-((A[i]*U[day])/Q[i]), sigma); // likelihood for daily uptake
    }
    for (j in 1:day){
    U[j]~ normal(theta*U[j-1]+b+w[j],sigma_day);//likelihood for daily change
    }

    U~normal(0,150); //priors
    intercept~normal(0,0.29);
    sigma~normal(0,1);
    }
    "
    ,fill=TRUE)
sink()

#run the MCMC

fit2 <- stan("model2.stan", data = datain2,  iter = 1000, chains = 4,control = list(max_treedepth = 15))

print(fit2)
traceplot(fit2,pars="U")
traceplot(fit2)





###model 3: daily uptake + time series assessment of change in daily uptake
data2<-read.csv("eco models 2.csv")
head(data2)
datain2 <- list(N = length(data2$conc), 
               conc =data2$conc, 
               Q = data2$Q, 
               A = data2$A, 
               t= data2$t,
               day=data2$day,
               day=data2$day)
datain2


sink("model3.stan")

cat("
    
    data {
    int <lower = 1> N;
    int <lower = 1> day;
    vector[N] Q;
    vector[N] A;
    vector[N] conc;
    vector[N] t;
    }
    
    parameters {
    real intercept;                // intercept
    real<lower = 0> sigma; //  standard deviation
    real int_day;
    real theta;
    real<lower = 0> sigma_day;
    vector[day] U;
    vector[day] y;
    vector[day] w;
     }
    
    model {
    for (i in 1:N){
    deltaconc[i]~ normal(intercept-((A[i]*U[day[i]])/Q[i]), sigma); // likelihood
    }
    for (j in 1:day){
    deltaconc[j]~ normal(theta*y[j-1]+theta*y[j-2]+w[j]+int_day,sigma_day);
    }
    U~normal(0,150); //priors
    intercept~normal(0,0.29);
    sigma~normal(0,1);
    theta~normal(0,1);
    int_day~normal(0,0.30);
    sigma_day~normal(0,1);
    w~normal(0,1);
    }
    "
    ,fill=TRUE)
sink()

#run the MCMC

fit3 <- stan("model3.stan", data = datain2,  iter = 1000, chains = 4,control = list(max_treedepth = 15))

print(fit3)
traceplot(fit3,pars="U")
traceplot(fit3)



###################predictors of nitrate uptake######################################3

daily<-read.csv("daily.csv")
head(daily)
daily <- list(ER =daily$ER, 
                GPP = daily$GPP, 
                nitrate = daily$nitrate, 
                day= daily$day,
                N=length(daily$N))
daily


sink("npred.stan")
cat("
    
    data {
    int <lower = 1> N;
    vector[N]GPP;
    vector[N] nitrate;
    }
    
    parameters {
    real b;                // intercept
    real m;                // slope
    real<lower = 0> sigma; //  standard deviation
    }
    
    model {
    nitrate ~ normal(b + m *GPP , sigma); // likelihood
    
    b~normal(0,10); //priors
    m~normal(0,100);
    
    }
    "
    ,fill=TRUE)
sink()

#run the MCMC
npredfit <- stan("npred.stan", data = daily,  iter = 1000, chains = 4)#chains always equal to 4
#with stan, 1000 iterations usually is enough
outputp<-print(npredfit)
traceplot(npredfit)


save(npredfit, file = "npredfit.RData")

npredfit2<-load("npredfit.RData")
fit_extract<-extract(npredfit)
fit_extract
#extract and observe the posterior distributions
plot(fit_extract$m, fit_extract$b, pch=16, col= 'black', cex=0.9, xlab="slope", ylab="intercept")
#plot model fit
par(mar = c(5,4,1,1))
par(mar = c(5,5,1,1))
plot(daily$GPP,daily$nitrate, xlab = expression("Daily GPP" ~ (g ~ O[2]~m^{-2} ~d^{-1})),ylab=expression("Daily nitrate uptake" ~ (mg ~ N~m^{-2} ~d^{-1})),cex.axis=1.5, cex.lab=1.7,las=1,pch=16,cex=2,col='black')
lines(daily$GPP, 38.24+ 163.28*daily$GPP,col="black",pch=18, lwd=3)




#plot the credible interval by extracting predicted slope and intercepts from 1000 iterations
x<-daily$GPP
y<-daily$nitrate
plot(x,y, pch= 16, col='blue')
lines(x, 4.07+ 1.03*x)

for (i in 1:1000){
  lines(x, fit_extract$m[i] + fit_extract$b[i]*x, col='light gray')
}
points(x,y, pch= 16, col='blue')



linearmod<-lm(daily$ER~daily$nitrate)
summary(linearmod)

abline(linearmod)
