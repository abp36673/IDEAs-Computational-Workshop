#Alexandria Purcell
#May 16 2019
#Estimation 

load('data.RData')     #load the data and plot flu cases
plot(flu,type='b',log='y',main='Epidemic in a British boarding school', cex.main=0.85,
     xlab='Day', ylab='Active influenza cases')

##Exercise 1: this equation shows the important one to one relationship
  #between R0 and the final epidemic size. Plot the relationship between 
  #the total epidemic size andR0 for the complete range of values between 0 and 1. 

epi.size<- seq(0.001,0.999, by= 0.001)

R0<-c()
for(ZN in epi.size){
  R0<-c(R0, log10(1-ZN)/-ZN)
  
}

plot(epi.size,R0)

model<-lm(log(flu[1:4])~day[1:4],data=flu);  #fit a linear model
summary(model)         #summary statistics for fit model
slope<-coef(model)[2]  #extract slope parameter
slope                 #print to screen

##EXercise 2: Supposing admission happened within 24 hours of the onset
  # of symptoms, how does this affect our estimate of R0? Twelve hours? 
gamma24=1
gamma12=1/0.5
R024<-(slope/gamma24)+1
R012<-(slope/gamma12)+1

#As time of admission is decreased from 2.5 days to 1 day the R0 drops from 
  #3.7 to 2.09. When admission time drops to 12 hours R0= 1.55

##Exercise 3: use this method to obtain estimates of R0 for measles 
  #from the first community assuming that the infectious period
  #is approximately two weeks or 0.0384 years

plot(niamey$V1,type='b',log='y',main='Measles Cases in Niamey V1', cex.main=0.85,
     xlab='Biweek', ylab='Active measles cases')



weeks<-seq(1,2*nrow(niamey),by=2)
niamey$weeks<-weeks

biweek<-seq(1,17,2)
correlation<-c()
for(t in biweek){
  model<-lm(log(niamey$V1[1:t])~biweek[1:t],data=niamey);  #fit a linear model
  print(summary(model))
}
# the going out to row 6 has the best R value of 0.951. 

model<-lm(log(niamey$V1[1:6])~biweek[1:6],data=niamey);  #fit a linear model
summary(model)
slope<-coef(model)[2]  #extract slope parameter
slope 

R0V1=(slope/0.5)+1
#R0V1=1.43


##Exercise 4: Plot the estimates of R0 obtained from n=3,4,5 data points 
  #against the standard error of the slope from the regrassion analysis
  #to show this tradeoff 

summary(model)$coefficients[2,2]

summary(model)[['coefficients']][,2]

biweek<-seq(1,31,2)
alues<-seq(3,9,1)
stan.error<-c()
rvals<-c()
for(x in alues){
  nimodel<-lm(log(V1[1:x])~biweek[1:x],data=niamey);  #fit a linear model
  SE<-coef(summary(nimodel))[2,2]
  slope<-coef(summary(nimodel))[1,2]
  RN<-slope/0.5+1
  stan.error<-c(stan.error,slope)
  rvals<-c(rvals,RN)
}
plot(rvals,stan.error)

load('data.RData')
niamey[5,3]<-0  #replace a "NA"
niamey<-data.frame(biweek=rep(seq(1,16),3),site=c(rep(1,16),rep(2,16),rep(3,16)),
                   cases=c(niamey[,1],niamey[,2],niamey[,3]))

plot(niamey$biweek,niamey$cases,type='p',col=niamey$site,xlab='Biweek',ylab='Cases')
lines(niamey$biweek[niamey$site==1],niamey$cases[niamey$site==1])
lines(niamey$biweek[niamey$site==2],niamey$cases[niamey$site==2],col=2)
lines(niamey$biweek[niamey$site==3],niamey$cases[niamey$site==3],col=3)

closed.sir.model <- function (t, x, params) {  #SIR model equations
  S <- x[1]
  I <- x[2]
  beta <- params
  dS <- -beta*S*I
  dI <- beta*S*I-(365/13)*I
  list(c(dS,dI))
}

sse.sir <- function(params0,data,site){  #function to calculate squared errors
  data<-data[data$site==site,]    #working dataset, based on site
  t <- data[,1]*gamma            #time in biweeks
  cases <- data[,3]               #number of cases
  beta <- exp(params0[1])            #parameter beta
  S0 <- exp(params0[2])           #initial susceptibles
  I0 <- exp(params0[3])           #initial infected        
  out <- as.data.frame(ode(c(S=S0,I=I0),times=t,closed.sir.model,beta,hmax=1/120))
  sse<-sum((out$I-cases)^2)       #sum of squared errors
}


library(deSolve)   #differential equation library
params0<-c(-3.2,7.3,-2.6)  #initial guess

fit1 <- optim(params0,sse.sir,data=niamey,site=1) #fit
exp(fit1$par)  #back-transform parameters
fit2 <- optim(params0,sse.sir,data=niamey,site=2) #fit
exp(fit2$par)  #back-transform parameters
fit3 <- optim(params0,sse.sir,data=niamey,site=3) #fit
exp(fit3$par)

par(mfrow=c(3,1))   #set up plotting area for multiple panels
plot(cases~biweek,data=subset(niamey,site==1),type='b',col='blue', pch=21) #plot site 1
t <- subset(niamey,site==1)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit1$par[2]),I=exp(fit1$par[3])),times=t,
                            closed.sir.model,exp(fit1$par[1]),hmax=1/120))
#obtain model predictions
lines(mod.pred$I~subset(niamey,site==1)[,1]) #and plot as a line

plot(cases~biweek,data=subset(niamey,site==2),type='b',col=site) #site 2
t <- subset(niamey,site==2)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit2$par[2]),I=exp(fit2$par[3])),times=t,
                            closed.sir.model,exp(fit2$par[1]),hmax=1/120))
lines(mod.pred$I~subset(niamey,site==2)[,1])


plot(cases~biweek,data=subset(niamey,site==3),type='b',col=site) #site 3
t <- subset(niamey,site==3)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit3$par[2]),I=exp(fit3$par[3])),times=t,
                            closed.sir.model,exp(fit3$par[1]),hmax=1/120))
lines(mod.pred$I~subset(niamey,site==3)[,1])

##Exercise 5: Modify the code aboe to estimate y and B simultaneously

closed.sir.model <- function (t, x, params) {  #SIR model equations
  S <- x[1]
  I <- x[2]
  beta <- params[1]
  gamma<-params[2]
  dS <- -beta*S*I
  dI <- beta*S*I-gamma*I
  list(c(dS,dI))
}

sse.sir <- function(params0,data,site){  #function to calculate squared errors
  data<-data[data$site==site,]    #working dataset, based on site
  t <- data[,1]*gamma          #time in biweeks
  cases <- data[,3]               #number of cases
  beta <- exp(params0[1])            #parameter beta
  S0 <- exp(params0[2])           #initial susceptibles
  I0 <- exp(params0[3]) #initial infected 
  gamma<-exp(params0[4])
  out <- as.data.frame(ode(c(S=S0,I=I0),times=t,closed.sir.model,parms=c(beta=beta, gamma=gamma), hmax=1/120))
  sse<-sum((out$I-cases)^2)       #sum of squared errors
}

library(deSolve)   #differential equation library
params0<-c(-3.2,7.3,-2.6,log(13/365))  #initial guess
lsoda()
fit1 <- optim(params0,sse.sir,data=niamey,site=1) #fit
exp(fit1$par)  #back-transform parameters
fit2 <- optim(params0,sse.sir,data=niamey,site=2) #fit
exp(fit2$par)  #back-transform parameters
fit3 <- optim(params0,sse.sir,data=niamey,site=3) #fit
exp(fit3$par)


par(mfrow=c(2,2))   #set up plotting area for multiple panels
plot(cases~biweek,data=subset(niamey,site==1),type='p', pch=21) #plot site 1
t <- subset(niamey,site==1)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit1$par[2]),I=exp(fit1$par[3])),times=t,
                            closed.sir.model,exp(fit1$par[1]),hmax=1/120))
#obtain model predictions
lines(mod.pred$I~subset(niamey,site==1)[,1]) #and plot as a line

plot(cases~biweek,data=subset(niamey,site==2),type='b',col=site) #site 2
t <- subset(niamey,site==2)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit2$par[2]),I=exp(fit2$par[3])),times=t,
                            closed.sir.model,exp(fit2$par[1]),hmax=1/120))
lines(mod.pred$I~subset(niamey,site==2)[,1])


plot(cases~biweek,data=subset(niamey,site==3),type='b',col=site) #site 3
t <- subset(niamey,site==3)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit3$par[2]),I=exp(fit3$par[3])),times=t,
                            closed.sir.model,exp(fit3$par[1]),hmax=1/120))
lines(mod.pred$I~subset(niamey,site==3)[,1])