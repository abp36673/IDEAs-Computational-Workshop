#Alexandria Purcell
#May 16 2019
# Pulse vaccination 

require(deSolve)


sir.model.open <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
    as.list(params),                   #this argument to "with" lets us use the variable names
    {                                  #the system of rate equations
      dS <- mu*(S+I+R) - beta*S*I - mu*S
      dI <- beta*S*I - gamma*I - mu*I
      dR <- gamma*I - mu*R
      dx <- c(dS,dI,dR)                #combine results into a single vector dx
      list(dx)                         #return result as a list
    }
  )
}

R0 <- 10
N <-  1					                       #population size
mu <- 0.02                             #per capita birth/death rate
gamma <- 365/10  			                 #recovery rate (in years)
beta <- R0*(gamma+mu)/N 	             #transmission rate
xstart <- c(S=0.2, I=0.001, R=1-0.2-0.001)	 #initial conditions, must sum to one
Tmax <- 120                            #integrate for 200 years after transients
params <- c(beta=beta, gamma=gamma, mu=mu)           #parameter vector
tau <- 0.1                             #size of time step
times <- seq(0, Tmax, by=tau)

out <- ode(xstart,times,sir.model.open,params, method='ode45', rtol=1e-7)

op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))                      #set graphical parameters
plot(I~time,data=out, type='l', lwd=2)                          #plot the I variable against time
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                      #re-set graphical parameters
plot(I~S,data=out,log='xy',yaxt='n',xlab='S', type='l', lwd=2)  #plot phase portrait
par(op) 

xstart<-out[which(out[,1]==50),2:4]

#adding pulse vaccinations:defining additional parameters

pv<-0.1
Tv<-4
vacc.events<-floor(Tmax/Tv)

data<-data.frame(S=out[which(out[,1]==50),2],
                 I=out[which(out[,1]==50),3],
                 R=out[which(out[,1]==50),4])
for(i in 1:vacc.events){
  out<-ode(xstart,seq(tau,Tv,by=tau),sir.model.open,params,method = 'ode45', rtol=1e-7)
  xstart<-out[dim(out)[1],2:4]
  xstart[1]<-(1-pv)*(tail(out,1)[2])
  xstart[3]<-xstart[3]+(pv)*(tail(out,1)[2])
  data<-rbind(data,out[,2:4])
}

data$time <- seq(50, Tmax+50, by=tau)
par(mar=c(5,4,4,4)+0.1)
plot(data$time[1:500], data$I[1:500], type='l', xlab='Time', ylab='', col='red', axes=FALSE)
axis(2, col.axis='red')
mtext(side=2, line=2.5, 'Infected', col='red')
box()
axis(1)
par(new=TRUE)
plot(data$time[1:500], data$S[1:500], type='l', xlab='', ylab='', axes=FALSE, col='black')
axis(4)
mtext(side=4, line=2.5, 'Susceptibles')

#Exercise 1: why might this periodicity be of interest? 
  #what might the consequences of "peaks" and troughs for public health 

#the periodicity is of interest because the fluctuations in the susceptibles over time is representing the replenishment 
  #of the susceptible population through births and this can be used to determine where a pulse should occur
#the peaks and troughs are important for public health because you can follow the fluctuations to determine when would be the 
  #the best times to pulse vaccinate. 

#Exeercise 2: by modifying the value of pv, see if you can locate the vaccination
  #threshold. does this agree with the analytically predicted threshold 

e<-exp(1)
value<- ((mu*Tv-pv)*(e^(mu*Tv)-1)+(mu*pv*Tv))/((mu*Tv)*(pv-1+e^(mu*Tv)))

#The analytically predicted threshold would be 1/R0 which is 0.1 but the current threshold that we 
  # have is 0.4 which is greater than 0.1. 
#To reach the predicted threshold of 0.1 you would want to have a pv value of at least 0.54 to bring you to a value of 0.09
  # any value above 0.54 would give you a threshold less than 0.1 

ex2<-function(pv){

sir.model.open <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
    as.list(params),                   #this argument to "with" lets us use the variable names
    {                                  #the system of rate equations
      dS <- mu*(S+I+R) - beta*S*I - mu*S
      dI <- beta*S*I - gamma*I - mu*I
      dR <- gamma*I - mu*R
      dx <- c(dS,dI,dR)                #combine results into a single vector dx
      list(dx)                         #return result as a list
    }
  )
}

R0 <- 10
N <-  1					                       #population size
mu <- 0.02                             #per capita birth/death rate
gamma <- 365/10  			                 #recovery rate (in years)
beta <- R0*(gamma+mu)/N 	             #transmission rate
xstart <- c(S=0.2, I=0.001, R=1-0.2-0.001)	 #initial conditions, must sum to one
Tmax <- 120                            #integrate for 200 years after transients
params <- c(beta=beta, gamma=gamma, mu=mu)           #parameter vector
tau <- 0.1                             #size of time step
times <- seq(0, Tmax, by=tau)

out <- ode(xstart,times,sir.model.open,params, method='ode45', rtol=1e-7)

#op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))                      #set graphical parameters
#plot(I~time,data=out, type='l', lwd=2)                          #plot the I variable against time
#par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                      #re-set graphical parameters
#plot(I~S,data=out,log='xy',yaxt='n',xlab='S', type='l', lwd=2)  #plot phase portrait
#par(op) 

xstart<-out[which(out[,1]==50),2:4]

#adding pulse vaccinations:defining additional parameters

#pv<-0.1
Tv<-4
vacc.events<-floor(Tmax/Tv)

data<-data.frame(S=out[which(out[,1]==50),2],
                 I=out[which(out[,1]==50),3],
                 R=out[which(out[,1]==50),4])
for(i in 1:vacc.events){
  out<-ode(xstart,seq(tau,Tv,by=tau),sir.model.open,params,method = 'ode45', rtol=1e-7)
  xstart<-out[dim(out)[1],2:4]
  xstart[1]<-(1-pv)*(tail(out,1)[2])
  xstart[3]<-xstart[3]+(pv)*(tail(out,1)[2])
  data<-rbind(data,out[,2:4])
}

data$time <- seq(50, Tmax+50, by=tau)
par(mar=c(5,4,4,4)+0.1)
plot(data$time[1:500], data$I[1:500], type='l', xlab='Time', ylab='', col='red', axes=FALSE)
axis(2, col.axis='red')
mtext(side=2, line=2.5, 'Infected', col='red')
box()
axis(1)
par(new=TRUE)
plot(data$time[1:500], data$S[1:500], type='l', xlab='', ylab='', axes=FALSE, col='black')
axis(4)
mtext(side=4, line=2.5, 'Susceptibles')
}

for(i in seq(0.1,0.95,0.5)){
  ex2(i)
}
#when run in R markdown the varying pv values are evaluated and plotted 

#Exercise 3: One thing we might be interested in knowing is how our vaccination 
#strategy affects mean disease prevalence. Devise a strategy to study this 
#question and produce a plot illustrating the result

ex3<-function(pv){
  
  sir.model.open <- function (t, x, params) {    #here we begin a function with three arguments
    S <- x[1]                               #create local variable S, the first element of x
    I <- x[2]                               #create local variable I
    R <- x[3]                               #create local variable R
    with(                                   #we can simplify code using "with"
      as.list(params),                   #this argument to "with" lets us use the variable names
      {                                  #the system of rate equations
        dS <- mu*(S+I+R) - beta*S*I - mu*S
        dI <- beta*S*I - gamma*I - mu*I
        dR <- gamma*I - mu*R
        dx <- c(dS,dI,dR)                #combine results into a single vector dx
        list(dx)                         #return result as a list
      }
    )
  }
  
  R0 <- 10
  N <-  1					                       #population size
  mu <- 0.02                             #per capita birth/death rate
  gamma <- 365/10  			                 #recovery rate (in years)
  beta <- R0*(gamma+mu)/N 	             #transmission rate
  xstart <- c(S=0.2, I=0.001, R=1-0.2-0.001)	 #initial conditions, must sum to one
  Tmax <- 120                            #integrate for 200 years after transients
  params <- c(beta=beta, gamma=gamma, mu=mu)           #parameter vector
  tau <- 0.1                             #size of time step
  times <- seq(0, Tmax, by=tau)

  out <- ode(xstart,times,sir.model.open,params, method='ode45', rtol=1e-7)
  xstart<-out[which(out[,1]==50),2:4]

  pv<-0.1
  Tv<-4
  vacc.events<-floor(Tmax/Tv)
  
  data<-data.frame(S=out[which(out[,1]==50),2],
                   I=out[which(out[,1]==50),3],
                   R=out[which(out[,1]==50),4]) 
  
  
for(i in 1:vacc.events){
  out<-ode(xstart,seq(tau,Tv,by=tau),sir.model.open,params,method = 'ode45', rtol=1e-7)
  xstart<-out[dim(out)[1],2:4]
  xstart[1]<-(1-pv)*(tail(out,1)[2])
  xstart[3]<-xstart[3]+(pv)*(tail(out,1)[2])
  data<-rbind(data,out[,2:4])
meanD<-mean(data$I)
list<-(list,meanD) 
}
}



