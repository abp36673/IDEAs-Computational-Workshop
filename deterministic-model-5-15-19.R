#Alexandria Purcell
# May-15-2019
#Numerical solution of deterministic epidemiological models 

#Boundary Value problems:
  #What is the trajectory of the quantity in the time of interest
  #phase portrait: has time being implicit and the two variables graphed
    #will be plotted in a relationship 
  #numerical integration: approximation, very good with small error
    #to obtain an approximate solution by solving a sequence of tractable
    #linear approximations at smaller and smaller step sizes until a 
    #specified tolerance is acheived 
  #in R numerical integration is readily performed by using a package 
    #called (deSolve)
  #particularly the function ode is useful since it automatically selects
    #the optimal solving algorithm based on numerical performance 


library(deSolve)
require(deSolve)
sir.model.closed<-function(t,x,params){
  S<-x[1]
  I<-x[2]
  R<-x[3]
  with(
    as.list(params),
    {
      dS<--beta*S*I
      dI<-beta*S*I-gamma*I
      dR<-gamma*I
      dx<-c(dS,dI,dR)
      list(dx)
    }
    )
}

times<-seq(0,120,by=5)
params<-c(beta=0.3,gamma=1/7)
xstart<-c(S=9999/10000,I=1/10000,R=0)

out<-as.data.frame(ode(xstart,times,sir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#Exercise 1: Explore the dynamics of the system for differnt values of
  #the B and y parameters by simulating and plotting trajectories as time 
  #series in phase and space (e.g.I vs S)

#changed the B parameter to 0.8 and kept gamma the same 
params<-c(beta=0.8,gamma=1/7)
xstart<-c(S=9999/10000,I=1/10000,R=0)

out<-as.data.frame(ode(xstart,times,sir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#changed the B parameter to 0.01 and kept gamma the same 
params<-c(beta=0.01,gamma=1/7)
xstart<-c(S=9999/10000,I=1/10000,R=0)

out<-as.data.frame(ode(xstart,times,sir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#change gamma to 1 and keep beta the same: faster recovery phase 
params<-c(beta=0.3,gamma=1)
xstart<-c(S=9999/10000,I=1/10000,R=0)

out<-as.data.frame(ode(xstart,times,sir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#change gamma to 1/100 and keep beta the same: slower recovery phase 
params<-c(beta=0.3,gamma=1/100)
xstart<-c(S=9999/10000,I=1/10000,R=0)

out<-as.data.frame(ode(xstart,times,sir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#Exercise 2: explore the dynamics of the system for one set of B and y at different
  #initial conditions. What happens if there is pre existing immunit in the population 

#change the number of recovered individuals to represent the pre existing immunity individuals 
xstart<-c(S=999/10000,I=1/10000,R=9000)

out<-as.data.frame(ode(xstart,times,sir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#change the number of infected individuals in the population 
xstart<-c(S=9000/10000,I=1000/10000,R=0)

out<-as.data.frame(ode(xstart,times,sir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#change the number of susceptibles in the population 
xstart<-c(S=10000/10000,I=0,R=0)

out<-as.data.frame(ode(xstart,times,sir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)


#Exercise 3: modify the codes given to study the dynamics of a 
  #demographically open SIR model 

sir.model.open<-function(t,x,params){
  S<-x[1]
  I<-x[2]
  R<-x[3]
  with(
    as.list(params),
    {
      dS<-b-beta*S*I-u*S
      dI<-beta*S*I-gamma*I-u*I
      dR<-gamma*I-u*R
      dx<-c(dS,dI,dR)
      list(dx)
    }
  )
}

times<-seq(0,120,by=5)
params<-c(beta=0.3,gamma=1/7,b=1/(75*365), u=1/(75*365))
xstart<-c(S=9999/10000,I=1/10000,R=0)
out<-as.data.frame(ode(xstart,times,sir.model.open,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)
#the open model reaches equilibrium, while in the closed model the 
  #epidemic comes to an end

#Exercise 4: modify the codes given to study the dynamics of an SEIR model

seir.model.closed<-function(t,x,params){
  S<-x[1]
  E<-x[2]
  I<-x[3]
  R<-x[4]
  with(
    as.list(params),
    {
      dS<--beta*S*I
      dE<-beta*S*I-sigma*E
      dI<-beta*S*I-gamma*I
      dR<-gamma*I
      dx<-c(dS,dE,dI,dR)
      list(dx)
    }
  )
}

times<-seq(0,120,by=5)
params<-c(beta=0.3,gamma=1/7,sigma=0.5)
xstart<-c(S=9899/10000,E=100/10000,I=1/10000,R=0)

out<-as.data.frame(ode(xstart,times,seir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#increased the number of exposed individuals 
seir.model.closed<-function(t,x,params){
  S<-x[1]
  E<-x[2]
  I<-x[3]
  R<-x[4]
  with(
    as.list(params),
    {
      dS<--beta*S*I
      dE<-beta*S*I-sigma*E
      dI<-sigma*E-gamma*I
      dR<-gamma*I
      dx<-c(dS,dE,dI,dR)
      list(dx)
    }
  )
}

times<-seq(0,120,by=5)
params<-c(beta=0.3,gamma=1/7,sigma=0.5)
xstart<-c(S=7999/10000,E=2000/10000,I=1/10000,R=0)

out<-as.data.frame(ode(xstart,times,seir.model.closed,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)

#using an open SEIR model 
seir.model.open<-function(t,x,params){
  S<-x[1]
  E<-x[2]
  I<-x[3]
  R<-x[4]
  with(
    as.list(params),
    {
      dS<-b-beta*S*I-u*S
      dE<-beta*S*I-sigma*E
      dI<-sigma*E-gamma*I-u*I
      dR<-gamma*I-u*R
      dx<-c(dS,dE,dI,dR)
      list(dx)
    }
  )
}

times<-seq(0,120,by=5)
params<-c(beta=0.3,gamma=1/7,b=1/(365*75), u=1/(365*75), sigma=0.5)
xstart<-c(S=9899/10000,E=100/10000,I=1/10000,R=0)
out<-as.data.frame(ode(xstart,times,seir.model.open,params))

op<-par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='b')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='b',yaxt='n',xlab='S')
par(op)