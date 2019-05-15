library(tidyverse)
library(magrittr)
library(GGally)
ld.data<-read_csv("lymediseasesummarydata.csv")

ld.data%>%ggpairs(column=c("prcp","avtemp","size",'cases'))

ld.data%<>%mutate(log10size=log10(size))
ld.data%<>%mutate(log10cases=log10(cases +1))
#you have to add 1 because the log10 of 1 is 0. 

#re-doing ggplots with log values. 
ld.data%>%ggpairs(column=c("prcp",'avtemp','log10size','log10cases'))

set.seed(222)
ld.data.sample<-ld.data%>%sample_n(100)

myplot<-ggplot(ld.data.sample)+geom_point(aes(prcp,avtemp))+geom_smooth(aes(prcp,avtemp),method = "lm")
myplot

mymodel<-lm(avtemp~prcp, data=ld.data.sample)
mymodel

#Task 7: 
#this is the slope of the line, the slope is 0.006720321
summary(mymodel)$coefficients[2,1]

#this is the p-value/significance, the p value is 3.190341e-06. The slope is significantly different from 0.
summary(mymodel)$coefficients[2,4]

# this is a table for the US population for each year
USpop<-ld.data%>%group_by(year)%>%summarize(sum(size))

#Task 8:
USpop%<>%rename(populatin='sum(size)')
USpop%<>%rename(population='populatin')

USpopplot<-ggplot(USpop)+geom_point(mapping=aes(x=year,y=population))+geom_smooth(aes(year,population), method = "lm")
USpopplot

us.pop.by.state<-ld.data%>%group_by(State)
us.pop.by.state%<>%nest
us.pop.by.state$data[[10]]

library(modelr)
#linear model for predicting population size by year 
lingrowth.model<-function(df2){model<-lm(size~year, data=df2)
return(model)
}

#Task 13:creating by_state column in the us.pop.by.state data set to assign each state a linear model 
us.pop.by.state%<>%mutate(model=purrr::map(data,lingrowth.model))

us.pop.by.state%<>%mutate(resids=purrr::map2(data,model,add_residuals))

#Task 14: the residual values are all negative, it appears that the model predicts to high. 
us.pop.by.state$resids[[10]]

#Task 15:
calc.total.resid<-function(x){
sum(abs(x$resid))
}

us.pop.by.state%<>%mutate(totalResid=purrr::map(resids,calc.total.resid))

#Task 16:
pop.slope<-function(model){
  model$coefficient[2]
}
#Added a column for population slopes for each state
us.pop.by.state%<>%mutate(PopSlope=purrr::map(model,pop.slope))
