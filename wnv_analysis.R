wnv<-read.csv('wnv2.csv')
head(wnv)
class(wnv$State)
class(wnv$Year)
class(wnv$EncephMen)
library(ggplot2)
class(wnv$Total)

mean<-function(x){
  s<-sum(x)
  n<-length(x)
  mean<-s/n
}

SE<-function(x){
  sd(x)/sqrt(length(x))
}

ggplot(data=wnv)+geom_histogram(mapping=aes(x=wnv$Year, fill=wnv$State)) 
ggplot(data=wnv, mapping=aes(x=wnv$Year, y=wnv$Total))+ geom_point(mapping=aes(color=wnv$State))+facet_wrap(~wnv$State)+
  scale_y_continuous(limits=c(0,3000))+ labs(x='Year', y='Count')
#taking the log of the total counts
log.cases<-log10(wnv$Total)
class(log.cases)
ggplot(data=wnv, mapping=aes(x=wnv$Year, y=log.cases))+ geom_point(mapping=aes(color=wnv$State))+facet_wrap(~wnv$State)+
  scale_y_continuous(limits=c(0,5))+ labs(x='Year', y='Count')
#assigned case fatality rate as deaths divided by total cases 
cfr<-(wnv$Fatal/wnv$Total)
ggplot(data=wnv, mapping=aes(x=wnv$Year, y=cfr))+ geom_point(mapping=aes(color=wnv$State))+facet_wrap(~wnv$State)+
  scale_y_continuous(limits=c(0,1))+ labs(x='Year', y='CFR')

#sum of febrile cases, neuroinvasive cases, and other cases to determine if it is equal to the Total. 
wnv$Total==(wnv$Fever+wnv$EncephMen+wnv$Other)

#calculating the annual case count for each state rounded to the nearest dozen. 
rounded.total<-wnv$Total%/%12
full.rounded<-rounded.total*12
rounded.total*12
rounded.error<-wnv$Total-full.rounded
wnv$Total-full.rounded

ndr<-function(state='Colorado', years=1999:2007){
x<-wnv[wnv$State %in% state & wnv$Year %in% years,]
y<-data.frame(state=x$State, ndr=x$EncephMen/x$Total)
m<-aggregate(y$ndr, by=list(y$state),FUN=mean)
se<-aggregate(y$ndr, by=list(y$state), FUN=function(x) sd(x)/sqrt(length(x)))
out<-merge(m,se,by='Group.1')
names(out)<-c('state','mean.ndr','se.ndr')
return(out)
}

disease<-ndr(state=c('California', 'Colorado', 'New York'))

ggplot(disease, aes(x=state, y=mean.ndr, fill=state))+geom_bar(stat='identity')+ geom_errorbar(aes(ymin=mean.ndr-se.ndr, ymax=mean.ndr-se.ndr))+labs(x='State', y='Neuroinvasive disease rate')

ndr<-function(state='colorado', years=1999:2007){

}


wnv %>% filter(State %in% c("x") & Year %in% 1999:2007) %>% group_by(State) %>% summarize(ndr=)




  