mers<-read.csv('cases.csv')
head(mers)
class(mers$onset)
mers$hospitalized[890]<-c('2015-02-20')
mers<-mers[-471,]
library(lubridate)
mers$onset2<-ymd(mers$onset)
mers$hospitalized2<-ymd(mers$hospitalized)
class(mers$onset2)
day0<-min(na.omit(mers$onset2))
#Question 1: You use the function na.omit to remove the data points that contain NA, meaning they are missing values. It will find the smallest value excluding the values missing data, since those would be considered the smallest value if they were included. 
mers$epi.day<-as.numeric(mers$onset2-day0)
#Question 2: as. numeric is going to assign the day0 a numeric value instead of a year-month-day Date. 

library(ggplot2)

ggplot(data=mers)+geom_bar(mapping=aes(x=epi.day))+labs(x='Epidemic day',y='Case count',title= 'Global count of MERS cases by date of symptom onset', caption= "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")

#Question 3: If you don't use + symbol in the previous code then the type of graph, such as bar graph, and the labels would not be added. 

ggplot(data=mers)+geom_bar(mapping=aes(x=epi.day,fill=country))+labs(x='Epidemic day',y='Case count',title= 'Global count of MERS cases by date of symptom onset', caption= "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")

ggplot(data=mers)+geom_bar(mapping=aes(x=epi.day,fill=gender))+labs(x='Epidemic day',y='Case count',title= 'Global count of MERS cases by date of symptom onset', caption= "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")
#Question 4: If you modify the argument fill position to other factors such as gender or outcome it will modify the graph to include the genders associated with each case count. Or organize the case counts for each epidemic day for the different outcomes assocated with the cases. 


ggplot(data=mers)+geom_bar(mapping=aes(x=epi.day,fill=country))+labs(x='Epidemic day',y='Case count',title= 'Global count of MERS cases by date of symptom onset', caption= "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")+coord_flip()+coord_polar()
#question 5: flipping the coordinates just flips the X and Y axis. The posittion in the circle is the epidmeic day and how far out the data point goes is the number of cases counted. 
mers$infectous.period<-mers$hospitalized2-mers$onset2
class(mers$infectious.period)
mers$infectious.period<-as.numeric(mers$infectous.period, unit="days")
ggplot(data=mers)+geom_histogram(aes(x=infectious.period))+labs(x='Infectious period', y='Frequency', title='Distribution of calculated MERS infectious period', caption="Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")
mers$infectious.period2<-ifelse(mers$infectious.period<0,0,mers$infectious.period)
ggplot(data=mers)+geom_histogram(aes(x=infectious.period2))+labs(x='Infectious period', y='Frequency', title= 'Distribution of calculated MERS infectious period (positive values only)', caption= "data")
ggplot(data=mers)+ geom_density(mapping=aes(x=infectious.period2))+ labs(x='Infectous period', y='Frequency', title ='probability density for MERS infectious period (positive values only)', caption='data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv')
ggplot(data=mers)+geom_area(stat='bin', mapping=aes(x=infectious.period2))+labs(x='Infectious period', y='Frequency', title='Area plot for MERS infectious period (positive values only)', caption= "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")
ggplot(data=mers)+geom_dotplot(stat='bin', mapping=aes(x=infectious.period2))+labs(x='Infectious period', y='Frequency', title='Area plot for MERS infectious period (positive values only)', caption= "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")
ggplot(data=mers)+geom_area(stat='bin', mapping=aes(x=infectious.period2))+labs(x='Infectious period', y='Frequency', title='Area plot for MERS infectious period (positive values only)', caption= "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")
ggplot(data=mers)+geom_line(stat='bin', mapping=aes(x=infectious.period2))+labs(x='Infectious period', y='Frequency', title='Area plot for MERS infectious period (positive values only)', caption= "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")
ggplot(data=mers)+geom_jitter(aes(x=epi.day, color=country, y=infectious.period2))+labs(x='Epidemic Day', y='Infectious Period')
ggplot(data=mers, mapping=aes(x=epi.day, y=infectious.period2))+ geom_point(mapping=aes(color=country))+ facet_wrap(~country)+scale_y_continuous(limits=c(0,50))+labs(x='Epidemic day', y='Infectious period', title='MERS infectious period (positive values only) over time')
ggplot(data=subset(mers, gender %in% c('M','F')& country %in% c('KSA', 'Oman', 'Iran', 'Jordan', 'Qatar', 'South Korea', 'UAE')), mapping=aes(x=epi.day, y=infectious.period2))+geom_point(mapping=aes(color=country))+facet_grid(gender~country)+scale_y_continuous(limits=c(0,50))+labs(x='Epidemic day', y='infectious period')
library(plotly)
epi.curve<-ggplot(data=mers)+geom_bar(mapping = aes(x=epi.day))+labs(x='Epidemic Day', y='Case count', title='Global count of MERS cases by date of symptom onset', caption="Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv")
ggplotly(epi.curve)



                                                                    