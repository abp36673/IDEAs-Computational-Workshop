library(dplyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(GGally)
library(maptools)
library(ggmap)
library(maps)

#Task-1
prism<-read_csv('climate.csv')
ld<-read_csv('lyme.csv')
pop<-read_csv('pop.csv')

#Task-2: It does not fit in the tidy data format because the population size variable does not have its own column,
#it is divided into multiple columns for each year. 

#Task-3

pop %<>% select(fips,starts_with("pop2"))
#comment: selecting the fips column in the pop data
pop %<>% gather(starts_with("pop2"),key="str_year",value="size") %>% na.omit
#gather fips, start year, and size columns, and omits any empty cell. 
pop %<>% mutate(year=str_replace_all(str_year,"pop",""))
#produced a year column and replaced fills it with start year data, and removes the pop so that it is only the year and not the pop year
pop %<>% mutate(year=as.integer(year))
#running year as an integer 
pop %<>% mutate(fips=str_replace_all(fips,"^0",""))
# found the fips starting with zero and replaced with nothing (removed the zeros)
pop %<>% mutate(fips=as.integer(fips))
#running fips as an integer 

#comment: you could remove the state wide summaries by deleting the fip codes that are associated with each state

#comment: you could remove the start year data by selecting the subset and -start_with('str_year') to choose everything but the start year and send that back to the data frame with %<>%

#gathering the data and creating a new column called cases and the cells are filled with the value of cases
ld%<>%gather(starts_with('Cases'),key="str_year", value="cases")
#produces a year column and fills the cells with the values from the str_year values that were associated with cases. 
ld%<>%mutate(year=str_replace_all(str_year, "Cases",""))
#converting the year values to be integers
ld%<>%mutate(year=as.integer(year))
#renaming state and county codes
ld%<>%rename(State=STNAME, County=CTYNAME)

function(state.code,county.code){
  if(str_length(county.code)==3){
    fips<-paste(as.character(state.code), as.character(county.code), sep="")%>% as.integer
  }
  else if(str_length(county.code)==2){
    fips<-paste(as.character(state.code),"0", as.character(county.code), sep="") %>% as.integer
  }
  else{fips<-paste(as.character(state.code),"00", as.character(county.code), sep='') %>% as.integer}
  return(fips)
}
ld%<>%rowwise()%<>%mutate(fips=fips.builder(STCODE,CTYCODE))
ld%<>%select(-c(STCODE,CTYCODE,str_year))

# to merge we are going to want to join based on the FIP codes 
ldprism<-inner_join(ld,prism)
#merging the ldprism data set with the population data set
summary.data<-inner_join(ldprism,pop)

#yearly cases arranged in descending order 
yearly.cases<-ld%>%ungroup%>%group_by(year)%>%summarize(total=sum(cases))%>%arrange(desc(total))

#average cases in each state 
state.by.state<-ld%>%ungroup%>%group_by(State)%>%summarize(total=mean(cases))%>%arrange(desc(total))

#saving the file
save(summary.data,file='get_summary.data.Rda')
save(yearly.cases,file='get_yearly.cases.Rda')
save(state.by.state,file='get_state.by.state.Rda')
save(climate,file='get_climate.Rda')
save(ld,file='get_ld.Rda')
save(pop,file='get_pop.Rda')
save(ldprism,file='get_ldprism.Rda')

write_csv(summary.data,path="lymediseasesummarydata.csv")


