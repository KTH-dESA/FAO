library(ggplot2)
library(dplyr)
library(reshape2)

CAPEX <- c(1,3) #write the levels that want to plot inisde c()
Fuel <- c(1,3) #write the levels that want to plot inisde c()
Environmental <- c(1,3) #write the levels that want to plot inisde c()

#Next define the path to the folder containing the results .csv files and the common part of the name in all files
#this path can be of any scenario, as the behaviour of the LCOE will be the same regardless of the scenario, in here I have scenario1
folder <- '/Users/camo/Desktop/KTH Work/NWSAS project/Phase II/Test-20190609/phase2_scenario1_part3 - results/'

#every file is called phase2_scenario1_part3_CAPEX_X_FEnv_XX so the common part will be phase2_scenario1_part3_CAPEX_
common_part <- 'phase2_scenario1_part3_CAPEX_' 

#Now save the file and hit the source botton, pdf file for each results file will be saved in the same folder

for (c in CAPEX){
  for (f in Fuel){
    for (e in Environmental){
      path <- paste(paste0(folder,common_part),c,'_FEnv_',f,e,'.csv',sep='')
      df <- read.csv(path)
      data <- df[c('country','diesel_lcoe','Epricelow','diesel_gasoil_lcoe','pv_lcoe')]
      names(data) = c('Country','Diesel','Grid','Diesel&Gasoil','PV')
      data[data$Country == 'Algeria',c('Grid','Diesel&Gasoil')] = NA
      data[data$Country == 'Tunisia',c('Grid','Diesel')] = NA
      data[data$Country == 'Libya',c('Diesel','Diesel&Gasoil')] = NA
      data = melt(data)
      
      p <- ggplot(data) + geom_boxplot(aes(x=variable, y=value, color=variable)) + facet_grid(.~Country) +
        labs(title='LCOE value per country',subtitle=paste('PV CAPEX level ',c,', fuel cost level ',f,' and environmental cost level ',e, sep=''),
             x='Technology',y='$/kWh',color='Technology') + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())#+ ylim(0.4, 1.6)
      
      #print(p)
      
      ggsave(paste(strsplit(path,'.csv'),'.pdf',sep=''),p,width = 7, height = 5, units='in')
    }
  }
}