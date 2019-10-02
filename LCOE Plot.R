library(ggplot2)
library(wesanderson)
library(reshape2)
library(gganimate)

CAPEX <- c(1,2,3) #write the levels that want to plot inside c()
Fuel <- c(1,2,3) #write the levels that want to plot inside c()
Environmental <- c(1) #write the levels that want to plot inisde c()

#Next define the path to the folder containing the results .csv files and the common part of the name in all files
#this path can be of any scenario, as the behaviour of the LCOE will be the same regardless of the scenario, in here I have scenario1
folder <- 'C:/Users/camilorg/Box Sync/NWSAS/NWSAS_Phase II/NewRun20190801/surf45_TDS2k_energyDes - results/'

#every file is called phase2_scenario1_part3_CAPEX_X_FEnv_XX so the common part will be phase2_scenario1_part3_CAPEX_
common_part <- 'surf45_TDS2k_energyDes_CAPEX_' 

#Now save the file and hit the source botton, pdf file for each results file will be saved in the same folder

all_data = data.frame('Country'=NULL,'Diesel'=NULL,'Elec'=NULL,
                      'PV'=NULL,'CAPEX'=NULL,
                      'Fuel'=NULL)

for (c in CAPEX){
  for (f in Fuel){
    for (e in Environmental){
      path <- paste(paste0(folder,common_part),c,'_FEnv_',f,e,'.csv',sep='')
      df <- read.csv(path)
      data <- df[c('country','diesel_lcoe','elec_lcoe','pv_lcoe')]
      names(data) = c('Country','Diesel','Elec','PV')
      data[data$Country == 'Algeria',c('Elec')] = NA
      data[data$Country == 'Tunisia',c('Elec')] = NA
      data[data$Country == 'Libya',c('Diesel')] = NA
      data = melt(data)
      data = na.omit(data)
      data['CAPEX'] = c
      data['Fuel'] = f
      
      all_data = rbind(all_data,data)
      
      # p <- ggplot(data) + geom_boxplot(aes(x=variable, y=value, color=variable)) + facet_grid(.~Country,scales='free_x') +
      #   labs(title='LCOE value per country',subtitle=paste('PV CAPEX level ',c,', fuel cost level ',f, sep=''),
      #        x='Technology',y='$/kWh',color='Technology') + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())#+ ylim(0.4, 1.6)
      # 
      # # 
      # ggsave(paste(strsplit(path,'.csv'),'.pdf',sep=''),p,width = 5, height = 5, units='in')
      
    }
  }
}

p <- ggplot(all_data) + geom_boxplot(aes(x=Country, y=value, fill=variable, color=variable),alpha=0.6,size=0.5) + facet_grid(paste('CAPEX',CAPEX)~paste('Fuel',Fuel),scales='free_x') +
     labs(title='LCOE value per country',subtitle=paste('PV CAPEX level 1 to ',c,', fuel cost level 1 to ',f, sep=''),x='Country',y='$/kWh',fill='Technology',color='Technology') + 
     theme_bw() + scale_color_manual(values=wes_palette("Darjeeling1", n = 3)) + scale_fill_manual(values=wes_palette("Darjeeling1", n = 3))
    # names(wes_palettes)
    # print(p)
# 

df <- all_data

CAPEX <-  as.matrix(df['CAPEX'])
Fuel <- as.matrix(df['Fuel'])
df['titles'] <- as.data.frame(matrix(paste0('PV CAPEX level ',CAPEX, ', Fuel level ', Fuel)))

Titles <- unique(df['titles'])


p_animate <- p <- ggplot(df) + geom_boxplot(aes(x=Country, y=value, fill=variable, color=variable),alpha=0.6,size=0.5) + #facet_grid(paste('CAPEX',CAPEX)~paste('Fuel',Fuel),scales='free_x') +
  theme_bw() + scale_color_manual(values=wes_palette("Darjeeling1", n = 3)) + scale_fill_manual(values=wes_palette("Darjeeling1", n = 3)) +
  transition_states(
    titles,
    transition_length = 2,
    state_length = 1) +
  labs(title='LCOE value per country',subtitle='{closest_state}',x='Country',y='$/kWh',fill='Technology',color='Technology')

anim_save('LCOE.gif',p_animate,width = 5, height = 5, units='in', res = 300)
ggsave(paste('LOCE_Grid_Des','.png',sep=''),p,width = 8, height = 8, units='in', dpi=300)
ggsave(paste('LOCE_Grid_Des','.pdf',sep=''),p,width = 8, height = 8, units='in')


anim <- ggplot(iris, aes(Sepal.Width, Petal.Width)) +
  geom_point() +
  labs(title = "{closest_state}") +
  transition_states(Species, transition_length = 3, state_length = 1)

head(iris)
