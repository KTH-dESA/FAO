library(ggplot2)

# there is no need to change the level of the sensitivity variables, as the water demand and energy demand wont change with those
CAPEX <- 1
Fuel <- 1
Environmental <- 1

# define the path to the folder containing all scenario results
folder <- 'C:/Users/camilorg/Box Sync/NWSAS/NWSAS_Phase II/NewRun20190801/'

# define the path of each scenario folder and the common part of the name in all files
#every file is called phase2_scenario1_part3_CAPEX_X_FEnv_XX so the common part will be phase2_scenarioX_part3X_CAPEX_
path1 <- 'surf45_TDS2k_energy - results/surf45_TDS2k_energy_CAPEX_'
path2 <- 'surf65_TDS2k_energy - results/surf65_TDS2k_energy_CAPEX_'
path3 <- 'drip85_TDS2k_energy - results/drip85_TDS2k_energy_CAPEX_'
path4 <- 'surf45_TDS2.5k_energyDes - results/surf45_TDS2.5k_energyDes_CAPEX_'
path5 <- 'surf65_TDS2.5k_energyDes - results/surf65_TDS2.5k_energyDes_CAPEX_'
path6 <- 'drip85_TDS2.5k_energyDes - results/drip85_TDS2.5k_energyDes_CAPEX_'
path7 <- 'surf45_TDS2k_energyDes - results/surf45_TDS2k_energyDes_CAPEX_'
path8 <- 'surf65_TDS2k_energyDes - results/surf65_TDS2k_energyDes_CAPEX_'
path9 <- 'drip85_TDS2k_energyDes - results/drip85_TDS2k_energyDes_CAPEX_'
path10 <- 'surf45_TDS3k_energyDes - results/surf45_TDS3k_energyDes_CAPEX_'
path11 <- 'surf65_TDS3k_energyDes - results/surf65_TDS3k_energyDes_CAPEX_'
path12 <- 'drip85_TDS3k_energyDes - results/drip85_TDS3k_energyDes_CAPEX_'

#save the file and hit the source botton
# 
# path1 <- paste(folder,path1,'.csv',sep='')
# path2 <- paste(folder,path2,'.csv',sep='')
# path3 <- paste(folder,path3,'.csv',sep='')
# path4 <- paste(folder,path4,'.csv',sep='')
# path5 <- paste(folder,path5,'.csv',sep='')
# path6 <- paste(folder,path6,'.csv',sep='')
path1 <- paste(folder,path1,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path2 <- paste(folder,path2,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path3 <- paste(folder,path3,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path4 <- paste(folder,path4,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path5 <- paste(folder,path5,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path6 <- paste(folder,path6,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path7 <- paste(folder,path7,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path8 <- paste(folder,path8,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path9 <- paste(folder,path9,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path10 <- paste(folder,path10,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path11 <- paste(folder,path11,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')
path12 <- paste(folder,path12,CAPEX,'_FEnv_',Fuel,Environmental,'.csv',sep='')

sep = ','
df1 <- read.csv(path1, sep=sep)
df2 <- read.csv(path2, sep=sep)
df3 <- read.csv(path3, sep=sep)
df4 <- read.csv(path4, sep=sep)
df5 <- read.csv(path5, sep=sep)
df6 <- read.csv(path6, sep=sep)
df7 <- read.csv(path7, sep=sep)
df8 <- read.csv(path8, sep=sep)
df9 <- read.csv(path9, sep=sep)
df10 <- read.csv(path10, sep=sep)
df11 <- read.csv(path11, sep=sep)
df12 <- read.csv(path12, sep=sep)

water_data = data.frame('water'=NULL,'NAME_1'=NULL,'total_demand'=NULL,
                       'max_capacity_diesel'=NULL,'max_capacity_pv'=NULL,
                       'least_cost_technology'=NULL,'scenario'=NULL,'country'=NULL)
scenarios <- c('Scenario 1: Surface irrigation','Scenario 2: Improved surface irrigation','Scenario 3: Drip irrigation')
i <- 1
for (df in list(df1,df2,df3)){
    data <- data.frame('water'=rowSums(df[names(df)[grep('SSWD_',names(df))]]))
    # data <- cbind(data,df[c('NAME_1','total_demand','max_capacity_diesel','max_capacity_pv','least_cost_technology','country')])'area_ha'
    data <- cbind(data,df[c('NAME_1','total_demand','country','area_ha')])

    # data[data['least_cost_technology'] == 'diesel_lcoe','capacity'] <- data[data['least_cost_technology'] == 'diesel_lcoe','max_capacity_diesel']
    # data[data['least_cost_technology'] == 'pv_lcoe','capacity'] <- data[data['least_cost_technology'] == 'pv_lcoe','max_capacity_pv']
    data['scenario'] <- scenarios[i]
    water_data = rbind(water_data,data)
    i = i + 1
}

get_energy <- function(df_list,energy_type,scenarios,desal=FALSE){
  i <- 1
  for (df in df_list){
    data <- data.frame('energy_desal'=rowSums(df[names(df)[grep('Edesal_GWh_',names(df))]]))
    data <- cbind(data,df[c('NAME_1','total_demand','max_capacity_diesel','max_capacity_pv','least_cost_technology','country')])
    
    if (desal){
      data['total_demand'] <- data['energy_desal'] * 1000000
    }
    
    data[data['least_cost_technology'] == 'diesel_lcoe','capacity'] <- data[data['least_cost_technology'] == 'diesel_lcoe','max_capacity_diesel']
    data[data['least_cost_technology'] == 'pv_lcoe','capacity'] <- data[data['least_cost_technology'] == 'pv_lcoe','max_capacity_pv']
    data['scenario'] <- scenarios[i]
    data['energy_type'] <- energy_type
    energy_data = rbind(energy_data,data)
    i <- i + 1
    # if (i==4){
    #   i <- 1
    # }
  }
  return(energy_data)
}

# i <- 1
# j <- 1

energy_data = data.frame('water'=NULL,'NAME_1'=NULL,'total_demand'=NULL,
                         'max_capacity_diesel'=NULL,'max_capacity_pv'=NULL,
                         'least_cost_technology'=NULL,'scenario'=NULL,'desalination'=NULL,'country'=NULL)

energy_type <- 'Pumping energy'
df_list <- list(df1,df2,df3)
energy_data <- get_energy(df_list,energy_type,scenarios)

energy_type <- 'Desalination energy'
df_list <- list(df4,df5,df6)
energy_data <- get_energy(df_list,energy_type,scenarios,TRUE)

water = aggregate(cbind(water,area_ha)~NAME_1+country+scenario, water_data, sum)
# area = aggregate(area_ha~NAME_1+country+scenario, water_data, sum)
energy = aggregate(cbind(total_demand,energy_desal)~NAME_1+country+scenario+energy_type, energy_data, sum)
# cost = aggregate(total_demand~NAME_1+scenario+desalination, energy_data, sum)

energy_data = data.frame('water'=NULL,'NAME_1'=NULL,'total_demand'=NULL,
                         'max_capacity_diesel'=NULL,'max_capacity_pv'=NULL,
                         'least_cost_technology'=NULL,'scenario'=NULL,'desalination'=NULL,'country'=NULL)

energy_type <- '2k'
df_list <- list(df7,df8,df9)
energy_data <- get_energy(df_list,energy_type,scenarios,TRUE)

energy_type <- '3k'
df_list <- list(df10,df11,df12)
energy_data <- get_energy(df_list,energy_type,scenarios,TRUE)

energy_sens = aggregate(cbind(total_demand,energy_desal)~NAME_1+country+scenario+energy_type, energy_data, sum)

# energy['sensitivity_low'] <- 0
# energy['sensitivity_high'] <- 0
energy[energy['energy_type']=='Desalination energy','sensitivity_low'] <- energy_sens[energy_sens['energy_type']=='2k','total_demand']
energy[energy['energy_type']=='Desalination energy','sensitivity_high'] <- energy_sens[energy_sens['energy_type']=='3k','total_demand']

# p_water <- ggplot(main_data) + geom_bar(aes(x=prov,y=water/1000000), stat = 'identity', fill='#2980B9') +
#          labs(x='Province',y='Water demand (Mm3/yr)') + 
#          scale_y_continuous(expand = c(0, 0)) + facet_grid(.~scenario) +
#          theme_bw()
p_water <- ggplot(water) + geom_bar(aes(x=reorder(NAME_1,-water),y=water/1000000,fill=scenario), color='white', size=0.2, stat = 'identity', position = 'dodge') +
  # labs(title='Irrigation water demand per province',x=element_blank(),y='Water demand (Mm3/yr)') + 
  labs(x=element_blank(),y='Water demand (Mm3/yr)') +
  scale_y_continuous(expand = c(0, 0))  + scale_fill_manual(values=c("#1b667e", "#548ea3", "#45b4d4")) + #scale_fill_brewer(palette="Set2",direction=-1) + 
  facet_grid(.~country,scales='free_x',space='free_x') +
  theme_minimal() + theme(legend.position="bottom",legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) #+ coord_flip()

p_water_int <- ggplot(water) + geom_bar(aes(x=reorder(NAME_1,-water/area_ha),y=water/area_ha,fill=scenario), color='white', size=0.2, stat = 'identity', position = 'dodge') +
  # labs(title='Irrigation water intensity per province',x=element_blank(),y='Water intensity (m3/ha)') + 
  labs(x=element_blank(),y='Water intensity (m3/ha)') + 
  scale_y_continuous(expand = c(0, 0))  + scale_fill_manual(values=c( "#007065", "#00a79d","#f5c181")) + #scale_fill_brewer(palette="Blues",direction=-1) + 
  facet_grid(.~country,scales='free_x',space='free_x') +
  theme_minimal() + theme(legend.position="bottom",legend.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1)) #+ coord_flip()


p_energy <- ggplot(energy, aes(x=reorder(NAME_1,-total_demand),fill=scenario)) + geom_bar(aes(y=total_demand/1000000), color='white', size=0.2, stat = 'identity', position = 'dodge') +
  # labs(title='Energy demand per province',x=element_blank(),y='Energy demand (GWh/yr)') +
  geom_errorbar(aes(ymin=sensitivity_high/1000000, ymax=sensitivity_low/1000000), width=0.5,
                position=position_dodge(0.9), color="#f6d55c", size = 0.8) +
  labs(x=element_blank(),y='Energy demand (GWh/yr)') + 
  scale_y_continuous(expand = c(0, 0))   + scale_fill_manual(values=c("#4d342f", "#dc5353", "#ff935c")) + #scale_fill_brewer(palette="Oranges",direction=-1) + 
  facet_grid(energy_type~country,scales='free_x',space='free_x') +
  theme_minimal() + theme(legend.position="bottom",legend.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1)) #+ coord_flip()


# ggsave('/Users/camo/Desktop/KTH Work/NWSAS project/water_demand.tiff',p_water,width = 7, height = 5, units='in',dpi=300)
ggsave(paste0(folder,'water_demand.png'),p_water,width = 8, height = 4, units='in', dpi = 300)
ggsave(paste0(folder,'water_intensity.png'),p_water_int,width = 8, height = 4, units='in', dpi = 300)
ggsave(paste0(folder,'energy_demand.png'),p_energy,width = 8, height = 6, units='in', dpi = 300)

ggsave(paste0(folder,'water_demand.pdf'),p_water,width = 8, height = 5, units='in')
ggsave(paste0(folder,'water_intensity.pdf'),p_water_int,width = 8, height = 5, units='in')
ggsave(paste0(folder,'energy_demand.pdf'),p_energy,width = 8, height = 7, units='in')


####### sttatistics calculations ###########
algeria_water = sum(subset(water, scenario == 'Scenario 1: Surface irrigation' & country == 'Algeria')['water'])/1000000
libya_water = sum(subset(water, scenario == 'Scenario 1: Surface irrigation' & country == 'Libya')['water'])/1000000
tunisia_water = sum(subset(water, scenario == 'Scenario 1: Surface irrigation' & country == 'Tunisia')['water'])/1000000
nwsas_water = algeria_water + libya_water + tunisia_water
algeria_prov_water = subset(water, scenario == 'Scenario 1: Surface irrigation' & (NAME_1 == 'Adrar' | NAME_1 == 'El Oued' | NAME_1 == 'Ghardaia' | NAME_1 == 'Ouargla'))

algeria_energy = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & country == 'Algeria' & energy_type == 'Pumping energy')['total_demand'])/1000000
libya_energy = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & country == 'Libya' & energy_type == 'Pumping energy')['total_demand'])/1000000
tunisia_energy = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & country == 'Tunisia' & energy_type == 'Pumping energy')['total_demand'])/1000000
nwsas_energy = algeria_energy + libya_energy + tunisia_energy
algeria_prov_energy = subset(energy, scenario == 'Scenario 1: Surface irrigation' & energy_type == 'Pumping energy' & (NAME_1 == 'Adrar' | NAME_1 == 'El Oued' | NAME_1 == 'Ghardaia' | NAME_1 == 'Ouargla'))

algeria_energy = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & country == 'Algeria' & energy_type == 'Desalination energy')['total_demand'])/1000000
eloeud_energy_desal = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & NAME_1 == 'El Oued' & energy_type == 'Desalination energy')['total_demand'])/1000000
libya_energy_desal = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & country == 'Libya' & energy_type == 'Desalination energy')['total_demand'])/1000000
tunisia_energy_desal = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & country == 'Tunisia' & energy_type == 'Desalination energy')['total_demand'])/1000000
nwsas_energy_desal = algeria_energy + libya_energy + tunisia_energy
algeria_prov_energy_desal = subset(energy, scenario == 'Scenario 1: Surface irrigation' & energy_type == 'Pumping energy' & (NAME_1 == 'Adrar' | NAME_1 == 'El Oued' | NAME_1 == 'Ghardaia' | NAME_1 == 'Ouargla'))

eloeud_energy_2k = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & NAME_1 == 'El Oued' & energy_type == 'Desalination energy')['sensitivity_low'])/1000000
eloeud_energy_3k = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & NAME_1 == 'El Oued' & energy_type == 'Desalination energy')['sensitivity_high'])/1000000
ouargla_energy_2k = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & NAME_1 == 'Ouargla' & energy_type == 'Desalination energy')['sensitivity_low'])/1000000
ouargla_energy_3k = sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & NAME_1 == 'Ouargla' & energy_type == 'Desalination energy')['sensitivity_high'])/1000000


energy_desal_2k <- sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & energy_type == 'Desalination energy')['sensitivity_low'])/1000000
energy_desal_3k <- sum(subset(energy, scenario == 'Scenario 1: Surface irrigation' & energy_type == 'Desalination energy')['sensitivity_high'])/1000000
desal_percentage_reduc <- (energy_desal_2k-energy_desal_3k)/energy_desal_2k




