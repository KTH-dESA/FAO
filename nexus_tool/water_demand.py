#Standard library imports
import pandas as pd
import numpy as np
import math

#Related third party imports
import pyeto

math.exp = np.exp
math.pow = np.power
math.sqrt = np.sqrt

#### default values:
eto = 'ETo_'
lat = 'lat'
elevation = 'elevation'
wind = 'wind'
srad = 'srad'
tmin = 'tmin'
tmax = 'tmax'
tavg = 'tavg'
eff = 'eff_'
prec = 'prec'
kc = 'kc_'
crop_share = 'crop_share'
crop_area = 'crop_area'
start = '_start'
end = '_end'
crop_column = 'crop'
acwr = 'ACWR_' 
pcwr = 'PCWR_'
pwd = 'PWD_'
sswd = 'SSWD_'
    
def set_cropland_share(df, crop_var, geo_boundary = 'global', 
                       boundary_name = None, crop_share = crop_share):
    if type(crop_var)!= dict:
        df[crop_share] = [{}]*df.shape[0]
        crop_list = set(df[crop_var])
        crop_dic = {i: 0 for i in crop_list}
        crop_dic_list = []
        for crop in crop_list:
            temp_dic = crop_dic.copy()
            temp_dic[crop] = 1
            df.loc[df[crop_var] == crop, crop_share] = [temp_dic]
    else:
        if geo_boundary == 'global':
            df[crop_share] = [crop_var]*df.shape[0]
        else:
            if boundary_name != None:
                if crop_share not in df.columns:
                    df[crop_share] = [{}]*df.shape[0]

                df.loc[df[geo_boundary]==boundary_name, crop_share] = [crop_var]
            else:
                print('Please provide a geo_boundary (e.g. "Province") and a boundary_name (e.g. "Name of province")')
    return df
  
def get_ky_list(df, crop_share = crop_share):
    temp_dic = {i: 0 for i in df.loc[0,crop_share].keys()}
    return temp_dic
    
def get_kc_list(df, crop_share = crop_share):
    temp_dic = {i: [0, 0, 0, 0] for i in df.loc[0,crop_share].keys()}
    return temp_dic

def get_evap_i(lat,elev,wind,srad,tmin,tmax,tavg,month):
    if month ==1:
        J = 15
    else:
        J = 15 + (month-1)*30
        
    latitude = pyeto.deg2rad(lat)
    atmosphericVapourPressure = pyeto.avp_from_tmin(tmin)
    saturationVapourPressure = pyeto.svp_from_t(tavg)
    ird = pyeto.inv_rel_dist_earth_sun(J)
    solarDeclination = pyeto.sol_dec(J)
    sha = [pyeto.sunset_hour_angle(l, solarDeclination) for l in latitude]
    extraterrestrialRad = [pyeto.et_rad(x, solarDeclination,y,ird) for 
                           x, y in zip(latitude,sha)]
    clearSkyRad = pyeto.cs_rad(elev,extraterrestrialRad)
    netInSolRadnet = pyeto.net_in_sol_rad(srad*0.001, albedo=0.23)
    netOutSolRadnet = pyeto.net_out_lw_rad(tmin, tmax, srad*0.001, clearSkyRad, 
                                           atmosphericVapourPressure)
    netRadiation = pyeto.net_rad(netInSolRadnet,netOutSolRadnet)
    tempKelvin = pyeto.celsius2kelvin(tavg)
    windSpeed2m = wind
    slopeSvp = pyeto.delta_svp(tavg)
    atmPressure = pyeto.atm_pressure(elev)
    psyConstant = pyeto.psy_const(atmPressure)
    
    return pyeto.fao56_penman_monteith(netRadiation, tempKelvin, windSpeed2m, 
                                       saturationVapourPressure, 
                                       atmosphericVapourPressure,
                                       slopeSvp, psyConstant, shf=0.0)

def get_eto(df, eto = eto, lat = lat, elevation = elevation, 
        wind = wind, srad = srad, tmin = tmin, tmax = tmax,
        tavg = tavg):
    '''
    calculate ETo for each row for each month 
    '''
    for i in range(1,13):
        df['{}{}'.format(eto, i)]=0
        df['{}{}'.format(eto, i)] = get_evap_i(df[lat],
                                               df[elevation],
                                               df['{}{}'.format(wind, i)],
                                               df['{}{}'.format(srad, i)],
                                               df['{}{}'.format(tmin, i)],
                                               df['{}{}'.format(tmax, i)],
                                               df['{}{}'.format(tavg, i)],
                                               i)
    return df

def get_eff_rainfall_i(prec,eto):
    return (1.253*((prec**0.824)-2.935))*10**(0.001*eto)
    
def get_effective_rainfall(df, eff = eff, prec = prec, eto = eto):
    for i in range(1,13):
        df['{}{}'.format(eff, i)]=0
        df.loc[df['{}{}'.format(prec, i)] < 12.5, '{}{}'.format(eff, i)] = df['{}{}'.format(prec, i)]/30
        df.loc[df['{}{}'.format(prec, i)] >= 12.5, '{}{}'.format(eff, i)] = get_eff_rainfall_i(df['{}{}'.format(prec, i)],df['{}{}'.format(eto, i)])/30 
    return df
    
def get_season_days(crop_calendar, season, start = start, end = end):
    season_start = pd.to_datetime(crop_calendar["".join([season, start])], 
                                  format='%d/%m')
    season_end = pd.to_datetime(crop_calendar["".join([season, end])], 
                                format='%d/%m')
    crop_calendar["_".join([season, 'days'])] = ((season_end - season_start).dt.days+1) % 365
    return crop_calendar

def get_calendar_days(crop_calendar, seasons, start = start, end = end):
    for season in seasons:
        crop_calendar = get_season_days(crop_calendar, season, start = start, end = end)
    return crop_calendar

def get_kc_i(plantation,Li,Ld,Lm,Le,kci,kcd,kcm,kce,isodate):
    """
    Each crop goes through four growing stages: 
    initial - development - mid-season and end-season 
    (check FAO-56 chapter 6 for more details)

    Inputs:
    Plantation = plantation datetime 
    Li = length of the initial stage (in days)
    Ld = length of the development stage (in days)
    Lm = length of the mid-season stage (in days)
    Le = length of the end-season stage (in days)

    kci = crop coefficient 'kc' at the initial stage. In this stage the ckc 
          value is constant and equal to kci
    kcm = crop coefficient 'kc' at the mid-season stage.
          In this stage the ckc value is constant and equal to kcm
    kce = crop coefficient 'kc' at the end-season stage. 
          In this stege the ckc value varies linearly between kce and kcm. 
          check equation 66 (page 132, FAO56). 
    isodate = current date (optional)

    Outputs: 
    * ckc : current crop coefficient, which is constant in the initial and 
            mid-season stages and varies linearly in the development   
            (increasing) and end-season (declining) stages. 

    Some Examples:
    Kc(plantation="2014-01-01",Li=25,Ld=25,Lm=30,Le=20,Kci=0.15,
       Kcm=1.19,Kce=0.35,isodate="2014-01-20") >>> 0.15
     
    Kc(plantation="2014-01-01",Li=25,Ld=25,Lm=30,Le=20,Kci=0.15,
       Kcm=1.19,Kce=0.35,isodate="2014-02-10") >>> 0.774
     
    Kc(plantation="2014-01-01",Li=25,Ld=25,Lm=30,Le=20,Kci=0.15,
       Kcm=1.19,Kce=0.35,isodate="2014-03-12") >>> 1.19
     
    Kc(plantation="2014-01-01",Li=25,Ld=25,Lm=30,Le=20,Kci=0.15,
       Kcm=1.19,Kce=0.35,isodate="2014-04-06") >>> 0.559
    """
    #step 1: 
    
    plantation = pd.to_datetime(plantation, format='%d/%m') #converting the plantation input info to data time
    isodate = pd.to_datetime(isodate , format='%d/%m')  #converting the current date input info to data time
    test = ((isodate-plantation).days+1)%365   #The difference in days between the current day and the plantation day.
    
    # Setting the plantation date and the current date (this is not used)
    Jc = test   
    Jp = 0
    J = (Jc - Jp)%365  # %365 means the remaing days of the year
    
    #Step 2: Calculating the day of the year when each crop stage ends placing the date in the number of days year betweem 0 (1/jan) and 365 (31/Jan)
    JLi = Jp + Li    #end of initial stage = plantation date + lenght of initial stage
#     JLi2 = JLi1 + Li2
    JLd = JLi + Ld   #end of development stage = end of initial stage + length of development stage
    JLm = JLd + Lm   #end of mid-season stage = end of development stage + length of mid-season stage
    JLe = JLm + Le   #end of end-season stage = end of mid-season stage + length of end-season stage

    #step 3: calculating ckc based on the end of each stage date

    if Jc > Jp and Jc < JLe:   #if the current date is greater than the plantation date and it is greater than the end of end-season stage
        if Jc <= JLi:    
            ckc = kci  #if the current date is before the end of initial stage then ckc = kci the coefficient of the initial stege
#         elif Jc > JLi1 and Jc <=JLi2: #New: to account for two init stages
#             ckc = kci2
        elif Jc > JLi and Jc <=JLd:  #if the current date is betweeen the end of the intial stege and the end of the development stage, then ckc is computed based on equation 66 (page 132.FAO56)
            ckc = kci + ((Jc-JLi)/Ld * (kcd-kci))
        elif Jc > JLd and Jc <= JLm: 
            ckc = kcm
        elif Jc > JLm and Jc <= JLe:
            ckc = kcm + ((Jc-JLm)/Le * (kce-kcm))       
    else:
        ckc = 0
    
    return ckc

def get_kc_values(crop_calendar, seasons, kc_dict, crop_column = crop_column, 
                  start = start, end = end, kc = kc):
    for i in range(1,13):
        crop_calendar['{}{}'.format(kc, i)]=0
        
    for index,row in crop_calendar.iterrows():
        crop = row[crop_column]
        for i in range(0,12):
            init_start = pd.to_datetime(crop_calendar[''.join([seasons[0], start])].iloc[index], 
                                         format='%d/%m') #read the plant start date from excel. 
            day_start= (init_start.day+1-31)%31   #what does this represent??   
            
            if (init_start.day-1==30):
                month_start = (init_start.month+1-12)%12  #next month
            else:
                month_start = (init_start.month-12)%12  #the current month
                
            month_start = (month_start+i)%12
            if (month_start==0):
                month_start = 12
            crop_calendar.loc[index,'{}{}'.format(kc, month_start)] = \
                get_kc_i(crop_calendar[''.join([seasons[0], start])].iloc[index],
                crop_calendar['_'.join([seasons[0], 'days'])].iloc[index],
                crop_calendar['_'.join([seasons[1], 'days'])].iloc[index],
                crop_calendar['_'.join([seasons[2], 'days'])].iloc[index],
                crop_calendar['_'.join([seasons[3], 'days'])].iloc[index],
                kc_dict[crop][0],
                kc_dict[crop][1],
                kc_dict[crop][2],
                kc_dict[crop][3],
                '{}/{}'.format(day_start,month_start))
    return crop_calendar

def get_harvest_fraction(i, crop_calendar, crop, init, late, 
                         crop_column, start_name = start, end_name = end):
    if i != 12:
        current_date = pd.to_datetime((i+1),format='%m')
    else:
        current_date = pd.to_datetime(1,format='%m')
    start = pd.to_datetime(crop_calendar.loc[crop_calendar[crop_column]==crop,''.join([init,start_name])], format='%d/%m') #defining the plant start date from excel and setting the correct month and days sequence to read.
    length = crop_calendar.loc[crop_calendar[crop_column]==crop,'_'.join([init,'days'])].iloc[0]
    days = ((current_date - start).iloc[0].days) % 365
    late_end = pd.to_datetime(crop_calendar.loc[crop_calendar[crop_column]==crop,''.join([late,end_name])], format='%d/%m').iloc[0]
    all_days = ((late_end - start).iloc[0].days+1) % 365
    if all_days == 0:
        all_days = 365
   
    if days == 0:
        start = pd.to_datetime(crop_calendar.loc[crop_calendar[crop_column]==crop,''.join([late,start_name])], format='%d/%m') #defining the plant start date from excel and setting the correct month and days sequence to read.
        length = crop_calendar.loc[crop_calendar[crop_column]==crop,'_'.join([late,'days'])].iloc[0]
        days = ((current_date - start).iloc[0].days) % 365
        if days <= length:
            return 1 #- days / length
        else:
            return 0
    elif days <= length:
        return days / length
    elif days <= all_days:
        return 1
    else:
        return 0
        
def get_water_demand(df, crop_calendar, ky_dict, crop_column, aeff, deff, 
                     init_season, late_season, pumping_hours_per_day, 
                     crop_area = crop_area, _eto = eto, _kc = kc, _eff = eff, 
                     _acwr = acwr, _pcwr = pcwr, _pwd = pwd, _sswd = sswd, 
                     start = start, end = end, crop_share = crop_share):
    for i in range (1,13):
        df['{}{}'.format(_pcwr,i)]=0 #PCWR: Peak Crop Water Requirement (l/s/ha) or "Duty", Previously PDWR
        df['{}{}'.format(_pwd,i)]=0  #PWD: Peak Water Demand in (l/s)
        df['{}{}'.format(_sswd,i)]=0 #SSWD: Seasonal Scheme Water Demand in (m3)
        
    #STEP 1: Compute the ACWR from ETc - check FAO1992- page 43-

    #acwr=row['ETo_{}'.format(i)]*30*row['kc_{}'.format(i)] - row['eff_{}'.format(i)]*30 - row['awc']/12 ))
    #once the available water content layer is obtained, the last past should be added to the equation

    for crop in crop_calendar[crop_column]:
        for i in range(1,13):
            eto = f'{_eto}{i}'
            kc = f'{_kc}{i}'
            eff = f'{_eff}{i}'
            acwr = f'{_acwr}{i}_{crop}'
            pcwr = f'{_pcwr}{i}'
            pwd = f'{_pwd}{i}'
            sswd = f'{_sswd}{i}'
            df[kc+'_'+crop] = float(crop_calendar.loc[crop_calendar[crop_column]==crop,kc])
            ky=ky_dict[crop] #Yield response factor coeff = 0.8 for date palms, source TABLE 53-FAO: http://www.fao.org/3/y4360e/y4360e0b.htm 
            df[acwr] = (df[eto]*30*df[kc+'_'+crop]*ky - df[eff]*30 - (0.12*df[eff])*30) #Assumption: awc=12% effective rainfall
            df.loc[df[acwr]<0,acwr] = 0
            df[pcwr] += ((df[acwr]*10)/30)*2*0.012
            
            df[f'harvest_{i}_'+crop] = df[crop_area] * \
                    np.array([x[crop] for x in df[crop_share]]) * \
                    get_harvest_fraction(i,crop_calendar, crop, init_season, 
                                         late_season, crop_column, 
                                         start_name = start, end_name = end)
            df[pwd] += (df[pcwr] *(df[f'harvest_{i}_'+crop]*24))/(pumping_hours_per_day*aeff*deff)
            df[sswd] += (df[acwr]*10*(df[f'harvest_{i}_'+crop])/(aeff*deff))
           