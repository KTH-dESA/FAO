def evap_i(lat,elev,wind,srad,tmin,tmax,tavg,month):
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
    extraterrestrialRad = [pyeto.et_rad(x, solarDeclination,y,ird) for x, y in zip(latitude,sha)]
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
                                       saturationVapourPressure, atmosphericVapourPressure,
                                       slopeSvp, psyConstant, shf=0.0)
    