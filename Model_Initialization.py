### Import Packages ###
import Calculate_ITP_Ice_Ocean_Stress as io_stress
import netCDF4 as nc
import seawater as sw
# from netCDF4 import Dataset, MFDataset, num2date
import numpy as np
import pandas as pd
from itp.itp_query import ItpQuery
from datetime import datetime
import matplotlib.dates as mdates
from scipy.interpolate import interp1d

def initialize(start_itp, end_itp, year, start, end, summer=False, io_drag_increased=False, io_drag_vary=False):

    itp = start_itp
    movicemot = False
    stationary=False
    intersect=False

    if summer == True:
        itpname1 = itp
    else:
        itpname1 = '{itp}winter'.format(itp=itp)

    if stationary==True:
        stat = 'stat'
    else:
        stat = ''

    # Set time range for ITP data
    start_date = start
    end_date = end

    time_range = [start_date, end_date]

    print_startdate = start_date.strftime('%Y_%m_%d')
    print_enddate = end_date.strftime('%Y_%m_%d')

    print(print_startdate, print_enddate)

    thermal = pd.read_csv("D:\ERA5-data\export\Surface_net_thermal_radiation_{itp}_{start}_{end}{stat}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    thermal.rename(columns={"0": "index", "0": "thermal"}, inplace=True)

    solar = pd.read_csv("D:\ERA5-data\export\Surface_solar_radiation_downwards_{itp}_{start}_{end}{stat}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    solar.rename(columns={"0": "index", "0": "solar"}, inplace=True)

    sens = pd.read_csv("D:\ERA5-data\export\Surface_sensible_heat_flux_{itp}_{start}_{end}{stat}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    sens.rename(columns={"0": "index", "0": "sens"}, inplace=True)

    latent = pd.read_csv("D:\ERA5-data\export\Surface_latent_heat_flux_{itp}_{start}_{end}{stat}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    latent.rename(columns={"0": "index", "0": "latent"}, inplace=True)

    seaice = pd.read_csv("D:\ERA5-data\export\Sea_ice_cover_{itp}_{start}_{end}{stat}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    seaice.rename(columns={"0": "index", "0": "ice"}, inplace=True)

    time = pd.read_csv("D:\ERA5-data\export\Surface_solar_radiation_downwards_dates{itp}_{start}_{end}{stat}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    time.rename(columns={"0": "index", "0": "time"}, inplace=True)

    if stationary == True:
        ### Ice motion data from NSIDC ###
        taux = pd.read_csv(r"D:\ERA5-data\export\ice_motion_u_{itp}_{start}_{end}{stat}.csv"
                                   .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
        taux.rename(columns={"0": "index", "0": "taux"}, inplace=True)

        tauy = pd.read_csv(r"D:\ERA5-data\export\ice_motion_v_{itp}_{start}_{end}{stat}.csv"
                                    .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
        tauy.rename(columns={"0": "index", "0": "tauy"}, inplace=True)

        geox = pd.read_csv(r"D:\ERA5-data\export\geostrophic_u_{itp}_{start}_{end}{stat}.csv"
                                   .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
        geox.rename(columns={"0": "index", "0": "geox"}, inplace=True)

        geoy = pd.read_csv(r"D:\ERA5-data\export\geostrophic_v_{itp}_{start}_{end}{stat}.csv"
                                   .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
        geoy.rename(columns={"0": "index", "0": "geoy"}, inplace=True)

        drift_dates = pd.read_csv("D:\ERA5-data\export\ice_motion_u_dates{itp}_{start}_{end}{stat}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
        drift_dates.rename(columns={"0": "index", "0": "drift_dates"}, inplace=True)
    #     drift_dates.drop([0], axis=0, inplace=True);
    elif movicemot == True:
         ### Ice motion data from NSIDC ###
        taux = pd.read_csv(r"D:\ERA5-data\export\mov_ice_motion_u_{itp}_{start}_{end}{stat}.csv"
                                   .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
        taux.rename(columns={"0": "index", "0": "taux"}, inplace=True)

        tauy = pd.read_csv(r"D:\ERA5-data\export\mov_ice_motion_v_{itp}_{start}_{end}{stat}.csv"
                                    .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
        tauy.rename(columns={"0": "index", "0": "tauy"}, inplace=True)

        drift_dates = pd.read_csv("D:\ERA5-data\export\mov_ice_motion_dates{itp}_{start}_{end}{stat}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
        drift_dates.rename(columns={"0": "index", "0": "drift_dates"}, inplace=True)
    #     drift_dates.drop([0], axis=0, inplace=True);
    else:
        ### Calculated ice-ocean stress from the ITP drift speeds ###
        io_stress.io_stress(start_itp, end_itp, year, start, end, summer=summer, stationary=stationary, intersect=intersect,
                            io_drag_increased=False, io_drag_vary=io_drag_vary)
        taux = pd.read_csv("D:\ERA5-data\export\wind_stress_x_{itp}_{start}_{end}.csv"
                                     .format(itp=itpname1, start = print_startdate, end = print_enddate))
        taux.rename(columns={"0": "index", "0": "taux"}, inplace=True)

        tauy = pd.read_csv("D:\ERA5-data\export\wind_stress_y_{itp}_{start}_{end}.csv"
                                     .format(itp=itpname1, start = print_startdate, end = print_enddate))
        tauy.rename(columns={"0": "index", "0": "tauy"}, inplace=True)

        drift_dates = pd.read_csv("D:\ERA5-data\export\drift_date_{itp}_{start}_{end}.csv"
                                 .format(itp=itpname1, start = print_startdate, end = print_enddate))
        drift_dates.rename(columns={"0": "index", "0": "drift_dates"}, inplace=True)
        drift_dates.drop([0], axis=0, inplace=True);

    air = pd.read_csv(r"D:\ERA5-data\export\air_temperature_{itp}_{start}_{end}{stat}.csv"
                      .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    air.rename(columns={"0": "index", "0": "air"}, inplace=True)

    uwind = pd.read_csv(r"D:\ERA5-data\export\u_component_of_wind_{itp}_{start}_{end}{stat}.csv"
                        .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    uwind.rename(columns={"0": "index", "0": "uwind"}, inplace=True)

    vwind = pd.read_csv(r"D:\ERA5-data\export\v_component_of_wind_{itp}_{start}_{end}{stat}.csv"
                        .format(itp=itpname1, start = print_startdate, end = print_enddate, stat=stat))
    vwind.rename(columns={"0": "index", "0": "vwind"}, inplace=True)
    
    if io_drag_vary==True:
        iodrag = pd.read_csv('D:\ERA5-data\export\iodrag_{itp}_{start}_{end}.csv'
                             .format(itp=itp, start = print_startdate, end =print_enddate))
        iodrag.rename(columns={"0": "index", "0": "drag"}, inplace=True)
        dragcf = iodrag['drag']
    elif io_drag_increased==True:
        dragcf = 4*5.5e-3
    else:
        dragcf = 5.5e-3

    tlen = len(time['time'])
    rain = np.zeros(tlen)

    thermal = thermal/3600
    solar = solar/3600
    sens = sens/3600
    latent = latent/3600

    hour = 1
    timestep = hour/24
    timearr = np.arange(0,300,timestep)
    t = timearr[0:tlen]


    #Calculate air-ocean stress
    u = uwind['uwind']
    v = vwind['vwind']

    pa = 1.2
    c = .0015
    pac = pa*c

    usq = (u)**2
    vsq = (v)**2

    utao = pac*usq
    vtao = pac*vsq

    def interpolate(inp, fi):
        k, f = int(fi // 1), fi % 1  # Split floating-point index into whole & fractional parts.
        j = k+1 if f > 0 else k  # Avoid index error.
        return (1-f) * inp[k] + f * inp[j]


    d = drift_dates['drift_dates']
    dens = 1027 #Density of seawater
    ### Calculate Ice-Ocean Stress ###
    if stationary == True and io_drag_vary==False:
    #     taux['taux'][:24] = taux['taux'][:24]/100 - geox['geox'][0]
    #     taux['taux'][24:] = taux['taux'][24:]/100 - geox['geox'][1]
    #     tauy['tauy'][:24] = tauy['tauy'][:24]/100 - geoy['geoy'][0]
    #     tauy['tauy'][24:] = tauy['tauy'][24:]/100 - geoy['geoy'][1]
    #     u = taux['taux']
    #     v = tauy['tauy']
        u = (taux['taux']/100)#-geox['geox'][0]
        v = (tauy['tauy']/100)#-geoy['geoy'][0]
        u_sq = u**2
        tx = np.array(dens*dragcf*u_sq)
        v_sq = v**2
        ty = np.array(dens*dragcf*v_sq)
        dnew = np.arange(0,len(tx),1)
        fx = interp1d(dnew, tx, axis=0, fill_value="extrapolate")
        tiox = fx(t)
        fy = interp1d(dnew, ty, axis=0, fill_value="extrapolate")
        tioy = fy(t)
    elif stationary == True and io_drag_vary==True:
        u = (taux['taux']/100)#-geox['geox'][0]
        v = (tauy['tauy']/100)#-geoy['geoy'][0]
        u_sq = u**2
        tx = np.array(dens*dragcf*u_sq)
        v_sq = v**2
        ty = np.array(dens*dragcf*v_sq)
        dnew = np.arange(0,len(tx),1)
        fx = interp1d(dnew, tx, axis=0, fill_value="extrapolate")
        tiox = fx(t)
        fy = interp1d(dnew, ty, axis=0, fill_value="extrapolate")
        tioy = fy(t)
    else:
        tx = taux['taux']
        ty = tauy['tauy']
#         dlen = len(d)
#         print(d)
#         di = datetime.fromisoformat(d[1])
#         df = datetime.fromisoformat(d[dlen-1])
#         ddiff = df-di
#         dstep = int(ddiff.total_seconds())/dlen
#         dnew = np.arange(0,int(ddiff.total_seconds()),dstep)
        new_len = len(t)
        delta = (len(tx)-1) / (new_len-1)
        tiox = [interpolate(tx, k*delta) for k in range(new_len)]
        delta = (len(ty)-1) / (new_len-1)
        tioy = [interpolate(ty, k*delta) for k in range(new_len)]

    f = nc.Dataset('input_data/forcings.nc','w', format='NETCDF4') #'w' stands for write

    f.createDimension('time', None)

    sw = f.createVariable('sw', 'f4', ('time'))
    lw = f.createVariable('lw', 'f4', ('time'))
    qlat = f.createVariable('qlat', 'f4', ('time'))
    qsens = f.createVariable('qsens', 'f4', ('time'))
    tix = f.createVariable('tix', 'f4', ('time'))
    tiy = f.createVariable('tiy', 'f4', ('time'))
    tax = f.createVariable('tax', 'f4', ('time'))
    tay = f.createVariable('tay', 'f4', ('time'))
    precip = f.createVariable('precip', 'f4', ('time'))
    ice = f.createVariable('ice', 'f4', ('time'))
    time = f.createVariable('time', 'f4', ('time'))

    time[:] = t
    sw[:] = solar['solar']
    lw[:] = thermal['thermal']
    qlat[:] = latent['latent']
    qsens[:] = sens['sens']
    tix[:] = tiox
    tiy[:] = tioy
    tax[:] = utao
    tay[:] = vtao
    precip[:] = air['air']
    ice[:] = seaice['ice']
    f.close()
    
    ### Create a NC file for model initial profiles ###
    path = r"C:\Users\eliza\OneDrive\Documents\Yale\itp_final_2022_10_21.db"
    query = ItpQuery(path, system=[itp], pressure=[0, 200], date_time=time_range)
    results = query.fetch()
    sal_data = results[0].absolute_salinity()
    temp_data = results[0].potential_temperature(p_ref=0)
    latitude = results[0].latitude
    depth = results[0].depth()
    lats = [p.latitude for p in results]

    avglat = sum(lats) / len(lats)
    f2 = nc.Dataset('input_data/initial_prof.nc','w', format='NETCDF4') #'w' stands for write
    f2.createDimension('depth', len(depth))
    lat = f2.createVariable('lat', 'f4')  
    z = f2.createVariable('z', 'f4', 'depth')
    t = f2.createVariable('t', 'f4', ('depth'))
    s = f2.createVariable('s', 'f4', ('depth'))
    lat[:] = avglat
    z[:] = depth
    t[:] = temp_data
    s[:] = sal_data
    #Add local attributes to variable instances
    lat.units = 'degrees north'
    t.units = 'Kelvin'
    s.units = 'psu'
    z.units = 'meters'
    f2.close()
    
    print('Done Initializing')
