import sys, datetime
import numpy as np
import matplotlib.pyplot as pp
sys.path.append('/home/liamk/seshat/glowcond/glow098release/GLOW/')
import pyglow098
from geospacepy import special_datetime,dmspcdf_tools,dmsp_spectrogram, satplottools
#from ovationpyme import ovation_prime

#op_cond_estimator = ovation_prime.ConductanceEstimator(datetime.datetime(2010,5,27,12),datetime.datetime(2010,5,30,12))

year,month,day = 2010,5,29
sat = 16

cdf = dmspcdf_tools.get_cdf(sat,year,month,day,'ssj')
cdfm = dmspcdf_tools.get_cdf(sat,year,month,day,'ssm')
cdfies = dmspcdf_tools.get_cdf(sat,year,month,day,'ssies')

dts = cdf['Epoch'][:]
uts = special_datetime.datetimearr2sod(cdf['Epoch'][:])
auroral_region = cdf['AURORAL_REGION'][:].flatten()
chen = cdf['CHANNEL_ENERGIES'][:]
glats = cdf['SC_GEOCENTRIC_LAT'][:]
glons = cdf['SC_GEOCENTRIC_LON'][:]
mlats = cdf['SC_APEX_LAT'][:]
mlts = cdf['SC_APEX_MLT'][:]
oi = cdf['ORBIT_INDEX'][:]

auroral = np.logical_and(auroral_region > 1,glats>0.)

auroral_eflux = cdf['ELE_DIFF_ENERGY_FLUX'][:][auroral,:]*1.6e-12
auroral_flux = cdf['ELE_DIFF_ENERGY_FLUX'][:][auroral,:]/chen
auroral_avg_energy = cdf['ELE_AVG_ENERGY'][:][auroral]
auroral_total_eflux = cdf['ELE_TOTAL_ENERGY_FLUX'][:][auroral]*1.6e-12
deltaB = cdfm['DELTA_B_APX'][:][auroral,:]
Vy = cdfies['VEL_ION_HORIZ_IDM'][:][auroral]

#subset = np.arange(1000,1900)
subset = np.arange(6700,7900)
orbit = 3
subset = np.flatnonzero(oi[auroral]==orbit)
#subset = np.arange(4000,4800)
#subset = np.arange(1800,1920)
#subset = np.arange(1801,1821)
#subset = np.arange(np.count_nonzero(auroral))


#Format input argument vectors
dts = dts[auroral][subset]
uts = uts[auroral][subset]
glats = glats[auroral][subset]
glons = glons[auroral][subset]
mlats = mlats[auroral][subset]
mlts = mlts[auroral][subset]
dBe,dBn = deltaB[subset,0],deltaB[subset,1]
idate,f107a,f107,f107p,ap = 16355,70,70,70,4
ncond = int(len(subset))
#glats = np.ones((ncond,))*70.
#glons = np.ones((ncond,))*0.
idates = np.ones((ncond,))*idate
f107as = np.ones((ncond,))*f107a
f107s = np.ones((ncond,))*f107
f107ps = np.ones((ncond,))*f107p
aps = np.ones((ncond,))*ap 
efs,ecs = auroral_total_eflux[subset].flatten(),auroral_avg_energy[subset].flatten()
phij = auroral_flux[subset,:]


print(len(idates),ncond)
#subroutine glowssjcond(idates,uts,glats,glongs,f107as,f107s,f107ps,aps,efs,ecs,xncond,phij)
usephij = 1 #Use SSJ fluxes or only Maxwellian
utsin = uts[0]*np.ones((ncond,))
efsin=np.ones_like(efs)*efs[0]
ecsin=np.ones_like(ecs)*ecs[0]
pyglow098.glowssjcond(idates,uts,glats,glons,f107as,f107s,f107ps,aps,efs,ecs,ncond,phij,usephij)
pedcond = pyglow098.cglow.pedcond
hallcond = pyglow098.cglow.hallcond
phitop = pyglow098.cglow.phitop
ener = pyglow098.cglow.ener
z = pyglow098.cglow.zz/1.0e5 #cm->km

print(hallcond.shape)
#pyglow098.glowssjcond(idates,uts,glats,glons,f107as,f107s,f107ps,aps,efs,ecs,ncond,phij,0)
#pedcondmax = pyglow098.cglow.pedcond
#hallcondmax = pyglow098.cglow.hallcond

#print np.nanemean(pedcond - pedcondmax)

#Predict Ovation Prime Conductance
"""
ovation_mlatgrid,ovation_mltgrid,pedgrid,hallgrid = op_cond_estimator.get_conductance(dts[len(dts)/2],hemi='N',
			auroral=True,solar=False,background_p=4.,background_h=4.)

ped_interpolator = ovation_prime.LatLocaltimeInterpolator(ovation_mlatgrid,ovation_mltgrid,pedgrid)
cond_ped_op = ped_interpolator.interpolate(mlats,mlts)

hall_interpolator = ovation_prime.LatLocaltimeInterpolator(ovation_mlatgrid,ovation_mltgrid,hallgrid)
cond_hall_op = hall_interpolator.interpolate(mlats,mlts)
"""
f = pp.figure()
#a0 = f.add_subplot(511)
#a05 = f.add_subplot(512)
a1 = f.add_subplot(311)
a2 = f.add_subplot(312)
a3 = f.add_subplot(313)
"""
z = np.array([80.0,81.5,83.0,84.5,86.0,87.5,89.0,90.5,92.0,93.5,
     95.0,97.5,99.0,100.5,102.0,103.5,105.0,107.5,109.,111.,
     113.,115.,117.5,120.,122.5,125.,128.,131.,135.,139.,
     143.,148.,153.,158.,164.,170.,176.,183.,190.,197.,
     205.,213.,221.,229.,237.,245.,254.,263.,272.,281.,
     290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,
     390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,
     490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,
     590.,600.,610.,620.,630.,640.])
"""
"""
z = np.array([80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
                90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
                100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,
                110.,111.5,113.,114.5,116.,118.,120.,122.,124.,126.,
                128.,130.,132.,134.,136.,138.,140.,142.,144.,146.,
                148.,150.,153.,156.,159.,162.,165.,168.,172.,176.,
                180.,185.,190.,195.,200.,205.,211.,217.,223.,230.,
                237.,244.,252.,260.,268.,276.,284.,292.,300.,309.,
                318.,327.,336.,345.,355.,365.,375.,385.,395.,406.,
                417.,428.,440.,453.,467.,482.,498.,515.,533.,551.,
                570.,590.,610.,630.,650.,670.,690.,710.,730.,750.,
                770.,790.,810.,830.,850.,870.,890.,910.,930.,950.])
"""
from matplotlib.colors import LogNorm

ped = pedcond[:,z<200.]
hall = hallcond[:,z<200.]
z = z[z<200.]

ped[np.logical_not(np.isfinite(ped))]=np.nanmedian(ped.flatten())
hall[np.logical_not(np.isfinite(hall))]=np.nanmedian(hall.flatten())

intped = np.trapz(ped,x=z*1000.) # Conductivity units of S/m?
inthall = np.trapz(hall,x=z*1000.) # Conductivity units of S/m?

"""
#Lazily calculate FAC
mu0 = np.pi*4*1e-7
#Velocity: 7493.082 +- 4.941 m/s for F16, 
#Velocity: 7493.685 +- 5.096 m/s for F17
vperp_jsheet = 7493

#Assume we are completely perpendicular
nTms_to_mAm2 = 1.0e-9*1.0e3 #Unit conversion
ddb = np.diff(dBe)
ddb = np.concatenate((ddb,ddb[-1:]),axis=0) #Make sure we are getting 
#the same number of points as the the original dB
#I removed the -1 from this because we were getting a sign error
jfac = 1./(mu0*vperp_jsheet)*ddb*nTms_to_mAm2

dVy = np.diff(Vy[subset])
dVy = np.concatenate((dVy,dVy[-1:]),axis=0)

#a0.plot(uts,dBe,'g.-')
a0.plot(uts,jfac,'m.-',label='FAC Density')
a0.legend()
a0.set_ylabel('[mA/m^2]')
a05.plot(uts,jfac**2/intped,'k.-',label='Joule Heat (J^2/SigmaP)')
a05.set_ylabel('[mW]')
a05.legend()

#a0.plot(ener,phitop)
#jj=9
#print idates[jj],uts[jj],glats[jj],glons[jj],f107as[jj],f107s[jj],f107ps[jj],aps[jj],efs[jj],ecs[jj] 
#a0.plot(ped[9,:],z,'b.')
#a0.legend()
#a0.set_xscale('log')
#a0.set_yscale('log')
#a0.plot(uts,dVy/np.nanmax(dVy),'g.-')
"""

dmsp_spectrogram.dmsp_spectrogram(uts,auroral_eflux[subset,:],
	chen,datalabel=None,cblims=[1e-7,1e-2],
	ax=a1,fluxunits='Electron\nEnergy Flux\n[mW/m^2]')


T,Z = np.meshgrid(uts,z)

print(np.nanmax(ped),np.nanmin(ped))
mappable = a3.pcolor(T, Z, ped.T, norm=LogNorm(vmin=np.nanmin(ped), vmax=np.nanmax(ped)), cmap='plasma')
cb = pp.colorbar(mappable,ax=a3)
cb.ax.set_ylabel('Pedersen\nConductivity\n [S/m]')
a3.set_ylabel('Altitude\n[km]')

#a2.plot(uts[:-1],np.diff(efs))

#a2.plot(uts,dVy)

a2.plot(uts,intped,'b-',label='DMSP+GLOW Ped')
#a2.plot(uts,cond_ped_op,'b--',label='OvationPyme Ped')
#a2.plot(uts[:-1],np.diff(inthall),label='dPed')
a2.plot(uts,inthall,'r-',label='DMSP+GLOW Hall')
#a2.plot(uts,cond_hall_op,'r--',label='OvationPyme Hall')
#a2.set_ylim([0,120])
#a2.set_ylim([0,1])
a2.legend(ncol=2,loc=0)
a2.set_ylabel('Conductance\n[S]')

cb2 = pp.colorbar(mappable,ax=a2)
cb2.ax.set_visible(False)
#pp.colorbar(mappable,ax=a05)
#pp.colorbar(mappable,ax=a0)

for a in [a1,a2]:
	a.xaxis.set_ticklabels('')
a3.set_xlabel('UT Second\n Apex Lat\nMLT')

satplottools.timepos_ticklabels(a3,uts,mlats,mlts,fs=10)
pp.subplots_adjust(bottom=.2)

f.suptitle('DMSP SSJ driving GLOW 0.98 \n F16, May 29, 2010 %s-%s' % (dts[0].strftime('%H:%M'),dts[-1].strftime('%H:%M')))
f.savefig('/home/liamk/Desktop/glowbeta_ssj_%d%d%d_orbit%d.png' % (year,month,day,orbit),dpi=300.)

pp.show()