import sys, datetime
import numpy as np
import matplotlib.pyplot as pp
sys.path.append('/home/liamk/seshat/glowcond/glow098beta/')
import pyglow098beta
from geospacepy import special_datetime,dmspcdf_tools,dmsp_spectrogram, satplottools
from ovationpyme import ovation_prime

op_cond_estimator = ovation_prime.ConductanceEstimator(datetime.datetime(2010,5,27,12),datetime.datetime(2010,5,30,12))

year,month,day = 2010,5,29
sat = 16

cdf = dmspcdf_tools.get_cdf(sat,year,month,day,'ssj')
cdfm = dmspcdf_tools.get_cdf(sat,year,month,day,'ssm')
cdfies = dmspcdf_tools.get_cdf(sat,year,month,day,'ssies')

for orbit in range(14):

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
	#subset = np.arange(6700,7900)
	subset = np.flatnonzero(oi[auroral]==orbit)
	
	if len(subset) < 10:
		print('No data for orbit %d' % (orbit))
		continue
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
	idate,f107a,f107,f107p,ap = 10149,70,70,70,20
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


	print len(idates),ncond
	#subroutine glowssjcond(idates,uts,glats,glongs,f107as,f107s,f107ps,aps,efs,ecs,xncond,phij)
	usephij = 1 #Use SSJ fluxes or only Maxwellian
	utsin = uts[0]*np.ones((ncond,))
	efsin=np.ones_like(efs)*efs[0]
	ecsin=np.ones_like(ecs)*ecs[0]
	pyglow098beta.glowssjcond(idates,uts,glats,glons,f107as,f107s,f107ps,aps,efs,ecs,ncond,phij,usephij)
	pedcond = pyglow098beta.cglow.pedcond
	hallcond = pyglow098beta.cglow.hallcond
	phitop = pyglow098beta.cglow.phitop
	ener = pyglow098beta.cglow.ener

	print hallcond.shape
	#pyglow098beta.glowssjcond(idates,uts,glats,glons,f107as,f107s,f107ps,aps,efs,ecs,ncond,phij,0)
	#pedcondmax = pyglow098beta.cglow.pedcond
	#hallcondmax = pyglow098beta.cglow.hallcond

	#print np.nanemean(pedcond - pedcondmax)

	#Predict Ovation Prime Conductance

	ovation_mlatgrid,ovation_mltgrid,pedgrid,hallgrid = op_cond_estimator.get_conductance(dts[len(dts)/2],hemi='N',
				auroral=True,solar=False,background_p=4.,background_h=4.)

	ped_interpolator = ovation_prime.LatLocaltimeInterpolator(ovation_mlatgrid,ovation_mltgrid,pedgrid)
	cond_ped_op = ped_interpolator.interpolate(mlats,mlts)

	hall_interpolator = ovation_prime.LatLocaltimeInterpolator(ovation_mlatgrid,ovation_mltgrid,hallgrid)
	cond_hall_op = hall_interpolator.interpolate(mlats,mlts)

	#Prepare Plots
	#There will be 3 line plots
	import matplotlib.gridspec as gridspec
	from matplotlib.colors import LogNorm

	gs = gridspec.GridSpec(4, 11)
	f = pp.figure()
	#a0 = f.add_subplot(511)
	#a05 = f.add_subplot(512)
	split = 9
	cbwidth = 1
	a1 = pp.subplot(gs[0,:split])
	a11 = pp.subplot(gs[0,split:split+cbwidth])
	a2 = pp.subplot(gs[1,:split])
	a22 = pp.subplot(gs[1:2,split:])
	a3 = pp.subplot(gs[2,:split])
	#a33 = pp.subplot(gs[2,split:])
	a4 = pp.subplot(gs[3,:split])
	a44 = pp.subplot(gs[3,split:split+cbwidth])

	# Heights from GLOW (will make an output in the future (maybe already is part of cglow?))
	z = np.array([80.0,81.5,83.0,84.5,86.0,87.5,89.0,90.5,92.0,93.5,
	     95.0,97.5,99.0,100.5,102.0,103.5,105.0,107.5,109.,111.,
	     113.,115.,117.5,120.,122.5,125.,128.,131.,135.,139.,
	     143.,148.,153.,158.,164.,170.,176.,183.,190.,197.,
	     205.,213.,221.,229.,237.,245.,254.,263.,272.,281.,
	     290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,
	     390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,
	     490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,
	     590.,600.,610.,620.,630.,640.])

	#Limit to heights at which GLOW is trustworthy
	ped = pedcond[:,z<200.]
	hall = hallcond[:,z<200.]
	z = z[z<200.]

	#Fill in errors
	ped[np.logical_not(np.isfinite(ped))]=np.nanmedian(ped.flatten())
	hall[np.logical_not(np.isfinite(hall))]=np.nanmedian(hall.flatten())

	#Integrate conductance
	intped = np.trapz(ped,x=z*1000.) # Conductivity units of S/m?
	inthall = np.trapz(hall,x=z*1000.) # Conductivity units of S/m?

	#Plot Electron energy flux
	dmsp_spectrogram.dmsp_spectrogram(uts,auroral_eflux[subset,:],
		chen,datalabel=None,cblims=[1e-7,1e-2],
		ax=a1,ax_cb=a11,fluxunits='Electron\nEnergy Flux\n[mW/m^2]')

	#Plot Pedersen Conductance
	a2.plot(uts,intped,'r.-',label='DMSP+GLOW')
	a2.plot(uts,cond_ped_op,'b.-',label='OP2010+Robinson')
	#a2.plot(uts[:-1],np.diff(inthall),label='dPed')
	#a2.set_ylim([0,120])
	#a2.set_ylim([0,1])
	a2.legend(ncol=2,loc=0)
	a2.set_ylabel('Pedersen\n Conductance\n[S]')

	satplottools.draw_dialplot(a22)
	x,y = satplottools.latlt2cart(mlats,mlts,'N')
	a22.plot(x,y,'k.')
	a22.text(x[-1],y[-1],'End')

	#Plot Hall Conductance
	a3.plot(uts,inthall,'r.-',label='DMSP+GLOW')
	a3.plot(uts,cond_hall_op,'b.-',label='OP2010+Robinson')
	a3.legend(ncol=2,loc=0)
	a3.set_ylabel('Hall\n Conductance\n[S]')

	#Draw conductivity as a pcolor plot with log-scaled color scale
	T,Z = np.meshgrid(uts,z)
	mappable = a4.pcolor(T, Z, ped.T, norm=LogNorm(vmin=np.nanmin(ped), vmax=np.nanmax(ped)), cmap='plasma')
	cb = pp.colorbar(mappable,cax=a44)
	cb.ax.set_ylabel('Pedersen\nConductivity\n [S/m]')
	a4.set_ylabel('Altitude\n[km]')

	#a2.plot(uts[:-1],np.diff(efs))

	#a2.plot(uts,dVy)

	#pp.colorbar(mappable,ax=a05)
	#pp.colorbar(mappable,ax=a0)

	for a in [a1,a2,a3]:
		a.xaxis.set_ticklabels('')
	a4.set_xlabel('UT Second\n Apex Lat\nMLT')

	satplottools.timepos_ticklabels(a4,uts,mlats,mlts,fs=10)
	pp.subplots_adjust(bottom=.2)

	f.suptitle('DMSP SSJ driving GLOW 0.98 \n F16, May 29, 2010 %s-%s' % (dts[0].strftime('%H:%M'),dts[-1].strftime('%H:%M')))
	f.savefig('/home/liamk/Desktop/glowbeta_ssj_ovation_%d%d%d_orbit%d.png' % (year,month,day,orbit),dpi=300.)

	pp.pause(1)