import sys, datetime, os, shutil
import numpy as np
import matplotlib.pyplot as pp
sys.path.append('/home/liamk/seshat/glowcond/glow098release/GLOW/')
from geospacepy import special_datetime,dmspcdf_tools,dmsp_spectrogram, satplottools, omnireader
from spacepy import pycdf

import multiprocessing

def mappable_glow(inputs):
	year,month,day,uts,glats,glons,ele_chan_nfluxes,ele_total_fluxes,ele_avg_energies,test,maxwellian = inputs[:]
	import pyglow098

	dt = datetime.datetime(year,month,day)

	doy = special_datetime.datetime2doy(dt)
	jd = special_datetime.datetime2jd(dt)

	f107,ap,f107p,f107a = get_f107_ap(dt)

	idoy = np.floor(doy)
	idate = (year-2000 if year >= 2000 else year-1900)*1000+idoy

	uts = uts.flatten()
	glats = glats.flatten()
	glons = glons.flatten()
	ncond = int(len(uts))
	idates = np.ones((ncond,))*idate
	f107as = np.ones((ncond,))*f107a
	f107s = np.ones((ncond,))*f107
	f107ps = np.ones((ncond,))*f107p
	aps = np.ones((ncond,))*ap
	efs,ecs = ele_total_fluxes.flatten(),ele_avg_energies.flatten()
	phij = ele_chan_nfluxes    

	print("Start GLOW run for %d SSJ points %.1f-%1.f" % (ncond,uts[0]/3600.,uts[-1]/3600.))

	usephij = 0 if maxwellian else 1 #Use SSJ fluxes (1) or only Maxwellian (0)
	#Test means test repeatability of GLOW (i.e. look effect of bug with save and intent(out))
	
	pyglow098.glowssjcond(idates,uts,glats,glons,f107as,f107s,f107ps,aps,efs,ecs,ncond,phij,usephij)

	print("End GLOW run for %d SSJ points %.1f-%1.f" % (ncond,uts[0]/3600.,uts[-1]/3600.))

	pedcond = pyglow098.cglow.pedcond
	hallcond = pyglow098.cglow.hallcond
	phitop = pyglow098.cglow.phitop
	ener = pyglow098.cglow.ener
	z = pyglow098.cglow.zz/1.0e5 #cm->km

	#Apply 200 km ceiling (above which GLOW is not so reliable)
	ped = pedcond[:,z<200.]
	hall = hallcond[:,z<200.]
	z = z[z<200.]

	#ped[np.logical_not(np.isfinite(ped))]=np.nanmedian(ped.flatten())
	#hall[np.logical_not(np.isfinite(hall))]=np.nanmedian(hall.flatten())

	intped = np.trapz(ped,x=z*1000.) # Conductivity units of S/m?
	inthall = np.trapz(hall,x=z*1000.) # Conductivity units of S/m?

	return z,pedcond,hallcond,intped,inthall

#from ovationpyme import ovation_prime

#op_cond_estimator = ovation_prime.ConductanceEstimator(datetime.datetime(2010,5,27,12),datetime.datetime(2010,5,30,12))
def get_f107_ap(dt,silent=False):
	"""
	Get F10.7 for day of interest (noon), and day previous and 81-day centered average
	Also get AP index for day of interest (noon)
	"""
	jd = special_datetime.datetime2jd(dt)

	oi = omnireader.omni_interval(dt-datetime.timedelta(days=41),dt+datetime.timedelta(days=41),'hourly')

	odt = oi['Epoch']
	ojd = special_datetime.datetimearr2jd(odt).flatten()
	f107_81 = oi['F10_INDEX']
	ap_81 = oi['AP_INDEX']
	f107a = np.nanmean(f107_81)
	today = np.logical_and(ojd > np.floor(jd),ojd < np.ceil(jd))
	yesterday = np.logical_and(ojd > np.floor(jd-1),ojd < np.ceil(jd-1))
	f107p = np.nanmean(f107_81[yesterday]) # Previous day f10.7
	f107 = np.nanmean(f107_81[today])
	ap = np.nanmean(ap_81[today])
	if not silent:
		print(dt.strftime('%Y%m%d %H:%M'))
		print('Yesterday F107 %.2f' % (f107p))
		print('Today F107 %.2f'  % (f107))
		print('81-day avg F107: %.2f'  % (f107a))
	return f107,ap,f107p,f107a

def _unused():
	cdfm = dmspcdf_tools.get_cdf(sat,year,month,day,'ssm')
	cdfies = dmspcdf_tools.get_cdf(sat,year,month,day,'ssies')


def get_cond_cdffn(ssj_cdffn,cond_cdf_dir='/home/liamk/code/condcdf'):
	if cond_cdf_dir is None:
		cond_cdffn = ssj_cdfn.splitext()[0]+'_GLOWcond.cdf'
	else:
		if not os.path.exists(cond_cdf_dir):
			os.makedirs(cond_cdf_dir)
			print("Made %s" %(cond_cdf_dir))

		ssj_cdfn_leaf = os.path.split(ssj_cdfn)[-1]
		cond_cdffn = os.path.join(cond_cdf_dir,ssj_cdffn.splitext()[0]+'_GLOWcond.cdf')
	return cond_cdffn

def run_ssj_glow(sat,year,month,day,minlat=50.,create_conductance_cdf=False,clobber=True,silent=False):
	"""
	Calculate the conductivities along DMSP satellite track
	for a particular day
	"""

	n_processers = 1
	
	test = False # Test for intent(out) bug
	maxwellian = False # Use maxwellian spectrum instead of DMSP SSJ

	dt = datetime.datetime(year,month,day,12,0,0)

	cdfn = dmspcdf_tools.get_cdf(sat,year,month,day,'ssj',return_file=True)
	with pycdf.CDF(cdfn) as cdf:

		dts = cdf['Epoch'][:]
		uts = special_datetime.datetimearr2sod(cdf['Epoch'][:]).flatten()
		n_times = uts.shape[0]
		auroral_region = cdf['AURORAL_REGION'][:].flatten()
		chen = cdf['CHANNEL_ENERGIES'][:]
		glats = cdf['SC_GEOCENTRIC_LAT'][:].flatten()
		mlats = cdf['SC_APEX_LAT'][:].flatten()
		glons = cdf['SC_GEOCENTRIC_LON'][:].flatten()

		nflux = cdf['ELE_DIFF_ENERGY_FLUX'][:]/chen
		avg_energy = cdf['ELE_AVG_ENERGY'][:].flatten()
		total_eflux = cdf['ELE_TOTAL_ENERGY_FLUX'][:].flatten()

		total_eflux *= 1.6e-12 #eV/cm/s/sr -> mW/m^2
		
		min_mlat = 50.
		
		#
		#Pack up the inputs
		#
		inputs,results,subset_masks = [],[],[]
		for i_subset in range(n_processers):
			subset_length = n_times/n_processers
			uts_start,uts_end = i_subset*subset_length,(i_subset+1)*subset_length
			in_lats = np.abs(mlats) > min_mlat 
			subset_mask = np.logical_and(uts >= uts_start,uts < uts_end)
			subset_masks.append(subset_mask)
			uts_in = uts[subset_mask]
			glats_in = glats[subset_mask]
			glons_in = glons[subset_mask]
			nflux_in = nflux[subset_mask,:]
			avg_energy_in = avg_energy[subset_mask]
			total_eflux_in = total_eflux[subset_mask]
			inputs.append((year,month,day,uts_in,glats_in,glons_in,nflux_in,avg_energy_in,total_eflux_in,test,maxwellian))
		
		#
		#Pass to the pool
		#

		if n_processers > 1:
			pool = multiprocessing.Pool(n_processers)
			results = pool.map(mappable_glow,inputs)
		else:
			results = [mappable_glow(inputs[0])]

		#
		#Unpack the results
		#
		for i_subset in range(n_processers):
			subset_mask = subset_masks[i_subset]
			z,subset_pedcond,subset_hallcond,subset_intped,subset_inthall = result[i_subset][:]
			if i_subset == 0:
				#Initialize the things
				pedcond,hallcond = np.zeros((n_times,len(z.flatten()))),np.zeros((n_times,len(z.flatten())))
				pedcond.fill(np.nan)
				hallcond.fill(np.nan)
				intped,inthall = np.zeros_like(glats),np.zeros_like(glats)
				intped.fill(np.nan)
				inthall.fill(np.nan)

			#Store the conductances
			pedcond[subset_mask,:]=subset_pedcond
			hallcond[subset_mask,:]=subset_pedcond
			intped[subset_mask] = subset_intped
			inthall[subset_mask] = subset_inthall

		if create_conductance_cdf:
			cdfn_cond = get_cond_cdffn(cdfn)
			if os.path.exists(cdfn_cond):
				if clobber:
					os.remove(cdfn_cond)
					shutil.copyfile(cdfn,cdfn_cond)
				else:
					raise IOError('Conductance CDF %s exists and clobber is False' % (cdfn_cond))
			#Copy the SSJ file
			shutil.copyfile(cdfn,cdfn_cond) 

			with pycdf.CDF(cdfn_cond) as cdf_cond:
				cdf_cond.readonly(False) # Modify the file
				cdf_cond['CONDUCTIVITY_ALTITUDES']=z
				cdf_cond['CONDUCTIVITY_ALTITUDES'].attrs['UNITS']='km'
				
				cdf_cond['PEDERSEN_CONDUCTANCE']=intped
				cdf_cond['PEDERSEN_CONDUCTANCE'].attrs['UNITS']='S'
				
				cdf_cond['HALL_CONDUCTANCE']=inthall
				cdf_cond['HALL_CONDUCTANCE'].attrs['UNITS']='S'
				
				cdf_cond['PEDERSEN_CONDUCTIVITY']=pedcond
				cdf_cond['PEDERSEN_CONDUCTIVITY'].attrs['UNITS']='S/m'
				
				cdf_cond['HALL_CONDUCTIVITY']=hallcond
				cdf_cond['HALL_CONDUCTIVITY'].attrs['UNITS']='S/m'
				
	 
def plot_orbit_cond(sat,year,month,day,orbit,minlat=50.,plotdir=None):

	#    
	#Prepare Plots
	#
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

	#
	# Get conductance CDF
	#
	cdfn = dmspcdf_tools.get_cdf(sat,year,month,day,'ssj',return_file=True)
	cdfn_cond = get_cond_cdffn(cdfn)

	#
	#   Figure Filename
	#
	if plotdir is None:
		plotdir = '/home/liamk/code/glowcond/%s' % (cdfn_cond_leaf)
	if not os.path.exists(plotdir):
		os.makedirs(plotdir)

	figfn = os.path.join(plotdir,'%s_%s_%d.png' % (cdfn_cond_leaf,hemi,np.abs(orbit)))
	hemi = 'N' if np.sign(orbit)==1 else 'S'

	with pycdf.CDF(cdfn_cond) as cdf:

		orbit_index = cdf['ORBIT_INDEX'][:].flatten()
		mlats = cdf['SC_APEX_LAT'][:].flatten()
		
		in_orbit = orbit_index == orbit
		subset = np.logical_and(in_orbit,np.abs(mlats)>minlat)
   
		dts = cdf['Epoch'][:][subset]
		uts = special_datetime.datetimearr2sod(cdf['Epoch'][:]).flatten()[subset]
		auroral_region = cdf['AURORAL_REGION'][:].flatten()[subset]
		chen = cdf['CHANNEL_ENERGIES'][:]
		glats = cdf['SC_GEOCENTRIC_LAT'][:].flatten()[subset]
		glons = cdf['SC_GEOCENTRIC_LON'][:].flatten()[subset]
		mlats = cdf['SC_APEX_LAT'][:].flatten()[subset]
		mlts = cdf['SC_APEX_MLT'][:].flatten()[subset]
		
		eflux = cdf['ELE_DIFF_ENERGY_FLUX'][:][subset,:]
		avg_energy = cdf['ELE_AVG_ENERGY'][:].flatten()[subset]
		total_eflux = cdf['ELE_TOTAL_ENERGY_FLUX'][:].flatten()[subset]

		total_eflux *= 1.6e-12 #eV/cm/s/sr -> mW/m^2
		eflux *= 1.6e-12

		z = cdf['CONDUCTIVITY_ALTITUDES'][:]
		ped = cdf['PEDERSEN_CONDCUCTIVITY'][:][subset]
		hall = cdf['HALL_CONDCUCTIVITY'][:][subset]
		intped = cdf['PEDERSEN_CONDUCTANCE'][:][subset]
		inthall = cdf['HALL_CONDUCTANCE'][:][subset]
		
		#Plot Electron energy flux
		dmsp_spectrogram.dmsp_spectrogram(uts,eflux,
			chen,datalabel=None,cblims=[1e-7,1e-2],
			ax=a1,ax_cb=a11,fluxunits='Electron\nEnergy Flux\n[mW/m^2]')

		#Plot Pedersen Conductance
		a2.plot(uts,intped,'r.-',label='DMSP+GLOW')
		#a2.plot(uts[:-1],np.diff(inthall),label='dPed')
		#a2.set_ylim([0,120])
		#a2.set_ylim([0,1])
		a2.legend(ncol=2,loc=0)
		a2.set_ylabel('Pedersen\n Conductance\n[S]')

		satplottools.draw_dialplot(a22)
		x,y = satplottools.latlt2cart(mlats,mlts,hemi)
		a22.plot(x,y,'k.')
		a22.text(x[-1],y[-1],'End')

		#Plot Hall Conductance
		a3.plot(uts,inthall,'r.-',label='DMSP+GLOW')
		a3.legend(ncol=2,loc=0)
		a3.set_ylabel('Hall\n Conductance\n[S]')

		#Draw conductivity as a pcolor plot with log-scaled color scale
		T,Z = np.meshgrid(uts,z)
		mappable = a4.pcolor(T, Z, ped.T, norm=LogNorm(vmin=np.nanmin(ped), vmax=np.nanmax(ped)), cmap='plasma')
		cb = pp.colorbar(mappable,cax=a44)
		cb.ax.set_ylabel('Pedersen\nConductivity\n [S/m]')
		a4.set_ylabel('Altitude\n[km]')

		f.savefig(figfn)

if __name__ == '__main__':
	
	sat,year,month,day = 16,2011,5,28
	run_ssj_glow(sat,year,month,day,create_conductance_cdf=True,clobber=True,silent=False)
	for orbit in range(1,15):
		plot_orbit_cond(sat,year,month,day,1*orbit,minlat=50.,plotdir=None)
		plot_orbit_cond(sat,year,month,day,-1*orbit,minlat=50.,plotdir=None)

