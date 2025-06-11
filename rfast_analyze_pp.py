# import statements
import emcee
import corner
import h5py
import sys
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
from astropy.table       import Table, Column, MaskedColumn
from astropy.io          import ascii
from rfast_routines      import spectral_grid
from rfast_routines      import gen_spec
from rfast_routines      import kernel_convol
from rfast_routines      import gen_spec_grid
from rfast_routines      import inputs
from rfast_routines      import init
from rfast_routines      import init_3d
from rfast_atm_routines  import set_gas_info
from rfast_atm_routines  import setup_atm
from rfast_atm_routines  import mmr2vmr
from rfast_atm_routines  import vmr2mmr
from rfast_opac_routines import opacities_info
from rfast_user_models   import surfalb
from rfast_user_models   import cloud_optprops
from rfast_user_models   import cloud_struct

# simple routine for importing emcee chain from h5 file
def reademceeh5(fn,nburn,thin,flatten=False):

  # open file, important data
  hf       = h5py.File(fn,'r')
  grps     = [item for item in hf['mcmc'].values()]

  # extract samples chain and log-likelihood, remove burn-in
  if (nburn >= 0):
    samples  = grps[1][nburn:,:,:]
    lnprob   = grps[2][nburn:,:]
  else:
    samples  = grps[1][nburn:,:,:]
    lnprob   = grps[2][nburn:,:]

  # thin
  samples  = samples[0::thin,:,:]
  lnprob   = lnprob[0::thin,:]

  # flatten
  if flatten:
    samples  = samples.reshape(samples.shape[0]*samples.shape[1],samples.shape[2])
    lnprob   = lnprob.reshape(lnprob.shape[0]*lnprob.shape[1])

  # close h5 file
  hf.close()

  return samples,lnprob

# get input script filename
if len(sys.argv) >= 2:
  filename_scr = sys.argv[1]    # if script name provided at command line
else:
  filename_scr = input("rfast inputs script filename: ") # otherwise ask for filename
  sys.argv.append(filename_scr) # poor practice, but prevents bug when importing Fx

# obtain input parameters from script
fnr,fnn,fns,dirout,Nlev,pmin,pmax,bg,\
species_r,f0,rdgas,fnatm,skpatm,colr,colpr,psclr,mmri,\
tpars,rdtmp,fntmp,skptmp,colt,colpt,psclt,\
species_l,species_c,\
lams,laml,res,regrid,smpl,opdir,\
Rp,Mp,gp,a,Apars,em,\
cld,phfc,opars,cpars,lamc0,fc,\
ray,ref,sct,fixp,pf,fixt,tf,p10,fp10,\
src,\
alpha,ntg,\
Ts,Rs,\
ntype,snr0,lam0,rnd,\
clr,fmin,mmrr,nwalkers,nstep,nburn,thin,restart,progress = inputs(filename_scr)


# input data filename
fn_dat = fnn + '.dat'

# set info for all radiatively active gases
Ngas,gasid,mmw0,ray0,nu0,mb,rayb = set_gas_info(bg)

# get initial gas mixing ratios
p,t,t0,z,grav,f,fb,m,nu = setup_atm(Nlev,Ngas,gasid,mmw0,pmin,pmax,
                                    tpars,rdtmp,fntmp,skptmp,colt,colpt,psclt,
                                    species_r,f0,rdgas,fnatm,skpatm,colr,colpr,psclr,
                                    mmri,mb,Mp,Rp,p10,fp10,src,ref,nu0)

# convert between mmr and vmr, if needed
#if (mmri != mmrr):
 # if mmri: # convert input mmr to vmr
 #   f,fb,f0 = mmr2vmr(mmw0,gasid,species_r,m,mb,f0,fb,f)
  #else: # otherwise convert input vmr to mmr
  #  f,fb,f0 = vmr2mmr(mmw0,gasid,species_r,m,mb,f0,fb,f)


# read input data
data        = ascii.read(dirout+fn_dat,data_start=1,delimiter='|') #same data as retrieve
lam         = data['col2'][:]
dlam        = data['col3'][:]
albedo      = data['col4'][:] ##added january 11th
dat         = data['col6'][:] #flux ratio
err         = data['col7'][:] #should be albedo error now

# save input radius for thermal emission case
Rpi = Rp

# generate wavelength grids
Nres             = 3 # no. of res elements to extend beyond grid edges to avoid edge sensitivity (2--4)
if regrid:
  lam_hr,dlam_hr = gen_spec_grid(lams,laml,np.float_(res)*smpl,Nres=np.rint(Nres*smpl))
else:
  x_low = min(0.01,min(lam)-dlam[0]*Nres) # note: prevent min wavelength of 0 um
  x_hgh = max(lam)+dlam[-1]*Nres
  lam_hr,dlam_hr = spectral_grid(x_low,x_hgh,res=lam/dlam*smpl,lamr=lam)

# assign photometric vs. spectroscopic points
mode           = np.copy(lam_hr)
mode[:]        = 1

# inform user of key opacities information
opacities_info(opdir)

# initialize opacities and convolution kernels
sigma_interp,cia_interp,ncia,ciaid,kern = init(lam,dlam,lam_hr,species_l,species_c,opdir,pf,tf,mode=mode)

# surface albedo model
As = surfalb(Apars,lam_hr)[0]

# cloud optical properties: asymmetry parameter, single scattering albedo, extinction efficiency
gc,wc,Qc = cloud_optprops(opars,cld,opdir,lam_hr)

# cloud vertical structure model
dtauc0,atm = cloud_struct(cpars,cld,p,t,z,grav,f,fb,m)
p,t,z,grav,f,fb,m = atm

# initialize disk integration quantities
threeD   = init_3d(src,ntg)

# package parameters for user-defined routines
tiso = tpars[0]
A0 = Apars[0]
A1 = Apars[1]
A2 = Apars[2]
A3 = Apars[3]
A4 = Apars[4]
# no parameters for cloud optical properties model
pt,dpc,tauc0 = cpars

# parameter names
lfN2,lfO2,lfH2O,lfO3,lfCO2,lfCO,lfCH4 = np.log10(f0[species_r=='n2'])[0],np.log10(f0[species_r=='o2'])[0],np.log10(f0[species_r=='h2o'])[0],np.log10(f0[species_r=='o3'])[0],np.log10(f0[species_r=='co2'])[0],np.log10(f0[species_r=='co'])[0],np.log10(f0[species_r=='ch4'])[0]

##added january 11th from  the other code

lpN2  = np.log10((10**lfN2)*pmax)
lpO2  = np.log10((10**lfO2)*pmax)
lpH2O = np.log10((10**lfH2O)*pmax)
lpO3  = np.log10((10**lfO3)*pmax)
lpCO2 = np.log10((10**lfCO2)*pmax)
lpCO  = np.log10((10**lfCO)*pmax)
lpCH4 = np.log10((10**lfCH4)*pmax)
pmax  = 10**lpN2+10**lpO2+10**lpH2O+10**lpO3+10**lpCO2+10**lpCO+10**lpCH4
lpmax = np.log10(pmax)
lfN2  = np.log10(10**lpN2/pmax)
lfO2  = np.log10(10**lpO2/pmax)
lfH2O = np.log10(10**lpH2O/pmax)
lfO3  = np.log10(10**lpO3/pmax)
lfCO2 = np.log10(10**lpCO2/pmax)
lfCO  = np.log10(10**lpCO/pmax)
lfCH4 = np.log10(10**lpCH4/pmax)

fN2,fO2,fH2O,fO3,fCO2,fCO,fCH4 = 10**lfN2,10**lfO2,10**lfH2O,10**lfO3,10**lfCO2,10**lfCO,10**lfCH4


lRp,lMp,lfc,lA0,lA1,lA2,lA3,lA4,ldpc,lpt,ltauc0 = np.log10(Rp),np.log10(Mp),np.log10(fc),np.log10(A0),np.log10(A1),np.log10(A2), np.log10(A3),np.log10(A4),np.log10(dpc),np.log10(pt),np.log10(tauc0)

pp_names  = [r"$\log\,$"+r"$p_{\rm N2}$",r"$\log\,$"+r"$p_{\rm O2}$",r"$\log\,$"+r"$p_{\rm H2O}$",r"$\log\,$"+r"$p_{\rm O3}$",r"$\log\,$"+r"$p_{\rm CO2}$",r"$\log\,$"+r"$p_{\rm CO}$",r"$\log\,$"+r"$p_{\rm CH4}$",r"$\log\,$"+r"$f_{Wbasalt}$",r"$\log\,$"+r"$f_{Kaolinite}$",r"$\log\,$"+r"$f_{Wgranite}$", r"$\log\,$"+r"$f_{Desert}$", r"$\log\,$"+r"$f_{Ice}$",r"$\log\,$"+r"$R_{\rm p}$",r"$\log\,$"+r"$M_{\rm p}$",r"$\log\,$"+r"$\Delta p_{\rm c}$",r"$\log\,$"+r"$p_{\rm t}$",r"$\log\,$"+r"$\tau_{\rm c}$",r"$\log\,$"+r"$f_{\rm c}$"]
vmr_names  = [r"$\log\,$"+r"$f_{\rm N2}$",r"$\log\,$"+r"$f_{\rm O2}$",r"$\log\,$"+r"$f_{\rm H2O}$",r"$\log\,$"+r"$f_{\rm O3}$",r"$\log\,$"+r"$f_{\rm CO2}$",r"$\log\,$"+r"$f_{\rm CO}$",r"$\log\,$"+r"$f_{\rm CH4}$",r"$\log\,$"+r"$A_{\rm 0}$",r"$\log\,$"+r"$R_{\rm p}$",r"$\log\,$"+r"$M_{\rm p}$",r"$\log\,$"+r"$\Delta p_{\rm c}$",r"$\log\,$"+r"$p_{\rm t}$",r"$\log\,$"+r"$\tau_{\rm c}$",r"$\log\,$"+r"$f_{\rm c}$",r"$\log\,$"+r"$p_{0}$"]


pp_truths = [lpN2,lpO2,lpH2O,lpO3,lpCO2,lpCO,lpCH4,lA0,lA1,lA2,lA3,lA4,lRp,lMp,ldpc,lpt,ltauc0,lfc]
vmr_truths = [lfN2,lfO2,lfH2O,lfO3,lfCO2,lfCO,lfCH4,lA0,lRp,lMp,ldpc,lpt,ltauc0,lfc,lpmax] #lpmax is listed as p0 in vmr_names


ndim   = len(pp_names)

# import chain data
samples,lnprob = reademceeh5(dirout+fnr+'.h5',nburn,thin)

#import chain data WITHOUT burnin/thin
samples_nbt,lnprob = reademceeh5(dirout+fnr+'.h5',0,1)
		

# print reduced chi-squared
lnp_max = np.amax(lnprob)
pos_max = np.where(lnprob == lnp_max)
print("Reduced chi-squared: ",-2*lnp_max/(albedo.shape[0]-ndim)) #change to albedo?

# relevant sizes
nstep    = samples.shape[0]
nwalkers = samples.shape[1]
ndim     = samples.shape[2]

# if doing center-log ratio, transform back to mixing ratio
if clr:
  gind = []
  xi   = samples[:,:,gind]
  clrs = np.sum(np.exp(xi),axis=2) + np.exp(-np.sum(xi,axis=2))
  clrs = np.repeat(clrs[:,:,np.newaxis], len(gind), axis=2)
  samples[:,:,gind] = np.log10(np.divide(np.exp(samples[:,:,gind]),clrs))

#plotting scheme requires co-located folders called 'walkers', 'corner_plots', 'Land_fractions', and 'spectra'

# plot the walker positions in each step
fig, axes = plt.subplots(ndim, 1, figsize=(8, 4 * ndim), tight_layout=True)
for i in range(ndim):
  for j in range(0,nwalkers):
    axes[i].plot(samples[:,j,i],color="black",linewidth=0.5)
    axes[i].set_ylabel(str(pp_names[i]))
    axes[i].set_xlabel('Step')
plt.savefig('walkers/walkers.png',format='png')
plt.close()

# plot the walker positions in each step WITHOUT burnin/thin added feb 6
fig, axes = plt.subplots(ndim, 1, figsize=(8, 4 * ndim), tight_layout=True)
for i in range(ndim):
  for j in range(0,nwalkers):
    axes[i].plot(samples_nbt[:,j,i],color="black",linewidth=0.5)
    axes[i].set_ylabel(str(pp_names[i]))
    axes[i].set_xlabel('Step')
plt.savefig('walkers/walkers_no_burnin_thin.png',format='png', bbox_inches='tight')
plt.close()


#ocean and total land walkers

fig,ax = plt.subplots()
for r in range(nwalkers):
	ocean_walk = 1 - 10**samples_nbt[:,r,7] - 10**samples_nbt[:,r,8] - 10**samples_nbt[:,r,9] - 10**samples_nbt[:,r,10] - 10**samples_nbt[:,r,11]
	ax.plot(ocean_walk, color = 'royalblue', linewidth = 0.5)
	ax.set_ylabel('Ocean fraction')
	ax.set_xlabel('Step')
plt.savefig('walkers/ocean_walkers.png', bbox_inches = 'tight', dpi = 300)
plt.close()	 

fig,ax = plt.subplots()
for r in range(nwalkers):
	land_walk = 10**samples_nbt[:,r,7] + 10**samples_nbt[:,r,8] + 10**samples_nbt[:,r,9] + 10**samples_nbt[:,r,10]
	ax.plot(land_walk, color = 'mediumseagreen', linewidth = 0.5)
	ax.set_ylabel('Land fraction')
	ax.set_xlabel('Step')
plt.savefig('walkers/land_walkers.png', bbox_inches = 'tight', dpi = 300)
plt.close()	 


# plot the corner plot
fig = corner.corner(samples.reshape((-1,ndim)), quantiles=[0.16, 0.5, 0.84], show_titles=True,
                    color='black', labels=pp_names, truths=pp_truths)
fig.savefig('corner_plots/corner_pp.png',format='png',bbox_inches='tight')
plt.close(fig)

#extract just the land fraction parameter distribution(s)

flat = samples.reshape((-1,ndim))

frac_label = pp_names[7]
frac_truth = pp_truths[7]

land1 = 10**(flat[:,7]).reshape(-1,1)
land2 = 10**(flat[:,8]).reshape(-1,1)
land3 = 10**(flat[:,9]).reshape(-1,1)
land4 = 10**(flat[:,10]).reshape(-1,1)
ice = 10**(flat[:,11]).reshape(-1,1)

fig = corner.corner(land1, quantiles = [0.16, 0.5, 0.84], show_titles = True,
                    color = 'black', labels = ['W basalt'], truths = [0.05])
fig.savefig('Land_fractions/w_basalt.png',format='png',bbox_inches='tight', dpi = 300)
fig.tight_layout()
plt.close(fig)

#for the second land fraction parameter, as a test
fig = corner.corner(land2, quantiles = [0.16, 0.5, 0.84], show_titles = True,
                    color = 'black', labels = ['Kaolinite'], truths = [0.05])
fig.savefig('Land_fractions/kaolinite.png',format='png',bbox_inches='tight', dpi = 300)
fig.tight_layout()
plt.close(fig)

fig = corner.corner(land3, quantiles = [0.16, 0.5, 0.84], show_titles = True,
                    color = 'black', labels = ['W granite'], truths = [0.1])
fig.savefig('Land_fractions/w_granite.png',format='png',bbox_inches='tight', dpi = 300)
fig.tight_layout()
plt.close(fig)

fig = corner.corner(land4, quantiles = [0.16, 0.5, 0.84], show_titles = True,
                    color = 'black', labels = ['Desert'], truths = [0.1])
fig.savefig('Land_fractions/desert.png',format='png',bbox_inches='tight', dpi = 300)
fig.tight_layout()
plt.close(fig)

#total land fraction
summ = land1 + land2 + land3 + land4

fig = corner.corner(summ, quantiles = [0.16, 0.5, 0.84], show_titles = True, color = 'black', labels = ['Total land'], truths = [0.3])
fig.savefig('Land_fractions/total_land.png',format='png',bbox_inches='tight', dpi = 300)
fig.tight_layout()
plt.close(fig)

#for ice
fig = corner.corner(ice, quantiles = [0.16, 0.5, 0.84], show_titles = True,
                    color = 'black', labels = ['Ice'], truths = [0.1])
fig.savefig('Land_fractions/ice.png',format='png',bbox_inches='tight', dpi = 300)
fig.tight_layout()
plt.close(fig)

ocean = 1 - summ - ice

data_out = Table([summ[:,0],ice[:,0],ocean[:,0]], names = ['land','ice','ocean'])
ascii.write(data_out, '../1.8um.csv', format = 'fixed_width', overwrite = True)

surfaces = []
for a,b,c,d,e in zip(land1,land2,land3,land4,ice):
	pars = a[0],b[0],c[0],d[0],e[0]

	surfaces.append(surfalb(pars,lam)[0])
	plt.plot(lam, surfalb(pars,lam)[0])

plt.ylabel('Albedo')
plt.xlabel('Wavelength ($\mu$m)')
plt.grid(alpha = 0.5)

plt.legend(facecolor = 'black', labelcolor = 'linecolor')
plt.savefig('surfaces.png', dpi = 300)
plt.close()

import scipy.stats

surfaces = np.array(surfaces)
surfaces = surfaces.T


ci_bot = []
ci_lower = []
ci_mid = []
ci_upper = []
ci_max = []
for t in surfaces:
    ci_bot.append(scipy.stats.scoreatpercentile(t, 0, interpolation_method = 'fraction', axis = 0))
    ci_lower.append(scipy.stats.scoreatpercentile(t, 2.5, interpolation_method = 'fraction', axis = 0))
    ci_mid.append(scipy.stats.scoreatpercentile(t, 50, interpolation_method = 'fraction', axis = 0))
    ci_upper.append(scipy.stats.scoreatpercentile(t, 97.5, interpolation_method = 'fraction', axis = 0))
    ci_max.append(scipy.stats.scoreatpercentile(t, 100, interpolation_method = 'fraction', axis = 0))

fig, ax = plt.subplots()
ax.fill_between(lam, ci_upper, ci_lower, alpha = 0.6, label = '95% CI')
ax.plot(lam, ci_mid, color = 'blue')
ax.set_ylabel('Albedo')
ax.set_xlabel('Wavelength ($\mu$m)')
ax.grid(alpha = 0.5)
ax.legend()
plt.savefig('95_ci.png',dpi = 300)
plt.close()  

data_out = Table([lam, ci_mid, ci_upper, ci_lower], names = ['lam', 'ci_mid', 'ci_upper', 'ci_lower'])
ascii.write(data_out, 'surface.csv', format = 'fixed_width', overwrite = True)

# plot best-fit model and residuals
gp = -1 # reverts to using Mp if gp not retrieved

# get best-fit parameters

lfN2,lfO2,lfH2O,lfO3,lfCO2,lfCO,lfCH4,lA0,lA1,lA2,lA3,lA4,lRp,lMp,dpc,lpt,ltauc0,lfc = samples_nbt[pos_max][0]
fN2,fO2,fH2O,fO3,fCO2,fCO,fCH4,Rp,Mp,fc,A0,A1,A2,A3,A4,dpc,pt,tauc0 = 10**(lfN2),10**(lfO2),10**(lfH2O),10**(lfO3),10**(lfCO2),10**(lfCO),10**(lfCH4),10**(lRp),10**(lMp),10**(lfc),10**(lA0),10**(lA1),10**(lA2),10**(lA3),10**(lA4),10**(ldpc),10**(lpt),10**(ltauc0)
f0[species_r=='N2'],f0[species_r=='o2'],f0[species_r=='h2o'],f0[species_r=='o3'],f0[species_r=='co2'],f0[species_r=='co'],f0[species_r=='ch4'] = fN2,fO2,fH2O,fO3,fCO2,fCO,fCH4

# package parameters for user-defined routines
tpars = tiso
Apars = A0,A1,A2,A3,A4
# no parameters for cloud optical properties model
cpars = pt,dpc,tauc0

x0 = f0,pmax,Rp,Mp,gp,Apars,tpars,opars,cpars,fc,a,alpha,mb,rayb #f0 here may not be mixing ratios...
#is samples coming from partial pressures in retrieve?

y  = lam_hr,pmin,Ngas,mmw0,nu0,em,Ts,Rs,ray0,gasid,ncia,ciaid,species_l,species_c,\
     Rpi,sigma_interp,cia_interp,kern,ray,ref,cld,sct,phfc,fixp,pf,fixt,tf,mmrr,\
     p10,fp10,src,threeD,fntmp,skptmp,colt,colpt,psclt,species_r,fnatm,skpatm,\
     colr,colpr,psclr,Nlev

# determine correct label for y axis
if (src == 'diff' or src == 'scnd' or src == 'cmbn' or src == 'phas'):
  ylab = 'Albedo'
if (src == 'thrm'):
  ylab = r'Specific flux (W/m$^2$/${\rm \mu}$m)'
if (src == 'trns'):
  ylab = r'Transit depth'

# best-fit model
from rfast_retrieve_pp import Fx
plt.scatter(lam, albedo, s = 12, c = 'palevioletred') ##modified for albedo january 11th
plt.errorbar(lam, albedo, yerr = err, fmt = '.k')
plt.plot(lam, Fx(x0,y), color = 'black')
plt.ylabel(ylab)
plt.xlabel(r'Wavelength (' + u'\u03bc' + 'm)')
plt.savefig(dirout+fnr+'_bestfit.png',format='png',bbox_inches='tight')
plt.close()

# residuals
plt.errorbar(lam, albedo-Fx(x0,y), yerr=err, fmt=".k")
plt.ylabel(ylab)
plt.xlabel(r'Wavelength (' + u'\u03bc' + 'm)')
plt.savefig(dirout+fnr+'_residuals.png',format='png',bbox_inches='tight')
plt.close()

# compute & print parameters, truths, mean inferred, and 16/84 percentile (credit: arnaud)
mean = np.zeros(len(pp_names))
std  = np.zeros([2,len(pp_names)])
for i in range(len(pp_names)):
  prcnt    = np.percentile(samples[:,:,i], [16, 50, 84])
  mean[i]  = prcnt[1]
  std[0,i] = np.diff(prcnt)[0]
  std[1,i] = np.diff(prcnt)[1]
colnames = ['Parameter','Input','Mean','- sig','+ sig']
data_out = Table([pp_names,pp_truths,mean,std[0,:],std[1,:]],names=colnames)
ascii.write(data_out,dirout+fnr+'.tab',format='fixed_width',overwrite=True)

#print out the parameters used in this run
print('Apars = ', Apars)
print('nwalkers = ', nwalkers)
print('nstep = ', nstep)







