# import statements
import time
import shutil
import sys
import numpy             as np
import matplotlib.pyplot as plt
from astropy.table       import Table, Column, MaskedColumn
from astropy.io          import ascii
from rfast_routines      import noise
from rfast_routines      import inputs

# get input script filename
if len(sys.argv) >= 2:
  filename_scr = sys.argv[1] # if script name provided at command line
else:
  filename_scr = input("rfast inputs script filename: ") # otherwise ask for filename

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
fn_dat = fns + '.raw'

# read input data
data        = ascii.read(dirout+fn_dat,data_start=1,delimiter='|')
lam         = data['col2'][:]
dlam        = data['col3'][:]
F1          = data['col4'][:] #albedo
F2          = data['col5'][:] #flux ratio

# snr0 constant w/wavelength case
if( len(snr0) == 1 ): #CHANGED F2 TO F1 HERE SO ITS CALCULATING ERROR OFF ALBEDO
  if (ntype != 'cppm'):
    err = noise(lam0,snr0,lam,dlam,F1,Ts,ntype)
  else:
    err    = np.zeros(F1.shape[0])
    err[:] = 1/snr0
else: # otherwise snr0 is bandpass dependent
  err = np.zeros(len(lam))
  for i in range(0,len(snr0)):
    ilam = np.where(np.logical_and(lam >= lams[i], lam <= laml[i]))
    if (len(lam0) == 1): # lam0 may be bandpass dependent
      lam0i = lam0
    else:
      lam0i = lam0[i]
    if (ntype != 'cppm'):
      erri      = noise(lam0i,snr0[i],lam,dlam,F1,Ts,ntype)
      err[ilam] = erri[ilam]
    else:
      err[ilam] = 1/snr0[i]

# generate faux spectrum, with random noise if requested
data = np.copy(F1)
if rnd:
  for k in range(0,len(lam)):
    data[k]  = np.random.normal(F1[k], err[k], 1)
    if data[k] < 0:
      data[k] = 0.

# write data file
if (src == 'diff' or src == 'cmbn'):
  names = ['wavelength','d_wavelength','albedo','flux_ratio','data','uncertainty']
if (src == 'thrm'):
  names = ['wavelength (um)','d wavelength (um)','Tb (K)','flux (W/m**2/um)','data','uncertainty']
if (src == 'scnd'):
  names = ['wavelength (um)','d wavelength (um)','Tb (K)','flux ratio','data','uncertainty']
if (src == 'trns'):
  names = ['wavelength (um)','d wavelength (um)','zeff (m)','transit depth','data','uncertainty']
if (src == 'phas'):
  names = ['wavelength (um)','d wavelength (um)','reflect','flux ratio','data','uncertainty']
data_out = Table([lam,dlam,F1,F2,data,err], names=names)
ascii.write(data_out,dirout+fnn+'.dat',format='fixed_width',overwrite=True)

from rfast_genspec_alb import name1, name2

#output as .csv for spectra_gen
ascii.write(data_out,'data_error.csv',overwrite=True)


# document parameters to file
shutil.copy(filename_scr,dirout+fnn+'.log')

# plot faux data
if (src == 'diff' or src == 'scnd' or src == 'cmbn' or src == 'phas'):
  ylab = 'Albedo'
if (src == 'thrm'):
  ylab = r'Specific flux (W/m$^2$/${\rm \mu}$m)'
if (src == 'trns'):
  ylab = r'Transit depth'
plt.errorbar(lam, data, yerr=err, fmt=".k") #data should now be albedo
plt.ylabel(ylab)
plt.grid(alpha = 0.5)
plt.xlabel(r'Wavelength (' + u'\u03bc' + 'm)')
plt.savefig('forward_model_noise.png',format='png',bbox_inches='tight')
plt.close()
