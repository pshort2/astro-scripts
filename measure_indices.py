import numpy as np
import matplotlib.pyplot as plt

from lmfit.models import GaussianModel

import targets as t
import index_functions as f

path = '../index_plots_results/' #output path


#Import spectrum of desired object
objects = ['14ae','14li','15oi','16fnl','ahk','azh','dsg','dyb', 'fyk','hyz','qiz'] 


def load_data(object):

	spec = np.genfromtxt('../spec/'+object+'_s_bin.txt',comments='#')

	return spec

##########################
############Data###########
##########################
d4, d4_err = [], []
hd, hd_err = [], []
ha, ha_err = [], []
for target in objects:

	#Load object dictionary (from targets.py)
	o = t.load_obj(target)
	obj = o['obj']
	z = o['z'] #redshift
	
	#Load spectrum
	#Infilling removed for ahk
	if obj == 'ahk':
		s = np.genfromtxt('../spec/ahk_s_bin_infill_removed.txt')
	else:
		s = load_data(obj)

	if obj == '14ae' or obj == 'qiz':
		s_nii = s
	elif  obj == 'hyz':
		s_nii = np.genfromtxt('../spec/nii_removed/'+obj+'_no_ha.txt') #tde h-alpha emission removed
		s_nii[:,0] /= (1+z)
	else:
		s_nii = np.genfromtxt('../spec/nii_removed/'+obj+'_no_nii.txt')
		s_nii[:,0] /= (1+z)

	#Correct for redshift
	s[:,0] /= (1+z)

	print (obj)

	#Use index functions to calculate indices
	hdelta = 	f.hdelta(s,target,True)
	d4000 = f.d4000(s,target,False)
	halpha = f.halpha(s_nii,target,False)

	d4.append(d4000[0])
	hd.append(hdelta[0])
	ha.append(halpha[0])

	d4_err.append(d4000[1])
	hd_err.append(hdelta[1])
	ha_err.append(halpha[1])

	#print (d4000[0])
	print (round(hdelta[0],5))


d4_out = np.zeros((len(d4),2))
d4_out[:,0] = d4
d4_out[:,1] = d4_err

hd_out = np.zeros((len(hd),2))
hd_out[:,0] = hd
hd_out[:,1] = hd_err

ha_out = np.zeros((len(ha),2))
ha_out[:,0] = ha
ha_out[:,1] = ha_err

np.savetxt(path+'d4000_20220418.txt', d4_out)
np.savetxt(path+'hd_20220418.txt', hd_out)
np.savetxt(path+'ha_20220418.txt', ha_out)


'''
###########################
###########Models###########
###########################
ages = np.arange(0.5, 13, 0.5)
tz = [(0.5,0.2), (0.5,1), (0.5,2.5), (1,0.2), (2,1)]

for x in tz:

	tau = x[0]
	Z = x[1]

	#def EW(ageects, redshift, phot, 
	i=0
	d4_mod = []
	hd_mod= []
	ha_mod=[]
	for age in ages:

		#Import spectrum
		s = np.genfromtxt('../../consistent_sfh_models/age_'+str(age)+'tau'+str(tau)+'_Z'+str(Z))
			
		#Correct for redshift 
		z = 0.015
		s[:,0] /= (1+z)

		hdd = f.hdelta(s,'s',False)
		print (round(hdd,5))

		hd_mod.append(hdd)
		ha_mod.append(f.halpha(s,'s',False))
		d4_mod.append(f.d4000(s,'s',False))
		
	np.savetxt(path+'D4000_tau'+str(tau)+'_Z'+str(Z)+'.txt',d4_mod)
	np.savetxt(path+'Hd_tau'+str(tau)+'_Z'+str(Z)+'.txt',hd_mod)
	np.savetxt(path+'Ha_tau'+str(tau)+'_Z'+str(Z)+'.txt',ha_mod)

###############################
###########Fitted spectra###########
###############################
def load_fits(run):

	spec = np.genfromtxt('../bp_spectra_fits/' + run + '.txt',comments='#')

	return spec

d4 = []
hd = []
#ha = []
for target in objects:

	#Load object dictionary (from targets.py)
	o = t.load_obj(target)
	obj = o['obj']
	z = o['z'] #redshift
	
	#Load spectrum
	end = "_uvb_psb_bin_20sigphoterr_lines20211220"
	run = obj+ end
	s = load_fits(run)

	#Correct for redshift
	s[:,0] /= (1+z)

	print (obj)

	#Use index functions to calculate indices
	hdelta = 	f.hdelta(s,target,False)
	d4000 = f.d4000(s,target,False)
	#halpha = f.halpha(s,target,False)

	d4.append(d4000)
	hd.append(hdelta)
	#ha.append(halpha)

	print (round(hdelta,5))

np.savetxt(path+'fit_d4000_20211220.txt', d4)
np.savetxt(path+'fit_hd_20211220.txt', hd)
#np.savetxt(path+'fit_ha_20211220.txt', ha)
'''

