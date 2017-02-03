# Generates spotted light curve 
# Uses the macula model by Kipping et al. 2011
# Written by Eunkyu Han


import numpy as np
import numpy.random as rand
from pymacula import MaculaModel
from pymacula import macula
import pdb
import george
from george import kernels
import numpy as np
import os
import pdb
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from astropy.io import fits
from astropy.io import ascii
from statsmodels.tsa.stattools import adfuller
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
import mpfit
from astropy.stats import LombScargle
import astropy.units as u


# Read in the McQuillan rotation period file
kid_table, Teff, logg, Mass, Prot, Prot_err, amp, Flag = np.genfromtxt('mcquillan.txt', comments='#', usecols=(0,1,2,4,5,6,7,8), unpack=True)

# Target KIC number
target = 893033				

    
'''
# Specify the data path 
# and read-in the files
path = '/Users/eunkyuh/Dropbox/BU/Classes/12.S680/Project/Data/public_Q6_long_1/'    # Path on ROSALIND
file = os.listdir(path)


# Match the target 
kicid = [] #np.empty(len(file))                 # Q6 data files
kic_name = []
for i in range(len(file)):
    kic_name = str.split(file[i], '-')       # Split the file name into useful portion
    chars = list(kic_name[0])                # Split the kplrXXXXXXXXXX into characters
    kicid.append(int(''.join(chars[4:])))       # Join them back so that they are complete numbers 
'''

kicid = []
path = '/home/eunkyu/Dropbox/BU/Classes/12.S680/Project/codes/Data/'
file = os.listdir(path) #'kplr001026146-2010265121752_llc.fits'
kic_name = str.split(file[0], '-')
chars = list(kic_name[0])
kicid.append(int(''.join(chars[4:])))




       
file_name = file #kic_name[0]+'-'+kic_name[1]
print 'File name is: ', file_name



# A function to read in the data
def read_data(file_name):
	index = np.where(target == kid_table)
	obs = fits.open(file_name)
	tbdata = obs[1].data
	head = obs[1].header
	ut_raw = tbdata['time']
	ut0 = head['BJDREFI']
	flux = tbdata['pdcsap_flux']
	fix = np.where(np.isnan(flux) == 0)             # Removing NaN values
	flux = flux[fix]
	ut_raw = ut_raw[fix]
	flux_err = tbdata['pdcsap_flux_err']
	flux_err = flux_err[fix]
	p_rotation = Prot[index]   

	t = ut0 + ut_raw 				# To feed in the GP

	# Removing the outliers by making 2.5*sigma cut
	std = np.std(flux)
	avg = np.mean(flux)
	out = np.where(flux > avg + 2.5*std)
	flux[out] = np.nan
	fix = np.where(np.isnan(flux) == 0)
	flux = flux[fix]
	ut_raw = ut_raw[fix]
	flux_err = flux_err[fix]
	

	flux_norm = flux/max(flux)
	
	max_index = np.where(max(flux) == flux)
	flux_norm_err = np.sqrt((flux_err/max(flux))**2 + (flux * flux_err[max_index] / max(flux)**2)**2)

	return ut_raw, flux, flux_err, flux_norm, flux_norm_err, p_rotation


'''
# Start reading in the data
# If I want to automate everything....
# Use this portion 

ut = []
flux_norm = []
flux_norm_err = []

for i in range(len(file_name)-1):
	
	q = read_data(path+file_name[i])

	ut = np.append(ut, q[0])
	flux_norm = np.append(flux_norm, q[3])
	flux_norm_err = np.append(flux_norm_err, q[4])
'''


# Read in the Q6 data
ut, flux, flux_err, flux_norm, flux_norm_err, p_rotation = read_data(path+file_name[0])			


# Binning, if needed
# Bin down by 3 points (~ 1.5 hours)
for i in range(np.round(len(ut)/3)):
	ut_binned.append(np.median(ut[i*3:i*3 + 3]))
	flux_norm_binned.append(np.median(flux_norm[i*3:i*3 + 3]))
	flux_norm_err_binned.append(np.median(flux_norm_err[i*3:i*3 + 3]))


delta_ut = ut[1]-ut[0]
ts = np.arange(ut[0], ut[len(ut)-1], 0.02)
#ts = np.arange(539, 719, 0.02)


ndata = 1e4
Nspot = np.array([3])
mmax = 1
derivatives = True
temporal = True
TdeltaV = True



# Defining the model parameters for the star
# Not varying these for fitting
# Except the rotation period will be changed for different star
Istar = np.pi/2			# Inclination of the star [rads]
Peq = p_rotation		# Rotation period of the star
kappa2 = 0.2			# Quadratic differential rotation coeff
kappa4 = 0.2			# Quartic differential rotation coeff
c1 = 0.4				# 1st of four-coeff stellar LD terms
c2 = 0.4269				# 2nd of four-coeff stellar LD terms
c3 = -0.0227 			# 3rd of four-coeff stellar LD terms
c4 = -0.0839			# 4th of four-coeff stellar LD terms
d1 = 0.4 				# 1st of four-coeff spot LD terms
d2 = 0.4269				# 2nd of four-coeff spot LD terms
d3 =-0.0227 			# 3rd of four-coeff spot LD terms
d4 = -0.0839			# 4th of four-coeff spot LD terms

theta_star = [Istar, Peq, kappa2, kappa4, c1, c2, c3, c4, d1, d2, d3, d4]

theta_inst = np.array([1.,1.])




Tstart = ut[0]
Tend = ut[len(ut)-1]

# Using the random distribution of the spots
p0_spot = np.zeros((8,Nspot))
theta_spot = []
for k in range(Nspot):
	p0_spot[0,k] = rand.random()*np.pi 						# longitude
	p0_spot[1,k] = rand.random()*np.pi/2 					# latitude
	p0_spot[2,k] = rand.random()*10*np.pi/180 				# alpha_max (size of spot)
	p0_spot[3,k] = rand.random()*0.5 						# fspot (flux contrast)
	p0_spot[4,k] = rand.random()*(Tend-Tstart) + Tstart 	# Time at which spot k is largest	
	p0_spot[5,k] = rand.random()*(Tend-Tstart)				# Lifetime of spot k (FWFM) [days]
	p0_spot[6,k] = rand.random()*(Tend-Tstart)				# Ingress duration of spot k [days]
	p0_spot[7,k] = rand.random()*(Tend-Tstart)				# Egress duration of spot k [days]
	
	index = np.argwhere(flux_norm_binned ==np.min(flux_norm_binned))

	for i in range(8):
		theta_spot.append(p0_spot[i, k])



# Set up to run the mpfit
# Fitting for the spots
parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'parname':'name'} for i in range(Nspot*8)]

################ 
# Spot 1
parinfo[0]['value'] = theta_spot[0]
parinfo[0]['fixed'] = 0
parinfo[0]['limited'][0] = 1
parinfo[0]['limited'][1] = 1
parinfo[0]['limits'][0] = 0
parinfo[0]['limits'][1] = 2*np.pi
parinfo[0]['parname'] = 'longitude of spot 1'

parinfo[1]['value'] = theta_spot[1]
parinfo[1]['fixed'] = 0
parinfo[1]['limited'][0] = 1
parinfo[1]['limited'][1] = 1
parinfo[1]['limits'][0] = -1.* np.pi
parinfo[1]['limits'][1] = np.pi
parinfo[1]['parname'] = 'latitude of spot 1'

parinfo[2]['value'] = theta_spot[2]
parinfo[2]['fixed'] = 0
parinfo[2]['limited'][0] = 1
parinfo[2]['limits'][0] = 0
parinfo[2]['parname'] = 'angluar size of spot 1'

parinfo[3]['value'] = theta_spot[3]
parinfo[3]['fixed'] = 0
parinfo[3]['limited'][0] = 1
parinfo[3]['limits'][0] = 0.0
parinfo[3]['parname'] = 'flux contrast'

parinfo[4]['value'] = theta_spot[4]
parinfo[4]['fixed'] = 0
parinfo[4]['limited'][0] = 1
parinfo[4]['limited'][1] = 1
parinfo[4]['limits'][0] = np.min(ut)
parinfo[4]['limits'][1] = np.max(ut)
parinfo[4]['parname'] = 'time at which spot is largest'

parinfo[5]['value'] = theta_spot[5]
parinfo[5]['fixed'] = 0
parinfo[5]['limited'][0] = 1
parinfo[5]['limits'][0] = 0
parinfo[5]['parname'] = 'life time of spot 1'

parinfo[6]['value'] = theta_spot[6]
parinfo[6]['fixed'] = 0
parinfo[6]['limited'][0] = 1
parinfo[6]['limits'][0] = 0
parinfo[6]['parname'] = 'ingress duration of spot 1'

parinfo[7]['value'] = theta_spot[7]
parinfo[7]['fixed'] = 0
parinfo[7]['limited'][0] = 1
parinfo[7]['limits'][0] = 0
parinfo[7]['parname'] = 'egress duration of spot 1'

################ 
# Spot 2
parinfo[8]['value'] = theta_spot[8]
parinfo[8]['fixed'] = 0
parinfo[8]['limited'][0] = 1
parinfo[8]['limited'][1] = 1
parinfo[8]['limits'][0] = 0
parinfo[8]['limits'][1] = 2*np.pi
parinfo[8]['parname'] = 'longitude of spot 2'

parinfo[9]['value'] = theta_spot[9]
parinfo[9]['fixed'] = 0
parinfo[9]['limited'][0] = 1
parinfo[9]['limited'][1] = 1
parinfo[9]['limits'][0] = -1.* np.pi
parinfo[9]['limits'][1] = np.pi
parinfo[9]['parname'] = 'latitude of spot 2'

parinfo[10]['value'] = theta_spot[10]
parinfo[10]['fixed'] = 0
parinfo[10]['limited'][0] = 1
parinfo[10]['limits'][0] = 0
parinfo[10]['parname'] = 'size of spot 2'

parinfo[11]['value'] = theta_spot[11]
parinfo[11]['fixed'] = 0
parinfo[11]['limited'][0] = 1
parinfo[11]['limits'][0] = 0
parinfo[11]['parname'] = 'flux contrast'

parinfo[12]['value'] = theta_spot[12]
parinfo[12]['fixed'] = 0
parinfo[12]['limited'][0] = 1
parinfo[12]['limited'][1] = 1
parinfo[12]['limits'][0] = np.min(ut)
parinfo[12]['limits'][1] = np.max(ut)
parinfo[12]['parname'] = 'time at which spot is largest'

parinfo[13]['value'] = theta_spot[13]
parinfo[13]['fixed'] = 0
parinfo[13]['limited'][0] = 1
parinfo[13]['limits'][0] = 0
parinfo[13]['parname'] = 'life time of spot 2'

parinfo[14]['value'] = theta_spot[14]
parinfo[14]['fixed'] = 0
parinfo[14]['limited'][0] = 1
parinfo[14]['limits'][0] = 0
parinfo[14]['parname'] = 'ingress duration of spot 2'

parinfo[15]['value'] = theta_spot[15]
parinfo[15]['fixed'] = 0
parinfo[15]['limited'][0] = 1
parinfo[15]['limits'][0] = 0
parinfo[15]['parname'] = 'egress duration of spot 2'

################ 
# Spot 3
parinfo[16]['value'] = theta_spot[16]
parinfo[16]['fixed'] = 0
parinfo[16]['limited'][0] = 1
parinfo[16]['limited'][1] = 1
parinfo[16]['limits'][0] = 0
parinfo[16]['limits'][1] = 2*np.pi
parinfo[16]['parname'] = 'longitude of spot 3'

parinfo[17]['value'] = theta_spot[17]
parinfo[17]['fixed'] = 0
parinfo[17]['limited'][0] = 1
parinfo[17]['limited'][1] = 1
parinfo[17]['limits'][0] = -1.* np.pi
parinfo[17]['limits'][1] = np.pi
parinfo[17]['parname'] = 'latitude of spot 3'

parinfo[18]['value'] = theta_spot[18]
parinfo[18]['fixed'] = 0
parinfo[18]['limited'][0] = 1
parinfo[18]['limits'][0] = 0
parinfo[18]['parname'] = 'size of spot 3'

parinfo[19]['value'] = theta_spot[19]
parinfo[19]['fixed'] = 0
parinfo[19]['limited'][0] = 1
parinfo[19]['limits'][0] = 0
parinfo[19]['parname'] = 'flux contrast'

parinfo[20]['value'] = theta_spot[20]
parinfo[20]['fixed'] = 0
parinfo[20]['limited'][0] = 1
parinfo[20]['limited'][1] = 1
parinfo[20]['limits'][0] = np.min(ut)
parinfo[20]['limits'][1] = np.max(ut)
parinfo[20]['parname'] = 'time at which spot is largest'

parinfo[21]['value'] = theta_spot[21]
parinfo[21]['fixed'] = 0
parinfo[21]['limited'][0] = 1
parinfo[21]['limits'][0] = 0
parinfo[21]['parname'] = 'life time of spot 3'

parinfo[22]['value'] = theta_spot[22]
parinfo[22]['fixed'] = 0
parinfo[22]['limited'][0] = 1
parinfo[22]['limits'][0] = 0
parinfo[22]['parname'] = 'ingress duration of spot 3'

parinfo[23]['value'] = theta_spot[23]
parinfo[23]['fixed'] = 0
parinfo[23]['limited'][0] = 1
parinfo[23]['limits'][0] = 0
parinfo[23]['parname'] = 'egress duration of spot 3'


def make_model(theta_spot, x, y, err):
	model_t = np.linspace(x[0]-0.001, x[len(x)-1]+0.001, 10000)
	model_t = np.linspace(x[0]-0.1, x[len(x)-1]+0.1, 10000)		# For the binned data

	# Re-formatting the parameters to feed into the model
	theta_spots = np.ndarray([8, Nspot])
	for i in range(Nspot):
		for j in range(8):
			theta_spots[j][i] = theta_spot[8*i+j]

	'''
	# Using the random distribution of the spots
	theta_spot = np.zeros((8,Nspot))
	for k in range(Nspot):
		theta_spot[0,k] = rand.random()*np.pi 						# longitude
		theta_spot[1,k] = rand.random()*np.pi/2 					# latitude
		theta_spot[2,k] = rand.random()*10*np.pi/180 				# alpha_max (size of spot)
		theta_spot[3,k] = rand.random()*0.5 						# fspot (flux contrast)
		theta_spot[4,k] = rand.random()*(Tend-Tstart) + Tstart 		# tmax
		theta_spot[5,k] = rand.random()*(Tend-Tstart)				# Time at which spot k is largest
		theta_spot[6,k] = rand.random()*(Tend-Tstart)				# Lifetime of spot k (FWFM) [days]
		theta_spot[7,k] = rand.random()*(Tend-Tstart)				# Ingress duration of spot k [days]
		#theta_spot[8,k]	= rand.random()*(Tend-Tstart)				# Egress duration of spot k  [days]
	'''


    # Models are made with the parameters
	model_y = macula(model_t, theta_star, theta_spots, theta_inst, derivatives, temporal)	
	return model_y, model_t


def myfunct(theta_spot, fjac=None, x=None, y=None, err=None):

	model_y, model_t = make_model(theta_spot, x, y, err)

	int_time = 1626.0/(3600.0*24.0) / p_rotation #period # sec
	deltaT_in_model = model_t[1] - model_t[0]
	smooth=scipy.ndimage.filters.uniform_filter(model_y, size = int_time/deltaT_in_model, mode = 'wrap')

	
	interp = interpolate.interp1d(model_t, smooth) 
	model_sm_t = x
	model_sm_y = interp(model_sm_t)
	

	sort = np.argsort(x)
	model_sm_y = model_sm_y[sort]
	#x = x[sort]
	#y = y[sort]
	#err = err[sort]
	
	status = 0
    
	plt.figure(1)
	plt.cla()
	plt.plot(x, y, 'bo', markersize = 4, label = 'Kepler Data')
	plt.plot(x, model_sm_y, 'r', label = 'Macula Model')
	plt.xlabel('BJD - 2454833', fontsize = 17)
	plt.ylabel('Relative Flux', fontsize = 17)
	plt.title(str(target))
	plt.legend(loc=0)


	return ([status, ((y - model_sm_y)/ err)])


# Structure that contains observational data
fa = {'x':ut, 'y':flux_norm, 'err':flux_norm_err}
#fa = {'x':ut_binned, 'y':flux_norm_binned, 'err':flux_norm_err_binned}

print '=============================================='
print ''
print 'Running mpfit'
print ''
print '=============================================='
initial_fit_parameters = mpfit.mpfit(myfunct, theta_spot, functkw=fa, parinfo=parinfo)


print "status: ", initial_fit_parameters.status
if (initial_fit_parameters.status <= 0): 
   print 'error message = ', initial_fit_parameters.errmsg
else:
   print "Iterations: ", initial_fit_parameters.niter
   print "Fitted pars: ", np.double(initial_fit_parameters.params)
   print "Uncertainties: ", np.double(initial_fit_parameters.perror)


# Check the trained model
model_y, model_t = make_model(initial_fit_parameters.params, ut, flux_norm, flux_norm_err)
model_yerr = np.std(model_y)


# Checking for the validity of the model
# For the Q6 data
frequency, power = LombScargle(ut, flux_norm, flux_norm_err).autopower()
model_f, model_power = LombScargle(model_t, model_y, model_yerr).autopower()


plt.figure(3)
plt.cla()
plt.plot(1/frequency, power, label = 'Kepler Data') 
plt.plot(1/model_f, model_power, color = 'r', label = 'Macula Model')
plt.axvline(x = p_rotation, color = 'k', linestyle = '--', label = 'Rotation Period')
plt.legend(loc = 0)
plt.xlim([0., 130.])
plt.title(str(target))
plt.xlabel('Period [day]', fontsize=17)
plt.ylabel('Power', fontsize=17)
plt.savefig(str(target)+'_Q6_kp_vs_pm_power.png')

plt.figure(1)
plt.savefig(str(target)+'_Q6_fit.png')

pdb.set_trace()
#########################
# Next quarter data
#########################

# Read-in the Q7 data
t, f, f_err, f_norm, f_norm_err, prot = read_data(path+file_name[1])
# Make a model to match the time
final_model_y, final_model_t = make_model(initial_fit_parameters.params, t, f_norm, f_norm_err)

# Check the validity of the model
frequency2, power2 = LombScargle(t, f_norm, f_norm_err).autopower()
final_model_f, final_model_power = LombScargle(final_model_t, final_model_y, model_yerr).autopower()


# Plot to see the result 
plt.figure(5)
plt.cla()
plt.plot(1/frequency2, power2, 'b', label= 'Kepler Data')
plt.plot(1/final_model_f, final_model_power, 'r', label = 'Macula Model')
plt.axvline(x = p_rotation, color = 'k', linestyle = '--', label = 'Rotation Period')
plt.legend(loc=0)
plt.xlim([0., 130.])
plt.title(str(target))
plt.xlabel('Period [day]', fontsize=17)
plt.ylabel('Power', fontsize=17)
plt.savefig(str(target)+'_Q7_kp_vs_pm_power.png')


plt.figure(2)
plt.cla()
plt.plot(t, f_norm, 'o', label = 'Kepler Data')
plt.plot(final_model_t, final_model_y*(np.median(f_norm)/np.median(final_model_y)), 'r', label = 'Macula Model')
plt.xlabel('BJD - 2454833', fontsize = 17)
plt.ylabel('Relative Flux', fontsize = 17)
plt.title(str(target))
plt.legend(loc=0)
plt.savefig(str(target)+'_Q7_fit.png')

