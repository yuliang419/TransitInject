import urllib
import numpy as np

kic_table, Teff, logg, Mass, Prot, Prot_err, amp, Flag = np.genfromtxt('mcquillan.txt', comments='#', usecols=(0,1,2,4,5,6,7,8), unpack=True)

for kic in kic_table:
	kic = str(kic)
	while len(kic) < 8:
		kic = '0'+kic
	print 'Downloading other '+kic
	getfile = urllib.URLopener()
	# download Q1 data (1 month)
	try:
		getfile.retrieve('https://archive.stsci.edu/pub/kepler/lightcurves/0'+kic[:3]+'/0'+kic+'/kplr0'+kic+'-2009259160929_llc.fits', 'mdwarfs/'+kic+'.fits')
	except IOError:
		continue