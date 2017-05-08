import urllib
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

# pd.compat.PY3 = True

# data = pd.read_csv('K2Campaign3targets.csv', header=0, index_col=False, usecols=[0,4], \
# 	na_values='null')
# mdwarf_props = 'GO1002|GO1006|GO1036|GO1037|GO1047|GO1050|GO1051|GO1053|GO1062|GO1075|GO2106|GO2073|GO2070|GO2069|GO2045|GO2031|GO2013|GO2011|\
# GO2011|GO3001|GO3011|GO301|GO3031|GO3034|GO3045|GO3069|GO3073|GO3094|GO3106|GO3107|GO3114'
# mdwarfs = data[data['Investigation IDs'].str.contains(mdwarf_props)]
campaign = '03'
mdwarfs = np.loadtxt('k2mdwarfs/mdwarfs.ls', dtype=str)

def download(epic):
	epic = str(epic)

	print 'Downloading M dwarf '+epic
	getfile = urllib.URLopener()
	# download campaign data
	try:
		getfile.retrieve('http://archive.stsci.edu/missions/hlsp/k2sff/c'+campaign+'/'+epic[:4]+'00000/'+epic[-5:]+\
			'/hlsp_k2sff_k2_lightcurve_'+epic+'-c'+campaign+'_kepler_v1_llc-default-aper.txt', 'k2mdwarfs/'+epic+'.txt')
	except IOError:
		print 'Download failed'
		pass

Parallel(n_jobs=2)(delayed(download)(epic) for epic in mdwarfs)
