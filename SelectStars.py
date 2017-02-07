import numpy as np
import pandas as pd
import glob

varcat = pd.read_table('varcat.txt', comment='#', header=0)
variables = varcat[(varcat.Class != 'Noise ') & (varcat.Class != 'OTHPER')].values

all_stars = pd.read_csv('exofop_campaign1.csv', comment='\\', header=1)
planets = all_stars[all_stars['Poss_pc'] == 'Y'].values

targlist = glob.glob('K2Campaign*targets.csv')
mdwarf_props = 'GO1002|GO1006|GO1036|GO1037|GO1047|GO1050|GO1051|GO1053|GO1062|GO1075|GO2106|GO2073|GO2070|GO2069|GO2045|GO2031|GO2013|GO2011|\
GO2011|GO3001|GO3011|GO301|GO3031|GO3034|GO3045|GO3069|GO3073|GO3094|GO3106|GO3107|GO3114'
mdwarfs = np.array([])
for camp in targlist:
	data = pd.read_csv(camp, header=0, index_col=False, usecols=[0,4], \
		na_values='null')

	md = data[data['Investigation IDs'].str.contains(mdwarf_props)].values
	mdwarfs = np.append(mdwarfs, md)

