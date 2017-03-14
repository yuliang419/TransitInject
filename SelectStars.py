import numpy as np
import pandas as pd
import glob

varcat = pd.read_table('varcat.txt', comment='#', header=0, skip_blank_lines=True)
all_bad = varcat[(varcat.Class != 'Noise ') & (varcat.Class != 'OTHPER')].ID.values
print 'No. of variables & EBs:', len(all_bad)

camps = glob.glob('exofop*csv')
for camp in camps:
	all_stars = pd.read_csv(camp, comment='\\', header=0)
	planets = all_stars[all_stars['Poss_pc'] == 'Y']
	planets = planets['EPIC ID'].values
	print 'No. of planet hosts:', len(planets)
	all_bad = np.append(all_bad, planets)

targlist = glob.glob('k2/K2Campaign*targets.csv')
mdwarf_props = 'GO3001|GO3011|GO301|GO3031|GO3034|GO3045|GO3069|GO3073|GO3094|GO3106|GO3107|GO3114'
mdwarfs = np.array([])
for camp in targlist:
	data = pd.read_csv(camp, header=0, index_col=False, usecols=[0,4], \
		na_values='null')

	md = data[data['Investigation IDs'].str.contains(mdwarf_props)]
	good_md = md[~md['EPIC ID'].isin(all_bad)]
	good_md = good_md['EPIC ID'].values.astype('str')

	mdwarfs = np.append(mdwarfs, good_md)

np.savetxt('k2/k2mdwarfs/mdwarfs.ls', mdwarfs, fmt='%10s')


