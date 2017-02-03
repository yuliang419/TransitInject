import numpy as np
import matplotlib.pyplot as plt
import glob


files = glob.glob('k2mdwarfs/*ltf.lc')

for f in files:
	name = f.split('/')[1]
	name = name.split('_')[0]
	print name
	t, f, newf = np.genfromtxt(f, unpack=True, usecols=(0,1,2))
	plt.close('all')
	fig = plt.figure(figsize=(10,3))
	plt.annotate(name, xy=(0.8, 0.9), xycoords='axes fraction')
	plt.plot(t,f, 'b.')
	plt.plot(t,newf,'r.')
	plt.savefig('k2mdwarfs/figs/'+name+'.png', dpi=150, bbox_inches='tight')
