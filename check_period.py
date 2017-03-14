import numpy as np
import glob
from joblib import Parallel, delayed

with open('success.txt', 'w') as outfile:
	print>>outfile, '# epic 	p 	p_true 	Rp 	bls		SNR 	search type'
with open('fail.txt', 'w') as outfile:
	print>>outfile, '# epic p Rp b'

epic_list = np.loadtxt('epic_list.txt')

def check(epic):
	files = glob.glob('bls/'+epic+'*blsanal')
	if len(files) == 0:
		return

	for f in files:
		fname = f.split('ltf.')
		pars = f.split('_')
		epic = pars[0]
		p_true = float(pars[2])
		rp_true = float(pars[3])
		b_true = float(pars[4])
		try:
			p, tc, snr = np.loadtxt(f, usecols=(1,2,20), unpack=True)
			p2, tc2, snr2 = np.loadtxt('bls_short/'+fname+'ltf.blsanal', usecols=(1,2,20), unpack=True)
		except:
			print 'Error reading BLS file'
			with open('blserr_log.txt', 'a') as errlog:
				print>>errlog, f
				continue

		for i in range(len(p)):
			if snr[i] > snr2[i]:
				mode = 'long'
				sn_opt = snr[i]
				p_opt = p[i]
			else:
				mode = 'short'
				sn_opt = snr2[i]
				p_opt = p2[i]

			if (abs(p_opt - p_true) < 0.1) and sn_opt > 9:
				with open('success.txt', 'a') as outfile:
					print>>outfile, epic, p_opt, p_true, rp_true, b_true, sn_opt, mode
			else:
				with open('fail.txt', 'a') as outfile:
					print>>outfile, epic, p_true, rp_true, b_true

Parallel(n_jobs=4)(delayed(check)(epic) for epic in epic_list)
