import numpy as np
from joblib import Parallel, delayed

with open('success.txt', 'w') as outfile:
    print>> outfile, '# epic 	p 	p_true 	Rp 	bls		SNR 	search type'
with open('fail.txt', 'w') as outfile:
    print>> outfile, '# epic p Rp b'

epic_list = np.loadtxt('test_epic_list.txt')


class InjectedLC:
    def __init__(self, name, indir='bls/', indir_short='bls_short/'):
        self.name = name
        self.indir = indir
        self.indir_short = indir_short
        fname = name.split('ltf.')[0]
        pars = fname.split('_')
        epic = pars[0]
        self.p_true = float(pars[2])
        self.rp_true = float(pars[3])
        self.b_true = float(pars[4])
        self.epic = epic
        self.p_all = []
        self.snr_all = []

    def __repr__(self):
        return '%s, %s, %s, %s' % (self.epic, self.p_true, self.rp_true, self.b_true)

    def get_peaks(self):
        fname = self.name.split('ltf.')[0]
        try:
            p, tc, snr = np.loadtxt(self.indir + self.name, usecols=(1, 2, 20), unpack=True)
            p2, tc2, snr2 = np.loadtxt(self.indir_short + fname + 'ltf.blsanal', usecols=(1, 2, 20), unpack=True)
        except:
            print 'Error reading BLS file'
            with open('blserr_log.txt', 'a') as errlog:
                print>> errlog, self.name
            return

        p_all = np.concatenate((p, p2))
        snr_all = np.concatenate((snr, snr2))
        sig = np.where(snr_all >= 11)
        p_all = p_all[sig]
        snr_all = snr_all[sig]

        inds = sorted(range(len(snr_all)), key=lambda k: snr_all[k])
        self.p_all = p_all[inds]
        self.snr_all = snr_all[inds]

    def check_period(self):
        correct = 0
        for i in range(len(self.p_all)):
            if abs(self.p_all[i] - self.p_true) < 0.1:
                correct = 1
                with open('success.txt', 'a') as outfile:
                    print>> outfile, self, self.p_all[i], self.snr_all[i]
            else:
                with open('fp.txt', 'a') as outfile:
                    print>> outfile, self, self.p_all[i], self.snr_all[i]
        if correct == 0:
            with open('failed.txt', 'a') as outfile:
                print>> outfile, self


def main(name):
    print name
    lc = InjectedLC(name)
    lc.get_peaks()
    lc.check_period()


if __name__ == '__main__':
    inlist = np.loadtxt('box/flattened-1.5/mdwarfs.ls')
    Parallel(n_jobs=6)(delayed(main)(name) for name in inlist)
