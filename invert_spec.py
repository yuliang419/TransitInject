import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


def ifft(t, nu, re, im, N, m):
    ts = np.zeros(len(t))
    mean_t = np.mean(t)
    for s in range(len(nu)):
        term = N/m * (re[s] * np.cos(2 * np.pi * nu[s] * (-t+mean_t)) - im[s] * np.sin(2 * np.pi * nu[s] * (-t+mean_t)))
        ts += term

    return ts


def get_spec_peaks(name, path='dft/'):
    """
    Use only most prominent peaks in DFT spectrum to reconstruct light curve.
    :param name: name of light curve.
    :param path: path containing DFT files.
    :return:
    """
    dspec = pd.read_csv(path + name + '.dftclean.dspec', sep=' ', dtype=np.float64, skiprows=1, names=['freq', 'pow', 're', 'im'])
    m = float(len(dspec))

    dftanal = pd.read_csv(path+name.split('.lc')[0]+'.dftanal', sep=' ', dtype=np.float64, index_col=0, skiprows=1, names=['freq', 'pow', 'snr'])
    nu = dftanal['freq'].values

    re = []
    im = []
    power = []
    for nu_s in nu:
        match = dspec[abs(dspec['freq']- nu_s) <1e-3].index[0]
        re.append(dspec.iloc[match]['re'])
        im.append(dspec.iloc[match]['im'])
        power.append(dspec.iloc[match]['pow'])

    return nu, power, re, im, m


def get_spec(name, path='dft/', fmax=None):
    """
    Use all frequencies < 0.5 to reconstruct light curve.
    :param name: name of light curve.
    :param path: path containing DFT files.
    :param fmax: max frequency to use in reconstruction (1/day). Set to None to automatically determine max freq.
    :return:
    """
    dspec = pd.read_csv(path + name + '.dftclean.dspec', sep=' ', dtype=np.float64, skiprows=1,
                        names=['freq', 'pow', 're', 'im'])
    m = float(len(dspec))
    all_pow = dspec['pow'].values

    if fmax is None:
        ind = np.argmax(all_pow)
        fmax = dspec.iloc[ind, 0]
        if fmax <= 0:
            fmax = -fmax
        if fmax < 0.5:
            fmax = 0.5
        fmax += 0.5

    print 'Max frequency:', fmax

    sel = dspec[abs(dspec['freq']) < fmax]
    nu = sel['freq'].values

    print "using %d peaks" % len(nu)
    re = []
    im = []
    power = []
    for nu_s in nu:
        match = dspec[abs(dspec['freq'] - nu_s) < 1e-3].index[0]
        re.append(dspec.iloc[match]['re'])
        im.append(dspec.iloc[match]['im'])
        power.append(dspec.iloc[match]['pow'])

    return nu, power, re, im, m


def main(lc, dftpath='box/dft/', lcpath='box/mags-30/'):
    nu, power, re, im, m = get_spec(lc, path=dftpath)
    t, mag = np.loadtxt(lcpath + lc, unpack=True, usecols=(0, 2))
    N = float(len(t))
    ts = ifft(t, nu, re, im, N, m)

    mag -= np.mean(mag)
    sig = np.std(mag)
    sigmod = np.std(ts)
    ts *= sig / sigmod
    return t, mag, ts


if __name__ == '__main__':
    lc = sys.argv[1]
    t, mag, ts = main(lc)
    plt.plot(t, mag, 'b.')
    plt.plot(t, ts, 'r')
    plt.show()