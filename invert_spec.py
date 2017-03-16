import numpy as np
import matplotlib.pyplot as plt


def ifft(t, nu, re, im, N, m):
    ts = np.zeros(len(t))
    for s in range(len(nu)):
        term = N*2 * (re[s] * np.cos(2 * np.pi * nu[s] * t) + im[s] * np.sin(2 * np.pi * nu[s] * t))
        ts += term

    return ts


def get_spec(name, path='dft/'):
    nu_all, power, re, im = np.loadtxt(path + name + '.dftclean.dspec', unpack=True)
    m = float(len(nu_all))
    # nu = np.loadtxt()
    nu = nu_all

    return nu, power, re, im, m


if __name__ == '__main__':
    lc = '206005460_2178.63_36.42_2.24_0.49_ltf.lc'
    nu, power, re, im, m = get_spec(lc)
    t, mag = np.loadtxt('mags/' + lc, unpack=True, usecols=(0, 1))
    N = float(len(t))
    ts = ifft(t, nu, re, im, N, m)

    plt.plot(t, mag, 'b.')
    plt.plot(t, ts, 'r')
    plt.show()
