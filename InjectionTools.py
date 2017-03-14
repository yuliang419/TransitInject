import numpy as np
import urllib
from lxml import etree
import batman
import ldtk
from scipy.interpolate import interp1d
from joblib import Parallel, delayed
import random
import glob
import matplotlib.pyplot as plt
import os


def planet_gen(target, outdir, model='mandelagol', return_lc=False):
    G = 2957.4493  # R_sun^3/M_sun/day^2
    t0 = random.uniform(min(target.time), max(target.time))
    P = 10.**random.uniform(np.log10(0.5), np.log10(80))
    rad = random.uniform(0.5, 4)

    # oversample
    cad = 29.4 / 60. / 24.  # cadence in days
    n_ev = 25
    dt_new = np.array([])

    for i, this_t in enumerate(target.time):
        temp = this_t + 1.0 / n_ev * cad * (np.arange(n_ev) - np.ceil(n_ev / 2))
        dt_new = np.append(dt_new, temp)

    if model == 'mandelagol':
        b = random.uniform(0.01, 1)

        params = batman.TransitParams()
        params.t0 = t0
        params.per = P
        params.rp = rad / 109. / target.radius
        params.a = (P ** 2 * G * target.mass / (4 * np.pi ** 2)) ** (1. / 3.) / target.radius  # in units of stellar rad
        params.inc = (180. / np.pi) * np.arccos(1. / params.a * b)
        params.ecc = 0
        params.w = 90
        params.limb_dark = 'quadratic'
        params.u = [target.u1, target.u2]
        m = batman.TransitModel(params, dt_new)
        model_lc = m.light_curve(params)

    elif model == 'box':
        a_Rs = (P ** 2 * G * target.mass / (4 * np.pi ** 2)) ** (1. / 3.) / target.radius
        depth = (rad / 109. / target.radius)**2
        dur = P / np.pi / a_Rs  # characteristic transit duration
        start = -np.floor((t0-target.time[0])/P)
        end = np.floor((target.time[-1]-t0)/P)
        midpts = t0 + np.arange(start, end)*P

        model_lc = np.ones(len(dt_new))

        for pt in midpts:
            now = np.where((dt_new > pt-dur/2) & (dt_new < pt+dur/2))
            model_lc[now] -= depth
    else:
        raise ValueError('Invalid transit model')

    model_oversamp = []
    for i in range(0, len(target.time)):
        mean_f = np.mean(model_lc[i * n_ev:i * n_ev + n_ev - 1])
        model_oversamp.append(mean_f)

    mod_lc = target.rawflux * np.array(model_oversamp)
    if np.isnan(np.mean(mod_lc)):
        print 'failed: ', str(target.epic)
        with open('errlog.txt', 'a') as f:
            print>>f, target.epic
    else:
        if model == 'mandelagol':
            np.savetxt(outdir + target.epic + '_%.2f_%.2f_%.2f_%.2f.txt' % (t0, P, rad, b),
                   np.transpose((target.time, mod_lc)),
                   fmt='%.6f')
        elif model == 'box':
            np.savetxt(outdir + target.epic + '_%.2f_%.2f_%.2f.txt' % (t0, P, rad),
                       np.transpose((target.time, mod_lc)),
                       fmt='%.6f')
    if return_lc:
        return mod_lc

class Target:
    wave, through = np.loadtxt('kepler_response_lowres1.txt', unpack=True, skiprows=9)
    through *= 100.

    def __init__(self, epic, indir):
        self.epic = epic
        try:
            t, f = np.genfromtxt(indir+epic+'.txt', unpack=True, usecols=(0, 1), comments='#')
        except IOError:
            print 'Error: Target does not exist'
            raise

        # numpts = 2*24*30  # 30-day campaign
        self.time = t
        self.rawflux = f
        self.radius = 0.5
        self.Teff = 4000
        self.logg = 4
        self.z = 0
        self.mass = 0.5
        self.u1 = 0.3
        self.u2 = 0.3
        self.get_info()
        self.get_ld()


    def get_info(self):
        url = "https://exofop.ipac.caltech.edu/k2/edit_target.php?id="+self.epic+"#files"
        s = urllib.urlopen(url).read()
        tree = etree.HTML(s)
        tb = tree.xpath('//table[@id="myTable2"]')
        td_content = [[td.text for td in tr.xpath('td')] for tr in tb[0].xpath('//tr')]

        try:
            match = [s for s in range(len(td_content)) if len(td_content[s]) > 0 and 'Stellar Parameters' in
                     td_content[s][0]][0]
            # params are [Teff, logg, radius, vsini, spectral type,
            # log R'HK, [Fe/H], dist, mass, density, Prot ...]

            params = [w.strip() for w in td_content[match+2]]

        except IndexError:
            print 'Stellar parameters not available on EXOFOP'
            raise

        self.radius = float(params[2])
        self.logg = float(params[1])
        if self.logg > 5:
            self.logg = 4.99  # ldtk can't handle logg >= 5
        self.Teff = float(params[0])
        self.z = float(params[6])
        self.mass = float(params[8])

    def get_ld(self):
        print self.Teff, self.logg, self.z
        f = interp1d(self.wave, self.through)
        wave_hires = np.linspace(min(self.wave), max(self.wave), 500, endpoint=True)
        through_hires = f(wave_hires)

        filters = [ldtk.TabulatedFilter('stis', wave_hires, through_hires)]
        # teff, logg, z
        sc = ldtk.LDPSetCreator([self.Teff, 100], [self.logg, 0.05], [self.z, 0.05], filters)
        ps = sc.create_profiles(nsamples=500)
        u1, u2 = ps.coeffs_qd(do_mc=False)
        self.u1 = u1[0][0]
        self.u2 = u2[0][0]

    def inject_transit(self, outdir, multi=True, model='mandelagol'):
        """
        Inject 2000 transits per star, with randomly selected parameters.
        :param outdir: directory to save output.
        :param multi: if True, loops will be executed in parallel.
        :return:
        """

        if multi:
            Parallel(n_jobs=16)(delayed(planet_gen)(self, outdir, model) for _ in range(2000))
        else:
            for _ in range(2000):
                planet_gen(self, outdir, model)


def test(epic, indir='k2/k2mdwarfs/', outdir='k2/box/injected/', model='box'):
    try:
        target = Target(epic, indir)
    except ValueError:
        print 'Invalid stellar parameters'
        return
    flux = planet_gen(target, outdir, model, return_lc=True)
    plt.plot(target.time, flux, 'b.')
    plt.show()


def main(epic, indir='k2/k2mdwarfs/', outdir='k2/injected/', multi=True, model='mandelagol'):
    try:
        target = Target(epic, indir)
    except (ValueError, IndexError):
        print 'Invalid stellar parameters'
        with open('errlog.txt', 'a') as f:
            print>>f, epic
        return
    except TypeError:
        print 'Unable to retrieve limb darkening parameters'
        with open('errlog.txt', 'a') as f:
            print>>f, epic
        return
    except IOError:
        return
    target.inject_transit(outdir=outdir, multi=multi, model=model)


if __name__ == '__main__':
    indir = 'k2/k2mdwarfs/'
    outdir = 'k2/box/injected/'
    stars = np.loadtxt(indir + 'mdwarfs.ls', dtype='S9')

    for epic in stars:
        old = glob.glob(outdir + epic + '*.txt')
        for f in old:
            os.remove(f)
        print 'Working on star', epic
        main(epic, indir, outdir, model='box')


