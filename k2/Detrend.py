#!/usr/bin/env python
# this file is to solve long trend filting by least square method using python's own routine.
import math
import numpy as np
import scipy as sp
import copy
from scipy import linalg
from scipy.optimize import leastsq
from scipy import interpolate
from scipy import signal
from scipy import stats
import ConfigParser
import matplotlib
import matplotlib.pyplot as plt
import logging
from gatspy.periodic import LombScargleFast
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel
from scipy.optimize import minimize

# create logger
module_logger = logging.getLogger('Lightcurve.detrend')

class DetrendFunc(object):

    def __init__(self):
        self.logger = logging.getLogger(DetrendFunc.__name__)
        self.logger.addHandler(logging.NullHandler())
        self.required_cols = ['jd', 'rlc']
        return
    
    def __getstate__(self):
        d = self.__dict__.copy()
        if 'logger' in d.keys():
            d['logger'] = d['logger'].name
        return d

    def __setstate__(self, d):
        if 'logger' in d.keys():
            d['logger'] = logging.getLogger(d['logger'])
        self.__dict__.update(d)

   
    @staticmethod
    def create(method, cfgfile='example.cfg'): 
        for c in [PolyDetrend, COSDetrend, EPDDetrend, FocusDetrend]:
            if c.istype(method):
                return c(cfgfile=cfgfile)

        raise NotImplementedError, 'detrending method %s is not recognized' % method 

    def detrend(self, data):
        ltflc = data['rlc']
        return ltflc


class PolyDetrend(DetrendFunc):
    # Detrend with polynomial fit
    METHOD = 'poly'
    def __init__(self, cfgfile='example.cfg'):
        self.logger = logging.getLogger(PolyDetrend.__name__)
        self.logger.addHandler(logging.NullHandler())
        self.required_cols = ['jd', 'rlc']
        self.required_keys = ['wn', 'norder']
        self.wn = 13
        self.norder = 5
        if not cfgfile == '':
            self.config(cfgfile)

    @staticmethod
    def istype(method):
        return method == PolyDetrend.METHOD

    def config(self, cfgfile):
        p = ConfigParser.RawConfigParser()
        p.read(cfgfile)
        self.logger.debug("%s", str(self.required_keys))
        for name in self.required_keys:
            try:
                # print name
                var = eval(p.get('LC', name))
            except ConfigParser.MissingSectionHeaderError:
                self.logger.error('The format of cfg %s file is not correct, ' \
                      'excute set_parse to see example.cfg, ' \
                      'missing section Header LC\n', cfgfile)
                raise
            except ConfigParser.NoOptionError:
                self.logger.warning('The option of cfg file %s is not correct, ' \
                      'excute set_parse to see example.cfg, ' \
                      'no option %s is given, ' \
                      'use default value' \
                      , cfgfile, name + ' ' + str(getattr(self, name)))
                continue
                # raise
            setattr(self, name, var)
        return

    def detrend(self, data, noplot=True, intran=[]):
        otime = data['jd']
        if not intran:
            intran = otime < 0
        oflux = data['rlc']
        order = self.norder
        length = len(oflux)
        ctime = otime[-intran]
        cflux = sp.signal.medfilt(oflux[-intran] - np.mean(oflux[-intran]), self.wn)
        fitflux, c = lspolyordern(ctime, cflux, order)
        rflux = np.zeros(length)
        for i in range(order + 1):
            rflux += c[i]*otime**(order-i)

        if not noplot:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(otime, oflux, '.')
            plt.show()

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(otime, rflux, '.', ctime, cflux, 'x')
            plt.show()
        dflux = oflux - rflux
        dflux -= np.mean(dflux) - np.mean(oflux)

        return dflux


class GPDetrend(DetrendFunc):
    METHOD = 'GP'
    def __init__(self):
        self.logger = logging.getLogger(GPDetrend.__name__)
        self.logger.addHandler(logging.NullHandler())
        self.required_cols = ['jd', 'ltflc']

    @staticmethod
    def istype(method):
        return method == GPDetrend.METHOD

    @staticmethod
    def get_period(data):
        time = data['jd']
        flux = data['ltflc']
        model = LombScargleFast(fit_offset=True).fit(time, flux)
        _, _ = model.periodogram_auto(oversampling=100)
        model.optimizer.period_range = (0.1, 30)
        print 'Best-fit period:', model.best_period
        return model.best_period

    @staticmethod
    def neg_log_like(params, y, gp):
        gp.kernel[:] = params
        ll = gp.lnlikelihood(y, quiet=True)
        # The scipy optimizer doesn't play well with infinities.
        return -ll if np.isfinite(ll) else 1e25

    def detrend(self, data):
        time = data['jd']
        flux = data['ltflc']
        period = self.get_period(time, flux)

        sigma = np.std(flux[:50])
        sorted = np.sort(flux)
        amp = sorted[-3] - sorted[3]
        a, l, G, P = amp / 2, 5, 0.1, period
        kernel = a * ExpSquaredKernel(l) * ExpSine2Kernel(G, P)
        gp = george.GP(kernel, mean=np.mean(flux))
        gp.compute(time, sigma)
        print 'Initial ln likelihood:', gp.lnlikelihood(flux)

        bounds = ((None, None), (np.log(period), None), (None, None), (np.log(period / 2), None))
        initial_params = gp.kernel.vector
        r = minimize(self.neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(flux, gp, 'george'))
        print r
        gp.kernel[:] = r.x
        pred_mean, var = gp.predict(flux, time)

        sigma = np.std(flux - pred_mean)
        good = np.where(abs(flux - pred_mean) < 3 * sigma)
        gp.compute(time[good], sigma)
        initial_params = gp.kernel.vector
        r = minimize(self.neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(flux[good], gp, 'george'))
        gp.kernel[:] = r.x
        pred_mean, var = gp.predict(flux[good], time)
        print 'Final ln likelihood:', gp.lnlikelihood(flux[good])

        return pred_mean


class COSDetrend(DetrendFunc):
    # Detrend with high pass filter
    METHOD = 'cos'
    def __init__(self, cfgfile='example.cfg'):
        self.logger = logging.getLogger(COSDetrend.__name__)
        self.logger.addHandler(logging.NullHandler())
        self.required_cols = ['jd', 'rlc']
        self.required_keys = ['tmin', 'wn', 'ncomp']
        self.tmin = 1.0
        self.wn = 13
        self.ncomp = 30
        if not cfgfile == '':
            self.config(cfgfile)
        self.logger.info("USE a high pass filter to detrend, it is configured with %s", self.__str__)
    
    def __str__(self):
        parastr=''
        for key in self.required_keys:
            parastr+="%s=%s, " % (key, getattr(self,key))
        return parastr

    @staticmethod
    def istype(method):
        return method == COSDetrend.METHOD

    def config(self, cfgfile):
        p = ConfigParser.RawConfigParser()
        p.read(cfgfile)
        self.logger.debug("%s", str(self.required_keys))
        for name in self.required_keys:
            try:
                # print name
                var = eval(p.get('LC', name))
            except ConfigParser.MissingSectionHeaderError:
                self.logger.error('The format of cfg %s file is not correct, ' \
                      'excute set_parse to see example.cfg, ' \
                      'missing section Header LC\n', cfgfile)
                raise
            except ConfigParser.NoOptionError:
                self.logger.warning('The option of cfg file %s is not correct, ' \
                      'excute set_parse to see example.cfg, ' \
                      'no option %s is given,use default value' \
                      , cfgfile, name + ' ' + str(getattr(self, name)))
                continue
                # raise
            setattr(self, name, var)
        return 

    def detrend(self, data, noplot=True, intran=[]):
        otime = data['jd'] 
        oflux = data['rlc']
        length = len(oflux)
        index = (~np.isnan(oflux))
        ctime = otime[index]
        cflux = oflux[index] - np.nanmean(oflux[index])
        cflux = sp.signal.medfilt(cflux, self.wn)
        
        k = (cflux[-1] - cflux[0]) / (ctime[-1] - ctime[0])
        b = cflux[0] - k * ctime[0]
        cflux -= (k * ctime + b)
        e0 = min(ctime)
        timespan = max(ctime) - min(ctime)
        self.logger.debug("E0=%f, len(cflux)=%d, len(ctime)=%d", e0, len(cflux),len(ctime))
        amatrix = matrixgen(ctime-e0, self.ncomp, timespan)
        # print A.shape,n
        c, resid, rank, sigma = linalg.lstsq(amatrix, cflux)
        n = self.ncomp
        # print resid
        e = np.rollaxis(np.sin(np.array(np.r_['c', 0: n]
                                        * (otime[np.arange(length)]-e0))
                               * math.pi/timespan), 1, 0)
        rflux = np.dot(e, c)
        eflux = k * otime + b
        # print rflux.shape,eflux.shape
        rflux += eflux
        if not noplot:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(otime, rflux, '.', otime, eflux, '+', ctime, cflux, 'x')
            plt.show()
        dflux = oflux - rflux
        dflux -= np.nanmean(dflux)
        dflux += np.nanmean(oflux)
        return dflux


class EPDDetrend(DetrendFunc):
    # Detrend with external parameters
    METHOD = 'epd'
    def __init__(self):
        self.logger = logging.getLogger(EPDDetrend.__name__)
        self.logger.addHandler(logging.NullHandler())
        self.required_cols = ['jd', 'rlc', 'x', 'y', 'bg', 'bgerr']
        return
    
    @staticmethod
    def istype(method):
        return method == EPDDetrend.METHOD


    def detrend(self, data):
        centx = data['x']
        centy = data['y']
        bg = data['bg']
        bgerr = data['bgerr']
        mag = data['ltf']

        def e(c, x, y, bg, bgdev, mag):
            return c[0] + c[1] * np.sin(2 * np.pi * x) \
                   + c[2] * np.cos(2 * np.pi * x) \
                   + c[3] * np.sin(2 * np.pi * y) \
                   + c[4] * np.cos(2 * np.pi * y) \
                   + c[5] * np.sin(4 * np.pi * x) \
                   + c[6] * np.cos(4 * np.pi * x) \
                   + c[7] * np.sin(4 * np.pi * y) \
                   + c[8] * np.cos(4 * np.pi * y) \
                   + c[9] * bg \
                   + c[10] * bgdev \
                   - mag
        # Kepmag = eval(result)
        c0 = np.zeros(11) + 1.
        c1, success = leastsq(e, c0, args=(centx, centy, bg, bgerr, magarr[a, :]))
        if (c1 == 1).all():
            dmag = mag * 1.0
            self.logger.warning("epd failed") 
        else:
            dmag = -e(c1, centx, centy, bg, bgerr, mag)
        dmag = dmag - sp.stats.nanmedian(dmag) + sp.stats.nanmedian(mag)
        return dmag


class FocusDetrend(DetrendFunc):
    # Detrend with Focus series only
    METHOD = 'focus'
    def __init__(self, cfgfile='example.cfg'):
        self.logger = logging.getLogger(FocusDetrend.__name__)
        self.logger.addHandler(logging.NullHandler())
        self.required_cols = ['cadence', 'rlc']
        self.required_keys = ['wn', 'norder']
        self.norder = 3
        self.wn = 13
        if not cfgfile == '':
            self.config(cfgfile)
        return
    
    @staticmethod
    def istype(method):
        return method == FocusDetrend.METHOD



    def config(self, cfgfile):
        p = ConfigParser.RawConfigParser()
        p.read(cfgfile)
        # print self.required_keys
        for name in self.required_keys:
            try:
                # print name
                var = eval(p.get('LC', name))
            except ConfigParser.MissingSectionHeaderError:
                self.logger.error('The format of cfg %s file is not correct, ' \
                      'excute set_parse to see example.cfg, ' \
                      'missing section Header LC\n', cfgfile)
                raise
            except ConfigParser.NoOptionError:
                self.logger.warning('#Warning: the option of cfg file %s is not correct, ' \
                      'excute set_parse to see example.cfg, no option %s is given,' \
                      'use default value', cfgfile, name + ' '
                                             + str(getattr(self, name)))
                continue
                # raise
            setattr(self, name, var)
        return 

    def detrend(self, data):
        xf = ((data['cadence'] / 48. / 13.7 + 0.5) % 1) - 0.5
        mag = data['rlc']
        focus = 0.0 + 10.0 * np.exp(-0.5 * (xf / 0.1)**2.)
        dmag = np.zeros(len(mag)) + mag
        noflux = sp.signal.medfilt(mag, self.wn)
        n = self.norder
        seg1 = data['cadence'] < 624
        x = data['cadence'][seg1]
        xi = sp.zeros(len(x)) + 1
        xii = sp.zeros(len(x)) + focus[seg1]
        amatrix = (x**n)[:, np.newaxis]
        for i in xrange(1, n):
            amatrix = np.hstack((amatrix, (x**(n-i))[:, np.newaxis]))
        amatrix = np.c_[amatrix, xi[:, np.newaxis]]
        amatrix = np.c_[amatrix, xii[:, np.newaxis]]
        c, resid, rank, sigma = sp.linalg.lstsq(amatrix, noflux[seg1])
        z = sp.zeros(len(x))
        for i in xrange(n+1):
            z += c[i] * x**(n-i)
        z += c[-1] * xii
        dmagseg1 = mag[seg1] - z + np.median(mag)
        
        seg2 = data['cadence'] > 624
        x = data['cadence'][seg2]
        xi = sp.zeros(len(x)) + 1
        xii = sp.zeros(len(x)) + focus[seg2]
        amatrix = (x**n)[:, np.newaxis]
        for i in xrange(1, n):
            amatrix = np.hstack((amatrix, (x**(n-i))[:, np.newaxis]))
        amatrix = np.c_[amatrix, xi[:, np.newaxis]]
        amatrix = np.c_[amatrix, xii[:, np.newaxis]]
        c, resid, rank, sigma = sp.linalg.lstsq(amatrix, noflux[seg2])
        z = sp.zeros(len(x))
        for i in xrange(n+1):
            z += c[i] * x**(n-i)
        z += c[-1] * xii
        dmagseg2 = mag[seg2] - z + np.median(mag)
        dmag[seg1] = dmagseg1
        dmag[seg2] = dmagseg2
        return dmag


def lspolyordern(x, y, n):
    x = sp.array(x)
    y = sp.array(y)
    xi = sp.zeros(len(x)) + 1
    amatrix = (x**n)[:, np.newaxis]
    for i in range(1, n):
        amatrix = np.hstack((amatrix, (x**(n-i))[:, np.newaxis]))

    amatrix = np.c_[amatrix, xi[:, np.newaxis]]
    c, resid, rank, sigma = sp.linalg.lstsq(amatrix, y)
    z = sp.zeros(len(x))
    for i in range(n+1):
        z += c[i]*x**(n-i)
    return [z, c]


def matrixgen(time, n, timespan):
    # generate least square fitting matrix with n cos filter,
    # t is the total time span, formulism refer to Huang and Bacos (2012) eq [1]'''
    tn = len(time)
    a = np.rollaxis(np.sin(np.array(np.r_['c', 0:n] * time[np.arange(tn)])
                           * math.pi / timespan), 1, 0)
    return a


def flatfunc(c, ctime, timespan):
    n = len(c) - 1
    b = np.rollaxis(np.sin(np.array(np.r_['c', 0:n] * ctime) * math.pi / timespan),
                    1, 0)
    rflux = np.dot(b, c[0:n])
    rflux -= np.mean(rflux)
    return rflux


def tranfunc(c, ctime, intran, timespan):
    n = len(c)-1
    b = np.rollaxis(np.sin(np.array(np.r_['c', 0: n] * ctime) * math.pi / timespan),
                    1, 0)
    rflux = np.dot(b, c[0: n])
    try:
        rflux[intran] += c[n]
    except TypeError:
        print c
        raise
    rflux -= np.mean(rflux)
    return rflux


def lsfitrecon(otime, oflux, intran, n=30, noplot=True, dipguess=0.0001):
    length = len(oflux)
    ctime = np.array(otime)
    epoch0 = min(ctime)
    timespan = max(ctime)-min(ctime)
    cflux = np.array(oflux)-np.mean(oflux)

    def e(c, time, index, flux, timespan):
        return tranfunc(c, time, index, timespan) - flux

    c0 = list(np.zeros(n))
    c0.append(dipguess)
    c, success = leastsq(e, c0, args=(ctime - epoch0, intran, cflux, timespan), maxfev=10000)
    b = np.rollaxis(np.sin(np.array(np.r_['c', 0:n]
                                    * (ctime[np.arange(length)]-E0))
                           * math.pi/timespan), 1, 0)
    rflux = np.dot(b, c[0:n])
    tran = sp.zeros(len(ctime))
    tran[intran] += c[n]
    tran -= np.mean(tran)
    # print c[n],len(tran[intran])
    if not noplot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(ctime, cflux, 'x', ctime, rflux, '.', ctime, tran, '+')
        plt.show()
        plt.close(fig)
    # rl = len(rflux)
    # nn = len(otime)
    dflux = oflux - rflux + np.mean(oflux)
    return dflux
